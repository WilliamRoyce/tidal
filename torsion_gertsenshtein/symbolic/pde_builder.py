"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid
from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from numpy.typing import NDArray
    from pde.grids.base import GridBase
    from pde.pdes.base import TState

    from torsion_gertsenshtein.utils import BCDescriptor

    NumericArray = NDArray[np.float64]


# ---------------------------------------------------------------------------
# Operator registry: maps operator name -> (handler, min_grid_dimension)
# Each handler takes (field: ScalarField, bc: BCDescriptor) -> ScalarField
# ---------------------------------------------------------------------------


def _op_laplacian(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    return field.laplace(bc=bc)


def _op_identity(field: ScalarField, _bc: BCDescriptor) -> ScalarField:
    return field.copy()


def _op_gradient(axis: int) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a gradient handler for a specific axis."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)
        component = grad[axis]
        assert isinstance(component, ScalarField)
        return component

    return _handler


def _op_directional_laplacian(
    axis: int,
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a directional Laplacian handler (∂²/∂x_i²)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)[axis]
        assert isinstance(grad, ScalarField)
        d2 = grad.gradient(bc=bc)[axis]
        assert isinstance(d2, ScalarField)
        return d2

    return _handler


def _op_cross_derivative(
    axis1: int, axis2: int
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a cross derivative handler (∂²/∂x_i ∂x_j)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad_j = field.gradient(bc=bc)[axis2]
        assert isinstance(grad_j, ScalarField)
        grad_ij = grad_j.gradient(bc=bc)[axis1]
        assert isinstance(grad_ij, ScalarField)
        return grad_ij

    return _handler


def _op_biharmonic(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    """Biharmonic operator: ∇⁴f = ∇²(∇²f)."""
    lap = field.laplace(bc=bc)
    assert isinstance(lap, ScalarField)
    bilap = lap.laplace(bc=bc)
    assert isinstance(bilap, ScalarField)
    return bilap


def _op_nth_derivative(
    axis: int, order: int
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a handler for Nth-order derivative along a single axis.

    Applies gradient in the given axis direction ``order`` times.
    E.g., order=3, axis=0 computes ∂³f/∂x³.
    """

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        result: ScalarField = field
        for _ in range(order):
            grad = result.gradient(bc=bc)
            component = grad[axis]
            assert isinstance(component, ScalarField)
            result = component
        return result

    return _handler


#: Map axis letter to axis index.
_AXIS_INDEX: dict[str, int] = {"x": 0, "y": 1, "z": 2}

#: Minimum grid dimension required for each axis.
_AXIS_MIN_DIM: dict[str, int] = {"x": 1, "y": 2, "z": 3}

#: Regex for parsing generic single-axis derivative names.
_GENERIC_SINGLE_RE = re.compile(r"^derivative_(\d+)_([xyz])$")


#: Registry mapping operator names to (handler, min_dimension) pairs.
#: To add a new operator, simply add an entry here.
_OPERATOR_REGISTRY: dict[str, tuple[Any, int]] = {
    "laplacian": (_op_laplacian, 1),
    "identity": (_op_identity, 1),
    "gradient_x": (_op_gradient(0), 1),
    "gradient_y": (_op_gradient(1), 2),
    "gradient_z": (_op_gradient(2), 3),
    "laplacian_x": (_op_directional_laplacian(0), 1),
    "laplacian_y": (_op_directional_laplacian(1), 2),
    "laplacian_z": (_op_directional_laplacian(2), 3),
    "cross_derivative_xy": (_op_cross_derivative(0, 1), 2),
    "cross_derivative_xz": (_op_cross_derivative(0, 2), 3),
    "cross_derivative_yz": (_op_cross_derivative(1, 2), 3),
    "biharmonic": (_op_biharmonic, 1),
    # Note: first_derivative_t is handled specially in _compute_rhs_for_component
}


@dataclass(frozen=True)
class ParsedFieldName:
    """Parsed field name components.

    Supports multiple field naming conventions:
    - standard: A_0, phi_1 (base_index)
    - tensor: stress_xy_0, u_x_1 (base_component_index)
    - compact: phi0, A1 (base+digits)
    - simple: phi, psi (letters only, index defaults to 0)
    """

    base: str
    index: int
    format: str

    @classmethod
    def parse(cls, name: str) -> ParsedFieldName:
        """Parse field name auto-detecting format.

        Parameters
        ----------
        name : str
            Field name to parse.

        Returns
        -------
        ParsedFieldName
            Parsed components with base, index, and format.
        """
        # Standard format: A_0, phi_1
        match = re.match(r"^([a-zA-Z]+)_([0-9]+)$", name)
        if match:
            return cls(
                base=match.group(1), index=int(match.group(2)), format="standard"
            )

        # Tensor format: stress_xy_0, u_x_1 (greedy match for base)
        match = re.match(r"^(.+)_([0-9]+)$", name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="tensor")

        # Compact format: phi0, A1
        match = re.match(r"^([a-zA-Z]+)([0-9]+)$", name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="compact")

        # Simple format: phi, psi (no index, defaults to 0)
        if re.match(r"^[a-zA-Z]+$", name):
            return cls(base=name, index=0, format="simple")

        # Fallback
        return cls(base=name, index=0, format="unknown")

    def to_momentum_name(self) -> str:
        """Convert to momentum field name."""
        return f"pi_{self.index}"


def parse_momentum_field_name(field_name: str) -> int | None:
    """Parse momentum field name and return index.

    Supports both pi_N and piN formats.

    Parameters
    ----------
    field_name : str
        Momentum field name like "pi_0", "pi0", "pi_1", "pi1".

    Returns
    -------
    int | None
        Index if valid momentum field name, None otherwise.
    """
    # Standard format: pi_N
    match = re.match(r"^pi_([0-9]+)$", field_name)
    if match:
        return int(match.group(1))

    # Compact format: piN
    match = re.match(r"^pi([0-9]+)$", field_name)
    if match:
        return int(match.group(1))

    return None


class PDEFromSpec(PDEBase):
    """Generic PDE class built from JSON equation specification.

    This class dynamically constructs the evolution equations from a parsed
    specification. NO physics is hardcoded - all equation structure comes
    from the EquationSystem that was derived from a Lagrangian.

    The state layout is determined by each component's time derivative order:
    - Second-order (wave, time_order>=2): [field_i, momentum_i] pair
    - First-order (heat/diffusion, time_order=1): [field_i] only
    - Constraint (elliptic, time_order=0): [field_i] only (no evolution)

    For second-order components:
        d/dt field_i = momentum_i
        d/dt momentum_i = RHS from specification
    For first-order components:
        d/dt field_i = RHS from specification
    For constraint components:
        d/dt field_i = 0

    Parameters
    ----------
    spec : EquationSystem
        The equation specification loaded from JSON.

    Attributes
    ----------
    spec : EquationSystem
        The equation specification.
    n_components : int
        Number of field components.
    explicit_time_dependence : bool
        True to support time-dependent coefficients (e.g., de Sitter spacetime).

    Examples
    --------
    >>> from torsion_gertsenshtein.symbolic import load_equation_system
    >>> from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec
    >>> spec = load_equation_system("examples/data/em_1d.json")
    >>> pde = PDEFromSpec(spec)
    >>> # pde can now be used with py-pde solvers
    """

    # Enable time-dependent coefficients for curved spacetime (e.g., de Sitter)
    # This allows evolution_rate to receive the current time t
    explicit_time_dependence = True

    def __init__(
        self,
        spec: EquationSystem,
        parameters: dict[str, float] | None = None,
    ) -> None:
        """Initialize PDE from equation specification.

        Parameters
        ----------
        spec : EquationSystem
            The equation specification loaded from JSON.
        parameters : dict[str, float] | None
            Optional parameter values to override symbolic coefficients.
            Keys are symbolic names (e.g., "dSH", "dSm2", "kappa"), values are numeric.
            When a term has a coefficient_symbolic that matches a key in this
            dict, the parameter value is used instead of the numeric coefficient.

            For time-dependent coefficients in curved spacetime, all symbols
            appearing in the coefficient expression must be provided here.
            The expressions are evaluated by substituting these values.
        """
        super().__init__()
        self.spec = spec
        self.n_components = spec.n_components
        self._component_name_to_index = {
            name: i for i, name in enumerate(spec.component_names)
        }
        self._parameters = parameters or {}

        # Build slot maps from state_layout for mixed time-order support
        self._field_slot_map: dict[str, int] = {}
        self._momentum_slot_map: dict[str, int] = {}
        for slot_idx, (name, slot_type) in enumerate(spec.state_layout):
            if slot_type == "field":
                self._field_slot_map[name] = slot_idx
            else:
                self._momentum_slot_map[name] = slot_idx

    def _resolve_coefficient(self, term: OperatorTerm) -> float:
        """Resolve the effective coefficient for a term.

        If the term has a symbolic coefficient name and that name (or its
        negation) is in the parameters dict, use the parameter value.
        Otherwise use the numeric coefficient from the JSON.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.

        Returns
        -------
        float
            The effective coefficient value.
        """
        if term.coefficient_symbolic is not None:
            sym = term.coefficient_symbolic

            # Check for negated symbol like "-m2"
            if sym.startswith("-") and sym[1:] in self._parameters:
                return -self._parameters[sym[1:]]
            if sym in self._parameters:
                return self._parameters[sym]

            # Symbolic coefficient present but not resolvable from parameters

            warnings.warn(
                f"Symbolic coefficient '{sym[:80]}' could not be resolved from "
                f"parameters {sorted(self._parameters.keys())}. "
                f"Falling back to numeric coefficient {term.coefficient}. "
                f"This may indicate unevaluated expressions in the JSON.",
                stacklevel=2,
            )

        # Default: use numeric coefficient from JSON
        return term.coefficient

    def _mathematica_to_python(self, expr: str) -> str:
        """Convert Mathematica InputForm expression to evaluable Python.

        Handles common Mathematica syntax:
        - ``E^(...)`` to ``exp(...)`` (Euler's number)
        - ``Sin[x]`` to ``sin(x)``, ``Cos[x]`` to ``cos(x)``, etc.
        - ``t[]`` to ``t`` (xCoba coordinate symbols, using actual coordinate names)
        - Mathematica brackets ``[``, ``]`` to Python parens ``(``, ``)``
        - Multiplication, grouping, negation are already Python-compatible

        Parameters
        ----------
        expr : str
            Mathematica InputForm expression string.

        Returns
        -------
        str
            Python-evaluable expression string.
        """
        result = expr
        # E^(...) -> exp(...) — Mathematica's Euler number raised to a power
        result = re.sub(r"\bE\^", "exp", result)
        # Common math functions: Sin[x] -> sin(x), etc.
        for mma_func, py_func in [
            ("Sin", "sin"),
            ("Cos", "cos"),
            ("Tan", "tan"),
            ("Log", "log"),
            ("Sqrt", "sqrt"),
            ("Abs", "abs"),
        ]:
            result = re.sub(rf"\b{mma_func}\b", py_func, result)
        # Mathematica brackets to Python parens (after function renaming)
        result = result.replace("[", "(").replace("]", ")")
        # Mathematica ^ to Python ** (AFTER E^ -> exp to avoid double-conversion)
        result = result.replace("^", "**")
        # xCoba coordinate symbols appear as zero-arg function calls: t() -> t, x() -> x, etc.
        # Uses actual coordinate names from the equation system (not hardcoded x/y/z).
        for coord in self.spec.effective_coordinates:
            result = result.replace(f"{coord}()", coord)
        return result

    def _resolve_coefficient_at_point(
        self,
        term: OperatorTerm,
        t: float,
        grid: GridBase | None = None,
    ) -> float | NumericArray:
        """Resolve a potentially coordinate-dependent coefficient.

        Handles constant, time-dependent, and position-dependent coefficients
        in a single unified code path. Uses numpy functions throughout so that
        the same evaluation works for both scalar (time-only) and array
        (position-dependent) results.

        Returns a scalar ``float`` for constant or time-only coefficients, or a
        ``numpy.ndarray`` (same shape as the grid) for position-dependent ones.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.
        t : float
            Current simulation time.
        grid : GridBase | None
            Simulation grid. Required when the term is position-dependent.

        Returns
        -------
        float | NumericArray
            Scalar or grid-shaped array of coefficient values.

        Raises
        ------
        ValueError
            If required parameters/grid are missing or expression cannot be evaluated.
        """
        # Fast path: no coordinate dependence → use simple parameter lookup
        if not term.time_dependent and not term.position_dependent:
            return self._resolve_coefficient(term)

        # Position-dependent: need grid for spatial coordinates
        if term.position_dependent and grid is None:
            msg = (
                f"Position-dependent coefficient '{term.coefficient_symbolic}' "
                f"requires grid info but no grid was provided."
            )
            raise ValueError(msg)

        sym = term.coefficient_symbolic or ""
        py_expr = self._mathematica_to_python(sym)

        # Build evaluation namespace with numpy functions (work on scalars AND arrays)
        namespace: dict[str, Any] = dict(self._parameters)
        namespace["t"] = t
        namespace["exp"] = np.exp
        namespace["sin"] = np.sin
        namespace["cos"] = np.cos
        namespace["tan"] = np.tan
        namespace["log"] = np.log
        namespace["sqrt"] = np.sqrt
        namespace["abs"] = np.abs

        # Inject spatial coordinates from grid when position-dependent
        if term.position_dependent:
            spatial_coords = self.spec.spatial_coordinates
            coords = grid.cell_coords  # type: ignore[union-attr]
            for i, name in enumerate(spatial_coords[: grid.num_axes]):  # type: ignore[union-attr]
                namespace[name] = np.asarray(coords[..., i])  # pyright: ignore[reportUnknownArgumentType]

        # Validate all symbols can be resolved
        identifiers = set(re.findall(r"\b[a-zA-Z_]\w*\b", py_expr))
        builtin_names = {"exp", "sin", "cos", "tan", "log", "sqrt", "abs"}
        coord_vars = set(self.spec.effective_coordinates)
        identifiers -= builtin_names | coord_vars
        missing = identifiers - set(self._parameters.keys())
        if missing:
            msg = (
                f"Parameters {sorted(missing)} are required for "
                f"coordinate-dependent coefficient '{sym}'. "
                f"Pass them via parameters={{...}} to PDEFromSpec or build_pde_from_json."
            )
            raise ValueError(msg)

        try:
            result = eval(py_expr, {"__builtins__": {}}, namespace)  # noqa: S307
            if isinstance(result, np.ndarray):
                return np.asarray(result, dtype=np.float64)
            return float(result)
        except Exception as e:
            msg = (
                f"Cannot evaluate coordinate-dependent coefficient '{sym}' "
                f"(Python form: '{py_expr}') at t={t}: {e}"
            )
            raise ValueError(msg) from e

    @staticmethod
    def _get_operator(
        operator_name: str, field: ScalarField, bc: BCDescriptor
    ) -> ScalarField:
        """Apply a named operator to a field.

        Uses the module-level ``_OPERATOR_REGISTRY`` for dispatch.
        Each operator specifies a handler function and minimum grid dimension.

        Parameters
        ----------
        operator_name : str
            Name of the operator ("laplacian", "identity", "gradient_x", etc.)
        field : ScalarField
            The field to operate on.
        bc : BCDescriptor
            Boundary condition specification.

        Returns
        -------
        ScalarField
            Result of applying the operator.

        Raises
        ------
        ValueError
            If the operator is not recognized or the grid dimension is too low.
        """
        entry = _OPERATOR_REGISTRY.get(operator_name)
        if entry is None:
            # Try dynamic resolution for generic Nth-order derivatives
            m = _GENERIC_SINGLE_RE.match(operator_name)
            if m:
                order = int(m.group(1))
                axis_letter = m.group(2)
                axis = _AXIS_INDEX[axis_letter]
                min_dim = _AXIS_MIN_DIM[axis_letter]
                entry = (_op_nth_derivative(axis, order), min_dim)
            else:
                msg = (
                    f"Unknown operator: '{operator_name}'. "
                    f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}. "
                    f"Dynamic patterns: derivative_N_x (N=integer, x/y/z=axis)."
                )
                raise ValueError(msg)

        handler, min_dim = entry
        if field.grid.dim < min_dim:
            msg = (
                f"Operator '{operator_name}' requires at least {min_dim}D grid, "
                f"but got {field.grid.dim}D grid."
            )
            raise ValueError(msg)

        return handler(field, bc)

    def _get_field_from_state(
        self,
        state: FieldCollection,
        field_name: str,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> ScalarField:
        """Get a field from state by name, supporting both field and momentum names.

        For mixed time-space derivatives like d_t d_x A, the Wolfram pipeline
        expresses these as gradients of momentum fields (e.g., gradient_x(pi_0)).
        This is valid because d_t A = pi, so d_x(d_t A) = d_x(pi).

        State layout is determined by ``spec.state_layout``:
        - Second-order components have [field, momentum] pairs
        - First-order components have [field] only (no momentum in state)

        For first-order components, ``pi_N`` references are resolved via the
        ``virtual_momenta`` dict, which contains the computed RHS of those
        first-order equations (since d_t(field) = RHS for first-order PDEs).

        Parameters
        ----------
        state : FieldCollection
            Current state with fields ordered by state_layout.
        field_name : str
            Name of the field to retrieve. Can be:
            - Regular field name like "A_0", "phi_0", "phi0"
            - Momentum field name like "pi_0", "pi_1", "pi0", "pi1"
        virtual_momenta : dict[str, ScalarField] | None
            Pre-computed RHS fields for first-order components, keyed by
            component name. Used when ``pi_N`` references a first-order
            component that has no momentum state variable.

        Returns
        -------
        ScalarField
            The requested field from the state.

        Raises
        ------
        ValueError
            If the field name is not recognized or momentum is unavailable.
        """
        # Check if this is a momentum field reference (pi_0, pi0, etc.)
        if field_name.startswith("pi"):
            momentum_idx = parse_momentum_field_name(field_name)
            if momentum_idx is not None:
                if not (0 <= momentum_idx < self.n_components):
                    msg = (
                        f"Momentum field index {momentum_idx} out of range. "
                        f"This system has {self.n_components} components "
                        f"(valid indices: 0 to {self.n_components - 1}). "
                        f"Field reference: '{field_name}'."
                    )
                    raise ValueError(msg)

                # Look up the component name for this index, then find its momentum slot
                comp_name = self.spec.component_names[momentum_idx]
                slot = self._momentum_slot_map.get(comp_name)
                if slot is not None:
                    # Second-order component: momentum is a state variable
                    momentum = state[slot]
                    assert isinstance(momentum, ScalarField)
                    return momentum

                # First-order component: check virtual_momenta
                if virtual_momenta is not None and comp_name in virtual_momenta:
                    return virtual_momenta[comp_name]

                eq = self.spec.equations[momentum_idx]
                msg = (
                    f"Momentum field '{field_name}' referenced but component "
                    f"'{comp_name}' has time_derivative_order={eq.time_derivative_order} "
                    f"(no momentum in state, and no virtual momentum computed). "
                    f"This may indicate a circular dependency between first-order components."
                )
                raise ValueError(msg)
            # If it looks like momentum but couldn't be parsed, raise clear error
            msg = (
                f"Invalid momentum field format: '{field_name}'. "
                f"Expected 'pi_N' or 'piN' where N is a numeric index (e.g., 'pi_0', 'pi0')."
            )
            raise ValueError(msg)

        # Regular field lookup via slot map
        slot = self._field_slot_map.get(field_name)
        if slot is not None:
            field = state[slot]
            assert isinstance(field, ScalarField)
            return field

        msg = f"Unknown field name: {field_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float = 0.0,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> ScalarField:
        """Compute the RHS for a single component's equation.

        This method evaluates all terms in the component's equation specification
        and sums them together.

        Parameters
        ----------
        component_idx : int
            Index of the component.
        state : FieldCollection
            Current state (all fields and momenta).
        bc : BCDescriptor
            Boundary condition specification.
        t : float
            Current simulation time (for time-dependent coefficients in curved spacetime).
        virtual_momenta : dict[str, ScalarField] | None
            Pre-computed RHS fields for first-order components. When a term
            references ``first_derivative_t`` or ``pi_N`` for a first-order
            component, the virtual momentum is used instead of a state variable.

        Returns
        -------
        ScalarField
            The computed RHS.

        Raises
        ------
        ValueError
            If a field or operator cannot be resolved.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term in eq.rhs_terms:
            target_field_name = term.field

            # Handle first_derivative_t operator specially
            # d_t(field) = momentum (second-order) or RHS (first-order via virtual_momenta)
            if term.operator == "first_derivative_t":
                target_idx = self._component_name_to_index.get(target_field_name)
                if target_idx is None:
                    msg = f"Unknown field for first_derivative_t: {target_field_name}"
                    raise ValueError(msg)

                comp_name = self.spec.component_names[target_idx]

                # Try state momentum slot first (second-order fields)
                slot = self._momentum_slot_map.get(comp_name)
                if slot is not None:
                    momentum = state[slot]
                    assert isinstance(momentum, ScalarField)
                    operated = momentum.copy()
                elif virtual_momenta is not None and comp_name in virtual_momenta:
                    operated = virtual_momenta[comp_name].copy()
                else:
                    eq_target = self.spec.equations[target_idx]
                    msg = (
                        f"first_derivative_t references '{target_field_name}' "
                        f"but it has time_derivative_order={eq_target.time_derivative_order} "
                        f"and no virtual momentum was computed."
                    )
                    raise ValueError(msg)
            else:
                # Standard operator handling (passes virtual_momenta for pi_N resolution)
                target_field = self._get_field_from_state(
                    state, target_field_name, virtual_momenta
                )
                operated = self._get_operator(term.operator, target_field, bc)

            # Resolve coefficient (supports time- and position-dependent)
            coefficient = self._resolve_coefficient_at_point(term, t, grid)

            # Add coefficient * operated to result
            if isinstance(coefficient, np.ndarray):
                contribution = ScalarField(grid, data=coefficient * operated.data)
            else:
                contribution = coefficient * operated
            result += contribution

        return result

    @override
    def evolution_rate(
        self,
        state: TState,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the time derivatives for all fields.

        Supports mixed time-derivative orders:
        - Second-order (wave): d/dt field = momentum, d/dt momentum = RHS
        - First-order (heat/diffusion): d/dt field = RHS
        - Constraint (elliptic, order=0): d/dt field = 0

        For first-order components, the computed RHS is stored as a "virtual
        momentum" so that second-order equations referencing ``pi_N`` or
        ``first_derivative_t`` of a first-order component can access it.

        Parameters
        ----------
        state : FieldCollection
            Current state with ``spec.state_size`` fields.
        t : float
            Current time. Used for time-dependent coefficients in curved spacetime
            (e.g., Hubble friction in de Sitter expansion).

        Returns
        -------
        FieldCollection
            Time derivatives for all fields.

        Raises
        ------
        ValueError
            If the state size or grid dimension does not match the spec.
        """
        assert isinstance(state, FieldCollection)
        expected_fields = self.spec.state_size
        if len(state) != expected_fields:
            msg = f"Expected {expected_fields} fields, got {len(state)}"
            raise ValueError(msg)

        # Validate grid dimension matches spec
        grid_dim = state.grid.dim
        expected_dim = self.spec.spatial_dimension
        if grid_dim != expected_dim:
            msg = (
                f"Grid dimension {grid_dim} does not match spec "
                f"spatial_dimension {expected_dim}. "
                f"The equation system expects a {expected_dim}D spatial grid "
                f"(from {self.spec.dimension}D spacetime)."
            )
            raise ValueError(msg)

        bc = infer_bc_from_grid(state.grid)
        grid = state.grid

        # Pass 1: Compute RHS for first-order components.
        # These become "virtual momenta" that second-order equations can reference
        # via pi_N or first_derivative_t operators.
        virtual_momenta: dict[str, ScalarField] = {}
        for i, eq in enumerate(self.spec.equations):
            if eq.time_derivative_order == 1:
                virtual_momenta[eq.field_name] = self._compute_rhs_for_component(
                    i, state, bc, t
                )

        # Pass 2: Build the full rates array using slot maps
        rates: list[ScalarField | None] = [None] * expected_fields

        for i, eq in enumerate(self.spec.equations):
            field_slot = self._field_slot_map[eq.field_name]

            if eq.time_derivative_order >= 2:  # noqa: PLR2004
                # Second-order: d/dt field = momentum, d/dt momentum = RHS
                momentum_slot = self._momentum_slot_map[eq.field_name]
                momentum = state[momentum_slot]
                assert isinstance(momentum, ScalarField)
                rates[field_slot] = momentum.copy()
                rates[momentum_slot] = self._compute_rhs_for_component(
                    i, state, bc, t, virtual_momenta
                )
            elif eq.time_derivative_order == 1:
                # First-order: d/dt field = RHS (already computed in Pass 1)
                rates[field_slot] = virtual_momenta[eq.field_name]
            else:
                # Constraint (order=0): no time evolution
                rates[field_slot] = ScalarField(grid, data=0.0)

        return FieldCollection(rates)  # type: ignore[arg-type]

    def _cache_key(self) -> dict[str, Any]:
        """Return a cache key for this PDE.

        The key includes the specification metadata to ensure different
        equation systems don't share cached operators.
        """
        return {
            "n_components": self.n_components,
            "component_names": self.spec.component_names,
            "metadata_hash": hash(frozenset(self.spec.metadata.items())),
        }


def build_pde_from_json(
    json_path: Path | str,
    parameters: dict[str, float] | None = None,
) -> PDEFromSpec:
    """Build a PDE from a JSON equation specification file.

    This is the main entry point for the Lagrangian-to-PDE pipeline on the
    Python side. Given a JSON file exported from Mathematica/xAct, this
    function creates a py-pde compatible PDE class.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file containing the equation specification.
    parameters : dict[str, float] | None
        Optional parameter values to override symbolic coefficients.
        Keys are symbolic names (e.g., "m2", "kappa"), values are numeric.
        Example: {"m2": 0.5, "kappa": 1.0}

    Returns
    -------
    PDEFromSpec
        A PDE instance ready for use with py-pde solvers.

    Examples
    --------
    >>> pde = build_pde_from_json("examples/data/em_1d.json")
    >>> # Create initial state and run simulation
    >>> from pde import CartesianGrid, ScalarField, FieldCollection
    >>> grid = CartesianGrid([(0, 100)], 256, periodic=True)
    >>> # ... create initial conditions and solve

    >>> # With custom parameter values:
    >>> pde = build_pde_from_json("examples/data/proca_1d.json", parameters={"m2": 2.0})
    """
    spec = load_equation_system(json_path)
    return PDEFromSpec(spec, parameters=parameters)


def create_initial_state(
    grid: GridBase,
    spec: EquationSystem,
    field_data: dict[str, NDArray[np.float64]] | None = None,
    momentum_data: dict[str, NDArray[np.float64]] | None = None,
) -> FieldCollection:
    """Create initial state for a PDEFromSpec simulation.

    The state layout is determined by ``spec.state_layout``:
    - Second-order components get [field, momentum] pairs
    - First-order/constraint components get [field] only

    Parameters
    ----------
    grid : GridBase
        The simulation grid.
    spec : EquationSystem
        The equation specification.
    field_data : dict[str, NDArray] | None
        Initial data for each field component. Keys are component names.
        Components not specified default to zero.
    momentum_data : dict[str, NDArray] | None
        Initial data for each momentum component. Keys are component names.
        Components not specified default to zero. Ignored for first-order
        and constraint components (which have no momentum slot).

    Returns
    -------
    FieldCollection
        Initial state with ``spec.state_size`` fields.
    """
    field_data = field_data or {}
    momentum_data = momentum_data or {}

    fields: list[ScalarField] = []

    for name, slot_type in spec.state_layout:
        if slot_type == "field":
            if name in field_data:
                fields.append(ScalarField(grid, data=field_data[name]))
            else:
                fields.append(ScalarField(grid, data=0.0))
        # Momentum slot (only for second-order components)
        elif name in momentum_data:
            fields.append(ScalarField(grid, data=momentum_data[name]))
        else:
            fields.append(ScalarField(grid, data=0.0))

    return FieldCollection(fields)
