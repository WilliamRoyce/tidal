"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid
from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)

if TYPE_CHECKING:
    from pathlib import Path

    import numpy as np
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


def _op_identity(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    return field.copy()


def _op_gradient(axis: int):
    """Create a gradient handler for a specific axis."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)
        component = grad[axis]
        assert isinstance(component, ScalarField)
        return component

    return _handler


def _op_directional_laplacian(axis: int):
    """Create a directional Laplacian handler (∂²/∂x_i²)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)[axis]
        assert isinstance(grad, ScalarField)
        d2 = grad.gradient(bc=bc)[axis]
        assert isinstance(d2, ScalarField)
        return d2

    return _handler


def _op_cross_derivative(axis1: int, axis2: int):
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

    The state is organized as a FieldCollection with 2 * n_components fields:
        [field_0, momentum_0, field_1, momentum_1, ...]

    For each component i, the evolution equations are:
        d/dt field_i = momentum_i
        d/dt momentum_i = sum of RHS terms from specification

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

        # Default: use numeric coefficient from JSON
        return term.coefficient

    @staticmethod
    def _mathematica_to_python(expr: str) -> str:
        """Convert Mathematica InputForm expression to evaluable Python.

        Handles common Mathematica syntax:
        - ``E^(...)`` to ``exp(...)`` (Euler's number)
        - ``Sin[x]`` to ``sin(x)``, ``Cos[x]`` to ``cos(x)``, etc.
        - ``t[]`` to ``t`` (xCoba coordinate symbols)
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
        # t() -> t — xCoba coordinate symbols appear as zero-arg function calls
        return result.replace("t()", "t")

    def _resolve_coefficient_at_time(self, term: OperatorTerm, t: float) -> float:
        """Resolve a potentially time-dependent coefficient at time t.

        Converts the Mathematica symbolic coefficient to a Python expression,
        substitutes parameter values and current time, then evaluates.
        This is general and works for any symbolic coefficient expression.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.
        t : float
            Current simulation time.

        Returns
        -------
        float
            The effective coefficient value at time t.

        Raises
        ------
        ValueError
            If required parameters are missing or expression cannot be evaluated.
        """
        # For non-time-dependent terms, use standard resolution
        if not term.time_dependent:
            return self._resolve_coefficient(term)

        sym = term.coefficient_symbolic or ""
        py_expr = self._mathematica_to_python(sym)

        # Build evaluation namespace: parameters + time + math functions
        namespace: dict[str, Any] = dict(self._parameters)
        namespace["t"] = t
        namespace["exp"] = math.exp
        namespace["sin"] = math.sin
        namespace["cos"] = math.cos
        namespace["tan"] = math.tan
        namespace["log"] = math.log
        namespace["sqrt"] = math.sqrt
        namespace["abs"] = abs

        # Check that all symbols in the expression can be resolved
        identifiers = set(re.findall(r"\b[a-zA-Z_]\w*\b", py_expr))
        identifiers -= {"exp", "sin", "cos", "tan", "log", "sqrt", "abs", "t"}
        missing = identifiers - set(self._parameters.keys())
        if missing:
            msg = (
                f"Parameters {sorted(missing)} are required for time-dependent "
                f"coefficient '{sym}'. Pass them via "
                f"parameters={{{', '.join(repr(k) + ': value' for k in sorted(missing))}}} "
                f"to PDEFromSpec or build_pde_from_json."
            )
            raise ValueError(msg)

        try:
            result = eval(py_expr, {"__builtins__": {}}, namespace)  # noqa: S307
            return float(result)
        except Exception as e:
            msg = (
                f"Cannot evaluate time-dependent coefficient '{sym}' "
                f"(Python form: '{py_expr}') with parameters {self._parameters} "
                f"at t={t}: {e}"
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
            msg = (
                f"Unknown operator: '{operator_name}'. "
                f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}"
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
        self, state: FieldCollection, field_name: str
    ) -> ScalarField:
        """Get a field from state by name, supporting both field and momentum names.

        For mixed time-space derivatives like d_t d_x A, the Wolfram pipeline
        expresses these as gradients of momentum fields (e.g., gradient_x(pi_0)).
        This is valid because d_t A = pi, so d_x(d_t A) = d_x(pi).

        The state is organized as: [field_0, pi_0, field_1, pi_1, ...]
        - Regular fields (A_0, phi_0, phi0, etc.) are at even indices: state[2*i]
        - Momentum fields (pi_0, pi0, pi_1, pi1, etc.) are at odd indices: state[2*i + 1]

        Supports flexible field naming conventions:
        - Standard: A_0, phi_1
        - Compact: A0, phi1, pi0, pi1
        - Tensor: stress_xy_0
        - Simple: phi (index defaults to 0)

        Parameters
        ----------
        state : FieldCollection
            Current state with interleaved field/momentum pairs.
        field_name : str
            Name of the field to retrieve. Can be:
            - Regular field name like "A_0", "phi_0", "phi0"
            - Momentum field name like "pi_0", "pi_1", "pi0", "pi1"

        Returns
        -------
        ScalarField
            The requested field from the state.

        Raises
        ------
        ValueError
            If the field name is not recognized.
        """
        # Check if this is a momentum field reference (pi_0, pi0, etc.)
        # First check if it looks like a momentum field (starts with "pi")
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

                # Momentum is at odd indices: state[2*i + 1]
                momentum = state[2 * momentum_idx + 1]
                assert isinstance(momentum, ScalarField)
                return momentum
            # If it looks like momentum but couldn't be parsed, raise clear error
            msg = (
                f"Invalid momentum field format: '{field_name}'. "
                f"Expected 'pi_N' or 'piN' where N is a numeric index (e.g., 'pi_0', 'pi0')."
            )
            raise ValueError(msg)

        # Regular field lookup - first try direct name match
        target_idx = self._component_name_to_index.get(field_name)
        if target_idx is not None:
            # Field is at even indices: state[2*i]
            field = state[2 * target_idx]
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
    ) -> ScalarField:
        """Compute the RHS for a single component's momentum equation.

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

        Returns
        -------
        ScalarField
            The computed RHS for d/dt momentum_i.

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
            # Find which field this term operates on
            # Supports both regular fields (A_0, phi_0) and momentum fields (pi_0, pi_1)
            # Momentum field support enables mixed time-space derivative handling
            target_field_name = term.field

            # Handle first_derivative_t operator specially
            # For d_t(field_i) = momentum_i, so first_derivative_t(field) returns momentum
            if term.operator == "first_derivative_t":
                # Get the field index for this term's target field
                target_idx = self._component_name_to_index.get(target_field_name)
                if target_idx is not None:
                    # Momentum is at odd indices: state[2*i + 1]
                    momentum = state[2 * target_idx + 1]
                    assert isinstance(momentum, ScalarField)
                    operated = momentum.copy()
                else:
                    msg = f"Unknown field for first_derivative_t: {target_field_name}"
                    raise ValueError(msg)
            else:
                # Standard operator handling
                target_field = self._get_field_from_state(state, target_field_name)
                operated = self._get_operator(term.operator, target_field, bc)

            # Resolve coefficient (use time-dependent resolution if needed)
            coefficient = self._resolve_coefficient_at_time(term, t)

            # Add coefficient * operated to result
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

        For each component i:
            d/dt field_i = momentum_i
            d/dt momentum_i = RHS from specification

        Parameters
        ----------
        state : FieldCollection
            Current state with 2 * n_components fields.
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
            If the state does not have exactly 2 * n_components fields.
        """
        assert isinstance(state, FieldCollection)
        expected_fields = 2 * self.n_components
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
        rates: list[ScalarField] = []

        for i in range(self.n_components):
            field_i = state[2 * i]
            momentum_i = state[2 * i + 1]

            assert isinstance(field_i, ScalarField)
            assert isinstance(momentum_i, ScalarField)

            # d/dt field_i = momentum_i
            d_field_dt = momentum_i.copy()

            # d/dt momentum_i = RHS from specification
            # Pass current time for time-dependent coefficients (curved spacetime)
            d_momentum_dt = self._compute_rhs_for_component(i, state, bc, t)

            rates.extend((d_field_dt, d_momentum_dt))

        return FieldCollection(rates)

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
        Components not specified default to zero.

    Returns
    -------
    FieldCollection
        Initial state with 2 * n_components fields.
    """
    field_data = field_data or {}
    momentum_data = momentum_data or {}

    fields: list[ScalarField] = []

    for name in spec.component_names:
        # Field
        if name in field_data:
            field = ScalarField(grid, data=field_data[name])
        else:
            field = ScalarField(grid, data=0.0)

        # Momentum
        if name in momentum_data:
            momentum = ScalarField(grid, data=momentum_data[name])
        else:
            momentum = ScalarField(grid, data=0.0)

        fields.extend((field, momentum))

    return FieldCollection(fields)
