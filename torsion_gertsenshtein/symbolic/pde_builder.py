"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

import math
import operator
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, SupportsFloat, cast

import numpy as np
from pde import FieldCollection, PDEBase, ScalarField
from scipy import special  # type: ignore[reportMissingTypeStubs]
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid
from torsion_gertsenshtein.symbolic.json_loader import (
    _CUSTOM_OPERATORS,  # type: ignore[reportPrivateUsage]
    AXIS_LETTERS,
    BoundaryCondition,
    ConstraintSolverConfig,
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
# Pre-compiled regex patterns for field name parsing (avoid re-compilation
# on every call to ParsedFieldName.parse / parse_momentum_field_name)
# ---------------------------------------------------------------------------
_STANDARD_FORMAT_RE = re.compile(r"^([a-zA-Z]+)_([0-9]+)$")
_TENSOR_FORMAT_RE = re.compile(r"^(.+)_([0-9]+)$")
_COMPACT_FORMAT_RE = re.compile(r"^([a-zA-Z]+)([0-9]+)$")
_SIMPLE_FORMAT_RE = re.compile(r"^[a-zA-Z]+$")
_MOMENTUM_STANDARD_RE = re.compile(r"^pi_([0-9]+)$")
_MOMENTUM_COMPACT_RE = re.compile(r"^pi([0-9]+)$")


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
        if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
            raise TypeError(msg)
        return component

    return _handler


def _op_directional_laplacian(
    axis: int,
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a directional Laplacian handler (∂²/∂x_i²)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)[axis]
        if not isinstance(grad, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad).__name__}"
            raise TypeError(msg)
        d2 = grad.gradient(bc=bc)[axis]
        if not isinstance(d2, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(d2).__name__}"
            raise TypeError(msg)
        return d2

    return _handler


def _op_cross_derivative(
    axis1: int, axis2: int
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a cross derivative handler (∂²/∂x_i ∂x_j)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad_j = field.gradient(bc=bc)[axis2]
        if not isinstance(grad_j, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad_j).__name__}"
            raise TypeError(msg)
        grad_ij = grad_j.gradient(bc=bc)[axis1]
        if not isinstance(grad_ij, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad_ij).__name__}"
            raise TypeError(msg)
        return grad_ij

    return _handler


def _op_biharmonic(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    """Biharmonic operator: ∇⁴f = ∇²(∇²f).

    Raises
    ------
    TypeError
        If intermediate results are not ``ScalarField``.
    """
    lap = field.laplace(bc=bc)
    if not isinstance(lap, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = f"Expected ScalarField from laplace, got {type(lap).__name__}"
        raise TypeError(msg)
    bilap = lap.laplace(bc=bc)
    if not isinstance(bilap, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = f"Expected ScalarField from laplace, got {type(bilap).__name__}"
        raise TypeError(msg)
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
            if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
                msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
                raise TypeError(msg)
            result = component
        return result

    return _handler


def _parse_multi_axis_spec(spec: str) -> list[tuple[int, int]]:
    """Parse a multi-axis derivative spec into (axis_index, order) pairs.

    Parameters
    ----------
    spec : str
        Multi-axis spec like ``"2x_1y"``, ``"3x_2z"``, ``"1x_1y_1z"``.

    Returns
    -------
    list[tuple[int, int]]
        List of (axis_index, order) pairs, sorted by axis index.

    Raises
    ------
    ValueError
        If the spec cannot be parsed or contains invalid axis names.
    """
    axis_pattern = re.compile(r"^(\d+)(" + _AXIS_RE_CLASS + r")$")
    parts = spec.split("_")
    result: list[tuple[int, int]] = []
    for part in parts:
        match = axis_pattern.match(part)
        if not match:
            msg = (
                f"Invalid multi-axis derivative part: '{part}'. "
                f"Expected format like '2x', '1y', '3z'."
            )
            raise ValueError(msg)
        order = int(match.group(1))
        axis_letter = match.group(2)
        axis = _AXIS_INDEX[axis_letter]
        result.append((axis, order))

    # Detect duplicate axes — e.g. "2x_1x_1y" has x twice
    seen_axes: set[int] = set()
    for axis, _order in result:
        if axis in seen_axes:
            axis_name = _AXIS_LETTER[axis]
            msg = (
                f"Duplicate axis '{axis_name}' in multi-axis derivative spec: '{spec}'. "
                f"Combine orders for the same axis into a single entry."
            )
            raise ValueError(msg)
        seen_axes.add(axis)

    return sorted(result, key=operator.itemgetter(0))


def _op_multi_axis_derivative(
    axes_and_orders: list[tuple[int, int]],
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a handler for multi-axis mixed spatial derivatives.

    Applies gradients sequentially for each (axis, order) pair.
    E.g., for ``derivative_2x_1y``: applies gradient in y once, then
    gradient in x twice, giving d^2/dx^2 d/dy.
    """

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        result: ScalarField = field
        for axis, order in axes_and_orders:
            for _ in range(order):
                grad = result.gradient(bc=bc)
                component = grad[axis]
                if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
                    msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
                    raise TypeError(msg)
                result = component
        return result

    return _handler


#: Map axis letter to axis index (derived from AXIS_LETTERS in json_loader).
_AXIS_INDEX: dict[str, int] = {letter: i for i, letter in enumerate(AXIS_LETTERS)}

#: Minimum grid dimension required for each axis.
_AXIS_MIN_DIM: dict[str, int] = {letter: i + 1 for i, letter in enumerate(AXIS_LETTERS)}

#: Map axis index back to letter.
_AXIS_LETTER: dict[int, str] = dict(enumerate(AXIS_LETTERS))

#: Character class matching all known axis letters.
_AXIS_RE_CLASS = "[" + "".join(AXIS_LETTERS) + "]"

#: Regex for parsing generic single-axis derivative names.
_GENERIC_SINGLE_RE = re.compile(
    r"^derivative_(\d+)_(" + _AXIS_RE_CLASS + r")$"
)

#: Regex for parsing generic multi-axis derivative names.
_GENERIC_MULTI_RE = re.compile(
    r"^derivative_(\d+" + _AXIS_RE_CLASS + r"(?:_\d+" + _AXIS_RE_CLASS + r")*)$"
)


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
    # first_derivative_t is handled specially in _compute_rhs_for_component.
    # This sentinel entry ensures _get_operator gives a clear error if called directly.
    "first_derivative_t": (None, 1),
}


def register_operator(
    name: str,
    handler: Callable[[ScalarField, BCDescriptor], ScalarField],
    min_dim: int = 1,
) -> None:
    """Register a custom spatial operator for use in the PDE pipeline.

    The operator becomes available both for JSON validation (``is_known_operator``)
    and for runtime evaluation (``PDEFromSpec._get_operator``).

    Parameters
    ----------
    name : str
        Operator name as it appears in JSON (e.g., ``"my_diffusion"``).
    handler : callable
        Function with signature ``(field: ScalarField, bc) -> ScalarField``.
    min_dim : int
        Minimum spatial grid dimension required for this operator.

    Raises
    ------
    ValueError
        If *name* would shadow a built-in operator.
    """
    if name in _OPERATOR_REGISTRY:
        msg = (
            f"Cannot register operator '{name}': "
            f"it shadows a built-in operator. Choose a different name."
        )
        raise ValueError(msg)
    _OPERATOR_REGISTRY[name] = (handler, min_dim)
    _CUSTOM_OPERATORS.add(name)


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
        match = _STANDARD_FORMAT_RE.match(name)
        if match:
            return cls(
                base=match.group(1), index=int(match.group(2)), format="standard"
            )

        # Tensor format: stress_xy_0, u_x_1 (greedy match for base)
        match = _TENSOR_FORMAT_RE.match(name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="tensor")

        # Compact format: phi0, A1
        match = _COMPACT_FORMAT_RE.match(name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="compact")

        # Simple format: phi, psi (no index, defaults to 0)
        if _SIMPLE_FORMAT_RE.match(name):
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
    match = _MOMENTUM_STANDARD_RE.match(field_name)
    if match:
        return int(match.group(1))

    # Compact format: piN
    match = _MOMENTUM_COMPACT_RE.match(field_name)
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
        *,
        constraint_eps: float = 1e-14,
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
        constraint_eps : float
            Tolerance for determining whether a Laplacian coefficient is
            effectively zero in the constraint solver.  Default ``1e-14``.
        """
        super().__init__()
        self.spec = spec
        self.n_components = spec.n_components
        self._component_name_to_index = {
            name: i for i, name in enumerate(spec.component_names)
        }
        self._parameters = parameters or {}
        self._constraint_eps = constraint_eps

        # Build slot maps from state_layout for mixed time-order support
        self._field_slot_map: dict[str, int] = {}
        self._momentum_slot_map: dict[str, int] = {}
        for slot_idx, (name, slot_type) in enumerate(spec.state_layout):
            if slot_type == "field":
                self._field_slot_map[name] = slot_idx
            else:
                self._momentum_slot_map[name] = slot_idx

        # B1: Cache for _mathematica_to_python() results (same symbolic expr → same output)
        self._expr_cache: dict[str, str] = {}

        # B2: Pre-build the static part of the coefficient namespace (math functions,
        # parameters). Only t and grid coordinates change per-call.
        self._base_namespace: dict[str, Any] = self._build_base_namespace()

        # B4: Pre-resolve constant coefficients (not position- or time-dependent).
        # Only pre-resolve when the symbolic coefficient is absent or resolvable
        # from the provided parameters — avoids emitting premature warnings for
        # terms whose symbolic cannot be resolved (the warning fires at runtime).
        self._preresolved: dict[tuple[int, int], float] = {}
        for eq_idx, eq in enumerate(spec.equations):
            for term_idx, term in enumerate(eq.rhs_terms):
                if (
                    not term.time_dependent
                    and not term.position_dependent
                    and self._is_resolvable(term)
                ):
                    self._preresolved[eq_idx, term_idx] = self._resolve_coefficient(
                        term
                    )

        # B5: Cache boundary conditions and grid coordinates (populated on first call)
        self._cached_bc: BCDescriptor | None = None
        self._cached_grid_id: int | None = None

        # Validate operator dimension requirements against spec at construction time
        self._validate_operator_dimensions()

    @staticmethod
    def _operator_min_dim(operator_name: str) -> int:
        """Return minimum spatial grid dimension required by an operator.

        Uses the same resolution logic as ``_get_operator`` (static registry
        then dynamic regex patterns) but without needing a grid or field.

        Parameters
        ----------
        operator_name : str
            Operator name (e.g. ``"gradient_y"``, ``"derivative_3_z"``).

        Returns
        -------
        int
            Minimum spatial dimension (1, 2, or 3).

        Raises
        ------
        ValueError
            If the operator name is not recognized.
        """
        entry = _OPERATOR_REGISTRY.get(operator_name)
        if entry is not None:
            return entry[1]
        m = _GENERIC_SINGLE_RE.match(operator_name)
        if m:
            return _AXIS_MIN_DIM[m.group(2)]
        m_multi = _GENERIC_MULTI_RE.match(operator_name)
        if m_multi:
            axes_and_orders = _parse_multi_axis_spec(m_multi.group(1))
            max_axis = max(axis for axis, _ in axes_and_orders)
            return _AXIS_MIN_DIM[_AXIS_LETTER[max_axis]]
        msg = (
            f"Unknown operator: '{operator_name}'. "
            f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}. "
            f"Dynamic patterns: derivative_N_x, derivative_Nx_My."
        )
        raise ValueError(msg)

    def _validate_operator_dimensions(self) -> None:
        """Validate all operators are compatible with the spec's spatial dimension.

        Called during ``__init__`` to fail fast when a JSON spec contains
        operators that require more spatial dimensions than the spec provides
        (e.g. ``gradient_z`` in a 2+1D spacetime).

        Raises
        ------
        ValueError
            If any operator requires more dimensions than ``spec.spatial_dimension``.
        """
        spatial_dim = self.spec.spatial_dimension
        for eq in self.spec.equations:
            for term in eq.rhs_terms:
                min_dim = self._operator_min_dim(term.operator)
                if min_dim > spatial_dim:
                    msg = (
                        f"Operator '{term.operator}' in equation for "
                        f"'{eq.field_name}' requires at least {min_dim}D "
                        f"spatial grid, but the spec has "
                        f"spatial_dimension={spatial_dim} "
                        f"(from {self.spec.dimension}D spacetime)."
                    )
                    raise ValueError(msg)

    def _is_resolvable(self, term: OperatorTerm) -> bool:
        """Check whether a term's coefficient can be resolved without warnings.

        Returns True if the term has no symbolic coefficient, or if its symbolic
        coefficient (possibly negated) matches a key in the parameters dict.
        """
        sym = term.coefficient_symbolic
        if sym is None:
            return True
        if sym.startswith("-") and sym[1:] in self._parameters:
            return True
        return sym in self._parameters

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

        Raises
        ------
        ValueError
            If a symbolic coefficient cannot be resolved from parameters.
        """
        if term.coefficient_symbolic is not None:
            sym = term.coefficient_symbolic

            # Check for negated symbol like "-m2"
            if sym.startswith("-") and sym[1:] in self._parameters:
                return -self._parameters[sym[1:]]
            if sym in self._parameters:
                return self._parameters[sym]

            # Symbolic coefficient present but not resolvable from parameters.
            # Fail fast: wrong physics that looks right is worse than a crash.
            msg = (
                f"Symbolic coefficient '{sym[:80]}' could not be resolved from "
                f"parameters {sorted(self._parameters.keys())}. "
                f"Pass the required parameter via "
                f"parameters={{'{sym.lstrip('-')}': <value>}} to "
                f"PDEFromSpec or build_pde_from_json."
            )
            raise ValueError(msg)

        # Default: use numeric coefficient from JSON (no symbolic → purely numeric term)
        return term.coefficient

    def _build_base_namespace(self) -> dict[str, Any]:
        """Build the static part of the coefficient evaluation namespace.

        Contains numpy/scipy math functions and user-provided parameters.
        The time variable ``t`` and spatial grid coordinates are injected
        per-call in ``_resolve_coefficient_at_point``.
        """
        ns: dict[str, Any] = dict(self._parameters)
        ns["exp"] = np.exp
        # Basic trig
        ns["sin"] = np.sin
        ns["cos"] = np.cos
        ns["tan"] = np.tan
        # Reciprocal trig (no direct numpy equivalents)
        ns["cot"] = lambda x: np.cos(x) / np.sin(x)  # type: ignore[reportUnknownLambdaType]
        ns["sec"] = lambda x: 1.0 / np.cos(x)  # type: ignore[reportUnknownLambdaType]
        ns["csc"] = lambda x: 1.0 / np.sin(x)  # type: ignore[reportUnknownLambdaType]
        # Inverse trig
        ns["arcsin"] = np.arcsin
        ns["arccos"] = np.arccos
        ns["arctan"] = np.arctan
        ns["arctan2"] = np.arctan2
        # Hyperbolic
        ns["sinh"] = np.sinh
        ns["cosh"] = np.cosh
        ns["tanh"] = np.tanh
        # Inverse hyperbolic
        ns["arcsinh"] = np.arcsinh
        ns["arccosh"] = np.arccosh
        ns["arctanh"] = np.arctanh
        # Other
        ns["log"] = np.log
        ns["sqrt"] = np.sqrt
        ns["abs"] = np.abs
        # Special functions (scipy.special)
        ns["erf"] = special.erf
        ns["jv"] = special.jv  # BesselJ
        ns["yv"] = special.yv  # BesselY
        return ns

    @staticmethod
    def _convert_power_function(expr: str) -> str:
        """Convert Power[base, exponent] to (base)**(exponent).

        Handles nested expressions in both arguments. Must run before bracket
        conversion ([] → ()).

        Parameters
        ----------
        expr : str
            Expression potentially containing Power[...] syntax.

        Returns
        -------
        str
            Expression with Power[...] converted to (...)**(...).
        """
        # Pattern: Power[<arg1>, <arg2>] where args may contain nested brackets
        # (?:[^[\]]|\[[^\]]*\]) matches either non-bracket chars or [...] pairs
        pattern = r"Power\[((?:[^[\]]|\[[^\]]*\])*),\s*((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            base = match.group(1).strip()
            exp = match.group(2).strip()
            return f"({base})**({exp})"

        # Multiple passes for nested Power calls
        prev = None
        result = expr
        while prev != result:
            prev = result
            result = re.sub(pattern, replacer, result)
        return result

    @staticmethod
    def _convert_arctan2(expr: str) -> str:
        """Convert ArcTan[x, y] to arctan2(y, x) with argument swap.

        Mathematica's ArcTan[x, y] computes atan2(y, x), so arguments must
        be swapped during conversion to match NumPy's arctan2(y, x) signature.

        Parameters
        ----------
        expr : str
            Expression potentially containing ArcTan[x, y] syntax.

        Returns
        -------
        str
            Expression with ArcTan[x, y] converted to arctan2(y, x).
        """
        # Pattern: ArcTan[<arg1>, <arg2>] - must handle before generic function conversion
        pattern = r"ArcTan\[((?:[^[\],]|\[[^\]]*\])*),\s*((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            x = match.group(1).strip()
            y = match.group(2).strip()
            return f"arctan2({y}, {x})"  # Swap x and y!

        return re.sub(pattern, replacer, expr)

    def _mathematica_to_python(self, expr: str) -> str:
        """Convert Mathematica InputForm expression to evaluable Python.

        Handles common Mathematica syntax:
        - ``E^(...)`` to ``exp(...)`` (Euler's number)
        - ``Power[x,y]`` to ``(x)**(y)`` (function form of exponentiation)
        - ``Sin[x]`` to ``sin(x)``, ``Cos[x]`` to ``cos(x)``, ``Tan[x]`` to ``tan(x)``
        - ``Cot[x]`` to ``cot(x)``, ``Sec[x]`` to ``sec(x)``, ``Csc[x]`` to ``csc(x)``
        - ``ArcSin[x]`` to ``arcsin(x)``, ``ArcCos[x]`` to ``arccos(x)``, etc.
        - ``ArcTan[x, y]`` to ``arctan2(y, x)`` (note argument order swap!)
        - ``Sinh[x]`` to ``sinh(x)``, ``Cosh[x]`` to ``cosh(x)``, etc.
        - ``ArcSinh[x]`` to ``arcsinh(x)``, ``ArcCosh[x]`` to ``arccosh(x)``, etc.
        - ``Erf[x]`` to ``erf(x)`` (scipy.special)
        - ``BesselJ[n, x]`` to ``jv(n, x)``, ``BesselY[n, x]`` to ``yv(n, x)``
        - ``Rational[p, q]`` to ``(p)/(q)`` (exact fractions from xAct InputForm)
        - ``Pi`` to ``np.pi`` (mathematical constant)
        - ``Sign[x]`` to ``np.sign(x)``, ``Max[a,b]`` to ``np.maximum(a,b)``
        - ``t[]`` to ``t`` (xCoba coordinate symbols, using actual coordinate names)
        - Mathematica brackets ``[``, ``]`` to Python parens ``(``, ``)``
        - Mathematica ``^`` to Python ``**``

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

        # Step 1: E^(...) → exp(...) — Mathematica's Euler number
        result = re.sub(r"\bE\^", "exp", result)

        # Step 2: Power[x, y] → (x)**(y) — must handle before bracket conversion
        result = self._convert_power_function(result)

        # Step 3: ArcTan[x, y] → arctan2(y, x) — special 2-arg case with swap
        # Must handle before generic function conversion
        result = self._convert_arctan2(result)

        # Step 3.5: Rational[p, q] → (p)/(q) — Mathematica exact fractions
        # xAct outputs Rational[1,2] for 1/2 in InputForm. Must handle before
        # bracket conversion since Rational uses Mathematica brackets.
        result = re.sub(
            r"Rational\[([^,\]]+),\s*([^,\]]+)\]", r"(\1)/(\2)", result
        )

        # Step 4: Function name conversions (batch)
        function_map = [
            # Basic trig
            ("Sin", "sin"),
            ("Cos", "cos"),
            ("Tan", "tan"),
            # Reciprocal trig
            ("Cot", "cot"),
            ("Sec", "sec"),
            ("Csc", "csc"),
            # Inverse trig (1-arg)
            ("ArcSin", "arcsin"),
            ("ArcCos", "arccos"),
            ("ArcTan", "arctan"),  # 1-arg version only (2-arg handled above)
            # Hyperbolic
            ("Sinh", "sinh"),
            ("Cosh", "cosh"),
            ("Tanh", "tanh"),
            # Inverse hyperbolic
            ("ArcSinh", "arcsinh"),
            ("ArcCosh", "arccosh"),
            ("ArcTanh", "arctanh"),
            # Other
            ("Log", "log"),
            ("Sqrt", "sqrt"),
            ("Abs", "abs"),
            ("Sign", "np.sign"),
            ("Max", "np.maximum"),
            ("Min", "np.minimum"),
            # Special functions (scipy.special)
            ("Erf", "erf"),
            ("BesselJ", "jv"),
            ("BesselY", "yv"),
        ]
        for mma_func, py_func in function_map:
            result = re.sub(rf"\b{mma_func}\b", py_func, result)

        # Step 4.5: Pi → np.pi — Mathematica constant (before bracket conversion)
        result = re.sub(r"\bPi\b", "np.pi", result)

        # Step 5: Mathematica brackets to Python parens (after function renaming)
        result = result.replace("[", "(").replace("]", ")")

        # Step 6: Mathematica ^ to Python ** (AFTER E^ → exp to avoid double-conversion)
        result = result.replace("^", "**")

        # Step 7: xCoba coordinate symbols: t() → t, x() → x, etc.
        # Uses actual coordinate names from equation system (not hardcoded x/y/z).
        for coord in self.spec.effective_coordinates:
            result = result.replace(f"{coord}()", coord)

        return result

    def _resolve_coefficient_at_point(
        self,
        term: OperatorTerm,
        t: float,
        grid: GridBase | None = None,
        coord_arrays: dict[str, NumericArray] | None = None,
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

        # B1: Use cached Mathematica→Python conversion
        if sym not in self._expr_cache:
            self._expr_cache[sym] = self._mathematica_to_python(sym)
        py_expr = self._expr_cache[sym]

        # B2: Clone pre-built base namespace and inject dynamic variables
        namespace: dict[str, Any] = dict(self._base_namespace)
        namespace["t"] = t

        # C2: Inject spatial coordinates — use pre-extracted arrays if available
        if term.position_dependent:
            if coord_arrays is not None:
                namespace.update(coord_arrays)
            else:
                spatial_coords = self.spec.spatial_coordinates
                coords = grid.cell_coords  # type: ignore[union-attr]
                for i, name in enumerate(spatial_coords[: grid.num_axes]):  # type: ignore[union-attr]
                    namespace[name] = np.asarray(coords[..., i], dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]

        # Validate all symbols can be resolved
        identifiers = set(re.findall(r"\b[a-zA-Z_]\w*\b", py_expr))
        # Derive builtin names from namespace (all math functions, excluding parameters and t)
        builtin_names = set(namespace.keys()) - set(self._parameters.keys()) - {"t"}
        # Exclude coordinate variables if position-dependent
        if term.position_dependent:
            builtin_names -= set(self.spec.spatial_coordinates)
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
        except Exception as e:
            msg = (
                f"Cannot evaluate coordinate-dependent coefficient '{sym}' "
                f"(Python form: '{py_expr}') at t={t}: {e}"
            )
            raise ValueError(msg) from e
        return self._validate_eval_result(result, sym, py_expr)

    @staticmethod
    def _validate_eval_result(
        result: object, sym: str, py_expr: str
    ) -> float | NumericArray:
        """Validate and coerce an eval() result to float or ndarray.

        Raises ValueError for complex, NaN, or Inf results with clear
        diagnostic messages pointing to the source expression.

        Raises
        ------
        TypeError
            If the result is complex.
        ValueError
            If the result is NaN or Inf.
        """
        if isinstance(result, complex):
            msg = (
                f"Coefficient '{sym}' evaluated to complex number {result} "
                f"(from '{py_expr}'). Only real-valued coefficients are supported."
            )
            raise TypeError(msg)

        if isinstance(result, np.ndarray):
            arr = np.asarray(result, dtype=np.float64)
            if np.any(np.isnan(arr)):
                msg = (
                    f"Coefficient '{sym}' produced NaN values "
                    f"(from '{py_expr}'). Check for 0/0 or invalid operations."
                )
                raise ValueError(msg)
            if np.any(np.isinf(arr)):
                msg = (
                    f"Coefficient '{sym}' produced Inf values "
                    f"(from '{py_expr}'). Check for division by zero."
                )
                raise ValueError(msg)
            return arr

        scalar = float(cast("SupportsFloat", result))
        if math.isnan(scalar):
            msg = (
                f"Coefficient '{sym}' evaluated to NaN "
                f"(from '{py_expr}'). Check for 0/0 or invalid operations."
            )
            raise ValueError(msg)
        if math.isinf(scalar):
            msg = (
                f"Coefficient '{sym}' evaluated to Inf "
                f"(from '{py_expr}'). Check for division by zero."
            )
            raise ValueError(msg)
        return scalar

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
        RuntimeError
            If the operator is in the registry but cannot be applied as a spatial operator.
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
                # Try multi-axis pattern: derivative_2x_1y, derivative_1x_1y_1z
                m_multi = _GENERIC_MULTI_RE.match(operator_name)
                if m_multi:
                    axes_and_orders = _parse_multi_axis_spec(m_multi.group(1))
                    max_axis = max(axis for axis, _ in axes_and_orders)
                    min_dim = _AXIS_MIN_DIM[_AXIS_LETTER[max_axis]]
                    entry = (_op_multi_axis_derivative(axes_and_orders), min_dim)
                else:
                    msg = (
                        f"Unknown operator: '{operator_name}'. "
                        f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}. "
                        f"Dynamic patterns: derivative_N_x, derivative_Nx_My."
                    )
                    raise ValueError(msg)

        handler, min_dim = entry

        # Sentinel check: some operators (first_derivative_t) are in the registry
        # but handled specially elsewhere — they cannot be applied as spatial operators.
        if handler is None:
            msg = (
                f"Operator '{operator_name}' cannot be applied as a spatial operator. "
                f"It is handled specially in _compute_rhs_for_component."
            )
            raise RuntimeError(msg)

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
        TypeError
            If a state element is not a ``ScalarField``.
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
                    if not isinstance(momentum, ScalarField):
                        msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                        raise TypeError(msg)
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
            if not isinstance(field, ScalarField):
                msg = f"Expected ScalarField, got {type(field).__name__}"
                raise TypeError(msg)
            return field

        msg = f"Unknown field name: {field_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(  # noqa: C901, PLR0912, PLR0914
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
        TypeError
            If a state element is not a ``ScalarField``.
        ValueError
            If a field or operator cannot be resolved.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # C2: Pre-extract spatial coordinate arrays once for all position-dependent terms
        coord_arrays: dict[str, NumericArray] | None = None
        if any(term.position_dependent for term in eq.rhs_terms):
            spatial_coords = self.spec.spatial_coordinates
            raw_coords = grid.cell_coords  # pyright: ignore[reportUnknownVariableType]
            coord_arrays = {
                name: np.asarray(raw_coords[..., i], dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]
                for i, name in enumerate(spatial_coords[: grid.num_axes])
            }

        # C1: Per-timestep coefficient cache — same symbolic expression at
        # same (t, grid) produces the same value, so deduplicate eval() calls
        coeff_cache: dict[str, float | NumericArray] = {}

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term_idx, term in enumerate(eq.rhs_terms):
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
                    if not isinstance(momentum, ScalarField):
                        msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                        raise TypeError(msg)
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

            # Resolve coefficient: B4 preresolved → C1 timestep cache → full eval
            preresolved = self._preresolved.get((component_idx, term_idx))
            if preresolved is not None:
                coefficient: float | NumericArray = preresolved
            else:
                cache_key = term.coefficient_symbolic
                if cache_key is not None and cache_key in coeff_cache:
                    coefficient = coeff_cache[cache_key]
                else:
                    coefficient = self._resolve_coefficient_at_point(
                        term, t, grid, coord_arrays=coord_arrays
                    )
                    if cache_key is not None:
                        coeff_cache[cache_key] = coefficient

            # Add coefficient * operated to result
            if isinstance(coefficient, np.ndarray):
                contribution = ScalarField(grid, data=coefficient * operated.data)
            else:
                contribution = coefficient * operated
            result += contribution

        return result

    def _build_constraint_bc(
        self,
        config: ConstraintSolverConfig,
        grid: GridBase,
    ) -> Any:  # noqa: ANN401
        """Convert a ConstraintSolverConfig to a py-pde boundary condition.

        Parameters
        ----------
        config : ConstraintSolverConfig
            Constraint solver configuration with per-axis BCs.
        grid : GridBase
            The simulation grid (used to check periodicity).

        Returns
        -------
        str | dict
            py-pde boundary condition descriptor.
        """
        spatial_coords = self.spec.spatial_coordinates

        # If all BCs are periodic and grid is periodic, use shorthand
        all_periodic = all(
            config.boundary_conditions.get(coord, BoundaryCondition("periodic")).type
            == "periodic"
            for coord in spatial_coords
        )
        if all_periodic and hasattr(grid, "periodic") and all(grid.periodic):
            return "auto_periodic_neumann"

        # Build explicit per-axis BC dict using coordinate names
        # py-pde expects: {"x": bc_x, "y": bc_y} where bc is "periodic"
        # or {"value": V} or {"derivative": D}
        bc_dict: dict[str, Any] = {}
        for i, coord in enumerate(spatial_coords):
            bc_config = config.boundary_conditions.get(coord)
            if bc_config is None:
                # Default: periodic if grid is periodic on this axis, else Neumann
                if hasattr(grid, "periodic") and grid.periodic[i]:
                    bc_dict[coord] = "periodic"
                else:
                    bc_dict[coord] = {"derivative": 0.0}
            elif bc_config.type == "periodic":
                bc_dict[coord] = "periodic"
            elif bc_config.type == "dirichlet":
                bc_dict[coord] = {
                    "value": bc_config.value if bc_config.value is not None else 0.0
                }
            elif bc_config.type == "neumann":
                bc_dict[coord] = {
                    "derivative": bc_config.derivative
                    if bc_config.derivative is not None
                    else 0.0
                }

        return bc_dict

    def _solve_constraint_equation(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
    ) -> FieldCollection:
        """Solve an elliptic constraint equation and update the state.

        For constraint equations in the form::

            0 = laplacian_coeff * laplacian(field) + source_terms

        This rearranges to standard Poisson form::

            nabla^2 field = -source_terms / laplacian_coeff

        and solves using py-pde's ``solve_poisson_equation``.

        Parameters
        ----------
        component_idx : int
            Index of the constraint equation in ``spec.equations``.
        state : FieldCollection
            Current state. A new FieldCollection is returned with the
            constraint field replaced by the solution.
        bc : BCDescriptor
            Boundary conditions for evaluating source-term operators.
        t : float
            Current time (for time-dependent source coefficients).

        Returns
        -------
        FieldCollection
            Updated state with the constraint field solved.

        Raises
        ------
        TypeError
            If the Poisson RHS is not a ``ScalarField``.
        ValueError
            If the equation lacks a ``laplacian(field)`` term or the
            Poisson solver fails.
        """
        from pde import solve_poisson_equation  # noqa: PLC0415, I001  # type: ignore[reportUnknownVariableType]

        eq = self.spec.equations[component_idx]
        grid = state.grid
        field_slot = self._field_slot_map[eq.field_name]

        # Separate RHS into the laplacian-of-self term and source terms
        laplacian_coeff: float | None = None
        source_terms: list[OperatorTerm] = []

        for term in eq.rhs_terms:
            if term.operator == "laplacian" and term.field == eq.field_name:
                if laplacian_coeff is not None:
                    msg = (
                        f"Multiple laplacian({eq.field_name}) terms in constraint "
                        f"equation. Expected exactly one."
                    )
                    raise ValueError(msg)
                laplacian_coeff = self._resolve_coefficient(term)
            else:
                source_terms.append(term)

        # Validate equation structure
        if laplacian_coeff is None:
            msg = (
                f"Constraint equation for {eq.field_name} lacks a "
                f"laplacian({eq.field_name}) term. "
                f"The elliptic solver requires the form: "
                f"laplacian(field) + source = 0."
            )
            raise ValueError(msg)

        if abs(laplacian_coeff) < self._constraint_eps:
            msg = (
                f"Laplacian coefficient for {eq.field_name} is effectively "
                f"zero ({laplacian_coeff}). Cannot solve elliptic equation."
            )
            raise ValueError(msg)

        # Compute source: S = sum(coeff_i * operator_i(field_i))
        rhs_source = ScalarField(grid, data=0.0)
        for term in source_terms:
            target_field = self._get_field_from_state(state, term.field)
            operated = self._get_operator(term.operator, target_field, bc)
            coefficient = self._resolve_coefficient_at_point(term, t, grid)

            if isinstance(coefficient, np.ndarray):
                contribution = ScalarField(grid, data=coefficient * operated.data)
            else:
                contribution = coefficient * operated
            rhs_source += contribution

        # Rearrange 0 = laplacian_coeff * nabla^2(phi) + S into the standard
        # Poisson form nabla^2(phi) = rhs, giving rhs = -S / laplacian_coeff.
        poisson_rhs = -rhs_source / laplacian_coeff
        if not isinstance(poisson_rhs, ScalarField):
            msg = f"Expected ScalarField for Poisson RHS, got {type(poisson_rhs).__name__}"
            raise TypeError(msg)

        # Build boundary conditions for the Poisson solver
        solver_bc = self._build_constraint_bc(eq.constraint_solver, grid)

        try:
            solution = solve_poisson_equation(
                rhs=poisson_rhs,
                bc=solver_bc,
                label=eq.field_name,
            )
        except Exception as e:
            rhs_max = float(np.max(np.abs(poisson_rhs.data)))
            msg = (
                f"Poisson solver failed for constraint {eq.field_name}:\n"
                f"  RHS max |f|: {rhs_max:.3e}\n"
                f"  BC: {solver_bc}\n"
                f"  Error: {e}"
            )
            raise ValueError(msg) from e

        # Update state in-place so dynamical equations (Pass 2) see the
        # solved constraint field via cross-field references.
        state[field_slot].data[:] = solution.data
        return state

    def _evolve_constraints(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
    ) -> tuple[TState, dict[str, ScalarField]]:
        """Solve constraint equations (time_derivative_order=0).

        Constraints are solved elliptically when enabled, updating the state
        in-place so that dynamical equations see the resolved fields.

        Returns
        -------
        tuple[TState, dict[str, ScalarField]]
            Updated state and virtual momenta dict with zero entries for
            each constraint field.

        Raises
        ------
        TypeError
            If state is not a FieldCollection.
        """
        # Type narrowing: constraints require FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "Constraint solving requires FieldCollection state"
            raise TypeError(msg)

        grid = state.grid
        virtual_momenta: dict[str, ScalarField] = {}
        for i, eq in enumerate(self.spec.equations):
            if eq.time_derivative_order == 0:
                if eq.constraint_solver.enabled:
                    state = self._solve_constraint_equation(i, state, bc, t)
                virtual_momenta[eq.field_name] = ScalarField(grid, data=0.0)
        return state, virtual_momenta

    def _evolve_first_order(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField],
    ) -> dict[str, ScalarField]:
        """Compute virtual momenta for first-order components.

        The RHS of each first-order equation becomes a "virtual momentum"
        that second-order equations can reference via ``pi_N`` or
        ``first_derivative_t``.  This runs after constraints are solved
        so that first-order equations see updated constraint fields.

        Returns
        -------
        dict[str, ScalarField]
            Updated virtual momenta dict including first-order entries.

        Raises
        ------
        TypeError
            If state is not a FieldCollection.
        """
        # Type narrowing: RHS computation requires FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "First-order evolution requires FieldCollection state"
            raise TypeError(msg)

        for i, eq in enumerate(self.spec.equations):
            if eq.time_derivative_order == 1:
                virtual_momenta[eq.field_name] = self._compute_rhs_for_component(
                    i, state, bc, t
                )
        return virtual_momenta

    def _evolve_second_order(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField],
        rates: list[ScalarField | None],
    ) -> list[ScalarField | None]:
        """Compute rates for all components from the virtual momenta.

        - Second-order: d/dt field = momentum, d/dt momentum = RHS
        - First-order: d/dt field = virtual momentum (already computed)
        - Constraint: d/dt field = 0

        Returns
        -------
        list[ScalarField | None]
            Populated rates array.

        Raises
        ------
        TypeError
            If momentum field is not a ScalarField or state is not FieldCollection.
        """
        # Type narrowing: second-order evolution requires FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "Second-order evolution requires FieldCollection state"
            raise TypeError(msg)

        grid = state.grid
        for i, eq in enumerate(self.spec.equations):
            field_slot = self._field_slot_map[eq.field_name]

            if eq.time_derivative_order >= 2:  # noqa: PLR2004
                momentum_slot = self._momentum_slot_map[eq.field_name]
                momentum = state[momentum_slot]
                if not isinstance(momentum, ScalarField):
                    msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                    raise TypeError(msg)
                rates[field_slot] = momentum.copy()
                rates[momentum_slot] = self._compute_rhs_for_component(
                    i, state, bc, t, virtual_momenta
                )
            elif eq.time_derivative_order == 1:
                rates[field_slot] = virtual_momenta[eq.field_name]
            else:
                rates[field_slot] = ScalarField(grid, data=0.0)

        return rates

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
        TypeError
            If ``state`` is not a ``FieldCollection``.
        ValueError
            If the state size or grid dimension does not match the spec.
        """
        if not isinstance(state, FieldCollection):
            msg = f"Expected FieldCollection, got {type(state).__name__}"
            raise TypeError(msg)
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

        # B5: Cache boundary conditions (same grid always produces same BCs)
        grid = state.grid
        grid_id = id(grid)
        if self._cached_bc is None or self._cached_grid_id != grid_id:
            self._cached_bc = infer_bc_from_grid(grid)
            self._cached_grid_id = grid_id
        bc = self._cached_bc

        # Pass 1a: Solve constraints
        state, virtual_momenta = self._evolve_constraints(state, bc, t)

        # Pass 1b: Compute first-order virtual momenta
        virtual_momenta = self._evolve_first_order(state, bc, t, virtual_momenta)

        # Pass 2: Build the full rates array
        rates: list[ScalarField | None] = [None] * expected_fields
        rates = self._evolve_second_order(state, bc, t, virtual_momenta, rates)

        return FieldCollection(rates)  # type: ignore[arg-type]

    def check_stability(self, dt: float, grid: GridBase) -> list[str]:
        """Check CFL / stability conditions for explicit time-stepping.

        Estimates maximum wave speeds and diffusivities from the equation
        coefficients and checks whether ``dt`` satisfies the CFL condition.

        Parameters
        ----------
        dt : float
            Proposed time step.
        grid : GridBase
            The spatial grid (used for cell spacing).

        Returns
        -------
        list[str]
            Warning messages for any violated conditions. Empty if stable.
        """
        warnings: list[str] = []
        dx_min = min(grid.discretization)

        for eq in self.spec.equations:
            if eq.time_derivative_order < 2:  # noqa: PLR2004
                continue

            # Estimate max wave speed from laplacian coefficients
            max_laplacian_coeff = 0.0
            max_diffusive_coeff = 0.0
            for term in eq.rhs_terms:
                coeff_abs = abs(term.coefficient)
                if (term.operator == "laplacian" and term.field == eq.field_name) or term.operator.startswith("laplacian_"):
                    max_laplacian_coeff = max(max_laplacian_coeff, coeff_abs)
                elif term.operator == "biharmonic":
                    max_diffusive_coeff = max(max_diffusive_coeff, coeff_abs)

            # CFL for wave equation: dt < dx / c where c = sqrt(laplacian_coeff)
            if max_laplacian_coeff > 0:
                c_max = math.sqrt(max_laplacian_coeff)
                cfl_limit = dx_min / c_max
                if dt > cfl_limit:
                    warnings.append(
                        f"CFL violated for {eq.field_name}: "
                        f"dt={dt:.3e} > dx/c={cfl_limit:.3e} "
                        f"(c={c_max:.3e}, dx={dx_min:.3e})"
                    )

            # Stability for biharmonic: dt < dx^4 / (2 * coeff)
            if max_diffusive_coeff > 0:
                biharm_limit = dx_min**4 / (2 * max_diffusive_coeff)
                if dt > biharm_limit:
                    warnings.append(
                        f"Biharmonic stability violated for {eq.field_name}: "
                        f"dt={dt:.3e} > dx^4/(2D)={biharm_limit:.3e}"
                    )

        return warnings

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
    *,
    constraint_eps: float = 1e-14,
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
    constraint_eps : float
        Tolerance for the constraint solver's Laplacian coefficient check.
        Default ``1e-14``.

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
    # Auto-load metadata parameter defaults when none explicitly provided
    if parameters is None:
        meta_params = spec.metadata.get("parameters")
        if isinstance(meta_params, dict):
            parameters = {k: float(v) for k, v in meta_params.items() if isinstance(v, (int, float))}
    return PDEFromSpec(spec, parameters=parameters, constraint_eps=constraint_eps)


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
