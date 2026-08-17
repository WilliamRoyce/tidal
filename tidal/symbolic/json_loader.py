"""Load and validate JSON equation specifications from Mathematica/xAct export.

This module provides the data structures and parsing logic for loading
field equations that were derived symbolically from Lagrangians.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import math
import re
import warnings
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

from tidal.symbolic._eval_utils import evaluate_coefficient
from tidal.symbolic._kinetic_eval import (
    KineticEvalError,
    evaluate_with_substitutions,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from tidal.solver.operators import SideBCSpec

# --- Constants & operator validation ---

#: Canonical spatial axis letters, ordered by dimension index.
#: Supports up to 6 spatial dimensions (sufficient for all foreseeable physics).
AXIS_LETTERS: tuple[str, ...] = ("x", "y", "z", "w", "v", "u")

#: Character class matching all known axis letters (for regex construction).
_AXIS_RE_CLASS = "[" + "".join(AXIS_LETTERS) + "]"

#: Pattern matching a coordinate call like ``x[]``, ``y[]`` in symbolic expressions.
#: Used to auto-detect position-dependent coefficients in hamiltonian_terms that
#: were exported before the ``coordinate_dependent`` field was added.
_COORD_CALL_RE = re.compile(r"\b[xyzwvut]\s*\[\s*\]")

logger = logging.getLogger(__name__)

#: Set of static operators supported by the pipeline.
#: Validated at JSON load time to catch typos early.
_STATIC_OPERATORS: frozenset[str] = frozenset(
    {
        "identity",
        "laplacian",
        "laplacian_x",
        "laplacian_y",
        "laplacian_z",
        "gradient_x",
        "gradient_y",
        "gradient_z",
        "cross_derivative_xy",
        "cross_derivative_xz",
        "cross_derivative_yz",
        "first_derivative_t",
        "biharmonic",
    },
)

#: Pattern for generic single-axis Nth-order derivatives: derivative_3_x, derivative_5_y, etc.
_GENERIC_SINGLE_AXIS_RE = re.compile(r"^derivative_(\d+)_(" + _AXIS_RE_CLASS + r")$")

#: Pattern for generic multi-axis derivatives: derivative_2x_1y, derivative_3x_2z, etc.
_GENERIC_MULTI_AXIS_RE = re.compile(
    r"^derivative_(\d+" + _AXIS_RE_CLASS + r"(?:_\d+" + _AXIS_RE_CLASS + r")*)$",
)


#: Mutable set of user-registered custom operators.
#: Currently unused; reserved as an extensibility point for future
#: operator registration by downstream code.
_CUSTOM_OPERATORS: set[str] = set()


def is_known_operator(name: str) -> bool:
    """Check whether an operator name is recognized.

    Accepts static operators (identity, laplacian, gradient_x, ...),
    user-registered custom operators (via ``register_operator``),
    dynamic patterns for generic Nth-order derivatives
    (derivative_3_x, derivative_5_y, derivative_2x_1y, ...),
    mixed time-space derivative operators (mixed_T2_S2x, mixed_T_S1x, ...),
    and pure higher-order time operators on RHS (d2_t, d3_t, d4_t).
    """
    if name in _STATIC_OPERATORS or name in _CUSTOM_OPERATORS:
        return True
    if _GENERIC_SINGLE_AXIS_RE.match(name) or _GENERIC_MULTI_AXIS_RE.match(name):
        return True
    # Mixed time-space operators:
    #   Wolfram format: mixed_T[n]_S[n][axis]_... (e.g., mixed_T2_S2x)
    #   Numeric format: mixed_[time]_[s1]_[s2]_... (e.g., mixed_1_0_0_1)
    if name.startswith("mixed_T"):
        return True
    if re.match(r"^mixed_\d+(_\d+)+$", name):
        return True
    # Pure time derivative operators on RHS: d2_t, d3_t, d4_t
    return bool(re.match(r"^d\d+_t$", name))


# --- LHS structure ---


@dataclass(frozen=True)
class LHSStructure:
    """Structure describing the left-hand side of a PDE.

    Supports different PDE types:
    - Elliptic (time_order=0): ∇²φ = f (Poisson, Laplace)
    - Parabolic (time_order=1): ∂_t φ = ... (heat, diffusion)
    - Hyperbolic (time_order=2): ∂²_t φ = ... (wave)
    - Higher order: ∂^n_t φ = ...

    Attributes
    ----------
    expression : str
        String representation (e.g., "d2_t(phi_0)", "d_t(phi)", "phi").
    time_order : int
        Order of time derivative on LHS (0, 1, 2, or higher).
    space_order : int
        Order of space derivative on LHS (usually 0).
    """

    expression: str
    time_order: int
    space_order: int = 0
    kinetic_coefficient_symbolic: str | None = None
    """Symbolic expression for the kinetic coefficient when ExportJSON left RHS unnormalized.

    When non-None, the RHS terms are stored WITHOUT the kinetic coefficient divisor.
    Python must multiply RHS terms by ``1/kinetic_coefficient`` at runtime.
    If this evaluates to zero for the given parameters, the field is a constraint
    (the kinetic term vanishes and the EOM becomes algebraic).

    Example: for a torsion theory with ``L ⊃ ½ξ(∂T)²``, the torsion EOM kinetic
    coefficient is ``xi``.  At ``xi=0`` torsion becomes algebraically constrained.
    """

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LHSStructure:
        """Create LHSStructure from structured JSON data.

        Expected format: {"expression": "...", "order": {"time": N, "space": 0}}
        Optional: ``"kinetic_coefficient_symbolic"`` — see class docstring.

        Parameters
        ----------
        data : Mapping[str, Any]
            The structured LHS data from JSON.

        Returns
        -------
        LHSStructure
            Parsed LHS structure.

        Raises
        ------
        ValueError
            If ``order`` dict does not contain a ``'time'`` key.
        """
        expression = str(data.get("expression", ""))
        order = data.get("order", {})
        if "time" not in order:
            msg = (
                "lhs.order must specify 'time' "
                "(0 for constraint, 1 for parabolic, >=2 for hyperbolic)"
            )
            raise ValueError(msg)
        time_order = int(order["time"])
        space_order = int(order.get("space", 0))
        kc_sym = data.get("kinetic_coefficient_symbolic") or None
        return cls(
            expression=expression,
            time_order=time_order,
            space_order=space_order,
            kinetic_coefficient_symbolic=kc_sym,
        )


# --- RHS operator terms ---


@dataclass(frozen=True)
class OperatorTerm:
    """A single term in the RHS of a field equation.

    Represents: coefficient * operator(field)

    Attributes
    ----------
    coefficient : float
        Numeric coefficient for this term.
    operator : str
        Name of the differential operator ("laplacian", "identity", "gradient_x", etc.)
    field : str
        Name of the field this operator acts on.
    coefficient_symbolic : str | None
        Optional symbolic name for the coefficient (e.g., "m2", "-kappa").
        When present, the coefficient can be overridden at runtime by passing
        a parameters dict to the PDE constructor.
    time_dependent : bool
        Whether the coefficient depends on time. For curved spacetime, terms
        like -2H∂_t φ (Hubble friction) have time-dependent coefficients when
        the conformal factor Ω(t) varies with time. Default False for flat spacetime.
    coordinate_dependent : tuple[str, ...]
        Coordinate names the coefficient depends on (e.g., ("x", "y") for
        position-dependent coefficients on curved spatial surfaces, or ("t",)
        for time-dependent). Empty tuple for constant coefficients.
    order_in_eps : int
        Order of this term in the small parameters configured by the theory's
        ``[perturbation]`` section. 0 for the base (unperturbed) theory; 1 for
        a first-order correction, etc. Set by the Wolfram ``ComputeOrderInEps``
        helper and consumed by :class:`EquationSystem.filter_by_order` to
        separate base from correction terms for the iterative perturbative
        solver (v6 plan, Stage 2). Defaults to 0 for backward compatibility
        with JSON files predating the field.
    """

    coefficient: float
    operator: str
    field: str
    coefficient_symbolic: str | None = None
    time_dependent: bool = False
    coordinate_dependent: tuple[str, ...] = ()
    order_in_eps: int = 0

    @property
    def position_dependent(self) -> bool:
        """Whether the coefficient depends on spatial coordinates.

        Returns ``True`` when ``coordinate_dependent`` is non-empty (explicit
        declaration), or when ``coefficient_symbolic`` contains a spatial
        coordinate call pattern such as ``x[]`` or ``y[]`` (auto-detection for
        JSON exports that predate the ``coordinate_dependent`` field).
        Time-only dependence (``t[]``) returns ``False``.
        """
        if self.coordinate_dependent:
            return bool(set(self.coordinate_dependent) - {"t"})
        if self.coefficient_symbolic is not None:
            # findall returns strings like "x[]", "t[]"; exclude time coord
            return any(
                m[0] != "t" for m in _COORD_CALL_RE.findall(self.coefficient_symbolic)
            )
        return False

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> OperatorTerm:
        """Create an OperatorTerm from a dictionary.

        Raises
        ------
        ValueError
            If required keys are missing or operator is unknown.
        """
        required_keys = {"coefficient", "operator", "field"}
        missing = required_keys - set(data.keys())
        if missing:
            msg = (
                f"RHS term missing required keys: {sorted(missing)}. Got: {dict(data)}"
            )
            raise ValueError(msg)

        operator = str(data["operator"])
        if not is_known_operator(operator):
            msg = (
                f"Unknown operator '{operator}'. "
                f"Known static operators: {sorted(_STATIC_OPERATORS)}. "
                f"Dynamic patterns: derivative_N_x, derivative_Nx_My (N,M=integers, x,y,z=axes)."
            )
            raise ValueError(msg)

        # Check for nonlinear coefficient warning from Wolfram export
        if data.get("warning") == "nonlinear_coefficient":
            import warnings  # noqa: PLC0415

            warnings.warn(
                f"Term '{operator}' on field '{data['field']}' has a "
                f"field-dependent (nonlinear) coefficient. The linear PDE "
                f"solver will treat it as a constant, producing wrong physics.",
                stacklevel=2,
            )

        return cls(
            coefficient=float(data["coefficient"]),
            operator=operator,
            field=str(data["field"]),
            coefficient_symbolic=data.get("coefficient_symbolic"),
            time_dependent=bool(data.get("time_dependent", False)),
            coordinate_dependent=tuple(data.get("coordinate_dependent", ())),
            order_in_eps=int(data.get("order_in_eps", 0)),
        )


# --- Boundary conditions ---

_VALID_BC_TYPES: frozenset[str] = frozenset(
    {"periodic", "dirichlet", "neumann", "robin"},
)


@dataclass(frozen=True)
class BoundaryCondition:
    """Boundary condition for one spatial axis.

    Attributes
    ----------
    type : str
        One of "periodic", "dirichlet", "neumann", or "robin".
    value : float | None
        Fixed value for Dirichlet BCs, or Robin beta.
    derivative : float | None
        Fixed normal derivative for Neumann BCs.
    gamma : float | None
        Robin coefficient gamma in d_n f + gamma*f = beta.
    """

    type: str
    value: float | None = None
    derivative: float | None = None
    gamma: float | None = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BoundaryCondition:
        """Create a BoundaryCondition from a dictionary.

        Raises
        ------
        ValueError
            If the BC type is not recognized.
        """
        bc_type = str(data["type"])
        if bc_type not in _VALID_BC_TYPES:
            msg = (
                f"Unknown BC type: {bc_type!r}. Valid types: {sorted(_VALID_BC_TYPES)}"
            )
            raise ValueError(msg)
        # Warn on irrelevant fields
        irrelevant_fields: dict[str, set[str]] = {
            "periodic": {"value", "derivative", "gamma"},
            "dirichlet": {"derivative", "gamma"},
            "neumann": {"gamma"},
        }
        extra = {
            k for k in irrelevant_fields.get(bc_type, set()) if data.get(k) is not None
        }
        if extra:
            logger.warning(
                "%s BC ignores field(s): %s",
                bc_type,
                ", ".join(sorted(extra)),
            )
        return cls(
            type=bc_type,
            value=data.get("value"),
            derivative=data.get("derivative"),
            gamma=data.get("gamma"),
        )

    def to_side_bc(self) -> SideBCSpec:
        """Convert to a ``SideBCSpec`` for the operator layer.

        Raises
        ------
        ValueError
            If the BC type is "periodic" (not representable as a side BC).
        """
        from tidal.solver.operators import SideBCSpec  # noqa: PLC0415

        if self.type == "periodic":
            msg = "Cannot convert periodic BC to SideBCSpec"
            raise ValueError(msg)
        return SideBCSpec(
            kind=self.type,
            value=self.value if self.value is not None else 0.0,
            derivative=self.derivative if self.derivative is not None else 0.0,
            gamma=self.gamma if self.gamma is not None else 0.0,
        )


# --- Constraint solver configuration ---

_VALID_CONSTRAINT_METHODS = frozenset({"auto", "fft", "matrix", "poisson"})


@dataclass(frozen=True)
class ConstraintSolverConfig:
    """Configuration for elliptic constraint solving.

    When ``enabled`` is True, the constraint equation is solved at each
    timestep rather than remaining frozen at its initial value.

    Solver Methods
    --------------
    - **"auto"** (default): Automatically selects the best solver.
      Uses FFT for fully periodic grids (O(N log N)), sparse matrix
      for non-periodic grids (O(N) via LU). Handles all constraint
      types: Poisson, Helmholtz, algebraic, anisotropic, etc.

    - **"fft"**: Force FFT solver. Requires fully periodic grid.

    - **"matrix"**: Force sparse matrix solver. Works with any BCs.

    - **"poisson"**: Original py-pde Poisson solver. Requires exactly
      one ``laplacian(self_field)`` term with no other self-referencing
      operators. Will warn if non-laplacian self-terms are present.
      Backward compatible.

    Coupled Constraint Parameters
    ----------------------------
    When multiple constraints reference each other's fields, the solver
    iterates using Gauss-Seidel until convergence or ``max_iterations``.
    For fully periodic grids, coupled constraints are solved exactly
    via Fourier-space block solve (no iteration needed).

    Attributes
    ----------
    enabled : bool
        Whether to solve the constraint elliptically. Default False
        preserves existing frozen-constraint behavior.
    method : str
        Solver method: "auto", "fft", "matrix", or "poisson".
    boundary_conditions : dict[str, BoundaryCondition]
        Per-axis boundary conditions (e.g., ``{"x": ..., "y": ...}``).
    max_iterations : int
        Maximum Gauss-Seidel iterations for coupled constraints.
        Must be >= 1.
    tolerance : float
        Convergence threshold for coupled constraint iteration.
        Iteration stops when max|field_new - field_old| < tolerance
        (scaled by field magnitude for robustness). Must be > 0.
    """

    enabled: bool = False
    method: str = "auto"
    boundary_conditions: dict[str, BoundaryCondition] = dataclass_field(
        default_factory=lambda: {},  # noqa: PIE807  # type: dict[str, BoundaryCondition]
    )
    max_iterations: int = 20
    tolerance: float = 1e-8

    @classmethod
    def from_dict(cls, data: Mapping[str, Any] | None) -> ConstraintSolverConfig:
        """Create from a dictionary or return default (disabled).

        Parameters
        ----------
        data : Mapping[str, Any] | None
            Parsed ``constraint_solver`` block from JSON, or None.

        Returns
        -------
        ConstraintSolverConfig
            Configuration instance.

        Raises
        ------
        ValueError
            If ``method`` is not one of the recognized solver methods.
        """
        if data is None:
            return cls()

        enabled = bool(data.get("enabled", False))
        method = str(data.get("method", "auto"))

        if method not in _VALID_CONSTRAINT_METHODS:
            msg = (
                f"Unknown constraint solver method '{method}'. "
                f"Must be one of: {', '.join(sorted(_VALID_CONSTRAINT_METHODS))}."
            )
            raise ValueError(msg)

        bc_data = data.get("boundary_conditions", {})
        boundary_conditions = {
            axis: BoundaryCondition.from_dict(bc_dict)
            for axis, bc_dict in bc_data.items()
        }

        max_iterations = int(data.get("max_iterations", 20))
        if max_iterations < 1:
            msg = f"max_iterations must be >= 1, got {max_iterations}"
            raise ValueError(msg)

        tolerance = float(data.get("tolerance", 1e-8))
        if tolerance <= 0:
            msg = f"tolerance must be > 0, got {tolerance}"
            raise ValueError(msg)

        return cls(
            enabled=enabled,
            method=method,
            boundary_conditions=boundary_conditions,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )


# === Phase K: Canonical Momentum Structures ===


@dataclass(frozen=True)
class HamiltonianFactor:
    """One factor in a quadratic Hamiltonian term.

    Represents a field (or its time/spatial derivative) that appears as a
    multiplicative factor in a term of the component-form Hamiltonian density.

    Attributes
    ----------
    field : str
        Component field name (e.g., "A_1", "phi_0").
    operator : str
        Differential operator: "identity" (bare field), "time_derivative"
        (∂_t field, maps to canonical momentum at evaluation), or a spatial
        operator ("gradient_x", "laplacian", etc.).
    """

    field: str
    operator: str

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HamiltonianFactor:
        """Parse from JSON dict."""
        return cls(field=str(data["field"]), operator=str(data["operator"]))


def _base_field(field_name: str) -> str:
    """Strip velocity prefix ``v_`` to get the base field name.

    E-L velocity form uses ``v_{field}`` for velocities. This extracts the
    underlying field name: ``"v_phi_0"`` → ``"phi_0"``, ``"phi_0"`` → ``"phi_0"``.
    """
    if field_name.startswith("v_"):
        return field_name[2:]
    return field_name


@dataclass(frozen=True)
class HamiltonianTerm:
    """A single quadratic term in the Hamiltonian density.

    H = Σ coefficient * factor_a * factor_b

    Attributes
    ----------
    coefficient : float
        Numeric coefficient.
    factor_a : HamiltonianFactor
        First field factor.
    factor_b : HamiltonianFactor
        Second field factor (may equal factor_a for squared terms).
    coefficient_symbolic : str | None
        Symbolic coefficient expression (for parameter override).
    coordinate_dependent : tuple[str, ...]
        Spatial axes the coefficient depends on (e.g. ``("x", "y")``).
        When non-empty, the coefficient must be evaluated on the grid.
        Older JSON exports omit this field; auto-detection via
        ``position_dependent`` covers those cases.
    term_class : str
        Classification: ``"self"`` (both factors reference the same base
        field) or ``"interaction"`` (cross-field coupling). Defaults to
        ``"unknown"`` for older JSONs; use ``is_self_energy`` property
        which auto-classifies by comparing factor field names.
    order_in_eps : int
        Perturbative order of this term: total exponent in declared
        ``small_parameters``.  Computed by Wolfram's ``ComputeOrderInEps``
        and emitted explicitly by ``ParseSingleHamiltonianTerm``.  Defaults
        to 0 for non-perturbative theories and direct hand-construction;
        legacy JSONs without the explicit field are tagged via heuristic
        (``1 if coefficient_symbolic else 0``) and should be re-derived.
    """

    coefficient: float
    factor_a: HamiltonianFactor
    factor_b: HamiltonianFactor
    coefficient_symbolic: str | None = None
    coordinate_dependent: tuple[str, ...] = ()
    term_class: str = "unknown"  # "self" or "interaction"
    # Perturbative order computed by Wolfram via ComputeOrderInEps; matches the
    # equation-side OperatorTerm.order_in_eps convention.  Default 0 covers
    # non-perturbative theories and direct hand-construction in tests.
    order_in_eps: int = 0

    @property
    def position_dependent(self) -> bool:
        """True if the coefficient is a function of spatial coordinates.

        Returns ``True`` when ``coordinate_dependent`` is non-empty (explicit
        declaration), or when ``coefficient_symbolic`` contains a coordinate
        call pattern such as ``x[]`` or ``y[]`` (auto-detection for JSON
        exports that predate the ``coordinate_dependent`` field).
        """
        if self.coordinate_dependent:
            return True
        if self.coefficient_symbolic is not None:
            return bool(_COORD_CALL_RE.search(self.coefficient_symbolic))
        return False

    @property
    def is_self_energy(self) -> bool:
        """True if both factors reference the same base field (self-energy).

        Uses ``term_class`` when available (from Wolfram export), otherwise
        classifies by comparing base field names (stripping ``v_`` prefix).
        """
        if self.term_class != "unknown":
            return self.term_class == "self"
        return _base_field(self.factor_a.field) == _base_field(self.factor_b.field)

    @property
    def base_field_a(self) -> str:
        """Base field name for factor_a (strips ``v_`` velocity prefix)."""
        return _base_field(self.factor_a.field)

    @property
    def base_field_b(self) -> str:
        """Base field name for factor_b (strips ``v_`` velocity prefix)."""
        return _base_field(self.factor_b.field)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HamiltonianTerm:
        """Parse from JSON dict.

        Reads the explicit ``order_in_eps`` field emitted by Wolfram's
        ``ComputeOrderInEps`` (single source of truth, matching the equation
        side at ``OperatorTerm.order_in_eps``).  For legacy JSONs without the
        explicit field, falls back to the heuristic ``1 if coefficient_symbolic
        is not None else 0`` — correct for current EH but coincidental in
        general.  The strict consistency check lives in ``tidal validate``.
        """
        coeff_sym = data.get("coefficient_symbolic")
        if "order_in_eps" in data:
            order = int(data["order_in_eps"])
        else:
            # Legacy JSON: infer via heuristic.  Coincidentally correct when
            # every coefficient_symbolic contains a small parameter (true for
            # Euler-Heisenberg) but not robust in general.  Re-derive the JSON
            # via `uv run tidal derive <toml>` to obtain Wolfram-emitted
            # rigorous tags.
            order = 1 if coeff_sym is not None else 0
        return cls(
            coefficient=float(data["coefficient"]),
            factor_a=HamiltonianFactor.from_dict(data["factor_a"]),
            factor_b=HamiltonianFactor.from_dict(data["factor_b"]),
            coefficient_symbolic=coeff_sym,
            coordinate_dependent=tuple(data.get("coordinate_dependent", [])),
            term_class=str(data.get("term_class", "unknown")),
            order_in_eps=order,
        )


@dataclass(frozen=True)
class CanonicalStructure:
    """Canonical structure: Hamiltonian terms for energy measurement.

    The Hamiltonian's bilinear terms are used by energy measurement to
    compute H(q, v).  The E-L equations are stored in the ``equations``
    array directly.

    Attributes
    ----------
    hamiltonian_terms : tuple[HamiltonianTerm, ...]
        Quadratic terms in the component-form Hamiltonian density.
        Used by energy measurement to compute H(q, v).
    volume_element : str or None
        Symbolic expression for ``sqrt|det(g_spatial)|``, the spatial
        volume element.  ``None`` for flat (Minkowski) spacetimes where
        the volume element is 1.  Used by energy measurement to weight
        the Hamiltonian density before spatial integration.
    """

    hamiltonian_terms: tuple[HamiltonianTerm, ...]
    volume_element: str | None = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CanonicalStructure:
        """Parse from JSON ``canonical`` section."""
        h_terms = tuple(HamiltonianTerm.from_dict(t) for t in data["hamiltonian_terms"])
        vol_elem = data.get("volume_element")  # None for flat spacetimes
        return cls(
            hamiltonian_terms=h_terms,
            volume_element=vol_elem,
        )


# --- Component equations ---


@dataclass(frozen=True)
class ComponentEquation:
    """Equation of motion for a single field component.

    For a wave-type equation:
        d^2/dt^2 field = sum of OperatorTerms

    Attributes
    ----------
    field_name : str
        Name of the field component (e.g., "A_0", "phi").
    field_index : int
        Index of this component in the field array.
    time_derivative_order : int
        Order of the time derivative on the LHS (2 for wave equations).
    rhs_terms : tuple[OperatorTerm, ...]
        Terms on the RHS of the equation.
    constraint_solver : ConstraintSolverConfig
        Configuration for elliptic constraint solving. Only meaningful
        when ``time_derivative_order == 0``.
    """

    field_name: str
    field_index: int
    time_derivative_order: int
    rhs_terms: tuple[OperatorTerm, ...]
    constraint_solver: ConstraintSolverConfig = dataclass_field(
        default_factory=ConstraintSolverConfig,
    )
    kinetic_coefficient_symbolic: str | None = None
    """Symbolic kinetic coefficient when ExportJSON left RHS unnormalized.

    Non-None only for equations derived with a parameter-based kinetic coefficient
    (e.g., ``xi`` in the dark photon torsion model).  The RHS terms are stored
    WITHOUT the 1/kinetic_coefficient divisor.  Use
    ``normalize_kinetic_coefficients(spec, params)`` to apply the normalization
    before passing the spec to a solver.  If this evaluates to zero at the given
    parameters, the field becomes a constraint (kinetic term vanishes).
    """

    def __post_init__(self) -> None:
        """Validate constraint_solver is only enabled for time_order=0.

        Raises
        ------
        ValueError
            If constraint_solver.enabled is True but time_derivative_order is not 0.
        """
        if self.constraint_solver.enabled and self.time_derivative_order != 0:
            msg = (
                f"constraint_solver.enabled=true is only valid for time_order=0 "
                f"(constraint equations), but {self.field_name} has "
                f"time_order={self.time_derivative_order}"
            )
            raise ValueError(msg)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        fields_lookup: dict[str, int],
    ) -> ComponentEquation:
        """Create a ComponentEquation from a dictionary.

        Raises
        ------
        ValueError
            If the RHS type is not "linear_combination", or if
            constraint_solver is enabled for a non-constraint equation.
        """
        field_name = str(data["field"])

        # Parse LHS to determine time derivative order
        # Requires structured format: {"expression": "...", "order": {"time": N, "space": 0}}
        lhs_data = data["lhs"]  # Required field
        lhs_structure = LHSStructure.from_dict(lhs_data)
        time_derivative_order = lhs_structure.time_order

        # Parse RHS terms
        rhs_data = data["rhs"]
        if rhs_data["type"] != "linear_combination":
            msg = f"Unsupported RHS type: {rhs_data['type']}"
            raise ValueError(msg)

        rhs_terms = tuple(
            OperatorTerm.from_dict(term_data) for term_data in rhs_data["terms"]
        )

        # Validate field_name exists in fields_lookup - no silent fallback to 0
        if field_name not in fields_lookup:
            valid_fields = list(fields_lookup.keys())
            msg = (
                f"Unknown field '{field_name}' in equation. "
                f"Valid fields are: {valid_fields}"
            )
            raise ValueError(msg)

        # Parse constraint solver config
        constraint_solver = ConstraintSolverConfig.from_dict(
            data.get("constraint_solver"),
        )

        return cls(
            field_name=field_name,
            field_index=fields_lookup[field_name],
            time_derivative_order=time_derivative_order,
            rhs_terms=rhs_terms,
            constraint_solver=constraint_solver,
            kinetic_coefficient_symbolic=lhs_structure.kinetic_coefficient_symbolic,
        )


# --- Symbolic coefficient resolution ---


def _accumulate_symbolic(existing: str | None, addition: str) -> str:
    """Sum a symbolic matrix contribution onto whatever is already there.

    A component can carry several ``identity`` terms acting on the same field —
    a base mass plus a background correction, say — and the matrix entry is
    their **sum**.  The numeric matrices have always accumulated; the symbolic
    ones used to overwrite, silently keeping only the last term (GH #403).

    Parameters
    ----------
    existing : str | None
        Expression accumulated so far, or ``None`` for the first contribution.
    addition : str
        The term's ``coefficient_symbolic``.

    Returns
    -------
    str
        ``addition`` when this is the first contribution, else ``"(a) + (b)"``.
        Operands are parenthesised so that summing e.g. ``-1 + 2*B0^2*rho``
        with ``-xi`` cannot re-associate into something else.
    """
    if existing is None:
        return addition
    return f"({existing}) + ({addition})"


class _UnresolvableCoefficientError(Exception):
    """Internal signal that a coefficient could not be reduced to a number."""


def _resolve_symbolic_coeff(sym: str, parameters: Mapping[str, float]) -> float | None:
    """Resolve a symbolic coefficient string with parameter values.

    Handles simple names (``"m2"``), negated names (``"-m2"``), and
    compound expressions (``"-2*m2"``).  Returns *None* if the expression
    cannot be evaluated with the given parameters, so the caller can fall
    back to the raw numeric coefficient.

    Two results are preserved deliberately, because callers depend on the
    ``None`` fallback rather than on an exception:

    * a division by zero (``"m2/0"``) returns ``None``, though the underlying
      evaluator raises;
    * a non-finite result (``"1e309"``) returns ``None``, since mass and
      coupling entries must be finite reals.

    Raises
    ------
    _UnresolvableCoefficientError
        Never propagates; raised internally to funnel an unresolved free
        parameter into the same ``None`` fallback as an evaluation failure.
    """
    # Fast path: negated parameter name
    if sym.startswith("-") and sym[1:] in parameters:
        return -parameters[sym[1:]]
    # Fast path: direct parameter name
    if sym in parameters:
        return parameters[sym]
    # Compound expression: e.g., "-2*m2", "3*lambda".
    #
    # Evaluated through the restricted-AST evaluator rather than eval(): these
    # strings come from JSON spec files, and `{"__builtins__": {}}` narrows the
    # namespace without restricting the grammar (GH #406).
    try:
        resolved = evaluate_with_substitutions(sym, dict(parameters))
        if resolved is None:  # a free parameter survived
            raise _UnresolvableCoefficientError  # noqa: TRY301
        value = float(resolved)
    except (KineticEvalError, _UnresolvableCoefficientError, TypeError, ValueError):
        logger.debug(
            "Could not resolve symbolic coefficient '%s' with parameters %s; "
            "matrix entry will use raw numeric coefficient",
            sym,
            sorted(parameters.keys()),
        )
        return None
    # Mass/coupling entries must be finite real numbers
    if not math.isfinite(value):
        return None
    return value


# --- Equation system ---


@dataclass(frozen=True)
class EquationSystem:
    """Complete system of field equations derived from a Lagrangian.

    Attributes
    ----------
    n_components : int
        Number of field components.
    dimension : int
        Spacetime dimension (e.g., 2 for 1+1D).
    spatial_dimension : int
        Number of spatial dimensions (dimension - 1).
    component_names : tuple[str, ...]
        Names of field components in order.
    equations : tuple[ComponentEquation, ...]
        Equations for each component.
    mass_matrix : tuple[tuple[float, ...], ...]
        Mass matrix M^2_ij for coupled systems.
    coupling_matrix : tuple[tuple[float, ...], ...]
        Coupling matrix for field interactions.
    metadata : dict[str, Any]
        Additional metadata (source, gauge, etc.)
    coordinates : tuple[str, ...]
        Coordinate names from JSON spacetime.coordinates (e.g., ("t", "x", "y")).
        Defaults to empty tuple; use ``effective_coordinates`` for a guaranteed
        non-empty result that infers names from dimension when not set.
    canonical : CanonicalStructure | None
        Canonical momentum and Hamiltonian structure from Legendre transform.
        Present when the JSON spec includes a ``"canonical"`` section (generated
        by ``tidal derive`` for non-linearization theories). None for legacy
        specs or linearization theories.
    """

    n_components: int
    dimension: int
    spatial_dimension: int
    component_names: tuple[str, ...]
    equations: tuple[ComponentEquation, ...]
    mass_matrix: tuple[tuple[float, ...], ...]
    coupling_matrix: tuple[tuple[float, ...], ...]
    metadata: dict[str, Any]
    coordinates: tuple[str, ...] = ()
    signature: tuple[int, ...] = ()
    mass_matrix_symbolic: tuple[tuple[str | None, ...], ...] = ()
    coupling_matrix_symbolic: tuple[tuple[str | None, ...], ...] = ()
    canonical: CanonicalStructure | None = None

    def __post_init__(self) -> None:
        """Validate the equation system.

        Raises
        ------
        ValueError
            If validation checks fail (invalid n_components, mismatched lengths,
            or invalid matrix dimensions).
        """
        if self.n_components < 1:
            msg = "n_components must be at least 1"
            raise ValueError(msg)

        if len(self.component_names) != self.n_components:
            msg = f"component_names length {len(self.component_names)} != n_components {self.n_components}"
            raise ValueError(msg)

        if len(self.equations) != self.n_components:
            msg = f"equations length {len(self.equations)} != n_components {self.n_components}"
            raise ValueError(msg)

        self._validate_matrix_dimensions()

        # Validate field references in equation terms
        self._validate_field_references()

        # Verify mass/coupling matrix consistency with identity operator terms.
        # Warns (does not error) if manually constructed EquationSystem has
        # matrices that don't match the convention: matrix[i][j] = -(identity coeff).
        raw_params = self.metadata.get("parameters", {})
        check_params: dict[str, float] = {
            k: float(v) for k, v in raw_params.items() if isinstance(v, (int, float))
        }
        expected_mass, expected_coupling, _, _ = self._compute_matrices_from_terms(
            self.equations,
            self.component_names,
            parameters=check_params or None,
        )
        if (
            self.mass_matrix != expected_mass
            or self.coupling_matrix != expected_coupling
        ):
            import warnings  # noqa: PLC0415

            warnings.warn(
                "mass_matrix/coupling_matrix inconsistent with identity operator terms. "
                f"Expected mass={expected_mass}, coupling={expected_coupling}; "
                f"got mass={self.mass_matrix}, coupling={self.coupling_matrix}. "
                "Use EquationSystem.from_dict() for auto-computation.",
                stacklevel=2,
            )

    def _validate_matrix_dimensions(self) -> None:
        """Validate mass and coupling matrix dimensions match n_components.

        Raises
        ------
        ValueError
            If any matrix has wrong number of rows or columns.
        """
        for name, matrix in [
            ("mass_matrix", self.mass_matrix),
            ("coupling_matrix", self.coupling_matrix),
        ]:
            if len(matrix) != self.n_components:
                msg = f"{name} rows {len(matrix)} != n_components {self.n_components}"
                raise ValueError(msg)
            for i, row in enumerate(matrix):
                if len(row) != self.n_components:
                    msg = f"{name} row {i} length {len(row)} != n_components {self.n_components}"
                    raise ValueError(msg)

    def _validate_field_references(self) -> None:
        """Validate that all field references in equation terms are valid.

        Field references can be:
        - Regular field names (e.g., ``"A_0"``, ``"A_1"``, ``"phi"``)
        - Velocity names: ``"v_field_name"`` (e.g., ``"v_A_1"``)

        Raises
        ------
        ValueError
            If a field reference is invalid.
        """
        valid_fields = set(self.component_names)
        # Build valid velocity references: v_field_name
        valid_velocities = {f"v_{name}" for name in self.component_names}

        for eq in self.equations:
            for term in eq.rhs_terms:
                field_ref = term.field

                # Check regular field names first (a field named "v_0" is
                # a valid field, not a velocity reference)
                if field_ref in valid_fields:
                    continue
                if field_ref in valid_velocities:
                    continue
                if field_ref.startswith("v_"):
                    msg = (
                        f"Invalid velocity reference '{field_ref}' "
                        f"in equation for {eq.field_name}. "
                        f"Valid: {sorted(valid_velocities)}."
                    )
                    raise ValueError(msg)
                msg = (
                    f"Unknown field reference '{field_ref}' "
                    f"in equation for {eq.field_name}. "
                    f"Valid fields: {sorted(valid_fields)}."
                )
                raise ValueError(msg)

    @staticmethod
    def _compute_matrices_from_terms(
        equations: tuple[ComponentEquation, ...],
        component_names: tuple[str, ...],
        parameters: Mapping[str, float] | None = None,
    ) -> tuple[
        tuple[tuple[float, ...], ...],
        tuple[tuple[float, ...], ...],
        tuple[tuple[str | None, ...], ...],
        tuple[tuple[str | None, ...], ...],
    ]:
        """Extract mass and coupling matrices from identity operator terms.

        Scans each equation's RHS terms for ``identity`` operators acting on
        known field names (not velocity references like ``v_N``).

        **This is the only implementation.** The matrices are *derived* from
        ``equations`` and are never stored: ``ExportJSON.wl`` used to emit a
        ``coupling`` section computed the same way in Wolfram, and the loader
        stopped reading it in c407240d (2026-02-07) without the export stopping.
        Two copies of one calculation then drifted apart unnoticed for six
        months — both truncating multi-term coefficients, and disagreeing on
        sign (GH #403, #404). The export was removed rather than repaired,
        because derived data does not belong in the file.

        The numeric and symbolic matrices describe the same quantity under
        different conventions, which is deliberate and load-bearing:

        * ``mass``/``coupling`` (numeric) — ``matrix[i][j] = -(coefficient)`` of
          the ``identity(field_j)`` term in equation *i*. The negation makes
          mass² positive for the standard Lagrangian sign convention
          ``∂²_t φ = … - m² φ``.
        * ``mass_sym``/``coupling_sym`` (symbolic) — the term's
          ``coefficient_symbolic`` **verbatim, un-negated**, so it can be
          resolved at runtime against real parameter values.

        So an identity term ``B0^2/8`` gives ``mass = -0.125`` and
        ``mass_sym = 'B0^2/8'``. A caller applying the numeric convention to the
        symbolic matrix gets the sign wrong; negate it explicitly.

        Neither matrix is normalised by the LHS ``kinetic_coefficient_symbolic``.
        For an effective mass², divide by it; see
        :func:`tidal.symbolic.spec_query.effective_coefficient`, which is the
        one implementation of that (#237, #258, #302).

        Contributions **accumulate**: a component may carry several ``identity``
        terms acting on the same field, and the entry is their sum, symbolic
        entries included (GH #403).

        When *parameters* are provided and a term has ``coefficient_symbolic``,
        the symbolic expression is resolved with the parameter values to
        produce the correct numeric matrix entry. Without parameters, the
        raw numeric coefficient (a shape-factor like ±1.0) is used.

        Parameters
        ----------
        equations : tuple[ComponentEquation, ...]
            Parsed component equations.
        component_names : tuple[str, ...]
            Names of all field components.
        parameters : Mapping[str, float] | None
            Default parameter values for resolving symbolic coefficients.

        Returns
        -------
        tuple
            ``(mass_matrix, coupling_matrix, mass_symbolic, coupling_symbolic)``
            Numeric matrices as nested float tuples; symbolic matrices as
            nested ``str | None`` tuples (empty tuple if no symbolic entries).
        """
        n = len(component_names)
        mass: list[list[float]] = [[0.0] * n for _ in range(n)]
        coupling: list[list[float]] = [[0.0] * n for _ in range(n)]
        mass_sym: list[list[str | None]] = [[None] * n for _ in range(n)]
        coupling_sym: list[list[str | None]] = [[None] * n for _ in range(n)]
        name_to_idx = {name: i for i, name in enumerate(component_names)}
        has_symbolic = False

        for i, eq in enumerate(equations):
            for term in eq.rhs_terms:
                if term.operator == "identity" and term.field in name_to_idx:
                    j = name_to_idx[term.field]
                    # Prefer symbolic + params for numeric value; fall back
                    # to the raw numeric coefficient (shape factor).
                    effective_coeff = term.coefficient
                    if term.coefficient_symbolic is not None and parameters:
                        resolved = _resolve_symbolic_coeff(
                            term.coefficient_symbolic,
                            parameters,
                        )
                        if resolved is not None:
                            effective_coeff = resolved
                    neg_coeff = -effective_coeff
                    if i == j:
                        mass[i][j] += neg_coeff
                        if term.coefficient_symbolic is not None:
                            mass_sym[i][j] = _accumulate_symbolic(
                                mass_sym[i][j],
                                term.coefficient_symbolic,
                            )
                            has_symbolic = True
                    else:
                        coupling[i][j] += neg_coeff
                        if term.coefficient_symbolic is not None:
                            coupling_sym[i][j] = _accumulate_symbolic(
                                coupling_sym[i][j],
                                term.coefficient_symbolic,
                            )
                            has_symbolic = True

        mass_sym_t: tuple[tuple[str | None, ...], ...] = ()
        coupling_sym_t: tuple[tuple[str | None, ...], ...] = ()
        if has_symbolic:
            mass_sym_t = tuple(tuple(row) for row in mass_sym)
            coupling_sym_t = tuple(tuple(row) for row in coupling_sym)

        return (
            tuple(tuple(row) for row in mass),
            tuple(tuple(row) for row in coupling),
            mass_sym_t,
            coupling_sym_t,
        )

    @property
    def time_orders(self) -> tuple[int, ...]:
        """Per-component time derivative orders."""
        return tuple(eq.time_derivative_order for eq in self.equations)

    @property
    def state_size(self) -> int:
        """Total number of state fields.

        Second-order components contribute 2 slots (field + momentum).
        First-order and constraint components contribute 1 slot (field only).
        """
        return sum(2 if t >= 2 else 1 for t in self.time_orders)  # noqa: PLR2004

    @property
    def state_layout(self) -> tuple[tuple[str, str], ...]:
        """State vector layout as (field_name, slot_type) tuples.

        slot_type is "field" or "momentum". Second-order components produce
        two entries (field, momentum); first-order/constraint produce one (field).
        """
        layout: list[tuple[str, str]] = []
        for eq in self.equations:
            layout.append((eq.field_name, "field"))
            if eq.time_derivative_order >= 2:  # noqa: PLR2004
                layout.append((eq.field_name, "momentum"))
        return tuple(layout)

    @property
    def effective_coordinates(self) -> tuple[str, ...]:
        """Coordinate names, inferred from dimension if not set explicitly."""
        if self.coordinates:
            return self.coordinates
        return ("t", *("x", "y", "z")[: self.spatial_dimension])

    @property
    def spatial_coordinates(self) -> tuple[str, ...]:
        """Spatial coordinate names (all except first, which is time)."""
        return self.effective_coordinates[1:]

    @cached_property
    def has_constraint_velocity_terms(self) -> bool:
        """True if the Hamiltonian has kinetic coupling between constraint fields.

        Detects ``time_derivative(C_i) x time_derivative(C_j)`` terms where
        both C_i and C_j are constraint fields (time_derivative_order == 0).
        These indicate the naive Hamiltonian H = sum(pi * v) - L is not the
        correct conserved quantity for this theory (Dirac-Bergmann theory).

        See GitHub issue #178 for details.
        """
        if self.canonical is None:
            return False
        constraint_fields = frozenset(
            eq.field_name for eq in self.equations if eq.time_derivative_order == 0
        )
        if not constraint_fields:
            return False
        for term in self.canonical.hamiltonian_terms:
            a, b = term.factor_a, term.factor_b
            # Kinetic constraint term: time_derivative(C_i) x time_derivative(C_j)
            # where both fields are constraints.  These arise from π_c · v_c
            # in the Legendre transform (spurious).  Cross-terms like
            # time_derivative(C) x identity(D) come from -L and are valid.
            if (
                a.operator == "time_derivative"
                and b.operator == "time_derivative"
                and a.field in constraint_fields
                and b.field in constraint_fields
            ):
                return True
        return False

    @cached_property
    def equation_map(self) -> dict[str, int]:
        """Map from field name to equation index. Cached on frozen dataclass."""
        return {eq.field_name: i for i, eq in enumerate(self.equations)}

    # ------------------------------------------------------------------ #
    # Perturbative-order filtering (v6 plan, Stage 2)                    #
    # ------------------------------------------------------------------ #
    # These methods partition an EquationSystem by the ``order_in_eps``  #
    # tag on each OperatorTerm, which is emitted by ExportJSON.wl when a #
    # theory has a ``[perturbation]`` section. Used by                   #
    # :class:`tidal.solver.perturbative_driver.PerturbativeSolver` to    #
    # drive Pass 0 (base) and Pass 1+ (correction) solves from a single #
    # derived theory file.                                               #

    def filter_by_order(self, n: int) -> EquationSystem:
        """Return a copy retaining only RHS terms with ``order_in_eps == n``.

        The equations (and their LHS structures) are preserved unchanged;
        only the ``rhs_terms`` tuple is filtered. Equations whose filtered
        RHS becomes empty are kept — callers may need them for state
        layout purposes (e.g., an evolution equation without a source in
        this order still requires its slot in the integrator).

        Use ``filter_by_order(0)`` for the Pass 0 base equations and
        ``filter_by_order(1)`` for the Pass 1 source terms of a linear
        perturbative expansion.
        """
        new_eqs = tuple(
            dataclasses.replace(
                eq,
                rhs_terms=tuple(t for t in eq.rhs_terms if t.order_in_eps == n),
            )
            for eq in self.equations
        )
        # Recompute mass/coupling matrices from the filtered equations
        # (#274). Filtering RHS terms changes the identity-
        # operator structure from which the matrices are derived; the
        # cached matrices on self reflect the full-spec structure and
        # would trigger __post_init__'s inconsistency UserWarning if
        # passed through unchanged.
        mass, coupling, mass_sym, coupling_sym = self._compute_matrices_from_terms(
            new_eqs,
            self.component_names,
            parameters=None,
        )
        return dataclasses.replace(
            self,
            equations=new_eqs,
            mass_matrix=mass,
            coupling_matrix=coupling,
            mass_matrix_symbolic=mass_sym,
            coupling_matrix_symbolic=coupling_sym,
        )

    def max_order(self) -> int:
        """Return the maximum ``order_in_eps`` across all RHS terms.

        Returns 0 for baseline theories (no ``[perturbation]`` section or
        no terms with non-zero order). Use to gate ``--perturbative-order``
        validation and to size the Pass loop in the driver.
        """
        orders = [t.order_in_eps for eq in self.equations for t in eq.rhs_terms]
        return max(orders) if orders else 0

    def has_corrections(self) -> bool:
        """Return True if any RHS term has ``order_in_eps > 0``.

        Cheap check used by the CLI to decide whether
        ``PerturbativeSolver`` is needed or the plain modal path
        suffices.
        """
        return any(t.order_in_eps > 0 for eq in self.equations for t in eq.rhs_terms)

    def canonicalize_kinetic_for_perturbation(
        self,
        small_parameters: Sequence[str],
    ) -> EquationSystem:
        """Translate small-parameter kinetic dependence into order-1 RHS terms.

        Uses the perturbative identity (see #301 Phase 3 / #303)::

            (M₀ + εM₁)⁻¹ (K₀ + εK₁) ≈ M₀⁻¹ K₀ + ε M₀⁻¹ (K₁ − M₁·M₀⁻¹·K₀)

        so the perturbative hierarchy stays clean even when a small parameter
        ε enters the LHS kinetic coefficient. For each equation whose
        ``kinetic_coefficient_symbolic`` mentions a small parameter, splits
        M(ε) into M₀ (parameter-free) and per-parameter corrections c_p·p
        via :func:`split_small_parameter_kinetic`, stores M₀ as the new
        kinetic coefficient, and appends synthesized order-1 RHS terms
        ``-c_p·p·K₀/M₀`` for each base (order-0) RHS term K₀. Pass 0 then
        sees a truly ε=0 baseline; Pass 1's existing Duhamel kernel
        integrates both the original K₁ terms and the synthesized
        corrections identically.

        No-op for equations with:

        * ``kinetic_coefficient_symbolic is None`` (M = 1 implicitly)
        * No small parameter present in the kinetic
        * Empty ``small_parameters`` argument

        Parameters
        ----------
        small_parameters
            The names declared in ``[perturbation].small_parameters`` of
            the TOML.

        Returns
        -------
        A new :class:`EquationSystem` with the canonicalized equations. The
        returned spec is idempotent under this transform; calling it twice
        produces the same result. ``self`` is unchanged.

        Raises
        ------
        tidal.symbolic._kinetic_eval.KineticEvalError
            If any kinetic coefficient has structure outside the
            perturbative contract (bilinear in two small parameters,
            quadratic in one, small-parameter denominator, parenthesized
            sub-sum). See :func:`split_small_parameter_kinetic`.
        """
        from tidal.symbolic._kinetic_eval import (  # noqa: PLC0415
            split_small_parameter_kinetic,
        )

        if not small_parameters:
            return self
        small_params_tuple = tuple(small_parameters)

        new_eqs: list[ComponentEquation] = []
        for eq in self.equations:
            kin_sym = eq.kinetic_coefficient_symbolic
            if kin_sym is None:
                new_eqs.append(eq)
                continue

            m0_expr, c_map = split_small_parameter_kinetic(
                kin_sym,
                small_params_tuple,
            )
            # No small-param dependence in the kinetic → pass through.
            nontrivial = any(c != "0" for c in c_map.values())
            if not nontrivial:
                new_eqs.append(eq)
                continue
            # M₀ collapses to zero (kinetic contains ONLY small parameters,
            # e.g., "b5" with b5 ∈ small_parameters). Such equations are
            # Gap-B demotion candidates: base_spec will turn them into
            # algebraic constraints via lhs_collapses_to_zero. Skip
            # synthesis here — dividing by M₀ = 0 would raise — and pass
            # the original equation through for base_spec to handle.
            from tidal.symbolic._kinetic_eval import (  # noqa: PLC0415
                evaluate_at_zero,
            )

            m0_at_zero = evaluate_at_zero(m0_expr, frozenset())
            if m0_at_zero is not None and m0_at_zero == 0.0:
                new_eqs.append(eq)
                continue

            # Synthesize Pass-1 RHS corrections from M₁·M₀⁻¹·K₀.
            base_rhs = tuple(t for t in eq.rhs_terms if t.order_in_eps == 0)
            synthesized: list[OperatorTerm] = []
            for param_name, c_expr in c_map.items():
                if c_expr == "0":
                    continue
                for base_term in base_rhs:
                    base_sym = base_term.coefficient_symbolic or str(
                        base_term.coefficient,
                    )
                    # coefficient_symbolic = -(c_p) * p * (a_k) / (M₀)
                    #   which is the order-1 source -M₁·M₀⁻¹·K₀ after Wolfram
                    #   sign-preservation conventions. Explicit parentheses
                    #   preserve operator precedence under evaluate_coefficient.
                    synth_sym = (
                        f"-(({c_expr}) * ({param_name}) * ({base_sym})) / ({m0_expr})"
                    )
                    # Numeric coefficient: the value of synth_sym with every
                    # parameter set to 1.0. The solver re-evaluates
                    # coefficient_symbolic at real parameters; the numeric
                    # value is a convenience for unit-test assertions and
                    # validation plots, mirroring the convention in existing
                    # Wolfram-emitted RHS terms.
                    synth_coef = _evaluate_synth_coefficient(
                        synth_sym,
                        base_term,
                    )
                    # GH #380: when split_small_parameter_kinetic translated
                    # `x[]` → `x` in c_expr / m0_expr (for ast.parse), the
                    # resulting synth_sym carries bare coord names. The
                    # downstream OperatorTerm.position_dependent auto-detect
                    # only matches `x[]`-style xCoba calls, so we must set
                    # coordinate_dependent explicitly — otherwise the modal
                    # source-evaluation call will not pass coord_arrays and
                    # eval() raises NameError on the bare coord.
                    detected_coords: list[str] = [
                        coord_name
                        for coord_name in ("t", "x", "y", "z")
                        if re.search(rf"\b{coord_name}\b", synth_sym)
                    ]
                    inherited = tuple(base_term.coordinate_dependent or ())
                    coord_dep = tuple(dict.fromkeys((*inherited, *detected_coords)))
                    synthesized.append(
                        OperatorTerm(
                            coefficient=synth_coef,
                            operator=base_term.operator,
                            field=base_term.field,
                            coefficient_symbolic=synth_sym,
                            coordinate_dependent=coord_dep,
                            order_in_eps=1,
                        ),
                    )

            new_eqs.append(
                dataclasses.replace(
                    eq,
                    kinetic_coefficient_symbolic=m0_expr,
                    rhs_terms=(*eq.rhs_terms, *synthesized),
                ),
            )

        return dataclasses.replace(self, equations=tuple(new_eqs))

    def base_spec(
        self,
        small_parameters: Sequence[str] | None = None,
    ) -> EquationSystem:
        """Return the Pass 0 base spec with LHS demoted where required.

        This is the correct entry point for the iterative perturbative
        driver (v6 plan, Gap B). It extends ``filter_by_order(0)`` with
        a check on each equation's ``kinetic_coefficient_symbolic``:

        * If the kinetic coefficient evaluates to literal zero when
          every small parameter is set to zero, the LHS is promoted by
          the correction. Demote it to an algebraic constraint:
          ``time_derivative_order = 0``, ``kinetic_coefficient_symbolic
          = None``. An ``identity`` self-term is prepended to the RHS
          if one is not already present, so Schur elimination detects
          the field as a proper constraint.
        * Otherwise keep the LHS as-is (the kinetic coefficient
          survives at ε=0 — this is a normal dynamical field or a
          kinetic term that only depends on non-perturbative
          parameters).

        Parameters
        ----------
        small_parameters : sequence of str, optional
            Small-parameter names. When absent, read from
            ``self.metadata.get("perturbation", {}).get(
            "small_parameters", [])``.

        Returns
        -------
        EquationSystem
            A new spec whose equations are the ε=0 base system. Every
            surviving equation has ``time_derivative_order <= 2`` by
            construction; a post-check raises ``ValueError`` if any
            residual has order > 2, indicating the ``[perturbation]``
            config misses some higher-derivative term.

        Raises
        ------
        ValueError
            If any ε=0 base equation still has ``time_derivative_order
            > 2`` after demotion (corresponds to a third- or
            higher-order term that isn't proportional to any declared
            small parameter — spec / config mismatch).
        """
        from tidal.symbolic._kinetic_eval import (  # noqa: PLC0415
            lhs_collapses_to_zero,
        )

        if small_parameters is None:
            pert_meta: dict[str, Any] = self.metadata.get("perturbation", {}) or {}
            small_parameters = list(pert_meta.get("small_parameters", []))

        new_eqs: list[ComponentEquation] = []
        for eq in self.equations:
            # Step 1: filter RHS to order-0 terms (existing logic).
            filtered_rhs = tuple(t for t in eq.rhs_terms if t.order_in_eps == 0)

            # Step 2: check whether the LHS itself is a correction.
            demote = lhs_collapses_to_zero(
                eq.kinetic_coefficient_symbolic,
                small_parameters,
            )

            if demote:
                # Ensure an identity self-term exists on the RHS so
                # downstream Schur elimination recognizes this as a
                # proper algebraic constraint (``1 * field = ...``).
                has_self_identity = any(
                    t.operator == "identity" and t.field == eq.field_name
                    for t in filtered_rhs
                )
                if not has_self_identity:
                    filtered_rhs = (
                        OperatorTerm(
                            coefficient=1.0,
                            operator="identity",
                            field=eq.field_name,
                        ),
                        *filtered_rhs,
                    )
                new_eqs.append(
                    dataclasses.replace(
                        eq,
                        time_derivative_order=0,
                        rhs_terms=filtered_rhs,
                        kinetic_coefficient_symbolic=None,
                    ),
                )
            else:
                new_eqs.append(dataclasses.replace(eq, rhs_terms=filtered_rhs))

        # Post-check: the base system must be 2nd-order at most.
        high_order = [eq.field_name for eq in new_eqs if eq.time_derivative_order > 2]  # noqa: PLR2004
        if high_order:
            msg = (
                f"base_spec: after applying [perturbation]="
                f"{list(small_parameters) if small_parameters else []}, "
                f"fields {high_order} still have time_order > 2 in the "
                f"base theory. The [perturbation] section must name "
                f"every small parameter that appears in a higher-"
                f"derivative kinetic coefficient, so the LHS collapses "
                f"to zero when all small parameters vanish. Add the "
                f"missing parameter(s) to [perturbation].small_parameters "
                f"in your theory.toml and re-derive."
            )
            raise ValueError(msg)

        # Recompute mass/coupling matrices from the demoted equations
        # (#274). The cached matrices on self reflect the
        # full spec's identity-term structure; after demotion the RHS
        # has a new ``1 * identity(self)`` term on each demoted
        # constraint, so mass_matrix[i][i] changes from its full-spec
        # value to reflect the algebraic-constraint form. Without this
        # recompute, ``__post_init__`` would raise UserWarning on every
        # ``base_spec`` call and downstream consumers that read
        # spec.mass_matrix get stale values.
        new_eqs_t = tuple(new_eqs)
        mass, coupling, mass_sym, coupling_sym = self._compute_matrices_from_terms(
            new_eqs_t,
            self.component_names,
            parameters=None,
        )

        return dataclasses.replace(
            self,
            equations=new_eqs_t,
            mass_matrix=mass,
            coupling_matrix=coupling,
            mass_matrix_symbolic=mass_sym,
            coupling_matrix_symbolic=coupling_sym,
        )

    @classmethod
    def from_dict(  # noqa: PLR0914, PLR0912, C901, PLR0915
        cls,
        data: Mapping[str, Any],
        *,
        strict_v6: bool = True,
    ) -> EquationSystem:
        """Create an EquationSystem from a dictionary (parsed JSON).

        Raises
        ------
        ValueError
            If the JSON data is invalid or component references are inconsistent.
        TypeError
            If a metadata parameter has an unsupported type.
        """
        # Extract spacetime info
        spacetime = data["spacetime"]
        dimension = int(spacetime["dimension"])
        spatial_dimension = dimension - 1  # Assuming 1 time dimension

        # Extract fields
        fields_data = data["fields"]
        component_names = tuple(f["name"] for f in fields_data)
        n_components = len(component_names)

        # Validate field name uniqueness
        if len(component_names) != len(set(component_names)):
            seen: set[str] = set()
            duplicates: list[str] = []
            for name in component_names:
                if name in seen:
                    duplicates.append(name)
                seen.add(name)
            msg = f"Duplicate field names: {duplicates}"
            raise ValueError(msg)

        # Build field name -> index lookup
        fields_lookup = {f["name"]: f["index"] for f in fields_data}

        # Extract tensor metadata (optional, from enriched JSON export)
        tensor_metadata: dict[str, dict[str, Any]] = {}
        for f in fields_data:
            if "tensor_head" in f:
                tensor_metadata[f["name"]] = {
                    "tensor_head": f["tensor_head"],
                    "tensor_rank": f.get("tensor_rank", 0),
                    "tensor_indices": f.get("tensor_indices", []),
                }

        # Parse equations
        equations = tuple(
            ComponentEquation.from_dict(eq_data, fields_lookup)
            for eq_data in data["equations"]
        )

        # Extract metadata (needed early for parameter-aware matrix computation)
        metadata = dict(data.get("metadata", {}))

        # Inject tensor metadata into metadata dict for LaTeX rendering
        if tensor_metadata:
            metadata["tensor_metadata"] = tensor_metadata

        # Extract default parameters from metadata (if available) so that
        # numeric matrices reflect actual parameter values, not just ±1.0
        # shape factors from the Wolfram pipeline.
        raw_params = metadata.get("parameters", {})
        default_params: dict[str, float] = {}
        for k, v in raw_params.items():
            if isinstance(v, (int, float)):
                default_params[k] = float(v)
            elif isinstance(v, str):
                try:
                    default_params[k] = float(v)
                except ValueError:
                    msg = f"Parameter '{k}': cannot convert '{v}' to float"
                    raise ValueError(msg) from None
            else:
                msg = f"Parameter '{k}': unsupported type {type(v).__name__}, expected number or numeric string"
                raise TypeError(msg)

        # Auto-compute mass/coupling matrices from identity operator terms.
        # This is the authoritative source — JSON values are ignored in favor
        # of values derived from the equation terms themselves.
        (
            mass_matrix,
            coupling_matrix,
            mass_matrix_symbolic,
            coupling_matrix_symbolic,
        ) = cls._compute_matrices_from_terms(
            equations,
            component_names,
            parameters=default_params or None,
        )

        # Extract coordinate names and metric signature
        coordinates = tuple(str(c) for c in spacetime.get("coordinates", []))
        signature = tuple(int(s) for s in spacetime.get("signature", []))

        # Parse canonical structure (Phase K) — optional for backward compat
        canonical_data = data.get("canonical")
        canonical: CanonicalStructure | None = None
        if canonical_data is not None:
            canonical = CanonicalStructure.from_dict(canonical_data)

        spec = cls(
            n_components=n_components,
            dimension=dimension,
            spatial_dimension=spatial_dimension,
            component_names=component_names,
            equations=equations,
            mass_matrix=mass_matrix,
            coupling_matrix=coupling_matrix,
            metadata=metadata,
            coordinates=coordinates,
            signature=signature,
            mass_matrix_symbolic=mass_matrix_symbolic,
            coupling_matrix_symbolic=coupling_matrix_symbolic,
            canonical=canonical,
        )

        # v6 guard: time_order > 2 equations require a [perturbation]
        # section.  Ostrogradsky reduction was deleted in v6 because it
        # introduces ghost modes; the iterative PerturbativeSolver
        # handles higher-derivative corrections without that cost.  When
        # no [perturbation] is configured, the loader refuses with a
        # migration hint rather than silently producing ghost-ful
        # evolution.
        has_pert = bool(metadata.get("perturbation", {}).get("small_parameters"))
        if (
            any(eq.time_derivative_order > 2 for eq in spec.equations)  # noqa: PLR2004
            and not has_pert
        ):
            fields_with_high_order = [
                eq.field_name
                for eq in spec.equations
                if eq.time_derivative_order > 2  # noqa: PLR2004
            ]
            msg = (
                f"JSON has fields with time_order > 2 "
                f"({fields_with_high_order}) but no [perturbation] section "
                f"in theory.toml. Ostrogradsky reduction was removed in "
                f"v6 because it introduces ghost modes.\n"
                f"Migration: add\n"
                f'    [perturbation]\n    small_parameters = ["<param>"]\n'
                f"to your theory.toml and re-derive. Then `tidal simulate` "
                f"will automatically use `--perturbative-order 1` by "
                f"default.\nSee docs/PERTURBATIVE_REDUCTION_IMPLEMENTATION.md."
            )
            if strict_v6:
                raise ValueError(msg)
            # Inspect / LaTeX-render callers (strict_v6=False) only read the
            # spec and never evolve it, so the Ostrogradsky-ghost concern
            # does not apply. Emit a warning so the condition stays visible
            # in logs.
            warnings.warn(
                "[json_loader] time_order > 2 fields present without "
                "[perturbation]; loaded for read-only use. " + msg.split("\n", 1)[0],
                stacklevel=2,
            )

        return spec


# --- JSON schema validation ---


def _validate_spacetime(spacetime: dict[str, Any]) -> None:
    """Validate spacetime configuration.

    Raises
    ------
    TypeError
        If spacetime is not a dictionary or spacetime.dimension is not an integer.
    ValueError
        If spacetime.dimension is missing.
    """
    if "dimension" not in spacetime:
        msg = "spacetime.dimension is required"
        raise ValueError(msg)
    if not isinstance(spacetime["dimension"], int):
        msg = "spacetime.dimension must be an integer"
        raise TypeError(msg)
    # Validate signature contains only ±1 and matches dimension
    if "signature" in spacetime:
        sig: list[int] = spacetime["signature"]
        if not isinstance(sig, list) or not all(s in {-1, 1} for s in sig):  # type: ignore[reportUnnecessaryIsInstance]
            msg = "spacetime.signature must be a list of +1 or -1 values"
            raise ValueError(msg)
        if len(sig) != spacetime["dimension"]:
            msg = (
                f"spacetime.signature length ({len(sig)}) "
                f"must match dimension ({spacetime['dimension']})"
            )
            raise ValueError(msg)


def _validate_fields(fields: list[Any]) -> None:
    """Validate fields list.

    Raises
    ------
    ValueError
        If the fields list is empty or a field is missing the 'name' key.
    """
    if len(fields) == 0:
        msg = "fields must be non-empty"
        raise ValueError(msg)
    for i, field in enumerate(fields):
        if not isinstance(field, dict) or "name" not in field:
            msg = f"fields[{i}] must be a dict with 'name' key"
            raise ValueError(msg)
    # Validate field indices are unique (when present)
    indices = [f["index"] for f in fields if "index" in f]
    if len(indices) != len(set(indices)):
        msg = f"Field indices must be unique, got: {indices}"
        raise ValueError(msg)


def _validate_equations(equations: list[Any]) -> None:
    """Validate equations list.

    Raises
    ------
    ValueError
        If the equations list is empty or required fields are missing.
    """
    if len(equations) == 0:
        msg = "equations must be a non-empty list"
        raise ValueError(msg)
    for i, eq in enumerate(equations):
        if "field" not in eq:
            msg = f"equations[{i}].field is required"
            raise ValueError(msg)
        if "rhs" not in eq:
            msg = f"equations[{i}].rhs is required"
            raise ValueError(msg)


def validate_json_schema(data: Mapping[str, Any]) -> None:
    """Validate that the JSON data matches the expected schema.

    Raises
    ------
    ValueError
        If required fields are missing or have invalid types.
    """
    required_top_level = ["spacetime", "fields", "equations"]
    for field in required_top_level:
        if field not in data:
            msg = f"Missing required top-level field: {field}"
            raise ValueError(msg)

    _validate_spacetime(data["spacetime"])
    _validate_fields(data["fields"])
    _validate_equations(data["equations"])

    # Cross-validate: equation field references must exist in fields list
    field_names = {f["name"] for f in data["fields"]}
    for i, eq in enumerate(data["equations"]):
        if eq["field"] not in field_names:
            msg = (
                f"equations[{i}].field '{eq['field']}' "
                f"not found in fields list: {sorted(field_names)}"
            )
            raise ValueError(msg)


# --- Public loader ---


def load_equation_system(
    json_path: Path | str,
    *,
    strict_v6: bool = True,
) -> EquationSystem:
    """Load an equation system from a JSON file.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file exported from Mathematica.
    strict_v6 : bool, optional
        When True (default), the v6 ``time_order > 2`` guard raises if the
        JSON has higher-derivative fields without a ``[perturbation]``
        section. Set to False for read-only callers (LaTeX rendering,
        inspection) that do not evolve the system — the guard becomes a
        warning so the spec still loads.

    Returns
    -------
    EquationSystem
        The parsed equation system.

    Raises
    ------
    FileNotFoundError
        If the JSON file does not exist.
    """
    path = Path(json_path)
    if not path.exists():
        msg = f"JSON file not found: {path}"
        raise FileNotFoundError(msg)

    with path.open(encoding="utf-8") as f:
        data = json.load(f)

    validate_json_schema(data)
    return EquationSystem.from_dict(data, strict_v6=strict_v6)


def _evaluate_synth_coefficient(
    synth_sym: str,
    base_term: OperatorTerm,
) -> float:
    """Return the numeric coefficient for a synthesized Phase-3 correction term.

    Evaluates ``synth_sym`` with every free name substituted by 1.0, matching
    the Wolfram-emitter convention that a JSON term's numeric ``coefficient``
    equals ``evaluate(coefficient_symbolic, params={name: 1.0 for ...})``. The
    actual runtime value comes from ``CoefficientEvaluator`` using the real
    parameters; this numeric is a convenience for validation plots and for
    matching the dataclass contract.

    Falls back to ``-base_term.coefficient`` if evaluation raises — a safe
    default (the synthesized term is −(...)·base_coef so ``-base_coef`` is
    the leading-sign approximation when symbolic evaluation hits an
    unsupported node).
    """
    from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

    try:
        # Collect names by parsing; substitute each to 1.0.
        import ast  # noqa: PLC0415

        normalized = synth_sym.replace("^", "**")
        tree = ast.parse(normalized, mode="eval")
        names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
        params = dict.fromkeys(names, 1.0)
        value = evaluate_coefficient(synth_sym, params, ())
        if isinstance(value, float):
            return value
        return float(value)
    except Exception:  # noqa: BLE001
        return -float(base_term.coefficient)


def normalize_kinetic_coefficients(
    spec: EquationSystem,
    params: dict[str, float],
) -> EquationSystem:
    """Apply symbolic kinetic-coefficient normalization to an equation system.

    .. deprecated::
        All time-domain solver backends (cvode/ida/leapfrog/scipy) now consume
        ``kinetic_coefficient_symbolic`` directly via
        :func:`tidal.solver._kinetic.build_inverse_kinetic_diag`, matching the
        modal solver's existing behavior. The canonical spec form carries the
        kinetic coefficient on the LHS and the un-normalized RHS; every
        backend applies M⁻¹ at setup. This function is retained only for
        external callers that pre-normalized specs before the root fix
        landed (see #301, #304). It will be removed in a future release.

    ExportJSON.wl emits equations with a parameter-based kinetic coefficient
    (e.g. ``xi``) **un-divided** when the coefficient cannot be safely divided
    symbolically (to avoid Wolfram producing ``f/xi`` in the JSON which
    becomes a ZeroDivisionError at xi=0).  The kinetic coefficient is stored
    in ``ComponentEquation.kinetic_coefficient_symbolic``.

    This function evaluates those coefficients at ``params`` and:

    * **K ≠ 0** — divides each RHS term by K (numeric ``coefficient / K``;
      symbolic ``coefficient_symbolic`` wrapped as ``(expr) / (kc_sym)`` so
      that subsequent ``CoefficientEvaluator`` calls remain correct).
    * **K = 0** — the kinetic term vanishes; the field becomes a constraint
      (``time_derivative_order=0``, empty ``rhs_terms``).  This is the
      physical xi=0 limit where a formerly dynamical field loses its kinetic
      energy and its EOM degenerates to an algebraic condition.

    Equations without ``kinetic_coefficient_symbolic`` are returned unchanged.
    The function is idempotent: calling it on an already-normalized spec (all
    ``kinetic_coefficient_symbolic`` are None) is a no-op.

    Parameters
    ----------
    spec : EquationSystem
        The equation system as loaded from JSON (pre-normalization).
    params : dict[str, float]
        User-provided parameter values (e.g. ``{"xi": 0.1, "alpha": 0.5}``).

    Returns
    -------
    EquationSystem
        A new ``EquationSystem`` with all kinetic coefficients normalized.
    """
    warnings.warn(
        "normalize_kinetic_coefficients is deprecated. Time-domain solver "
        "backends (cvode/ida/leapfrog/scipy) now consume "
        "kinetic_coefficient_symbolic directly via "
        "tidal.solver._kinetic.build_inverse_kinetic_diag, matching modal. "
        "Drop the pre-normalization call from your code. See #301, #304.",
        DeprecationWarning,
        stacklevel=2,
    )
    new_eqs: list[ComponentEquation] = []
    changed = False

    for eq in spec.equations:
        kc_sym = eq.kinetic_coefficient_symbolic
        if kc_sym is None:
            new_eqs.append(eq)
            continue

        # Evaluate kinetic coefficient at current params
        try:
            kc_value = evaluate_coefficient(kc_sym, params, spec.effective_coordinates)
        except Exception:  # noqa: BLE001
            # If evaluation fails (e.g. symbol not in params), skip normalization
            new_eqs.append(eq)
            continue

        if not isinstance(kc_value, float):
            kc_value = float(kc_value)  # type: ignore[arg-type]

        changed = True

        if kc_value == 0.0:
            # Kinetic coefficient vanishes → field becomes a constraint
            # (time_derivative_order=0, no RHS — algebraic zero constraint)
            new_eq = dataclasses.replace(
                eq,
                time_derivative_order=0,
                rhs_terms=(),
                kinetic_coefficient_symbolic=None,
            )
        else:
            # Divide each RHS term by kc_value
            new_terms: list[OperatorTerm] = []
            for term in eq.rhs_terms:
                new_coeff = term.coefficient / kc_value
                if term.coefficient_symbolic is not None:
                    new_cs = f"({term.coefficient_symbolic}) / ({kc_sym})"
                else:
                    new_cs = term.coefficient_symbolic
                new_terms.append(
                    dataclasses.replace(
                        term,
                        coefficient=new_coeff,
                        coefficient_symbolic=new_cs,
                    ),
                )
            new_eq = dataclasses.replace(
                eq,
                rhs_terms=tuple(new_terms),
                kinetic_coefficient_symbolic=None,
            )

        new_eqs.append(new_eq)

    if not changed:
        return spec

    # Recompute mass/coupling matrices from the normalized equations
    # (#274). Dividing every RHS term by the kinetic
    # coefficient rescales the identity-operator coefficients, so the
    # cached matrices on spec no longer match.
    new_eqs_t = tuple(new_eqs)
    mass, coupling, mass_sym, coupling_sym = spec._compute_matrices_from_terms(  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        new_eqs_t,
        spec.component_names,
        parameters=None,
    )
    return dataclasses.replace(
        spec,
        equations=new_eqs_t,
        mass_matrix=mass,
        coupling_matrix=coupling,
        mass_matrix_symbolic=mass_sym,
        coupling_matrix_symbolic=coupling_sym,
    )
