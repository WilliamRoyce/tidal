"""Load and validate JSON equation specifications from Mathematica/xAct export.

This module provides the data structures and parsing logic for loading
field equations that were derived symbolically from Lagrangians.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping

#: Canonical spatial axis letters, ordered by dimension index.
#: Supports up to 6 spatial dimensions (sufficient for all foreseeable physics).
AXIS_LETTERS: tuple[str, ...] = ("x", "y", "z", "w", "v", "u")

#: Character class matching all known axis letters (for regex construction).
_AXIS_RE_CLASS = "[" + "".join(AXIS_LETTERS) + "]"

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
    }
)

#: Pattern for generic single-axis Nth-order derivatives: derivative_3_x, derivative_5_y, etc.
_GENERIC_SINGLE_AXIS_RE = re.compile(
    r"^derivative_(\d+)_(" + _AXIS_RE_CLASS + r")$"
)

#: Pattern for generic multi-axis derivatives: derivative_2x_1y, derivative_3x_2z, etc.
_GENERIC_MULTI_AXIS_RE = re.compile(
    r"^derivative_(\d+" + _AXIS_RE_CLASS + r"(?:_\d+" + _AXIS_RE_CLASS + r")*)$"
)

# Backward-compatible alias
KNOWN_OPERATORS: frozenset[str] = _STATIC_OPERATORS

#: Mutable set of user-registered custom operators.
#: Populated via ``register_operator()`` in pde_builder.py.
_CUSTOM_OPERATORS: set[str] = set()


def is_known_operator(name: str) -> bool:
    """Check whether an operator name is recognized.

    Accepts static operators (identity, laplacian, gradient_x, ...),
    user-registered custom operators (via ``register_operator``),
    and dynamic patterns for generic Nth-order derivatives
    (derivative_3_x, derivative_5_y, derivative_2x_1y, ...).
    """
    return (
        name in _STATIC_OPERATORS
        or name in _CUSTOM_OPERATORS
        or bool(_GENERIC_SINGLE_AXIS_RE.match(name))
        or bool(_GENERIC_MULTI_AXIS_RE.match(name))
    )


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

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LHSStructure:
        """Create LHSStructure from structured JSON data.

        Expected format: {"expression": "...", "order": {"time": N, "space": 0}}

        Parameters
        ----------
        data : Mapping[str, Any]
            The structured LHS data from JSON.

        Returns
        -------
        LHSStructure
            Parsed LHS structure.
        """
        expression = str(data.get("expression", ""))
        order = data.get("order", {})
        time_order = int(order.get("time", 2))  # Default to 2 for hyperbolic
        space_order = int(order.get("space", 0))
        return cls(
            expression=expression, time_order=time_order, space_order=space_order
        )


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
    """

    coefficient: float
    operator: str
    field: str
    coefficient_symbolic: str | None = None
    time_dependent: bool = False
    coordinate_dependent: tuple[str, ...] = ()

    @property
    def position_dependent(self) -> bool:
        """Whether the coefficient depends on spatial coordinates.

        A coordinate is spatial if it is not the time coordinate ``"t"``.
        This works for any coordinate naming convention (Cartesian, spherical, etc.).
        """
        return bool(set(self.coordinate_dependent) - {"t"})

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
        )


_VALID_BC_TYPES: frozenset[str] = frozenset({"periodic", "dirichlet", "neumann"})


@dataclass(frozen=True)
class BoundaryCondition:
    """Boundary condition for one spatial axis.

    Attributes
    ----------
    type : str
        One of "periodic", "dirichlet", or "neumann".
    value : float | None
        Fixed value for Dirichlet BCs.
    derivative : float | None
        Fixed normal derivative for Neumann BCs.
    """

    type: str
    value: float | None = None
    derivative: float | None = None

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
            msg = f"Unknown BC type: {bc_type!r}. Valid types: {sorted(_VALID_BC_TYPES)}"
            raise ValueError(msg)
        return cls(
            type=bc_type,
            value=data.get("value"),
            derivative=data.get("derivative"),
        )


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
        default_factory=lambda: {}  # noqa: PIE807  # type: dict[str, BoundaryCondition]
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
        default_factory=ConstraintSolverConfig
    )

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
        cls, data: Mapping[str, Any], fields_lookup: dict[str, int]
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
            data.get("constraint_solver")
        )

        return cls(
            field_name=field_name,
            field_index=fields_lookup[field_name],
            time_derivative_order=time_derivative_order,
            rhs_terms=rhs_terms,
            constraint_solver=constraint_solver,
        )


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
    mass_matrix_symbolic: tuple[tuple[str | None, ...], ...] = ()
    coupling_matrix_symbolic: tuple[tuple[str | None, ...], ...] = ()

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
        expected_mass, expected_coupling, _, _ = self._compute_matrices_from_terms(
            self.equations, self.component_names
        )
        if self.mass_matrix != expected_mass or self.coupling_matrix != expected_coupling:
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
        - Regular field names (e.g., "A_0", "A_1", "phi")
        - Momentum field names (e.g., "pi_0", "pi_1") for mixed time-space derivatives

        Raises
        ------
        ValueError
            If a field reference is invalid.
        """
        valid_fields = set(self.component_names)

        for eq in self.equations:
            for term in eq.rhs_terms:
                field_ref = term.field

                # Check for momentum field reference (pi_*)
                if field_ref.startswith("pi_"):
                    parts = field_ref.split("_")
                    if len(parts) != 2:  # noqa: PLR2004
                        msg = (
                            f"Invalid momentum field reference '{field_ref}' "
                            f"in equation for {eq.field_name}. "
                            f"Expected format 'pi_N' where N is a numeric index."
                        )
                        raise ValueError(msg)

                    idx_str = parts[1]
                    if not idx_str.isdigit():
                        msg = (
                            f"Invalid momentum field index in '{field_ref}' "
                            f"(equation for {eq.field_name}). "
                            f"Expected numeric index, got '{idx_str}'."
                        )
                        raise ValueError(msg)

                    idx = int(idx_str)
                    if not (0 <= idx < self.n_components):
                        msg = (
                            f"Momentum field index {idx} out of range in '{field_ref}' "
                            f"(equation for {eq.field_name}). "
                            f"Valid indices: 0 to {self.n_components - 1}."
                        )
                        raise ValueError(msg)
                # Regular field reference
                elif field_ref not in valid_fields:
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
    ) -> tuple[
        tuple[tuple[float, ...], ...],
        tuple[tuple[float, ...], ...],
        tuple[tuple[str | None, ...], ...],
        tuple[tuple[str | None, ...], ...],
    ]:
        """Extract mass and coupling matrices from identity operator terms.

        Scans each equation's RHS terms for ``identity`` operators acting on
        known field names (not momentum references like ``pi_N``).

        Convention: ``matrix[i][j] = -(coefficient)`` where ``coefficient``
        is the numeric coefficient of the ``identity(field_j)`` term in
        equation *i*.  This makes mass² positive for the standard Lagrangian
        sign convention ``∂²_t φ = … - m² φ``.

        When a term has ``coefficient_symbolic``, it is stored as the
        authoritative representation of that matrix entry. The symbolic
        expression is preserved as-is from the term (not evaluated) so that
        it can be resolved at runtime with actual parameter values.

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
                    neg_coeff = -term.coefficient
                    if i == j:
                        mass[i][j] += neg_coeff
                        if term.coefficient_symbolic is not None:
                            mass_sym[i][j] = term.coefficient_symbolic
                            has_symbolic = True
                    else:
                        coupling[i][j] += neg_coeff
                        if term.coefficient_symbolic is not None:
                            coupling_sym[i][j] = term.coefficient_symbolic
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

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EquationSystem:  # noqa: PLR0914
        """Create an EquationSystem from a dictionary (parsed JSON).

        Raises
        ------
        ValueError
            If the JSON data is invalid or component references are inconsistent.
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

        # Parse equations
        equations = tuple(
            ComponentEquation.from_dict(eq_data, fields_lookup)
            for eq_data in data["equations"]
        )

        # Auto-compute mass/coupling matrices from identity operator terms.
        # This is the authoritative source — JSON values are ignored in favour
        # of values derived from the equation terms themselves.
        # Symbolic expressions are extracted from per-term coefficient_symbolic
        # and preserved unevaluated — they are only resolved at runtime when
        # actual parameter values are supplied to the PDE solver.
        (
            mass_matrix,
            coupling_matrix,
            mass_matrix_symbolic,
            coupling_matrix_symbolic,
        ) = cls._compute_matrices_from_terms(equations, component_names)

        # Extract metadata
        metadata = dict(data.get("metadata", {}))

        # Extract coordinate names
        coordinates = tuple(str(c) for c in spacetime.get("coordinates", []))

        return cls(
            n_components=n_components,
            dimension=dimension,
            spatial_dimension=spatial_dimension,
            component_names=component_names,
            equations=equations,
            mass_matrix=mass_matrix,
            coupling_matrix=coupling_matrix,
            metadata=metadata,
            coordinates=coordinates,
            mass_matrix_symbolic=mass_matrix_symbolic,
            coupling_matrix_symbolic=coupling_matrix_symbolic,
        )


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


def load_equation_system(json_path: Path | str) -> EquationSystem:
    """Load an equation system from a JSON file.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file exported from Mathematica.

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
    return EquationSystem.from_dict(data)
