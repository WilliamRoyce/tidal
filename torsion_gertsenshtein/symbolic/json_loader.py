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
_GENERIC_SINGLE_AXIS_RE = re.compile(r"^derivative_(\d+)_([xyz])$")

#: Pattern for generic multi-axis derivatives: derivative_2x_1y, derivative_3x_2z, etc.
_GENERIC_MULTI_AXIS_RE = re.compile(r"^derivative_(\d+[xyz](?:_\d+[xyz])*)$")

# Backward-compatible alias
KNOWN_OPERATORS: frozenset[str] = _STATIC_OPERATORS


def is_known_operator(name: str) -> bool:
    """Check whether an operator name is recognized.

    Accepts both static operators (identity, laplacian, gradient_x, ...)
    and dynamic patterns for generic Nth-order derivatives
    (derivative_3_x, derivative_5_y, derivative_2x_1y, ...).
    """
    return (
        name in _STATIC_OPERATORS
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


@dataclass(frozen=True)
class ConstraintSolverConfig:
    """Configuration for elliptic constraint solving.

    When ``enabled`` is True, the constraint equation is solved at each
    timestep using py-pde's Poisson solver rather than remaining frozen.

    Attributes
    ----------
    enabled : bool
        Whether to solve the constraint elliptically. Default False
        preserves existing frozen-constraint behavior.
    method : str
        Solver method. Currently only ``"poisson"`` is supported.
    boundary_conditions : dict[str, BoundaryCondition]
        Per-axis boundary conditions (e.g., ``{"x": ..., "y": ...}``).
    """

    enabled: bool = False
    method: str = "poisson"
    boundary_conditions: dict[str, BoundaryCondition] = dataclass_field(
        default_factory=dict
    )

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
        """
        if data is None:
            return cls()

        enabled = bool(data.get("enabled", False))
        method = str(data.get("method", "poisson"))

        bc_data = data.get("boundary_conditions", {})
        boundary_conditions = {
            axis: BoundaryCondition.from_dict(bc_dict)
            for axis, bc_dict in bc_data.items()
        }

        return cls(
            enabled=enabled,
            method=method,
            boundary_conditions=boundary_conditions,
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
        """Validate constraint_solver is only enabled for time_order=0."""
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

        if len(self.mass_matrix) != self.n_components:
            msg = f"mass_matrix rows {len(self.mass_matrix)} != n_components {self.n_components}"
            raise ValueError(msg)

        for i, row in enumerate(self.mass_matrix):
            if len(row) != self.n_components:
                msg = f"mass_matrix row {i} length {len(row)} != n_components {self.n_components}"
                raise ValueError(msg)

        if len(self.coupling_matrix) != self.n_components:
            msg = f"coupling_matrix rows {len(self.coupling_matrix)} != n_components {self.n_components}"
            raise ValueError(msg)

        for i, row in enumerate(self.coupling_matrix):
            if len(row) != self.n_components:
                msg = f"coupling_matrix row {i} length {len(row)} != n_components {self.n_components}"
                raise ValueError(msg)

        # Validate field references in equation terms
        self._validate_field_references()

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

        # Parse coupling matrices
        # Note: Using list comprehension to avoid Python mutable aliasing bug
        # where [[0.0] * n] * n creates shared row references
        coupling_data = data.get("coupling", {})

        def _default_zero_matrix(n: int) -> list[list[float]]:
            """Create a proper zero matrix without shared row references."""
            return [[0.0 for _ in range(n)] for _ in range(n)]

        mass_matrix = tuple(
            tuple(float(x) for x in row)
            for row in coupling_data.get(
                "mass_matrix", _default_zero_matrix(n_components)
            )
        )
        coupling_matrix = tuple(
            tuple(float(x) for x in row)
            for row in coupling_data.get(
                "coupling_matrix", _default_zero_matrix(n_components)
            )
        )

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
