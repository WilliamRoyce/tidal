"""Load and validate JSON equation specifications from Mathematica/xAct export.

This module provides the data structures and parsing logic for loading
field equations that were derived symbolically from Lagrangians.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping


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
    """

    coefficient: float
    operator: str
    field: str

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> OperatorTerm:
        """Create an OperatorTerm from a dictionary."""
        return cls(
            coefficient=float(data["coefficient"]),
            operator=str(data["operator"]),
            field=str(data["field"]),
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
    """

    field_name: str
    field_index: int
    time_derivative_order: int
    rhs_terms: tuple[OperatorTerm, ...]

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any], fields_lookup: dict[str, int]
    ) -> ComponentEquation:
        """Create a ComponentEquation from a dictionary.

        Raises
        ------
        ValueError
            If the RHS type is not "linear_combination".
        """
        field_name = str(data["field"])

        # Parse LHS to determine time derivative order
        lhs = str(data["lhs"])
        if "d2_t" in lhs:
            time_derivative_order = 2
        elif "d_t" in lhs:
            time_derivative_order = 1
        else:
            time_derivative_order = 0

        # Parse RHS terms
        rhs_data = data["rhs"]
        if rhs_data["type"] != "linear_combination":
            msg = f"Unsupported RHS type: {rhs_data['type']}"
            raise ValueError(msg)

        rhs_terms = tuple(
            OperatorTerm.from_dict(term_data) for term_data in rhs_data["terms"]
        )

        return cls(
            field_name=field_name,
            field_index=fields_lookup.get(field_name, 0),
            time_derivative_order=time_derivative_order,
            rhs_terms=rhs_terms,
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
    """

    n_components: int
    dimension: int
    spatial_dimension: int
    component_names: tuple[str, ...]
    equations: tuple[ComponentEquation, ...]
    mass_matrix: tuple[tuple[float, ...], ...]
    coupling_matrix: tuple[tuple[float, ...], ...]
    metadata: dict[str, Any]

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

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EquationSystem:
        """Create an EquationSystem from a dictionary (parsed JSON)."""
        # Extract spacetime info
        spacetime = data["spacetime"]
        dimension = int(spacetime["dimension"])
        spatial_dimension = dimension - 1  # Assuming 1 time dimension

        # Extract fields
        fields_data = data["fields"]
        component_names = tuple(f["name"] for f in fields_data)
        n_components = len(component_names)

        # Build field name -> index lookup
        fields_lookup = {f["name"]: f["index"] for f in fields_data}

        # Parse equations
        equations = tuple(
            ComponentEquation.from_dict(eq_data, fields_lookup)
            for eq_data in data["equations"]
        )

        # Parse coupling matrices
        coupling_data = data.get("coupling", {})
        mass_matrix = tuple(
            tuple(float(x) for x in row)
            for row in coupling_data.get(
                "mass_matrix", [[0.0] * n_components] * n_components
            )
        )
        coupling_matrix = tuple(
            tuple(float(x) for x in row)
            for row in coupling_data.get(
                "coupling_matrix", [[0.0] * n_components] * n_components
            )
        )

        # Extract metadata
        metadata = dict(data.get("metadata", {}))

        return cls(
            n_components=n_components,
            dimension=dimension,
            spatial_dimension=spatial_dimension,
            component_names=component_names,
            equations=equations,
            mass_matrix=mass_matrix,
            coupling_matrix=coupling_matrix,
            metadata=metadata,
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
