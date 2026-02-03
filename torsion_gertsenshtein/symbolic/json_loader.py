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


def validate_json_schema(data: Mapping[str, Any]) -> None:
    """Validate that the JSON data matches the expected schema.

    Raises
    ------
    ValueError
        If required fields are missing or have invalid types.
    TypeError
        If fields have invalid types (e.g., spacetime.dimension is not an integer).
    """
    required_top_level = ["spacetime", "fields", "equations"]
    for field in required_top_level:
        if field not in data:
            msg = f"Missing required top-level field: {field}"
            raise ValueError(msg)

    # Validate spacetime
    spacetime = data["spacetime"]
    if "dimension" not in spacetime:
        msg = "spacetime.dimension is required"
        raise ValueError(msg)
    if not isinstance(spacetime["dimension"], int):
        msg = "spacetime.dimension must be an integer"
        raise TypeError(msg)

    # Validate fields
    fields = data["fields"]
    if not isinstance(fields, list) or len(fields) == 0:
        msg = "fields must be a non-empty list"
        raise ValueError(msg)
    for i, field in enumerate(fields):
        if "name" not in field:
            msg = f"fields[{i}].name is required"
            raise ValueError(msg)

    # Validate equations
    equations = data["equations"]
    if not isinstance(equations, list) or len(equations) == 0:
        msg = "equations must be a non-empty list"
        raise ValueError(msg)
    for i, eq in enumerate(equations):
        if "field" not in eq:
            msg = f"equations[{i}].field is required"
            raise ValueError(msg)
        if "rhs" not in eq:
            msg = f"equations[{i}].rhs is required"
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
