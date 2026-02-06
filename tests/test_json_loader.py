"""Tests for the JSON equation loader."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from torsion_gertsenshtein.symbolic.json_loader import (
    KNOWN_OPERATORS,
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
    validate_json_schema,
)
from torsion_gertsenshtein.symbolic.pde_builder import (
    _OPERATOR_REGISTRY,  # noqa: PLC2701
)

# === Fixtures ===


@pytest.fixture
def em_json_path() -> Path:
    """Path to the EM 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "em_1d.json"


@pytest.fixture
def kg_json_path() -> Path:
    """Path to the Klein-Gordon 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "klein_gordon_1d.json"


@pytest.fixture
def em_json_data() -> dict[str, Any]:
    """EM equation JSON data for testing."""
    return {
        "metadata": {
            "source": "xAct",
            "lagrangian_expr": "-1/4 F[-a,-b] F[a,b]",
            "gauge": "lorenz",
        },
        "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
        "fields": [
            {"name": "A_0", "index": 0, "is_dynamical": True},
            {"name": "A_1", "index": 1, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "A_0",
                "lhs": {"expression": "d2_t(A_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}
                    ],
                },
            },
            {
                "field": "A_1",
                "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "A_1"}
                    ],
                },
            },
        ],
        "coupling": {
            "mass_matrix": [[0.0, 0.0], [0.0, 0.0]],
            "coupling_matrix": [[0.0, 0.0], [0.0, 0.0]],
        },
    }


@pytest.fixture
def kg_json_data() -> dict[str, Any]:
    """Klein-Gordon equation JSON data for testing."""
    return {
        "metadata": {"source": "xAct", "lagrangian_expr": "KG Lagrangian"},
        "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
        "fields": [{"name": "phi", "index": 0, "is_dynamical": True}],
        "equations": [
            {
                "field": "phi",
                "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
                        {"coefficient": -1.0, "operator": "identity", "field": "phi"},
                    ],
                },
            }
        ],
        "coupling": {"mass_matrix": [[1.0]], "coupling_matrix": [[0.0]]},
    }


# === OperatorTerm Tests ===


class TestOperatorTerm:
    """Tests for OperatorTerm dataclass."""

    def test_from_dict(self) -> None:
        """Test creating OperatorTerm from dictionary."""
        data = {"coefficient": 1.5, "operator": "laplacian", "field": "phi"}
        term = OperatorTerm.from_dict(data)

        assert term.coefficient == 1.5  # noqa: PLR2004
        assert term.operator == "laplacian"
        assert term.field == "phi"

    def test_frozen(self) -> None:
        """Test that OperatorTerm is immutable."""
        term = OperatorTerm(coefficient=1.0, operator="identity", field="A_0")
        with pytest.raises(AttributeError):
            term.coefficient = 2.0  # type: ignore[misc]

    def test_numeric_conversion(self) -> None:
        """Test that coefficient is converted to float."""
        data = {"coefficient": 1, "operator": "laplacian", "field": "phi"}
        term = OperatorTerm.from_dict(data)
        assert isinstance(term.coefficient, float)
        assert term.coefficient == 1.0

    def test_from_dict_with_coefficient_symbolic(self) -> None:
        """Test parsing coefficient_symbolic from JSON."""
        data = {
            "coefficient": -1.0,
            "operator": "identity",
            "field": "phi",
            "coefficient_symbolic": "-m2",
        }
        term = OperatorTerm.from_dict(data)

        assert term.coefficient == -1.0
        assert term.operator == "identity"
        assert term.field == "phi"
        assert term.coefficient_symbolic == "-m2"

    def test_from_dict_without_coefficient_symbolic(self) -> None:
        """Test that coefficient_symbolic defaults to None."""
        data = {"coefficient": 1.5, "operator": "laplacian", "field": "phi"}
        term = OperatorTerm.from_dict(data)

        assert term.coefficient_symbolic is None

    def test_coefficient_symbolic_positive(self) -> None:
        """Test positive symbolic coefficient."""
        data = {
            "coefficient": 1.0,
            "operator": "identity",
            "field": "A_0",
            "coefficient_symbolic": "kappa",
        }
        term = OperatorTerm.from_dict(data)

        assert term.coefficient_symbolic == "kappa"


# === ComponentEquation Tests ===


class TestComponentEquation:
    """Tests for ComponentEquation dataclass."""

    def test_from_dict_wave_equation(self) -> None:
        """Test creating ComponentEquation for wave-type equation."""
        data = {
            "field": "A_0",
            "lhs": {"expression": "d2_t(A_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}
                ],
            },
        }
        fields_lookup = {"A_0": 0, "A_1": 1}

        eq = ComponentEquation.from_dict(data, fields_lookup)

        assert eq.field_name == "A_0"
        assert eq.field_index == 0
        time_derivative_order = 2
        assert eq.time_derivative_order == time_derivative_order
        assert len(eq.rhs_terms) == 1
        assert eq.rhs_terms[0].operator == "laplacian"

    def test_from_dict_with_mass_term(self) -> None:
        """Test equation with mass term."""
        data = {
            "field": "phi",
            "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
                    {"coefficient": -1.0, "operator": "identity", "field": "phi"},
                ],
            },
        }
        fields_lookup = {"phi": 0}

        eq = ComponentEquation.from_dict(data, fields_lookup)

        assert eq.field_name == "phi"
        time_derivative_order = 2
        assert eq.time_derivative_order == time_derivative_order
        term_count = 2
        assert len(eq.rhs_terms) == term_count

        # Check terms
        operators = {term.operator for term in eq.rhs_terms}
        assert "laplacian" in operators
        assert "identity" in operators

    def test_unsupported_rhs_type_raises(self) -> None:
        """Test that unsupported RHS type raises ValueError."""
        data: dict[str, Any] = {
            "field": "phi",
            "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
            "rhs": {"type": "nonlinear", "terms": []},
        }
        with pytest.raises(ValueError, match="Unsupported RHS type"):
            ComponentEquation.from_dict(data, {})


# === EquationSystem Tests ===


class TestEquationSystem:
    """Tests for EquationSystem dataclass."""

    def test_from_dict_em(self, em_json_data: dict[str, Any]) -> None:
        """Test creating EquationSystem from EM JSON data."""
        system = EquationSystem.from_dict(em_json_data)

        num_em_components = 2
        assert system.n_components == num_em_components
        dimension = 2
        assert system.dimension == dimension
        assert system.spatial_dimension == 1
        assert system.component_names == ("A_0", "A_1")
        num_equations = 2
        assert len(system.equations) == num_equations
        assert system.metadata["gauge"] == "lorenz"

    def test_from_dict_kg(self, kg_json_data: dict[str, Any]) -> None:
        """Test creating EquationSystem from Klein-Gordon JSON data."""
        system = EquationSystem.from_dict(kg_json_data)

        num_kg_components = 1
        assert system.n_components == num_kg_components
        dimension = 2
        assert system.dimension == dimension
        assert system.spatial_dimension == 1
        assert system.component_names == ("phi",)
        assert len(system.equations) == 1

        # Check mass matrix
        assert system.mass_matrix == ((1.0,),)

    def test_validation_n_components(self) -> None:
        """Test validation of n_components."""
        with pytest.raises(ValueError, match="n_components must be at least 1"):
            EquationSystem(
                n_components=0,
                dimension=2,
                spatial_dimension=1,
                component_names=(),
                equations=(),
                mass_matrix=(),
                coupling_matrix=(),
                metadata={},
            )

    def test_validation_component_names_length(self) -> None:
        """Test validation of component_names length."""
        with pytest.raises(ValueError, match="component_names length"):
            EquationSystem(
                n_components=2,
                dimension=2,
                spatial_dimension=1,
                component_names=("A_0",),  # Wrong length
                equations=(
                    ComponentEquation("A_0", 0, 2, ()),
                    ComponentEquation("A_1", 1, 2, ()),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )


# === Schema Validation Tests ===


class TestValidateJsonSchema:
    """Tests for JSON schema validation."""

    def test_valid_schema(self, em_json_data: dict[str, Any]) -> None:
        """Test that valid data passes validation."""
        validate_json_schema(em_json_data)  # Should not raise

    def test_missing_spacetime(self) -> None:
        """Test that missing spacetime raises ValueError."""
        data: dict[str, Any] = {"fields": [], "equations": []}
        with pytest.raises(
            ValueError, match="Missing required top-level field: spacetime"
        ):
            validate_json_schema(data)

    def test_missing_fields(self) -> None:
        """Test that missing fields raises ValueError."""
        data: dict[str, Any] = {"spacetime": {"dimension": 2}, "equations": []}
        with pytest.raises(
            ValueError, match="Missing required top-level field: fields"
        ):
            validate_json_schema(data)

    def test_missing_equations(self) -> None:
        """Test that missing equations raises ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2},
            "fields": [{"name": "phi"}],
        }
        with pytest.raises(
            ValueError, match="Missing required top-level field: equations"
        ):
            validate_json_schema(data)

    def test_empty_fields(self) -> None:
        """Test that empty fields list raises ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2},
            "fields": [],
            "equations": [],
        }
        with pytest.raises(ValueError, match="fields must be non-empty"):
            validate_json_schema(data)


# === File Loading Tests ===


class TestLoadEquationSystem:
    """Tests for loading equation systems from files."""

    def test_load_em_file(self, em_json_path: Path) -> None:
        """Test loading EM equations from file."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        system = load_equation_system(em_json_path)

        num_em_components = 2
        assert system.n_components == num_em_components
        assert "A_0" in system.component_names
        assert "A_1" in system.component_names

    def test_load_kg_file(self, kg_json_path: Path) -> None:
        """Test loading Klein-Gordon equations from file."""
        if not kg_json_path.exists():
            pytest.skip(f"Test file not found: {kg_json_path}")

        system = load_equation_system(kg_json_path)

        num_kg_components = 1
        assert system.n_components == num_kg_components
        assert system.component_names == ("phi_0",)

    def test_file_not_found(self, tmp_path: Path) -> None:
        """Test that FileNotFoundError is raised for missing file."""
        nonexistent = tmp_path / "nonexistent.json"
        with pytest.raises(FileNotFoundError):
            load_equation_system(nonexistent)

    def test_load_with_string_path(self, em_json_path: Path) -> None:
        """Test loading with string path."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        system = load_equation_system(str(em_json_path))
        num_em_components = 2
        assert system.n_components == num_em_components


# === Phase 4: Field Reference Validation Tests ===


class TestFieldReferenceValidation:
    """Tests for field reference validation in equation terms."""

    def test_valid_regular_field_references(self) -> None:
        """Test that valid regular field references pass validation."""
        system = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("A_0", "A_1"),
            equations=(
                ComponentEquation(
                    "A_0",
                    0,
                    2,
                    (
                        OperatorTerm(1.0, "laplacian", "A_0"),
                        OperatorTerm(0.5, "gradient_x", "A_1"),  # Cross-field ref
                    ),
                ),
                ComponentEquation("A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )
        # Should not raise - valid references
        assert system.n_components == 2  # noqa: PLR2004

    def test_valid_momentum_field_references(self) -> None:
        """Test that valid momentum field references (pi_*) pass validation."""
        system = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("A_0", "A_1"),
            equations=(
                ComponentEquation(
                    "A_0",
                    0,
                    2,
                    (
                        OperatorTerm(1.0, "laplacian", "A_0"),
                        OperatorTerm(0.5, "gradient_x", "pi_1"),  # Momentum reference
                    ),
                ),
                ComponentEquation("A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )
        # Should not raise - valid momentum reference
        assert system.n_components == 2  # noqa: PLR2004

    def test_invalid_regular_field_reference_raises(self) -> None:
        """Test that invalid regular field reference raises ValueError."""
        with pytest.raises(ValueError, match="Unknown field reference 'B_0'"):
            EquationSystem(
                n_components=2,
                dimension=2,
                spatial_dimension=1,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation(
                        "A_0",
                        0,
                        2,
                        (
                            OperatorTerm(1.0, "laplacian", "B_0"),  # Invalid field
                        ),
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )

    def test_momentum_reference_out_of_range_raises(self) -> None:
        """Test that out-of-range momentum reference raises ValueError."""
        with pytest.raises(ValueError, match="out of range"):
            EquationSystem(
                n_components=2,
                dimension=3,
                spatial_dimension=2,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation(
                        "A_0",
                        0,
                        2,
                        (
                            OperatorTerm(1.0, "gradient_x", "pi_5"),  # Out of range
                        ),
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )

    def test_malformed_momentum_reference_raises(self) -> None:
        """Test that malformed momentum reference raises ValueError."""
        # Non-numeric index
        with pytest.raises(ValueError, match="numeric index"):
            EquationSystem(
                n_components=2,
                dimension=3,
                spatial_dimension=2,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation(
                        "A_0",
                        0,
                        2,
                        (
                            OperatorTerm(1.0, "gradient_x", "pi_abc"),  # Non-numeric
                        ),
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )

        # Wrong format (too many underscores)
        with pytest.raises(ValueError, match="pi_N"):
            EquationSystem(
                n_components=2,
                dimension=3,
                spatial_dimension=2,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation(
                        "A_0",
                        0,
                        2,
                        (
                            OperatorTerm(1.0, "gradient_x", "pi_0_extra"),  # Wrong format
                        ),
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )


# === Phase 8A: Operator Registry Sync Test ===


class TestOperatorRegistrySync:
    """Validate KNOWN_OPERATORS and _OPERATOR_REGISTRY stay consistent."""

    def test_all_registry_operators_are_known(self) -> None:
        """Every operator in _OPERATOR_REGISTRY must be in KNOWN_OPERATORS."""
        unknown = set(_OPERATOR_REGISTRY.keys()) - KNOWN_OPERATORS
        assert unknown == set(), (
            f"Operators in _OPERATOR_REGISTRY but not in KNOWN_OPERATORS: {sorted(unknown)}. "
            "Add them to KNOWN_OPERATORS in json_loader.py."
        )

    def test_known_minus_registry_is_only_first_derivative_t(self) -> None:
        """The only KNOWN_OPERATOR not in _OPERATOR_REGISTRY is first_derivative_t."""
        special_cased = KNOWN_OPERATORS - set(_OPERATOR_REGISTRY.keys())
        assert special_cased == {"first_derivative_t"}, (
            f"Expected only {{'first_derivative_t'}} to be special-cased, "
            f"but got: {sorted(special_cased)}"
        )


# === Phase 8D: Python Validation Tests ===


class TestValidationErrors:
    """Tests for fail-fast validation in json_loader."""

    def test_duplicate_field_names_raises(self) -> None:
        """Duplicate field names in JSON should raise ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
            "fields": [
                {"name": "phi", "index": 0},
                {"name": "phi", "index": 1},  # Duplicate
            ],
            "equations": [
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": "phi"}],
                    },
                },
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": "phi"}],
                    },
                },
            ],
        }
        with pytest.raises(ValueError, match="Duplicate field names"):
            EquationSystem.from_dict(data)

    def test_unknown_operator_raises(self) -> None:
        """Unknown operator in RHS term should raise ValueError."""
        data = {"coefficient": 1.0, "operator": "teleportation", "field": "phi"}
        with pytest.raises(ValueError, match=r"Unknown operator.*teleportation"):
            OperatorTerm.from_dict(data)

    def test_missing_required_keys_raises(self) -> None:
        """Missing required keys in RHS term should raise ValueError."""
        # Missing 'operator'
        data: dict[str, Any] = {"coefficient": 1.0, "field": "phi"}
        with pytest.raises(ValueError, match="missing required keys"):
            OperatorTerm.from_dict(data)

        # Missing 'coefficient' and 'field'
        data2: dict[str, Any] = {"operator": "laplacian"}
        with pytest.raises(ValueError, match="missing required keys"):
            OperatorTerm.from_dict(data2)

    def test_mismatched_mass_matrix_rows_raises(self) -> None:
        """Mass matrix with wrong number of rows should raise ValueError."""
        with pytest.raises(ValueError, match="mass_matrix rows"):
            EquationSystem(
                n_components=2,
                dimension=2,
                spatial_dimension=1,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation("A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)),
                    ComponentEquation("A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)),
                ),
                mass_matrix=((0.0, 0.0),),  # 1 row instead of 2
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )

    def test_mismatched_mass_matrix_cols_raises(self) -> None:
        """Mass matrix with wrong number of columns should raise ValueError."""
        with pytest.raises(ValueError, match="mass_matrix row 0 length"):
            EquationSystem(
                n_components=2,
                dimension=2,
                spatial_dimension=1,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation("A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)),
                    ComponentEquation("A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)),
                ),
                mass_matrix=((0.0,), (0.0,)),  # 1 col instead of 2
                coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
                metadata={},
            )

    def test_mismatched_coupling_matrix_raises(self) -> None:
        """Coupling matrix with wrong dimensions should raise ValueError."""
        with pytest.raises(ValueError, match="coupling_matrix rows"):
            EquationSystem(
                n_components=2,
                dimension=2,
                spatial_dimension=1,
                component_names=("A_0", "A_1"),
                equations=(
                    ComponentEquation("A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)),
                    ComponentEquation("A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0), (0.0, 0.0)),  # 3 rows
                metadata={},
            )

    def test_unknown_field_in_equation_raises(self) -> None:
        """Field in equation that doesn't exist in fields_lookup should raise ValueError."""
        data: dict[str, Any] = {
            "field": "nonexistent",
            "lhs": {"expression": "d2_t(nonexistent)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": "phi"}],
            },
        }
        with pytest.raises(ValueError, match="Unknown field 'nonexistent'"):
            ComponentEquation.from_dict(data, {"phi": 0})
