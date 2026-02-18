"""Tests for the JSON equation loader."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

from tidal.symbolic.json_loader import (
    _STATIC_OPERATORS,
    ComponentEquation,
    EquationSystem,
    LHSStructure,
    OperatorTerm,
    _resolve_symbolic_coeff,
    load_equation_system,
    validate_json_schema,
)
from tidal.symbolic.pde_builder import (
    _OPERATOR_REGISTRY,
)

if TYPE_CHECKING:
    from pde import ScalarField

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

        assert term.coefficient == 1.5
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

    def test_invalid_signature_values(self) -> None:
        """Test that signature with non ±1 values raises ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 2]},
            "fields": [{"name": "phi"}],
            "equations": [
                {"field": "phi", "rhs": {"type": "linear_combination", "terms": []}}
            ],
        }
        with pytest.raises(ValueError, match="signature must be a list of"):
            validate_json_schema(data)

    def test_signature_dimension_mismatch(self) -> None:
        """Test that signature length != dimension raises ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1, 1]},
            "fields": [{"name": "phi"}],
            "equations": [
                {"field": "phi", "rhs": {"type": "linear_combination", "terms": []}}
            ],
        }
        with pytest.raises(ValueError, match=r"signature length.*must match dimension"):
            validate_json_schema(data)

    def test_duplicate_field_indices(self) -> None:
        """Test that duplicate field indices raise ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2},
            "fields": [{"name": "phi", "index": 0}, {"name": "chi", "index": 0}],
            "equations": [
                {"field": "phi", "rhs": {"type": "linear_combination", "terms": []}},
                {"field": "chi", "rhs": {"type": "linear_combination", "terms": []}},
            ],
        }
        with pytest.raises(ValueError, match="Field indices must be unique"):
            validate_json_schema(data)

    def test_equation_references_nonexistent_field(self) -> None:
        """Test that equation referencing unknown field raises ValueError."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2},
            "fields": [{"name": "phi"}],
            "equations": [
                {
                    "field": "nonexistent",
                    "rhs": {"type": "linear_combination", "terms": []},
                }
            ],
        }
        with pytest.raises(ValueError, match="not found in fields list"):
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
                ComponentEquation(
                    "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )
        # Should not raise - valid references
        assert system.n_components == 2

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
                ComponentEquation(
                    "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )
        # Should not raise - valid momentum reference
        assert system.n_components == 2

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
                            OperatorTerm(
                                1.0, "gradient_x", "pi_0_extra"
                            ),  # Wrong format
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
    """Validate _STATIC_OPERATORS and _OPERATOR_REGISTRY stay consistent."""

    def test_all_registry_operators_are_known(self) -> None:
        """Every operator in _OPERATOR_REGISTRY must be in _STATIC_OPERATORS."""
        unknown = set(_OPERATOR_REGISTRY.keys()) - _STATIC_OPERATORS
        assert unknown == set(), (
            f"Operators in _OPERATOR_REGISTRY but not in _STATIC_OPERATORS: {sorted(unknown)}. "
            "Add them to _STATIC_OPERATORS in json_loader.py."
        )

    def test_all_known_operators_in_registry(self) -> None:
        """Every _STATIC_OPERATORS entry has a handler in _OPERATOR_REGISTRY."""
        special_cased = _STATIC_OPERATORS - set(_OPERATOR_REGISTRY.keys())
        assert special_cased == set(), (
            f"_STATIC_OPERATORS not in _OPERATOR_REGISTRY: {sorted(special_cased)}. "
            "Add them to _OPERATOR_REGISTRY (use None handler for special-cased operators)."
        )


class TestRegisterOperator:
    """Tests for the custom operator registration API."""

    def test_register_and_use_custom_operator(self) -> None:
        """Registered operator is accepted by is_known_operator and usable in PDE."""
        from tidal.symbolic.json_loader import (
            _CUSTOM_OPERATORS,
            is_known_operator,
        )
        from tidal.symbolic.pde_builder import (
            _OPERATOR_REGISTRY,
            register_operator,
        )

        name = "_test_custom_op"

        def _handler(field: ScalarField, _bc: object) -> ScalarField:
            return field * 2.0  # type: ignore[return-value]

        try:
            register_operator(name, _handler, min_dim=1)
            assert is_known_operator(name)
            assert name in _OPERATOR_REGISTRY
        finally:
            # Clean up so other tests aren't affected
            _OPERATOR_REGISTRY.pop(name, None)
            _CUSTOM_OPERATORS.discard(name)

    def test_register_shadow_builtin_raises(self) -> None:
        """Registering an operator that shadows a built-in raises ValueError."""
        from tidal.symbolic.pde_builder import (
            register_operator,
        )

        with pytest.raises(ValueError, match="shadows a built-in"):
            register_operator("laplacian", lambda f, _bc: f, min_dim=1)


# === Phase 8D: Python Validation Tests ===


class TestValidationErrors:
    """Tests for fail-fast validation in json_loader."""

    def test_duplicate_field_names_raises(self) -> None:
        """Duplicate field names in JSON should raise ValueError."""
        data: dict[str, Any] = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [
                {"name": "phi", "index": 0},
                {"name": "phi", "index": 1},  # Duplicate
            ],
            "equations": [
                {
                    "field": "phi",
                    "lhs": {
                        "expression": "d2_t(phi)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi",
                            }
                        ],
                    },
                },
                {
                    "field": "phi",
                    "lhs": {
                        "expression": "d2_t(phi)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi",
                            }
                        ],
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
                    ComponentEquation(
                        "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
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
                    ComponentEquation(
                        "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
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
                    ComponentEquation(
                        "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                    ),
                    ComponentEquation(
                        "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                    ),
                ),
                mass_matrix=((0.0, 0.0), (0.0, 0.0)),
                coupling_matrix=((0.0, 0.0), (0.0, 0.0), (0.0, 0.0)),  # 3 rows
                metadata={},
            )

    def test_unknown_field_in_equation_raises(self) -> None:
        """Field in equation that doesn't exist in fields_lookup should raise ValueError."""
        data: dict[str, Any] = {
            "field": "nonexistent",
            "lhs": {
                "expression": "d2_t(nonexistent)",
                "order": {"time": 2, "space": 0},
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian", "field": "phi"}
                ],
            },
        }
        with pytest.raises(ValueError, match="Unknown field 'nonexistent'"):
            ComponentEquation.from_dict(data, {"phi": 0})


class TestCoordinateDependentTerms:
    """Tests for coordinate_dependent field on OperatorTerm."""

    def test_coordinate_dependent_default_empty(self) -> None:
        """OperatorTerm defaults to empty coordinate_dependent tuple."""
        term = OperatorTerm(1.0, "laplacian", "phi_0")
        assert term.coordinate_dependent == ()
        assert not term.position_dependent

    def test_coordinate_dependent_spatial(self) -> None:
        """OperatorTerm with spatial coordinate_dependent is position_dependent."""
        term = OperatorTerm(
            1.0,
            "laplacian_x",
            "phi_0",
            coefficient_symbolic="x()^2/(2*sphR^2)",
            coordinate_dependent=("x",),
        )
        assert term.coordinate_dependent == ("x",)
        assert term.position_dependent

    def test_coordinate_dependent_time_only(self) -> None:
        """OperatorTerm with only time dependence is NOT position_dependent."""
        term = OperatorTerm(
            1.0,
            "identity",
            "phi_0",
            coefficient_symbolic="E^(2*dSH*t())",
            time_dependent=True,
            coordinate_dependent=("t",),
        )
        assert term.coordinate_dependent == ("t",)
        assert not term.position_dependent
        assert term.time_dependent

    def test_coordinate_dependent_mixed(self) -> None:
        """OperatorTerm with both spatial and time is position_dependent."""
        term = OperatorTerm(
            1.0,
            "identity",
            "phi_0",
            coefficient_symbolic="x()*t()",
            time_dependent=True,
            coordinate_dependent=("x", "t"),
        )
        assert term.position_dependent
        assert term.time_dependent

    def test_from_dict_with_coordinate_dependent(self) -> None:
        """OperatorTerm.from_dict parses coordinate_dependent field."""
        data = {
            "coefficient": 1.0,
            "operator": "laplacian_x",
            "field": "phi_0",
            "coefficient_symbolic": "x()^2/(2*sphR^2)",
            "coordinate_dependent": ["x"],
        }
        term = OperatorTerm.from_dict(data)
        assert term.coordinate_dependent == ("x",)
        assert term.position_dependent

    def test_from_dict_without_coordinate_dependent(self) -> None:
        """OperatorTerm.from_dict defaults to empty when field absent."""
        data = {
            "coefficient": 1.0,
            "operator": "laplacian",
            "field": "phi_0",
        }
        term = OperatorTerm.from_dict(data)
        assert term.coordinate_dependent == ()
        assert not term.position_dependent

    def test_from_dict_multi_coordinate(self) -> None:
        """OperatorTerm.from_dict with multiple coordinate dependencies."""
        data = {
            "coefficient": 1.0,
            "operator": "laplacian_x",
            "field": "phi_0",
            "coefficient_symbolic": "(x()^2*y()^2)/(2*sphR^4)",
            "coordinate_dependent": ["x", "y"],
        }
        term = OperatorTerm.from_dict(data)
        assert term.coordinate_dependent == ("x", "y")
        assert term.position_dependent

    def test_position_dependent_non_cartesian(self) -> None:
        """position_dependent works for non-Cartesian coordinate names."""
        term = OperatorTerm(
            1.0,
            "laplacian",
            "phi_0",
            coordinate_dependent=("r", "theta"),
        )
        assert term.position_dependent

    def test_position_dependent_only_t(self) -> None:
        """A term depending only on 't' is NOT position_dependent."""
        term = OperatorTerm(
            1.0,
            "identity",
            "phi_0",
            coordinate_dependent=("t",),
            time_dependent=True,
        )
        assert not term.position_dependent


class TestEquationSystemCoordinates:
    """Tests for EquationSystem coordinate handling."""

    def _make_spec(self, **kwargs: Any) -> EquationSystem:  # noqa: ANN401 - test helper accepts any EquationSystem parameters
        """Create a minimal EquationSystem with customizable fields."""
        defaults: dict[str, Any] = {
            "n_components": 1,
            "dimension": 2,
            "spatial_dimension": 1,
            "component_names": ("phi",),
            "equations": (
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                ),
            ),
            "mass_matrix": ((0.0,),),
            "coupling_matrix": ((0.0,),),
            "metadata": {},
        }
        defaults.update(kwargs)
        return EquationSystem(**defaults)

    def test_effective_coordinates_inferred_1d(self) -> None:
        """1+1D spec without explicit coordinates infers ('t', 'x')."""
        spec = self._make_spec(dimension=2, spatial_dimension=1)
        assert spec.effective_coordinates == ("t", "x")
        assert spec.spatial_coordinates == ("x",)

    def test_effective_coordinates_inferred_2d(self) -> None:
        """2+1D spec without explicit coordinates infers ('t', 'x', 'y')."""
        spec = self._make_spec(dimension=3, spatial_dimension=2)
        assert spec.effective_coordinates == ("t", "x", "y")
        assert spec.spatial_coordinates == ("x", "y")

    def test_effective_coordinates_explicit(self) -> None:
        """Explicit coordinates are used as-is."""
        spec = self._make_spec(
            dimension=3,
            spatial_dimension=2,
            coordinates=("t", "r", "theta"),
        )
        assert spec.effective_coordinates == ("t", "r", "theta")
        assert spec.spatial_coordinates == ("r", "theta")

    def test_coordinates_from_json(self) -> None:
        """from_dict extracts coordinates from JSON spacetime."""
        data = {
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
            },
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            }
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.coordinates == ("t", "x", "y")
        assert spec.effective_coordinates == ("t", "x", "y")

    def test_coordinates_default_empty_from_json(self) -> None:
        """from_dict defaults to empty tuple when coordinates absent."""
        data = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
            },
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            }
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.coordinates == ()
        assert spec.effective_coordinates == ("t", "x")


class TestAutoComputedMatrices:
    """Tests for auto-computed mass_matrix and coupling_matrix from terms."""

    def test_single_field_mass_from_identity(self) -> None:
        """Identity term coefficient populates mass_matrix diagonal."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -3.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((3.0,),)
        assert spec.coupling_matrix == ((0.0,),)

    def test_multi_field_cross_coupling(self) -> None:
        """Cross-field identity terms populate coupling_matrix off-diagonal."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "chi_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -0.5,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "chi_0",
                            },
                            {
                                "coefficient": -4.0,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                            {
                                "coefficient": -0.5,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((1.0, 0.0), (0.0, 4.0))
        assert spec.coupling_matrix == ((0.0, 0.5), (0.5, 0.0))

    def test_massless_field_zero_matrix(self) -> None:
        """No identity terms produces zero matrices."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "A_0", "index": 0}],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "d2_t(A_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "A_0",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((0.0,),)
        assert spec.coupling_matrix == ((0.0,),)

    def test_json_mass_matrix_ignored(self) -> None:
        """Auto-computation overrides whatever mass_matrix JSON provides."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -2.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
            # JSON says 99.0 — auto-computation should override to 2.0
            "coupling": {"mass_matrix": [[99.0]], "coupling_matrix": [[0.0]]},
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((2.0,),)

    def test_momentum_refs_ignored_in_matrix(self) -> None:
        """Identity terms referencing pi_N should not populate matrices."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "A_0", "index": 0}],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "d2_t(A_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "pi_0",
                            },
                            {
                                "coefficient": -5.0,
                                "operator": "identity",
                                "field": "A_0",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((5.0,),)

    def test_symbolic_matrix_from_terms(self) -> None:
        """Symbolic mass matrix is auto-computed from term coefficient_symbolic."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-procaMassSquared",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((1.0,),)
        # Symbolic is auto-computed from term, not from JSON coupling section
        assert spec.mass_matrix_symbolic == (("-procaMassSquared",),)

    def test_symbolic_matrix_without_json_coupling(self) -> None:
        """Symbolic is available even when JSON has no mass_matrix_symbolic."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-(dSm2*E^(2*dSH*t[]))",
                            },
                        ],
                    },
                }
            ],
            # No mass_matrix_symbolic in coupling — auto-compute from terms
            "coupling": {"mass_matrix": [[0.0]], "coupling_matrix": [[0.0]]},
        }
        spec = EquationSystem.from_dict(data)
        # Numeric is unreliable for time-dependent terms, but follows convention
        assert spec.mass_matrix == ((-1.0,),)
        # Symbolic is the authoritative value, preserved from term
        assert spec.mass_matrix_symbolic == (("-(dSm2*E^(2*dSH*t[]))",),)

    def test_no_symbolic_for_constant_coefficients(self) -> None:
        """Constant identity terms produce no symbolic matrix."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -2.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        assert spec.mass_matrix == ((2.0,),)
        # No symbolic since term has no coefficient_symbolic
        assert spec.mass_matrix_symbolic == ()


# === Integration helpers for real JSON files (skip if absent) ===

_DATA_DIR = Path(__file__).resolve().parent.parent / "examples" / "data"


def _load_json(name: str) -> dict[str, Any]:
    """Load a JSON file from examples/data/, skip test if absent."""
    import json

    path = _DATA_DIR / name
    if not path.exists():
        pytest.skip(f"{name} not found (run tidal derive)")
    with path.open() as f:
        return json.load(f)


class TestParameterResolvedMatrices:
    """Verify that numeric matrices reflect actual parameter values."""

    def test_mass_matrix_resolved_with_parameters(self) -> None:
        """Numeric mass_matrix reflects parameter values, not ±1.0 shape factors."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "coordinates": ["t", "x"]},
            "fields": [{"name": "phi", "index": 0}],
            "metadata": {"parameters": {"m2": 5.0}},
            "equations": [
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi",
                                "coefficient_symbolic": "-m2",
                            },
                            {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        # With m2=5.0: -(-m2) = -(-5.0) = 5.0 (not 1.0 shape factor)
        assert spec.mass_matrix == ((5.0,),)
        assert spec.mass_matrix_symbolic == (("-m2",),)

    def test_coupling_matrix_resolved_with_parameters(self) -> None:
        """Numeric coupling_matrix reflects parameter values."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "coordinates": ["t", "x"]},
            "fields": [{"name": "phi", "index": 0}, {"name": "chi", "index": 1}],
            "metadata": {"parameters": {"m2": 1.0, "g": 0.3}},
            "equations": [
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi",
                                "coefficient_symbolic": "-m2",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "chi",
                                "coefficient_symbolic": "g",
                            },
                        ],
                    },
                },
                {
                    "field": "chi",
                    "lhs": {"expression": "d2_t(chi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi",
                                "coefficient_symbolic": "g",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "chi",
                                "coefficient_symbolic": "-m2",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        # mass: -(-m2) = -(-1.0) = 1.0
        assert spec.mass_matrix == ((1.0, 0.0), (0.0, 1.0))
        # coupling[0][1] = -(g) = -0.3 (phi eq, chi field)
        # coupling[1][0] = -(g) = -0.3 (chi eq, phi field)
        assert spec.coupling_matrix == ((0.0, -0.3), (-0.3, 0.0))

    def test_no_parameters_uses_shape_factor(self) -> None:
        """Without parameters, numeric matrix uses raw coefficient (shape factor)."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "coordinates": ["t", "x"]},
            "fields": [{"name": "phi", "index": 0}],
            "metadata": {},  # No parameters
            "equations": [
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi",
                                "coefficient_symbolic": "-m2",
                            },
                        ],
                    },
                }
            ],
        }
        spec = EquationSystem.from_dict(data)
        # No parameters → falls back to coefficient shape factor: -(-1.0) = 1.0
        assert spec.mass_matrix == ((1.0,),)

    def test_coupled_proca_correct_masses(self) -> None:
        """Coupled Proca: A-fields have mA2=1.0, B-fields have mB2=2.0."""
        spec = EquationSystem.from_dict(_load_json("coupled_proca_3d.json"))
        # A components: -(-mA2) with mA2=1.0 → 1.0
        assert spec.mass_matrix[0][0] == pytest.approx(1.0)
        assert spec.mass_matrix[1][1] == pytest.approx(1.0)
        assert spec.mass_matrix[2][2] == pytest.approx(1.0)
        # B components: -(-mB2) with mB2=2.0 → 2.0
        assert spec.mass_matrix[3][3] == pytest.approx(2.0)
        assert spec.mass_matrix[4][4] == pytest.approx(2.0)
        assert spec.mass_matrix[5][5] == pytest.approx(2.0)
        # Coupling: -(gcoup) with gcoup=0.5 → -0.5
        assert spec.coupling_matrix[0][3] == pytest.approx(-0.5)
        assert spec.coupling_matrix[3][0] == pytest.approx(-0.5)


class TestResolveSymbolicCoeff:
    """Direct tests for _resolve_symbolic_coeff edge cases."""

    def test_simple_parameter(self) -> None:
        assert _resolve_symbolic_coeff("m2", {"m2": 5.0}) == 5.0

    def test_negated_parameter(self) -> None:
        assert _resolve_symbolic_coeff("-m2", {"m2": 5.0}) == -5.0

    def test_compound_expression(self) -> None:
        assert _resolve_symbolic_coeff("-2*m2", {"m2": 3.0}) == -6.0

    def test_power_expression(self) -> None:
        # Mathematica ^ → Python **
        assert _resolve_symbolic_coeff("m2^2", {"m2": 3.0}) == 9.0

    def test_unresolvable_returns_none(self) -> None:
        assert _resolve_symbolic_coeff("unknown", {"m2": 1.0}) is None

    def test_coordinate_dependent_returns_none(self) -> None:
        # Mathematica coordinate syntax can't be resolved to float
        assert _resolve_symbolic_coeff("x[]^(-2)", {"m2": 1.0}) is None

    def test_division_by_zero_returns_none(self) -> None:
        assert _resolve_symbolic_coeff("m2/0", {"m2": 1.0}) is None

    def test_inf_result_returns_none(self) -> None:
        # math.isfinite rejects Inf
        assert _resolve_symbolic_coeff("1e309", {}) is None


class TestLHSStructureMissingTimeOrder:
    """Fail-fast: LHS order must include 'time' key."""

    def test_missing_time_raises(self) -> None:
        with pytest.raises(ValueError, match=r"lhs\.order must specify 'time'"):
            LHSStructure.from_dict({"expression": "phi", "order": {"space": 0}})

    def test_empty_order_raises(self) -> None:
        with pytest.raises(ValueError, match=r"lhs\.order must specify 'time'"):
            LHSStructure.from_dict({"expression": "phi", "order": {}})

    def test_valid_order_passes(self) -> None:
        lhs = LHSStructure.from_dict(
            {"expression": "phi", "order": {"time": 2, "space": 0}}
        )
        assert lhs.time_order == 2


class TestMetadataParameterParsing:
    """Fail-fast: metadata parameters must be numeric."""

    def _make_spec_dict(self, params: dict[str, object]) -> dict[str, object]:
        return {
            "metadata": {"parameters": params},
            "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            }
                        ],
                    },
                }
            ],
        }

    def test_string_numeric_param_accepted(self) -> None:
        spec = EquationSystem.from_dict(self._make_spec_dict({"m2": "1.5"}))
        assert spec is not None

    def test_string_non_numeric_param_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot convert"):
            EquationSystem.from_dict(self._make_spec_dict({"m2": "abc"}))
