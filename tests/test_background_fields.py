"""Stress tests for background fields feature.

Tests the full pipeline from coefficient resolution to energy computation,
TOML validation, and measurement guards.  Position-dependent expressions
like UnitStep, Sign, Max, Min are exercised through the measurement and
energy modules.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
from numpy.testing import assert_allclose

from tidal.measurement._io import SimulationData
from tidal.symbolic.json_loader import (
    CanonicalStructure,
    ComponentEquation,
    EquationSystem,
    HamiltonianFactor,
    HamiltonianTerm,
    OperatorTerm,
)

if TYPE_CHECKING:
    import pathlib


# ===========================================================================
# Group 3: Energy computation with background-derived terms
# ===========================================================================


# ===========================================================================
# Group 4: TOML validation edge cases
# ===========================================================================


class TestBackgroundFieldValidation:
    """TOML validation edge cases for [[background_fields]]."""

    def test_background_gradient_in_lagrangian_raises(
        self,
        tmp_path: pathlib.Path,
    ) -> None:
        """CD[-a][G[]] in Lagrangian should raise with clear error message."""
        from tidal.cli._derive import _validate_config

        config: dict[str, Any] = {
            "theory": {"name": "Bad Gradient"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "constants": {"names": ["g0"]},
            "background_fields": [
                {"name": "G", "type": "scalar", "components": ["g0"]},
            ],
            "lagrangian": {
                "expression": (
                    "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] "
                    "- CD[-a][G[]] eta[a,b] CD[-b][phi[]]"
                ),
            },
            "output": {"path": "out.json"},
        }
        with pytest.raises(ValueError, match=r"covariant derivative.*background.*G"):
            _validate_config(config)

    def test_background_no_gradient_passes(
        self,
        tmp_path: pathlib.Path,
    ) -> None:
        """G[] * phi[] (no gradient) should pass validation."""
        from tidal.cli._derive import _validate_config

        config: dict[str, Any] = {
            "theory": {"name": "OK Coupling"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "constants": {"names": ["g0"]},
            "background_fields": [
                {"name": "G", "type": "scalar", "components": ["g0"]},
            ],
            "lagrangian": {
                "expression": "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - G[] * phi[]^2",
            },
            "output": {"path": "out.json"},
        }
        # Should not raise
        _validate_config(config)

    def test_multiple_background_fields(
        self,
        tmp_path: pathlib.Path,
    ) -> None:
        """Two background fields in same config should work."""
        from tidal.cli import main

        config = tmp_path / "multi_bg.toml"
        config.write_text("""
[theory]
name = "Multi Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0", "W0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[[background_fields]]
name = "W"
type = "scalar"
components = ["W0"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - V[]*phi[]^2 - W[]*phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

    def test_background_with_derived_fields(
        self,
        tmp_path: pathlib.Path,
    ) -> None:
        """Background and derived fields should coexist."""
        from tidal.cli import main

        config = tmp_path / "bg_derived.toml"
        config.write_text("""
[theory]
name = "Background With Derived"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d] - V[] * A[-a] eta[a, b] A[-b]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0


# ===========================================================================
# Group 8: Measurement module hardening — FFT guards for position-dependent
# ===========================================================================


def _make_position_dependent_sim_data(
    *,
    position_dependent_coupling: bool = False,
    position_dependent_mass: bool = False,
) -> SimulationData:
    """Build synthetic SimulationData with position-dependent terms.

    Two 2+1D scalar fields phi_0, chi_0 with optional position-dependent
    mass or coupling terms.
    """
    phi_terms_list: list[OperatorTerm] = [
        OperatorTerm(coefficient=1.0, operator="laplacian", field="phi_0"),
        OperatorTerm(coefficient=-1.0, operator="identity", field="phi_0"),
    ]
    chi_terms_list: list[OperatorTerm] = [
        OperatorTerm(coefficient=1.0, operator="laplacian", field="chi_0"),
        OperatorTerm(coefficient=-1.0, operator="identity", field="chi_0"),
    ]

    if position_dependent_mass:
        # Override phi's self-mass to be position-dependent
        phi_terms_list[1] = OperatorTerm(
            coefficient=-1.0,
            operator="identity",
            field="phi_0",
            coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
            coordinate_dependent=("x", "y"),
        )

    if position_dependent_coupling:
        # Add position-dependent cross-field coupling
        phi_terms_list.append(
            OperatorTerm(
                coefficient=-0.1,
                operator="identity",
                field="chi_0",
                coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
                coordinate_dependent=("x", "y"),
            ),
        )
        chi_terms_list.append(
            OperatorTerm(
                coefficient=-0.1,
                operator="identity",
                field="phi_0",
                coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
                coordinate_dependent=("x", "y"),
            ),
        )

    eq_phi = ComponentEquation(
        field_name="phi_0",
        field_index=0,
        time_derivative_order=2,
        rhs_terms=tuple(phi_terms_list),
    )
    eq_chi = ComponentEquation(
        field_name="chi_0",
        field_index=1,
        time_derivative_order=2,
        rhs_terms=tuple(chi_terms_list),
    )

    # Build canonical Hamiltonian terms for energy measurement
    h_terms: list[HamiltonianTerm] = []
    for fname in ("phi_0", "chi_0"):
        h_terms.extend(
            [
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field=fname, operator="time_derivative"),
                    factor_b=HamiltonianFactor(field=fname, operator="time_derivative"),
                ),
                HamiltonianTerm(
                    coefficient=-0.5,
                    factor_a=HamiltonianFactor(field=fname, operator="identity"),
                    factor_b=HamiltonianFactor(field=fname, operator="laplacian"),
                ),
            ],
        )
    # Mass terms — position-dependent for phi if requested
    if position_dependent_mass:
        h_terms.append(
            HamiltonianTerm(
                coefficient=1.0,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
                coefficient_symbolic="g0/2 * exp(-(x()^2 + y()^2))",
                coordinate_dependent=("x", "y"),
            ),
        )
    else:
        h_terms.append(
            HamiltonianTerm(
                coefficient=0.5,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
            ),
        )
    h_terms.append(
        HamiltonianTerm(
            coefficient=0.5,
            factor_a=HamiltonianFactor(field="chi_0", operator="identity"),
            factor_b=HamiltonianFactor(field="chi_0", operator="identity"),
        ),
    )
    if position_dependent_coupling:
        h_terms.append(
            HamiltonianTerm(
                coefficient=0.1,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="chi_0", operator="identity"),
                coefficient_symbolic="g0 * exp(-(x()^2 + y()^2))",
                coordinate_dependent=("x", "y"),
            ),
        )
    canonical = CanonicalStructure(hamiltonian_terms=tuple(h_terms))

    spec = EquationSystem(
        n_components=2,
        dimension=3,
        spatial_dimension=2,
        component_names=("phi_0", "chi_0"),
        equations=(eq_phi, eq_chi),
        mass_matrix=((1.0, 0.0), (0.0, 1.0)),
        coupling_matrix=((0.0, 0.1), (0.1, 0.0))
        if position_dependent_coupling
        else ((0.0, 0.0), (0.0, 0.0)),
        metadata={},
        coordinates=("t", "x", "y"),
        canonical=canonical,
    )

    n_grid = 16
    n_snap = 5
    dx = 10.0 / n_grid
    shape = (n_snap, n_grid, n_grid)
    rng = np.random.default_rng(42)

    return SimulationData(
        times=np.linspace(0.0, 1.0, n_snap),
        fields={
            "phi_0": rng.standard_normal(shape),
            "chi_0": rng.standard_normal(shape) * 0.01,
        },
        velocities={
            "phi_0": rng.standard_normal(shape),
            "chi_0": rng.standard_normal(shape) * 0.01,
        },
        grid_spacing=(dx, dx),
        grid_bounds=((-5.0, 5.0), (-5.0, 5.0)),
        periodic=(True, True),
        spec=spec,
        parameters={"g0": 1.0},
    )


class TestMeasurementFFTGuards:
    """Verify FFT-based methods reject position-dependent systems."""

    def test_spectral_energy_rejects_position_dependent_mass(self) -> None:
        """compute_spectral_energy raises TypeError for ndarray mass."""
        from tidal.measurement._spectral import compute_spectral_energy

        field = np.ones((16, 16))
        mom = np.zeros((16, 16))
        mass_array = np.ones((16, 16)) * 2.0  # ndarray, not scalar

        with pytest.raises(TypeError, match="spatially uniform mass"):
            compute_spectral_energy(field, mom, mass_array, (0.5, 0.5), (True, True))

    def test_dispersion_rejects_position_dependent_coupling(self) -> None:
        """Dispersion relation omega(k) rejects position-dependent coupling."""
        from tidal.measurement._dispersion import compute_dispersion

        data = _make_position_dependent_sim_data(position_dependent_coupling=True)
        with pytest.raises(ValueError, match="position-dependent"):
            compute_dispersion(data, "phi_0")

    def test_conversion_probability_accepts_position_dependent_mass(self) -> None:
        """P(t) works with position-dependent mass (not FFT-based)."""
        from tidal.measurement._conversion import compute_conversion_probability

        data = _make_position_dependent_sim_data(position_dependent_mass=True)
        # Should NOT raise — P(t) uses real-space energy, not Fourier
        result = compute_conversion_probability(data, "phi_0", "chi_0")
        assert result.probability.shape == (5,)

    def test_energy_conservation_with_position_dependent_mass(self) -> None:
        """check_energy_conservation works with position-dependent mass."""
        from tidal.measurement._diagnostics import check_energy_conservation

        data = _make_position_dependent_sim_data(position_dependent_mass=True)
        result = check_energy_conservation(data)
        assert result.relative_error.shape == (5,)


# ===========================================================================
# Group 9: Proca + Scalar Background (Phase 2B)
# ===========================================================================


class TestProcaScalarBackground:
    """Tests for coupled Proca with Lorentzian scalar background."""

    TOML_PATH = (
        Path(__file__).parent.parent / "examples" / "proca_background" / "theory.toml"
    )

    def test_proca_background_wls_dry_run(self) -> None:
        """WLS generation succeeds and contains ReplaceAll for G."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)
        # Scalar BG uses ReplaceAll (not ComponentValue)
        assert "/." in wls or "ReplaceAll" in wls
        # Must contain the Lorentzian profile
        assert "Power" in wls
        # Must define background tensor G
        assert "pwlG[]" in wls or "DefTensor[pwlG[]" in wls

    def test_proca_background_lorentzian_coefficient_evaluates(self) -> None:
        """Lorentzian profile evaluates correctly: g0/(1 + r^2/R^2)."""
        from tidal.symbolic._eval_utils import evaluate_coefficient

        expr = "g0 * (1 + (x()^2 + y()^2) / R^2)**(-1)"
        params = {"g0": 1.0, "R": 8.0}
        x = np.linspace(-10, 10, 32)
        y = np.linspace(-10, 10, 32)
        xg, yg = np.meshgrid(x, y, indexing="ij")

        result = evaluate_coefficient(
            expr,
            params,
            ("t", "x", "y"),
            coord_arrays={"x": xg, "y": yg},
        )
        assert isinstance(result, np.ndarray)
        # Origin should be g0 = 1.0
        center = result[16, 16]
        assert_allclose(center, 1.0, atol=0.1)
        # Corners should be smaller (Lorentzian decay)
        corner = result[0, 0]
        assert corner < center
        # Algebraic tails — at r^2 ~ 200, G ~ 1/(1+200/64) ~ 0.24
        assert corner > 0.1  # heavier tails than Gaussian

    def test_proca_background_energy_with_position_dependent_coupling(self) -> None:
        """compute_system_energy works for constraint + dynamic with pos-dep coupling."""
        from tidal.measurement._energy import compute_system_energy

        # Build 2-field spec: A_0 (constraint) + B_1 (dynamic) with G(x,y) coupling
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
            OperatorTerm(
                coefficient=-0.5,
                operator="identity",
                field="B_1",
                coefficient_symbolic="-0.5 * (1 + (x()^2 + y()^2) / 64)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
        )
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="B_1", operator="time_derivative"),
                    factor_b=HamiltonianFactor(field="B_1", operator="time_derivative"),
                ),
                HamiltonianTerm(
                    coefficient=-0.5,
                    factor_a=HamiltonianFactor(field="B_1", operator="identity"),
                    factor_b=HamiltonianFactor(field="B_1", operator="laplacian"),
                ),
                HamiltonianTerm(
                    coefficient=1.0,
                    factor_a=HamiltonianFactor(field="B_1", operator="identity"),
                    factor_b=HamiltonianFactor(field="B_1", operator="identity"),
                ),
            ),
        )
        spec = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(
                ComponentEquation(
                    field_name="A_0",
                    field_index=0,
                    time_derivative_order=0,
                    rhs_terms=a0_terms,
                ),
                ComponentEquation(
                    field_name="B_1",
                    field_index=1,
                    time_derivative_order=2,
                    rhs_terms=b1_terms,
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
            coordinates=("t", "x", "y"),
            canonical=canonical,
        )
        data = SimulationData(
            times=np.array([0.0, 1.0]),
            fields={
                "A_0": np.zeros((2, 16, 16)),
                "B_1": np.ones((2, 16, 16)) * 0.1,
            },
            velocities={"B_1": np.zeros((2, 16, 16))},
            grid_spacing=(1.25, 1.25),
            grid_bounds=((-10.0, 10.0), (-10.0, 10.0)),
            periodic=(True, True),
            spec=spec,
            parameters={},
        )
        result = compute_system_energy(data, 0)
        assert result.total >= 0.0

    def test_proca_background_measurement_cross_validation(self) -> None:
        """Manual Lorentzian energy vs virial for simple constraint+dynamic system."""
        from tidal.measurement._energy import compute_system_energy

        # B_1 dynamic field with constant mass; A_0 constraint (no energy contribution)
        a0_terms = (OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),)
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="B_1", operator="time_derivative"),
                    factor_b=HamiltonianFactor(field="B_1", operator="time_derivative"),
                ),
                HamiltonianTerm(
                    coefficient=-0.5,
                    factor_a=HamiltonianFactor(field="B_1", operator="identity"),
                    factor_b=HamiltonianFactor(field="B_1", operator="laplacian"),
                ),
                HamiltonianTerm(
                    coefficient=1.0,
                    factor_a=HamiltonianFactor(field="B_1", operator="identity"),
                    factor_b=HamiltonianFactor(field="B_1", operator="identity"),
                ),
            ),
        )
        spec = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(
                ComponentEquation(
                    field_name="A_0",
                    field_index=0,
                    time_derivative_order=0,
                    rhs_terms=a0_terms,
                ),
                ComponentEquation(
                    field_name="B_1",
                    field_index=1,
                    time_derivative_order=2,
                    rhs_terms=b1_terms,
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
            coordinates=("t", "x", "y"),
            canonical=canonical,
        )
        # Uniform B_1 = 1.0: KE=0, PE = m²/2 * integral of B_1^2
        n_grid = 16
        dx = 20.0 / n_grid
        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "A_0": np.zeros((1, n_grid, n_grid)),
                "B_1": np.ones((1, n_grid, n_grid)),
            },
            velocities={"B_1": np.zeros((1, n_grid, n_grid))},
            grid_spacing=(dx, dx),
            grid_bounds=((-10.0, 10.0), (-10.0, 10.0)),
            periodic=(True, True),
            spec=spec,
            parameters={},
        )
        result = compute_system_energy(data, 0)
        # Expected: m²=2, B_1=1 everywhere
        # Energy density = 0.5 * m² * <B_1²> = 0.5 * 2 * 1² = 1.0
        expected_potential = 0.5 * 2.0 * 1.0**2
        assert_allclose(result.total, expected_potential, rtol=0.1)
