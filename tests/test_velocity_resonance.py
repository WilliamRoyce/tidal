"""Tests for velocity extraction and resonance analysis measurements.

Covers F3a (group/phase velocity) and F3b (resonance analysis).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tidal.measurement import SimulationData
from tidal.measurement._resonance import (
    ResonanceResult,
    compute_resonance_analysis,
)
from tidal.measurement._velocity import (
    VelocityMismatchResult,
    VelocityResult,
    _compute_group_velocity,
    compute_velocities,
    compute_velocity_mismatch,
)
from tidal.symbolic.json_loader import EquationSystem

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"


# ============================================================
# Helpers
# ============================================================


def _build_coupled_scalars_spec() -> EquationSystem:
    """Build coupled scalar spec using the conftest inline spec (identity coupling)."""
    from tests.conftest import _COUPLED_SCALARS_SPEC

    return EquationSystem.from_dict(_COUPLED_SCALARS_SPEC)  # type: ignore[arg-type]


def _make_two_field_wave_data(
    k0: float = 3.0,
    m2_phi: float = 1.0,
    m2_chi: float = 1.0,
    n_grid: int = 128,
    n_snapshots: int = 51,
    t_end: float = 5.0,
) -> SimulationData:
    """Build two-field data with travelling waves at different masses.

    Both fields have the same wavenumber k0 but different frequencies
    omega_phi = sqrt(k0^2 + m2_phi) and omega_chi = sqrt(k0^2 + m2_chi).
    """
    spec = _build_coupled_scalars_spec()

    domain = 20.0
    dx = domain / n_grid
    x = np.linspace(0, domain, n_grid, endpoint=False)
    times = np.linspace(0, t_end, n_snapshots)

    omega_phi = np.sqrt(k0**2 + m2_phi)
    omega_chi = np.sqrt(k0**2 + m2_chi)

    phi_fields: list[np.ndarray] = []
    v_phi_fields: list[np.ndarray] = []
    chi_fields: list[np.ndarray] = []
    v_chi_fields: list[np.ndarray] = []

    for t in times:
        phi = np.cos(k0 * x - omega_phi * t)
        v_phi = omega_phi * np.sin(k0 * x - omega_phi * t)
        chi = 0.5 * np.cos(k0 * x - omega_chi * t)
        v_chi = 0.5 * omega_chi * np.sin(k0 * x - omega_chi * t)

        phi_fields.append(phi)
        v_phi_fields.append(v_phi)
        chi_fields.append(chi)
        v_chi_fields.append(v_chi)

    params = {
        k: float(v)
        for k, v in spec.metadata.get("parameters", {}).items()
        if isinstance(v, (int, float))
    }

    return SimulationData(
        times=times,
        fields={
            "phi_0": np.stack(phi_fields),
            "chi_0": np.stack(chi_fields),
        },
        velocities={
            "phi_0": np.stack(v_phi_fields),
            "chi_0": np.stack(v_chi_fields),
        },
        grid_spacing=(dx,),
        grid_bounds=((0.0, domain),),
        periodic=(True,),
        spec=spec,
        parameters=params,
    )


# ============================================================
# _compute_group_velocity unit tests
# ============================================================


class TestComputeGroupVelocity:
    """Test the numerical derivative helper."""

    def test_linear_omega_gives_constant_vg(self) -> None:
        """For omega = c * k, group velocity should be c everywhere."""
        k = np.linspace(0.5, 5.0, 20)
        c = 0.8
        omega = c * k
        vg = _compute_group_velocity(k, omega)
        np.testing.assert_allclose(vg, c, atol=0.01)

    def test_quadratic_omega(self) -> None:
        """For omega = k^2, v_g = 2k."""
        k = np.linspace(1.0, 5.0, 50)
        omega = k**2
        vg = _compute_group_velocity(k, omega)
        expected = 2 * k
        # Interior points should be very close; edges less accurate
        np.testing.assert_allclose(vg[2:-2], expected[2:-2], rtol=0.05)

    def test_single_point(self) -> None:
        """Single point returns zero velocity."""
        k = np.array([1.0])
        omega = np.array([2.0])
        vg = _compute_group_velocity(k, omega)
        assert len(vg) == 1
        assert vg[0] == 0.0

    def test_two_points(self) -> None:
        """Two points use forward/backward differences."""
        k = np.array([1.0, 2.0])
        omega = np.array([3.0, 5.0])
        vg = _compute_group_velocity(k, omega)
        assert len(vg) == 2
        np.testing.assert_allclose(vg, [2.0, 2.0])


# ============================================================
# Velocity extraction tests
# ============================================================


class TestVelocities:
    """Test group and phase velocity extraction."""

    def test_basic_extraction(self) -> None:
        """Verify velocity extraction runs without error."""
        data = _make_two_field_wave_data()
        result = compute_velocities(data, "phi_0")
        assert isinstance(result, VelocityResult)
        assert result.n_active_modes > 0
        assert result.group_velocity_mean > 0
        assert result.phase_velocity_mean > 0

    def test_phase_velocity_positive(self) -> None:
        """Phase velocity should be positive for all active modes with k > 0."""
        data = _make_two_field_wave_data()
        result = compute_velocities(data, "phi_0")
        # Phase velocity should be non-negative everywhere
        assert np.all(result.phase_velocity >= 0.0)
        # Mean should be positive (massive field has v_p > 0)
        assert result.phase_velocity_mean > 0.0

    def test_field_name_preserved(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_velocities(data, "phi_0")
        assert result.field_name == "phi_0"

    def test_multi_field_group(self) -> None:
        """Can extract velocity for a field group."""
        data = _make_two_field_wave_data()
        result = compute_velocities(data, ["phi_0", "chi_0"])
        assert isinstance(result, VelocityResult)
        assert "phi_0" in result.field_name
        assert "chi_0" in result.field_name


# ============================================================
# Velocity mismatch tests
# ============================================================


class TestVelocityMismatch:
    """Test velocity mismatch between two field groups."""

    def test_same_mass_small_mismatch(self) -> None:
        """Same mass fields should have small velocity mismatch."""
        data = _make_two_field_wave_data(m2_phi=1.0, m2_chi=1.0)
        result = compute_velocity_mismatch(data, "phi_0", "chi_0")
        assert isinstance(result, VelocityMismatchResult)
        # Same mass -> same dispersion -> near-zero mismatch
        assert result.max_mismatch < 0.5

    def test_different_mass_larger_mismatch(self) -> None:
        """Different mass fields should have measurable velocity mismatch."""
        data = _make_two_field_wave_data(m2_phi=1.0, m2_chi=4.0)
        result = compute_velocity_mismatch(data, "phi_0", "chi_0")
        assert isinstance(result, VelocityMismatchResult)
        assert len(result.shared_wavenumbers) > 0
        # Different mass means different group velocities
        assert result.mean_mismatch >= 0.0

    def test_source_velocity_populated(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_velocity_mismatch(data, "phi_0", "chi_0")
        assert result.source_velocity.n_active_modes > 0
        assert result.target_velocity.n_active_modes > 0


# ============================================================
# Resonance analysis tests
# ============================================================


class TestResonanceAnalysis:
    """Test resonance mode identification."""

    def test_same_mass_many_resonant(self) -> None:
        """Equal masses should give many resonant modes."""
        data = _make_two_field_wave_data(m2_phi=1.0, m2_chi=1.0)
        result = compute_resonance_analysis(data, "phi_0", "chi_0")
        assert isinstance(result, ResonanceResult)
        # With same mass, all modes should be near-resonant
        assert result.n_resonant_modes > 0

    def test_different_mass_fewer_resonant(self) -> None:
        """Very different masses should give fewer resonant modes."""
        data = _make_two_field_wave_data(m2_phi=0.1, m2_chi=10.0)
        result = compute_resonance_analysis(
            data, "phi_0", "chi_0", resonance_threshold=0.05
        )
        # Very different mass -> large frequency mismatch -> fewer resonant
        assert isinstance(result, ResonanceResult)

    def test_peak_conversion_k_exists(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_resonance_analysis(data, "phi_0", "chi_0")
        assert result.peak_conversion_k >= 0.0

    def test_bandwidth_nonnegative(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_resonance_analysis(data, "phi_0", "chi_0")
        assert result.conversion_bandwidth >= 0.0

    def test_field_names_preserved(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_resonance_analysis(data, "phi_0", "chi_0")
        assert result.source_field == "phi_0"
        assert result.target_field == "chi_0"

    def test_mismatch_shape_matches_wavenumbers(self) -> None:
        data = _make_two_field_wave_data()
        result = compute_resonance_analysis(data, "phi_0", "chi_0")
        assert len(result.resonance_mismatch) == len(result.wavenumbers)
        assert len(result.resonant_modes) == len(result.wavenumbers)

    def test_no_shared_modes_raises(self) -> None:
        """If fields have no shared active modes, should raise ValueError."""
        spec = _build_coupled_scalars_spec()
        domain = 20.0
        dx = domain / 64
        x = np.linspace(0, domain, 64, endpoint=False)
        times = np.linspace(0, 2.0, 21)

        # phi_0 has signal, chi_0 is silent (no active modes)
        phi_fields = [np.cos(2 * x - 3 * t) for t in times]
        chi_fields = [np.zeros(64) for _ in times]
        v_phi = [3 * np.sin(2 * x - 3 * t) for t in times]
        v_chi = [np.zeros(64) for _ in times]

        params = {
            k: float(v)
            for k, v in spec.metadata.get("parameters", {}).items()
            if isinstance(v, (int, float))
        }

        data = SimulationData(
            times=times,
            fields={"phi_0": np.stack(phi_fields), "chi_0": np.stack(chi_fields)},
            velocities={"phi_0": np.stack(v_phi), "chi_0": np.stack(v_chi)},
            grid_spacing=(dx,),
            grid_bounds=((0.0, domain),),
            periodic=(True,),
            spec=spec,
            parameters=params,
        )

        with pytest.raises(ValueError, match="No shared active"):
            compute_resonance_analysis(data, "phi_0", "chi_0")
