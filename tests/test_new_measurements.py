"""Tests for the new sweep-oriented measurements.

Covers:
- effective_mass: m²_eff = ω² - k² from dispersion relation
- asymptotic conversion: P_final + directional P_transmitted / P_reflected
- peak_conversion: scalar summary from conversion probability
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tidal.measurement import SimulationData
from tidal.measurement._asymptotic import (
    AsymptoticConversionResult,
    _build_k_grids,
    _directional_split,
    _source_wavevector,
    compute_asymptotic_conversion,
)
from tidal.measurement._effective_mass import (
    EffectiveMassResult,
    compute_effective_mass,
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


def _make_traveling_wave_data(
    k0: float = 3.0,
    n_grid: int = 128,
    n_snapshots: int = 51,
) -> SimulationData:
    """Build two-field data with a right-traveling wave in phi_0.

    phi_0(x, 0) = cos(k0 * x) * gaussian(x)
    chi_0(x, 0) = 0

    Uses the coupled_scalars spec with weak coupling so chi_0 picks
    up a small amount of energy via conversion.
    """
    spec = _build_coupled_scalars_spec()

    domain = 20.0
    dx = domain / n_grid
    x = np.linspace(0, domain, n_grid, endpoint=False)
    t_end = 5.0
    times = np.linspace(0, t_end, n_snapshots)

    m2_phi = float(spec.mass_matrix[0][0])

    # Right-traveling wave packet
    omega = np.sqrt(k0**2 + m2_phi)
    sigma = 2.0
    envelope = np.exp(-((x - domain / 4) ** 2) / (2 * sigma**2))

    phi_fields: list[np.ndarray] = []
    v_phi_fields: list[np.ndarray] = []
    chi_fields: list[np.ndarray] = []
    v_chi_fields: list[np.ndarray] = []

    for t in times:
        phi = envelope * np.cos(k0 * x - omega * t)
        v_phi = envelope * omega * np.sin(k0 * x - omega * t)
        # Small conversion to chi (simplified)
        chi = 0.02 * envelope * np.cos(k0 * x - omega * t + 0.3)
        v_chi = 0.02 * envelope * omega * np.sin(k0 * x - omega * t + 0.3)

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


def _make_standing_wave_data(
    n_grid: int = 64,
    n_snapshots: int = 21,
) -> SimulationData:
    """Build two-field data with a standing wave (no net propagation direction)."""
    spec = _build_coupled_scalars_spec()

    domain = 10.0
    dx = domain / n_grid
    x = np.linspace(0, domain, n_grid, endpoint=False)
    k0 = 2.0 * np.pi / domain * 3  # 3 wavelengths
    m2_phi = float(spec.mass_matrix[0][0])
    omega = np.sqrt(k0**2 + m2_phi)
    times = np.linspace(0, 3.0, n_snapshots)

    phi_fields: list[np.ndarray] = []
    v_phi_fields: list[np.ndarray] = []
    chi_fields: list[np.ndarray] = []
    v_chi_fields: list[np.ndarray] = []

    for t in times:
        # Standing wave = cos(kx)cos(ωt) — no net propagation
        phi = np.cos(k0 * x) * np.cos(omega * t)
        v_phi = -np.cos(k0 * x) * omega * np.sin(omega * t)
        chi = np.zeros(n_grid)
        v_chi = np.zeros(n_grid)

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
# EffectiveMass tests
# ============================================================


class TestEffectiveMass:
    """Tests for compute_effective_mass()."""

    def test_basic_result_type(self) -> None:
        """compute_effective_mass returns EffectiveMassResult."""
        data = _make_traveling_wave_data()
        result = compute_effective_mass(data, ["phi_0"])
        assert isinstance(result, EffectiveMassResult)
        assert result.n_active_modes > 0
        assert result.field_name == "phi_0"

    def test_effective_mass_finite(self) -> None:
        """m²_eff should be a finite number."""
        data = _make_traveling_wave_data()
        result = compute_effective_mass(data, ["phi_0"])
        assert np.isfinite(result.m2_eff)
        assert np.isfinite(result.m2_eff_std)

    def test_arrays_consistent(self) -> None:
        """Per-mode arrays should have consistent shapes."""
        data = _make_traveling_wave_data()
        result = compute_effective_mass(data, ["phi_0"])
        assert len(result.m2_eff_values) == result.n_active_modes
        assert len(result.wavenumbers) == result.n_active_modes
        assert len(result.frequencies) == result.n_active_modes

    def test_group_fields(self) -> None:
        """Effective mass over a field group sums spectral power."""
        data = _make_traveling_wave_data()
        result = compute_effective_mass(data, ["phi_0", "chi_0"])
        assert isinstance(result, EffectiveMassResult)
        assert "phi_0" in result.field_name
        assert "chi_0" in result.field_name


# ============================================================
# AsymptoticConversion tests
# ============================================================


class TestAsymptoticConversion:
    """Tests for compute_asymptotic_conversion()."""

    def test_basic_result(self) -> None:
        """Returns AsymptoticConversionResult with expected fields."""
        data = _make_traveling_wave_data()
        result = compute_asymptotic_conversion(data, "phi_0", "chi_0")
        assert isinstance(result, AsymptoticConversionResult)
        assert result.P_final >= 0
        assert result.P_transmitted >= 0
        assert result.P_reflected >= 0
        assert result.source_field == "phi_0"
        assert result.target_field == "chi_0"

    def test_transmitted_reflected_sum(self) -> None:
        """P_transmitted + P_reflected ≈ P_final."""
        data = _make_traveling_wave_data()
        result = compute_asymptotic_conversion(data, "phi_0", "chi_0")
        assert result.P_transmitted + result.P_reflected == pytest.approx(
            result.P_final,
            rel=1e-2,
        )

    def test_source_wavevector_positive(self) -> None:
        """Source wavevector should have positive k for a right-traveling wave."""
        data = _make_traveling_wave_data(k0=3.0)
        result = compute_asymptotic_conversion(data, "phi_0", "chi_0")
        # 1D: source_wavevector is a 1-tuple
        assert len(result.source_wavevector) == 1
        assert result.source_wavevector[0] > 0  # Right-traveling

    def test_negative_wavevector(self) -> None:
        """Source wavevector should be negative for a left-traveling wave."""
        data = _make_traveling_wave_data(k0=-3.0)
        result = compute_asymptotic_conversion(data, "phi_0", "chi_0")
        assert result.source_wavevector[0] < 0  # Left-traveling

    def test_standing_wave_equal_split(self) -> None:
        """Standing wave has no net direction → equal forward/reflected split."""
        data = _make_standing_wave_data()
        # The chi field has no energy, but test the source wavevector
        k_source = _source_wavevector(data, ["phi_0"])
        # Standing wave: spectral centroid should be ~0 (symmetric positive and
        # negative k contributions cancel)
        assert abs(k_source[0]) < 0.1

    def test_auto_detect_target(self) -> None:
        """target_fields=None auto-detects remaining dynamical fields."""
        data = _make_traveling_wave_data()
        result = compute_asymptotic_conversion(data, "phi_0")
        assert result.target_field == "chi_0"

    def test_overlap_raises(self) -> None:
        """Overlapping source and target raises ValueError."""
        data = _make_traveling_wave_data()
        with pytest.raises(ValueError, match="overlap"):
            compute_asymptotic_conversion(data, "phi_0", "phi_0")

    def test_energies_positive(self) -> None:
        """Source and target energies are non-negative."""
        data = _make_traveling_wave_data()
        result = compute_asymptotic_conversion(data, "phi_0", "chi_0")
        assert result.E_source_initial > 0
        assert result.E_target_final >= 0


class TestBuildKGrids:
    """Tests for _build_k_grids helper."""

    def test_1d(self) -> None:
        """1D k-grid has correct shape."""
        grids = _build_k_grids((32,), (0.1,))
        assert len(grids) == 1
        assert grids[0].shape == (32,)

    def test_2d(self) -> None:
        """2D k-grid has correct shape."""
        grids = _build_k_grids((16, 32), (0.1, 0.2))
        assert len(grids) == 2
        assert grids[0].shape == (16, 32)
        assert grids[1].shape == (16, 32)


class TestSourceWavevector:
    """Tests for _source_wavevector."""

    def test_returns_correct_shape(self) -> None:
        """Source wavevector has ndim components."""
        data = _make_traveling_wave_data()
        k = _source_wavevector(data, ["phi_0"])
        assert k.shape == (1,)  # 1D problem

    def test_plane_wave_recovers_k0(self) -> None:
        """For a clean plane wave, spectral centroid ≈ k0."""
        k0 = 3.0
        data = _make_traveling_wave_data(k0=k0)
        k = _source_wavevector(data, ["phi_0"])
        # The spectral centroid should be close to k0
        # (not exact due to Gaussian envelope and grid discretization)
        assert abs(k[0] - k0) < 1.0


class TestDirectionalSplit:
    """Tests for _directional_split (analytic signal decomposition)."""

    def test_traveling_wave_mostly_forward(self) -> None:
        """Right-traveling wave cos(kx - ωt) → mostly forward energy."""
        data = _make_traveling_wave_data(k0=3.0)
        k_hat = np.array([1.0])  # Forward = positive k
        fwd, ref = _directional_split(data, ["phi_0"], 0, k_hat)
        # Analytic signal correctly assigns most energy to forward direction
        assert fwd > 0.8
        assert ref < 0.2

    def test_fractions_sum_to_one(self) -> None:
        """Forward + reflected should sum to ≈ 1."""
        data = _make_traveling_wave_data()
        k_hat = np.array([1.0])
        fwd, ref = _directional_split(data, ["phi_0"], 0, k_hat)
        assert fwd + ref == pytest.approx(1.0, abs=0.02)

    def test_standing_wave_equal_split(self) -> None:
        """Standing wave cos(kx)cos(ωt) → ~50/50 split."""
        data = _make_standing_wave_data()
        k_hat = np.array([1.0])
        fwd, ref = _directional_split(data, ["phi_0"], 0, k_hat)
        assert fwd == pytest.approx(0.5, abs=0.1)
        assert ref == pytest.approx(0.5, abs=0.1)


# ============================================================
# PeakConversion tests (via measure CLI runner)
# ============================================================


class TestPeakConversion:
    """Tests for _run_peak_conversion in the CLI measure module."""

    def test_returns_scalar_metrics(self) -> None:
        """peak_conversion produces P_max, P_max_time, P_final."""
        from tidal.cli._measure import (
            _run_peak_conversion,  # pyright: ignore[reportPrivateUsage]
        )

        data = _make_traveling_wave_data()
        result = _run_peak_conversion(data, ("phi_0",), ("chi_0",))
        assert "P_max" in result
        assert "P_max_time" in result
        assert "P_final" in result
        assert result["P_max"] >= result["P_final"]  # Peak ≥ final


# ============================================================
# CLI integration tests
# ============================================================


class TestMeasureCLINewTypes:
    """Test that new measurement types are valid --what options."""

    def test_valid_measurement_names(self) -> None:
        """New measurements are in _VALID_MEASUREMENTS."""
        from tidal.cli._measure import _VALID_MEASUREMENTS

        assert "effective_mass" in _VALID_MEASUREMENTS
        assert "asymptotic" in _VALID_MEASUREMENTS
        assert "peak_conversion" in _VALID_MEASUREMENTS

    def test_parse_new_measurements(self) -> None:
        """_parse_measurements accepts the new names."""
        from tidal.cli._measure import (
            _parse_measurements,  # pyright: ignore[reportPrivateUsage]
        )

        result = _parse_measurements("effective_mass,asymptotic,peak_conversion")
        assert result == {"effective_mass", "asymptotic", "peak_conversion"}
