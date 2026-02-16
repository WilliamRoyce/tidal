"""Tests for the tidal.measurement module.

Covers energy computation, conversion probability, spectral analysis,
and diagnostics for coupled field systems.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import cast

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from tidal.measurement import (
    ConversionResult,
    DispersionResult,
    EnergyDiagnostics,
    MixingResult,
    MixingSpectrum,
    SimulationData,
    SpectralConversion,
    SystemEnergy,
    check_energy_conservation,
    compute_conversion_probability,
    compute_dispersion,
    compute_energy_timeseries,
    compute_field_energy,
    compute_group_conversion,
    compute_group_spectral_conversion,
    compute_mixing_length,
    compute_mixing_spectrum,
    compute_mode_amplitudes,
    compute_spectral_conversion,
    compute_spectral_energy,
    compute_spectrum,
    compute_system_energy,
)
from tidal.measurement._energy import (
    _apply_spatial_operator,  # pyright: ignore[reportPrivateUsage]
    _compute_constraint_self_energy,  # pyright: ignore[reportPrivateUsage]
    _compute_virial_potential,  # pyright: ignore[reportPrivateUsage]
    _gradient_energy_density,  # pyright: ignore[reportPrivateUsage]
    _is_momentum_field,  # pyright: ignore[reportPrivateUsage]
    _resolve_term_target,  # pyright: ignore[reportPrivateUsage]
    _self_gradient_axes,  # pyright: ignore[reportPrivateUsage]
)
from tidal.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)

# ============================================================
# Helpers
# ============================================================

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"


def _build_coupled_scalars_spec() -> EquationSystem:
    """Load coupled_scalars.json spec."""
    path = DATA_DIR / "coupled_scalars.json"
    if not path.exists():
        pytest.skip("coupled_scalars.json not found")
    return load_equation_system(path)


def _make_sim_data_two_fields(
    n_grid: int = 32,
    n_snapshots: int = 11,
    amplitude: float = 1.0,
) -> SimulationData:
    """Build synthetic SimulationData for two coupled scalar fields.

    Creates a uniform (k=0) mode in phi and zero in chi, then
    computes the exact coupled-oscillator time evolution.
    """
    spec = _build_coupled_scalars_spec()

    dx = 10.0 / n_grid
    grid_spacing = (dx,)
    grid_bounds = ((0.0, 10.0),)
    periodic = (True,)

    # For uniform mode (k=0), the equations reduce to:
    # d²φ/dt² = -m²_φ φ - g χ
    # d²χ/dt² = -m²_χ χ - g φ
    # Use mass matrix from spec to get correct values
    m2_phi = float(spec.mass_matrix[0][0])
    m2_chi = float(spec.mass_matrix[1][1])
    g_val = float(spec.coupling_matrix[0][1])

    # Time range
    t_end = 10.0
    times = np.linspace(0.0, t_end, n_snapshots)

    # Exact solution for uniform mode: eigenfrequency analysis
    m_eff = np.array([[m2_phi, g_val], [g_val, m2_chi]])
    eigenvalues, eigenvectors = np.linalg.eigh(m_eff)
    omega = np.sqrt(np.maximum(eigenvalues, 0.0))

    # IC: phi(0) = amplitude (uniform), chi(0) = 0, pi_phi(0) = 0, pi_chi(0) = 0
    ic = np.array([amplitude, 0.0])
    # Decompose into eigenmodes
    c = eigenvectors.T @ ic  # coefficients in eigenmode basis

    fields_lists: dict[str, list[np.ndarray]] = {"phi_0": [], "chi_0": []}
    momenta_lists: dict[str, list[np.ndarray]] = {"phi_0": [], "chi_0": []}

    for t in times:
        # Mode amplitudes
        mode_vals = c * np.cos(omega * t)
        mode_dots = -c * omega * np.sin(omega * t)

        # Transform back to field basis
        field_vals = eigenvectors @ mode_vals
        mom_vals = eigenvectors @ mode_dots

        # Uniform field across grid
        phi_arr = np.full(n_grid, field_vals[0])
        chi_arr = np.full(n_grid, field_vals[1])
        pi_phi_arr = np.full(n_grid, mom_vals[0])
        pi_chi_arr = np.full(n_grid, mom_vals[1])

        fields_lists["phi_0"].append(phi_arr)
        fields_lists["chi_0"].append(chi_arr)
        momenta_lists["phi_0"].append(pi_phi_arr)
        momenta_lists["chi_0"].append(pi_chi_arr)

    fields_np = {k: np.stack(v) for k, v in fields_lists.items()}
    momenta_np = {k: np.stack(v) for k, v in momenta_lists.items()}

    # Extract metadata parameters so symbolic coefficients resolve correctly
    params = {
        k: float(v)
        for k, v in spec.metadata.get("parameters", {}).items()
        if isinstance(v, (int, float))
    }

    return SimulationData(
        times=times,
        fields=fields_np,
        momenta=momenta_np,
        grid_spacing=grid_spacing,
        grid_bounds=grid_bounds,
        periodic=periodic,
        spec=spec,
        parameters=params,
    )


def _make_single_field_data(
    n_grid: int = 64,
    n_snapshots: int = 5,
) -> SimulationData:
    """Build synthetic data for a single KG field (Gaussian pulse)."""
    path = DATA_DIR / "klein_gordon_1d.json"
    if not path.exists():
        pytest.skip("klein_gordon_1d.json not found")
    spec = load_equation_system(path)

    domain_len = 20.0
    dx = domain_len / n_grid
    x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)
    sigma = 2.0

    times = np.linspace(0.0, 5.0, n_snapshots)

    # Just repeat same Gaussian for simplicity (energy should be constant
    # for the t=0 snapshot test; exact evolution not needed for unit tests)
    phi = np.exp(-((x - domain_len / 2) ** 2) / (2 * sigma**2))
    pi_field = np.zeros(n_grid)

    fields = {"phi_0": np.stack([phi] * n_snapshots)}
    momenta_dict = {"phi_0": np.stack([pi_field] * n_snapshots)}

    return SimulationData(
        times=times,
        fields=fields,
        momenta=momenta_dict,
        grid_spacing=(dx,),
        grid_bounds=((0.0, domain_len),),
        periodic=(True,),
        spec=spec,
        parameters={},
    )


# ============================================================
# SimulationData tests
# ============================================================


class TestSimulationDataFromStorage:
    """Test SimulationData.from_storage."""

    def test_extracts_fields_and_momenta(self) -> None:
        """from_storage correctly maps state_layout to fields/momenta."""
        spec = _build_coupled_scalars_spec()
        grid = CartesianGrid([(0.0, 10.0)], 32, periodic=True)

        # Create synthetic storage with start_writing
        fields_init = [ScalarField.random_uniform(grid) for _ in range(spec.state_size)]
        state_init = FieldCollection(fields_init)
        storage = MemoryStorage()
        storage.start_writing(state_init)
        storage.append(state_init, 0.0)

        for t_val in [1.0, 2.0]:
            fields_t = [
                ScalarField.random_uniform(grid) for _ in range(spec.state_size)
            ]
            state_t = FieldCollection(fields_t)
            storage.append(state_t, t_val)

        data = SimulationData.from_storage(storage, spec, grid)

        assert data.n_snapshots == 3
        assert "phi_0" in data.fields
        assert "chi_0" in data.fields
        assert "phi_0" in data.momenta  # 2nd order → has momentum
        assert "chi_0" in data.momenta
        assert data.fields["phi_0"].shape == (3, 32)
        assert data.momenta["phi_0"].shape == (3, 32)

    def test_empty_storage_raises(self) -> None:
        """from_storage raises ValueError on empty storage."""
        spec = _build_coupled_scalars_spec()
        grid = CartesianGrid([(0, 10)], 32, periodic=True)
        storage = MemoryStorage()

        with pytest.raises(ValueError, match="empty"):
            SimulationData.from_storage(storage, spec, grid)


class TestSimulationDataFromNpz:
    """Test SimulationData.from_npz."""

    def test_round_trip(self, tmp_path: Path) -> None:
        """Save to NPZ and reload, verifying field data matches."""
        data_orig = _make_single_field_data()

        # Save manually
        npz_path = tmp_path / "test.npz"
        save_data: dict[str, np.ndarray] = {"times": data_orig.times}
        for name in data_orig.fields:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"{name}_t{t_idx}"] = data_orig.fields[name][t_idx]
        for name in data_orig.momenta:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"pi_{name}_t{t_idx}"] = data_orig.momenta[name][t_idx]
        np.savez(str(npz_path), **save_data)  # type: ignore[reportArgumentType]

        # Reload
        data_loaded = SimulationData.from_npz(
            npz_path,
            data_orig.spec,
            data_orig.grid_spacing,
            data_orig.grid_bounds,
            data_orig.periodic,
        )

        np.testing.assert_allclose(data_loaded.times, data_orig.times)
        np.testing.assert_allclose(
            data_loaded.fields["phi_0"],
            data_orig.fields["phi_0"],
        )

    def test_missing_file_raises(self) -> None:
        """from_npz raises FileNotFoundError for missing file."""
        spec = _build_coupled_scalars_spec()
        with pytest.raises(FileNotFoundError):
            SimulationData.from_npz(
                "/nonexistent/path.npz",
                spec,
                (0.1,),
                ((0.0, 10.0),),
                (True,),
            )


# ============================================================
# Energy tests
# ============================================================


class TestFieldEnergy:
    """Test compute_field_energy."""

    def test_gaussian_energy_positive(self) -> None:
        """Energy of a Gaussian pulse is positive."""
        data = _make_single_field_data()
        fe = compute_field_energy(
            data.fields["phi_0"][0],
            data.momenta["phi_0"][0],
            mass_squared=1.0,
            grid_spacing=data.grid_spacing,
            periodic=data.periodic,
        )
        assert fe.total > 0
        assert fe.kinetic == 0.0  # zero momentum
        assert fe.gradient > 0  # non-constant field
        assert fe.mass > 0  # nonzero mass and field

    def test_nan_raises(self) -> None:
        """NaN in field data raises ValueError."""
        bad = np.array([1.0, np.nan, 3.0])
        with pytest.raises(ValueError, match="NaN"):
            compute_field_energy(bad, None, 1.0, (0.1,), (True,))

    def test_inf_raises(self) -> None:
        """Inf in field data raises ValueError."""
        bad = np.array([1.0, np.inf, 3.0])
        with pytest.raises(ValueError, match="NaN or Inf"):
            compute_field_energy(bad, None, 1.0, (0.1,), (True,))

    def test_zero_field_zero_energy(self) -> None:
        """All-zero field has zero energy."""
        zero = np.zeros(32)
        fe = compute_field_energy(zero, zero, 1.0, (0.1,), (True,))
        assert fe.total == 0.0


class TestSystemEnergy:
    """Test compute_system_energy."""

    def test_coupled_system(self) -> None:
        """System energy includes interaction term."""
        data = _make_sim_data_two_fields()
        se = compute_system_energy(data, 0)

        assert isinstance(se, SystemEnergy)
        assert "phi_0" in se.per_field
        assert "chi_0" in se.per_field
        assert se.total > 0

    def test_out_of_range_raises(self) -> None:
        """Out-of-range t_idx raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="out of range"):
            compute_system_energy(data, 99)


class TestEnergyTimeseries:
    """Test compute_energy_timeseries."""

    def test_shapes(self) -> None:
        """Timeseries arrays have correct shape."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        times, per_field, _, total = compute_energy_timeseries(data)

        assert len(times) == 11
        assert len(total) == 11
        assert "phi_0" in per_field
        assert len(per_field["phi_0"]) == 11

    def test_coupled_oscillator_energy_conservation(self) -> None:
        """Total energy (field + interaction) is conserved for exact solution."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        _, _, _, total = compute_energy_timeseries(data)

        # Energy should be conserved to floating-point precision
        # (synthetic exact data, no solver error)
        relative_drift = np.abs(total - total[0]) / total[0]
        assert np.max(relative_drift) < 1e-10


# ============================================================
# Conversion probability tests
# ============================================================


class TestConversionProbability:
    """Test compute_conversion_probability."""

    def test_basic_conversion(self) -> None:
        """Conversion probability grows from zero for excited source."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        result = compute_conversion_probability(data, "phi_0", "chi_0")

        assert isinstance(result, ConversionResult)
        assert result.source_field == "phi_0"
        assert result.target_field == "chi_0"
        # P(0) should be ~0 (chi starts at zero)
        assert result.probability[0] < 1e-10
        # P should grow > 0 at some point
        assert np.max(result.probability) > 0.01

    def test_rabi_oscillation(self) -> None:
        """Conversion matches exact analytical P(t) for uniform-mode oscillators.

        For uniform (k=0) initial conditions in a coupled system,
        the exact conversion probability includes both kinetic and mass
        contributions:
            P(t) = [pi_chi(t)^2 + m2_chi * chi(t)^2] / [m2_phi * phi(0)^2]
        """
        data = _make_sim_data_two_fields(n_snapshots=201)

        m2_phi = float(data.spec.mass_matrix[0][0])
        m2_chi = float(data.spec.mass_matrix[1][1])
        g_val = float(data.spec.coupling_matrix[0][1])

        # Compute analytical P(t) at each snapshot time
        m_eff = np.array([[m2_phi, g_val], [g_val, m2_chi]])
        eigenvalues, eigenvectors = np.linalg.eigh(m_eff)
        omega = np.sqrt(np.maximum(eigenvalues, 0.0))

        ic = np.array([1.0, 0.0])  # phi=1, chi=0
        c = eigenvectors.T @ ic

        result = compute_conversion_probability(data, "phi_0", "chi_0")

        # E_phi(0) = 0.5 * L * m2_phi * 1^2 (uniform field, zero momentum)
        # E_chi(t) = 0.5 * L * [pi_chi^2 + m2_chi * chi^2]
        # P(t) = E_chi(t) / E_phi(0) = [pi_chi^2 + m2_chi * chi^2] / m2_phi
        p_expected = np.zeros(len(data.times))
        for i, t in enumerate(data.times):
            chi_t = c[0] * eigenvectors[1, 0] * np.cos(omega[0] * t) + c[
                1
            ] * eigenvectors[1, 1] * np.cos(omega[1] * t)
            pi_chi_t = -c[0] * eigenvectors[1, 0] * omega[0] * np.sin(omega[0] * t) - c[
                1
            ] * eigenvectors[1, 1] * omega[1] * np.sin(omega[1] * t)
            p_expected[i] = (pi_chi_t**2 + m2_chi * chi_t**2) / m2_phi

        # Pointwise comparison — exact data, should match to floating-point precision
        np.testing.assert_allclose(
            result.probability, p_expected, rtol=1e-10, atol=1e-15
        )

        assert result.probability.max() > 0.01  # nontrivial conversion

    def test_same_field_raises(self) -> None:
        """Same source and target raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="different"):
            compute_conversion_probability(data, "phi_0", "phi_0")

    def test_invalid_field_raises(self) -> None:
        """Invalid field name raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="not in spec"):
            compute_conversion_probability(data, "phi_0", "nonexistent")

    def test_zero_source_energy_raises(self) -> None:
        """Zero initial source energy raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5, amplitude=0.0)
        # Both fields zero → zero energy
        with pytest.raises(ValueError, match="zero initial energy"):
            compute_conversion_probability(data, "phi_0", "chi_0")


# ============================================================
# Spectral tests
# ============================================================


class TestSpectrum:
    """Test compute_spectrum."""

    def test_plane_wave_peak(self) -> None:
        """A plane wave produces a spectral peak at the correct k."""
        n = 128
        domain_len = 20.0
        dx = domain_len / n
        x = np.linspace(0, domain_len - dx, n)
        k0 = 2.0 * np.pi * 3 / domain_len  # 3 wavelengths across domain
        field = np.cos(k0 * x)

        snap = compute_spectrum(field, (dx,), (True,))

        # Peak should be near k0
        peak_k = snap.wavenumbers[np.argmax(snap.power_spectrum)]
        assert abs(peak_k - k0) < 0.5  # within half a bin width

    def test_non_periodic_warning(self) -> None:
        """Non-periodic grid emits windowing warning."""
        field = np.random.default_rng(42).standard_normal(32)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compute_spectrum(field, (0.1,), (False,))
            assert any("Hann window" in str(warning.message) for warning in w)


class TestModeAmplitudes:
    """Test compute_mode_amplitudes."""

    def test_shape(self) -> None:
        """Mode amplitudes have correct (n_snapshots, n_modes) shape."""
        data = _make_single_field_data(n_snapshots=5)
        times, wavenumbers, amplitudes = compute_mode_amplitudes(data, "phi_0")

        assert len(times) == 5
        assert len(wavenumbers) > 0
        assert amplitudes.shape[0] == 5
        assert amplitudes.shape[1] == len(wavenumbers)

    def test_invalid_field_raises(self) -> None:
        """Invalid field name raises ValueError."""
        data = _make_single_field_data(n_snapshots=3)
        with pytest.raises(ValueError, match="not in spec"):
            compute_mode_amplitudes(data, "nonexistent")


# ============================================================
# Diagnostics tests
# ============================================================


class TestEnergyConservation:
    """Test check_energy_conservation."""

    def test_conserved_system(self) -> None:
        """Exact coupled oscillator data passes conservation check."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        diag = check_energy_conservation(data, threshold=1e-6)

        assert isinstance(diag, EnergyDiagnostics)
        assert diag.is_conserved is True
        assert diag.max_relative_error < 1e-6

    def test_bad_threshold_raises(self) -> None:
        """Non-positive threshold raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="positive"):
            check_energy_conservation(data, threshold=0.0)


# ============================================================
# New tests from critical review
# ============================================================


class TestNonPeriodicGradient:
    """Test gradient energy with non-periodic boundary conditions."""

    def test_linear_field_gradient(self) -> None:
        """Linear field f(x) = x has constant gradient 1.

        Gradient energy = 0.5 * integral(1^2) * dV = 0.5 * L.
        """
        n = 64
        domain_len = 10.0
        dx = domain_len / n
        x = np.linspace(dx / 2, domain_len - dx / 2, n)
        field = x.copy()

        grad_sq = _gradient_energy_density(field, (dx,), (False,))

        # np.gradient of a linear field with uniform spacing gives 1.0 everywhere
        # (except possibly boundary cells, but close enough)
        np.testing.assert_allclose(grad_sq, 1.0, atol=1e-10)


class TestSpectrum2D:
    """Test spectral analysis in 2D."""

    def test_2d_plane_wave(self) -> None:
        """2D plane wave cos(kx*x + ky*y) peaks at |k| = sqrt(kx^2 + ky^2)."""
        nx, ny = 64, 64
        lx, ly = 10.0, 10.0
        dx, dy = lx / nx, ly / ny
        x = np.linspace(0, lx - dx, nx)
        y = np.linspace(0, ly - dy, ny)
        xx, yy = np.meshgrid(x, y, indexing="ij")

        # 2 wavelengths in x, 3 in y
        kx0 = 2.0 * np.pi * 2 / lx
        ky0 = 2.0 * np.pi * 3 / ly
        k_expected = np.sqrt(kx0**2 + ky0**2)

        field = np.cos(kx0 * xx + ky0 * yy)
        snap = compute_spectrum(field, (dx, dy), (True, True))

        peak_k = snap.wavenumbers[np.argmax(snap.power_spectrum)]
        assert abs(peak_k - k_expected) < 1.0  # within 1 bin width


class TestDecoupledFields:
    """Test that zero coupling gives zero conversion."""

    def test_zero_coupling_no_conversion(self) -> None:
        """With g=0 (decoupled), P(t) = 0 for all t."""
        spec = _build_coupled_scalars_spec()

        n_grid = 32
        n_snapshots = 11
        dx = 10.0 / n_grid
        times = np.linspace(0.0, 10.0, n_snapshots)

        # For g=0, phi just oscillates at its own frequency, chi stays at zero
        m2_phi = float(spec.mass_matrix[0][0])
        omega_phi = np.sqrt(max(m2_phi, 0.0))

        fields_np: dict[str, np.ndarray] = {}
        momenta_np: dict[str, np.ndarray] = {}

        phi_list: list[np.ndarray] = []
        chi_list: list[np.ndarray] = []
        pi_phi_list: list[np.ndarray] = []
        pi_chi_list: list[np.ndarray] = []

        for t in times:
            phi_list.append(np.full(n_grid, np.cos(omega_phi * t)))
            chi_list.append(np.zeros(n_grid))
            pi_phi_list.append(np.full(n_grid, -omega_phi * np.sin(omega_phi * t)))
            pi_chi_list.append(np.zeros(n_grid))

        fields_np["phi_0"] = np.stack(phi_list)
        fields_np["chi_0"] = np.stack(chi_list)
        momenta_np["phi_0"] = np.stack(pi_phi_list)
        momenta_np["chi_0"] = np.stack(pi_chi_list)

        data = SimulationData(
            times=times,
            fields=fields_np,
            momenta=momenta_np,
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        result = compute_conversion_probability(data, "phi_0", "chi_0")
        # chi energy is zero at all times → P(t) = 0
        np.testing.assert_allclose(result.probability, 0.0, atol=1e-12)


class TestNearZeroEnergyThreshold:
    """Test _ENERGY_FLOOR threshold for near-zero source energy."""

    def test_near_zero_source_raises(self) -> None:
        """Near-zero source energy (below _ENERGY_FLOOR) raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5, amplitude=1e-15)
        with pytest.raises(ValueError, match="zero initial energy"):
            compute_conversion_probability(data, "phi_0", "chi_0")

    def test_above_floor_succeeds(self) -> None:
        """Source energy above _ENERGY_FLOOR succeeds."""
        # amplitude=1.0 gives substantial energy — should not raise
        data = _make_sim_data_two_fields(n_snapshots=5, amplitude=1.0)
        result = compute_conversion_probability(data, "phi_0", "chi_0")
        assert result.probability is not None


class TestPositionDependentValidation:
    """Test that position-dependent mass/coupling is rejected."""

    def _make_spec_with_position_dependent_mass(self) -> EquationSystem:
        """Build a synthetic spec with position-dependent identity term."""
        base_spec = _build_coupled_scalars_spec()

        # Replace first equation's identity term with position-dependent one
        new_terms: list[OperatorTerm] = []
        for term in base_spec.equations[0].rhs_terms:
            if term.operator == "identity" and term.field == "phi_0":
                new_terms.append(
                    OperatorTerm(
                        coefficient=term.coefficient,
                        operator="identity",
                        field="phi_0",
                        coefficient_symbolic=term.coefficient_symbolic,
                        coordinate_dependent=("x",),
                    )
                )
            else:
                new_terms.append(term)

        new_eq = ComponentEquation(
            field_name=base_spec.equations[0].field_name,
            field_index=base_spec.equations[0].field_index,
            time_derivative_order=base_spec.equations[0].time_derivative_order,
            rhs_terms=tuple(new_terms),
        )
        new_equations = (new_eq, base_spec.equations[1])

        return EquationSystem(
            n_components=base_spec.n_components,
            dimension=base_spec.dimension,
            spatial_dimension=base_spec.spatial_dimension,
            equations=new_equations,
            component_names=base_spec.component_names,
            mass_matrix=base_spec.mass_matrix,
            coupling_matrix=base_spec.coupling_matrix,
            metadata=base_spec.metadata,
            coordinates=base_spec.coordinates,
            mass_matrix_symbolic=base_spec.mass_matrix_symbolic,
            coupling_matrix_symbolic=base_spec.coupling_matrix_symbolic,
        )

    def test_position_dependent_mass_raises(self) -> None:
        """Position-dependent mass term raises ValueError in system energy."""
        spec = self._make_spec_with_position_dependent_mass()
        data = SimulationData(
            times=np.array([0.0, 1.0]),
            fields={
                "phi_0": np.ones((2, 32)),
                "chi_0": np.zeros((2, 32)),
            },
            momenta={
                "phi_0": np.zeros((2, 32)),
                "chi_0": np.zeros((2, 32)),
            },
            grid_spacing=(0.3125,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        with pytest.raises(ValueError, match="Position-dependent"):
            compute_system_energy(data, 0)


class TestFromNpzAuto:
    """Test SimulationData.from_npz_auto."""

    def test_round_trip_with_parameters(self, tmp_path: Path) -> None:
        """Save with grid metadata + parameters, reload via from_npz_auto."""
        data_orig = _make_single_field_data()
        params = {"m2": 1.0, "g": 0.5}

        # Save with full metadata (mimicking _save_npz)
        npz_path = tmp_path / "test_auto.npz"
        save_data: dict[str, np.ndarray] = {"times": data_orig.times}
        for name in data_orig.fields:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"{name}_t{t_idx}"] = data_orig.fields[name][t_idx]
        for name in data_orig.momenta:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"pi_{name}_t{t_idx}"] = data_orig.momenta[name][t_idx]

        # Grid metadata
        save_data["grid_spacing"] = np.array(data_orig.grid_spacing)
        save_data["grid_bounds"] = np.array(data_orig.grid_bounds)
        save_data["grid_periodic"] = np.array(data_orig.periodic)

        # Parameters
        save_data["_param_names"] = np.array(list(params.keys()))
        save_data["_param_values"] = np.array(list(params.values()))

        np.savez(str(npz_path), **save_data)  # type: ignore[reportArgumentType]

        # Reload via from_npz_auto
        data_loaded = SimulationData.from_npz_auto(npz_path, data_orig.spec)

        np.testing.assert_allclose(data_loaded.times, data_orig.times)
        np.testing.assert_allclose(
            data_loaded.fields["phi_0"],
            data_orig.fields["phi_0"],
        )
        assert data_loaded.parameters == pytest.approx(params)
        assert data_loaded.grid_spacing == pytest.approx(data_orig.grid_spacing)
        assert data_loaded.periodic == data_orig.periodic

    def test_missing_metadata_raises(self, tmp_path: Path) -> None:
        """from_npz_auto raises ValueError when grid metadata is missing."""
        spec = _build_coupled_scalars_spec()
        npz_path = tmp_path / "no_meta.npz"
        np.savez(str(npz_path), times=np.array([0.0, 1.0]))

        with pytest.raises(ValueError, match=r"missing.*grid_spacing"):
            SimulationData.from_npz_auto(npz_path, spec)


# ============================================================
# Group conversion tests
# ============================================================


class TestGroupConversion:
    """Test compute_group_conversion for multi-field groups."""

    def test_single_source_single_target_matches_pairwise(self) -> None:
        """Single-field groups degenerate to pairwise conversion."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        pairwise = compute_conversion_probability(data, "phi_0", "chi_0")
        group = compute_group_conversion(data, "phi_0", "chi_0")
        np.testing.assert_allclose(group.probability, pairwise.probability)
        assert group.source_field == "phi_0"
        assert group.target_field == "chi_0"

    def test_none_target_uses_all_dynamical(self) -> None:
        """target_fields=None auto-selects all other dynamical fields."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        result = compute_group_conversion(data, "phi_0")
        assert result.target_field == "chi_0"

    def test_multi_target_explicit_equals_auto(self) -> None:
        """Explicit target list matches auto-target."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        explicit = compute_group_conversion(data, "phi_0", ["chi_0"])
        auto = compute_group_conversion(data, "phi_0")
        np.testing.assert_allclose(explicit.probability, auto.probability)

    def test_overlap_raises(self) -> None:
        """Source and target groups must not overlap."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="overlap"):
            compute_group_conversion(data, "phi_0", ["phi_0", "chi_0"])

    def test_invalid_field_raises(self) -> None:
        """Invalid field name raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="not in spec"):
            compute_group_conversion(data, "nonexistent")

    def test_all_fields_as_source_empty_target_raises(self) -> None:
        """All dynamical fields as source leaves empty target -> raises."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="empty"):
            compute_group_conversion(data, ["phi_0", "chi_0"])


# ============================================================
# Spatial operator tests
# ============================================================


class TestApplySpatialOperator:
    """Test _apply_spatial_operator dispatch."""

    def test_identity_returns_copy(self) -> None:
        """Identity operator returns a copy, not the same array."""
        field = np.array([1.0, 2.0, 3.0])
        result = _apply_spatial_operator("identity", field, (0.1,), (True,))
        np.testing.assert_array_equal(result, field)
        assert result is not field  # copy, not alias

    def test_gradient_x_linear(self) -> None:
        """Gradient of a linear field f(x) = 2x gives constant 2."""
        n = 64
        dx = 10.0 / n
        x = np.linspace(dx / 2, 10.0 - dx / 2, n)
        field = 2.0 * x  # f(x) = 2x

        result = _apply_spatial_operator("gradient_x", field, (dx,), (False,))
        np.testing.assert_allclose(result, 2.0, atol=1e-10)

    def test_laplacian_x_cosine_periodic(self) -> None:
        """Laplacian of cos(kx) = -k^2 cos(kx) for periodic FFT."""
        n = 128
        domain = 2.0 * np.pi
        dx = domain / n
        x = np.linspace(0, domain - dx, n)
        k = 3.0
        field = np.cos(k * x)

        result = _apply_spatial_operator("laplacian_x", field, (dx,), (True,))
        expected = -(k**2) * np.cos(k * x)
        np.testing.assert_allclose(result, expected, atol=1e-10)

    def test_laplacian_sum_of_axes(self) -> None:
        """Isotropic laplacian = laplacian_x + laplacian_y."""
        n = 32
        domain = 2.0 * np.pi
        dx = domain / n
        x = np.linspace(0, domain - dx, n)
        y = np.linspace(0, domain - dx, n)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        field = np.cos(xx) + np.cos(yy)

        spacing = (dx, dx)
        periodic = (True, True)

        lap_iso = _apply_spatial_operator("laplacian", field, spacing, periodic)
        lap_x = _apply_spatial_operator("laplacian_x", field, spacing, periodic)
        lap_y = _apply_spatial_operator("laplacian_y", field, spacing, periodic)

        np.testing.assert_allclose(lap_iso, lap_x + lap_y, atol=1e-10)

    def test_cross_derivative_xy(self) -> None:
        """Cross derivative of sin(x)*sin(y) = cos(x)*cos(y)."""
        n = 64
        domain = 2.0 * np.pi
        dx = domain / n
        x = np.linspace(0, domain - dx, n)
        y = np.linspace(0, domain - dx, n)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        field = np.sin(xx) * np.sin(yy)

        spacing = (dx, dx)
        periodic = (True, True)

        result = _apply_spatial_operator(
            "cross_derivative_xy", field, spacing, periodic
        )
        expected = np.cos(xx) * np.cos(yy)
        np.testing.assert_allclose(result, expected, atol=1e-8)

    def test_unknown_operator_raises(self) -> None:
        """Unknown operator name raises ValueError."""
        field = np.zeros(10)
        with pytest.raises(ValueError, match="Unknown spatial operator"):
            _apply_spatial_operator("unknown_op", field, (0.1,), (True,))


class TestIsMomentumField:
    """Test _is_momentum_field regex."""

    def test_pi_underscore(self) -> None:
        assert _is_momentum_field("pi_0") is True
        assert _is_momentum_field("pi_1") is True
        assert _is_momentum_field("pi_12") is True

    def test_pi_no_underscore(self) -> None:
        assert _is_momentum_field("pi0") is True
        assert _is_momentum_field("pi1") is True

    def test_regular_fields(self) -> None:
        assert _is_momentum_field("phi_0") is False
        assert _is_momentum_field("A_0") is False
        assert _is_momentum_field("pi_phi") is False  # non-numeric


# ============================================================
# Virial potential tests
# ============================================================


class TestVirialPotential:
    """Test _compute_virial_potential."""

    def test_single_kg_virial_matches_canonical(self) -> None:
        """For single KG field, virial = gradient + mass energy."""
        data = _make_single_field_data(n_snapshots=3)

        v_virial = _compute_virial_potential(data, 0)
        fe = compute_field_energy(
            data.fields["phi_0"][0],
            data.momenta["phi_0"][0],
            mass_squared=float(data.spec.mass_matrix[0][0]),
            grid_spacing=data.grid_spacing,
            periodic=data.periodic,
        )

        # Virial should equal gradient + mass (the total potential energy)
        np.testing.assert_allclose(v_virial, fe.gradient + fe.mass, rtol=1e-6)

    def test_coupled_oscillator_energy_still_conserved(self) -> None:
        """Existing exact coupled-oscillator data still conserves energy."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        _, _, _, total = compute_energy_timeseries(data)

        relative_drift = np.abs(total - total[0]) / total[0]
        assert np.max(relative_drift) < 1e-10

    def test_coupled_scalars_interaction_nonzero(self) -> None:
        """Coupled scalars should have nonzero interaction energy."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        # At some time step, chi should be nonzero → interaction nonzero
        se = compute_system_energy(data, 5)
        # Total = kinetic + virial + constraint. Interaction = virial - self_potentials.
        # For coupled uniform oscillators, the coupling term is nonzero.
        assert se.interaction != 0.0

    def test_single_field_zero_interaction(self) -> None:
        """Single KG field has zero interaction energy."""
        data = _make_single_field_data(n_snapshots=3)
        se = compute_system_energy(data, 0)
        np.testing.assert_allclose(se.interaction, 0.0, atol=1e-12)


# ============================================================
# Operator-aware gradient tests
# ============================================================


class TestSelfGradientAxes:
    """Tests for _self_gradient_axes helper."""

    def test_full_laplacian_returns_none(self) -> None:
        """Equation with isotropic 'laplacian' → None (all axes)."""
        terms = (OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),)
        eq = ComponentEquation(
            field_name="phi", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) is None

    def test_directional_laplacian_y(self) -> None:
        """Equation with only laplacian_y → [1]."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="A_1"),
            OperatorTerm(coefficient=0.5, operator="cross_derivative_xy", field="A_2"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_1"),
        )
        eq = ComponentEquation(
            field_name="A_1", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [1]

    def test_directional_laplacian_x(self) -> None:
        """Equation with only laplacian_x → [0]."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_2"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_2"),
        )
        eq = ComponentEquation(
            field_name="A_2", field_index=1,
            time_derivative_order=2, rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [0]

    def test_both_directional_laplacians(self) -> None:
        """Equation with laplacian_x + laplacian_y → [0, 1]."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi"),
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="phi"),
        )
        eq = ComponentEquation(
            field_name="phi", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [0, 1]

    def test_cross_field_laplacian_ignored(self) -> None:
        """Laplacian on a different field is not counted as self."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="A_1"),
            OperatorTerm(coefficient=0.5, operator="laplacian_x", field="A_2"),
        )
        eq = ComponentEquation(
            field_name="A_1", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        # Only laplacian_y(A_1) is self; laplacian_x(A_2) is cross-field
        assert _self_gradient_axes(eq) == [1]

    def test_no_laplacian_returns_empty(self) -> None:
        """Equation with no laplacian at all → empty list."""
        terms = (
            OperatorTerm(coefficient=-1.0, operator="identity", field="phi"),
        )
        eq = ComponentEquation(
            field_name="phi", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == []


class TestOperatorAwareGradient:
    """Tests for operator-aware per-field gradient in system energy."""

    def test_vector_single_component_zero_interaction(self) -> None:
        """Vector field with directional laplacian: interaction = 0 at t=0.

        Mimics Proca A_1 with laplacian_y only.  When only A_1 is excited,
        interaction should be zero because the per-field gradient only includes
        the y-axis (matching the virial).
        """
        # Build a 2D spec with two fields: A_1 (laplacian_y) + A_2 (laplacian_x)
        a1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="A_1"),
            OperatorTerm(coefficient=0.5, operator="cross_derivative_xy", field="A_2"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_1"),
        )
        a2_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_2"),
            OperatorTerm(coefficient=0.5, operator="cross_derivative_xy", field="A_1"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_2"),
        )
        eq_a1 = ComponentEquation(
            field_name="A_1", field_index=0,
            time_derivative_order=2, rhs_terms=a1_terms,
        )
        eq_a2 = ComponentEquation(
            field_name="A_2", field_index=1,
            time_derivative_order=2, rhs_terms=a2_terms,
        )

        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            equations=(eq_a1, eq_a2),
            component_names=("A_1", "A_2"),
            mass_matrix=((1.0, 0.0), (0.0, 1.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )

        # 2D grid
        grid = CartesianGrid(
            bounds=[(0, np.pi), (0, np.pi)], shape=[16, 16], periodic=True,
        )
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        y = cast("np.ndarray", grid.cell_coords[..., 1])
        gaussian = 0.5 * np.exp(
            -((x - np.pi / 2) ** 2 + (y - np.pi / 2) ** 2) / (2 * 0.5**2)
        )

        # State: A_1 field, A_1 momentum, A_2 field, A_2 momentum
        fields = [
            ScalarField(grid, data=gaussian),  # A_1 field
            ScalarField(grid, data=0.0),        # A_1 momentum
            ScalarField(grid, data=0.0),        # A_2 field
            ScalarField(grid, data=0.0),        # A_2 momentum
        ]
        state = FieldCollection(fields)
        storage = MemoryStorage()
        storage.start_writing(state)
        storage.append(state, time=0.0)

        data = SimulationData.from_storage(
            storage, spec, grid, parameters={"m2": 1.0},
        )

        se = compute_system_energy(data, 0)

        # Key assertion: interaction should be ~0 (not -0.196 as with isotropic)
        np.testing.assert_allclose(se.interaction, 0.0, atol=1e-10)

        # Per-field gradient for A_1 should only include y-axis
        # (roughly half of the full gradient)
        assert se.per_field["A_1"].gradient > 0.0
        assert se.per_field["A_2"].gradient == 0.0  # A_2 is zero everywhere

    def test_scalar_full_laplacian_unchanged(self) -> None:
        """Scalar field with full laplacian: gradient uses all axes (no change)."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="phi"),
        )
        eq = ComponentEquation(
            field_name="phi", field_index=0,
            time_derivative_order=2, rhs_terms=terms,
        )
        spec = EquationSystem(
            n_components=1, dimension=3, spatial_dimension=2,
            equations=(eq,),
            component_names=("phi",),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = CartesianGrid(
            bounds=[(0, np.pi), (0, np.pi)], shape=[16, 16], periodic=True,
        )
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        y = cast("np.ndarray", grid.cell_coords[..., 1])
        gaussian = 0.5 * np.exp(
            -((x - np.pi / 2) ** 2 + (y - np.pi / 2) ** 2) / (2 * 0.5**2)
        )

        fields = [
            ScalarField(grid, data=gaussian),  # phi field
            ScalarField(grid, data=0.0),        # phi momentum
        ]
        state = FieldCollection(fields)
        storage = MemoryStorage()
        storage.start_writing(state)
        storage.append(state, time=0.0)

        data = SimulationData.from_storage(
            storage, spec, grid, parameters={"m2": 1.0},
        )

        se = compute_system_energy(data, 0)

        # Interaction is zero: single scalar field with full laplacian
        np.testing.assert_allclose(se.interaction, 0.0, atol=1e-10)

        # Gradient should include both x and y axes
        # Compare with explicit full gradient computation
        fe_full = compute_field_energy(
            gaussian, None, 1.0,
            data.grid_spacing, data.periodic, gradient_axes=None,
        )
        np.testing.assert_allclose(
            se.per_field["phi"].gradient, fe_full.gradient, rtol=1e-12,
        )


# ============================================================
# Constraint self-energy tests
# ============================================================


def _make_constraint_spec() -> EquationSystem:
    """Build a synthetic spec with one constraint and one dynamical field.

    Mimics a simplified gauge theory: A_0 (constraint) + A_1 (dynamical).
    """
    # A_0: constraint (time_order=0), simplified gauge theory
    a0_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_0"),
        OperatorTerm(
            coefficient=0.5,
            operator="identity",
            field="A_0",
            coefficient_symbolic="Am2",
        ),
    )
    eq_a0 = ComponentEquation(
        field_name="A_0",
        field_index=0,
        time_derivative_order=0,
        rhs_terms=a0_terms,
    )

    # A_1 equation: dynamical (time_order = 2)
    # d2_t(A_1) = laplacian_x(A_1) - 0.5 * identity(A_1)
    a1_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_1"),
        OperatorTerm(
            coefficient=-0.5,
            operator="identity",
            field="A_1",
            coefficient_symbolic="-Am2",
        ),
    )
    eq_a1 = ComponentEquation(
        field_name="A_1",
        field_index=1,
        time_derivative_order=2,
        rhs_terms=a1_terms,
    )

    return EquationSystem(
        n_components=2,
        dimension=2,
        spatial_dimension=1,
        equations=(eq_a0, eq_a1),
        component_names=("A_0", "A_1"),
        mass_matrix=((-0.5, 0.0), (0.0, 0.5)),
        coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
        metadata={"parameters": {"Am2": 0.5}},
        coordinates=("t", "x"),
        mass_matrix_symbolic=(("-Am2", "0"), ("0", "Am2")),
        coupling_matrix_symbolic=(("0", "0"), ("0", "0")),
    )


class TestConstraintSelfEnergy:
    """Test _compute_constraint_self_energy."""

    def test_constraint_energy_negative(self) -> None:
        """Constraint field with nonzero gradient + mass gives negative energy."""
        spec = _make_constraint_spec()
        n = 64
        dx = 10.0 / n
        x = np.linspace(dx / 2, 10.0 - dx / 2, n)

        # A_0 = cos(kx), A_1 = 0
        k = 2.0 * np.pi / 10.0
        a0_field = np.cos(k * x)

        data = SimulationData(
            times=np.array([0.0]),
            fields={"A_0": a0_field[np.newaxis, :], "A_1": np.zeros((1, n))},
            momenta={"A_1": np.zeros((1, n))},
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"Am2": 0.5},
        )

        energy = _compute_constraint_self_energy(data, 0)

        # Should be negative: -1/2 |grad(A_0)|^2 - 1/2 m^2 A_0^2
        assert energy < 0.0

    def test_no_constraints_returns_zero(self) -> None:
        """System without constraints gives zero constraint self-energy."""
        data = _make_single_field_data(n_snapshots=3)
        energy = _compute_constraint_self_energy(data, 0)
        assert energy == 0.0

    def test_constraint_field_not_in_data_skipped(self) -> None:
        """Constraint field not stored in data is silently skipped."""
        spec = _make_constraint_spec()
        n = 32

        # Only A_1 in data, A_0 omitted
        data = SimulationData(
            times=np.array([0.0]),
            fields={"A_1": np.ones((1, n))},
            momenta={"A_1": np.zeros((1, n))},
            grid_spacing=(10.0 / n,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"Am2": 0.5},
        )

        energy = _compute_constraint_self_energy(data, 0)
        assert energy == 0.0


class TestVirialWithConstraints:
    """Test the full energy computation with constraint fields."""

    def test_total_includes_constraint_energy(self) -> None:
        """Total system energy includes constraint self-energy."""
        spec = _make_constraint_spec()
        n = 64
        dx = 10.0 / n
        x = np.linspace(dx / 2, 10.0 - dx / 2, n)

        k = 2.0 * np.pi * 2.0 / 10.0
        a0 = 0.5 * np.cos(k * x)
        a1 = np.cos(k * x)
        pi_a1 = np.zeros(n)

        data = SimulationData(
            times=np.array([0.0]),
            fields={"A_0": a0[np.newaxis, :], "A_1": a1[np.newaxis, :]},
            momenta={"A_1": pi_a1[np.newaxis, :]},
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"Am2": 0.5},
        )

        se = compute_system_energy(data, 0)

        # Constraint self-energy should be negative (g^{00} = -1 sign flip)
        constraint_e = _compute_constraint_self_energy(data, 0)
        assert constraint_e < 0.0

        # Verify total = kinetic + virial + constraint_self
        virial = _compute_virial_potential(data, 0)
        expected_total = se.per_field["A_1"].kinetic + virial + constraint_e
        np.testing.assert_allclose(se.total, expected_total, rtol=1e-10)


class TestScalarVectorEnergyConservation:
    """Regression test: scalar-vector coupling conserves energy after pi_N fix.

    The scalar_vector_coupling.json had wrong pi_N indices (parse index instead
    of global index) which caused spurious -∂_x(∂_tφ) forces in A_1/A_2 equations.
    With corrected indices, the simulation conserves the virial Hamiltonian.
    """

    def test_energy_conserved(self) -> None:
        """Short scalar-vector simulation conserves energy to < 5%."""
        from tidal.symbolic import build_pde_from_json, load_equation_system
        from tidal.utils import normalize_solve_result

        json_path = (
            Path(__file__).resolve().parent.parent
            / "examples"
            / "data"
            / "scalar_vector_coupling.json"
        )
        if not json_path.exists():
            pytest.skip("scalar_vector_coupling.json not found")

        params = {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2}
        spec = load_equation_system(str(json_path))
        pde = build_pde_from_json(str(json_path), parameters=params)

        grid = CartesianGrid([(0, 10), (0, 10)], [32, 32], periodic=True)

        # 7-slot state: phi_0, pi_phi, A_0, A_1, pi_A1, A_2, pi_A2
        x_coords = cast("np.ndarray", grid.cell_coords[..., 0])
        y_coords = cast("np.ndarray", grid.cell_coords[..., 1])
        gaussian = np.exp(-((x_coords - 5.0) ** 2 + (y_coords - 5.0) ** 2) / 2.0)

        state = FieldCollection(
            [
                ScalarField(grid, data=gaussian, label="phi_0"),
                ScalarField(grid, data=0.0, label="pi_phi"),
                ScalarField(grid, data=0.0, label="A_0"),
                ScalarField(grid, data=0.0, label="A_1"),
                ScalarField(grid, data=0.0, label="pi_A1"),
                ScalarField(grid, data=0.0, label="A_2"),
                ScalarField(grid, data=0.0, label="pi_A2"),
            ]
        )

        storage = MemoryStorage()
        result = pde.solve(
            state,
            t_range=2.0,
            dt=0.01,
            scheme="runge-kutta",
            tracker=storage.tracker(0.2),
        )
        normalize_solve_result(result)

        data = SimulationData.from_storage(storage, spec, grid, params)

        energies: list[float] = []
        for t_idx in range(data.n_snapshots):
            se = compute_system_energy(data, t_idx)
            energies.append(se.total)

        e0 = energies[0]
        max_drift = max(abs(e - e0) for e in energies) / max(abs(e0), 1e-12)

        assert max_drift < 0.05, (
            f"Energy drift {max_drift:.2%} exceeds 5% threshold. "
            f"This likely indicates incorrect pi_N indexing in the JSON."
        )


# ============================================================
# Mixing length tests
# ============================================================


def _make_conversion_result(
    times: np.ndarray,
    probability: np.ndarray,
) -> ConversionResult:
    """Build a ConversionResult with synthetic data for testing."""
    return ConversionResult(
        times=times,
        probability=probability,
        source_energy=1.0 - probability,
        target_energy=probability.copy(),
        total_energy=np.ones_like(times),
        relative_energy_error=np.zeros_like(times),
        source_field="source",
        target_field="target",
    )


class TestMixingLength:
    """Tests for compute_mixing_length — spectral peak detection."""

    def test_sin_squared_mixing_length(self) -> None:
        """sin²(2t) oscillates at ω=4 → L_mix = π/4."""
        times = np.linspace(0, 20, 1000)
        prob = np.sin(2.0 * times) ** 2
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv)

        # sin²(2t) = 0.5(1 - cos(4t)), so P(t) oscillates at ω=4
        assert isinstance(result, MixingResult)
        assert abs(result.dominant_frequency - 4.0) / 4.0 < 0.05
        expected_lmix = np.pi / 4.0
        assert abs(result.mixing_length - expected_lmix) / expected_lmix < 0.05
        assert result.max_conversion > 0.95

    def test_damped_oscillation(self) -> None:
        """Non-sin² form: damped oscillation still has a detectable peak."""
        times = np.linspace(0, 10, 500)
        # P(t) = (1 - e^{-t}) * 0.5 * (1 - cos(4t))
        prob = (1.0 - np.exp(-times)) * 0.5 * (1.0 - np.cos(4.0 * times))
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv)

        assert result.mixing_length > 0
        assert result.mixing_length_uncertainty > 0
        assert result.max_conversion > 0

    def test_delayed_onset_spectral(self) -> None:
        """Delayed sin²(2(t-5)) still has ω_dom ≈ 4, same as no delay."""
        times = np.linspace(0, 20, 1000)
        prob = np.where(
            times < 5.0,
            0.0,
            np.sin(2.0 * (times - 5.0)) ** 2,
        )
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv)

        # Spectral method correctly finds ω=4 regardless of delay
        assert abs(result.dominant_frequency - 4.0) / 4.0 < 0.10
        expected_lmix = np.pi / 4.0
        assert abs(result.mixing_length - expected_lmix) / expected_lmix < 0.10

    def test_monotonic_raises(self) -> None:
        """Constant P(t) has no spectral peaks → ValueError."""
        times = np.linspace(0, 10, 100)
        prob = np.full_like(times, 0.5)  # No oscillation → no peaks
        conv = _make_conversion_result(times, prob)

        with pytest.raises(ValueError, match="No spectral peaks"):
            compute_mixing_length(conv)

    def test_too_few_points_raises(self) -> None:
        """Fewer than 3 time points → ValueError."""
        times = np.array([0.0, 1.0])
        prob = np.array([0.0, 0.5])
        conv = _make_conversion_result(times, prob)

        with pytest.raises(ValueError, match="at least 3"):
            compute_mixing_length(conv)

    def test_invalid_min_prominence(self) -> None:
        """min_prominence must be in (0, 1)."""
        times = np.linspace(0, 20, 1000)
        prob = np.sin(2.0 * times) ** 2
        conv = _make_conversion_result(times, prob)

        with pytest.raises(ValueError, match="min_prominence"):
            compute_mixing_length(conv, min_prominence=0.0)
        with pytest.raises(ValueError, match="min_prominence"):
            compute_mixing_length(conv, min_prominence=1.0)

    def test_min_prominence_filtering(self) -> None:
        """Low prominence → more peaks; high prominence → fewer peaks."""
        omega_1, omega_2 = 2.0, 5.0
        times = np.linspace(0, 40, 4000)
        prob = 0.7 * np.cos(omega_1 * times) + 0.3 * np.cos(omega_2 * times)
        # Shift to positive
        prob = prob - prob.min() + 0.01
        conv = _make_conversion_result(times, prob)

        result_low = compute_mixing_length(conv, min_prominence=0.01)
        result_high = compute_mixing_length(conv, min_prominence=0.5)

        # Low prominence should detect more peaks than high prominence
        assert len(result_low.peaks) >= len(result_high.peaks)
        assert len(result_high.peaks) >= 1

    def test_uncertainty_from_hwhm(self) -> None:
        """Long timeseries with clean cos gives sharp peak → small uncertainty."""
        omega_0 = 3.0
        times = np.linspace(0, 100, 10000)
        prob = 0.5 + 0.5 * np.cos(omega_0 * times)
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv)

        # dL = (pi/omega^2) * HWHM where HWHM = FWHM/2
        expected_unc = (np.pi / (result.dominant_frequency**2)) * result.frequency_fwhm / 2.0
        assert abs(result.mixing_length_uncertainty - expected_unc) < 1e-10
        # Clean signal over long time → narrow peak → small uncertainty
        assert result.mixing_length_uncertainty < result.mixing_length * 0.1

    def test_multiple_peaks_detected(self) -> None:
        """Two-frequency signal detects both peaks."""
        omega_1, omega_2 = 2.0, 5.0
        times = np.linspace(0, 40, 4000)
        prob = 0.7 * np.cos(omega_1 * times) + 0.3 * np.cos(omega_2 * times)
        prob = prob - prob.min() + 0.01
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv, min_prominence=0.01)

        assert len(result.peaks) >= 2
        # Dominant peak should be near omega_1 (stronger component)
        assert abs(result.peaks[0].frequency - omega_1) / omega_1 < 0.05

    def test_peaks_sorted_by_power(self) -> None:
        """Detected peaks are sorted by power descending."""
        omega_1, omega_2 = 2.0, 5.0
        times = np.linspace(0, 40, 4000)
        prob = 0.7 * np.cos(omega_1 * times) + 0.3 * np.cos(omega_2 * times)
        prob = prob - prob.min() + 0.01
        conv = _make_conversion_result(times, prob)

        result = compute_mixing_length(conv, min_prominence=0.01)

        for i in range(len(result.peaks) - 1):
            assert result.peaks[i].power >= result.peaks[i + 1].power

    def test_rayleigh_resolution_stored(self) -> None:
        """MixingSpectrum stores Rayleigh resolution = 2*pi/T."""
        times = np.linspace(0, 20, 200)
        prob = 0.5 + 0.5 * np.cos(3.0 * times)
        conv = _make_conversion_result(times, prob)

        spectrum = compute_mixing_spectrum(conv)

        t_total = float(times[-1] - times[0])
        expected = 2.0 * np.pi / t_total
        assert abs(spectrum.rayleigh_resolution - expected) / expected < 0.01


class TestMixingSpectrum:
    """Tests for compute_mixing_spectrum — temporal FFT of P(t)."""

    def test_single_oscillation(self) -> None:
        """cos(ω₀ t) has dominant peak at ω₀."""
        omega_0 = 3.0
        times = np.linspace(0, 20, 2000)
        prob = 0.5 + 0.5 * np.cos(omega_0 * times)
        conv = _make_conversion_result(times, prob)

        spectrum = compute_mixing_spectrum(conv)

        assert isinstance(spectrum, MixingSpectrum)
        # Dominant frequency should be close to omega_0
        assert abs(spectrum.dominant_frequency - omega_0) / omega_0 < 0.05
        # Mixing length = π/ω
        expected_lmix = np.pi / omega_0
        assert abs(spectrum.dominant_mixing_length - expected_lmix) / expected_lmix < 0.05

    def test_two_frequencies(self) -> None:
        """Two-frequency signal shows two peaks, dominant is the stronger one."""
        omega_1, omega_2 = 2.0, 5.0
        times = np.linspace(0, 40, 4000)
        prob = 0.7 * np.cos(omega_1 * times) + 0.3 * np.cos(omega_2 * times)
        conv = _make_conversion_result(times, prob)

        spectrum = compute_mixing_spectrum(conv)

        # Dominant should be omega_1 (stronger component)
        assert abs(spectrum.dominant_frequency - omega_1) / omega_1 < 0.05

        # Both peaks should be visible in the spectrum
        # Find the bin closest to omega_2
        idx_2 = int(np.argmin(np.abs(spectrum.frequencies - omega_2)))
        # Power at omega_2 should be substantial (not zero)
        assert spectrum.power[idx_2] > 0.01 * spectrum.power.max()

    def test_dc_excluded(self) -> None:
        """DC bin is excluded — first frequency is positive."""
        times = np.linspace(0, 10, 200)
        prob = 0.5 + 0.5 * np.cos(2.0 * times)
        conv = _make_conversion_result(times, prob)

        spectrum = compute_mixing_spectrum(conv)

        assert spectrum.frequencies[0] > 0.0

    def test_non_uniform_timestep_raises(self) -> None:
        """Non-uniform time spacing → ValueError."""
        times = np.array([0.0, 1.0, 3.0, 6.0, 10.0])
        prob = np.array([0.0, 0.5, 0.3, 0.8, 0.2])
        conv = _make_conversion_result(times, prob)

        with pytest.raises(ValueError, match="Non-uniform"):
            compute_mixing_spectrum(conv)

    def test_consistent_with_mixing_length(self) -> None:
        """compute_mixing_length derives L_mix from the same spectrum."""
        omega = 2.0
        times = np.linspace(0, 20, 2000)
        prob = np.sin(omega * times) ** 2
        conv = _make_conversion_result(times, prob)

        mixing_result = compute_mixing_length(conv)
        spectrum = compute_mixing_spectrum(conv)

        # Both use the same FFT; mixing_length = π/ω_dom from the dominant peak
        assert abs(mixing_result.mixing_length - spectrum.dominant_mixing_length) < 0.01
        assert abs(mixing_result.dominant_frequency - spectrum.dominant_frequency) < 0.01


# ============================================================
# Spectral conversion tests
# ============================================================


class TestSpectralConversion:
    """Test compute_spectral_conversion for per-mode P(k,t)."""

    def test_shapes(self) -> None:
        """Output arrays have correct shapes."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        result = compute_spectral_conversion(data, "phi_0", "chi_0")

        assert isinstance(result, SpectralConversion)
        n_snap = data.n_snapshots
        n_modes = len(result.wavenumbers)
        assert result.probability.shape == (n_snap, n_modes)
        assert result.source_spectral_energy.shape == (n_snap, n_modes)
        assert result.target_spectral_energy.shape == (n_snap, n_modes)
        assert result.active_modes.shape == (n_modes,)
        assert len(result.times) == n_snap

    def test_field_names(self) -> None:
        """Source and target field names are stored correctly."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        result = compute_spectral_conversion(data, "phi_0", "chi_0")
        assert result.source_field == "phi_0"
        assert result.target_field == "chi_0"

    def test_same_field_raises(self) -> None:
        """Same source and target raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="different"):
            compute_spectral_conversion(data, "phi_0", "phi_0")

    def test_invalid_field_raises(self) -> None:
        """Invalid field name raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="not in spec"):
            compute_spectral_conversion(data, "phi_0", "nonexistent")

    def test_zero_source_energy_raises(self) -> None:
        """Zero source energy raises ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=5, amplitude=0.0)
        with pytest.raises(ValueError, match="no modes above energy floor"):
            compute_spectral_conversion(data, "phi_0", "chi_0")

    def test_zero_target_gives_zero_probability(self) -> None:
        """Target at zero → P(k,t=0) = 0 everywhere."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        result = compute_spectral_conversion(data, "phi_0", "chi_0")
        # chi starts at zero, so P(k, t=0) should be ~0
        np.testing.assert_allclose(result.probability[0], 0.0, atol=1e-10)

    def test_energy_floor_masks_inactive(self) -> None:
        """Inactive modes (below energy floor) have P=0 and active_modes=False."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        # Very high floor raises (all modes below threshold)
        with pytest.raises(ValueError, match="no modes above energy floor"):
            compute_spectral_conversion(
                data, "phi_0", "chi_0", energy_floor=1e10,
            )
        # Normal case: inactive modes have P=0
        result = compute_spectral_conversion(data, "phi_0", "chi_0")
        inactive = ~result.active_modes
        if np.any(inactive):
            assert np.all(result.probability[:, inactive] == 0.0)

    def test_p_nonzero_at_later_times(self) -> None:
        """P(k,t) should grow above zero as coupling transfers energy."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        result = compute_spectral_conversion(data, "phi_0", "chi_0")
        # At final time, there should be nonzero conversion in active modes
        active_p_final = result.probability[-1, result.active_modes]
        assert np.max(active_p_final) > 0.01


class TestGroupSpectralConversion:
    """Test compute_group_spectral_conversion."""

    def test_single_field_matches_pairwise(self) -> None:
        """Single-field groups should match pairwise spectral conversion."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        pairwise = compute_spectral_conversion(data, "phi_0", "chi_0")
        group = compute_group_spectral_conversion(data, "phi_0", "chi_0")
        np.testing.assert_allclose(
            group.probability, pairwise.probability, atol=1e-15,
        )

    def test_overlap_raises(self) -> None:
        """Source and target groups must not overlap."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="overlap"):
            compute_group_spectral_conversion(
                data, "phi_0", ["phi_0", "chi_0"],
            )

    def test_none_target_auto_selects(self) -> None:
        """target_fields=None auto-selects remaining dynamical fields."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        result = compute_group_spectral_conversion(data, "phi_0")
        assert result.target_field == "chi_0"


# ============================================================
# Review Pass 4 tests (A8)
# ============================================================


class TestResolveTermTarget:
    """Tests for _resolve_term_target fail-fast behavior (A3)."""

    def test_unknown_field_raises(self) -> None:
        """Unresolvable field reference should raise ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        with pytest.raises(ValueError, match="Unresolvable field reference"):
            _resolve_term_target(data, "nonexistent_field", 0)

    def test_constraint_returns_none(self) -> None:
        """Constraint field momentum (time_derivative_order == 0) should return None."""
        import dataclasses

        # Build data with a constraint field (order 0)
        spec = _build_coupled_scalars_spec()
        # Modify first equation to be a constraint (order 0)
        eq0 = spec.equations[0]
        constraint_eq = dataclasses.replace(eq0, time_derivative_order=0)
        constraint_spec = dataclasses.replace(
            spec, equations=(constraint_eq, spec.equations[1]),
        )

        data = _make_sim_data_two_fields(n_snapshots=3)
        # Replace spec with the constraint version
        constraint_data = dataclasses.replace(data, spec=constraint_spec)
        # pi_0 refers to the first field — which is now a constraint
        result = _resolve_term_target(constraint_data, "pi_0", 0)
        assert result is None

    def test_momentum_out_of_range_raises(self) -> None:
        """Momentum index beyond spec fields should raise ValueError."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        # pi_99 is out of range for a 2-field spec
        with pytest.raises(ValueError, match="resolves to index 99"):
            _resolve_term_target(data, "pi_99", 0)


class TestBincountRegression:
    """Verify np.bincount radial binning matches expected behavior (A4)."""

    def test_bincount_sums_correctly(self) -> None:
        """Radial binning should sum values into correct bins."""
        from tidal.measurement._spectral import (
            _radial_bin,  # pyright: ignore[reportPrivateUsage]
        )

        # Build a simple 2D k-magnitude grid and values
        grid_spacing = (1.0, 1.0)
        field_shape = (8, 8)
        k_arrays = [
            np.fft.fftfreq(8, d=1.0) * 2 * np.pi,
            np.fft.rfftfreq(8, d=1.0) * 2 * np.pi,
        ]
        k_grid = np.meshgrid(*k_arrays, indexing="ij")
        k_mag = np.sqrt(k_grid[0] ** 2 + k_grid[1] ** 2)

        # Uniform values -- sum in each bin should equal count * value
        values = np.ones_like(k_mag)
        centers, binned = _radial_bin(k_mag, values, grid_spacing, field_shape)

        # Total sum must be preserved
        np.testing.assert_allclose(binned.sum(), values.sum(), rtol=1e-12)
        # All bins non-negative
        assert np.all(binned >= 0)
        # Must have at least one bin
        assert len(centers) >= 1
        assert len(centers) == len(binned)


class TestSpectralEnergyPhysics:
    """Tests for spectral energy normalization (A6) and edge cases."""

    def test_spectral_energy_matches_field_energy_uniform(self) -> None:
        """Spectral energy should exactly match real-space energy for uniform field.

        For a uniform field (all power in k=0), rfftn returns the DC bin fully
        (no conjugate-symmetry issue), so the Parseval relationship holds exactly.
        After the dV fix (A6), sum(spectral_bins) == field_energy.
        """
        n_grid = 32
        dx = 10.0 / n_grid
        amplitude = 2.0
        field_arr = np.full(n_grid, amplitude)
        mom_arr = np.zeros(n_grid)
        m2 = 1.0  # arbitrary positive mass squared

        # Real-space energy: purely mass term (gradient=0 for uniform, kinetic=0)
        fe = compute_field_energy(field_arr, mom_arr, m2, (dx,), (True,))

        # Spectral energy
        _wn, se_bins = compute_spectral_energy(
            field_arr, mom_arr, m2, (dx,), (True,),
        )
        spectral_total = float(se_bins.sum())

        # Exact match for uniform field (k=0 only)
        np.testing.assert_allclose(spectral_total, fe.total, rtol=1e-10)

    def test_spectral_energy_negative_mass(self) -> None:
        """Tachyonic m² < 0 should not crash; energy can be negative per mode."""
        n_grid = 32
        dx = 1.0
        field = np.random.default_rng(42).standard_normal(n_grid)
        momentum = np.random.default_rng(43).standard_normal(n_grid)
        m2_tachyonic = -1.0

        wavenumbers, se = compute_spectral_energy(
            field, momentum, m2_tachyonic, (dx,), (True,),
        )
        # Should produce finite results (some bins may be negative)
        assert np.all(np.isfinite(se))
        assert len(wavenumbers) == len(se)


# ============================================================
# Dispersion relation tests (Part B)
# ============================================================


def _make_plane_wave_data(
    n_grid: int = 64,
    n_snapshots: int = 64,
    mode_index: int = 3,
) -> SimulationData:
    """Build synthetic data: single-mode plane wave with exact KG evolution.

    Creates phi(x, t) = A * cos(k0*x) * cos(omega0*t) where
    omega0 = sqrt(k0^2 + m^2) from the KG dispersion relation.
    """
    spec = load_equation_system(DATA_DIR / "klein_gordon_1d.json")
    m2 = float(spec.mass_matrix[0][0])

    domain_len = 20.0
    dx = domain_len / n_grid
    x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

    k0 = 2.0 * np.pi * mode_index / domain_len
    omega0 = float(np.sqrt(k0**2 + m2))

    # Time span: long enough for several oscillations
    t_end = 8.0 * np.pi / omega0  # ~4 full oscillations
    times = np.linspace(0.0, t_end, n_snapshots)

    amplitude = 1.0
    fields_list: list[np.ndarray] = []
    momenta_list: list[np.ndarray] = []

    for t in times:
        phi = amplitude * np.cos(k0 * x) * np.cos(omega0 * t)
        pi_field = -amplitude * omega0 * np.cos(k0 * x) * np.sin(omega0 * t)
        fields_list.append(phi)
        momenta_list.append(pi_field)

    return SimulationData(
        times=times,
        fields={"phi_0": np.stack(fields_list)},
        momenta={"phi_0": np.stack(momenta_list)},
        grid_spacing=(dx,),
        grid_bounds=((0.0, domain_len),),
        periodic=(True,),
        spec=spec,
        parameters={},
    )


class TestDispersionRelation:
    """Tests for compute_dispersion."""

    def test_shapes(self) -> None:
        """Output arrays should have consistent shapes."""
        data = _make_plane_wave_data()
        result = compute_dispersion(data, "phi_0")

        assert isinstance(result, DispersionResult)
        n_modes = len(result.wavenumbers)
        n_freq = len(result.frequencies)
        assert result.power.shape == (n_modes, n_freq)
        assert result.peak_frequencies.shape == (n_modes,)
        assert result.peak_powers.shape == (n_modes,)

    def test_field_name_stored(self) -> None:
        """field_name should round-trip."""
        data = _make_plane_wave_data()
        result = compute_dispersion(data, "phi_0")
        assert result.field_name == "phi_0"

    def test_invalid_field_raises(self) -> None:
        """Unknown field should raise ValueError."""
        data = _make_plane_wave_data()
        with pytest.raises(ValueError, match="not in spec"):
            compute_dispersion(data, "nonexistent")

    def test_too_few_snapshots_raises(self) -> None:
        """Fewer than 3 snapshots should raise ValueError."""
        data = _make_plane_wave_data(n_snapshots=2)
        with pytest.raises(ValueError, match="at least 3"):
            compute_dispersion(data, "phi_0")

    def test_non_uniform_timestep_raises(self) -> None:
        """Non-uniform timestep should raise ValueError."""
        import dataclasses

        data = _make_plane_wave_data(n_snapshots=10)
        # Make timestep non-uniform
        bad_times = data.times.copy()
        bad_times[-1] += 1.0
        bad_data = dataclasses.replace(data, times=bad_times)
        with pytest.raises(ValueError, match="Non-uniform"):
            compute_dispersion(bad_data, "phi_0")

    def test_peak_detection_plane_wave(self) -> None:
        """Monochromatic plane wave should peak at the injected (k0, omega0)."""
        mode_index = 3
        data = _make_plane_wave_data(mode_index=mode_index)
        spec = data.spec
        m2 = float(spec.mass_matrix[0][0])
        domain_len = float(data.grid_bounds[0][1] - data.grid_bounds[0][0])

        k0 = 2.0 * np.pi * mode_index / domain_len
        omega0 = float(np.sqrt(k0**2 + m2))

        result = compute_dispersion(data, "phi_0")

        # Find the k-bin closest to k0
        k_idx = int(np.argmin(np.abs(result.wavenumbers - k0)))
        detected_omega = result.peak_frequencies[k_idx]

        # Should match within Rayleigh resolution
        np.testing.assert_allclose(
            detected_omega, omega0, atol=2 * result.rayleigh_resolution,
        )

    def test_kg_dispersion_curve(self) -> None:
        """Free KG field: active modes should follow omega = sqrt(k^2 + m^2).

        Uses parameter variables from the spec, not hardcoded numbers.
        """
        data = _make_plane_wave_data(n_grid=64, n_snapshots=128)
        m2 = float(data.spec.mass_matrix[0][0])
        result = compute_dispersion(data, "phi_0")

        # For active modes, compare detected omega to analytical
        active = result.peak_frequencies > 0.0
        if np.any(active):
            k_active = result.wavenumbers[active]
            omega_detected = result.peak_frequencies[active]
            omega_expected = np.sqrt(k_active**2 + m2)

            # Within 2x Rayleigh resolution
            np.testing.assert_allclose(
                omega_detected, omega_expected,
                atol=2 * result.rayleigh_resolution,
            )

    def test_rayleigh_resolution(self) -> None:
        """Rayleigh resolution should be 2*pi / T."""
        data = _make_plane_wave_data(n_snapshots=32)
        result = compute_dispersion(data, "phi_0")

        t_span = float(data.times[-1] - data.times[0])
        expected = 2.0 * np.pi / t_span
        np.testing.assert_allclose(result.rayleigh_resolution, expected, rtol=1e-10)

    def test_inactive_modes_zeroed(self) -> None:
        """Modes below min_amplitude should have peak_frequencies = 0."""
        data = _make_plane_wave_data()
        result = compute_dispersion(data, "phi_0", min_amplitude=1e-12)

        # Plane wave only excites one mode — many modes should be inactive
        n_active = np.count_nonzero(result.peak_frequencies > 0.0)
        # At most a handful of modes should be active (the injected mode
        # plus possibly its neighbors due to spectral leakage)
        assert n_active < len(result.wavenumbers)
        # Inactive modes should have zero peak frequency
        inactive = result.peak_frequencies == 0.0
        assert np.all(result.peak_powers[inactive] == 0.0)

    def test_power_nonnegative(self) -> None:
        """All S(k, omega) values should be >= 0."""
        data = _make_plane_wave_data()
        result = compute_dispersion(data, "phi_0")
        assert np.all(result.power >= 0.0)

    def test_single_frequency_dominates(self) -> None:
        """Monochromatic wave should concentrate power at one frequency bin."""
        data = _make_plane_wave_data(mode_index=3)
        result = compute_dispersion(data, "phi_0")

        # Find mode with most power
        total_per_mode = result.power.sum(axis=1)
        dominant_k = int(np.argmax(total_per_mode))

        # The peak power should be > 50% of the total for that mode
        mode_power = result.power[dominant_k]
        assert float(np.max(mode_power)) > 0.5 * float(np.sum(mode_power))

    def test_frequencies_positive(self) -> None:
        """All frequency values should be > 0 (DC excluded)."""
        data = _make_plane_wave_data()
        result = compute_dispersion(data, "phi_0")
        assert np.all(result.frequencies > 0.0)


# ============================================================
# Group 12 — SnapshotWriter and disk-backed storage
# ============================================================

from tidal.measurement._writer import (  # noqa: E402
    SnapshotWriter,
    compute_snapshot_count,
)


class TestSnapshotWriter:
    """Tests for SnapshotWriter disk-backed streaming."""

    def test_round_trip(self, tmp_path: Path) -> None:
        """Write snapshots via SnapshotWriter, read via from_directory."""
        spec = _build_coupled_scalars_spec()
        n_grid = 16
        n_snapshots = 5
        grid_shape = (n_grid,)
        grid_spacing = (0.5,)
        grid_bounds = ((0.0, 8.0),)
        periodic = (True,)
        params = {"m_phi": 1.0, "m_chi": 2.0, "g": 0.1}

        # Generate synthetic data
        rng = np.random.default_rng(42)
        expected_times = np.linspace(0.0, 4.0, n_snapshots)
        expected_fields: dict[str, np.ndarray] = {
            "phi_0": rng.standard_normal((n_snapshots, n_grid)),
            "chi_0": rng.standard_normal((n_snapshots, n_grid)),
        }
        expected_momenta: dict[str, np.ndarray] = {
            "phi_0": rng.standard_normal((n_snapshots, n_grid)),
            "chi_0": rng.standard_normal((n_snapshots, n_grid)),
        }

        output_dir = tmp_path / "snap_out"
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=grid_shape,
            n_snapshots=n_snapshots,
            grid_spacing=grid_spacing,
            grid_bounds=grid_bounds,
            periodic=periodic,
            parameters=params,
        ) as writer:
            for i in range(n_snapshots):
                writer.append(
                    float(expected_times[i]),
                    {
                        "phi_0": expected_fields["phi_0"][i],
                        "chi_0": expected_fields["chi_0"][i],
                    },
                    {
                        "phi_0": expected_momenta["phi_0"][i],
                        "chi_0": expected_momenta["chi_0"][i],
                    },
                )

        # Load via from_directory
        data = SimulationData.from_directory(output_dir, spec)

        np.testing.assert_allclose(data.times, expected_times)
        for name in ("phi_0", "chi_0"):
            np.testing.assert_allclose(data.fields[name], expected_fields[name])
            np.testing.assert_allclose(data.momenta[name], expected_momenta[name])

        assert data.grid_spacing == grid_spacing
        assert data.grid_bounds == grid_bounds
        assert data.periodic == periodic
        assert data.parameters == params

    def test_mmap_readonly(self, tmp_path: Path) -> None:
        """from_directory returns memory-mapped arrays (read-only)."""
        spec = _build_coupled_scalars_spec()
        n_grid = 8
        n_snapshots = 3

        output_dir = tmp_path / "mmap_test"
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=n_snapshots,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 8.0),),
            periodic=(True,),
        ) as writer:
            for i in range(n_snapshots):
                writer.append(
                    float(i),
                    {"phi_0": np.zeros(n_grid), "chi_0": np.ones(n_grid)},
                    {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
                )

        data = SimulationData.from_directory(output_dir, spec)

        # Memory-mapped arrays should be np.memmap instances
        assert isinstance(data.fields["phi_0"], np.memmap)
        assert isinstance(data.momenta["phi_0"], np.memmap)
        assert isinstance(data.times, np.memmap)

    def test_count_property(self, tmp_path: Path) -> None:
        """Count tracks how many snapshots have been written."""
        output_dir = tmp_path / "count_test"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=3,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        assert writer.count == 0
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        assert writer.count == 1
        writer.append(1.0, {"phi_0": np.zeros(4)}, {})
        assert writer.count == 2
        writer.close()
        assert writer.count == 2

    def test_overflow_raises(self, tmp_path: Path) -> None:
        """Writing more snapshots than pre-allocated raises ValueError."""
        output_dir = tmp_path / "overflow_test"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=1,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        with pytest.raises(ValueError, match="Cannot write snapshot"):
            writer.append(1.0, {"phi_0": np.zeros(4)}, {})
        writer.close()

    def test_missing_field_raises(self, tmp_path: Path) -> None:
        """Append raises ValueError if a required field is missing."""
        output_dir = tmp_path / "missing_field"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=2,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        with pytest.raises(ValueError, match="Missing field 'chi_0'"):
            writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        writer.close()

    def test_write_after_close_raises(self, tmp_path: Path) -> None:
        """Writing to a closed writer raises ValueError."""
        output_dir = tmp_path / "closed_test"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=2,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        writer.close()
        with pytest.raises(ValueError, match="already closed"):
            writer.append(0.0, {"phi_0": np.zeros(4)}, {})

    def test_zero_snapshots_raises(self) -> None:
        """n_snapshots < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_snapshots must be >= 1"):
            SnapshotWriter(
                output_dir=Path("/tmp/unused"),
                field_names=["phi_0"],
                momentum_names=[],
                grid_shape=(4,),
                n_snapshots=0,
                grid_spacing=(1.0,),
                grid_bounds=((0.0, 4.0),),
                periodic=(False,),
            )

    def test_empty_field_names_raises(self) -> None:
        """Empty field_names raises ValueError."""
        with pytest.raises(ValueError, match="field_names must be non-empty"):
            SnapshotWriter(
                output_dir=Path("/tmp/unused"),
                field_names=[],
                momentum_names=[],
                grid_shape=(4,),
                n_snapshots=1,
                grid_spacing=(1.0,),
                grid_bounds=((0.0, 4.0),),
                periodic=(False,),
            )

    def test_single_snapshot(self, tmp_path: Path) -> None:
        """Single-snapshot round-trip works correctly."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "single"
        n_grid = 4
        field_data = np.array([1.0, 2.0, 3.0, 4.0])

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=1,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
        ) as writer:
            writer.append(
                0.0,
                {"phi_0": field_data, "chi_0": np.zeros(n_grid)},
                {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
            )

        data = SimulationData.from_directory(output_dir, spec)
        assert data.n_snapshots == 1
        np.testing.assert_allclose(data.fields["phi_0"][0], field_data)

    def test_partial_close_warns(self, tmp_path: Path) -> None:
        """Closing with fewer snapshots than pre-allocated emits a warning."""
        output_dir = tmp_path / "partial"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=5,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        writer.append(1.0, {"phi_0": np.ones(4)}, {})

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            writer.close()
            assert len(w) == 1
            assert "2/5" in str(w[0].message)

    def test_2d_grid(self, tmp_path: Path) -> None:
        """Round-trip with a 2D grid shape."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "grid2d"
        grid_shape = (8, 6)
        n_snapshots = 3

        rng = np.random.default_rng(123)
        expected = rng.standard_normal((n_snapshots, *grid_shape))

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=grid_shape,
            n_snapshots=n_snapshots,
            grid_spacing=(1.0, 1.0),
            grid_bounds=((0.0, 8.0), (0.0, 6.0)),
            periodic=(True, True),
        ) as writer:
            for i in range(n_snapshots):
                writer.append(
                    float(i),
                    {"phi_0": expected[i], "chi_0": np.zeros(grid_shape)},
                    {"phi_0": np.zeros(grid_shape), "chi_0": np.zeros(grid_shape)},
                )

        data = SimulationData.from_directory(output_dir, spec)
        np.testing.assert_allclose(data.fields["phi_0"], expected)

    def test_spec_path_recorded(self, tmp_path: Path) -> None:
        """spec_path is recorded in metadata.json."""
        import json

        output_dir = tmp_path / "specpath"
        spec_path = Path("examples/data/coupled_scalars.json")

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=1,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
            spec_path=spec_path,
        ) as writer:
            writer.append(0.0, {"phi_0": np.zeros(4)}, {})

        metadata = json.loads((output_dir / "metadata.json").read_text())
        assert metadata["spec_path"] == str(spec_path)
        assert metadata["version"] == 1


class TestCrashRecovery:
    """Test recovery from incomplete writes (no metadata.json)."""

    def test_recovery_without_metadata(self, tmp_path: Path) -> None:
        """from_directory recovers when metadata.json is missing (crash)."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "crash"
        n_grid = 4
        n_snapshots = 5

        # Write 3 of 5 snapshots, then DON'T close (simulates crash)
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=n_snapshots,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
        )
        for i in range(3):
            writer.append(
                float(i) + 0.1,  # non-zero times
                {"phi_0": np.full(n_grid, float(i)), "chi_0": np.zeros(n_grid)},
                {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
            )
        # Flush without writing metadata.json
        writer._flush_mmaps()  # pyright: ignore[reportPrivateUsage]
        # Delete metadata.json if it exists (it shouldn't, but be safe)
        meta_path = output_dir / "metadata.json"
        if meta_path.exists():
            meta_path.unlink()

        # from_directory should infer 3 snapshots from times.npy
        data = SimulationData.from_directory(output_dir, spec)
        assert data.n_snapshots == 3
        np.testing.assert_allclose(data.times[0], 0.1)
        np.testing.assert_allclose(data.times[2], 2.1)
        np.testing.assert_allclose(data.fields["phi_0"][1], np.full(n_grid, 1.0))

    def test_recovery_all_zero_times(self, tmp_path: Path) -> None:
        """Recovery with only t=0 written yields 1 snapshot."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "crash_t0"
        n_grid = 4
        n_snapshots = 10

        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=n_snapshots,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
        )
        # Write only the initial snapshot at t=0
        writer.append(
            0.0,
            {"phi_0": np.ones(n_grid), "chi_0": np.zeros(n_grid)},
            {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
        )
        writer._flush_mmaps()  # pyright: ignore[reportPrivateUsage]

        # Remove metadata.json
        meta_path = output_dir / "metadata.json"
        if meta_path.exists():
            meta_path.unlink()

        data = SimulationData.from_directory(output_dir, spec)
        assert data.n_snapshots == 1
        np.testing.assert_allclose(data.fields["phi_0"][0], np.ones(n_grid))


class TestSimulationDataLoad:
    """Test SimulationData.load() universal entry point."""

    def test_load_npz(self, tmp_path: Path) -> None:
        """load() correctly dispatches to from_npz_auto for .npz files."""
        data_orig = _make_single_field_data()

        # Save as NPZ with grid metadata
        npz_path = tmp_path / "test.npz"
        save_data: dict[str, np.ndarray] = {"times": data_orig.times}
        for name in data_orig.fields:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"{name}_t{t_idx}"] = data_orig.fields[name][t_idx]
        for name in data_orig.momenta:
            for t_idx in range(data_orig.n_snapshots):
                save_data[f"pi_{name}_t{t_idx}"] = data_orig.momenta[name][t_idx]

        # Grid metadata required by from_npz_auto
        save_data["grid_spacing"] = np.array(data_orig.grid_spacing)
        save_data["grid_bounds"] = np.array(data_orig.grid_bounds)
        save_data["grid_periodic"] = np.array(data_orig.periodic)
        np.savez(str(npz_path), **save_data)  # type: ignore[reportArgumentType]

        data_loaded = SimulationData.load(npz_path, data_orig.spec)
        np.testing.assert_allclose(data_loaded.times, data_orig.times)
        np.testing.assert_allclose(
            data_loaded.fields["phi_0"], data_orig.fields["phi_0"],
        )

    def test_load_directory(self, tmp_path: Path) -> None:
        """load() correctly dispatches to from_directory for directories."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "loaddir"
        n_grid = 4
        n_snapshots = 2

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=n_snapshots,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
        ) as writer:
            for i in range(n_snapshots):
                writer.append(
                    float(i),
                    {"phi_0": np.full(n_grid, float(i)), "chi_0": np.zeros(n_grid)},
                    {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
                )

        data = SimulationData.load(output_dir, spec)
        assert data.n_snapshots == n_snapshots
        np.testing.assert_allclose(data.fields["phi_0"][1], np.full(n_grid, 1.0))


class TestComputeSnapshotCount:
    """Tests for compute_snapshot_count utility."""

    def test_exact_count(self) -> None:
        """10.0 / 0.1 + 1 = 101 snapshots."""
        assert compute_snapshot_count(10.0, 0.1) == 101

    def test_non_divisible(self) -> None:
        """Non-divisible interval truncates via int()."""
        # 10.0 / 3.0 = 3.333... → int(3.333) + 1 = 4
        assert compute_snapshot_count(10.0, 3.0) == 4

    def test_single_snapshot_interval(self) -> None:
        """Interval equal to t_end → 2 snapshots (t=0 and t=t_end)."""
        assert compute_snapshot_count(5.0, 5.0) == 2

    def test_zero_t_end_raises(self) -> None:
        with pytest.raises(ValueError, match="t_end must be positive"):
            compute_snapshot_count(0.0, 1.0)

    def test_negative_interval_raises(self) -> None:
        with pytest.raises(ValueError, match="snapshot_interval must be positive"):
            compute_snapshot_count(10.0, -1.0)
