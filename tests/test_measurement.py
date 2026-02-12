"""Tests for the tidal.measurement module.

Covers energy computation, conversion probability, spectral analysis,
and diagnostics for coupled field systems.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from tidal.measurement import (
    ConversionResult,
    EnergyDiagnostics,
    SimulationData,
    SystemEnergy,
    check_energy_conservation,
    compute_conversion_probability,
    compute_energy_timeseries,
    compute_field_energy,
    compute_group_conversion,
    compute_mode_amplitudes,
    compute_spectrum,
    compute_system_energy,
)
from tidal.measurement._energy import (
    _apply_spatial_operator,  # pyright: ignore[reportPrivateUsage]
    _compute_constraint_self_energy,  # pyright: ignore[reportPrivateUsage]
    _compute_virial_potential,  # pyright: ignore[reportPrivateUsage]
    _gradient_energy_density,  # pyright: ignore[reportPrivateUsage]
    _is_momentum_field,  # pyright: ignore[reportPrivateUsage]
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

    return SimulationData(
        times=times,
        fields=fields_np,
        momenta=momenta_np,
        grid_spacing=grid_spacing,
        grid_bounds=grid_bounds,
        periodic=periodic,
        spec=spec,
        parameters={},
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

        result = _apply_spatial_operator("cross_derivative_xy", field, spacing, periodic)
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
# Constraint self-energy tests
# ============================================================


def _make_constraint_spec() -> EquationSystem:
    """Build a synthetic spec with one constraint and one dynamical field.

    Mimics a simplified gauge theory: A_0 (constraint) + A_1 (dynamical).
    """
    # A_0 equation: constraint (time_order = 0)
    # A_0 = laplacian_x(A_0) + 0.5 * identity(A_0)
    a0_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_0"),
        OperatorTerm(coefficient=0.5, operator="identity", field="A_0",
                     coefficient_symbolic="Am2"),
    )
    eq_a0 = ComponentEquation(
        field_name="A_0", field_index=0,
        time_derivative_order=0, rhs_terms=a0_terms,
    )

    # A_1 equation: dynamical (time_order = 2)
    # d2_t(A_1) = laplacian_x(A_1) - 0.5 * identity(A_1)
    a1_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_1"),
        OperatorTerm(coefficient=-0.5, operator="identity", field="A_1",
                     coefficient_symbolic="-Am2"),
    )
    eq_a1 = ComponentEquation(
        field_name="A_1", field_index=1,
        time_derivative_order=2, rhs_terms=a1_terms,
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

        # Total = kinetic + virial + constraint_self
        # The interaction captures the difference between virial total
        # and per-field self-potentials (including constraint contribution)
        a1_self = se.per_field["A_1"].gradient + se.per_field["A_1"].mass
        virial = _compute_virial_potential(data, 0)
        expected_total = se.per_field["A_1"].kinetic + virial + constraint_e
        np.testing.assert_allclose(se.total, expected_total, rtol=1e-10)
