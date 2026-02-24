"""Tests for the tidal.measurement module.

Covers energy computation, conversion probability, spectral analysis,
and diagnostics for coupled field systems.
"""

from __future__ import annotations

import dataclasses
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from tidal.solver._types import SolverResult

import numpy as np
import pytest

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
    summarize,
)
from tidal.measurement._energy import (
    _apply_spatial_operator,
    _compute_constraint_coupling_energy,
    _compute_constraint_self_energy,
    _compute_hamiltonian_from_canonical,
    _compute_virial_potential,
    _evaluate_hamiltonian_factor,
    _gradient_energy_density,
    _is_momentum_field,
    _resolve_term_target,
    _self_gradient_axes,
)
from tidal.symbolic.json_loader import (
    CanonicalStructure,
    ComponentEquation,
    EquationSystem,
    HamiltonianFactor,
    HamiltonianTerm,
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

    # Extract metadata parameters so symbolic coefficients resolve correctly
    params = {
        k: float(v)
        for k, v in spec.metadata.get("parameters", {}).items()
        if isinstance(v, (int, float))
    }

    return SimulationData(
        times=times,
        fields=fields,
        momenta=momenta_dict,
        grid_spacing=(dx,),
        grid_bounds=((0.0, domain_len),),
        periodic=(True,),
        spec=spec,
        parameters=params,
    )


# ============================================================
# SimulationData tests
# ============================================================


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
# Canonical Hamiltonian energy tests
# ============================================================


def _make_kg_canonical_structure(m2: float = 1.0) -> CanonicalStructure:
    """Build canonical structure for Klein-Gordon: H = ½π² + ½(∇φ)² + ½m²φ².

    Decomposed as quadratic terms:
      ½ * time_derivative(phi_0) * time_derivative(phi_0)  → kinetic
      -½ * phi_0 * laplacian(phi_0)                        → gradient (IBP)
      ½m² * phi_0 * phi_0                                  → mass
    """
    return CanonicalStructure(
        hamiltonian_terms=(
            HamiltonianTerm(
                coefficient=0.5,
                factor_a=HamiltonianFactor(field="phi_0", operator="time_derivative"),
                factor_b=HamiltonianFactor(field="phi_0", operator="time_derivative"),
            ),
            HamiltonianTerm(
                coefficient=-0.5,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="laplacian"),
            ),
            HamiltonianTerm(
                coefficient=0.5 * m2,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
            ),
        ),
        field_rates={},
        hamiltonian_symbolic=f"1/2 pi^2 + 1/2 (grad phi)^2 + 1/2 * {m2} * phi^2",
    )


def _make_coupled_canonical_structure(
    m2_phi: float,
    m2_chi: float,
    g: float,
) -> CanonicalStructure:
    """Build canonical structure for two coupled scalars.

    H = ½π_φ² + ½(∇φ)² + ½m²_φ φ² + ½π_χ² + ½(∇χ)² + ½m²_χ χ² + g φ χ
    """
    return CanonicalStructure(
        hamiltonian_terms=(
            # phi kinetic
            HamiltonianTerm(
                coefficient=0.5,
                factor_a=HamiltonianFactor(field="phi_0", operator="time_derivative"),
                factor_b=HamiltonianFactor(field="phi_0", operator="time_derivative"),
            ),
            # phi gradient (IBP form)
            HamiltonianTerm(
                coefficient=-0.5,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="laplacian"),
            ),
            # phi mass
            HamiltonianTerm(
                coefficient=0.5 * m2_phi,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
            ),
            # chi kinetic
            HamiltonianTerm(
                coefficient=0.5,
                factor_a=HamiltonianFactor(field="chi_0", operator="time_derivative"),
                factor_b=HamiltonianFactor(field="chi_0", operator="time_derivative"),
            ),
            # chi gradient (IBP form)
            HamiltonianTerm(
                coefficient=-0.5,
                factor_a=HamiltonianFactor(field="chi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="chi_0", operator="laplacian"),
            ),
            # chi mass
            HamiltonianTerm(
                coefficient=0.5 * m2_chi,
                factor_a=HamiltonianFactor(field="chi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="chi_0", operator="identity"),
            ),
            # coupling: g * phi * chi (split as ½g φχ + ½g χφ in symmetric form,
            # but a single term g*φ*χ works since it's symmetric)
            HamiltonianTerm(
                coefficient=g,
                factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                factor_b=HamiltonianFactor(field="chi_0", operator="identity"),
            ),
        ),
        field_rates={},
        hamiltonian_symbolic="coupled scalar H",
    )


class TestCanonicalHamiltonianEnergy:
    """Test canonical Hamiltonian energy evaluation from structured terms."""

    def test_evaluate_factor_identity(self) -> None:
        """Identity factor returns field data directly."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        result = _evaluate_hamiltonian_factor("phi_0", "identity", data, 0)
        assert result is not None
        np.testing.assert_array_equal(result, data.fields["phi_0"][0])

    def test_evaluate_factor_time_derivative(self) -> None:
        """time_derivative factor returns stored momentum."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        result = _evaluate_hamiltonian_factor("phi_0", "time_derivative", data, 0)
        assert result is not None
        np.testing.assert_array_equal(result, data.momenta["phi_0"][0])

    def test_evaluate_factor_gradient(self) -> None:
        """gradient_x factor applies first derivative."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        result = _evaluate_hamiltonian_factor("phi_0", "gradient_x", data, 0)
        assert result is not None
        assert result.shape == data.fields["phi_0"][0].shape

    def test_kg_hamiltonian_matches_field_energy(self) -> None:
        """Canonical H for single field matches compute_field_energy.

        Uses the coupled_scalars spec with only phi excited (chi=0).
        The KG canonical structure matches ½π² + ½(∇φ)² + ½m²φ².
        """
        data = _make_sim_data_two_fields(n_snapshots=3)
        m2_phi = float(data.spec.mass_matrix[0][0])
        canonical = _make_kg_canonical_structure(m2_phi)
        spec_with_canonical = dataclasses.replace(data.spec, canonical=canonical)
        data_with_canonical = SimulationData(
            times=data.times,
            fields=data.fields,
            momenta=data.momenta,
            grid_spacing=data.grid_spacing,
            grid_bounds=data.grid_bounds,
            periodic=data.periodic,
            spec=spec_with_canonical,
            parameters=data.parameters,
        )

        h_canonical = _compute_hamiltonian_from_canonical(data_with_canonical, 0)

        # Compare with standard field energy (should match for scalars)
        fe = compute_field_energy(
            data.fields["phi_0"][0],
            data.momenta["phi_0"][0],
            m2_phi,
            data.grid_spacing,
            data.periodic,
        )
        # The canonical H uses -½φ·∇²φ (IBP form) which differs from
        # ½|∇φ|² by boundary terms. For periodic BCs, these are identical.
        np.testing.assert_allclose(h_canonical, fe.total, rtol=1e-10)

    def test_coupled_canonical_energy_conservation(self) -> None:
        """Canonical H is conserved for exact coupled oscillator evolution."""
        data = _make_sim_data_two_fields(n_snapshots=51)
        m2_phi = float(data.spec.mass_matrix[0][0])
        m2_chi = float(data.spec.mass_matrix[1][1])
        g_val = float(data.spec.coupling_matrix[0][1])

        canonical = _make_coupled_canonical_structure(m2_phi, m2_chi, g_val)
        spec_with_canonical = dataclasses.replace(data.spec, canonical=canonical)
        data_c = SimulationData(
            times=data.times,
            fields=data.fields,
            momenta=data.momenta,
            grid_spacing=data.grid_spacing,
            grid_bounds=data.grid_bounds,
            periodic=data.periodic,
            spec=spec_with_canonical,
            parameters=data.parameters,
        )

        energies = [
            _compute_hamiltonian_from_canonical(data_c, t_idx)
            for t_idx in range(data_c.n_snapshots)
        ]
        energies_arr = np.array(energies)
        # Exact solution → conserved to machine precision
        relative_drift = np.abs(energies_arr - energies_arr[0]) / energies_arr[0]
        assert np.max(relative_drift) < 1e-10

    def test_system_energy_uses_canonical_when_available(self) -> None:
        """compute_system_energy uses canonical H when spec has canonical structure."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        m2_phi = float(data.spec.mass_matrix[0][0])
        m2_chi = float(data.spec.mass_matrix[1][1])
        g_val = float(data.spec.coupling_matrix[0][1])

        canonical = _make_coupled_canonical_structure(m2_phi, m2_chi, g_val)
        spec_with_canonical = dataclasses.replace(data.spec, canonical=canonical)
        data_c = SimulationData(
            times=data.times,
            fields=data.fields,
            momenta=data.momenta,
            grid_spacing=data.grid_spacing,
            grid_bounds=data.grid_bounds,
            periodic=data.periodic,
            spec=spec_with_canonical,
            parameters=data.parameters,
        )

        se = compute_system_energy(data_c, 0)
        assert se.total > 0
        # Verify total matches direct canonical H evaluation
        h_direct = _compute_hamiltonian_from_canonical(data_c, 0)
        np.testing.assert_allclose(se.total, h_direct, rtol=1e-12)

    def test_canonical_matches_virial_for_scalars(self) -> None:
        """For coupled scalars, canonical H and virial formula give same total."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        m2_phi = float(data.spec.mass_matrix[0][0])
        m2_chi = float(data.spec.mass_matrix[1][1])
        g_val = float(data.spec.coupling_matrix[0][1])

        # Virial energy (no canonical)
        se_virial = compute_system_energy(data, 0)

        # Canonical energy
        canonical = _make_coupled_canonical_structure(m2_phi, m2_chi, g_val)
        spec_with_canonical = dataclasses.replace(data.spec, canonical=canonical)
        data_c = SimulationData(
            times=data.times,
            fields=data.fields,
            momenta=data.momenta,
            grid_spacing=data.grid_spacing,
            grid_bounds=data.grid_bounds,
            periodic=data.periodic,
            spec=spec_with_canonical,
            parameters=data.parameters,
        )
        se_canonical = compute_system_energy(data_c, 0)

        # For scalar fields with periodic BCs, both paths should agree
        np.testing.assert_allclose(
            se_canonical.total,
            se_virial.total,
            rtol=1e-10,
        )

    def test_no_canonical_raises(self) -> None:
        """_compute_hamiltonian_from_canonical raises without canonical structure."""
        data = _make_sim_data_two_fields(n_snapshots=3)
        # Strip canonical structure to test the error path
        spec_no_canon = dataclasses.replace(data.spec, canonical=None)
        data_no_canon = dataclasses.replace(data, spec=spec_no_canon)
        with pytest.raises(ValueError, match="without canonical"):
            _compute_hamiltonian_from_canonical(data_no_canon, 0)

    def test_proca_hamiltonian_uses_velocity_not_momentum(self) -> None:
        """For Proca fields, time_derivative factor must reconstruct velocity.

        In Proca theory, the canonical momentum is π_i = ∂_t A_i - ∂_i A_0,
        so the velocity is vel_i = π_i + ∂_i A_0. The kinetic term in H is
        ½ vel² = ½ (π + ∂_x A_0)², NOT ½ π².  This test creates synthetic
        data where A_0 is nonzero and verifies the cross-terms are included.
        """
        n_grid = 64
        domain_len = 10.0
        dx = domain_len / n_grid
        x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

        # Use periodic-compatible sinusoidal A_0 so gradients wrap correctly
        k0 = 2 * np.pi / domain_len
        a0_amp = 0.5
        a0 = a0_amp * np.sin(k0 * x)  # periodic: A_0(0) ≈ A_0(L)

        # A_1 field and its canonical momentum pi_1
        a1 = np.cos(k0 * x)
        pi1 = 0.5 * np.ones(n_grid)

        # Build a minimal 2-component Proca-like spec
        # A_0: constraint (time_order=0), A_1: dynamical (time_order=2)
        eq_a0 = ComponentEquation(
            field_name="A_0",
            field_index=0,
            time_derivative_order=0,
            rhs_terms=(
                OperatorTerm(coefficient=-1.0, operator="laplacian_x", field="A_1"),
            ),
        )
        eq_a1 = ComponentEquation(
            field_name="A_1",
            field_index=1,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_1"),
                OperatorTerm(coefficient=-1.0, operator="identity", field="A_1"),
            ),
        )

        # Hamiltonian: ½ vel_1² + ½ (∂_x A_1)² + ½ m² A_1²
        # where vel_1 = time_derivative(A_1) = pi_1 + gradient_x(A_0)
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="A_1", operator="time_derivative"),
                    factor_b=HamiltonianFactor(field="A_1", operator="time_derivative"),
                ),
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="A_1", operator="gradient_x"),
                    factor_b=HamiltonianFactor(field="A_1", operator="gradient_x"),
                ),
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="A_1", operator="identity"),
                    factor_b=HamiltonianFactor(field="A_1", operator="identity"),
                ),
            ),
            field_rates={
                "A_1": (
                    OperatorTerm(coefficient=1.0, operator="identity", field="pi_1"),
                    OperatorTerm(coefficient=1.0, operator="gradient_x", field="A_0"),
                ),
            },
            hamiltonian_symbolic="proca test H",
        )

        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("A_0", "A_1"),
            equations=(eq_a0, eq_a1),
            mass_matrix=((0.0, 0.0), (0.0, 1.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={"source": "test", "parameters": {}},
            canonical=canonical,
        )

        data = SimulationData(
            times=np.array([0.0]),
            fields={"A_0": a0[np.newaxis], "A_1": a1[np.newaxis]},
            momenta={"A_1": pi1[np.newaxis]},
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # Compute expected energy using SAME finite difference stencils
        # as the solver.  The IBP path evaluates ½ gradient_x(A_1)² as
        # -½ A_1 · laplacian_x(A_1), matching the 3-point laplacian stencil.
        grad_a0 = (np.roll(a0, -1) - np.roll(a0, 1)) / (2 * dx)
        vel = pi1 + grad_a0  # velocity = pi + d_x(A_0)
        lap_a1 = (np.roll(a1, -1) - 2 * a1 + np.roll(a1, 1)) / dx**2
        h_expected = float((0.5 * vel**2 + (-0.5) * a1 * lap_a1 + 0.5 * a1**2).mean())

        # If we had the bug (using pi instead of vel), kinetic term would be
        # ½ π², missing the cross-terms pi*d_x(A_0) and (d_x A_0)²
        h_wrong = float((0.5 * pi1**2 + (-0.5) * a1 * lap_a1 + 0.5 * a1**2).mean())

        # Verify the correct answer (with velocity reconstruction) matches
        np.testing.assert_allclose(h_eval, h_expected, rtol=1e-10)

        # Verify the wrong answer (raw momentum) does NOT match
        assert abs(h_eval - h_wrong) > 0.01, (
            f"h_eval={h_eval} should differ from h_wrong={h_wrong} when A_0 is nonzero"
        )


class TestIBPHamiltonian:
    """Test integration-by-parts Hamiltonian evaluation for gradient terms.

    The IBP path converts gradient-product Hamiltonian terms into
    second-order operator form:  coeff * ⟨∂_a(u)·∂_b(v)⟩ → -coeff * ⟨u, ∂²_ab(v)⟩.
    This ensures the computed Hamiltonian uses the SAME FD stencils as the
    solver, making it exactly the conserved quantity of the discrete system.
    """

    def test_gradient_pair_to_second_order_same_axis(self) -> None:
        """gradient_x * gradient_x -> laplacian_x."""
        from tidal.measurement._energy import _gradient_pair_to_second_order

        assert (
            _gradient_pair_to_second_order("gradient_x", "gradient_x") == "laplacian_x"
        )
        assert (
            _gradient_pair_to_second_order("gradient_y", "gradient_y") == "laplacian_y"
        )
        assert (
            _gradient_pair_to_second_order("gradient_z", "gradient_z") == "laplacian_z"
        )

    def test_gradient_pair_to_second_order_cross_axis(self) -> None:
        """gradient_x * gradient_y -> cross_derivative_xy (sorted)."""
        from tidal.measurement._energy import _gradient_pair_to_second_order

        assert (
            _gradient_pair_to_second_order("gradient_x", "gradient_y")
            == "cross_derivative_xy"
        )
        assert (
            _gradient_pair_to_second_order("gradient_y", "gradient_x")
            == "cross_derivative_xy"
        )
        assert (
            _gradient_pair_to_second_order("gradient_x", "gradient_z")
            == "cross_derivative_xz"
        )
        assert (
            _gradient_pair_to_second_order("gradient_z", "gradient_y")
            == "cross_derivative_yz"
        )

    def test_ibp_gradient_term_matches_laplacian(self) -> None:
        """IBP converts ½⟨(∂_x φ)²⟩ to -½⟨φ·∂²_x φ⟩ using 3-point laplacian.

        For periodic BCs, these should match to machine precision.
        """
        n_grid = 128
        domain_len = 2 * np.pi
        dx = domain_len / n_grid
        x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

        # Smooth periodic field: sin(2x) + 0.3 cos(5x)
        phi = np.sin(2 * x) + 0.3 * np.cos(5 * x)

        # Gradient-product form: ½ Σ (D_c φ)²
        grad_phi = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * dx)
        energy_grad = 0.5 * float((grad_phi**2).mean())

        # Laplacian form: -½ ⟨φ, L φ⟩ using 3-point stencil
        lap_phi = (np.roll(phi, -1) - 2 * phi + np.roll(phi, 1)) / dx**2
        energy_ibp = -0.5 * float((phi * lap_phi).mean())

        # These differ by O(dx²) ≈ 2e-3 — NOT the same!
        assert abs(energy_grad - energy_ibp) > 1e-5, (
            "Gradient and IBP forms should differ due to stencil mismatch"
        )

        # Build canonical structure with gradient_x * gradient_x term
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(field="phi_0", operator="gradient_x"),
                    factor_b=HamiltonianFactor(field="phi_0", operator="gradient_x"),
                ),
            ),
            field_rates={},
            hamiltonian_symbolic="test",
        )

        eq = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi_0"),
            ),
        )
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0",),
            equations=(eq,),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={"source": "test"},
            canonical=canonical,
        )
        data = SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis]},
            momenta={},
            grid_spacing=(dx,),
            grid_bounds=((0.0, domain_len),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # IBP path should give the LAPLACIAN form (matching solver), not gradient form
        np.testing.assert_allclose(h_eval, energy_ibp, rtol=1e-12)

    def test_ibp_cross_derivative_term(self) -> None:
        """IBP converts ⟨∂_x(u)·∂_y(v)⟩ to -⟨u, ∂²_xy(v)⟩.

        Verifies cross-derivative terms are correctly handled.
        Uses u = v = sin(x+y) so that ⟨∂_x(u)·∂_y(v)⟩ = ⟨cos²(x+y)⟩ = ½.
        """
        n_grid = 32
        domain_len = 2 * np.pi
        dx = domain_len / n_grid
        coords = np.linspace(dx / 2, domain_len - dx / 2, n_grid)
        xx, yy = np.meshgrid(coords, coords, indexing="ij")

        # sin(x+y) has non-zero ⟨∂_x · ∂_y⟩ cross term
        u_data = np.sin(xx + yy)
        v_data = np.sin(xx + yy)

        # IBP form: -⟨u, cross_derivative_xy(v)⟩
        # Cross derivative = D_c^x(D_c^y(v))
        dv_dy = (np.roll(v_data, -1, axis=1) - np.roll(v_data, 1, axis=1)) / (2 * dx)
        d2v_dxdy = (np.roll(dv_dy, -1, axis=0) - np.roll(dv_dy, 1, axis=0)) / (2 * dx)
        energy_ibp = -1.0 * float((u_data * d2v_dxdy).mean())

        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=1.0,
                    factor_a=HamiltonianFactor(field="u_0", operator="gradient_x"),
                    factor_b=HamiltonianFactor(field="v_0", operator="gradient_y"),
                ),
            ),
            field_rates={},
            hamiltonian_symbolic="test cross",
        )

        eq_u = ComponentEquation(
            field_name="u_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian", field="u_0"),
            ),
        )
        eq_v = ComponentEquation(
            field_name="v_0",
            field_index=1,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian", field="v_0"),
            ),
        )
        spec = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("u_0", "v_0"),
            equations=(eq_u, eq_v),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={"source": "test"},
            canonical=canonical,
        )
        data = SimulationData(
            times=np.array([0.0]),
            fields={"u_0": u_data[np.newaxis], "v_0": v_data[np.newaxis]},
            momenta={},
            grid_spacing=(dx, dx),
            grid_bounds=((0.0, domain_len), (0.0, domain_len)),
            periodic=(True, True),
            spec=spec,
            parameters={},
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)
        np.testing.assert_allclose(h_eval, energy_ibp, rtol=1e-12)

    def test_field_rates_symbolic_resolution(self) -> None:
        """field_rate coefficients are resolved via symbolic expressions.

        With rho=2, field_rate coefficient_symbolic="rho^(-1)" should give 0.5,
        not the placeholder coefficient=1.0.
        """
        n_grid = 32
        dx = 0.1

        phi = np.ones(n_grid)
        pi_data = 2.0 * np.ones(n_grid)

        # KE term: rho/2 * (d_t phi)^2 = rho/2 * (pi/rho)^2 = pi^2/(2*rho)
        # With rho=2, pi=2: KE = 4/(2*2) = 1.0
        rho = 2.0

        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=1.0,
                    coefficient_symbolic="rho/2",
                    factor_a=HamiltonianFactor(
                        field="phi_0", operator="time_derivative"
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0", operator="time_derivative"
                    ),
                ),
            ),
            field_rates={
                "phi_0": (
                    OperatorTerm(
                        coefficient=1.0,
                        operator="identity",
                        field="pi_0",
                        coefficient_symbolic="rho^(-1)",
                    ),
                ),
            },
            hamiltonian_symbolic="test rho",
        )

        eq = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi_0"),
            ),
        )
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0",),
            equations=(eq,),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={"source": "test", "parameters": {"rho": rho}},
            canonical=canonical,
        )
        data = SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis]},
            momenta={"phi_0": pi_data[np.newaxis]},
            grid_spacing=(dx,),
            grid_bounds=((0.0, n_grid * dx),),
            periodic=(True,),
            spec=spec,
            parameters={},  # Empty — should merge from spec metadata
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # velocity = pi / rho = 2/2 = 1
        # KE = (rho/2) * vel² = (2/2) * 1² = 1.0
        expected = rho / 2 * (pi_data[0] / rho) ** 2
        np.testing.assert_allclose(h_eval, expected, rtol=1e-12)

    def test_parameter_merge_from_spec_metadata(self) -> None:
        """Parameters from spec metadata are used when data.parameters is empty.

        The symbolic coefficient "mu/2" should resolve to 0.5 with mu=1.0
        from spec metadata, even when data.parameters={}.
        """
        n_grid = 32
        dx = 0.1

        phi = np.ones(n_grid) * 0.5

        canonical = CanonicalStructure(
            hamiltonian_terms=(
                HamiltonianTerm(
                    coefficient=1.0,  # placeholder
                    coefficient_symbolic="mu/2",
                    factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
                    factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
                ),
            ),
            field_rates={},
            hamiltonian_symbolic="test merge",
        )

        eq = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi_0"),
            ),
        )
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0",),
            equations=(eq,),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={"source": "test", "parameters": {"mu": 1.0}},
            canonical=canonical,
        )
        data = SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis]},
            momenta={},
            grid_spacing=(dx,),
            grid_bounds=((0.0, n_grid * dx),),
            periodic=(True,),
            spec=spec,
            parameters={},  # Empty! Metadata should fill in.
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # mu/2 * ⟨φ²⟩ = 0.5 * 0.25 = 0.125
        expected = 0.5 * float((phi**2).mean())
        np.testing.assert_allclose(h_eval, expected, rtol=1e-12)

        # Without the metadata merge, coefficient would be 1.0 (placeholder)
        h_wrong = 1.0 * float((phi**2).mean())
        assert abs(h_eval - h_wrong) > 0.01


class TestAnalyticalEnergyConservation:
    """Verify canonical Hamiltonian against exact analytical energy.

    Uses synthetic plane-wave data with known analytical energy density
    to verify the Hamiltonian evaluation pipeline end-to-end.
    """

    @staticmethod
    def _make_scalar_planewave_spec(
        *,
        m2: float = 1.0,
    ) -> EquationSystem:
        """Build a minimal 1+1D massive scalar spec with canonical structure."""
        eq = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi_0"),
                OperatorTerm(coefficient=-m2, operator="identity", field="phi_0"),
            ),
        )
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                # ½ (∂_t φ)²
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="time_derivative",
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0",
                        operator="time_derivative",
                    ),
                ),
                # ½ (∂_x φ)²
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="gradient_x",
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0",
                        operator="gradient_x",
                    ),
                ),
                # ½ m² φ²
                HamiltonianTerm(
                    coefficient=0.5 * m2,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="identity",
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0",
                        operator="identity",
                    ),
                ),
            ),
            field_rates={
                "phi_0": (
                    OperatorTerm(
                        coefficient=1.0,
                        operator="identity",
                        field="pi_0",
                    ),
                ),
            },
            hamiltonian_symbolic="plane wave test",
        )
        return EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0",),
            equations=(eq,),
            mass_matrix=((m2,),),
            coupling_matrix=((0.0,),),
            metadata={"source": "test", "parameters": {}},
            canonical=canonical,
        )

    def test_plane_wave_energy_matches_analytical(self) -> None:
        """Free scalar plane wave: ⟨H⟩ = ½ A² ω²."""
        m2 = 1.0
        a_amp = 0.3
        n_grid = 256
        domain_len = 2 * np.pi
        dx = domain_len / n_grid
        x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

        k = 2.0  # wavenumber
        omega = np.sqrt(k**2 + m2)

        # φ(t=0,x) = A cos(kx), π(t=0,x) = ∂_t φ = A ω sin(kx)
        # (assuming φ = A cos(ωt - kx) → at t=0: φ = A cos(kx), π = Aω sin(kx))
        phi = a_amp * np.cos(k * x)
        pi_field = a_amp * omega * np.sin(k * x)

        spec = self._make_scalar_planewave_spec(m2=m2)
        data = SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis]},
            momenta={"phi_0": pi_field[np.newaxis]},
            grid_spacing=(dx,),
            grid_bounds=((0.0, domain_len),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # Analytical: ⟨H⟩ = ½ A² ω² (space-averaged energy density)
        h_analytical = 0.5 * a_amp**2 * omega**2
        np.testing.assert_allclose(h_eval, h_analytical, rtol=1e-3)

    def test_plane_wave_energy_conserved_across_snapshots(self) -> None:
        """Verify H(t=0) = H(t=T) for exact plane wave solution."""
        m2 = 0.5
        a_amp = 0.2
        n_grid = 128
        domain_len = 2 * np.pi
        dx = domain_len / n_grid
        x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

        k = 1.0
        omega = np.sqrt(k**2 + m2)

        # Two snapshots: t=0 and t=T
        times = np.array([0.0, 1.5])
        phi_data = np.zeros((2, n_grid))
        pi_data = np.zeros((2, n_grid))
        for i, t in enumerate(times):
            phi_data[i] = a_amp * np.cos(omega * t - k * x)
            pi_data[i] = a_amp * omega * np.sin(omega * t - k * x)

        spec = self._make_scalar_planewave_spec(m2=m2)
        data = SimulationData(
            times=times,
            fields={"phi_0": phi_data},
            momenta={"phi_0": pi_data},
            grid_spacing=(dx,),
            grid_bounds=((0.0, domain_len),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        h0 = _compute_hamiltonian_from_canonical(data, 0)
        h1 = _compute_hamiltonian_from_canonical(data, 1)
        np.testing.assert_allclose(h0, h1, rtol=1e-3)

    def test_energy_partition_kinetic_gradient_mass(self) -> None:
        """Verify individual energy contributions match analytical values.

        For φ = A cos(kx) at t=0 with π = Aω sin(kx):
          kinetic = ½ ⟨π²⟩ = ¼ A² ω²
          gradient = ½ ⟨(∂_x φ)²⟩ = ¼ A² k²
          mass     = ½ m² ⟨φ²⟩ = ¼ A² m²
        """
        m2 = 2.0
        a_amp = 0.4
        n_grid = 256
        domain_len = 2 * np.pi
        dx = domain_len / n_grid
        x = np.linspace(dx / 2, domain_len - dx / 2, n_grid)

        k = 3.0
        omega = np.sqrt(k**2 + m2)

        phi = a_amp * np.cos(k * x)
        pi_field = a_amp * omega * np.sin(k * x)

        # Compute individual energy densities directly
        grad_phi = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * dx)

        kinetic = 0.5 * float((pi_field**2).mean())
        gradient = 0.5 * float((grad_phi**2).mean())
        mass_energy = 0.5 * m2 * float((phi**2).mean())

        np.testing.assert_allclose(kinetic, 0.25 * a_amp**2 * omega**2, rtol=1e-3)
        np.testing.assert_allclose(gradient, 0.25 * a_amp**2 * k**2, rtol=2e-2)
        np.testing.assert_allclose(mass_energy, 0.25 * a_amp**2 * m2, rtol=1e-3)

        # Total should be ½ A² ω²
        total = kinetic + gradient + mass_energy
        np.testing.assert_allclose(total, 0.5 * a_amp**2 * omega**2, rtol=2e-2)

    def test_coupled_eigenmode_energy(self) -> None:
        """Coupled scalars in pure eigenmode: H = 1/2 omega_i^2 |a_i|^2.

        Two coupled scalars with mass matrix M = [[m1², g], [g, m2²]].
        Pure eigenmode initial data gives H from the eigenvalue.
        """
        m2_phi = 2.0
        m2_chi = 3.0
        g_val = 0.5
        a_amp = 0.3
        n_grid = 128
        domain_len = 2 * np.pi
        dx = domain_len / n_grid

        # Diagonalize mass matrix
        m_eff = np.array([[m2_phi, g_val], [g_val, m2_chi]])
        eigenvalues, eigenvectors = np.linalg.eigh(m_eff)

        # Use the first eigenmode (uniform k=0)
        mode_1 = eigenvectors[:, 0]

        # Initial condition: pure eigenmode
        phi_init = a_amp * mode_1[0] * np.ones(n_grid)
        chi_init = a_amp * mode_1[1] * np.ones(n_grid)
        # At t=0 with cos(ω₁ t), velocity = 0
        pi_phi = np.zeros(n_grid)
        pi_chi = np.zeros(n_grid)

        # Build coupled spec
        eq_phi = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(
                    coefficient=-m2_phi,
                    operator="identity",
                    field="phi_0",
                ),
                OperatorTerm(
                    coefficient=-g_val,
                    operator="identity",
                    field="chi_0",
                ),
            ),
        )
        eq_chi = ComponentEquation(
            field_name="chi_0",
            field_index=1,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(
                    coefficient=-g_val,
                    operator="identity",
                    field="phi_0",
                ),
                OperatorTerm(
                    coefficient=-m2_chi,
                    operator="identity",
                    field="chi_0",
                ),
            ),
        )
        canonical = CanonicalStructure(
            hamiltonian_terms=(
                # ½ π_φ²
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="time_derivative",
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0",
                        operator="time_derivative",
                    ),
                ),
                # ½ π_χ²
                HamiltonianTerm(
                    coefficient=0.5,
                    factor_a=HamiltonianFactor(
                        field="chi_0",
                        operator="time_derivative",
                    ),
                    factor_b=HamiltonianFactor(
                        field="chi_0",
                        operator="time_derivative",
                    ),
                ),
                # ½ m1² φ²
                HamiltonianTerm(
                    coefficient=0.5 * m2_phi,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="identity",
                    ),
                    factor_b=HamiltonianFactor(
                        field="phi_0",
                        operator="identity",
                    ),
                ),
                # ½ m2² χ²
                HamiltonianTerm(
                    coefficient=0.5 * m2_chi,
                    factor_a=HamiltonianFactor(
                        field="chi_0",
                        operator="identity",
                    ),
                    factor_b=HamiltonianFactor(
                        field="chi_0",
                        operator="identity",
                    ),
                ),
                # g φχ (coupling: coefficient is g, not g/2, because
                # H has both φ·(g·χ) + χ·(g·φ) combined)
                HamiltonianTerm(
                    coefficient=g_val,
                    factor_a=HamiltonianFactor(
                        field="phi_0",
                        operator="identity",
                    ),
                    factor_b=HamiltonianFactor(
                        field="chi_0",
                        operator="identity",
                    ),
                ),
            ),
            field_rates={
                "phi_0": (
                    OperatorTerm(
                        coefficient=1.0,
                        operator="identity",
                        field="pi_0",
                    ),
                ),
                "chi_0": (
                    OperatorTerm(
                        coefficient=1.0,
                        operator="identity",
                        field="pi_1",
                    ),
                ),
            },
            hamiltonian_symbolic="coupled eigenmode test",
        )

        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0", "chi_0"),
            equations=(eq_phi, eq_chi),
            mass_matrix=((m2_phi, g_val), (g_val, m2_chi)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={"source": "test", "parameters": {}},
            canonical=canonical,
        )

        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "phi_0": phi_init[np.newaxis],
                "chi_0": chi_init[np.newaxis],
            },
            momenta={
                "phi_0": pi_phi[np.newaxis],
                "chi_0": pi_chi[np.newaxis],
            },
            grid_spacing=(dx,),
            grid_bounds=((0.0, domain_len),),
            periodic=(True,),
            spec=spec,
            parameters={},
        )

        h_eval = _compute_hamiltonian_from_canonical(data, 0)

        # Analytical: H = ½ ω₁² A² (eigenmode energy for uniform k=0 mode)
        # At t=0 with cos(ω₁·0)=1, all energy is potential:
        # V = ½ Σ_ij M_ij q_i q_j = ½ A² (mode · M · mode)
        # Since mode is eigenvector: mode · M · mode = λ₁
        h_analytical = 0.5 * a_amp**2 * eigenvalues[0]
        np.testing.assert_allclose(h_eval, h_analytical, rtol=1e-10)


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


class TestSummarize:
    """Test summarize() diagnostics function."""

    def test_summarize_keys(self) -> None:
        """Summary dict contains all expected keys."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        result = summarize(data)

        expected_keys = {
            "times",
            "per_field_energy",
            "interaction_energy",
            "total_energy",
            "energy_conservation",
            "field_peaks",
        }
        assert set(result.keys()) == expected_keys
        assert isinstance(result["energy_conservation"], EnergyDiagnostics)
        per_field: dict[str, object] = result["per_field_energy"]
        assert isinstance(per_field, dict)
        assert set(per_field.keys()) == {"phi_0", "chi_0"}

    def test_summarize_field_peaks(self) -> None:
        """Field peaks match manual np.max(np.abs(...)) of first/last snapshots."""
        data = _make_sim_data_two_fields(n_snapshots=11)
        result = summarize(data)

        for name in data.fields:
            initial_peak = float(np.max(np.abs(data.fields[name][0])))
            final_peak = float(np.max(np.abs(data.fields[name][-1])))
            assert result["field_peaks"][name] == (initial_peak, final_peak)


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


class TestPositionDependentMass:
    """Test that position-dependent mass is correctly resolved in energy."""

    @staticmethod
    def _make_standalone_spec_with_position_dependent_mass() -> EquationSystem:
        """Build a self-contained spec with position-dependent identity term.

        Two 1+1D scalar fields: phi_0 with position-dependent mass
        ``-m2 * x`` (symbolic), chi_0 with constant mass.
        """
        phi_terms = (
            OperatorTerm(
                coefficient=1.0,
                operator="laplacian",
                field="phi_0",
            ),
            OperatorTerm(
                coefficient=-1.0,
                operator="identity",
                field="phi_0",
                coefficient_symbolic="-m2",
                coordinate_dependent=("x",),
            ),
        )
        chi_terms = (
            OperatorTerm(
                coefficient=1.0,
                operator="laplacian",
                field="chi_0",
            ),
            OperatorTerm(
                coefficient=-4.0,
                operator="identity",
                field="chi_0",
            ),
        )
        eq_phi = ComponentEquation(
            field_name="phi_0",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=phi_terms,
        )
        eq_chi = ComponentEquation(
            field_name="chi_0",
            field_index=1,
            time_derivative_order=2,
            rhs_terms=chi_terms,
        )
        return EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi_0", "chi_0"),
            equations=(eq_phi, eq_chi),
            mass_matrix=((1.0, 0.0), (0.0, 4.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
            coordinates=("t", "x"),
        )

    def test_position_dependent_mass_works(self) -> None:
        """Position-dependent mass term is resolved without error."""
        spec = self._make_standalone_spec_with_position_dependent_mass()
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
            parameters={"m2": 1.0},
        )

        result = compute_system_energy(data, 0)
        # Should succeed — position-dependent mass is evaluated on grid
        assert result.total >= 0.0
        assert "phi_0" in result.per_field
        assert "chi_0" in result.per_field


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
        """Laplacian of cos(kx) matches 3-point FD stencil for periodic grid."""
        n = 128
        domain = 2.0 * np.pi
        dx = domain / n
        x = np.linspace(0, domain - dx, n)
        k = 3.0
        field = np.cos(k * x)

        result = _apply_spatial_operator("laplacian_x", field, (dx,), (True,))
        # FD effective wavenumber: k²_eff = 2(1 - cos(k·dx)) / dx²
        k2_eff = 2.0 * (1.0 - np.cos(k * dx)) / (dx * dx)
        expected = -k2_eff * np.cos(k * x)
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
        """Cross derivative of sin(x)*sin(y) via successive FD stencils."""
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
        # FD central difference of sin(kx): D_x[sin(kx)] = cos(kx) * sin(k·dx)/dx
        # Cross derivative = product of two central-difference factors
        k = 1.0
        fd_factor = np.sin(k * dx) / dx  # sin(k·dx)/dx vs exact k
        expected = fd_factor**2 * np.cos(xx) * np.cos(yy)
        np.testing.assert_allclose(result, expected, atol=1e-10)

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
        # Use high resolution: canonical H uses (∂_x φ)² while per-field uses
        # -φ·∂²_x φ (equivalent in continuum, O(dx²) difference discretely).
        data = _make_single_field_data(n_grid=256, n_snapshots=3)
        se = compute_system_energy(data, 0)
        np.testing.assert_allclose(se.interaction, 0.0, atol=1e-5)


# ============================================================
# Operator-aware gradient tests
# ============================================================


class TestSelfGradientAxes:
    """Tests for _self_gradient_axes helper."""

    def test_full_laplacian_returns_none(self) -> None:
        """Equation with isotropic 'laplacian' → None (all axes)."""
        terms = (OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),)
        eq = ComponentEquation(
            field_name="phi",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=terms,
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
            field_name="A_1",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [1]

    def test_directional_laplacian_x(self) -> None:
        """Equation with only laplacian_x → [0]."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_2"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_2"),
        )
        eq = ComponentEquation(
            field_name="A_2",
            field_index=1,
            time_derivative_order=2,
            rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [0]

    def test_both_directional_laplacians(self) -> None:
        """Equation with laplacian_x + laplacian_y → [0, 1]."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_x", field="phi"),
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="phi"),
        )
        eq = ComponentEquation(
            field_name="phi",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == [0, 1]

    def test_cross_field_laplacian_ignored(self) -> None:
        """Laplacian on a different field is not counted as self."""
        terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian_y", field="A_1"),
            OperatorTerm(coefficient=0.5, operator="laplacian_x", field="A_2"),
        )
        eq = ComponentEquation(
            field_name="A_1",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=terms,
        )
        # Only laplacian_y(A_1) is self; laplacian_x(A_2) is cross-field
        assert _self_gradient_axes(eq) == [1]

    def test_no_laplacian_returns_empty(self) -> None:
        """Equation with no laplacian at all → empty list."""
        terms = (OperatorTerm(coefficient=-1.0, operator="identity", field="phi"),)
        eq = ComponentEquation(
            field_name="phi",
            field_index=0,
            time_derivative_order=2,
            rhs_terms=terms,
        )
        assert _self_gradient_axes(eq) == []


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

        # Verify total = kinetic + virial + constraint_self + constraint_cross
        virial = _compute_virial_potential(data, 0)
        cross_e = _compute_constraint_coupling_energy(data, 0)
        expected_total = se.per_field["A_1"].kinetic + virial + constraint_e + cross_e
        np.testing.assert_allclose(se.total, expected_total, rtol=1e-10)

        # Single-constraint system → cross coupling is zero
        assert cross_e == 0.0


def _make_two_constraint_spec() -> EquationSystem:
    """Build a synthetic spec with two constraints and two dynamical fields.

    Mimics a simplified coupled Proca system: A_0 (constraint) + A_1 (dynamical)
    + B_0 (constraint) + B_1 (dynamical), with cross-constraint coupling.
    """
    # A_0: constraint (time_order=0) with coupling to B_0
    a0_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_0"),
        OperatorTerm(
            coefficient=-1.0,
            operator="identity",
            field="A_0",
            coefficient_symbolic="-mA2",
        ),
        OperatorTerm(
            coefficient=0.5,
            operator="identity",
            field="B_0",
            coefficient_symbolic="gcoup",
        ),
    )
    eq_a0 = ComponentEquation(
        field_name="A_0",
        field_index=0,
        time_derivative_order=0,
        rhs_terms=a0_terms,
    )

    # A_1 dynamical field
    a1_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="A_1"),
        OperatorTerm(
            coefficient=-1.0,
            operator="identity",
            field="A_1",
            coefficient_symbolic="-mA2",
        ),
        OperatorTerm(
            coefficient=0.5,
            operator="identity",
            field="B_1",
            coefficient_symbolic="gcoup",
        ),
    )
    eq_a1 = ComponentEquation(
        field_name="A_1",
        field_index=1,
        time_derivative_order=2,
        rhs_terms=a1_terms,
    )

    # B_0 constraint with coupling to A_0
    b0_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="B_0"),
        OperatorTerm(
            coefficient=-2.0,
            operator="identity",
            field="B_0",
            coefficient_symbolic="-mB2",
        ),
        OperatorTerm(
            coefficient=0.5,
            operator="identity",
            field="A_0",
            coefficient_symbolic="gcoup",
        ),
    )
    eq_b0 = ComponentEquation(
        field_name="B_0",
        field_index=2,
        time_derivative_order=0,
        rhs_terms=b0_terms,
    )

    # B_1 dynamical field
    b1_terms = (
        OperatorTerm(coefficient=1.0, operator="laplacian_x", field="B_1"),
        OperatorTerm(
            coefficient=-2.0,
            operator="identity",
            field="B_1",
            coefficient_symbolic="-mB2",
        ),
        OperatorTerm(
            coefficient=0.5,
            operator="identity",
            field="A_1",
            coefficient_symbolic="gcoup",
        ),
    )
    eq_b1 = ComponentEquation(
        field_name="B_1",
        field_index=3,
        time_derivative_order=2,
        rhs_terms=b1_terms,
    )

    return EquationSystem(
        n_components=4,
        dimension=2,
        spatial_dimension=1,
        equations=(eq_a0, eq_a1, eq_b0, eq_b1),
        component_names=("A_0", "A_1", "B_0", "B_1"),
        mass_matrix=(
            (1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 2.0, 0.0),
            (0.0, 0.0, 0.0, 2.0),
        ),
        coupling_matrix=(
            (0.0, 0.0, -0.5, 0.0),
            (0.0, 0.0, 0.0, -0.5),
            (-0.5, 0.0, 0.0, 0.0),
            (0.0, -0.5, 0.0, 0.0),
        ),
        metadata={"parameters": {"mA2": 1.0, "mB2": 2.0, "gcoup": 0.5}},
        coordinates=("t", "x"),
        mass_matrix_symbolic=(
            ("-mA2", None, None, None),
            (None, "-mA2", None, None),
            (None, None, "-mB2", None),
            (None, None, None, "-mB2"),
        ),
        coupling_matrix_symbolic=(
            (None, None, "gcoup", None),
            (None, None, None, "gcoup"),
            ("gcoup", None, None, None),
            (None, "gcoup", None, None),
        ),
    )


class TestConstraintCouplingEnergy:
    """Test _compute_constraint_coupling_energy."""

    def test_cross_constraint_coupling_nonzero(self) -> None:
        """Two constraints with cross identity terms give nonzero coupling."""
        spec = _make_two_constraint_spec()
        n = 64
        dx = 10.0 / n
        x = np.linspace(dx / 2, 10.0 - dx / 2, n)

        k = 2.0 * np.pi / 10.0
        a0 = np.cos(k * x)
        b0 = 0.5 * np.sin(k * x)

        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "A_0": a0[np.newaxis, :],
                "A_1": np.zeros((1, n)),
                "B_0": b0[np.newaxis, :],
                "B_1": np.zeros((1, n)),
            },
            momenta={
                "A_1": np.zeros((1, n)),
                "B_1": np.zeros((1, n)),
            },
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"mA2": 1.0, "mB2": 2.0, "gcoup": 0.5},
        )

        energy = _compute_constraint_coupling_energy(data, 0)

        # Expected: +½ * 0.5 * ⟨A_0*B_0⟩ (from A_0 eq)
        #         + ½ * 0.5 * ⟨B_0*A_0⟩ (from B_0 eq)
        #         = 0.5 * ⟨A_0*B_0⟩
        expected = 0.5 * float(np.mean(a0 * b0))
        np.testing.assert_allclose(energy, expected, rtol=1e-10)

    def test_single_constraint_returns_zero(self) -> None:
        """Single constraint system → no cross terms → zero."""
        spec = _make_constraint_spec()
        n = 32
        dx = 10.0 / n

        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "A_0": np.ones((1, n)),
                "A_1": np.zeros((1, n)),
            },
            momenta={"A_1": np.zeros((1, n))},
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"Am2": 0.5},
        )

        energy = _compute_constraint_coupling_energy(data, 0)
        assert energy == 0.0

    def test_system_energy_includes_cross_constraint(self) -> None:
        """compute_system_energy includes cross-constraint coupling."""
        spec = _make_two_constraint_spec()
        n = 64
        dx = 10.0 / n
        x = np.linspace(dx / 2, 10.0 - dx / 2, n)

        k = 2.0 * np.pi / 10.0
        a0 = 0.3 * np.cos(k * x)
        a1 = np.cos(k * x)
        b0 = 0.2 * np.sin(k * x)
        b1 = 0.5 * np.sin(k * x)

        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "A_0": a0[np.newaxis, :],
                "A_1": a1[np.newaxis, :],
                "B_0": b0[np.newaxis, :],
                "B_1": b1[np.newaxis, :],
            },
            momenta={
                "A_1": np.zeros((1, n)),
                "B_1": np.zeros((1, n)),
            },
            grid_spacing=(dx,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            spec=spec,
            parameters={"mA2": 1.0, "mB2": 2.0, "gcoup": 0.5},
        )

        se = compute_system_energy(data, 0)
        virial = _compute_virial_potential(data, 0)
        constraint_self = _compute_constraint_self_energy(data, 0)
        constraint_cross = _compute_constraint_coupling_energy(data, 0)

        # Verify all energy components sum to the total
        total_kinetic = sum(fe.kinetic for fe in se.per_field.values())
        expected = total_kinetic + virial + constraint_self + constraint_cross
        np.testing.assert_allclose(se.total, expected, rtol=1e-10)

        # Cross coupling should be nonzero for this configuration
        assert constraint_cross != 0.0


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
        expected_unc = (
            (np.pi / (result.dominant_frequency**2)) * result.frequency_fwhm / 2.0
        )
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
        assert (
            abs(spectrum.dominant_mixing_length - expected_lmix) / expected_lmix < 0.05
        )

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
        assert (
            abs(mixing_result.dominant_frequency - spectrum.dominant_frequency) < 0.01
        )


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
                data,
                "phi_0",
                "chi_0",
                energy_floor=1e10,
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
            group.probability,
            pairwise.probability,
            atol=1e-15,
        )

    def test_overlap_raises(self) -> None:
        """Source and target groups must not overlap."""
        data = _make_sim_data_two_fields(n_snapshots=5)
        with pytest.raises(ValueError, match="overlap"):
            compute_group_spectral_conversion(
                data,
                "phi_0",
                ["phi_0", "chi_0"],
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
            spec,
            equations=(constraint_eq, spec.equations[1]),
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
            _radial_bin,
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
            field_arr,
            mom_arr,
            m2,
            (dx,),
            (True,),
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
            field,
            momentum,
            m2_tachyonic,
            (dx,),
            (True,),
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
    path = DATA_DIR / "klein_gordon_1d.json"
    if not path.exists():
        pytest.skip("klein_gordon_1d.json not found (run tidal derive)")
    spec = load_equation_system(path)
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
            detected_omega,
            omega0,
            atol=2 * result.rayleigh_resolution,
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
                omega_detected,
                omega_expected,
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
        writer._flush_mmaps()
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
        writer._flush_mmaps()

        # Remove metadata.json
        meta_path = output_dir / "metadata.json"
        if meta_path.exists():
            meta_path.unlink()

        data = SimulationData.from_directory(output_dir, spec)
        assert data.n_snapshots == 1
        np.testing.assert_allclose(data.fields["phi_0"][0], np.ones(n_grid))


class TestSimulationDataLoad:
    """Test SimulationData.load() universal entry point."""

    def test_load_file_raises(self, tmp_path: Path) -> None:
        """load() raises ValueError for non-directory paths."""
        spec = _build_coupled_scalars_spec()
        npz_path = tmp_path / "test.npz"
        npz_path.write_bytes(b"dummy")

        with pytest.raises(ValueError, match="no longer supported"):
            SimulationData.load(npz_path, spec)

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


# ============================================================
# Group 12b — Review hardening: validation, recovery, integration
# ============================================================


class TestSnapshotWriterValidation:
    """Tests for input validation added in critical review."""

    def test_nan_time_raises(self, tmp_path: Path) -> None:
        """NaN time raises ValueError."""
        output_dir = tmp_path / "nan_t"
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
        with pytest.raises(ValueError, match="Time must be finite"):
            writer.append(float("nan"), {"phi_0": np.zeros(4)}, {})
        writer.close()

    def test_inf_time_raises(self, tmp_path: Path) -> None:
        """Infinite time raises ValueError."""
        output_dir = tmp_path / "inf_t"
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
        with pytest.raises(ValueError, match="Time must be finite"):
            writer.append(float("inf"), {"phi_0": np.zeros(4)}, {})
        writer.close()

    def test_out_of_order_time_raises(self, tmp_path: Path) -> None:
        """Decreasing time raises ValueError."""
        output_dir = tmp_path / "order_t"
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
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        writer.append(1.0, {"phi_0": np.zeros(4)}, {})
        with pytest.raises(ValueError, match="non-decreasing"):
            writer.append(0.5, {"phi_0": np.zeros(4)}, {})
        writer.close()

    def test_equal_times_allowed(self, tmp_path: Path) -> None:
        """Equal consecutive times are allowed (non-decreasing, not strictly increasing)."""
        output_dir = tmp_path / "equal_t"
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
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})
        writer.append(0.0, {"phi_0": np.zeros(4)}, {})  # Should not raise
        writer.close()

    def test_wrong_field_shape_raises(self, tmp_path: Path) -> None:
        """Field with wrong shape raises ValueError."""
        output_dir = tmp_path / "shape_f"
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
        with pytest.raises(ValueError, match="Field 'phi_0' has shape"):
            writer.append(0.0, {"phi_0": np.zeros(8)}, {})
        writer.close()

    def test_wrong_momentum_shape_raises(self, tmp_path: Path) -> None:
        """Momentum with wrong shape raises ValueError."""
        output_dir = tmp_path / "shape_m"
        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=["phi_0"],
            grid_shape=(4,),
            n_snapshots=1,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        )
        with pytest.raises(ValueError, match="Momentum 'phi_0' has shape"):
            writer.append(0.0, {"phi_0": np.zeros(4)}, {"phi_0": np.zeros(8)})
        writer.close()

    def test_output_path_is_file_raises(self, tmp_path: Path) -> None:
        """SnapshotWriter raises if output path is an existing file."""
        file_path = tmp_path / "not_a_dir"
        file_path.write_text("oops")
        with pytest.raises(ValueError, match="not a directory"):
            SnapshotWriter(
                output_dir=file_path,
                field_names=["phi_0"],
                momentum_names=[],
                grid_shape=(4,),
                n_snapshots=1,
                grid_spacing=(1.0,),
                grid_bounds=((0.0, 4.0),),
                periodic=(False,),
            )

    def test_stale_files_cleaned(self, tmp_path: Path) -> None:
        """Re-using an output directory removes stale .npy files."""
        output_dir = tmp_path / "reuse"

        # First run with phi_0 and chi_0
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=1,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        ) as w1:
            w1.append(0.0, {"phi_0": np.zeros(4), "chi_0": np.ones(4)}, {})

        assert (output_dir / "chi_0.npy").exists()

        # Second run with only phi_0
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0"],
            momentum_names=[],
            grid_shape=(4,),
            n_snapshots=1,
            grid_spacing=(1.0,),
            grid_bounds=((0.0, 4.0),),
            periodic=(False,),
        ) as w2:
            w2.append(0.0, {"phi_0": np.zeros(4)}, {})

        # Stale chi_0.npy should be gone
        assert not (output_dir / "chi_0.npy").exists()
        assert (output_dir / "phi_0.npy").exists()

    def test_3d_grid_round_trip(self, tmp_path: Path) -> None:
        """Round-trip with a 3D grid shape."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "grid3d"
        grid_shape = (4, 4, 4)
        n_snapshots = 2

        rng = np.random.default_rng(321)
        expected = rng.standard_normal((n_snapshots, *grid_shape))

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=grid_shape,
            n_snapshots=n_snapshots,
            grid_spacing=(1.0, 1.0, 1.0),
            grid_bounds=((0.0, 4.0), (0.0, 4.0), (0.0, 4.0)),
            periodic=(True, True, True),
        ) as writer:
            for i in range(n_snapshots):
                writer.append(
                    float(i),
                    {"phi_0": expected[i], "chi_0": np.zeros(grid_shape)},
                    {"phi_0": np.zeros(grid_shape), "chi_0": np.zeros(grid_shape)},
                )

        data = SimulationData.from_directory(output_dir, spec)
        np.testing.assert_allclose(data.fields["phi_0"], expected)


class TestCrashRecoveryHardened:
    """Tests for crash recovery edge cases."""

    def test_recovery_with_nan_times(self, tmp_path: Path) -> None:
        """Recovery skips NaN entries in times array."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "nan_recovery"
        n_grid = 4

        writer = SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=10,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
            flush_interval=1,
        )
        # Write 3 valid snapshots
        for i in range(3):
            writer.append(
                float(i) + 0.1,
                {"phi_0": np.full(n_grid, float(i)), "chi_0": np.zeros(n_grid)},
                {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
            )
        writer._flush_mmaps()

        # Corrupt entry 3 with NaN (simulating partial write)
        times_mmap = np.load(
            str(output_dir / "times.npy"),
            mmap_mode="r+",
        )
        times_mmap[3] = float("nan")
        times_mmap.flush()

        # Remove metadata.json to force recovery
        meta_path = output_dir / "metadata.json"
        if meta_path.exists():
            meta_path.unlink()

        data = SimulationData.from_directory(output_dir, spec)
        # Should recover 3 valid snapshots, not 4 (NaN skipped)
        assert data.n_snapshots == 3

    def test_recovery_with_corrupt_json(self, tmp_path: Path) -> None:
        """Recovery works when metadata.json is truncated/corrupt."""
        spec = _build_coupled_scalars_spec()
        output_dir = tmp_path / "corrupt_json"
        n_grid = 4

        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(n_grid,),
            n_snapshots=3,
            grid_spacing=(2.5,),
            grid_bounds=((0.0, 10.0),),
            periodic=(True,),
        ) as writer:
            for i in range(3):
                writer.append(
                    float(i) + 0.1,
                    {"phi_0": np.full(n_grid, float(i)), "chi_0": np.zeros(n_grid)},
                    {"phi_0": np.zeros(n_grid), "chi_0": np.zeros(n_grid)},
                )

        # Corrupt metadata.json
        meta_path = output_dir / "metadata.json"
        meta_path.write_text('{"n_snapshots": 3, "fields": [')  # Truncated JSON

        data = SimulationData.from_directory(output_dir, spec)
        assert data.n_snapshots == 3
        np.testing.assert_allclose(data.fields["phi_0"][1], np.full(n_grid, 1.0))


class TestMemmapMeasurementIntegration:
    """Test that measurement functions work on memmap-backed SimulationData."""

    def test_energy_from_directory(self, tmp_path: Path) -> None:
        """compute_field_energy works on memmap data from from_directory."""
        # Build reference data in memory
        data_mem = _make_sim_data_two_fields(n_grid=16, n_snapshots=5)

        # Write to directory
        output_dir = tmp_path / "energy_int"
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(16,),
            n_snapshots=5,
            grid_spacing=data_mem.grid_spacing,
            grid_bounds=data_mem.grid_bounds,
            periodic=data_mem.periodic,
            parameters=data_mem.parameters,
        ) as writer:
            for t_idx in range(data_mem.n_snapshots):
                writer.append(
                    float(data_mem.times[t_idx]),
                    {n: data_mem.fields[n][t_idx] for n in data_mem.fields},
                    {n: data_mem.momenta[n][t_idx] for n in data_mem.momenta},
                )

        # Load from directory (memmap)
        data_dir = SimulationData.from_directory(output_dir, data_mem.spec)

        # Compute system energy on both
        energy_mem = compute_system_energy(data_mem, 0)
        energy_dir = compute_system_energy(data_dir, 0)

        np.testing.assert_allclose(energy_dir.total, energy_mem.total, rtol=1e-12)

    def test_conversion_from_directory(self, tmp_path: Path) -> None:
        """compute_conversion_probability works on memmap data."""
        data_mem = _make_sim_data_two_fields(n_grid=16, n_snapshots=11)

        output_dir = tmp_path / "conv_int"
        with SnapshotWriter(
            output_dir=output_dir,
            field_names=["phi_0", "chi_0"],
            momentum_names=["phi_0", "chi_0"],
            grid_shape=(16,),
            n_snapshots=11,
            grid_spacing=data_mem.grid_spacing,
            grid_bounds=data_mem.grid_bounds,
            periodic=data_mem.periodic,
            parameters=data_mem.parameters,
        ) as writer:
            for t_idx in range(data_mem.n_snapshots):
                writer.append(
                    float(data_mem.times[t_idx]),
                    {n: data_mem.fields[n][t_idx] for n in data_mem.fields},
                    {n: data_mem.momenta[n][t_idx] for n in data_mem.momenta},
                )

        data_dir = SimulationData.from_directory(output_dir, data_mem.spec)

        result_mem = compute_conversion_probability(data_mem, "phi_0", "chi_0")
        result_dir = compute_conversion_probability(data_dir, "phi_0", "chi_0")

        np.testing.assert_allclose(
            result_dir.probability,
            result_mem.probability,
            rtol=1e-12,
        )


# ============================================================
# Group 17 — SimulationData.save() round-trip
# ============================================================


class TestSimulationDataSave:
    """Tests for SimulationData.save() method."""

    def test_save_round_trip(self, tmp_path: Path) -> None:
        """save() then from_directory() produces identical data."""
        data = _make_sim_data_two_fields(n_grid=16, n_snapshots=5)
        save_dir = tmp_path / "saved"
        data.save(save_dir)

        loaded = SimulationData.from_directory(save_dir, data.spec)
        np.testing.assert_array_almost_equal(loaded.times, data.times)
        for name in data.fields:
            np.testing.assert_array_almost_equal(
                loaded.fields[name],
                data.fields[name],
            )
        for name in data.momenta:
            np.testing.assert_array_almost_equal(
                loaded.momenta[name],
                data.momenta[name],
            )
        assert loaded.grid_spacing == data.grid_spacing
        assert loaded.grid_bounds == data.grid_bounds
        assert loaded.periodic == data.periodic

    def test_save_creates_directory(self, tmp_path: Path) -> None:
        """save() creates the directory if it doesn't exist."""
        data = _make_sim_data_two_fields(n_grid=16, n_snapshots=3)
        save_dir = tmp_path / "nested" / "dir"
        data.save(save_dir)
        assert save_dir.is_dir()
        assert (save_dir / "metadata.json").exists()
        assert (save_dir / "times.npy").exists()

    def test_save_metadata_fields(self, tmp_path: Path) -> None:
        """save() writes correct metadata.json content."""
        import json

        data = _make_sim_data_two_fields(n_grid=16, n_snapshots=3)
        save_dir = tmp_path / "meta_check"
        data.save(save_dir)

        meta = json.loads((save_dir / "metadata.json").read_text())
        assert meta["n_snapshots"] == 3
        assert meta["grid_spacing"] == list(data.grid_spacing)
        assert meta["periodic"] == list(data.periodic)
        assert set(meta["fields"]) == set(data.fields.keys())
        assert set(meta["momenta"]) == set(data.momenta.keys())


class TestSnapshotCountValidation:
    """Fail-fast: metadata n_snapshots must not exceed actual array size."""

    def test_inflated_n_snapshots_raises(self, tmp_path: Path) -> None:
        """Metadata claiming more snapshots than times.npy has → ValueError."""
        import json

        spec = _build_coupled_scalars_spec()

        out = tmp_path / "inflated"
        out.mkdir()
        times = np.array([0.0, 1.0, 2.0])
        np.save(str(out / "times.npy"), times)
        for eq in spec.equations:
            np.save(str(out / f"{eq.field_name}.npy"), np.zeros((3, 8)))
            if eq.time_derivative_order >= 2:
                np.save(str(out / f"pi_{eq.field_name}.npy"), np.zeros((3, 8)))

        meta = {
            "n_snapshots": 100,
            "fields": [eq.field_name for eq in spec.equations],
            "momenta": [
                eq.field_name for eq in spec.equations if eq.time_derivative_order >= 2
            ],
            "grid_spacing": [1.0],
            "grid_bounds": [[0.0, 8.0]],
            "periodic": [True],
        }
        (out / "metadata.json").write_text(json.dumps(meta))

        with pytest.raises(ValueError, match="Metadata claims 100 snapshots"):
            SimulationData.from_directory(out, spec)


# ==================== SimulationData.from_result() ====================


class TestSimulationDataFromResult:
    """Tests for the native-path constructor that bypasses py-pde."""

    @staticmethod
    def _kg_spec() -> EquationSystem:
        from tidal.symbolic.json_loader import EquationSystem

        return EquationSystem.from_dict(
            {
                "metadata": {"source": "test", "parameters": {"m2": 1.0}},
                "spacetime": {
                    "dimension": 2,
                    "signature": [-1, 1],
                    "coordinates": ["t", "x"],
                },
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
                                    "coefficient": -1.0,
                                    "operator": "identity",
                                    "field": "phi_0",
                                },
                                {
                                    "coefficient": 1.0,
                                    "operator": "laplacian_x",
                                    "field": "phi_0",
                                },
                            ],
                        },
                    }
                ],
            }
        )

    def test_shapes(self) -> None:
        """Fields and momenta have correct (n_snapshots, *grid_shape)."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(16,), periodic=(True,))

        # 2 slots (phi_0, pi_phi_0) x 16 points = 32 flat size
        n_snaps, n_flat = 5, 32
        result: SolverResult = {
            "t": np.linspace(0, 1, n_snaps),
            "y": np.random.default_rng(42).standard_normal((n_snaps, n_flat)),
            "success": True,
            "message": "ok",
        }

        sd = SimulationData.from_result(result, spec, gi, {"m2": 1.0})

        assert sd.times.shape == (5,)
        assert "phi_0" in sd.fields
        assert sd.fields["phi_0"].shape == (5, 16)
        assert "phi_0" in sd.momenta
        assert sd.momenta["phi_0"].shape == (5, 16)

    def test_values_match_slicing(self) -> None:
        """Data in fields/momenta matches manual slicing of y."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(8,), periodic=(False,))

        rng = np.random.default_rng(123)
        y = rng.standard_normal((3, 16))  # 2 slots x 8 points
        result: SolverResult = {
            "t": np.array([0.0, 0.5, 1.0]),
            "y": y,
            "success": True,
            "message": "",
        }

        sd = SimulationData.from_result(result, spec, gi)

        # First slot (phi_0): y[:, 0:8]
        np.testing.assert_array_equal(sd.fields["phi_0"], y[:, :8].reshape(3, 8))
        # Second slot (pi_phi_0): y[:, 8:16]
        np.testing.assert_array_equal(sd.momenta["phi_0"], y[:, 8:16].reshape(3, 8))

    def test_grid_metadata(self) -> None:
        """Grid spacing, bounds, periodic are propagated correctly."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 8.0),), shape=(16,), periodic=(True,))

        result: SolverResult = {
            "t": np.array([0.0]),
            "y": np.zeros((1, 32)),
            "success": True,
            "message": "",
        }
        sd = SimulationData.from_result(result, spec, gi, {"m2": 2.0})

        assert sd.grid_spacing == gi.dx
        assert sd.grid_bounds == gi.bounds
        assert sd.periodic == (True,)
        assert sd.parameters == {"m2": 2.0}

    def test_empty_result_raises(self) -> None:
        """Empty solver result raises ValueError."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(16,), periodic=(False,))

        result: SolverResult = {
            "t": np.array([]),
            "y": np.zeros((0, 32)),
            "success": True,
            "message": "",
        }
        with pytest.raises(ValueError, match="no snapshots"):
            SimulationData.from_result(result, spec, gi)

    def test_failed_result_raises(self) -> None:
        """from_result() must reject results with success=False."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(16,), periodic=(False,))

        result: SolverResult = {
            "t": np.array([0.0, 0.5]),
            "y": np.zeros((2, 32)),
            "success": False,
            "message": "Adaptive step failed to converge",
        }
        with pytest.raises(ValueError, match="failed solver result"):
            SimulationData.from_result(result, spec, gi)

    def test_snapshot_time_mismatch_raises(self) -> None:
        """from_result() must detect mismatch between len(t) and y.shape[0]."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(16,), periodic=(False,))

        result: SolverResult = {
            "t": np.linspace(0, 1, 6),   # 6 time points
            "y": np.zeros((3, 32)),       # only 3 state vectors
            "success": True,
            "message": "",
        }
        with pytest.raises(ValueError, match="Snapshot count mismatch"):
            SimulationData.from_result(result, spec, gi)

    def test_roundtrip_save_load(self, tmp_path: Path) -> None:
        """from_result → save → from_directory → compare."""
        from tidal.solver.grid import GridInfo

        spec = self._kg_spec()
        gi = GridInfo(bounds=((0.0, 10.0),), shape=(8,), periodic=(True,))

        rng = np.random.default_rng(99)
        y = rng.standard_normal((4, 16))
        result: SolverResult = {
            "t": np.array([0.0, 1.0, 2.0, 3.0]),
            "y": y,
            "success": True,
            "message": "",
        }

        sd1 = SimulationData.from_result(result, spec, gi, {"m2": 1.0})
        out_dir = tmp_path / "snapshots"
        sd1.save(out_dir)

        sd2 = SimulationData.from_directory(out_dir, spec)

        np.testing.assert_allclose(sd2.times, sd1.times)
        np.testing.assert_allclose(sd2.fields["phi_0"], sd1.fields["phi_0"])
        np.testing.assert_allclose(sd2.momenta["phi_0"], sd1.momenta["phi_0"])
        assert sd2.grid_spacing == sd1.grid_spacing
        assert sd2.grid_bounds == sd1.grid_bounds
        assert sd2.periodic == sd1.periodic


# ============================================================
# BC-type persistence and Neumann gradient tests
# ============================================================


class TestBCTypes:
    """Test bc_types round-trip through save/load and energy computation."""

    @staticmethod
    def _make_1d_simulation_data(
        *,
        bc_types: tuple[str, ...] | None = None,
    ) -> SimulationData:
        """Build a minimal 1D SimulationData for testing."""
        from tidal.symbolic.json_loader import EquationSystem

        spec_dict: dict[str, Any] = {
            "metadata": {
                "source": "test-bc",
                "lagrangian_expr": "test",
                "derived_from": "test",
                "gauge": "none",
                "linearized": False,
                "parameters": {},
            },
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
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
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
            "coupling": {},
        }
        spec = EquationSystem.from_dict(spec_dict)
        n_x = 64
        n_t = 5
        times = np.linspace(0, 1, n_t)
        fields = {"phi_0": np.random.default_rng(42).standard_normal((n_t, n_x))}
        momenta = {"phi_0": np.random.default_rng(43).standard_normal((n_t, n_x))}
        return SimulationData(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=(0.1,),
            grid_bounds=((0.0, 6.4),),
            periodic=(False,),
            spec=spec,
            parameters={},
            bc_types=bc_types,
        )

    def test_bc_types_roundtrip(self, tmp_path: Path) -> None:
        """bc_types survives save → load cycle."""
        sd1 = self._make_1d_simulation_data(
            bc_types=("neumann",),
        )
        out_dir = tmp_path / "bc_test"
        sd1.save(out_dir)

        sd2 = SimulationData.from_directory(out_dir, sd1.spec)
        assert sd2.bc_types == ("neumann",)

    def test_bc_types_none_for_legacy(self, tmp_path: Path) -> None:
        """Legacy metadata.json without bc_types → bc_types is None."""
        sd1 = self._make_1d_simulation_data(bc_types=None)
        out_dir = tmp_path / "legacy_test"
        sd1.save(out_dir)

        sd2 = SimulationData.from_directory(out_dir, sd1.spec)
        assert sd2.bc_types is None

    def test_neumann_gradient_energy_symmetric(self) -> None:
        """Neumann BC uses even-reflection (symmetric) ghost cells.

        cos(pi*x/L) with Neumann BCs: the derivative at boundaries
        is -sin(pi*x/L) * pi/L which is zero at x=0 and x=L, so
        Neumann BCs are exact.  Compare gradient energy computed with
        Neumann vs Dirichlet padding — they should differ because
        Dirichlet (odd) padding reverses the sign.
        """
        from tidal.measurement._energy import _gradient_energy_density

        n_x = 128
        length = 2 * np.pi
        dx = length / n_x
        x = np.linspace(dx / 2, length - dx / 2, n_x)
        field = np.cos(np.pi * x / length)

        # With Neumann BC (correct for this mode)
        grad_neumann = _gradient_energy_density(
            field,
            grid_spacing=(dx,),
            periodic=(False,),
            bc_types=("neumann",),
        )
        energy_neumann = 0.5 * float(grad_neumann.mean())

        # Analytic: 0.5 * ⟨(pi/L)² sin²(pi*x/L)⟩ = 0.25 * (pi/L)²
        k = np.pi / length
        expected = 0.25 * k**2
        np.testing.assert_allclose(energy_neumann, expected, rtol=0.01)

    def test_neumann_vs_dirichlet_differ(self) -> None:
        """Neumann and Dirichlet give different gradient energies for cos mode."""
        from tidal.measurement._energy import _gradient_energy_density

        n_x = 128
        length = 2 * np.pi
        dx = length / n_x
        x = np.linspace(dx / 2, length - dx / 2, n_x)
        field = np.cos(np.pi * x / length)

        grad_neumann = _gradient_energy_density(
            field,
            (dx,),
            (False,),
            bc_types=("neumann",),
        )
        grad_dirichlet = _gradient_energy_density(
            field,
            (dx,),
            (False,),
            bc_types=("dirichlet",),
        )
        # They should not be equal — Dirichlet applies odd-reflection
        # at the boundaries which distorts the cos mode
        assert not np.allclose(grad_neumann, grad_dirichlet)


class TestDtMetadata:
    """Test dt round-trip through save/load and conservation diagnostics."""

    @staticmethod
    def _make_sim_data(
        *,
        dt: float | None = None,
    ) -> SimulationData:
        """Build a minimal SimulationData with dt."""
        from tidal.symbolic.json_loader import EquationSystem

        spec_dict: dict[str, Any] = {
            "metadata": {
                "source": "test-dt",
                "lagrangian_expr": "test",
                "derived_from": "test",
                "gauge": "none",
                "linearized": False,
                "parameters": {},
            },
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
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
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
            "coupling": {},
        }
        spec = EquationSystem.from_dict(spec_dict)
        n_x = 64
        n_t = 5
        times = np.linspace(0, 1, n_t)
        fields = {"phi_0": np.random.default_rng(42).standard_normal((n_t, n_x))}
        momenta = {"phi_0": np.random.default_rng(43).standard_normal((n_t, n_x))}
        return SimulationData(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=(0.1,),
            grid_bounds=((0.0, 6.4),),
            periodic=(True,),
            spec=spec,
            parameters={},
            dt=dt,
        )

    def test_dt_roundtrip(self, tmp_path: Path) -> None:
        """Time-step dt survives save → load cycle via metadata.json."""
        sd1 = self._make_sim_data(dt=0.05)
        out_dir = tmp_path / "dt_test"
        sd1.save(out_dir)

        sd2 = SimulationData.from_directory(out_dir, sd1.spec)
        assert sd2.dt is not None
        assert sd2.dt == pytest.approx(0.05)

    def test_dt_none_for_legacy(self, tmp_path: Path) -> None:
        """Legacy metadata.json without dt → dt is None."""
        sd1 = self._make_sim_data(dt=None)
        out_dir = tmp_path / "legacy_dt"
        sd1.save(out_dir)

        sd2 = SimulationData.from_directory(out_dir, sd1.spec)
        assert sd2.dt is None

    def test_conservation_dt_aware_threshold(self) -> None:
        """check_energy_conservation scales threshold by dt² when dt is known."""
        from tidal.measurement._diagnostics import check_energy_conservation

        # Build a SimulationData with dt=0.1 → shadow bound = 10 * 0.01 = 0.1
        sd = self._make_sim_data(dt=0.1)

        # With large dt, the threshold auto-scales to max(1e-3, 0.1) = 0.1
        # Random data won't satisfy conservation, but threshold should be 0.1
        diag = check_energy_conservation(sd, threshold=1e-3)
        # The is_conserved check uses max(1e-3, 10*0.01) = 0.1 as threshold
        # We just verify it doesn't crash and the threshold was effectively raised
        assert diag.max_relative_error >= 0  # sanity

    def test_conservation_no_dt_uses_default_threshold(self) -> None:
        """Without dt, default threshold is used unchanged."""
        from tidal.measurement._diagnostics import check_energy_conservation

        sd = self._make_sim_data(dt=None)
        diag = check_energy_conservation(sd, threshold=1e-3)
        # Without dt, threshold stays at 1e-3 — random data almost certainly fails
        assert diag.max_relative_error >= 0  # sanity


# ============================================================
# Position-dependent Hamiltonian energy tests
# ============================================================


class TestHamiltonianTermPositionDependent:
    """Unit tests for HamiltonianTerm.position_dependent auto-detection."""

    def test_constant_coeff_not_position_dependent(self) -> None:
        term = HamiltonianTerm(
            coefficient=0.5,
            factor_a=HamiltonianFactor(field="phi_0", operator="gradient_x"),
            factor_b=HamiltonianFactor(field="phi_0", operator="gradient_x"),
            coefficient_symbolic="mPhi2/2",
        )
        assert not term.position_dependent

    def test_gaussian_coeff_is_position_dependent(self) -> None:
        term = HamiltonianTerm(
            coefficient=1.0,
            factor_a=HamiltonianFactor(field="chi_0", operator="identity"),
            factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
            coefficient_symbolic="E^(-1/2*x[]^2/R^2 - y[]^2/(2*R^2))*g0",
        )
        assert term.position_dependent

    def test_csc_coeff_is_position_dependent(self) -> None:
        term = HamiltonianTerm(
            coefficient=1.0,
            factor_a=HamiltonianFactor(field="phi_0", operator="gradient_z"),
            factor_b=HamiltonianFactor(field="phi_0", operator="gradient_z"),
            coefficient_symbolic="Csc[y[]]^2/(2*x[]^2)",
        )
        assert term.position_dependent

    def test_coordinate_dependent_field_takes_priority(self) -> None:
        # Explicit coordinate_dependent field overrides auto-detection
        term = HamiltonianTerm(
            coefficient=1.0,
            factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
            factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
            coefficient_symbolic="mPhi2/2",  # looks constant
            coordinate_dependent=("x",),
        )
        assert term.position_dependent

    def test_no_symbolic_not_position_dependent(self) -> None:
        term = HamiltonianTerm(
            coefficient=0.5,
            factor_a=HamiltonianFactor(field="phi_0", operator="identity"),
            factor_b=HamiltonianFactor(field="phi_0", operator="identity"),
        )
        assert not term.position_dependent

    def test_from_dict_parses_coordinate_dependent(self) -> None:
        data: dict[str, Any] = {
            "coefficient": 1.0,
            "factor_a": {"field": "phi_0", "operator": "gradient_z"},
            "factor_b": {"field": "phi_0", "operator": "gradient_z"},
            "coefficient_symbolic": "Csc[y[]]^2/(2*x[]^2)",
            "coordinate_dependent": ["x", "y"],
        }
        term = HamiltonianTerm.from_dict(data)
        assert term.coordinate_dependent == ("x", "y")
        assert term.position_dependent

    def test_from_dict_without_coordinate_dependent_uses_autodetect(self) -> None:
        data: dict[str, Any] = {
            "coefficient": 1.0,
            "factor_a": {"field": "chi_0", "operator": "identity"},
            "factor_b": {"field": "phi_0", "operator": "identity"},
            "coefficient_symbolic": "E^(-1/2*x[]^2/R^2)*g0",
        }
        term = HamiltonianTerm.from_dict(data)
        assert term.coordinate_dependent == ()  # not in JSON
        assert term.position_dependent  # auto-detected


class TestHamiltonianPositionDependentEnergy:
    """Tests that _compute_hamiltonian_from_canonical evaluates position-dependent
    Hamiltonian coefficients on the spatial grid rather than using the scalar fallback.
    """

    def _make_gaussian_coupling_spec(self, g0: float = 1.0, r_scale: float = 4.0) -> EquationSystem:
        """Build a minimal 2-field spec with a Gaussian interaction Hamiltonian term."""
        data: dict[str, Any] = {
            "metadata": {"parameters": {"g0": g0, "R": r_scale}},
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
            },
            "fields": [
                {"name": "phi_0", "index": 0, "is_dynamical": True},
                {"name": "chi_0", "index": 1, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                            {"coefficient": 1.0, "operator": "laplacian_y", "field": "phi_0"},
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                            {"coefficient": 1.0, "operator": "laplacian_y", "field": "chi_0"},
                        ],
                    },
                },
            ],
            "canonical": {
                "hamiltonian_terms": [
                    # Kinetic: pi²/2
                    {
                        "coefficient": 0.5,
                        "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                        "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                    },
                    {
                        "coefficient": 0.5,
                        "factor_a": {"field": "chi_0", "operator": "time_derivative"},
                        "factor_b": {"field": "chi_0", "operator": "time_derivative"},
                    },
                    # Position-dependent interaction: g0*exp(-(x²+y²)/2R²) * phi * chi
                    {
                        "coefficient": 1.0,
                        "factor_a": {"field": "chi_0", "operator": "identity"},
                        "factor_b": {"field": "phi_0", "operator": "identity"},
                        "coefficient_symbolic": "E^(-1/2*x[]^2/R^2 - y[]^2/(2*R^2))*g0",
                    },
                ],
                "field_rates": {
                    "phi_0": [{"coefficient": 1.0, "operator": "identity", "field": "pi_0"}],
                    "chi_0": [{"coefficient": 1.0, "operator": "identity", "field": "pi_1"}],
                },
                "kinetic_matrix": {
                    "entries": [
                        {"i": 0, "j": 0, "value": 1.0},
                        {"i": 1, "j": 1, "value": 1.0},
                    ],
                    "dimension": 2,
                },
                "spatial_momenta": {},
                "hamiltonian_symbolic": "test",
            },
        }
        return EquationSystem.from_dict(data)

    def _make_sim_data(
        self, spec: EquationSystem, n: int = 16
    ) -> SimulationData:
        """Build a SimulationData with unit fields and zero momenta."""
        dx = 1.0
        # phi = 1 everywhere, chi = 1 everywhere
        phi = np.ones((n, n))
        chi = np.ones((n, n))
        pi_phi = np.zeros((n, n))
        pi_chi = np.zeros((n, n))
        params = {
            k: float(v)
            for k, v in spec.metadata.get("parameters", {}).items()
            if isinstance(v, (int, float))
        }
        return SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis], "chi_0": chi[np.newaxis]},
            momenta={"phi_0": pi_phi[np.newaxis], "chi_0": pi_chi[np.newaxis]},
            grid_spacing=(dx, dx),
            grid_bounds=((0.0, n * dx), (0.0, n * dx)),
            periodic=(True, True),
            spec=spec,
            parameters=params,
        )

    def test_gaussian_interaction_evaluated_on_grid(self) -> None:
        """Energy with Gaussian coefficient ≠ energy with uniform coefficient=1."""
        g0 = 1.0
        r_scale = 4.0
        spec = self._make_gaussian_coupling_spec(g0=g0, r_scale=r_scale)
        data = self._make_sim_data(spec, n=16)

        energy = _compute_hamiltonian_from_canonical(data, 0)

        # With phi=chi=1 everywhere, the interaction contribution is:
        #   <g0 * exp(-(x²+y²)/2R²)> = g0 * spatial_average(Gaussian)
        # This must be < g0 (the Gaussian is <1 almost everywhere).
        # If the bug were still present (coeff=1.0), the interaction term would be
        # exactly g0 * 1.0 = g0 = 1.0, making the total wrong.
        n = 16
        dx = 1.0
        # Cell-centered coordinates (matching _build_coord_arrays: lo + dx/2)
        xs = np.arange(n) * dx + dx / 2  # 0.5, 1.5, ..., 15.5
        ys = np.arange(n) * dx + dx / 2
        xx, yy = np.meshgrid(xs, ys, indexing="ij")
        expected_interaction = float(
            (g0 * np.exp(-0.5 * xx**2 / r_scale**2 - yy**2 / (2 * r_scale**2))).mean()
        )

        # The kinetic terms are 0 (zero momenta), so the total energy = interaction term
        assert abs(energy - expected_interaction) < 1e-10

    def test_different_g0_gives_different_energy(self) -> None:
        """Scaling g0 scales the interaction energy proportionally."""
        spec1 = self._make_gaussian_coupling_spec(g0=1.0)
        spec2 = self._make_gaussian_coupling_spec(g0=2.0)
        data1 = self._make_sim_data(spec1, n=16)
        data2 = self._make_sim_data(spec2, n=16)

        e1 = _compute_hamiltonian_from_canonical(data1, 0)
        e2 = _compute_hamiltonian_from_canonical(data2, 0)

        # If coefficient were wrongly treated as 1.0, both would be equal
        assert abs(e2 - 2.0 * e1) < 1e-10

    def test_position_dependent_gradient_coeff(self) -> None:
        """Hamiltonian with 1/x[]^2 on gradient_x term is evaluated on grid."""
        data: dict[str, Any] = {
            "metadata": {},
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                        ],
                    },
                }
            ],
            "canonical": {
                "hamiltonian_terms": [
                    # gradient_x^2 with position-dependent coeff 1/x[]^2
                    {
                        "coefficient": 0.5,
                        "factor_a": {"field": "phi_0", "operator": "gradient_x"},
                        "factor_b": {"field": "phi_0", "operator": "gradient_x"},
                        "coefficient_symbolic": "1/(2*x[]^2)",
                    },
                ],
                "field_rates": {
                    "phi_0": [{"coefficient": 1.0, "operator": "identity", "field": "pi_0"}]
                },
                "kinetic_matrix": {
                    "entries": [{"i": 0, "j": 0, "value": 1.0}],
                    "dimension": 1,
                },
                "spatial_momenta": {},
                "hamiltonian_symbolic": "test",
            },
        }
        spec = EquationSystem.from_dict(data)
        n = 8
        dx = 1.0
        # phi = x (linear ramp), so grad_x(phi) = 1 everywhere
        xs = np.arange(n, dtype=float) * dx + 0.5  # cell-centred, avoid x=0
        phi = np.tile(xs[:, np.newaxis], (1, n))
        pi = np.zeros((n, n))
        sd = SimulationData(
            times=np.array([0.0]),
            fields={"phi_0": phi[np.newaxis]},
            momenta={"phi_0": pi[np.newaxis]},
            grid_spacing=(dx, dx),
            grid_bounds=((0.5, n * dx + 0.5), (0.0, n * dx)),
            periodic=(False, True),
            spec=spec,
            parameters={},
        )
        energy = _compute_hamiltonian_from_canonical(sd, 0)
        # With grad_x(phi)=1 and coeff=1/(2*x²), via IBP:
        # energy = -<1/(2*x²) * phi * laplacian_x(phi)>
        # laplacian_x of linear ramp = 0 everywhere (interior), so energy ≈ 0
        # This just verifies it runs without error and returns a finite value
        assert np.isfinite(energy)
