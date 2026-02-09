"""Physics validation tests: energy conservation and analytical solutions.

These tests verify that the PDE pipeline produces physically correct results,
not just structurally valid output. They catch sign errors, wrong operator
implementations, and numerical issues that structural tests miss.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField, VectorField

from tidal.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
)
from tidal.symbolic.pde_builder import PDEFromSpec, create_initial_state

if TYPE_CHECKING:
    from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _wave_hamiltonian(state: FieldCollection, grid: CartesianGrid) -> float:
    """Compute H = ½ ∫ (π² + |∇φ|²) dx for the 1D wave equation.

    Parameters
    ----------
    state : FieldCollection
        [phi, pi] field collection.
    grid : CartesianGrid
        The simulation grid.

    Returns
    -------
    float
        Total energy (Hamiltonian).
    """
    phi = cast("ScalarField", state[0])
    pi = cast("ScalarField", state[1])
    dx = grid.discretization[0]

    # Kinetic: ½ π²
    kinetic = 0.5 * np.sum(pi.data**2) * dx

    # Gradient: ½ |∇φ|², using centered differences
    grad: VectorField = phi.gradient(bc="periodic")
    grad_sq = sum(np.sum(g.data**2) for g in grad) * dx  # type: ignore[union-attr]

    return float(kinetic + 0.5 * grad_sq)


def _kg_hamiltonian(state: FieldCollection, grid: CartesianGrid, m2: float) -> float:
    """Compute H = ½ ∫ (π² + |∇φ|² + m²φ²) dx for the Klein-Gordon equation."""
    phi = cast("ScalarField", state[0])
    pi = cast("ScalarField", state[1])
    dx = grid.discretization[0]

    kinetic = 0.5 * np.sum(pi.data**2) * dx
    grad: VectorField = phi.gradient(bc="periodic")
    gradient_energy = 0.5 * sum(np.sum(g.data**2) for g in grad) * dx  # type: ignore[union-attr]
    mass_energy = 0.5 * m2 * np.sum(phi.data**2) * dx

    return float(kinetic + gradient_energy + mass_energy)


# ---------------------------------------------------------------------------
# Wave equation spec (massless Klein-Gordon)
# ---------------------------------------------------------------------------


@pytest.fixture
def wave_spec() -> EquationSystem:
    """1D wave equation: ∂²φ/∂t² = ∂²φ/∂x² (Laplacian)."""
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"description": "1D wave equation"},
    )


@pytest.fixture
def kg_spec() -> EquationSystem:
    """1D Klein-Gordon: ∂²φ/∂t² = ∇²φ - m²φ."""
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(1.0, "laplacian", "phi"),
                    OperatorTerm(-1.0, "identity", "phi", "-m2"),
                ),
            ),
        ),
        mass_matrix=((1.0,),),
        coupling_matrix=((0.0,),),
        mass_matrix_symbolic=(("-m2",),),
        metadata={"description": "1D Klein-Gordon"},
    )


# ---------------------------------------------------------------------------
# #10: Energy conservation tests
# ---------------------------------------------------------------------------


class TestEnergyConservation:
    """Verify Hamiltonian conservation for conservative systems."""

    def test_wave_energy_conserved(self, wave_spec: EquationSystem) -> None:
        """Wave equation H = ½∫(π² + |∇φ|²)dx is conserved over 200 steps.

        Uses explicit Euler with small dt to keep integration error small.
        Tolerance allows for O(dt²) error accumulation.
        """
        pde = PDEFromSpec(wave_spec)
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        x = cast("np.ndarray", grid.cell_coords[..., 0])

        # Initial: traveling Gaussian pulse
        sigma = 0.5
        x0 = 5.0
        phi_data = np.exp(-((x - x0) ** 2) / (2 * sigma**2))
        pi_data = np.zeros_like(phi_data)

        state = FieldCollection(
            [
                ScalarField(grid, data=phi_data),
                ScalarField(grid, data=pi_data),
            ]
        )

        h_initial = _wave_hamiltonian(state, grid)
        assert h_initial > 0, "Initial energy should be positive"

        # Evolve with small dt for 200 steps
        dt = 0.001
        for _ in range(200):
            rates = pde.evolution_rate(state, t=0.0)
            new_fields: list[ScalarField] = []
            for field, rate in zip(state, rates, strict=False):
                new_data = field.data + dt * rate.data
                new_fields.append(ScalarField(grid, data=new_data))
            state = FieldCollection(new_fields)

        h_final = _wave_hamiltonian(state, grid)

        # Energy should be conserved within ~1% (Euler method, small dt)
        relative_change = abs(h_final - h_initial) / h_initial
        assert relative_change < 0.02, (
            f"Energy not conserved: H_0={h_initial:.6f}, H_f={h_final:.6f}, "
            f"relative change={relative_change:.4f}"
        )

    def test_kg_energy_conserved(self, kg_spec: EquationSystem) -> None:
        """Klein-Gordon H = ½∫(π² + |∇φ|² + m²φ²)dx is conserved."""
        m2 = 2.0
        pde = PDEFromSpec(kg_spec, parameters={"m2": m2})
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        x = cast("np.ndarray", grid.cell_coords[..., 0])

        # Gaussian initial data
        phi_data = np.exp(-((x - 5.0) ** 2) / (2 * 0.5**2))
        state = FieldCollection(
            [
                ScalarField(grid, data=phi_data),
                ScalarField(grid, data=0.0),
            ]
        )

        h_initial = _kg_hamiltonian(state, grid, m2)

        dt = 0.001
        for _ in range(200):
            rates = pde.evolution_rate(state, t=0.0)
            new_fields: list[ScalarField] = []
            for field, rate in zip(state, rates, strict=False):
                new_fields.append(ScalarField(grid, data=field.data + dt * rate.data))
            state = FieldCollection(new_fields)

        h_final = _kg_hamiltonian(state, grid, m2)
        relative_change = abs(h_final - h_initial) / h_initial
        assert relative_change < 0.02, (
            f"KG energy not conserved: relative change={relative_change:.4f}"
        )


# ---------------------------------------------------------------------------
# #11: Analytical solution comparison tests
# ---------------------------------------------------------------------------


class TestAnalyticalSolutions:
    """Compare numerical output against known analytical solutions."""

    def test_wave_dispersion_relation(self, wave_spec: EquationSystem) -> None:
        """For ∂²φ/∂t² = ∇²φ, a plane wave φ = cos(kx) has ∂²φ/∂t² = -k²cos(kx).

        This directly tests that the RHS = ∇²φ = -k²φ for the eigenfunction.
        """
        pde = PDEFromSpec(wave_spec)
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        k = 2 * np.pi / 10  # wavenumber fitting the domain

        # Initial: φ = cos(kx), π = 0
        phi_data = np.cos(k * x)
        state = FieldCollection(
            [
                ScalarField(grid, data=phi_data),
                ScalarField(grid, data=0.0),
            ]
        )

        rates = pde.evolution_rate(state, t=0.0)

        # d_t φ = π = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d_t π = ∇²φ = -k² cos(kx)
        expected_accel = -(k**2) * np.cos(k * x)
        assert_allclose(rates[1].data, expected_accel, rtol=1e-3)

    def test_kg_dispersion_relation(self, kg_spec: EquationSystem) -> None:
        """For ∂²φ/∂t² = ∇²φ - m²φ, plane wave has ω² = k² + m².

        The acceleration should be -(k² + m²) cos(kx).
        """
        m2 = 4.0
        pde = PDEFromSpec(kg_spec, parameters={"m2": m2})
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        k = 2 * np.pi / 10

        phi_data = np.cos(k * x)
        state = FieldCollection(
            [
                ScalarField(grid, data=phi_data),
                ScalarField(grid, data=0.0),
            ]
        )

        rates = pde.evolution_rate(state, t=0.0)

        # d_t π = ∇²φ - m²φ = -k²cos(kx) - m²cos(kx) = -(k²+m²)cos(kx)
        omega_sq = k**2 + m2
        expected = -omega_sq * np.cos(k * x)
        assert_allclose(rates[1].data, expected, rtol=1e-3)

    def test_standing_wave_quarter_period(self, wave_spec: EquationSystem) -> None:
        """A standing wave φ(0) = cos(kx), π(0) = 0 evolves to φ ≈ 0, π ≈ -ω sin(kx)
        after a quarter period T/4 = π/(2ω) where ω = k.

        Uses Störmer-Verlet (leapfrog) for better energy conservation.
        """
        pde = PDEFromSpec(wave_spec)
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        k = 2 * np.pi / 10
        omega = k  # dispersion: ω = k for wave equation

        phi = ScalarField(grid, data=np.cos(k * x))
        pi = ScalarField(grid, data=0.0)

        # Störmer-Verlet (symplectic) integration for T/4
        quarter_period = np.pi / (2 * omega)
        dt = 0.005
        n_steps = int(quarter_period / dt)

        for _ in range(n_steps):
            # Half-step momentum
            state = FieldCollection([phi, pi])
            rates = pde.evolution_rate(state, t=0.0)
            pi_half = ScalarField(grid, data=pi.data + 0.5 * dt * rates[1].data)

            # Full-step position
            phi_new = ScalarField(grid, data=phi.data + dt * pi_half.data)

            # Half-step momentum (with new position)
            state = FieldCollection([phi_new, pi_half])
            rates = pde.evolution_rate(state, t=0.0)
            pi_new = ScalarField(grid, data=pi_half.data + 0.5 * dt * rates[1].data)

            phi = phi_new
            pi = pi_new

        # Analytical: φ(t) = cos(kx)cos(ωt), π(t) = -ω cos(kx)sin(ωt)
        # At T/4: φ ≈ 0, π ≈ -ω cos(kx)
        assert np.max(np.abs(phi.data)) < 0.15, (
            f"φ amplitude should be ~0 at T/4, got max |φ|={np.max(np.abs(phi.data)):.4f}"
        )

        # π should match -ω cos(kx)
        expected_pi = -omega * np.cos(k * x)
        # Normalize to check shape agreement (amplitude may differ slightly)
        pi_norm = np.max(np.abs(pi.data))
        expected_norm = np.max(np.abs(expected_pi))
        assert_allclose(
            pi.data / pi_norm,
            expected_pi / expected_norm,
            atol=0.15,
        )


# ---------------------------------------------------------------------------
# #20: CFL / stability diagnostics
# ---------------------------------------------------------------------------


class TestCFLStability:
    """Verify CFL condition checker detects unstable configurations."""

    def test_stable_dt_returns_no_warnings(self, wave_spec: EquationSystem) -> None:
        """Small dt relative to dx produces no warnings."""
        pde = PDEFromSpec(wave_spec)
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        # dx = 10/128 ≈ 0.078, c=1 → CFL limit ≈ 0.078
        dt = 0.01
        warnings = pde.check_stability(dt, grid)
        assert warnings == []

    def test_unstable_dt_returns_warning(self, wave_spec: EquationSystem) -> None:
        """Large dt relative to dx triggers CFL warning."""
        pde = PDEFromSpec(wave_spec)
        grid = CartesianGrid([(0, 10)], 128, periodic=True)
        # dx ≈ 0.078, c=1 → CFL limit ≈ 0.078, dt=0.5 >> limit
        dt = 0.5
        warnings = pde.check_stability(dt, grid)
        assert len(warnings) == 1
        assert "CFL violated" in warnings[0]
        assert "phi" in warnings[0]

    def test_kg_cfl_includes_wave_speed(self, kg_spec: EquationSystem) -> None:
        """KG equation CFL check uses the Laplacian coefficient."""
        pde = PDEFromSpec(kg_spec, parameters={"m2": 4.0})
        grid = CartesianGrid([(0, 10)], 64, periodic=True)
        # dx = 10/64 ≈ 0.156, laplacian coeff = 1.0 → c = 1.0
        # CFL limit ≈ 0.156
        dt = 0.001
        warnings = pde.check_stability(dt, grid)
        assert warnings == []


# ---------------------------------------------------------------------------
# Massive 3-form: rank-3 antisymmetric tensor (symmetry reduction example)
# ---------------------------------------------------------------------------


class TestMassive3Form:
    """Validate the massive 3-form example (rank-3 antisymmetric tensor).

    This exercises the tensor symmetry reduction pipeline
    (EnumerateComponentTuples) and verifies that each independent component
    satisfies the Klein-Gordon equation independently.
    """

    @pytest.fixture
    def massive_3form_spec(self) -> EquationSystem:
        """Load the massive 3-form JSON generated by the Wolfram pipeline."""
        from pathlib import Path

        from tidal.symbolic.json_loader import (
            load_equation_system,
        )

        json_path = (
            Path(__file__).parent.parent / "examples" / "data" / "massive_3form.json"
        )
        if not json_path.exists():
            pytest.skip("massive_3form.json not found (run Wolfram script first)")
        return load_equation_system(json_path)

    def test_component_count(self, massive_3form_spec: EquationSystem) -> None:
        """Antisymmetric rank-3 in 4D: 64 raw → 4 independent components."""
        assert massive_3form_spec.n_components == 4
        assert massive_3form_spec.dimension == 4
        assert massive_3form_spec.spatial_dimension == 3
        assert len(massive_3form_spec.equations) == 4

    def test_kg_dispersion_all_components(
        self, massive_3form_spec: EquationSystem
    ) -> None:
        """Each component satisfies KG dispersion: ω² = k² + m².

        For a plane wave φ = cos(k·x), the acceleration should be
        -(k_x² + k_y² + k_z² + m²) cos(k·x) for each component.
        """
        m2 = 2.0
        pde = PDEFromSpec(massive_3form_spec, parameters={"m2": m2})
        grid = CartesianGrid([(0, 10), (0, 10), (0, 10)], 16, periodic=True)
        coords: np.ndarray = grid.cell_coords  # type: ignore[reportUnknownVariableType]
        kx = 2 * np.pi / 10  # one wavelength across the domain

        # State layout is interleaved: [C_0, pi_0, C_1, pi_1, C_2, pi_2, C_3, pi_3]
        # Set all C_i = cos(kx * x), all pi_i = 0
        cos_data: NDArray[np.float64] = np.cos(kx * coords[..., 0])  # type: ignore[reportUnknownArgumentType]
        names = massive_3form_spec.component_names
        field_data: dict[str, NDArray[np.float64]] = dict.fromkeys(names, cos_data)
        state = create_initial_state(grid, massive_3form_spec, field_data=field_data)

        rates = pde.evolution_rate(state, t=0.0)

        # For each component, check the interleaved rates:
        # rate[2*i]   = d_t C_i = pi_i = 0
        # rate[2*i+1] = d_t pi_i = ∇²C_i - m²C_i
        omega_sq = kx**2 + m2
        expected = -omega_sq * np.cos(kx * coords[..., 0])  # type: ignore[reportUnknownArgumentType]

        for i in range(len(names)):
            assert_allclose(rates[2 * i].data, 0.0, atol=1e-10)
            assert_allclose(rates[2 * i + 1].data, expected, rtol=0.1)

    def test_component_independence(self, massive_3form_spec: EquationSystem) -> None:
        """Components are decoupled: exciting C_0 only leaves C_1..C_3 unchanged.

        Since there's no cross-field coupling, the evolution rate for
        components 1-3 should be zero when only component 0 has initial data.
        """
        m2 = 1.0
        pde = PDEFromSpec(massive_3form_spec, parameters={"m2": m2})
        grid = CartesianGrid([(0, 10), (0, 10), (0, 10)], 16, periodic=True)
        coords: np.ndarray = grid.cell_coords  # type: ignore[reportUnknownVariableType]

        # State layout: [C_0, pi_0, C_1, pi_1, C_2, pi_2, C_3, pi_3]
        # Set only C_0 = Gaussian, everything else = 0
        names = massive_3form_spec.component_names
        gaussian: NDArray[np.float64] = np.exp(  # type: ignore[reportUnknownArgumentType]
            -np.sum((coords - 5.0) ** 2, axis=-1) / 2  # type: ignore[reportUnknownArgumentType]
        )
        field_data = {names[0]: gaussian}
        state = create_initial_state(grid, massive_3form_spec, field_data=field_data)

        rates = pde.evolution_rate(state, t=0.0)

        # rate[0] = d_t C_0 = pi_0 = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # rate[1] = d_t pi_0 should be non-zero (Gaussian has curvature)
        assert np.max(np.abs(rates[1].data)) > 0.01

        # C_1, C_2, C_3 and their momenta should all be exactly zero
        for i in range(1, len(names)):
            assert_allclose(
                rates[2 * i].data,
                0.0,
                atol=1e-14,
                err_msg=f"C_{i} rate should be zero",
            )
            assert_allclose(
                rates[2 * i + 1].data,
                0.0,
                atol=1e-14,
                err_msg=f"pi_{i} rate should be zero",
            )
