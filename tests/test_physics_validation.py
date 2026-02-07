"""Physics validation tests: energy conservation and analytical solutions.

These tests verify that the PDE pipeline produces physically correct results,
not just structurally valid output. They catch sign errors, wrong operator
implementations, and numerical issues that structural tests miss.
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
)
from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec


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
    phi = state[0]
    pi = state[1]
    dx = grid.discretization[0]

    # Kinetic: ½ π²
    kinetic = 0.5 * np.sum(pi.data**2) * dx

    # Gradient: ½ |∇φ|², using centered differences
    grad = phi.gradient(bc="periodic")
    grad_sq = sum(np.sum(g.data**2) for g in grad) * dx  # type: ignore[union-attr]

    return float(kinetic + 0.5 * grad_sq)


def _kg_hamiltonian(
    state: FieldCollection, grid: CartesianGrid, m2: float
) -> float:
    """Compute H = ½ ∫ (π² + |∇φ|² + m²φ²) dx for the Klein-Gordon equation."""
    phi = state[0]
    pi = state[1]
    dx = grid.discretization[0]

    kinetic = 0.5 * np.sum(pi.data**2) * dx
    grad = phi.gradient(bc="periodic")
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

        state = FieldCollection([
            ScalarField(grid, data=phi_data),
            ScalarField(grid, data=pi_data),
        ])

        h_initial = _wave_hamiltonian(state, grid)
        assert h_initial > 0, "Initial energy should be positive"

        # Evolve with small dt for 200 steps
        dt = 0.001
        for _ in range(200):
            rates = pde.evolution_rate(state, t=0.0)
            new_fields = []
            for field, rate in zip(state, rates):
                new_data = field.data + dt * rate.data
                new_fields.append(ScalarField(grid, data=new_data))
            state = FieldCollection(new_fields)

        h_final = _wave_hamiltonian(state, grid)

        # Energy should be conserved within ~1% (Euler method, small dt)
        relative_change = abs(h_final - h_initial) / h_initial
        assert relative_change < 0.02, (  # noqa: PLR2004
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
        state = FieldCollection([
            ScalarField(grid, data=phi_data),
            ScalarField(grid, data=0.0),
        ])

        h_initial = _kg_hamiltonian(state, grid, m2)

        dt = 0.001
        for _ in range(200):
            rates = pde.evolution_rate(state, t=0.0)
            new_fields = []
            for field, rate in zip(state, rates):
                new_fields.append(ScalarField(grid, data=field.data + dt * rate.data))
            state = FieldCollection(new_fields)

        h_final = _kg_hamiltonian(state, grid, m2)
        relative_change = abs(h_final - h_initial) / h_initial
        assert relative_change < 0.02, (  # noqa: PLR2004
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
        state = FieldCollection([
            ScalarField(grid, data=phi_data),
            ScalarField(grid, data=0.0),
        ])

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
        state = FieldCollection([
            ScalarField(grid, data=phi_data),
            ScalarField(grid, data=0.0),
        ])

        rates = pde.evolution_rate(state, t=0.0)

        # d_t π = ∇²φ - m²φ = -k²cos(kx) - m²cos(kx) = -(k²+m²)cos(kx)
        omega_sq = k**2 + m2
        expected = -omega_sq * np.cos(k * x)
        assert_allclose(rates[1].data, expected, rtol=1e-3)

    def test_standing_wave_quarter_period(
        self, wave_spec: EquationSystem
    ) -> None:
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
        assert np.max(np.abs(phi.data)) < 0.15, (  # noqa: PLR2004
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
