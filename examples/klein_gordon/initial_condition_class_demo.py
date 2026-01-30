"""Demonstration of class-based InitialCondition API.

This example shows how to use the new InitialCondition base class pattern
for creating initial conditions. This approach provides:

1. Better code reuse through inheritance
2. Elimination of coordinate handling duplication
3. Easier composition and customization
4. Type-safe parameter validation

Compare with older function-based examples to see the improvements.
"""

from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from torsion_gertsenshtein.kgsim import (
    GaussianPulse,
    GridConfig,
    InitialCondition,
    KGParameters,
    KleinGordonPDE,
    RingPulse2D,
    SimulationConfig,
    make_grid,
    run_with_snapshots,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


class DoubleGaussianPulse(InitialCondition):
    """Custom initial condition: two Gaussian pulses.

    This demonstrates how easy it is to create custom ICs by inheriting
    from the InitialCondition base class.
    """

    def __init__(  # noqa: PLR0913 - Demo class showing multiple parameters
        self,
        *,
        amplitude1: float,
        center1: tuple[float, ...],
        width1: float,
        amplitude2: float,
        center2: tuple[float, ...],
        width2: float,
    ) -> None:
        """Initialize two-pulse parameters."""
        self.amplitude1 = amplitude1
        self.center1 = center1
        self.width1 = width1
        self.amplitude2 = amplitude2
        self.center2 = center2
        self.width2 = width2

    def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
        """Compute sum of two Gaussians."""
        # First pulse
        r1 = self._compute_distances_from_center(grid, self.center1)
        pulse1 = self.amplitude1 * np.exp(-(r1**2) / (2.0 * self.width1**2))

        # Second pulse
        r2 = self._compute_distances_from_center(grid, self.center2)
        pulse2 = self.amplitude2 * np.exp(-(r2**2) / (2.0 * self.width2**2))

        return pulse1 + pulse2


def demo_gaussian_pulse_1d() -> None:
    """Demonstrate basic GaussianPulse class usage in 1D."""
    print("\n" + "=" * 70)
    print("DEMO 1: GaussianPulse class in 1D")
    print("=" * 70)

    # Setup
    grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0.0, 100.0),)))

    # Create initial condition using class-based API
    ic = GaussianPulse(amplitude=1.0, width=3.0, center=[50.0], initial_velocity=0.0)
    state = ic.build(grid)

    print("✓ Created GaussianPulse IC")
    print(f"  Amplitude: {ic.amplitude}")
    print(f"  Width: {ic.width}")
    print(f"  Center: {ic.center}")
    print(f"  Max φ: {state[0].data.max():.4f}")

    # Quick simulation
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=20.0, solver="scipy", backend="numpy", progress=False
    )
    result, _ = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=1.0
    )

    print(f"✓ Simulation complete, final max φ: {result[0].data.max():.4f}")


def demo_ring_pulse_2d() -> None:
    """Demonstrate RingPulse2D class usage."""
    print("\n" + "=" * 70)
    print("DEMO 2: RingPulse2D class in 2D")
    print("=" * 70)

    # Setup
    grid = make_grid(
        GridConfig(dim=2, shape=(64, 64), bounds=((-20.0, 20.0), (-20.0, 20.0)))
    )

    # Create ring initial condition
    ic = RingPulse2D(amplitude=1.0, initial_radius=8.0, sigma=1.5)
    state = ic.build(grid)

    print("✓ Created RingPulse2D IC")
    print(f"  Amplitude: {ic.amplitude}")
    print(f"  Ring radius: {ic.initial_radius}")
    print(f"  Width: {ic.width}")
    print(f"  Max φ: {state[0].data.max():.4f}")

    # Visualize
    _fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        state[0].data.T, origin="lower", extent=(-20, 20, -20, 20), cmap="RdBu_r"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("RingPulse2D Initial Condition")
    ax.set_aspect("equal")
    plt.colorbar(im, ax=ax, label="φ")

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out_path = "outputs/ring_pulse_2d_ic.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved {out_path}")


def demo_custom_double_pulse() -> None:
    """Demonstrate custom DoubleGaussianPulse IC."""
    print("\n" + "=" * 70)
    print("DEMO 3: Custom DoubleGaussianPulse (inheritance)")
    print("=" * 70)

    # Setup
    grid = make_grid(
        GridConfig(dim=2, shape=(64, 64), bounds=((-20.0, 20.0), (-20.0, 20.0)))
    )

    # Create custom two-pulse IC
    ic = DoubleGaussianPulse(
        amplitude1=1.0,
        center1=(-8.0, 0.0),
        width1=2.0,
        amplitude2=0.8,
        center2=(8.0, 0.0),
        width2=2.5,
    )
    state = ic.build(grid)

    print("✓ Created custom DoubleGaussianPulse IC")
    print(f"  Pulse 1: A={ic.amplitude1}, center={ic.center1}, width={ic.width1}")
    print(f"  Pulse 2: A={ic.amplitude2}, center={ic.center2}, width={ic.width2}")
    print(f"  Max φ: {state[0].data.max():.4f}")

    # Visualize
    _fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        state[0].data.T, origin="lower", extent=(-20, 20, -20, 20), cmap="viridis"
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Custom DoubleGaussianPulse Initial Condition")
    ax.set_aspect("equal")
    plt.colorbar(im, ax=ax, label="φ")

    out_path = "outputs/double_gaussian_ic.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved {out_path}")


def main() -> None:
    """Run all InitialCondition class demonstrations."""
    print("\n" + "=" * 70)
    print("INITIAL CONDITION CLASS-BASED API DEMONSTRATION")
    print("=" * 70)
    print("""
The new InitialCondition base class provides:

✓ Elimination of coordinate handling duplication
✓ Clean inheritance for custom ICs
✓ Type-safe parameter validation
✓ Consistent API across all IC types
✓ Easy composition and customization
    """)

    demo_gaussian_pulse_1d()
    demo_ring_pulse_2d()
    demo_custom_double_pulse()

    print("\n" + "=" * 70)
    print("KEY BENEFITS")
    print("=" * 70)
    print("""
1. CODE REUSE: The base class handles all coordinate extraction logic
   - No more duplicated coordinate handling across 10+ IC functions
   - Consistent behavior across different py-pde versions

2. EXTENSIBILITY: Create custom ICs by inheriting and implementing _compute_phi()
   - Only need to define the field profile, not infrastructure
   - Example: DoubleGaussianPulse above is only ~20 lines

3. TYPE SAFETY: Parameter validation in __init__
   - Catches errors at IC creation, not during simulation
   - Better error messages

4. BACKWARD COMPATIBILITY: Original function-based API still works
   - Gradual migration path
   - No breaking changes for existing code
    """)


if __name__ == "__main__":
    main()
