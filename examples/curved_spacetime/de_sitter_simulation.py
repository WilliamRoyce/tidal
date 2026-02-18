"""De Sitter Klein-Gordon 2+1D simulation - Hubble friction in expanding spacetime.

This script demonstrates a scalar field on 2+1D de Sitter (exponentially expanding)
spacetime, where Hubble friction genuinely appears.

Physics:
  Metric: ds^2 = Omega(t)^2 * (-dt^2 + dx^2 + dy^2)  where Omega(t) = exp(H*t)

  The conformal d'Alembertian in (n+1) dimensions with Omega = exp(Ht):
    Box_g phi = Omega^{-2}[-d2_t phi + nabla^2 phi - (n-1)*H*d_t phi]

  For 2+1D (n=2 spatial dimensions):
    d2_t phi = nabla^2 phi - H * d_t phi - m^2 * exp(2Ht) * phi

  Key observation:
  - Hubble friction: -(n-1)*H * d_t phi = -H * d_t phi for n=2
  - This is genuine physics - in 1+1D, friction = 0 due to conformal invariance
  - Energy is NOT conserved - waves lose energy to expansion
  - Amplitude decays due to Hubble damping

Verification:
  Compare with flat space KG - de Sitter should show:
  1. Amplitude decay (Hubble damping)
  2. Energy loss to cosmic expansion
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CallbackTracker, CartesianGrid, FieldCollection, ScalarField

from tidal.measurement import (
    SimulationData,
    compute_energy_timeseries,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
HUBBLE_PARAM = 0.1  # Hubble constant (expansion rate)
MASS_SQUARED = 1.0  # Mass parameter

# Grid
X_MIN = 0.0
X_MAX = 50.0
Y_MIN = 0.0
Y_MAX = 50.0
GRID_SHAPE = [64, 64]
GRID_PERIODIC: list[bool] = [True, True]

# Time integration
T_END = 20.0
SNAPSHOT_INTERVAL = 0.2

# Initial conditions
CENTER_X = 25.0
CENTER_Y = 25.0
PULSE_WIDTH = 3.0
PULSE_AMPLITUDE = 1.0

# Analysis
AMPLITUDE_DECAY_THRESHOLD = 0.9  # Threshold to detect Hubble friction effect

# Output
OUTPUT_FILENAME = "de_sitter_kg_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for de Sitter KG simulation results."""

    grid: CartesianGrid
    data: SimulationData
    x_coords: NumericArray
    y_coords: NumericArray


def _print_header() -> None:
    print("=" * 60)
    print("De Sitter Klein-Gordon 2+1D Simulation (Hubble Friction)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = exp(2Ht) * (-dt^2 + dx^2 + dy^2)")
    print(
        "Wave equation (2+1D): d2_t phi = nabla^2 phi - H * d_t phi - m^2 * exp(2Ht) * phi"
    )
    print("                                               ^^^ (n-1)H = H for n=2")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (2+1D)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        terms_str = " + ".join(
            f"{term.coefficient:.1f}*{term.operator}({term.field})"
            for term in eq.rhs_terms
        )
        print(f"    d2_t({eq.field_name}) = {terms_str}")
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(
        json_path, parameters={"dSH": HUBBLE_PARAM, "dSm2": MASS_SQUARED}
    )
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Hubble parameter: H = {HUBBLE_PARAM}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 2D simulation grid...")
    grid = CartesianGrid(
        bounds=[(X_MIN, X_MAX), (Y_MIN, Y_MAX)],
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  Domain: [{X_MIN:.0f}, {X_MAX:.0f}] x [{Y_MIN:.0f}, {Y_MAX:.0f}]")
    print(f"  Resolution: {grid.shape[0]}x{grid.shape[1]} points")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    # Initialize a 2D Gaussian pulse at center
    gaussian = PULSE_AMPLITUDE * np.exp(
        -((x - CENTER_X) ** 2 + (y - CENTER_Y) ** 2) / (2 * PULSE_WIDTH**2)
    )

    phi = ScalarField(grid, data=gaussian, label="dSphi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")
    state = FieldCollection([phi, pi])

    print(f"  2D Gaussian pulse at center ({CENTER_X}, {CENTER_Y})")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print("  Initial momentum: zero")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection, json_path: Path
) -> SimulationResult:
    print("Step 5: Running simulation...")

    spec = load_equation_system(json_path)
    params = {"dSH": HUBBLE_PARAM, "dSm2": MASS_SQUARED}
    output_data_dir = Path(__file__).parent.parent.parent / "outputs" / "de_sitter"

    writer, callback = create_snapshot_callback(
        output_dir=output_data_dir,
        spec=spec,
        grid=grid,
        t_end=T_END,
        snapshot_interval=SNAPSHOT_INTERVAL,
        parameters=params,
        spec_path=json_path,
    )
    tracker = CallbackTracker(callback, interrupts=SNAPSHOT_INTERVAL)
    pde.solve(
        state,
        t_range=T_END,
        solver="scipy",
        method="DOP853",
        tracker=tracker,
    )
    writer.close()

    data = SimulationData.from_directory(output_data_dir, spec)
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    print(f"  Duration: {T_END} time units")
    print(f"  {data.n_snapshots} snapshots saved to {output_data_dir}")
    print()

    return SimulationResult(
        grid=grid,
        data=data,
        x_coords=x,
        y_coords=y,
    )


def _analyze_results(result: SimulationResult) -> dict[str, float]:
    print("Step 6: Analyzing results...")

    data = result.data
    field_name = next(iter(data.fields))

    # Extract field amplitudes
    initial_max = float(np.max(np.abs(data.fields[field_name][0])))
    final_max = float(np.max(np.abs(data.fields[field_name][-1])))

    print(f"  Initial: max|phi| = {initial_max:.4f}")
    print(f"  Final:   max|phi| = {final_max:.4f}")

    # Calculate amplitude decay ratio
    decay_ratio = final_max / initial_max if initial_max > 0 else 0.0
    expected_decay = float(np.exp(-0.5 * HUBBLE_PARAM * T_END))

    print(f"  Amplitude ratio (final/initial): {decay_ratio:.4f}")
    print(f"  Reference decay exp(-H*t/2): {expected_decay:.4f}")

    if decay_ratio < AMPLITUDE_DECAY_THRESHOLD:
        print("  Hubble friction observed: amplitude has decayed")
    else:
        print("  Note: For small H or short time, damping may be subtle")

    # Energy evolution via measurement module (NOT conserved for de Sitter)
    print("Computing energy evolution (virial formula)...")
    _times, _per_field, _interaction, total = compute_energy_timeseries(data)
    e0 = total[0]
    e_final = total[-1]
    print(f"  Initial energy density (virial): {e0:.6f}")
    print(f"  Final energy density (virial):   {e_final:.6f}")
    if e0 > 0:
        print(f"  Energy loss: {100 * (1 - e_final / e0):.1f}%")
    print()

    return {"decay_ratio": decay_ratio}


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    data = result.data
    grid = result.grid
    field_name = next(iter(data.fields))
    n_snapshots = data.n_snapshots

    initial_phi = data.fields[field_name][0]
    final_phi = data.fields[field_name][-1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initial field (2D)
    ax = axes[0, 0]
    im = ax.imshow(
        initial_phi.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-PULSE_AMPLITUDE,
        vmax=PULSE_AMPLITUDE,
        extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Initial $\phi$ (t=0)")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    # Final field (2D)
    ax = axes[0, 1]
    im = ax.imshow(
        final_phi.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-PULSE_AMPLITUDE,
        vmax=PULSE_AMPLITUDE,
        extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(rf"Final $\phi$ (t={T_END:.0f})")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    # Amplitude and field-norm decay (monotonic indicators of Hubble damping)
    ax = axes[1, 0]
    time_values = list(data.times)
    amplitudes = [
        float(np.max(np.abs(data.fields[field_name][i])))
        for i in range(n_snapshots)
    ]
    field_norms = [
        float(np.sum(data.fields[field_name][i] ** 2))
        for i in range(n_snapshots)
    ]

    t_array = np.array(time_values)
    ax.plot(
        time_values,
        np.array(amplitudes) / amplitudes[0],
        "b-",
        linewidth=2,
        label=r"max$|\phi|$ / max$|\phi_0|$",
    )
    ax.plot(
        time_values,
        np.array(field_norms) / field_norms[0],
        "g-",
        linewidth=2,
        label=r"$\int\phi^2 / \int\phi_0^2$",
    )

    # Reference: amplitude decay exp(-Ht/2) for weakly damped oscillator
    ax.plot(
        time_values,
        np.exp(-0.5 * HUBBLE_PARAM * t_array),
        "r--",
        linewidth=1.5,
        label=r"$e^{-Ht/2}$",
    )

    ax.set_xlabel("Time")
    ax.set_ylabel("Normalized quantity")
    ax.set_title(f"Hubble damping (H={HUBBLE_PARAM})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Cross-section through center (y = center_y)
    ax = axes[1, 1]
    center_idx = grid.shape[1] // 2
    times_to_plot = [
        0,
        n_snapshots // 4,
        n_snapshots // 2,
        3 * n_snapshots // 4,
        n_snapshots - 1,
    ]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times_to_plot)))

    x_1d = cast("np.ndarray", grid.cell_coords[:, 0, 0])
    for i, t_idx in enumerate(times_to_plot):
        ax.plot(
            x_1d,
            data.fields[field_name][t_idx][:, center_idx],
            color=colors[i],
            label=f"t={time_values[t_idx]:.0f}",
            alpha=0.8,
        )
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"$\phi$ cross-section at y=25 (note amplitude decay)")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        rf"De Sitter KG 2+1D: $g_{{\mu\nu}} = e^{{2Ht}} \eta_{{\mu\nu}}$, H={HUBBLE_PARAM}, m$^2$={MASS_SQUARED}"
        "\nHubble friction causes energy loss to expansion",
        fontsize=12,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved plot to: {output_path}")
    plt.close()
    print()


def _print_footer(decay_ratio: float) -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print(f"  1. De Sitter metric: g_ab = exp(2Ht) * eta_ab with H={HUBBLE_PARAM}")
    print("  2. Hubble friction term: -(n-1)*H * d_t phi where n = spatial dimensions")
    print("     For this 2+1D simulation (n=2): -H * d_t phi")
    print("     Note: In 1+1D, friction = 0 (conformal invariance)")
    print("  3. Christoffel symbols computed from metric definition (not hardcoded)")
    print("  4. Energy is NOT conserved - lost to cosmic expansion")
    print(f"  5. Amplitude decayed by factor ~{decay_ratio:.3f} over t={T_END}")
    print()
    print("  Physical interpretation:")
    print("    Waves in expanding spacetime lose energy to expansion.")
    print("    The Hubble friction term acts like a velocity-dependent damping.")
    print("    This is why the CMB photons have redshifted since the Big Bang.")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the 2+1D de Sitter Klein-Gordon simulation.

    Raises
    ------
    FileNotFoundError
        If the JSON spec has not been generated yet.
    """
    json_path = Path(__file__).parent.parent / "data" / "de_sitter_kg.json"
    if not json_path.exists():
        msg = (
            f"{json_path} not found. "
            "Run 'tidal derive theory.toml' first (requires wolframscript)."
        )
        raise FileNotFoundError(msg)
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid)
    result = _run_simulation(pde, grid, state, json_path)
    analysis = _analyze_results(result)
    _plot_results(result)
    _print_footer(decay_ratio=analysis["decay_ratio"])


if __name__ == "__main__":
    main()
