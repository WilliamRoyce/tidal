"""Measurement example: scalar-to-vector conversion in a CS+divergence system.

Demonstrates the measurement module applied to the scalar-vector coupling:

    L = -1/2 (d phi)^2 - phim2/2 phi^2
        -1/4 F_ab F^ab - Am2/2 A^2
        + kCS/2 eps^abc A_a d_b A_c
        + gSV phi div(A)

A Gaussian pulse excites the scalar field phi.  The divergence coupling
(gSV) transfers energy to the vector components A_1, A_2.  The Chern-Simons
term (kCS) creates helical mixing between A components.

This example uses ``compute_group_conversion`` to measure the total
scalar-to-vector conversion P(t) = E_{A_1+A_2}(t) / E_phi(0), which is
the natural group-to-group measurement for mixed-rank field systems.

Usage:
    uv run python examples/scalar_vector_coupling/measure_conversion.py
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from tidal.measurement import (
    SimulationData,
    check_energy_conservation,
    compute_conversion_probability,
    compute_energy_timeseries,
    compute_group_conversion,
    compute_spectrum,
)
from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.utils import normalize_solve_result

if TYPE_CHECKING:
    from pde.trackers.base import TrackerBase

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

PARAMS: dict[str, float] = {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2}
T_END = 15.0
TRACKER_INTERVAL = 0.3
OUTPUT_FILENAME = "scalar_vector_measurement.png"


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def _run_simulation() -> SimulationData:
    """Run the scalar-vector simulation and return measurement-ready data."""
    json_path = Path(__file__).parent.parent / "data" / "scalar_vector_coupling.json"

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid([(0, 50), (0, 50)], [64, 64], periodic=True)

    # Build 7-slot initial state: phi, pi_phi, A_0, A_1, pi_A1, A_2, pi_A2
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    gaussian = np.exp(-((x - 25.0) ** 2 + (y - 25.0) ** 2) / (2 * 5.0**2))

    state = FieldCollection([
        ScalarField(grid, data=gaussian, label="phi_0"),
        ScalarField(grid, data=0.0, label="pi_phi"),
        ScalarField(grid, data=0.0, label="A_0"),
        ScalarField(grid, data=0.0, label="A_1"),
        ScalarField(grid, data=0.0, label="pi_A1"),
        ScalarField(grid, data=0.0, label="A_2"),
        ScalarField(grid, data=0.0, label="pi_A2"),
    ])

    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=TRACKER_INTERVAL)
    result = pde.solve(
        state, t_range=T_END, dt=0.01,
        scheme="runge-kutta", tracker=tracker,
    )
    normalize_solve_result(result)

    return SimulationData.from_storage(storage, spec, grid, PARAMS)


# ---------------------------------------------------------------------------
# Measurement + printing
# ---------------------------------------------------------------------------


def _print_summary(
    total: ConversionResult,
    r_a1: ConversionResult,
    r_a2: ConversionResult,
    diag: EnergyDiagnostics,
) -> None:
    """Print quantitative summary to stdout."""
    print("=" * 60)
    print("Scalar-Vector Coupling: Measurement Summary")
    print("=" * 60)
    print()

    # Parameters
    print(f"  phim2 = {PARAMS['phim2']},  Am2 = {PARAMS['Am2']},  "
          f"kCS = {PARAMS['kCS']},  gSV = {PARAMS['gSV']}")
    print()

    # Group conversion
    peak_idx = int(np.argmax(total.probability))
    print("  Total scalar->vector P(t):")
    print(f"    Peak P = {total.probability[peak_idx]:.6f}")
    print(f"    at t = {total.times[peak_idx]:.2f}")
    print(f"    Initial source energy E_phi(0) = {total.source_energy[0]:.4f}")
    print()

    # Per-component
    peak_a1 = int(np.argmax(r_a1.probability))
    peak_a2 = int(np.argmax(r_a2.probability))
    print("  Per-component breakdown:")
    print(f"    P(phi->A_1): peak = {r_a1.probability[peak_a1]:.6f} at t = {r_a1.times[peak_a1]:.2f}")
    print(f"    P(phi->A_2): peak = {r_a2.probability[peak_a2]:.6f} at t = {r_a2.times[peak_a2]:.2f}")
    print()

    # Conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print("    threshold  = 1e-2")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _plot_results(  # noqa: PLR0915
    data: SimulationData,
    total: ConversionResult,
    r_a1: ConversionResult,
    r_a2: ConversionResult,
    diag: EnergyDiagnostics,
) -> Path:
    """Generate 2x3 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        "Scalar-Vector Coupling: Measurement Analysis\n"
        r"$\mathcal{L} \supset g_{\mathrm{SV}}\,\varphi\,\nabla\!\cdot\!A"
        r" + \frac{k_{\mathrm{CS}}}{2}\epsilon^{abc}A_a\partial_b A_c$",
        fontsize=13,
    )

    peak_idx = int(np.argmax(total.probability))

    # [0,0] Total scalar->vector P(t)
    ax = axes[0, 0]
    ax.plot(total.times, total.probability, "b-", linewidth=1.5)
    ax.plot(
        total.times[peak_idx], total.probability[peak_idx],
        "ro", markersize=6,
    )
    ax.annotate(
        f"P = {total.probability[peak_idx]:.4f}\nt = {total.times[peak_idx]:.1f}",
        xy=(total.times[peak_idx], total.probability[peak_idx]),
        xytext=(15, -10), textcoords="offset points", fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "gray"},
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t) = E_{A_1+A_2}(t)\,/\,E_\varphi(0)$")
    ax.set_title(r"Total $\varphi \to \mathbf{A}$ Conversion")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Per-component P(t): A_1 and A_2
    ax = axes[0, 1]
    ax.plot(r_a1.times, r_a1.probability, "r-", label=r"$P(\varphi \to A_1)$", linewidth=1.2)
    ax.plot(r_a2.times, r_a2.probability, "g-", label=r"$P(\varphi \to A_2)$", linewidth=1.2)
    ax.plot(total.times, total.probability, "b--", label="Total", linewidth=1.0, alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t)$")
    ax.set_title("Per-Component Conversion")
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,2] Energy conservation
    ax = axes[0, 2]
    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(1e-2, color="r", linestyle="--", alpha=0.5, label="threshold (1e-2)")
    ax.axhline(-1e-2, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E\,/\,E_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,0] Per-field energy timeseries
    ax = axes[1, 0]
    times, per_field, interaction, energy_total = compute_energy_timeseries(data)
    ax.plot(times, per_field["phi_0"], "b-", label=r"$E_\varphi$", linewidth=1.2)
    if "A_1" in per_field:
        ax.plot(times, per_field["A_1"], "r-", label=r"$E_{A_1}$", linewidth=1.2)
    if "A_2" in per_field:
        ax.plot(times, per_field["A_2"], "g-", label=r"$E_{A_2}$", linewidth=1.2)
    ax.plot(times, interaction, "m--", label=r"$E_\mathrm{int}$", linewidth=1.0, alpha=0.7)
    ax.plot(times, energy_total, "k-", label=r"$E_{\mathrm{total}}$", linewidth=1.0, alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy")
    ax.set_title("Energy Decomposition")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)

    # [1,1] 2D spectral power: phi at t=0 vs A_1 at t=peak
    ax = axes[1, 1]
    snap_phi = compute_spectrum(data.fields["phi_0"][0], data.grid_spacing, data.periodic)
    snap_a1 = compute_spectrum(
        data.fields["A_1"][peak_idx], data.grid_spacing, data.periodic,
    )
    ax.semilogy(snap_phi.wavenumbers, snap_phi.power_spectrum, "b-", label=r"$\varphi(t=0)$")
    ax.semilogy(
        snap_a1.wavenumbers, snap_a1.power_spectrum,
        "r-", label=rf"$A_1(t={total.times[peak_idx]:.1f})$",
    )
    ax.set_xlabel(r"$|k|$")
    ax.set_ylabel(r"$|\hat{\phi}(k)|^2$")
    ax.set_title("Power Spectrum")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,2] Summary text
    ax = axes[1, 2]
    lines = [
        "Parameters:",
        f"  $m_\\varphi^2 = {PARAMS['phim2']}$,  $m_A^2 = {PARAMS['Am2']}$",
        f"  $k_{{CS}} = {PARAMS['kCS']}$,  $g_{{SV}} = {PARAMS['gSV']}$",
        "",
        "Measurement Results:",
        f"  Total $P(\\varphi \\to A) = {total.probability[peak_idx]:.6f}$",
        f"  at $t = {total.times[peak_idx]:.2f}$",
        "",
        "Per-Component Peaks:",
        f"  $P(\\varphi \\to A_1) = {r_a1.probability.max():.6f}$",
        f"  $P(\\varphi \\to A_2) = {r_a2.probability.max():.6f}$",
        "",
        f"  max $|\\Delta E / E_0| = {diag.max_relative_error:.2e}$",
        f"  Conservation: {'PASS' if diag.is_conserved else 'FAIL'}",
        "",
        "Initial condition:",
        "  Gaussian in $\\varphi$, $A = 0$",
    ]
    ax.text(
        0.05, 0.95, "\n".join(lines),
        transform=ax.transAxes, fontsize=10,
        verticalalignment="top", fontfamily="monospace",
    )
    ax.axis("off")

    plt.tight_layout()
    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Run simulation and perform full measurement analysis."""
    print("Running scalar-vector coupling simulation...")
    data = _run_simulation()
    print(f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]")

    print("Computing group conversion probability (phi -> {A_1, A_2})...")
    total = compute_group_conversion(data, "phi_0")

    print("Computing per-component conversion...")
    r_a1 = compute_conversion_probability(data, "phi_0", "A_1")
    r_a2 = compute_conversion_probability(data, "phi_0", "A_2")

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=1e-2)

    _print_summary(total, r_a1, r_a2, diag)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(data, total, r_a1, r_a2, diag)
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
