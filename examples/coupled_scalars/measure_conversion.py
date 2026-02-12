"""Measurement example: wave conversion in coupled Klein-Gordon fields.

Demonstrates the measurement module applied to the coupled scalar system:

    L = 1/2(d phi)^2 - 1/2 m_phi^2 phi^2
      + 1/2(d chi)^2 - 1/2 m_chi^2 chi^2 - g phi chi

A Gaussian pulse excites the lighter field (phi, m^2=1).  Coupling (g=0.5)
transfers energy to the heavier field (chi, m^2=4).  The measurement module
quantifies this conversion using canonical Hamiltonian energy, spectral
decomposition, and energy conservation diagnostics.

Usage:
    uv run python examples/coupled_scalars/measure_conversion.py
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, MemoryStorage

from tidal.measurement import (
    SimulationData,
    check_energy_conservation,
    compute_conversion_probability,
    compute_energy_timeseries,
    compute_mode_amplitudes,
    compute_spectrum,
)
from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.vectorfield import ComponentGaussianPulse

if TYPE_CHECKING:
    from pde.trackers.base import TrackerBase

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

PARAMS: dict[str, float] = {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5}
T_END = 30.0
TRACKER_INTERVAL = 0.2
OUTPUT_FILENAME = "coupled_scalars_measurement.png"


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def _run_simulation() -> tuple[SimulationData, dict[str, float]]:
    """Run the coupled KG simulation and return measurement-ready data."""
    json_path = Path(__file__).parent.parent / "data" / "coupled_scalars.json"

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid([(0, 100)], 256, periodic=True)

    pulse = ComponentGaussianPulse(
        center=(30.0,),
        width=5.0,
        amplitude=1.0,
        active_components={"phi_0": 1.0},
    )
    initial_state = pulse.create(grid, spec)

    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=TRACKER_INTERVAL)
    pde.solve(
        initial_state, t_range=T_END, dt=0.01,
        scheme="runge-kutta", tracker=tracker,
    )

    data = SimulationData.from_storage(storage, spec, grid, PARAMS)
    return data, PARAMS


# ---------------------------------------------------------------------------
# Measurement + printing
# ---------------------------------------------------------------------------


def _print_summary(
    result: ConversionResult,
    diag: EnergyDiagnostics,
    params: dict[str, float],
) -> None:
    """Print quantitative summary to stdout."""
    print("=" * 60)
    print("Coupled Klein-Gordon: Measurement Summary")
    print("=" * 60)
    print()

    # Parameters
    print(f"  m_phi^2 = {params['mPhi2']},  m_chi^2 = {params['mChi2']},  g = {params['gCpl']}")

    # Normal mode frequencies (uniform-mode, k=0)
    m_phi2, m_chi2, g = params["mPhi2"], params["mChi2"], params["gCpl"]
    avg = (m_phi2 + m_chi2) / 2.0
    delta = np.sqrt(((m_phi2 - m_chi2) / 2.0) ** 2 + g**2)
    omega_plus = np.sqrt(avg + delta)
    omega_minus = np.sqrt(avg - delta)
    print(f"  Normal modes (k=0): omega_+ = {omega_plus:.4f},  omega_- = {omega_minus:.4f}")
    print(f"  Beat frequency: delta_omega = {omega_plus - omega_minus:.4f}")
    print(f"  Rabi period: T_Rabi = {2 * np.pi / (omega_plus - omega_minus):.2f}")
    print()

    # Conversion
    peak_idx = int(np.argmax(result.probability))
    print(f"  Peak conversion P(t) = {result.probability[peak_idx]:.6f}")
    print(f"    at t = {result.times[peak_idx]:.2f}")
    print(f"  Initial source energy E_phi(0) = {result.source_energy[0]:.4f}")
    print(f"  Final  source energy E_phi(T) = {result.source_energy[-1]:.4f}")
    print(f"  Final  target energy E_chi(T) = {result.target_energy[-1]:.4f}")
    print()

    # Conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print("    threshold  = 1e-3")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _plot_results(  # noqa: PLR0914, PLR0915
    data: SimulationData,
    result: ConversionResult,
    diag: EnergyDiagnostics,
    params: dict[str, float],
) -> Path:
    """Generate 2x3 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        "Coupled Klein-Gordon: Measurement Analysis\n"
        r"$\mathcal{L} = \frac{1}{2}(\partial\varphi)^2"
        r" - \frac{1}{2}m_\varphi^2\varphi^2"
        r" + \frac{1}{2}(\partial\chi)^2"
        r" - \frac{1}{2}m_\chi^2\chi^2 - g\,\varphi\chi$",
        fontsize=13,
    )

    # [0,0] Conversion probability P(t)
    ax = axes[0, 0]
    ax.plot(result.times, result.probability, "b-", linewidth=1.5)
    peak_idx = int(np.argmax(result.probability))
    ax.plot(
        result.times[peak_idx], result.probability[peak_idx],
        "ro", markersize=6,
    )
    ax.annotate(
        f"P = {result.probability[peak_idx]:.4f}\nt = {result.times[peak_idx]:.1f}",
        xy=(result.times[peak_idx], result.probability[peak_idx]),
        xytext=(15, -10), textcoords="offset points", fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "gray"},
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t) = E_\chi(t)\,/\,E_\varphi(0)$")
    ax.set_title("Conversion Probability")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Per-field energy timeseries
    ax = axes[0, 1]
    times, per_field, interaction, total = compute_energy_timeseries(data)
    ax.plot(times, per_field["phi_0"], "b-", label=r"$E_\varphi$", linewidth=1.2)
    ax.plot(times, per_field["chi_0"], "r-", label=r"$E_\chi$", linewidth=1.2)
    ax.plot(times, interaction, "g--", label=r"$E_\mathrm{int}$", linewidth=1.0, alpha=0.7)
    ax.plot(times, total, "k-", label=r"$E_\mathrm{total}$", linewidth=1.0, alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy")
    ax.set_title("Energy Decomposition")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)

    # [0,2] Energy conservation (relative error)
    ax = axes[0, 2]
    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(1e-3, color="r", linestyle="--", alpha=0.5, label="threshold (1e-3)")
    ax.axhline(-1e-3, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E\,/\,E_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,0] Spectral power at t=0 and t=peak
    ax = axes[1, 0]
    snap0 = compute_spectrum(data.fields["phi_0"][0], data.grid_spacing, data.periodic)
    snap_peak = compute_spectrum(
        data.fields["chi_0"][peak_idx], data.grid_spacing, data.periodic,
    )
    ax.semilogy(snap0.wavenumbers, snap0.power_spectrum, "b-", label=r"$\varphi(t=0)$")
    ax.semilogy(
        snap_peak.wavenumbers, snap_peak.power_spectrum,
        "r-", label=rf"$\chi(t={result.times[peak_idx]:.1f})$",
    )
    ax.set_xlabel(r"$|k|$")
    ax.set_ylabel(r"$|\hat{\phi}(k)|^2$")
    ax.set_title("Power Spectrum")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,1] Mode amplitude waterfall for chi
    ax = axes[1, 1]
    mode_times, wavenumbers, amplitudes = compute_mode_amplitudes(data, "chi_0")
    if len(wavenumbers) > 1 and len(mode_times) > 1:
        n_modes = min(20, len(wavenumbers))
        extent = [float(mode_times[0]), float(mode_times[-1]),
                  float(wavenumbers[0]), float(wavenumbers[min(n_modes, len(wavenumbers) - 1)])]
        im = ax.imshow(
            amplitudes[:, :n_modes].T,
            aspect="auto", origin="lower", extent=extent,
            cmap="inferno",
        )
        plt.colorbar(im, ax=ax, label=r"$|\hat{\chi}(k)|$")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$|k|$")
    ax.set_title(r"$\chi$ Mode Amplitudes")

    # [1,2] Summary text
    ax = axes[1, 2]
    m_phi2, m_chi2, g = params["mPhi2"], params["mChi2"], params["gCpl"]
    avg = (m_phi2 + m_chi2) / 2.0
    delta = np.sqrt(((m_phi2 - m_chi2) / 2.0) ** 2 + g**2)
    omega_plus = np.sqrt(avg + delta)
    omega_minus = np.sqrt(avg - delta)
    lines = [
        "Parameters:",
        f"  $m_\\varphi^2 = {m_phi2}$,  $m_\\chi^2 = {m_chi2}$,  $g = {g}$",
        "",
        "Normal Modes (k=0):",
        f"  $\\omega_+ = {omega_plus:.4f}$,  $\\omega_- = {omega_minus:.4f}$",
        f"  $T_{{Rabi}} = {2 * np.pi / (omega_plus - omega_minus):.2f}$",
        "",
        "Measurement Results:",
        f"  Peak $P(t) = {result.probability[peak_idx]:.6f}$",
        f"  at $t = {result.times[peak_idx]:.2f}$",
        f"  max $|\\Delta E / E_0| = {diag.max_relative_error:.2e}$",
        f"  Conservation: {'PASS' if diag.is_conserved else 'FAIL'}",
        "",
        "Initial condition:",
        "  Gaussian in $\\varphi$, $\\chi = 0$",
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
    print("Running coupled Klein-Gordon simulation...")
    data, params = _run_simulation()
    print(f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]")

    print("Computing conversion probability...")
    result = compute_conversion_probability(data, "phi_0", "chi_0")

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=1e-3)

    _print_summary(result, diag, params)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(data, result, diag, params)
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
