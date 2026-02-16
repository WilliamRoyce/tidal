"""Measurement example: A->B energy transfer in coupled Proca system.

Demonstrates the measurement module applied to the coupled Proca system:

    L = -1/4 F^A F^A - 1/4 F^B F^B
        - mA2/2 A^2 - mB2/2 B^2 + gcoup A.B

A Gaussian pulse in A_1 transfers energy to the B sector via the gcoup
coupling.  The A_0 and B_0 constraint fields are solved via coupled FFT
(periodic boundary conditions).

This example demonstrates:
- **Group conversion** P(t) via ``compute_group_conversion``
- **Spectral conversion** P(k,t) via ``compute_group_spectral_conversion``
- **Dispersion relation** omega(k) via ``compute_dispersion``
- Per-component conversion breakdown

Usage:
    uv run python examples/coupled_proca/measure_conversion.py
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
    compute_dispersion,
    compute_energy_timeseries,
    compute_group_conversion,
    compute_group_spectral_conversion,
    compute_mixing_length,
    compute_mixing_spectrum,
)
from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.utils import normalize_solve_result

if TYPE_CHECKING:
    from pde.trackers.base import TrackerBase

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics
    from tidal.measurement._dispersion import DispersionResult
    from tidal.measurement._spectral_conversion import SpectralConversion

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

PARAMS: dict[str, float] = {"mA2": 1.0, "mB2": 5.0, "gcoup": 0.5}
T_END = 100.0
TRACKER_INTERVAL = 0.16
ENERGY_THRESHOLD = 0.002  # Depends strongly on resolution and BCs; periodic BCs give machine-precision conservation, while Dirichlet can have O(1%) errors.
OUTPUT_FILENAME = "coupled_proca_measurement.png"


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def _run_simulation() -> SimulationData:
    """Run the coupled Proca simulation and return measurement data."""
    json_path = Path(__file__).parent.parent / "data" / "coupled_proca_3d.json"

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid(
        bounds=[(0, np.pi), (0, np.pi)],
        shape=[96, 96],
        periodic=True,
    )

    # Build initial state from spec layout — Gaussian pulse in A_1
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    gaussian = 0.5 * np.exp(
        -((x - np.pi / 2) ** 2 + (y - np.pi / 2) ** 2) / (2 * 0.5**2)
    )

    fields: list[ScalarField] = []
    for name, slot_type in spec.state_layout:
        sf = ScalarField(grid, data=0.0, label=f"{name}_{slot_type}")
        if name == "A_1" and slot_type == "field":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)

    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=TRACKER_INTERVAL)
    result = pde.solve(
        state,
        t_range=T_END,
        dt=0.01,
        scheme="runge-kutta",
        tracker=tracker,
    )
    normalize_solve_result(result)

    return SimulationData.from_storage(storage, spec, grid, PARAMS)


# ---------------------------------------------------------------------------
# Measurement + printing
# ---------------------------------------------------------------------------


def _print_summary(  # noqa: PLR0913, PLR0915, PLR0917
    total: ConversionResult,
    r_a2: ConversionResult,
    r_b1: ConversionResult,
    r_b2: ConversionResult,
    diag: EnergyDiagnostics,
    spectral_conv: SpectralConversion | None,
    disp_a1: DispersionResult | None,
) -> None:
    """Print quantitative summary to stdout."""
    print("=" * 60)
    print("Coupled Proca: Measurement Summary")
    print("=" * 60)
    print()

    print(
        f"  mA2 = {PARAMS['mA2']},  mB2 = {PARAMS['mB2']},  gcoup = {PARAMS['gcoup']}"
    )
    print()

    # Group conversion
    peak_idx = int(np.argmax(total.probability))
    print("  Total A_1 -> {A_2, B_1, B_2} conversion P(t):")
    print(f"    Peak P = {total.probability[peak_idx]:.6f}")
    print(f"    at t = {total.times[peak_idx]:.2f}")
    print(f"    Initial source energy E_A1(0) = {total.source_energy[0]:.4f}")
    print()

    # Per-component
    peak_a2 = int(np.argmax(r_a2.probability))
    peak_b1 = int(np.argmax(r_b1.probability))
    peak_b2 = int(np.argmax(r_b2.probability))
    print("  Per-component breakdown:")
    print(
        f"    P(A_1->A_2): peak = {r_a2.probability[peak_a2]:.6f} at t = {r_a2.times[peak_a2]:.2f}"
    )
    print(
        f"    P(A_1->B_1): peak = {r_b1.probability[peak_b1]:.6f} at t = {r_b1.times[peak_b1]:.2f}"
    )
    print(
        f"    P(A_1->B_2): peak = {r_b2.probability[peak_b2]:.6f} at t = {r_b2.times[peak_b2]:.2f}"
    )
    print()

    # Mixing length
    try:
        mixing = compute_mixing_length(total)
        print("  Mixing length (spectral):")
        print(
            f"    L_mix      = {mixing.mixing_length:.4f} +/- {mixing.mixing_length_uncertainty:.4f}"
        )
        print(
            f"    omega_dom  = {mixing.dominant_frequency:.4f}  (FWHM = {mixing.frequency_fwhm:.4f})"
        )
        print(f"    max P(t)   = {mixing.max_conversion:.6f}")
        if len(mixing.peaks) > 1:
            print(f"    ({len(mixing.peaks)} spectral peaks detected)")
    except ValueError as e:
        print(f"  Mixing length: not extracted ({e})")
    print()

    # Mixing spectrum
    try:
        spectrum = compute_mixing_spectrum(total)
        print("  Mixing spectrum (temporal FFT of P(t)):")
        print(
            f"    dominant oscillation freq: omega = {spectrum.dominant_frequency:.4f}"
        )
        print(f"    dominant spectral L_mix:   {spectrum.dominant_mixing_length:.4f}")
        print(f"    frequency bins: {len(spectrum.frequencies)}")
    except ValueError as e:
        print(f"  Mixing spectrum: not computed ({e})")
    print()

    # Spectral conversion
    if spectral_conv is not None:
        n_active = int(spectral_conv.active_modes.sum())
        print("  Spectral conversion P(k,t):")
        print(
            f"    Active k-modes: {n_active} / {len(spectral_conv.wavenumbers)}"
        )
        if n_active > 0:
            peak_k_idx = int(np.argmax(spectral_conv.probability[-1]))
            print(
                f"    Peak P(k, t_final) at |k| = {spectral_conv.wavenumbers[peak_k_idx]:.4f}:"
                f"  P = {spectral_conv.probability[-1, peak_k_idx]:.6f}"
            )
    else:
        print("  Spectral conversion: not computed")
    print()

    # Dispersion
    if disp_a1 is not None:
        n_active = int(np.count_nonzero(disp_a1.peak_frequencies > 0.0))
        print("  Dispersion (A_1):")
        print(
            f"    Active k-modes: {n_active} / {len(disp_a1.wavenumbers)}"
        )
        print(
            f"    Rayleigh resolution: {disp_a1.rayleigh_resolution:.4f} rad/time"
        )
    else:
        print("  Dispersion (A_1): not computed")
    print()

    # Conservation — periodic BCs give machine-precision conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print(f"    threshold  = {ENERGY_THRESHOLD}")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _plot_results(  # noqa: PLR0913, PLR0914, PLR0915, PLR0917
    data: SimulationData,
    total: ConversionResult,
    r_a2: ConversionResult,
    r_b1: ConversionResult,
    r_b2: ConversionResult,
    diag: EnergyDiagnostics,
    spectral_conv: SpectralConversion | None,
    disp_a1: DispersionResult | None,
) -> Path:
    """Generate 2x3 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    fig.suptitle(
        "Coupled Proca 2+1D (Periodic): Measurement Analysis\n"
        r"$\mathcal{L} \supset -\frac{m_A^2}{2}A^2 - \frac{m_B^2}{2}B^2"
        r" + g\,A\!\cdot\!B$"
        f"\nmA2={PARAMS['mA2']}, mB2={PARAMS['mB2']}, gcoup={PARAMS['gcoup']}",
        fontsize=12,
    )

    peak_idx = int(np.argmax(total.probability))

    # [0,0] Total P(t) — group conversion from A_1 + mixing length
    ax = axes[0, 0]
    ax.plot(total.times, total.probability, "b-", linewidth=1.5)
    ax.plot(
        total.times[peak_idx],
        total.probability[peak_idx],
        "ro",
        markersize=6,
    )
    ax.annotate(
        f"P = {total.probability[peak_idx]:.4f}\nt = {total.times[peak_idx]:.1f}",
        xy=(total.times[peak_idx], total.probability[peak_idx]),
        xytext=(15, -10),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "gray"},
    )
    # Mixing length annotation
    try:
        mixing = compute_mixing_length(total)
        ax.annotate(
            f"$L_{{mix}}$ = {mixing.mixing_length:.2f} $\\pm$ {mixing.mixing_length_uncertainty:.2f}",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            ha="right",
            va="top",
            fontsize=9,
            color="green",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "alpha": 0.8},
        )
    except ValueError:
        pass
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t) = E_{\mathrm{target}}(t)\,/\,E_{A_1}(0)$")
    ax.set_title(r"Total $A_1 \to \{A_2, B_1, B_2\}$ Conversion")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Per-component P(t)
    ax = axes[0, 1]
    ax.plot(
        r_a2.times, r_a2.probability, "c-", label=r"$P(A_1 \to A_2)$", linewidth=1.2
    )
    ax.plot(
        r_b1.times, r_b1.probability, "r-", label=r"$P(A_1 \to B_1)$", linewidth=1.2
    )
    ax.plot(
        r_b2.times, r_b2.probability, "g-", label=r"$P(A_1 \to B_2)$", linewidth=1.2
    )
    ax.plot(
        total.times, total.probability, "b--", label="Total", linewidth=1.0, alpha=0.5
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t)$")
    ax.set_title("Per-Component Conversion")
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,2] Mixing spectrum (temporal FFT of P(t))
    ax = axes[0, 2]
    try:
        spectrum = compute_mixing_spectrum(total)
        ax.semilogy(
            spectrum.frequencies, spectrum.power, "b-", linewidth=0.5, alpha=0.7
        )
        ax.axvline(
            spectrum.dominant_frequency,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=rf"$\omega_{{dom}}$ = {spectrum.dominant_frequency:.2f}",
        )
        # Focus on the interesting frequency range (up to 10x dominant)
        x_max = min(10 * spectrum.dominant_frequency, spectrum.frequencies[-1])
        ax.set_xlim(0, x_max)
        ax.set_xlabel(r"Angular frequency $\omega$ (rad/time)")
        ax.set_ylabel(r"Power $|\hat{P}(\omega)|^2$")
        ax.set_title("Mixing Spectrum")
        ax.legend(fontsize=8)
    except ValueError as e:
        ax.text(
            0.5,
            0.5,
            f"Not computed:\n{e}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
        )
        ax.set_title("Mixing Spectrum")
    ax.grid(visible=True, alpha=0.3)

    # [1,0] Energy decomposition
    ax = axes[1, 0]
    times, per_field, interaction, energy_total = compute_energy_timeseries(data)
    for name, series in per_field.items():
        ax.plot(times, series, linewidth=1.2, label=rf"$E_{{{name}}}$")
    ax.plot(
        times, interaction, "m--", label=r"$E_\mathrm{int}$", linewidth=1.0, alpha=0.7
    )
    ax.plot(
        times,
        energy_total,
        "k-",
        label=r"$E_{\mathrm{total}}$",
        linewidth=1.0,
        alpha=0.5,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy")
    ax.set_title("Energy Decomposition")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(visible=True, alpha=0.3)

    # [1,1] Energy conservation
    ax = axes[1, 1]
    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(
        ENERGY_THRESHOLD,
        color="r",
        linestyle="--",
        alpha=0.5,
        label=f"threshold ({ENERGY_THRESHOLD})",
    )
    ax.axhline(-ENERGY_THRESHOLD, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E\,/\,E_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,2] Summary text panel
    ax = axes[1, 2]
    ax.axis("off")
    try:
        mixing = compute_mixing_length(total)
        summary_text = (
            "Mixing Length Summary\n"
            "---------------------\n"
            f"L_mix     = {mixing.mixing_length:.4f}\n"
            f"  +/-     = {mixing.mixing_length_uncertainty:.4f}\n"
            f"omega_dom = {mixing.dominant_frequency:.4f}\n"
            f"FWHM      = {mixing.frequency_fwhm:.4f}\n"
            f"max P     = {mixing.max_conversion:.6f}\n"
            f"peaks     = {len(mixing.peaks)}\n"
        )
    except ValueError:
        summary_text = "No mixing detected\n"
    if spectral_conv is not None:
        n_active = int(spectral_conv.active_modes.sum())
        summary_text += f"\nSpectral Conv: {n_active} active modes"
    if disp_a1 is not None:
        n_act = int(np.count_nonzero(disp_a1.peak_frequencies > 0.0))
        summary_text += f"\nDispersion: {n_act} active modes"
    ax.text(
        0.1,
        0.9,
        summary_text,
        transform=ax.transAxes,
        fontsize=10,
        family="monospace",
        verticalalignment="top",
    )

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
    print("Running coupled Proca simulation (periodic BCs)...")
    data = _run_simulation()
    print(
        f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]"
    )

    print("Computing group conversion probability (A_1 -> {A_2, B_1, B_2})...")
    total = compute_group_conversion(data, "A_1")

    print("Computing per-component conversion...")
    r_a2 = compute_conversion_probability(data, "A_1", "A_2")
    r_b1 = compute_conversion_probability(data, "A_1", "B_1")
    r_b2 = compute_conversion_probability(data, "A_1", "B_2")

    print("Computing spectral conversion P(k,t)...")
    spectral_conv: SpectralConversion | None = None
    try:
        spectral_conv = compute_group_spectral_conversion(
            data, ["A_1"], ["B_1", "B_2"],
        )
    except ValueError as e:
        print(f"  Spectral conversion: not computed ({e})")

    print("Computing dispersion relation...")
    disp_a1: DispersionResult | None = None
    try:
        disp_a1 = compute_dispersion(data, "A_1")
    except ValueError as e:
        print(f"  Dispersion (A_1): not computed ({e})")

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=ENERGY_THRESHOLD)

    _print_summary(total, r_a2, r_b1, r_b2, diag, spectral_conv, disp_a1)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data, total, r_a2, r_b1, r_b2, diag, spectral_conv, disp_a1,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
