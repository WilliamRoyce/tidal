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
from pde import CallbackTracker, CartesianGrid, FieldCollection, ScalarField

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
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.utils import normalize_solve_result

if TYPE_CHECKING:
    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics
    from tidal.measurement._dispersion import DispersionResult
    from tidal.measurement._mixing import MixingResult, MixingSpectrum
    from tidal.measurement._spectral_conversion import SpectralConversion

# ── Configuration ─────────────────────────────────────────────

# Physics
PARAMS: dict[str, float] = {"mA2": 1.0, "mB2": 5.0, "gcoup": 0.5}

# Grid
GRID_BOUNDS = [(0, 50), (0, 50)]
GRID_SHAPE = [64, 64]
GRID_PERIODIC = True

# Time integration
T_END = 10.0
DT = 0.01
TRACKER_INTERVAL = 0.2

# Initial conditions
PULSE_CENTER_X = 25.0
PULSE_CENTER_Y = 25.0
PULSE_WIDTH = 5.0

# Analysis
ENERGY_THRESHOLD = 0.002  # periodic BCs give machine-precision conservation

# Output
OUTPUT_FILENAME = "coupled_proca_measurement.png"

# ── Simulation ────────────────────────────────────────────────


def _run_simulation() -> tuple[SimulationData, Path]:
    """Run the coupled Proca simulation and return measurement data."""
    json_path = Path(__file__).parent.parent / "data" / "coupled_proca_3d.json"

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid(
        bounds=GRID_BOUNDS,
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )

    # Build initial state from spec layout — Gaussian pulse in A_1
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    gaussian = np.exp(
        -((x - PULSE_CENTER_X) ** 2 + (y - PULSE_CENTER_Y) ** 2)
        / (2 * PULSE_WIDTH**2)
    )

    fields: list[ScalarField] = []
    for name, slot_type in spec.state_layout:
        sf = ScalarField(grid, data=0.0, label=f"{name}_{slot_type}")
        if name == "A_1" and slot_type == "field":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)

    output_data_dir = Path(__file__).parent.parent / "data" / "coupled_proca_output"
    writer, callback = create_snapshot_callback(
        output_dir=output_data_dir,
        spec=spec,
        grid=grid,
        t_end=T_END,
        snapshot_interval=TRACKER_INTERVAL,
        parameters=PARAMS,
        spec_path=json_path,
    )
    tracker = CallbackTracker(callback, interrupts=TRACKER_INTERVAL)
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=tracker,
    )
    normalize_solve_result(result)
    writer.close()

    return SimulationData.from_directory(output_data_dir, spec), output_data_dir


# ── Measurement + printing ────────────────────────────────────


def _print_summary(  # noqa: PLR0913, PLR0915, PLR0917
    total: ConversionResult,
    r_a2: ConversionResult,
    r_b1: ConversionResult,
    r_b2: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
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
    if mixing is not None:
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
    else:
        print("  Mixing length: not extracted (no spectral peaks)")
    print()

    # Mixing spectrum
    if spectrum is not None:
        print("  Mixing spectrum (temporal FFT of P(t)):")
        print(
            f"    dominant oscillation freq: omega = {spectrum.dominant_frequency:.4f}"
        )
        print(f"    dominant spectral L_mix:   {spectrum.dominant_mixing_length:.4f}")
        print(f"    frequency bins: {len(spectrum.frequencies)}")
    else:
        print("  Mixing spectrum: not computed")
    print()

    # Spectral conversion
    if spectral_conv is not None:
        n_active = int(spectral_conv.active_modes.sum())
        print("  Spectral conversion P(k,t):")
        print(f"    Active k-modes: {n_active} / {len(spectral_conv.wavenumbers)}")
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
        print(f"    Active k-modes: {n_active} / {len(disp_a1.wavenumbers)}")
        print(f"    Rayleigh resolution: {disp_a1.rayleigh_resolution:.4f} rad/time")
    else:
        print("  Dispersion (A_1): not computed")
    print()

    # Conservation — periodic BCs give machine-precision conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print(f"    threshold  = {ENERGY_THRESHOLD}")
    print("=" * 60)


# ── Plotting ──────────────────────────────────────────────────


def _plot_results(  # noqa: C901, PLR0912, PLR0913, PLR0914, PLR0915, PLR0917
    data: SimulationData,
    total: ConversionResult,
    r_a2: ConversionResult,
    r_b1: ConversionResult,
    r_b2: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_a1: DispersionResult | None,
) -> Path:
    """Generate 2x4 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
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
    ax.plot(total.times[peak_idx], total.probability[peak_idx], "ro", markersize=6)
    ax.annotate(
        f"P = {total.probability[peak_idx]:.4f}\nt = {total.times[peak_idx]:.1f}",
        xy=(total.times[peak_idx], total.probability[peak_idx]),
        xytext=(15, -10),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "gray"},
    )
    if mixing is not None:
        ax.annotate(
            f"$L_{{mix}}$ = {mixing.mixing_length:.2f}"
            f" $\\pm$ {mixing.mixing_length_uncertainty:.2f}",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            ha="right",
            va="top",
            fontsize=9,
            color="green",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "alpha": 0.8},
        )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t) = E_{\mathrm{target}}(t)\,/\,E_{A_1}(0)$")
    ax.set_title(r"Total $A_1 \to \{A_2, B_1, B_2\}$ Conversion")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Per-component P(t)
    ax = axes[0, 1]
    ax.plot(
        r_a2.times,
        r_a2.probability,
        "c-",
        label=r"$P(A_1 \to A_2)$",
        linewidth=1.2,
    )
    ax.plot(
        r_b1.times,
        r_b1.probability,
        "r-",
        label=r"$P(A_1 \to B_1)$",
        linewidth=1.2,
    )
    ax.plot(
        r_b2.times,
        r_b2.probability,
        "g-",
        label=r"$P(A_1 \to B_2)$",
        linewidth=1.2,
    )
    ax.plot(
        total.times,
        total.probability,
        "b--",
        label="Total",
        linewidth=1.0,
        alpha=0.5,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t)$")
    ax.set_title("Per-Component Conversion")
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,2] Energy conservation
    ax = axes[0, 2]
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

    # [0,3] Mixing spectrum (temporal FFT of P(t))
    ax = axes[0, 3]
    if spectrum is not None:
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
        x_max = min(10 * spectrum.dominant_frequency, spectrum.frequencies[-1])
        ax.set_xlim(0, x_max)
        ax.set_xlabel(r"$\omega$ (rad/time)")
        ax.set_ylabel(r"$|\hat{P}(\omega)|^2$")
        ax.legend(fontsize=8)
    else:
        ax.text(
            0.5,
            0.5,
            "Not computed",
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
        times,
        interaction,
        "m--",
        label=r"$E_\mathrm{int}$",
        linewidth=1.0,
        alpha=0.7,
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

    # [1,1] Spectral conversion P(k,t) heatmap
    ax = axes[1, 1]
    if spectral_conv is not None and np.any(spectral_conv.active_modes):
        mesh = ax.pcolormesh(
            spectral_conv.wavenumbers,
            spectral_conv.times,
            spectral_conv.probability,
            shading="nearest",
            cmap="inferno",
        )
        fig.colorbar(mesh, ax=ax, label=r"$P(k,t)$", pad=0.02)
        # Crop x-axis to active k-range
        k_active_max = float(
            spectral_conv.wavenumbers[spectral_conv.active_modes].max()
        )
        ax.set_xlim(0, k_active_max * 1.15)
        ax.set_xlabel(r"$|k|$")
        ax.set_ylabel("Time")
        ax.set_title(r"Spectral Conversion $P(k,t)$")
    else:
        ax.text(
            0.5,
            0.5,
            "No spectral\nconversion data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
        )
        ax.set_title("Spectral Conversion")
    ax.grid(visible=True, alpha=0.3)

    # [1,2] Dispersion relation S(k, omega) heatmap
    ax = axes[1, 2]
    if disp_a1 is not None and np.any(disp_a1.peak_frequencies > 0.0):
        log_power = np.log10(np.maximum(disp_a1.power, 1e-30))
        log_max = float(log_power.max())
        mesh = ax.pcolormesh(
            disp_a1.wavenumbers,
            disp_a1.frequencies,
            log_power.T,
            shading="nearest",
            cmap="viridis",
            vmin=log_max - 20,
            vmax=log_max,
        )
        fig.colorbar(mesh, ax=ax, label=r"$\log_{10} S(k, \omega)$", pad=0.02)
        active = disp_a1.peak_frequencies > 0.0
        ax.plot(
            disp_a1.wavenumbers[active],
            disp_a1.peak_frequencies[active],
            "w--",
            linewidth=1.5,
            alpha=0.9,
            label=r"$\omega(k)$ peak",
        )
        # Data-driven axis cropping
        # k-axis: match spectral conversion's active_modes range when available
        if spectral_conv is not None and np.any(spectral_conv.active_modes):
            ax.set_xlim(
                0,
                float(spectral_conv.wavenumbers[spectral_conv.active_modes].max())
                * 1.15,
            )
        elif np.any(active):
            max_peak = float(np.max(disp_a1.peak_powers))
            if max_peak > 0:
                strong = disp_a1.peak_powers >= max_peak * 1e-3
                if np.any(strong):
                    ax.set_xlim(0, float(disp_a1.wavenumbers[strong].max()) * 1.1)
        # ω-axis: power threshold (limits to actual oscillation frequencies)
        peak_pwr = float(np.maximum(np.max(disp_a1.power), 1e-30))
        sig_f = np.max(disp_a1.power, axis=0) >= peak_pwr * 1e-6
        if np.any(sig_f):
            ax.set_ylim(0, float(disp_a1.frequencies[sig_f].max()) * 1.1)
        ax.set_xlabel(r"$|k|$")
        ax.set_ylabel(r"$\omega$ (rad/time)")
        ax.set_title(r"Dispersion $\omega(k)$ ($A_1$)")
        ax.legend(fontsize=8, loc="upper left")
    else:
        ax.text(
            0.5,
            0.5,
            "No dispersion data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
        )
        ax.set_title("Dispersion")
    ax.grid(visible=True, alpha=0.3)

    # [1,3] Summary text panel
    ax = axes[1, 3]
    ax.axis("off")
    lines = [
        "Parameters:",
        f"  $m_A^2 = {PARAMS['mA2']}$,  $m_B^2 = {PARAMS['mB2']}$",
        f"  $g_{{coup}} = {PARAMS['gcoup']}$",
        "",
        "Results:",
        f"  Peak $P(t) = {total.probability[peak_idx]:.6f}$",
        f"  at $t = {total.times[peak_idx]:.2f}$",
        "",
        "Per-Component Peaks:",
        f"  $P(A_1 \\to A_2) = {r_a2.probability.max():.6f}$",
        f"  $P(A_1 \\to B_1) = {r_b1.probability.max():.6f}$",
        f"  $P(A_1 \\to B_2) = {r_b2.probability.max():.6f}$",
    ]
    if mixing is not None:
        lines += [
            "",
            "Mixing Length:",
            f"  $L_{{mix}} = {mixing.mixing_length:.4f}$",
            f"  $\\omega_{{dom}} = {mixing.dominant_frequency:.4f}$",
        ]
    lines += [
        "",
        f"  max $|\\Delta E / E_0| = {diag.max_relative_error:.2e}$",
        f"  Conservation: {'PASS' if diag.is_conserved else 'FAIL'}",
        "",
        "Initial condition:",
        "  Gaussian in $A_1$, all others $= 0$",
    ]
    ax.text(
        0.05,
        0.9,
        "\n".join(lines),
        transform=ax.transAxes,
        fontsize=9,
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


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run simulation and perform full measurement analysis."""
    print("Running coupled Proca simulation (periodic BCs)...")
    data, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(
        f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]"
    )

    print("Computing group conversion probability (A_1 -> {A_2, B_1, B_2})...")
    total = compute_group_conversion(data, "A_1")

    print("Computing per-component conversion...")
    r_a2 = compute_conversion_probability(data, "A_1", "A_2")
    r_b1 = compute_conversion_probability(data, "A_1", "B_1")
    r_b2 = compute_conversion_probability(data, "A_1", "B_2")

    print("Computing mixing length and spectrum...")
    mixing: MixingResult | None = None
    spectrum: MixingSpectrum | None = None
    try:
        mixing = compute_mixing_length(total)
    except ValueError as e:
        print(f"  Mixing length: not extracted ({e})")
    try:
        spectrum = compute_mixing_spectrum(total)
    except ValueError as e:
        print(f"  Mixing spectrum: not computed ({e})")

    print("Computing spectral conversion P(k,t)...")
    spectral_conv: SpectralConversion | None = None
    try:
        spectral_conv = compute_group_spectral_conversion(
            data,
            ["A_1"],
            ["B_1", "B_2"],
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

    _print_summary(
        total,
        r_a2,
        r_b1,
        r_b2,
        diag,
        mixing,
        spectrum,
        spectral_conv,
        disp_a1,
    )

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        total,
        r_a2,
        r_b1,
        r_b2,
        diag,
        mixing,
        spectrum,
        spectral_conv,
        disp_a1,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
