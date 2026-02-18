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
scalar-to-vector conversion P(t) = E_{A_1+A_2}(t) / E_phi(0), plus:

- **Spectral conversion** P(k,t) via ``compute_group_spectral_conversion``
- **Dispersion relation** omega(k) via ``compute_dispersion``
- **Mixing length** and **mixing spectrum**

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
from pde import CallbackTracker, CartesianGrid, FieldCollection, ScalarField

from tidal.measurement import (
    SimulationData,
    check_energy_conservation,
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
PARAMS: dict[str, float] = {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2}

# Grid
GRID_BOUNDS = [(0, 50), (0, 50)]
GRID_SHAPE = [96, 96]
GRID_PERIODIC = True

# Time integration
T_END = 200.0
DT = 0.01
TRACKER_INTERVAL = 0.2

# Initial conditions
PULSE_CENTER_X = 25.0
PULSE_CENTER_Y = 25.0
PULSE_WIDTH = 5.0

# Output
OUTPUT_FILENAME = "scalar_vector_measurement.png"

# ── Simulation ────────────────────────────────────────────────


def _run_simulation() -> tuple[SimulationData, Path]:
    """Run the scalar-vector simulation and return measurement-ready data."""
    json_path = Path(__file__).parent.parent / "data" / "scalar_vector_coupling.json"

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid(GRID_BOUNDS, GRID_SHAPE, periodic=GRID_PERIODIC)

    # Build 7-slot initial state: phi, pi_phi, A_0, A_1, pi_A1, A_2, pi_A2
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    gaussian = np.exp(
        -((x - PULSE_CENTER_X) ** 2 + (y - PULSE_CENTER_Y) ** 2)
        / (2 * PULSE_WIDTH**2)
    )

    state = FieldCollection(
        [
            ScalarField(grid, data=gaussian, label="phi_0"),
            ScalarField(grid, data=0.0, label="pi_phi"),
            ScalarField(grid, data=0.0, label="A_0"),
            ScalarField(grid, data=0.0, label="A_1"),
            ScalarField(grid, data=0.0, label="pi_A1"),
            ScalarField(grid, data=0.0, label="A_2"),
            ScalarField(grid, data=0.0, label="pi_A2"),
        ]
    )

    output_data_dir = (
        Path(__file__).parent.parent / "data" / "scalar_vector_coupling_output"
    )
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


def _print_summary(  # noqa: PLR0915
    total: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
) -> None:
    """Print quantitative summary to stdout."""
    print("=" * 60)
    print("Scalar-Vector Coupling: Measurement Summary")
    print("=" * 60)
    print()

    # Parameters
    print(
        f"  phim2 = {PARAMS['phim2']},  Am2 = {PARAMS['Am2']},  "
        f"kCS = {PARAMS['kCS']},  gSV = {PARAMS['gSV']}"
    )
    print()

    # Group conversion
    peak_idx = int(np.argmax(total.probability))
    print("  Total scalar->vector P(t):")
    print(f"    Peak P = {total.probability[peak_idx]:.6f}")
    print(f"    at t = {total.times[peak_idx]:.2f}")
    print(f"    Initial source energy E_phi(0) = {total.source_energy[0]:.4f}")
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
    if disp_phi is not None:
        n_active = int(np.count_nonzero(disp_phi.peak_frequencies > 0.0))
        print("  Dispersion (phi):")
        print(f"    Active k-modes: {n_active} / {len(disp_phi.wavenumbers)}")
        print(f"    Rayleigh resolution: {disp_phi.rayleigh_resolution:.4f} rad/time")
    else:
        print("  Dispersion (phi): not computed")
    print()

    # Conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print("    threshold  = 1e-3")
    print("=" * 60)


# ── Plotting ──────────────────────────────────────────────────


def _plot_results(  # noqa: C901, PLR0912, PLR0914, PLR0915
    data: SimulationData,
    total: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
) -> Path:
    """Generate 2x4 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    fig.suptitle(
        "Scalar-Vector Coupling: Measurement Analysis\n"
        r"$\mathcal{L} \supset g_{\mathrm{SV}}\,\varphi\,\nabla\!\cdot\!A"
        r" + \frac{k_{\mathrm{CS}}}{2}\epsilon^{abc}A_a\partial_b A_c$",
        fontsize=13,
    )

    peak_idx = int(np.argmax(total.probability))

    # [0,0] Total scalar->vector P(t) + mixing length annotation
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
    ax.set_ylabel(r"$P(t) = E_{A_1+A_2}(t)\,/\,E_\varphi(0)$")
    ax.set_title(r"Total $\varphi \to \mathbf{A}$ Conversion")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Energy conservation
    ax = axes[0, 1]
    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(1e-3, color="r", linestyle="--", alpha=0.5, label="threshold (1e-3)")
    ax.axhline(-1e-3, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E\,/\,E_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [0,2] Mixing spectrum (temporal FFT of P(t))
    ax = axes[0, 2]
    if spectrum is not None:
        ax.semilogy(
            spectrum.frequencies, spectrum.power, "b-", linewidth=0.5, alpha=0.7
        )
        ax.axvline(
            spectrum.dominant_frequency,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=rf"$\omega_{{\mathrm{{dom}}}}$ = {spectrum.dominant_frequency:.2f}",
        )
        x_max = min(10 * spectrum.dominant_frequency, spectrum.frequencies[-1])
        ax.set_xlim(0, x_max)
        ax.set_xlabel(r"$\omega$ (rad/time)")
        ax.set_ylabel(r"$|\hat{P}(\omega)|^2$")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "Not computed", transform=ax.transAxes, ha="center")
    ax.set_title("Mixing Spectrum")
    ax.grid(visible=True, alpha=0.3)

    # [0,3] Per-field energy timeseries
    ax = axes[0, 3]
    times, per_field, interaction, energy_total = compute_energy_timeseries(data)
    ax.plot(times, per_field["phi_0"], "b-", label=r"$E_\varphi$", linewidth=1.2)
    if "A_1" in per_field:
        ax.plot(times, per_field["A_1"], "r-", label=r"$E_{A_1}$", linewidth=1.2)
    if "A_2" in per_field:
        ax.plot(times, per_field["A_2"], "g-", label=r"$E_{A_2}$", linewidth=1.2)
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
    ax.set_ylabel("Energy density")
    ax.set_title("Energy Decomposition")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)

    # [1,0] Spectral conversion P(k,t) heatmap
    ax = axes[1, 0]
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

    # [1,1] Dispersion relation S(k, omega) heatmap
    ax = axes[1, 1]
    if disp_phi is not None and np.any(disp_phi.peak_frequencies > 0.0):
        log_power = np.log10(np.maximum(disp_phi.power, 1e-30))
        log_max = float(log_power.max())
        mesh = ax.pcolormesh(
            disp_phi.wavenumbers,
            disp_phi.frequencies,
            log_power.T,
            shading="nearest",
            cmap="viridis",
            vmin=log_max - 20,
            vmax=log_max,
        )
        fig.colorbar(mesh, ax=ax, label=r"$\log_{10} S(k, \omega)$", pad=0.02)
        active = disp_phi.peak_frequencies > 0.0
        ax.plot(
            disp_phi.wavenumbers[active],
            disp_phi.peak_frequencies[active],
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
            max_peak = float(np.max(disp_phi.peak_powers))
            if max_peak > 0:
                strong = disp_phi.peak_powers >= max_peak * 1e-3
                if np.any(strong):
                    ax.set_xlim(0, float(disp_phi.wavenumbers[strong].max()) * 1.1)
        # omega-axis: power threshold (limits to actual oscillation frequencies)
        peak_pwr = float(np.maximum(np.max(disp_phi.power), 1e-30))
        sig_f = np.max(disp_phi.power, axis=0) >= peak_pwr * 1e-6
        if np.any(sig_f):
            ax.set_ylim(0, float(disp_phi.frequencies[sig_f].max()) * 1.1)
        ax.set_xlabel(r"$|k|$")
        ax.set_ylabel(r"$\omega$ (rad/time)")
        ax.set_title(r"Dispersion $\omega(k)$ ($\varphi$)")
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
        "  Gaussian in $\\varphi$, $A = 0$",
    ]
    ax.text(
        0.05,
        0.95,
        "\n".join(lines),
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        fontfamily="monospace",
    )
    ax.axis("off")

    # [1,3] blank
    axes[1, 3].axis("off")

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
    print("Running scalar-vector coupling simulation...")
    data, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(
        f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]"
    )

    print("Computing group conversion probability (phi -> {A_1, A_2})...")
    total = compute_group_conversion(data, "phi_0")

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
            "phi_0",
            ["A_1", "A_2"],
        )
    except ValueError as e:
        print(f"  Spectral conversion: not computed ({e})")

    print("Computing dispersion relation...")
    disp_phi: DispersionResult | None = None
    try:
        disp_phi = compute_dispersion(data, "phi_0")
    except ValueError as e:
        print(f"  Dispersion (phi): not computed ({e})")

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=1e-3)

    _print_summary(total, diag, mixing, spectrum, spectral_conv, disp_phi)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        total,
        diag,
        mixing,
        spectrum,
        spectral_conv,
        disp_phi,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
