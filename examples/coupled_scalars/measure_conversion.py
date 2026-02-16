"""Measurement example: wave conversion in coupled Klein-Gordon fields.

Demonstrates the measurement module applied to the coupled scalar system:

    L = 1/2(d phi)^2 - 1/2 m_phi^2 phi^2
      + 1/2(d chi)^2 - 1/2 m_chi^2 chi^2 - g phi chi

A Gaussian pulse excites the lighter field (phi, m^2=1).  Coupling (g=0.5)
transfers energy to the heavier field (chi, m^2=4).  The measurement module
quantifies this conversion using:

- **Conversion probability** P(t) = E_chi(t) / E_phi(0)
- **Mixing length** L_mix = pi / omega_dom (half-period of dominant oscillation)
- **Spectral conversion** P(k,t) — per-mode conversion probability
- **Dispersion relation** omega(k) — extracted via spacetime FFT
- **Energy conservation** diagnostics

Usage:
    uv run python examples/coupled_scalars/measure_conversion.py
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from pde import CallbackTracker, CartesianGrid

from tidal.measurement import (
    SimulationData,
    check_energy_conservation,
    compute_conversion_probability,
    compute_dispersion,
    compute_energy_timeseries,
    compute_mixing_length,
    compute_mixing_spectrum,
    compute_spectral_conversion,
    compute_spectrum,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.vectorfield import ComponentGaussianPulse

if TYPE_CHECKING:
    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics
    from tidal.measurement._dispersion import DispersionResult
    from tidal.measurement._mixing import MixingResult, MixingSpectrum
    from tidal.measurement._spectral_conversion import SpectralConversion

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

PARAMS: dict[str, float] = {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5}
T_END = 200.0
TRACKER_INTERVAL = 0.2
OUTPUT_FILENAME = "coupled_scalars_measurement.png"


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def _run_simulation() -> tuple[SimulationData, dict[str, float], Path]:
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

    output_data_dir = Path(__file__).parent.parent / "data" / "coupled_scalars_output"
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
    pde.solve(
        initial_state,
        t_range=T_END,
        dt=0.01,
        scheme="runge-kutta",
        tracker=tracker,
    )
    writer.close()

    data = SimulationData.from_directory(output_data_dir, spec)
    return data, PARAMS, output_data_dir


# ---------------------------------------------------------------------------
# Measurement + printing
# ---------------------------------------------------------------------------


def _print_summary(  # noqa: PLR0912, PLR0913, PLR0915, PLR0917
    result: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
    disp_chi: DispersionResult | None,
    params: dict[str, float],
) -> None:
    """Print quantitative summary to stdout."""
    print("=" * 60)
    print("Coupled Klein-Gordon: Measurement Summary")
    print("=" * 60)
    print()

    # Parameters
    print(
        f"  m_phi^2 = {params['mPhi2']},  m_chi^2 = {params['mChi2']},  g = {params['gCpl']}"
    )

    # Normal mode frequencies (uniform-mode, k=0)
    m_phi2, m_chi2, g = params["mPhi2"], params["mChi2"], params["gCpl"]
    avg = (m_phi2 + m_chi2) / 2.0
    delta = np.sqrt(((m_phi2 - m_chi2) / 2.0) ** 2 + g**2)
    omega_plus = np.sqrt(avg + delta)
    omega_minus = np.sqrt(avg - delta)
    print(
        f"  Normal modes (k=0): omega_+ = {omega_plus:.4f},  omega_- = {omega_minus:.4f}"
    )
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
        if n_active > 0:
            # Compare measured omega(k=0) with predicted
            omega_measured = disp_phi.peak_frequencies[0]
            print(f"    omega(k~0) measured: {omega_measured:.4f}")
    else:
        print("  Dispersion (phi): not computed")
    if disp_chi is not None:
        n_active = int(np.count_nonzero(disp_chi.peak_frequencies > 0.0))
        print("  Dispersion (chi):")
        print(f"    Active k-modes: {n_active} / {len(disp_chi.wavenumbers)}")
    else:
        print("  Dispersion (chi): not computed")
    print()

    # Conservation
    print(f"  Energy conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print(f"    max |dE/E| = {diag.max_relative_error:.2e}")
    print("    threshold  = 1e-3")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _plot_results(  # noqa: C901, PLR0912, PLR0913, PLR0914, PLR0915, PLR0917
    data: SimulationData,
    result: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
    params: dict[str, float],
) -> Path:
    """Generate 2x4 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    fig.suptitle(
        "Coupled Klein-Gordon: Measurement Analysis\n"
        r"$\mathcal{L} = \frac{1}{2}(\partial\varphi)^2"
        r" - \frac{1}{2}m_\varphi^2\varphi^2"
        r" + \frac{1}{2}(\partial\chi)^2"
        r" - \frac{1}{2}m_\chi^2\chi^2 - g\,\varphi\chi$",
        fontsize=13,
    )

    peak_idx = int(np.argmax(result.probability))

    # [0,0] Conversion probability P(t) + mixing length annotation
    ax = axes[0, 0]
    ax.plot(result.times, result.probability, "b-", linewidth=1.5)
    ax.plot(result.times[peak_idx], result.probability[peak_idx], "ro", markersize=6)
    ax.annotate(
        f"P = {result.probability[peak_idx]:.4f}\nt = {result.times[peak_idx]:.1f}",
        xy=(result.times[peak_idx], result.probability[peak_idx]),
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
    ax.set_ylabel(r"$P(t) = E_\chi(t)\,/\,E_\varphi(0)$")
    ax.set_title("Conversion Probability")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)

    # [0,1] Per-field energy timeseries
    ax = axes[0, 1]
    times, per_field, interaction, total = compute_energy_timeseries(data)
    ax.plot(times, per_field["phi_0"], "b-", label=r"$E_\varphi$", linewidth=1.2)
    ax.plot(times, per_field["chi_0"], "r-", label=r"$E_\chi$", linewidth=1.2)
    ax.plot(
        times,
        interaction,
        "g--",
        label=r"$E_\mathrm{int}$",
        linewidth=1.0,
        alpha=0.7,
    )
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
        # ω-axis: power threshold (limits to actual oscillation frequencies)
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

    # [1,2] Power spectrum (phi at t=0, chi at t=peak)
    ax = axes[1, 2]
    snap0 = compute_spectrum(data.fields["phi_0"][0], data.grid_spacing, data.periodic)
    snap_peak = compute_spectrum(
        data.fields["chi_0"][peak_idx],
        data.grid_spacing,
        data.periodic,
    )
    ax.semilogy(snap0.wavenumbers, snap0.power_spectrum, "b-", label=r"$\varphi(t=0)$")
    ax.semilogy(
        snap_peak.wavenumbers,
        snap_peak.power_spectrum,
        "r-",
        label=rf"$\chi(t={result.times[peak_idx]:.1f})$",
    )
    ax.set_xlabel(r"$|k|$")
    ax.set_ylabel(r"$|\hat{\phi}(k)|^2$")
    ax.set_title("Power Spectrum")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # [1,3] Summary text
    ax = axes[1, 3]
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
        "Results:",
        f"  Peak $P(t) = {result.probability[peak_idx]:.6f}$",
        f"  at $t = {result.times[peak_idx]:.2f}$",
    ]
    if mixing is not None:
        lines += [
            f"  $L_{{mix}} = {mixing.mixing_length:.4f}$",
            f"  $\\omega_{{dom}} = {mixing.dominant_frequency:.4f}$",
        ]
    lines += [
        "",
        f"  max $|\\Delta E / E_0| = {diag.max_relative_error:.2e}$",
        f"  Conservation: {'PASS' if diag.is_conserved else 'FAIL'}",
        "",
        "Initial condition:",
        "  Gaussian in $\\varphi$, $\\chi = 0$",
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
    data, params, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(
        f"  {data.n_snapshots} snapshots collected over t=[0, {float(data.times[-1]):.1f}]"
    )

    print("Computing conversion probability...")
    result = compute_conversion_probability(data, "phi_0", "chi_0")

    print("Computing mixing length and spectrum...")
    mixing: MixingResult | None = None
    spectrum: MixingSpectrum | None = None
    try:
        mixing = compute_mixing_length(result)
    except ValueError as e:
        print(f"  Mixing length: not extracted ({e})")
    try:
        spectrum = compute_mixing_spectrum(result)
    except ValueError as e:
        print(f"  Mixing spectrum: not computed ({e})")

    print("Computing spectral conversion P(k,t)...")
    spectral_conv: SpectralConversion | None = None
    try:
        spectral_conv = compute_spectral_conversion(data, "phi_0", "chi_0")
    except ValueError as e:
        print(f"  Spectral conversion: not computed ({e})")

    print("Computing dispersion relations...")
    disp_phi: DispersionResult | None = None
    disp_chi: DispersionResult | None = None
    try:
        disp_phi = compute_dispersion(data, "phi_0")
    except ValueError as e:
        print(f"  Dispersion (phi): not computed ({e})")
    try:
        disp_chi = compute_dispersion(data, "chi_0")
    except ValueError as e:
        print(f"  Dispersion (chi): not computed ({e})")

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=1e-3)

    _print_summary(
        result,
        diag,
        mixing,
        spectrum,
        spectral_conv,
        disp_phi,
        disp_chi,
        params,
    )

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        result,
        diag,
        mixing,
        spectrum,
        spectral_conv,
        disp_phi,
        params,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
