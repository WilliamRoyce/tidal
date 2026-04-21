"""Measurement visualization for ``tidal measure --output``.

Generates a 2x3 measurement figure:

.. code-block::

    [0,0] Conversion P(t)     [0,1] Per-field energy     [0,2] Energy conservation
    [1,0] Power spectrum       [1,1] Mixing spectrum      [1,2] Summary text

Each panel is guarded by ``if key in results`` — gracefully skips
missing measurements.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pathlib import Path

    from matplotlib.axes import Axes

    from tidal.measurement._io import SimulationData

# ---------------------------------------------------------------------------
# Panel renderers
# ---------------------------------------------------------------------------

_ENERGY_FLOOR = 1e-30


def _plot_conversion(ax: Axes, results: dict[str, Any]) -> None:
    """[0,0] Conversion probability P(t)."""
    import numpy as np

    if "conversion" not in results or "error" in results["conversion"]:
        ax.text(0.5, 0.5, "No conversion data", transform=ax.transAxes, ha="center")
        return

    conv = results["conversion"]
    conv_result = conv.get("_result_obj")
    if conv_result is None:
        return

    ax.plot(conv_result.times, conv_result.probability, "b-", linewidth=1.5)
    peak_idx = int(np.argmax(conv_result.probability))
    ax.plot(
        conv_result.times[peak_idx],
        conv_result.probability[peak_idx],
        "ro",
        markersize=6,
    )
    ax.set_ylim(bottom=0)
    src, tgt = ", ".join(conv["source"]), ", ".join(conv["target"])
    ax.set_ylabel(rf"$P(t) = E_{{\mathrm{{{tgt}}}}} / E_{{\mathrm{{{src}}}}}(0)$")

    if "mixing" in results and "error" not in results["mixing"]:
        mix = results["mixing"]
        ax.annotate(
            f"$L_{{mix}}$ = {mix['mixing_length']:.2f}"
            f" $\\pm$ {mix['mixing_length_uncertainty']:.2f}",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            ha="right",
            va="top",
            fontsize=9,
            color="green",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "alpha": 0.8},
        )


def _plot_energy(ax: Axes, results: dict[str, Any]) -> None:
    """[0,1] Per-field energy timeseries."""
    if "energy" not in results or "error" in results["energy"]:
        ax.text(0.5, 0.5, "No energy data", transform=ax.transAxes, ha="center")
        return

    eng = results["energy"]
    times = eng["times"]
    for name, series in eng["per_field"].items():
        ax.plot(times, series, label=name, linewidth=1.2)
    ax.plot(
        times, eng["interaction"], "--", label="interaction", linewidth=1.0, alpha=0.7,
    )
    ax.plot(times, eng["total"], "k-", label="total", linewidth=1.0, alpha=0.5)
    ax.legend(fontsize=7, ncol=2)


def _plot_conservation(ax: Axes, results: dict[str, Any]) -> None:
    """[0,2] Energy conservation (relative error)."""
    import numpy as np

    if "energy" not in results or "error" in results["energy"]:
        ax.text(0.5, 0.5, "No conservation data", transform=ax.transAxes, ha="center")
        return

    eng = results["energy"]
    total = np.array(eng["total"])
    e0 = total[0]
    if abs(e0) > _ENERGY_FLOOR:
        rel_err = (total - e0) / e0
        ax.plot(eng["times"], rel_err, "k-", linewidth=1.0)

    if "conservation" in results and "error" not in results["conservation"]:
        threshold = results["conservation"]["threshold"]
        ax.axhline(
            threshold,
            color="r",
            linestyle="--",
            alpha=0.5,
            label=f"threshold ({threshold:.0e})",
        )
        ax.axhline(-threshold, color="r", linestyle="--", alpha=0.5)
        ax.legend(fontsize=8)


def _plot_spectrum(ax: Axes, results: dict[str, Any]) -> None:
    """[1,0] Spatial power spectrum (initial + final)."""
    if "spectrum" not in results or "error" in results["spectrum"]:
        ax.set_visible(False)
        return

    for name, snap in results["spectrum"].items():
        ax.semilogy(
            snap["initial"]["wavenumbers"],
            snap["initial"]["power"],
            label=f"{name}(t=0)",
            alpha=0.7,
        )
        ax.semilogy(
            snap["final"]["wavenumbers"],
            snap["final"]["power"],
            "--",
            label=f"{name}(t=final)",
            alpha=0.7,
        )
    ax.legend(fontsize=7)


def _plot_dispersion(ax: Axes, results: dict[str, Any]) -> None:
    """[1,0] Dispersion relation S(k, omega) heatmap with peak curve overlay."""
    import numpy as np

    if "dispersion" not in results or "error" in results["dispersion"]:
        ax.set_visible(False)
        return

    disp = results["dispersion"]
    result = disp.get("_result_obj")
    if result is None:
        ax.set_visible(False)
        return

    wn = result.wavenumbers
    freqs = result.frequencies
    power = result.power

    if power.size == 0 or len(freqs) == 0:
        ax.text(0.5, 0.5, "No frequency bins", transform=ax.transAxes, ha="center")
        return

    # Log-scale power with dynamic range clipping
    log_power = np.log10(np.maximum(power, 1e-30))
    log_max = float(log_power.max())

    mesh = ax.pcolormesh(
        wn,
        freqs,
        log_power.T,
        shading="nearest",
        cmap="viridis",
        vmin=log_max - 20,
        vmax=log_max,
    )
    ax.figure.colorbar(mesh, ax=ax, label=r"$\log_{10} S(k, \omega)$", pad=0.02)  # pyright: ignore[reportOptionalMemberAccess]

    # Overlay peak dispersion curve (only active modes)
    active = result.peak_frequencies > 0.0
    if np.any(active):
        ax.plot(
            wn[active],
            result.peak_frequencies[active],
            "w--",
            linewidth=1.5,
            alpha=0.9,
            label=r"$\omega(k)$ peak",
        )
        ax.legend(fontsize=8, loc="upper left")

    # Data-driven axis cropping
    if np.any(active):
        max_peak = float(np.max(result.peak_powers))
        if max_peak > 0:
            strong = result.peak_powers >= max_peak * 1e-3
            if np.any(strong):
                ax.set_xlim(0, float(wn[strong].max()) * 1.1)
    # ω-axis: power threshold (limits to actual oscillation frequencies)
    peak_pwr = float(np.maximum(np.max(power), 1e-30))
    sig_f = np.max(power, axis=0) >= peak_pwr * 1e-6
    if np.any(sig_f):
        ax.set_ylim(0, float(freqs[sig_f].max()) * 1.1)

    ax.set_xlabel(r"$|k|$")
    ax.set_ylabel(r"$\omega$ (rad/time)")
    ax.set_title(f"Dispersion ({disp['field']})", fontsize=10)


def _plot_mixing_spectrum(ax: Axes, results: dict[str, Any]) -> None:
    """[1,1] Mixing spectrum (temporal FFT of P(t))."""
    # Prefer cached spectrum from _run_mixing; fall back to recomputation
    spectrum = None
    if "mixing" in results and "error" not in results["mixing"]:
        spectrum = results["mixing"].get("_spectrum_obj")

    if spectrum is None:
        # Recompute from conversion result as fallback
        if "conversion" not in results or "error" in results["conversion"]:
            ax.text(0.5, 0.5, "No mixing data", transform=ax.transAxes, ha="center")
            return

        conv_result = results["conversion"].get("_result_obj")
        if conv_result is None:
            return

        from tidal.measurement import compute_mixing_spectrum

        try:
            spectrum = compute_mixing_spectrum(conv_result)
        except ValueError:
            ax.text(0.5, 0.5, "Not computed", transform=ax.transAxes, ha="center")
            return

    ax.semilogy(spectrum.frequencies, spectrum.power, "b-", linewidth=0.5, alpha=0.7)
    ax.axvline(
        spectrum.dominant_frequency,
        color="red",
        linestyle="--",
        alpha=0.7,
        label=rf"$\omega_{{\mathrm{{dom}}}}$ = {spectrum.dominant_frequency:.2f}",
    )
    x_max = min(10 * spectrum.dominant_frequency, spectrum.frequencies[-1])
    ax.set_xlim(0, x_max)
    ax.legend(fontsize=8)


def _plot_summary_text(ax: Axes, data: SimulationData, results: dict[str, Any]) -> None:
    """[1,2] Summary text panel."""
    lines: list[str] = []

    if data.parameters:
        lines.append("Parameters:")
        for k, v in data.parameters.items():
            lines.append(f"  {k} = {v}")
        lines.append("")

    if "conservation" in results and "error" not in results["conservation"]:
        cons = results["conservation"]
        lines.extend(
            [
                f"Conservation: {'PASS' if cons['is_conserved'] else 'FAIL'}",
                f"  max |dE/E| = {cons['max_relative_error']:.2e}",
                "",
            ],
        )

    if "conversion" in results and "error" not in results["conversion"]:
        conv = results["conversion"]
        lines.extend(
            [
                f"Peak P(t) = {conv['peak_probability']:.6f}",
                f"  at t = {conv['peak_time']:.2f}",
                "",
            ],
        )

    if "mixing" in results and "error" not in results["mixing"]:
        mix = results["mixing"]
        lines.extend(
            [
                f"L_mix = {mix['mixing_length']:.4f}",
                f"  +/- {mix['mixing_length_uncertainty']:.4f}",
                f"omega_dom = {mix['dominant_frequency']:.4f}",
                "",
            ],
        )

    ax.text(
        0.05,
        0.95,
        "\n".join(lines),
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        fontfamily="monospace",
    )
    ax.axis("off")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def save_measurement_plot(
    path: Path,
    data: SimulationData,
    results: dict[str, Any],
) -> None:
    """Generate and save a 2x3 measurement figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("TIDAL Measurement Analysis", fontsize=13)

    # Row 0
    _plot_conversion(axes[0, 0], results)
    axes[0, 0].set_xlabel("Time")
    axes[0, 0].set_title("Conversion Probability")
    axes[0, 0].grid(visible=True, alpha=0.3)

    _plot_energy(axes[0, 1], results)
    axes[0, 1].set_xlabel("Time")
    axes[0, 1].set_ylabel("Energy")
    axes[0, 1].set_title("Energy Decomposition")
    axes[0, 1].grid(visible=True, alpha=0.3)

    _plot_conservation(axes[0, 2], results)
    axes[0, 2].set_xlabel("Time")
    axes[0, 2].set_ylabel(r"$\Delta E / E_0$")
    axes[0, 2].set_title("Energy Conservation")
    axes[0, 2].grid(visible=True, alpha=0.3)

    # Row 1 — priority: dispersion > spectrum
    if "dispersion" in results and "error" not in results.get("dispersion", {}):
        _plot_dispersion(axes[1, 0], results)
    else:
        _plot_spectrum(axes[1, 0], results)
        axes[1, 0].set_xlabel(r"$|k|$")
        axes[1, 0].set_ylabel(r"$|\hat{\phi}(k)|^2$")
        axes[1, 0].set_title("Power Spectrum")
    axes[1, 0].grid(visible=True, alpha=0.3)

    _plot_mixing_spectrum(axes[1, 1], results)
    axes[1, 1].set_xlabel(r"$\omega$ (rad/time)")
    axes[1, 1].set_ylabel(r"$|\hat{P}(\omega)|^2$")
    axes[1, 1].set_title("Mixing Spectrum")
    axes[1, 1].grid(visible=True, alpha=0.3)

    _plot_summary_text(axes[1, 2], data, results)

    plt.tight_layout()
    fig.savefig(str(path), dpi=150, bbox_inches="tight")
    plt.close(fig)
