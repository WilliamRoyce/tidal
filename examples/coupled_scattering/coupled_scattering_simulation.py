"""Measurement example: wave conversion in coupled scalar scattering (2+1D).

Demonstrates the measurement module applied to a position-dependent coupling:

    L = -1/2 (d phi)^2 - 1/2 (d chi)^2
        - mPhi2/2 phi^2 - mChi2/2 chi^2 - G(x,y) phi chi

where G(x,y) = g0 * exp(-(x^2+y^2) / (2*R^2)) is a localized Gaussian
background field centered at the origin.

A rightward-propagating phi wave packet is incident on the coupling region.
Chi starts at zero — any chi that appears is purely from conversion.
The measurement module quantifies this conversion using:

- **Conversion probability** P(t) = E_chi(t) / E_phi(0)
- **Mixing length** L_mix = pi / omega_dom  (half-period of dominant oscillation)
- **Spectral conversion** P(k,t) — per-mode conversion probability
- **Dispersion relation** omega(k) — extracted via spacetime FFT
- **Hamiltonian energy** decomposition (manual, cross-validated against
  virial energy from the measurement module)

Analogy to Gertsenshtein effect:
  phi <-> photon,  chi <-> graviton,  G(x,y) <-> background B_0(x,y)

Usage:
    uv run python examples/coupled_scattering/coupled_scattering_simulation.py
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
    compute_conversion_probability,
    compute_dispersion,
    compute_mixing_length,
    compute_mixing_spectrum,
    compute_spectral_conversion,
    compute_system_energy,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._dispersion import DispersionResult
    from tidal.measurement._mixing import MixingResult, MixingSpectrum
    from tidal.measurement._spectral_conversion import SpectralConversion

# NOTE: The measurement module's virial energy (compute_system_energy) now
# supports position-dependent coupling G(x,y).  We keep the manual
# Hamiltonian decomposition as a cross-validation — it gives a physics-
# specific breakdown (phi energy, chi energy, coupling energy) that the
# generic virial formula does not.

# ── Configuration ─────────────────────────────────────────────

# Physics
MPHI2 = 1.0
MCHI2 = 4.0
G0 = 1.0
COUPLING_RADIUS = 8.0

PARAMS: dict[str, float] = {
    "mPhi2": MPHI2,
    "mChi2": MCHI2,
    "g0": G0,
    "R": COUPLING_RADIUS,
}

# Grid (2D, periodic)
DOMAIN = (-50.0, 50.0)
N_CELLS = 128

# Time integration
T_END = 1000.0
DT = 0.02
TRACKER_INTERVAL = 10.0

# Initial conditions — rightward-propagating phi wave packet
X0 = -25.0  # center of wave packet (left of coupling region)
SIGMA = 4.0  # Gaussian envelope width
K0 = 3.0  # wavevector (k0 > sqrt(mPhi2) for propagation)
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "coupled_scattering_measurement.png"


# ── Hamiltonian energy (manual) ───────────────────────────────


@dataclass(frozen=True)
class HamiltonianDecomposition:
    """Hamiltonian energy decomposition over time.

    H = sum_i [1/2 pi_i^2 + 1/2 |grad phi_i|^2 + m_i^2/2 phi_i^2]
        + integral G(x,y) phi chi dA
    """

    times: NDArray[np.float64]
    phi_energy: NDArray[np.float64]
    chi_energy: NDArray[np.float64]
    coupling_energy: NDArray[np.float64]
    total_energy: NDArray[np.float64]


def _gradient_energy_2d(
    field: NDArray[np.float64],
    dx: float,
    dy: float,
) -> float:
    """Compute 1/2 integral |grad phi|^2 dA using periodic central differences."""
    dv = dx * dy
    gx = (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2 * dx)
    gy = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2 * dy)
    return 0.5 * float(np.sum(gx**2 + gy**2)) * dv


def _compute_coupling_field(data: SimulationData) -> NDArray[np.float64]:
    """Compute G(x,y) = g0 * exp(-r^2 / (2*R^2)) on the simulation grid."""
    shape = data.fields["phi_0"].shape[1:]  # spatial (nx, ny)
    dx, dy = data.grid_spacing
    x_1d = np.linspace(
        data.grid_bounds[0][0] + dx / 2,
        data.grid_bounds[0][1] - dx / 2,
        shape[0],
    )
    y_1d = np.linspace(
        data.grid_bounds[1][0] + dy / 2,
        data.grid_bounds[1][1] - dy / 2,
        shape[1],
    )
    x, y = np.meshgrid(x_1d, y_1d, indexing="ij")
    g0 = data.parameters.get("g0", G0)
    r_param = data.parameters.get("R", COUPLING_RADIUS)
    return g0 * np.exp(-(x**2 + y**2) / (2 * r_param**2))


def _compute_hamiltonian_energies(
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> HamiltonianDecomposition:
    """Compute full Hamiltonian decomposition over all snapshots.

    Provides a physics-specific breakdown (phi energy, chi energy,
    coupling energy) as a cross-validation against the measurement
    module's virial energy.
    """
    dx, dy = data.grid_spacing
    dv = dx * dy
    n = data.n_snapshots
    phi_e = np.empty(n)
    chi_e = np.empty(n)
    cpl_e = np.empty(n)

    m_phi2 = data.parameters.get("mPhi2", MPHI2)
    m_chi2 = data.parameters.get("mChi2", MCHI2)

    for i in range(n):
        phi = data.fields["phi_0"][i]
        pi_phi = data.momenta["phi_0"][i]
        chi = data.fields["chi_0"][i]
        pi_chi = data.momenta["chi_0"][i]

        phi_e[i] = (
            0.5 * float(np.sum(pi_phi**2)) * dv
            + _gradient_energy_2d(phi, dx, dy)
            + 0.5 * m_phi2 * float(np.sum(phi**2)) * dv
        )
        chi_e[i] = (
            0.5 * float(np.sum(pi_chi**2)) * dv
            + _gradient_energy_2d(chi, dx, dy)
            + 0.5 * m_chi2 * float(np.sum(chi**2)) * dv
        )
        cpl_e[i] = float(np.sum(coupling_field * phi * chi)) * dv

    return HamiltonianDecomposition(
        times=data.times,
        phi_energy=phi_e,
        chi_energy=chi_e,
        coupling_energy=cpl_e,
        total_energy=phi_e + chi_e + cpl_e,
    )


# ── Simulation ────────────────────────────────────────────────


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    """Create ICs: traveling phi wave packet, chi = 0 (pure conversion)."""
    coords = cast("np.ndarray", grid.cell_coords)
    x = np.asarray(coords[..., 0], dtype=float)
    y = np.asarray(coords[..., 1], dtype=float)

    envelope = PULSE_AMPLITUDE * np.exp(
        -((x - X0) ** 2) / (2 * SIGMA**2) - y**2 / (2 * SIGMA**2)
    )
    phi_data = envelope * np.cos(K0 * (x - X0))
    omega = np.sqrt(K0**2 + MPHI2)
    pi_phi_data = omega * envelope * np.sin(K0 * (x - X0))

    return FieldCollection(
        [
            ScalarField(grid, data=phi_data, label="phi_0"),
            ScalarField(grid, data=pi_phi_data, label="pi_phi_0"),
            ScalarField(grid, data=np.zeros_like(x), label="chi_0"),
            ScalarField(grid, data=np.zeros_like(x), label="pi_chi_0"),
        ]
    )


def _run_simulation() -> tuple[SimulationData, dict[str, float], Path]:
    """Run the coupled scattering simulation and return measurement-ready data.

    Raises
    ------
    FileNotFoundError
        If the JSON spec has not been generated yet.
    """
    json_path = Path(__file__).parent.parent / "data" / "coupled_scattering.json"
    if not json_path.exists():
        msg = (
            f"{json_path} not found. "
            "Run 'tidal derive theory.toml' first (requires wolframscript)."
        )
        raise FileNotFoundError(msg)

    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters=PARAMS)
    grid = CartesianGrid([DOMAIN, DOMAIN], N_CELLS, periodic=True)
    initial_state = _create_initial_state(grid)

    output_data_dir = (
        Path(__file__).parent.parent / "data" / "coupled_scattering_output"
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
    pde.solve(
        initial_state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=tracker,
    )
    writer.close()

    data = SimulationData.from_directory(output_data_dir, spec)
    return data, PARAMS, output_data_dir


# ── Summary ───────────────────────────────────────────────────


def _print_summary(  # noqa: PLR0913, PLR0915, PLR0917
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
    params: dict[str, float],
) -> None:
    """Print quantitative measurement summary to stdout."""
    print("=" * 65)
    print("Coupled Scalar Scattering: Measurement Summary")
    print("=" * 65)
    print()

    m_phi2 = params["mPhi2"]
    m_chi2 = params["mChi2"]
    g0 = params["g0"]
    r_param = params["R"]
    print(f"  m_phi^2 = {m_phi2},  m_chi^2 = {m_chi2},  g0 = {g0},  R = {r_param}")

    omega_k = np.sqrt(K0**2 + m_phi2)
    v_group = K0 / omega_k
    print(f"  Wave packet: k0 = {K0},  omega = {omega_k:.4f},  v_group = {v_group:.4f}")
    print(f"  Stability: det(M) = {m_phi2 * m_chi2 - g0**2:.2f} > 0")
    print()

    peak_idx = int(np.argmax(result.probability))
    print(f"  Peak conversion P(t) = {result.probability[peak_idx]:.6f}")
    print(f"    at t = {result.times[peak_idx]:.2f}")
    print(f"  Initial source energy E_phi(0) = {result.source_energy[0]:.4f}")
    print(f"  Final  target energy E_chi(T) = {result.target_energy[-1]:.4f}")
    print()

    if mixing is not None:
        print("  Mixing length (spectral):")
        print(
            f"    L_mix      = {mixing.mixing_length:.4f}"
            f" +/- {mixing.mixing_length_uncertainty:.4f}"
        )
        print(
            f"    omega_dom  = {mixing.dominant_frequency:.4f}"
            f"  (FWHM = {mixing.frequency_fwhm:.4f})"
        )
    else:
        print("  Mixing length: not extracted (no spectral peaks)")
    print()

    if spectrum is not None:
        print("  Mixing spectrum (temporal FFT of P(t)):")
        print(f"    dominant frequency: {spectrum.dominant_frequency:.4f}")
        print(f"    dominant L_mix:     {spectrum.dominant_mixing_length:.4f}")
    else:
        print("  Mixing spectrum: not computed")
    print()

    if spectral_conv is not None:
        n_active = int(spectral_conv.active_modes.sum())
        print(f"  Spectral conversion: {n_active} active k-modes")
    else:
        print("  Spectral conversion: not computed")

    if disp_phi is not None:
        n_active = int(np.count_nonzero(disp_phi.peak_frequencies > 0.0))
        print(f"  Dispersion (phi): {n_active} active k-modes")
    else:
        print("  Dispersion: not computed")
    print()

    h0 = hamiltonian.total_energy[0]
    h_final = hamiltonian.total_energy[-1]
    max_drift = float(
        np.max(np.abs((hamiltonian.total_energy - h0) / max(abs(h0), 1e-30)))
    )
    print(f"  Hamiltonian H(0) = {h0:.4f},  H(T) = {h_final:.4f}")
    print(f"  max |dH/H| = {max_drift:.2e}")
    print("  (Manual decomposition cross-validated against virial energy)")
    print("=" * 65)


# ── Plotting ──────────────────────────────────────────────────


def _plot_results(  # noqa: PLR0913, PLR0917
    data: SimulationData,
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    spectral_conv: SpectralConversion | None,
    disp_phi: DispersionResult | None,
    coupling_field: NDArray[np.float64],
    params: dict[str, float],
) -> Path:
    """Generate 2x4 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    fig.suptitle(
        "Coupled Scalar Scattering (2+1D): Measurement Analysis\n"
        r"$\mathcal{L} = -\frac{1}{2}(\partial\phi)^2"
        r" - \frac{1}{2}(\partial\chi)^2"
        r" - \frac{m_\phi^2}{2}\phi^2"
        r" - \frac{m_\chi^2}{2}\chi^2"
        r" - G(x,y)\,\phi\,\chi$",
        fontsize=13,
    )

    peak_idx = int(np.argmax(result.probability))

    # [0,0] Conversion probability P(t)
    _plot_conversion(axes[0, 0], result, peak_idx, mixing)

    # [0,1] Hamiltonian energy decomposition
    _plot_hamiltonian(axes[0, 1], hamiltonian)

    # [0,2] Energy conservation |dH/H0|
    _plot_energy_conservation(axes[0, 2], hamiltonian)

    # [0,3] Mixing spectrum
    _plot_mixing_spectrum(axes[0, 3], spectrum)

    # [1,0] Spectral conversion P(k,t)
    _plot_spectral_conversion(axes[1, 0], fig, spectral_conv)

    # [1,1] Dispersion omega(k)
    _plot_dispersion(axes[1, 1], fig, disp_phi, spectral_conv)

    # [1,2] |chi| at final time + coupling contours
    _plot_chi_heatmap(axes[1, 2], data, coupling_field)

    # [1,3] Summary text
    _plot_summary_text(axes[1, 3], result, hamiltonian, mixing, params, peak_idx)

    plt.tight_layout()
    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _plot_conversion(
    ax: Axes,
    result: ConversionResult,
    peak_idx: int,
    mixing: MixingResult | None,
) -> None:
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
    ax.set_ylabel(r"$P(t) = E_\chi(t)\,/\,E_\phi(0)$")
    ax.set_title("Conversion Probability")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)


def _plot_hamiltonian(ax: Axes, h: HamiltonianDecomposition) -> None:
    ax.plot(h.times, h.phi_energy, "b-", label=r"$E_\phi$", linewidth=1.2)
    ax.plot(h.times, h.chi_energy, "r-", label=r"$E_\chi$", linewidth=1.2)
    ax.plot(
        h.times,
        h.coupling_energy,
        "g--",
        label=r"$\int G\phi\chi\,dA$",
        linewidth=1.0,
        alpha=0.7,
    )
    ax.plot(
        h.times,
        h.total_energy,
        "k-",
        label=r"$H$ (total)",
        linewidth=1.0,
        alpha=0.5,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy")
    ax.set_title("Hamiltonian Energy")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)


def _plot_energy_conservation(ax: Axes, h: HamiltonianDecomposition) -> None:
    h0 = h.total_energy[0]
    relative = (h.total_energy - h0) / max(abs(h0), 1e-30)
    ax.plot(h.times, relative, "k-", linewidth=1.0)
    ax.axhline(1e-3, color="r", linestyle="--", alpha=0.5, label="threshold (1e-3)")
    ax.axhline(-1e-3, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta H\,/\,H_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def _plot_mixing_spectrum(ax: Axes, spectrum: MixingSpectrum | None) -> None:
    if spectrum is not None:
        ax.semilogy(
            spectrum.frequencies,
            spectrum.power,
            "b-",
            linewidth=0.5,
            alpha=0.7,
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


def _plot_spectral_conversion(
    ax: Axes,
    fig: Figure,
    spectral_conv: SpectralConversion | None,
) -> None:
    if spectral_conv is not None and np.any(spectral_conv.active_modes):
        prob = spectral_conv.probability
        p_active = prob[:, spectral_conv.active_modes]
        vmax = float(np.percentile(p_active, 99)) if p_active.size > 0 else 1.0
        mesh = ax.pcolormesh(
            spectral_conv.wavenumbers,
            spectral_conv.times,
            prob,
            shading="nearest",
            cmap="inferno",
            vmax=max(vmax, 1e-6),
        )
        fig.colorbar(mesh, ax=ax, label=r"$P(k,t)$", pad=0.02)
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


def _plot_dispersion(
    ax: Axes,
    fig: Figure,
    disp_phi: DispersionResult | None,
    spectral_conv: SpectralConversion | None,
) -> None:
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
        fig.colorbar(mesh, ax=ax, label=r"$\log_{10} S(k,\omega)$", pad=0.02)
        active = disp_phi.peak_frequencies > 0.0
        ax.plot(
            disp_phi.wavenumbers[active],
            disp_phi.peak_frequencies[active],
            "w--",
            linewidth=1.5,
            alpha=0.9,
            label=r"$\omega(k)$ peak",
        )
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
        peak_pwr = float(np.maximum(np.max(disp_phi.power), 1e-30))
        sig_f = np.max(disp_phi.power, axis=0) >= peak_pwr * 1e-6
        if np.any(sig_f):
            ax.set_ylim(0, float(disp_phi.frequencies[sig_f].max()) * 1.1)
        ax.set_xlabel(r"$|k|$")
        ax.set_ylabel(r"$\omega$ (rad/time)")
        ax.set_title(r"Dispersion $\omega(k)$ ($\phi$)")
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


def _plot_chi_heatmap(
    ax: Axes,
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> None:
    """Plot |chi| at final time with coupling G(x,y) contours."""
    final_chi = data.fields["chi_0"][-1]
    bounds = data.grid_bounds
    extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    im = ax.imshow(
        np.abs(final_chi).T,
        aspect="equal",
        origin="lower",
        extent=extent,
        cmap="inferno",
    )
    # Coupling contours
    shape = final_chi.shape
    dx, dy = data.grid_spacing
    x_1d = np.linspace(bounds[0][0] + dx / 2, bounds[0][1] - dx / 2, shape[0])
    y_1d = np.linspace(bounds[1][0] + dy / 2, bounds[1][1] - dy / 2, shape[1])
    xg, yg = np.meshgrid(x_1d, y_1d, indexing="ij")
    ax.contour(
        xg.T,
        yg.T,
        coupling_field.T,
        levels=[G0 * 0.1, G0 * 0.5],
        colors="cyan",
        linewidths=0.8,
        linestyles="--",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    t_final = float(data.times[-1])
    ax.set_title(rf"$|\chi|$ at $t={t_final:.0f}$ + coupling $G$ contours")
    plt.colorbar(im, ax=ax)


def _plot_summary_text(  # noqa: PLR0913, PLR0917
    ax: Axes,
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    mixing: MixingResult | None,
    params: dict[str, float],
    peak_idx: int,
) -> None:
    m_phi2 = params["mPhi2"]
    m_chi2 = params["mChi2"]
    g0 = params["g0"]
    r_param = params["R"]
    omega_k = np.sqrt(K0**2 + m_phi2)
    h0 = hamiltonian.total_energy[0]
    max_drift = float(
        np.max(np.abs((hamiltonian.total_energy - h0) / max(abs(h0), 1e-30)))
    )
    lines = [
        "Parameters:",
        f"  $m_\\phi^2 = {m_phi2}$,  $m_\\chi^2 = {m_chi2}$",
        f"  $g_0 = {g0}$,  $R = {r_param}$",
        "",
        f"Wave packet: $k_0 = {K0}$",
        f"  $\\omega = {omega_k:.4f}$,  $v_g = {K0 / omega_k:.4f}$",
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
        f"  max $|\\Delta H / H_0| = {max_drift:.2e}$",
        "",
        "Initial condition:",
        "  Gaussian $\\phi$ wave packet, $\\chi = 0$",
        "  (pure conversion scenario)",
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


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run simulation and perform full measurement analysis."""
    print("Running coupled scattering simulation (2+1D)...")
    data, params, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(f"  {data.n_snapshots} snapshots over t=[0, {float(data.times[-1]):.1f}]")

    print("Computing coupling field G(x,y)...")
    coupling_field = _compute_coupling_field(data)

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
        # Higher floor than default 1e-12: in 2D, radial binning spreads
        # a directional wave packet across circular shells, leaving most bins
        # with near-zero initial energy.  1e-4 keeps only meaningful modes.
        spectral_conv = compute_spectral_conversion(
            data,
            "phi_0",
            "chi_0",
            energy_floor=1e-4,
        )
    except ValueError as e:
        print(f"  Spectral conversion: not computed ({e})")

    print("Computing dispersion relation...")
    disp_phi: DispersionResult | None = None
    try:
        disp_phi = compute_dispersion(data, "phi_0")
    except ValueError as e:
        print(f"  Dispersion (phi): not computed ({e})")

    print("Computing Hamiltonian energy (manual decomposition)...")
    hamiltonian = _compute_hamiltonian_energies(data, coupling_field)

    # Cross-validate with measurement module's virial energy
    print("Cross-validating with virial energy (measurement module)...")
    virial = compute_system_energy(data, 0)
    h_manual = hamiltonian.total_energy[0]
    h_virial = virial.total
    rel_diff = abs(h_manual - h_virial) / max(abs(h_manual), 1e-30)
    print(f"  Manual H(0) = {h_manual:.6f}")
    print(f"  Virial H(0) = {h_virial:.6f}")
    print(f"  Relative difference: {rel_diff:.2e}")

    _print_summary(
        result,
        hamiltonian,
        mixing,
        spectrum,
        spectral_conv,
        disp_phi,
        params,
    )

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        result,
        hamiltonian,
        mixing,
        spectrum,
        spectral_conv,
        disp_phi,
        coupling_field,
        params,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
