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
- **Energy decomposition** via ``compute_energy_timeseries`` (measurement module)
- **Energy conservation** via ``check_energy_conservation``

Note: Spectral conversion P(k,t) and dispersion omega(k) are NOT available for
this system because the position-dependent Gaussian coupling G(x,y) breaks
spatial translation invariance required by FFT-based analysis.

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
    EnergyDiagnostics,
    SimulationData,
    check_energy_conservation,
    compute_conversion_probability,
    compute_energy_timeseries,
    compute_mixing_length,
    compute_mixing_spectrum,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._mixing import MixingResult, MixingSpectrum

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
N_CELLS = 96

# Time integration
T_END = 500.0
DT = 0.02
TRACKER_INTERVAL = 5.0

# Initial conditions — rightward-propagating phi wave packet
X0 = -25.0  # center of wave packet (left of coupling region)
SIGMA = 4.0  # Gaussian envelope width
K0 = 3.0  # wavevector (k0 > sqrt(mPhi2) for propagation)
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "coupled_scattering_measurement.png"


# ── Energy decomposition (via measurement module) ────────────


@dataclass(frozen=True)
class EnergyDecomposition:
    """Energy decomposition by field group via the generic virial formula.

    Computed by ``compute_energy_timeseries``, which reads operators directly
    from the JSON spec.  Automatically handles position-dependent coefficients
    like the Gaussian coupling G(x,y).
    """

    times: NDArray[np.float64]
    phi_energy: NDArray[np.float64]
    chi_energy: NDArray[np.float64]
    coupling_energy: NDArray[np.float64]
    total_energy: NDArray[np.float64]


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


def _compute_energy_decomposition(data: SimulationData) -> EnergyDecomposition:
    """Compute energy decomposition grouped by field sector.

    Uses the generic virial formula via ``compute_energy_timeseries``,
    which reads operators from the JSON spec.  This correctly handles
    position-dependent coefficients and all spatial operators.
    """
    times, per_field, interaction, total = compute_energy_timeseries(data)

    phi_energy = per_field["phi_0"]
    chi_energy = per_field["chi_0"]

    return EnergyDecomposition(
        times=times,
        phi_energy=phi_energy,
        chi_energy=chi_energy,
        coupling_energy=interaction,
        total_energy=total,
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


def _print_summary(
    result: ConversionResult,
    energy: EnergyDecomposition,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
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

    h0 = energy.total_energy[0]
    h_final = energy.total_energy[-1]
    print(f"  Energy density ⟨ε⟩(0) = {h0:.6f},  ⟨ε⟩(T) = {h_final:.6f}")
    print(f"  max |dE/E0| = {diag.max_relative_error:.2e}")
    print(f"  Conserved: {diag.is_conserved}")
    print("=" * 65)


# ── Plotting ──────────────────────────────────────────────────


def _plot_results(
    data: SimulationData,
    result: ConversionResult,
    energy: EnergyDecomposition,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    coupling_field: NDArray[np.float64],
    params: dict[str, float],
) -> Path:
    """Generate 2x3 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
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

    # [0,1] Energy decomposition
    _plot_hamiltonian(axes[0, 1], energy)

    # [0,2] Energy conservation via check_energy_conservation
    _plot_energy_conservation(axes[0, 2], diag)

    # [1,0] Mixing spectrum
    _plot_mixing_spectrum(axes[1, 0], spectrum)

    # [1,1] |chi| at final time + coupling contours
    _plot_chi_heatmap(axes[1, 1], data, coupling_field)

    # [1,2] Summary text
    _plot_summary_text(axes[1, 2], result, energy, diag, mixing, params, peak_idx)

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


def _plot_hamiltonian(ax: Axes, h: EnergyDecomposition) -> None:
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
    ax.set_ylabel("Energy density")
    ax.set_title("Energy Decomposition")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)


def _plot_energy_conservation(ax: Axes, diag: EnergyDiagnostics) -> None:
    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(1e-3, color="r", linestyle="--", alpha=0.5, label="threshold (1e-3)")
    ax.axhline(-1e-3, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E\,/\,E_0$")
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


def _plot_summary_text(
    ax: Axes,
    result: ConversionResult,
    energy: EnergyDecomposition,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    params: dict[str, float],
    peak_idx: int,
) -> None:
    m_phi2 = params["mPhi2"]
    m_chi2 = params["mChi2"]
    g0 = params["g0"]
    r_param = params["R"]
    omega_k = np.sqrt(K0**2 + m_phi2)
    _ = energy  # used for type consistency; total energy printed via diag
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
        f"  max $|\\Delta E / E_0| = {diag.max_relative_error:.2e}$",
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

    print("Computing energy decomposition (virial formula)...")
    energy = _compute_energy_decomposition(data)

    print("Checking energy conservation...")
    diag = check_energy_conservation(data)

    _print_summary(
        result,
        energy,
        diag,
        mixing,
        spectrum,
        params,
    )

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        result,
        energy,
        diag,
        mixing,
        spectrum,
        coupling_field,
        params,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
