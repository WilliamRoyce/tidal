"""Measurement example: vector conversion with Lorentzian background (2+1D).

Demonstrates the measurement module applied to a position-dependent coupling
between two massive Proca vector fields:

    L = -1/4 F^A_{ab} F^{A,ab} - 1/4 F^B_{ab} F^{B,ab}
        - mA2/2 A_a A^a - mB2/2 B_a B^a
        + gcoup * G(x,y) * A_a B^a

where G(x,y) = g0 / (1 + (x^2+y^2)/R^2) is a Lorentzian profile centered
at the origin.  Unlike the Gaussian in coupled_scattering, the Lorentzian
has algebraic (1/r^2) tails rather than exponential decay.

A Gaussian pulse in A_1 is placed at the center of the coupling region.
B starts at zero -- any B that appears is purely from conversion.
The measurement module quantifies this conversion using:

- **Conversion probability** P(t) = E_B(t) / E_A(0) via group conversion
- **Mixing length** and **mixing spectrum** (temporal FFT of P(t))
- **Energy conservation** via ``check_energy_conservation``
- **Hamiltonian energy** decomposition via the generic virial formula
  (from the measurement module's ``compute_energy_timeseries``)

Usage:
    uv run python examples/proca_background/proca_background_simulation.py
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
    check_energy_conservation,
    compute_energy_timeseries,
    compute_group_conversion,
    compute_mixing_length,
    compute_mixing_spectrum,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._diagnostics import EnergyDiagnostics
    from tidal.measurement._mixing import MixingResult, MixingSpectrum

# -- Configuration --------------------------------------------------------

# Physics
MA2 = 1.0
MB2 = 2.0
GCOUP = 0.5
G0 = 1.0
COUPLING_RADIUS = 8.0

PARAMS: dict[str, float] = {
    "mA2": MA2,
    "mB2": MB2,
    "gcoup": GCOUP,
    "g0": G0,
    "R": COUPLING_RADIUS,
}

# Grid (2D, periodic, centered at origin for the Lorentzian)
DOMAIN = (-30.0, 30.0)
N_CELLS = 64

# Time integration
T_END = 20.0
DT = 0.02
TRACKER_INTERVAL = 0.2

# Initial conditions -- Gaussian pulse in A_1 at center
PULSE_AMPLITUDE = 0.5
PULSE_WIDTH = 3.0

# Analysis
ENERGY_THRESHOLD = 0.001

# Output
OUTPUT_FILENAME = "proca_background_measurement.png"


# -- Energy decomposition (via measurement module) ------------------------


@dataclass(frozen=True)
class EnergyDecomposition:
    """Energy decomposition by field group via the generic virial formula.

    Computed by ``compute_energy_timeseries``, which reads operators directly
    from the JSON spec.  Automatically handles directional laplacians (curl
    energy for Proca vectors) and position-dependent coefficients.
    """

    times: NDArray[np.float64]
    a_energy: NDArray[np.float64]
    b_energy: NDArray[np.float64]
    coupling_energy: NDArray[np.float64]
    total_energy: NDArray[np.float64]


def _compute_coupling_field(data: SimulationData) -> NDArray[np.float64]:
    """Compute G(x,y) = g0 / (1 + (x^2+y^2)/R^2) on the simulation grid."""
    shape = data.fields["A_1"].shape[1:]  # spatial (nx, ny)
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
    return g0 / (1.0 + (x**2 + y**2) / r_param**2)


def _compute_energy_decomposition(data: SimulationData) -> EnergyDecomposition:
    """Compute energy decomposition grouped by field sector.

    Uses the generic virial formula via ``compute_energy_timeseries``,
    which reads operators from the JSON spec.  This correctly handles
    directional laplacians (Proca curl energy), cross-derivatives,
    constraint self-energy, and position-dependent coefficients.
    """
    times, per_field, interaction, total = compute_energy_timeseries(data)

    a_fields = [n for n in per_field if n.startswith("A_")]
    b_fields = [n for n in per_field if n.startswith("B_")]

    a_energy = cast("NDArray[np.float64]", sum(per_field[n] for n in a_fields))
    b_energy = cast("NDArray[np.float64]", sum(per_field[n] for n in b_fields))

    return EnergyDecomposition(
        times=times,
        a_energy=a_energy,
        b_energy=b_energy,
        coupling_energy=interaction,
        total_energy=total,
    )


# -- Simulation -----------------------------------------------------------


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    """Create ICs: Gaussian A_1 pulse at center, all B = 0."""
    coords = cast("np.ndarray", grid.cell_coords)
    x = np.asarray(coords[..., 0], dtype=float)
    y = np.asarray(coords[..., 1], dtype=float)

    gaussian = PULSE_AMPLITUDE * np.exp(-(x**2 + y**2) / (2 * PULSE_WIDTH**2))

    # State layout: A_0, A_1, pi_A1, A_2, pi_A2, B_0, B_1, pi_B1, B_2, pi_B2
    return FieldCollection(
        [
            ScalarField(grid, data=0.0, label="A_0_field"),  # A_0 constraint
            ScalarField(grid, data=gaussian, label="A_1_field"),  # A_1 pulse
            ScalarField(grid, data=0.0, label="A_1_momentum"),  # pi_A1
            ScalarField(grid, data=0.0, label="A_2_field"),  # A_2
            ScalarField(grid, data=0.0, label="A_2_momentum"),  # pi_A2
            ScalarField(grid, data=0.0, label="B_0_field"),  # B_0 constraint
            ScalarField(grid, data=0.0, label="B_1_field"),  # B_1
            ScalarField(grid, data=0.0, label="B_1_momentum"),  # pi_B1
            ScalarField(grid, data=0.0, label="B_2_field"),  # B_2
            ScalarField(grid, data=0.0, label="B_2_momentum"),  # pi_B2
        ]
    )


def _run_simulation() -> tuple[SimulationData, dict[str, float], Path]:
    """Run the coupled Proca + Lorentzian background simulation.

    Raises
    ------
    FileNotFoundError
        If the JSON spec has not been generated yet.
    """
    json_path = Path(__file__).parent.parent / "data" / "proca_background.json"
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

    output_data_dir = Path(__file__).parent.parent.parent / "outputs" / "proca_background"
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


# -- Summary --------------------------------------------------------------


def _print_summary(
    result: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    params: dict[str, float],
) -> None:
    """Print quantitative measurement summary to stdout."""
    print("=" * 65)
    print("Proca + Lorentzian Background: Measurement Summary")
    print("=" * 65)
    print()

    m_a2 = params["mA2"]
    m_b2 = params["mB2"]
    gcoup = params["gcoup"]
    g0 = params["g0"]
    r_param = params["R"]
    print(f"  mA2={m_a2}, mB2={m_b2}, gcoup={gcoup}, g0={g0}, R={r_param}")
    print(
        f"  Stability: mA2*mB2 = {m_a2 * m_b2:.1f} > (g0*gcoup)^2 = {(g0 * gcoup) ** 2:.2f}"
    )
    print()

    peak_idx = int(np.argmax(result.probability))
    print(f"  Peak conversion P(t) = {result.probability[peak_idx]:.6f}")
    print(f"    at t = {result.times[peak_idx]:.2f}")
    print(f"  Initial source energy E_A(0) = {result.source_energy[0]:.4f}")
    print(f"  Final  target energy E_B(T) = {result.target_energy[-1]:.4f}")
    print()

    if mixing is not None:
        print("  Mixing length:")
        print(f"    L_mix = {mixing.mixing_length:.4f}")
        print(f"    omega_dom = {mixing.dominant_frequency:.4f}")
        print(f"    uncertainty = {mixing.mixing_length_uncertainty:.4f}")
    else:
        print("  Mixing length: NOT EXTRACTED")
    print()

    print(f"  Energy conservation: max |dE/E0| = {diag.max_relative_error:.2e}")
    print(f"  Conservation: {'PASS' if diag.is_conserved else 'FAIL'}")
    print("  (virial formula from measurement module)")
    print("=" * 65)


# -- Plotting -------------------------------------------------------------


def _plot_results(
    data: SimulationData,
    result: ConversionResult,
    hamiltonian: EnergyDecomposition,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    spectrum: MixingSpectrum | None,
    coupling_field: NDArray[np.float64],
    params: dict[str, float],
) -> Path:
    """Generate 2x4 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    fig.suptitle(
        "Proca + Lorentzian Background (2+1D): Measurement Analysis\n"
        r"$\mathcal{L} = -\frac{1}{4}F_A^2 - \frac{1}{4}F_B^2"
        r" - \frac{m_A^2}{2}A^2 - \frac{m_B^2}{2}B^2"
        r" + g\,G(x,y)\,A\!\cdot\!B$",
        fontsize=13,
    )

    peak_idx = int(np.argmax(result.probability))

    # [0,0] Total P(t) + mixing length annotation
    _plot_conversion(axes[0, 0], result, mixing)

    # [0,1] Energy conservation via check_energy_conservation
    _plot_energy_conservation(axes[0, 1], diag)

    # [0,2] Mixing spectrum (temporal FFT of P(t))
    _plot_mixing_spectrum(axes[0, 2], spectrum)

    # [0,3] Energy decomposition (E_A, E_B, coupling, total)
    _plot_hamiltonian(axes[0, 3], hamiltonian)

    # [1,0] Coupling G(x,y) heatmap
    _plot_coupling_field(axes[1, 0], fig, data, coupling_field)

    # [1,1] |B_1| at peak conversion time
    _plot_field_heatmap(axes[1, 1], data, "B_1", coupling_field, snapshot_idx=peak_idx)

    # [1,2] Summary text
    _plot_summary_text(axes[1, 2], result, diag, mixing, params)

    # [1,3] blank
    axes[1, 3].axis("off")

    plt.tight_layout()
    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _plot_conversion(
    ax: Axes, result: ConversionResult, mixing: MixingResult | None
) -> None:
    peak_idx = int(np.argmax(result.probability))
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
    ax.set_ylabel(r"$P(t) = E_B(t)\,/\,E_A(0)$")
    ax.set_title("Conversion Probability (A -> B)")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)


def _plot_energy_conservation(ax: Axes, diag: EnergyDiagnostics) -> None:
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


def _plot_mixing_spectrum(ax: Axes, spectrum: MixingSpectrum | None) -> None:
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


def _plot_hamiltonian(ax: Axes, h: EnergyDecomposition) -> None:
    ax.plot(h.times, h.a_energy, "b-", label=r"$E_A$", linewidth=1.2)
    ax.plot(h.times, h.b_energy, "r-", label=r"$E_B$", linewidth=1.2)
    ax.plot(
        h.times,
        h.coupling_energy,
        "g--",
        label=r"$E_\mathrm{int}$",
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


def _plot_coupling_field(
    ax: Axes,
    fig: Figure,
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> None:
    bounds = data.grid_bounds
    extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    im = ax.imshow(
        coupling_field.T,
        aspect="equal",
        origin="lower",
        extent=extent,
        cmap="YlOrRd",
    )
    fig.colorbar(im, ax=ax, label=r"$G(x,y)$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Coupling $G(x,y) = g_0 / (1 + r^2/R^2)$")


def _plot_field_heatmap(
    ax: Axes,
    data: SimulationData,
    field_name: str,
    coupling_field: NDArray[np.float64],
    snapshot_idx: int = -1,
) -> None:
    """Plot |field| at given snapshot with coupling G(x,y) contours."""
    final = data.fields[field_name][snapshot_idx]
    bounds = data.grid_bounds
    extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    im = ax.imshow(
        np.abs(final).T,
        aspect="equal",
        origin="lower",
        extent=extent,
        cmap="inferno",
    )
    # Coupling contours
    shape = final.shape
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
    t_snap = float(data.times[snapshot_idx])
    ax.set_title(rf"$|{field_name}|$ at $t={t_snap:.0f}$ (peak) + coupling contours")
    plt.colorbar(im, ax=ax)


def _plot_summary_text(
    ax: Axes,
    result: ConversionResult,
    diag: EnergyDiagnostics,
    mixing: MixingResult | None,
    params: dict[str, float],
) -> None:
    m_a2 = params["mA2"]
    m_b2 = params["mB2"]
    gcoup = params["gcoup"]
    g0 = params["g0"]
    r_param = params["R"]
    peak_idx = int(np.argmax(result.probability))
    lines = [
        "Parameters:",
        f"  $m_A^2 = {m_a2}$,  $m_B^2 = {m_b2}$",
        f"  $g_{{coup}} = {gcoup}$,  $g_0 = {g0}$,  $R = {r_param}$",
        "",
        f"Stability: $m_A^2 m_B^2 = {m_a2 * m_b2:.1f}$",
        f"  $> (g_0 g_{{coup}})^2 = {(g0 * gcoup) ** 2:.2f}$",
        "",
        "Results:",
        f"  Peak $P(t) = {result.probability[peak_idx]:.6f}$",
        f"  at $t = {result.times[peak_idx]:.2f}$",
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
        "Background: Lorentzian $1/(1+r^2/R^2)$",
        "  (algebraic tails, non-compact)",
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


# -- Entry point ----------------------------------------------------------


def main() -> None:
    """Run simulation and perform measurement analysis."""
    print("Running Proca + Lorentzian background simulation (2+1D)...")
    data, params, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(f"  {data.n_snapshots} snapshots over t=[0, {float(data.times[-1]):.1f}]")

    print("Computing coupling field G(x,y)...")
    coupling_field = _compute_coupling_field(data)

    print("Computing group conversion probability (A -> B)...")
    result = compute_group_conversion(
        data,
        source_fields=["A_1", "A_2"],
        target_fields=["B_1", "B_2"],
    )

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
    hamiltonian = _compute_energy_decomposition(data)

    print("Checking energy conservation...")
    diag = check_energy_conservation(data, threshold=ENERGY_THRESHOLD)

    _print_summary(result, diag, mixing, params)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(
        data,
        result,
        hamiltonian,
        diag,
        mixing,
        spectrum,
        coupling_field,
        params,
    )
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
