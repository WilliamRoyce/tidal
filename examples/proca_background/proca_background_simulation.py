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
- **Hamiltonian energy** decomposition (manual, cross-validated against
  virial energy from the measurement module)

Spectral measurements (P(k,t), omega(k), mixing length) are NOT available
for this system because the position-dependent Lorentzian background breaks
translation invariance.  Only real-space measurements are physical.

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
    compute_group_conversion,
    compute_system_energy,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult

# NOTE: The measurement module's virial energy (compute_system_energy) now
# supports position-dependent coupling G(x,y).  We keep the manual
# Hamiltonian decomposition as a cross-validation -- it gives a physics-
# specific breakdown (A energy, B energy, coupling energy) that the
# generic virial formula does not.

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
TRACKER_INTERVAL = 1.0

# Initial conditions -- Gaussian pulse in A_1 at center
PULSE_AMPLITUDE = 0.5
PULSE_WIDTH = 3.0

# Output
OUTPUT_FILENAME = "proca_background_measurement.png"


# -- Hamiltonian energy (manual) ------------------------------------------


@dataclass(frozen=True)
class HamiltonianDecomposition:
    """Hamiltonian energy decomposition over time.

    H = sum_i [1/2 pi_Ai^2 + 1/2 |grad A_i|^2 + mA2/2 A_i^2]
      + sum_i [1/2 pi_Bi^2 + 1/2 |grad B_i|^2 + mB2/2 B_i^2]
      - gcoup * integral G(x,y) (A_1*B_1 + A_2*B_2) dA

    Note: A_0 and B_0 are constraints (t_order=0) and have no kinetic
    energy.  Only spatial components (A_1, A_2, B_1, B_2) contribute.
    """

    times: NDArray[np.float64]
    a_energy: NDArray[np.float64]
    b_energy: NDArray[np.float64]
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


def _compute_hamiltonian_energies(  # noqa: PLR0914
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> HamiltonianDecomposition:
    """Compute full Hamiltonian decomposition over all snapshots.

    Only spatial components (A_1, A_2, B_1, B_2) contribute to kinetic
    and mass energy.  Constraints A_0, B_0 have no conjugate momenta.
    """
    dx, dy = data.grid_spacing
    dv = dx * dy
    n = data.n_snapshots
    a_e = np.empty(n)
    b_e = np.empty(n)
    cpl_e = np.empty(n)

    m_a2 = data.parameters.get("mA2", MA2)
    m_b2 = data.parameters.get("mB2", MB2)
    gcoup = data.parameters.get("gcoup", GCOUP)

    for i in range(n):
        # A sector: A_1 and A_2 (spatial evolution components)
        a1 = data.fields["A_1"][i]
        a2 = data.fields["A_2"][i]
        pi_a1 = data.momenta["A_1"][i]
        pi_a2 = data.momenta["A_2"][i]

        a_e[i] = (
            0.5 * float(np.sum(pi_a1**2 + pi_a2**2)) * dv
            + _gradient_energy_2d(a1, dx, dy)
            + _gradient_energy_2d(a2, dx, dy)
            + 0.5 * m_a2 * float(np.sum(a1**2 + a2**2)) * dv
        )

        # B sector: B_1 and B_2
        b1 = data.fields["B_1"][i]
        b2 = data.fields["B_2"][i]
        pi_b1 = data.momenta["B_1"][i]
        pi_b2 = data.momenta["B_2"][i]

        b_e[i] = (
            0.5 * float(np.sum(pi_b1**2 + pi_b2**2)) * dv
            + _gradient_energy_2d(b1, dx, dy)
            + _gradient_energy_2d(b2, dx, dy)
            + 0.5 * m_b2 * float(np.sum(b1**2 + b2**2)) * dv
        )

        # Coupling: -gcoup * integral G(x,y) * (A_1*B_1 + A_2*B_2) dA
        # Sign: Lagrangian has +gcoup*G*A.B, Hamiltonian has -gcoup*G*A.B
        # But in Minkowski with (-,+,+): A.B = -A_0*B_0 + A_1*B_1 + A_2*B_2
        # For the potential energy term, the spatial dot product contributes.
        cpl_e[i] = -gcoup * float(
            np.sum(coupling_field * (a1 * b1 + a2 * b2))
        ) * dv

    return HamiltonianDecomposition(
        times=data.times,
        a_energy=a_e,
        b_energy=b_e,
        coupling_energy=cpl_e,
        total_energy=a_e + b_e + cpl_e,
    )


# -- Simulation -----------------------------------------------------------


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    """Create ICs: Gaussian A_1 pulse at center, all B = 0."""
    coords = cast("np.ndarray", grid.cell_coords)
    x = np.asarray(coords[..., 0], dtype=float)
    y = np.asarray(coords[..., 1], dtype=float)

    gaussian = PULSE_AMPLITUDE * np.exp(
        -(x**2 + y**2) / (2 * PULSE_WIDTH**2)
    )

    # State layout: A_0, A_1, pi_A1, A_2, pi_A2, B_0, B_1, pi_B1, B_2, pi_B2
    return FieldCollection(
        [
            ScalarField(grid, data=0.0, label="A_0_field"),       # A_0 constraint
            ScalarField(grid, data=gaussian, label="A_1_field"),   # A_1 pulse
            ScalarField(grid, data=0.0, label="A_1_momentum"),     # pi_A1
            ScalarField(grid, data=0.0, label="A_2_field"),        # A_2
            ScalarField(grid, data=0.0, label="A_2_momentum"),     # pi_A2
            ScalarField(grid, data=0.0, label="B_0_field"),        # B_0 constraint
            ScalarField(grid, data=0.0, label="B_1_field"),        # B_1
            ScalarField(grid, data=0.0, label="B_1_momentum"),     # pi_B1
            ScalarField(grid, data=0.0, label="B_2_field"),        # B_2
            ScalarField(grid, data=0.0, label="B_2_momentum"),     # pi_B2
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

    output_data_dir = (
        Path(__file__).parent.parent / "data" / "proca_background_output"
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


# -- Summary --------------------------------------------------------------


def _print_summary(
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
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
    print(f"  Stability: mA2*mB2 = {m_a2 * m_b2:.1f} > (g0*gcoup)^2 = {(g0 * gcoup)**2:.2f}")
    print()

    peak_idx = int(np.argmax(result.probability))
    print(f"  Peak conversion P(t) = {result.probability[peak_idx]:.6f}")
    print(f"    at t = {result.times[peak_idx]:.2f}")
    print(f"  Initial source energy E_A(0) = {result.source_energy[0]:.4f}")
    print(f"  Final  target energy E_B(T) = {result.target_energy[-1]:.4f}")
    print()

    # NOTE: Spectral measurements are not available for this system.
    print("  Spectral measurements: NOT AVAILABLE")
    print("    (Lorentzian background breaks translation invariance)")
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


# -- Plotting -------------------------------------------------------------


def _plot_results(
    data: SimulationData,
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    coupling_field: NDArray[np.float64],
    params: dict[str, float],
) -> Path:
    """Generate 2x3 measurement figure. Returns path to saved PNG."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(
        "Proca + Lorentzian Background (2+1D): Measurement Analysis\n"
        r"$\mathcal{L} = -\frac{1}{4}F_A^2 - \frac{1}{4}F_B^2"
        r" - \frac{m_A^2}{2}A^2 - \frac{m_B^2}{2}B^2"
        r" + g\,G(x,y)\,A\!\cdot\!B$",
        fontsize=13,
    )

    # [0,0] Conversion probability P(t)
    _plot_conversion(axes[0, 0], result)

    # [0,1] Hamiltonian energy decomposition
    _plot_hamiltonian(axes[0, 1], hamiltonian)

    # [0,2] Energy conservation |dH/H0|
    _plot_energy_conservation(axes[0, 2], hamiltonian)

    # [1,0] Lorentzian coupling G(x,y)
    _plot_coupling_field(axes[1, 0], fig, data, coupling_field)

    # [1,1] |B_1| at final time
    _plot_field_heatmap(axes[1, 1], data, "B_1", coupling_field)

    # [1,2] Summary text
    _plot_summary_text(axes[1, 2], result, hamiltonian, params)

    plt.tight_layout()
    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _plot_conversion(ax: Axes, result: ConversionResult) -> None:
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
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P(t) = E_B(t)\,/\,E_A(0)$")
    ax.set_title("Conversion Probability (A -> B)")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)


def _plot_hamiltonian(ax: Axes, h: HamiltonianDecomposition) -> None:
    ax.plot(h.times, h.a_energy, "b-", label=r"$E_A$", linewidth=1.2)
    ax.plot(h.times, h.b_energy, "r-", label=r"$E_B$", linewidth=1.2)
    ax.plot(
        h.times,
        h.coupling_energy,
        "g--",
        label=r"$-g\!\int\! G\,A\!\cdot\!B\,dA$",
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


def _plot_coupling_field(
    ax: Axes,
    fig: plt.Figure,
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
) -> None:
    """Plot |field| at final time with coupling G(x,y) contours."""
    final = data.fields[field_name][-1]
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
    t_final = float(data.times[-1])
    ax.set_title(rf"$|{field_name}|$ at $t={t_final:.0f}$ + coupling contours")
    plt.colorbar(im, ax=ax)


def _plot_summary_text(
    ax: Axes,
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    params: dict[str, float],
) -> None:
    m_a2 = params["mA2"]
    m_b2 = params["mB2"]
    gcoup = params["gcoup"]
    g0 = params["g0"]
    r_param = params["R"]
    peak_idx = int(np.argmax(result.probability))
    h0 = hamiltonian.total_energy[0]
    max_drift = float(
        np.max(np.abs((hamiltonian.total_energy - h0) / max(abs(h0), 1e-30)))
    )
    lines = [
        "Parameters:",
        f"  $m_A^2 = {m_a2}$,  $m_B^2 = {m_b2}$",
        f"  $g_{{coup}} = {gcoup}$,  $g_0 = {g0}$,  $R = {r_param}$",
        "",
        f"Stability: $m_A^2 m_B^2 = {m_a2 * m_b2:.1f}$",
        f"  $> (g_0 g_{{coup}})^2 = {(g0 * gcoup)**2:.2f}$",
        "",
        "Results:",
        f"  Peak $P(t) = {result.probability[peak_idx]:.6f}$",
        f"  at $t = {result.times[peak_idx]:.2f}$",
        f"  max $|\\Delta H / H_0| = {max_drift:.2e}$",
        "",
        "Background: Lorentzian $1/(1+r^2/R^2)$",
        "  (algebraic tails, non-compact)",
        "",
        "Spectral methods: N/A",
        "  (breaks translation invariance)",
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
        source=["A_1", "A_2"],
        target=["B_1", "B_2"],
    )

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

    _print_summary(result, hamiltonian, params)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(data, result, hamiltonian, coupling_field, params)
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
