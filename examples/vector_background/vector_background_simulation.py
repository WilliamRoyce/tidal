"""Measurement example: scalar-vector conversion with tanh domain wall (2+1D).

Demonstrates the measurement module applied to a vector background field:

    L = -1/2 (d phi)^2 - mPhi2/2 phi^2
        -1/2 (d_a A_b)(d^a A^b) + 1/2 (d_a A_b)(d^b A^a)
        - mA2/2 A_a A^a + gBV * B^a A_a phi

where B_a(x,y) = (0, 0, B0 * tanh(x/W) * exp(-y^2/(2R^2))) is an external
vector background with a tanh domain wall in x, localized in y by a Gaussian
envelope.  This is qualitatively different from both:
  - Gaussian (coupled_scattering): symmetric, exponential tails
  - Lorentzian (proca_background): symmetric, algebraic tails
  - Tanh wall: antisymmetric, transitions from -B0 to +B0 across x=0

Only the B_2 component is nonzero, so after VarD:
  - phi eq: ... + gBV * B_2(x,y) * A_2  (position-dependent identity)
  - A_2 eq: ... + gBV * B_2(x,y) * phi  (position-dependent identity)
  - A_0, A_1 eqs: coupling vanishes (B_0 = B_1 = 0)

A rightward-propagating phi wave packet starts at x=-15 (negative coupling
region) and crosses the domain wall at x=0.  The sign-changing background
creates asymmetric conversion behavior.

The measurement module quantifies conversion using:
- **Conversion probability** P(t) = E_A2(t) / E_phi(0)
- **Hamiltonian energy** decomposition (manual, cross-validated against virial)

Spectral measurements (P(k,t), omega(k), mixing length) are NOT available
because the position-dependent background breaks translation invariance.

Usage:
    uv run python examples/vector_background/vector_background_simulation.py
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
    compute_system_energy,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult

# NOTE: Only A_2 couples to phi (B_0=B_1=0).  A_1 should remain at zero
# throughout the simulation.  This is a key validation check.

# -- Configuration --------------------------------------------------------

# Physics
MPHI2 = 1.0
MA2 = 2.0
GBV = 0.5
B0 = 1.0
WALL_WIDTH = 3.0
COUPLING_RADIUS = 8.0

PARAMS: dict[str, float] = {
    "mPhi2": MPHI2,
    "mA2": MA2,
    "gBV": GBV,
    "B0": B0,
    "W": WALL_WIDTH,
    "R": COUPLING_RADIUS,
}

# Grid (2D, periodic, centered at origin)
DOMAIN = (-30.0, 30.0)
N_CELLS = 64

# Time integration
T_END = 20.0
DT = 0.02
TRACKER_INTERVAL = 1.0

# Initial conditions -- rightward-propagating phi wave packet
X0 = -15.0  # left of domain wall
SIGMA = 3.0
K0 = 2.0  # wavevector (k0 > sqrt(mPhi2) for propagation)
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "vector_background_measurement.png"


# -- Hamiltonian energy (manual) ------------------------------------------


@dataclass(frozen=True)
class HamiltonianDecomposition:
    """Hamiltonian energy decomposition over time.

    H = 1/2 pi_phi^2 + 1/2 |grad phi|^2 + mPhi2/2 phi^2
      + sum_{i=1,2} [1/2 pi_Ai^2 + 1/2 |grad A_i|^2 + mA2/2 A_i^2]
      - gBV * integral B_2(x,y) * A_2 * phi dA

    Note: A_0 is a constraint and has no kinetic energy.
    Only B_2 is nonzero so coupling reduces to a single integral.
    """

    times: NDArray[np.float64]
    phi_energy: NDArray[np.float64]
    a_energy: NDArray[np.float64]
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
    """Compute B_2(x,y) = B0 * tanh(x/W) * exp(-y^2/(2R^2)) on the grid."""
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
    b0 = data.parameters.get("B0", B0)
    w = data.parameters.get("W", WALL_WIDTH)
    r_param = data.parameters.get("R", COUPLING_RADIUS)
    return b0 * np.tanh(x / w) * np.exp(-(y**2) / (2 * r_param**2))


def _compute_hamiltonian_energies(  # noqa: PLR0914
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> HamiltonianDecomposition:
    """Compute full Hamiltonian decomposition over all snapshots.

    A_0 is a constraint (no kinetic energy).  A_1 has no background
    coupling (B_1=0), but we still include its free energy.
    """
    dx, dy = data.grid_spacing
    dv = dx * dy
    n = data.n_snapshots
    phi_e = np.empty(n)
    a_e = np.empty(n)
    cpl_e = np.empty(n)

    m_phi2 = data.parameters.get("mPhi2", MPHI2)
    m_a2 = data.parameters.get("mA2", MA2)
    gbv = data.parameters.get("gBV", GBV)

    for i in range(n):
        # Phi sector
        phi = data.fields["phi_0"][i]
        pi_phi = data.momenta["phi_0"][i]

        phi_e[i] = (
            0.5 * float(np.sum(pi_phi**2)) * dv
            + _gradient_energy_2d(phi, dx, dy)
            + 0.5 * m_phi2 * float(np.sum(phi**2)) * dv
        )

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

        # Coupling: -gBV * integral B_2(x,y) * A_2 * phi dA
        # Sign: Lagrangian has +gBV*B^a*A_a*phi, with B^2 = B_2 (upper index
        # in Minkowski flat spatial).  Hamiltonian flips sign of potential.
        cpl_e[i] = -gbv * float(np.sum(coupling_field * a2 * phi)) * dv

    return HamiltonianDecomposition(
        times=data.times,
        phi_energy=phi_e,
        a_energy=a_e,
        coupling_energy=cpl_e,
        total_energy=phi_e + a_e + cpl_e,
    )


# -- Simulation -----------------------------------------------------------


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    """Create ICs: traveling phi wave packet, A = 0 (pure conversion)."""
    coords = cast("np.ndarray", grid.cell_coords)
    x = np.asarray(coords[..., 0], dtype=float)
    y = np.asarray(coords[..., 1], dtype=float)

    envelope = PULSE_AMPLITUDE * np.exp(
        -((x - X0) ** 2) / (2 * SIGMA**2) - y**2 / (2 * SIGMA**2)
    )
    phi_data = envelope * np.cos(K0 * (x - X0))
    omega = np.sqrt(K0**2 + MPHI2)
    pi_phi_data = omega * envelope * np.sin(K0 * (x - X0))

    # State layout: phi_0, pi_phi_0, A_0 (constraint), A_1, pi_A1, A_2, pi_A2
    return FieldCollection(
        [
            ScalarField(grid, data=phi_data, label="phi_0"),
            ScalarField(grid, data=pi_phi_data, label="pi_phi_0"),
            ScalarField(grid, data=np.zeros_like(x), label="A_0"),
            ScalarField(grid, data=np.zeros_like(x), label="A_1"),
            ScalarField(grid, data=np.zeros_like(x), label="pi_A_1"),
            ScalarField(grid, data=np.zeros_like(x), label="A_2"),
            ScalarField(grid, data=np.zeros_like(x), label="pi_A_2"),
        ]
    )


def _run_simulation() -> tuple[SimulationData, dict[str, float], Path]:
    """Run the scalar-vector + tanh domain wall simulation.

    Raises
    ------
    FileNotFoundError
        If the JSON spec has not been generated yet.
    """
    json_path = Path(__file__).parent.parent / "data" / "vector_background.json"
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
        Path(__file__).parent.parent / "data" / "vector_background_output"
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
    print("Vector Background Domain Wall: Measurement Summary")
    print("=" * 65)
    print()

    m_phi2 = params["mPhi2"]
    m_a2 = params["mA2"]
    gbv = params["gBV"]
    b0 = params["B0"]
    w = params["W"]
    r_param = params["R"]
    print(f"  mPhi2={m_phi2}, mA2={m_a2}, gBV={gbv}, B0={b0}, W={w}, R={r_param}")
    max_coupling = gbv * b0
    print(f"  Max coupling (x>>W): gBV*B0 = {max_coupling:.2f}")
    print(f"  Stability: mPhi2*mA2 = {m_phi2 * m_a2:.1f} > (gBV*B0)^2 = {max_coupling**2:.2f}")
    print()

    omega_k = np.sqrt(K0**2 + m_phi2)
    v_group = K0 / omega_k
    print(f"  Wave packet: k0 = {K0}, omega = {omega_k:.4f}, v_group = {v_group:.4f}")
    print()

    peak_idx = int(np.argmax(result.probability))
    print(f"  Peak conversion P(t) = {result.probability[peak_idx]:.6f}")
    print(f"    at t = {result.times[peak_idx]:.2f}")
    print(f"  Initial source energy E_phi(0) = {result.source_energy[0]:.4f}")
    print(f"  Final  target energy E_A2(T) = {result.target_energy[-1]:.4f}")
    print()

    print("  Spectral measurements: NOT AVAILABLE")
    print("    (tanh domain wall breaks translation invariance)")
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
        "Scalar-Vector + Tanh Domain Wall Background (2+1D)\n"
        r"$\mathcal{L} = -\frac{1}{2}(\partial\phi)^2"
        r" - \frac{m_\phi^2}{2}\phi^2"
        r" - \frac{1}{2}(\partial_a A_b)(\partial^a A^b)"
        r" + \ldots"
        r" + g_{BV}\,B^a A_a\,\phi$",
        fontsize=13,
    )

    # [0,0] Conversion probability P(t)
    _plot_conversion(axes[0, 0], result)

    # [0,1] Hamiltonian energy decomposition
    _plot_hamiltonian(axes[0, 1], hamiltonian)

    # [0,2] B_2(x,y) domain wall coupling field
    _plot_coupling_field(axes[0, 2], fig, data, coupling_field)

    # [1,0] A_2 at final time (coupled component)
    _plot_field_heatmap(axes[1, 0], data, "A_2", coupling_field, "A_2 (coupled)")

    # [1,1] A_1 at final time (uncoupled -- should stay ~0)
    _plot_field_heatmap(axes[1, 1], data, "A_1", coupling_field, "A_1 (uncoupled)")

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
    ax.set_ylabel(r"$P(t) = E_{A_2}(t)\,/\,E_\phi(0)$")
    ax.set_title(r"Conversion Probability ($\phi \to A_2$)")
    ax.set_ylim(bottom=0)
    ax.grid(visible=True, alpha=0.3)


def _plot_hamiltonian(ax: Axes, h: HamiltonianDecomposition) -> None:
    ax.plot(h.times, h.phi_energy, "b-", label=r"$E_\phi$", linewidth=1.2)
    ax.plot(h.times, h.a_energy, "r-", label=r"$E_A$", linewidth=1.2)
    ax.plot(
        h.times,
        h.coupling_energy,
        "g--",
        label=r"$-g\!\int\! B_2\,A_2\,\phi\,dA$",
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


def _plot_coupling_field(
    ax: Axes,
    fig: plt.Figure,
    data: SimulationData,
    coupling_field: NDArray[np.float64],
) -> None:
    """Plot B_2(x,y) with antisymmetric contours."""
    bounds = data.grid_bounds
    extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    vmax = float(np.max(np.abs(coupling_field)))
    im = ax.imshow(
        coupling_field.T,
        aspect="equal",
        origin="lower",
        extent=extent,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    fig.colorbar(im, ax=ax, label=r"$B_2(x,y)$")
    ax.axvline(0, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Domain wall $B_2 = B_0\tanh(x/W)\,e^{-y^2/(2R^2)}$")


def _plot_field_heatmap(
    ax: Axes,
    data: SimulationData,
    field_name: str,
    coupling_field: NDArray[np.float64],
    title_label: str,
) -> None:
    """Plot |field| at final time with coupling B_2(x,y) contours."""
    final = data.fields[field_name][-1]
    bounds = data.grid_bounds
    extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
    vmax = max(float(np.max(np.abs(final))), 1e-10)
    im = ax.imshow(
        np.abs(final).T,
        aspect="equal",
        origin="lower",
        extent=extent,
        cmap="inferno",
        vmin=0,
        vmax=vmax,
    )
    # Domain wall contours
    shape = final.shape
    dx, dy = data.grid_spacing
    x_1d = np.linspace(bounds[0][0] + dx / 2, bounds[0][1] - dx / 2, shape[0])
    y_1d = np.linspace(bounds[1][0] + dy / 2, bounds[1][1] - dy / 2, shape[1])
    xg, yg = np.meshgrid(x_1d, y_1d, indexing="ij")
    bg_max = float(np.max(np.abs(coupling_field)))
    if bg_max > 0:
        ax.contour(
            xg.T,
            yg.T,
            coupling_field.T,
            levels=[-bg_max * 0.5, 0, bg_max * 0.5],
            colors="cyan",
            linewidths=0.8,
            linestyles="--",
        )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    t_final = float(data.times[-1])
    ax.set_title(rf"$|{title_label}|$ at $t={t_final:.0f}$")
    plt.colorbar(im, ax=ax)


def _plot_summary_text(
    ax: Axes,
    result: ConversionResult,
    hamiltonian: HamiltonianDecomposition,
    params: dict[str, float],
) -> None:
    m_phi2 = params["mPhi2"]
    m_a2 = params["mA2"]
    gbv = params["gBV"]
    b0 = params["B0"]
    w = params["W"]
    r_param = params["R"]
    omega_k = np.sqrt(K0**2 + m_phi2)
    peak_idx = int(np.argmax(result.probability))
    h0 = hamiltonian.total_energy[0]
    max_drift = float(
        np.max(np.abs((hamiltonian.total_energy - h0) / max(abs(h0), 1e-30)))
    )
    lines = [
        "Parameters:",
        f"  $m_\\phi^2 = {m_phi2}$,  $m_A^2 = {m_a2}$",
        f"  $g_{{BV}} = {gbv}$,  $B_0 = {b0}$",
        f"  $W = {w}$,  $R = {r_param}$",
        "",
        f"Wave packet: $k_0 = {K0}$",
        f"  $\\omega = {omega_k:.4f}$,  $v_g = {K0 / omega_k:.4f}$",
        "",
        "Results:",
        f"  Peak $P(t) = {result.probability[peak_idx]:.6f}$",
        f"  at $t = {result.times[peak_idx]:.2f}$",
        f"  max $|\\Delta H / H_0| = {max_drift:.2e}$",
        "",
        "Domain wall: tanh(x/W)",
        "  coupling changes sign at x=0",
        "  A_1 should stay $\\approx 0$ (B_1=0)",
        "",
        "Spectral methods: N/A",
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
    print("Running scalar-vector + tanh domain wall simulation (2+1D)...")
    data, params, data_dir = _run_simulation()
    print(f"  Data saved to: {data_dir}")
    print(f"  {data.n_snapshots} snapshots over t=[0, {float(data.times[-1]):.1f}]")

    print("Computing coupling field B_2(x,y)...")
    coupling_field = _compute_coupling_field(data)

    print("Computing conversion probability (phi -> A_2)...")
    result = compute_conversion_probability(data, source="phi_0", target="A_2")

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

    # Check that A_1 stays at zero (B_1=0, uncoupled)
    a1_max = float(np.max(np.abs(data.fields["A_1"])))
    a2_max = float(np.max(np.abs(data.fields["A_2"])))
    print(f"  A_1 max (uncoupled, should be ~0): {a1_max:.2e}")
    print(f"  A_2 max (coupled): {a2_max:.4f}")

    _print_summary(result, hamiltonian, params)

    print()
    print("Generating measurement plots...")
    output_path = _plot_results(data, result, hamiltonian, coupling_field, params)
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    main()
