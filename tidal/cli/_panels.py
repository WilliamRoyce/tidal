"""Render functions for ``tidal plot`` — one plot type per function.

Each function takes a matplotlib ``Axes`` and a ``SimulationData`` object,
producing a single focused visualization.  The ``tidal plot`` command
calls exactly one of these per invocation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData

# -- Shared constants -------------------------------------------------------

DPI = 150
VMAX_FLOOR = 0.01  # Minimum symmetric colour-scale magnitude
_LINE_COLORS = (
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:olive",
    "tab:cyan",
)


# -- Helpers ----------------------------------------------------------------


def field_names(data: SimulationData, requested: list[str] | None) -> list[str]:
    """Resolve field names, defaulting to all fields in spec order.

    Raises
    ------
    ValueError
        If any requested field name is not found in data.
    """
    all_names = list(data.fields.keys())
    if requested is None:
        return all_names
    unknown = [f for f in requested if f not in data.fields]
    if unknown:
        msg = (
            f"Unknown field(s): {', '.join(unknown)}. Available: {', '.join(all_names)}"
        )
        raise ValueError(msg)
    return requested


def single_field(data: SimulationData, requested: str | None) -> str:
    """Resolve a single field name, defaulting to first dynamical field.

    When *requested* is ``None``, returns the first field that has
    momentum data (i.e. a dynamical, time-evolved field) rather than a
    constraint field.  Falls back to the absolute first field if every
    field is a constraint.

    Raises
    ------
    ValueError
        If the requested field name is not found in data.
    """
    all_names = list(data.fields.keys())
    if requested is None:
        # Prefer dynamical fields over constraints.  When multiple
        # dynamical fields exist, pick the one with the largest peak
        # amplitude so the overview shows the most interesting field
        # (e.g. A_2 in EM instead of zero-amplitude A_1).
        dyn = [n for n in all_names if n in set(data.dynamical_fields)]
        if not dyn:
            return all_names[0]
        if len(dyn) == 1:
            return dyn[0]
        best, best_amp = dyn[0], 0.0
        for name in dyn:
            amp = float(np.max(np.abs(data.fields[name])))
            if amp > best_amp:
                best, best_amp = name, amp
        return best
    if requested not in data.fields:
        msg = f"Unknown field '{requested}'. Available: {', '.join(all_names)}"
        raise ValueError(msg)
    return requested


def _spatial_dim(data: SimulationData) -> int:
    """Return the number of spatial dimensions (1, 2, or 3)."""
    first = next(iter(data.fields.values()))
    return first.ndim - 1  # (n_snapshots, *spatial_shape)


def _z_profile_1d(
    snap: np.ndarray, data: SimulationData,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract z-profile at the centre of x, y for a 3D snapshot.

    Returns ``(z_coords, profile)`` where *z_coords* are cell centres
    along the last spatial axis and *profile* is ``snap[ix, iy, :]``.
    """
    ix = snap.shape[0] // 2
    iy = snap.shape[1] // 2
    z_min, z_max = data.grid_bounds[2]
    n_z = snap.shape[2]
    z = np.linspace(
        z_min + data.grid_spacing[2] / 2, z_max - data.grid_spacing[2] / 2, n_z,
    )
    return z, snap[ix, iy, :]


def resolve_time_indices(
    data: SimulationData,
    raw: list[int] | None,
    *,
    n_default: int = 5,
) -> list[int]:
    """Resolve time indices, supporting negative indexing.

    Raises
    ------
    ValueError
        If any resolved index is out of range.
    """
    n = data.n_snapshots
    if raw is not None:
        resolved = [i if i >= 0 else n + i for i in raw]
        for idx in resolved:
            if idx < 0 or idx >= n:
                msg = f"Time index {idx} out of range [0, {n - 1}]"
                raise ValueError(msg)
        return resolved
    # Default: evenly spaced including first and last
    if n <= n_default:
        return list(range(n))
    step = (n - 1) / (n_default - 1)
    return [round(i * step) for i in range(n_default)]


def _apply_cross_section(
    field_data: NDArray[np.float64],
    data: SimulationData,
    cross_section: tuple[str, float],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Slice 2D field data along a cross-section, returning (x_coords, sliced_1d).

    Parameters
    ----------
    field_data
        Array of shape ``(*spatial_shape)`` for a single snapshot.
    data
        Simulation data (for grid bounds and spacing).
    cross_section
        Tuple of ``(axis_name, value)`` e.g. ``("y", 25.0)``.

    Returns
    -------
    x_coords
        1D coordinate array along the remaining axis.
    sliced
        1D field values along the remaining axis.

    Raises
    ------
    ValueError
        If the data is not 2D or the axis name is unknown.
    """
    axis_name, value = cross_section
    bounds = data.grid_bounds
    spacing = data.grid_spacing

    if _spatial_dim(data) != 2:  # noqa: PLR2004
        msg = f"--cross-section requires 2D spatial data, got {_spatial_dim(data)}D"
        raise ValueError(msg)

    # axis_name must be "x" or "y"
    if axis_name == "y":
        # Slice along y → return x-axis data
        y_min, y_max = bounds[1]
        n_y = field_data.shape[1]
        y_coords = np.linspace(y_min + spacing[1] / 2, y_max - spacing[1] / 2, n_y)
        idx = int(np.argmin(np.abs(y_coords - value)))
        sliced = field_data[:, idx]  # shape (n_x,)
        x_min, x_max = bounds[0]
        n_x = field_data.shape[0]
        x_coords = np.linspace(x_min + spacing[0] / 2, x_max - spacing[0] / 2, n_x)
        return x_coords, sliced
    if axis_name == "x":
        # Slice along x → return y-axis data
        x_min, x_max = bounds[0]
        n_x = field_data.shape[0]
        x_coords_full = np.linspace(x_min + spacing[0] / 2, x_max - spacing[0] / 2, n_x)
        idx = int(np.argmin(np.abs(x_coords_full - value)))
        sliced = field_data[idx, :]  # shape (n_y,)
        y_min, y_max = bounds[1]
        n_y = field_data.shape[1]
        remaining = np.linspace(y_min + spacing[1] / 2, y_max - spacing[1] / 2, n_y)
        return remaining, sliced

    msg = f"Unknown cross-section axis '{axis_name}', expected 'x' or 'y'"
    raise ValueError(msg)


# ===========================================================================
# Render functions — one per plot type
# ===========================================================================


def render_heatmap(
    ax: Axes,
    data: SimulationData,
    field: str,
    *,
    cmap: str = "RdBu_r",
) -> None:
    """Render spacetime heatmap φ(x, t) for 1D spatial data.

    Raises
    ------
    ValueError
        If the data is not 1D spatial.
    """
    if _spatial_dim(data) != 1:
        msg = f"Heatmap requires 1D spatial data, got {_spatial_dim(data)}D"
        raise ValueError(msg)

    arr = data.fields[field]  # (n_t, n_x)
    vmax = max(float(np.max(np.abs(arr))), VMAX_FLOOR)

    x_min, x_max = data.grid_bounds[0]
    t_min = float(data.times[0])
    t_max = float(data.times[-1])

    ax.imshow(
        arr,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
        extent=(x_min, x_max, t_min, t_max),
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(field)


def render_snapshot(
    ax: Axes,
    data: SimulationData,
    field: str,
    time_index: int,
    *,
    cmap: str = "RdBu_r",
) -> None:
    """Render field snapshot at a single time — 2D imshow or 1D line plot.

    Raises
    ------
    ValueError
        If time_index is out of range or dimension is unsupported.
    """
    idx = time_index if time_index >= 0 else data.n_snapshots + time_index
    if idx < 0 or idx >= data.n_snapshots:
        msg = f"Time index {time_index} out of range [0, {data.n_snapshots - 1}]"
        raise ValueError(msg)

    snap = data.fields[field][idx]
    t_val = float(data.times[idx])
    dim = _spatial_dim(data)

    if dim == 1:
        x_min, x_max = data.grid_bounds[0]
        x = np.linspace(
            x_min + data.grid_spacing[0] / 2,
            x_max - data.grid_spacing[0] / 2,
            snap.shape[0],
        )
        ax.plot(x, snap, "b-", linewidth=1.5)
        ax.set_xlabel("x")
        ax.set_ylabel(field)
    elif dim == 2:  # noqa: PLR2004
        vmax = max(float(np.max(np.abs(snap))), VMAX_FLOOR)
        (x_min, x_max), (y_min, y_max) = data.grid_bounds[:2]
        im = ax.imshow(
            snap.T,
            origin="lower",
            cmap=cmap,
            vmin=-vmax,
            vmax=vmax,
            extent=(x_min, x_max, y_min, y_max),
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.figure.colorbar(im, ax=ax, label=field)  # type: ignore[union-attr]
    elif dim == 3:  # noqa: PLR2004
        z, z_profile = _z_profile_1d(snap, data)
        ax.plot(z, z_profile, "b-", linewidth=1.5)
        ax.set_xlabel("z")
        ax.set_ylabel(field)
    else:
        msg = f"Snapshot not yet supported for {dim}D data"
        raise ValueError(msg)

    ax.set_title(f"{field} (t={t_val:.2f})")


def render_amplitude(
    ax: Axes,
    data: SimulationData,
    fields: list[str],
    *,
    overlay: str | None = None,
) -> None:
    """Peak amplitude max|field| vs time for one or more fields."""
    times = data.times

    for i, name in enumerate(fields):
        arr = data.fields[name]  # (n_t, *spatial)
        peaks = np.array(
            [float(np.max(np.abs(arr[t]))) for t in range(data.n_snapshots)],
        )
        color = _LINE_COLORS[i % len(_LINE_COLORS)]
        ax.plot(times, peaks, color=color, linewidth=1.5, label=name)

    if overlay is not None:
        _plot_overlay(ax, data, overlay)

    ax.set_xlabel("Time")
    ax.set_ylabel(r"max $|$field$|$")
    ax.set_title("Peak Amplitude")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def _plot_overlay(ax: Axes, data: SimulationData, formula: str) -> None:
    """Evaluate and plot an analytic overlay formula as a dashed reference.

    Raises
    ------
    ValueError
        If the formula evaluation fails.
    """
    from tidal.cli._simulate import FORMULA_NAMESPACE

    t = data.times
    ns = {**FORMULA_NAMESPACE, "t": t}
    try:
        values = eval(formula, {"__builtins__": {}}, ns)  # noqa: S307
    except Exception as exc:
        msg = f"Error evaluating overlay formula '{formula}': {exc}"
        raise ValueError(msg) from exc

    ax.plot(t, values, "k--", linewidth=1.2, alpha=0.7, label=f"ref: {formula}")


def render_energy(
    ax: Axes,
    data: SimulationData,
    fields: list[str],
) -> None:
    """Per-field L² energy (∫|field|² dx) vs time."""
    times = data.times
    vol = data.volume_element

    total = np.zeros(data.n_snapshots)
    for i, name in enumerate(fields):
        arr = data.fields[name]  # (n_t, *spatial)
        energy = np.array(
            [float(np.sum(arr[t] ** 2)) * vol for t in range(data.n_snapshots)],
        )
        color = _LINE_COLORS[i % len(_LINE_COLORS)]
        ax.plot(times, energy, color=color, linewidth=1.5, label=name)
        total += energy

    if len(fields) > 1:
        ax.plot(times, total, "k--", linewidth=1.0, alpha=0.5, label="total")

    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\int |$field$|^2 \, dx$")
    ax.set_title("Field Energy")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def render_profile(
    ax: Axes,
    data: SimulationData,
    field: str,
    time_indices: list[int],
    *,
    cross_section: tuple[str, float] | None = None,
) -> None:
    """Render 1D field profile at multiple times, overlaid with colour gradient.

    Raises
    ------
    ValueError
        If data is 2D without cross_section, or dimension is unsupported.
    """
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    n_lines = len(time_indices)
    colors = cmap(np.linspace(0.2, 0.8, n_lines))

    dim = _spatial_dim(data)
    for i, t_idx in enumerate(time_indices):
        idx = t_idx if t_idx >= 0 else data.n_snapshots + t_idx
        snap = data.fields[field][idx]
        t_val = float(data.times[idx])

        if dim == 1:
            x_min, x_max = data.grid_bounds[0]
            n_x = snap.shape[0]
            x = np.linspace(
                x_min + data.grid_spacing[0] / 2, x_max - data.grid_spacing[0] / 2, n_x,
            )
            ax.plot(
                x,
                snap,
                color=colors[i],
                linewidth=1.2,
                label=f"t={t_val:.1f}",
                alpha=0.8,
            )
        elif dim == 2 and cross_section is not None:  # noqa: PLR2004
            x, sliced = _apply_cross_section(snap, data, cross_section)
            ax.plot(
                x,
                sliced,
                color=colors[i],
                linewidth=1.2,
                label=f"t={t_val:.1f}",
                alpha=0.8,
            )
        elif dim == 3:  # noqa: PLR2004
            z, z_profile = _z_profile_1d(snap, data)
            ax.plot(
                z,
                z_profile,
                color=colors[i],
                linewidth=1.2,
                label=f"t={t_val:.1f}",
                alpha=0.8,
            )
        else:
            msg = (
                "Profile requires 1D, 2D (with --cross-section), or 3D data. "
                f"Got {dim}D data."
            )
            raise ValueError(msg)

    if dim == 3:  # noqa: PLR2004
        ax.set_xlabel("z")
    elif cross_section is None or cross_section[0] == "y":
        ax.set_xlabel("x")
    else:
        ax.set_xlabel("y")
    ax.set_ylabel(field)
    ax.set_title(f"{field} profile")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)


def render_compare(
    ax: Axes,
    data: SimulationData,
    fields: list[str],
    *,
    cross_section: tuple[str, float] | None = None,
) -> None:
    """Compare initial vs final field profiles for one or more fields.

    Raises
    ------
    ValueError
        If data is 2D without cross_section, or dimension is unsupported.
    """
    dim = _spatial_dim(data)

    for i, name in enumerate(fields):
        initial = data.fields[name][0]
        final = data.fields[name][-1]
        color = _LINE_COLORS[i % len(_LINE_COLORS)]

        if dim == 1:
            x_min, x_max = data.grid_bounds[0]
            x = np.linspace(
                x_min + data.grid_spacing[0] / 2,
                x_max - data.grid_spacing[0] / 2,
                initial.shape[0],
            )
            ax.plot(x, initial, color=color, linewidth=1.5, label=f"{name} (t=0)")
            ax.plot(
                x,
                final,
                color=color,
                linewidth=1.5,
                linestyle="--",
                label=f"{name} (t={float(data.times[-1]):.1f})",
            )
        elif dim == 2 and cross_section is not None:  # noqa: PLR2004
            x_i, sliced_i = _apply_cross_section(initial, data, cross_section)
            x_f, sliced_f = _apply_cross_section(final, data, cross_section)
            ax.plot(x_i, sliced_i, color=color, linewidth=1.5, label=f"{name} (t=0)")
            ax.plot(
                x_f,
                sliced_f,
                color=color,
                linewidth=1.5,
                linestyle="--",
                label=f"{name} (t={float(data.times[-1]):.1f})",
            )
        elif dim == 2:  # noqa: PLR2004
            msg = "Compare for 2D data requires --cross-section for line overlay"
            raise ValueError(msg)
        elif dim == 3:  # noqa: PLR2004
            z, init_profile = _z_profile_1d(initial, data)
            ax.plot(z, init_profile, color=color, linewidth=1.5, label=f"{name} (t=0)")
            ax.plot(
                z,
                _z_profile_1d(final, data)[1],
                color=color,
                linewidth=1.5,
                linestyle="--",
                label=f"{name} (t={float(data.times[-1]):.1f})",
            )
        else:
            msg = f"Compare not supported for {dim}D data"
            raise ValueError(msg)

    if dim == 3:  # noqa: PLR2004
        ax.set_xlabel("z")
    elif cross_section is None or cross_section[0] == "y":
        ax.set_xlabel("x")
    else:
        ax.set_xlabel("y")
    ax.set_ylabel("field value")
    ax.set_title("Initial vs Final")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)


def render_hamiltonian(
    ax: Axes,
    data: SimulationData,
    fields: list[str],
) -> None:
    """Hamiltonian energy density decomposition vs time.

    Shows per-field energy, interaction energy, and total Hamiltonian
    energy density over the simulation time.  Uses the canonical
    Hamiltonian (Legendre transform) when available, otherwise falls
    back to the virial formula.
    """
    from tidal.measurement._energy import compute_energy_timeseries

    times, per_field, interaction, total = compute_energy_timeseries(data)

    for i, name in enumerate(fields):
        if name in per_field:
            color = _LINE_COLORS[i % len(_LINE_COLORS)]
            ax.plot(times, per_field[name], color=color, linewidth=1.5, label=name)

    if len(fields) > 1:
        ax.plot(
            times,
            interaction,
            color="gray",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            label="interaction",
        )

    ax.plot(times, total, "k-", linewidth=1.0, alpha=0.5, label="total")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\langle \mathcal{H} \rangle$")
    ax.set_title("Hamiltonian Energy")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def render_conservation(
    ax: Axes,
    data: SimulationData,
    *,
    threshold: float = 1e-3,
) -> None:
    r"""Relative energy drift :math:`\Delta E / E_0` vs time.

    Annotates the plot with PASS/FAIL based on whether the maximum
    relative error stays below *threshold*.  When ``data.dt`` is
    available, the threshold is automatically scaled by the shadow-
    Hamiltonian bound (``10 * dt²``), and the effective value is shown.
    """
    from tidal.measurement._diagnostics import check_energy_conservation

    diag = check_energy_conservation(data, threshold=threshold)

    # Compute effective threshold (may differ from input after dt² scaling)
    effective_threshold = threshold
    if data.dt is not None:
        effective_threshold = max(threshold, 10.0 * data.dt**2)

    ax.plot(diag.times, diag.relative_error, "k-", linewidth=1.0)
    ax.axhline(
        effective_threshold,
        color="r",
        linestyle="--",
        alpha=0.5,
        label=f"threshold ({effective_threshold:.1e})",
    )
    ax.axhline(-effective_threshold, color="r", linestyle="--", alpha=0.5)

    status = "PASS" if diag.is_conserved else "FAIL"
    color = "green" if diag.is_conserved else "red"
    dt_note = f"  (dt={data.dt:.4f})" if data.dt is not None else ""
    ax.annotate(
        f"{status}\nmax |dE/E| = {diag.max_relative_error:.2e}{dt_note}",
        xy=(0.95, 0.95),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=color,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "alpha": 0.8},
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta E / E_0$")
    ax.set_title("Energy Conservation")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)
