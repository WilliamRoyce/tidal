"""``tg simulate`` — Run PDE simulation from a JSON specification."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

    from matplotlib.axes import Axes
    from pde import CartesianGrid, FieldCollection, MemoryStorage

    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

# Default grid shapes per spatial dimension
_DEFAULT_SHAPES: dict[int, list[int]] = {
    1: [64],
    2: [32, 32],
    3: [16, 16, 16],
}

_DEFAULT_BOUND = (0.0, 10.0)
_SPATIAL_DIM_2D = 2

# CFL safety factor for auto-dt computation
_CFL_FACTOR = 0.5

# Visualization defaults
_DPI = 150
_VMAX_FLOOR = 0.01

# Curated namespace for --ic-formula eval().
# Includes np for backward compatibility (e.g. np.exp(...) in formulas)
# plus named math functions for convenience.
_FORMULA_NAMESPACE: dict[str, object] = {
    "np": np,
    "pi": np.pi,
    "e": np.e,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "exp": np.exp,
    "log": np.log,
    "log10": np.log10,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "tanh": np.tanh,
    "cosh": np.cosh,
    "sinh": np.sinh,
    "arctan": np.arctan,
    "arctan2": np.arctan2,
    "heaviside": np.heaviside,
    "where": np.where,
}


@dataclass(frozen=True)
class _PlotContext:
    """Bundled context for plot generation."""

    spec: EquationSystem
    storage: MemoryStorage
    grid: CartesianGrid
    initial_state: FieldCollection
    params: dict[str, float]


def _field_slots(spec: EquationSystem) -> dict[str, int]:
    """Map component names to their field slot indices in the state vector.

    Uses ``spec.state_layout`` which accounts for mixed time-orders:
    second-order fields have [field, momentum] pairs, while first-order
    and constraint fields have only [field].
    """
    return {
        name: idx
        for idx, (name, slot_type) in enumerate(spec.state_layout)
        if slot_type == "field"
    }


def _parse_params(raw: list[str], spec: EquationSystem) -> dict[str, float]:
    """Parse --param KEY=VAL arguments into a dict.

    Also merges default parameters from metadata when not overridden.

    Raises
    ------
    ValueError
        If parameter format is invalid or value is non-numeric.
    """
    params: dict[str, float] = {}

    # Start with metadata defaults
    meta_params: object = spec.metadata.get("parameters", {})
    if isinstance(meta_params, dict):
        meta_dict = cast("dict[str, object]", meta_params)
        for key, val in meta_dict.items():
            try:
                params[key] = float(val)  # type: ignore[arg-type]
            except (ValueError, TypeError):
                print(
                    f"  Warning: metadata parameter '{key}' has non-numeric value {val!r}, skipping",
                    file=sys.stderr,
                )

    # Override with CLI params
    for item in raw:
        if "=" not in item:
            msg = f"Invalid --param format: '{item}'. Expected KEY=VALUE (e.g. --param m2=1.0)"
            raise ValueError(msg)
        key, val_str = item.split("=", 1)
        try:
            params[key.strip()] = float(val_str.strip())
        except ValueError:
            msg = f"Invalid parameter value: '{val_str}' for key '{key}'. Must be a number."
            raise ValueError(msg) from None

    return params


def _parse_grid_shape(raw: str | None, spatial_dim: int) -> list[int]:
    """Parse --grid-shape argument.

    Raises
    ------
    ValueError
        If the number of values doesn't match spatial dimension.
    """
    if raw is None:
        return _DEFAULT_SHAPES.get(spatial_dim, [16] * spatial_dim)

    parts = [s.strip() for s in raw.split(",")]
    if len(parts) == 1:
        # Uniform: "32" → [32, 32, 32]
        n = int(parts[0])
        return [n] * spatial_dim
    if len(parts) == spatial_dim:
        return [int(p) for p in parts]

    msg = (
        f"--grid-shape expects 1 or {spatial_dim} values "
        f"(got {len(parts)}). Example: --grid-shape 32 or --grid-shape 32,32,32"
    )
    raise ValueError(msg)


def _parse_bounds(raw: str | None, spatial_dim: int) -> list[tuple[float, float]]:
    """Parse --bounds argument.

    Raises
    ------
    ValueError
        If the number of values doesn't match spatial dimension or format is invalid.
    """
    if raw is None:
        return [_DEFAULT_BOUND] * spatial_dim

    parts = [s.strip() for s in raw.split(",")]
    if len(parts) == 1:
        lo, hi = _parse_single_bound(parts[0])
        return [(lo, hi)] * spatial_dim
    if len(parts) == spatial_dim:
        return [_parse_single_bound(p) for p in parts]

    msg = (
        f"--bounds expects 1 or {spatial_dim} values "
        f"(got {len(parts)}). Example: --bounds 0:20 or --bounds 0:20,0:10,0:10"
    )
    raise ValueError(msg)


def _parse_single_bound(s: str) -> tuple[float, float]:
    """Parse 'LO:HI' → (float, float).

    Raises
    ------
    ValueError
        If format is not LO:HI.
    """
    if ":" not in s:
        msg = f"Invalid bound format: '{s}'. Expected LO:HI (e.g. 0:20)"
        raise ValueError(msg)
    lo_str, hi_str = s.split(":", 1)
    lo, hi = float(lo_str), float(hi_str)
    if lo >= hi:
        msg = f"Invalid bound: lower ({lo}) must be less than upper ({hi})"
        raise ValueError(msg)
    return lo, hi


_VALID_BC_TYPES = {"periodic", "neumann", "dirichlet"}


def _parse_bc(raw: str | None, *, periodic: bool, spatial_dim: int) -> bool | list[bool]:
    """Parse --bc argument into periodic specification for CartesianGrid.

    Parameters
    ----------
    raw : str | None
        Raw --bc argument (e.g. "neumann,periodic").
    periodic : bool
        Value from --periodic flag (used when --bc is None).
    spatial_dim : int
        Number of spatial dimensions.

    Returns
    -------
    bool | list[bool]
        Single bool or per-axis list for CartesianGrid periodic parameter.

    Raises
    ------
    ValueError
        If BC count doesn't match dimension or BC type is unknown.
    """
    if raw is None:
        return periodic

    bc_list = [b.strip().lower() for b in raw.split(",")]

    if len(bc_list) == 1:
        bc_list *= spatial_dim
    elif len(bc_list) != spatial_dim:
        msg = (
            f"--bc expects 1 or {spatial_dim} values "
            f"(got {len(bc_list)}). Example: --bc neumann,periodic"
        )
        raise ValueError(msg)

    for bc in bc_list:
        if bc not in _VALID_BC_TYPES:
            msg = f"Invalid boundary condition: '{bc}'. Must be one of: {', '.join(sorted(_VALID_BC_TYPES))}"
            raise ValueError(msg)

    return [bc == "periodic" for bc in bc_list]


def _build_grid(
    args: Namespace,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> CartesianGrid:
    """Build a CartesianGrid from CLI arguments and pre-parsed bounds."""
    from pde import CartesianGrid

    shape = _parse_grid_shape(args.grid_shape, spec.spatial_dimension)
    periodic = _parse_bc(args.bc, periodic=args.periodic, spatial_dim=spec.spatial_dimension)

    return CartesianGrid(
        bounds=bounds,
        shape=shape,
        periodic=periodic,
    )


def _apply_formula_ic(
    args: Namespace,
    grid: CartesianGrid,
    spec: EquationSystem,
    component: str,
) -> FieldCollection:
    """Apply formula-based initial condition.

    Raises
    ------
    ValueError
        If --ic-formula expression is not provided.
    """
    from torsion_gertsenshtein.symbolic import create_initial_state

    if args.ic_formula is None:
        msg = "--ic=formula requires --ic-formula=EXPR"
        raise ValueError(msg)

    coords = spec.spatial_coordinates
    namespace = dict(_FORMULA_NAMESPACE)
    for i, name in enumerate(coords):
        namespace[name] = grid.cell_coords[..., i]

    field_arr = eval(args.ic_formula, {"__builtins__": {}}, namespace)  # noqa: S307
    field_arr = np.asarray(field_arr, dtype=float)

    # Broadcast scalar results to grid shape
    if field_arr.shape == ():
        field_arr = np.full(grid.shape, float(field_arr))

    return create_initial_state(grid, spec, field_data={component: field_arr})


def _build_initial_state(
    args: Namespace,
    grid: CartesianGrid,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> FieldCollection:
    """Build initial state from CLI IC arguments.

    Uses ``create_initial_state`` which respects ``state_layout`` for
    mixed time-order systems (constraint + dynamical fields).

    Raises
    ------
    ValueError
        If component name is unknown or IC type is invalid.
    """
    from torsion_gertsenshtein.symbolic import create_initial_state
    from torsion_gertsenshtein.vectorfield.initial_conditions import (
        ComponentGaussianPulse,
        ComponentPlaneWave,
    )

    ic_type = args.ic
    component = args.ic_component or spec.component_names[0]

    if component not in spec.component_names:
        msg = (
            f"Unknown component '{component}'. "
            f"Available: {', '.join(spec.component_names)}"
        )
        raise ValueError(msg)

    if ic_type == "zero":
        if args.ic_component is not None:
            print(
                f"  Note: --ic-component '{args.ic_component}' is ignored for zero IC"
            )
        return create_initial_state(grid, spec)

    if ic_type == "gaussian":
        if args.ic_center is not None:
            center = tuple(float(c) for c in args.ic_center.split(","))
        else:
            center = tuple((lo + hi) / 2.0 for lo, hi in bounds)
        domain_size = min(hi - lo for lo, hi in bounds)
        width = args.ic_width if args.ic_width is not None else domain_size / 10.0

        pulse = ComponentGaussianPulse(
            center=center,
            width=width,
            amplitude=args.ic_amplitude,
            active_components={component: 1.0},
        )
        # Use create_initial_state with field_data for layout-aware IC
        field_arr = pulse.compute_gaussian(grid)
        return create_initial_state(grid, spec, field_data={component: field_arr})

    if ic_type == "plane-wave":
        if args.ic_wavevector is not None:
            kvec = tuple(float(k) for k in args.ic_wavevector.split(","))
        else:
            # Default: one wavelength across domain in x-direction
            lx = bounds[0][1] - bounds[0][0]
            kvec = tuple(
                2.0 * math.pi / lx if i == 0 else 0.0
                for i in range(spec.spatial_dimension)
            )

        wave = ComponentPlaneWave(
            wavevector=kvec,
            amplitude=args.ic_amplitude,
            active_components={component: 1.0},
        )
        # Use create_initial_state for layout-aware IC
        field_arr, momentum_arr = wave.compute_plane_wave(grid)
        return create_initial_state(
            grid,
            spec,
            field_data={component: field_arr},
            momentum_data={component: momentum_arr},
        )

    if ic_type == "formula":
        return _apply_formula_ic(args, grid, spec, component)

    msg = f"Unknown IC type: {ic_type}"
    raise ValueError(msg)


def _infer_output_format(args: Namespace) -> str:
    """Determine output format from --format or file extension."""
    if args.no_plot:
        return "summary"
    if args.output_format is not None:
        return args.output_format
    if args.output is not None:
        ext = Path(args.output).suffix.lower()
        if ext == ".npz":
            return "npz"
        if ext in {".png", ".pdf", ".jpg", ".svg"}:
            return "png"
    return "png"


def _generate_output(args: Namespace, ctx: _PlotContext) -> None:
    """Generate output based on format selection."""
    fmt = _infer_output_format(args)

    final = cast("FieldCollection", ctx.storage[-1])
    times = list(ctx.storage.times)

    # Always print summary
    _print_summary(ctx.spec, ctx.initial_state, final, times, ctx.params)

    if fmt == "summary":
        return

    # Determine output path (default: next to JSON spec file)
    if args.output is not None:
        output_path = Path(args.output)
    else:
        json_file = Path(args.json_path).resolve()
        output_path = json_file.parent / f"{json_file.stem}_output.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "npz":
        _save_npz(output_path, ctx.spec, ctx.storage)
    else:
        _save_plot(output_path, ctx)


def _print_constraint_summary(
    spec: EquationSystem,
    initial: FieldCollection,
    solved: FieldCollection,
    params: dict[str, float],
) -> None:
    """Print constraint-solve summary to stdout."""
    print()
    print("Results (constraint solve):")
    print(f"  Parameters: {params}")
    print()

    slots = _field_slots(spec)
    for name in spec.component_names:
        slot = slots[name]
        init_peak = float(np.max(np.abs(initial[slot].data)))
        solved_peak = float(np.max(np.abs(solved[slot].data)))
        print(f"  {name}: peak {init_peak:.4f} → {solved_peak:.4f}")


def _print_summary(
    spec: EquationSystem,
    initial: FieldCollection,
    final: FieldCollection,
    times: list[float],
    params: dict[str, float],
) -> None:
    """Print simulation summary to stdout."""
    print()
    print("Results:")
    print(f"  Time range: 0.0 → {times[-1]:.2f} ({len(times)} snapshots)")
    print(f"  Parameters: {params}")
    print()

    slots = _field_slots(spec)
    for name in spec.component_names:
        slot = slots[name]
        init_peak = float(np.max(np.abs(initial[slot].data)))
        final_peak = float(np.max(np.abs(final[slot].data)))
        if init_peak > 0:
            ratio = final_peak / init_peak
            print(
                f"  {name}: peak {init_peak:.4f} → {final_peak:.4f} (ratio: {ratio:.4f})"
            )
        else:
            print(f"  {name}: peak {init_peak:.4f} → {final_peak:.4f}")


def _save_npz(
    path: Path,
    spec: EquationSystem,
    storage: MemoryStorage,
) -> None:
    """Save simulation data as .npz file."""
    times = np.array(storage.times)
    field_data: dict[str, np.ndarray] = {"times": times}

    slots = _field_slots(spec)
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        for name in spec.component_names:
            key = f"{name}_t{t_idx}"
            field_data[key] = snapshot[slots[name]].data

    np.savez(str(path), **field_data)  # type: ignore[reportArgumentType]
    print(f"  Saved data to: {path}")


def _save_constraint_output(
    args: Namespace,
    spec: EquationSystem,
    state: FieldCollection,
) -> None:
    """Save constraint-solve output to NPZ if requested."""
    if args.no_plot or args.output is None:
        return
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix == ".npz":
        slots = _field_slots(spec)
        field_data: dict[str, np.ndarray] = {}
        for name in spec.component_names:
            field_data[name] = state[slots[name]].data
        np.savez(str(output_path), **field_data)  # type: ignore[reportArgumentType]
        print(f"  Saved data to: {output_path}")


def _save_plot(path: Path, ctx: _PlotContext) -> None:
    """Generate and save visualization."""
    import matplotlib as mpl

    mpl.use("Agg")

    spatial_dim = ctx.spec.spatial_dimension

    if spatial_dim == 1:
        _plot_1d(path, ctx)
    elif spatial_dim == _SPATIAL_DIM_2D:
        _plot_2d(path, ctx)
    else:
        _plot_3d(path, ctx)


# --- Shared plot helpers ---


def _plot_amplitude_decay(
    ax: Axes,
    ctx: _PlotContext,
    slot: int,
    name: str,
    times: list[float],
) -> None:
    """Plot peak amplitude decay over time (shared by 1D/2D/3D)."""
    peaks = [float(np.max(np.abs(ctx.initial_state[slot].data)))]
    for t_idx in range(len(ctx.storage)):
        snap = cast("FieldCollection", ctx.storage[t_idx])
        peaks.append(float(np.max(np.abs(snap[slot].data))))
    ax.plot([0.0, *times], peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)


def _plot_1d(path: Path, ctx: _PlotContext) -> None:
    """1D visualization: spacetime heatmap + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    times = list(ctx.storage.times)
    x_coords = cast("np.ndarray", ctx.grid.cell_coords[..., 0])
    slots = _field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    # Panel 1: spacetime heatmap for first field
    ax = axes[0]
    heatmap_data: list[np.ndarray] = []
    for t_idx in range(len(ctx.storage)):
        snap = cast("FieldCollection", ctx.storage[t_idx])
        heatmap_data.append(snap[slot].data)

    heatmap = np.array(heatmap_data)
    ax.imshow(
        heatmap,
        aspect="auto",
        origin="lower",
        extent=[float(x_coords[0]), float(x_coords[-1]), 0, times[-1]],
        cmap="RdBu_r",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(f"{name} spacetime evolution")

    # Panel 2: amplitude decay
    _plot_amplitude_decay(axes[1], ctx, slot, name, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(
        f"{ctx.spec.component_names[0]} ({param_str})"
        if ctx.params
        else ctx.spec.component_names[0]
    )
    plt.tight_layout()
    plt.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


def _plot_2d(path: Path, ctx: _PlotContext) -> None:
    """2D visualization: initial + final snapshots + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    times = list(ctx.storage.times)
    final = cast("FieldCollection", ctx.storage[-1])
    slots = _field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    bounds = ctx.grid.axes_bounds

    # Panel 1: initial x-y
    ax = axes[0]
    init_data = ctx.initial_state[slot].data
    vmax = max(float(np.max(np.abs(init_data))), _VMAX_FLOOR)
    ax.imshow(
        init_data.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} (t=0)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 2: final x-y
    ax = axes[1]
    final_data = final[slot].data
    vmax_f = max(float(np.max(np.abs(final_data))), _VMAX_FLOOR)
    ax.imshow(
        final_data.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax_f,
        vmax=vmax_f,
    )
    ax.set_title(f"{name} (t={times[-1]:.1f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 3: amplitude decay
    _plot_amplitude_decay(axes[2], ctx, slot, name, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(f"{name} ({param_str})" if ctx.params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


# --- 3D plot sub-helpers ---


def _plot_z_profile(
    ax: Axes,
    ctx: _PlotContext,
    slot: int,
    name: str,
    times: list[float],
) -> None:
    """Plot z-profile at center of x-y plane."""
    n = ctx.grid.shape[0]
    ic = n // 2
    z_1d = cast("np.ndarray", ctx.grid.cell_coords[ic, ic, :, 2])
    final = cast("FieldCollection", ctx.storage[-1])

    ax.plot(z_1d, ctx.initial_state[slot].data[ic, ic, :], "b-", linewidth=2, label="t=0")
    ax.plot(
        z_1d, final[slot].data[ic, ic, :], "r-", linewidth=2, label=f"t={times[-1]:.1f}"
    )
    ax.set_xlabel("z")
    ax.set_ylabel(name)
    ax.set_title(f"{name} z-profile (x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def _plot_xy_slice(
    ax: Axes,
    ctx: _PlotContext,
    slot: int,
    name: str,
) -> None:
    """Plot x-y slice at z=center."""
    n = ctx.grid.shape[0]
    ic = n // 2
    bounds = ctx.grid.axes_bounds

    init_slice = ctx.initial_state[slot].data[:, :, ic]
    vmax = max(float(np.max(np.abs(init_slice))), _VMAX_FLOOR)
    ax.imshow(
        init_slice.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} x-y (t=0, z=center)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def _plot_component_independence(
    ax: Axes,
    ctx: _PlotContext,
    slots: dict[str, int],
    times: list[float],
) -> None:
    """Plot amplitude of other components over time."""
    colors = ["red", "green", "purple", "orange", "brown"]
    for i in range(1, ctx.spec.n_components):
        comp_name = ctx.spec.component_names[i]
        comp_slot = slots[comp_name]
        comp_peaks = [0.0]
        for t_idx in range(len(ctx.storage)):
            snap = cast("FieldCollection", ctx.storage[t_idx])
            comp_peaks.append(float(np.max(np.abs(snap[comp_slot].data))))
        ax.plot(
            [0.0, *times],
            comp_peaks,
            color=colors[(i - 1) % len(colors)],
            linewidth=2,
            label=comp_name,
        )
    ax.set_xlabel("Time")
    ax.set_ylabel("max |field|")
    ax.set_title("Other components")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))


def _plot_3d(path: Path, ctx: _PlotContext) -> None:
    """3D visualization: z-profile + x-y slice + amplitude decay + component check."""
    import matplotlib.pyplot as plt

    times = list(ctx.storage.times)
    slots = _field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    n_panels = 4 if ctx.spec.n_components > 1 else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

    # Panel 1: z-profile
    _plot_z_profile(axes[0], ctx, slot, name, times)

    # Panel 2: x-y slice initial
    _plot_xy_slice(axes[1], ctx, slot, name)

    # Panel 3: amplitude decay
    _plot_amplitude_decay(axes[2], ctx, slot, name, times)

    # Panel 4: component independence (multi-field only)
    if ctx.spec.n_components > 1:
        _plot_component_independence(axes[3], ctx, slots, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(f"{name} ({param_str})" if ctx.params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


# --- Command entry point ---


def _check_cfl_stability(
    pde: object, dt: float, grid: CartesianGrid
) -> None:
    """Print CFL stability warnings to stderr."""
    warnings_list = cast("list[str]", pde.check_stability(dt, grid))  # type: ignore[attr-defined]
    sys.stderr.writelines(f"  Warning: {w}\n" for w in warnings_list)


def simulate_command(args: Namespace) -> int:
    """Execute the simulate command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system

    json_path = Path(args.json_path)
    if not json_path.exists():
        print(f"Error: file not found: {json_path}", file=sys.stderr)
        return 1

    # Step 1: Load spec
    print("Loading equation specification...")
    spec = load_equation_system(json_path)
    print(
        f"  {spec.n_components} component(s), {spec.dimension}D ({spec.spatial_dimension}+1D)"
    )

    # Step 2: Parse parameters
    params = _parse_params(args.param, spec)
    if params:
        print(f"  Parameters: {params}")

    # Step 3: Build PDE
    print("Building PDE...")
    pde = build_pde_from_json(json_path, parameters=params)

    # Step 4: Create grid
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    grid = _build_grid(args, spec, bounds)
    print(f"  Grid: {'x'.join(str(s) for s in grid.shape)}, bounds: {grid.axes_bounds}")

    # Step 5: Initial conditions
    state = _build_initial_state(args, grid, spec, bounds)
    initial_state = state.copy()
    print(f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}")

    # Constraint-only mode: single constraint solve, no time evolution
    if args.mode == "constraint":
        print("Solving constraints...")
        pde.evolution_rate(state, t=0.0)  # type: ignore[attr-defined]
        print("  Constraint solve complete.")
        _print_constraint_summary(spec, initial_state, state, params)
        _save_constraint_output(args, spec, state)
        return 0

    # Step 6: Determine dt
    from pde import MemoryStorage

    from torsion_gertsenshtein.utils import normalize_solve_result

    dt = args.dt
    if dt is None:
        # Auto dt from grid spacing
        dt = _CFL_FACTOR * min(
            (b[1] - b[0]) / s for b, s in zip(grid.axes_bounds, grid.shape, strict=True)
        )

    _check_cfl_stability(pde, dt, grid)

    # Step 7: Run simulation
    print(
        f"Running simulation (t=0 → {args.t_end}, dt={dt:.4f}, scheme={args.scheme})..."
    )
    storage = MemoryStorage()
    tracker = storage.tracker(
        args.snapshots if args.snapshots is not None else args.t_end / 100.0
    )
    if args.scheme == "scipy":
        # py-pde's ScipySolver is a separate solver class, not a scheme of ExplicitSolver
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]
                state, t_range=args.t_end, dt=dt, solver="scipy", tracker=tracker,
            )
        )
    else:
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]
                state, t_range=args.t_end, dt=dt, scheme=args.scheme, tracker=tracker,
            )
        )
    print(f"  {len(storage)} snapshots stored")

    # Step 8: Output
    ctx = _PlotContext(
        spec=spec,
        storage=storage,
        grid=grid,
        initial_state=initial_state,
        params=params,
    )
    _generate_output(args, ctx)

    return 0
