"""``tg simulate`` — Run PDE simulation from a JSON specification."""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

    from pde import CartesianGrid, FieldCollection, MemoryStorage

    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

# Default grid shapes per spatial dimension
_DEFAULT_SHAPES: dict[int, list[int]] = {
    1: [64],
    2: [32, 32],
    3: [16, 16, 16],
}

_DEFAULT_BOUND = (0.0, 10.0)


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
    """
    params: dict[str, float] = {}

    # Start with metadata defaults
    meta_params = spec.metadata.get("parameters", {})
    if isinstance(meta_params, dict):
        for k, v in meta_params.items():
            try:
                params[k] = float(v)
            except (ValueError, TypeError):
                print(
                    f"  Warning: metadata parameter '{k}' has non-numeric value {v!r}, skipping"
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
    """Parse --grid-shape argument."""
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
    """Parse --bounds argument."""
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
    """Parse 'LO:HI' → (float, float)."""
    if ":" not in s:
        msg = f"Invalid bound format: '{s}'. Expected LO:HI (e.g. 0:20)"
        raise ValueError(msg)
    lo_str, hi_str = s.split(":", 1)
    return float(lo_str), float(hi_str)


_VALID_BC_TYPES = {"periodic", "neumann", "dirichlet"}


def _parse_bc(raw: str | None, periodic: bool, spatial_dim: int) -> bool | list[bool]:
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
    """
    if raw is None:
        return periodic

    bc_list = [b.strip().lower() for b in raw.split(",")]

    if len(bc_list) == 1:
        bc_list = bc_list * spatial_dim
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


def _build_grid(args: Namespace, spec: EquationSystem) -> CartesianGrid:
    """Build a CartesianGrid from CLI arguments."""
    from pde import CartesianGrid

    shape = _parse_grid_shape(args.grid_shape, spec.spatial_dimension)
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    periodic = _parse_bc(args.bc, args.periodic, spec.spatial_dimension)

    return CartesianGrid(
        bounds=bounds,
        shape=shape,
        periodic=periodic,
    )


def _build_initial_state(
    args: Namespace,
    grid: CartesianGrid,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> FieldCollection:
    """Build initial state from CLI IC arguments.

    Uses ``create_initial_state`` which respects ``state_layout`` for
    mixed time-order systems (constraint + dynamical fields).
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
            print(f"  Note: --ic-component '{args.ic_component}' is ignored for zero IC")
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
        field_arr = pulse._compute_gaussian(grid)  # noqa: SLF001
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
        field_arr, momentum_arr = wave._compute_plane_wave(grid)  # noqa: SLF001
        return create_initial_state(
            grid, spec,
            field_data={component: field_arr},
            momentum_data={component: momentum_arr},
        )

    if ic_type == "formula":
        if args.ic_formula is None:
            msg = "--ic=formula requires --ic-formula=EXPR"
            raise ValueError(msg)

        coords = spec.spatial_coordinates
        namespace: dict[str, object] = {"np": np, "pi": np.pi}
        for i, name in enumerate(coords):
            namespace[name] = grid.cell_coords[..., i]

        field_arr = eval(args.ic_formula, {"__builtins__": {}}, namespace)  # noqa: S307
        field_arr = np.asarray(field_arr, dtype=float)

        # Broadcast scalar results to grid shape
        if field_arr.shape == ():
            field_arr = np.full(grid.shape, float(field_arr))

        return create_initial_state(grid, spec, field_data={component: field_arr})

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


def _generate_output(
    args: Namespace,
    spec: EquationSystem,
    storage: MemoryStorage,
    grid: CartesianGrid,
    initial_state: FieldCollection,
    params: dict[str, float],
) -> None:
    """Generate output based on format selection."""
    fmt = _infer_output_format(args)

    final = cast("FieldCollection", storage[-1])
    times = list(storage.times)

    # Always print summary
    _print_summary(spec, initial_state, final, times, params)

    if fmt == "summary":
        return

    # Determine output path
    if args.output is not None:
        output_path = Path(args.output)
    else:
        stem = Path(args.json_path).stem
        output_path = Path("outputs") / f"{stem}_output.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "npz":
        _save_npz(output_path, spec, storage, grid)
    else:
        _save_plot(output_path, spec, storage, grid, initial_state, params)


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
    grid: CartesianGrid,
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

    np.savez(path, **field_data)
    print(f"  Saved data to: {path}")


def _save_plot(
    path: Path,
    spec: EquationSystem,
    storage: MemoryStorage,
    grid: CartesianGrid,
    initial_state: FieldCollection,
    params: dict[str, float],
) -> None:
    """Generate and save visualization."""
    import matplotlib as mpl

    mpl.use("Agg")

    spatial_dim = spec.spatial_dimension

    if spatial_dim == 1:
        _plot_1d(path, spec, storage, grid, initial_state, params)
    elif spatial_dim == 2:
        _plot_2d(path, spec, storage, grid, initial_state, params)
    else:
        _plot_3d(path, spec, storage, grid, initial_state, params)


def _plot_1d(
    path: Path,
    spec: EquationSystem,
    storage: MemoryStorage,
    grid: CartesianGrid,
    initial_state: FieldCollection,
    params: dict[str, float],
) -> None:
    """1D visualization: spacetime heatmap + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    times = list(storage.times)
    x_coords = cast("np.ndarray", grid.cell_coords[..., 0])
    slots = _field_slots(spec)
    name = spec.component_names[0]
    slot = slots[name]

    # Panel 1: spacetime heatmap for first field
    ax = axes[0]
    field_data = []
    for t_idx in range(len(storage)):
        snap = cast("FieldCollection", storage[t_idx])
        field_data.append(snap[slot].data)

    heatmap = np.array(field_data)
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
    ax = axes[1]
    peaks = [float(np.max(np.abs(initial_state[slot].data)))]
    for t_idx in range(len(storage)):
        snap = cast("FieldCollection", storage[t_idx])
        peaks.append(float(np.max(np.abs(snap[slot].data))))
    ax.plot([0.0, *times], peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    param_str = ", ".join(f"{k}={v}" for k, v in params.items())
    fig.suptitle(
        f"{spec.component_names[0]} ({param_str})"
        if params
        else spec.component_names[0]
    )
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


def _plot_2d(
    path: Path,
    spec: EquationSystem,
    storage: MemoryStorage,
    grid: CartesianGrid,
    initial_state: FieldCollection,
    params: dict[str, float],
) -> None:
    """2D visualization: initial + final snapshots + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    times = list(storage.times)
    final = cast("FieldCollection", storage[-1])
    slots = _field_slots(spec)
    name = spec.component_names[0]
    slot = slots[name]

    bounds = grid.axes_bounds

    # Panel 1: initial x-y
    ax = axes[0]
    init_data = initial_state[slot].data
    vmax = max(float(np.max(np.abs(init_data))), 0.01)
    ax.imshow(
        init_data.T,
        origin="lower",
        extent=[bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]],
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
    vmax_f = max(float(np.max(np.abs(final_data))), 0.01)
    ax.imshow(
        final_data.T,
        origin="lower",
        extent=[bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]],
        cmap="RdBu_r",
        vmin=-vmax_f,
        vmax=vmax_f,
    )
    ax.set_title(f"{name} (t={times[-1]:.1f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 3: amplitude decay
    ax = axes[2]
    peaks = [float(np.max(np.abs(initial_state[slot].data)))]
    for t_idx in range(len(storage)):
        snap = cast("FieldCollection", storage[t_idx])
        peaks.append(float(np.max(np.abs(snap[slot].data))))
    ax.plot([0.0, *times], peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    param_str = ", ".join(f"{k}={v}" for k, v in params.items())
    fig.suptitle(f"{name} ({param_str})" if params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


def _plot_3d(
    path: Path,
    spec: EquationSystem,
    storage: MemoryStorage,
    grid: CartesianGrid,
    initial_state: FieldCollection,
    params: dict[str, float],
) -> None:
    """3D visualization: z-profile + x-y slice + amplitude decay + component check."""
    import matplotlib.pyplot as plt

    n = grid.shape[0]
    ic = n // 2  # center indices
    times = list(storage.times)
    final = cast("FieldCollection", storage[-1])
    slots = _field_slots(spec)
    name = spec.component_names[0]
    slot = slots[name]

    z_1d = cast("np.ndarray", grid.cell_coords[ic, ic, :, 2])
    bounds = grid.axes_bounds

    n_panels = 4 if spec.n_components > 1 else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

    # Panel 1: z-profile
    ax = axes[0]
    ax.plot(z_1d, initial_state[slot].data[ic, ic, :], "b-", linewidth=2, label="t=0")
    ax.plot(
        z_1d, final[slot].data[ic, ic, :], "r-", linewidth=2, label=f"t={times[-1]:.1f}"
    )
    ax.set_xlabel("z")
    ax.set_ylabel(name)
    ax.set_title(f"{name} z-profile (x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # Panel 2: x-y slice initial
    ax = axes[1]
    init_slice = initial_state[slot].data[:, :, ic]
    vmax = max(float(np.max(np.abs(init_slice))), 0.01)
    ax.imshow(
        init_slice.T,
        origin="lower",
        extent=[bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} x-y (t=0, z=center)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 3: amplitude decay
    ax = axes[2]
    peaks = [float(np.max(np.abs(initial_state[slot].data)))]
    for t_idx in range(len(storage)):
        snap = cast("FieldCollection", storage[t_idx])
        peaks.append(float(np.max(np.abs(snap[slot].data))))
    ax.plot([0.0, *times], peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    # Panel 4: component independence (multi-field only)
    if spec.n_components > 1:
        ax = axes[3]
        colors = ["red", "green", "purple", "orange", "brown"]
        for i in range(1, spec.n_components):
            comp_name = spec.component_names[i]
            comp_slot = slots[comp_name]
            comp_peaks = [0.0]
            for t_idx in range(len(storage)):
                snap = cast("FieldCollection", storage[t_idx])
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

    param_str = ", ".join(f"{k}={v}" for k, v in params.items())
    fig.suptitle(f"{name} ({param_str})" if params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


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
    from pde import MemoryStorage

    from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system
    from torsion_gertsenshtein.utils import normalize_solve_result

    json_path = Path(args.json_path)
    if not json_path.exists():
        print(f"Error: file not found: {json_path}")
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
    grid = _build_grid(args, spec)
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    print(f"  Grid: {'x'.join(str(s) for s in grid.shape)}, bounds: {grid.axes_bounds}")

    # Step 5: Initial conditions
    state = _build_initial_state(args, grid, spec, bounds)
    initial_state = state.copy()
    print(f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}")

    # Constraint-only mode: single constraint solve, no time evolution
    if args.mode == "constraint":
        print("Solving constraints...")
        pde.evolution_rate(state, t=0.0)
        print("  Constraint solve complete.")

        # Print result summary (no time evolution → no storage/snapshots)
        _print_constraint_summary(spec, initial_state, state, params)

        if not args.no_plot and args.output is not None:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            if output_path.suffix == ".npz":
                slots = _field_slots(spec)
                field_data: dict[str, np.ndarray] = {}
                for name in spec.component_names:
                    field_data[name] = state[slots[name]].data
                np.savez(output_path, **field_data)
                print(f"  Saved data to: {output_path}")

        return 0

    # Step 6: Determine dt
    dt = args.dt
    if dt is None:
        # Auto dt from grid spacing
        dx = min(
            (b[1] - b[0]) / s for b, s in zip(grid.axes_bounds, grid.shape, strict=True)
        )
        dt = 0.5 * dx  # CFL-safe default

    # CFL stability check
    stability_warnings = pde.check_stability(dt, grid)
    if stability_warnings:
        import sys

        for warning in stability_warnings:
            print(f"  Warning: {warning}", file=sys.stderr)

    snapshot_interval = (
        args.snapshots if args.snapshots is not None else args.t_end / 20.0
    )

    # Step 7: Run simulation
    print(
        f"Running simulation (t=0 → {args.t_end}, dt={dt:.4f}, scheme={args.scheme})..."
    )
    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=args.t_end,
        dt=dt,
        scheme=args.scheme,
        tracker=storage.tracker(snapshot_interval),
    )
    normalize_solve_result(result)
    print(f"  {len(storage)} snapshots stored")

    # Step 8: Output
    _generate_output(args, spec, storage, grid, initial_state, params)

    return 0
