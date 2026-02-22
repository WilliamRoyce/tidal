"""``tidal simulate`` — Run PDE simulation from a JSON specification."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Callable

    from pde import CartesianGrid, FieldCollection, MemoryStorage

    from tidal.measurement._io import SimulationData
    from tidal.measurement._writer import SnapshotWriter
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem

# Default grid shapes per spatial dimension
_DEFAULT_SHAPES: dict[int, list[int]] = {
    1: [64],
    2: [32, 32],
    3: [16, 16, 16],
}

_DEFAULT_BOUND = (0.0, 10.0)
SPATIAL_DIM_2D = 2

# CFL safety factor for auto-dt computation
_CFL_FACTOR = 0.5

# Threshold for zero-evolution diagnostic (effectively machine epsilon)
_ZERO_RATE_THRESHOLD = 1e-14

# Visualization defaults
DPI = 150
VMAX_FLOOR = 0.01

# Curated namespace for --ic-formula eval().
# Includes np for backward compatibility (e.g. np.exp(...) in formulas)
# plus named math functions for convenience.
FORMULA_NAMESPACE: dict[str, object] = {
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
class PlotContext:
    """Bundled context for plot generation."""

    spec: EquationSystem
    storage: MemoryStorage
    grid: CartesianGrid
    initial_state: FieldCollection
    params: dict[str, float]

    def to_simulation_data(self) -> SimulationData:
        """Convert to a :class:`~tidal.measurement.SimulationData` for measurement."""
        from tidal.measurement import SimulationData

        return SimulationData.from_storage(
            self.storage,
            self.spec,
            self.grid,
            self.params,
        )


def field_slots(spec: EquationSystem) -> dict[str, int]:
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


def _parse_params(raw: list[str], spec: EquationSystem) -> dict[str, float]:  # noqa: C901
    """Parse --param KEY=VAL arguments into a dict.

    Also merges default parameters from metadata when not overridden.
    Warns on CLI parameters not found in the equation spec.

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
    cli_keys: set[str] = set()
    for item in raw:
        if "=" not in item:
            msg = f"Invalid --param format: '{item}'. Expected KEY=VALUE (e.g. --param m2=1.0)"
            raise ValueError(msg)
        key, val_str = item.split("=", 1)
        key = key.strip()
        cli_keys.add(key)
        try:
            params[key] = float(val_str.strip())
        except ValueError:
            msg = f"Invalid parameter value: '{val_str}' for key '{key}'. Must be a number."
            raise ValueError(msg) from None

    # Warn on unknown CLI params
    if cli_keys:
        from tidal.cli._inspect import discover_parameters

        try:
            known: set[str] = set(discover_parameters(spec).keys())
        except TypeError:
            known = set()
        if isinstance(meta_params, dict):
            known |= set(cast("dict[str, object]", meta_params).keys())
        for key in sorted(cli_keys - known):
            print(
                f"  Warning: parameter '{key}' not found in equation spec. Possible typo?",
                file=sys.stderr,
            )

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


_VALID_BC_TYPES = {"periodic", "neumann"}

# Native path supports Dirichlet (our apply_operator handles it; py-pde doesn't)
_NATIVE_BC_TYPES = frozenset({"periodic", "neumann", "dirichlet"})


def _parse_bc(
    raw: str | None, *, periodic: bool, spatial_dim: int
) -> bool | list[bool]:
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
        if bc == "dirichlet":
            msg = (
                "Dirichlet boundary conditions are not supported by py-pde's CartesianGrid. "
                "Use 'neumann' (zero-flux) or 'periodic' instead."
            )
            raise ValueError(msg)
        if bc not in _VALID_BC_TYPES:
            msg = f"Invalid boundary condition: '{bc}'. Must be one of: {', '.join(sorted(_VALID_BC_TYPES))}"
            raise ValueError(msg)

    return [bc == "periodic" for bc in bc_list]


def _parse_periodic(
    bc_str: str | None,
    *,
    periodic: bool,
    spatial_dim: int,
) -> tuple[bool, ...]:
    """Parse boundary spec into periodic flags for GridInfo (native path).

    Unlike ``_parse_bc`` (which returns ``bool | list[bool]`` for
    CartesianGrid), this always returns a fixed-length tuple and accepts
    ``"dirichlet"`` as a valid BC type.

    Parameters
    ----------
    bc_str : str | None
        Raw ``--bc`` argument (e.g. ``"neumann,periodic"``).
    periodic : bool
        Value from ``--periodic`` flag (used when ``--bc`` is None).
    spatial_dim : int
        Number of spatial dimensions.

    Returns
    -------
    tuple[bool, ...]
        Per-axis periodic flags, length ``spatial_dim``.

    Raises
    ------
    ValueError
        If BC count doesn't match dimension or BC type is unknown.
    """
    if bc_str is None:
        return tuple(periodic for _ in range(spatial_dim))

    bc_list = [b.strip().lower() for b in bc_str.split(",")]

    if len(bc_list) == 1:
        bc_list *= spatial_dim
    elif len(bc_list) != spatial_dim:
        msg = (
            f"--bc expects 1 or {spatial_dim} values "
            f"(got {len(bc_list)}). Example: --bc neumann,periodic"
        )
        raise ValueError(msg)

    for bc in bc_list:
        if bc not in _NATIVE_BC_TYPES:
            msg = (
                f"Invalid boundary condition: '{bc}'. "
                f"Must be one of: {', '.join(sorted(_NATIVE_BC_TYPES))}"
            )
            raise ValueError(msg)

    return tuple(bc == "periodic" for bc in bc_list)


def _build_grid_info(
    args: Namespace,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> GridInfo:
    """Build GridInfo directly from CLI arguments (no py-pde).

    This is the native-path equivalent of ``_build_grid()``.
    """
    from tidal.solver.grid import GridInfo

    shape = _parse_grid_shape(args.grid_shape, spec.spatial_dimension)
    periodic = _parse_periodic(
        args.bc, periodic=args.periodic, spatial_dim=spec.spatial_dimension
    )
    return GridInfo(
        bounds=tuple(bounds),
        shape=tuple(shape),
        periodic=periodic,
    )


def _build_grid(
    args: Namespace,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> CartesianGrid:
    """Build a CartesianGrid from CLI arguments and pre-parsed bounds."""
    from pde import CartesianGrid

    shape = _parse_grid_shape(args.grid_shape, spec.spatial_dimension)
    periodic = _parse_bc(
        args.bc, periodic=args.periodic, spatial_dim=spec.spatial_dimension
    )

    return CartesianGrid(
        bounds=bounds,
        shape=shape,
        periodic=periodic,
    )


def _validate_formula_ast(expr: str, allowed_names: set[str]) -> None:
    """Validate a formula expression using AST analysis.

    Only allows safe math constructs: literals, names from the allowed set,
    binary/unary ops, comparisons, function calls (by name only), subscripts,
    and ternary expressions. Rejects attribute access, imports, assignments,
    lambda, comprehensions, and other potentially unsafe constructs.

    Raises
    ------
    ValueError
        If the expression contains disallowed names or syntax errors.
    TypeError
        If the expression contains disallowed AST node types.
    """
    import ast

    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        msg = f"Invalid formula syntax: {exc}"
        raise ValueError(msg) from exc

    safe_nodes = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Compare,
        ast.BoolOp,
        ast.IfExp,
        ast.Call,
        ast.Subscript,
        ast.Slice,
        ast.Tuple,
        ast.List,
        # Operators
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.FloorDiv,
        ast.Mod,
        ast.Pow,
        ast.USub,
        ast.UAdd,
        ast.Not,
        ast.And,
        ast.Or,
        # Comparisons
        ast.Eq,
        ast.NotEq,
        ast.Lt,
        ast.LtE,
        ast.Gt,
        ast.GtE,
        # Literals
        ast.Constant,
        ast.Starred,
        # Context nodes (internal AST markers)
        ast.Load,
        ast.Store,
        ast.Del,
    )

    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            if node.id not in allowed_names:
                msg = (
                    f"Disallowed name '{node.id}' in formula. "
                    f"Allowed: {', '.join(sorted(allowed_names))}"
                )
                raise ValueError(msg)
        elif isinstance(node, ast.Attribute):
            # Allow np.something (one level only)
            if isinstance(node.value, ast.Name) and node.value.id == "np":
                continue
            msg = f"Attribute access not allowed in formula: '.{node.attr}'"
            raise ValueError(msg)
        elif not isinstance(node, safe_nodes):
            msg = f"Disallowed construct in formula: {type(node).__name__}"
            raise TypeError(msg)


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
        If --ic-formula expression is not provided or contains unsafe constructs.
    """
    from tidal.symbolic import create_initial_state

    if args.ic_formula is None:
        msg = "--ic=formula requires --ic-formula=EXPR"
        raise ValueError(msg)

    coords = spec.spatial_coordinates
    namespace = dict(FORMULA_NAMESPACE)
    for i, name in enumerate(coords):
        namespace[name] = grid.cell_coords[..., i]

    allowed_names = set(namespace.keys())
    _validate_formula_ast(args.ic_formula, allowed_names)

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
    from tidal.symbolic import create_initial_state
    from tidal.vectorfield.initial_conditions import (
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
            if len(center) != spec.spatial_dimension:
                msg = (
                    f"--ic-center has {len(center)} values but spatial dimension is "
                    f"{spec.spatial_dimension}. Expected {spec.spatial_dimension} comma-separated values."
                )
                raise ValueError(msg)
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


# --- Native initial condition helpers (no py-pde) ---


def _gaussian_y0(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
    component: str,
) -> np.ndarray:
    """Compute Gaussian IC as flat state vector (native path).

    Raises
    ------
    ValueError
        If ``--ic-center`` dimension count doesn't match spatial dimension.
    """
    from tidal.solver.fields import FieldSet
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)

    if args.ic_center is not None:
        center = tuple(float(c) for c in args.ic_center.split(","))
        if len(center) != spec.spatial_dimension:
            msg = (
                f"--ic-center has {len(center)} values but spatial dimension is "
                f"{spec.spatial_dimension}. Expected {spec.spatial_dimension} "
                f"comma-separated values."
            )
            raise ValueError(msg)
    else:
        center = tuple((lo + hi) / 2.0 for lo, hi in bounds)

    domain_size = min(hi - lo for lo, hi in bounds)
    width = args.ic_width if args.ic_width is not None else domain_size / 10.0

    coords = grid_info.cell_coords  # (*grid_shape, ndim)
    dist_sq = np.zeros(grid_info.shape, dtype=np.float64)
    for dim in range(grid_info.ndim):
        if dim < len(center):
            dist_sq += (coords[..., dim] - center[dim]) ** 2

    field_arr = args.ic_amplitude * np.exp(-dist_sq / (2 * width**2))
    return FieldSet.from_dict(layout, grid_info.shape, {component: field_arr}).flat.copy()


def _plane_wave_y0(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
    component: str,
) -> np.ndarray:
    """Compute plane-wave IC as flat state vector (native path)."""
    from tidal.solver.fields import FieldSet
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)

    if args.ic_wavevector is not None:
        kvec = tuple(float(k) for k in args.ic_wavevector.split(","))
    else:
        lx = bounds[0][1] - bounds[0][0]
        kvec = tuple(
            2.0 * math.pi / lx if i == 0 else 0.0
            for i in range(spec.spatial_dimension)
        )

    coords = grid_info.cell_coords
    k_dot_x = np.zeros(grid_info.shape, dtype=np.float64)
    for dim in range(min(grid_info.ndim, len(kvec))):
        k_dot_x += kvec[dim] * coords[..., dim]

    k_mag = float(np.sqrt(sum(k**2 for k in kvec)))
    amplitude = args.ic_amplitude

    field_arr = amplitude * np.cos(k_dot_x)
    momentum_arr = -amplitude * k_mag * np.sin(k_dot_x)

    slot_data: dict[str, np.ndarray] = {component: field_arr}
    mom_name = f"pi_{component}"
    if mom_name in FieldSet(layout, grid_info.shape):
        slot_data[mom_name] = momentum_arr

    return FieldSet.from_dict(layout, grid_info.shape, slot_data).flat.copy()


def _formula_y0(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    component: str,
) -> np.ndarray:
    """Compute formula-based IC as flat state vector (native path).

    Raises
    ------
    ValueError
        If ``--ic-formula`` is not provided or contains unsafe constructs.
    """
    from tidal.solver.fields import FieldSet
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)

    if args.ic_formula is None:
        msg = "--ic=formula requires --ic-formula=EXPR"
        raise ValueError(msg)

    namespace = dict(FORMULA_NAMESPACE)
    for i, name in enumerate(spec.spatial_coordinates):
        namespace[name] = grid_info.cell_coords[..., i]

    allowed_names = set(namespace.keys())
    _validate_formula_ast(args.ic_formula, allowed_names)

    field_arr = eval(args.ic_formula, {"__builtins__": {}}, namespace)  # noqa: S307
    field_arr = np.asarray(field_arr, dtype=float)

    if field_arr.shape == ():
        field_arr = np.full(grid_info.shape, float(field_arr))

    return FieldSet.from_dict(layout, grid_info.shape, {component: field_arr}).flat.copy()


def _build_initial_y0(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
) -> np.ndarray:
    """Build initial state as flat numpy vector (native path, no py-pde).

    This is the native-path equivalent of ``_build_initial_state()``.
    Uses ``GridInfo.cell_coords`` for coordinate arrays and packs results
    via ``FieldSet.from_dict`` + ``StateLayout``.

    Returns
    -------
    np.ndarray
        Flat state vector of length ``StateLayout.total_size``.

    Raises
    ------
    ValueError
        If component name is unknown or IC type is invalid.
    """
    from tidal.solver.fields import FieldSet
    from tidal.solver.state import StateLayout

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
        layout = StateLayout.from_spec(spec, grid_info.num_points)
        return FieldSet.zeros(layout, grid_info.shape).flat.copy()

    if ic_type == "gaussian":
        return _gaussian_y0(args, spec, grid_info, bounds, component)

    if ic_type == "plane-wave":
        return _plane_wave_y0(args, spec, grid_info, bounds, component)

    if ic_type == "formula":
        return _formula_y0(args, spec, grid_info, component)

    msg = f"Unknown IC type: {ic_type}"
    raise ValueError(msg)


# --- Native output pipeline (no py-pde) ---


def _print_summary_native(sim_data: SimulationData) -> None:
    """Print simulation summary using SimulationData (no py-pde types)."""
    times = sim_data.times
    print()
    print("Results:")
    print(f"  Time range: {float(times[0]):.2f} → {float(times[-1]):.2f} ({len(times)} snapshots)")
    print(f"  Parameters: {sim_data.parameters}")
    print()

    for name in sim_data.spec.component_names:
        if name not in sim_data.fields:
            continue
        init_peak = float(np.max(np.abs(sim_data.fields[name][0])))
        final_peak = float(np.max(np.abs(sim_data.fields[name][-1])))
        if init_peak > 0:
            ratio = final_peak / init_peak
            print(
                f"  {name}: peak {init_peak:.4f} → {final_peak:.4f} (ratio: {ratio:.4f})"
            )
        else:
            print(f"  {name}: peak {init_peak:.4f} → {final_peak:.4f}")


def _generate_output_native(
    args: Namespace,
    sim_data: SimulationData,
    grid_info: GridInfo,
) -> None:
    """Generate output for the native solver path (no py-pde types)."""
    fmt = _infer_output_format(args)

    if fmt in {"summary", "directory"}:
        if sim_data.n_snapshots > 0:
            _print_summary_native(sim_data)
        return

    _print_summary_native(sim_data)

    if args.output is not None:
        output_path = Path(args.output)
    else:
        json_file = Path(args.json_path).resolve()
        output_path = json_file.parent / f"{json_file.stem}_output.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    from tidal.cli._plot import save_plot_native

    save_plot_native(output_path, sim_data, grid_info)


def _infer_output_format(args: Namespace) -> str:
    """Determine output format from --format or file extension.

    Raises
    ------
    ValueError
        If ``--output`` has a ``.npz`` extension (no longer supported).
    """
    if args.no_plot:
        return "summary"
    if args.output_format is not None:
        return args.output_format
    if args.output is not None:
        ext = Path(args.output).suffix.lower()
        if ext == ".npz":
            msg = (
                "NPZ format is no longer supported. "
                "Use a directory path (no extension) for disk-backed output."
            )
            raise ValueError(msg)
        if ext in {".png", ".pdf", ".jpg", ".svg"}:
            return ext.lstrip(".")
        # No extension → directory format (disk-backed streaming)
        if not ext:
            return "directory"
    return "png"


def _generate_output(args: Namespace, ctx: PlotContext) -> None:
    """Generate output based on format selection.

    For ``"directory"`` format, this is a no-op — data was already streamed
    to disk by ``SnapshotWriter`` before this function is called.
    """
    fmt = _infer_output_format(args)

    if fmt in {"summary", "directory"}:
        # Directory/summary: data already on disk; print summary if storage
        # captured at least one snapshot (may be empty if solver overshot).
        if len(ctx.storage) > 0:
            final = cast("FieldCollection", ctx.storage[-1])
            times = list(ctx.storage.times)
            _print_summary(ctx.spec, ctx.initial_state, final, times, ctx.params)
        return

    final = cast("FieldCollection", ctx.storage[-1])
    times = list(ctx.storage.times)

    # Always print summary
    _print_summary(ctx.spec, ctx.initial_state, final, times, ctx.params)

    # Determine output path (default: next to JSON spec file)
    if args.output is not None:
        output_path = Path(args.output)
    else:
        json_file = Path(args.json_path).resolve()
        output_path = json_file.parent / f"{json_file.stem}_output.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    from tidal.cli._plot import save_plot

    save_plot(output_path, ctx)


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

    slots = field_slots(spec)
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

    slots = field_slots(spec)
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


def _save_constraint_output(
    args: Namespace,
    spec: EquationSystem,
    state: FieldCollection,
) -> None:
    """Save constraint-solve output to a directory if requested.

    Raises
    ------
    ValueError
        If ``--output`` has a ``.npz`` extension (no longer supported).
    """
    if args.no_plot or args.output is None:
        return
    output_path = Path(args.output)
    if output_path.suffix == ".npz":
        msg = (
            "NPZ format is no longer supported. "
            "Use a directory path (no extension) for output."
        )
        raise ValueError(msg)
    output_path.mkdir(parents=True, exist_ok=True)
    slots = field_slots(spec)
    for name in spec.component_names:
        np.save(str(output_path / f"{name}.npy"), state[slots[name]].data)
    print(f"  Saved data to: {output_path}")


# --- BC string derivation (used by both native and legacy paths) ---


def _bc_from_args(
    args: Namespace, periodic_flags: tuple[bool, ...]
) -> str | tuple[str, ...]:
    """Determine BC string for spatial operators from CLI args."""
    if args.bc:
        bc_parts = [b.strip().lower() for b in args.bc.split(",")]
        return tuple(bc_parts) if len(bc_parts) > 1 else bc_parts[0]
    if all(periodic_flags):
        return "periodic"
    if not any(periodic_flags):
        return "neumann"
    return tuple("periodic" if p else "neumann" for p in periodic_flags)


# --- Native simulation path (no py-pde) ---


def _warn_zero_evolution_native(
    spec: EquationSystem,
    grid_info: GridInfo,
    y0: np.ndarray,
    params: dict[str, float],
    bc: str | tuple[str, ...] | None,
) -> None:
    """Warn if all RHS rates are zero at t=0 (native path)."""
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.fields import FieldSet
    from tidal.solver.rhs import RHSEvaluator
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid_info, params)
    rhs_eval = RHSEvaluator(spec, grid_info, coeff_eval, bc=bc)
    fieldset = FieldSet.from_flat(layout, grid_info.shape, y0)

    max_rate = 0.0
    for eq_idx in range(len(spec.equations)):
        rhs = rhs_eval.evaluate(eq_idx, fieldset, t=0.0)
        max_rate = max(max_rate, float(np.max(np.abs(rhs))))

    if max_rate < _ZERO_RATE_THRESHOLD:
        print(
            "  Warning: all evolution rates are zero at t=0. "
            "The initial condition may be a static configuration. "
            "For gauge fields, try --ic plane-wave to provide "
            "non-zero conjugate momentum.",
            file=sys.stderr,
        )


def _setup_disk_writer_native(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    params: dict[str, float],
    snapshot_interval: float,
) -> tuple[SnapshotWriter, Callable[[float, np.ndarray], None]]:
    """Set up disk-backed SnapshotWriter using StateLayout (no py-pde).

    Returns (writer, snapshot_callback).
    """
    from tidal.measurement._writer import SnapshotWriter, compute_snapshot_count
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    output_dir = Path(args.output) if args.output else Path("output")
    n_snaps = compute_snapshot_count(args.t_end, snapshot_interval)

    field_names = [s.name for s in layout.slots if s.kind != "momentum"]
    momentum_names = [s.field_name for s in layout.slots if s.kind == "momentum"]

    writer = SnapshotWriter(
        output_dir=output_dir,
        field_names=field_names,
        momentum_names=momentum_names,
        grid_shape=grid_info.shape,
        n_snapshots=n_snaps,
        grid_spacing=tuple(float(d) for d in grid_info.dx),
        grid_bounds=grid_info.bounds,
        periodic=grid_info.periodic,
        parameters=params,
        spec_path=Path(args.json_path),
    )

    # Build slot index maps from layout
    n_pts = grid_info.num_points
    shape = grid_info.shape
    field_set = set(field_names)
    momentum_set = set(momentum_names)

    field_slots_map: dict[str, int] = {}
    momentum_slots_map: dict[str, int] = {}
    for i, slot in enumerate(layout.slots):
        if slot.kind == "momentum":
            momentum_slots_map[slot.field_name] = i
        elif slot.name in field_set:
            field_slots_map[slot.name] = i

    def _disk_callback(t: float, y_flat: np.ndarray) -> None:
        fields_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in field_slots_map.items()
        }
        moms_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in momentum_slots_map.items()
            if name in momentum_set
        }
        writer.append(t, fields_d, moms_d)

    return writer, _disk_callback


def _simulate_native(
    args: Namespace,
    spec: EquationSystem,
    params: dict[str, float],
) -> int:
    """Run simulation via native TIDAL solver (no py-pde).

    Self-contained flow: GridInfo -> IC -> solve -> SimulationData -> output.
    Handles both IDA and leapfrog schemes.
    """
    from tidal.measurement._io import SimulationData

    # Progress printer
    def _noop(*_a: object, **_kw: object) -> None:
        pass

    log = _noop if args.quiet else print

    # 1. Grid
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    grid_info = _build_grid_info(args, spec, bounds)
    log(f"  Grid: {'x'.join(str(s) for s in grid_info.shape)}, bounds: {grid_info.bounds}")

    # 2. BC
    bc = _bc_from_args(args, grid_info.periodic)

    # 3. Initial conditions
    y0 = _build_initial_y0(args, spec, grid_info, bounds)
    log(f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}")

    # 4. Diagnostics
    _validate_solver_params(args)
    if args.mode != "constraint":
        _warn_zero_evolution_native(spec, grid_info, y0, params, bc)

    # 5. Snapshot configuration
    snapshot_interval = (
        args.snapshots if args.snapshots is not None else args.t_end / 100.0
    )

    # 6. Disk writer (if directory output)
    fmt = _infer_output_format(args)
    writer: SnapshotWriter | None = None
    snapshot_cb: Callable[[float, np.ndarray], None] | None = None

    if fmt == "directory":
        writer, snapshot_cb = _setup_disk_writer_native(
            args, spec, grid_info, params, snapshot_interval,
        )

    # 7. Solve
    if args.scheme == "ida":
        from tidal.solver.ida import solve_ida

        num_snapshots = max(int(args.t_end / snapshot_interval) + 1, 2)
        log(f"Running IDA solver (t=0 → {args.t_end}, {num_snapshots} snapshots)...")
        result = solve_ida(
            spec, grid_info, y0,
            t_span=(0.0, args.t_end),
            bc=bc, parameters=params,
            num_snapshots=num_snapshots,
            snapshot_callback=snapshot_cb,
        )
    else:  # leapfrog
        from tidal.solver.leapfrog import solve_leapfrog

        dt = args.dt
        if dt is None:
            dt = _CFL_FACTOR * min(float(d) for d in grid_info.dx)
        log(f"Running leapfrog solver (t=0 → {args.t_end}, dt={dt:.4f})...")
        result = solve_leapfrog(
            spec, grid_info, y0,
            t_span=(0.0, args.t_end), dt=dt,
            bc=bc, parameters=params,
            snapshot_interval=snapshot_interval,
            snapshot_callback=snapshot_cb,
        )

    if not result["success"]:
        print(f"Error: solver failed: {result['message']}", file=sys.stderr)

    if writer is not None:
        writer.close()
        log(f"  {writer.count} snapshots streamed to: {writer.output_dir}")

    # 8. Build SimulationData
    sim_data = SimulationData.from_result(result, spec, grid_info, params)
    log(f"  {sim_data.n_snapshots} snapshots stored")

    # 9. Output
    _generate_output_native(args, sim_data, grid_info)
    return 0


# --- Command entry point ---


def _check_cfl_stability(pde: object, dt: float, grid: CartesianGrid) -> None:
    """Print CFL stability warnings to stderr."""
    warnings_list = cast("list[str]", pde.check_stability(dt, grid))  # type: ignore[attr-defined]
    sys.stderr.writelines(f"  Warning: {w}\n" for w in warnings_list)


def _validate_solver_params(args: Namespace) -> None:
    """Validate solver-related CLI arguments.

    Raises
    ------
    ValueError
        If ``--t-end``, ``--dt``, or ``--snapshots`` are non-positive.
    """
    if args.t_end <= 0:
        msg = f"--t-end must be positive, got {args.t_end}"
        raise ValueError(msg)
    if args.dt is not None and args.dt <= 0:
        msg = f"--dt must be positive, got {args.dt}"
        raise ValueError(msg)
    if args.snapshots is not None and args.snapshots <= 0:
        msg = f"--snapshots must be positive, got {args.snapshots}"
        raise ValueError(msg)


def _warn_zero_evolution(pde: object, state: FieldCollection) -> None:
    """Warn if all evolution rates are zero at t=0.

    Catches initial conditions that produce static configurations, e.g.
    pure-gauge EM with Gaussian A_1 and zero conjugate momentum.
    """
    rates = pde.evolution_rate(state, t=0.0)  # type: ignore[attr-defined]
    max_rate = max(float(np.max(np.abs(f.data))) for f in rates)
    if max_rate < _ZERO_RATE_THRESHOLD:
        print(
            "  Warning: all evolution rates are zero at t=0. "
            "The initial condition may be a static configuration. "
            "For gauge fields, try --ic plane-wave to provide "
            "non-zero conjugate momentum.",
            file=sys.stderr,
        )


def simulate_command(args: Namespace) -> int:  # noqa: C901, PLR0914, PLR0915
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
    from tidal.symbolic import load_equation_system

    json_path = Path(args.json_path)
    if not json_path.exists():
        print(f"Error: file not found: {json_path}", file=sys.stderr)
        return 1

    # Progress printer — suppressed by --quiet
    def _noop(*_a: object, **_kw: object) -> None:
        pass

    log = _noop if args.quiet else print

    # Step 1: Load spec
    log("Loading equation specification...")
    spec = load_equation_system(json_path)
    log(
        f"  {spec.n_components} component(s), {spec.dimension}D ({spec.spatial_dimension}+1D)"
    )

    # Step 2: Parse parameters
    params = _parse_params(args.param, spec)
    if params:
        log(f"  Parameters: {params}")

    # ─── Native path (IDA / leapfrog) ───────────────────────
    # Self-contained: no py-pde imports, no type conversions.
    if args.scheme in {"ida", "leapfrog"}:
        return _simulate_native(args, spec, params)

    # ─── Legacy py-pde path (scipy / runge-kutta) ──────────
    from tidal.symbolic import build_pde_from_json

    # Step 3: Build PDE
    log("Building PDE...")
    pde = build_pde_from_json(json_path, parameters=params)

    # Step 4: Create grid
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    grid = _build_grid(args, spec, bounds)
    log(f"  Grid: {'x'.join(str(s) for s in grid.shape)}, bounds: {grid.axes_bounds}")

    # Step 5: Initial conditions
    state = _build_initial_state(args, grid, spec, bounds)
    initial_state = state.copy()
    log(f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}")

    # Safety net: warn if all evolution rates are zero at t=0
    # (e.g., pure-gauge EM with Gaussian IC and π=0)
    if args.mode != "constraint":
        _warn_zero_evolution(pde, state)

    # Constraint-only mode: single constraint solve, no time evolution
    if args.mode == "constraint":
        log("Solving constraints...")
        pde.evolution_rate(state, t=0.0)  # type: ignore[attr-defined]
        log("  Constraint solve complete.")
        _print_constraint_summary(spec, initial_state, state, params)
        _save_constraint_output(args, spec, state)
        return 0

    # Step 6: Validate solver parameters and determine dt
    _validate_solver_params(args)

    # py-pde path: runge-kutta or scipy
    from tidal.utils import normalize_solve_result

    dt = args.dt
    if dt is None:
        # Auto dt from grid spacing
        dt = _CFL_FACTOR * min(
            (b[1] - b[0]) / s for b, s in zip(grid.axes_bounds, grid.shape, strict=True)
        )

    _check_cfl_stability(pde, dt, grid)

    # Step 7: Run simulation
    snapshot_interval = (
        args.snapshots if args.snapshots is not None else args.t_end / 100.0
    )
    log(
        f"Running simulation (t=0 → {args.t_end}, dt={dt:.4f}, scheme={args.scheme})..."
    )

    fmt = _infer_output_format(args)
    use_directory = fmt == "directory"
    writer_pypde: SnapshotWriter | None = None

    if use_directory:
        # Disk-backed streaming: O(1) memory regardless of snapshot count
        storage, tracker, writer_pypde = _setup_disk_backed(
            args,
            spec,
            grid,
            snapshot_interval,
            params,
        )
    else:
        # In-memory path (for plot output)
        from pde import MemoryStorage

        storage = MemoryStorage()
        tracker = storage.tracker(snapshot_interval)

    # Solve (tracker type varies: StorageTracker for legacy, list for disk-backed)
    if args.scheme == "scipy":
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]
                state,
                t_range=args.t_end,
                dt=dt,
                solver="scipy",
                tracker=tracker,  # pyright: ignore[reportArgumentType]
            )
        )
    else:
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]
                state,
                t_range=args.t_end,
                dt=dt,
                scheme=args.scheme,
                tracker=tracker,  # pyright: ignore[reportArgumentType]
            )
        )

    if writer_pypde is not None:
        writer_pypde.close()
        log(f"  {writer_pypde.count} snapshots streamed to: {writer_pypde.output_dir}")
    else:
        log(f"  {len(storage)} snapshots stored")

    # Step 8: Output
    ctx = PlotContext(
        spec=spec,
        storage=storage,
        grid=grid,
        initial_state=initial_state,
        params=params,
    )
    _generate_output(args, ctx)

    return 0


def _setup_disk_backed(
    args: Namespace,
    spec: EquationSystem,
    grid: CartesianGrid,
    snapshot_interval: float,
    params: dict[str, float],
) -> tuple[MemoryStorage, list[object], SnapshotWriter]:
    """Set up SnapshotWriter + CallbackTracker for disk-backed simulation.

    Also creates a minimal MemoryStorage that only stores the final snapshot
    (for PlotContext summary compatibility).

    Returns (storage, tracker_list, writer).
    """
    from pde import CallbackTracker, MemoryStorage

    from tidal.measurement._writer import create_snapshot_callback

    output_dir = Path(args.output) if args.output else Path("output")

    writer, callback = create_snapshot_callback(
        output_dir=output_dir,
        spec=spec,
        grid=grid,
        t_end=args.t_end,
        snapshot_interval=snapshot_interval,
        parameters=params,
        spec_path=Path(args.json_path),
    )

    snapshot_tracker = CallbackTracker(callback, interrupts=snapshot_interval)

    # Minimal MemoryStorage for PlotContext summary (only the final snapshot)
    storage = MemoryStorage()
    summary_tracker = storage.tracker(args.t_end)

    return storage, [snapshot_tracker, summary_tracker], writer
