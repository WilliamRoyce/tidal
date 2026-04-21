"""``tidal simulate`` — Run PDE simulation from a JSON specification."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from tidal.cli._console import debug as _cdebug
from tidal.cli._console import error as _cerror
from tidal.cli._console import error_with_hint as _cerror_hint
from tidal.cli._console import log as _clog
from tidal.cli._console import warn as _cwarn

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Callable

    from tidal.measurement._io import SimulationData
    from tidal.measurement._writer import SnapshotWriter
    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import AxisBCSpec, BCSpec
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

# Laplacian-like operators that contribute to the wave speed
_LAPLACIAN_OPS = frozenset(
    {
        "laplacian",
        "laplacian_x",
        "laplacian_y",
        "laplacian_z",
    }
)

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


_NATIVE_BC_TYPES = frozenset({"periodic", "neumann", "dirichlet", "robin"})


def _parse_bc_entry(entry: str) -> tuple[str, dict[str, float]]:
    """Parse a single BC entry like ``"dirichlet:1.0"`` or ``"robin:gamma=1:beta=0"``.

    Returns ``(bc_type, params)`` where params is a dict of keyword arguments.

    Supported formats:
    - ``"periodic"`` -- simple BC type
    - ``"dirichlet:1.0"`` -- BC type with positional value
    - ``"neumann:deriv=0.5"`` -- BC type with keyword argument
    - ``"robin:gamma=1:beta=0"`` -- BC type with multiple keyword arguments

    Raises
    ------
    ValueError
        If the BC type is unknown or parameters are malformed.
    """
    parts = entry.strip().lower().split(":")
    bc_type = parts[0]
    if bc_type not in _NATIVE_BC_TYPES:
        msg = (
            f"Invalid boundary condition: '{bc_type}'. "
            f"Must be one of: {', '.join(sorted(_NATIVE_BC_TYPES))}"
        )
        raise ValueError(msg)

    params: dict[str, float] = {}
    for part in parts[1:]:
        if "=" in part:
            key, val_str = part.split("=", 1)
            try:
                params[key.strip()] = float(val_str.strip())
            except ValueError:
                msg = f"Invalid BC parameter value: '{part}' (expected key=number)"
                raise ValueError(msg) from None
        else:
            # Positional value: "dirichlet:1.0" or "neumann:0.5"
            try:
                val = float(part.strip())
            except ValueError:
                msg = f"Invalid BC parameter: '{part}' (expected number or key=number)"
                raise ValueError(msg) from None
            if bc_type == "neumann":
                params["deriv"] = val
            else:
                params["value"] = val

    return bc_type, params


def _bc_entry_to_axis_bc(bc_type: str, params: dict[str, float]) -> AxisBCSpec:
    """Convert a parsed BC entry to an AxisBCSpec.

    Raises
    ------
    ValueError
        If the BC type is unknown or parameters are invalid.
    """
    from tidal.solver.operators import AxisBCSpec, SideBCSpec

    if bc_type == "periodic":
        if params:
            msg = "Periodic BC does not accept parameters"
            raise ValueError(msg)
        return AxisBCSpec(periodic=True)

    if bc_type == "dirichlet":
        side = SideBCSpec(kind="dirichlet", value=params.get("value", 0.0))
    elif bc_type == "neumann":
        side = SideBCSpec(kind="neumann", derivative=params.get("deriv", 0.0))
    elif bc_type == "robin":
        side = SideBCSpec(
            kind="robin",
            value=params.get("beta", 0.0),
            gamma=params.get("gamma", 0.0),
        )
    else:
        msg = f"Unknown BC type: {bc_type!r}"
        raise ValueError(msg)

    return AxisBCSpec(periodic=False, low=side, high=side)


def _parse_periodic(
    bc_str: str | None,
    *,
    periodic: bool,
    spatial_dim: int,
) -> tuple[bool, ...]:
    """Parse boundary spec into periodic flags for GridInfo.

    Returns a fixed-length tuple of booleans.  Accepts simple types
    (``"periodic"``, ``"neumann"``, ``"dirichlet"``, ``"robin"``)
    and extended syntax (``"dirichlet:1.0"``, ``"robin:gamma=1:beta=0"``).

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

    result: list[bool] = []
    for entry in bc_list:
        bc_type, _params = _parse_bc_entry(entry)
        result.append(bc_type == "periodic")

    return tuple(result)


def _parse_axis_bcs(
    bc_str: str | None,
    *,
    spatial_dim: int,
) -> tuple[AxisBCSpec, ...] | None:
    """Parse ``--bc`` string into structured AxisBCSpec objects.

    Returns ``None`` when ``--bc`` is not specified (use default inference).
    Returns a tuple of ``AxisBCSpec`` when explicit BCs are given.

    Raises
    ------
    ValueError
        If BC count doesn't match dimension or parameters are invalid.
    """
    if bc_str is None:
        return None

    bc_list = [b.strip().lower() for b in bc_str.split(",")]
    if len(bc_list) == 1:
        bc_list *= spatial_dim
    elif len(bc_list) != spatial_dim:
        msg = (
            f"--bc expects 1 or {spatial_dim} values "
            f"(got {len(bc_list)}). Example: --bc neumann,periodic"
        )
        raise ValueError(msg)

    specs: list[AxisBCSpec] = []
    for entry in bc_list:
        bc_type, params = _parse_bc_entry(entry)
        specs.append(_bc_entry_to_axis_bc(bc_type, params))

    return tuple(specs)


def _build_grid_info(
    args: Namespace,
    spec: EquationSystem,
    bounds: list[tuple[float, float]],
) -> GridInfo:
    """Build GridInfo from CLI arguments."""
    from tidal.solver.grid import GridInfo

    shape = _parse_grid_shape(args.grid_shape, spec.spatial_dimension)
    periodic = _parse_periodic(
        args.bc, periodic=args.periodic, spatial_dim=spec.spatial_dimension
    )
    axis_bcs = _parse_axis_bcs(args.bc, spatial_dim=spec.spatial_dimension)

    # Legacy string BC tuple for backward compat with GridInfo.bc
    bc: tuple[str, ...] | None = None
    if args.bc:
        bc_list = [b.strip().lower().split(":")[0] for b in args.bc.split(",")]
        if len(bc_list) == 1:
            bc = tuple(bc_list[0] for _ in range(spec.spatial_dimension))
        else:
            bc = tuple(bc_list)

    return GridInfo(
        bounds=tuple(bounds),
        shape=tuple(shape),
        periodic=periodic,
        bc=bc,
        axis_bcs=axis_bcs,
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


# --- Initial condition helpers ---
#
# Architecture: Each IC type has a builder function that returns
# ``dict[str, np.ndarray]`` (slot_data) mapping slot names to arrays.
# ``_build_initial_y0()`` orchestrates the pipeline:
#   1. Validate component name
#   2. Call the type-specific builder → slot_data dict
#   3. Apply ``--ic-field`` per-field formula overrides
#   4. Pack via ``FieldSet.from_dict()`` → flat state vector
#
# The constraint IC solver (``ensure_consistent_ic``) runs later inside
# the solver functions, so constraint fields may be overwritten there.


def _eval_formula_expr(
    expr: str,
    spec: EquationSystem,
    grid_info: GridInfo,
) -> np.ndarray:
    """Evaluate a formula expression against spatial coordinates.

    Validates the expression via AST analysis, then evaluates it in a
    sandboxed namespace containing ``FORMULA_NAMESPACE`` + spatial coords.
    """
    namespace = dict(FORMULA_NAMESPACE)
    for i, name in enumerate(spec.spatial_coordinates):
        namespace[name] = grid_info.cell_coords[..., i]

    allowed_names = set(namespace.keys())
    _validate_formula_ast(expr, allowed_names)

    result = eval(expr, {"__builtins__": {}}, namespace)  # noqa: S307
    result = np.asarray(result, dtype=float)

    if result.shape == ():
        result = np.full(grid_info.shape, float(result))

    return result


def _has_velocity_slot(
    layout: object,
    component: str,
) -> bool:
    """Check whether a field has a velocity slot in the state layout."""
    from tidal.solver.state import StateLayout

    assert isinstance(layout, StateLayout)
    return component in layout.velocity_slot_map


def _gaussian_slots(  # noqa: PLR0913, PLR0917
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
    component: str,
    layout: object,
) -> dict[str, np.ndarray]:
    """Compute Gaussian IC as slot_data dict.

    When ``--ic-wavevector`` is provided, creates a travelling wave packet:
    a Gaussian envelope modulated by a carrier wave, with matching velocity
    for unidirectional propagation.  Positive wavevector gives a right-mover.

    Raises
    ------
    ValueError
        If ``--ic-center`` dimension count doesn't match spatial dimension.
    """
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

    envelope = args.ic_amplitude * np.exp(-dist_sq / (2 * width**2))

    slot_data: dict[str, np.ndarray] = {}

    if args.ic_wavevector is not None:
        kvec = tuple(float(k) for k in args.ic_wavevector.split(","))
        k_dot_x = np.zeros(grid_info.shape, dtype=np.float64)
        for dim in range(min(grid_info.ndim, len(kvec))):
            k_dot_x += kvec[dim] * coords[..., dim]
        k_mag = float(np.sqrt(sum(k**2 for k in kvec)))
        slot_data[component] = envelope * np.cos(k_dot_x)
        if _has_velocity_slot(layout, component):
            slot_data[f"v_{component}"] = envelope * k_mag * np.sin(k_dot_x)
    else:
        slot_data[component] = envelope

    return slot_data


def _plane_wave_slots(  # noqa: PLR0913, PLR0917
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
    component: str,
    layout: object,
) -> dict[str, np.ndarray]:
    """Compute plane-wave IC as slot_data dict.

    Uses ``cos(k·x)`` for field and ``+|k|·sin(k·x)`` for velocity
    (right-mover for positive k).
    """
    if args.ic_wavevector is not None:
        kvec = tuple(float(k) for k in args.ic_wavevector.split(","))
    else:
        lx = bounds[0][1] - bounds[0][0]
        kvec = tuple(
            2.0 * math.pi / lx if i == 0 else 0.0 for i in range(spec.spatial_dimension)
        )

    coords = grid_info.cell_coords
    k_dot_x = np.zeros(grid_info.shape, dtype=np.float64)
    for dim in range(min(grid_info.ndim, len(kvec))):
        k_dot_x += kvec[dim] * coords[..., dim]

    k_mag = float(np.sqrt(sum(k**2 for k in kvec)))
    amplitude = args.ic_amplitude

    slot_data: dict[str, np.ndarray] = {
        component: amplitude * np.cos(k_dot_x),
    }
    if _has_velocity_slot(layout, component):
        slot_data[f"v_{component}"] = amplitude * k_mag * np.sin(k_dot_x)

    return slot_data


def _formula_slots(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    component: str,
    layout: object,
) -> dict[str, np.ndarray]:
    """Compute formula-based IC as slot_data dict.

    Supports optional ``--ic-formula-velocity`` for setting the velocity
    slot alongside the field, enabling custom travelling wave ICs.

    Raises
    ------
    ValueError
        If ``--ic-formula`` is not provided.
    """
    if args.ic_formula is None:
        msg = "--ic=formula requires --ic-formula=EXPR"
        raise ValueError(msg)

    slot_data: dict[str, np.ndarray] = {
        component: _eval_formula_expr(args.ic_formula, spec, grid_info),
    }

    vel_expr = getattr(args, "ic_formula_velocity", None)
    if vel_expr is not None:
        if _has_velocity_slot(layout, component):
            slot_data[f"v_{component}"] = _eval_formula_expr(vel_expr, spec, grid_info)
        else:
            print(
                f"  Warning: --ic-formula-velocity ignored for '{component}' "
                f"(no velocity slot — first-order or constraint field)",
                file=sys.stderr,
            )

    return slot_data


def _file_slots_npy(
    ic_path: Path,
    spec: EquationSystem,
    grid_info: GridInfo,
) -> dict[str, np.ndarray]:
    """Load flat state vector from .npy and unpack into slot_data.

    Raises
    ------
    ValueError
        If .npy size doesn't match state vector layout.
    """
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    flat = np.load(ic_path)
    if flat.shape != (layout.total_size,):
        msg = (
            f"IC .npy has {flat.size} elements but state vector needs "
            f"{layout.total_size} (fields={len(spec.equations)}, "
            f"grid={grid_info.shape})"
        )
        raise ValueError(msg)

    slot_data: dict[str, np.ndarray] = {}
    for slot in layout.slots:
        idx = layout.field_slot_map.get(slot.name)
        if idx is None and slot.kind == "velocity":
            idx = layout.velocity_slot_map.get(slot.field_name)
        if idx is not None:
            start = idx * layout.num_points
            end = start + layout.num_points
            slot_data[slot.name] = flat[start:end].reshape(grid_info.shape).copy()
    return slot_data


def _file_slots_dir(
    ic_path: Path,
    spec: EquationSystem,
    grid_info: GridInfo,
) -> dict[str, np.ndarray]:
    """Load final snapshot from simulation output directory.

    Raises
    ------
    ValueError
        If saved grid shape doesn't match current grid.
    """
    import json

    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)

    meta_path = ic_path / "metadata.json"
    if meta_path.exists():
        with meta_path.open(encoding="utf-8") as f:
            meta = json.load(f)
        saved_shape = tuple(meta.get("grid_shape", []))
        if saved_shape and saved_shape != grid_info.shape:
            msg = (
                f"Grid shape mismatch: saved {saved_shape} vs current {grid_info.shape}"
            )
            raise ValueError(msg)

    slot_data: dict[str, np.ndarray] = {}
    for slot in layout.slots:
        npy_file = ic_path / f"{slot.name}.npy"
        if npy_file.exists():
            arr = np.load(npy_file)
            # Take final snapshot if multi-snapshot (n_snapshots, *grid_shape)
            if arr.ndim > len(grid_info.shape):
                arr = arr[-1]
            slot_data[slot.name] = arr
    return slot_data


def _file_slots(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
) -> dict[str, np.ndarray]:
    """Load IC from a .npy file or simulation output directory.

    Raises
    ------
    ValueError
        If ``--ic-file`` is not provided, path doesn't exist, or format unknown.
    """
    ic_path_str = getattr(args, "ic_file", None)
    if ic_path_str is None:
        msg = "--ic=file requires --ic-file=PATH"
        raise ValueError(msg)

    ic_path = Path(ic_path_str)
    if not ic_path.exists():
        msg = f"IC file not found: {ic_path}"
        raise ValueError(msg)

    if ic_path.suffix == ".npy":
        return _file_slots_npy(ic_path, spec, grid_info)

    if ic_path.is_dir():
        return _file_slots_dir(ic_path, spec, grid_info)

    msg = f"IC path must be a .npy file or directory, got: {ic_path}"
    raise ValueError(msg)


class ResumeState:  # noqa: B903
    """Checkpoint state loaded from a snapshot directory for simulation resume."""

    __slots__ = (
        "bc_types",
        "dt",
        "grid_bounds",
        "grid_shape",
        "parameters",
        "periodic",
        "snapshot_index",
        "t_start",
        "y0",
    )

    def __init__(  # noqa: PLR0913
        self,
        *,
        y0: np.ndarray,
        t_start: float,
        parameters: dict[str, float],
        grid_shape: tuple[int, ...],
        grid_bounds: tuple[tuple[float, float], ...],
        periodic: tuple[bool, ...],
        bc_types: tuple[str, ...] | None,
        dt: float | None,
        snapshot_index: int,
    ) -> None:
        self.y0 = y0
        self.t_start = t_start
        self.parameters = parameters
        self.grid_shape = grid_shape
        self.grid_bounds = grid_bounds
        self.periodic = periodic
        self.bc_types = bc_types
        self.dt = dt
        self.snapshot_index = snapshot_index


def _load_resume_state(  # noqa: PLR0914
    resume_dir: Path,
    spec: EquationSystem,
    snapshot_index: int | None = None,
) -> ResumeState:
    """Load checkpoint state from a snapshot directory for resume.

    Parameters
    ----------
    resume_dir : Path
        Path to a simulation output directory (containing metadata.json,
        times.npy, and per-field .npy files).
    spec : EquationSystem
        Equation system specification (used to build StateLayout and
        validate field compatibility).
    snapshot_index : int | None
        Which snapshot to load (0-based). ``None`` loads the last snapshot.

    Returns
    -------
    ResumeState
        Checkpoint state ready for solver initialization.

    Raises
    ------
    FileNotFoundError
        If the resume directory or required files are missing.
    ValueError
        If field names don't match or snapshot index is out of range.
    """
    import json

    from tidal.solver.fields import FieldSet
    from tidal.solver.state import StateLayout

    # 1. Load metadata
    meta_path = resume_dir / "metadata.json"
    if not meta_path.exists():
        msg = f"Resume directory missing metadata.json: {resume_dir}"
        raise FileNotFoundError(msg)

    with meta_path.open(encoding="utf-8") as f:
        meta = json.load(f)

    saved_fields: list[str] = meta.get("fields", [])
    grid_shape = tuple(meta.get("grid_shape", []))
    grid_bounds = tuple(tuple(b) for b in meta.get("grid_bounds", []))
    periodic = tuple(meta.get("periodic", []))
    bc_types_raw = meta.get("bc_types")
    bc_types = tuple(bc_types_raw) if bc_types_raw is not None else None
    parameters: dict[str, float] = meta.get("parameters", {})
    dt = meta.get("dt")

    # 2. Validate field compatibility
    spec_fields = set(spec.component_names)
    saved_field_set = set(saved_fields)
    missing = spec_fields - saved_field_set
    if missing:
        msg = (
            f"Checkpoint missing fields required by spec: {sorted(missing)}. "
            f"Saved fields: {sorted(saved_field_set)}"
        )
        raise ValueError(msg)

    # 3. Load times and determine snapshot index
    times_path = resume_dir / "times.npy"
    if not times_path.exists():
        msg = f"Resume directory missing times.npy: {resume_dir}"
        raise FileNotFoundError(msg)

    times = np.load(times_path)
    n_snapshots = len(times)

    if snapshot_index is None:
        snapshot_index = n_snapshots - 1
    elif snapshot_index < 0 or snapshot_index >= n_snapshots:
        msg = f"Snapshot index {snapshot_index} out of range (0..{n_snapshots - 1})"
        raise ValueError(msg)

    t_start = float(times[snapshot_index])

    # 4. Load field and velocity data at the chosen snapshot
    layout = StateLayout.from_spec(spec, int(np.prod(grid_shape)))
    slot_data: dict[str, np.ndarray] = {}

    for slot in layout.slots:
        npy_file = resume_dir / f"{slot.name}.npy"
        if npy_file.exists():
            arr = np.load(npy_file)
            # Extract specific snapshot if multi-snapshot
            if arr.ndim > len(grid_shape):
                arr = arr[snapshot_index]
            slot_data[slot.name] = arr
        elif slot.kind in {"field", "velocity"}:
            import warnings

            warnings.warn(
                f"Resume: '{slot.name}' not found in checkpoint "
                f"({resume_dir}) — defaulting to zero",
                stacklevel=2,
            )

    # 5. Pack into flat state vector
    y0 = FieldSet.from_dict(layout, grid_shape, slot_data).flat.copy()

    return ResumeState(
        y0=y0,
        t_start=t_start,
        parameters=parameters,
        grid_shape=grid_shape,
        grid_bounds=grid_bounds,
        periodic=periodic,
        bc_types=bc_types,
        dt=dt,
        snapshot_index=snapshot_index,
    )


def _validate_resume_grid(resume: ResumeState, grid_info: GridInfo) -> None:
    """Validate that the resume checkpoint is compatible with the current grid.

    Raises
    ------
    ValueError
        If grid shape or bounds don't match.
    """
    if resume.grid_shape != grid_info.shape:
        msg = (
            f"Grid shape mismatch: checkpoint has {resume.grid_shape} "
            f"but current grid is {grid_info.shape}"
        )
        raise ValueError(msg)
    # Compare bounds with tolerance for float rounding
    bounds_tol = 1e-10
    for i, (saved, current) in enumerate(
        zip(resume.grid_bounds, grid_info.bounds, strict=True)
    ):
        if (
            abs(saved[0] - current[0]) > bounds_tol
            or abs(saved[1] - current[1]) > bounds_tol
        ):
            msg = (
                f"Grid bounds mismatch on axis {i}: "
                f"checkpoint has {saved} but current grid is {current}"
            )
            raise ValueError(msg)


def _noise_slots(
    args: Namespace,
    grid_info: GridInfo,
    component: str,
) -> dict[str, np.ndarray]:
    """Generate white Gaussian noise IC on a single field.

    Reproducible with ``--ic-noise-seed``. Velocity is zero.
    """
    seed = getattr(args, "ic_noise_seed", None)
    rng = np.random.default_rng(seed)

    return {
        component: args.ic_amplitude * rng.standard_normal(grid_info.shape),
    }


def _apply_ic_perturbation(
    y0: np.ndarray,
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
) -> np.ndarray:
    """Add small Gaussian noise to field slots (not velocities) for ensemble variation.

    The perturbation scale is relative to ``--ic-amplitude``:
    ``perturbation = scale * amplitude * N(0,1)``.

    This is analogous to ensemble weather forecasting where perturbed ICs
    generate ensemble members for uncertainty estimation.

    References
    ----------
    Palmer, T.N. et al. (1993) "Ensemble prediction", ECMWF Tech Memo 188.

    Parameters
    ----------
    y0 : ndarray
        Initial state vector (modified in-place and returned).
    args : Namespace
        CLI args; uses ``ic_perturbation`` (scale) and ``ic_perturbation_seed``.
    spec : EquationSystem
        Equation system for layout construction.
    grid_info : GridInfo
        Grid information for layout construction.

    Returns
    -------
    ndarray
        The perturbed state vector (same object as *y0*).
    """
    from tidal.solver.state import StateLayout

    scale = getattr(args, "ic_perturbation", None)
    if scale is None or scale == 0.0:
        return y0

    seed = getattr(args, "ic_perturbation_seed", None)
    rng = np.random.default_rng(seed)
    amplitude = getattr(args, "ic_amplitude", 1.0) or 1.0
    layout = StateLayout.from_spec(spec, grid_info.num_points)
    n = layout.num_points

    for i, slot in enumerate(layout.slots):
        if slot.kind == "field":
            start = i * n
            y0[start : start + n] += scale * amplitude * rng.standard_normal(n)

    return y0


def _apply_ic_field_overrides(
    slot_data: dict[str, np.ndarray],
    ic_field_args: list[str],
    spec: EquationSystem,
    grid_info: GridInfo,
    layout: object,
) -> None:
    """Apply ``--ic-field`` per-field formula overrides to slot_data.

    Each entry has the format ``FIELD:EXPR`` (sets field slot) or
    ``FIELD:velocity:EXPR`` (sets velocity slot). Modifies slot_data
    in-place.

    Raises
    ------
    ValueError
        If field name is unknown or format is invalid.
    """
    for entry in ic_field_args:
        parts = entry.split(":", maxsplit=2)

        if len(parts) == 2:  # noqa: PLR2004
            field_name, expr = parts
            is_velocity = False
        elif len(parts) == 3 and parts[1] == "velocity":  # noqa: PLR2004
            field_name, _, expr = parts
            is_velocity = True
        else:
            msg = (
                f"Invalid --ic-field format: '{entry}'. "
                f"Expected FIELD:EXPR or FIELD:velocity:EXPR"
            )
            raise ValueError(msg)

        if field_name not in spec.component_names:
            msg = (
                f"Unknown field '{field_name}' in --ic-field. "
                f"Available: {', '.join(spec.component_names)}"
            )
            raise ValueError(msg)

        arr = _eval_formula_expr(expr, spec, grid_info)

        if is_velocity:
            if _has_velocity_slot(layout, field_name):
                slot_data[f"v_{field_name}"] = arr
            else:
                print(
                    f"  Warning: velocity override for '{field_name}' ignored "
                    f"(no velocity slot)",
                    file=sys.stderr,
                )
        else:
            slot_data[field_name] = arr


def _build_initial_y0(
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    bounds: list[tuple[float, float]],
) -> np.ndarray:
    """Build initial state as flat numpy vector.

    Pipeline:
      1. Validate ``--ic-component``
      2. Dispatch to type-specific builder → ``slot_data`` dict
      3. Apply ``--ic-field`` per-field formula overrides
      4. Pack via ``FieldSet.from_dict()`` → flat state vector

    The constraint IC solver (``ensure_consistent_ic``) runs later inside
    the solver, so constraint fields set here may be adjusted for consistency.

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

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    ic_type = args.ic
    component = args.ic_component or spec.component_names[0]

    if component not in spec.component_names:
        msg = (
            f"Unknown component '{component}'. "
            f"Available: {', '.join(spec.component_names)}"
        )
        raise ValueError(msg)

    # Step 1: Build base IC as slot_data dict
    if ic_type == "zero":
        if args.ic_component is not None:
            print(
                f"  Note: --ic-component '{args.ic_component}' is ignored for zero IC"
            )
        slot_data: dict[str, np.ndarray] = {}

    elif ic_type == "gaussian":
        slot_data = _gaussian_slots(args, spec, grid_info, bounds, component, layout)

    elif ic_type == "plane-wave":
        slot_data = _plane_wave_slots(args, spec, grid_info, bounds, component, layout)

    elif ic_type == "formula":
        slot_data = _formula_slots(args, spec, grid_info, component, layout)

    elif ic_type == "file":
        slot_data = _file_slots(args, spec, grid_info)

    elif ic_type == "noise":
        slot_data = _noise_slots(args, grid_info, component)

    else:
        msg = f"Unknown IC type: {ic_type}"
        raise ValueError(msg)

    # Step 2: Apply per-field formula overrides
    ic_field_list: list[str] = getattr(args, "ic_field", None) or []
    if ic_field_list:
        _apply_ic_field_overrides(slot_data, ic_field_list, spec, grid_info, layout)

    # Step 3: Pack into flat state vector
    return FieldSet.from_dict(layout, grid_info.shape, slot_data).flat.copy()


# --- Native output pipeline (no py-pde) ---


def _print_summary(sim_data: SimulationData) -> None:
    """Print simulation summary using SimulationData (no py-pde types)."""
    times = sim_data.times
    print()
    print("Results:")
    print(
        f"  Time range: {float(times[0]):.2f} → {float(times[-1]):.2f} ({len(times)} snapshots)"
    )
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


def _generate_output(
    args: Namespace,
    sim_data: SimulationData,
    grid_info: GridInfo,
) -> None:
    """Generate output for the native solver path (no py-pde types)."""
    fmt = _infer_output_format(args)

    if fmt in {"summary", "directory"}:
        if sim_data.n_snapshots > 0:
            _print_summary(sim_data)
        if fmt == "directory" and args.output:
            from tidal.cli._plot import save_plot

            overview = Path(args.output) / "overview.png"
            save_plot(overview, sim_data, grid_info)
        return

    _print_summary(sim_data)

    if args.output is not None:
        output_path = Path(args.output)
    else:
        json_file = Path(args.json_path).resolve()
        output_path = json_file.parent / f"{json_file.stem}_output.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    from tidal.cli._plot import save_plot

    save_plot(output_path, sim_data, grid_info)


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


# --- Native simulation path (no py-pde) ---


def _warn_zero_evolution(
    spec: EquationSystem,
    grid_info: GridInfo,
    y0: np.ndarray,
    params: dict[str, float],
    bc: BCSpec | None,
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

    # Inject zero constraint velocities so v_field_name refs resolve.
    # At t=0 the actual velocities are unknown (IDA computes them), but
    # zero is the best estimate for this diagnostic check.
    for eq in spec.equations:
        if eq.time_derivative_order == 0:
            fieldset.set_aux(
                f"v_{eq.field_name}",
                np.zeros(grid_info.shape),
            )

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


def _check_result_finite(result: SolverResult) -> None:
    """Raise SimulationDivergedError if the final state contains NaN or Inf.

    This is a single post-simulation check — zero per-step overhead.

    Raises
    ------
    SimulationDivergedError
        If the final state contains non-finite values.
    """
    from tidal.solver._exceptions import SimulationDivergedError

    y = result["y"]
    if len(y) == 0:
        return
    final = y[-1]
    if not np.isfinite(final).all():
        msg = (
            "Simulation produced non-finite values (NaN or Inf). "
            "The system is likely physically unstable."
        )
        raise SimulationDivergedError(msg)


def _setup_disk_writer_native(  # noqa: PLR0913, PLR0917
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    params: dict[str, float],
    snapshot_interval: float,
    dt: float | None = None,
    num_snapshots: int | None = None,
) -> tuple[SnapshotWriter, Callable[..., None]]:
    """Set up disk-backed SnapshotWriter using StateLayout (no py-pde).

    Returns (writer, snapshot_callback).  The callback accepts
    ``(t, y_flat)`` or ``(t, y_flat, yp_flat)`` — IDA passes ``yp``
    for constraint velocity extraction.
    """
    from tidal.measurement._writer import SnapshotWriter, compute_snapshot_count
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    output_dir = Path(args.output) if args.output else Path("output")
    n_snaps = num_snapshots or compute_snapshot_count(args.t_end, snapshot_interval)

    field_names = [s.name for s in layout.slots if s.kind != "velocity"]
    velocity_names = [s.field_name for s in layout.slots if s.kind == "velocity"]

    # Identify constraint fields whose velocities should be written to disk.
    # IDA provides exact ∂_t(constraint) via yp — needed for energy measurement.
    constraint_names = [
        eq.field_name for eq in spec.equations if eq.time_derivative_order == 0
    ]
    constraint_slot_map: dict[str, int] = {}
    for cname in constraint_names:
        if cname in layout.field_slot_map:
            constraint_slot_map[cname] = layout.field_slot_map[cname]

    # Include constraint velocity names alongside dynamical velocities
    all_velocity_names = velocity_names + list(constraint_slot_map.keys())

    writer = SnapshotWriter(
        output_dir=output_dir,
        field_names=field_names,
        velocity_names=all_velocity_names,
        grid_shape=grid_info.shape,
        n_snapshots=n_snaps,
        grid_spacing=tuple(float(d) for d in grid_info.dx),
        grid_bounds=grid_info.bounds,
        periodic=grid_info.periodic,
        parameters=params,
        spec_path=Path(args.json_path),
        bc_types=grid_info.bc_types,
        dt=dt,
    )

    # Build slot index maps from layout
    n_pts = grid_info.num_points
    shape = grid_info.shape
    field_set = set(field_names)
    velocity_set = set(velocity_names)

    field_slots_map: dict[str, int] = {}
    velocity_slots_map: dict[str, int] = {}
    for i, slot in enumerate(layout.slots):
        if slot.kind == "velocity":
            velocity_slots_map[slot.field_name] = i
        elif slot.name in field_set:
            field_slots_map[slot.name] = i

    def _disk_callback(
        t: float,
        y_flat: np.ndarray,
        yp_flat: np.ndarray | None = None,
    ) -> None:
        fields_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in field_slots_map.items()
        }
        vels_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in velocity_slots_map.items()
            if name in velocity_set
        }
        # Extract constraint velocities from IDA's yp vector.
        # For modal solver (no yp), provide zeros — constraint velocities
        # are not directly available from the eigendecomposition output.
        for cname, slot_idx in constraint_slot_map.items():
            if yp_flat is not None:
                vels_d[cname] = yp_flat[
                    slot_idx * n_pts : (slot_idx + 1) * n_pts
                ].reshape(shape)
            else:
                vels_d[cname] = np.zeros(shape)
        writer.append(t, fields_d, vels_d)

    return writer, _disk_callback


def _setup_memory_accumulator_native(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    grid_info: GridInfo,
    params: dict[str, float],
    snapshot_interval: float,
    dt: float | None = None,
    num_snapshots: int | None = None,
) -> tuple[Any, Callable[..., None]]:
    """In-memory twin of :func:`_setup_disk_writer_native`.

    Produces an :class:`InMemoryAccumulator` and a ``(t, y_flat, yp_flat)``
    callback that writes to it.  Intended for the inference likelihood
    path: the accumulator is materialised into a ``SimulationData`` at
    the end of ``_simulate`` and never touches disk.
    """
    from tidal.measurement._writer import InMemoryAccumulator, compute_snapshot_count
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid_info.num_points)
    n_snaps = num_snapshots or compute_snapshot_count(
        snapshot_interval * 10,  # placeholder; caller already computed num_snapshots
        snapshot_interval,
    )
    if num_snapshots is not None:
        n_snaps = num_snapshots

    field_names = [s.name for s in layout.slots if s.kind != "velocity"]
    velocity_names = [s.field_name for s in layout.slots if s.kind == "velocity"]

    # Constraint fields whose velocities we may want in the output.
    constraint_names = [
        eq.field_name for eq in spec.equations if eq.time_derivative_order == 0
    ]
    constraint_slot_map: dict[str, int] = {
        cname: layout.field_slot_map[cname]
        for cname in constraint_names
        if cname in layout.field_slot_map
    }

    all_velocity_names = velocity_names + list(constraint_slot_map.keys())

    accumulator = InMemoryAccumulator(
        field_names=field_names,
        velocity_names=all_velocity_names,
        grid_shape=grid_info.shape,
        n_snapshots=n_snaps,
        grid_spacing=tuple(float(d) for d in grid_info.dx),
        grid_bounds=grid_info.bounds,
        periodic=grid_info.periodic,
        parameters=params,
        bc_types=grid_info.bc_types,
        dt=dt,
    )

    n_pts = grid_info.num_points
    shape = grid_info.shape
    field_set = set(field_names)
    velocity_set = set(velocity_names)

    field_slots_map: dict[str, int] = {}
    velocity_slots_map: dict[str, int] = {}
    for i, slot in enumerate(layout.slots):
        if slot.kind == "velocity":
            velocity_slots_map[slot.field_name] = i
        elif slot.name in field_set:
            field_slots_map[slot.name] = i

    def _memory_callback(
        t: float,
        y_flat: np.ndarray,
        yp_flat: np.ndarray | None = None,
    ) -> None:
        fields_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in field_slots_map.items()
        }
        vels_d = {
            name: y_flat[idx * n_pts : (idx + 1) * n_pts].reshape(shape)
            for name, idx in velocity_slots_map.items()
            if name in velocity_set
        }
        for cname, slot_idx in constraint_slot_map.items():
            if yp_flat is not None:
                vels_d[cname] = yp_flat[
                    slot_idx * n_pts : (slot_idx + 1) * n_pts
                ].reshape(shape)
            else:
                vels_d[cname] = np.zeros(shape)
        accumulator.append(t, fields_d, vels_d)

    return accumulator, _memory_callback


def _extract_constraint_bc(
    spec: EquationSystem,
) -> tuple[str, ...] | None:
    """Extract BCs from constraint_solver blocks in the JSON spec.

    Returns a per-axis BC tuple if any constraint equation defines BCs,
    otherwise None (use global BC).
    """
    for eq in spec.equations:
        if eq.time_derivative_order != 0:
            continue
        cs = eq.constraint_solver
        if cs and cs.enabled and cs.boundary_conditions:
            bc_list: list[str] = []
            for coord in spec.spatial_coordinates:
                bc_entry = cs.boundary_conditions.get(coord)
                if bc_entry is not None:
                    bc_list.append(bc_entry.type)
                else:
                    bc_list.append("neumann")
            return tuple(bc_list)
    return None


def _constraint_mode(  # noqa: PLR0913, PLR0917
    args: Namespace,
    spec: EquationSystem,
    grid_info: GridInfo,
    y0: np.ndarray,
    params: dict[str, float],
    bc: BCSpec | None,
    log: object,
) -> int:
    """Solve constraints via IDA's algebraic initial condition solver.

    BCs are taken from the constraint_solver block in the JSON spec when
    available, overriding the global ``--bc`` argument.
    """
    from tidal.measurement._io import SimulationData
    from tidal.solver.ida import solve_ida

    log_fn = cast("Callable[..., None]", log)

    # Use constraint-specific BCs if defined in the JSON spec
    constraint_bc = _extract_constraint_bc(spec)
    if constraint_bc is not None:
        log_fn(f"  Using constraint BCs from spec: {constraint_bc}")
        bc = constraint_bc

    log_fn("Solving constraints via IDA...")

    # IDA's calc_initcond="yp0" adjusts yp0 to satisfy F(t0, y0, yp0)=0.
    # A short time span suffices — we only need IDA to find consistent
    # initial conditions, not evolve.
    result = solve_ida(
        spec,
        grid_info,
        y0,
        t_span=(0.0, 0.01),
        bc=bc,
        parameters=params,
        num_snapshots=2,
        allow_inconsistent_ic=getattr(args, "allow_inconsistent_ic", False),
    )

    if not result["success"]:
        _cerror_hint(
            f"constraint solve failed: {result['message']}",
            [
                "Check `--bc` for consistency with constraint equations",
                "Try `--allow-inconsistent-ic`",
            ],
        )
        return 1

    log_fn("  Constraint solve complete.")

    sim_data = SimulationData.from_result(result, spec, grid_info, params)
    _print_summary(sim_data)

    # Save constraint output if directory requested
    if args.output is not None:
        output_path = Path(args.output)
        output_path.mkdir(parents=True, exist_ok=True)
        for name, arr in sim_data.fields.items():
            np.save(str(output_path / f"{name}.npy"), arr[-1])
        log_fn(f"  Saved data to: {output_path}")

    return 0


def _compute_cfl_dt(
    spec: EquationSystem,
    grid: GridInfo,
    params: dict[str, float],
) -> float:
    """Compute CFL-safe dt accounting for position-dependent coefficients.

    For each dynamical equation (time_order >= 2), sums the absolute values
    of all laplacian-type coefficient magnitudes.  For position-dependent
    coefficients, evaluates on the grid and takes the spatial maximum.

    Returns ``CFL_FACTOR * dx_min / sqrt(max_c2)`` where ``max_c2`` is the
    maximum summed laplacian coefficient across all equations and grid points.
    Falls back to ``CFL_FACTOR * dx_min`` when no laplacian terms exist.
    """
    from tidal.symbolic._eval_utils import evaluate_coefficient

    dx_min = min(float(d) for d in grid.dx)
    coords = spec.coordinates[1:]  # spatial only
    coord_arrays = dict(zip(coords, grid.coord_arrays(), strict=False))

    max_c2 = 0.0

    for eq in spec.equations:
        if eq.time_derivative_order < 2:  # noqa: PLR2004
            continue

        # Sum |coefficient| for all laplacian-type self-terms
        c2_grid: float | np.ndarray = 0.0
        for term in eq.rhs_terms:
            if term.operator not in _LAPLACIAN_OPS:
                continue
            if term.field != eq.field_name:
                continue  # cross-field laplacian — skip for CFL estimate

            if term.coefficient_symbolic is not None and term.position_dependent:
                val = evaluate_coefficient(
                    term.coefficient_symbolic,
                    params,
                    coords,
                    coord_arrays=coord_arrays,
                )
                c2_grid = c2_grid + np.abs(val)  # noqa: PLR6104
            else:
                c2_grid = c2_grid + abs(term.coefficient)  # noqa: PLR6104

        if isinstance(c2_grid, np.ndarray):
            eq_max = float(np.max(c2_grid))
        else:
            eq_max = float(c2_grid)

        max_c2 = max(max_c2, eq_max)

    if max_c2 > 0:
        return float(_CFL_FACTOR * dx_min / math.sqrt(max_c2))

    return float(_CFL_FACTOR * dx_min)


def _has_dissipation(spec: EquationSystem) -> bool:
    """Check if the system has dissipative terms (``first_derivative_t``).

    Dissipation breaks symplecticity, so Yoshida 4th-order leapfrog should
    not be used.  Also drives IDA auto-selection in ``_resolve_scheme()``.
    """
    return any(
        term.operator == "first_derivative_t"
        for eq in spec.equations
        for term in eq.rhs_terms
    )


def _has_time_dependent_coeffs(spec: EquationSystem) -> bool:
    """Check if any RHS term has a time-dependent coefficient.

    Time-dependent coefficients mean the Hamiltonian is not conserved,
    and Yoshida's negative middle sub-step (w₂ < 0) can introduce
    artifacts when the system is non-autonomous.
    """
    return any(term.time_dependent for eq in spec.equations for term in eq.rhs_terms)


def _resolve_scheme(  # noqa: C901
    scheme: str,
    spec: EquationSystem,
    grid: GridInfo | None = None,
    bc: BCSpec | None = None,
) -> str:
    """Resolve ``'auto'`` to the best solver for the equation system.

    Auto-selection algorithm (checked in order):

    1. Modal eligible (flat metric + all-periodic + time-independent +
       supported operators) → modal. Handles constraints via Fourier-space
       Schur complement elimination (Hairer & Wanner 1996, Ch. VII).
    2. Constraint equations (time_order=0) not modal-eligible → IDA
       (algebraic constraints need DAE residual form).
    3. First-order (time_order=1) equations → IDA (diffusion/transport needs
       implicit time integration; eligible periodic first-order caught by 1).
    4. Dissipation (``first_derivative_t`` operator in any RHS) → IDA (breaks
       symplecticity, BDF handles well).
    5. No canonical Hamiltonian structure for second-order wave equations
       → warning + CVODE fallback (may indicate hand-crafted or legacy JSON).
    6. Pure wave (all second-order, Hamiltonian) → CVODE BDF (adaptive,
       tolerance-controlled, same SUNDIALS ecosystem as IDA).

    Raises
    ------
    RuntimeError
        If ``--scheme modal`` is explicitly requested but the system is not
        eligible (wrong metric, BCs, or operators).

    References
    ----------
        Hindmarsh et al., "SUNDIALS", ACM TOMS, 2005.
        Moler & Van Loan (2003), SIAM Review 45(1):3-49 (modal solver).
    """
    # When a perturbation-tagged spec will be routed through the
    # PerturbativeSolver, the relevant eligibility check is on the
    # base_spec (post-Gap-B demotion), not the full spec. The correction
    # RHS terms carry higher-order operators (d3_t, d4_t, mixed_T3_S1x)
    # that Pass 0 never sees — they only appear as Duhamel source
    # coefficients. Build the base spec once here for both explicit and
    # auto paths.
    # If the user's [perturbation] config is malformed (missing small
    # parameters, 3rd-order residual after demotion), base_spec's
    # ValueError carries an actionable message — propagate it instead
    # of silently falling back to the full spec, which would make the
    # modal-eligibility check fail for an unrelated reason and mislead
    # the user with "auto-selected: CVODE" (#277).
    eligibility_spec = spec.base_spec() if spec.has_corrections() else spec

    if scheme != "auto":
        if scheme == "modal" and grid is not None:
            # Validate modal eligibility when explicitly requested
            from tidal.solver.modal import can_use_modal

            if not can_use_modal(eligibility_spec, grid, bc):
                msg = (
                    "--scheme modal requested but system is not eligible. "
                    "Modal solver requires: flat metric, all-periodic BCs, "
                    "time-independent coefficients, and supported spatial "
                    "operators.  Use 'auto' or another solver."
                )
                raise RuntimeError(msg)
        return scheme

    # 1-2. Modal solver (handles constraints via Schur complement if eligible)
    #      Flat metric, all-periodic, time-independent, supported operators.
    #      Constraints are Fourier-eliminable if their operators have exact
    #      Fourier multipliers.
    if grid is not None:
        from tidal.solver.modal import can_use_modal

        if can_use_modal(eligibility_spec, grid, bc):
            return "modal"

    # 1b. Constraint equations not modal-eligible → IDA (DAE solver required)
    for eq in spec.equations:
        if eq.time_derivative_order == 0:
            return "ida"

    # 3. First-order (diffusion/transport) equations → IDA
    for eq in spec.equations:
        if eq.time_derivative_order == 1:
            return "ida"

    # 4. Dissipation (first_derivative_t in any RHS term) → IDA
    if _has_dissipation(spec):
        return "ida"

    # 5. No canonical Hamiltonian structure → fall through to CVODE with
    #    a warning.  Missing canonical indicates hand-crafted or legacy JSON;
    #    pipeline-derived specs should always include it.
    if spec.canonical is None:
        has_wave = any(eq.time_derivative_order >= 2 for eq in spec.equations)  # noqa: PLR2004
        if has_wave:
            import warnings

            warnings.warn(
                "Second-order wave equations missing canonical Hamiltonian "
                "structure in JSON spec.  This may indicate a hand-crafted "
                "or legacy JSON — consider using 'tidal derive' to generate "
                "specs from the Wolfram pipeline.",
                UserWarning,
                stacklevel=2,
            )

    # 6. Pure wave, Hamiltonian → CVODE BDF (adaptive, tolerance-controlled)
    return "cvode"


def _simulate(  # noqa: C901, PLR0911, PLR0912, PLR0914, PLR0915
    args: Namespace,
    spec: EquationSystem,
    params: dict[str, float],
    *,
    in_memory_out: list[Any] | None = None,
) -> int:
    """Run simulation via native TIDAL solver.

    Self-contained flow: GridInfo -> IC -> solve -> SimulationData -> output.
    Handles IDA, CVODE, scipy, and leapfrog schemes.

    Parameters
    ----------
    in_memory_out : list[SimulationData] | None
        If provided, skip all disk output (writer, plots, constraint .npy
        files, HTML report) and append the resulting ``SimulationData`` to
        the list instead.  Used by the Bayesian-inference likelihood path
        to avoid the disk round-trip that dominates per-evaluation wall
        time (issue #269).
    """
    from tidal.measurement._io import SimulationData
    from tidal.solver.operators import set_fd_order, set_spectral

    log = _clog

    # 0a. FD order — must be set before any operator evaluation.
    # CLI default is 4 (5-point Fornberg stencil); module default is 2
    # for backward compatibility with library/test code.
    fd_order: int = getattr(args, "fd_order", 4)
    set_fd_order(fd_order)
    if fd_order != 4:  # noqa: PLR2004
        log(f"  FD order: {fd_order}")

    # 1. Grid
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    grid_info = _build_grid_info(args, spec, bounds)
    log(
        f"  Grid: {'x'.join(str(s) for s in grid_info.shape)}, bounds: {grid_info.bounds}"
    )

    _cdebug(f"periodic={grid_info.periodic}, dx={grid_info.dx}")

    # 2. BC (stored in GridInfo, derive tuple for solver calls)
    bc = grid_info.effective_bc

    # 2a. Validate grid size vs FD order — a stencil of width (fd_order + 1)
    # requires at least that many grid points on each axis.
    min_n = min(grid_info.shape)
    required_n = fd_order + 1
    if min_n < required_n:
        fd_order_explicit = getattr(args, "fd_order", 4) != 4  # noqa: PLR2004
        if fd_order_explicit:
            msg = (
                f"Grid too small for --fd-order {fd_order}: minimum axis has "
                f"{min_n} points but stencil width requires >= {required_n}."
            )
            _cerror_hint(
                msg,
                [
                    f"Increase `--grid-shape` (need >= {required_n} points per axis)",
                    "Or use `--fd-order 2` for coarse grids",
                ],
            )
            return 1
        # Default fd-order 4 on a tiny grid — fall back to order 2
        fd_order = 2
        set_fd_order(fd_order)
        log(f"  FD order: reduced to {fd_order} (grid too small for order 4)")

    # 2b. Spectral mode — auto-detect or validate.
    # Three states: None (auto-detect), True (force on), False (force off).
    # Auto: enabled when ALL BCs are periodic (spectral requires periodicity).
    # Ref: Burns et al. (2020), Phys. Rev. Research 2:023068.
    spectral_arg = getattr(args, "spectral", None)
    all_periodic = all(grid_info.periodic)

    if spectral_arg is None:
        # Auto-detect: enable spectral when all BCs are periodic
        use_spectral = all_periodic
        if use_spectral:
            log("  Auto-selected: spectral operators (all BCs periodic)")
    elif spectral_arg:
        # User explicitly requested --spectral: validate periodic BCs
        use_spectral = True
        if not all_periodic:
            non_periodic = [i for i, p in enumerate(grid_info.periodic) if not p]
            msg = (
                f"--spectral requires all boundary conditions to be periodic. "
                f"Non-periodic axes: {non_periodic}."
            )
            _cerror_hint(
                msg,
                [
                    "Use `--periodic` for all axes",
                    "Or remove `--spectral` to use FD stencils",
                ],
            )
            return 1
    else:
        # User explicitly passed --no-spectral: force FD stencils
        use_spectral = False

    set_spectral(use_spectral)
    if use_spectral and spectral_arg is True:
        log("  Operators: spectral (FFT)")

    # 3. Initial conditions (or resume from checkpoint)
    resume_state: ResumeState | None = None
    t_start = 0.0

    resume_path = getattr(args, "resume", None)
    # Distinguish from sweep's boolean --resume (which means "resume sweep")
    if isinstance(resume_path, str):
        resume_state = _load_resume_state(
            Path(resume_path), spec, getattr(args, "snapshot", None)
        )
        _validate_resume_grid(resume_state, grid_info)
        y0 = resume_state.y0
        t_start = resume_state.t_start
        log(
            f"  IC: resume from {resume_path} "
            f"(snapshot {resume_state.snapshot_index}, t={t_start:.4f})"
        )
    else:
        y0 = _build_initial_y0(args, spec, grid_info, bounds)
        y0 = _apply_ic_perturbation(y0, args, spec, grid_info)
        ic_desc = f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}"
        ic_field_list: list[str] = getattr(args, "ic_field", None) or []
        if ic_field_list:
            ic_desc += f" + {len(ic_field_list)} field override(s)"
        log(ic_desc)

    # Handle --t-additional (only with --resume)
    if getattr(args, "t_additional", None) is not None:
        if resume_state is None:
            print(
                "Warning: --t-additional without --resume; ignored",
                file=sys.stderr,
            )
        else:
            args.t_end = t_start + args.t_additional

    # Validate t_end > t_start for resumed simulations
    if resume_state is not None and args.t_end <= t_start:
        _cerror_hint(
            f"--t-end ({args.t_end}) must be greater than checkpoint time ({t_start})",
            [
                f"Checkpoint is at t={t_start}. Use a larger `--t-end`",
                "Or use `--t-additional T` to extend by T time units",
            ],
        )
        return 1

    # 4. Diagnostics
    _validate_solver_params(args)

    # Resolve solver scheme (auto-select based on equation operators)
    scheme = _resolve_scheme(args.scheme, spec, grid_info, bc)
    if args.scheme == "auto":
        log(f"  Auto-selected solver: {scheme}")
    _cdebug(f"solver={scheme}, fd_order={fd_order}, spectral={use_spectral}")

    # Spectral + IDA incompatibility: spectral operators produce dense
    # coupling (every grid point depends on every other), incompatible
    # with IDA's sparse Jacobian infrastructure.
    if use_spectral and scheme == "ida":
        if spectral_arg is None:
            # Auto-detected spectral + IDA needed → silently disable spectral
            use_spectral = False
            set_spectral(False)
            log("  Note: IDA requires FD operators; spectral auto-disabled")
        elif args.scheme == "auto":
            # User explicitly requested --spectral, scheme auto-selected IDA
            # → switch to CVODE to honour spectral request
            scheme = "cvode"
            log("  Note: --spectral incompatible with IDA; switching to CVODE")
        else:
            # User explicitly requested both IDA + spectral — error
            msg = (
                "--spectral is incompatible with --scheme ida "
                "(spectral operators produce dense coupling). "
                "Use --scheme cvode, scipy, or leapfrog instead."
            )
            _cerror(msg)
            return 1
    # Modal solver operates in pure k-space — spectral operators are not
    # used during time evolution.  However, keep spectral=True so that
    # energy measurements use FFT operators matching the modal solver's
    # conserved Hamiltonian.  Without this, energy measurement uses FD
    # operators that differ from the exact Fourier Hamiltonian, producing
    # spurious conservation errors that increase with grid size.
    if use_spectral and scheme == "modal":
        log("  Note: modal solver uses k-space natively; spectral auto-disabled")
    log(f"  Scheme: {scheme}")

    # --- Dry-run: preview setup and exit ---
    if getattr(args, "dry_run", False):
        n_fields = spec.n_components
        n_eqs = len(spec.equations)
        grid_pts = 1
        for s in grid_info.shape:
            grid_pts *= s
        # Rough memory estimate: state vector + snapshots
        n_snapshots = max(int((args.t_end - t_start) / (args.t_end / 20.0)) + 1, 2)
        mem_bytes = (
            grid_pts * n_fields * 2 * 8 * n_snapshots
        )  # fields + velocities, float64
        if mem_bytes < 1024 * 1024:
            mem_str = f"{mem_bytes / 1024:.0f} KB"
        else:
            mem_str = f"{mem_bytes / (1024 * 1024):.1f} MB"
        spec_name = getattr(args, "json_path", "unknown")
        print(
            f"  Spec:     {Path(spec_name).name} ({n_fields} fields, {n_eqs} equations)"
        )
        print(
            f"  Grid:     {'x'.join(str(s) for s in grid_info.shape)} points, "
            f"bounds {grid_info.bounds}, "
            f"{'periodic' if all(grid_info.periodic) else 'mixed BCs'}"
        )
        print(
            f"  Solver:   {scheme} (auto-selected)"
            if args.scheme == "auto"
            else f"  Solver:   {scheme}"
        )
        print(f"  FD order: {fd_order}")
        print(f"  Steps:    ~{n_snapshots} snapshots, t={t_start}→{args.t_end}")
        print(f"  Est. memory: ~{mem_str}")
        return 0

    # Constraint-only mode: solve algebraic equations via IDA, no time evolution
    if args.mode == "constraint":
        return _constraint_mode(args, spec, grid_info, y0, params, bc, log)

    import contextlib

    # Diagnostic warning that builds its own CoefficientEvaluator + runs
    # one full RHS evaluation.  Useful for catching IC bugs in the
    # `tidal simulate` user flow but pure overhead (~3ms/call) in the
    # nested-sampling likelihood path — skip it when we're in-memory
    # mode (see #269, #291).
    if in_memory_out is None:
        with contextlib.suppress(ValueError):
            # May raise ValueError for systems with time-derivative operators
            # (d2_t, mixed_T2_S1x, etc.) that the physical-space RHS evaluator
            # cannot handle. These are simulated by the modal solver's
            # generalized mass-matrix path which works directly in Fourier space.
            _warn_zero_evolution(spec, grid_info, y0, params, bc)
    # Note: mass stability pre-check removed — the modal solver's eigenvalue
    # pre-check (in _evolve_per_mode) provides more accurate instability
    # detection using the full evolution matrix, not a simplified mass proxy.

    # 5. Compute dt for leapfrog (needed before snapshot configuration)
    dt: float | None = None
    lf_order_arg: int | None = getattr(args, "leapfrog_order", None)
    lf_order: int = lf_order_arg if lf_order_arg is not None else 2
    if scheme == "leapfrog":
        # 5a. Auto-detect leapfrog order when not explicitly specified.
        # Yoshida 4th-order is preferred for time-independent, non-dissipative
        # systems: O(dt⁴) accuracy allows ~2x larger dt → net speedup.
        # Ref: Yoshida (1990), Physics Letters A 150(5-7), pp. 262-268.
        dissipative = _has_dissipation(spec)
        time_dep = _has_time_dependent_coeffs(spec)
        if lf_order_arg is None:
            if not dissipative and not time_dep:
                lf_order = 4
                log(
                    "  Auto-selected: Yoshida 4th-order leapfrog "
                    "(time-independent, non-dissipative system)"
                )
            else:
                lf_order = 2
                reasons: list[str] = []
                if dissipative:
                    reasons.append("dissipative terms")
                if time_dep:
                    reasons.append("time-dependent coefficients")
                log(
                    f"  Auto-selected: Störmer-Verlet 2nd-order leapfrog "
                    f"({', '.join(reasons)} detected)"
                )
        elif lf_order_arg == 4:  # noqa: PLR2004
            # User explicitly requested Yoshida — warn if inappropriate.
            if dissipative:
                import warnings

                warnings.warn(
                    "Yoshida 4th-order leapfrog is not recommended for "
                    "dissipative systems (first_derivative_t terms detected). "
                    "Dissipation breaks symplecticity. Consider --scheme ida "
                    "or --leapfrog-order 2.",
                    stacklevel=2,
                )
            if time_dep:
                import warnings

                warnings.warn(
                    "Yoshida 4th-order leapfrog: time-dependent coefficients "
                    "detected. The negative middle sub-step (w₂ ≈ -1.70) "
                    "evolves backward in time, which may introduce artifacts "
                    "for non-autonomous systems.",
                    stacklevel=2,
                )

        dt = args.dt
        if dt is None:
            dt = _compute_cfl_dt(spec, grid_info, params)
            _cdebug(f"CFL dt={dt:.6e} (safety={_CFL_FACTOR})")
        # Yoshida CFL correction: the effective CFL limit is reduced by
        # max(|wᵢ|) ≈ 1.70 because the middle sub-step has |w₂| > 1.
        # Must happen before snapshot configuration so the writer
        # pre-allocates the correct number of snapshots.
        if lf_order == 4:  # noqa: PLR2004
            from tidal.solver.leapfrog import YOSHIDA_WEIGHTS

            cfl_factor = max(abs(w) for w in YOSHIDA_WEIGHTS)
            dt /= cfl_factor

    # 6. Snapshot configuration — clamp interval to dt for leapfrog,
    # since the solver can't save more often than once per timestep.
    duration = args.t_end - t_start
    # Default to 20 snapshots across the duration. Empirically this resolves
    # P_max and other peak-finding measurements to well under 0.1% for the
    # physics of interest (e.g. dark photon conversion: P_max matches 100-
    # snapshot runs to 4 significant figures).  Use --snapshots to specify
    # a smaller interval when finer time resolution is needed.
    snapshot_interval = (
        args.snapshots if args.snapshots is not None else duration / 20.0
    )
    if dt is not None and snapshot_interval < dt:
        log(
            f"  Note: snapshot interval {snapshot_interval:.4f} < dt {dt:.4f}; "
            f"saving every step"
        )
        snapshot_interval = dt

    num_snapshots = max(int(duration / snapshot_interval) + 1, 2)
    # When snapshot_interval ≈ dt (i.e., saving every step), ceil() in the
    # leapfrog step count can exceed floor() in snapshot count by 1.  Add
    # a safety margin to prevent writer overflow.
    if dt is not None:
        n_steps_est = max(1, math.ceil(duration / dt - 1e-10))
        num_snapshots = max(num_snapshots, n_steps_est + 2)

    # 7. Disk writer (if directory output) or in-memory accumulator
    # (inference path: skip disk entirely, see issue #269).
    fmt = _infer_output_format(args)
    writer: SnapshotWriter | None = None
    accumulator: Any = None
    snapshot_cb: Callable[[float, np.ndarray], None] | None = None

    if in_memory_out is not None:
        accumulator, snapshot_cb = _setup_memory_accumulator_native(
            spec,
            grid_info,
            params,
            snapshot_interval,
            dt=dt,
            num_snapshots=num_snapshots,
        )
    elif fmt == "directory":
        writer, snapshot_cb = _setup_disk_writer_native(
            args,
            spec,
            grid_info,
            params,
            snapshot_interval,
            dt=dt,
            num_snapshots=num_snapshots,
        )

    # 8. Progress bar (suppressed by --quiet and non-TTY stderr)
    from tidal.solver.progress import SimulationProgress

    progress: SimulationProgress | None = None
    if not args.quiet:
        solver_labels = {
            "ida": "IDA",
            "cvode": "CVODE",
            "scipy": "scipy",
            "leapfrog": "leapfrog",
        }
        progress = SimulationProgress(
            t_start,
            args.t_end,
            solver_name=solver_labels.get(scheme, scheme),
        )

    # 9. Solve
    if scheme == "ida":
        from tidal.solver.ida import solve_ida

        log(
            f"Running IDA solver (t={t_start} → {args.t_end}, {num_snapshots} snapshots, "
            f"rtol={args.rtol:.0e}, atol={args.atol:.0e})..."
        )
        # Skip constraint IC solving when resuming (state already consistent)
        allow_inconsistent = getattr(args, "allow_inconsistent_ic", False)
        if resume_state is not None:
            allow_inconsistent = True
        result = solve_ida(
            spec,
            grid_info,
            y0,
            t_span=(t_start, args.t_end),
            bc=bc,
            parameters=params,
            num_snapshots=num_snapshots,
            rtol=args.rtol,
            atol=args.atol,
            snapshot_callback=snapshot_cb,
            allow_inconsistent_ic=allow_inconsistent,
            progress=progress,
        )
    elif scheme == "cvode":
        from tidal.solver.cvode import solve_cvode

        method = args.method or "BDF"
        max_step = args.max_step or 0.0
        log(
            f"Running CVODE solver ({method}, t={t_start} → {args.t_end}, "
            f"rtol={args.rtol:.0e}, atol={args.atol:.0e})..."
        )
        result = solve_cvode(
            spec,
            grid_info,
            y0,
            t_span=(t_start, args.t_end),
            bc=bc,
            parameters=params,
            method=method,
            rtol=args.rtol,
            atol=args.atol,
            max_step=max_step,
            num_snapshots=num_snapshots,
            snapshot_callback=snapshot_cb,
            progress=progress,
        )
    elif scheme == "scipy":
        from tidal.solver.scipy_solver import solve_scipy

        method = args.method or "DOP853"
        cfl_dt = _compute_cfl_dt(spec, grid_info, params)
        max_step = args.max_step if args.max_step is not None else cfl_dt
        log(
            f"Running scipy solver ({method}, t={t_start} → {args.t_end}, "
            f"rtol={args.rtol:.0e}, atol={args.atol:.0e})..."
        )
        result = solve_scipy(
            spec,
            grid_info,
            y0,
            t_span=(t_start, args.t_end),
            bc=bc,
            parameters=params,
            method=method,
            rtol=args.rtol,
            atol=args.atol,
            max_step=max_step,
            num_snapshots=num_snapshots,
            snapshot_callback=snapshot_cb,
            progress=progress,
        )
    elif scheme == "modal":
        # Default --perturbative-order: 1 when the JSON declares a
        # perturbation block, 0 otherwise. An explicit --perturbative-order
        # flag on the CLI overrides this default.  Use getattr so callers
        # built from leaner parsers (e.g. `tidal sample`) still work.
        pert_meta: dict[str, Any] = spec.metadata.get("perturbation") or {}
        pert_order_arg = getattr(args, "perturbative_order", None)
        if pert_order_arg is not None:
            pert_order = int(pert_order_arg)
        else:
            pert_order = 1 if pert_meta.get("small_parameters") else 0

        if pert_order > 0 and spec.has_corrections():
            from tidal.solver.modal import solve_modal
            from tidal.solver.perturbative_driver import (
                PerturbativeSolver,
            )

            log(
                f"Running perturbative modal solver (order={pert_order}, "
                f"t={t_start} → {args.t_end}, {num_snapshots} snapshots)..."
            )
            pert_solver = PerturbativeSolver(spec)
            pert_result = pert_solver.solve(
                y0,
                grid_info,
                t_span=(t_start, args.t_end),
                order=pert_order,
                parameters=params,
                num_snapshots=num_snapshots,
                small_parameters=list(pert_meta.get("small_parameters") or []),
            )
            result = cast("SolverResult", pert_result.total)
            # Pass 0 / Pass 1 outputs live in base_spec's layout (h_4/h_7/h_9
            # demoted to algebraic constraints at ε=0). The solver now
            # carries this layout on the result itself (#276), so
            # downstream SimulationData / measurement callers pair the
            # state vector with the correct spec without the prior
            # `spec = pert_solver.base_spec` rebind workaround. The full
            # original spec stays on pert_result.full_spec.
            assert pert_result.spec is not None, (
                "PerturbativeResult.spec must be populated. See #276."
            )
            spec = pert_result.spec
            log(
                f"Perturbative layout swap: full-spec fields={len(pert_result.full_spec.equations) if pert_result.full_spec else '?'} "
                f"→ base-spec fields={len(spec.equations)} "
                f"(demoted constraints: "
                f"{sum(1 for e in spec.equations if e.time_derivative_order == 0)})"
            )
            validity = pert_result.validity
            corr_level = validity.get("correction_level", validity.get("warn_level"))
            base_level = validity.get("base_level", "ok")
            if corr_level == "warn":
                _cwarn(
                    "Perturbative truncation validity: "
                    f"ε·ω·t_end ≈ {validity['validity_param']:.2g} > 0.1 "
                    f"(dominant parameter: {validity['dominant_parameter']!r}). "
                    "O(ε²) truncation error may exceed 10%. Consider a "
                    "smaller small-parameter, shorter t_end, or coarser grid.",
                )
            elif corr_level == "error":
                _cwarn(
                    "Perturbative truncation validity: "
                    f"ε·ω·t_end ≈ {validity['validity_param']:.2g} > 1.0 — "
                    "the EFT regime is violated for this configuration. "
                    "Results should NOT be trusted. Reduce the small "
                    "parameter or shrink t_end / k_max.",
                )
            # #285: separate base-theory stability diagnostic.
            if base_level == "warn":
                _cwarn(
                    "Base-theory stability: "
                    f"max Re(λ)·t_end ≈ {validity['base_stability_param']:.2g} > 0 "
                    f"(dominant: {validity.get('dominant_tachyon_field')!r}). "
                    "The Pass 0 solution contains an exponentially growing "
                    "mode (tachyon / Jeans / ghost). Pass 1 corrections on "
                    "top of a diverging Pass 0 are not physically meaningful.",
                )
            elif base_level == "error":
                _cwarn(
                    "Base-theory stability: "
                    f"max Re(λ)·t_end ≈ {validity['base_stability_param']:.2g} > 30 "
                    f"(dominant: {validity.get('dominant_tachyon_field')!r}). "
                    "Pass 0 will overflow double precision within the "
                    "simulation window. Shorten t_end or fix the base-"
                    "theory mass spectrum before re-running.",
                )
        else:
            from tidal.solver.modal import solve_modal

            log(
                f"Running modal solver (t={t_start} → {args.t_end}, "
                f"{num_snapshots} snapshots)..."
            )
            result = solve_modal(
                spec,
                grid_info,
                y0,
                t_span=(t_start, args.t_end),
                bc=bc,
                parameters=params,
                rtol=args.rtol,
                atol=args.atol,
                num_snapshots=num_snapshots,
                snapshot_callback=snapshot_cb,
                progress=progress,
            )
    else:  # leapfrog
        from tidal.solver.leapfrog import solve_leapfrog

        assert dt is not None  # computed in step 5
        if lf_order == 4:  # noqa: PLR2004
            log(
                f"Running Yoshida 4th-order leapfrog "
                f"(t={t_start} → {args.t_end}, dt={dt:.4f})..."
            )
        else:
            log(f"Running leapfrog solver (t={t_start} → {args.t_end}, dt={dt:.4f})...")
        result = solve_leapfrog(
            spec,
            grid_info,
            y0,
            t_span=(t_start, args.t_end),
            dt=dt,
            bc=bc,
            parameters=params,
            snapshot_interval=snapshot_interval,
            snapshot_callback=snapshot_cb,
            progress=progress,
            order=lf_order,
        )

    if not result["success"]:
        _cerror_hint(
            f"solver failed: {result['message']}",
            [
                "Try coarser grid (`--grid-shape 32`)",
                "Increase tolerances (`--rtol 1e-5`)",
                "Try different solver (`--scheme scipy`)",
            ],
        )
        return 1

    # Post-simulation divergence check: verify final state is finite.
    # This is a single check at the end — zero per-step overhead.
    _check_result_finite(result)

    if writer is not None:
        writer.close()
        log(f"  {writer.count} snapshots streamed to: {writer.output_dir.resolve()}")

        # Save constraint velocities from modal solver (exact ∂_t for constraints).
        # "Constraint" is a solver concept — constraint fields have physical
        # velocities determined by coupling to dynamical fields.
        cv = result.get("constraint_velocities")
        if cv:
            for c_name, c_vel_arr in cv.items():
                np.save(
                    str(writer.output_dir / f"v_{c_name}.npy"),
                    np.asarray(c_vel_arr, dtype=np.float64),
                )
            log(f"  {len(cv)} constraint velocity arrays saved")

    # 8. Build SimulationData — use memory-mapped directory reader when
    # snapshots were already streamed to disk (avoids double-buffering).
    if writer is not None:
        sim_data = SimulationData.from_directory(writer.output_dir, spec)
    elif accumulator is not None:
        accumulator.close()
        # Inject constraint velocities from modal solver (same semantics as
        # the disk path's `np.save` block above, but into the in-memory
        # SimulationData.velocities dict).
        cv = result.get("constraint_velocities")
        if cv:
            for c_name, c_vel_arr in cv.items():
                accumulator.set_velocity(
                    c_name, np.asarray(c_vel_arr, dtype=np.float64)
                )
        sim_data = accumulator.to_sim_data(spec)
    else:
        sim_data = SimulationData.from_result(result, spec, grid_info, params, dt=dt)
    log(f"  {sim_data.n_snapshots} snapshots stored")

    # 9. Output — in-memory mode skips plotting, report generation, and
    # the output dispatch entirely.  Caller consumes sim_data via
    # in_memory_out.
    if in_memory_out is not None:
        in_memory_out.append(sim_data)
        return 0

    _generate_output(args, sim_data, grid_info)

    # 10. HTML report (optional)
    report_path = getattr(args, "report", None)
    if report_path:
        from tidal.cli._report import generate_report

        generate_report(
            sim_data=sim_data,
            spec=spec,
            params=params,
            grid_info=grid_info,
            scheme=scheme,
            report_path=report_path,
        )
        log(f"  Report saved to: {Path(report_path).resolve()}")

    return 0


# --- Command entry point ---


def _validate_solver_params(args: Namespace) -> None:
    """Validate solver-related CLI arguments.

    Raises
    ------
    ValueError
        If ``--t-end``, ``--dt``, ``--snapshots``, ``--rtol``, ``--atol``,
        or ``--max-step`` are non-positive, or if tolerance flags are used
        with leapfrog.
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
    if args.rtol <= 0:
        msg = f"--rtol must be positive, got {args.rtol}"
        raise ValueError(msg)
    if args.atol <= 0:
        msg = f"--atol must be positive, got {args.atol}"
        raise ValueError(msg)
    if args.max_step is not None and args.max_step <= 0:
        msg = f"--max-step must be positive, got {args.max_step}"
        raise ValueError(msg)


def simulate_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0915
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
    if getattr(args, "list_schemes", False):
        print("Available solver schemes:")
        print("  auto      Auto-select based on equation structure (default)")
        print("  modal     Fourier modal solver (periodic, time-independent)")
        print("  cvode     SUNDIALS CVODE (adaptive ODE, tolerance-controlled)")
        print("  ida       SUNDIALS IDA (DAE, algebraic constraints)")
        print("  leapfrog  Symplectic leapfrog (exact energy conservation)")
        print("  scipy     SciPy solve_ivp (DOP853, Radau, BDF)")
        return 0

    if args.json_path is None:
        _cerror("json_path is required")
        return 1

    from tidal.symbolic import load_equation_system

    json_path = Path(args.json_path)
    if not json_path.exists():
        _cerror_hint(
            f"file not found: {json_path}",
            [
                "Run 'tidal list' to see available specifications.",
                "Run 'tidal derive <theory.toml>' to generate from a Lagrangian.",
            ],
        )
        return 1

    log = _clog

    # Step 1: Load spec
    log("Loading equation specification...")
    spec = load_equation_system(json_path)
    log(
        f"  {spec.n_components} component(s), {spec.dimension}D ({spec.spatial_dimension}+1D)"
    )

    # Step 2: Validate resume args and inherit config from checkpoint
    resume_dir: Path | None = None
    resume_meta: dict[str, Any] | None = None
    if args.resume is not None:
        import json as _json

        resume_dir = Path(args.resume)
        if not resume_dir.is_dir():
            _cerror_hint(
                f"resume directory not found: {resume_dir}",
                [
                    "Create output with `tidal simulate ... --output DIR`, then `--resume DIR`",
                ],
            )
            return 1

        # --resume and --ic are mutually exclusive.  argparse default is
        # "gaussian", so a non-gaussian value means the user explicitly set --ic.
        if args.ic != "gaussian":
            _cerror_hint(
                "--resume and --ic cannot be used together",
                ["`--resume` loads IC from checkpoint. Omit `--ic` to resume."],
            )
            return 1

        meta_path = resume_dir / "metadata.json"
        if meta_path.exists():
            with meta_path.open(encoding="utf-8") as f:
                resume_meta = _json.load(f)

            # Inherit grid config if not explicitly provided
            meta: dict[str, Any] = resume_meta  # type: ignore[assignment]
            if args.grid_shape is None:
                grid_shape = cast("list[int]", meta["grid_shape"])
                args.grid_shape = ",".join(str(s) for s in grid_shape)
            if args.bounds is None:
                grid_bounds = cast("list[list[float]]", meta["grid_bounds"])
                args.bounds = ",".join(f"{b[0]}:{b[1]}" for b in grid_bounds)
            if args.bc is None and "bc_types" in meta:
                bc_types = cast("list[str]", meta["bc_types"])
                args.bc = ",".join(bc_types)
            log(f"  Resuming from: {resume_dir}")

    if args.snapshot is not None and args.resume is None:
        _cerror_hint(
            "--snapshot requires --resume",
            ["Usage: `--resume DIR --snapshot 3`"],
        )
        return 1

    # Step 3: Parse parameters (merge with checkpoint metadata if resuming)
    params = _parse_params(args.param, spec)
    if resume_meta is not None:
        saved_params: dict[str, float] = resume_meta.get("parameters", {})  # type: ignore[assignment]
        # Saved params as defaults, CLI --param overrides
        params = {**saved_params, **params}
    if params:
        log(f"  Parameters: {params}")

    # Check output directory collision
    if args.output and Path(args.output).exists() and not getattr(args, "force", False):
        _cerror_hint(
            f"output directory already exists: {args.output}",
            ["Use --force to overwrite", "Or choose a different --output path"],
        )
        return 1

    # Step 3b: previously called `normalize_kinetic_coefficients` to work
    # around the modal solver's hardcoded `M_mat[fi,fi]=1.0` assumption.
    # That root cause is now fixed in `_build_evolution_matrices`, which
    # reads `kinetic_coefficient_symbolic` directly into the M diagonal
    # and uses SVD to handle asymmetric M.  The normalize step is
    # therefore redundant and has been removed — keeping it would
    # produce an equivalent but different-looking equation system that
    # some downstream code (energy measurement, constraint IC solver)
    # may not handle correctly.

    # All simulation goes through the native IDA/leapfrog path
    try:
        return _simulate(args, spec, params)
    except Exception as exc:
        from tidal.solver._exceptions import SimulationDivergedError

        if isinstance(exc, SimulationDivergedError):
            _cerror_hint(
                str(exc),
                [
                    "Reduce `--ic-amplitude` or check `--param` values",
                    "Try smaller timestep with `--dt`",
                ],
            )
            return 1
        raise
