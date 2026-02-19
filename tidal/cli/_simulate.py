"""``tidal simulate`` — Run PDE simulation from a JSON specification."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

    from pde import CartesianGrid, FieldCollection, MemoryStorage

    from tidal.measurement._io import SimulationData
    from tidal.measurement._writer import SnapshotWriter
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


def parse_params(raw: list[str], spec: EquationSystem) -> dict[str, float]:  # noqa: C901
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


def validate_formula_ast(expr: str, allowed_names: set[str]) -> None:
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
    validate_formula_ast(args.ic_formula, allowed_names)

    field_arr = eval(args.ic_formula, {"__builtins__": {}}, namespace)  # noqa: S307
    field_arr = np.asarray(field_arr, dtype=float)

    # Broadcast scalar results to grid shape
    if field_arr.shape == ():
        field_arr = np.full(grid.shape, float(field_arr))

    return create_initial_state(grid, spec, field_data={component: field_arr})


def _validate_ic_args(
    ic_type: str,
    component: str,
    spec: EquationSystem,
    args: Namespace,
) -> None:
    """Emit warnings for problematic IC argument combinations.

    Warns when IC targets a constraint field or when arguments are
    silently ignored for a given IC type.
    """
    import warnings

    # Warn if IC targets a constraint field (time_derivative_order == 0)
    if ic_type != "zero":
        eq = next((e for e in spec.equations if e.field_name == component), None)
        if eq is not None and eq.time_derivative_order == 0:
            dynamical = [
                e.field_name
                for e in spec.equations
                if e.time_derivative_order >= 2  # noqa: PLR2004
            ]
            suggestion = dynamical[0] if dynamical else "a dynamical field"
            warnings.warn(
                f"IC applied to constraint field '{component}' "
                f"(time_derivative_order=0). "
                f"Constraint solver will overwrite this IC. "
                f"Consider --ic-component {suggestion} instead.",
                UserWarning,
                stacklevel=3,
            )

    # Warn about silently ignored --ic-wavevector with gaussian IC
    if ic_type == "gaussian" and args.ic_wavevector is not None:
        warnings.warn(
            "--ic-wavevector is ignored for --ic=gaussian "
            "(only applies to --ic=plane-wave)",
            UserWarning,
            stacklevel=3,
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

    _validate_ic_args(ic_type, component, spec, args)

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
    # Default: disk-backed directory output so simulation data always persists.
    # Users can override with --output foo.png for in-memory image output.
    return "directory"


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


# --- Adaptive solver helpers ---

_IMPLICIT_METHODS = frozenset({"Radau", "BDF"})


def _format_solver_log(args: Namespace, dt: float, max_step: float | None) -> str:
    """Format the solver startup log message."""
    parts = [f"Running simulation (t=0 → {args.t_end}"]
    if args.scheme == "scipy":
        method = args.method or "DOP853"
        parts.append(f"method={method}")
        if args.rtol is not None:
            parts.append(f"rtol={args.rtol:.0e}")
        if args.atol is not None:
            parts.append(f"atol={args.atol:.0e}")
        parts.append(f"first_step={dt:.4f}")
    elif args.adaptive:
        tol = args.tolerance if args.tolerance is not None else 1e-4
        parts.extend((f"adaptive RK, tolerance={tol:.0e}", f"initial dt={dt:.4f}"))
    else:
        parts.append(f"dt={dt:.4f}, scheme={args.scheme}")
    if max_step is not None:
        parts.append(f"max_step={max_step:.4f}")
    return ", ".join(parts) + ")..."


def _build_scipy_kwargs(
    args: Namespace,
    max_step: float | None,
    pde: object,
    grid: CartesianGrid,
) -> dict[str, Any]:
    """Build kwargs to forward to ``solve_ivp`` via py-pde's ScipySolver."""
    kwargs: dict[str, Any] = {"method": args.method or "DOP853"}
    if args.rtol is not None:
        kwargs["rtol"] = args.rtol
    if args.atol is not None:
        kwargs["atol"] = args.atol
    if max_step is not None:
        kwargs["max_step"] = max_step
    # For implicit methods, compute Jacobian sparsity for performance
    method = args.method or "DOP853"
    if method in _IMPLICIT_METHODS:
        jac_sp = pde.jacobian_sparsity(grid)  # type: ignore[attr-defined]
        if jac_sp is not None:
            kwargs["jac_sparsity"] = jac_sp
        elif grid.dim > 1:
            import warnings as _warnings

            _warnings.warn(
                f"Implicit method {method} on {grid.dim}D grid: "
                f"Jacobian sparsity not available, falling back to dense "
                f"estimation. This may be very slow for large grids.",
                stacklevel=2,
            )
    return kwargs


def _build_energy_monitor(
    threshold: float,
    snapshot_interval: float,
) -> object:
    """Build a CallbackTracker that monitors energy conservation.

    Uses the L2 norm of the full state vector as an energy proxy.
    This is conserved for Hamiltonian systems (up to discretization error)
    and does not require knowledge of the equation structure.

    Raises RuntimeError during simulation if ``|dE/E0|`` exceeds *threshold*.
    """
    from pde import CallbackTracker

    energy_epsilon = 1e-30
    e0_holder: list[float | None] = [None]  # mutable container for closure

    def _energy_callback(state_now: object, t: float) -> None:
        data = state_now.data  # type: ignore[union-attr]
        e_total = float(np.sum(data * data))  # pyright: ignore[reportUnknownArgumentType]
        if e0_holder[0] is None:
            e0_holder[0] = e_total
            return
        e0 = e0_holder[0]
        if abs(e0) < energy_epsilon:
            return  # can't compute relative drift from ~zero energy
        drift = abs((e_total - e0) / e0)
        if drift > threshold:
            msg = (
                f"Energy conservation violated at t={t:.4f}: "
                f"|dE/E0|={drift:.3e} > threshold={threshold:.3e}"
            )
            raise RuntimeError(msg)

    return CallbackTracker(_energy_callback, interrupts=snapshot_interval)


# --- Command entry point ---


def _check_cfl_stability(pde: object, dt: float, grid: CartesianGrid) -> None:
    """Print CFL stability warnings to stderr."""
    warnings_list = cast("list[str]", pde.check_stability(dt, grid))  # type: ignore[attr-defined]
    sys.stderr.writelines(f"  Warning: {w}\n" for w in warnings_list)


def _validate_solver_params(args: Namespace) -> None:  # noqa: C901, PLR0912
    """Validate solver-related CLI arguments.

    Raises
    ------
    ValueError
        If ``--t-end``, ``--dt``, or ``--snapshots`` are non-positive,
        or if adaptive-solver flags are combined incorrectly.
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

    # Adaptive-solver flag validation
    if args.method is not None and args.scheme != "scipy":
        msg = "--method requires --scheme scipy"
        raise ValueError(msg)
    if args.rtol is not None and args.scheme != "scipy":
        msg = "--rtol requires --scheme scipy"
        raise ValueError(msg)
    if args.atol is not None and args.scheme != "scipy":
        msg = "--atol requires --scheme scipy"
        raise ValueError(msg)
    if args.tolerance is not None and (
        args.scheme != "runge-kutta" or not args.adaptive
    ):
        msg = "--tolerance requires --scheme runge-kutta --adaptive"
        raise ValueError(msg)
    if args.rtol is not None and args.rtol <= 0:
        msg = f"--rtol must be positive, got {args.rtol}"
        raise ValueError(msg)
    if args.atol is not None and args.atol <= 0:
        msg = f"--atol must be positive, got {args.atol}"
        raise ValueError(msg)
    if args.tolerance is not None and args.tolerance <= 0:
        msg = f"--tolerance must be positive, got {args.tolerance}"
        raise ValueError(msg)
    if args.max_step is not None and args.max_step <= 0:
        msg = f"--max-step must be positive, got {args.max_step}"
        raise ValueError(msg)
    if args.max_step is not None and args.scheme != "scipy":
        msg = "--max-step requires --scheme scipy (py-pde's ExplicitSolver does not support max_step)"
        raise ValueError(msg)
    if args.energy_monitor is not None and args.energy_monitor <= 0:
        msg = f"--energy-monitor must be positive, got {args.energy_monitor}"
        raise ValueError(msg)


def _stiffness_advisory(  # noqa: C901, PLR0912
    spec: EquationSystem,
    grid: CartesianGrid,
    args: Namespace,
    log: object,
) -> None:
    """Log an advisory if the system appears stiff and an explicit method is chosen."""
    method = args.method or ("DOP853" if args.scheme == "scipy" else "runge-kutta")
    if method in _IMPLICIT_METHODS:
        return  # already using an implicit method

    stiffness_threshold = 100
    dx_min = min(grid.discretization)

    # --- Mass-based stiffness check ---
    mass = spec.mass_matrix
    if mass:
        try:
            max_mass_sq = max(
                abs(mass[i][j]) for i in range(len(mass)) for j in range(len(mass[0]))
            )
        except (TypeError, IndexError):
            max_mass_sq = 0.0
        if max_mass_sq > 0:
            max_c_sq = 0.0
            for eq in spec.equations:
                if eq.time_derivative_order < 2:  # noqa: PLR2004
                    continue
                for term in eq.rhs_terms:
                    if term.operator == "laplacian" or term.operator.startswith(
                        "laplacian_"
                    ):
                        max_c_sq = max(max_c_sq, abs(term.coefficient))
            if max_c_sq > 0:
                stiffness_ratio = max_mass_sq * dx_min**2 / max_c_sq
                if stiffness_ratio > stiffness_threshold:
                    msg = (
                        f"  Note: system may be stiff "
                        f"(m²·dx²/c²={stiffness_ratio:.0f}). "
                        f"Consider --scheme scipy --method Radau "
                        f"for better performance."
                    )
                    cast("object", log)(msg)  # type: ignore[operator]

    # --- Anisotropic Laplacian coefficient spread ---
    # If directional Laplacians (laplacian_x, laplacian_y, ...) have
    # very different coefficients, the fastest axis dominates CFL while
    # the slowest axis wastes steps — a form of stiffness.
    dir_coeffs: list[float] = []
    for eq in spec.equations:
        if eq.time_derivative_order < 2:  # noqa: PLR2004
            continue
        dir_coeffs.extend(
            abs(term.coefficient)
            for term in eq.rhs_terms
            if term.operator.startswith("laplacian_")
        )
    if len(dir_coeffs) >= 2:  # noqa: PLR2004
        min_c = min(dir_coeffs)
        max_c = max(dir_coeffs)
        if min_c > 0 and max_c / min_c > stiffness_threshold:
            msg = (
                f"  Note: anisotropic Laplacian stiffness "
                f"(max/min={max_c / min_c:.0f}). "
                f"Consider --scheme scipy --method Radau."
            )
            cast("object", log)(msg)  # type: ignore[operator]


def simulate_command(args: Namespace) -> int:  # noqa: C901, PLR0912, PLR0914, PLR0915
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
    from tidal.symbolic import build_pde_from_json, load_equation_system

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
    params = parse_params(args.param, spec)
    if params:
        log(f"  Parameters: {params}")

    # Step 3: Build PDE
    log("Building PDE...")
    pde = build_pde_from_json(json_path, parameters=params)

    # Step 4: Create grid
    bounds = _parse_bounds(args.bounds, spec.spatial_dimension)
    grid = _build_grid(args, spec, bounds)
    log(f"  Grid: {'x'.join(str(s) for s in grid.shape)}, bounds: {grid.axes_bounds}")

    # Stiffness advisory: warn if system appears stiff and explicit method is chosen
    _stiffness_advisory(spec, grid, args, log)

    # Step 5: Initial conditions
    state = _build_initial_state(args, grid, spec, bounds)
    initial_state = state.copy()
    log(f"  IC: {args.ic} on {args.ic_component or spec.component_names[0]}")

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

    from tidal.utils import normalize_solve_result

    # Compute CFL-based dt hint and max_step guard
    cfl_dt = _CFL_FACTOR * min(
        (b[1] - b[0]) / s for b, s in zip(grid.axes_bounds, grid.shape, strict=True)
    )
    cfl_max_step = pde.cfl_limit(grid)  # type: ignore[attr-defined]

    # Determine dt and max_step based on solver path
    if args.scheme == "scipy":
        # Scipy: dt is just a first_step hint; tolerances drive step selection
        dt: float = args.dt if args.dt is not None else cfl_dt
        max_step: float | None = (
            args.max_step if args.max_step is not None else cfl_max_step
        )
    elif args.adaptive:
        # Adaptive explicit RK: dt is initial guess, adapted by tolerance
        dt = args.dt if args.dt is not None else cfl_dt
        max_step = None  # py-pde ExplicitSolver has no max_step support
    else:
        # Fixed-step explicit RK (original behavior)
        dt = args.dt if args.dt is not None else cfl_dt
        max_step = None
        _check_cfl_stability(pde, dt, grid)

    # Step 7: Run simulation
    snapshot_interval = (
        args.snapshots if args.snapshots is not None else args.t_end / 100.0
    )
    log(_format_solver_log(args, dt, max_step))

    fmt = _infer_output_format(args)
    use_directory = fmt == "directory"
    writer: SnapshotWriter | None = None

    if use_directory:
        # Disk-backed streaming: O(1) memory regardless of snapshot count
        storage, tracker, writer = _setup_disk_backed(
            args,
            spec,
            grid,
            snapshot_interval,
            params,
        )
        log(f"  Output directory: {writer.output_dir}")
    else:
        # In-memory path (for plot output)
        from pde import MemoryStorage

        storage = MemoryStorage()
        tracker = storage.tracker(snapshot_interval)

    # Add blow-up detection callback
    initial_max = float(np.max(np.abs(initial_state.data)))
    blowup_threshold = min(max(initial_max * 1e6, 1e6), 1e15)

    def _check_blowup(state_now: object, t: float) -> None:
        current_max = float(np.max(np.abs(state_now.data)))  # type: ignore[union-attr]
        if current_max > blowup_threshold:
            msg = (
                f"Blow-up detected at t={t:.4f}: max|state|={current_max:.3e} "
                f"exceeds threshold {blowup_threshold:.3e}"
            )
            raise RuntimeError(msg)

    from pde import CallbackTracker

    blowup_tracker = CallbackTracker(_check_blowup, interrupts=snapshot_interval)
    if isinstance(tracker, list):
        tracker.append(blowup_tracker)
    else:
        tracker = [tracker, blowup_tracker]

    # Add energy monitor callback if requested
    if args.energy_monitor is not None:
        energy_tracker = _build_energy_monitor(
            args.energy_monitor,
            snapshot_interval,
        )
        tracker.append(energy_tracker)  # pyright: ignore[reportArgumentType]
        # Implicit methods (BDF, Radau) introduce numerical dissipation,
        # which causes the L2 norm to drift even for stable systems.
        # Warn the user to use a generous threshold.
        method = args.method or ("DOP853" if args.scheme == "scipy" else "")
        if method in _IMPLICIT_METHODS:
            log(
                f"  Note: --energy-monitor with implicit method {method}: "
                f"expect L2 norm drift from numerical dissipation. "
                f"Use threshold >= 0.1 to avoid false alarms."
            )

    import time as _time

    t_wall_start = _time.perf_counter()

    # Solve — three paths: fixed explicit, adaptive explicit, scipy adaptive
    if args.scheme == "scipy":
        solver_kwargs = _build_scipy_kwargs(args, max_step, pde, grid)
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]  # pyright: ignore[reportCallIssue]
                state,
                t_range=args.t_end,
                dt=dt,
                solver="scipy",
                tracker=tracker,  # pyright: ignore[reportArgumentType]
                **solver_kwargs,
            )
        )
    elif args.adaptive:
        adaptive_kwargs: dict[str, Any] = {"adaptive": True}
        if args.tolerance is not None:
            adaptive_kwargs["tolerance"] = args.tolerance
        normalize_solve_result(
            pde.solve(  # type: ignore[attr-defined]  # pyright: ignore[reportCallIssue]
                state,
                t_range=args.t_end,
                dt=dt,
                scheme=args.scheme,
                tracker=tracker,  # pyright: ignore[reportArgumentType]
                **adaptive_kwargs,
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

    t_wall_elapsed = _time.perf_counter() - t_wall_start
    log(f"  Completed in {t_wall_elapsed:.2f}s")

    if writer is not None:
        writer.close()
        log(f"  {writer.count} snapshots streamed to: {writer.output_dir}")
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

    if args.output:
        output_dir = Path(args.output)
    else:
        # Auto-generate from spec name: examples/data/foo.json → examples/data/foo_output/
        json_file = Path(args.json_path).resolve()
        output_dir = json_file.parent / f"{json_file.stem}_output"

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
