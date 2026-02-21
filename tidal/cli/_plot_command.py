"""``tidal plot`` — Generate individual plots from simulation output.

Reads disk-backed simulation output (from ``tidal simulate --output``)
and produces a single focused plot per invocation.  Users compose
what they need via multiple calls in shell scripts.

Plot types: heatmap, snapshot, amplitude, energy, profile, compare.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

_VALID_TYPES = frozenset(
    {"heatmap", "snapshot", "amplitude", "energy", "profile", "compare"}
)

DPI_DEFAULT = 150


# ------------------------------------------------------------------
# Argument helpers
# ------------------------------------------------------------------


def _parse_fields(raw: str | None) -> list[str] | None:
    """Parse comma-separated ``--fields`` into a list, or *None* if absent."""
    if raw is None:
        return None
    return [s.strip() for s in raw.split(",") if s.strip()]


def _parse_time_indices(raw: str | None) -> list[int] | None:
    """Parse comma-separated ``--time-indices`` into a list of ints.

    Raises
    ------
    ValueError
        If any element is not a valid integer.
    """
    if raw is None:
        return None
    try:
        return [int(s.strip()) for s in raw.split(",")]
    except ValueError as exc:
        msg = f"Invalid --time-indices: {raw!r} (must be comma-separated integers)"
        raise ValueError(msg) from exc


def _parse_cross_section(raw: str | None) -> tuple[str, float] | None:
    """Parse ``--cross-section AXIS=VAL`` into (axis_name, value).

    Raises
    ------
    ValueError
        If the format is invalid or the axis is not x, y, or z.
    """
    if raw is None:
        return None
    if "=" not in raw:
        msg = f"Invalid --cross-section: {raw!r} (expected AXIS=VAL, e.g. y=25.0)"
        raise ValueError(msg)
    axis, val_str = raw.split("=", 1)
    axis = axis.strip().lower()
    if axis not in {"x", "y", "z"}:
        msg = f"Invalid cross-section axis '{axis}', expected x, y, or z"
        raise ValueError(msg)
    try:
        val = float(val_str.strip())
    except ValueError as exc:
        msg = f"Invalid cross-section value: {val_str!r} (must be a number)"
        raise ValueError(msg) from exc
    return (axis, val)


def _parse_figsize(raw: str | None) -> tuple[float, float] | None:
    """Parse ``--figsize W,H`` into (width, height).

    Raises
    ------
    ValueError
        If the format is invalid or values are not numbers.
    """
    if raw is None:
        return None
    parts = raw.split(",")
    if len(parts) != 2:  # noqa: PLR2004
        msg = f"Invalid --figsize: {raw!r} (expected W,H)"
        raise ValueError(msg)
    try:
        return (float(parts[0].strip()), float(parts[1].strip()))
    except ValueError as exc:
        msg = f"Invalid --figsize values: {raw!r} (must be numbers)"
        raise ValueError(msg) from exc


def _validate_overlay(formula: str) -> None:
    """Validate an overlay formula using AST analysis."""
    from tidal.cli._simulate import FORMULA_NAMESPACE, _validate_formula_ast

    allowed = set(FORMULA_NAMESPACE.keys()) | {"t"}
    _validate_formula_ast(formula, allowed)


# ------------------------------------------------------------------
# Default output filename
# ------------------------------------------------------------------

_SINGLE_FIELD_TYPES = frozenset({"heatmap", "snapshot", "profile"})


def _default_filename(plot_type: str, args: Namespace) -> str:
    """Build a collision-resistant default filename from plot args.

    Encodes the field name and time index into the filename so that
    multiple ``tidal plot`` invocations with different parameters do
    not overwrite each other.

    Examples: ``heatmap_phi_0.png``, ``snapshot_A_0_t0.png``,
    ``snapshot_A_1_final.png``, ``amplitude.png``.
    """
    parts: list[str] = [plot_type]

    # Include field name for single-field plot types
    field: str | None = getattr(args, "field", None)
    if field is not None and plot_type in _SINGLE_FIELD_TYPES:
        parts.append(field)

    # Include time index for snapshot
    if plot_type == "snapshot":
        time_index: int | None = getattr(args, "time_index", None)
        if time_index is not None and time_index == -1:
            parts.append("final")
        elif time_index is not None:
            parts.append(f"t{time_index}")
        else:
            parts.append("final")

    return "_".join(parts) + ".png"


# ------------------------------------------------------------------
# Main command
# ------------------------------------------------------------------


def plot_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0914, PLR0915
    """Execute the ``tidal plot`` subcommand."""
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    from tidal.cli._measure import _load_data, _resolve_spec_path
    from tidal.cli._panels import (
        field_names,
        render_amplitude,
        render_compare,
        render_energy,
        render_heatmap,
        render_profile,
        render_snapshot,
        resolve_time_indices,
        single_field,
    )

    data_path = Path(args.data_dir)
    if not data_path.is_dir():
        print(f"Error: '{data_path}' is not a directory", file=sys.stderr)
        return 1

    plot_type: str = args.type
    if plot_type not in _VALID_TYPES:
        print(
            f"Error: unknown plot type '{plot_type}'. "
            f"Valid: {', '.join(sorted(_VALID_TYPES))}",
            file=sys.stderr,
        )
        return 1

    # Parse options
    try:
        fields_list = _parse_fields(args.fields)
        time_indices = _parse_time_indices(args.time_indices)
        cross_section = _parse_cross_section(args.cross_section)
        figsize = _parse_figsize(args.figsize)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    # Validate overlay formula if provided
    overlay: str | None = args.overlay
    if overlay is not None:
        try:
            _validate_overlay(overlay)
        except (ValueError, TypeError) as exc:
            print(f"Error in --overlay formula: {exc}", file=sys.stderr)
            return 1

    # Load data
    try:
        spec_path = _resolve_spec_path(data_path, args.spec)
        param_overrides: list[str] = args.param or []
        data = _load_data(data_path, spec_path, param_overrides)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"Error loading data: {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(
            f"Loaded {data.n_snapshots} snapshots, "
            f"fields: {', '.join(data.fields.keys())}"
        )

    # Output path
    output_path = (
        Path(args.output)
        if args.output
        else data_path / _default_filename(plot_type, args)
    )

    # Create figure
    dpi = args.dpi or DPI_DEFAULT
    cmap = args.cmap or "RdBu_r"

    if figsize is not None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig, ax = plt.subplots(1, 1)

    # Dispatch to render function
    try:
        if plot_type == "heatmap":
            field = single_field(data, args.field)
            render_heatmap(ax, data, field, cmap=cmap)

        elif plot_type == "snapshot":
            field = single_field(data, args.field)
            t_idx = args.time_index if args.time_index is not None else -1
            render_snapshot(ax, data, field, t_idx, cmap=cmap)

        elif plot_type == "amplitude":
            fields = field_names(data, fields_list)
            render_amplitude(ax, data, fields, overlay=overlay)

        elif plot_type == "energy":
            fields = field_names(data, fields_list)
            render_energy(ax, data, fields)

        elif plot_type == "profile":
            field = single_field(data, args.field)
            indices = resolve_time_indices(data, time_indices)
            render_profile(ax, data, field, indices, cross_section=cross_section)

        elif plot_type == "compare":
            fields = field_names(data, fields_list)
            render_compare(ax, data, fields, cross_section=cross_section)

    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        plt.close(fig)
        return 1

    # Apply custom title
    if args.title:
        fig.suptitle(args.title)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    if not args.quiet:
        print(f"Saved: {output_path}")

    return 0
