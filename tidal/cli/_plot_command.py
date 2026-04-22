"""``tidal plot`` — Generate individual plots from simulation output.

Reads disk-backed simulation output (from ``tidal simulate --output``)
and produces a single focused plot per invocation.  Users compose
what they need via multiple calls in shell scripts.

Plot types: heatmap, snapshot, amplitude, energy, profile, compare.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

_VALID_TYPES = frozenset(
    {
        "heatmap",
        "snapshot",
        "amplitude",
        "energy",
        "profile",
        "compare",
        "hamiltonian",
        "conservation",
        "sweep",
        "sweep-compare",
        "sweep-parallel",
        "sweep-tornado",
        "sweep-scatter",
        "sweep-grouped",
        "sweep-histogram",
        "sweep-divergence",
        "campaign",
        "convergence",
        "replicate-convergence",
        "corner",
    },
)

_SWEEP_TYPES = frozenset(
    {
        "sweep",
        "sweep-compare",
        "sweep-parallel",
        "sweep-tornado",
        "sweep-scatter",
        "sweep-grouped",
        "sweep-histogram",
        "sweep-divergence",
        "campaign",
        "convergence",
        "replicate-convergence",
    },
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
    from tidal.cli._simulate import (
        FORMULA_NAMESPACE,
        _validate_formula_ast,  # pyright: ignore[reportPrivateUsage]
    )

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


def _corner_plot(data_path: Path, args: Namespace) -> int:
    """Render a post-hoc corner plot from a pulled ``tidal sample`` output.

    Loads via :meth:`InferenceResult.from_directory` (reads
    ``inference.json`` + ``results.csv``) and delegates to
    :func:`plot_corner` which handles anesthetic version drift.
    Resolves issue #288.
    """
    from tidal.cli._console import error_with_hint
    from tidal.inference._results import InferenceResult
    from tidal.inference._visualize import plot_corner

    try:
        result = InferenceResult.from_directory(data_path)
    except (FileNotFoundError, ValueError, OSError) as exc:
        error_with_hint(
            f"loading inference result from {data_path}: {exc}",
            [
                "Directory must contain inference.json + results.csv "
                "(produced by `tidal sample`)",
                "For raw PolyChord chain dirs, use `anesthetic.read_chains` directly",
            ],
        )
        return 1

    # Inject priors supplied via --priors into result.metadata so the
    # corner plotter can overlay them. This lets us replot chains pulled
    # before priors were persisted in inference.json (see #308/#309).
    prior_specs = list(getattr(args, "priors", []) or [])
    if prior_specs:
        from tidal.inference._prior import parse_prior

        parsed = []
        for spec in prior_specs:
            try:
                p = parse_prior(spec)
            except ValueError as exc:
                error_with_hint(
                    f"--priors: {exc}",
                    ["Syntax: NAME=DIST:LO:HI (e.g. mA2=log_uniform:0.001:1.0)"],
                )
                return 1
            parsed.append(
                {
                    "name": p.name,
                    "distribution": p.distribution,
                    "low": p.low,
                    "high": p.high,
                },
            )
        result.metadata["priors"] = parsed

    out_path = Path(args.output) if args.output else data_path / "corner.png"
    plot_corner(result, out_path)
    if not args.quiet:
        print(f"Saved corner plot to: {out_path}")
    return 0


# ------------------------------------------------------------------
# Main command
# ------------------------------------------------------------------


def plot_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0914, PLR0915
    """Execute the ``tidal plot`` subcommand."""
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    from tidal.cli._console import error as _cerror
    from tidal.cli._console import error_with_hint
    from tidal.cli._measure import (
        _load_data,  # pyright: ignore[reportPrivateUsage]
        _resolve_spec_path,  # pyright: ignore[reportPrivateUsage]
    )
    from tidal.cli._panels import (
        field_names,
        render_amplitude,
        render_compare,
        render_conservation,
        render_energy,
        render_hamiltonian,
        render_heatmap,
        render_profile,
        render_snapshot,
        resolve_time_indices,
        single_field,
    )

    data_path = Path(args.data_dir)
    if not data_path.is_dir():
        error_with_hint(
            f"'{data_path}' is not a directory",
            ["Use simulation output directory: `tidal plot output/ --type heatmap`"],
        )
        return 1

    plot_type: str = args.type
    if plot_type not in _VALID_TYPES:
        error_with_hint(
            f"unknown plot type '{plot_type}'. Valid: {', '.join(sorted(_VALID_TYPES))}",
            [
                "Valid: heatmap, snapshot, amplitude, energy, profile, compare, hamiltonian",
            ],
        )
        return 1

    # Campaign: multi-directory composition (handled separately)
    if plot_type == "campaign":
        from tidal.cli._campaign import render_campaign

        return render_campaign(args)

    # Sweep plot types: dispatch to sweep-specific handler
    if plot_type in _SWEEP_TYPES:
        return _sweep_plot(args, data_path, plot_type)

    # Corner plot: load an InferenceResult (not SimulationData) and
    # render via anesthetic.  Accepts the pulled HPC artefact dir
    # directly (the one that contains inference.json + results.csv +
    # _chains/).
    if plot_type == "corner":
        return _corner_plot(data_path, args)

    # Parse options
    try:
        fields_list = _parse_fields(args.fields)
        time_indices = _parse_time_indices(args.time_indices)
        cross_section = _parse_cross_section(args.cross_section)
        figsize = _parse_figsize(args.figsize)
    except ValueError as exc:
        _cerror(str(exc))
        return 1

    # Validate overlay formula if provided
    overlay: str | None = args.overlay
    if overlay is not None:
        try:
            _validate_overlay(overlay)
        except (ValueError, TypeError) as exc:
            error_with_hint(
                f"in --overlay formula: {exc}",
                ["Check syntax. Example: `--overlay 'sin(x)*t'`"],
            )
            return 1

    # Load data
    try:
        spec_path = _resolve_spec_path(data_path, args.spec)
        param_overrides: list[str] = args.param or []
        data = _load_data(data_path, spec_path, param_overrides)
    except (FileNotFoundError, ValueError, OSError) as exc:
        _cerror(f"loading data: {exc}")
        return 1

    if not args.quiet:
        print(
            f"Loaded {data.n_snapshots} snapshots, "
            f"fields: {', '.join(data.fields.keys())}",
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

        elif plot_type == "hamiltonian":
            fields = field_names(data, fields_list)
            render_hamiltonian(ax, data, fields)

        elif plot_type == "conservation":
            threshold = args.threshold if args.threshold is not None else 1e-3
            render_conservation(ax, data, threshold=threshold)

    except ValueError as exc:
        _cerror(str(exc))
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
        print(f"Saved: {output_path.resolve()}")

    return 0


# ------------------------------------------------------------------
# Sweep plot dispatch
# ------------------------------------------------------------------


def _sweep_plot(args: Namespace, data_path: Path, plot_type: str) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0914, PLR0915
    """Handle sweep-specific plot types.

    Loads ``SweepResults`` from *data_path* and dispatches to the
    appropriate render function in ``_sweep_panels``.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    from tidal.cli._console import error as _cerror
    from tidal.cli._console import error_with_hint
    from tidal.cli._sweep_panels import (
        render_convergence,
        render_replicate_convergence,
        render_sweep_1d,
        render_sweep_1d_grouped,
        render_sweep_1d_multi,
        render_sweep_2d,
        render_sweep_2d_with_overlay,
        render_sweep_compare,
        render_sweep_divergence,
        render_sweep_histogram,
        render_sweep_parallel,
        render_sweep_scatter,
        render_sweep_tornado,
    )
    from tidal.measurement._sweep_results import SweepResults

    # Load sweep data
    sweep_json = data_path / "sweep.json"
    if not sweep_json.exists():
        error_with_hint(
            f"'{data_path}' is not a sweep directory (no sweep.json)",
            ["Use `tidal sweep --output` directory for sweep plots"],
        )
        return 1

    try:
        results = SweepResults.from_directory(data_path)
    except (FileNotFoundError, ValueError) as exc:
        _cerror(f"loading sweep data: {exc}")
        return 1

    # Load paired baseline sweep if requested
    baseline_sweep_dir = getattr(args, "baseline_sweep", None)
    if baseline_sweep_dir is not None:
        try:
            baseline_results = SweepResults.from_directory(Path(baseline_sweep_dir))
            results.add_paired_baseline(baseline_results)
        except (FileNotFoundError, ValueError) as exc:
            _cerror(f"loading baseline sweep: {exc}")
            return 1

    # Auto-derive amplification columns when a derived metric is requested
    raw_metric: str | None = getattr(args, "metric", None)
    paired_metrics = {"A_paired", "log10_A_paired", "P_baseline"}
    derived_metrics = {"A", "log10_A", "C0", "P_EM", "P_valid"}
    if raw_metric is not None:
        requested = {s.strip() for s in raw_metric.split(",")}
        if requested & derived_metrics:
            baseline = getattr(args, "baseline_formula", None)
            results.add_derived_columns(
                baseline_formula=baseline or "sin(kappa * B0 * t_end / 2)**2",
            )
        if (requested & paired_metrics) and baseline_sweep_dir is None:
            error_with_hint(
                f"Metric '{raw_metric}' requires --baseline-sweep",
                [
                    "Run a baseline sweep with deltam=0 on the same parameter grid",
                    "Then: tidal plot <signal> --type sweep --metric A_paired "
                    "--baseline-sweep <baseline>",
                ],
            )
            return 1
    if raw_metric is None and plot_type in {
        "sweep",
        "convergence",
        "sweep-parallel",
        "sweep-tornado",
        "sweep-scatter",
        "sweep-grouped",
        "sweep-histogram",
        "replicate-convergence",
    }:
        # Try to auto-detect a sensible metric
        from tidal.measurement._sweep_results import DEFAULT_METRIC_CANDIDATES

        for candidate in DEFAULT_METRIC_CANDIDATES:
            if results.rows and candidate in results.rows[0]:
                raw_metric = candidate
                break
        if raw_metric is None:
            error_with_hint(
                "--metric is required for sweep plots. "
                f"Available: {', '.join(results.metric_names)}",
                ["Example: `--metric P_max`"],
            )
            return 1

    figsize = _parse_figsize(getattr(args, "figsize", None))
    dpi = getattr(args, "dpi", None) or DPI_DEFAULT

    # Colormap / divergent-center options (used by sweep, sweep-parallel, sweep-scatter)
    dc: float | None = getattr(args, "divergent_center", None)
    cmap_name: str = getattr(args, "cmap", None) or (
        "RdBu_r" if dc is not None else "viridis"
    )

    try:
        if plot_type == "sweep":
            metrics = [s.strip() for s in raw_metric.split(",")]  # type: ignore[union-attr]
            n_swept = len(results.swept_params)
            overlay: str | None = getattr(args, "overlay", None)

            log_scale: bool = getattr(args, "log_scale", False)
            log_y: bool = getattr(args, "log_y", False)
            thresholds: list[str] = getattr(args, "hline", []) or []

            if n_swept == 1:
                if len(metrics) == 1:
                    fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
                    use_scatter: bool = getattr(args, "scatter", False)
                    ann_extremes: bool = getattr(args, "annotate_extremes", False)
                    use_arctan: bool = getattr(args, "arctan_axes", False)
                    try:
                        render_sweep_1d(
                            ax,
                            results,
                            metrics[0],
                            overlay=overlay,
                            log_y=log_y,
                            thresholds=thresholds,
                            scatter=use_scatter,
                            annotate_extremes=ann_extremes,
                            arctan_axes=use_arctan,
                        )
                    except ValueError as exc:
                        error_with_hint(
                            f"in --overlay formula: {exc}",
                            ["Check syntax. Example: `--overlay 'sin(x)*t'`"],
                        )
                        return 1
                else:
                    fig = plt.figure(figsize=figsize or (8, 3 * len(metrics)))
                    render_sweep_1d_multi(fig, results, metrics, log_y=log_y)
            elif n_swept == 2:  # noqa: PLR2004
                if overlay:
                    fig = plt.figure(figsize=figsize or (15, 5))
                    try:
                        render_sweep_2d_with_overlay(
                            fig,
                            results,
                            metrics[0],
                            overlay,
                            log_scale=log_scale,
                        )
                    except ValueError as exc:
                        error_with_hint(
                            f"in --overlay formula: {exc}",
                            ["Check syntax. Example: `--overlay 'sin(x)*t'`"],
                        )
                        return 1
                else:
                    clamp_raw: str | None = getattr(args, "clamp_color", None)
                    clamp: tuple[float, float] | None = None
                    if clamp_raw:
                        parts = clamp_raw.split(":")
                        clamp = (float(parts[0]), float(parts[1]))
                    fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 6))
                    render_sweep_2d(
                        ax,
                        results,
                        metrics[0],
                        log_scale=log_scale,
                        divergent_center=dc,
                        cmap_name=cmap_name,
                        clamp_color=clamp,
                    )
            else:
                error_with_hint(
                    f"sweep plot supports 1 or 2 swept parameters, got {n_swept}",
                    ["Use `--type sweep-parallel` for 3+ parameters"],
                )
                return 1

        elif plot_type == "sweep-compare":
            measurement = raw_metric or "conversion"
            fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
            spec_override = getattr(args, "spec", None)
            render_sweep_compare(ax, results, measurement, spec_path=spec_override)

        elif plot_type == "convergence":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for convergence plots",
                    ["Example: `--metric P_max`"],
                )
                return 1
            fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
            render_convergence(ax, results, raw_metric)

        elif plot_type == "sweep-parallel":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for sweep-parallel plots",
                    ["Example: `--metric P_max`"],
                )
                return 1
            fig, ax = plt.subplots(1, 1, figsize=figsize or (10, 6))
            render_sweep_parallel(
                ax,
                results,
                raw_metric,
                cmap_name=cmap_name,
                divergent_center=dc,
            )

        elif plot_type == "sweep-tornado":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for sweep-tornado plots",
                    ["Example: `--metric P_max`"],
                )
                return 1
            fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
            render_sweep_tornado(ax, results, raw_metric)

        elif plot_type == "sweep-scatter":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for sweep-scatter plots",
                    ["Example: `--metric P_max`"],
                )
                return 1
            scatter_params_raw: str | None = getattr(args, "params", None)
            scatter_params: list[str] | None = (
                [s.strip() for s in scatter_params_raw.split(",")]
                if scatter_params_raw
                else None
            )
            log_diag: bool = getattr(args, "log_diagonal", False)
            n_params = len(scatter_params or list(results.swept_params.keys()))
            # constrained_layout handles the colorbar placement cleanly,
            # avoiding the overlap we used to get with tight_layout.
            fig = plt.figure(
                figsize=figsize or (3 * n_params, 3 * n_params),
                constrained_layout=True,
            )
            scatter_clamp_raw: str | None = getattr(args, "clamp_color", None)
            scatter_clamp: tuple[float, float] | None = None
            if scatter_clamp_raw:
                sc_parts = scatter_clamp_raw.split(":")
                scatter_clamp = (float(sc_parts[0]), float(sc_parts[1]))
            render_sweep_scatter(
                fig,
                results,
                raw_metric,
                cmap_name=cmap_name,
                divergent_center=dc,
                log_diagonal=log_diag,
                params=scatter_params,
                clamp_color=scatter_clamp,
            )

        elif plot_type == "sweep-grouped":
            x_param = getattr(args, "x_param", None)
            group_by = getattr(args, "group_by", None)
            if x_param is None or group_by is None:
                error_with_hint(
                    "--x-param and --group-by are required for sweep-grouped plots",
                    [
                        "Example: `--type sweep-grouped --metric C0 "
                        "--x-param mA2 --group-by B0`",
                    ],
                )
                return 1
            fig, ax = plt.subplots(
                1,
                1,
                figsize=figsize or (9, 6),
                constrained_layout=True,
            )
            render_sweep_1d_grouped(
                ax,
                results,
                raw_metric or "",
                x_param=x_param,
                group_param=group_by,
                overlay=getattr(args, "overlay", None),
                log_y=getattr(args, "log_y", False),
                log_x=getattr(args, "log_x", False),
            )

        elif plot_type == "sweep-histogram":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for sweep-histogram plots",
                    ["Example: `--metric log10_A`"],
                )
                return 1
            fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
            render_sweep_histogram(ax, results, raw_metric)

        elif plot_type == "sweep-divergence":
            fig, ax = plt.subplots(1, 1, figsize=figsize or (10, 6))
            render_sweep_divergence(ax, results)

        elif plot_type == "replicate-convergence":
            if raw_metric is None:
                error_with_hint(
                    "--metric is required for replicate-convergence plots",
                    ["Example: `--metric P_max`"],
                )
                return 1
            if not results.has_replicates:
                error_with_hint(
                    "replicate-convergence requires ensemble data "
                    "(use --n-replicates in sweep)",
                    ["Re-run sweep with `--n-replicates 10`"],
                )
                return 1
            fig, ax = plt.subplots(1, 1, figsize=figsize or (8, 5))
            render_replicate_convergence(ax, results, raw_metric)

        else:
            _cerror(f"unknown sweep plot type '{plot_type}'")
            return 1

    except ValueError as exc:
        _cerror(str(exc))
        return 1

    # Apply custom title
    if args.title:
        fig.suptitle(args.title)

    # Output path
    output_path = Path(args.output) if args.output else data_path / f"{plot_type}.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # sweep-scatter and sweep-grouped manage their own layout via
    # constrained_layout; calling tight_layout on them produces warnings
    # and overlapping colorbars. Skip those cases.
    manual_layout = {"sweep-scatter", "sweep-grouped"}
    if plot_type not in manual_layout:
        plt.tight_layout()
    fig.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    if not args.quiet:
        print(f"Saved: {output_path.resolve()}")

    return 0
