"""``tidal plot campaign`` — multi-directory sweep composition.

Reads multiple sweep subdirectories and composes a multi-panel figure
where each subplot is a standard sweep panel (1D or 2D).

Usage::

    tidal plot campaign /tmp/sweeps/phase1 --metric log10_A --layout 3x3

This replaces standalone scripts like ``analyze_phase1.py`` and much of
``visualize_sweep_campaign.py`` by composing CLI-standard panels.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace


def render_campaign(args: Namespace) -> int:  # noqa: C901, PLR0912, PLR0915
    """Render a multi-panel figure from multiple sweep subdirectories.

    Returns 0 on success, 1 on error.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    from tidal.cli._console import error as _cerror
    from tidal.cli._console import error_with_hint
    from tidal.cli._sweep_panels import (
        _greek_label,  # pyright: ignore[reportPrivateUsage]
        render_sweep_1d,
        render_sweep_2d,
    )
    from tidal.measurement._sweep_results import SweepResults

    campaign_dir = Path(args.data_path)
    if not campaign_dir.is_dir():
        _cerror(f"'{campaign_dir}' is not a directory")
        return 1

    # Discover sweep subdirectories (each must have sweep.json)
    subdirs = sorted(
        d for d in campaign_dir.iterdir() if d.is_dir() and (d / "sweep.json").exists()
    )

    if not subdirs:
        error_with_hint(
            f"No sweep subdirectories found in '{campaign_dir}'",
            ["Each subdirectory must contain sweep.json + results.json"],
        )
        return 1

    # Parse layout
    layout_raw: str | None = getattr(args, "layout", None)
    n_panels = len(subdirs)
    if layout_raw:
        parts = layout_raw.lower().split("x")
        nrows, ncols = int(parts[0]), int(parts[1])
    else:
        # Auto layout
        import math

        ncols = min(3, n_panels)
        nrows = math.ceil(n_panels / ncols)

    raw_metric: str | None = getattr(args, "metric", None)
    if raw_metric is None:
        error_with_hint(
            "--metric is required for campaign plots",
            ["Example: `--metric log10_A` or `--metric P_max`"],
        )
        return 1

    figsize_raw: str | None = getattr(args, "figsize", None)
    if figsize_raw:
        w, h = (float(x) for x in figsize_raw.split(","))
        figsize = (w, h)
    else:
        figsize = (5 * ncols, 4 * nrows)

    dpi = getattr(args, "dpi", None) or 150
    use_scatter: bool = getattr(args, "scatter", False)
    ann_extremes: bool = getattr(args, "annotate_extremes", False)
    use_arctan: bool = getattr(args, "arctan_axes", False)
    log_y: bool = getattr(args, "log_y", False)

    # Classify panels if requested
    classify_raw: str | None = getattr(args, "classify", None)
    classify_threshold: float | None = None
    if classify_raw:
        # Format: "active:A>1.01" — extract the threshold
        parts = classify_raw.split(":")
        if len(parts) >= 2:  # noqa: PLR2004
            import re

            match = re.search(r"[><=]+([\d.]+)", parts[1])
            if match:
                classify_threshold = float(match.group(1))

    # Derived metrics auto-derive
    derived_metrics = {"A", "log10_A", "C0", "P_EM", "P_valid"}
    need_derive = raw_metric in derived_metrics
    baseline = getattr(args, "baseline_formula", None)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    fig.suptitle(
        getattr(args, "title", None) or f"Campaign: {raw_metric}",
        fontsize=14,
        fontweight="bold",
    )

    for idx, subdir in enumerate(subdirs):
        row_idx = idx // ncols
        col_idx = idx % ncols
        if row_idx >= nrows:
            break
        ax = axes[row_idx][col_idx]

        try:
            results = SweepResults.from_directory(subdir)
        except (FileNotFoundError, ValueError) as exc:
            ax.text(0.5, 0.5, f"Error: {exc}", transform=ax.transAxes, ha="center")
            ax.set_title(subdir.name)
            continue

        if need_derive:
            results.add_derived_columns(
                baseline_formula=baseline or "sin(kappa * B0 * t_end / 2)**2",
            )

        n_swept = len(results.swept_params)
        param_name = next(iter(results.swept_params.keys())) if n_swept >= 1 else "?"

        # Classification coloring
        if classify_threshold is not None and raw_metric in results.metric_names:
            metric_vals = results.column(raw_metric)
            valid_vals = metric_vals[np.isfinite(metric_vals)]
            if len(valid_vals) > 0:
                max_val = float(np.max(np.abs(valid_vals)))
                is_active = max_val > classify_threshold
                if is_active:
                    ax.set_facecolor("#FFF8E1")
                else:
                    ax.set_facecolor("#F5F5F5")

        if n_swept == 1:
            thresholds: list[str] = getattr(args, "hline", []) or []
            render_sweep_1d(
                ax,
                results,
                raw_metric,
                log_y=log_y,
                scatter=use_scatter,
                annotate_extremes=ann_extremes,
                arctan_axes=use_arctan,
                thresholds=thresholds,
            )
        elif n_swept == 2:  # noqa: PLR2004
            dc: float | None = getattr(args, "divergent_center", None)
            cmap_name: str = getattr(args, "cmap", None) or (
                "RdBu_r" if dc is not None else "viridis"
            )
            clamp_raw: str | None = getattr(args, "clamp_color", None)
            clamp: tuple[float, float] | None = None
            if clamp_raw:
                parts = clamp_raw.split(":")
                clamp = (float(parts[0]), float(parts[1]))
            render_sweep_2d(
                ax,
                results,
                raw_metric,
                divergent_center=dc,
                cmap_name=cmap_name,
                clamp_color=clamp,
            )
        else:
            ax.text(
                0.5,
                0.5,
                f"{n_swept}D sweep\n(unsupported in campaign)",
                transform=ax.transAxes,
                ha="center",
            )

        # Title: show all swept param names
        if n_swept == 1:
            title_label = _greek_label(param_name)
        elif n_swept >= 2:  # noqa: PLR2004
            pnames = list(results.swept_params.keys())
            title_label = " x ".join(_greek_label(p) for p in pnames)
        else:
            title_label = subdir.name
        ax.set_title(
            title_label,
            fontsize=11,
            fontweight="bold"
            if ax.get_facecolor()[:3] != (0.96, 0.96, 0.96)
            else "normal",
        )

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        row_idx = idx // ncols
        col_idx = idx % ncols
        if row_idx < nrows:
            axes[row_idx][col_idx].set_visible(False)

    plt.tight_layout(rect=(0, 0, 1, 0.95))

    # Output
    output: str | None = getattr(args, "output", None)
    out_path = Path(output) if output else campaign_dir / "campaign.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return 0
