"""Cubed-sphere atlas plot for posteriors over a coupling sphere.

Produces a 2N-face panel grid showing the per-face posterior over the
face-local cube coordinates ``chi_i^{k+/-}``.  See the module
:mod:`tidal.inference._sphere` for the underlying chart and
:class:`tidal.inference._prior.RadialAngularPrior` for the prior.

A survey directory contains per-tile output subdirectories named
``<face_label>_tile<sub_tile>/``, each written by ``tidal sample
--joint-prior`` for one cell of the cubed-sphere.  The atlas pools
tiles belonging to the same face and renders one filled-contour KDE
per face, arranged in a grid (3x4 for N=6 -> 12 faces).

Public API
----------
``plot_atlas(survey_dir, output_path=None, ...) -> Path``
"""

from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from tidal.inference._sphere import (
    cell_label,
    face_label_math,
    face_to_axis_sign,
    n_faces,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

logger = logging.getLogger(__name__)


_CELL_DIR_RE = re.compile(r"^(\d{2})([pm])_tile(\d+(?:_\d+)*)$")


def _parse_cell_dir_name(name: str) -> tuple[int, tuple[int, ...]] | None:
    """Parse ``"01p_tile2_3"`` into ``(face_idx, sub_tile)``.

    Returns ``None`` if the name does not match the cell-dir convention.
    """
    m = _CELL_DIR_RE.match(name)
    if m is None:
        return None
    axis_one_indexed = int(m.group(1))
    sign_char = m.group(2)
    sub_tile = tuple(int(s) for s in m.group(3).split("_"))
    face_idx = 2 * (axis_one_indexed - 1) + (1 if sign_char == "p" else 2)
    return face_idx, sub_tile


def _load_cell_metadata(cell_dir: Path) -> dict[str, object] | None:
    """Load ``inference.json`` from one cell's output dir.

    Returns ``None`` if the file does not exist or is unreadable.
    """
    meta_path = cell_dir / "inference.json"
    if not meta_path.exists():
        return None
    try:
        with meta_path.open() as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def _find_joint_prior_record(
    meta: dict[str, object],
) -> dict[str, object] | None:
    """Return the radial-angular prior record from a saved metadata dict."""
    priors = meta.get("priors", [])
    if not isinstance(priors, list):
        return None
    for p in priors:
        if isinstance(p, dict) and p.get("kind") == "radial_angular":
            return p
    return None


def _load_cell_samples(
    cell_dir: Path,
    param_names: Sequence[str],
) -> NDArray[np.float64] | None:
    """Load ``results.csv`` columns for the given ``param_names``.

    Returns an ``(n_samples, n_params)`` matrix, or ``None`` if the file
    is missing / empty.
    """
    csv_path = cell_dir / "results.csv"
    if not csv_path.exists():
        return None
    import csv

    rows: list[dict[str, str]] = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        return None
    samples = np.zeros((len(rows), len(param_names)), dtype=np.float64)
    for i, row in enumerate(rows):
        for j, name in enumerate(param_names):
            samples[i, j] = float(row[name])
    return samples


def _load_cell_weights(cell_dir: Path) -> NDArray[np.float64] | None:
    """Load the posterior ``weight`` column from one cell's ``results.csv``.

    Nested-sampling chains have a ``weight`` column for posterior-weighted
    plotting; MC chains do not.  Returns ``None`` when the column is
    absent so the caller can fall back to uniform weights.
    """
    csv_path = cell_dir / "results.csv"
    if not csv_path.exists():
        return None
    import csv

    with csv_path.open() as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or "weight" not in reader.fieldnames:
            return None
        weights = [float(row["weight"]) for row in reader]
    if not weights:
        return None
    return np.asarray(weights, dtype=np.float64)


def _physical_to_face_local(
    c: NDArray[np.float64],
    face_idx: int,
    Q: NDArray[np.float64],  # noqa: N803
) -> NDArray[np.float64]:
    """Invert the cubed-sphere projection.

    Given an ``(n_samples, N)`` matrix of physical coupling vectors
    (samples drawn from ``RadialAngularPrior``), recover the face-local
    cube coordinates ``chi`` of shape ``(n_samples, N - 1)`` for the
    given ``face_idx``.

    Derivation:  ``c = r * Q^T (axis_vec + embed @ chi) / |.|``.
    Set ``v_unit = Q c / |c|``; then ``v_unit = v / |v|`` for the
    cube-vector ``v`` whose dominant component is ``v_{k-1} = s``.
    Therefore ``chi_j = v_unit_j / |v_unit_{k-1}|`` for the (N - 1)
    non-dominant slots.
    """
    q_arr = np.asarray(Q)
    n = q_arr.shape[0]
    if c.shape[1] != n:
        msg = f"c has {c.shape[1]} cols, expected {n}"
        raise ValueError(msg)

    # v_unit[i] = (Q @ c[i]) / |c[i]|; vectorize across samples.
    norms = np.linalg.norm(c, axis=1, keepdims=True)
    safe = np.where(norms > 0, norms, 1.0)
    v_unit = c @ q_arr.T / safe  # shape (n_samples, N)

    k, _s = face_to_axis_sign(face_idx)
    dom = np.abs(v_unit[:, k - 1])
    safe_dom = np.where(dom > 0, dom, 1.0)
    # Drop the dominant axis column to get chi of shape (n_samples, N - 1).
    others = np.delete(v_unit, k - 1, axis=1)
    return others / safe_dom[:, None]


# ----------------------------------------------------------------------
# Layout
# ----------------------------------------------------------------------


def _grid_shape(n_dims: int, cols: int | None = None) -> tuple[int, int]:
    """``(rows, cols)`` for the atlas layout.

    Default (``cols=None``) gives the ``2 x n_dims`` layout: top row
    positive faces, bottom row negative faces, column ``k - 1`` = axis
    ``k``.  Used for survey theories where ``2 * n_dims`` panels fit
    naturally side-by-side (N=2..6).

    With ``cols=K`` the layout becomes ``2 * ceil(n_dims / K) x K``: the
    axes are split into ``ceil(n_dims / K)`` blocks of ``K`` consecutive
    axes, each block occupying two rows (up over down).  Used for high-N
    theories where ``2 x n_dims`` is too wide for a single page; e.g.
    N=8 with ``cols=4`` becomes a 4x4 square layout.

    Raises ``ValueError`` if ``n_dims < 2`` or ``cols < 1``.
    """
    if n_dims < 2:
        msg = f"n_dims must be >= 2; got {n_dims}"
        raise ValueError(msg)
    if cols is None:
        return 2, n_dims
    if cols < 1:
        msg = f"cols must be >= 1; got {cols}"
        raise ValueError(msg)
    blocks = math.ceil(n_dims / cols)
    return 2 * blocks, cols


def _face_to_slot(face_idx: int, cols: int | None = None) -> tuple[int, int]:
    """``face_idx`` (1-indexed) -> ``(row, col)`` in the atlas layout.

    Default (``cols=None``) is the 2 x N rule: face ``2k - 1`` is the
    ``+k`` ("up") face -> ``(0, k - 1)``; face ``2k`` is the ``-k``
    ("down") face -> ``(1, k - 1)``.

    With ``cols=K`` (block-pair layout): axis ``k`` sits in block
    ``(k - 1) // K`` at column ``(k - 1) % K``; up-faces go to the top
    row of the block, down-faces to the bottom row.  Final row =
    ``2 * block + (0 if up else 1)``.
    """
    axis_one_indexed, sign = face_to_axis_sign(face_idx)
    if cols is None:
        row = 0 if sign > 0 else 1
        col = axis_one_indexed - 1
        return row, col
    if cols < 1:
        msg = f"cols must be >= 1; got {cols}"
        raise ValueError(msg)
    block = (axis_one_indexed - 1) // cols
    col = (axis_one_indexed - 1) % cols
    row = 2 * block + (0 if sign > 0 else 1)
    return row, col


# ----------------------------------------------------------------------
# Per-face render
# ----------------------------------------------------------------------


def _render_face_panel(
    fig: object,  # matplotlib Figure | SubFigure (untyped to keep deps light)
    face_idx: int,
    chi: NDArray[np.float64] | None,
    n_dims: int,
    weights: NDArray[np.float64] | None = None,
    *,
    show_xlabels: bool = True,
    show_ylabels: bool = True,
    fill_colors: tuple[str, str] | None = None,
    color: str | None = None,
    color_alpha: float = 0.5,
) -> None:
    r"""Render one face panel via anesthetic — same visual style as the standalone corner plot.

    Each panel is an ``(N-1) x (N-1)`` lower-triangle 2D KDE corner plot
    of the face-local cube coordinates ``chi_1, ..., chi_{N-1}``, with
    a 1D KDE marginal on the diagonal — two-tone fills (default
    Planck-blue ``#aac8e9`` / ``#3877b8``; override via ``fill_colors``).
    Numeric tick marks are stripped (``ticks=None``).  The
    ``\\chi_1, \\chi_2, ..., \\chi_{N-1}`` LaTeX axis labels are kept
    on the leftmost column / bottom row of each panel BUT only when
    ``show_xlabels`` / ``show_ylabels`` are True — the atlas driver
    sets these so labels appear only on the outer edges of the whole
    atlas grid (bottom row of the outer grid for x-labels, leftmost
    column for y-labels).  A ``Face k↑`` / ``Face k↓`` text annotation
    is drawn in the panel's upper-right whitespace.

    Composability is via :func:`anesthetic.make_2d_axes`'s
    ``subplot_spec`` argument, exposed through
    :func:`tidal.inference._visualize._render_anesthetic_corner_into`.
    """
    from tidal.inference._visualize import _render_anesthetic_corner_into

    n_chi = n_dims - 1
    fm = face_label_math(face_idx)

    # Empty / sparse cells: draw an empty placeholder + face label and bail.
    # anesthetic's KDE requires a non-trivial sample (>= ~5 finite points
    # spanning at least 2D), and a face with no live tile would otherwise
    # explode.
    if chi is None or len(chi) < 5:
        ax = fig.subplots(1, 1)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.5,
            0.5,
            "(no samples)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=8,
            color="0.5",
        )
        fig.text(
            0.72,
            0.78,
            rf"$\mathrm{{Face}}\ {fm}$",
            ha="center",
            va="center",
            fontsize=9,
        )
        return

    # Build an anesthetic Samples object over the face-local chi columns.
    # Column names are kept legal-identifier (chi_1, ..., chi_{N-1}); the
    # plain-chi LaTeX labels are passed via the ``labels`` dict.  We
    # deliberately drop the face superscript: at 4x4 atlas size
    # \chi_i^{k\uparrow/\downarrow} is unreadable and the face identity
    # is already carried by the panel-title annotation below.  The
    # caption explains that chi_i means a different physical coupling
    # per face.
    from anesthetic import MCMCSamples

    col_names = [f"theta_{i + 1}" for i in range(n_chi)]
    latex_labels = {col_names[i]: rf"$\theta_{{{i + 1}}}$" for i in range(n_chi)}
    data = {col_names[i]: chi[:, i] for i in range(n_chi)}
    if weights is not None and len(weights) == len(chi):
        # Drop columns to constructor + set weights post-hoc; anesthetic's
        # MCMCSamples accepts a ``weights`` kwarg directly.
        samples = MCMCSamples(data=data, weights=weights)
    else:
        samples = MCMCSamples(data=data)

    # Render into the SubFigure passed in by the atlas driver.
    # ``ticks=None`` strips numeric tick marks; the axis labels live in
    # the ``labels`` mapping.
    _, axes_df = _render_anesthetic_corner_into(
        samples,
        col_names,
        target_fig=fig,
        labels=latex_labels,
        ticks=None,
        show_diagonal=True,
        fill_two_tone=True,
        fill_colors=fill_colors,
        color=color,
        color_alpha=color_alpha,
    )

    # Tighten the SubFigure's internal layout so the corner-plot grid
    # fills nearly the whole panel: anesthetic's default GridSpec leaves
    # generous left/right/top/bottom margins which, with 16 panels in
    # the YM--PGT atlas, compounds into a lot of dead whitespace.  Leave
    # enough bottom/left margin for outer-edge axis labels (handled by
    # the show_xlabels / show_ylabels flags above) but otherwise push
    # the grid flush to the SubFigure boundary.
    fig.subplots_adjust(
        left=0.10,
        right=0.98,
        bottom=0.10,
        top=0.98,
        wspace=0.0,
        hspace=0.0,
    )

    # Suppress non-edge labels.  ``set_visible(False)`` skips layout
    # allocation entirely; setting the text to "" would leave a blank
    # padding band of approximately the label height.
    if not show_xlabels or not show_ylabels:
        for ax in axes_df.values.flatten():  # type: ignore[attr-defined]
            if ax is None:
                continue
            if not show_xlabels:
                ax.xaxis.label.set_visible(False)
            if not show_ylabels:
                ax.yaxis.label.set_visible(False)

    fig.text(
        0.72,
        0.78,
        rf"$\mathrm{{Face}}\ {fm}$",
        ha="center",
        va="center",
        fontsize=9,
    )


def _render_axis_panel(
    fig: object,  # matplotlib Figure | SubFigure
    axis_idx: int,
    chi_up: NDArray[np.float64] | None,
    chi_down: NDArray[np.float64] | None,
    n_dims: int,
    *,
    weights_up: NDArray[np.float64] | None = None,
    weights_down: NDArray[np.float64] | None = None,
    show_xlabels: bool = True,
    show_ylabels: bool = True,
    color_up: str,
    color_down: str,
    color_alpha: float = 0.5,
) -> None:
    r"""Render one axis panel via anesthetic — overlay of up- and down-face posteriors.

    Parallel to :func:`_render_face_panel`, but pools the positive-sign
    (``+x_k`` = face ``2k - 1``) and negative-sign (``-x_k`` = face
    ``2k``) face-local chi data into a single panel.  Up posterior is
    drawn first in ``color_up`` (typically IBM magenta ``#dc267f``);
    down is overlaid second in ``color_down`` (typically IBM yellow
    ``#ffb000``) — both at ``color_alpha = 0.5`` so overlapping support
    reads as a darker mix.

    Panel title is ``Axis k ↑↓`` with the ↑ glyph in ``color_up`` and
    the ↓ glyph in ``color_down`` so the title doubles as a per-panel
    legend.
    """
    n_chi = n_dims - 1

    # Empty / sparse panel: at least one of (chi_up, chi_down) must carry
    # enough samples for anesthetic's KDE.  If both are empty, draw a
    # placeholder and bail.
    have_up = chi_up is not None and len(chi_up) >= 5
    have_down = chi_down is not None and len(chi_down) >= 5
    if not have_up and not have_down:
        ax = fig.subplots(1, 1)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.5,
            0.5,
            "(no samples)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=8,
            color="0.5",
        )
        _draw_axis_title(fig, axis_idx, color_up, color_down)
        return

    from anesthetic import MCMCSamples, make_2d_axes

    col_names = [f"theta_{i + 1}" for i in range(n_chi)]
    latex_labels = {col_names[i]: rf"$\theta_{{{i + 1}}}$" for i in range(n_chi)}

    # Build the axes grid once for both overlays.
    _, axes_df = make_2d_axes(
        col_names,
        labels=latex_labels,
        ticks=None,
        lower=True,
        diagonal=True,
        upper=False,
        fig=fig,
    )

    # Overlay each posterior on the same axes_df.  anesthetic's plot_2d
    # accepts an explicit AxesDataFrame to render into; calling it twice
    # paints both posteriors on the same panel.  Order matters
    # cosmetically only: down drawn second sits on top of up where they
    # overlap, which mixes the two alphas into a darker tone.
    for chi, weights, color, present in (
        (chi_up, weights_up, color_up, have_up),
        (chi_down, weights_down, color_down, have_down),
    ):
        if not present or chi is None:
            continue
        data = {col_names[i]: chi[:, i] for i in range(n_chi)}
        if weights is not None and len(weights) == len(chi):
            samples = MCMCSamples(data=data, weights=weights)
        else:
            samples = MCMCSamples(data=data)
        samples.plot_2d(axes_df, color=color, alpha=color_alpha)

    # ``plot_2d`` overwrites the LaTeX labels with the column names;
    # re-apply via the public ``set_labels`` method.
    axes_df.set_labels(latex_labels)  # type: ignore[attr-defined]

    # Anesthetic auto-scales each panel to its sample range, which
    # gives the visual impression of a hard cutoff when the posterior
    # falls short of the cube edge.  The cube coordinates are by
    # construction in [-1, 1], so force every panel axis to show the
    # full extent.
    for ax in axes_df.values.flatten():  # type: ignore[attr-defined]
        if ax is None:
            continue
        ax.set_xlim(-1.0, 1.0)
        if hasattr(ax, "position") and ax.position == "diagonal":
            # Diagonal cells carry 1D KDEs; only x is bounded by [-1, 1],
            # the y-axis is density-scaled so leave it alone.
            continue
        ax.set_ylim(-1.0, 1.0)

    # Tighten the SubFigure's internal layout (same constants as the
    # single-face renderer; the 2x4 outer grid gives each panel ~3.5x
    # the area of the previous 4x4 case so margins are not pinched).
    fig.subplots_adjust(
        left=0.10,
        right=0.98,
        bottom=0.10,
        top=0.98,
        wspace=0.0,
        hspace=0.0,
    )

    # Suppress non-edge labels.
    if not show_xlabels or not show_ylabels:
        for ax in axes_df.values.flatten():  # type: ignore[attr-defined]
            if ax is None:
                continue
            if not show_xlabels:
                ax.xaxis.label.set_visible(False)
            if not show_ylabels:
                ax.yaxis.label.set_visible(False)

    _draw_axis_title(fig, axis_idx, color_up, color_down)


def _draw_axis_title(
    fig: object,
    axis_idx: int,
    color_up: str,
    color_down: str,
) -> None:
    r"""Draw ``Axis k ↑↓`` in the panel's upper-right whitespace.

    The ↑ glyph renders in ``color_up`` and the ↓ glyph in ``color_down``,
    so the title doubles as a per-panel legend for the two-posterior
    overlay below.

    Implemented as three matplotlib ``fig.text`` calls (axis number in
    black, ↑ in color_up, ↓ in color_down) rather than one LaTeX
    ``\textcolor`` call because matplotlib's mathtext renderer does
    not honor ``\textcolor`` without a full ``text.usetex=True``
    config; per-text ``color=`` kwargs work in every backend.
    """
    # Layout: three text elements packed left-to-right inside the empty
    # upper-right triangle of the lower-triangle corner plot.  All three
    # use ``ha="left"`` from successive anchor x-positions so the glyphs
    # do not overlap (matplotlib's mathtext ``ha="center"`` was tight
    # enough at small panel size that ↑↓ collided with the trailing axis
    # digit).  Single consistent fontsize keeps the arrows in line with
    # the axis number.
    fig.text(
        0.60,
        0.78,
        rf"$\mathrm{{Axis}}\ {axis_idx}$",
        ha="left",
        va="center",
        fontsize=9,
    )
    fig.text(
        0.85,
        0.78,
        r"$\uparrow$",
        ha="left",
        va="center",
        fontsize=10,
        color=color_up,
    )
    fig.text(
        0.92,
        0.78,
        r"$\downarrow$",
        ha="left",
        va="center",
        fontsize=10,
        color=color_down,
    )


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def plot_atlas(
    survey_dir: Path,
    output_path: Path | None = None,
    *,
    show: bool = False,
    layout_cols: int | None = None,
    fill_colors: tuple[str, str] | None = None,
    color: str | None = None,
    color_alpha: float = 0.5,
    overlay_up_down: bool = False,
    color_up: str | None = None,
    color_down: str | None = None,
) -> Path:
    """Render the cubed-sphere atlas for a survey directory.

    The survey directory holds per-tile subdirectories named
    ``<face_label>_tile<sub_tile>/``, each written by ``tidal sample
    --joint-prior``.  This routine pools tiles by face, recovers the
    face-local cube coordinates ``chi`` for each sample, and renders a
    lower-triangle 2D-histogram panel per face.  Diagonal panels are
    hidden; the upper triangle is blank; the face label is drawn in the
    top-right corner of each panel.

    Layout: default ``layout_cols=None`` produces the ``2 x N`` grid
    (top row +k faces, bottom row -k faces, column = axis).  Pass
    ``layout_cols=K`` to use the ``2 * ceil(N / K) x K`` block-pair
    layout — used for high-N theories where the default 2 x N strip is
    too wide for a single page (e.g. N=8 with ``layout_cols=4`` gives a
    4 x 4 square).

    Parameters
    ----------
    survey_dir : Path
        Directory containing one ``<face_label>_tile<sub>/`` subdir per
        cell (output of ``tidal sample --joint-prior``).
    output_path : Path or None
        Where to write the PDF.  Defaults to ``survey_dir / atlas.pdf``.
    show : bool
        If True, display the figure interactively (blocks).
    layout_cols : int or None
        Column count for the block-pair layout (see above).  ``None``
        (default) preserves the ``2 x n_dims`` layout.
    fill_colors : (str, str) or None
        ``(outer_95, inner_68)`` hex pair overriding the default
        two-tone fills.  Mutually exclusive with ``color`` (the
        single-tone path wins; ``fill_colors`` is ignored when
        ``color`` is set).  ``None`` (default) preserves the existing
        Planck-blue palette.
    color : str or None
        Single hex color applied uniformly to every artist (1D KDE
        curves, 2D fills, contour outlines) via anesthetic's
        ``samples.plot_2d(color=...)``.  Matches the manuscript's
        results-section corner-plot style — see
        ``scripts/figures/_corner_style.py``.  ``None`` (default) keeps
        anesthetic's default palette.
    color_alpha : float
        Alpha for single-tone ``color`` and for ``overlay_up_down``
        overlays.  Defaults to 0.5, matching the manuscript's
        ``OVERLAY_ALPHA`` from ``_corner_style.py``.
    overlay_up_down : bool
        When True, render N panels (one per coupling axis) with the
        positive-axis face overlay in ``color_up`` and the
        negative-axis face overlay in ``color_down``.  Halves the panel
        count vs. the default single-face-per-panel mode and is the
        rendering used by the bespoke K.2 manuscript figure.  When
        False (default), retains the per-face rendering and
        ``color_up`` / ``color_down`` are ignored.  ``layout_cols`` is
        ignored in overlay mode (layout fixed at ``2 x ceil(N/2)``).
        ``color`` and ``fill_colors`` are also ignored in overlay mode
        (the up/down pair carries the coloring).
    color_up : str or None
        Hex color for the positive-axis (``+x_k``, "up") face overlay.
        Used only when ``overlay_up_down=True``.  Recommend IBM
        colorblind magenta ``#dc267f`` (``_palette.AMP_COLOR``).
    color_down : str or None
        Hex color for the negative-axis (``-x_k``, "down") face
        overlay.  Used only when ``overlay_up_down=True``.  Recommend
        IBM colorblind yellow ``#ffb000`` (``_palette.SUP_COLOR``).

    Returns
    -------
    Path
        The output PDF path.
    """
    import matplotlib as mpl

    if not show:
        mpl.use("Agg")
    import matplotlib.pyplot as plt

    survey_dir = Path(survey_dir)
    if output_path is None:
        output_path = survey_dir / "atlas.pdf"
    output_path = Path(output_path)

    # Discover cells.
    cells: list[tuple[int, tuple[int, ...], Path]] = []
    for child in sorted(survey_dir.iterdir()):
        if not child.is_dir():
            continue
        parsed = _parse_cell_dir_name(child.name)
        if parsed is None:
            continue
        face_idx, sub_tile = parsed
        cells.append((face_idx, sub_tile, child))

    if not cells:
        msg = (
            f"no <face_label>_tile<sub>/ subdirectories found in "
            f"{survey_dir}; nothing to plot"
        )
        raise ValueError(msg)

    # Pool by face.  Determine N from the first cell's metadata.  Parallel
    # ``by_face_weights`` carries the per-sample posterior weight column
    # (when the chain has one) so the panel KDE is posterior-weighted —
    # an unweighted KDE on a nested-sampling chain misrepresents the
    # posterior because dead points carry exp(logZ) weights that vary by
    # many orders of magnitude.
    by_face: dict[int, list[NDArray[np.float64]]] = {}
    by_face_weights: dict[int, list[NDArray[np.float64]]] = {}
    n_dims: int | None = None
    for face_idx, _sub_tile, cell_dir in cells:
        meta = _load_cell_metadata(cell_dir)
        if meta is None:
            logger.warning("cell %s: no inference.json, skipping", cell_dir.name)
            continue
        record = _find_joint_prior_record(meta)
        if record is None:
            logger.warning(
                "cell %s: no radial_angular prior in metadata, skipping",
                cell_dir.name,
            )
            continue
        names = record["names"]
        Q = np.asarray(record["Q"], dtype=np.float64)  # noqa: N806
        if Q.shape != (len(names), len(names)):
            logger.warning(
                "cell %s: Q shape %s does not match N=%d",
                cell_dir.name,
                Q.shape,
                len(names),
            )
            continue
        if n_dims is None:
            n_dims = len(names)
        elif len(names) != n_dims:
            logger.warning(
                "cell %s: N=%d differs from atlas N=%d, skipping",
                cell_dir.name,
                len(names),
                n_dims,
            )
            continue

        if face_idx != record.get("face_idx"):
            logger.warning(
                "cell %s: dir face_idx=%d differs from metadata face_idx=%d",
                cell_dir.name,
                face_idx,
                record.get("face_idx"),
            )
            continue

        c = _load_cell_samples(cell_dir, names)
        if c is None or len(c) == 0:
            logger.warning("cell %s: no samples in results.csv", cell_dir.name)
            continue

        chi = _physical_to_face_local(c, face_idx, Q)
        by_face.setdefault(face_idx, []).append(chi)
        w = _load_cell_weights(cell_dir)
        if w is not None and len(w) == len(c):
            by_face_weights.setdefault(face_idx, []).append(w)

    if n_dims is None:
        msg = (
            f"no usable cell metadata found in {survey_dir}; "
            f"each cell must have inference.json + results.csv"
        )
        raise ValueError(msg)

    n_faces_total = n_faces(n_dims)
    if overlay_up_down:
        if color_up is None or color_down is None:
            msg = (
                "overlay_up_down=True requires color_up and color_down "
                "(IBM colorblind magenta/yellow recommended); got "
                f"color_up={color_up!r}, color_down={color_down!r}"
            )
            raise ValueError(msg)
        # 2 x ceil(N/2) grid: top row axes 1..ceil(N/2), bottom row
        # axes ceil(N/2)+1..N.  For N=8 this is 2x4.
        rows = 2
        cols = math.ceil(n_dims / 2)
    else:
        rows, cols = _grid_shape(n_dims, layout_cols)

    # KDE sample-compression in anesthetic draws from NumPy's legacy
    # global RNG, so renders are non-deterministic without a seed.  The
    # per-panel renders funnel through _render_anesthetic_corner_into,
    # which seeds and *restores* the global state around each plot_2d
    # (see _deterministic_render).  A bare np.random.seed() here used to
    # leave the process in a fixed state for whatever ran next; that is
    # no longer needed and was itself a leak.  See issue #388.

    # Layout figure.
    panel_in = 1.6  # inches per face panel; gives ~5x6 inch atlas for N=6
    fig = plt.figure(
        figsize=(cols * panel_in, rows * panel_in), constrained_layout=False
    )
    subfigs = fig.subfigures(rows, cols, wspace=0.0, hspace=0.0)
    subfigs = np.atleast_2d(subfigs)

    if overlay_up_down:
        # One panel per axis k: pool face 2k-1 (up) and face 2k (down).
        for axis_idx in range(1, n_dims + 1):
            # 2 x ceil(N/2) layout: top row holds axes 1..ceil(N/2),
            # bottom row holds the rest.
            r = (axis_idx - 1) // cols
            c_idx = (axis_idx - 1) % cols
            sub = subfigs[r, c_idx]
            face_up = 2 * axis_idx - 1
            face_down = 2 * axis_idx
            chunks_up = by_face.get(face_up)
            chunks_down = by_face.get(face_down)
            chi_up: NDArray[np.float64] | None = (
                np.concatenate(chunks_up, axis=0) if chunks_up else None
            )
            chi_down: NDArray[np.float64] | None = (
                np.concatenate(chunks_down, axis=0) if chunks_down else None
            )
            w_up_chunks = by_face_weights.get(face_up)
            w_down_chunks = by_face_weights.get(face_down)
            w_up: NDArray[np.float64] | None = (
                np.concatenate(w_up_chunks, axis=0) if w_up_chunks else None
            )
            w_down: NDArray[np.float64] | None = (
                np.concatenate(w_down_chunks, axis=0) if w_down_chunks else None
            )
            show_xlabels = r == rows - 1
            show_ylabels = c_idx == 0
            _render_axis_panel(
                sub,
                axis_idx,
                chi_up,
                chi_down,
                n_dims,
                weights_up=w_up,
                weights_down=w_down,
                show_xlabels=show_xlabels,
                show_ylabels=show_ylabels,
                color_up=color_up,
                color_down=color_down,
                color_alpha=color_alpha,
            )
    else:
        for face_idx in range(1, n_faces_total + 1):
            # 2 x N up/down layout: face 2k - 1 -> (row=0, col=k-1) (positive
            # axis); face 2k -> (row=1, col=k-1) (negative axis).  When
            # layout_cols is set, the block-pair rule kicks in — see
            # _face_to_slot for details.
            r, c_idx = _face_to_slot(face_idx, layout_cols)
            sub = subfigs[r, c_idx]
            chunks = by_face.get(face_idx)
            chi: NDArray[np.float64] | None
            chi = np.concatenate(chunks, axis=0) if chunks else None
            weights = by_face_weights.get(face_idx)
            w_arr: NDArray[np.float64] | None = (
                np.concatenate(weights, axis=0) if weights else None
            )
            # Show x-axis labels only on the bottom row of the outer grid
            # and y-axis labels only on the leftmost column.  At 4x4 atlas
            # size, per-panel labels become visual noise; outer-edge-only
            # labeling preserves identification without repetition.
            show_xlabels = r == rows - 1
            show_ylabels = c_idx == 0
            _render_face_panel(
                sub,
                face_idx,
                chi,
                n_dims,
                weights=w_arr,
                show_xlabels=show_xlabels,
                show_ylabels=show_ylabels,
                fill_colors=fill_colors,
                color=color,
                color_alpha=color_alpha,
            )

    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return output_path


__all__ = ["cell_label", "plot_atlas"]
