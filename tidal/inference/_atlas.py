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

    # v_unit[i] = (Q @ c[i]) / |c[i]|; vectorise across samples.
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
) -> None:
    """Render one face panel via anesthetic — same visual style as the standalone corner plot.

    Each panel is an ``(N-1) x (N-1)`` lower-triangle 2D KDE corner plot
    of the face-local cube coordinates ``chi_1, ..., chi_{N-1}``, with
    a 1D KDE marginal on the diagonal — two-tone Planck-blue fills
    (``#aac8e9`` / ``#3877b8``) identical to the standalone per-face
    corner plots.  Numeric tick marks are stripped (``ticks=None``) but
    the LaTeX axis labels ``chi_i^{face↑/↓}`` are kept on the leftmost
    column / bottom row to identify the face-local coordinate.  A
    ``Face k↑`` / ``Face k↓`` text annotation is drawn in the panel's
    upper-right whitespace.

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
            0.99,
            0.99,
            rf"$\mathrm{{Face}}\ {fm}$",
            ha="right",
            va="top",
            fontsize=10,
        )
        return

    # Build an anesthetic Samples object over the face-local chi columns.
    # Column names are kept legal-identifier (chi_1, ..., chi_{N-1}); the
    # face-aware LaTeX labels are passed via the ``labels`` dict so they
    # appear on the axes without polluting the underlying column index.
    from anesthetic import MCMCSamples

    col_names = [f"chi_{i + 1}" for i in range(n_chi)]
    latex_labels = {col_names[i]: rf"$\chi_{{{i + 1}}}^{{{fm}}}$" for i in range(n_chi)}
    data = {col_names[i]: chi[:, i] for i in range(n_chi)}
    if weights is not None and len(weights) == len(chi):
        # Drop columns to constructor + set weights post-hoc; anesthetic's
        # MCMCSamples accepts a ``weights`` kwarg directly.
        samples = MCMCSamples(data=data, weights=weights)
    else:
        samples = MCMCSamples(data=data)

    # Render into the SubFigure passed in by the atlas driver.
    # ``ticks=None`` strips numeric tick marks; the axis labels live in
    # the ``labels`` mapping (face-aware LaTeX from face_label_math).
    _render_anesthetic_corner_into(
        samples,
        col_names,
        target_fig=fig,
        labels=latex_labels,
        ticks=None,
        show_diagonal=True,
        fill_two_tone=True,
    )

    fig.text(
        0.99,
        0.99,
        rf"$\mathrm{{Face}}\ {fm}$",
        ha="right",
        va="top",
        fontsize=10,
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
    rows, cols = _grid_shape(n_dims, layout_cols)

    # Seed numpy's legacy RNG so KDE sample-compression in anesthetic
    # (anesthetic.utils.triangular_sample_compression_2d uses bare
    # `np.random.choice` without a Generator) is deterministic across
    # renders.  Without this, repeated atlas renders of the same survey
    # produce pixel-different output, which is bad for review diffs.
    # The legacy global API is required here because anesthetic does not
    # accept an RNG / Generator handle to thread through; NPY002 cannot
    # be addressed at this layer.
    np.random.seed(0)  # noqa: NPY002

    # Layout figure.
    panel_in = 1.6  # inches per face panel; gives ~5x6 inch atlas for N=6
    fig = plt.figure(
        figsize=(cols * panel_in, rows * panel_in), constrained_layout=False
    )
    subfigs = fig.subfigures(rows, cols, wspace=0.05, hspace=0.05)
    subfigs = np.atleast_2d(subfigs)

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
        _render_face_panel(sub, face_idx, chi, n_dims, weights=w_arr)

    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return output_path


__all__ = ["cell_label", "plot_atlas"]
