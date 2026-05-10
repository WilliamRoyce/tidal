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


def _grid_shape(faces: int) -> tuple[int, int]:
    """Pick (rows, cols) close to square with rows * cols >= faces.

    Hard-codes the supervisor's example layout for the most common cases.
    """
    layouts: dict[int, tuple[int, int]] = {
        4: (2, 2),
        6: (2, 3),
        8: (2, 4),
        10: (2, 5),
        12: (3, 4),
        14: (2, 7),
        16: (4, 4),
    }
    if faces in layouts:
        return layouts[faces]
    rows = max(1, math.floor(math.sqrt(faces)))
    cols = math.ceil(faces / rows)
    return rows, cols


# ----------------------------------------------------------------------
# Per-face render
# ----------------------------------------------------------------------


def _render_face_panel(
    fig: object,  # matplotlib Figure | SubFigure (untyped to keep deps light)
    face_idx: int,
    chi: NDArray[np.float64] | None,
    n_dims: int,
) -> None:
    """Render one face panel: lower-triangle 2D histograms.

    Diagonal hidden, upper triangle blank, face label top-right.
    """
    n_chi = n_dims - 1
    axes = fig.subplots(n_chi, n_chi, sharex=False, sharey=False)
    # Collapse axes to 2-D array for uniform indexing even when n_chi == 1.
    axes = np.atleast_2d(axes)

    fm = face_label_math(face_idx)
    fig.text(
        0.99,
        0.99,
        rf"$\mathrm{{Face}}\ {fm}$",
        ha="right",
        va="top",
        fontsize=10,
    )

    for i in range(n_chi):
        for j in range(n_chi):
            ax = axes[i, j]
            if j > i or i == j:
                ax.set_visible(False)
                continue
            if chi is None or len(chi) < 5:
                ax.set_xlim(-1, 1)
                ax.set_ylim(-1, 1)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            # 2D histogram on the [-1, 1]^2 face cube
            hist, _xe, _ye = np.histogram2d(
                chi[:, j],
                chi[:, i],
                bins=32,
                range=((-1.0, 1.0), (-1.0, 1.0)),
            )
            ax.imshow(
                hist.T,
                origin="lower",
                extent=(-1.0, 1.0, -1.0, 1.0),
                aspect="auto",
                cmap="Blues",
                interpolation="bilinear",
            )
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
            ax.set_xticks([])
            ax.set_yticks([])

    # LaTeX axis labels on the leftmost column and bottom row.
    for i in range(n_chi):
        ax = axes[i, 0]
        if not ax.get_visible():
            continue
        ax.set_ylabel(rf"$\chi_{{{i + 2}}}^{{{fm}}}$", fontsize=8)
    for j in range(n_chi):
        ax = axes[n_chi - 1, j]
        if not ax.get_visible():
            continue
        ax.set_xlabel(rf"$\chi_{{{j + 1}}}^{{{fm}}}$", fontsize=8)


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def plot_atlas(
    survey_dir: Path,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> Path:
    """Render the cubed-sphere atlas for a survey directory.

    The survey directory holds per-tile subdirectories named
    ``<face_label>_tile<sub_tile>/``, each written by ``tidal sample
    --joint-prior``.  This routine pools tiles by face, recovers the
    face-local cube coordinates ``chi`` for each sample, and renders a
    lower-triangle 2D-histogram panel per face.  Diagonal panels are
    hidden; the upper triangle is blank; the face label is drawn in the
    top-right corner of each panel.  Panels are arranged in an
    auto-sized grid (3 rows x 4 cols for N=6).

    Parameters
    ----------
    survey_dir : Path
        Directory containing one ``<face_label>_tile<sub>/`` subdir per
        cell (output of ``tidal sample --joint-prior``).
    output_path : Path or None
        Where to write the PDF.  Defaults to ``survey_dir / atlas.pdf``.
    show : bool
        If True, display the figure interactively (blocks).

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

    # Pool by face.  Determine N from the first cell's metadata.
    by_face: dict[int, list[NDArray[np.float64]]] = {}
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

    if n_dims is None:
        msg = (
            f"no usable cell metadata found in {survey_dir}; "
            f"each cell must have inference.json + results.csv"
        )
        raise ValueError(msg)

    n_faces_total = n_faces(n_dims)
    rows, cols = _grid_shape(n_faces_total)

    # Layout figure.
    panel_in = 1.6  # inches per face panel; gives ~5x6 inch atlas for N=6
    fig = plt.figure(
        figsize=(cols * panel_in, rows * panel_in), constrained_layout=False
    )
    subfigs = fig.subfigures(rows, cols, wspace=0.05, hspace=0.05)
    subfigs = np.atleast_2d(subfigs)

    for face_idx in range(1, n_faces_total + 1):
        slot = face_idx - 1  # 0-indexed slot in the grid
        r = slot // cols
        c_idx = slot % cols
        sub = subfigs[r, c_idx]
        chunks = by_face.get(face_idx)
        chi: NDArray[np.float64] | None
        chi = np.concatenate(chunks, axis=0) if chunks else None
        _render_face_panel(sub, face_idx, chi, n_dims)

    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return output_path


__all__ = ["cell_label", "plot_atlas"]
