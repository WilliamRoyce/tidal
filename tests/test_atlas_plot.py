"""Smoke tests for :func:`tidal.inference._atlas.plot_atlas`."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

import numpy as np
import pytest

from tidal.inference._atlas import (
    _grid_shape,
    _parse_cell_dir_name,
    _physical_to_face_local,
    plot_atlas,
)
from tidal.inference._sphere import (
    cell_label,
    enumerate_cells,
    face_to_direction,
    n_faces,
)

if TYPE_CHECKING:
    from pathlib import Path

# ----------------------------------------------------------------------
# Cell-directory name parsing
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("01p_tile1", (1, (1,))),
        ("01p_tile2_3", (1, (2, 3))),
        ("01m_tile1_1_1", (2, (1, 1, 1))),
        ("06p_tile2_2_2_2_2", (11, (2, 2, 2, 2, 2))),
        ("06m_tile1_1_1_1_1", (12, (1, 1, 1, 1, 1))),
    ],
)
def test_parse_cell_dir_name(name: str, expected: tuple[int, tuple[int, ...]]) -> None:
    assert _parse_cell_dir_name(name) == expected


@pytest.mark.parametrize(
    "name",
    [
        "results.csv",
        "_chains",
        "01x_tile1",
        "1p_tile1",  # not zero-padded
        "01p",  # missing tile suffix
        "01p_tile",  # missing tile index
    ],
)
def test_parse_cell_dir_name_rejects_invalid(name: str) -> None:
    assert _parse_cell_dir_name(name) is None


def test_parse_cell_dir_name_round_trip() -> None:
    """``_parse_cell_dir_name(cell_label(face, sub))`` recovers (face, sub)."""
    for face_idx, sub_tile in enumerate_cells(n_dims=3, M=2):
        name = cell_label(face_idx, sub_tile)
        parsed = _parse_cell_dir_name(name)
        assert parsed == (face_idx, sub_tile)


# ----------------------------------------------------------------------
# Inverse projection: physical c -> face-local chi
# ----------------------------------------------------------------------


def test_physical_to_face_local_inverts_face_to_direction() -> None:
    n = 4
    Q = np.eye(n)
    rng = np.random.default_rng(1234)
    for face_idx in range(1, n_faces(n) + 1):
        for _ in range(20):
            chi_in = rng.uniform(-0.95, 0.95, size=n - 1)
            theta_hat = face_to_direction(face_idx, chi_in, Q)
            # Use a non-trivial magnitude r:
            r = float(rng.uniform(0.1, 10.0))
            c = (r * theta_hat)[None, :]  # shape (1, N)
            chi_out = _physical_to_face_local(c, face_idx, Q)
            np.testing.assert_allclose(chi_out[0], chi_in, atol=1e-12)


def test_physical_to_face_local_handles_random_q() -> None:
    """Non-identity Q must invert correctly too."""
    from tidal.inference._sphere import random_rotation

    n = 5
    Q = random_rotation(n, seed=42)
    rng = np.random.default_rng(7)
    chi_in = rng.uniform(-0.9, 0.9, size=n - 1)
    theta_hat = face_to_direction(3, chi_in, Q)
    c = (2.5 * theta_hat)[None, :]
    chi_out = _physical_to_face_local(c, 3, Q)
    np.testing.assert_allclose(chi_out[0], chi_in, atol=1e-12)


# ----------------------------------------------------------------------
# Layout helper
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("n_dims", "expected"),
    [(2, (2, 2)), (3, (2, 3)), (4, (2, 4)), (5, (2, 5)), (6, (2, 6)), (7, (2, 7))],
)
def test_grid_shape_is_2_by_n(n_dims: int, expected: tuple[int, int]) -> None:
    """2 x N up/down layout — top row positive faces, bottom row negative."""
    assert _grid_shape(n_dims) == expected


def test_grid_shape_rejects_n_dims_below_2() -> None:
    with pytest.raises(ValueError, match="n_dims must be >= 2"):
        _grid_shape(1)


def test_face_to_slot_up_down_columns() -> None:
    """Face 2k - 1 (positive axis k) -> (0, k - 1); face 2k -> (1, k - 1)."""
    from tidal.inference._atlas import _face_to_slot

    for k in range(1, 7):
        assert _face_to_slot(2 * k - 1) == (0, k - 1)
        assert _face_to_slot(2 * k) == (1, k - 1)


# ----------------------------------------------------------------------
# layout_cols block-pair rule
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("n_dims", "cols", "expected"),
    [
        (8, 4, (4, 4)),  # YM-PGT 4x4 case
        (6, 3, (4, 3)),  # 2 blocks of 3
        (6, 4, (4, 4)),  # 6 over 4 cols -> 2 blocks (with trailing empties)
        (5, 5, (2, 5)),  # 1 block of 5, equivalent to default
        (10, 5, (4, 5)),  # 2 blocks of 5
        (12, 6, (4, 6)),  # 2 blocks of 6
        (16, 4, (8, 4)),  # 4 blocks of 4 — tall layout
    ],
)
def test_grid_shape_layout_cols(
    n_dims: int, cols: int, expected: tuple[int, int]
) -> None:
    """Block-pair layout: 2 * ceil(n_dims/cols) rows x cols columns."""
    assert _grid_shape(n_dims, cols=cols) == expected


def test_grid_shape_layout_cols_none_preserves_2_by_n() -> None:
    """Default behaviour (cols=None) is unchanged from the 2 x N rule."""
    assert _grid_shape(5, cols=None) == (2, 5)
    assert _grid_shape(8, cols=None) == (2, 8)


def test_grid_shape_rejects_cols_below_1() -> None:
    with pytest.raises(ValueError, match="cols must be >= 1"):
        _grid_shape(8, cols=0)


def test_face_to_slot_layout_cols_4_for_n8() -> None:
    """The YM-PGT case: N=8 in a 4x4 block-pair layout.

    Top two rows hold axes 1-4 (axis k = (0, k-1) up, (1, k-1) down).
    Bottom two rows hold axes 5-8 (axis k = (2, k-5) up, (3, k-5) down).
    """
    from tidal.inference._atlas import _face_to_slot

    # Axes 1-4: block 0, rows 0 (up) / 1 (down)
    assert _face_to_slot(1, cols=4) == (0, 0)  # face 1 = +x_1
    assert _face_to_slot(2, cols=4) == (1, 0)  # face 2 = -x_1
    assert _face_to_slot(3, cols=4) == (0, 1)  # face 3 = +x_2
    assert _face_to_slot(7, cols=4) == (0, 3)  # face 7 = +x_4
    assert _face_to_slot(8, cols=4) == (1, 3)  # face 8 = -x_4

    # Axes 5-8: block 1, rows 2 (up) / 3 (down)
    assert _face_to_slot(9, cols=4) == (2, 0)  # face 9 = +x_5
    assert _face_to_slot(10, cols=4) == (3, 0)  # face 10 = -x_5
    assert _face_to_slot(15, cols=4) == (2, 3)  # face 15 = +x_8
    assert _face_to_slot(16, cols=4) == (3, 3)  # face 16 = -x_8


def test_face_to_slot_layout_cols_default_matches_no_cols() -> None:
    """cols=None and the no-cols call agree on every face for the 2 x N rule."""
    from tidal.inference._atlas import _face_to_slot

    for k in range(1, 9):
        for sign_idx in (2 * k - 1, 2 * k):
            assert _face_to_slot(sign_idx) == _face_to_slot(sign_idx, cols=None)


def test_face_to_slot_layout_cols_trailing_empty_block() -> None:
    """N=6, cols=4 -> first block holds axes 1-4, second block holds axes 5-6
    with the last two columns of the second block empty.
    """
    from tidal.inference._atlas import _face_to_slot

    # Axes 1-4: block 0
    assert _face_to_slot(1, cols=4) == (0, 0)
    assert _face_to_slot(8, cols=4) == (1, 3)
    # Axes 5-6: block 1, only first two columns populated
    assert _face_to_slot(9, cols=4) == (2, 0)  # +x_5
    assert _face_to_slot(10, cols=4) == (3, 0)  # -x_5
    assert _face_to_slot(11, cols=4) == (2, 1)  # +x_6
    assert _face_to_slot(12, cols=4) == (3, 1)  # -x_6


def test_face_to_slot_layout_cols_rejects_zero() -> None:
    from tidal.inference._atlas import _face_to_slot

    with pytest.raises(ValueError, match="cols must be >= 1"):
        _face_to_slot(1, cols=0)


def test_plot_atlas_layout_cols_smoke_render(tmp_path: Path) -> None:
    """End-to-end: synthetic N=4 survey with layout_cols=2 produces a 4 x 2
    grid PDF without exception.  Asserts that plot_atlas accepts the
    layout_cols kwarg and writes a non-empty PDF.
    """
    n = 4
    names = ("a", "b", "c", "d")
    Q = np.eye(n)
    survey_dir = tmp_path / "survey"
    survey_dir.mkdir()
    for face_idx, sub_tile in enumerate_cells(n_dims=n, M=1):
        cd = survey_dir / cell_label(face_idx, sub_tile)
        _make_synthetic_cell(cd, names, face_idx, sub_tile, M=1, Q=Q)

    out_pdf = survey_dir / "atlas_2cols.pdf"
    result = plot_atlas(survey_dir, out_pdf, layout_cols=2)
    assert result == out_pdf
    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


# ----------------------------------------------------------------------
# Synthetic survey -> atlas PDF
# ----------------------------------------------------------------------


def _make_synthetic_cell(
    cell_dir: Path,
    names: tuple[str, ...],
    face_idx: int,
    sub_tile: tuple[int, ...],
    M: int,
    Q: np.ndarray,
    n_samples: int = 200,
    seed: int = 0,
) -> None:
    """Write inference.json + results.csv into ``cell_dir`` for testing."""
    cell_dir.mkdir(parents=True, exist_ok=True)
    n = len(names)
    rng = np.random.default_rng(seed)

    # Synthetic samples on the sub-tile of the sphere, scaled by random r.
    from tidal.inference._sphere import tile_bounds

    u_lo, u_hi = tile_bounds(sub_tile, M)
    samples = np.zeros((n_samples, n))
    for i in range(n_samples):
        u_face = u_lo + rng.uniform(0.0, 1.0, size=n - 1) * (u_hi - u_lo)
        theta_hat = face_to_direction(face_idx, u_face, Q)
        r = float(np.exp(rng.uniform(np.log(1e-2), np.log(1e2))))
        samples[i, :] = r * theta_hat

    # inference.json
    meta = {
        "param_names": list(names),
        "priors": [
            {
                "kind": "radial_angular",
                "names": list(names),
                "r_lo": 1e-3,
                "r_hi": 1e3,
                "face_idx": face_idx,
                "sub_tile": list(sub_tile),
                "M": M,
                "Q": Q.tolist(),
            }
        ],
    }
    (cell_dir / "inference.json").write_text(json.dumps(meta))

    # results.csv
    with (cell_dir / "results.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([*names, "log_likelihood", "log_prior", "weight"])
        for i in range(n_samples):
            writer.writerow([*samples[i, :].tolist(), 0.0, 0.0, 1.0])


def test_plot_atlas_renders_pdf(tmp_path: Path) -> None:
    """End-to-end smoke: synthetic cells -> atlas.pdf with 2N panels."""
    n = 3
    names = ("a", "b", "c")
    Q = np.eye(n)
    survey_dir = tmp_path / "survey"
    survey_dir.mkdir()

    # Populate one tile per face, M=1 (12 tiles total for N=3 -> 6 faces).
    for face_idx, sub_tile in enumerate_cells(n_dims=n, M=1):
        cd = survey_dir / cell_label(face_idx, sub_tile)
        _make_synthetic_cell(cd, names, face_idx, sub_tile, M=1, Q=Q)

    out_pdf = survey_dir / "atlas.pdf"
    result = plot_atlas(survey_dir, out_pdf)
    assert result == out_pdf
    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


def test_plot_atlas_default_output_path(tmp_path: Path) -> None:
    n = 3
    names = ("a", "b", "c")
    Q = np.eye(n)
    survey_dir = tmp_path / "survey"
    survey_dir.mkdir()
    for face_idx, sub_tile in enumerate_cells(n_dims=n, M=1):
        cd = survey_dir / cell_label(face_idx, sub_tile)
        _make_synthetic_cell(cd, names, face_idx, sub_tile, M=1, Q=Q)
    out = plot_atlas(survey_dir)
    assert out == survey_dir / "atlas.pdf"
    assert out.exists()


def test_plot_atlas_handles_partial_coverage(tmp_path: Path) -> None:
    """Not every face needs a populated tile; missing faces stay blank."""
    n = 3
    names = ("a", "b", "c")
    Q = np.eye(n)
    survey_dir = tmp_path / "survey"
    survey_dir.mkdir()
    # Populate only faces 1 and 2 (out of 6).
    for face_idx in (1, 2):
        cd = survey_dir / cell_label(face_idx, (1,))
        _make_synthetic_cell(cd, names, face_idx, (1,), M=1, Q=Q)

    out_pdf = survey_dir / "atlas.pdf"
    result = plot_atlas(survey_dir, out_pdf)
    assert result == out_pdf
    assert out_pdf.exists()


def test_render_face_panel_uses_no_matplotlib_histogram(tmp_path: Path) -> None:
    """Regression: the panel renderer must not fall back to the
    pre-refactor ``np.histogram2d`` + ``contourf`` path.  Atlas now goes
    through ``_render_anesthetic_corner_into`` so panels match the
    standalone corner plot style (two-tone fills + diagonal marginals).
    """
    import re
    from pathlib import Path as _Path

    src = _Path(
        "/workspaces/torsion-gertsenshtein/tidal/inference/_atlas.py",
    ).read_text(encoding="utf-8")
    # Strip docstrings + comments so the check looks at actual call sites only.
    code_only = re.sub(r'""".*?"""', "", src, flags=re.DOTALL)
    code_only = re.sub(r"#.*", "", code_only)
    assert "histogram2d" not in code_only, (
        "_atlas.py contains a histogram2d call — refactor regressed to "
        "the manual matplotlib path"
    )
    assert "contourf" not in code_only, (
        "_atlas.py contains a contourf call — atlas should render via "
        "anesthetic's KDE pipeline, not raw matplotlib contours"
    )
    assert "gaussian_filter" not in code_only, (
        "_atlas.py contains a scipy.ndimage.gaussian_filter call — "
        "atlas should defer KDE smoothing to anesthetic"
    )


def test_plot_atlas_calls_helper_with_ticks_none(tmp_path: Path) -> None:
    """Atlas panels must call the composable helper with ``ticks=None``
    so numeric tick labels are suppressed (the face-local LaTeX axis
    labels remain).  We mock-patch the helper to capture its kwargs.
    """
    n = 3
    names = ("a", "b", "c")
    Q = np.eye(n)
    survey_dir = tmp_path / "survey"
    survey_dir.mkdir()
    for face_idx, sub_tile in enumerate_cells(n_dims=n, M=1):
        cd = survey_dir / cell_label(face_idx, sub_tile)
        _make_synthetic_cell(cd, names, face_idx, sub_tile, M=1, Q=Q)

    import tidal.inference._visualize as _viz

    real_helper = _viz._render_anesthetic_corner_into
    seen_kwargs: list[dict[str, object]] = []

    def _spy(*args: object, **kwargs: object) -> object:
        seen_kwargs.append(dict(kwargs))
        return real_helper(*args, **kwargs)  # type: ignore[arg-type]

    _viz._render_anesthetic_corner_into = _spy  # type: ignore[attr-defined]
    try:
        plot_atlas(survey_dir, survey_dir / "atlas.pdf")
    finally:
        _viz._render_anesthetic_corner_into = real_helper  # type: ignore[attr-defined]

    # The helper should have been called at least once per face (6 for N=3).
    assert len(seen_kwargs) >= n_faces(n)
    # Every call must have ticks=None and target_fig set (composable path).
    for kw in seen_kwargs:
        assert kw["ticks"] is None
        assert kw["target_fig"] is not None
        # Face-aware LaTeX labels must be wired in.
        labels = kw["labels"]
        assert isinstance(labels, dict)
        for v in labels.values():
            assert isinstance(v, str)
            assert "\\chi" in v
            assert "\\uparrow" in v or "\\downarrow" in v


def test_plot_atlas_rejects_empty_dir(tmp_path: Path) -> None:
    survey_dir = tmp_path / "empty"
    survey_dir.mkdir()
    with pytest.raises(ValueError, match="no <face_label>"):
        plot_atlas(survey_dir, survey_dir / "atlas.pdf")


def test_plot_atlas_rejects_dir_without_metadata(tmp_path: Path) -> None:
    """A subdir matching the cell-name regex but missing inference.json
    is skipped; if no usable cells remain, raise.
    """
    survey_dir = tmp_path / "broken"
    survey_dir.mkdir()
    (survey_dir / "01p_tile1").mkdir()
    (survey_dir / "01p_tile1" / "results.csv").write_text("a,b,c\n0.1,0.2,0.3\n")
    with pytest.raises(ValueError, match="no usable cell metadata"):
        plot_atlas(survey_dir, survey_dir / "atlas.pdf")
