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
    ("faces", "expected"),
    [(4, (2, 2)), (6, (2, 3)), (8, (2, 4)), (10, (2, 5)), (12, (3, 4))],
)
def test_grid_shape_known_layouts(faces: int, expected: tuple[int, int]) -> None:
    assert _grid_shape(faces) == expected


def test_grid_shape_generic_covers_all_faces() -> None:
    for f in (2, 18, 20, 22):
        rows, cols = _grid_shape(f)
        assert rows * cols >= f


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
