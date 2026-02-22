"""Tests for tidal.solver.grid.GridInfo."""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.grid import GridInfo


class TestGridInfoConstruction:
    def test_1d(self) -> None:
        g = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        assert g.ndim == 1
        assert g.num_points == 64
        assert g.dx == pytest.approx((10 / 64,))

    def test_2d(self) -> None:
        g = GridInfo(bounds=((0, 2), (0, 3)), shape=(4, 6), periodic=(True, False))
        assert g.ndim == 2
        assert g.num_points == 24
        assert g.dx == pytest.approx((0.5, 0.5))

    def test_3d(self) -> None:
        g = GridInfo(bounds=((0, 1), (0, 2), (0, 3)), shape=(4, 5, 6), periodic=(False, True, False))
        assert g.ndim == 3
        assert g.num_points == 120
        assert g.dx == pytest.approx((0.25, 0.4, 0.5))


class TestGridInfoValidation:
    def test_bounds_shape_mismatch(self) -> None:
        with pytest.raises(ValueError, match=r"bounds has 2.*shape has 1"):
            GridInfo(bounds=((0, 1), (0, 2)), shape=(10,), periodic=(True, True))

    def test_bounds_periodic_mismatch(self) -> None:
        with pytest.raises(ValueError, match=r"bounds has 1.*periodic has 2"):
            GridInfo(bounds=((0, 1),), shape=(10,), periodic=(True, True))

    def test_upper_not_greater(self) -> None:
        with pytest.raises(ValueError, match="upper must exceed lower"):
            GridInfo(bounds=((5, 3),), shape=(10,), periodic=(True,))

    def test_shape_too_small(self) -> None:
        with pytest.raises(ValueError, match="need at least 2"):
            GridInfo(bounds=((0, 1),), shape=(1,), periodic=(True,))


class TestGridInfoBC:
    def test_default_bc_none(self) -> None:
        g = GridInfo(bounds=((0, 10),), shape=(20,), periodic=(True,))
        assert g.bc is None

    def test_explicit_bc(self) -> None:
        g = GridInfo(bounds=((0, 10),), shape=(20,), periodic=(True,), bc=("periodic",))
        assert g.bc == ("periodic",)

    def test_effective_bc_from_periodic(self) -> None:
        g = GridInfo(bounds=((0, 1), (0, 2)), shape=(5, 10), periodic=(True, False))
        assert g.effective_bc == ("periodic", "neumann")

    def test_effective_bc_from_explicit(self) -> None:
        g = GridInfo(
            bounds=((0, 1), (0, 2)),
            shape=(5, 10),
            periodic=(False, False),
            bc=("dirichlet", "neumann"),
        )
        assert g.effective_bc == ("dirichlet", "neumann")

    def test_bc_length_mismatch(self) -> None:
        with pytest.raises(ValueError, match="bc has 1 entries"):
            GridInfo(bounds=((0, 1), (0, 2)), shape=(5, 10), periodic=(True, True), bc=("periodic",))

    def test_bc_invalid_type(self) -> None:
        with pytest.raises(ValueError, match="must be one of"):
            GridInfo(bounds=((0, 1),), shape=(10,), periodic=(True,), bc=("invalid",))

    def test_bc_periodic_inconsistency(self) -> None:
        with pytest.raises(ValueError, match="must be consistent"):
            GridInfo(bounds=((0, 1),), shape=(10,), periodic=(False,), bc=("periodic",))


class TestGridInfoCoords:
    def test_1d_cell_coords(self) -> None:
        g = GridInfo(bounds=((0, 10),), shape=(4,), periodic=(True,))
        cc = g.cell_coords
        assert cc.shape == (4, 1)
        np.testing.assert_allclose(cc[:, 0], [1.25, 3.75, 6.25, 8.75])

    def test_2d_cell_coords_shape(self) -> None:
        g = GridInfo(bounds=((0, 2), (0, 3)), shape=(4, 3), periodic=(True, True))
        cc = g.cell_coords
        assert cc.shape == (4, 3, 2)

    def test_axes_coords(self) -> None:
        g = GridInfo(bounds=((0, 2), (0, 3)), shape=(4, 3), periodic=(True, True))
        x = g.axes_coords(0)
        y = g.axes_coords(1)
        assert x.shape == (4,)
        assert y.shape == (3,)
        np.testing.assert_allclose(x, [0.25, 0.75, 1.25, 1.75])
        np.testing.assert_allclose(y, [0.5, 1.5, 2.5])

    def test_axes_coords_out_of_range(self) -> None:
        g = GridInfo(bounds=((0, 1),), shape=(10,), periodic=(True,))
        with pytest.raises(IndexError, match="axis 1 out of range"):
            g.axes_coords(1)

    def test_coord_arrays(self) -> None:
        g = GridInfo(bounds=((0, 2), (0, 3)), shape=(4, 3), periodic=(True, True))
        xs, ys = g.coord_arrays()
        assert xs.shape == (4, 3)
        assert ys.shape == (4, 3)
        # x constant along columns
        np.testing.assert_allclose(xs[:, 0], xs[:, 1])
        # y constant along rows
        np.testing.assert_allclose(ys[0, :], ys[1, :])
        # coord_arrays consistent with cell_coords
        cc = g.cell_coords
        np.testing.assert_allclose(cc[..., 0], xs)
        np.testing.assert_allclose(cc[..., 1], ys)

    def test_3d_coord_arrays(self) -> None:
        g = GridInfo(bounds=((0, 1), (0, 2), (0, 3)), shape=(4, 5, 6), periodic=(False, True, False))
        xs, _ys, zs = g.coord_arrays()
        assert xs.shape == (4, 5, 6)
        # z constant along axes 0, 1
        np.testing.assert_allclose(zs[0, 0, :], zs[1, 2, :])


class TestGridInfoFrozen:
    def test_immutable(self) -> None:
        g = GridInfo(bounds=((0, 1),), shape=(10,), periodic=(True,))
        with pytest.raises(AttributeError):
            g.shape = (20,)  # type: ignore[misc]
