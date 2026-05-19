"""Tests for :class:`RadialAngularPrior` + mixed-prior dispatch in
:mod:`tidal.inference._prior`.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tidal.inference._prior import (
    Prior,
    RadialAngularPrior,
    build_prior_transform,
    parse_joint_prior,
    prior_param_names,
    total_prior_ndim,
    validate_prior_names,
)
from tidal.inference._sphere import face_to_direction, random_rotation, tile_bounds


def _make_prior(names: tuple[str, ...] = ("a", "b", "c", "d")) -> RadialAngularPrior:
    return RadialAngularPrior(
        names=names,
        r_lo=1e-3,
        r_hi=1e3,
        face_idx=1,
        sub_tile=(1,) * (len(names) - 1),
        M=1,
        Q=np.eye(len(names)),
    )


# ----------------------------------------------------------------------
# Construction & validation
# ----------------------------------------------------------------------


def test_post_init_rejects_too_few_names() -> None:
    with pytest.raises(ValueError, match="at least 2 couplings"):
        RadialAngularPrior(
            names=("a",),
            r_lo=1e-3,
            r_hi=1e3,
            face_idx=1,
            sub_tile=(),
            M=1,
            Q=np.eye(1),
        )


def test_post_init_rejects_negative_r() -> None:
    with pytest.raises(ValueError, match="r bounds must be positive"):
        RadialAngularPrior(
            names=("a", "b"),
            r_lo=-1.0,
            r_hi=1e3,
            face_idx=1,
            sub_tile=(1,),
            M=1,
            Q=np.eye(2),
        )


def test_post_init_rejects_swapped_r_bounds() -> None:
    """r_lo strictly greater than r_hi is rejected (equal bounds are allowed —
    that's fixed-radius mode; see :func:`test_post_init_accepts_equal_r_bounds`).
    """
    with pytest.raises(ValueError, match="r_lo must be <= r_hi"):
        RadialAngularPrior(
            names=("a", "b"),
            r_lo=10.0,
            r_hi=1.0,
            face_idx=1,
            sub_tile=(1,),
            M=1,
            Q=np.eye(2),
        )


def test_post_init_accepts_equal_r_bounds() -> None:
    """``r_lo == r_hi`` is accepted as fixed-radius mode."""
    p = RadialAngularPrior(
        names=("a", "b"),
        r_lo=1.5,
        r_hi=1.5,
        face_idx=1,
        sub_tile=(1,),
        M=1,
        Q=np.eye(2),
    )
    assert p.is_fixed_radius
    assert p.r_lo == 1.5
    assert p.r_hi == 1.5


def test_is_fixed_radius_false_for_variable_r() -> None:
    p = _make_prior()
    assert not p.is_fixed_radius


def test_transform_fixed_radius_ignores_u0() -> None:
    """In fixed-radius mode the magnitude is independent of ``u[0]``."""
    p = RadialAngularPrior(
        names=("a", "b", "c"),
        r_lo=0.01,
        r_hi=0.01,
        face_idx=1,
        sub_tile=(1, 1),
        M=1,
        Q=np.eye(3),
    )
    c1 = p.transform(np.array([0.1, 0.5, 0.5]))
    c2 = p.transform(np.array([0.9, 0.5, 0.5]))
    assert math.isclose(float(np.linalg.norm(c1)), 0.01, rel_tol=1e-12)
    assert math.isclose(float(np.linalg.norm(c2)), 0.01, rel_tol=1e-12)
    # u[0] is discarded entirely — angular outputs must match bit-for-bit.
    np.testing.assert_array_equal(c1, c2)


def test_parse_joint_prior_fixed_radius() -> None:
    """The CLI parser accepts ``r_lo == r_hi`` and produces a fixed-radius prior."""
    spec = "names=a,b,c;type=cubed_sphere;M=1;face=1;sub=1_1;r_lo=1.5;r_hi=1.5"
    p = parse_joint_prior(spec)
    assert isinstance(p, RadialAngularPrior)
    assert p.is_fixed_radius
    assert p.r_lo == 1.5
    assert p.r_hi == 1.5
    # Round-trip: transform delivers magnitude exactly r_lo regardless of u[0].
    c = p.transform(np.array([0.3, 0.5, 0.5]))
    assert math.isclose(float(np.linalg.norm(c)), 1.5, rel_tol=1e-12)


def test_post_init_rejects_wrong_subtile_length() -> None:
    with pytest.raises(ValueError, match="sub_tile must have length"):
        RadialAngularPrior(
            names=("a", "b", "c"),
            r_lo=1e-3,
            r_hi=1e3,
            face_idx=1,
            sub_tile=(1,),
            M=1,
            Q=np.eye(3),
        )


def test_post_init_rejects_wrong_q_shape() -> None:
    with pytest.raises(ValueError, match="Q must have shape"):
        RadialAngularPrior(
            names=("a", "b", "c"),
            r_lo=1e-3,
            r_hi=1e3,
            face_idx=1,
            sub_tile=(1, 1),
            M=1,
            Q=np.eye(2),
        )


def test_n_dims_matches_names_length() -> None:
    p = _make_prior(("a", "b", "c", "d", "e", "f"))
    assert p.n_dims == 6


# ----------------------------------------------------------------------
# Transform: cube -> physical coupling
# ----------------------------------------------------------------------


def test_transform_returns_correct_shape() -> None:
    p = _make_prior(("a", "b", "c", "d"))
    c = p.transform(np.full(4, 0.5))
    assert c.shape == (4,)


def test_transform_geometric_centre() -> None:
    """``u = (0.5, 0.5, ..., 0.5)`` on M=1 face 1+ projects to face centre.

    At u[0] = 0.5 the magnitude is the geometric mean of r_lo and r_hi
    (log_uniform midpoint).  At u[1:] = 0.5 the face-cube position is
    ``u_face = -1 + 1 * 1 = 0`` -- i.e. the face centre, which projects
    to ``+e_1`` for face 1+.  So c = r_mid * +e_1.
    """
    p = _make_prior(("a", "b", "c"))  # N = 3
    c = p.transform(np.array([0.5, 0.5, 0.5]))
    r_mid = math.sqrt(p.r_lo * p.r_hi)  # geometric mean
    expected = np.array([r_mid, 0.0, 0.0])
    np.testing.assert_allclose(c, expected, atol=1e-12)


def test_transform_magnitude_scales_with_u0() -> None:
    """Varying u[0] only rescales the whole vector (direction unchanged)."""
    p = _make_prior(("a", "b", "c"))
    u_dir = np.array([0.3, 0.7])  # fixed direction part
    c0 = p.transform(np.array([0.0, *u_dir]))
    c1 = p.transform(np.array([1.0, *u_dir]))
    # c1 / c0 should be uniform across components (same direction, different r)
    ratio = c1 / c0
    np.testing.assert_allclose(ratio, ratio[0] * np.ones(3), atol=1e-12)
    # Magnitude ratio = r_hi / r_lo
    np.testing.assert_allclose(ratio[0], p.r_hi / p.r_lo, atol=1e-12)


def test_transform_unit_direction_independent_of_r() -> None:
    """``c / |c|`` should depend only on u[1:], not on u[0]."""
    p = _make_prior(("a", "b", "c", "d"))
    u_dir = np.array([0.2, 0.4, 0.6])
    dirs = []
    for u0 in (0.0, 0.3, 0.7, 1.0):
        c = p.transform(np.array([u0, *u_dir]))
        dirs.append(c / np.linalg.norm(c))
    for d in dirs[1:]:
        np.testing.assert_allclose(d, dirs[0], atol=1e-12)


def test_transform_lands_in_correct_tile() -> None:
    """A sample on face k+ at sub-tile s with M-subdivision must lie inside
    that tile's gnomonic image on the sphere.
    """
    n = 4
    p = RadialAngularPrior(
        names=("a", "b", "c", "d"),
        r_lo=0.1,
        r_hi=10.0,
        face_idx=3,  # +x_2
        sub_tile=(2, 1, 2),
        M=2,
        Q=np.eye(n),
    )
    rng = np.random.default_rng(13)
    for _ in range(20):
        u = rng.uniform(0.0, 1.0, size=n)
        c = p.transform(u)
        c_unit = c / np.linalg.norm(c)
        # Expected face-local cube coords of the sphere point:
        u_lo, u_hi = tile_bounds((2, 1, 2), 2)
        u_face_expected = u_lo + u[1:] * (u_hi - u_lo)
        # Recompute via face_to_direction and check the unit direction matches
        c_unit_expected = face_to_direction(3, u_face_expected, np.eye(n))
        np.testing.assert_allclose(c_unit, c_unit_expected, atol=1e-12)


def test_transform_with_random_q_orthogonal() -> None:
    """Pre-rotation Q only changes the basis; magnitudes preserved."""
    n = 5
    Q = random_rotation(n, seed=99)
    p = RadialAngularPrior(
        names=tuple("abcde"),
        r_lo=0.5,
        r_hi=2.0,
        face_idx=2,
        sub_tile=(1,) * (n - 1),
        M=1,
        Q=Q,
    )
    rng = np.random.default_rng(0)
    u = rng.uniform(0.0, 1.0, size=n)
    c = p.transform(u)
    # Magnitude lies in [r_lo, r_hi]
    r = float(np.linalg.norm(c))
    assert p.r_lo <= r <= p.r_hi


# ----------------------------------------------------------------------
# CLI parsing
# ----------------------------------------------------------------------


def test_parse_joint_prior_minimal() -> None:
    spec = (
        "names=alpha1,alpha2,alpha3,delta1;"
        "type=cubed_sphere;M=2;face=1;sub=1_1_1;"
        "r_lo=1e-3;r_hi=1e3"
    )
    p = parse_joint_prior(spec)
    assert p.names == ("alpha1", "alpha2", "alpha3", "delta1")
    assert p.r_lo == 1e-3
    assert p.r_hi == 1e3
    assert p.face_idx == 1
    assert p.sub_tile == (1, 1, 1)
    assert p.M == 2
    np.testing.assert_array_equal(p.Q, np.eye(4))


def test_parse_joint_prior_random_q() -> None:
    spec = (
        "names=a,b,c;type=cubed_sphere;M=1;face=1;sub=1_1;"
        "r_lo=0.01;r_hi=100;Q=random:42"
    )
    p = parse_joint_prior(spec)
    np.testing.assert_allclose(p.Q, random_rotation(3, seed=42), atol=1e-12)


def test_parse_joint_prior_explicit_identity_q() -> None:
    spec = "names=a,b;type=cubed_sphere;M=1;face=1;sub=1;r_lo=0.01;r_hi=100;Q=identity"
    p = parse_joint_prior(spec)
    np.testing.assert_array_equal(p.Q, np.eye(2))


def test_parse_joint_prior_rejects_missing_field() -> None:
    spec = "names=a,b;type=cubed_sphere;M=1;face=1;sub=1"  # no r_lo / r_hi
    with pytest.raises(ValueError, match="missing required field"):
        parse_joint_prior(spec)


def test_parse_joint_prior_rejects_unknown_type() -> None:
    spec = "names=a,b;type=spherical_polar;M=1;face=1;sub=1;r_lo=0.01;r_hi=100"
    with pytest.raises(ValueError, match="only type=cubed_sphere"):
        parse_joint_prior(spec)


def test_parse_joint_prior_rejects_bad_q() -> None:
    spec = "names=a,b;type=cubed_sphere;M=1;face=1;sub=1;r_lo=0.01;r_hi=100;Q=fancy"
    with pytest.raises(ValueError, match="Q must be"):
        parse_joint_prior(spec)


def test_parse_joint_prior_rejects_no_kv_field() -> None:
    spec = "names=a,b;type=cubed_sphere;M=1;face=1;sub=1;notakv;r_lo=0.01;r_hi=1"
    with pytest.raises(ValueError, match="key=value"):
        parse_joint_prior(spec)


# ----------------------------------------------------------------------
# Mixed-prior dispatch
# ----------------------------------------------------------------------


def test_total_prior_ndim_per_param_only() -> None:
    priors = [
        Prior(name="alpha", distribution="uniform", low=0.0, high=1.0),
        Prior(name="beta", distribution="log_uniform", low=0.1, high=10.0),
    ]
    assert total_prior_ndim(priors) == 2


def test_total_prior_ndim_joint_only() -> None:
    p = _make_prior(("a", "b", "c", "d"))
    assert total_prior_ndim([p]) == 4


def test_total_prior_ndim_mixed() -> None:
    p = _make_prior(("a", "b", "c"))
    extra = Prior(name="kappa", distribution="log_uniform", low=0.5, high=2.0)
    assert total_prior_ndim([p, extra]) == 4
    assert total_prior_ndim([extra, p]) == 4


def test_prior_param_names_mixed_order_preserved() -> None:
    p = _make_prior(("alpha", "beta", "gamma"))
    extra = Prior(name="kappa", distribution="log_uniform", low=0.5, high=2.0)
    assert prior_param_names([p, extra]) == ["alpha", "beta", "gamma", "kappa"]
    assert prior_param_names([extra, p]) == ["kappa", "alpha", "beta", "gamma"]


def test_validate_prior_names_rejects_duplicates() -> None:
    p = _make_prior(("alpha", "beta", "gamma"))
    extra = Prior(name="alpha", distribution="uniform", low=0.0, high=1.0)
    with pytest.raises(ValueError, match="duplicate"):
        validate_prior_names([p, extra])


def test_validate_prior_names_accepts_disjoint() -> None:
    p = _make_prior(("alpha", "beta", "gamma"))
    extra = Prior(name="kappa", distribution="log_uniform", low=0.5, high=2.0)
    validate_prior_names([p, extra])  # no raise


def test_build_prior_transform_per_param_only() -> None:
    """Legacy code path: per-parameter Priors only — output unchanged."""
    priors = [
        Prior(name="a", distribution="uniform", low=0.0, high=10.0),
        Prior(name="b", distribution="uniform", low=-1.0, high=1.0),
    ]
    transform = build_prior_transform(priors)
    out = transform(np.array([0.3, 0.5]))
    np.testing.assert_allclose(out, [3.0, 0.0])


def test_build_prior_transform_joint_only() -> None:
    p = _make_prior(("a", "b", "c"))
    transform = build_prior_transform([p])
    u = np.array([0.5, 0.5, 0.5])
    out = transform(u)
    expected = p.transform(u)
    np.testing.assert_allclose(out, expected)
    assert out.shape == (3,)


def test_build_prior_transform_mixed_order_independent() -> None:
    """Mixed: a per-param prior plus a joint prior; verify column ordering."""
    p = _make_prior(("alpha", "beta", "gamma"))
    extra = Prior(name="kappa", distribution="uniform", low=0.0, high=10.0)

    # joint first, then per-param
    transform_jp = build_prior_transform([p, extra])
    u_jp = np.array([0.5, 0.5, 0.5, 0.7])  # u[0:3] -> joint, u[3] -> kappa
    out_jp = transform_jp(u_jp)
    assert out_jp.shape == (4,)
    np.testing.assert_allclose(out_jp[:3], p.transform(u_jp[:3]))
    np.testing.assert_allclose(out_jp[3], 7.0)

    # per-param first, then joint
    transform_pj = build_prior_transform([extra, p])
    u_pj = np.array([0.7, 0.5, 0.5, 0.5])
    out_pj = transform_pj(u_pj)
    assert out_pj.shape == (4,)
    np.testing.assert_allclose(out_pj[0], 7.0)
    np.testing.assert_allclose(out_pj[1:], p.transform(u_pj[1:]))


def test_build_prior_transform_unchanged_for_legacy_callers() -> None:
    """Existing per-parameter Prior callers must keep working unchanged."""
    priors_legacy = [
        Prior(name="alpha", distribution="uniform", low=-1.0, high=1.0),
        Prior(name="beta", distribution="log_uniform", low=1e-2, high=1e2),
    ]
    transform = build_prior_transform(priors_legacy)
    out = transform(np.array([0.5, 0.5]))
    # Reproduces existing behaviour exactly
    np.testing.assert_allclose(out[0], 0.0, atol=1e-12)
    np.testing.assert_allclose(out[1], 1.0, atol=1e-12)
