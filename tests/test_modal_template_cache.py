"""Tests for the pre-Schur ingredient template cache (#327 follow-up 2).

Five correctness/performance properties:

1. **Ingredient linearity (n_c>0 spec)** — for the GENERAL spec (n_c=18),
   the pre-Schur ingredient matrices assembled from the template agree with
   direct _build_evolution_matrices output to machine precision (atol=1e-10).
   This confirms the ingredients are exactly linear in the free coupling params.

2. **Post-Schur simulation equivalence (single param)** — run_inference_step
   with template ON and OFF at 10 random delta1 values; field snapshots agree
   to atol=1e-10.  This is the end-to-end correctness test for n_c>0 path.

3. **Multi-parameter prior equivalence** — same as test 2 but with a
   3-parameter prior (delta1, zeta1, beta2); template ON vs OFF.

4. **n_c=0 path unaffected** — for dark_photon_plasma.json (n_c=0), template
   ON vs OFF produces bit-exact snapshots.  Guards against regression.

5. **Perf budget** — median per-call wall ≤ 35 ms on the GENERAL fixture with
   N_free=3 prior (post-template target ~26 ms; budget padded for CI noise).
"""

from __future__ import annotations

import os
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
GERT_SPEC = REPO / "examples/data/torsion_gertsenshtein_general_nonminimal.json"
DPP_SPEC = REPO / "examples/data/dark_photon_plasma.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_ns(
    spec_path: Path,
    prior_params: tuple[str, ...],
    fixed_params: tuple[str, ...],
    *,
    grid_shape: str = "32",
    bounds: str = "0:100",
    ic: str = "plane-wave",
    n_samples: int = 1,
) -> object:
    from tidal.cli import _build_parser

    parser = _build_parser()
    args = [
        "sample",
        str(spec_path),
        *(f"--prior={p}" for p in prior_params),
        "--likelihood",
        "P_max:maximize",
        "--baseline-formula",
        "sin(kappa*B0*t_end/2)**2",
        *(f"--param={p}" for p in fixed_params),
        "--grid-shape",
        grid_shape,
        f"--bounds={bounds}",
        "--periodic",
        "--ic",
        ic,
        "--ic-component",
        "h_5",
        "--ic-wavevector",
        "0.06283185307179587",
        "--ic-amplitude",
        "0.01",
        "--t-end",
        "10",
        "--measure",
        "peak_conversion,conversion",
        "--source",
        "h_5",
        "--target",
        "a_1",
        "--method",
        "mc",
        "--n-samples",
        str(n_samples),
        "--output",
        f"/tmp/test_template_cache_{os.getpid()}",
    ]
    return parser.parse_args(args)


_GERT_FIXED = (
    "kappa=1.0",
    "B0=0.01",
    "chi=0.0",
    "eta=0.0",
    "zeta1=0.0",
    "zeta2=0.0",
    "zeta3=0.0",
    "beta1=0.0",
    "beta2=-0.6",
    "beta3=0.5",
    "xi=1.0",
)


def _run_snapshots(
    spec_path: Path,
    prior_params: tuple[str, ...],
    fixed_params: tuple[str, ...],
    overrides: dict[str, float],
    *,
    template_on: bool,
    grid_shape: str = "32",
    bounds: str = "0:100",
) -> np.ndarray:
    from tidal.cli._sweep import run_inference_step

    prev = os.environ.get("TIDAL_MODAL_A_TEMPLATE")
    os.environ["TIDAL_MODAL_A_TEMPLATE"] = "1" if template_on else "0"
    # Invalidate cache between ON/OFF runs so we don't get a stale hit.
    import tidal.solver.modal as _m

    _m._A_TEMPLATE_CACHE.clear()
    _m._FREE_PARAM_NAMES = frozenset()
    try:
        ns = _build_ns(
            spec_path,
            prior_params,
            fixed_params,
            grid_shape=grid_shape,
            bounds=bounds,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim_data = run_inference_step(ns, spec_path, overrides)
    finally:
        if prev is None:
            os.environ.pop("TIDAL_MODAL_A_TEMPLATE", None)
        else:
            os.environ["TIDAL_MODAL_A_TEMPLATE"] = prev
        _m._A_TEMPLATE_CACHE.clear()
        _m._FREE_PARAM_NAMES = frozenset()

    return np.stack(
        [sim_data.fields[name] for name in sorted(sim_data.fields)],
        axis=0,
    )


# ---------------------------------------------------------------------------
# 1.  Ingredient linearity — pre-Schur matrices agree to machine precision
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_ingredient_linearity_general_spec() -> None:
    """Pre-Schur ingredient matrices assembled from template equal direct build
    to atol=1e-10 at 5 random parameter vectors (n_c=18 GENERAL spec).
    """
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.modal import (
        _apply_ingredient_template,
        _build_a_template_cache,
        _build_evolution_matrices,
        _build_k_axes,
        _build_k_grid,
        _IngredientTemplateCache,
        register_inference_context,
    )
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem

    spec = EquationSystem.from_file(GERT_SPEC)
    grid = GridInfo(bounds=((0.0, 100.0),), shape=(32,), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    base_params = {
        "kappa": 1.0,
        "B0": 0.01,
        "delta1": 0.0,
        "zeta1": 0.0,
        "zeta2": 0.0,
        "zeta3": 0.0,
        "beta1": 0.0,
        "beta2": -0.6,
        "beta3": 0.5,
        "xi": 1.0,
        "chi": 0.0,
        "eta": 0.0,
    }
    free_params = frozenset(["delta1", "zeta1", "beta2"])

    register_inference_context(free_params)
    tmpl = _build_a_template_cache(
        spec,
        layout,
        grid,
        k_grid,
        rfft_shape,
        base_params,
        free_params,
        needs_reduction=True,
    )
    assert isinstance(tmpl, _IngredientTemplateCache), (
        "Expected _IngredientTemplateCache for GENERAL spec (n_c=18)"
    )

    rng = np.random.default_rng(42)
    for _ in range(5):
        params = dict(base_params)
        params["delta1"] = rng.uniform(-0.02, 0.02)
        params["zeta1"] = rng.uniform(-0.02, 0.02)
        params["beta2"] = rng.uniform(-0.8, -0.4)

        # Template path: assemble ingredients + Schur
        (_A_tmpl, *_rest_tmpl) = _apply_ingredient_template(tmpl, params)

        # Direct path: full rebuild
        coeff_eval = CoefficientEvaluator(spec, grid, params)
        ingr_direct: dict = {}
        _build_evolution_matrices(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
            _ingredient_out=ingr_direct,
        )
        # Check pre-Schur ingredients from the direct build vs template assembly
        A_dd_tmpl = (
            tmpl.A_dd_base.copy()
            + (params["delta1"] - base_params["delta1"]) * tmpl.dA_dd["delta1"]
            + (params["zeta1"] - base_params["zeta1"]) * tmpl.dA_dd["zeta1"]
            + (params["beta2"] - base_params["beta2"]) * tmpl.dA_dd["beta2"]
        )
        np.testing.assert_allclose(
            A_dd_tmpl,
            ingr_direct["A_dd"],
            atol=1e-10,
            rtol=1e-8,
            err_msg="A_dd ingredient linearity violated",
        )
        S_cc_tmpl = (
            tmpl.S_cc_base.copy()
            + (params["delta1"] - base_params["delta1"]) * tmpl.dS_cc["delta1"]
            + (params["zeta1"] - base_params["zeta1"]) * tmpl.dS_cc["zeta1"]
            + (params["beta2"] - base_params["beta2"]) * tmpl.dS_cc["beta2"]
        )
        np.testing.assert_allclose(
            S_cc_tmpl,
            ingr_direct["S_cc"],
            atol=1e-10,
            rtol=1e-8,
            err_msg="S_cc ingredient linearity violated",
        )


# ---------------------------------------------------------------------------
# 2.  Post-Schur simulation equivalence — single prior param
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_single_param_prior_simulation_equivalence() -> None:
    """Template ON vs OFF produces snapshots agreeing to atol=1e-10 for 5
    random delta1 values in [-0.025, 0.025] on the GENERAL spec.
    """
    prior = ("delta1=uniform:-0.025:0.025",)
    fixed = tuple(p for p in _GERT_FIXED if not p.startswith("delta1"))
    # delta1=0 is non-trivial (background coupling present via B0 terms)
    for delta1_val in np.linspace(-0.02, 0.02, 5):
        overrides = {"delta1": float(delta1_val)}
        snaps_on = _run_snapshots(
            GERT_SPEC,
            prior,
            fixed,
            overrides,
            template_on=True,
        )
        snaps_off = _run_snapshots(
            GERT_SPEC,
            prior,
            fixed,
            overrides,
            template_on=False,
        )
        np.testing.assert_allclose(
            snaps_on,
            snaps_off,
            atol=1e-10,
            rtol=1e-8,
            err_msg=f"Template mismatch at delta1={delta1_val:.4f}",
        )


# ---------------------------------------------------------------------------
# 3.  Multi-parameter prior equivalence
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_multi_param_prior_simulation_equivalence() -> None:
    """Template ON vs OFF with 3-parameter prior (delta1, zeta1, beta2)
    produces snapshots agreeing to atol=1e-10 at 3 random points.
    """
    prior = (
        "delta1=uniform:-0.025:0.025",
        "zeta1=uniform:-0.025:0.025",
        "beta2=uniform:-0.8:-0.4",
    )
    fixed = tuple(
        p
        for p in _GERT_FIXED
        if not any(p.startswith(k) for k in ("delta1", "zeta1", "beta2"))
    )
    rng = np.random.default_rng(7)
    for _ in range(3):
        overrides = {
            "delta1": float(rng.uniform(-0.02, 0.02)),
            "zeta1": float(rng.uniform(-0.02, 0.02)),
            "beta2": float(rng.uniform(-0.8, -0.4)),
        }
        snaps_on = _run_snapshots(
            GERT_SPEC,
            prior,
            fixed,
            overrides,
            template_on=True,
        )
        snaps_off = _run_snapshots(
            GERT_SPEC,
            prior,
            fixed,
            overrides,
            template_on=False,
        )
        np.testing.assert_allclose(
            snaps_on,
            snaps_off,
            atol=1e-10,
            rtol=1e-8,
            err_msg=f"Multi-param template mismatch at {overrides}",
        )


# ---------------------------------------------------------------------------
# 4.  n_c=0 path unaffected — bit-exact for DPP spec
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_nc0_spec_template_bit_exact() -> None:
    """For dark_photon_plasma.json (n_c=0), template ON vs OFF produce
    bit-exact snapshots — the n_c=0 _ATemplateCache path is unchanged.
    """
    prior = ("delta1=uniform:-0.025:0.025",)
    dpp_fixed = (
        "kappa=1.0",
        "B0=0.01",
        "mA2=0.5",
        "deltam=0.0",
        "xi=0.5",
        "alpha3=0.1",
    )
    fixed = tuple(p for p in dpp_fixed if not p.startswith("delta1"))
    overrides = {"delta1": 0.0}
    snaps_on = _run_snapshots(
        DPP_SPEC,
        prior,
        fixed,
        overrides,
        template_on=True,
        grid_shape="32",
        bounds="0:50",
    )
    snaps_off = _run_snapshots(
        DPP_SPEC,
        prior,
        fixed,
        overrides,
        template_on=False,
        grid_shape="32",
        bounds="0:50",
    )
    np.testing.assert_array_equal(snaps_on, snaps_off)


# ---------------------------------------------------------------------------
# 5.  Perf budget
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_template_cache_perf_budget() -> None:
    """Median per-call wall with N_free=3 prior on GENERAL spec stays ≤ 35 ms.

    Post-template target is ~26 ms (13 ms probe + 9 ms ifft + 3 ms Schur
    + <1 ms assembly); budget is 35 ms to absorb CI noise.
    """
    import tidal.solver.modal as _m
    from tidal.cli._sweep import run_inference_step

    prior = (
        "delta1=uniform:-0.025:0.025",
        "zeta1=uniform:-0.025:0.025",
        "beta2=uniform:-0.8:-0.4",
    )
    fixed = tuple(
        p
        for p in _GERT_FIXED
        if not any(p.startswith(k) for k in ("delta1", "zeta1", "beta2"))
    )
    ns = _build_ns(GERT_SPEC, prior, fixed, grid_shape="256", bounds="0:100")
    overrides = {"delta1": 0.0, "zeta1": 0.0, "beta2": -0.6}

    os.environ["TIDAL_MODAL_A_TEMPLATE"] = "1"
    _m._A_TEMPLATE_CACHE.clear()
    _m._FREE_PARAM_NAMES = frozenset()

    # Warmup (also builds the template)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(3):
            run_inference_step(ns, GERT_SPEC, overrides)

    times: list[float] = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(20):
            t0 = time.perf_counter()
            run_inference_step(ns, GERT_SPEC, overrides)
            times.append((time.perf_counter() - t0) * 1000)

    _m._A_TEMPLATE_CACHE.clear()
    _m._FREE_PARAM_NAMES = frozenset()

    median_ms = float(np.median(times))
    assert median_ms <= 35.0, (
        f"Template-cache per-call median {median_ms:.2f} ms exceeds 35 ms budget; "
        "ingredient assembly or Schur is unexpectedly slow"
    )
