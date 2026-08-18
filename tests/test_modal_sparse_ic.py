"""Regression tests for the sparse-IC mode-skip optimization in
``_evolve_per_mode_pade`` (#327 follow-up).

Pins three correctness properties and one perf budget:

1. **Bit-exact equivalence on dense IC** — Gaussian IC has every
   mode active by construction, so the sparse-IC code path must
   produce byte-identical output to the legacy all-modes path.
2. **Equivalence within sigfig precision on plane-wave IC** — the
   only difference vs legacy is that inactive modes' contribution
   to ifft is exactly zero (sparse-IC) instead of ~1e-14·growth
   (legacy via FFT-floor noise); for non-tachyonic params this
   difference is below ifft normalization.
3. **Divergence detection preserved** — the eigenvalue-based pre-
   check must still raise ``SimulationDivergedError`` on the
   noise-amplification regime documented in
   ``docs/tex/stability_probe.tex`` §sec:ieee-floor-verification.
4. **Perf budget** — likelihood-evaluation median wall on the
   Gertsenshtein fixture under plane-wave IC stays within budget.

The kill-switch ``TIDAL_MODAL_SPARSE_IC=0`` recovers legacy
behavior bit-exactly; tests use it to compare paths.
"""

from __future__ import annotations

import os
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

from tidal.solver._exceptions import SimulationDivergedError

REPO = Path(__file__).resolve().parents[1]
GERT_SPEC = REPO / "examples/data/torsion_gertsenshtein_general_nonminimal.json"
DPP_SPEC = REPO / "examples/data/dark_photon_plasma.json"
D1_SPEC = REPO / "examples/data/torsion_gertsenshtein_nonminimal.json"


def _build_inference_args(  # internal test helper
    *,
    n_samples: int = 1,
    grid_shape: str = "256",
    bounds: str = "0:100",
    ic: str = "plane-wave",
    spec_path: Path = GERT_SPEC,
    fixed_params: tuple[str, ...] = (
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
    ),
    extra: tuple[str, ...] = (),
):
    """Construct an argparse Namespace for ``run_inference_step``."""
    from tidal.cli import _build_parser

    parser = _build_parser()
    args = [
        "sample",
        str(spec_path),
        "--prior",
        "delta1=uniform:-0.025:0.025",
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
        f"/tmp/test_modal_sparse_ic_{os.getpid()}",
        *extra,
    ]
    return parser.parse_args(args)


def _run_sim_snapshots(  # internal test helper, returns NDArray
    *,
    sparse_ic: bool,
    overrides: dict[str, float],
    spec_path: Path = GERT_SPEC,
    grid_shape: str = "256",
    bounds: str = "0:100",
    ic: str = "plane-wave",
    fixed_params: tuple[str, ...] = (
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
    ),
):
    """Run one in-memory simulation and return the snapshot array."""
    from tidal.cli._sweep import run_inference_step

    prev = os.environ.get("TIDAL_MODAL_SPARSE_IC")
    os.environ["TIDAL_MODAL_SPARSE_IC"] = "1" if sparse_ic else "0"
    try:
        ns = _build_inference_args(
            spec_path=spec_path,
            grid_shape=grid_shape,
            bounds=bounds,
            ic=ic,
            fixed_params=fixed_params,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim_data = run_inference_step(ns, spec_path, overrides)
    finally:
        if prev is None:
            os.environ.pop("TIDAL_MODAL_SPARSE_IC", None)
        else:
            os.environ["TIDAL_MODAL_SPARSE_IC"] = prev
    # Stack field arrays in deterministic name order so the two
    # paths' outputs are directly comparable element-by-element.
    return np.stack(
        [sim_data.fields[name] for name in sorted(sim_data.fields)],
        axis=0,
    )


# ---------------------------------------------------------------------
# 1.  Dense IC — bit-exact across paths
# ---------------------------------------------------------------------


@pytest.mark.slow
def test_dense_gaussian_ic_bit_exact() -> None:
    """For Gaussian IC every Fourier mode is active; sparse and legacy
    paths must produce byte-identical output (rtol=0, atol=0).
    """
    overrides = {"delta1": 0.0}
    snaps_sparse = _run_sim_snapshots(
        sparse_ic=True,
        overrides=overrides,
        spec_path=DPP_SPEC,
        grid_shape="32",
        bounds="0:50",
        ic="gaussian",
        fixed_params=(
            "kappa=1.0",
            "B0=0.01",
            "mA2=0.5",
            "deltam=0.0",
            "xi=0.5",
            "alpha3=0.1",
        ),
    )
    snaps_legacy = _run_sim_snapshots(
        sparse_ic=False,
        overrides=overrides,
        spec_path=DPP_SPEC,
        grid_shape="32",
        bounds="0:50",
        ic="gaussian",
        fixed_params=(
            "kappa=1.0",
            "B0=0.01",
            "mA2=0.5",
            "deltam=0.0",
            "xi=0.5",
            "alpha3=0.1",
        ),
    )
    np.testing.assert_array_equal(snaps_sparse, snaps_legacy)


# ---------------------------------------------------------------------
# 2.  Plane-wave IC — sigfig equivalence on non-tachyonic params
# ---------------------------------------------------------------------


@pytest.mark.slow
def test_plane_wave_ic_sigfig_equivalence() -> None:
    """For plane-wave IC at the campaign working point, sparse and
    legacy paths agree to sigfig precision: the only difference is
    the IEEE 754 FFT floor at non-fundamental bins (~ 1e-14)
    propagated through the matvec.  Tolerance set to ``atol=1e-12``
    to absorb this without false alarms.
    """
    overrides = {"delta1": 0.0}
    snaps_sparse = _run_sim_snapshots(
        sparse_ic=True,
        overrides=overrides,
    )
    snaps_legacy = _run_sim_snapshots(
        sparse_ic=False,
        overrides=overrides,
    )
    # Shape matches.
    assert snaps_sparse.shape == snaps_legacy.shape
    # Numerical agreement to sigfig precision.
    np.testing.assert_allclose(snaps_sparse, snaps_legacy, atol=1e-12, rtol=1e-9)


# ---------------------------------------------------------------------
# 3.  Divergence detection — eigenvalue pre-check still trips
# ---------------------------------------------------------------------


@pytest.mark.slow
def test_eigenvalue_pre_check_catches_noise_amplified_divergence() -> None:
    """Sample 2860 from the Stage C truth table (γ_v1 = 17.11) is
    the documented IEEE-floor-amplification regime
    (``docs/tex/stability_probe.tex`` §sec:ieee-floor-verification).
    With γ ≫ log(100)/t_end ≈ 0.46, the noise floor at non-IC
    bins amplifies through ``exp(γ·t)`` to non-perturbative
    amplitude well before t_end.  The sparse-IC eigenvalue pre-check
    must raise ``SimulationDivergedError`` before the time loop —
    this is the load-bearing safety property the new pre-check
    inherits from the all-modes expm pre-check.
    """
    fixed = (
        "kappa=1.0",
        "B0=0.01",
        "alpha1=-0.354871297172597",
        "alpha2=-1.44478374604432",
        "alpha3=0.0794345149654468",
    )
    overrides = {"delta1": 0.077796221228164}
    with pytest.raises(SimulationDivergedError):
        _run_sim_snapshots(
            sparse_ic=True,
            overrides=overrides,
            spec_path=D1_SPEC,
            grid_shape="256",
            bounds="0:100",
            ic="plane-wave",
            fixed_params=fixed,
        )


# ---------------------------------------------------------------------
# 4.  Perf budget
# ---------------------------------------------------------------------


@pytest.mark.slow
def test_inference_eval_perf() -> None:
    """Likelihood-eval median wall on the Gertsenshtein fixture
    under sparse-IC plane-wave IC stays under budget.

    The budget is set conservatively to absorb CI noise: the
    optimization's measured win on this fixture is ~5 % wall
    (most of the per-call cost is matrix construction +
    coefficient resolution, which sparse-IC does not address).
    A 200-ms budget here is ~3× the typical 75 ms median, leaving
    plenty of headroom for slow CI workers without flakiness.
    """
    from tidal.cli._sweep import run_inference_step

    os.environ["TIDAL_MODAL_SPARSE_IC"] = "1"
    overrides = {"delta1": 0.0}
    ns = _build_inference_args()

    # Warmup
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

    median_ms = float(np.median(times))
    assert median_ms <= 200.0, (
        f"likelihood-eval median wall {median_ms:.2f} ms exceeds 200 ms budget; "
        "investigate matrix-construction or coefficient-resolution paths"
    )
