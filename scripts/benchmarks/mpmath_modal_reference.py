r"""mpmath reference oracle for the GH #367 modal-solver discretization.

Question this script answers: is the Fourier-spectral discretization of the
pos-dep + periodic problem fundamentally broken, or is it correct in
principle and only float64 precision is the issue?

The diagnostic:

- exp(A·t)·y0 computed by scipy.sparse.linalg.expm_multiply at float64
- exp(A·t)·y0 computed by scipy.linalg.expm at float64 (dense Padé)
- exp(A·t)·y0 computed by mpmath.expm at dps=50 (arbitrary-precision oracle)
- exp(A·t)·y0 computed by solve_cvode in physical space (truth reference)

Interpretation:

- mpmath ≈ CVODE  → the spectral formulation is correct in principle; the
  float64 issue is precision/conditioning. Hou-Li filter / better
  discretization should work.
- mpmath ≈ scipy (still wildly wrong) → the discretization itself is
  fundamentally broken at finite N; we'd need mixed real/Fourier (path 4).

N=32 is used so the matrix (n_total≈198) is tractable in mpmath at dps=50
in a few minutes. The spurious eigenvalue at N=32 is k_max=1.005, so at
t_end=20 the artifact factor exp(k_max·t) ≈ exp(20) ≈ 5e8 — still firmly
in the broken regime, so the diagnostic question is preserved.

Usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
        uv run python scripts/benchmarks/mpmath_modal_reference.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.sparse.linalg

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.coefficients import CoefficientEvaluator  # noqa: E402
from tidal.solver.cvode import solve_cvode  # noqa: E402
from tidal.solver.grid import GridInfo  # noqa: E402
from tidal.solver.modal import (  # noqa: E402
    _build_convolution_matrix,
    _build_k_axes,
    _build_k_grid,
    _fft_slots,
    _ifft_slots,
)
from tidal.solver.state import StateLayout  # noqa: E402
from tidal.symbolic.json_loader import EquationSystem  # noqa: E402

THEORY_JSON = REPO_ROOT / "examples/data/gertsenshtein_e0_dual_gaussian.json"
PARAMS: dict[str, float] = {
    "kappa": 1.0,
    "Bpeak": 0.01,
    "sigB": 5.0,
    "zc1": 25.0,
    "zc2": 75.0,
}
IC_AMPLITUDE = 1e-2
IC_WIDTH = 5.0
IC_CENTER = 25.0
DOMAIN_LENGTH = 100.0
IC_SLOT = "h_5"

N = 32  # small enough for mpmath; large enough to expose the artifact
T_END = 20.0  # firmly in the broken regime (k_max*t = 20.1)
MPMATH_DPS = 50  # ~50 decimal digits of precision
DIAGNOSTIC_THRESHOLD = (
    0.01  # 1% tolerance for the mpmath-vs-cvode / mpmath-vs-scipy classifier
)


def build_y0_physical(layout: StateLayout, n_points: int) -> np.ndarray:
    y0 = np.zeros(layout.total_size, dtype=float)
    x_axis = np.linspace(0.0, DOMAIN_LENGTH, n_points, endpoint=False)
    gauss = IC_AMPLITUDE * np.exp(-((x_axis - IC_CENTER) ** 2) / (2 * IC_WIDTH**2))
    slot_idx = layout.slot_name_to_idx[IC_SLOT]
    y0[layout.slot_slice(slot_idx)] = gauss
    return y0


def h5_peak_from_flat(y_flat: np.ndarray, layout: StateLayout) -> float:
    h5_slice = layout.slot_slice(layout.slot_name_to_idx[IC_SLOT])
    return float(np.max(np.abs(y_flat[h5_slice])))


def main() -> None:
    print(f"Loading {THEORY_JSON.relative_to(REPO_ROOT)}")
    with Path(THEORY_JSON).open(encoding="utf-8") as f:
        spec = EquationSystem.from_dict(json.load(f))

    grid = GridInfo(shape=(N,), bounds=[(0.0, DOMAIN_LENGTH)], periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    y0 = build_y0_physical(layout, N)
    print(f"  N={N}, t_end={T_END}, n_slots={layout.num_slots}")
    print(f"  initial h5_peak = {h5_peak_from_flat(y0, layout):.6e}")

    # --- Physical-space truth via CVODE -------------------------------------
    print("\n[1/4] solve_cvode (physical-space truth)")
    t0 = time.perf_counter()
    cvode_result = solve_cvode(
        spec, grid, y0, t_span=(0.0, T_END), parameters=PARAMS, num_snapshots=2
    )
    cvode_time = time.perf_counter() - t0
    cvode_final = np.asarray(cvode_result["y"][-1])
    cvode_h5 = h5_peak_from_flat(cvode_final, layout)
    print(f"  cvode h5_peak = {cvode_h5:.6e}   wall = {cvode_time:.2f}s")

    # --- Build the convolution matrix A_full ---------------------------------
    print("\n[2/4] _build_convolution_matrix")
    rfft_shape = (N // 2 + 1,)
    n_modes = int(np.prod(rfft_shape))
    print(f"  rfft_shape = {rfft_shape}, n_modes = {n_modes}")

    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    coeff_eval = CoefficientEvaluator(spec, grid, PARAMS)
    t0 = time.perf_counter()
    A_full = _build_convolution_matrix(
        spec, layout, grid, coeff_eval, k_grid, rfft_shape
    )
    build_time = time.perf_counter() - t0
    n_total = A_full.shape[0]
    print(f"  A_full shape = {A_full.shape}   build wall = {build_time:.2f}s")
    eigs = np.linalg.eigvals(A_full)
    max_re = float(np.max(eigs.real))
    print(f"  max Re(λ) = {max_re:.4f}   (matches handoff table: N=32 → 1.005)")

    # FFT y0 to k-space, flatten
    y0_hat = _fft_slots(y0, layout, grid)
    y0_flat = y0_hat.ravel()
    print(f"  y0_flat shape = {y0_flat.shape}   norm = {np.linalg.norm(y0_flat):.4e}")

    # --- float64 references --------------------------------------------------
    print("\n[3/4] float64 matrix-exponential references")
    t0 = time.perf_counter()
    y_expm_mult = np.asarray(
        scipy.sparse.linalg.expm_multiply(
            A_full, y0_flat, start=0.0, stop=T_END, num=2
        )[-1],
        dtype=np.complex128,
    )
    em_time = time.perf_counter() - t0
    y_expm_mult_phys = _ifft_slots(
        y_expm_mult.reshape(layout.num_slots, n_modes), layout, grid
    )
    em_h5 = h5_peak_from_flat(y_expm_mult_phys, layout)
    print(f"  expm_multiply h5_peak = {em_h5:.6e}   wall = {em_time:.2f}s")

    t0 = time.perf_counter()
    y_dense = scipy.linalg.expm(A_full * T_END) @ y0_flat
    dp_time = time.perf_counter() - t0
    y_dense_phys = _ifft_slots(y_dense.reshape(layout.num_slots, n_modes), layout, grid)
    dp_h5 = h5_peak_from_flat(y_dense_phys, layout)
    print(f"  scipy.linalg.expm h5_peak = {dp_h5:.6e}   wall = {dp_time:.2f}s")

    # --- mpmath oracle -------------------------------------------------------
    print(f"\n[4/4] mpmath.expm oracle (dps={MPMATH_DPS}, may take minutes)")
    import mpmath

    mpmath.mp.dps = MPMATH_DPS
    t0 = time.perf_counter()
    A_mp = mpmath.matrix(n_total, n_total)
    for i in range(n_total):
        for j in range(n_total):
            v = A_full[i, j]
            A_mp[i, j] = mpmath.mpc(v.real, v.imag)
    convert_time = time.perf_counter() - t0
    print(f"  converted A → mpmath in {convert_time:.1f}s")

    t0 = time.perf_counter()
    expAt_mp = mpmath.expm(A_mp * mpmath.mpf(T_END))
    expm_time = time.perf_counter() - t0
    print(f"  mpmath.expm(A·t) in {expm_time:.1f}s")

    t0 = time.perf_counter()
    y0_mp = mpmath.matrix(n_total, 1)
    for i in range(n_total):
        v = y0_flat[i]
        y0_mp[i, 0] = mpmath.mpc(v.real, v.imag)
    y_mp = expAt_mp * y0_mp
    matvec_time = time.perf_counter() - t0
    print(f"  matvec in {matvec_time:.1f}s")

    y_mp_np = np.zeros(n_total, dtype=np.complex128)
    for i in range(n_total):
        y_mp_np[i] = complex(y_mp[i, 0])
    y_mp_phys = _ifft_slots(y_mp_np.reshape(layout.num_slots, n_modes), layout, grid)
    mp_h5 = h5_peak_from_flat(y_mp_phys, layout)
    print(f"  mpmath h5_peak = {mp_h5:.6e}")

    # --- Summary -------------------------------------------------------------
    print("\n=== Summary ===")
    print(f"  CVODE (physical truth)  h5_peak = {cvode_h5:.6e}")
    print(f"  scipy.expm_multiply     h5_peak = {em_h5:.6e}")
    print(f"  scipy.linalg.expm       h5_peak = {dp_h5:.6e}")
    print(f"  mpmath.expm (dps={MPMATH_DPS}) h5_peak = {mp_h5:.6e}")
    print()
    print(f"  mpmath / CVODE ratio    = {mp_h5 / cvode_h5:.3e}")
    print(f"  scipy / CVODE ratio     = {em_h5 / cvode_h5:.3e}")
    print()
    if abs(mp_h5 - cvode_h5) / cvode_h5 < DIAGNOSTIC_THRESHOLD:
        print(
            "DIAGNOSTIC RESULT: mpmath ≈ CVODE → spectral formulation is correct "
            "in principle. The float64 issue is precision/conditioning. "
            "Filter / better representation should work."
        )
    elif abs(mp_h5 - em_h5) / max(em_h5, 1e-300) < DIAGNOSTIC_THRESHOLD:
        print(
            "DIAGNOSTIC RESULT: mpmath ≈ scipy (both wildly wrong) → the "
            "discretization itself is fundamentally broken at finite N. "
            "Hou-Li filter alone cannot fix it; mixed real/Fourier needed."
        )
    else:
        print(
            "DIAGNOSTIC RESULT: mpmath differs from both — unexpected. "
            "Inspect the per-mode breakdown."
        )

    # --- Save fixture --------------------------------------------------------
    fixture_dir = REPO_ROOT / "tests/fixtures"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    fixture_path = fixture_dir / "gh367_mpmath_h5_n32.npz"
    np.savez(
        fixture_path,
        N=N,
        t_end=T_END,
        mpmath_dps=MPMATH_DPS,
        y0_flat=y0_flat,
        y_mpmath_kspace=y_mp_np,
        y_mpmath_physical=y_mp_phys,
        cvode_h5_peak=cvode_h5,
        mpmath_h5_peak=mp_h5,
        scipy_expm_multiply_h5_peak=em_h5,
        scipy_dense_expm_h5_peak=dp_h5,
        max_real_eigenvalue=max_re,
    )
    print(f"\nWrote {fixture_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
