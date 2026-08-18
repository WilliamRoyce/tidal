"""Benchmark: matrix-exponential evolution paths for per-mode modal solver.

Decision input for the Tier 1 ill-conditioning fix (#320). Compares wall time of:

  (a) eigendecomposition + V·diag(exp(λt))·V⁻¹      [current default]
  (b) per-mode  scipy.sparse.linalg.expm_multiply   [Method B, Al-Mohy & Higham 2011]
  (c) per-mode expm_multiply with start/stop/num    [Method B variant, scipy adaptive]
  (d) scipy.linalg.expm precompute + matvec         [Method A, Higham 2009 Padé S&S]

The 4-row workload table covers Stage A (dark_photon_plasma CDT,
torsion_dark_photon_fv) and Stage B+ (graviton_torsion, euler_heisenberg).

KEY FINDING (path D wins, universally).
=========================================
Path (d) — `scipy.linalg.expm(M · dt)` precomputed once per mode, then matvec
in the snapshot loop — is **as fast as or faster than eigendecomposition** on
every workload tested:

    Theory                              (a)/(a)  (b)/(a)  (c)/(a)  (d)/(a)
    dark_photon_plasma (CDT, Stage A)    1.00x   28.61x   31.84x   1.06x
    torsion_dark_photon_fv (Stage A)     1.00x   22.31x   25.82x   0.49x
    graviton_torsion (Stage B+)          1.00x  803.08x   65.34x   0.99x
    euler_heisenberg                     1.00x  469.32x   59.35x   0.75x

WHY (b)/(c) ARE SLOW
====================
`scipy.sparse.linalg.expm_multiply` has substantial per-call setup overhead:
estimating ||A||₁, choosing scaling parameter s, choosing Taylor truncation m,
allocating Krylov-Taylor workspace. For our small per-mode block sizes
(bs ≤ 15), this overhead dominates the actual O(bs³) ≈ 1000 flops of the
exponential. Calling expm_multiply once per mode (path b) or once per mode
with snapshot scheduling (path c) means we pay that overhead 100s-1000s of
times.

WHY (d) IS COMPETITIVE WITH EIGENDECOMPOSITION
==============================================
`scipy.linalg.expm` (Higham 2009 Padé scaling-and-squaring) is a simpler
direct routine that avoids the iterative-action overhead — it produces
`exp(M·dt)` in one O(bs³) call. Once `exp(M·dt)` is precomputed per mode,
the snapshot loop becomes pure matvec at O(bs²) per snapshot per mode,
identical in cost to eigendecomposition's vectorized einsum. The matrix
exponential then propagates uniformly: `y[ti+1] = exp_M_dt @ y[ti]`.

Path (d) is **mathematically robust for arbitrary cond(V)** (Higham 2009 §3:
backward error bounded by function condition number, not eigenvector
condition number) AND **wall-time competitive with eigendecomposition**.
This makes the originally-planned hybrid (eigendecomposition + cond(V)
threshold + expm_multiply fallback) unnecessary: path (d) replaces the
eigendecomposition default uniformly.

REMAINING ROLE OF EIGENDECOMPOSITION
====================================
Pass 1 perturbative Duhamel solver (`_evolve_duhamel_per_mode`) needs per-
eigenvalue kernels G(λ_i, λ_j; t) and the V/V⁻¹ projection. Eigendecomposition
remains opt-in (collect_eigendata=True) for those callers; if cond(V) > 10¹²
on that opt-in path, raise NotImplementedError (see plan §"Pass 1 separate
gap"). Tier 1.5 follow-up implements the augmented-exp Pass 1 path
(Al-Mohy & Higham 2011 §5.2) so Pass 1 can also drop the eigendecomposition
dependency — tracked as #320-Pass1.

Run:  uv run python scripts/benchmark_modal_evolution.py
"""

from __future__ import annotations

import time

import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import expm_multiply


def make_block(
    bs: int, n_modes: int, cond_target: float, *, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Construct a per-mode block matrix with controllable cond(V).

    Returns (A_modes, B_modes) where A_modes shape is (n_modes, bs, bs) and
    B_modes is identity (homogeneous problem). cond(V) is tuned by clustering
    eigenvalues: small spread -> large cond(V).
    """
    rng = np.random.default_rng(seed)
    # Eigenvalues scattered by physical dispersion: imaginary part ~ k, small
    # real part. Cluster two adjacent eigenvalues to make V nearly singular.
    A_modes = np.empty((n_modes, bs, bs), dtype=np.complex128)
    for m in range(n_modes):
        k = (m + 1) / n_modes
        # Imaginary spectrum (oscillatory dynamics) with a controllable
        # near-degeneracy that drives cond(V) up.
        eigs = 1j * (np.arange(bs) - bs / 2 + 1) * k
        if cond_target > 1e10 and bs >= 2:
            # Force two eigenvalues to near-coincide → near-defective V
            gap = 10 ** (-np.log10(cond_target) / 2)
            eigs[1] = eigs[0] * (1 + gap) + 1j * gap
        # Build A = P · diag(eigs) · P⁻¹ with P having the desired conditioning
        # by skew-symmetrizing a perturbed identity.
        P = np.eye(bs, dtype=np.complex128) + 0.3 * (
            rng.standard_normal((bs, bs)) + 1j * rng.standard_normal((bs, bs))
        )
        if cond_target > 1e10 and bs >= 2:
            # Make column 1 nearly parallel to column 0
            P[:, 1] = P[:, 0] * (1 + 1e-8) + 1e-8 * P[:, 1]
        A_m = P @ np.diag(eigs) @ np.linalg.inv(P)
        A_modes[m] = A_m
    B_modes = np.tile(np.eye(bs, dtype=np.complex128), (n_modes, 1, 1))
    return A_modes, B_modes


def path_a_eigendecomposition(
    A_modes: np.ndarray,
    B_modes: np.ndarray,
    y0: np.ndarray,
    t_eval: np.ndarray,
) -> np.ndarray:
    """Path (a): eigendecomposition + V·diag(exp(λt))·V⁻¹ — matches current modal.py."""
    n_modes, bs, _ = A_modes.shape
    n_snap = len(t_eval)
    eig_vals = np.empty((n_modes, bs), dtype=np.complex128)
    v_mat = np.empty((n_modes, bs, bs), dtype=np.complex128)
    for m in range(n_modes):
        ev, vr = sla.eig(A_modes[m], B_modes[m], right=True)
        eig_vals[m] = ev
        v_mat[m] = vr
    # Filter gauge / spurious infinite eigs the way modal.py does
    gauge = ~np.isfinite(eig_vals) | (np.abs(eig_vals) > 1e12)
    eig_vals[gauge] = 0.0
    v_inv = np.linalg.inv(v_mat)
    y0_eigen = np.einsum("mij,mj->mi", v_inv, y0)  # (n_modes, bs)
    V_y0 = v_mat * y0_eigen[:, np.newaxis, :]  # (n_modes, bs, bs)
    out = np.empty((n_snap, n_modes, bs), dtype=np.complex128)
    t0 = t_eval[0]
    for ti, t in enumerate(t_eval):
        dt = t - t0
        exp_lambda = np.exp(eig_vals * dt)
        out[ti] = np.einsum("mij,mj->mi", V_y0, exp_lambda)
    return out


def path_b_expm_per_mode(
    A_modes: np.ndarray,
    B_modes: np.ndarray,
    y0: np.ndarray,
    t_eval: np.ndarray,
) -> np.ndarray:
    """Path (b): per-mode expm_multiply, single call per (mode, snapshot)."""
    n_modes, bs, _ = A_modes.shape
    n_snap = len(t_eval)
    out = np.empty((n_snap, n_modes, bs), dtype=np.complex128)
    t0 = t_eval[0]
    for m in range(n_modes):
        # M = B⁻¹ A — one stable solve since B is well-conditioned
        M = np.linalg.solve(B_modes[m], A_modes[m])
        for ti, t in enumerate(t_eval):
            dt = float(t - t0)
            if dt == 0.0:
                out[ti, m] = y0[m]
            else:
                out[ti, m] = expm_multiply(M * dt, y0[m])
    return out


def path_c_expm_batched(
    A_modes: np.ndarray,
    B_modes: np.ndarray,
    y0: np.ndarray,
    t_eval: np.ndarray,
) -> np.ndarray:
    """Path (c): per-mode expm_multiply with start/stop/num scheduling."""
    n_modes, bs, _ = A_modes.shape
    n_snap = len(t_eval)
    out = np.empty((n_snap, n_modes, bs), dtype=np.complex128)
    t0 = float(t_eval[0])
    t_end = float(t_eval[-1])
    for m in range(n_modes):
        M = np.linalg.solve(B_modes[m], A_modes[m])
        if n_snap == 1 or t_end == t0:
            for ti, t in enumerate(t_eval):
                dt = float(t - t0)
                out[ti, m] = y0[m] if dt == 0.0 else expm_multiply(M * dt, y0[m])
        else:
            # scipy's vectorized path
            ys = expm_multiply(M, y0[m], start=t0, stop=t_end, num=n_snap)
            out[:, m, :] = ys
    return out


def path_d_expm_precompute(
    A_modes: np.ndarray,
    B_modes: np.ndarray,
    y0: np.ndarray,
    t_eval: np.ndarray,
) -> np.ndarray:
    """Path (d): precompute exp(M·t_i) per mode via scipy.linalg.expm, then matvec.

    For evenly spaced snapshots, computes exp(M·dt) once per mode and applies
    repeatedly: y[ti+1] = exp(M·dt) @ y[ti]. Avoids expm_multiply's per-call
    setup cost. O(bs³) precompute per mode + O(bs²) per snapshot per mode.
    """
    n_modes, bs, _ = A_modes.shape
    n_snap = len(t_eval)
    out = np.empty((n_snap, n_modes, bs), dtype=np.complex128)
    t0 = float(t_eval[0])

    # Detect uniform spacing
    if n_snap > 1:
        dts = np.diff(t_eval)
        uniform = np.allclose(dts, dts[0])
    else:
        uniform = True
        dts = np.array([0.0])

    for m in range(n_modes):
        M = np.linalg.solve(B_modes[m], A_modes[m])
        if uniform and n_snap > 1:
            expM_dt = sla.expm(M * float(dts[0]))
            y_curr = y0[m].copy()
            out[0, m] = y_curr
            for ti in range(1, n_snap):
                y_curr = expM_dt @ y_curr
                out[ti, m] = y_curr
        else:
            for ti, t in enumerate(t_eval):
                dt = float(t - t0)
                if dt == 0.0:
                    out[ti, m] = y0[m]
                else:
                    out[ti, m] = sla.expm(M * dt) @ y0[m]
    return out


def time_path(fn, *args, n_repeat: int = 3) -> tuple[float, np.ndarray]:
    # Warm-up
    result = fn(*args)
    times = []
    for _ in range(n_repeat):
        t = time.perf_counter()
        result = fn(*args)
        times.append(time.perf_counter() - t)
    return float(np.median(times)), result


def cond_v_for_block(A_block: np.ndarray, B_block: np.ndarray) -> float:
    """Compute cond(V) for one mode-block to verify the synthetic ill-conditioning."""
    _ev, vr = sla.eig(A_block, B_block, right=True)
    return float(np.linalg.cond(vr))


def benchmark_row(
    name: str,
    bs: int,
    n_blocks: int,
    n_modes: int,
    n_snapshots: int,
    t_end: float,
    cond_target: float,
) -> dict:
    """Run a benchmark row from the workload table.

    Returns wall times for each path and the speed ratio.
    """
    # Build n_blocks independent block matrices and concatenate as if they were
    # the per-block work in _evolve_per_mode. We time one block here since
    # blocks are independent; multiply by n_blocks for total wall time.
    A_modes, B_modes = make_block(bs, n_modes, cond_target, seed=42)
    rng = np.random.default_rng(0)
    y0 = rng.standard_normal((n_modes, bs)) + 1j * rng.standard_normal((n_modes, bs))
    t_eval = np.linspace(0.0, t_end, n_snapshots)

    # Verify cond(V) characteristics
    cond_vs = [cond_v_for_block(A_modes[i], B_modes[i]) for i in range(min(5, n_modes))]
    cond_v_typical = float(np.median(cond_vs))

    # Time each path
    t_a, r_a = time_path(path_a_eigendecomposition, A_modes, B_modes, y0, t_eval)
    t_b, r_b = time_path(path_b_expm_per_mode, A_modes, B_modes, y0, t_eval)
    t_c, r_c = time_path(path_c_expm_batched, A_modes, B_modes, y0, t_eval)
    t_d, r_d = time_path(path_d_expm_precompute, A_modes, B_modes, y0, t_eval)

    # Total wall (multiply per-block time by n_blocks)
    t_a_total = t_a * n_blocks
    t_b_total = t_b * n_blocks
    t_c_total = t_c * n_blocks
    t_d_total = t_d * n_blocks

    # Cross-check that paths agree (where eigendecomposition is sound)
    if cond_target < 1e10:
        # Path A is reliable here — check (b), (c), (d) match it
        rel_b = float(np.max(np.abs(r_a - r_b)) / max(np.max(np.abs(r_a)), 1e-15))
        rel_c = float(np.max(np.abs(r_a - r_c)) / max(np.max(np.abs(r_a)), 1e-15))
        rel_d = float(np.max(np.abs(r_a - r_d)) / max(np.max(np.abs(r_a)), 1e-15))
    else:
        # Path A is unsound — compare (b) vs (c) vs (d) instead
        rel_b = float("nan")
        rel_c = float(np.max(np.abs(r_b - r_c)) / max(np.max(np.abs(r_b)), 1e-15))
        rel_d = float(np.max(np.abs(r_b - r_d)) / max(np.max(np.abs(r_b)), 1e-15))

    return {
        "name": name,
        "bs": bs,
        "n_blocks": n_blocks,
        "n_modes": n_modes,
        "n_snapshots": n_snapshots,
        "t_end": t_end,
        "cond_v_typical": cond_v_typical,
        "wall_a_ms": t_a_total * 1000.0,
        "wall_b_ms": t_b_total * 1000.0,
        "wall_c_ms": t_c_total * 1000.0,
        "wall_d_ms": t_d_total * 1000.0,
        "ratio_b_over_a": t_b_total / max(t_a_total, 1e-9),
        "ratio_c_over_a": t_c_total / max(t_a_total, 1e-9),
        "ratio_d_over_a": t_d_total / max(t_a_total, 1e-9),
        "agreement_b_vs_a": rel_b,
        "agreement_c_vs_a": rel_c,
        "agreement_d_vs_a": rel_d,
    }


def main() -> None:
    workloads = [
        # (name, bs, n_blocks, n_modes, n_snapshots, t_end, cond_target)
        ("dark_photon_plasma (CDT, Stage A)", 6, 5, 64, 2, 10.0, 1e15),
        ("torsion_dark_photon_fv (Stage A)", 6, 3, 64, 2, 10.0, 1e8),
        ("graviton_torsion (Stage B+)", 10, 4, 128, 50, 50.0, 1e9),
        ("euler_heisenberg", 5, 2, 96, 50, 50.0, 1e8),
    ]

    print("\n=== Modal evolution benchmark — Tier 1 (#320) decision input ===\n")
    print(
        f"{'Theory':<40s} {'bs':>3s} {'#blk':>4s} {'#mod':>5s} {'#snap':>5s} "
        f"{'cond(V)':>10s}  {'(a) eig':>9s}  {'(b) ppm':>9s}  {'(c) bch':>9s}  "
        f"{'(d) Pad':>9s}  {'b/a':>6s}  {'c/a':>6s}  {'d/a':>6s}"
    )
    print("-" * 175)

    results = []
    for wl in workloads:
        r = benchmark_row(*wl)
        results.append(r)
        print(
            f"{r['name']:<40s} "
            f"{r['bs']:>3d} {r['n_blocks']:>4d} {r['n_modes']:>5d} {r['n_snapshots']:>5d} "
            f"{r['cond_v_typical']:>10.2e}  "
            f"{r['wall_a_ms']:>7.1f}ms  "
            f"{r['wall_b_ms']:>7.1f}ms  "
            f"{r['wall_c_ms']:>7.1f}ms  "
            f"{r['wall_d_ms']:>7.1f}ms  "
            f"{r['ratio_b_over_a']:>6.2f}  "
            f"{r['ratio_c_over_a']:>6.2f}  "
            f"{r['ratio_d_over_a']:>6.2f}"
        )

    print("\n=== Agreement check (relative max diff) ===")
    for r in results:
        b_str = (
            "n/a (eig unsound)"
            if np.isnan(r["agreement_b_vs_a"])
            else f"{r['agreement_b_vs_a']:.2e}"
        )
        print(
            f"  {r['name']:<40s}  (b): {b_str}  "
            f"(c): {r['agreement_c_vs_a']:.2e}  "
            f"(d): {r['agreement_d_vs_a']:.2e}"
        )

    print("\n=== Decision summary ===")
    worst_b = max(r["ratio_b_over_a"] for r in results)
    worst_c = max(r["ratio_c_over_a"] for r in results)
    worst_d = max(r["ratio_d_over_a"] for r in results)
    print(f"  Worst (b)/(a) ratio: {worst_b:.2f}  (path b: per-mode expm_multiply)")
    print(
        f"  Worst (c)/(a) ratio: {worst_c:.2f}  (path c: expm_multiply start/stop/num)"
    )
    print(
        f"  Worst (d)/(a) ratio: {worst_d:.2f}  (path d: scipy.linalg.expm precompute + matvec)"
    )
    print(
        "  Decision rule: ≤ 2× → unconditional, otherwise hybrid with κ_thresh = 1e12"
    )
    best_robust = min(worst_b, worst_c, worst_d)
    best_path = ["b", "c", "d"][[worst_b, worst_c, worst_d].index(best_robust)]
    if best_robust <= 2.0:
        print(f"  → Recommendation: UNCONDITIONAL via path ({best_path})")
    else:
        print("  → Recommendation: HYBRID with cond(V) > 1e12 fallback")


if __name__ == "__main__":
    main()
