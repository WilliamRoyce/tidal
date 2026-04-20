#!/usr/bin/env python3
"""Schur complement analysis of torsion-mediated Gertsenshtein amplification.

Extracts the effective coupling mu_eff and mass m2_eff for the h5<->a1 channel
after eliminating all torsion constraint fields via the Fourier-space Schur
complement. Computes the amplification factor A = |mu_eff / mu_GR|^2 and
maps the instability boundary and suppression valley.

Physics: The nonminimal coupling delta1 * R_tilde_{[mu nu]} F^{mu nu} creates
a torsion-mediated feedback loop. Algebraic torsion constraints are eliminated
via the Schur complement, modifying both the effective h5<->a1 coupling (mu_eff)
and the effective photon mass (m2_eff). The coupling reverses sign at a critical
delta1_crit, producing a suppression valley (zero crossing) and amplification
(|mu_eff| > |mu_GR|) for delta1 > delta1_crit.

References
----------
  - Plan: docs/AMPLIFICATION_INVESTIGATION.md
  - Physics: docs/tex/amplification_mechanism.tex
  - Sweep data: examples/data/nonminimal_heatmap_d1_a2_hires/

Usage:
  python schur_complement_analysis.py [--spec JSON] [--output DIR]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import operator

from tidal.solver import _build_evolution_matrices
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import (
    load_equation_system,
    normalize_kinetic_coefficients,
)


def build_grid(N: int = 256, L: float = 100.0) -> GridInfo:
    """Standard 1D periodic grid matching sweep parameters."""
    return GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))


def get_k_grid(grid: GridInfo) -> tuple[list[np.ndarray], tuple[int, ...]]:
    """Build Fourier wavenumber grid."""
    N = grid.shape[0]
    dx = grid.dx[0]
    k_vals = 2 * np.pi * np.fft.rfftfreq(N, dx)
    return [k_vals], (N // 2 + 1,)


def extract_h5_a1_block(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    params: dict[str, float],
) -> np.ndarray:
    """Extract the 4x4 {h5, vh5, a1, va1} block from the Schur complement.

    Returns
    -------
        Block matrix of shape (n_modes, 4, 4) with rows/cols [h5, vh5, a1, va1].
    """
    coeff_eval = CoefficientEvaluator(spec, grid, params)
    A_reduced, _, _, _, _, mapping = _build_evolution_matrices(
        spec, layout, grid, coeff_eval, k_grid, rfft_shape
    )

    # Map original slots to reduced slots
    h5_red = mapping[layout.field_slot_map["h_5"]]
    vh5_red = mapping[layout.velocity_slot_map["h_5"]]
    a1_red = mapping[layout.field_slot_map["a_1"]]
    va1_red = mapping[layout.velocity_slot_map["a_1"]]

    idx = np.array([h5_red, vh5_red, a1_red, va1_red])
    return A_reduced[:, idx[:, None], idx[None, :]]


def analyze_point(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    base_params: dict[str, float],
    delta1: float,
    alpha2: float,
    k_idx: int = 8,
) -> dict:
    """Analyze the Schur complement at a single (delta1, alpha2) point.

    Returns dict with:
        mu_GR: baseline h5<->a1 coupling (imaginary part)
        mu_eff: effective coupling with torsion feedback
        mu_ratio: mu_eff / mu_GR (amplification of coupling)
        m2_a1_GR: baseline a1 effective mass^2
        m2_a1_eff: effective a1 mass^2 with torsion feedback
        m2_shift: mass shift from torsion feedback
        max_re_eig: maximum real eigenvalue of the 4x4 block
        stable: whether max_re_eig < growth threshold
        A_coupling: |mu_ratio|^2 (amplification factor from coupling alone)
    """
    k_vals = k_grid[0]

    # Baseline (delta1=0)
    params_gr = {**base_params, "delta1": 0.0, "alpha2": alpha2}
    B_gr = extract_h5_a1_block(spec, layout, grid, k_grid, rfft_shape, params_gr)

    # With torsion coupling
    params = {**base_params, "delta1": delta1, "alpha2": alpha2}
    B_eff = extract_h5_a1_block(spec, layout, grid, k_grid, rfft_shape, params)

    # Extract coupling: B[1,2] = v_h5 <- a1 (off-diagonal)
    mu_gr = B_gr[k_idx, 1, 2].imag
    mu_eff = B_eff[k_idx, 1, 2].imag

    # Extract mass: B[3,2] = v_a1 <- a1 (diagonal, real part)
    m2_a1_gr = B_gr[k_idx, 3, 2].real
    m2_a1_eff = B_eff[k_idx, 3, 2].real

    # Eigenvalues of full 4x4 block (all modes)
    all_eigs = np.linalg.eigvals(B_eff)
    all_eigs_gr = np.linalg.eigvals(B_gr)
    max_re = np.max(np.real(all_eigs), axis=1)
    max_re_gr = np.max(np.real(all_eigs_gr), axis=1)

    # Excess growth: how much faster do torsion-modified modes grow
    # compared to the GR baseline? The baseline already has Re(eig)=k
    # from wave propagation, so only the EXCESS matters.
    excess = max_re - max_re_gr

    # Growth threshold: the modal solver diverges when
    # exp(Re(lambda)*t_end) > 1e8, i.e., Re(lambda) > 0.37
    # But we want the excess over GR baseline
    t_end = base_params.get("t_end", 50.0)
    growth_threshold = np.log(1e8) / t_end  # ~0.37

    # Find most unstable mode (by excess over baseline)
    most_unstable_k = np.argmax(excess)
    max_re_eig = max_re[most_unstable_k]
    max_excess = float(excess[most_unstable_k])

    mu_ratio = mu_eff / mu_gr if abs(mu_gr) > 1e-30 else float("nan")

    return {
        "delta1": delta1,
        "alpha2": alpha2,
        "k": k_vals[k_idx],
        "mu_GR": mu_gr,
        "mu_eff": mu_eff,
        "mu_ratio": mu_ratio,
        "m2_a1_GR": m2_a1_gr,
        "m2_a1_eff": m2_a1_eff,
        "m2_shift": m2_a1_eff - m2_a1_gr,
        "max_re_eig": max_re_eig,
        "max_excess": max_excess,
        "most_unstable_k": k_vals[most_unstable_k],
        "stable": bool(max_excess < growth_threshold),
        "A_coupling": mu_ratio**2,
    }


def sweep_coupling(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    base_params: dict[str, float],
    delta1_values: np.ndarray,
    alpha2_values: np.ndarray,
    k_idx: int = 8,
) -> list[dict]:
    """Sweep (delta1, alpha2) parameter space."""
    results = []
    total = len(delta1_values) * len(alpha2_values)
    for i, d1 in enumerate(delta1_values):
        for j, a2 in enumerate(alpha2_values):
            n = i * len(alpha2_values) + j + 1
            if n % 50 == 0 or n == total:
                print(f"  Point {n}/{total}: delta1={d1:.3f}, alpha2={a2:.3f}")
            result = analyze_point(
                spec,
                layout,
                grid,
                k_grid,
                rfft_shape,
                base_params,
                d1,
                a2,
                k_idx,
            )
            results.append(result)
    return results


def find_instability_boundary_from_sweep(
    sweep_csv: str | Path,
) -> list[dict]:
    """Extract instability boundary from sweep results.

    The boundary is defined empirically: for each positive delta1 value,
    the highest alpha2 where the sweep diverged. This uses the modal
    solver's IC-projection-aware divergence guard, which is the correct
    criterion (simple eigenvalue thresholds fail because the GR baseline
    already has large real eigenvalues from wave propagation).

    Returns list of dicts with delta1, alpha2_crit (upper boundary),
    alpha2_lower (lower boundary if exists), and alpha2_first_valid.
    """
    sweep_csv = Path(sweep_csv)
    with sweep_csv.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    from collections import defaultdict

    by_d1: dict[float, list[tuple[float, str, str]]] = defaultdict(list)
    for r in rows:
        d1 = float(r["delta1"])
        a2 = float(r["alpha2"])
        status = r.get("run_status", "success")
        pmax = r.get("P_max", "")
        by_d1[d1].append((a2, status, pmax))

    results = []
    for d1 in sorted(by_d1.keys()):
        if d1 < 0:
            continue  # use delta1 -> -delta1 symmetry
        pts = sorted(by_d1[d1], key=operator.itemgetter(0))

        # Classify each alpha2 as diverged or valid
        diverged_a2 = []
        valid_a2 = []
        for a2, st, pmax in pts:
            if st == "diverged":
                diverged_a2.append(a2)
            elif pmax:
                try:
                    p = float(pmax)
                    if 0 < p < 0.1:
                        valid_a2.append(a2)
                except (ValueError, TypeError):
                    pass

        # Upper boundary: highest diverged alpha2 that is below a valid point
        # Lower boundary: lowest diverged alpha2 that is above a valid point
        upper_crit = float("nan")
        lower_crit = float("nan")
        first_valid = valid_a2[0] if valid_a2 else float("nan")

        if diverged_a2 and valid_a2:
            # Find upper boundary (stable window upper edge)
            candidates = [a for a in diverged_a2 if a > min(valid_a2)]
            if candidates:
                upper_crit = min(candidates)

            # Find lower boundary (stable window lower edge)
            candidates = [a for a in diverged_a2 if a < max(valid_a2)]
            if candidates:
                lower_crit = max(candidates)

        results.append(
            {
                "delta1": d1,
                "alpha2_upper": upper_crit,
                "alpha2_lower": lower_crit,
                "alpha2_first_valid": first_valid,
            }
        )
    return results


def find_zero_crossing(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    base_params: dict[str, float],
    alpha2: float,
    delta1_range: tuple[float, float] = (0.0, 1.5),
    n_delta1: int = 100,
    k_idx: int = 8,
) -> float | None:
    """Find delta1 where mu_eff crosses zero (suppression valley)."""
    d1_vals = np.linspace(*delta1_range, n_delta1)
    prev_mu = None
    for d1 in d1_vals:
        r = analyze_point(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            d1,
            alpha2,
            k_idx,
        )
        mu = r["mu_eff"]
        if prev_mu is not None and prev_mu * mu < 0:
            # Linear interpolation for zero crossing
            d1_prev = d1 - (d1_vals[1] - d1_vals[0])
            d1_crit = d1_prev - prev_mu * (d1 - d1_prev) / (mu - prev_mu)
            return float(d1_crit)
        prev_mu = mu
    return None


def mini_solver_amplification(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    base_params: dict[str, float],
    delta1: float,
    alpha2: float,
    *,
    ic_type: str = "gaussian",
    ic_amplitude: float = 0.1,
    ic_width: float = 5.0,
    ic_center: float = 50.0,
    t_end: float = 50.0,
    n_snapshots: int = 201,
) -> dict:
    """Compute amplification factor A by evolving the 4x4 h5+a1 block analytically.

    Matches the actual modal solver evolution:
    1. Build A_reduced for all k modes
    2. Extract 4x4 {a1, va1, h5, vh5} block at each mode
    3. Eigendecompose each block
    4. Construct IC matching the sweep (Gaussian or plane-wave)
    5. Evolve via exp(lambda*t)
    6. Compute P(t) = E_a1(t) / E_h5(0), take P_max
    7. Compare delta1>0 vs delta1=0 baseline

    Returns dict with P_max_torsion, P_max_GR, A = P_max_torsion/P_max_GR.
    """
    k_vals = k_grid[0]
    n_modes = len(k_vals)
    grid.dx[0]
    N = grid.shape[0]

    # h5=slot 8, vh5=slot 9, a1=slot 2, va1=slot 3 in reduced system
    np.array([2, 3, 8, 9])  # a1, va1, h5, vh5

    def get_full_system(params):
        """Get the full 14x14 reduced system and slot mapping."""
        coeff_eval = CoefficientEvaluator(spec, grid, params)
        A_red, _, _, _, _, mapping = _build_evolution_matrices(
            spec, layout, grid, coeff_eval, k_grid, rfft_shape
        )
        return A_red, mapping

    # Build IC in Fourier space for the FULL 14-slot system
    x = np.linspace(0, grid.bounds[0][1], N, endpoint=False)
    if ic_type == "gaussian":
        h5_x = ic_amplitude * np.exp(-((x - ic_center) ** 2) / (2 * ic_width**2))
        vh5_x = np.zeros(N)
    elif ic_type == "plane-wave":
        k0 = 2 * np.pi / (grid.bounds[0][1] - grid.bounds[0][0])
        h5_x = ic_amplitude * np.cos(k0 * x)
        vh5_x = ic_amplitude * k0 * np.sin(k0 * x)
    else:
        msg = f"Unknown IC type: {ic_type}"
        raise ValueError(msg)

    h5_hat = np.fft.rfft(h5_x) / N
    vh5_hat = np.fft.rfft(vh5_x) / N

    def evolve_and_measure(params):
        """Evolve full 14x14 reduced system, measure h5->a1 conversion."""
        A_red, mapping = get_full_system(params)
        n_slots = A_red.shape[1]  # 14
        t_eval = np.linspace(0, t_end, n_snapshots)

        # Map field names to reduced slots
        h5_r = mapping[layout.field_slot_map["h_5"]]
        vh5_r = mapping[layout.velocity_slot_map["h_5"]]
        a1_r = mapping[layout.field_slot_map["a_1"]]

        # IC per mode: only h5 and vh5 are nonzero
        ic_per_mode = np.zeros((n_modes, n_slots), dtype=complex)
        ic_per_mode[:, h5_r] = h5_hat
        ic_per_mode[:, vh5_r] = vh5_hat

        # Initial h5 energy (Parseval)
        e_h5_0 = np.sum(np.abs(h5_hat) ** 2) * N

        # Eigendecompose all modes with tachyonic suppression
        all_V_y0 = []
        all_eig_vals = []
        for m in range(n_modes):
            eigs, V = np.linalg.eig(A_red[m])
            Vinv = np.linalg.inv(V)
            c = Vinv @ ic_per_mode[m]

            # Tachyonic suppression
            max_c = np.max(np.abs(c))
            noise_floor = 1e-12 * max_c if max_c > 0 else 1e-30
            for j in range(n_slots):
                if eigs[j].real > 1e-8 and abs(c[j]) < noise_floor:
                    eigs[j] = 0.0

            all_eig_vals.append(eigs)
            all_V_y0.append(V * c[np.newaxis, :])

        all_eig_vals = np.array(all_eig_vals)  # (n_modes, n_slots)

        # Track a1 energy over time
        p_max = 0.0
        for t in t_eval:
            exp_lt = np.exp(all_eig_vals * t)  # (n_modes, n_slots)
            a1_hat_t = np.zeros(n_modes, dtype=complex)
            for m in range(n_modes):
                y_t = all_V_y0[m] @ exp_lt[m]
                a1_hat_t[m] = y_t[a1_r]

            e_a1_t = np.sum(np.abs(a1_hat_t) ** 2) * N
            p_t = e_a1_t / e_h5_0 if e_h5_0 > 1e-30 else 0.0
            p_max = max(p_max, p_t)

        return p_max

    # Run for both GR and torsion
    params_gr = {**base_params, "delta1": 0.0, "alpha2": alpha2}
    params_t = {**base_params, "delta1": delta1, "alpha2": alpha2}

    p_max_gr = evolve_and_measure(params_gr)
    p_max_t = evolve_and_measure(params_t)

    A = p_max_t / p_max_gr if p_max_gr > 1e-30 else float("nan")

    return {
        "delta1": delta1,
        "alpha2": alpha2,
        "P_max_GR": p_max_gr,
        "P_max_torsion": p_max_t,
        "A": A,
        "log10_A": np.log10(A) if A > 0 else float("nan"),
        "ic_type": ic_type,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spec",
        default="examples/data/torsion_gertsenshtein_nonminimal.json",
        help="Path to JSON spec",
    )
    parser.add_argument("--output", default=None, help="Output directory")
    parser.add_argument("--k-idx", type=int, default=8, help="Fourier mode index")
    parser.add_argument(
        "--mode",
        choices=[
            "point",
            "sweep",
            "boundary",
            "zero-crossing",
            "validate",
            "mini-solver",
        ],
        default="sweep",
        help="Analysis mode",
    )
    parser.add_argument("--delta1", type=float, default=1.0)
    parser.add_argument("--alpha2", type=float, default=-0.6)
    args = parser.parse_args()

    base_params = {
        "B0": 0.0001,
        "kappa": 1.0,
        "alpha1": 0.0,
        "alpha3": 1.0,
    }

    print("Loading spec...")
    spec_raw = load_equation_system(args.spec)
    # CRITICAL: normalize kinetic coefficients (divide RHS by LHS kinetic coeff).
    # Without this, h_5 appears to have wrong-sign mass (+k^2 instead of -k^2).
    # The simulation/sweep code does this automatically; we must do it explicitly.
    spec = normalize_kinetic_coefficients(spec_raw, base_params)
    grid = build_grid()
    layout = StateLayout.from_spec(spec, grid.num_points)
    k_grid, rfft_shape = get_k_grid(grid)
    k_vals = k_grid[0]

    if args.mode == "point":
        print(
            f"\nAnalyzing point delta1={args.delta1}, alpha2={args.alpha2}, k={k_vals[args.k_idx]:.4f}"
        )
        r = analyze_point(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            args.delta1,
            args.alpha2,
            args.k_idx,
        )
        for key, val in r.items():
            print(f"  {key}: {val}")

    elif args.mode == "sweep":
        print("\nSweeping (delta1, alpha2) parameter space...")
        d1_vals = np.linspace(0.0, 1.0, 21)
        a2_vals = np.linspace(-1.5, -0.3, 25)
        results = sweep_coupling(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            d1_vals,
            a2_vals,
            args.k_idx,
        )
        # Save results
        out_dir = (
            Path(args.output)
            if args.output
            else Path("examples/data/schur_complement_sweep")
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        fieldnames = list(results[0].keys())
        with (out_dir / "schur_results.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(results)
        print(f"\nSaved {len(results)} points to {out_dir / 'schur_results.csv'}")

    elif args.mode == "boundary":
        print("\nExtracting instability boundary from sweep data...")
        sweep_csv = Path("examples/data/nonminimal_heatmap_d1_a2_hires/results.csv")
        if not sweep_csv.exists():
            print(f"  Sweep data not found: {sweep_csv}")
            return
        results = find_instability_boundary_from_sweep(sweep_csv)
        print("\nInstability boundary (stable window edges):")
        print(f"  {'delta1':>7} {'upper':>8} {'lower':>8} {'first_valid':>12}")
        for r in results:
            u = (
                f"{r['alpha2_upper']:.3f}"
                if np.isfinite(r["alpha2_upper"])
                else "  none"
            )
            lo = (
                f"{r['alpha2_lower']:.3f}"
                if np.isfinite(r["alpha2_lower"])
                else "  none"
            )
            fv = (
                f"{r['alpha2_first_valid']:.3f}"
                if np.isfinite(r["alpha2_first_valid"])
                else "  none"
            )
            print(f"  {r['delta1']:+7.3f} {u:>8} {lo:>8} {fv:>12}")

    elif args.mode == "zero-crossing":
        print(f"\nFinding coupling zero-crossing at alpha2={args.alpha2}...")
        d1_crit = find_zero_crossing(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            args.alpha2,
            k_idx=args.k_idx,
        )
        if d1_crit is not None:
            print(f"  Zero-crossing at delta1 = {d1_crit:.4f}")
        else:
            print("  No zero-crossing found")

    elif args.mode == "validate":
        print("\nValidating against sweep data...")
        csv_path = Path("examples/data/nonminimal_heatmap_d1_a2_hires/results.csv")
        if not csv_path.exists():
            print(f"  Sweep data not found: {csv_path}")
            return
        with csv_path.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

        n_compared = 0
        errors = []
        for row in rows:
            log10_A_str = row.get("log10_A", "")
            if not log10_A_str:
                continue
            d1 = float(row["delta1"])
            a2 = float(row["alpha2"])
            A_sweep = 10 ** float(log10_A_str)

            r = analyze_point(
                spec,
                layout,
                grid,
                k_grid,
                rfft_shape,
                base_params,
                d1,
                a2,
                args.k_idx,
            )
            A_schur = abs(r["A_coupling"])
            if A_schur > 0 and A_sweep > 0:
                rel_err = abs(A_schur - A_sweep) / A_sweep
                errors.append(rel_err)
                n_compared += 1
                if n_compared <= 10 or rel_err > 0.5:
                    print(
                        f"  d1={d1:+.3f} a2={a2:.3f}: "
                        f"A_sweep={A_sweep:.2e} A_schur={A_schur:.2e} "
                        f"err={rel_err:.1%}"
                    )

        if errors:
            print(f"\n  Compared {n_compared} points")
            print(f"  Median relative error: {np.median(errors):.1%}")
            print(f"  Mean relative error: {np.mean(errors):.1%}")
            print(f"  Max relative error: {max(errors):.1%}")

    elif args.mode == "mini-solver":
        print(f"\nMini-solver: delta1={args.delta1}, alpha2={args.alpha2}")
        r = mini_solver_amplification(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            args.delta1,
            args.alpha2,
        )
        for key, val in r.items():
            print(f"  {key}: {val}")


if __name__ == "__main__":
    main()
