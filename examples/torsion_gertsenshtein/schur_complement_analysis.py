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

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import _build_constraint_eliminated_matrices
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import load_equation_system


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
    A_reduced, _, _, _, mapping = _build_constraint_eliminated_matrices(
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
    max_re = np.max(np.real(all_eigs), axis=1)

    # Growth threshold for stability (matching modal solver)
    t_end = base_params.get("t_end", 50.0)
    growth_threshold = np.log(1e8) / t_end  # ~0.37

    # Find most unstable mode
    most_unstable_k = np.argmax(max_re)
    max_re_eig = max_re[most_unstable_k]

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
        "most_unstable_k": k_vals[most_unstable_k],
        "stable": bool(max_re_eig < growth_threshold),
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


def find_instability_boundary(
    spec,
    layout: StateLayout,
    grid: GridInfo,
    k_grid: list[np.ndarray],
    rfft_shape: tuple[int, ...],
    base_params: dict[str, float],
    delta1_values: np.ndarray,
    alpha2_range: tuple[float, float] = (-3.0, -0.3),
    n_alpha2: int = 100,
    k_idx: int = 8,
) -> list[dict]:
    """Find the instability boundary alpha2_crit(delta1)."""
    alpha2_vals = np.linspace(*alpha2_range, n_alpha2)
    results = []
    for d1 in delta1_values:
        prev_stable = None
        boundary_a2 = None
        for a2 in alpha2_vals:
            r = analyze_point(
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
            if prev_stable is not None and not prev_stable and r["stable"]:
                boundary_a2 = a2
                break
            prev_stable = r["stable"]

        results.append(
            {
                "delta1": d1,
                "alpha2_crit": boundary_a2 if boundary_a2 is not None else float("nan"),
            }
        )
        print(f"  delta1={d1:.3f}: boundary at alpha2={boundary_a2}")
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
        choices=["point", "sweep", "boundary", "zero-crossing", "validate"],
        default="sweep",
        help="Analysis mode",
    )
    parser.add_argument("--delta1", type=float, default=1.0)
    parser.add_argument("--alpha2", type=float, default=-0.6)
    args = parser.parse_args()

    print("Loading spec...")
    spec = load_equation_system(args.spec)
    grid = build_grid()
    layout = StateLayout.from_spec(spec, grid.num_points)
    k_grid, rfft_shape = get_k_grid(grid)
    k_vals = k_grid[0]

    base_params = {
        "B0": 0.0001,
        "kappa": 1.0,
        "alpha1": 0.0,
        "alpha3": 1.0,
    }

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
        print("\nFinding instability boundary alpha2_crit(delta1)...")
        d1_vals = np.linspace(0.05, 1.0, 20)
        results = find_instability_boundary(
            spec,
            layout,
            grid,
            k_grid,
            rfft_shape,
            base_params,
            d1_vals,
            k_idx=args.k_idx,
        )
        print("\nInstability boundary:")
        for r in results:
            print(f"  delta1={r['delta1']:.3f}: alpha2_crit={r['alpha2_crit']:.4f}")

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


if __name__ == "__main__":
    main()
