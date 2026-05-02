r"""Random 9-D joint pre-flight check for the D2.3 composed prior.

Samples ``N=2000`` random draws from the full D2.3 prior block specified in
``docs/tex/pgt_stability_priors.tex`` §"D2.3 (Full T5 YM--PGT, all enabled)",
applies the D2.0 soft constraints, and runs ``check_conversion_stability``
against the canonical T5 spec at the same probe settings used by the
single/pair sensitivity scans.

Goal: estimate the joint prior survival rate (post-constraint × post-probe)
and confirm that no new 9-D-only joint instability exists. Single-parameter
and pairwise scans (D2.0/D2.1/D2.2) are necessary but not sufficient; this
script samples the full corner cube to catch any joint coupling missed by the
2-D recipe.

Run:  ``uv run python examples/data/d2_prior_scan/d23_joint_check/_analyze.py``

Outputs:
  d23_joint_survival.csv  — per-sample row (9 params + β-pass + γ_eff + verdict)
  SUMMARY.md              — survival rates and decision verdict (proceed/warn/stop)
"""

from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

from tidal.measurement._stability import check_conversion_stability
from tidal.solver.grid import GridInfo
from tidal.symbolic._spec_cache import load_spec_cached

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SPEC_PATH = REPO / "examples/data/torsion_gertsenshtein_general_nonminimal.json"

GRID = GridInfo(shape=(256,), bounds=((0.0, 100.0),), periodic=(True,))
IC_K = 2 * math.pi / 100.0  # 0.06283185307179587
KAPPA = 1.0
B0 = 0.01

N_SAMPLES = 2000
RNG_SEED = 20260502  # date-keyed for reproducibility

# Decision-gate thresholds (effective survival = post-constraint AND post-probe)
GATE_PROCEED = 0.05
GATE_WARN = 0.01

# ---------------------------------------------------------------------------
# Prior bounds — exact D2.3 composed block from docs/tex/pgt_stability_priors.tex
# ---------------------------------------------------------------------------

PRIOR_BOUNDS = {
    "beta1": ("uniform", -1.5, 1.0),
    "beta2": ("uniform", -3.0, -0.3),
    "beta3": ("uniform", -1.0, 1.0),
    "xi": ("log_uniform", 0.01, 10.0),
    "delta1": ("uniform", -0.025, 0.025),
    "chi": ("uniform", -0.009, 0.009),
    "zeta1": ("uniform", -0.050, 0.050),
    "zeta2": ("uniform", -0.021, 0.021),
    "zeta3": ("uniform", -0.050, 0.050),
}

PARAM_NAMES = list(PRIOR_BOUNDS.keys())


def _draw(rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
    """Draw n samples from the composed prior — vectorised."""
    out: dict[str, np.ndarray] = {}
    for name, (kind, lo, hi) in PRIOR_BOUNDS.items():
        if kind == "uniform":
            out[name] = rng.uniform(lo, hi, size=n)
        elif kind == "log_uniform":
            out[name] = np.exp(rng.uniform(math.log(lo), math.log(hi), size=n))
        else:
            msg = f"unknown prior kind {kind}"
            raise ValueError(msg)
    return out


# ---------------------------------------------------------------------------
# Soft constraints (PolyChord-side: --constraint clauses from D2.0)
# ---------------------------------------------------------------------------


def _beta_constraints_pass(b1: float, b2: float, b3: float, xi: float) -> bool:
    """Same four soft constraints PolyChord sees in D2.0/D2.1/D2.2/D2.3."""
    return (
        (b1 + b2 < -0.5)
        and (b1 - b2 > 0.25)
        and (b3 + (2 * b1 + b2) / 3 < 1.333)
        and (xi > 0.01)
    )


# ---------------------------------------------------------------------------
# Stability probe — same settings as zeta_scan/_analyze.py
# ---------------------------------------------------------------------------


def _probe(spec: object, params: dict[str, float]) -> tuple[float, bool]:
    result = check_conversion_stability(
        spec,  # type: ignore[arg-type]
        GRID,
        params,
        source="h_5",
        target="a_1",
        ic_wavevector=IC_K,
    )
    return float(result.max_excess), bool(result.stable)


def _params_dict(row: dict[str, float]) -> dict[str, float]:
    return {
        "kappa": KAPPA,
        "B0": B0,
        "eta": 0.0,
        **{k: float(row[k]) for k in PARAM_NAMES},
    }


def _fmt(x: float) -> str:
    if not math.isfinite(x):
        return "inf"
    return f"{x:.6g}"


# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------


def run_scan(spec: object, n_samples: int = N_SAMPLES) -> Path:
    rng = np.random.default_rng(RNG_SEED)
    draws = _draw(rng, n_samples)

    out = HERE / "d23_joint_survival.csv"
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "sample",
                *PARAM_NAMES,
                "beta_constraint_pass",
                "max_excess",
                "probe_pass",
            ]
        )

        for i in range(n_samples):
            row = {k: float(draws[k][i]) for k in PARAM_NAMES}
            beta_pass = _beta_constraints_pass(
                row["beta1"], row["beta2"], row["beta3"], row["xi"]
            )
            if beta_pass:
                params = _params_dict(row)
                exc, stable = _probe(spec, params)
            else:
                exc, stable = float("nan"), False
            w.writerow(
                [
                    i,
                    *(f"{row[k]:.6g}" for k in PARAM_NAMES),
                    int(beta_pass),
                    _fmt(exc),
                    int(stable),
                ]
            )

            if (i + 1) % 200 == 0:
                print(f"  [{i + 1}/{n_samples}] processed")

    return out


# ---------------------------------------------------------------------------
# Analysis + SUMMARY.md
# ---------------------------------------------------------------------------


def _decision_verdict(eff_survival: float) -> tuple[str, str]:
    if eff_survival >= GATE_PROCEED:
        return (
            "PROCEED",
            f"effective survival ≥ {GATE_PROCEED:.0%} → D2.3-B unchanged",
        )
    if eff_survival >= GATE_WARN:
        return (
            "PROCEED-WITH-WARN",
            f"effective survival in [{GATE_WARN:.0%}, {GATE_PROCEED:.0%}) → submit, but "
            "warn in slurm output that the chain will be slower than D2.0/D2.1/D2.2",
        )
    return (
        "STOP",
        f"effective survival < {GATE_WARN:.0%} → likely a 9-D-only joint coupling missed "
        "by 2-D scans; tighten priors or add a constraint before submitting",
    )


def _analyse(csv_path: Path) -> dict[str, object]:
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    n_total = len(rows)
    n_post_beta = sum(1 for r in rows if int(r["beta_constraint_pass"]) == 1)
    n_post_probe = sum(
        1
        for r in rows
        if int(r["beta_constraint_pass"]) == 1 and int(r["probe_pass"]) == 1
    )

    rate_beta = n_post_beta / n_total
    rate_probe_given_beta = (n_post_probe / n_post_beta) if n_post_beta else 0.0
    eff_survival = n_post_probe / n_total

    # Top-5 highest finite max_excess (closest to the rejection threshold)
    # among samples that passed beta-constraints. Tells us which corners of the
    # 9-D prior cube are nearly tachyonic.
    surv_rows = [r for r in rows if int(r["beta_constraint_pass"]) == 1]
    finite = [
        r for r in surv_rows if r["max_excess"] != "inf" and r["max_excess"] != "nan"
    ]
    finite.sort(key=lambda r: float(r["max_excess"]), reverse=True)
    top5 = finite[:5]

    return {
        "n_total": n_total,
        "n_post_beta": n_post_beta,
        "n_post_probe": n_post_probe,
        "rate_beta": rate_beta,
        "rate_probe_given_beta": rate_probe_given_beta,
        "eff_survival": eff_survival,
        "top5": top5,
    }


def _write_summary_md(csv_path: Path, stats: dict[str, object]) -> Path:
    out = HERE / "SUMMARY.md"
    eff = float(stats["eff_survival"])  # type: ignore[arg-type]
    verdict, justification = _decision_verdict(eff)

    lines: list[str] = [
        "# D2.3 Joint Pre-Flight Check — Summary\n",
        "Scan: random 9-D draws from the composed D2.3 prior, with D2.0 soft "
        "constraints applied, then probed via `check_conversion_stability` at "
        "the canonical settings (κ=1, B₀=0.01, h_5→a_1, "
        "ic_wavevector=2π/L=0.0628…).\n",
        f"Source: `{csv_path.relative_to(REPO)}`  ",
        "Recipe: `docs/tex/pgt_stability_priors.tex` §D2.3  ",
        "Template: `examples/data/d2_prior_scan/zeta_scan/_analyze.py`\n",
        "## Survival Rates\n",
        "| Stage | Count | Rate |",
        "|-------|------:|-----:|",
        f"| Total draws | {stats['n_total']} | 100.00% |",
        f"| Post β-constraints | {stats['n_post_beta']} | "
        f"{float(stats['rate_beta']):.2%} |",
        f"| Post stability probe (given β) | {stats['n_post_probe']} | "
        f"{float(stats['rate_probe_given_beta']):.2%} |",
        f"| **Effective survival** | **{stats['n_post_probe']}** | **{eff:.2%}** |",
        "",
        "Reference rates: D2.0/D2.1 prior stability sweep ≈ 60% rejected "
        "(40% effective survival); D2.2 sweep 59.1% rejected (40.9% effective).\n",
        "## Top-5 Surviving Samples by γ_eff\n",
        "These are the samples that passed both the β-constraints and the "
        "stability probe but with the largest max_excess values (closest to the "
        "γ_eff > 0.3 rejection threshold). If the top-5 cluster on a particular "
        "(param₁, param₂) corner, that signals a near-instability worth a "
        "tightening or constraint.\n",
        "| sample | β₁ | β₂ | β₃ | ξ | δ₁ | χ | ζ₁ | ζ₂ | ζ₃ | γ_eff |",
        "|--------|----|----|----|---|----|---|----|----|----|-------|",
    ]
    top5 = stats["top5"]
    assert isinstance(top5, list)
    lines.extend(
        f"| {r['sample']} | {r['beta1']} | {r['beta2']} | {r['beta3']} | "
        f"{r['xi']} | {r['delta1']} | {r['chi']} | {r['zeta1']} | "
        f"{r['zeta2']} | {r['zeta3']} | {r['max_excess']} |"
        for r in top5
    )

    lines.extend(
        (
            "",
            "## Decision Verdict\n",
            f"**Verdict: {verdict}**\n",
            f"{justification}\n",
            "Gate thresholds: PROCEED ≥ 5%, PROCEED-WITH-WARN ∈ [1%, 5%), "
            "STOP < 1% effective survival.\n",
        )
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> None:
    spec = load_spec_cached(SPEC_PATH)
    print(f"Loaded spec: {SPEC_PATH.name}")
    print(f"Drawing {N_SAMPLES} random samples from the composed D2.3 prior")
    print(f"Probe settings: κ=1, B₀=0.01, h_5→a_1, ic_wavevector={IC_K:.10g}\n")

    csv_path = run_scan(spec)
    n_rows = sum(1 for _ in csv_path.open()) - 1
    print(f"\n  wrote {csv_path.name} ({n_rows} rows)")

    stats = _analyse(csv_path)
    print("\n=== D2.3 joint pre-flight survival ===\n")
    print(f"  Total draws:                {stats['n_total']}")
    print(
        f"  Post β-constraints:         {stats['n_post_beta']} "
        f"({float(stats['rate_beta']):.2%})"
    )
    print(
        f"  Post stability probe:       {stats['n_post_probe']} "
        f"({float(stats['rate_probe_given_beta']):.2%} of post-β)"
    )
    eff = float(stats["eff_survival"])  # type: ignore[arg-type]
    print(f"  Effective survival:         {eff:.2%}")

    verdict, justification = _decision_verdict(eff)
    print(f"\n  Verdict: {verdict}")
    print(f"  {justification}")

    md = _write_summary_md(csv_path, stats)
    print(f"\n  wrote {md.name}")


if __name__ == "__main__":
    main()
