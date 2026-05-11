#!/usr/bin/env python3
r"""Generate a v2-vs-v3 comparison table from two ``hpc_results/<jobid>/`` directories.

Reads ``inference.json`` from each (PolyChord nested-sampling output written
by ``tidal sample``) and emits a markdown row per parameter capturing log Z,
ESS, MAP shift in σ-units, per-coupling marginal D_KL, and prior compactness.

Usage
-----
    uv run python scripts/v3_v2_comparison.py \\
        --v2 hpc_results/28896653 --v3 hpc_results/29149987 \\
        --label "D1 amp" --output docs/comparison/d1_amp_v2_v3.md

Notes
-----
* The MAP-shift σ-unit metric uses v3's posterior standard deviation (from
  the 95% credible interval / 3.92) as the denominator — v3 is the wider
  prior and therefore the more conservative reference for "did the MAP
  meaningfully move?"
* Parameter names must match between v2 and v3 ``param_names``.  Extra v3
  params (e.g. if a Lagrangian was de-pruned) are reported as "(v3-only)";
  missing v3 params are flagged as warnings.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


def _load_inference(run_dir: Path) -> dict[str, Any]:
    inference_path = run_dir / "inference.json"
    if not inference_path.exists():
        msg = f"inference.json not found in {run_dir}"
        raise FileNotFoundError(msg)
    return json.loads(inference_path.read_text())


def _ci_to_sigma(ci_lo: float, ci_hi: float) -> float:
    """Convert a 95% credible interval width to a Gaussian-equivalent sigma.

    For an approximately-Gaussian marginal, the 95% CI spans ~3.92σ
    (1.96σ either side).  Returns NaN if the CI is degenerate.
    """
    width = ci_hi - ci_lo
    if not math.isfinite(width) or width <= 0:
        return float("nan")
    return width / 3.92


def _prior_compactness(prior: dict[str, Any]) -> str:
    """Single-line summary of a prior's effective coupling-space coverage."""
    dist = prior.get("distribution", "?")
    lo = prior.get("low", float("nan"))
    hi = prior.get("high", float("nan"))
    if dist == "arctan_uniform":
        # Cauchy-tail: tan(89°) ~ 57; effectively full real line
        effective = math.tan(math.radians(hi)) if hi > 0 else float("inf")
        return f"arctan[{lo:.0f}°..{hi:.0f}°] → ±{effective:.1f}"
    if dist == "log_uniform":
        return f"log[{lo:.0e}..{hi:.0e}]"
    if dist == "uniform":
        return f"uniform[{lo:g}..{hi:g}]"
    return f"{dist}[{lo}..{hi}]"


def render_table(  # noqa: PLR0914
    v2: dict[str, Any],
    v3: dict[str, Any],
    *,
    label: str,
    v2_jobid: str,
    v3_jobid: str,
) -> str:
    """Format the comparison as GitHub-flavoured markdown."""
    lines: list[str] = []
    lines.extend(
        (
            f"# {label} — v2 vs v3 comparison",
            "",
            f"- **v2 reference**: `hpc_results/{v2_jobid}/`",
            f"- **v3 chain**: `hpc_results/{v3_jobid}/`",
            "",
            "## Headline",
            "",
            "| Metric | v2 | v3 | Δ |",
            "| --- | --- | --- | --- |",
        )
    )
    lz2 = v2.get("log_evidence", float("nan"))
    lz3 = v3.get("log_evidence", float("nan"))
    lze2 = v2.get("log_evidence_err", float("nan"))
    lze3 = v3.get("log_evidence_err", float("nan"))
    ess2 = v2.get("effective_sample_size", float("nan"))
    ess3 = v3.get("effective_sample_size", float("nan"))
    pi2 = v2.get("parameter_importance", {})
    pi3 = v3.get("parameter_importance", {})
    dkl2 = pi2.get("d_kl", float("nan"))
    dkl3 = pi3.get("d_kl", float("nan"))
    lines.extend(
        (
            f"| log Z | {lz2:+.3f} ± {lze2:.3f} | {lz3:+.3f} ± {lze3:.3f} | {lz3 - lz2:+.2f} nats |",
            f"| ESS | {ess2:.0f} | {ess3:.0f} | {ess3 - ess2:+.0f} |",
            f"| Joint D_KL | {dkl2:.2f} nats | {dkl3:.2f} nats | {dkl3 - dkl2:+.2f} nats |",
            f"| n_samples | {v2.get('n_samples', '?')} | {v3.get('n_samples', '?')} | — |",
            "",
            "## Per-coupling",
            "",
            "| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |",  # noqa: RUF001
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        )
    )

    v2_priors = {p["name"]: p for p in v2.get("priors", [])}
    v3_priors = {p["name"]: p for p in v3.get("priors", [])}
    v2_map = v2.get("map_estimate", {})
    v3_map = v3.get("map_estimate", {})
    v2.get("credible_interval_95", {})
    v3_ci = v3.get("credible_interval_95", {})
    v2_dkl = pi2.get("marginal_d_kl", {})
    v3_dkl = pi3.get("marginal_d_kl", {})

    all_params = sorted(set(v2_priors) | set(v3_priors))
    for p in all_params:
        v2p = _prior_compactness(v2_priors[p]) if p in v2_priors else "(v3-only)"
        v3p = _prior_compactness(v3_priors[p]) if p in v3_priors else "(v2-only)"
        m2 = v2_map.get(p, float("nan"))
        m3 = v3_map.get(p, float("nan"))
        ci3 = v3_ci.get(p, [float("nan"), float("nan")])
        sigma3 = _ci_to_sigma(ci3[0], ci3[1])
        shift = (
            (m3 - m2) / sigma3 if math.isfinite(sigma3) and sigma3 > 0 else float("nan")
        )
        d2 = v2_dkl.get(p, float("nan"))
        d3 = v3_dkl.get(p, float("nan"))
        lines.append(
            f"| {p} | {v2p} | {v3p} | {m2:+.3g} | {m3:+.3g} | {shift:+.2f} | {d2:.2f} | {d3:.2f} |"
        )

    lines.extend(
        (
            "",
            "## Notes",
            "",
            "* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.",
            "* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.",  # noqa: RUF001
            "* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--v2", type=Path, required=True, help="v2 hpc_results/<jobid>/ directory"
    )
    parser.add_argument(
        "--v3", type=Path, required=True, help="v3 hpc_results/<jobid>/ directory"
    )
    parser.add_argument(
        "--label", type=str, default="comparison", help="Campaign label for header"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output markdown path (defaults to stdout)",
    )
    args = parser.parse_args()

    v2 = _load_inference(args.v2)
    v3 = _load_inference(args.v3)
    out = render_table(
        v2,
        v3,
        label=args.label,
        v2_jobid=args.v2.name,
        v3_jobid=args.v3.name,
    )

    if args.output is None:
        sys.stdout.write(out)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(out)
        print(f"Wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
