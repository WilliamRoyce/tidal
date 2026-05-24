#!/usr/bin/env python3
"""Phase E pre-flight safety check.

For the dual-Gaussian geometry, the leading-order Gertsenshtein conversion
amplitude per Gaussian peak is

    P_GR ~= sin^2( arg ),   arg = kappa * Bpeak * sigB * sqrt(2*pi) / 2.

If `arg` is too small the signal sits at numerical noise; too large and we
leave the perturbative regime where A = P_theory / P_GR is well-defined; and
near n*pi we are at a sin^2 node where P_GR vanishes (a geometric "null" that
masquerades as a physics result, per the user's catch).

This script reads `scripts/hpc_submit_drafts/v3e_localised/_geometry.env`,
computes `arg`, and asserts:

  * `arg in [1e-3, 0.3]`
  * `min_n |arg - n*pi|       > 0.3`  (clear of zero crossings of sin^2)
  * `min_n |arg - (n+0.5)*pi| > 0.3`  (clear of maxima — kept symmetric so any
                                        change that pushes us into a high-P
                                        regime also trips the check and is
                                        flagged before submit)

Exit 0 on PASS. Exit non-zero with a clear message on FAIL.

Usage:
    python scripts/v3e_boccaletti_preflight.py                # default env path
    python scripts/v3e_boccaletti_preflight.py <env-file>     # override
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

DEFAULT_ENV = (
    Path(__file__).resolve().parent
    / "hpc_submit_drafts"
    / "v3e_localised"
    / "_geometry.env"
)


def _parse_env(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        line = line.removeprefix("export ")
        key, _, value = line.partition("=")
        out[key.strip()] = value.split("#", 1)[0].strip().strip('"').strip("'")
    return out


def evaluate(env: dict[str, str], *, kappa: float = 1.0) -> dict[str, float]:
    bpeak = float(env["BPEAK"])
    sigb = float(env["SIGB"])
    arg = kappa * bpeak * sigb * math.sqrt(2 * math.pi) / 2.0
    # Distance to nearest sin^2 extremum (zeros at n*pi, maxima at (n+0.5)*pi).
    # The n=0 zero is at arg=0 — that *is* the perturbative regime we want
    # (P_GR ~ arg^2). The geometric-null concern is the next-and-beyond zeros
    # at n*pi for n >= 1. Same for maxima: (n+0.5)*pi for n >= 0 is the first
    # max at pi/2, which is non-perturbative full-conversion territory.
    n_zero = max(1, round(arg / math.pi))
    n_max = max(0, round((arg / math.pi) - 0.5))
    dist_zero = abs(arg - n_zero * math.pi)
    dist_max = abs(arg - (n_max + 0.5) * math.pi)
    return {
        "kappa": kappa,
        "Bpeak": bpeak,
        "sigB": sigb,
        "arg": arg,
        "P_GR_leading": math.sin(arg) ** 2,
        "dist_to_nearest_higher_zero": dist_zero,
        "dist_to_nearest_max": dist_max,
    }


def check(env_path: Path = DEFAULT_ENV) -> int:
    env = _parse_env(env_path)
    m = evaluate(env)
    failures: list[str] = []
    if not (1e-3 <= m["arg"] <= 0.3):
        failures.append(f"arg = {m['arg']:.4g} outside perturbative window [1e-3, 0.3]")
    if m["dist_to_nearest_higher_zero"] <= 0.3:
        failures.append(
            f"arg = {m['arg']:.4g} within 0.3 of sin^2 zero at n*pi (n>=1) "
            f"(distance {m['dist_to_nearest_higher_zero']:.4g})"
        )
    if m["dist_to_nearest_max"] <= 0.3:
        failures.append(
            f"arg = {m['arg']:.4g} within 0.3 of sin^2 maximum "
            f"(distance {m['dist_to_nearest_max']:.4g})"
        )
    print(f"Phase E Boccaletti preflight: env={env_path}")
    print(f"  Bpeak={m['Bpeak']:.4g}  sigB={m['sigB']:.4g}  kappa={m['kappa']:.4g}")
    print(
        f"  arg = kappa*Bpeak*sigB*sqrt(2*pi)/2 = {m['arg']:.4g} "
        f"(leading-order P_GR = {m['P_GR_leading']:.4g})"
    )
    if failures:
        print("FAIL:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("PASS")
    return 0


def main(argv: list[str]) -> int:
    env_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_ENV
    if not env_path.exists():
        print(f"ERROR: env file not found: {env_path}", file=sys.stderr)
        return 2
    return check(env_path)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv))
