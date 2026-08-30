"""F1-a: chain identities on the vacuum (Bpeak=0) pencil.

Gate: every diffeomorphism chain must annihilate the pure-EH operator at
machine precision for every gauge-profile mode -- this is the E[Phibar]=0
limit of the linearized Noether identity, and it isolates candidate-SHAPE
correctness from background breaking.  The U(1) chain is the NEGATIVE
control: Lorenz gauge fixing is in the action, so chi must FAIL at O(1)
even at Bpeak=0 (breaking_source = gauge_fixing, not background).

usage: f1a_exactness.py [N ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from orbit import GENERATORS, capture, chain, equilibrate, sigma  # noqa: E402

NS = [int(a) for a in sys.argv[1:]] or [16, 24, 32]
TOL = 1e-12

status = 0
for n in NS:
    ctx = capture(n, bpeak=0.0)
    A, B, _ = equilibrate(ctx.A, ctx.B)
    scale = max(np.linalg.norm(A, 1), np.linalg.norm(B, 1))
    print(f"\nN={n}  dim={ctx.dim}  ||(A,B)||_1={scale:.3e}   [Bpeak=0 vacuum pencil]")
    for gen in GENERATORS:
        worst, worst_mode = 0.0, -1
        for m in range(n):
            R = chain(ctx, gen, m)
            s = sigma(A, B, R) / scale
            if s > worst:
                worst, worst_mode = s, m
        expect_exact = gen != "chi_u1"
        ok = (worst <= TOL) if expect_exact else (worst > 1e-3)
        label = "EXACT" if expect_exact else "must-BREAK (Lorenz control)"
        print(
            f"  {gen:8s} depth={len(chain(ctx, gen, 1)):d}  "
            f"worst rel-defect={worst:.2e} @mode {worst_mode:2d}   "
            f"{label:26s} {'PASS' if ok else 'FAIL'}"
        )
        status |= 0 if ok else 1

print("\nGATE F1-a:", "PASS" if status == 0 else "FAIL")
sys.exit(status)
