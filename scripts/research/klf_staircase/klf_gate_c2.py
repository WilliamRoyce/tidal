"""Gate (c) stage-2: graded recursion — deep-tau staircase, then re-staircase
the finite block at shallow tau2 to pin the weakly-determined sector.

usage: klf_gate_c2.py N tau1 tau2
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import klf_gate_c as g
from klf_port import compose, klf_right

TAU1 = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10
TAU2 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6

A0, B0 = g.A0, g.B0
n_tot = A0.shape[0]
mask = np.zeros(n_tot, dtype=bool)
mask[g.obs_idx] = True

s = np.maximum(np.max(np.abs(A0), axis=1), np.max(np.abs(B0), axis=1))
s[s == 0.0] = 1.0
An, Bn = A0 / s[:, None], B0 / s[:, None]

t0 = time.perf_counter()
dec1 = klf_right(An, Bn, TAU1)
nf1 = dec1["nf"]
roff1, coff1 = dec1["roff"], dec1["coff"]
fc = slice(coff1, coff1 + nf1)
fr = slice(roff1, roff1 + nf1)
W1 = dec1["Z"][:, fc]
Af = dec1["M"][fr, fc]
Ef = dec1["N"][fr, fc]
assert roff1 == 0, f"stage-1 right structure unexpected: roff={roff1}"

# stage 2 on the clean finite pencil (equilibrate rows again)
s2 = np.maximum(np.max(np.abs(Af), axis=1), np.max(np.abs(Ef), axis=1))
s2[s2 == 0.0] = 1.0
Af2, Ef2 = Af / s2[:, None], Ef / s2[:, None]
dec2 = klf_right(Af2, Ef2, TAU2)
A2_eff, proj2, pinned2, info2 = compose(dec2, TAU2)
dt = time.perf_counter() - t0

# overall generator on original coordinates
A_tot = W1 @ A2_eff @ W1.conj().T
proj_tot = W1 @ proj2 @ W1.conj().T
pinned_tot = W1 @ pinned2

# contract against ALL original (equilibrated) rows
rng = np.random.default_rng(0)
Y = proj_tot @ (rng.standard_normal((n_tot, 8)) + 1j * rng.standard_normal((n_tot, 8)))
contract = float(
    np.linalg.norm(Bn @ (A_tot @ Y) - An @ Y) / max(np.linalg.norm(An @ Y), 1e-300)
)

blk = A_tot[np.ix_(mask, mask)]
blk_rel = np.linalg.norm(blk - A0[np.ix_(mask, mask)]) / np.linalg.norm(
    A0[np.ix_(mask, mask)]
)
cross = max(
    np.linalg.norm(A_tot[np.ix_(mask, ~mask)]),
    np.linalg.norm(A_tot[np.ix_(~mask, mask)]),
)
ev = np.linalg.eigvals(A2_eff)
max_re = float(np.max(ev.real)) if ev.size else 0.0
# pin overlap with the obs sector must be ~0 (pins are far-field gauge)
pin_obs = float(np.linalg.norm(pinned_tot[mask, :])) if pinned_tot.shape[1] else 0.0

print(
    f"N={g.N} dim={n_tot} tau1={TAU1:.0e} tau2={TAU2:.0e}: nf1={nf1} "
    f"nf2={info2['n_finite']} pinned={pinned_tot.shape[1]} contract={contract:.2e}  "
    f"maxRe={max_re:+.3e} |Im|max={float(np.max(np.abs(ev.imag))) if ev.size else 0:.3f}  "
    f"obs rel={blk_rel:.2e} cross={cross:.2e} pin-obs={pin_obs:.2e}  [{dt:.1f}s]"
)
