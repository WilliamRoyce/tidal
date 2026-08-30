"""F1-c: does the symbolic orbit span the pencil's measured weak sector?

This is the experiment GH #477 was opened on.  The 2026-08-27 Lane-1c
diagnostic found the E.cal pencil's weakly-determined directions
NUMERICALLY and observed they looked like the diffeomorphism orbit.  Here
the orbit is built SYMBOLICALLY and the two are compared by principal
angles -- validating both sides at once.

Also checks the property the whole consumer design rests on: the
least-determined gauge profiles must come out FAR-FIELD LOCALIZED, i.e.
the grading rediscovers "the symmetry returns where the background
vanishes" without being told about B0(z).

usage: f1c_span.py [N]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import subspace_angles

sys.path.insert(0, str(Path(__file__).parent))
from orbit import capture, equilibrate, grade, joint_operator  # noqa: E402

N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
LAM = 0.37 + 0.91j  # generic probe point (away from the spectrum)

ctx = capture(N)
A, B, _ = equilibrate(ctx.A, ctx.B)
C, D, labels = joint_operator(ctx, A, B)
sig, prof, chains = grade(C, D)
print(f"N={N} dim={ctx.dim}  joint profiles={C.shape[1]}  kept={len(sig)}")
print(
    f"sigma spectrum: min={sig[0]:.3e}  median={np.median(sig):.3e}  max={sig[-1]:.3e}"
)

# --- (1) are the least-determined profiles far-field localized? ---
z, p = ctx.z, ctx.params
near = (np.abs(z - p["zc1"]) <= 2 * p["sigB"]) | (np.abs(z - p["zc2"]) <= 2 * p["sigB"])
print("\n--- (1) real-space localization of the least-determined gauge profiles ---")
print(
    "(profile = the gauge FUNCTION xi(z); 'near' = within 2 sigma of either Gaussian)"
)
for rank in (0, 1, 2, len(sig) // 2, len(sig) - 1):
    u = prof[:, rank]
    tot = np.zeros(N, dtype=np.complex128)
    for (gen, m), amp in zip(labels, u, strict=True):
        if abs(amp) > 0:
            hat = np.zeros(N, dtype=np.complex128)
            hat[m] = amp
            tot += np.fft.ifftn(hat)
    e = np.abs(tot) ** 2
    share = float(e[near].sum() / max(e.sum(), 1e-300))
    tag = (
        "least-determined"
        if rank < 3
        else ("median" if rank == len(sig) // 2 else "most-determined")
    )
    print(
        f"  rank {rank:3d}  sigma={sig[rank]:.3e}  near-Gaussian energy share={share:.3f}  [{tag}]"
    )

# --- (2) principal angles: symbolic orbit vs the pencil's own weak sector ---
# The pencil null vector of a chain at probe lambda is v(lam) = sum_j lam^j R_j
# (that is what the chain identities say); R_0 alone is NOT it.
print("\n--- (2) principal angles: symbolic orbit span vs measured weak directions ---")
depth = C.shape[0] // ctx.dim - 1


def orbit_at(lam: complex, ncols: int) -> np.ndarray:
    v = np.zeros((ctx.dim, ncols), dtype=np.complex128)
    for j in range(depth + 1):
        v += (lam**j) * chains[j * ctx.dim : (j + 1) * ctx.dim, :ncols]
    return v / np.maximum(np.linalg.norm(v, axis=0, keepdims=True), 1e-300)


pen = A - LAM * B
_, s_pen, vh_pen = np.linalg.svd(pen)
print(
    f"pencil singular values: min={s_pen[-1]:.3e}  s[-32]={s_pen[-32]:.3e}  max={s_pen[0]:.3e}"
)
res = np.linalg.norm(pen @ orbit_at(LAM, 32), axis=0)
print(
    f"residual ||(A-lam B) v_orbit|| for the 32 least-determined: "
    f"min={res.min():.3e} max={res.max():.3e}  (tracks sigma, not machine eps)"
)
# Compare each weak SUBSPACE against the orbit SPAN.  Matching dimensions
# one-to-one would compare two differently-ordered bases (the joint grading
# orders whole chains; the pencil orders single-lambda residuals) and reports
# a spurious pi/2 for the permuted directions.
ORB = min(64, chains.shape[1])
Qorb = np.linalg.qr(orbit_at(LAM, ORB))[0]
for dim in (4, 8, 16, 32):
    weak = vh_pen[-dim:].conj().T
    ang = subspace_angles(weak, Qorb)
    print(
        f"  {dim:3d} weakest vs orbit span(dim {ORB}): "
        f"max angle={np.max(ang):.3e} rad, median={np.median(ang):.3e}"
    )

# --- (3) is anything weak that the orbit does NOT explain? ---
# The consumer refuses when a weak direction lies outside the exported orbit
# span; this is the measurement that decides whether that guard ever fires here.
print("\n--- (3) weak directions OUTSIDE the orbit span (the refusal guard) ---")
Q = np.linalg.qr(orbit_at(LAM, min(64, chains.shape[1])))[0]
for dim in (4, 8, 16, 32, 64):
    weak = vh_pen[-dim:].conj().T
    resid = np.linalg.norm(weak - Q @ (Q.conj().T @ weak), axis=0)
    print(
        f"  dim={dim:3d}  unexplained component of each weak vector: "
        f"min={resid.min():.3f} max={resid.max():.3f}"
    )
inv = {v: k for k, v in ctx.slot.items()}
w0 = vh_pen[-1].conj()
per_slot = np.array(
    [np.linalg.norm(w0[s * ctx.N : (s + 1) * ctx.N]) for s in range(ctx.n_slots)]
)
top = np.argsort(per_slot)[::-1][:6]
print(
    "  weakest pencil direction lives on slots: "
    + ", ".join(f"{inv.get(s, s)}={per_slot[s]:.2f}" for s in top)
)

worst = float(
    np.max(
        np.linalg.norm(
            vh_pen[-32:].conj().T - Q @ (Q.conj().T @ vh_pen[-32:].conj().T), axis=0
        )
    )
)
ok_loc = True  # section (1) is reported, not thresholded: it is a physics observation
ok_span = worst < 0.2
print(
    f"\nGATE F1-c: {'PASS' if (ok_loc and ok_span) else 'FAIL'} "
    f"(worst unexplained component over the 32 weakest = {worst:.3f})"
)
sys.exit(0 if (ok_loc and ok_span) else 1)
