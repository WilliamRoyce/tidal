"""F1-b: is the chain defect the background-EOM violation?

The linearized Noether identity fixes the defect to be exactly linear in
the background-EOM residual E[Phibar].  For E.cal that residual has two
pieces with DIFFERENT powers of the background amplitude:

    Einstein violation   -kappa^2 T^EM ~ B0^2   ->  h-row defect ~ Bpeak^2
    Maxwell violation    J = curl B  ~ d_z B0   ->  a-row defect ~ Bpeak^1

so the identity makes a sharp, falsifiable prediction about SCALING that
no other origin of the defect (wrong chain shape, discretization
artifact, convention error) would reproduce.  It also predicts the defect
is spatially LOCALIZED where the background lives -- the statement that
the gauge symmetry returns as B0 -> 0, which is the whole reason the
far-field directions are weakly determined.

usage: f1b_breaking.py [N] [mode]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from orbit import (  # noqa: E402
    BASE_PARAMS,
    GENERATORS,
    capture,
    chain,
    defect,
    equilibrate,
)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
MODE = int(sys.argv[2]) if len(sys.argv) > 2 else 1
BPEAKS = [0.0025, 0.005, 0.01, 0.02]
FAR = 4.0  # far field := |z - zc| > FAR*sigB for both Gaussians


def row_mask(ctx, family: str) -> np.ndarray:
    """Row mask for equations belonging to a field family ('a' or 'h')."""
    m = np.zeros(ctx.dim, dtype=bool)
    for name, s in ctx.slot.items():
        base = name.removeprefix("v_")
        if base.startswith(family):
            m[s * ctx.N : (s + 1) * ctx.N] = True
    return m


def split_defect(ctx, A, B, R):
    """Per-row-class defect norms, normalized by the chain norm."""
    d = defect(A, B, R)
    c = max(np.linalg.norm(np.concatenate(R)), 1e-300)
    nblk = len(d) // ctx.dim
    blocks = d.reshape(nblk, ctx.dim)
    out = {}
    for fam in ("a", "h"):
        out[fam] = float(np.linalg.norm(blocks[:, row_mask(ctx, fam)]) / c)
    return out, blocks / c


print(f"N={N}  gauge profile = Fourier mode {MODE}")
print("\n--- (1) SCALING: d log||defect|| / d log Bpeak  (identity predicts a:1, h:2) ---")
records: dict[str, dict[str, list[float]]] = {g: {"a": [], "h": []} for g in GENERATORS}
ctxs = {}
for bp in BPEAKS:
    ctx = capture(N, bpeak=bp)
    ctxs[bp] = ctx
    A, B, _ = equilibrate(ctx.A, ctx.B)
    for g in GENERATORS:
        parts, _ = split_defect(ctx, A, B, chain(ctx, g, MODE))
        for fam in ("a", "h"):
            records[g][fam].append(parts[fam])

lb = np.log(BPEAKS)
ZERO = 1e-12  # a row class the orbit does not touch at all


def _slope(y: list[float]) -> tuple[float, str]:
    """Log-log slope, or 'exact' when the class is identically zero."""
    arr = np.array(y)
    if np.max(arr) < ZERO:
        return float("nan"), "exact-0"
    return float(np.polyfit(lb, np.log(np.maximum(arr, 1e-300)), 1)[0]), "slope"


print(f"{'generator':10s} {'a-row':>14s} {'h-row':>14s}   verdict")
status = 0
for g in GENERATORS:
    got = {fam: _slope(records[g][fam]) for fam in ("a", "h")}
    txt = {
        fam: ("identically 0" if kind == "exact-0" else f"{val:.3f}")
        for fam, (val, kind) in got.items()
    }
    if g == "chi_u1":
        ok, verdict = True, "(U(1): broken by the ACTION, background-independent)"
    else:
        # The identity predicts a:1 (Maxwell violation J ~ d_z B0) and h:2
        # (Einstein violation T^EM ~ B0^2).  A row class the orbit never
        # reaches has NO defect to scale -- that is exact satisfaction, the
        # strongest possible outcome, not a failed fit.
        a_ok = got["a"][1] == "exact-0" or abs(got["a"][0] - 1.0) < 0.05
        h_ok = got["h"][1] == "exact-0" or abs(got["h"][0] - 2.0) < 0.05
        ok = a_ok and h_ok
        verdict = "MATCHES Noether prediction" if ok else "DOES NOT match (a:1, h:2)"
        status |= 0 if ok else 1
    print(f"{g:10s} {txt['a']:>14s} {txt['h']:>14s}   {verdict}")

print(
    "  note: xi_x reaches NO photon row (the a-rows read h_2, h_5, h_0/h_4/h_7/h_9;\n"
    "        the xi_x orbit is h_1,h_6,v_h_6) -- so its Maxwell-violation defect is\n"
    "        exactly zero rather than O(Bpeak)."
)

print("\n--- (2) LOCALIZATION: fraction of defect energy in the far field ---")
ctx = ctxs[BASE_PARAMS["Bpeak"]]
A, B, _ = equilibrate(ctx.A, ctx.B)
z, p = ctx.z, ctx.params
far = (np.abs(z - p["zc1"]) > FAR * p["sigB"]) & (np.abs(z - p["zc2"]) > FAR * p["sigB"])
print(f"far field = {far.sum()}/{ctx.N} grid points (|z-zc| > {FAR:g} sigma)")
for g in GENERATORS:
    _, blocks = split_defect(ctx, A, B, chain(ctx, g, MODE))
    tot = frac = 0.0
    for blk in blocks:
        for s in range(ctx.n_slots):
            prof = np.fft.ifftn(blk[s * ctx.N : (s + 1) * ctx.N])
            tot += float(np.sum(np.abs(prof) ** 2))
            frac += float(np.sum(np.abs(prof[far]) ** 2))
    print(f"  {g:10s} far-field share of defect energy = {frac / max(tot, 1e-300):.3e}")

print("\n--- (3) xi_z: why is it absent from the measured weak set? ---")
print("relative determination sigma (larger = better determined):")
print("H1 = xi_z too determined via delta a_2 ~ B0 ; H2 = depth-3 (L_2) chain shape")
for g in GENERATORS:
    sig = [
        float(np.linalg.norm(defect(A, B, chain(ctx, g, m))))
        / max(float(np.linalg.norm(np.concatenate(chain(ctx, g, m)))), 1e-300)
        for m in range(ctx.N)
    ]
    print(f"  {g:10s} median={np.median(sig):.3e}  min={np.min(sig):.3e}  max={np.max(sig):.3e}")

print("\n--- (4) REPRESENTATIVE (gate F1-f): full Lie vs covariant photon action ---")
print("delta a_mu = L_xi Abar_mu   vs   delta a_mu = xi^nu Fbar_{nu mu}")
print("(they differ by a U(1) with chi = xi.Abar, which the Lorenz term breaks)")
for g in ("xi_y", "xi_z"):
    for cov in (False, True):
        R = chain(ctx, g, MODE, covariant=cov)
        parts, blocks = split_defect(ctx, A, B, R)
        tot = frac = 0.0
        for blk in blocks:
            for s in range(ctx.n_slots):
                prof = np.fft.ifftn(blk[s * ctx.N : (s + 1) * ctx.N])
                tot += float(np.sum(np.abs(prof) ** 2))
                frac += float(np.sum(np.abs(prof[far]) ** 2))
        tag = "covariant" if cov else "full Lie "
        print(
            f"  {g:5s} {tag}: a-row={parts['a']:.3e}  h-row={parts['h']:.3e}  "
            f"far-field share={frac / max(tot, 1e-300):.3e}"
        )
print(
    "  xi_z is the CONTROL: Abar_z = 0 so the two representatives coincide and the\n"
    "  numbers must match exactly."
)

print("\nGATE F1-b:", "PASS" if status == 0 else "FAIL")
sys.exit(status)
