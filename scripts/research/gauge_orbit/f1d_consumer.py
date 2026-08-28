"""F1-d: the symbolic-orbit quotient, end to end (the arc's CHECKPOINT gate).

Construction (arc-F plan steps 3-5, CORRECTED 2026-08-26 by measurement):

  * grade the JOINT chain span; pin profiles with sigma < tau_pin*||(A,B)||;
  * the reduced dynamics lives on the QUOTIENT, not on a slice.  Given a
    state y in the gauge slice, the true ydot may DRIFT along the pinned
    span W; that drift is projected out.  So the generator solves

        B (Cb xdot + W w) = A y          for (xdot, w)

    and keeps xdot.  Forbidding the drift instead (demanding ydot itself
    lie in the slice) over-constrains the system and fails at O(0.2) --
    measured, and it fails EVEN FOR AN EXACT SYMMETRY (Bpeak=0), which is
    how the error was caught.  No equation is ever dropped, so the GH #474
    information-discarding defect is structurally impossible here.

Gates: contract on the constraint manifold; the {h_5,a_1} physical block
exact; cross-coupling into physics; and max Re lambda N- AND tau_pin-
independent (the anti-pincer gate -- GH #473 failed exactly there).

usage: f1d_consumer.py [N ...] [--tau=1e-6,1e-5,1e-4]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

sys.path.insert(0, str(Path(__file__).parent))
from orbit import capture, equilibrate, grade, joint_operator  # noqa: E402

OBS = ("h_5", "a_1")


def build(ctx, A, B, tau_pin: float):
    """Symbolic-orbit quotient -> generator on the slice-manifold."""
    Cmap, Dmap, _ = joint_operator(ctx, A, B)
    sig, _, chains = grade(Cmap, Dmap)
    n = ctx.dim
    scale = max(np.linalg.norm(A, 1), np.linalg.norm(B, 1))
    depth = Cmap.shape[0] // n - 1
    pinned = np.where(sig < tau_pin * scale)[0]
    if len(pinned) == 0:
        return None

    cols = []
    for i in pinned:
        for j in range(depth + 1):
            v = chains[j * n : (j + 1) * n, i]
            nv = np.linalg.norm(v)
            if nv > 1e-12:
                cols.append(v / nv)
    U, s, _ = np.linalg.svd(np.array(cols).T, full_matrices=False)
    W = U[:, s > 1e-10 * s[0]]
    p = W.shape[1]
    Uw, _, _ = np.linalg.svd(W, full_matrices=True)
    Cb = Uw[:, p:]

    # algebraic constraints: rows with no B content demand (A y)_alg = 0
    Ub, sb, _ = np.linalg.svd(B)
    rk = int((sb > 1e-10 * sb[0]).sum())
    Malg = Ub[:, rk:].conj().T @ A
    slice_man = null_space(np.vstack([Malg, W.conj().T]), rcond=1e-10)
    if slice_man.shape[1] == 0:
        return None

    # generator: solve B(Cb xdot + W w) = A y, keep xdot, express in slice coords
    M = np.hstack([B @ Cb, B @ W])
    rhs = A @ slice_man
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    ydot = Cb @ sol[: Cb.shape[1], :]
    G = slice_man.conj().T @ ydot
    lsq_resid = float(
        np.linalg.norm(M @ sol - rhs) / max(np.linalg.norm(rhs), 1e-300)
    )
    return {
        "G": G, "S": slice_man, "W": W, "p": p, "sigma": sig,
        "n_profiles": len(pinned), "lsq_resid": lsq_resid,
        "sigma_pinned_max": float(sig[pinned].max()),
        "sigma_first_retained": float(sig[len(pinned)]) if len(pinned) < len(sig) else float("nan"),
    }


def run(n: int, taus: list[float]) -> int:
    ctx = capture(n)
    A, B, _ = equilibrate(ctx.A, ctx.B)
    mask = np.zeros(ctx.dim, dtype=bool)
    for f in OBS:
        for nm in (f, "v_" + f):
            s = ctx.slot[nm]
            mask[s * ctx.N : (s + 1) * ctx.N] = True
    rng = np.random.default_rng(0)
    print(f"\nN={n}  dim={ctx.dim}")
    status = 0
    for tau in taus:
        out = build(ctx, A, B, tau)
        if out is None:
            print(f"  tau={tau:.0e}: nothing pinned / empty slice")
            continue
        G, S, W, p = out["G"], out["S"], out["W"], out["p"]
        x = rng.standard_normal((S.shape[1], 8)) + 1j * rng.standard_normal((S.shape[1], 8))
        y, ydot = S @ x, S @ (G @ x)
        contract = float(
            np.linalg.norm(B @ ydot - A @ y) / max(np.linalg.norm(A @ y), 1e-300)
        )
        A_eff = S @ G @ S.conj().T
        blk = float(
            np.linalg.norm(A_eff[np.ix_(mask, mask)] - A[np.ix_(mask, mask)])
            / np.linalg.norm(A[np.ix_(mask, mask)])
        )
        cross = max(
            float(np.linalg.norm(A_eff[np.ix_(mask, ~mask)])),
            float(np.linalg.norm(A_eff[np.ix_(~mask, mask)])),
        )
        ev = np.linalg.eigvals(G)
        pin_obs = float(np.linalg.norm(W[mask, :]))
        ok = contract < 1e-10
        status |= 0 if ok else 1
        print(
            f"  tau={tau:.0e}: p={p:4d} slice={S.shape[1]:4d} contract={contract:.2e} "
            f"obs-blk={blk:.2e} cross={cross:.2e} maxRe={np.max(ev.real):+.3e} "
            f"|Im|max={np.max(np.abs(ev.imag)):.2f} pin-on-obs={pin_obs:.1e}"
        )
    return status


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tau_arg = next((a for a in sys.argv[1:] if a.startswith("--tau")), None)
    taus = [float(x) for x in tau_arg.split("=")[1].split(",")] if tau_arg else [1e-6, 1e-5, 1e-4]
    rc = 0
    for nn in ([int(a) for a in args] or [16, 24]):
        rc |= run(nn, taus)
    print("\nGATE F1-d:", "PASS" if rc == 0 else "FAIL")
    sys.exit(rc)
