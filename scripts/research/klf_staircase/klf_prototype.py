"""Stage-1 prototype: Kronecker-like staircase deflation (GH #473).

Two-phase Van Dooren staircase on the row-equilibrated pencil B·ẏ = A·y
(here written M − λN with M := A, N := B):

Phase 1 (left + infinite structure → trailing corner): repeatedly find
rows with no derivative content (left null space of N), column-compress
their algebraic content, move the block to the bottom-right, recurse on
the rest. Trailing columns are algebraically RECOVERED from the leading
state; surplus rows are constraint chains (kept, checked by the full-row
contract — never dropped).

Phase 2 (right structure → leading corner): on the remaining middle,
repeatedly find columns with no derivative content (null of N),
row-compress their equation content; columns with no equations at all
are FREE (gauge) — pinned; determined ones join the leading corner with
their defining rows. What survives has N invertible: the REGULAR part,
whose generator is a solve against N.

Outputs compose exactly like the engine expects: A_eff (evolution on the
manifold), proj (on-manifold projector incl. recovered trailing
columns), pinned basis (for the certificate). Acceptance is the
deflation contract over ALL original rows on manifold states.
"""

from __future__ import annotations

import numpy as np

Arr = np.ndarray


def _null_cols(mat: Arr, tol: float) -> tuple[Arr, Arr]:
    """(null-space basis, complement basis) of mat's columns at tol·smax."""
    if mat.shape[0] == 0:
        n = mat.shape[1]
        return np.eye(n, dtype=complex), np.zeros((n, 0), dtype=complex)
    _u, s, vh = np.linalg.svd(mat, full_matrices=True)
    smax = s[0] if s.size else 0.0
    r = int(np.sum(s > tol * max(smax, 1.0)))
    v = vh.conj().T
    return v[:, r:], v[:, :r]


def _range_rows(mat: Arr, tol: float) -> tuple[Arr, Arr, int]:
    """(range basis U_r, null-complement U_n, rank) of mat's rows."""
    if mat.shape[1] == 0 or mat.shape[0] == 0:
        m = mat.shape[0]
        return np.zeros((m, 0), dtype=complex), np.eye(m, dtype=complex), 0
    u, s, _ = np.linalg.svd(mat, full_matrices=True)
    smax = s[0] if s.size else 0.0
    r = int(np.sum(s > tol * max(smax, 1.0)))
    return u[:, :r], u[:, r:], r


class Staircase:
    """Result of the two-phase reduction, in the transformed basis.

    Column order (of Z): [pinned (right/free) | finite | trailing
    (recovered)]. Row order (of Q): [consumed-right rows | finite rows |
    trailing algebraic rows]. All of Q, Z unitary by construction.
    """

    def __init__(self, A: Arr, B: Arr, tol: float) -> None:
        n = A.shape[0]
        self.n = n
        self.tol = tol
        Q = np.eye(n, dtype=complex)
        Z = np.eye(n, dtype=complex)
        M = A.astype(complex).copy()
        N = B.astype(complex).copy()

        # ---- Phase 1: peel algebraic rows (left + infinite) to the bottom
        r_hi = n  # rows [0, r_hi) are "active"; [r_hi, n) trailing
        c_hi = n  # cols [0, c_hi) active; [c_hi, n) trailing (recovered)
        while True:
            Nw = N[:r_hi, :c_hi]
            M[:r_hi, :c_hi]
            # rows with no derivative content: left null space of Nw
            u_r, u_null, _ = _range_rows(Nw, tol)
            p = u_null.shape[1]
            if p == 0:
                break
            # rotate rows: derivative rows first, algebraic rows last
            Urot = np.hstack([u_r, u_null])
            M[:r_hi, :] = Urot.conj().T @ M[:r_hi, :]
            N[:r_hi, :] = Urot.conj().T @ N[:r_hi, :]
            Q[:, :r_hi] @= Urot
            alg = slice(r_hi - p, r_hi)
            # column-compress the algebraic rows' content: determined cols
            v_free, v_det = _null_cols(M[alg, :c_hi], tol)
            q = v_det.shape[1]
            Vrot = np.hstack([v_free, v_det])  # determined cols last
            M[:, :c_hi] @= Vrot
            N[:, :c_hi] @= Vrot
            Z[:, :c_hi] @= Vrot
            r_hi -= p
            c_hi -= min(q, p)  # at most p columns are determined by p rows
            if q > p:
                # over-determined: q can exceed p only through rank
                # interplay; retreat columns by q anyway is wrong — the
                # extra rows are constraint chains on the leading state.
                c_hi -= q - p
            if p == 0 and q == 0:
                break

        self.r_fin_rows = r_hi
        self.c_lead = c_hi  # columns not consumed by phase 1

        # ---- Phase 2: peel free/right columns to the front
        c_lo = 0
        r_lo = 0
        while True:
            Nw = N[r_lo:r_hi, c_lo:c_hi]
            M[r_lo:r_hi, c_lo:c_hi]
            v_null, v_reg = _null_cols(Nw, tol)
            m = v_null.shape[1]
            if m == 0:
                break
            Vrot = np.hstack([v_null, v_reg])  # N-null cols first
            M[:, c_lo:c_hi] @= Vrot
            N[:, c_lo:c_hi] @= Vrot
            Z[:, c_lo:c_hi] @= Vrot
            head = slice(c_lo, c_lo + m)
            # rows acting on the N-null columns
            u_r, u_rest, rho = _range_rows(M[r_lo:r_hi, head], tol)
            Urot = np.hstack([u_r, u_rest])
            M[r_lo:r_hi, :] = Urot.conj().T @ M[r_lo:r_hi, :]
            N[r_lo:r_hi, :] = Urot.conj().T @ N[r_lo:r_hi, :]
            Q[:, r_lo:r_hi] @= Urot
            # the m columns move to the leading corner with their rho rows;
            # (m - rho) of them are free at this level (gauge)
            c_lo += m
            r_lo += rho
            if m == 0 and rho == 0:
                break

        self.c_pin = c_lo  # leading columns = right structure (incl. free)
        self.r_used = r_lo
        self.M, self.N, self.Q, self.Z = M, N, Q, Z
        # finite (regular) block
        self.fin_rows = slice(r_lo, r_hi)
        self.fin_cols = slice(c_lo, c_hi)
        self.trail_cols = slice(c_hi, n)
        self.trail_rows = slice(r_hi, n)

    def compose(self) -> tuple[Arr, Arr, Arr, dict]:
        """(A_eff, proj, pinned_basis, info) in the ORIGINAL basis.

        Semantics: right-structure columns pinned to zero; finite block
        evolves by G = E⁻¹A_f; trailing columns recovered algebraically
        from the finite state via the trailing rows (least squares, full
        column rank); manifold = span(finite ⊕ recovered-trailing).
        """
        M, N, Z = self.M, self.N, self.Z
        fr, fc = self.fin_rows, self.fin_cols
        tr, tc = self.trail_rows, self.trail_cols
        E = N[fr, fc]
        Af = M[fr, fc]
        nf = E.shape[0]
        info = {
            "n_pinned": self.c_pin,
            "n_finite": nf,
            "n_trail": self.n - (self.c_hi if hasattr(self, "c_hi") else 0),
        }
        # trailing recovery: algebraic rows: M[tr, tc] z_t + M[tr, fc] z_f = 0
        T = M[tr, tc]
        S = M[tr, fc]
        if T.shape[1] > 0:
            R_t = -np.linalg.lstsq(T, S, rcond=None)[0]  # z_t = R_t z_f
        else:
            R_t = np.zeros((0, nf), dtype=complex)
        # finite dynamics: rows fr: N[fr,fc] ż_f + N[fr,tc] ż_t = M[fr,fc] z_f + M[fr,tc] z_t
        # (couplings to pinned cols vanish: z_pin ≡ 0)
        E_eff = E + N[fr, tc] @ R_t
        A_rows = Af + M[fr, tc] @ R_t
        G = np.linalg.solve(E_eff, A_rows)  # ż_f = G z_f
        # state map: y = Z_f z_f + Z_t z_t = (Z_f + Z_t R_t) z_f
        Zf = Z[:, fc]
        Zt = Z[:, tc]
        X = Zf + Zt @ R_t  # (n, nf) manifold parametrization
        # A_eff y = d/dt y = X ż_f = X G z_f ; z_f from y via pseudo-inverse
        Xp = np.linalg.pinv(X)
        A_eff = X @ G @ Xp
        proj = X @ Xp
        pinned = Z[:, : self.c_pin]
        return A_eff, proj, pinned, info


def klf_deflate(A: Arr, B: Arr, tau: float) -> tuple[Arr, Arr, Arr, dict]:
    """Row-equilibrate, staircase, compose. Returns (A_eff, proj, pinned, info)."""
    s = np.maximum(np.max(np.abs(A), axis=1), np.max(np.abs(B), axis=1))
    s[s == 0.0] = 1.0
    An, Bn = A / s[:, None], B / s[:, None]
    st = Staircase(An, Bn, tau)
    st.c_hi = st.fin_cols.stop
    A_eff, proj, pinned, info = st.compose()
    # full-row contract on manifold states (the acceptance)
    rng = np.random.default_rng(0)
    Y = proj @ (
        rng.standard_normal((A.shape[0], 8)) + 1j * rng.standard_normal((A.shape[0], 8))
    )
    denom = max(np.linalg.norm(An @ Y), 1e-300)
    info["contract"] = float(np.linalg.norm(Bn @ (A_eff @ Y) - An @ Y) / denom)
    info["max_re"] = float(np.max(np.linalg.eigvals(A_eff).real)) if A.shape[0] else 0.0
    return A_eff, proj, pinned, info


# ---------------------------------------------------------------------------
# Gate (a): exact synthetics
# ---------------------------------------------------------------------------


def gate_a() -> None:
    ok = True

    # (a1) invertible B → A_eff == B^-1 A, proj == I, nothing pinned
    rng = np.random.default_rng(3)
    A = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    B = np.diag(np.array([2.0, -0.5, 1.5, 1.0], dtype=complex))
    A_eff, proj, pinned, info = klf_deflate(A, B, 1e-10)
    e1 = np.linalg.norm(A_eff - np.linalg.solve(B, A))
    print(
        f"a1 invertible-B: err={e1:.1e} pinned={pinned.shape[1]} contract={info['contract']:.1e}"
    )
    ok &= e1 < 1e-10 and pinned.shape[1] == 0

    # (a2) L1 orphan pair (qφ,vφ,qc,vc): row for vc missing → (qc,vc) chain free
    A = np.zeros((4, 4), dtype=complex)
    B = np.zeros((4, 4), dtype=complex)
    B[0, 0] = 1.0
    A[0, 1] = 1.0
    B[1, 1] = 1.0
    A[1, 0] = -1.0
    B[2, 2] = 1.0
    A[2, 3] = 1.0
    A_eff, proj, pinned, info = klf_deflate(A, B, 1e-10)
    y = np.array([0.7, -0.2, 0.0, 0.0], dtype=complex)
    dy = A_eff @ y
    e2 = abs(dy[0] + 0.2) + abs(dy[1] + 0.7)
    row_pin = np.sum(np.abs(pinned) ** 2, axis=1)
    print(
        f"a2 orphan pair: osc-err={e2:.1e} pinned={pinned.shape[1]} pin-rows={np.round(row_pin, 3)} contract={info['contract']:.1e}"
    )
    ok &= (
        e2 < 1e-9 and pinned.shape[1] == 2 and row_pin[0] < 1e-12 and row_pin[1] < 1e-12
    )

    # (a3) #465 duplicate rows: same equation twice → chain q_diff + λ v_diff free
    #     q̇1 = v1 ; v̇1-row duplicated; fields (q1,v1,q2,v2) with q2 an exact copy's partner
    #     Simplest faithful synthetic: two identical algebraic rows on (q1,q2):
    #     rows: q̇1=v1 ; v̇1=-q1 ; alg: q1 - q2 = 0 ; alg dup: q1 - q2 = 0
    A = np.zeros((4, 4), dtype=complex)
    B = np.zeros((4, 4), dtype=complex)
    B[0, 0] = 1.0
    A[0, 1] = 1.0
    B[1, 1] = 1.0
    A[1, 0] = -1.0
    A[2, 0] = 1.0
    A[2, 2] = -1.0
    A[3, 0] = 1.0
    A[3, 2] = -1.0
    # column 3 (v2) appears nowhere → free; q2 recovered from q1
    A_eff, proj, pinned, info = klf_deflate(A, B, 1e-10)
    y = proj @ np.array([0.5, 0.1, 0.5, 0.0], dtype=complex)
    dy = A_eff @ y
    e3 = abs(dy[0] - y[1]) + abs(dy[1] + y[0]) + abs(dy[2] - dy[0])
    print(
        f"a3 dup rows: err={e3:.1e} pinned={pinned.shape[1]} contract={info['contract']:.1e} maxRe={info['max_re']:+.1e}"
    )
    ok &= e3 < 1e-9 and info["contract"] < 1e-9

    # (a4) index-2 chain (hand-solved, from test_pencil_engine)
    w2, a, b, d, g = 1.7, 0.9, 0.35, -0.6, 0.4j
    A = np.zeros((4, 4), dtype=complex)
    B = np.zeros((4, 4), dtype=complex)
    B[0, 0] = 1.0
    A[0, 1] = 1.0
    B[1, 1] = 1.0
    A[1, 0] = -w2
    A[1, 3] = g
    B[2, 2] = 1.0
    A[2, 3] = 1.0
    A[3, 2] = a
    A[3, 1] = b
    A[3, 0] = d
    A_eff, proj, pinned, info = klf_deflate(A, B, 1e-10)
    qphi, vphi = 0.8, -0.3
    qc = -(d * qphi + b * vphi) / a
    vc = (b * w2 * qphi - d * vphi) / (a + g * b)
    y = np.array([qphi, vphi, qc, vc], dtype=complex)
    dy = A_eff @ y
    vdot = (-w2 * a / (a + g * b)) * qphi + (-g * d / (a + g * b)) * vphi
    e4 = abs(dy[0] - vphi) + abs(dy[1] - vdot)
    e4p = np.linalg.norm(proj @ y - y)
    print(
        f"a4 index-2: dyn-err={e4:.1e} proj-err={e4p:.1e} pinned={pinned.shape[1]} contract={info['contract']:.1e}"
    )
    ok &= e4 < 1e-8 and e4p < 1e-8

    # (a5) left chain: constraint + its time derivative (over-determined, consistent)
    #     q̇=v ; v̇=-q ; alg: q + s = 0 (defines s) ; alg2: v + w = 0 (defines w) ;
    #     extra row: derivative consistency s + ... keep simple: fields (q,v,s)
    #     rows: q̇=v ; v̇=-q ; 0 = q + s ; 0 = v + ṡ  (ṡ has B-content: left chain)
    A = np.zeros((4, 3), dtype=complex)  # rectangular not supported: embed square
    A = np.zeros((4, 4), dtype=complex)
    B = np.zeros((4, 4), dtype=complex)
    B[0, 0] = 1.0
    A[0, 1] = 1.0
    B[1, 1] = 1.0
    A[1, 0] = -1.0
    A[2, 0] = 1.0
    A[2, 2] = 1.0  # 0 = q + s
    B[3, 2] = -1.0
    A[3, 1] = 1.0  # ṡ = v  (consistent with s = -q)
    # column 3 unused → free direction; s recovered; the two s-rows consistent
    A_eff, proj, pinned, info = klf_deflate(A, B, 1e-10)
    y = proj @ np.array([0.4, -0.7, -0.4, 0.0], dtype=complex)
    dy = A_eff @ y
    e5 = abs(dy[0] - y[1]) + abs(dy[1] + y[0]) + abs(dy[2] + dy[0])
    print(
        f"a5 left chain: err={e5:.1e} pinned={pinned.shape[1]} contract={info['contract']:.1e} maxRe={info['max_re']:+.1e}"
    )
    ok &= e5 < 1e-8 and info["contract"] < 1e-8

    print("GATE A:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    gate_a()
