"""Faithful numpy port of the GUPTRI-style KLF staircase (GH #473 stage 1).

Ported from MatrixPencils.jl (Andreas Varga): `_preduceBF!`, `_preduce2!`,
`_preduce3!`, `_preduce4!` (SVD variants) and the `klf_right!` driver —
SVD rank decisions, Givens choreography preserving E's triangular form,
orthogonal accumulation in Q, Z. Convention: pencil M − λN over
y-columns; TIDAL's B·ẏ = A·y enters as M := A, N := B.

Layout after `klf_right` (finite_infinite = true path of `klf`):

        [ Mr-λNr |   *    |    *     ]   rows: sum(νr) | nf | sum(ν)
    Q'(M-λN)Z =  [   0    | Mf-λNf |  *  ]   cols: sum(μr) | nf | sum(μ)
        [   0    |   0    | Mil-λNil ]

Mr: right (gauge) staircase, full row rank, more columns than rows;
Mf − λNf: regular part, Nf invertible; Mil: infinite + left staircase,
full column rank → its states vanish identically for the homogeneous
problem. Composition: z_trail ≡ 0; the right rows are solved min-norm
for their determined part (nilpotent fixed point over the chain depth);
the free part is pinned (gauge, certified downstream).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import qr as _qr

Arr = np.ndarray
C128 = np.complex128


def _giv(f: complex, g: complex) -> tuple[float, complex, complex]:
    """Givens (c real, s complex, r): c*f + s*g = r, -conj(s)*f + c*g = 0."""
    if g == 0:
        return 1.0, 0.0 + 0.0j, f
    if f == 0:
        d = abs(g)
        return 0.0, np.conj(g) / d, d
    d = np.hypot(abs(f), abs(g))
    fs = f / abs(f)
    c = abs(f) / d
    s = fs * np.conj(g) / d
    return float(c), complex(s), complex(fs * d)


def _rot_rows(
    mats: list[tuple[Arr, slice]], i1: int, i2: int, c: float, s: complex
) -> None:
    """Apply G to rows i1, i2 (global indices) of each (mat, col-slice)."""
    for m, js in mats:
        r1 = m[i1, js].copy()
        r2 = m[i2, js].copy()
        m[i1, js] = c * r1 + s * r2
        m[i2, js] = -np.conj(s) * r1 + c * r2


def _rot_cols(
    mats: list[tuple[Arr, slice]], j1: int, j2: int, c: float, s: complex
) -> None:
    """Apply G' on the right to columns j1, j2 of each (mat, row-slice)."""
    for m, isl in mats:
        c1 = m[isl, j1].copy()
        c2 = m[isl, j2].copy()
        m[isl, j1] = c * c1 + np.conj(s) * c2
        m[isl, j2] = -s * c1 + c * c2


def preduce_bf(
    M: Arr,
    N: Arr,
    Q: Arr,
    Z: Arr,
    tol_n: float,
    roff: int = 0,
    coff: int = 0,
    rtrail: int = 0,
    ctrail: int = 0,
) -> tuple[int, int, int]:
    """SVD variant of `_preduceBF!`: N22 → [0 | E] with E = diag, nonsingular."""
    mM, nM = M.shape
    npp = mM - rtrail - roff
    npm = nM - ctrail - coff
    if npp == 0 or npm == 0:
        return 0, npm, npp
    i22 = slice(roff, roff + npp)
    j22 = slice(coff, coff + npm)
    i12 = slice(0, roff + npp)
    j23 = slice(coff, nM)
    jN23 = slice(coff + npm, nM)
    U, S, Vt = np.linalg.svd(N[i22, j22])
    n = int(np.sum(tol_n < S)) if S.size else 0  # absolute tau (equilibrated rows)
    m = npm - n
    p = npp - n
    V = Vt.conj().T
    # column permutation: null cols first, range cols last
    perm = np.concatenate([np.arange(n, npm), np.arange(0, n)])
    Vp = V[:, perm]
    M[i12, j22] @= Vp
    N[0:roff, j22] @= Vp
    M[i22, j23] = U.conj().T @ M[i22, j23]
    N[i22, jN23] = U.conj().T @ N[i22, jN23]
    Q[:, i22] @= U
    Z[:, j22] @= Vp
    block = np.zeros((npp, npm), dtype=C128)
    block[:n, m:] = np.diag(S[:n])
    N[i22, j22] = block
    return n, m, p


def preduce4(
    n: int,
    m: int,
    p: int,
    M: Arr,
    N: Arr,
    Q: Arr,
    Z: Arr,
    tol: float,
    roff: int,
    coff: int,
    rtrail: int,
    ctrail: int,
) -> int:
    """SVD variant of `_preduce4!`: column-compress C, keep E triangular."""
    n + m
    M.shape[0]
    nM = M.shape[1]
    ie0 = roff
    ic0 = roff + n
    jc0 = coff + m
    if n > p:
        for i in range(1, p + 1):
            ii = ic0 + (p - i)  # global row of C being swept
            for jj in range(n - i):  # local col indices in C
                j1g = jc0 + jj
                j2g = jc0 + jj + 1
                cg, sg, r = _giv(np.conj(M[ii, j2g]), np.conj(M[ii, j1g]))
                # right rotation on columns (j2g, j1g)
                _rot_cols(
                    [(M, slice(0, ii)), (N, slice(0, roff + jj + 2))], j2g, j1g, cg, sg
                )
                _rot_cols([(Z, slice(0, Z.shape[0]))], j2g, j1g, cg, sg)
                M[ii, j2g] = np.conj(r)
                M[ii, j1g] = 0.0
                # restore E triangular: row rotation (jj, jj+1) within E rows
                i1g = ie0 + jj
                i2g = ie0 + jj + 1
                cg2, sg2, r2 = _giv(N[i1g, j1g], N[i2g, j1g])
                _rot_rows([(N, slice(j1g + 1, nM))], i1g, i2g, cg2, sg2)
                _rot_rows([(M, slice(coff, nM))], i1g, i2g, cg2, sg2)
                N[i1g, j1g] = r2
                N[i2g, j1g] = 0.0
                _rot_cols([(Q, slice(0, Q.shape[0]))], i1g, i2g, cg2, sg2)
    pn = min(n, p)
    if pn == 0:
        return 0
    jcs = slice(jc0 + n - pn, jc0 + n)
    Cblk = M[ic0 : ic0 + p, jcs]
    U, S, Vt = np.linalg.svd(Cblk)
    rho = int(np.sum(tol < S))
    if rho == pn:
        return rho  # C already full column rank — leave untouched (ref line ~1365)
    Q1 = U[:, ::-1]
    newC = np.zeros((p, pn), dtype=C128)
    if rho > 0:
        newC[:, pn - rho :] = Q1[:, p - rho :] @ np.diag(S[:rho][::-1])
    M[ic0 : ic0 + p, jcs] = newC
    if rho == 0:
        return rho
    Z1 = Vt.conj().T[:, ::-1]
    M[0 : roff + n, jcs] @= Z1
    N[0 : roff + n, jcs] @= Z1
    Z[:, jcs] @= Z1
    # re-triangularize E22 = N[it, jcs]
    it = slice(roff + n - pn, roff + n)
    E22 = N[it, jcs]
    Qf, Rf = np.linalg.qr(E22)
    N[it, jcs] = Rf
    jt1 = slice(jc0 + n, nM)
    N[it, jt1] = Qf.conj().T @ N[it, jt1]
    M[it, coff:nM] = Qf.conj().T @ M[it, coff:nM]
    Q[:, it] @= Qf
    return rho


def preduce2(
    n: int,
    m: int,
    p: int,
    M: Arr,
    N: Arr,
    Q: Arr,
    Z: Arr,
    tol: float,
    roff: int,
    coff: int,
    rtrail: int,
    ctrail: int,
) -> tuple[int, int]:
    """SVD variant of `_preduce2!`: compress [D C] columns; returns (tau, rho)."""
    npp = n + p
    npm = n + m
    M.shape[1]
    ia = slice(0, roff + n)
    ic = slice(roff + n, roff + npp)
    jb = slice(coff, coff + m)
    if m > 0:
        D = M[ic, jb]
        tau = (
            int(np.sum(np.linalg.svd(D, compute_uv=False) > tol)) if min(D.shape) else 0
        )
        Qf, Rf, piv = _qr(D.conj().T, mode="full", pivoting=True)
        # BE ← BE·Qf, then reverse columns (of the m-block)
        M[ia, jb] = (M[ia, jb] @ Qf)[:, ::-1]
        N[0:roff, jb] = (N[0:roff, jb] @ Qf)[:, ::-1]
        Z[:, jb] = (Z[:, jb] @ Qf)[:, ::-1]
        newD = np.zeros((p, m), dtype=C128)
        if tau > 0:
            newD[:, m - tau :] = Rf[:tau, ::-1][::-1, :].conj().T
        M[ic, jb] = newD
        # C rows permuted by piv then reversed
        jc = slice(coff + m, coff + npm)
        M[ic, jc] = M[ic, jc][piv, :][::-1, :]
        # Step 2: Givens sweep zeroing the sub-antidiagonal of the tau rows
        k = 1
        for ig in range(roff + npp - 1, roff + npp - 1 - tau, -1):
            for jj in range(coff + m - k, coff + m + n - k):
                j2 = jj + 1
                cg, sg, r = _giv(np.conj(M[ig, j2]), np.conj(M[ig, jj]))
                _rot_cols([(M, slice(0, ig)), (N, slice(0, ig + 1))], j2, jj, cg, sg)
                _rot_cols([(Z, slice(0, Z.shape[0]))], j2, jj, cg, sg)
                M[ig, j2] = np.conj(r)
                M[ig, jj] = 0.0
            k += 1
    else:
        tau = 0
        piv = np.arange(0)
    rho = preduce4(
        n, m - tau, p - tau, M, N, Q, Z, tol, roff, coff, rtrail + tau, ctrail + tau
    )
    if m > 0:
        jrt = slice(coff + npm - (tau + rho), coff + npm)
        blk = M[ic, jrt][::-1, :]
        inv = np.empty_like(piv)
        inv[piv] = np.arange(len(piv))
        M[ic, jrt] = blk[inv, :]
    return tau, rho


def preduce3(
    n: int,
    m: int,
    M: Arr,
    N: Arr,
    Q: Arr,
    Z: Arr,
    tol: float,
    roff: int,
    coff: int,
    rtrail: int,
    ctrail: int,
) -> int:
    """SVD variant of `_preduce3!`: row-compress B, keep E triangular."""
    nM = M.shape[1]
    ib0 = roff
    je0 = coff + m
    if n > m:
        for j in range(m):
            jg = coff + j
            for iiL in range(n - 1, j, -1):
                i1 = ib0 + iiL - 1
                i2 = ib0 + iiL
                cg, sg, r = _giv(M[i1, jg], M[i2, jg])
                _rot_rows(
                    [(M, slice(jg + 1, nM)), (N, slice(je0 + iiL - 1, nM))],
                    i1,
                    i2,
                    cg,
                    sg,
                )
                M[i1, jg] = r
                M[i2, jg] = 0.0
                _rot_cols([(Q, slice(0, Q.shape[0]))], i1, i2, cg, sg)
                # restore E: column rotation on (je0+iiL, je0+iiL-1)
                j1 = je0 + iiL
                j2 = je0 + iiL - 1
                cg2, sg2, r2 = _giv(np.conj(N[i2, j1]), np.conj(N[i2, j2]))
                _rot_cols(
                    [(N, slice(0, i2)), (M, slice(0, M.shape[0]))], j1, j2, cg2, sg2
                )
                _rot_cols([(Z, slice(0, Z.shape[0]))], j1, j2, cg2, sg2)
                N[i2, j1] = np.conj(r2)
                N[i2, j2] = 0.0
    mn = min(n, m)
    if mn == 0:
        return 0
    B = M[ib0 : ib0 + mn, coff : coff + m]
    U, S, Vt = np.linalg.svd(B)
    rho = int(np.sum(tol < S))
    newB = np.zeros((mn, m), dtype=C128)
    if rho > 0:
        newB[:rho, :] = np.diag(S[:rho]) @ Vt[:rho, :]
    if rho == mn:
        return rho
    M[ib0 : ib0 + mn, coff : coff + m] = newB
    if rho == 0:
        return 0
    ibt = slice(ib0, ib0 + mn)
    jt = slice(je0, nM)
    Q[:, ibt] @= U
    N[ibt, jt] = U.conj().T @ N[ibt, jt]
    M[ibt, jt] = U.conj().T @ M[ibt, jt]
    # re-triangularize E11 = N[ibt, je0:je0+mn] via reversed-QR
    jt1 = slice(je0, je0 + mn)
    E11 = N[ibt, jt1]
    Qf, _ = np.linalg.qr(E11[::-1, :].conj().T)
    for mat, isl in (
        (M, slice(0, M.shape[0])),
        (N, slice(0, ib0 + mn)),
        (Z, slice(0, Z.shape[0])),
    ):
        mat[isl, jt1] = (mat[isl, jt1] @ Qf)[:, ::-1]
    E11v = N[ibt, jt1]
    E11v[:] = np.triu(E11v)
    return rho


def _invariants(
    M: Arr,
    N: Arr,
    Q: Arr,
    Z: Arr,
    A0: Arr,
    B0: Arr,
    roff: int,
    coff: int,
    rtrail: int,
    ctrail: int,
    label: str,
) -> None:
    mM, nM = M.shape
    eq_m = np.linalg.norm(Q.conj().T @ A0 @ Z - M)
    eq_n = np.linalg.norm(Q.conj().T @ B0 @ Z - N)
    tr_m = np.max(np.abs(M[mM - rtrail :, : nM - ctrail])) if rtrail else 0.0
    tr_n = np.max(np.abs(N[mM - rtrail :, : nM - ctrail])) if rtrail else 0.0
    ld_m = np.max(np.abs(M[roff:, :coff])) if coff else 0.0
    ld_n = np.max(np.abs(N[roff:, :coff])) if coff else 0.0
    nz_n = np.max(np.abs(N[mM - rtrail :, nM - ctrail :])) if rtrail and ctrail else 0.0
    print(
        f"    [{label}] equiv M/N: {eq_m:.1e}/{eq_n:.1e}  "
        f"trail-zero M/N: {tr_m:.1e}/{tr_n:.1e}  lead-zero M/N: {ld_m:.1e}/{ld_n:.1e}  "
        f"N-trail-block: {nz_n:.1e}"
    )


def klf_right(A: Arr, B: Arr, tol_rel: float, debug: bool = False) -> dict:
    """PreduceBF + klf_right! driver. Returns transformed pencil + layout."""
    nsz = A.shape[0]
    M = A.astype(C128).copy()
    N = B.astype(C128).copy()
    Q = np.eye(nsz, dtype=C128)
    Z = np.eye(nsz, dtype=C128)
    tol_n = tol_rel
    n, m, p = preduce_bf(M, N, Q, Z, tol_n)
    if debug:
        print(f"  preduce_bf: n={n} m={m} p={p}")
        _invariants(M, N, Q, Z, A, B, 0, 0, 0, 0, "bf")
    # Absolute tolerance on the row-equilibrated pencil (rows O(1)): the
    # determination-floor semantics of the existing engine. The reference's
    # rtol*opnorm(M,1) scaling would inflate tau by ~norm1(M) and truncate
    # genuine B0^2-scale determination content.
    tol = tol_rel
    roff = coff = rtrail = ctrail = 0
    nu, mu = [], []
    while p > 0:
        tau, rho = preduce2(n, m, p, M, N, Q, Z, tol, roff, coff, rtrail, ctrail)
        nu.append(p)
        mu.append(rho + tau)
        ctrail += tau + rho
        rtrail += p
        n -= rho
        p = rho
        m -= tau
        if debug:
            print(
                f"  preduce2: tau={tau} rho={rho} -> n={n} m={m} p={p} rtrail={rtrail} ctrail={ctrail}"
            )
            _invariants(M, N, Q, Z, A, B, roff, coff, rtrail, ctrail, "p2")
    nur, mur = [], []
    while m > 0:
        rho = preduce3(n, m, M, N, Q, Z, tol, roff, coff, rtrail, ctrail)
        nur.append(rho)
        mur.append(m)
        roff += rho
        coff += m
        n -= rho
        m = rho
        if debug:
            print(f"  preduce3: rho={rho} -> n={n} m={m} roff={roff} coff={coff}")
            _invariants(M, N, Q, Z, A, B, roff, coff, rtrail, ctrail, "p3")
    return {
        "M": M,
        "N": N,
        "Q": Q,
        "Z": Z,
        "roff": roff,
        "coff": coff,
        "nf": n,
        "rtrail": rtrail,
        "ctrail": ctrail,
        "nur": nur,
        "mur": mur,
        "nu": nu,
        "mu": mu,
        "tol": tol,
    }


def compose(dec: dict, tol: float) -> tuple[Arr, Arr, Arr, dict]:
    """(A_eff, proj, pinned, info) from the klf_right layout."""
    M, N, Z = dec["M"], dec["N"], dec["Z"]
    nsz = M.shape[0]
    roff, coff, nf = dec["roff"], dec["coff"], dec["nf"]
    fr = slice(roff, roff + nf)
    fc = slice(coff, coff + nf)
    Ef = N[fr, fc]
    Af = M[fr, fc]
    G = np.linalg.solve(Ef, Af) if nf else np.zeros((0, 0), dtype=C128)
    # right rows: solve min-norm for the determined part of z_r
    Rr = slice(0, roff)
    Mr = M[Rr, 0:coff]
    Nr = N[Rr, 0:coff]
    CM = M[Rr, fc]
    CN = N[Rr, fc]
    if roff > 0 and nf > 0:
        P = np.linalg.pinv(Mr, rcond=1e-10)
        X = -P @ (CM + CN @ G)
        for _ in range(60):
            Xn = -P @ (CM + CN @ G - Nr @ X @ G)
            if np.linalg.norm(Xn - X) <= 1e-14 * max(np.linalg.norm(Xn), 1.0):
                X = Xn
                break
            X = Xn
    else:
        X = np.zeros((coff, nf), dtype=C128)
    W = Z[:, fc] + Z[:, 0:coff] @ X  # manifold parametrization
    Wp = np.linalg.pinv(W)
    A_eff = W @ G @ Wp
    proj = W @ Wp
    # pinned = span of the right block's polynomial null family: the gauge
    # freedom of the underdetermined rows. Union of null(Mr − λNr) at two
    # generic probes (chains of degree ≤ 1 exact; higher degrees enter via
    # the union's span in practice — validated by the gates).
    if coff > 0:
        if roff > 0:
            vecs = []
            for lam0 in (0.5488 + 0.7152j, -1.3113 + 0.4227j):
                pen = Mr - lam0 * Nr
                _u, s, vh = np.linalg.svd(pen)
                smax = s[0] if s.size else 1.0
                r = int(np.sum(s > tol * max(smax, 1.0)))
                vecs.append(vh[r:].conj().T)
            stack = np.hstack(vecs) if vecs else np.zeros((coff, 0), dtype=C128)
            u2, s2, _ = np.linalg.svd(stack, full_matrices=False)
            kdim = int(np.sum(s2 > 1e-10 * max(s2[0] if s2.size else 1.0, 1.0)))
            ker = u2[:, :kdim]
        else:
            ker = np.eye(coff, dtype=C128)
        pinned = Z[:, 0:coff] @ ker
    else:
        pinned = np.zeros((nsz, 0), dtype=C128)
    info = {
        "n_pinned": pinned.shape[1],
        "n_finite": nf,
        "n_right_cols": coff,
        "n_trail_cols": dec["ctrail"],
    }
    return A_eff, proj, pinned, info


def klf_deflate(A: Arr, B: Arr, tau: float) -> tuple[Arr, Arr, Arr, dict]:
    s = np.maximum(np.max(np.abs(A), axis=1), np.max(np.abs(B), axis=1))
    s[s == 0.0] = 1.0
    An, Bn = A / s[:, None], B / s[:, None]
    dec = klf_right(An, Bn, tau)
    A_eff, proj, pinned, info = compose(dec, tau)
    rng = np.random.default_rng(0)
    Y = proj @ (
        rng.standard_normal((A.shape[0], 8)) + 1j * rng.standard_normal((A.shape[0], 8))
    )
    denom = max(np.linalg.norm(An @ Y), 1e-300)
    info["contract"] = float(np.linalg.norm(Bn @ (A_eff @ Y) - An @ Y) / denom)
    ev = np.linalg.eigvals(A_eff)
    info["max_re"] = float(np.max(ev.real)) if ev.size else 0.0
    return A_eff, proj, pinned, info
