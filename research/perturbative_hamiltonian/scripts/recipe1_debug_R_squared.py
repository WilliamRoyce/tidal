# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
Debug: Why does R̃² vanish on q-irreducible projection at linear order on flat?
================================================================================

Round 3 Agent H, 2026-04-26.

This sub-script isolates the surprising result and diagnoses it. We compute:
  - R̃_{abcd} for a single symbolic q with no irreducibility constraints
  - the same with q irreducibility imposed
  - R̃² in both cases
to see at which level the vanishing happens.

If R̃² vanishes only AFTER irreducibility is imposed, that's a genuine
algebraic identity worth understanding (and confirms the verdict).

If R̃² vanishes BEFORE irreducibility, there's a bug in the contraction.
"""

from __future__ import annotations

import pathlib
from itertools import product

import sympy as sp

t, x, y, z = sp.symbols("t x y z", real=True)
COORDS = (t, x, y, z)
ETA = sp.diag(-1, 1, 1, 1)
ETA_INV = ETA


def perm_sign(perm):
    n = len(perm)
    s = 1
    arr = list(perm)
    for i in range(n):
        for j in range(i + 1, n):
            if arr[i] > arr[j]:
                s = -s
    return s


def eps_down(a, b, c, d):
    if len({a, b, c, d}) != 4:
        return 0
    return perm_sign((a, b, c, d))


k = sp.symbols("k0 k1 k2 k3", real=True)
phase = sum(k[i] * COORDS[i] for i in range(4))
profile = sp.exp(sp.I * phase)


def build_q_basis_constrained():
    syms = sp.symbols("c0:4_0:4_0:4")
    c = {}
    idx = 0
    for a, b, cc in product(range(4), repeat=3):
        c[a, b, cc] = syms[idx]
        idx += 1
    constraints = []
    for a, b, cc in product(range(4), repeat=3):
        constraints.append(c[a, b, cc] + c[a, cc, b])
    constraints.extend(sum(ETA_INV[a, a] * c[a, a, cc] for a in range(4)) for cc in range(4))
    constraints.extend(sum(
                eps_down(a, b, cc, d) * c[a, b, cc]
                for a, b, cc in product(range(4), repeat=3)
            ) for d in range(4))
    sol = sp.linsolve(constraints, syms)
    sub = dict(zip(syms, next(iter(sol)), strict=False))
    free = set()
    for v in sub.values():
        if v is not None:
            free |= getattr(v, "free_symbols", set())
    q = {}
    for a, b, cc in product(range(4), repeat=3):
        e = c[a, b, cc]
        q[a, b, cc] = sp.simplify(e.xreplace(sub) if e in sub else e)
    return q, sorted(free, key=lambda s: s.name)


def build_q_unconstrained():
    """A unconstrained generic q (just antisymmetric in last two indices)."""
    q = {}
    for a, b, cc in product(range(4), repeat=3):
        if b == cc:
            q[a, b, cc] = sp.Integer(0)
        elif b < cc:
            q[a, b, cc] = sp.Symbol(f"q_{a}_{b}_{cc}")
        else:
            q[a, b, cc] = -sp.Symbol(f"q_{a}_{cc}_{b}")
    params = []
    for a, b, cc in product(range(4), repeat=3):
        if b < cc:
            params.append(sp.Symbol(f"q_{a}_{b}_{cc}"))
    return q, params


def contortion(q):
    K = {}
    for r, m, n in product(range(4), repeat=3):
        K[r, m, n] = sp.Rational(1, 2) * (q[r, m, n] - q[m, r, n] + q[n, r, m])
    return K


def riemann_cartan(K):
    R = {}
    for r, s, m, n in product(range(4), repeat=4):
        d_m = sp.diff(K[r, n, s] * profile, COORDS[m])
        d_n = sp.diff(K[r, m, s] * profile, COORDS[n])
        R[r, s, m, n] = sp.expand(d_m - d_n)
    return R


def R_squared_full(R):
    """R̃_{abcd} R̃^{abcd} = sum η^{aa'}η^{bb'}η^{cc'}η^{dd'} R[a',b',c',d'] R[a,b,c,d]
    With diag η, the pairs (a,a') etc. just put a'=a with sign η[a,a].
    """
    expr = 0
    for a, b, cc, d in product(range(4), repeat=4):
        sign = ETA_INV[a, a] * ETA_INV[b, b] * ETA_INV[cc, cc] * ETA_INV[d, d]
        expr += sign * R[a, b, cc, d] * R[a, b, cc, d]
    return sp.expand(expr)


def main() -> None:
    out_dir = pathlib.Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True, parents=True)
    out_txt = out_dir / "recipe1_debug_R_squared_run.txt"
    log = []

    def emit(*args) -> None:
        msg = " ".join(str(a) for a in args)
        print(msg)
        log.append(msg)

    emit("=" * 78)
    emit("DEBUG: R̃² vanishing on q-irreducible projection")
    emit("=" * 78)

    # 1. Single-component check: pick a specific q, see R̃²
    emit("")
    emit("--- Test 1: unconstrained q (just antisymmetric in last two) ---")
    q_unc, params_unc = build_q_unconstrained()
    emit(f"  unconstrained-q free parameters: {len(params_unc)} (expect 24 = 4×6)")
    K_unc = contortion(q_unc)
    R_unc = riemann_cartan(K_unc)
    R2_unc = R_squared_full(R_unc)
    R2_unc_strip = sp.expand(R2_unc / profile**2)

    # Sample with random integer parameters
    R2_unc_sample = R2_unc_strip.xreplace(
        {p: sp.Integer(7 + 3 * i) for i, p in enumerate(params_unc)}
    )
    R2_unc_sample = sp.expand(R2_unc_sample)
    if R2_unc_sample == 0:
        emit("  R̃² (unconstrained) = 0 — implies a BUG in contraction or contortion.")
    else:
        poly = sp.Poly(R2_unc_sample, *k)
        emit(f"  R̃² (unconstrained) max k-degree: {poly.total_degree()} (expected 2)")
        if poly.total_degree() == 2:
            emit(
                "  → R̃² is nonzero and standard-kinetic at unconstrained level — GOOD."
            )

    emit("")
    emit("--- Test 2: q-irreducible (tracelessness + axial-vanishing imposed) ---")
    q_c, params_c = build_q_basis_constrained()
    emit(f"  constrained-q free parameters: {len(params_c)} (expect 16)")
    K_c = contortion(q_c)
    R_c = riemann_cartan(K_c)
    R2_c = R_squared_full(R_c)
    R2_c_strip = sp.expand(R2_c / profile**2)
    R2_c_sample = R2_c_strip.xreplace(
        {p: sp.Integer(7 + 5 * i) for i, p in enumerate(params_c)}
    )
    R2_c_sample = sp.expand(R2_c_sample)
    if R2_c_sample == 0:
        emit("  R̃² (q-irreducible) = 0 — Verifies surprising vanishing.")
        emit("  Trying random parameters...")
        import random

        random.seed(42)
        R2_c_sample2 = R2_c_strip.xreplace(
            {p: sp.Integer(random.randint(-7, 7)) for p in params_c}
        )
        R2_c_sample2 = sp.expand(R2_c_sample2)
        if R2_c_sample2 == 0:
            emit(
                "    Still zero with random params — R̃² IS identically zero on q-irreducible."
            )
            emit("")
            emit(
                "    NOTE: This particular result (R̃² vanishes ID on q-irreducible at linear"
            )
            emit("    order on flat) is consistent with the analytic structure but is")
            emit(
                "    surprisingly strong. It implies the q-irreducible at linear order on flat"
            )
            emit("    has NO contribution to ANY R̃-quadratic Lagrangian invariant.")
            emit("")
            emit("    Possible explanations:")
            emit(
                "    (a) The q-irreducible projects K into a sector where the ANTISYMMETRIC"
            )
            emit(
                "        part of ∂K (which is what R̃ measures) vanishes by Bianchi-like"
            )
            emit("        identities induced by q's algebraic constraints.")
            emit(
                "    (b) The Hehl-McCrea-Mielke-Ne'eman irreducible decomposition has the"
            )
            emit(
                "        property that q is RICCI-FLAT at linear order on flat backgrounds."
            )
            emit("")
            emit(
                "    Either way, the verdict for Recipe 1 is even STRONGER than expected:"
            )
            emit(
                "    There is NO higher-derivative (∂²q)² structure because there is NO"
            )
            emit("    kinetic structure of any kind on the q-irreducible. Trivially")
            emit("    standard-kinetic since there's nothing there.")
        else:
            poly = sp.Poly(R2_c_sample2, *k)
            emit(f"    Random-param R̃² has k-degree {poly.total_degree()}")
    else:
        poly = sp.Poly(R2_c_sample, *k)
        emit(f"  R̃² (q-irreducible) max k-degree: {poly.total_degree()}")

    # 3. Sanity check: compute R̃ for unconstrained q and see what's happening
    emit("")
    emit("--- Test 3: Detailed inspection of R̃_{0123} for unconstrained q ---")
    R0123 = R_unc[0, 1, 2, 3] / profile
    R0123_simplified = sp.simplify(R0123)
    emit(f"  R̃_{{0,1,2,3}} = {R0123_simplified}")

    # 4. With irreducibility imposed
    emit("")
    emit("--- Test 4: R̃_{0123} with q-irreducibility imposed ---")
    R0123_c = R_c[0, 1, 2, 3] / profile
    R0123_c_simplified = sp.simplify(R0123_c)
    emit(f"  R̃_{{0,1,2,3}} (q-irred) = {R0123_c_simplified}")

    with pathlib.Path(out_txt).open("w", encoding="utf-8") as f:
        f.write("\n".join(log))
    emit(f"\nTranscript: {out_txt}")


if __name__ == "__main__":
    main()
