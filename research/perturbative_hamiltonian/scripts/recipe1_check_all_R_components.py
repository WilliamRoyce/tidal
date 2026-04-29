# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""Check all R̃ components on q-irreducible to see if they all vanish."""

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


def build_q_basis():
    syms = sp.symbols("c0:4_0:4_0:4")
    c = {}
    idx = 0
    for a, b, cc in product(range(4), repeat=3):
        c[a, b, cc] = syms[idx]
        idx += 1
    constraints = []
    for a, b, cc in product(range(4), repeat=3):
        constraints.append(c[a, b, cc] + c[a, cc, b])
    constraints.extend(
        sum(ETA_INV[a, a] * c[a, a, cc] for a in range(4)) for cc in range(4)
    )
    constraints.extend(
        sum(
            eps_down(a, b, cc, d) * c[a, b, cc]
            for a, b, cc in product(range(4), repeat=3)
        )
        for d in range(4)
    )
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


def contortion(q):
    K = {}
    for r, m, n in product(range(4), repeat=3):
        K[r, m, n] = sp.Rational(1, 2) * (q[r, m, n] + q[m, r, n] + q[n, r, m])
    return K


def riemann_cartan(K):
    R = {}
    for r, s, m, n in product(range(4), repeat=4):
        d_m = sp.diff(K[r, n, s] * profile, COORDS[m])
        d_n = sp.diff(K[r, m, s] * profile, COORDS[n])
        R[r, s, m, n] = sp.expand(d_m - d_n)
    return R


def main() -> None:
    out_dir = pathlib.Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True, parents=True)
    out_txt = out_dir / "recipe1_R_components_run.txt"
    log = []

    def emit(*args) -> None:
        msg = " ".join(str(a) for a in args)
        print(msg)
        log.append(msg)

    q, params = build_q_basis()
    K = contortion(q)
    R = riemann_cartan(K)

    # How many R-components are nonzero?
    n_nonzero = 0
    nonzero_examples = []
    for rsmn in product(range(4), repeat=4):
        e = sp.expand(R[rsmn] / profile)
        if e != 0:
            n_nonzero += 1
            if len(nonzero_examples) < 5:
                nonzero_examples.append((rsmn, e))

    emit(f"Number of nonzero R̃-components on q-irreducible: {n_nonzero} out of 256")
    emit("First 5 nonzero examples:")
    for rsmn, e in nonzero_examples:
        emit(f"  R̃[{rsmn}] = {e}")

    # Now check K-components: how many are nonzero?
    n_K_nonzero = 0
    for rmn in product(range(4), repeat=3):
        if K[rmn] != 0:
            n_K_nonzero += 1
    emit(f"\nNumber of nonzero K-components: {n_K_nonzero} out of 64")

    # Most important: is K=0 identically?
    K_total = sum(sp.Abs(K[rmn]) for rmn in product(range(4), repeat=3))
    K_total_sub = K_total.xreplace(
        {p: sp.Integer(7 + 5 * i) for i, p in enumerate(params)}
    )
    emit(f"\nSum of |K| with sample params: {K_total_sub}")

    with pathlib.Path(out_txt).open("w", encoding="utf-8") as f:
        f.write("\n".join(log))
    emit(f"\nTranscript: {out_txt}")


if __name__ == "__main__":
    main()
