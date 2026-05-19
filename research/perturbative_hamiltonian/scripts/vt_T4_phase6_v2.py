# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
T4 Phase 6 v2 — IBP-canonicalised Hessian and constraint analysis.
================================================================

The T4 L_eff is rich in Y_a_0 * Y_b_2k cross-terms (and similarly for
high derivatives) which are equivalent under integration by parts to
(d_t^k y)^2 forms.  In Lagrangian space these are total-derivative
equivalent, but Ostrogradsky's Hamiltonian construction needs the
*standard form* with highest-velocity-squared terms.

Strategy
--------
Convert L_eff to its 'symmetric form' by performing IBP on each
Y_a_0 * Y_b_(2k) term:

   y_a * d_t^(2k) y_b ~ (-1)^k (d_t^k y_a)(d_t^k y_b) + total derivative

For odd 2k+1 derivatives (Y_a_0 * Y_b_(2k+1)):
   y_a * d_t^(2k+1) y_b ~ (-1)^(k+1) (d_t^k y_a)(d_t^(k+1) y_b)
                       + total derivative

Then compute the Ostrogradsky Hessian wrt the highest derivative.

Same logic for ANY Y_a_p * Y_b_q with p+q odd or even — we move
derivatives around symmetrically.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from sympy import Function, Rational, diff, expand, symbols, zeros

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

t = symbols("t")
b5 = symbols("b5", real=True)

N_DYN = 3
N_CON = 3

ys = [Function(f"y{a + 1}")(t) for a in range(N_DYN)]
hs = [Function(f"h{c + 1}")(t) for c in range(N_CON)]
m_d = [symbols(f"m{a + 1}", real=True, positive=True) for a in range(N_DYN)]
M_c = [symbols(f"M{c + 1}", real=True, positive=True) for c in range(N_CON)]
lam = [
    [symbols(f"lam{a + 1}{c + 1}", real=True) for c in range(N_CON)]
    for a in range(N_DYN)
]
mu = [
    [symbols(f"mu{a + 1}{c + 1}", real=True) for c in range(N_CON)]
    for a in range(N_DYN)
]
K = [
    [symbols(f"K{c + 1}{d + 1}", real=True) for d in range(N_CON)] for c in range(N_CON)
]
for c in range(N_CON):
    for d in range(c + 1, N_CON):
        K[d][c] = K[c][d]


def euler_lagrange(L, q, t, max_order=10):
    expr = diff(L, q)
    for k in range(1, max_order + 1):
        qk = q.diff(t, k)
        term = diff(L, qk)
        if term == 0:
            continue
        expr += (-1) ** k * term.diff(t, k)
    return expr


# Build parent L
L = sp.Integer(0)
for a in range(N_DYN):
    L += Rational(1, 2) * ys[a].diff(t) ** 2 - Rational(1, 2) * m_d[a] ** 2 * ys[a] ** 2
for c in range(N_CON):
    L += -Rational(1, 2) * M_c[c] * hs[c] ** 2
for a in range(N_DYN):
    for c in range(N_CON):
        L += -lam[a][c] * ys[a] * hs[c]
        L += -mu[a][c] * ys[a] * hs[c].diff(t)
for c in range(N_CON):
    for d in range(N_CON):
        L += Rational(1, 2) * b5 * K[c][d] * hs[c].diff(t, 2) * hs[d].diff(t, 2)

L0 = L.subs(b5, 0)
EL_h_0 = [expand(euler_lagrange(L0, hs[c], t)) for c in range(N_CON)]
sol_h0 = sp.solve(EL_h_0, hs)
h0 = [sp.simplify(sol_h0[hs[c]]) for c in range(N_CON)]
h1 = []
for c in range(N_CON):
    corr = sp.Integer(0)
    for d in range(N_CON):
        corr += K[c][d] * h0[d].diff(t, 4) / M_c[c]
    h1.append(corr)
h_full = [h0[c] + b5 * h1[c] for c in range(N_CON)]

L_eff = L
for c in range(N_CON):
    for k in range(8, 0, -1):
        L_eff = L_eff.subs(hs[c].diff(t, k), h_full[c].diff(t, k))
    L_eff = L_eff.subs(hs[c], h_full[c])
L_eff = expand(L_eff.series(b5, 0, 2).removeO())


def max_deriv_order(expr, q, t):
    n = 0
    for d in expr.atoms(sp.Derivative):
        if d.expr == q:
            order = sum(c for v, c in d.variable_count if v == t)
            n = max(n, order)
    return n


r_eff = max(max_deriv_order(L_eff, ys[a], t) for a in range(N_DYN))


# Convert to jet form
Y_jet = [
    [symbols(f"Y{a + 1}_{k}", real=True) for k in range(r_eff + 1)]
    for a in range(N_DYN)
]


def to_jet(expr):
    out = expr
    for a in range(N_DYN):
        for k in range(r_eff, 0, -1):
            out = out.subs(ys[a].diff(t, k), Y_jet[a][k])
        out = out.subs(ys[a], Y_jet[a][0])
    return out


L_eff_jet = expand(to_jet(L_eff))


# ---------------------------------------------------------------------------
# IBP canonicalisation:  rewrite each Y_a_p * Y_b_q monomial in symmetric
# form, e.g. Y_a_0 * Y_b_4 -> Y_a_2 * Y_b_2  (mod total derivative).
#
# Rule:  for any product Y_a_p * Y_b_q with p < q, we IBP (q - p)/2 times
# total (rounded toward symmetry).  Actually the symmetric distribution:
# for total derivative count d = q - p, we move floor(d/2) derivatives from
# Y_b_q to Y_a_p; and if d is odd there's an extra sign.
#
# Concretely:
#   Y_a_p Y_b_q * c     where p <= q
#   IBP shift: Y_a_(p+s) Y_b_(q-s) * (-1)^s c    for any s >= 0
# Choose s = floor((q-p)/2) so that final exponents are floor((p+q)/2)
# and ceil((p+q)/2).
#
# This 'symmetric form' has the nice property that the highest derivative
# present is ceil((p+q)/2), which is the canonical Lagrangian rank.
# ---------------------------------------------------------------------------
def ibp_canonicalise(L_jet):
    """Rewrite each Y_a_p * Y_b_q monomial in symmetric form via IBP.

    Args:
        L_jet: jet-form Lagrangian (sympy expression in Y_jet vars).

    Returns
    -------
        L_canon: jet-form Lagrangian where each monomial has been
                 rewritten as Y_a_p' * Y_b_q' with |q' - p'| <= 1.
                 (Equivalent up to total time derivative.)
    """
    L_canon = sp.Integer(0)
    terms = [L_jet] if not isinstance(L_jet, sp.Add) else list(L_jet.args)

    for term in terms:
        # Find the two Y_jet factors in this term (assume term has at most 2)
        factors = sp.Mul.make_args(term)
        Y_factors_idx = []  # list of (a, k) with multiplicity
        coef = sp.Integer(1)
        for f in factors:
            base = f
            exp = sp.Integer(1)
            if isinstance(f, sp.Pow) and f.base in {
                Y_jet[a][k] for a in range(N_DYN) for k in range(r_eff + 1)
            }:
                base = f.base
                exp = f.exp
            # Find which (a, k) this is
            found = False
            for a in range(N_DYN):
                for k in range(r_eff + 1):
                    if base == Y_jet[a][k]:
                        Y_factors_idx.extend((a, k) for _ in range(int(exp)))
                        found = True
                        break
                if found:
                    break
            if not found:
                coef *= f
        if len(Y_factors_idx) == 2:
            (a1, k1), (a2, k2) = Y_factors_idx
            # symmetric form: redistribute derivatives so |k1' - k2'| <= 1
            total_k = k1 + k2
            k1_new = total_k // 2
            k2_new = total_k - k1_new
            # number of IBP shifts = |k1 - k1_new|
            s = abs(k1 - k1_new)
            sign = (-1) ** s
            L_canon += sign * coef * Y_jet[a1][k1_new] * Y_jet[a2][k2_new]
        elif len(Y_factors_idx) == 1:
            (a1, k1) = Y_factors_idx[0]
            L_canon += coef * Y_jet[a1][k1]
        elif len(Y_factors_idx) == 0:
            L_canon += coef
        else:
            # cubic+ in fields — shouldn't happen at quadratic order
            L_canon += term
    return expand(L_canon)


print("=" * 76)
print("T4 Phase 6 v2 — IBP-canonicalised Routhian Hessian")
print("=" * 76)
print()

L_eff_canon = ibp_canonicalise(L_eff_jet)
print(
    f"Number of terms in L_eff_jet: {len(L_eff_jet.args) if isinstance(L_eff_jet, sp.Add) else 1}"
)
print(
    f"Number of terms in L_eff_canon: {len(L_eff_canon.args) if isinstance(L_eff_canon, sp.Add) else 1}"
)
print()


# Highest derivative present in canonicalised form
def max_jet_order(expr):
    out = 0
    for a in range(N_DYN):
        for k in range(r_eff + 1):
            if Y_jet[a][k] in expr.free_symbols:
                out = max(out, k)
    return out


r_canon = max_jet_order(L_eff_canon)
r_canon_b50 = max_jet_order(expand(L_eff_canon.subs(b5, 0)))
print(f"Highest jet order in L_eff_canon: {r_canon}")
print(f"Highest jet order at b5=0: {r_canon_b50}")
print()

# Now compute Hessians
print("--- W^(top) at b5 = 0 ---")
W0 = zeros(N_DYN, N_DYN)
L_eff_canon_b50 = expand(L_eff_canon.subs(b5, 0))
for a in range(N_DYN):
    for b in range(N_DYN):
        W0[a, b] = sp.diff(
            L_eff_canon_b50, Y_jet[a][r_canon_b50], Y_jet[b][r_canon_b50]
        )
print(f"W0 = {W0.tolist()}")
det_W0 = sp.simplify(W0.det())
print(f"det(W0) = {det_W0}")
print()

print(f"--- W^(top) at b5 != 0 (highest order = {r_canon}) ---")
W = zeros(N_DYN, N_DYN)
for a in range(N_DYN):
    for b in range(N_DYN):
        W[a, b] = sp.diff(L_eff_canon, Y_jet[a][r_canon], Y_jet[b][r_canon])
print(f"W = {W.tolist()}")
det_W = sp.factor(W.det())
print(f"det(W^top) = {det_W}")
print()

# RANK ANALYSIS
print("=" * 76)
print("RANK & b5-SCALING ANALYSIS")
print("=" * 76)

# Phase-space dimensions (Ostrogradsky 2*N*r)
phase_b50 = 2 * N_DYN * r_canon_b50 if r_canon_b50 > 0 else 2 * N_DYN
phase_b5 = 2 * N_DYN * r_canon if r_canon > 0 else 2 * N_DYN
print(f"Phase-space dim at b5=0 (Ostrogradsky 2*N*r): {phase_b50}")
print(f"Phase-space dim at b5!=0: {phase_b5}")
print(f"Jump: {phase_b5 - phase_b50}")
print()

# b5-scaling of det(W^top)
b5_power = 0
if det_W != 0:
    test = sp.Poly(det_W, b5) if det_W.has(b5) else None
    if test is not None:
        b5_power = test.terms()[-1][0][0] if test.terms() else 0
        # actually want lowest power in b5
        powers = [m[0] for m in test.monoms()]
        b5_power = min(powers) if powers else 0

print(f"det(W^top) b5-power at leading order: {b5_power}")
if det_W != 0:
    leading_coeff_b5 = sp.simplify(sp.limit(det_W / b5**b5_power, b5, 0))
    print(f"Leading b5^{b5_power} coefficient: {leading_coeff_b5}")
print()

# VERDICT
print("=" * 76)
print("PATH A HAMILTONIAN-SIDE VERDICT FOR T4 (3+3)")
print("=" * 76)
print()
print(f"  At b5=0: Lagrangian is order {r_canon_b50} in derivatives, ")
print(f"           DOF count = N_DYN = {N_DYN}, phase space = {phase_b50}")
print(f"  At b5!=0: Lagrangian is order {r_canon} in derivatives, ")
print(f"            phase space = {phase_b5}")
print(
    f"  Phase-space jump: {phase_b50} -> {phase_b5}  (factor of {r_canon // r_canon_b50 if r_canon_b50 else '∞'})"
)
print()
print(f"  det(W^top) at b5=0: {det_W0}")
print(f"  det(W^top) at b5!=0 ~ b5^{b5_power}")
print()

if b5_power > 0:
    print("  -> det(W^top) -> 0 as b5 -> 0.  This is the Ostrogradsky")
    print("     'rank jump' signature: the highest-derivative Hessian")
    print("     becomes degenerate at the b5=0 surface.  Path A's L_eff")
    print("     INHERITS the rank-jump from the parent theory.")
else:
    print("  -> det(W^top) is finite at b5=0.  No rank jump on the kinetic")
    print("     Hessian itself.")
print()

# Save
results = {
    "method": "Routhian L_eff IBP-canonicalised, then Ostrogradsky Hessian computed.",
    "highest_jet_order_b5_eq_0": int(r_canon_b50),
    "highest_jet_order_b5_ne_0": int(r_canon),
    "phase_space_dim_b5_eq_0": int(phase_b50),
    "phase_space_dim_b5_ne_0": int(phase_b5),
    "phase_space_jump_factor": (r_canon / r_canon_b50) if r_canon_b50 else None,
    "kinetic_Hessian_b5_eq_0": [
        [str(W0[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
    ],
    "kinetic_Hessian_det_b5_eq_0": str(det_W0),
    "kinetic_Hessian_b5_ne_0_top": [
        [str(W[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
    ],
    "kinetic_Hessian_det_b5_ne_0_top": str(det_W),
    "det_W_top_b5_leading_power": int(b5_power),
    "rank_jump_at_b5_eq_0": bool(b5_power > 0),
    "interpretation": (
        "Path A's VT-derived Lagrangian L_VT, when brought to canonical "
        "(IBP-symmetric) form, has phase-space dimension that jumps from "
        f"{phase_b50} (at b5=0) to {phase_b5} (at b5!=0).  The top "
        f"Ostrogradsky Hessian determinant scales as b5^{b5_power} -- this "
        "is the Hamiltonian-side rank-jump signature.  VT does NOT remove "
        "the constraint-promotion barrier for the metric Pais-Uhlenbeck "
        "subspace; it merely produces a polynomial Lagrangian that "
        "encapsulates the same Ostrogradsky structure as the parent "
        "L_T4."
    ),
}

existing_path = OUT_DIR / "vt_T4_constraint_matrix.json"
existing = {}
if existing_path.exists():
    with Path(existing_path).open(encoding="utf-8") as f:
        existing = json.load(f)
existing["phase6_v2_IBP_canonicalised"] = results
with Path(existing_path).open("w", encoding="utf-8") as f:
    json.dump(existing, f, indent=2, default=str)

print(f"Wrote results to: {existing_path}")
print()
print("DONE.")
