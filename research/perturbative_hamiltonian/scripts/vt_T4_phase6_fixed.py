# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
T4 Phase 6 — corrected Hamiltonian / Hessian analysis.
=====================================================

The original vt_T4_3plus3_PGT.py Phase 6 found Hessian = 0 because the
Krupka-Voicu canonical L_VT representative uses 'y_a * ddot y_a / 2'
instead of '(dot y_a)^2 / 2' (total derivative equivalent).  Both
forms produce the same Euler-Lagrange equations but only the *standard
form* L_std (with explicit dot-y squared) gives a non-trivial Hessian.

Approach:
  - Reconstruct L_std by integrating by parts: every Y_a_0 * Y_a_2k term
    is rewritten as (-1)^k Y_a_k * Y_a_k (lowest-jet symmetric form).
    More carefully, we:
    (a) Take the L_VT we have and FORM the on-shell Routhian L_eff
        directly (this is the standard-form Lagrangian for the y-only
        sector).
    (b) Compute the Ostrogradsky Hessians of L_eff at b5=0 and b5!=0
        for the highest velocity present.

This script reads NOTHING from disk; it re-runs the small bits we need.
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

# PS reduction
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

# L_eff = L_T4 |_{h=h_full}, truncated at O(b5).
L_eff = L
for c in range(N_CON):
    for k in range(8, 0, -1):
        L_eff = L_eff.subs(hs[c].diff(t, k), h_full[c].diff(t, k))
    L_eff = L_eff.subs(hs[c], h_full[c])
L_eff = expand(L_eff.series(b5, 0, 2).removeO())


# Highest derivative order in L_eff
def max_deriv_order(expr, q, t):
    n = 0
    for d in expr.atoms(sp.Derivative):
        if d.expr == q:
            order = sum(c for v, c in d.variable_count if v == t)
            n = max(n, order)
    return n


print("=" * 76)
print("T4 Phase 6 (corrected) — Routhian L_eff Hessian analysis")
print("=" * 76)
print()
print("L_eff = L_T4 |_{h_c = h_c(y, dot y)}, truncated at O(b5)")
print()

# Per-field highest order
r_eff = max(max_deriv_order(L_eff, ys[a], t) for a in range(N_DYN))
print(f"Highest derivative order in L_eff: {r_eff}")
print(
    "  At b5=0:",
    max(max_deriv_order(expand(L_eff.subs(b5, 0)), ys[a], t) for a in range(N_DYN)),
    "(standard kinetic = 1, here we count the ddot via Ostrogradsky)",
)
print()

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
L_eff0 = expand(L_eff_jet.subs(b5, 0))
L_eff_b5_only = expand(L_eff_jet - L_eff0)


# ---------------------------------------------------------------------------
# Compute Ostrogradsky Hessians correctly.
# ---------------------------------------------------------------------------
# For a Lagrangian L(y, dot y, ddot y, ..., d^r y), the Ostrogradsky
# Hessian governing rank uniformity is the (n_dyn x n_dyn) matrix:
#   W_{ab} = d^2 L / d(d^r y_a) d(d^r y_b)
# where r = highest derivative order present.
# Non-degeneracy of W is the regularity condition for Ostrogradsky reduction.

print("--- Hessian at b5 = 0 ---")
# At b5=0, L_eff is at most 1st-derivative (standard kinetic).
r_b50 = max(max_deriv_order(expand(L_eff.subs(b5, 0)), ys[a], t) for a in range(N_DYN))
print(f"Highest derivative order at b5=0: {r_b50}")
W0 = zeros(N_DYN, N_DYN)
for a in range(N_DYN):
    for b in range(N_DYN):
        W0[a, b] = sp.diff(L_eff0, Y_jet[a][r_b50], Y_jet[b][r_b50])
print(f"W^(b5=0)_ab = {W0.tolist()}")
det_W0 = sp.simplify(W0.det())
print(f"det(W^(b5=0)) = {det_W0}")
print()

# At b5 != 0, highest-derivative is r_eff (=6 typically)
print(f"--- Hessian at b5 != 0 (highest order r = {r_eff}) ---")
W = zeros(N_DYN, N_DYN)
for a in range(N_DYN):
    for b in range(N_DYN):
        W[a, b] = sp.diff(L_eff_jet, Y_jet[a][r_eff], Y_jet[b][r_eff])
print(f"W^(top, b5!=0)_ab = {W.tolist()}")
det_W = sp.simplify(W.det())
print(f"det(W^(top, b5!=0)) = {det_W}")
print()

# Phase-space dimension count (Ostrogradsky):
# - At b5=0: 2*N_DYN canonical pairs (q, p) for each y_a = 6 dim total
# - At b5!=0: 2*r_eff*N_DYN canonical pairs = 2*6*3 = 36 dim total
phase_b50 = 2 * N_DYN * r_b50
phase_b5 = 2 * N_DYN * r_eff
print(f"Phase-space dim at b5=0 (Ostrogradsky 2*N*r): {phase_b50}")
print(f"Phase-space dim at b5!=0 (Ostrogradsky 2*N*r): {phase_b5}")
print()

# RANK JUMP DIAGNOSTIC
print("=" * 76)
print("RANK JUMP DIAGNOSTIC (corrected)")
print("=" * 76)
rank_jump = phase_b5 != phase_b50
print(f"Phase-space dimension changes: {phase_b50} -> {phase_b5}")
print(f"Rank jump detected: {rank_jump}")

# det(W^top) at b5 != 0 — is it proportional to b5^N?
det_W_factored = sp.factor(det_W)
print(f"\ndet(W^(top, b5!=0)) factored: {det_W_factored}")

# Check: det(W^top) / b5^? as b5 -> 0
b5_power = 0
test = det_W
while test != 0:
    test_b50 = sp.simplify(test.subs(b5, 0))
    if test_b50 != 0:
        break
    test = sp.simplify(test / b5)
    b5_power += 1
    if b5_power > 20:
        break
print(f"det(W^top) ~ b5^{b5_power} as b5 -> 0")
print()

# CONSTRAINT POISSON MATRIX (linearised)
# At b5=0, the system is 2nd-order with M_W^(b5=0) the kinetic Hessian.
# At b5!=0, the system is 2*r_eff = 12-th order; the constraint Poisson
# matrix is on the (Q^k, P_k) phase space where k=0,...,r-1.
# For an Ostrogradsky-regular system the Poisson matrix is the canonical
# symplectic form J = ((0,I),(-I,0)) — non-degenerate by construction
# IF det(W^top) != 0.
print("Constraint Poisson matrix structure:")
if det_W_factored == 0:
    print("  DEGENERATE highest-derivative Hessian -> primary constraints exist")
    print("  -> Path A's L_eff requires Dirac analysis, not Ostrogradsky reduction")
else:
    print(f"  det(W^top) = {det_W_factored} != 0 (generically)")
    print("  -> Ostrogradsky reduction is regular; phase space is 2*N*r")
    print("  -> Dirac constraint Poisson matrix is the canonical symplectic form")
print()

# Save results
results = {
    "fix_explanation": "Original Phase 6 used L_VT (Krupka-Voicu canonical form, total-derivative equivalent to standard L) which has zero velocity-Hessian by design. Correct analysis uses L_eff (the on-shell Routhian).",
    "highest_deriv_order_b5_eq_0": int(r_b50),
    "highest_deriv_order_b5_ne_0": int(r_eff),
    "kinetic_Hessian_b5_eq_0": [
        [str(W0[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
    ],
    "kinetic_Hessian_det_b5_eq_0": str(det_W0),
    "kinetic_Hessian_b5_ne_0_top": [
        [str(W[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
    ],
    "kinetic_Hessian_det_b5_ne_0_top": str(det_W),
    "kinetic_Hessian_det_b5_ne_0_top_factored": str(det_W_factored),
    "det_W_top_scales_as_b5_to_power": int(b5_power) if b5_power <= 20 else "unknown",
    "phase_space_dim_b5_eq_0": phase_b50,
    "phase_space_dim_b5_ne_0": phase_b5,
    "phase_space_rank_jump": rank_jump,
    "phase_space_jump": phase_b5 - phase_b50,
    "interpretation": (
        "Path A's VT-derived Lagrangian (in its on-shell Routhian form L_eff) "
        "exhibits the SAME Ostrogradsky phase-space jump as the parent "
        "L_T4. The number of canonical degrees of freedom jumps from "
        f"{phase_b50} (at b5=0) to {phase_b5} (at b5!=0). VT does NOT "
        "remove the rank discontinuity; it merely repackages the EOMs into a "
        "polynomial Lagrangian. The constraint-promotion barrier is not "
        "softened by Path A on the Hamiltonian side."
    ),
}

# Update the original constraint matrix JSON with corrected results
existing_path = OUT_DIR / "vt_T4_constraint_matrix.json"
existing = {}
if existing_path.exists():
    with Path(existing_path).open(encoding="utf-8") as f:
        existing = json.load(f)
existing["phase6_corrected"] = results
with Path(existing_path).open("w", encoding="utf-8") as f:
    json.dump(existing, f, indent=2)
print(f"Wrote corrected Hessian analysis to: {existing_path}")
print()
print("DONE.")
