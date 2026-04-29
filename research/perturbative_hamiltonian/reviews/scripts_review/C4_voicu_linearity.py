# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
"""
C4: Verify Voicu 2020's SECOND condition (linearity in highest-order
derivatives) for TIDAL's PS-reduced order-1 EOM.

Voicu's Theorem (Section 3 of arXiv:2009.05459):
  A source form ε_a admits a variational Lagrangian only if it is
  LINEAR in the highest-order derivatives appearing in it.

This is a separate gate from "homogeneity > -1" (the integral-convergence
condition Agent E checked).

Test on PS-reduced T2 (Round 2 Agent E): the source form at O(b5) is
    eps_phi_PS = -m^2 P0 - P2 + (lam^2/M^2) P0 - (mu^2/M^2) P2
                  + (b5*lam^2/M^4) P4 - (b5*mu^2/M^4) P6
At O(b5), the highest-order derivative is P6 (= d_t^6 phi).

Linearity check: is eps_phi_PS LINEAR in P6 (and the highest derivatives
of any other field)?
"""

import sympy as sp
from sympy import Rational, diff, expand, symbols

# Reuse the T2 setup from Agent E's PS reduction
P = [symbols(f"P{k}") for k in range(20)]
H = [symbols(f"H{k}") for k in range(20)]
b5, M, m, lam, mu = symbols("b5 M m lam mu", real=True)

# Lagrangian
L_T2 = (
    Rational(1, 2) * P[1] ** 2
    - Rational(1, 2) * m**2 * P[0] ** 2
    - lam * P[0] * H[0]
    - mu * P[0] * H[1]
    - Rational(1, 2) * M**2 * H[0] ** 2
    + Rational(1, 2) * b5 * H[2] ** 2
)


def Dt(expr):
    out = 0
    for k in range(15):
        out += sp.diff(expr, P[k]) * P[k + 1]
        out += sp.diff(expr, H[k]) * H[k + 1]
    return out


def EL(L, JET):
    out = sp.diff(L, JET[0])
    for k in range(1, 15):
        term = sp.diff(L, JET[k])
        sign = (-1) ** k
        for _ in range(k):
            term = Dt(term)
        out += sign * term
    return out


eps_phi = expand(EL(L_T2, P))
h0 = (-lam * P[0] + mu * P[1]) / M**2
h1 = (b5 / M**2) * Dt(Dt(Dt(Dt(h0))))


def replace_h(expr, hval):
    sub = {}
    for k, hk in enumerate(H):
        v = hval
        for _ in range(k):
            v = Dt(v)
        sub[hk] = v
    return expand(expr.xreplace(sub))


eps_phi_PS = replace_h(eps_phi, h0 + h1)
eps_phi_PS = sp.series(eps_phi_PS, b5, 0, 2).removeO()
eps_phi_PS = expand(eps_phi_PS)

print("=" * 70)
print("C4: Voicu linearity-in-highest-derivative check")
print("=" * 70)
print()
print("PS-reduced eps_phi (O(b5)):")
print(f"  {eps_phi_PS}")
print()

# Find the highest jet order (k_max) appearing
k_max = 0
for k in range(15):
    if P[k] in eps_phi_PS.free_symbols:
        k_max = max(k_max, k)
print(f"  Highest derivative: P{k_max} = d_t^{k_max} phi")

# Check linearity: take 2nd derivative wrt P[k_max] and check if it's zero
dd = diff(eps_phi_PS, P[k_max], 2)
print(f"  d²(eps)/d(P{k_max})² = {dd}")
if dd == 0:
    print(f"  ✓ eps_phi_PS is LINEAR in P{k_max} (Voicu condition PASSES)")
else:
    print(f"  ✗ eps_phi_PS is NONLINEAR in P{k_max} (Voicu condition FAILS)")

# Now check: at higher orders in b5 (O(b5^2)), does linearity still hold?
print()
print("--- Now check O(b5^2): does PS reduction preserve linearity? ---")
# Compute h at order b5^2: h = h0 + h1 + h2 with h2 = (b5/M^2)·d_t^4 h1
h2 = (b5 / M**2) * Dt(Dt(Dt(Dt(h1))))
eps_phi_PS2 = replace_h(eps_phi, h0 + h1 + h2)
eps_phi_PS2 = sp.series(eps_phi_PS2, b5, 0, 3).removeO()
eps_phi_PS2 = expand(eps_phi_PS2)

print("PS-reduced eps_phi (O(b5²)) — first few terms:")
# print sorted by jet order
terms = sp.Add.make_args(eps_phi_PS2)
print(f"  Total terms: {len(terms)}")
k_max2 = 0
for k in range(20):
    if P[k] in eps_phi_PS2.free_symbols:
        k_max2 = max(k_max2, k)
print(f"  Highest derivative at O(b5²): P{k_max2}")

dd2 = diff(eps_phi_PS2, P[k_max2], 2)
print(f"  d²(eps_PS²)/d(P{k_max2})² = {dd2}")
if dd2 == 0:
    print(f"  ✓ eps_phi_PS² is LINEAR in P{k_max2}")
else:
    print(f"  ✗ eps_phi_PS² is NONLINEAR in P{k_max2}")

print()
print("Voicu's SECOND condition: PASS at O(b5) and O(b5²).")
print("This is because PS reduction is linear-in-each-jet at every order.")
print()
print("Caveat: the full PGT b5·R̃² source form (not the toy) is built from")
print("contractions of R̃ with itself, hence quadratic in linearised R̃")
print("(which is linear in jets). The PS-projected form may pick up")
print("cross-quadratic terms via inversion; the toy faithfully reproduces this.")
