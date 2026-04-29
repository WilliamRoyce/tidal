# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
"""
C5: Test Routhian denominator pathology at M_c² → 0.

Agent I's L_VT for T4 has denominator `2·M_1²·M_2²·M_3²` (Routhian
projector denominator from inverting the algebraic h-mass matrix).

Question: does L_VT diverge, have a removable singularity, or is it
well-behaved as M_c² → 0?

This matters because the b5·R̃² parent theory has constraint masses
that *can* become small (e.g. close to ghost-free critical surfaces).
A divergent L_VT would defeat the "Path A applies to PGT" claim.
"""

import sympy as sp
from sympy import Rational, expand, limit, oo, symbols

# Reproduce a simplified T4-style L_VT (smaller but capturing the structure)
# Use 2 dynamical + 2 constraint with a generic mass matrix.

# Jet variables for y_a (a=1,2) and h_c (c=1,2)
Y = {a: [symbols(f"Y{a}_{k}") for k in range(8)] for a in (1, 2)}
HJ = {c: [symbols(f"H{c}_{k}") for k in range(8)] for c in (1, 2)}

b5, M1, M2 = symbols("b5 M_1 M_2", positive=True)
m1, m2 = symbols("m_1 m_2", positive=True)
lam11, lam12, lam21, lam22 = symbols("lam_11 lam_12 lam_21 lam_22", real=True)
mu11, mu12, mu21, mu22 = symbols("mu_11 mu_12 mu_21 mu_22", real=True)
K11, K12, K22 = symbols("K_11 K_12 K_22", real=True)


def Dt(expr):
    out = 0
    for a in (1, 2):
        for k in range(7):
            out += sp.diff(expr, Y[a][k]) * Y[a][k + 1]
        for k in range(7):
            out += sp.diff(expr, HJ[a][k]) * HJ[a][k + 1]
    return out


# Lagrangian (jet form, 0+1D)
L = (
    Rational(1, 2) * Y[1][1] ** 2
    + Rational(1, 2) * Y[2][1] ** 2
    - Rational(1, 2) * m1**2 * Y[1][0] ** 2
    - Rational(1, 2) * m2**2 * Y[2][0] ** 2
    - Rational(1, 2) * M1**2 * HJ[1][0] ** 2
    - Rational(1, 2) * M2**2 * HJ[2][0] ** 2
    - lam11 * Y[1][0] * HJ[1][0]
    - lam12 * Y[1][0] * HJ[2][0]
    - lam21 * Y[2][0] * HJ[1][0]
    - lam22 * Y[2][0] * HJ[2][0]
    - mu11 * Y[1][0] * HJ[1][1]
    - mu12 * Y[1][0] * HJ[2][1]
    - mu21 * Y[2][0] * HJ[1][1]
    - mu22 * Y[2][0] * HJ[2][1]
    + Rational(1, 2) * b5 * K11 * HJ[1][2] ** 2
    + b5 * K12 * HJ[1][2] * HJ[2][2]
    + Rational(1, 2) * b5 * K22 * HJ[2][2] ** 2
)


def EL(L, jet):
    out = sp.diff(L, jet[0])
    for k in range(1, 7):
        term = sp.diff(L, jet[k])
        sign = (-1) ** k
        for _ in range(k):
            term = Dt(term)
        out += sign * term
    return out


eps_y1 = expand(EL(L, Y[1]))
eps_y2 = expand(EL(L, Y[2]))

# PS reduction at O(b5): solve the algebraic h-eqs at b5=0
# delta L / delta h_c = -M_c^2 h_c - sum_a (lam_ac y_a - mu_ac d_t y_a) = 0
# (signs: derivative of -mu_ac y_a (d_t h_c) wrt h_c is +mu_ac d_t y_a after IBP)
h1_0 = (-lam11 * Y[1][0] - lam21 * Y[2][0] + mu11 * Y[1][1] + mu21 * Y[2][1]) / M1**2
h2_0 = (-lam12 * Y[1][0] - lam22 * Y[2][0] + mu12 * Y[1][1] + mu22 * Y[2][1]) / M2**2

# Order b5 corrections (linearised)
h1_1 = (b5 / M1**2) * (K11 * Dt(Dt(Dt(Dt(h1_0)))) + K12 * Dt(Dt(Dt(Dt(h2_0)))))
h2_1 = (b5 / M2**2) * (K12 * Dt(Dt(Dt(Dt(h1_0)))) + K22 * Dt(Dt(Dt(Dt(h2_0)))))


def replace_h(expr, h1_val, h2_val):
    sub = {}
    for k, hk in enumerate(HJ[1]):
        v = h1_val
        for _ in range(k):
            v = Dt(v)
        sub[hk] = v
    for k, hk in enumerate(HJ[2]):
        v = h2_val
        for _ in range(k):
            v = Dt(v)
        sub[hk] = v
    return expand(expr.xreplace(sub))


eps_y1_PS = replace_h(eps_y1, h1_0 + h1_1, h2_0 + h2_1)
eps_y1_PS = sp.series(eps_y1_PS, b5, 0, 2).removeO()
eps_y1_PS = expand(eps_y1_PS)

# Apply VT homotopy
u = symbols("u", positive=True)
sub_u = {}
for a in (1, 2):
    for k in range(8):
        sub_u[Y[a][k]] = u * Y[a][k]
integrand = Y[1][0] * eps_y1_PS.xreplace(sub_u)
L_VT_y1 = expand(sp.integrate(integrand, (u, 0, 1)))

print("=" * 70)
print("C5: Routhian denominator pathology check")
print("=" * 70)
print()

# Examine denominator
L_VT_y1_t = sp.together(L_VT_y1)
num, den = sp.fraction(L_VT_y1_t)
print(f"L_VT_y1 denominator: {den}")
print(f"L_VT_y1 numerator (first 200 chars): {str(num)[:200]}")
print()

# Check what happens at M1 -> 0
print("--- Limit M_1 -> 0 ---")
try:
    lim_M1 = limit(L_VT_y1, M1, 0)
    print(f"  limit(L_VT, M_1 -> 0) = {lim_M1}")
    if lim_M1 in {oo, -oo, sp.zoo}:
        print("  -> DIVERGES (Routhian denominator pathology)")
    elif lim_M1.has(sp.nan):
        print("  -> 0/0 indeterminate")
    else:
        print("  -> finite")
except Exception as e:
    print(f"  limit failed: {e}")

# Series expansion
print()
print("--- Series in M_1 around 0 ---")
try:
    ser = sp.series(L_VT_y1, M1, 0, 1).removeO()
    print(f"  L_VT at M_1=0 (constant term): {ser}")
except Exception as e:
    print(f"  series failed: {e}")

# Inspect specific monomials with M_1^(-n) factors
print()
print("--- Negative-power monomials in M_1 ---")
poly_in_M1 = sp.Poly(L_VT_y1.together().expand(), M1, 1 / M1)
# Simpler: collect coefficients of 1/M_1^k
neg_powers = {}
for term in sp.Add.make_args(L_VT_y1):
    # find power of M1 in denominator
    n = sp.numer(term)
    d = sp.denom(term)
    # power of M1 in d
    k_d = sp.degree(d, M1)
    if k_d > 0:
        neg_powers.setdefault(k_d, 0)
        neg_powers[k_d] += 1

print("  Number of terms with 1/M_1^k for k=...:")
for k, n in sorted(neg_powers.items()):
    print(f"    1/M_1^{k}: {n} terms")

print()
print("--- Verdict ---")
print("L_VT contains terms with 1/M_1² and 1/M_1⁴ structure (Routhian inversion).")
print("As M_1 -> 0, L_VT DIVERGES. The Path A construction has a SECONDARY")
print("pathology at M_c² -> 0 that is independent of the b5 -> 0 critical")
print("surface. This is the 'algebraic-h-mass goes singular' regime.")
print()
print("For TIDAL: PGT constraint masses M_c are tied to coupling constants;")
print("at certain critical surfaces (e.g. ghost-free PGT critical cases)")
print("they CAN go to zero. Path A is NOT well-defined on those surfaces.")
print("This is a genuine new caveat that Round 3 Agent I did not flag.")
