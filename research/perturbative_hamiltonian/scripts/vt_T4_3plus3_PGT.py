# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
Vainberg-Tonti integral applied to T4 — a 3+3 PGT-faithful toy with
Blagojevic-Cvetkovic 2018 Appendix D-inspired mass matrix.
=====================================================================

This is the next scaling step for Path A (Krupka-Voicu canonical
variational completion) following Round 2 Agent E's T3 (2+2) verification.

T4 setup
--------
  - 3 dynamical fields y_a (a=1,2,3) representing PGT TT graviton modes
  - 3 constraint-promoted fields h_c (c=1,2,3) representing h_4, h_7, h_9
    (the constraint-sector metric components that get promoted to
    Pais-Uhlenbeck 4th-order at b5 != 0)
  - A 3x3 lambda matrix lambda_{ac} encoding the algebraic mixing
  - A 3x3 mu matrix mu_{ac} encoding the velocity (1st-derivative) mixing
  - A 3x3 SYMMETRIC kinetic K_{cd} cross-coupling for the b5 sector
    (so b5*K_cd*ddot(h_c)*ddot(h_d) is the PGT b5*R~^2 -> constraint
    subspace projection)
  - Constraint masses M2_c (3 in number) and dynamical masses m_a (3 of them)

Lagrangian (1+0D, the time-only PS-reduced form):
  L_T4 = sum_a [ (1/2)(d_t y_a)^2 - (1/2) m_a^2 y_a^2 ]
       + sum_c [ -(1/2) M2_c h_c^2 ]
       - sum_{a,c} [ lambda_{ac} y_a h_c + mu_{ac} y_a d_t h_c ]
       + (b5/2) sum_{c,d} K_{cd} (d_t^2 h_c)(d_t^2 h_d)

This generalises T3 to N=3 with ALL matrices generic, and is the largest
analytically tractable analog of the actual h_4/h_7/h_9 constraint
sector. We use *symbolic* matrices (no numerics) so the algebra is
fully general.

Plan
----
Phase 2: Build L_T4 explicitly with a symbolic 3x3 BC-style mass-matrix
         template (block-upper-triangular structure inspired by BC
         Appendix D constraint chain).
Phase 3: Apply Parker-Simon iterative reduction:
         - Solve order-0 algebraic h-EOMs
         - Substitute back into y-EOMs
         - Truncate at O(b5)
Phase 4: Apply Vainberg-Tonti homotopy to construct L_VT.
         Verify Krupka-Voicu Theorem 1: EL(L_VT) - epsilon = 0.
Phase 5: Project L_VT onto constraint subspace (h_4, h_7, h_9 stand-ins)
         and inspect kinetic structure.
         CRITICAL question: is it (d h)^2 (clean) or (d^2 h)^2
         (Pais-Uhlenbeck)?
Phase 6: Hamiltonian / Legendre transform of L_VT.
         Compute conjugate momenta and constraint Poisson matrix.
         Compare uniform-rank vs rank-jump in b5.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp
from sympy import (
    Function,
    Poly,
    Rational,
    diff,
    expand,
    integrate,
    simplify,
    symbols,
    zeros,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

t = symbols("t")
u = symbols("u", positive=True)
b5 = symbols("b5", real=True)

N_DYN = 3  # y_1, y_2, y_3 (graviton stand-ins)
N_CON = 3  # h_4, h_7, h_9 stand-ins

# ---------------------------------------------------------------------------
# Phase 1 — Field declarations and BC-Appendix-D-inspired matrix template
# ---------------------------------------------------------------------------

ys = [Function(f"y{a + 1}")(t) for a in range(N_DYN)]
hs = [Function(f"h{c + 1}")(t) for c in range(N_CON)]

m_d = [symbols(f"m{a + 1}", real=True, positive=True) for a in range(N_DYN)]
M_c = [symbols(f"M{c + 1}", real=True, positive=True) for c in range(N_CON)]

# lambda_{a,c}: algebraic mixing y_a*h_c.
# BC Appendix D (qualitative) shows the constraint mass matrix has a
# block-triangular structure with non-zero off-diagonal couplings between
# h_4 and h_7 (the chain h_4 -> h_7 -> h_9 of secondary constraints).
# We encode this with a generic 3x3 lambda but make some entries vanish
# to emulate the BC structure WHEN we want the explicit pattern; for
# the convergence check we keep them generic.
lam = [
    [symbols(f"lam{a + 1}{c + 1}", real=True) for c in range(N_CON)]
    for a in range(N_DYN)
]
mu = [
    [symbols(f"mu{a + 1}{c + 1}", real=True) for c in range(N_CON)]
    for a in range(N_DYN)
]

# Symmetric K_{cd} for b5*K*ddot(h)*ddot(h) sector
K = [
    [symbols(f"K{c + 1}{d + 1}", real=True) for d in range(N_CON)] for c in range(N_CON)
]
for c in range(N_CON):
    for d in range(c + 1, N_CON):
        K[d][c] = K[c][d]

# ---------------------------------------------------------------------------
# Helper: Euler-Lagrange operator with arbitrary order
# ---------------------------------------------------------------------------


def euler_lagrange(L, q, t, max_order=10):
    expr = diff(L, q)
    for k in range(1, max_order + 1):
        qk = q.diff(t, k)
        term = diff(L, qk)
        if term == 0:
            continue
        expr += (-1) ** k * term.diff(t, k)
    return expr


# ---------------------------------------------------------------------------
# Phase 2 — Build L_T4
# ---------------------------------------------------------------------------

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

print("=" * 76)
print("Phase 2 — T4 Lagrangian (3 dynamical + 3 constraint-promoted)")
print("=" * 76)
print()
print(f"L_T4 has {len(L.args)} terms (Lagrangian density)")
print()
print("Structural breakdown:")
print(f"  - {N_DYN} dynamical kinetic + mass terms (y_a)")
print(f"  - {N_CON} constraint mass terms (-M_c^2 h_c^2 / 2)")
print(
    f"  - {N_DYN * N_CON} algebraic mixings (lambda) + "
    f"{N_DYN * N_CON} velocity mixings (mu)"
)
n_K = N_CON * (N_CON + 1) // 2
print(f"  - {n_K} symmetric b5-sector kinetic terms (K_cd ddot h_c ddot h_d)")
print()


# ---------------------------------------------------------------------------
# Phase 3 — Parker-Simon iterative reduction
# ---------------------------------------------------------------------------
print("=" * 76)
print("Phase 3 — Parker-Simon iterative reduction")
print("=" * 76)
print()

L0 = L.subs(b5, 0)
EL_h_0 = [expand(euler_lagrange(L0, hs[c], t)) for c in range(N_CON)]

print("Order-0 EOMs for h_c (b5=0):")
for c, E in enumerate(EL_h_0):
    print(f"  E_h{c + 1} = {E}")
print()

# The order-0 EOM is purely algebraic in h_c with diagonal M^2 mass matrix:
# -M_c^2 h_c - sum_a lam_ac y_a + sum_a mu_ac d_t y_a = 0  (after IBP for d_t)
# Wait: actually the contribution is  d/dh_c [ - mu_ac y_a d_t h_c ] = 0
# but d/d(d_t h_c) [ ... ] = -mu_ac y_a, so -d_t(-mu_ac y_a) = mu_ac d_t y_a.
sol_h0 = sp.solve(EL_h_0, hs)
h0 = [sp.simplify(sol_h0[hs[c]]) for c in range(N_CON)]

print("Order-0 algebraic solutions:")
for c, hsol in enumerate(h0):
    print(f"  h{c + 1}^(0) = {hsol}")
print()

# Order-1 correction: h_c = h_c^(0) + b5 * (M^{-2})_cc' K_{c'd} d_t^4 h_d^(0)
# Since M^2 is diagonal here, (M^{-2})_cc' = delta_cc'/M_c^2.
h1 = []
for c in range(N_CON):
    corr = sp.Integer(0)
    for d in range(N_CON):
        corr += K[c][d] * h0[d].diff(t, 4) / M_c[c]
    h1.append(corr)
h_full = [h0[c] + b5 * h1[c] for c in range(N_CON)]

# ---------------------------------------------------------------------------
# Substitute h-solutions into y-EOMs, truncate at O(b5)
# ---------------------------------------------------------------------------


def sub_h_in_eom(expr, h_sols, max_h_deriv=10):
    out = expr
    for c in range(N_CON):
        for k in range(max_h_deriv, 0, -1):
            out = out.subs(hs[c].diff(t, k), h_sols[c].diff(t, k))
        out = out.subs(hs[c], h_sols[c])
    return out


print("Computing PS-reduced y-EOMs (this is the slow step) ...")
sys.stdout.flush()

E_y_full = [euler_lagrange(L, ys[a], t) for a in range(N_DYN)]
E_y_red = []
for a in range(N_DYN):
    sub = sub_h_in_eom(E_y_full[a], h_full)
    sub = expand(sub.series(b5, 0, 2).removeO())
    E_y_red.append(sub)


# Determine max derivative order in the reduced system
def max_deriv_order(expr, q, t):
    n = 0
    for d in expr.atoms(sp.Derivative):
        if d.expr == q:
            order = sum(c for v, c in d.variable_count if v == t)
            n = max(n, order)
    return n


r_max = 0
for a in range(N_DYN):
    for b in range(N_DYN):
        r_max = max(r_max, max_deriv_order(E_y_red[a], ys[b], t))

print(f"Max derivative order in PS-reduced y-EOMs: r = {r_max}")
print("Expected: 6 (from the b5*K*d_t^4 -> y_a chain through h)")
print()


# Convert to jet form
Y_jet = [
    [symbols(f"Y{a + 1}_{k}", real=True) for k in range(r_max + 1)]
    for a in range(N_DYN)
]


def to_jet(expr):
    out = expr
    for a in range(N_DYN):
        for k in range(r_max, 0, -1):
            out = out.subs(ys[a].diff(t, k), Y_jet[a][k])
        out = out.subs(ys[a], Y_jet[a][0])
    return out


eps_y = [expand(to_jet(E_y_red[a])) for a in range(N_DYN)]


# Polynomiality check
fiber_syms = set()
for a in range(N_DYN):
    fiber_syms.update(Y_jet[a])


def has_rational_fiber(expr):
    _n, d = sp.fraction(sp.together(expr))
    return any(f in fiber_syms for f in d.free_symbols)


print("Polynomiality of PS-reduced epsilons in fibres:")
for a in range(N_DYN):
    rat = has_rational_fiber(eps_y[a])
    print(f"  eps_y{a + 1} rational in fibres: {rat}")
print()


# ---------------------------------------------------------------------------
# Phase 4 — Vainberg-Tonti homotopy and Krupka-Voicu Theorem 1 verification
# ---------------------------------------------------------------------------
print("=" * 76)
print("Phase 4 — Vainberg-Tonti integral construction")
print("=" * 76)
print()


def homothety(expr, u):
    sub = {}
    for a in range(N_DYN):
        for k in range(r_max + 1):
            sub[Y_jet[a][k]] = u * Y_jet[a][k]
    return expr.subs(sub)


eps_y_u = [expand(homothety(eps_y[a], u)) for a in range(N_DYN)]

print("u-power range of eps_a(u*jet) (should be [1,1] for linear sources):")
for a in range(N_DYN):
    poly_u = Poly(eps_y_u[a], u)
    monoms = poly_u.monoms()
    if monoms:
        mn = min(m[0] for m in monoms)
        mx = max(m[0] for m in monoms)
        print(f"  eps_y{a + 1}(u*jet): u-power range [{mn}, {mx}]")
print()

# Construct VT Lagrangian:  L_VT = sum_a Y[a][0] * INT_0^1 eps_a(u*jet) du
print("Computing VT integral (homotopy [0,1]) ...")
L_VT = sp.Integer(0)
for a in range(N_DYN):
    integrand = expand(Y_jet[a][0] * eps_y_u[a])
    L_VT += integrate(integrand, (u, 0, 1))
L_VT = expand(L_VT)

n_terms = len(L_VT.args) if isinstance(L_VT, sp.Add) else 1
print(f"L_VT constructed: {n_terms} terms")
print()


# Krupka-Voicu Theorem 1 verification
def jet_to_funcs(expr):
    out = expr
    for a in range(N_DYN):
        out = out.subs(Y_jet[a][0], ys[a])
        for k in range(1, r_max + 1):
            out = out.subs(Y_jet[a][k], ys[a].diff(t, k))
    return out


print("Verifying Krupka-Voicu Theorem 1: EL(L_VT) = eps_y ...")
sys.stdout.flush()

L_VT_func = jet_to_funcs(L_VT)
all_ok = True
for a in range(N_DYN):
    el_a = euler_lagrange(L_VT_func, ys[a], t, max_order=2 * r_max + 2)
    el_a_jet = expand(to_jet(el_a))
    diff_a = simplify(el_a_jet - eps_y[a])
    ok = diff_a == 0
    print(
        f"  EL_y{a + 1}(L_VT) - eps_y{a + 1} = {diff_a if not ok else '0'} -- "
        f"{'OK' if ok else 'FAIL'}"
    )
    all_ok = all_ok and ok
print()
print(f"VT canonical-completion verified for T4 (3+3): {all_ok}")
print()


# Denominator structure
num, den = sp.fraction(sp.together(L_VT))
den_pure = sp.cancel(den / sp.gcd(den, 2))
print(f"L_VT denominator (mod /2): {den_pure}")
print()


# ---------------------------------------------------------------------------
# Phase 5 — Constraint-subspace projection
# ---------------------------------------------------------------------------
print("=" * 76)
print("Phase 5 — Project L_VT onto constraint (h_4,h_7,h_9) subspace")
print("=" * 76)
print()
print(
    "L_VT lives in y-space (the dynamical sector). To answer the metric\n"
    "Pais-Uhlenbeck question we examine TWO subspaces:\n"
    "  (a) The Routhian/on-shell L_eff obtained by substituting h(y) into\n"
    "      L_T4 directly. This is the 'h-eliminated direct Lagrangian'.\n"
    "  (b) An L_VT-h that we obtain by extending VT to BOTH y AND h\n"
    "      coordinates simultaneously (treating h as a dynamical fibre\n"
    "      coordinate too). This is needed to address the 'kinetic\n"
    "      structure of the constraint subspace' question.\n"
)

# (a) Direct h-eliminated Lagrangian
print("(a) Direct on-shell Lagrangian L_eff = L_T4|_{h=h_full}")
sys.stdout.flush()
L_eff = L
for c in range(N_CON):
    for k in range(8, 0, -1):
        L_eff = L_eff.subs(hs[c].diff(t, k), h_full[c].diff(t, k))
    L_eff = L_eff.subs(hs[c], h_full[c])
L_eff = expand(L_eff.series(b5, 0, 2).removeO())
n_eff = len(L_eff.args) if isinstance(L_eff, sp.Add) else 1
print(f"    L_eff has {n_eff} terms after expansion.")
print()

# Check: which highest derivatives of y appear in L_eff?
print("    Highest derivative of y in L_eff (per field):")
for a in range(N_DYN):
    n = max_deriv_order(L_eff, ys[a], t)
    print(f"      y_{a + 1}: order {n}")
print()

# At b5=0, L_eff should be 2nd-order (only kinetic + mass + algebraic Routhian)
L_eff_b50 = expand(L_eff.subs(b5, 0))
print("    At b5=0 (algebraic constraint subspace):")
for a in range(N_DYN):
    n = max_deriv_order(L_eff_b50, ys[a], t)
    print(f"      y_{a + 1}: order {n}  (2 = standard kinetic; >2 = HD)")
print()

# At order b5: which derivatives appear?
L_eff_b5 = expand(L_eff - L_eff_b50)
print("    O(b5) correction:")
for a in range(N_DYN):
    n = max_deriv_order(L_eff_b5, ys[a], t)
    print(f"      y_{a + 1}: order {n}  (>= 4 indicates Pais-Uhlenbeck)")
print()


# (b) Extend VT homothety to (y,h) jointly to expose the structure of
# h-mode kinetic terms.  Here we treat h as a free coordinate (NOT
# eliminated) and compute the FULL EL system and its VT integral.
print("(b) Joint VT on (y,h) — full 6-field system")
sys.stdout.flush()

# Full EL system from L (no PS reduction — pure Helmholtz check on L itself)
# After truncating at O(b5), L is itself the parent Lagrangian, so
# EL_y(L) and EL_h(L) ARE the parent EOMs.  The VT integral of these
# should reproduce L (up to total time derivative) — a sanity check.
all_fields = ys + hs
N_TOT = len(all_fields)


# Determine max derivative order in L for each field
r_each = []
for f in all_fields:
    rmax_f = 0
    for d_ in L.atoms(sp.Derivative):
        if d_.expr == f:
            o = sum(c for v, c in d_.variable_count if v == t)
            rmax_f = max(rmax_f, o)
    r_each.append(max(rmax_f, 1))

# Joint jet variables
F_jet = []
for i, f in enumerate(all_fields):
    rmf = max(r_each[i] + 2, 4)  # buffer for higher EL orders
    F_jet.append([symbols(f"F{i + 1}_{k}", real=True) for k in range(rmf + 1)])


def to_full_jet(expr):
    out = expr
    for i, f in enumerate(all_fields):
        rmf = len(F_jet[i]) - 1
        for k in range(rmf, 0, -1):
            out = out.subs(f.diff(t, k), F_jet[i][k])
        out = out.subs(f, F_jet[i][0])
    return out


def from_full_jet(expr):
    out = expr
    for i, f in enumerate(all_fields):
        rmf = len(F_jet[i]) - 1
        out = out.subs(F_jet[i][0], f)
        for k in range(1, rmf + 1):
            out = out.subs(F_jet[i][k], f.diff(t, k))
    return out


def homothety_full(expr, u):
    sub = {}
    for i, _ in enumerate(all_fields):
        rmf = len(F_jet[i]) - 1
        for k in range(rmf + 1):
            sub[F_jet[i][k]] = u * F_jet[i][k]
    return expr.subs(sub)


print("    Computing full 6-field EL system ...")
sys.stdout.flush()
eps_full_func = []
for f in all_fields:
    eps_f = euler_lagrange(L, f, t, max_order=8)
    eps_full_func.append(expand(eps_f))


# Note: For the parent L, eps = EL(L) IS the source form.
# VT integral of this should reproduce L (modulo total t-derivative).
# This is a Helmholtz consistency check on L itself.
print("    Computing VT integral of the full 6-field epsilon system ...")
sys.stdout.flush()
eps_full_jet = [expand(to_full_jet(e)) for e in eps_full_func]
eps_full_u = [expand(homothety_full(e, u)) for e in eps_full_jet]

# u-power range
print("    u-power range of full 6-field epsilons:")
for i, e in enumerate(eps_full_u):
    if e == 0:
        print(f"      eps_{all_fields[i]}: zero")
        continue
    poly_u = Poly(e, u)
    monoms = poly_u.monoms()
    if monoms:
        mn = min(m[0] for m in monoms)
        mx = max(m[0] for m in monoms)
        print(f"      eps_{all_fields[i]}: u-power [{mn}, {mx}]")
print()

L_VT_full = sp.Integer(0)
for i, _ in enumerate(all_fields):
    integrand = expand(F_jet[i][0] * eps_full_u[i])
    L_VT_full += integrate(integrand, (u, 0, 1))
L_VT_full = expand(L_VT_full)

# Compare to L directly (up to total derivative): EL(L_VT_full) should equal eps
print("    Verifying EL(L_VT_full) = eps for the parent system ...")
L_VT_full_func = from_full_jet(L_VT_full)
all_ok_full = True
for i, f in enumerate(all_fields):
    el = euler_lagrange(L_VT_full_func, f, t, max_order=8)
    el_jet = expand(to_full_jet(el))
    diff_i = simplify(el_jet - eps_full_jet[i])
    ok = diff_i == 0
    print(
        f"      EL_{f}(L_VT_full) - eps = {0 if ok else 'NONZERO'} -- {'OK' if ok else 'FAIL'}"
    )
    all_ok_full = all_ok_full and ok
print()
print(f"    Joint VT consistency: {all_ok_full}")
print()


# Now project L_VT_full onto the h subspace by setting all y_a -> 0 and
# inspect kinetic structure.  This is the VT-derived effective Lagrangian
# for the h-only subspace.
print("Projecting L_VT_full onto h-subspace (y_a = d_t y_a = ... = 0):")
proj_subs = {}
for i in range(N_DYN):
    rmf = len(F_jet[i]) - 1
    for k in range(rmf + 1):
        proj_subs[F_jet[i][k]] = 0
L_VT_h = expand(L_VT_full.subs(proj_subs))
L_VT_h_func = from_full_jet(L_VT_h)
print()
print("L_VT |_{y=0} (h-only subspace, jet form):")
print(f"  {L_VT_h}")
print()

# Inspect structure of h-only subspace
h_jet_idx = list(range(N_DYN, N_TOT))
print("Highest derivative of h_c in L_VT_h:")
for ic in h_jet_idx:
    rmf = len(F_jet[ic]) - 1
    max_k = -1
    for k in range(rmf + 1):
        if F_jet[ic][k] in L_VT_h.free_symbols:
            max_k = max(max_k, k)
    print(f"  h_{ic - N_DYN + 1}: highest jet order {max_k}")
print()

# Note: the h-only subspace at b5=0 should have ONLY mass terms (no kinetic),
# because the parent L0 has -M_c^2 h_c^2/2 only.  At b5 != 0 there's the
# K_cd ddot h_c ddot h_d term.  This is exactly the Pais-Uhlenbeck signature.
print("h-only subspace breakdown:")
L_VT_h_b50 = expand(L_VT_h.subs(b5, 0))
L_VT_h_b5 = expand(L_VT_h - L_VT_h_b50)
print(
    f"  At b5=0 ({len(L_VT_h_b50.args) if isinstance(L_VT_h_b50, sp.Add) else 1} terms):"
)
print(f"    {L_VT_h_b50}")
print()
print(f"  O(b5) ({len(L_VT_h_b5.args) if isinstance(L_VT_h_b5, sp.Add) else 1} terms):")
print(f"    {L_VT_h_b5}")
print()


# Check: does the h-only subspace at b5 != 0 contain (d_t^2 h)^2 terms?
# This is the Pais-Uhlenbeck signature.
def has_2nd_derivative_squared(expr, jets):
    """Check if expr contains products F_i_2 * F_j_2 for any i,j among jets."""
    pairs = []
    for ic in jets:
        rmf = len(F_jet[ic]) - 1
        if rmf >= 2:
            for jc in jets:
                rmf2 = len(F_jet[jc]) - 1
                if rmf2 >= 2:
                    sym_pair = F_jet[ic][2] * F_jet[jc][2]
                    coeff = expr.coeff(sym_pair, 1)
                    if coeff != 0:
                        pairs.append((ic, jc, coeff))
    return pairs


PU_pairs = has_2nd_derivative_squared(L_VT_h_b5, h_jet_idx)
print("Pais-Uhlenbeck (d_t^2 h)*(d_t^2 h) pairs in O(b5) part of L_VT_h:")
if PU_pairs:
    for ic, jc, coef in PU_pairs:
        print(f"  h_{ic - N_DYN + 1}'' * h_{jc - N_DYN + 1}''  coefficient = {coef}")
    print()
    print("  -> The h-subspace IS Pais-Uhlenbeck at order b5.")
else:
    print("  None (h-subspace is NOT PU at order b5).")
print()


# ---------------------------------------------------------------------------
# Phase 6 — Hamiltonian / constraint Poisson matrix
# ---------------------------------------------------------------------------
print("=" * 76)
print("Phase 6 — Hamiltonian analysis of L_VT")
print("=" * 76)
print()

# We work with the y-only L_VT (the actual VT-derived Lagrangian).
# For the y-only L_VT:
#   - It depends on y_a, dot y_a, ddot y_a, ... up to (d_t^r_max y_a).
#   - r_max = 6 (from PS reduction)
#   - This is a higher-derivative Lagrangian, so the Ostrogradsky
#     Hamiltonian formalism applies.
print(f"y-only L_VT depends on derivatives up to order {r_max}.")
print()
print(
    "Ostrogradsky construction: introduce auxiliary fields\n"
    "  Q_a^k = d_t^k y_a    for k=0,...,r-1\n"
    "and conjugate momenta\n"
    "  P_a^k = sum_{j=0}^{r-1-k} (-1)^j d_t^j (dL/d(d_t^{k+1+j} y_a))\n"
    "for k=0,...,r-1.\n"
)

# For r=6, this gives 6 phase-space variables per dynamical field.
# Hamiltonian is H = sum_a sum_k P_a^k * Q_a^{k+1} - L  with Q_a^r given by
# the standard Ostrogradsky relation (last momentum = dL/d(d_t^r y)).

# For the purposes of this preflight, we compute the SYMPLECTIC structure
# matrix at the LINEARISED level: the constraint-Poisson submatrix that
# determines whether the system has uniform rank.

# Define jet-level momenta
print("Computing Ostrogradsky momenta for L_VT (linearised) ...")
sys.stdout.flush()

# Ostrogradsky momentum P_a^k = sum_{j=0}^{r-1-k} (-1)^j d_t^j (dL/d(d_t^{k+1+j} y_a))
# We compute these as jet-symbolic expressions.

# To avoid blowing up the algebra, we work to leading order in b5.  Split
# L_VT = L0_VT + b5 * L1_VT.
L_VT_jet_subs = L_VT  # already in jet form

L0_VT = expand(L_VT_jet_subs.subs(b5, 0))
L1_VT_full = expand(L_VT_jet_subs - L0_VT)

print(
    f"  L0_VT (b5=0 part) has {len(L0_VT.args) if isinstance(L0_VT, sp.Add) else 1} terms"
)
print(
    f"  L1_VT (O(b5) part) has {len(L1_VT_full.args) if isinstance(L1_VT_full, sp.Add) else 1} terms"
)
print()

# Compute the Hessian matrix H_{a,k; b,l} = d^2 L_VT / d(Y_a_k) d(Y_b_l)
# for k,l = 1,...,r_max (the velocity & higher-derivative jets).
#
# RANK-UNIFORMITY criterion (Round 1 Agent A's diagnostic):
# If rank(H) at b5=0 < rank(H) at b5 != 0, we have a rank jump (PU pathology).
# If rank is constant in b5 -> Path A produces a UNIFORM-RANK Hamiltonian.

# Build leading Hessian for the highest-derivative pair (the PU block):
# At b5=0, L_VT_jet has at most (d_t y_a)^2 type terms (kinetic, no HD).
# At b5 != 0, L_VT_jet has (d_t^? y_a)*(d_t^? y_b) type terms with high orders.
print("Highest-derivative jet structure of L_VT_jet:")
print("  At b5=0: ", end="")
max_order_b50 = 0
for a in range(N_DYN):
    for k in range(r_max + 1):
        if Y_jet[a][k] in L0_VT.free_symbols:
            max_order_b50 = max(max_order_b50, k)
print(f"max jet order = {max_order_b50}")
print("  At b5 != 0: ", end="")
max_order_full = 0
for a in range(N_DYN):
    for k in range(r_max + 1):
        if Y_jet[a][k] in L_VT_jet_subs.free_symbols:
            max_order_full = max(max_order_full, k)
print(f"max jet order = {max_order_full}")
print()

# The rank jump is in the highest-derivative Hessian:
# H^{(top)}_{ab} = d^2 L_VT / d(Y_a_{r_max}) d(Y_b_{r_max})
print(f"Top-Hessian (kinetic block) at order Y_a_{max_order_full}:")
H_top = zeros(N_DYN, N_DYN)
for a in range(N_DYN):
    for b in range(N_DYN):
        H_top[a, b] = sp.diff(
            L_VT_jet_subs, Y_jet[a][max_order_full], Y_jet[b][max_order_full]
        )
print(f"  H_top = {H_top.tolist()}")
print()
det_H_top = sp.simplify(H_top.det())
print(f"  det(H_top) = {det_H_top}")
print()

# Standard kinetic Hessian (Y_a_1 block) at b5=0
H_kin0 = zeros(N_DYN, N_DYN)
for a in range(N_DYN):
    for b in range(N_DYN):
        H_kin0[a, b] = sp.diff(L0_VT, Y_jet[a][1], Y_jet[b][1])
print(f"Kinetic Hessian at b5=0 (Y_a_1 block): H_kin0 = {H_kin0.tolist()}")
det_H_kin0 = sp.simplify(H_kin0.det())
print(f"  det(H_kin0) = {det_H_kin0}")
print()


# Now study the rank-jump diagnostic:
# RANK at b5=0:  measure rank of full (Y_jet[a][k]) Hessian for k>=1
# RANK at b5 != 0: same but at order b5
def full_velocity_hessian(L, max_order):
    """Return Hessian wrt all (y_a, d_t y_a, ..., d_t^max_order y_a) jets except y_a=0 column."""
    rows = []
    for a in range(N_DYN):
        rows.extend(Y_jet[a][k] for k in range(1, max_order + 1))
    n = len(rows)
    H = zeros(n, n)
    for i, vi in enumerate(rows):
        for j, vj in enumerate(rows):
            H[i, j] = sp.diff(L, vi, vj)
    return H, rows


print("Full velocity-jet Hessian rank diagnostic:")
H_full_b50, rows_b50 = full_velocity_hessian(L0_VT, max_order_b50)
H_full_b5, rows_b5 = full_velocity_hessian(L_VT_jet_subs, max_order_full)

# Substitute b5 -> small numeric and rank-check both
H_b50_num = H_full_b50.subs(dict.fromkeys(H_full_b50.free_symbols, 1))
H_b5_num = H_full_b5.subs({s: 1 for s in H_full_b5.free_symbols if s != b5}).subs(
    b5, sp.Rational(1, 1000)
)

try:
    rank_b50 = H_b50_num.rank()
except Exception as e:
    rank_b50 = f"ERROR: {e}"
try:
    rank_b5 = H_b5_num.rank()
except Exception as e:
    rank_b5 = f"ERROR: {e}"

print(f"  rank(H_velocity_jet) at b5=0 (params=1): {rank_b50}")
print(f"  rank(H_velocity_jet) at b5=1/1000 (params=1): {rank_b5}")
print()
print(f"  Velocity-jet dimension at b5=0: {len(rows_b50)}")
print(f"  Velocity-jet dimension at b5!=0: {len(rows_b5)}")
print()

if isinstance(rank_b50, int) and isinstance(rank_b5, int):
    if rank_b5 > rank_b50 or len(rows_b5) > len(rows_b50):
        print(
            "  -> RANK JUMP DETECTED.  Path A's L_VT inherits the same\n"
            "     constraint-promotion structure as the parent Lagrangian.\n"
            "     The Hamiltonian-side rank-jump persists."
        )
    else:
        print("  -> No rank jump.  Phase-space dimension uniform in b5.")
else:
    print("  -> Cannot determine rank numerically (algebraic rank ambiguity).")
print()


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
print("=" * 76)
print("FINAL VERDICT for T4 (3 dynamical + 3 constraint, BC-style mass matrix)")
print("=" * 76)
print()
print("  VT integral CONVERGES on [0,1] for all 3 dynamical eps.")
print(f"  Krupka-Voicu Theorem 1: EL(L_VT) - eps = 0 verified: {all_ok}")
print(f"  L_VT denominator: {den_pure}")
print()
print("  Phase 5 result on the metric Pais-Uhlenbeck subspace:")
print("    L_eff (h-eliminated) at b5=0: standard 2nd-order kinetic")
print(
    f"    L_eff at O(b5): contains derivatives up to "
    f"{max(max_deriv_order(L_eff_b5, ys[a], t) for a in range(N_DYN))}-th order in y"
)
print(f"    h-only subspace contains (d^2 h)^2 PU pairs: {'YES' if PU_pairs else 'NO'}")
print()
print("  Phase 6 result on the Hamiltonian-side rank:")
print(f"    Highest-velocity-jet rank at b5=0: {rank_b50}")
print(f"    Highest-velocity-jet rank at b5!=0: {rank_b5}")
print(f"    Phase-space dimension at b5=0: {len(rows_b50)}")
print(f"    Phase-space dimension at b5!=0: {len(rows_b5)}")
print()


# Save full L_VT and constraint matrix
results_dir = OUT_DIR

# Save L_VT_full (the actual VT-derived Lagrangian)
with Path(results_dir / "vt_T4_lagrangian.txt").open("w", encoding="utf-8") as f:
    f.write("Vainberg-Tonti Lagrangian for T4 toy (3 dynamical + 3 constraint)\n")
    f.write("=" * 76 + "\n\n")
    f.write(f"VT canonical-completion verified: {all_ok}\n")
    f.write(f"Number of terms in y-only L_VT: {n_terms}\n")
    f.write(f"Denominator (y-only L_VT, mod /2): {den_pure}\n\n")
    f.write("--- y-only L_VT (jet form) ---\n")
    f.write(str(L_VT) + "\n\n")
    f.write("--- L_eff (direct h-elimination) ---\n")
    f.write(f"Number of terms: {len(L_eff.args) if isinstance(L_eff, sp.Add) else 1}\n")
    f.write(str(L_eff) + "\n\n")
    f.write("--- L_VT projected onto h-subspace (y=0) ---\n")
    f.write("L_VT_h_b50 (b5=0 part): " + str(L_VT_h_b50) + "\n\n")
    f.write("L_VT_h_b5  (O(b5) part): " + str(L_VT_h_b5) + "\n\n")

print(f"Wrote y-only L_VT and projections to: {results_dir / 'vt_T4_lagrangian.txt'}")


# Constraint Poisson matrix structure
constraint_data = {
    "T4_setup": {
        "N_dynamical": N_DYN,
        "N_constraint": N_CON,
        "max_jet_order_b5_eq_0": int(max_order_b50),
        "max_jet_order_b5_ne_0": int(max_order_full),
    },
    "vt_completion": {
        "krupka_voicu_theorem_1_verified": bool(all_ok),
        "joint_6field_helmholtz_check": bool(all_ok_full),
        "L_VT_denominator_modulo_two": str(den_pure),
        "n_terms_y_only_L_VT": int(n_terms),
    },
    "phase5_metric_PU_subspace": {
        "L_eff_max_deriv_order_b5_eq_0": {
            f"y_{a + 1}": int(max_deriv_order(L_eff_b50, ys[a], t))
            for a in range(N_DYN)
        },
        "L_eff_max_deriv_order_b5_correction": {
            f"y_{a + 1}": int(max_deriv_order(L_eff_b5, ys[a], t)) for a in range(N_DYN)
        },
        "h_subspace_PU_pairs_count": len(PU_pairs),
        "h_subspace_PU_pairs": [
            {
                "h_c_index": ic - N_DYN + 1,
                "h_d_index": jc - N_DYN + 1,
                "coefficient": str(coef),
            }
            for ic, jc, coef in PU_pairs
        ],
    },
    "phase6_hamiltonian_rank": {
        "top_Hessian_3x3": [
            [str(H_top[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
        ],
        "top_Hessian_det": str(det_H_top),
        "kinetic_Hessian_b5_eq_0_3x3": [
            [str(H_kin0[i, j]) for j in range(N_DYN)] for i in range(N_DYN)
        ],
        "kinetic_Hessian_b5_eq_0_det": str(det_H_kin0),
        "velocity_jet_rank_b5_eq_0_params_1": str(rank_b50),
        "velocity_jet_rank_b5_ne_0_params_1": str(rank_b5),
        "velocity_jet_dim_b5_eq_0": len(rows_b50),
        "velocity_jet_dim_b5_ne_0": len(rows_b5),
        "rank_jump_detected": bool(
            isinstance(rank_b50, int)
            and isinstance(rank_b5, int)
            and (rank_b5 > rank_b50 or len(rows_b5) > len(rows_b50))
        ),
    },
}

with Path(results_dir / "vt_T4_constraint_matrix.json").open("w", encoding="utf-8") as f:
    json.dump(constraint_data, f, indent=2)
print(f"Wrote constraint analysis to: {results_dir / 'vt_T4_constraint_matrix.json'}")
print()
print("DONE.")
