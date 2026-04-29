# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
Vainberg-Tonti integral convergence test for T3 toy.
====================================================

T3 is a 3+3 PGT-faithful toy:

  - 3 dynamical fields q_a (a=1,2,3) representing TT graviton modes
  - 3 constraint-promoted fields h_b (b=1,2,3) representing h_4, h_7, h_9
  - Cross-couplings via lambda_{ab} q_a h_b   (algebraic mixing)
  - Velocity-mixing via mu_{ab} q_a d_t h_b   (derivative mixing)
  - PGT b5 R_tilde^2 -> b5 K_{ab} (d_t^2 h_a)(d_t^2 h_b)  with K_{ab}
    a non-trivial cross-coupling matrix (not just diagonal).

Lagrangian:
  L_T3 = sum_a (1/2)(d_t q_a)^2 - (1/2) m_a^2 q_a^2
       - sum_{a,b} ( lambda_{ab} q_a h_b + mu_{ab} q_a d_t h_b )
       + sum_a (-1/2) M_a^2 h_a^2
       + (b5/2) sum_{a,b} K_{ab} (d_t^2 h_a)(d_t^2 h_b)

This activates the constraint-field kinetic structure ONLY via b5 (so b5=0
is the algebraic-constraint surface).  Mass terms M_a^2 give h_a^(0)
algebraic solutions.

We follow the same PS-reduction recipe as T2:
  h_a^(0) = (M^2)^{-1}_{ab} ( - lambda_{cb} q_c + mu_{cb} d_t q_c )
  h_a    = h_a^(0) + b5 * sum_b (M^2)^{-1}_{ab} K_{bc} d_t^4 h_c^(0)

Substitute into d/dq_a EOMs and truncate at O(b5).  Compute VT integral.

For analytical tractability, take 2 dynamical + 2 constraint fields with
generic 2x2 lambda, mu, M^2, K matrices.  This already exercises the
matrix structure that motivated b5 R_tilde^2 cross-couplings in PGT.
"""

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
)

t = symbols("t")
u = symbols("u", positive=True)
b5 = symbols("b5", real=True)

# --- 2 dynamical, 2 constraint version (analytically tractable) ---
N = 2  # try 2; can bump to 3 by parameter

# Dynamical fields q_a, constraint fields h_b
qs = [Function(f"q{a + 1}")(t) for a in range(N)]
hs = [Function(f"h{a + 1}")(t) for a in range(N)]

# Parameters: masses, couplings
m = [symbols(f"m{a + 1}", real=True, positive=True) for a in range(N)]
M2 = [symbols(f"M2_{a + 1}", real=True, positive=True) for a in range(N)]
lam_ab = [[symbols(f"L{a + 1}{b + 1}", real=True) for b in range(N)] for a in range(N)]
mu_ab = [[symbols(f"mu{a + 1}{b + 1}", real=True) for b in range(N)] for a in range(N)]
# Symmetric kinetic-cross-coupling matrix K_{ab} (b5 sector)
K_ab = [[symbols(f"K{a + 1}{b + 1}", real=True) for b in range(N)] for a in range(N)]
# enforce K symmetric
for a in range(N):
    for b in range(a + 1, N):
        K_ab[b][a] = K_ab[a][b]


def euler_lagrange(L, q, t, max_order=8):
    expr = diff(L, q)
    for k in range(1, max_order + 1):
        qk = q.diff(t, k)
        term = diff(L, qk)
        if term == 0:
            continue
        expr += (-1) ** k * term.diff(t, k)
    return expr


# --- Build Lagrangian ---
L = 0
for a in range(N):
    L += Rational(1, 2) * qs[a].diff(t) ** 2 - Rational(1, 2) * m[a] ** 2 * qs[a] ** 2
for a in range(N):
    L += -Rational(1, 2) * M2[a] * hs[a] ** 2
for a in range(N):
    for b in range(N):
        L += -lam_ab[a][b] * qs[a] * hs[b]
        L += -mu_ab[a][b] * qs[a] * hs[b].diff(t)
for a in range(N):
    for b in range(N):
        L += Rational(1, 2) * b5 * K_ab[a][b] * hs[a].diff(t, 2) * hs[b].diff(t, 2)

print("=" * 72)
print("Step 1: T3 Lagrangian and order-0 algebraic constraint")
print("=" * 72)
print(f"L_T3 (N={N} dyn + {N} constr):")
print("  ", expand(L))
print()

# --- Order-0: solve algebraic h-EOMs at b5=0 ---
L0 = L.subs(b5, 0)
EL_h_0 = [expand(euler_lagrange(L0, hs[a], t)) for a in range(N)]
print("Order-0 EOMs for h_a (b5=0):")
for a, E in enumerate(EL_h_0):
    print(f"  E_h{a + 1} =", E)
print()

# Solve for hs[0], hs[1] as linear system
# (the system is purely algebraic at b5=0 because no h-derivatives appear
# without a b5 prefactor: derivatives appear only via mu_{ab} q_a d_t h_b
# which after E-L gives -d_t(mu_{ab} q_a) * delta-like structure... wait,
# d/dh_b of mu*q*d_th_b = 0; d/d(d_th_b) of that = mu*q.  So the E-L for h_b
# is: -M2*h_b - sum_a lambda_{ab} q_a + d_t(sum_a mu_{ab} q_a) + b5 * higher.
# Hence d/dt mu*q is the contribution -- still algebraic in h.)
sol_h0 = sp.solve(EL_h_0, hs)
print("Order-0 solutions h_a^(0):")
h0_list = [sp.simplify(sol_h0[hs[a]]) for a in range(N)]
for a, hsol in enumerate(h0_list):
    print(f"  h{a + 1}^(0) =", hsol)
print()

# --- Order-1 correction: h_a = h_a^(0) + b5 * sum_b (M2)^{-1}_{ab} K_{bc} d_t^4 h_c^(0)
# (Diagonal M^2 -> M2_{ab}^{-1} = delta_{ab}/M2_a)
h1_list = []
for a in range(N):
    corr = 0
    for c in range(N):
        corr += K_ab[a][c] * h0_list[c].diff(t, 4) / M2[a]
    h1_list.append(corr)
h_sol = [h0_list[a] + b5 * h1_list[a] for a in range(N)]
print("Iterative solutions h_a = h_a^(0) + b5 * h_a^(1) + O(b5^2):")
for a, hs_ in enumerate(h_sol):
    print(f"  h{a + 1} =", expand(hs_))
print()


# --- Substitute into q_a EOMs, truncate at O(b5) ---
def sub_h_in_eom(expr, h_sols, max_h_deriv=8):
    out = expr
    for a in range(N):
        for k in range(max_h_deriv, 0, -1):
            out = out.subs(hs[a].diff(t, k), h_sols[a].diff(t, k))
        out = out.subs(hs[a], h_sols[a])
    return out


E_q_full = [euler_lagrange(L, qs[a], t) for a in range(N)]
E_q_red = []
for a in range(N):
    sub = sub_h_in_eom(E_q_full[a], h_sol)
    sub = expand(sub.series(b5, 0, 2).removeO())
    E_q_red.append(sub)

print("=" * 72)
print("Step 2: PS-reduced EOMs E^(1)_qa truncated at O(b5)")
print("=" * 72)
for a in range(N):
    print(f"E^(1)_q{a + 1} =", E_q_red[a])
    print()


# --- Convert to jet form and check polynomiality ---
def max_deriv_order(expr, q, t):
    n = 0
    for d in expr.atoms(sp.Derivative):
        if d.expr == q:
            order = sum(c for v, c in d.variable_count if v == t)
            n = max(n, order)
    return n


r_max = 0
for a in range(N):
    for b in range(N):
        r_max = max(r_max, max_deriv_order(E_q_red[a], qs[b], t))
print(f"Max derivative order in PS-reduced EOMs: r = {r_max}")
print()

Q_jet = [
    [symbols(f"Q{a + 1}_{k}", real=True) for k in range(r_max + 1)] for a in range(N)
]


def to_jet(expr):
    out = expr
    for a in range(N):
        for k in range(r_max, 0, -1):
            out = out.subs(qs[a].diff(t, k), Q_jet[a][k])
        out = out.subs(qs[a], Q_jet[a][0])
    return out


E_q_jet = [expand(to_jet(E_q_red[a])) for a in range(N)]
print("=" * 72)
print("Step 3: jet-form EOMs")
print("=" * 72)
for a in range(N):
    print(f"eps_q{a + 1}(jet) =", E_q_jet[a])
    print()


# Check polynomiality in fiber variables
def has_rational_fiber(expr, fiber_syms) -> bool:
    _n, d = sp.fraction(sp.together(expr))
    return any(f in fiber_syms for f in d.free_symbols)


fiber_syms = set()
for a in range(N):
    fiber_syms.update(Q_jet[a][k] for k in range(r_max + 1))

print("Polynomiality in fibers:")
for a in range(N):
    print(
        f"  eps_q{a + 1} rational in fibers: "
        f"{has_rational_fiber(E_q_jet[a], fiber_syms)}"
    )
print()


# --- VT homothety and integral ---
def homothety(expr, u):
    sub = {}
    for a in range(N):
        for k in range(r_max + 1):
            sub[Q_jet[a][k]] = u * Q_jet[a][k]
    return expr.subs(sub)


eps_q_u = [expand(homothety(E_q_jet[a], u)) for a in range(N)]
print("=" * 72)
print("Step 4: VT integrand u-power analysis")
print("=" * 72)
for a in range(N):
    poly_u = Poly(eps_q_u[a], u)
    monoms = poly_u.monoms()
    if monoms:
        min_u_pow = min(m[0] for m in monoms)
        max_u_pow = max(m[0] for m in monoms)
        print(f"  eps_q{a + 1}(u*jet): u-power range [{min_u_pow}, {max_u_pow}]")
    else:
        print(f"  eps_q{a + 1}(u*jet): zero")
print()


# Compute VT Lagrangian: L_VT = sum_a Q[a][0] * INT_0^1 eps_qa(u*jet) du
L_VT = 0
for a in range(N):
    integrand = expand(Q_jet[a][0] * eps_q_u[a])
    L_VT += integrate(integrand, (u, 0, 1))
L_VT = expand(L_VT)

print("=" * 72)
print("Step 5: VT Lagrangian (truncated print)")
print("=" * 72)
print(f"L_VT has {len(L_VT.args) if isinstance(L_VT, sp.Add) else 1} terms")
print()


# --- Verify EL(L_VT) = E^(1)_qa ---
def jet_to_funcs(expr):
    out = expr
    for a in range(N):
        out = out.subs(Q_jet[a][0], qs[a])
        for k in range(1, r_max + 1):
            out = out.subs(Q_jet[a][k], qs[a].diff(t, k))
    return out


L_VT_func = jet_to_funcs(L_VT)
EL_VT_jet = []
for a in range(N):
    el = euler_lagrange(L_VT_func, qs[a], t, max_order=2 * r_max + 2)
    el_jet = expand(to_jet(el))
    EL_VT_jet.append(el_jet)

print("=" * 72)
print("Step 6: consistency check  EL(L_VT) == eps_qa")
print("=" * 72)
all_ok = True
for a in range(N):
    diff_a = simplify(EL_VT_jet[a] - E_q_jet[a])
    ok = diff_a == 0
    print(f"  EL_q{a + 1}(L_VT) - eps_q{a + 1} = {diff_a}  -- {'OK' if ok else 'FAIL'}")
    all_ok = all_ok and ok
print()
print(f"Variational completion verified: {all_ok}")
print()


# --- Pathology checks ---
print("=" * 72)
print("Step 7: parameter-limit / denominator check")
print("=" * 72)
num, den = sp.fraction(sp.together(L_VT))
print(f"L_VT denominator: {den}")
# Strip the trivial /2 factors from polynomial integration
den_pure = sp.cancel(den / sp.gcd(den, 2))
print(f"L_VT denominator (mod /2): {den_pure}")
free = den_pure.free_symbols if hasattr(den_pure, "free_symbols") else set()
if not free:
    print("  -> L_VT is REGULAR (polynomial after /2 normalization)")
else:
    print(f"  -> L_VT has rational dependence on: {free}")
print()


# limits: b5 -> 0, M2_a -> 0, m_a -> 0
print("Limit b5 -> 0:")
L_b50 = expand(L_VT.subs(b5, 0))
print(f"  Number of terms: {len(L_b50.args) if isinstance(L_b50, sp.Add) else 1}")
print()

# Check denominator after b5 -> 0 limit (does it remain finite?)
num0, den0 = sp.fraction(sp.together(L_b50))
print(f"  Denominator at b5=0: {den0}")
print()

# Test M2 -> 0 (singular limit): M2_a in denominator should appear from
# the matrix-inverse h^(0).
print("Limit M2_1 -> 0:")
try:
    L_M0 = expand(L_VT.subs(M2[0], 0))
    print(
        f"  Result: {'finite' if L_M0.is_finite is None or L_M0.is_finite else 'INFINITE'}"
    )
except Exception as e:
    print(f"  EXCEPTION: {e}")

# Check denominator structure
num_full, den_full = sp.fraction(sp.together(L_VT))
free_M = [M2_ for M2_ in M2 if M2_ in den_full.free_symbols]
print(f"M2 appearing in L_VT denominator: {free_M}")
print()


# --- Final verdict ---
print("=" * 72)
print(f"FINAL VERDICT for T3 (N={N}+{N}):")
print("=" * 72)

divergent = False
reasons = []
for a in range(N):
    if has_rational_fiber(E_q_jet[a], fiber_syms):
        divergent = True
        reasons.append(f"eps_q{a + 1} rational in fibers")
    poly_u = Poly(eps_q_u[a], u)
    monoms = poly_u.monoms()
    if monoms:
        min_pow = min(m[0] for m in monoms)
        if min_pow < 0:
            divergent = True
            reasons.append(f"eps_q{a + 1}(u*jet) has u^{min_pow} singular term")

if divergent:
    print("  VT integral DIVERGES.  Reasons:")
    for r in reasons:
        print(f"    - {r}")
else:
    print("  VT integral CONVERGES on the homotopy [0,1].")
    print("  All eps_qa(u*jet) are polynomial in u.")
    print(f"  Variational completion EL match: {all_ok}")
    print("  L_VT denominators are polynomial in (M2, m, lambda, mu, K, b5).")
    print()
    print("  The expected M2_a^{-1} factors come from inverting the")
    print("  algebraic h-mass matrix (Routhian projector); these are NOT")
    print("  homotopy-integral pathologies.  M2 -> 0 is a separate physical")
    print("  limit (massless constraint field, not the b5 -> 0 critical")
    print("  surface).")
