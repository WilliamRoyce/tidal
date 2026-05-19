# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
"""
Vainberg-Tonti integral convergence test for T2 toy.
====================================================

Toy model (Agent B Round 1):

    L_T2 = (1/2)(d_t phi)^2 + (1/2)(d_t chi)^2
         - (1/2) m_phi^2 phi^2 - (1/2) m_chi^2 chi^2
         - lambda * phi * h - lambda_chi * chi * d_t h
         - (1/2) h^2 + (1/2) b5 (d_t^2 h)^2

The constraint field h is purely auxiliary at b5=0 (algebraic EOM).
Parker-Simon reduction substitutes h -> h_sol = h^(0) + b5*h^(1) + O(b5^2),
yielding reduced EOMs E^(1)_phi, E^(1)_chi on the (phi, chi) jet bundle.

Agent B confirmed Helmholtz residues vanish, i.e., the reduced source form
is variational.  Per Krupka-Voicu (arXiv:1406.6646, eq. 11) we can then
construct the Vainberg-Tonti Lagrangian as

    L_VT(jet) = sum_sigma  y^sigma * INT_0^1  eps_sigma( u * jet )  du

where the homothety chi_u scales all jet variables uniformly by u.

Voicu 2020 (arXiv:2009.05459) shows that for 4D-Gauss-Bonnet this integral
DIVERGES because the integrand behaves like u^{-1} near u=0 (the
truncated A_munu transforms with negative degree of homogeneity in g).
The pole appears as a (D-4)^{-1} prefactor when regulated.

This script:
  1. Reproduces Agent B's PS-reduced EOMs in jet form.
  2. Computes the VT integrand and the u-integral symbolically.
  3. Verifies the Euler-Lagrange of L_VT reproduces E^(1)_phi, E^(1)_chi.
  4. Tests parameter-limit pathologies: b5 -> 0, m_phi -> 0, m_chi -> 0,
     lambda -> 0, lambda_chi -> 0.

Convergence criterion: VT integrand must be polynomial (or at most
polynomial * non-singular) in u across the entire integration domain
[0,1].  Any 1/u, log(u), or 1/(u-c) factor with c in [0,1] -> divergent.
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
b5, mphi, mchi, lam, lamc = symbols("b5 m_phi m_chi lambda lambda_chi", real=True)

phi = Function("phi")(t)
chi = Function("chi")(t)
h = Function("h")(t)


# ---------------------------------------------------------------------------
# Step 1: build T2 Lagrangian and reproduce PS reduction (Agent B)
# ---------------------------------------------------------------------------

L_T2 = (
    Rational(1, 2) * phi.diff(t) ** 2
    + Rational(1, 2) * chi.diff(t) ** 2
    - Rational(1, 2) * mphi**2 * phi**2
    - Rational(1, 2) * mchi**2 * chi**2
    - lam * phi * h
    - lamc * chi * h.diff(t)
    - Rational(1, 2) * h**2
    + Rational(1, 2) * b5 * h.diff(t, 2) ** 2
)


def euler_lagrange(L, q, t, max_order=8):
    expr = diff(L, q)
    for k in range(1, max_order + 1):
        qk = q.diff(t, k)
        term = diff(L, qk)
        if term == 0:
            continue
        expr += (-1) ** k * term.diff(t, k)
    return expr


# Order-0 algebraic constraint for h:
#   E_h = -lambda*phi + lambda_chi*d_t(chi) - h + b5*d_t^4(h) = 0
# at b5 = 0 -> h^(0) = -lambda*phi + lambda_chi*d_t(chi)
h0 = -lam * phi + lamc * chi.diff(t)

# Order-1 correction comes from substituting h^(0) into the b5*d_t^4 h piece:
#   h^(1) = d_t^4 h^(0)
h_sol = h0 + b5 * h0.diff(t, 4)

# Substitute h(t) and its derivatives into E_phi, E_chi and truncate at O(b5)
E_phi_full = euler_lagrange(L_T2, phi, t)
E_chi_full = euler_lagrange(L_T2, chi, t)


def sub_h_in_eom(expr, h_sol_expr, max_h_deriv=8):
    out = expr
    for k in range(max_h_deriv, 0, -1):
        out = out.subs(h.diff(t, k), h_sol_expr.diff(t, k))
    return out.subs(h, h_sol_expr)


E_phi_red = expand(sub_h_in_eom(E_phi_full, h_sol).series(b5, 0, 2).removeO())
E_chi_red = expand(sub_h_in_eom(E_chi_full, h_sol).series(b5, 0, 2).removeO())

print("=" * 72)
print("Step 1: Parker-Simon-reduced EOMs (truncated at O(b5))")
print("=" * 72)
print("E^(1)_phi =", E_phi_red)
print()
print("E^(1)_chi =", E_chi_red)
print()


# ---------------------------------------------------------------------------
# Step 2: convert EOMs to jet form (replace d_t^k phi -> Phi_k symbol)
# ---------------------------------------------------------------------------

# Maximum derivative order appearing in PS-reduced EOMs.
# We need to scan E_phi_red / E_chi_red.


def max_deriv_order(expr, q, t):
    n = 0
    for d in expr.atoms(sp.Derivative):
        if d.expr == q:
            order = sum(c for v, c in d.variable_count if v == t)
            n = max(n, order)
    return n


r_phi = max_deriv_order(E_phi_red, phi, t)
r_chi_in_phi = max_deriv_order(E_phi_red, chi, t)
r_phi_in_chi = max_deriv_order(E_chi_red, phi, t)
r_chi = max_deriv_order(E_chi_red, chi, t)
r_max = max(r_phi, r_chi_in_phi, r_phi_in_chi, r_chi)
print(f"Max derivative order in PS-reduced EOMs: r = {r_max}")
print(f"  phi appears in E_phi_red up to order {r_phi}")
print(f"  chi appears in E_phi_red up to order {r_chi_in_phi}")
print(f"  phi appears in E_chi_red up to order {r_phi_in_chi}")
print(f"  chi appears in E_chi_red up to order {r_chi}")
print()


# Jet variables: Phi_k = d_t^k phi, Chi_k = d_t^k chi for k = 0..r_max
Phi = [symbols(f"Phi_{k}", real=True) for k in range(r_max + 1)]
Chi = [symbols(f"Chi_{k}", real=True) for k in range(r_max + 1)]


def to_jet(expr):
    """Replace d_t^k phi(t) -> Phi[k], d_t^k chi(t) -> Chi[k]."""
    out = expr
    for k in range(r_max, 0, -1):
        out = out.subs(phi.diff(t, k), Phi[k])
        out = out.subs(chi.diff(t, k), Chi[k])
    out = out.subs(phi, Phi[0])
    return out.subs(chi, Chi[0])


E_phi_jet = expand(to_jet(E_phi_red))
E_chi_jet = expand(to_jet(E_chi_red))
print("=" * 72)
print("Step 2: jet-form EOMs")
print("=" * 72)
print("eps_phi(jet) =", E_phi_jet)
print()
print("eps_chi(jet) =", E_chi_jet)
print()


# ---------------------------------------------------------------------------
# Step 3: build VT integrand and integrate over u in [0, 1]
# ---------------------------------------------------------------------------
# VT formula (Krupka-Voicu eq. 11):
#   L_VT(jet) = sum_sigma y^sigma * INT_0^1 eps_sigma(u * jet) du
#             = Phi[0] * INT_0^1 eps_phi(u * jet) du
#             + Chi[0] * INT_0^1 eps_chi(u * jet) du
#
# The homothety chi_u: every Phi[k] -> u*Phi[k], every Chi[k] -> u*Chi[k].
# Note: u also appears NOT scaling the parameters (b5, m, lambda, lambda_chi).


def homothety(expr, u):
    sub = {}
    for k in range(r_max + 1):
        sub[Phi[k]] = u * Phi[k]
        sub[Chi[k]] = u * Chi[k]
    return expr.subs(sub)


eps_phi_u = expand(homothety(E_phi_jet, u))
eps_chi_u = expand(homothety(E_chi_jet, u))

print("=" * 72)
print("Step 3: VT integrand structure")
print("=" * 72)
print("eps_phi(u*jet) =", eps_phi_u)
print()
print("eps_chi(u*jet) =", eps_chi_u)
print()

# Check that eps_phi(u*jet) is polynomial in u (sufficient for convergence
# of the homotopy integral over [0,1]).
poly_phi = Poly(eps_phi_u, u)
poly_chi = Poly(eps_chi_u, u)
print("Polynomial degree in u:")
print(f"  eps_phi(u*jet): degree {poly_phi.degree()}")
print(f"  eps_chi(u*jet): degree {poly_chi.degree()}")

# Check leading u-power: if any coefficient at u^0 is non-zero AND it
# multiplies Phi[0] or Chi[0] (degree-0 in fiber), no divergence
# from the t->0 endpoint.  The general convergence test is whether
# Phi[0]*eps_phi(u*jet) is bounded as u -> 0.
integrand_phi = expand(Phi[0] * eps_phi_u)
integrand_chi = expand(Chi[0] * eps_chi_u)
poly_int_phi = Poly(integrand_phi, u)
poly_int_chi = Poly(integrand_chi, u)
print()
print("VT integrand polynomial degree in u (must be >= 0 polynomial):")
print(
    f"  Phi[0]*eps_phi(u*jet): degree {poly_int_phi.degree()},"
    f" lowest power {poly_int_phi.monoms()[-1] if poly_int_phi.monoms() else 'NA'}"
)
print(
    f"  Chi[0]*eps_chi(u*jet): degree {poly_int_chi.degree()},"
    f" lowest power {poly_int_chi.monoms()[-1] if poly_int_chi.monoms() else 'NA'}"
)
print()

# Integrate over u in [0,1]
I_phi = integrate(integrand_phi, (u, 0, 1))
I_chi = integrate(integrand_chi, (u, 0, 1))
L_VT = expand(I_phi + I_chi)

print("=" * 72)
print("Step 4: VT Lagrangian (integrated over u in [0,1])")
print("=" * 72)
print("INT[Phi*eps_phi du] =", expand(I_phi))
print()
print("INT[Chi*eps_chi du] =", expand(I_chi))
print()
print("L_VT(jet) =", L_VT)
print()


# ---------------------------------------------------------------------------
# Step 5: verify L_VT reproduces the original EOMs via Euler-Lagrange
# ---------------------------------------------------------------------------
# This is the canonical-completion consistency check (Krupka-Voicu Thm 1
# / Voicu 2020 eq. 17): if the source form IS variational, EL of L_VT
# returns it exactly.


def jet_to_funcs(expr):
    """Replace Phi[k] -> d_t^k phi(t), Chi[k] -> d_t^k chi(t)."""
    out = expr
    out = out.subs(Phi[0], phi)
    out = out.subs(Chi[0], chi)
    for k in range(1, r_max + 1):
        out = out.subs(Phi[k], phi.diff(t, k))
        out = out.subs(Chi[k], chi.diff(t, k))
    return out


L_VT_func = jet_to_funcs(L_VT)
EL_phi_VT = expand(euler_lagrange(L_VT_func, phi, t, max_order=2 * r_max + 2))
EL_chi_VT = expand(euler_lagrange(L_VT_func, chi, t, max_order=2 * r_max + 2))

# Reduce to jet form for comparison
EL_phi_VT_jet = expand(to_jet(EL_phi_VT))
EL_chi_VT_jet = expand(to_jet(EL_chi_VT))

diff_phi = expand(EL_phi_VT_jet - E_phi_jet)
diff_chi = expand(EL_chi_VT_jet - E_chi_jet)

print("=" * 72)
print("Step 5: consistency check  EL(L_VT) == E^(1)_sigma ?")
print("=" * 72)
print("EL_phi(L_VT) =", EL_phi_VT_jet)
print()
print("EL_chi(L_VT) =", EL_chi_VT_jet)
print()
print("EL_phi(L_VT) - eps_phi =", simplify(diff_phi))
print("EL_chi(L_VT) - eps_chi =", simplify(diff_chi))
print()

verified = (simplify(diff_phi) == 0) and (simplify(diff_chi) == 0)
print(f"Variational completion verified: {verified}")
print()


# ---------------------------------------------------------------------------
# Step 6: parameter-limit pathology check
# ---------------------------------------------------------------------------
# Voicu 2020's GB-divergence appeared as a (D-4)^{-1} pole.  Here we don't
# have a dimension parameter, but we do have b5, m, lambda, lambda_chi.
# Test L_VT in each limit:

print("=" * 72)
print("Step 6: parameter-limit pathology check")
print("=" * 72)

limits = [
    ("b5 -> 0", b5, 0),
    ("m_phi -> 0", mphi, 0),
    ("m_chi -> 0", mchi, 0),
    ("lambda -> 0", lam, 0),
    ("lambda_chi -> 0", lamc, 0),
]

for name, sym, val in limits:
    L_lim = expand(L_VT.subs(sym, val))
    print(f"  {name}: L_VT = {L_lim}")
print()

# Also: extreme-value rationality check.  Are there ANY 1/(parameters)
# factors in L_VT?  Examine denominators.
num, den = sp.fraction(sp.together(L_VT))
print(f"  L_VT denominator: {den}")
print(
    f"  -> {'REGULAR (polynomial in params)' if den == 1 else 'has rational structure'}"
)
print()


# ---------------------------------------------------------------------------
# Step 7: explicit comparison with the Agent-B-direct effective Lagrangian
# ---------------------------------------------------------------------------
# Agent B also computed L_eff = L_T2 |_{h = h_sol(b5)} and showed its EL
# equations equal the PS-reduced EOMs.  L_eff and L_VT should be EQUAL up
# to a total time derivative.  Compute their difference.

L_eff = (
    L_T2.subs(h.diff(t, 2), h_sol.diff(t, 2))
    .subs(h.diff(t), h_sol.diff(t))
    .subs(h, h_sol)
)
L_eff = expand(L_eff.series(b5, 0, 2).removeO())
L_eff_jet = expand(to_jet(L_eff))

# A total-time-derivative term in jet variables has the form
# d_t F = sum_k Phi[k+1] dF/dPhi[k] + sum_k Chi[k+1] dF/dChi[k].
# The difference (L_eff_jet - L_VT) being a total derivative is hard
# to verify directly in jet form, but they should at least agree on EOMs.
EL_phi_eff_jet = expand(to_jet(euler_lagrange(L_eff, phi, t, max_order=8)))
EL_chi_eff_jet = expand(to_jet(euler_lagrange(L_eff, chi, t, max_order=8)))

print("=" * 72)
print("Step 7: comparison L_eff (direct on-shell Routhian) vs L_VT")
print("=" * 72)
print("L_eff(jet) =", L_eff_jet)
print()
print("L_VT (jet) =", L_VT)
print()
print("L_eff - L_VT =", expand(L_eff_jet - L_VT))
print()
print("EL_phi(L_eff) - EL_phi(L_VT) =", simplify(EL_phi_eff_jet - EL_phi_VT_jet))
print("EL_chi(L_eff) - EL_chi(L_VT) =", simplify(EL_chi_eff_jet - EL_chi_VT_jet))
print()


# ---------------------------------------------------------------------------
# Step 8: explicit homogeneity-of-source-form analysis (Voicu 2020 key)
# ---------------------------------------------------------------------------
# Voicu 2020 identifies divergence by NEGATIVE homogeneity of E_munu under
# uniform fiber rescaling.  Under chi_u: jet -> u*jet, the integrand
# u * eps_sigma(u*jet) must be integrable on (0,1].  If eps_sigma is a
# polynomial of degree d_min >= 0 in fibers, the integrand is u^(1+d_min)
# which is fine.  If eps_sigma contains 1/Phi_k or 1/Chi_k structures,
# eps_sigma(u*jet) ~ u^(-1) and the integrand is u^0 (bounded) but a
# higher-order singularity gives a pole.

# For T2, eps_phi and eps_chi are POLYNOMIAL in jet variables (we built
# them by E-L of a polynomial Lagrangian).  Confirm:
ej_phi_atoms = E_phi_jet.atoms(sp.Symbol)
ej_chi_atoms = E_chi_jet.atoms(sp.Symbol)


# Check that no 1/Phi or 1/Chi appears (i.e. no rational dependence on
# fiber variables).
def has_rational_fiber(expr) -> bool:
    """Returns True if expr has 1/(Phi_k) or 1/(Chi_k) factors."""
    _n, d = sp.fraction(sp.together(expr))
    fiber_syms = set(Phi) | set(Chi)
    return any(f in fiber_syms for f in d.free_symbols)


print("=" * 72)
print("Step 8: source-form homogeneity analysis")
print("=" * 72)
print(f"eps_phi rational in fibers: {has_rational_fiber(E_phi_jet)}")
print(f"eps_chi rational in fibers: {has_rational_fiber(E_chi_jet)}")
print()
# Compute degree of homogeneity of eps_sigma(u*jet) in u (i.e., the
# minimal power of u appearing).
poly_eps_phi_u = Poly(eps_phi_u, u)
poly_eps_chi_u = Poly(eps_chi_u, u)
phi_min_u = (
    min(m[0] for m in poly_eps_phi_u.monoms()) if poly_eps_phi_u.monoms() else None
)
chi_min_u = (
    min(m[0] for m in poly_eps_chi_u.monoms()) if poly_eps_chi_u.monoms() else None
)
print(f"Min power of u in eps_phi(u*jet): {phi_min_u}")
print(f"Min power of u in eps_chi(u*jet): {chi_min_u}")
print()
print(
    "VT integrand u * Phi[0] * eps_phi(u*jet) has min u-power"
    f" = {phi_min_u}  (must be >= 0 for convergence)"
)
# Note: VT formula in Krupka-Voicu eq. 11 is L_VT = y * INT eps(u*y) du,
# WITHOUT a leading u factor.  The leading u in eq. 12 is for the
# correction-term tau.  For the Lagrangian itself, integrand is just
# y * eps(u*y).
print()


# ---------------------------------------------------------------------------
# Final verdict
# ---------------------------------------------------------------------------
print("=" * 72)
print("FINAL VERDICT for T2:")
print("=" * 72)
divergent = False
reasons = []
if has_rational_fiber(E_phi_jet) or has_rational_fiber(E_chi_jet):
    divergent = True
    reasons.append("source form is rational in fibers")
if phi_min_u is not None and phi_min_u < 0:
    divergent = True
    reasons.append(f"eps_phi(u*jet) has u^{phi_min_u} singular term")
if chi_min_u is not None and chi_min_u < 0:
    divergent = True
    reasons.append(f"eps_chi(u*jet) has u^{chi_min_u} singular term")

if divergent:
    print("  VT integral DIVERGES.  Reasons:")
    for r in reasons:
        print(f"    - {r}")
else:
    print("  VT integral CONVERGES.")
    print("  Source form is polynomial in fibers => integrand is")
    print("  polynomial in u => integral over [0,1] is finite.")
    print(f"  Reproduced original EOMs via EL: {verified}")
    if den == 1:
        print("  L_VT is regular (polynomial) in all parameters")
        print("  (b5, m_phi, m_chi, lambda, lambda_chi).")
    print()
    print("  L_VT =", L_VT)
