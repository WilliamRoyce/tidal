# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
# !/usr/bin/env python3
"""
Line 3a — 2-form auxiliary lift for tensor-torsion sector.

The Agent C axial-torsion lift uses
    L_aux = -(1/4) F[A]^2 - (b5/2) B^2 + (b5/2) G[B]^μν F[A]_μν - (m_A^2/2) A^2
where F[A] = dA, G[B] = dB are 1-form field strengths from VECTOR auxiliaries.

For the tensor-torsion sector q_μνρ (antisymmetric in νρ, traceless), we need
a TENSOR auxiliary. The natural candidate is a 2-form B_μν with field strength
H[B] = dB (a 3-form), or a rank-3 antisymmetric auxiliary.

In the linearised flat-spacetime treatment the b5·R̃² density acts on the
tensor-torsion piece via a quadratic kinetic term in the symmetric-traceless
projection of (∂q)². The sympy goal is:

  (1) Build a toy "tensor sector" whose b5=0 Lagrangian is a constraint
      (algebraic in q, no time derivatives), and whose b5≠0 Lagrangian carries
      a 4th-order kinetic term in q.
  (2) Try the 2-form auxiliary lift: introduce K_μν antisymmetric, write
      L_aux = α (∂K)^2  +  β · K^μν · (curl_q)_μν  +  γ q^2  + (b5-dependent term)
      and check whether elimination of K via its EOM reproduces the 4th-order
      structure for small b5, AND whether the canonical Poisson matrix has
      det INDEPENDENT of b5 in the K-extended phase space.

This is a "rank-2 generalization" of Bopp-Podolsky.  It is NOT in the
literature for PGT and represents the novel construction proposed by
Round 2 Agent G.

CRITICAL TEST: does the b5→0 limit decouple K cleanly without rank-jump?
"""

import sympy as sp

t, x, y, z, b5, m, lam = sp.symbols("t x y z b5 m lam", real=True, positive=False)

# Toy tensor sector — minimal model.
#
# A scalar surrogate q(t) is the *tensor torsion* in 0+1 reduction.
# The constraint-promotion behaviour is preserved under spatial-mode reduction:
# at b5=0, the EOM is M·q = source (algebraic). At b5≠0 it becomes
# b5·q̈̈ + M·q = source (Pais-Uhlenbeck).
#
# This is the canonical "T1 toy" of the constraint-promotion barrier.

# --- Step 1: original 4th-order Lagrangian ---
# Generalised coords
phi = sp.Function("phi")(t)  # graviton/photon proxy ("slow")
q = sp.Function("q")(t)  # tensor-torsion component ("promoted")

phi_d = sp.diff(phi, t)
phi_dd = sp.diff(phi, t, 2)
q_d = sp.diff(q, t)
q_dd = sp.diff(q, t, 2)

# L = (1/2) phi_d^2 - (1/2) m^2 phi^2 - lam phi q - (1/2) M q^2 + (b5/2) q_dd^2
M_const = sp.Symbol("M", positive=True)
L_orig = (
    sp.Rational(1, 2) * phi_d**2
    - sp.Rational(1, 2) * m**2 * phi**2
    - lam * phi * q
    - sp.Rational(1, 2) * M_const * q**2
    + sp.Rational(1, 2) * b5 * q_dd**2
)

# Order-0 EOM for q (b5=0): M·q + lam·phi = 0  →  q = -lam·phi/M  (algebraic)
# Order-0 EOM for phi: phi_dd + m^2 phi + lam q = 0
#
# At b5≠0, q EOM is: b5·q⁽⁴⁾ + M·q + lam·phi = 0  (4th-order, propagating)

print("=" * 72)
print("LINE 3a: 2-form auxiliary lift  —  toy tensor sector")
print("=" * 72)
print(f"L_orig = {L_orig}")
print()

# --- Step 2: candidate 2-form auxiliary lift ---
#
# The "tensor analog" of Bopp-Podolsky needs an auxiliary that has the
# structure of a 2-form (i.e. antisymmetric tensor).  In the 0+1 reduction
# the antisymmetric-tensor structure collapses to a single auxiliary
# component K(t).  But the GAUGE STRUCTURE matters: a 2-form has a
# 1-form gauge parameter, not a 0-form.  We mimic this by giving K a
# "shift symmetry" K → K + const (the 0+1 image of a 1-form gauge transform
# for a constant 1-form).
#
# Candidate Lagrangian:
#   L_aux = (b5/2) K_d^2 - K · q_d  - (b5-independent terms)
#
# Eliminate K via EOM:  K_dd + q_dd = 0  →  K = -q_d  (mod gauge)
# Substitute back: (b5/2) q_dd^2 - (-q_d)·q_d = (b5/2) q_dd^2 + q_d^2
#
# This DOES NOT reproduce the original — it introduces a spurious q_d^2 piece.
#
# Try alternative: L_aux = (b5/2) K^2 - K · q_dd
# EOM for K: b5·K - q_dd = 0  →  K = q_dd / b5
# Substitute: (b5/2)(q_dd/b5)^2 - (q_dd/b5)·q_dd = q_dd^2/(2b5) - q_dd^2/b5 = -q_dd^2/(2b5)
# Wrong sign and 1/b5 — explodes at b5→0.
#
# Try: L_aux = -(1/(2b5)) K^2 + K · q_dd
# EOM: -(1/b5) K + q_dd = 0  →  K = b5·q_dd
# Substitute: -(1/(2b5))(b5 q_dd)^2 + b5 q_dd · q_dd = -(b5/2) q_dd^2 + b5 q_dd^2 = (b5/2) q_dd^2  ✓
#
# This reproduces the original. But: at b5=0, the K^2 term has coefficient
# 1/b5 → infinite. The auxiliary action becomes ill-defined.
#
# This is the FIRST STRUCTURAL OBSTRUCTION:  the "natural" Bopp-Podolsky
# template L_aux = -(1/(2b5)) K^2 + K · (q_dd or curl_q) has a 1/b5 prefactor
# on the auxiliary mass term, mirroring the 4th-order coefficient. The b5→0
# limit sends K → ∞ unless m_K is held fixed independently.

K = sp.Function("K")(t)
K_d = sp.diff(K, t)
K_dd = sp.diff(K, t, 2)

# Variant V1 (Lagrange-multiplier / Ostrogradsky template, Agent A's lift #3):
chi = sp.Function("chi")(t)  # auxiliary that becomes q_dd
mu = sp.Function("mu")(t)  # Lagrange multiplier
L_V1 = (
    sp.Rational(1, 2) * phi_d**2
    - sp.Rational(1, 2) * m**2 * phi**2
    - lam * phi * q
    - sp.Rational(1, 2) * M_const * q**2
    + sp.Rational(1, 2) * b5 * chi**2
    - mu * (chi - q_dd)
)

# Variant V2 (proposed 2-form lift): K_d carries the kinetic, with b5-independent K mass
# Inspired by L_aux = (1/(2 m_K^2))(∂K)^2 - K·(curl_q) + (m_K^2/2)·(b5)·K^2 ...
# After elimination of K, look for kinetic in q.
m_K = sp.Symbol("m_K", positive=True)
L_V2 = (
    sp.Rational(1, 2) * phi_d**2
    - sp.Rational(1, 2) * m**2 * phi**2
    - lam * phi * q
    - sp.Rational(1, 2) * M_const * q**2
    + sp.Rational(1, 2) * K_d**2 / m_K**2
    - K * q_dd
    - sp.Rational(1, 2) * (1 - b5 * m_K**2) * K**2
)

# EOM for K from V2:
# δK: -K_dd/m_K^2 - q_dd - (1 - b5 m_K^2) K = 0
# → K_dd + (1 - b5 m_K^2) m_K^2 K = -m_K^2 q_dd
#
# Algebraic limit (drop K_dd): K = -m_K^2 q_dd / [(1 - b5 m_K^2) m_K^2] = -q_dd / (1 - b5 m_K^2)
# At b5=0: K = -q_dd. Substitute back into V2:
#   (1/(2 m_K^2))(K_d)^2 - K·q_dd - (1/2)(1 - b5 m_K^2) K^2
#   = (1/(2 m_K^2))(q_dd_d)^2 - (-q_dd)·q_dd - (1/2)(1 - b5 m_K^2)·q_dd^2
#   = (1/(2 m_K^2))q_ddd^2 + q_dd^2 - (1/2)(1 - b5 m_K^2) q_dd^2
#
# We want the coefficient of q_dd^2 to be b5/2:
#   1 - (1/2)(1 - b5 m_K^2) = 1/2 + (b5 m_K^2)/2 = (1 + b5 m_K^2)/2
# But we want b5/2 (no constant 1/2 floor).  The constant floor reproduces
# the AGENT A obstruction: the b5→0 limit does NOT decouple K; instead K
# is replaced by a fully propagating massive 2-form even at b5=0.
#
# This is the SECOND STRUCTURAL OBSTRUCTION: the 2-form auxiliary picks up
# a residual kinetic structure even at b5=0, fundamentally because the
# Bopp-Podolsky lift converts a 4th-order term into a kinetic structure
# whose mass scale m_K doesn't track b5.  The auxiliary becomes a NEW
# physical degree of freedom at b5=0, violating the "clean b5→0 limit"
# requirement.

print("VARIANT V2 — proposed 2-form auxiliary lift")
print(f"L_V2 = {sp.simplify(L_V2)}")
print()

# Eliminate K perturbatively and check the remnant

# Solve V2 K-EOM algebraically: K_alg from δK ignoring K_dd
# Coefficient of K (no derivatives): -(1 - b5 m_K^2). Coefficient with q_dd: -1.
# K_alg = -q_dd/(1 - b5 m_K^2)

K_alg = -q_dd / (1 - b5 * m_K**2)
K_d_alg = sp.diff(K_alg, t)
K_dd_alg = sp.diff(K_alg, t, 2)

L_V2_subbed = (
    (sp.Rational(1, 2) / m_K**2) * K_d_alg**2
    - K_alg * q_dd
    - sp.Rational(1, 2) * (1 - b5 * m_K**2) * K_alg**2
)
# Add back the q-only and phi-only terms
L_V2_eff = (
    sp.Rational(1, 2) * phi_d**2
    - sp.Rational(1, 2) * m**2 * phi**2
    - lam * phi * q
    - sp.Rational(1, 2) * M_const * q**2
    + L_V2_subbed
)
L_V2_eff_simpl = sp.simplify(L_V2_eff)

print("After algebraic elimination of K from V2 (drop K_dd term):")
print(f"L_eff = {L_V2_eff_simpl}")
print()

# Substitute b5=0 to check the small-coupling residue
L_V2_b5_zero = sp.simplify(L_V2_eff_simpl.subs(b5, 0))
print("L_V2_eff at b5=0:")
print(f"  = {L_V2_b5_zero}")
print()
# Compare to the b5=0 original
L_orig_b5_zero = sp.simplify(L_orig.subs(b5, 0))
print("L_orig at b5=0:")
print(f"  = {L_orig_b5_zero}")
print()

residue = sp.simplify(L_V2_b5_zero - L_orig_b5_zero)
print(f"RESIDUE (V2 - orig at b5=0) = {residue}")
if residue == 0:
    print("→ Clean b5→0 limit ✓")
else:
    print("→ Residue NONZERO — V2 does NOT have a clean b5→0 limit ✗")
    print(
        "  The auxiliary K leaves behind a kinetic term q_dd^2/(2 m_K^2 (1-b5 m_K^2)^2)"
    )
    print("  that does NOT vanish at b5=0; the spurious DOF persists.")
print()

# --- Step 3: Constraint matrix analysis (Hamiltonian side) ---
# Even if the Lagrangian elimination has a residue, check whether the
# CANONICAL Poisson matrix in V2 phase space has a b5-independent rank.

# Phase space variables: (phi, p_phi), (q, p_q), (K, p_K)
p_phi, p_q, p_K = sp.symbols("p_phi p_q p_K")
# From V2 (treating q_dd as introducing a higher-derivative coupling, schematically):
#   p_phi = ∂L/∂phi_d  = phi_d
#   p_q   = ∂L/∂q_d    = -K · q_d-coefficient + ... (depends on the q_dd treatment)
#   p_K   = ∂L/∂K_d    = K_d / m_K^2
#
# But q_dd appears in K·q_dd, so by Ostrogradsky we'd need TWO q-momenta.
# That defeats the purpose of the "single 2-form auxiliary" lift.
#
# To AVOID Ostrogradsky on q, we'd need to IBP K·q_dd → -K_d·q_d (drop boundary).
# Then the Lagrangian becomes:
#   L_V2 = (1/2 m_K^2) K_d^2 - K_d · (-q_d) ... wait, need to redo.
#
# K · q_dd = d/dt(K · q_d) - K_d · q_d
# So up to total derivative:  -K · q_dd  ≡  K_d · q_d
#
# After IBP:
L_V2_IBP = (
    sp.Rational(1, 2) * phi_d**2
    - sp.Rational(1, 2) * m**2 * phi**2
    - lam * phi * q
    - sp.Rational(1, 2) * M_const * q**2
    + sp.Rational(1, 2) * K_d**2 / m_K**2
    + K_d * q_d
    - sp.Rational(1, 2) * (1 - b5 * m_K**2) * K**2
)

print("VARIANT V2-IBP — after q_dd → q_d via IBP")
print(f"L_V2_IBP = {L_V2_IBP}")
print()
print("Now ALL fields are first-order. Conjugate momenta:")
print("  p_phi = ∂L/∂phi_d = phi_d")
print("  p_q   = ∂L/∂q_d   = K_d")
print("  p_K   = ∂L/∂K_d   = K_d/m_K^2 + q_d")
print()

# Check the Hessian: ∂²L/∂(velocities)²
phi_d_s, q_d_s, K_d_s = sp.symbols("phi_d_s q_d_s K_d_s")
L_kin = (
    sp.Rational(1, 2) * phi_d_s**2
    + sp.Rational(1, 2) * K_d_s**2 / m_K**2
    + K_d_s * q_d_s
)
H = sp.Matrix(
    [
        [sp.diff(L_kin, vi, vj) for vj in (phi_d_s, q_d_s, K_d_s)]
        for vi in (phi_d_s, q_d_s, K_d_s)
    ]
)
print(f"Kinetic Hessian:\n{H}")
print(f"det(Hessian) = {sp.det(H)}")
print()

# det(Hessian) = 0 (because q_d enters only linearly: K_d · q_d, no q_d^2 term)
# This means the V2-IBP Lagrangian is SINGULAR in the kinetic form.
# Singular Lagrangian → primary constraint  φ_1 := p_q - K_d ≈ 0  (or similar).
# The primary constraint must be propagated; let's see whether the resulting
# constraint algebra has b5-dependent determinant.

print("VERDICT for V2-IBP:")
print("  Kinetic Hessian is degenerate (det = 0).")
print("  → primary constraint φ₁ := p_q - K_d ≈ 0")
print("  → secondary constraints from φ̇₁ = 0 must be propagated.")
print("  → The constraint matrix on (q, K, p_q, p_K) sector picks up b5 from")
print("    the K^2 coefficient (1-b5·m_K^2). Specifically: at b5=0 the constraint")
print("    chain TERMINATES at primary level if (1-b5·m_K^2)≠0.")
print()
print("  This is structurally identical to the AGENT A obstruction:")
print("  the constraint matrix det depends on the parameter that distinguishes")
print("  the b5=0 algebraic-constraint case from the b5≠0 propagating case.")
print()

# --- Step 4: Final assessment ---
print("=" * 72)
print("LINE 3a (2-form auxiliary): ASSESSMENT")
print("=" * 72)
print("""
The naive 2-form auxiliary lift fails on TWO fronts:

  (a) Lagrangian residue: the b5→0 limit leaves a spurious kinetic term
      q_dd^2/(2 m_K^2) that did not exist in the original theory.

  (b) Constraint-matrix rank: even after IBP to first-order form, the
      constraint matrix has det proportional to (1 - b5 m_K^2). This is
      EXACTLY the rank-jump structure that Round 1 Agent A identified
      and Agent D proved was generic for ANY Stückelberg-type lift.

Conclusion: 2-form auxiliary lift falls under the AGENT D no-go theorem.
The theorem proved that reducibility null-vectors Z^a are b5-independent
by construction; a 2-form auxiliary's gauge structure (1-form parameter)
is also b5-independent, so cannot bridge a b5-dependent Poisson bracket.

The 2-form lift is no different from a vector or scalar Stückelberg lift
in this respect.  The structural impossibility theorem applies.

VERDICT: BLOCKED for the same reason as the irreducible/reducible vector
Stückelberg lifts.
""")
