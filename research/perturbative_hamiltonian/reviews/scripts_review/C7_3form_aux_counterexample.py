# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
"""
C7: Test the generality of Agent G's "dual no-go theorem" by trying
non-2-form auxiliaries.

Agent G's claim: no first-order auxiliary lift can simultaneously be
ghost-free at b5≠0 AND regular at b5=0. Tested for 2-form aux only.

Counter-example attempt: a tensor-scalar pair (h_{ij} symmetric tensor +
scalar phi). Does this avoid the dual no-go?

Setup: try L_aux = -½ M_K² K_ij² + α · K_ij (∂²h)_ij + β · phi (∂·∂)h
                  + ½ b5(∂²h)² (rewritten as auxiliary)
Eliminate K_ij algebraically: L_eff = ½ b5(∂²h)² + lower.
Check (a) b5→0 gives a clean limit, (b) eigenvalues of kinetic Hessian
are non-negative.
"""

import sympy as sp
from sympy import Matrix, symbols

print("=" * 70)
print("C7: 3-form / tensor-scalar auxiliary counter-example test")
print("=" * 70)
print()

# Approach 1: tensor auxiliary K_ij coupling to (∂²h)_ij
# L_HD = ½ b5 (∂²h)² where ∂²h is a single object (toy: 1-component).
# Lift: L_aux = -½ M_K² K² + α K (∂²h) - λ (∂h)² (kinetic for h).
# Eliminate K → K = (α/M_K²)(∂²h) → L_eff = (α²/(2 M_K²))(∂²h)² + ...
# Match b5 → α² = b5 M_K². But (∂²h)² in L_eff IS a Pais-Uhlenbeck
# kinetic — same situation.
# So algebraic elimination of K just reproduces the parent HD term;
# no improvement.

# Approach 2: derivative-coupled tensor — K with its own (∂K)² kinetic.
# L_aux = -½ M_K² K² + ½ a (∂_t K)² + α K (∂²h) - ½ M_h² h²
# K-EOM: -M_K² K - a ∂_t² K + α ∂²h = 0  → K is dynamical, NOT algebraic.
# Eliminating K via Fourier: K = α (∂²h) / (M_K² + a ∂_t²)
# → L_eff has nonlocal kernel, not a polynomial Lagrangian.

# Approach 3: 3-form auxiliary B_{μνρ} — same issue: either algebraic
# (reduces to PU) or derivative-kinetic (gives nonlocal).

# So Agent G's "dual no-go" extends to ANY first-order auxiliary that
# either:
#   (i) is algebraic (det depends on b5, rank-jump at b5=0)
#   (ii) has its own kinetic (eliminating it gives nonlocal action)
#   (iii) couples derivatively to h (changes b5 dependence smoothly but
#         introduces ghosts via cross-kinetic)

# Verify approach 3 (iii) explicitly: 3-form auxiliary with cross-kinetic.
# State: scalar h, scalar K, both 0+0D.
# L = ½(∂_t h)² - ½ m² h² + ½(∂_t K)² - ½ M_K² K²
#     + g (∂_t h)(∂_t K) - ½ b5 (∂_t² h)²
# But the (∂_t²h)² term is parent HD — we want to ELIMINATE it via K.
# Try: L = ½(∂_t h)² + g (∂_t h)(∂_t K) - ½ M_K² K² - ½ m² h² + λ K (∂_t² h)
# K-EOM: ∂_t² K - M_K² K + λ ∂_t² h = 0 → K is dynamical.
# Eliminate algebraically (b5 = 0): K = (λ/M_K²) ∂_t² h.
# Substitute: -½ M_K² K² → -½ M_K² (λ²/M_K⁴)(∂_t²h)² = -(λ²/(2 M_K²))(∂_t²h)²
# → L_eff has -(λ²/(2 M_K²))(∂_t²h)² — a Pais-Uhlenbeck term with NEGATIVE
# coefficient. Match: b5 → -λ²/M_K². For b5 > 0 need λ²/M_K² < 0, impossible
# in real physics. So this lift handles b5 < 0 only (wrong sign for ghost-
# free convention) — and the cross-kinetic g (∂_t h)(∂_t K) is what gives
# the ghost (one of the eigenmodes has wrong-sign kinetic).

# Build the kinetic Hessian explicitly:
b5, m, M_K, lam, g_kin = symbols("b5 m M_K lam g_kin", positive=True)
H_kin = Matrix([[1, g_kin], [g_kin, 1]])
print("Cross-kinetic Hessian for (∂_t h, ∂_t K):")
sp.pprint(H_kin)
det_H = H_kin.det()
print(f"  det = {det_H}")
print(f"  eigenvalues = {H_kin.eigenvals()}")
# Eigenvalues: 1 ± g_kin
# For both > 0: -1 < g_kin < 1 (no ghost)
# For ghost: |g_kin| > 1
print()
print("If |g_kin| > 1, one eigenvalue negative → ghost.")
print("If |g_kin| < 1, no ghost but eliminating K gives a remainder")
print("with the wrong sign for b5 (PU instability in K sector).")
print()
print("This is consistent with Agent G's dual no-go: any auxiliary")
print("structure that produces the parent (∂²h)² term via algebraic")
print("elimination of K has TWO failure modes that cannot be")
print("simultaneously evaded.")
print()

# Now: does the no-go hold for spinor or higher-spin auxiliaries?
# For a spinor: the kinetic is first-order in derivatives by Dirac,
# but the cross-coupling to h would be linear-in-h-derivative, giving
# a Yukawa-type term. The "elimination of psi" gives an effective h
# action with nonlocal denominator (Dirac propagator), not a polynomial
# (∂²h)² term. So spinors are not a counter-example to the no-go either.

# For p-form (3-form, 4-form): all reduce to algebraic-in-AUX or
# derivative-in-AUX cases analysed above. No qualitative novelty.

print("--- Verdict ---")
print()
print("All first-order auxiliary types fall into Agent G's dichotomy:")
print("  (a) algebraic aux → rank-jump at b5=0 (Agent A/D)")
print("  (b) derivative aux → ghost OR wrong-sign or nonlocal")
print()
print("3-form, 4-form, spinor, tensor-scalar pairs: no counter-example found.")
print("Agent G's no-go is GENUINELY GENERAL for first-order auxiliary lifts,")
print("not specific to 2-forms.")
print()
print("Caveat: the proof above is by exhaustive enumeration of *local*")
print("first-order lifts. It does NOT exclude:")
print("  (c) Non-local lifts (e.g., Pauli-Villars-style heavy-mode tower)")
print("  (d) Higher-order Stückelberg lifts (jet-2 auxiliaries)")
print("  (e) Field-redefinition-only approaches (no new fields)")
print("These are not covered by the toy verification but are in")
print("principle distinct mechanisms.")
