# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
"""
C8: Diff Agent F's axial constraint Poisson matrix vs Agent J's
Curtright-q kinetic Hessian to quantify the actual structural similarity.

Round 3 synthesis claim: "three independent paths converge". Are they
really independent and structurally parallel?

Agent F's M_axial:
  M_A_block       = ((0, -1), (1, 0)) * m_A^2 -> det = m_A^4 (b-independent)
  M_aux_block     = ((0, -1), (1, 0)) * (1/b) -> det = 1/b^2 (DIVERGES at b=0)

Agent J's H_kin (Curtright):
  4x4 kinetic Hessian after canonical rescaling, det = 1 - lam_a^2
  (λ_a = mass-rescaling parameter, NOT b5)
  -> b5-independent (under canonical rescaling).

DIFF:
  - Agent F's "uniformity" applies ONLY to the A-sector block; the
    auxiliary block has a singular b -> 0 limit (det -> infinity).
  - Agent J's "uniformity" applies to the full kinetic Hessian after
    rescaling.
  - Agent F's mechanism: aux phi has mass = 1/b; b -> 0 sends phi to
    infinite mass, decoupling.  This is the SAME singular limit as the
    parent theory (Lyakhovich rank-jump on aux block).
  - Agent J's mechanism: q has mass = M/b5; b5 -> 0 sends q to infinite
    mass, decoupling.  Auxiliaries h, b, a have finite kinetic.
"""

import sympy as sp
from sympy import Matrix, symbols

print("=" * 70)
print("C8: Structural diff — Agent F (axial) vs Agent J (Curtright-q)")
print("=" * 70)
print()

# Reproduce both matrices for direct comparison
b, m_A = symbols("b m_A", positive=True)
lam_a = symbols("lam_a", real=True)
b5 = symbols("b5", positive=True)

# Agent F: full 4x4 constraint Poisson matrix
M_F_full = Matrix(
    [
        [0, 0, -(m_A**2), 0],
        [0, 0, 0, -1 / b],
        [m_A**2, 0, 0, 0],
        [0, 1 / b, 0, 0],
    ]
)
print("Agent F's full constraint Poisson matrix (4x4):")
sp.pprint(M_F_full)
print(f"  det = {M_F_full.det()}")
print(f"  det(b * M) = {(b * M_F_full).det()}  -> vanishes at b=0  -> RANK-JUMP")
print()

# Agent J's 4x4 kinetic Hessian (canonical-rescaled)
H_J = Matrix(
    [
        [1, lam_a, 0, 0],
        [lam_a, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
)
print("Agent J's kinetic Hessian (4x4, canonical-rescaled):")
sp.pprint(H_J)
print(f"  det = {H_J.det()}")
print("  -> b5-independent")
print()

print("--- STRUCTURAL DIFF ---")
print()
print("Property                | Agent F (axial)        | Agent J (Curtright-q)")
print("-" * 75)
print("Object                  | Constraint Poisson M    | Kinetic Hessian H")
print("Dimension               | 4x4                     | 4x4 (1+1D toy)")
print("b/b5 dependence in det  | det ~ 1/b^2             | det = 1 - lam_a^2")
print("Behavior at b->0        | DIVERGES                | UNCHANGED")
print("Rank-jump?              | YES (det vanishes when  | NO")
print("                        |  multiplied by b^2)     |")
print("Decoupling mechanism    | aux mass -> infty as    | q mass -> infty as")
print("                        | b -> 0                  | b5 -> 0")
print()

print("--- STRUCTURAL PARALLELISM CLAIM ---")
print()
print("Round 3 synthesis claims: 'Both decouple the heavy field via")
print("infinite-mass; tensor case needs 3 auxiliaries, axial 1.'")
print()
print("This is TRUE at the qualitative level: both lifts have a heavy")
print("field that decouples in the b->0 limit.")
print()
print("But the QUANTITATIVE structures are DIFFERENT:")
print()
print("  Agent F: Constraint Poisson DIVERGES at b=0. The rank-uniformity")
print("    claim is restricted to the A-sector block only. The auxiliary")
print("    block has the same Lyakhovich rank-jump as the parent theory")
print("    -- it has been MOVED, not REMOVED.")
print()
print("  Agent J: Kinetic Hessian is FINITE and b5-independent at b5=0.")
print("    The mass-term divergence is the only b5-singular content,")
print("    and it confines q to zero amplitude (decoupling).")
print()
print("So Agent J's lift is structurally STRONGER than Agent F's: Agent F")
print("only re-shuffles the rank-jump; Agent J genuinely uniforms the")
print("kinetic structure.")
print()
print("The 'cross-validation' framing in Round 3 synthesis equates these")
print("two as if they were the same mechanism. They are NOT. Agent F's")
print("is a Lyakhovich-style lift with Poisson rank-jump pushed into")
print("the auxiliary; Agent J's is a Stueckelberg-style lift with")
print("kinetic uniformity. The latter is what would be needed for the")
print("metric h_4/h_7/h_9 sector, and Agent J flags it does NOT apply")
print("there.")
print()
print("This is consistent with Review 2's finding that 'cross-validation'")
print("inflation is overstated. The two lifts are qualitatively similar")
print("(both decouple a heavy field at b->0) but quantitatively distinct")
print("(rank-shuffling vs kinetic uniformization).")
