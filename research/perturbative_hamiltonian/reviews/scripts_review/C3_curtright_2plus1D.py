# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
"""
C3: Re-do Agent J's gauge-invariance verification in 2+1D (and 3+1D
sketch), to check if the Stückelberg gauge invariance survives beyond
the 1+1D toy.

In D dimensions, a (2,1) Curtright tensor T_{μν|ρ} (antisymmetric in
[μν], no full antisymmetry) has D^2(D-1)/2 independent components
before cyclic constraint, and (D-1)·D·(D+1)/3 after cyclic.

  D=2: 2 components before cyclic, 1 after.   (the 1+1D toy)
  D=3: 9 before, 8 after.                     (2+1D)
  D=4: 24 before, 16 after.                   (3+1D, full PGT)

Stückelberg auxiliaries: h symmetric (D(D+1)/2), b antisymmetric
(D(D-1)/2), a vector (D).
  D=2: h=3, b=1, a=2. Total = 6 auxiliaries for 1 physical.
  D=3: h=6, b=3, a=3. Total = 12 auxiliaries for 8 physical.
  D=4: h=10, b=6, a=4. Total = 20 auxiliaries for 16 physical.

Test: in D=3, build a momentum-space (k0, k1, k2) plane-wave ansatz,
verify δF̊ = 0 under the (s, β, α) Stückelberg shifts.
"""

import sys

import sympy as sp
from sympy import I, Rational, symbols

D = int(sys.argv[1]) if len(sys.argv) > 1 else 3  # spacetime dimension
print("=" * 70)
print(f"C3: Curtright Stückelberg gauge invariance in D={D}")
print("=" * 70)
print()

# Plane-wave ansatz: k_μ
k = [symbols(f"k{i}", real=True) for i in range(D)]

# T[mu, nu, rho]: amplitude of the Curtright tensor.
# Antisymmetric in [mu, nu]. Free amplitudes T_amp[mu][nu][rho] for mu<nu.
T_amp = {}
for mu in range(D):
    for nu in range(mu + 1, D):
        for rho in range(D):
            T_amp[mu, nu, rho] = symbols(f"T_{mu}{nu}_{rho}", real=False)


def T(mu, nu, rho):
    if mu == nu:
        return 0
    if mu < nu:
        return T_amp[mu, nu, rho]
    return -T_amp[nu, mu, rho]


# Auxiliary h_{μν} (symmetric)
h_amp = {}
for mu in range(D):
    for nu in range(mu, D):
        h_amp[mu, nu] = symbols(f"h_{mu}{nu}", real=False)


def h(mu, nu):
    if mu <= nu:
        return h_amp[mu, nu]
    return h_amp[nu, mu]


# Auxiliary b_{μν} (antisymmetric)
b_amp = {}
for mu in range(D):
    for nu in range(mu + 1, D):
        b_amp[mu, nu] = symbols(f"b_{mu}{nu}", real=False)


def bf(mu, nu):
    if mu == nu:
        return 0
    if mu < nu:
        return b_amp[mu, nu]
    return -b_amp[nu, mu]


# Auxiliary a_μ (vector)
a_amp = [symbols(f"a_{mu}", real=False) for mu in range(D)]


# Field strength F̊_{μν|ρ} = T_{μν|ρ} - 2 ∂_[μ h_{ν]ρ} - 2 ∂_[μ b_{ν]ρ} + 2 ∂_ρ b_{μν} - 2 ∂_ρ ∂_[μ a_{ν]}
# In momentum: ∂_μ → i k_μ.
def Fring(mu, nu, rho):
    res = T(mu, nu, rho)
    # -2 ∂_[μ h_{ν]ρ}  = -(∂_μ h_{νρ} - ∂_ν h_{μρ})  = -i(k_μ h_{νρ} - k_ν h_{μρ})
    res += -I * (k[mu] * h(nu, rho) - k[nu] * h(mu, rho))
    # -2 ∂_[μ b_{ν]ρ}  similarly
    res += -I * (k[mu] * bf(nu, rho) - k[nu] * bf(mu, rho))
    # +2 ∂_ρ b_{μν}    (note: full coefficient already in front of D[ρ] b_{μν})
    res += 2 * I * k[rho] * bf(mu, nu)
    # -2 ∂_ρ ∂_[μ a_{ν]} = -2·∂_ρ · (1/2)(∂_μ a_ν - ∂_ν a_μ) = -(i k_ρ)·(i k_μ a_ν - i k_ν a_μ)
    res += -((I * k[rho]) * (I * k[mu] * a_amp[nu] - I * k[nu] * a_amp[mu]))
    return sp.expand(res)


# Stückelberg gauge parameters: s_{μν} symmetric, β_{μν} antisymmetric, α_μ vector
s_amp = {}
for mu in range(D):
    for nu in range(mu, D):
        s_amp[mu, nu] = symbols(f"s_{mu}{nu}", real=False)


def s(mu, nu):
    if mu <= nu:
        return s_amp[mu, nu]
    return s_amp[nu, mu]


beta_amp = {}
for mu in range(D):
    for nu in range(mu + 1, D):
        beta_amp[mu, nu] = symbols(f"beta_{mu}{nu}", real=False)


def beta(mu, nu):
    if mu == nu:
        return 0
    if mu < nu:
        return beta_amp[mu, nu]
    return -beta_amp[nu, mu]


alpha_amp = [symbols(f"alpha_{mu}", real=False) for mu in range(D)]


# Apply Stückelberg shifts:
#   δT_{μν|ρ} = 2 ∂_[μ s_{ν]ρ} + 2 ∂_[μ β_{ν]ρ} - 2 ∂_ρ β_{μν} + 2 ∂_ρ ∂_[μ α_{ν]}
#   δh_{μν} = s_{μν}
#   δb_{μν} = β_{μν}
#   δa_μ = α_μ
def deltaT(mu, nu, rho):
    return sp.expand(
        I * (k[mu] * s(nu, rho) - k[nu] * s(mu, rho))
        + I * (k[mu] * beta(nu, rho) - k[nu] * beta(mu, rho))
        + (-2) * I * k[rho] * beta(mu, nu)
        + 2
        * (I * k[rho])
        * Rational(1, 2)
        * (I * k[mu] * alpha_amp[nu] - I * k[nu] * alpha_amp[mu])
    )


# δF̊_{μν|ρ}: substitute T → T + δT, h → h + s, b → b + β, a → a + α
def deltaFring(mu, nu, rho):
    Fbase = Fring(mu, nu, rho)
    # Substitute one at a time and subtract the base
    subs = {}
    # h shift
    for (mm, nn), hsym in h_amp.items():
        if mm == nn:
            subs[hsym] = hsym + s(mm, nn)
        else:
            subs[hsym] = hsym + s(mm, nn)
    # b shift
    for (mm, nn), bsym in b_amp.items():
        subs[bsym] = bsym + beta(mm, nn)
    # a shift
    for mu0 in range(D):
        subs[a_amp[mu0]] = a_amp[mu0] + alpha_amp[mu0]
    # T shift (apply to T_amp directly using the explicit formula)
    for (mm, nn, rr), tsym in T_amp.items():
        subs[tsym] = tsym + deltaT(mm, nn, rr)
    Fnew = sp.expand(Fbase.xreplace(subs))
    return sp.expand(Fnew - Fbase)


# Check all (μ < ν, ρ) component variations
print("Checking δF̊_{μν|ρ} for all antisymmetric pairs and ρ values...")
all_zero = True
nonzero_components = []
for mu in range(D):
    for nu in range(mu + 1, D):
        for rho in range(D):
            df = deltaFring(mu, nu, rho)
            if df != 0:
                # try simplify
                df_s = sp.simplify(df)
                if df_s != 0:
                    all_zero = False
                    nonzero_components.append((mu, nu, rho, df_s))

if all_zero:
    print(f"  ✓ ALL {D * (D - 1) // 2 * D} components of δF̊ vanish identically.")
    print(f"  -> Curtright Stückelberg gauge invariance holds in D={D}.")
else:
    print(f"  ✗ FAIL: {len(nonzero_components)} components have nonzero variation:")
    for mu, nu, rho, df in nonzero_components[:5]:
        print(f"    δF̊_{{{mu}{nu}|{rho}}} = {df}")

# Component count check
print()
print(f"Component count check (D={D}):")
print(
    f"  T_{{μν|ρ}} amplitudes (mu<nu, all rho): {len(T_amp)}  (expected D²(D-1)/2 = {D * D * (D - 1) // 2})"
)
print(f"  Auxiliary h: {len(h_amp)}  (expected D(D+1)/2 = {D * (D + 1) // 2})")
print(f"  Auxiliary b: {len(b_amp)}  (expected D(D-1)/2 = {D * (D - 1) // 2})")
print(f"  Auxiliary a: {len(a_amp)}  (expected D = {D})")
