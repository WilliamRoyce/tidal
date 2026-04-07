# Amplification Investigation Tracker

Investigation of amplification and suppression of the Gertsenshtein effect in nonminimal PGT+EM models. See `docs/tex/amplification_mechanism.tex` for the physics documentation.

## Status Key
- [ ] Not started
- [~] In progress
- [x] Complete
- [!] Blocked / needs discussion

## Priority 1 -- Foundation

- [x] Numerical Schur complement extraction (Phase A)
  - [x] Write `examples/torsion_gertsenshtein/schur_complement_analysis.py`
  - [x] Extract mu_eff and m2_eff at grid of (delta1, alpha2) points
  - [x] Map instability boundary alpha2_crit(delta1) from sweep data
  - [x] Extract coupling zero-crossing delta1_crit(alpha2) = 0.62
  - [!] Mini-solver validation: blocked by persistent tachyonic eigenvalues in h5-a1 block

- [x] Wolfram symbolic Schur complement (Phase A2)
  - [x] Write `examples/torsion_gertsenshtein/schur_complement_symbolic.wl`
  - [x] Build Block A (3x3) and Block B (5x5) symbolically
  - [x] Compute Inverse[S_cc] for each block
  - [x] Derive closed-form mu_eff and mu_ratio formulas
  - [x] Apply kinetic normalization (divide h5 row by -1/kappa^2)
  - [ ] Full numerical verification against Python (partial: signs match, magnitudes differ due to lambda-dependence)

- [!] Validate vs sweep data (Phase B) -- BLOCKED
  - [!] A_coupling = |mu_ratio|^2 at single k does NOT predict A_sweep
  - [!] Mini-solver approach blocked by tachyonic modes in reduced system
  - [!] Root cause: h5-a1 block has real eigenvalues +-lambda even after normalization
  - [!] The actual modal solver handles this via block-aware evolution + Hamiltonian energy
  - [ ] NEXT: investigate how Hamiltonian energy makes P bounded despite tachyonic modes

- [x] Documentation: `docs/tex/amplification_mechanism.tex`
  - [x] Coupling modification mechanism (Schur complement derivation)
  - [x] Distinction table vs Cembranos/Berlin/Lella
  - [x] Constraint block structure (Block A + Block B)
  - [x] delta1 symmetry explanation
  - [x] B0 scaling verification
  - [x] Observational implications (physical units, constraints, spectral signature)
  - [x] Full literature citations

## Priority 2 -- Precision & Scaling

- [x] B0 scaling check (Exp 3)
  - [x] Sweep B0 at amplification peak: C0 = 30,204 constant to 4 sig figs
  - [x] A = C0/C0_EM = 48.3 (plane-wave IC)

- [~] High-res boundary zoom (Exp 2)
  - [~] sweep_alpha2_hires.sh: 200 points, delta1=1.0, alpha2=-1.2:-0.5, B0=0.0001
  - [ ] Post-process and plot log10(A) vs alpha2
  - [ ] Extract scaling exponent near boundary

- [ ] Reparametrize plots
  - [x] Instability boundary extracted from sweep data
  - [ ] Re-plot heatmap with boundary-distance axis

- [x] Suppression valley phase analysis (Exp 5)
  - [x] Zero-crossing at delta1 = 0.62 (confirmed numerically)
  - [x] Coupling sign reversal verified

## Priority 3 -- Extensions

- [ ] Frequency dependence (Exp 4)
  - [x] Coupling ratio vs k table computed (1.55x to 7.48x)
  - [ ] Systematic k sweep via mini-solver or simulations
  - [ ] Characterize scaling exponent

- [ ] Three-sector survey (Exp 6)
  - [ ] 3D LHS over (alpha1, alpha2, alpha3)
  - [ ] Identify which sectors contribute

- [x] Observational implications
  - [x] Physical unit translation in amplification_mechanism.tex
  - [x] kappa_eff = sqrt(A)*kappa ~ 7*kappa
  - [x] Lella et al. bounds strengthened by factor sqrt(A)
  - [x] Blue-tilted spectral signature identified as unique prediction

## Key Discoveries

### 2026-04-07: Coupling reversal mechanism identified

Numerical Schur complement analysis of the h5-a1 block reveals:

1. **Not light-mediator 1/m^2**: The torsion mass matrix M_torsion is well-conditioned
   (smallest eigenvalue -0.5 to -3.5) at all amplified points. The M_torsion singularity
   is at alpha2=-0.25, far from the instability boundary (alpha2 ~ -1.1 to -1.8).

2. **Coupling reversal**: The primary amplification mechanism is SIGN REVERSAL of the
   effective h5-a1 coupling. At alpha2=-0.95, the coupling crosses zero at delta1 ~ 0.63.
   Below: slight suppression. Above: coupling reverses sign and grows.

3. **Quantitative coupling ratios** (alpha2=-0.82, delta1=1.0):
   - k=0.063: mu_eff/mu_GR = -1.55
   - k=0.503: mu_eff/mu_GR = -1.92
   - k=2.011: mu_eff/mu_GR = -7.48
   Enhancement is k-dependent (grows with wavenumber).

4. **Mass shift**: Torsion feedback also shifts the photon effective mass:
   m2_eff(a1) = -k^2 + shift(delta1^2, alpha_i, k). When shift > k^2,
   the photon becomes tachyonic (instability). The instability boundary
   shifts with delta1 (from alpha2=-1.78 at small delta1 to -1.14 at delta1=1).

5. **Two independent constraint blocks mediate the coupling**:
   - Block A (3x3): {t_0, t_15, t_22} -- spatial coupling
   - Block B (5x5): {t_3, t_4, t_7, t_8, t_18} -- time-derivative coupling

6. **Suppression valley** = exact destructive interference where mu_eff = 0.
   Zero-crossing mapped: at alpha2=-0.95, delta1_crit = 0.63.

### 2026-04-07: Exponential h5 growth and measurement normalization

The h5 graviton equation in the PGT nonminimal model has POSITIVE effective mass:
  d2_t h5 = +k^2/kappa^2 * h5 + B0^2/2 * h5 + coupling*a1
This gives exponential growth h5(t) ~ cosh(kt), NOT oscillatory wave propagation.
The 4x4 h5+a1 block from A_reduced shows P growing to ~0.56 at t=10, while the
actual simulation gives P=0.005.

The resolution: the conversion measurement normalizes P(t) = E_target(t)/E_source(0),
which uses the initial source energy, not the instantaneous. The modal solver correctly
handles the exponentially growing modes. The amplification factor A = P_torsion/P_GR
is well-defined because BOTH numerator and denominator have the same h5 growth.

**Consequence**: The Schur complement coupling ratio mu_eff/mu_GR does NOT directly
give the amplification factor A. The full eigenvalue structure (including growth rates)
must be accounted for. The correct approach to predict A is to evolve the 4x4 block
with identical ICs for both GR and torsion cases, then take the ratio of P_max values.

### 2026-04-07: B0 scaling verification PASSED

C0 = P/B0^2 is constant to 4 significant figures across B0 = 1e-6 to 1e-4:
  C0 = 30,204 (torsion-modified, delta1=1.0, alpha2=-0.82)
  C0_EM = 625 (Gertsenshtein baseline)
  A = C0/C0_EM = 48.3 (amplification factor at this parameter point)

Deviations at B0=1e-3 (+0.03%) and B0=1e-2 (-3.3%) from nonlinear backreaction.
This confirms the amplification is a coupling modification, NOT a nonlinear effect.

### 2026-04-07: Kinetic normalization bug found and fixed

The schur_complement_analysis.py was missing normalize_kinetic_coefficients(), causing
h_5 to appear with wrong-sign mass (+k^2 instead of -k^2). After fixing:
- GR baseline eigenvalues are purely imaginary (oscillatory)
- No more spurious exponential growth
- Coupling ratios and zero-crossing location unchanged (sign-independent)

### 2026-04-07: Two instability regions confirmed

Eigenvalue analysis with baseline-relative criterion reveals two instability regions:
- Lower boundary (Block A denominator pole): alpha2 -> -7/(4*kappa^2) = -1.75 at kappa=1
- Upper boundary (tachyonic photon): alpha2 varies with delta1 (~-0.80 at delta1=1)
- Stable window at delta1=1: approximately -1.14 < alpha2 < -0.76
