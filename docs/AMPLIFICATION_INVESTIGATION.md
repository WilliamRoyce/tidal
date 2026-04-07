# Amplification Investigation Tracker

Investigation of amplification and suppression of the Gertsenshtein effect in nonminimal PGT+EM models. See `docs/tex/amplification_mechanism.tex` for the physics documentation.

## Status Key
- [ ] Not started
- [~] In progress
- [x] Complete
- [!] Blocked / needs discussion

## Priority 1 -- Foundation

- [~] Numerical Schur complement extraction (Phase A)
  - [ ] Write `examples/torsion_gertsenshtein/schur_complement_analysis.py`
  - [ ] Extract mu_eff and m2_eff at grid of (delta1, alpha2) points
  - [ ] Map instability boundary alpha2_crit(delta1) from eigenvalue crossing
  - [ ] Extract coupling zero-crossing delta1_crit(alpha2) (suppression valley)
  - [ ] Compute A_analytic = |mu_eff/mu_GR|^2 for validation

- [ ] Wolfram symbolic Schur complement (Phase A2)
  - [ ] Write `examples/torsion_gertsenshtein/schur_complement_symbolic.wl`
  - [ ] Import JSON, parse equations to symbolic form
  - [ ] Build Block A (3x3: t_0, t_15, t_22) mass matrix symbolically
  - [ ] Build Block B (5x5: t_3, t_4, t_7, t_8, t_18) symbolically
  - [ ] Compute Inverse[S_cc] for each block
  - [ ] Derive mu_eff(delta1, alpha_i, kappa, k) closed form
  - [ ] Derive m2_eff(delta1, alpha_i, kappa, k) closed form
  - [ ] Verify numerically against Python Schur complement
  - [ ] Export to TeX via TeXForm[]

- [ ] Validate vs sweep data (Phase B)
  - [ ] Load hires heatmap results.csv
  - [ ] Compare A_analytic vs A_numerical at all valid points
  - [ ] Scatter plot + residual analysis

- [ ] Documentation: `docs/tex/amplification_mechanism.tex`
  - [ ] Coupling modification mechanism (Schur complement derivation)
  - [ ] Distinction table vs Cembranos/Berlin/Lella
  - [ ] Mass spectrum connection and constraint block structure
  - [ ] delta1 symmetry and suppression valley explanation
  - [ ] Full literature citations

## Priority 2 -- Precision & Scaling

- [ ] B0 scaling check (Exp 3)
  - [ ] Sweep B0 at amplification peak
  - [ ] Verify C0 = P/B0^2 is constant

- [ ] High-res boundary zoom (Exp 2)
  - [ ] delta1=1.0, sweep alpha2 from -1.2 to -0.5, 200 points, B0=1e-6
  - [ ] Plot log10(A) vs (alpha2 - alpha2_crit)
  - [ ] Extract scaling exponent

- [ ] Reparametrize plots
  - [ ] Compute instability boundary alpha2_crit(delta1) for boundary-distance axis
  - [ ] Add sector mass transform (for literature comparison)
  - [ ] Re-plot heatmap with boundary-distance axis

- [ ] Suppression valley phase analysis (Exp 5)
  - [ ] Extract effective coupling sign at delta1 ~ 0.63
  - [ ] Confirm destructive interference (zero-crossing)

## Priority 3 -- Extensions

- [ ] Frequency dependence (Exp 4)
  - [ ] Sweep k at amplification peak
  - [ ] Characterize spectral signature

- [ ] Three-sector survey (Exp 6)
  - [ ] 3D LHS over (alpha1, alpha2, alpha3)
  - [ ] Identify which sectors contribute

- [ ] Observational implications
  - [ ] Physical unit translation
  - [ ] Mass bounds from non-observation
  - [ ] Spectral signature prediction

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
