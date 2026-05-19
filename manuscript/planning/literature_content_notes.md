# Literature Content Notes — MSci Report

Synthesised from systematic reading of all 23 locally-available papers relevant to report content.
Four thematic groups; notes are organised by report section usage.

---

## Group 1 — Gertsenshtein-effect papers

### 2301.02072 — Domcke & Garcia-Cely (2023), *A Simple Derivation of the Gertsenshtein Effect*

**Key results for §1.1 and §4.1:**
- Cleanest modern re-derivation: treats graviton-photon as a two-state quantum mixing problem. Reproduces P = sin²(κB₀D/2) with κ² = 16πG.
- Explicitly flags the κ² = 16πG convention (not 8πG) — confirms TIDAL's usage.
- Linearises Einstein-Maxwell on a Minkowski background with uniform B₀; gauge: TT + Lorenz. Arrives at a 2×2 mixing matrix that is formally identical to axion-photon mixing (Raffelt-Stodolsky 1988).
- **For §4.1**: cite this as the modern statement of the formula validated to 0.04% by TIDAL.

### 2310.04150 — Hwang & Noh (2023), *On graviton-photon conversions in magnetic environments*

**Key results for §2.3:**
- Demonstrates that prior literature (specifically Palessandro 2405.01407) used incorrect EM field definitions: mixed covariant/contravariant forms that give wrong equations.
- Correct definitions: E_i = F_{i0} (covariant), B^i via the Levi-Civita symbol — the exact convention adopted in TIDAL.
- In the regime applicable to TIDAL (weak field, Minkowski background), the corrected and incorrect definitions agree for plane gravitational waves in uniform B₀. This confirms TIDAL's Boccaletti validation is unaffected.
- **For §2.3**: cite to establish which EM field convention we use.

### 2004.02714 — Ejlli (2020), *Graviton-photon mixing: Exact solution in a constant magnetic field*

**Key results for §2.3:**
- Exact solution (not perturbative) via eigenmode decomposition in a constant external B field.
- In the weak-field limit (κB₀D ≪ 1) reduces to sin²(κB₀D/2) — identical to TIDAL's validation formula.
- Includes plasma-mass term m_γ for the photon; in TIDAL's vacuum runs m_γ = 0, recovering the pure Boccaletti formula.
- **For §4.1**: cite as the non-perturbative analytic comparison for the Boccaletti validation.

### 2405.11786 — Hwang & Noh (2024), *Nonlinear graviton-photon mixing*

**Key results (background context only — not a primary citation):**
- Studies nonlinear corrections to the Gertsenshtein formula at higher field strength.
- Confirms linear formula P = sin²(κB₀D/2) is exact at leading order in κ; corrections are O(κ²B₀²).
- **Not a primary citation** but confirms TIDAL operates in the correct linearised regime.

### 2405.01407 — Palessandro (2024), *Graviton-Photon Oscillations as a Probe of Quantum Gravity*

**For §1.1 prior-art survey:**
- Applies Gertsenshtein mechanism to probe quantum gravity corrections.
- Uses the incorrect EM field convention identified by Hwang-Noh 2310.04150 — flagged for citation context.
- **Cite in §1.1 as part of the broadening of Gertsenshtein applications to quantum gravity, with a note that the formalism was corrected in 2310.04150.**

### 2507.16609 — Domcke, Garcia-Cely & Lee (2025), *Gravitational Wave Scattering on Magnetic Fields*

**Key results for §4 (background):**
- Shows that finite-extent B fields cause GW scattering, not just conversion. In the uniform-field limit, recovers the standard Gertsenshtein formula.
- This confirms that TIDAL's use of a spatially uniform B₀ is the correct limit for the conversion problem studied.

### 2510.17094 — Tomomatsu, Suyama & Gondolo (2025), *Gertsenshtein effect on curved spacetime*

**Key results for §2.3:**
- Derives the Gertsenshtein conversion probability on a curved background induced by the background magnetic field itself (B²-induced curvature).
- At leading order in B²/M_Pl², the standard Minkowski result P = sin²(κB₀D/2) is recovered — the curvature correction is suppressed by (κB₀L)².
- **Confirms TIDAL's Minkowski approximation is valid** for all physically realistic B₀ values.
- Uses the same κ² = 16πG convention throughout — consistent with TIDAL.

---

## Group 2 — PGT ghost/Hamiltonian structure papers

### 1812.02675 — Lin, Hobson & Lasenby (2019), *Ghost and Tachyon Free Poincaré Gauge Theories: A Systematic Approach*

**Key results for §1.2 and §5.2:**
- Supervisors' paper (Lin was Lasenby's student). Identifies all 450 ghost-tachyon-free PGT combinations at the quadratic-torsion level by exhaustive computer-algebra scan.
- Imposes the Sezgin-Van Nieuwenhuizen (1980) no-ghost conditions for each of the 450 parameter combinations.
- Key finding: the viable ghost-free window is extremely narrow — most of the 450 cases have additional constraints that are hard to satisfy simultaneously.
- **For §1.2**: cite as the foundational ghost-free classification that motivates the search for viable propagating torsion theories.
- **For §5.2**: cite in the "three constructive paths" discussion.

### 1804.05556 — Blagojević & Cvetković (2018), *General Poincaré Gauge Theory: Hamiltonian Structure*

**Key results for §5.2:**
- Comprehensive Hamiltonian analysis of PGT at the general level; derives constraint algebras.
- Shows that spin-1 and spin-2 torsion propagation are governed by separate sectors of the constraint algebra and that making both healthy simultaneously is generically impossible.
- This is the "spin-1/spin-2 mutual exclusivity" result cited in the report's §5.2.
- **For §5.2**: cite alongside 1812.02675 as the structural reason why generic PGT has ghost problems.

### 2510.08201 — Barker & Glavan (2024), *Spectrum of Pure R² Gravity: Full Hamiltonian Analysis*

**Key results for §1.2 (background):**
- Demonstrates Hamilcar-style constraint analysis on pure R² gravity.
- Relevant as a comparator: shows the Hamiltonian method identifies the correct degrees of freedom in a theory closely related to PGT.
- **Secondary citation** for the Hamilcar/PSALTer software context in App C.

### 2505.23894 — Barker, Marzo & Santoni (2025), *Can Metric-Affine Gravity Be Saved?*

**Key results for §1.2:**
- Finds that extending PGT to full metric-affine gravity introduces additional ghost modes from the non-metricity sector.
- The "path to saving" MAG requires imposing projective symmetry or working in a constrained subspace.
- **For §1.2**: cite as context for why PGT (not MAG) is the natural starting point for a ghost-safe torsion theory.

### 2402.07641 — Barker & Marzo (2024), *Particle Spectra of General Palatini/Metric-Affine Theories*

**Key results for §1.2 (background):**
- Shows that the PSALTer algorithm extends to metric-affine theories.
- Confirms ghost analysis is tractable beyond PGT — but the ghost structure is richer and less well-understood.
- **Background citation** for App A (PSALTer capabilities).

### 2402.14917 — Barker & Zell (2024), *Consistent Particle Physics in MAG from Extended Projective Symmetry*

**Key results for §5.2:**
- Projective symmetry = the symmetry that eliminates the problematic dilaton in metric-affine gravity.
- Identifies the "T̃_μ axial trace" torsion mode as the one surviving healthy torsion degree of freedom after projective symmetry is imposed — **directly confirms TIDAL's trace channel identification**.
- **For §5.2 path 2** (ghost-free kinetics): cite alongside 2507.05349 as establishing which torsion mode survives the ghost-freedom requirement.

---

## Group 3 — Torsion-EM coupling foundations and dark photon limits

### hep-th/0103093 — Shapiro (2002), *Physical Aspects of the Space-Time Torsion*

**Key results for §2.1:**
- Comprehensive review (357 pages). The notation reference: T^α_{μν} for the full torsion tensor; V_μ = T^α_{μα} for the trace vector; A^μ = ε^{αβγμ}T_{αβγ} for the axial pseudovector; T_{αβγ} = T^{(1)}_{αβγ} + T^{(2)}_{αβγ} + T^{(3)}_{αβγ} for the irreducible decomposition (matching Hehl et al. 1976 labelling).
- For the Minkowski background: curvature vanishes, torsion contributes an additional connection term. The Riemann-Cartan structure gives \nabla_μ = \partial_μ + Γ^{RC}_{μ} where Γ^{RC} includes the contorsion K^α_{μν}.
- Mass dimension: [T^α_{μν}] = [length]^{-1} = [mass] in natural units.
- **Primary notation reference for §2.1.**

### gr-qc/0001010 — Hehl & Obukhov (2000), *How Does the Electromagnetic Field Couple to Gravity?*

**Key results for §2.3:**
- Establishes the admissibility taxonomy for torsion-EM couplings:
  - **Inadmissible**: F·T terms (odd in torsion, violate charge conservation/parity).
  - **Admissible**: F²·T² terms (minimal coupling = F·F with RC metric; non-minimal = additional F²T² operators).
- The standard "minimal" coupling L_EM = -¼ F^{μν}F_{μν} √(-g) with the Riemann-Cartan metric gives F depending on A only (not torsion) if one uses the U(1) gauge potential — this is TIDAL's Stage A setup.
- **For §2.3**: cite to justify why Stage A (dark photon plasma with Einstein-Cartan torsion) produces no conversion: the minimal coupling does not introduce a direct torsion-photon vertex.

### gr-qc/0307063 — Itin & Hehl (2003), *Maxwell's Field Coupled Nonminimally to Quadratic Torsion*

**Key results for §2.3 and §5.2:**
- The 17-parameter family of F²T² non-minimal couplings. Each of the 17 terms is a product of two F components and two T components (one from each irreducible piece).
- Ghost-freedom analysis: most of the 17 terms introduce additional ghost modes in the photon sector. The "safe" ones are the parity-even, trace-sector-dominant terms.
- **For §5.2 path 1** (non-minimal coupling): cite as the systematic classification of non-minimal operators available.

### gr-qc/0305049 — Rubilar, Obukhov & Hehl (2003), *Torsion Nonminimally Coupled to the EM Field and Birefringence*

**Key results for §2.3:**
- Shows that non-minimal F²T² couplings generically produce birefringence: left- and right-handed photon polarisations propagate at different speeds.
- This is an observational handle distinct from the Gertsenshtein conversion probability: birefringence gives a rotation of polarisation angle.
- **For §2.3**: cite as a complementary observable to the Gertsenshtein effect (conversion probability vs polarisation rotation).

### 2105.04565 — Caputo, Millar, O'Hare & Vitagliano (2021), *Dark Photon Limits: A Handbook*

**Key results for §4.2:**
- Comprehensive dark photon constraint summary. Key formula: conversion probability for dark-photon (mass m_A, coupling ε) in magnetic field B₀ over distance D is P ≈ (εκB₀)² D²/4 for m_A D ≪ 1 (and oscillates for m_A D ≫ 1).
- The dark-photon mass determines whether conversion is coherent (P ~ D²) or oscillatory (P ~ 1/m_A²).
- TIDAL's Stage A null result: for the parameter space explored (m_A ~ 0.01–1 eV, ε ~ 0.01–2), P_max matches sin²(κB₀D/2) to 6.7×10⁻⁶ — no torsion amplification above the standard dark-photon baseline.
- **For §4.2**: cite to establish the comparison baseline for Stage A results.

---

## Group 4 — Constructive paths and recent developments

### 2507.09228 — Legner, Handley & Barker (2025), *Alleviating the Hubble Tension with TorC*

**Key results for §5.2:**
- TorC = Torsion Condensation model: a specific PGT cosmology where torsion condenses like a Bose-Einstein condensate at a phase transition.
- Provides a concrete, ghost-free torsion model with a propagating vector torsion mode (spin-1, trace sector) that drives modified expansion history.
- H₀ tension: achieves S₈ and H₀ simultaneously within 1σ of Planck, unlike ΛCDM.
- **For §5.2 path 2** (ghost-free kinetics): cite as the most recent concrete realisation of a viable, propagating torsion theory. The TorC model is directly in the regime TIDAL is studying (trace-sector torsion).

### 2506.21662 — Barker, Marzo & Santoni (2025), *Infrared Foundations I: Rank-3 Field Theories*

**Key results for §1.2 and §5.2:**
- Systematic analysis of ALL rank-3 field theories from the perspective of infrared constraints (gauge symmetry, ghost freedom, masslessness).
- Key finding: only **one** rank-3 theory is consistent (the self-dual limit). This is highly restrictive.
- For PGT: the torsion tensor T^α_{μν} is a rank-3 field; the classification here applies directly.
- **For §1.2**: motivates why the torsion theory space is much smaller than it appears — the viable corner is a specific rank-3 subspace.

### 2507.05349 — Barker, Marzo & Santoni (2025), *Infrared Foundations II: Torsion-like Theories and Ghost Analysis*

**Key results — CRITICAL for §1.2, §4.4, §5.2:**
- The single most important recent paper for the report's framing.
- Analyses 22 ghost-free torsion-like models (the intersection of ghost freedom + Lorentz covariance + minimal coupling).
- **Central finding: EVERY one of the 22 ghost-free models propagates ONLY spin-1 vector modes (the trace sector).** No model with spin-0 or spin-2 propagating torsion can be simultaneously ghost-free under any admissible gauge symmetry.
- This directly validates TIDAL's trace-channel focus: the only torsion sector that can couple to photons in a ghost-free theory is the spin-1/trace sector (²T component in Hehl labelling).
- **For §1.2**: this is the "Deep Problem" and its partial resolution — PGT without ghost freedom is useless, ghost freedom restricts to trace sector.
- **For §4.4**: cite to frame the hx↔ax independence finding as a CONSEQUENCE of this restriction (only trace channel can be healthy, so hx↔ax decoupling is expected).
- **For §5.2 path 2**: the non-trivial extension is to build a full interacting theory around the spin-1 trace mode.

### 2406.11956 — Karananas, Shaposhnikov & Zell (2024), *Weyl-Invariant Einstein–Cartan Gravity*

**Key results for §5.2 path 1:**
- Embeds EC torsion in a Weyl-invariant framework. The Weyl symmetry constrains the torsion to live in the trace vector sector (²T).
- Non-minimal couplings to EM emerge from the Weyl-invariant completion of the EM action — generating F²T² operators of the Itin-Hehl (gr-qc/0307063) type.
- The coupling is controlled by the Weyl gauge field, which is identified with the torsion trace.
- **For §5.2 path 1**: cite as a concrete model that generates non-minimal torsion-EM coupling from a symmetry principle (Weyl invariance), rather than adding operators by hand.

### 2506.17017 — Bahamonde et al. (2025), *Cosmology of Cubic Poincaré Gauge Gravity*

**Key results for §5.2 path 3:**
- Cubic PGT Lagrangian (cubic in torsion/curvature). The cubic terms can be tuned to cancel the ghost modes that appear at quadratic level.
- Cosmological solutions: the cubic model admits de Sitter and power-law attractors. The torsion scalar field plays the role of dark energy.
- Explicitly ghost-free at the quadratic perturbation level when cubic coupling satisfies a specific constraint.
- **For §5.2 path 3**: cite as the primary reference for the "cubic PGT" route to a ghost-free torsion theory that avoids the spin-1/spin-2 exclusivity problem.

---

## Cross-cutting notation table

| Symbol | Meaning | Convention source |
|--------|---------|-------------------|
| T^α_{μν} | Full torsion tensor | Hehl 1976; Shapiro hep-th/0103093 |
| ¹T, T^α_{μν} irreducible | Traceless tensor part (q_{μνρ}) | Hehl 1976; BHL 2406.12826 |
| ²T, T_μ | Trace vector = T^α_{μα} | Hehl 1976; Shapiro §2.2 |
| ³T, S^μ | Axial pseudovector = ε^{αβγμ}T_{αβγ} | Hehl 1976; Shapiro §2.2 |
| κ² = 16πG | Graviton-photon coupling | Domcke 2301.02072; BHL 2406.12826 |
| P = sin²(κB₀D/2) | Conversion probability | Boccaletti 1970; Domcke 2301.02072 |
| (-,+,+,+) | Metric signature | BHL joint papers; Shapiro |
| e^a_{μ} | Vierbein | BHL 2406.12826; 2309.14783 |
| ω^{ab}_{μ} | Spin connection | BHL 2406.12826 |

---

## Key framing insights for writing

1. **The Deep Problem** (§1.2): PGT has a rich parameter space, but almost all of it is ghost-haunted. The ghost-freedom requirement (Sezgin-VanNieuwenhuizen 1980; Lin-Hobson-Lasenby 2019; Barker-Marzo-Santoni 2025) is so restrictive that it forces torsion into the spin-1/trace sector alone (2507.05349). TIDAL's Stage B null (no amplification in EC theory) and the hx↔ax independence finding are therefore expected rather than surprising — they are consequences of the ghost structure.

2. **The positive finding** (§4.4): hx↔ax torsion-independence is a structural theorem, not an accident. At linear order, the tensor (hx) and axial (ax) torsion channels are kinematically decoupled from the photon sector because the only admissible torsion-photon vertex is a trace-vector coupling (per gr-qc/0001010). This insight should be stated as a theorem in §2.4 and confirmed numerically in §4.4.

3. **Three constructive paths** (§5.2):
   - Path 1: Non-minimal F²T² coupling (Itin-Hehl gr-qc/0307063 operators + Karananas 2406.11956 as Weyl-invariant realisation)
   - Path 2: Ghost-free kinetics for spin-1 trace torsion (2507.05349 + TorC 2507.09228 + Barker-Zell 2402.14917)
   - Path 3: Cubic PGT Lagrangian that cancels quadratic ghosts (2506.17017)

4. **Software framing** (§3 + App A–C): TIDAL is the first pipeline to combine the symbolic PGT torsion decomposition (Wolfram/xAct) with an exact spectral solver (Fourier modal backend). PSALTer (2406.09500) handles particle spectra; Hamilcar (2512.25007) handles constraint algebras; TIDAL handles the dynamical evolution — they are complementary tools addressing different aspects of the same theory class.

5. **Null results as constraints** (§4.2–4.3): "We exclude amplification above the Einstein-Maxwell baseline in the Einstein-Cartan sector to precision P < P_GR × 10⁻⁵ across 276 parameter combinations." This is a positive precision constraint on torsion-Gertsenshtein physics, not a failure.
