# Research References for TIDAL

These references inform design decisions and should be cited where appropriate in reports and documentation.

## Scientific Codebases (Architecture & Feature Influence)

| Reference                                                                                                                                                                                                          | Relevance to TIDAL                                                                                        |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------- |
| **Dedalus Project** — Burns et al., "Dedalus: A Flexible Framework for Numerical Simulations with Spectral Methods", Phys. Rev. Research 2, 023068 (2020). [arXiv:1905.10388](https://arxiv.org/abs/1905.10388)    | Spectral spatial discretization, native eigenvalue problems, HDF5 analysis output, on-the-fly diagnostics |
| **FEniCS / FEniCSx** — Baratta et al., "DOLFINx: The next generation FEniCS problem solving environment" (2023). [fenicsproject.org](https://fenicsproject.org/)                                                   | Variational form DSL, automated dimensional analysis for PDEs, adaptive mesh refinement                   |
| **MEEP** — Oskooi et al., "MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method", Computer Physics Communications (2010). [meep.readthedocs.io](https://meep.readthedocs.io/) | PML absorbing boundaries, mode decomposition (EME), Poynting flux monitors, convergence testing           |
| **Cadabra** — Peeters, "Cadabra: A field-theory motivated approach to computer algebra", Computer Physics Communications (2007). [cadabra.science](https://cadabra.science/)                                       | Tensor algebra symbolic computation, field theory notation                                                |
| **xAct** — Martín-García et al., "xAct: Efficient tensor computer algebra for Mathematica". [xact.es](https://xact.es/articles.html)                                                                               | Core symbolic tensor algebra (already used by TIDAL)                                                      |
| **xPert** — Brizuela et al., "xPert: Computer algebra for metric perturbation theory", General Relativity and Gravitation 41 (2009). [ResearchGate](https://www.researchgate.net/publication/1740524)              | Metric perturbation linearization (already used by TIDAL)                                                 |
| **py-pde** — Zwicker, "py-pde: A Python package for solving partial differential equations", Journal of Open Source Software 5(48), 2158 (2020). [GitHub](https://github.com/zwicker-group/py-pde)                 | Original PDE backend; FD stencil conventions retained in TIDAL's native operators                          |

## Gertsenshtein Effect & Wave Conversion Physics

| Reference                                                                                                                                                                                                            | Relevance                                                                                                     |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| **Gertsenshtein (1962)** — "Wave resonance of light and gravitational waves", JETP 14, 84                                                                                                                            | Original prediction of graviton-photon conversion in B-field                                                  |
| **Palessandro & Rothman (2023)** — "A simple derivation of the Gertsenshtein effect", *Phys. Dark Univ.* **40**, 101187, [arXiv:2301.02072](https://arxiv.org/abs/2301.02072)                                        | ⚠️ **Formula is wrong by √(4π)**: they use a non-canonical `√G·B₀` normalization, giving `P = sin²(√G B₀ D)`. Do **not** copy it — see `docs/tex/gertsenshtein_formula.tex` §4 for TIDAL's verified form. (Previously mis-attributed here to Domcke & Garcia-Cely; corrected 2026-08-30.) |
| **Hwang & Noh (2023)** — "On graviton-photon conversions in magnetic environments", [arXiv:2310.04150](https://arxiv.org/abs/2310.04150)                                                                              | Proper EM field definitions in curved spacetime; formulation-dependent graviton mass; critique of ad-hoc derivations |
| **Berlin et al. (2024)** — "Numerical analysis of resonant axion-photon mixing", [arXiv:2405.08865](https://arxiv.org/abs/2405.08865)                                                                                | Eq.(309): Schrödinger-like two-state mixing `2ik∂ₓẼ ≃ -(mₐ²-ωₚ²)Ẽ - coupling`; Eq.(335): resonant conversion probability (= Raffelt-Stodolsky 1988 on resonance). Directly analogous to TIDAL's plasma-Gertsenshtein approach. |
| **Domcke, Garcia-Cely & Lee (2025)** — "Gravitational wave conversion into photons in magnetized plasmas", [arXiv:2507.16609](https://arxiv.org/abs/2507.16609)                                                      | Eq.(216): massless GW-sourced photon wave eq `□Aₕ = -j_eff` (j_eff defined Eq.(187)). Eq.(1484) §"Towards including Medium Effects": `[□-μ²(x)]Aₕ = -j_eff` with plasma mass — EOM-level counterpart to TIDAL's perturbation-Proca term. |
| **Ejlli et al. (2019)** — "Upper limits on the amplitude of ultra-high-frequency gravitational waves", Eur. Phys. J. C 79, 1032 (2019). [Springer](https://link.springer.com/article/10.1140/epjc/s10052-019-7542-5) | Experimental bounds on Gertsenshtein conversion                                                               |
| **Ejlli (2020)** — "Graviton-photon mixing: exact solution in a constant magnetic field", [arXiv:2004.02714](https://arxiv.org/abs/2004.02714) | First exact (non-perturbative) Gertsenshtein solution |
| **Palessandro (2024)** — "Graviton-Photon Oscillations as a Probe of Quantum Gravity", [arXiv:2405.01407](https://arxiv.org/abs/2405.01407) | Extended Gertsenshtein treatment (**contains normalization error** corrected by TIDAL) |
| **Hwang & Noh (2024)** — "Graviton-photon conversions in Euler-Heisenberg nonlinear electrodynamics", [arXiv:2405.11786](https://arxiv.org/abs/2405.11786) | Nonlinear QED corrections to Gertsenshtein; chiral GW propagation |
| **Adler (1971)** — "Photon splitting and photon dispersion in a strong magnetic field", Ann. Phys. 67, 599. [DOI](https://doi.org/10.1016/0003-4916(71)90154-0) | QED vacuum birefringence (Heaviside-Lorentz, natural units): `n_∥ − 1 = (14α²/45)(B/B_c)²`, `n_⊥ − 1 = (8α²/45)(B/B_c)²`. Ratio 7/4 + individual coefficients analytically reproduced by TIDAL's EH Lagrangian under the convention mapping `ρ = 8c₁ = 4α²/45`, `σ = 2c₂ = 7α²/180` (so σ/ρ = 7/16). See [perturbative_reduction.tex §Validation](tex/perturbative_reduction.tex). |
| **Dunne (2004)** — "Heisenberg–Euler effective Lagrangians: Basics and extensions", [arXiv:hep-th/0406216](https://arxiv.org/abs/hep-th/0406216) | Pedagogical review. Eq. "lightlight" (line 332 of local TeX source) gives `S^(1) = (e⁴/(360π²m⁴))∫[(E²-B²)² + 7(E·B)²]`. Uses halved-dual `F̃ = ½ε·F`; TIDAL uses un-halved `ε·F·F` (a factor 4 on squaring), so Dunne `c₂` → TIDAL `σ = 2c₂`. Heaviside-Lorentz `e⁴ = 16π²α²` cancels the π² in Dunne's prefactor to give `c₁ = α²/(90m⁴)` (no π² — earlier "convention gap" claim was a derivation error). |
| **Graviton-photon oscillation in general modified gravity (2023)** — [arXiv:2302.08186](https://arxiv.org/abs/2302.08186) | Gertsenshtein in Horndeski, massive gravity, cosmic backgrounds |

## Primary Classical References (pre-arXiv)

These papers establish the foundational theory used by TIDAL's Gertsenshtein
and plasma-Gertsenshtein examples. They predate arXiv and are not stored in
`Literature/`; the key formulas are transcribed here with cross-references
to the TIDAL files that implement them.

### Gertsenshtein 1962 — original graviton-photon conversion

**Citation:** M. E. Gertsenshtein, "Wave resonance of light and gravitational
waves", Sov. Phys. JETP **14**, 84 (1962). Russian original: ZhETF **41**, 113.

**Key result:** A gravitational wave propagating through a static magnetic
field `B₀` coherently converts into an electromagnetic wave (and vice versa)
with probability

```
P(D) = sin²(κ · B₀ · D / 2)     [κ = √(8πG), D = propagation distance]
```

**Implemented in:** [examples/gertsenshtein/theory.toml](examples/gertsenshtein/theory.toml)
and validated by [sweep_B0.sh](examples/gertsenshtein/sweep_B0.sh) (RMS < 0.015
vs analytical).

**Canonical normalization:** The √(8πG) factor is often mis-stated in the
literature (see Palessandro 2024 critique in arXiv:2405.01407). Dandoy et al.
2024 ([Literature/2406.17853](Literature/2406.17853)) reproduces the correct
normalization independently, which TIDAL uses.

### Boccaletti et al. 1970 — localized B-field integral form

**Citation:** D. Boccaletti, V. De Sabbata, P. Fortini, C. Gualdi, "Conversion
of photons into gravitons and vice versa in a static electromagnetic field",
Nuovo Cimento B **70**, 129 (1970).

**Key result:** For a magnetic field profile `B(z)` localized in space, the
massless-vacuum conversion probability generalizes from the uniform-field
`sin²(κB₀D/2)` to

```
P = sin²((κ/2) · ∫ B(z) dz)     [integral over propagation path]
```

**Implemented in:** [docs/tex/gertsenshtein_localized.tex](docs/tex/gertsenshtein_localized.tex)
and [examples/gertsenshtein/run_localized.sh](examples/gertsenshtein/run_localized.sh).
Cross-reference: Domcke, Garcia-Cely & Lee 2025 ([Literature/2507.16609](Literature/2507.16609))
extends the formula to full 3D scattering with Born approximation.

### Raffelt & Stodolsky 1988 — two-state mixing & Lorentzian resonance

**Citation:** G. Raffelt & L. Stodolsky, "Mixing of the photon with low-mass
particles", Phys. Rev. D **37**, 1237 (1988).

**Key result:** For a photon mixing with a second state carrying effective
mass `mA²`, the coherent conversion probability has a Lorentzian resonance
in the detuning `Δm² = mA² - m_g²_eff`:

```
P_max = (2·κ·B₀·k)² / [(2·κ·B₀·k)² + Δm²²]     [HWHM at Δm² = 2·κ·B₀·k]
```

where `m_g²_eff = κ²B₀²/2` is the graviton effective mass from background
EM stress-energy. On resonance (`Δm² = 0`) the mixing is maximal; off by
the HWHM `2·κ·B₀·k` the amplitude is suppressed by a factor of 2.

**Implemented in:** [examples/gertsenshtein_proca/theory.toml](examples/gertsenshtein_proca/theory.toml)
which adds a perturbation-level Proca mass `-(mA²/2) a·a` to Einstein-Maxwell.
Validated end-to-end by [sweep_resonance_1d.sh](examples/gertsenshtein_proca/sweep_resonance_1d.sh)
(observed HWHM = 0.405 vs theory 0.402 at B₀=0.10, commit `1526d77`) and
[sweep_resonance.sh](examples/gertsenshtein_proca/sweep_resonance.sh) (2D
resonance map, HWHM tracks theory within 9% across B₀ ∈ [0.09, 0.25]).

**Why perturbation-level Proca, not full-field:** effective photon mass in
plasma is a dispersion modification, not a Lorentz-invariant Lagrangian mass.
Placing the mass on the full field `A_μ A^μ` contaminates the graviton EOM
with spurious position-dependent `|Ā|²` terms (issue #142). By writing it on
`a[LI[1]]` only, xPert correctly injects the mass into the O(ε²) Lagrangian
without polluting the background. This is equivalent to Domcke, Garcia-Cely
& Lee 2025's EOM-level treatment `[☐ - μ²]A_h = -j_eff` via variational
calculus, and uses the curved-spacetime Maxwell conventions of Hwang & Noh
2023 ([Literature/2310.04150](Literature/2310.04150)).

## Torsion in Gauge Gravity

| Reference | Relevance |
| --- | --- |
| **Blagojević & Hehl (2013)** — *Gauge Theories of Gravitation*. Imperial College Press | Comprehensive PGT textbook |
| **Shapiro (2002)** — "Physical Aspects of the Space-Time Torsion", [arXiv:hep-th/0103093](https://arxiv.org/abs/hep-th/0103093) | Torsion phenomenology, propagation constraints, ghost-tachyon conditions |
| **Hehl et al. (1976)** — "General Relativity with Spin and Torsion", Rev. Mod. Phys. 48:393 | Foundational Einstein-Cartan review |
| **Sezgin & van Nieuwenhuizen (1980)** — "New ghost-free gravity Lagrangians", Phys. Rev. D 21:3269 | Linearized PGT propagating modes, ghost-free parameter conditions |
| **Nikiforova et al. (2009)** — "Stability of the Massive Torsion Modes", [arXiv:0905.4007](https://arxiv.org/abs/0905.4007) | Ghost-free windows, dispersion relations for torsion |
| **Barker (2024)** — "Every Poincaré gauge theory is conformal", [arXiv:2406.12826](https://arxiv.org/abs/2406.12826) | Particle spectra, no-ghost condition β₃ > 0 |
| **Obukhov & Trukhanova (2024)** — "Electrodynamics in Poincaré Gauge Theory", [arXiv:2410.01355](https://arxiv.org/abs/2410.01355) | Spin-torsion coupling to electrodynamics |
| **Blagojević & Cvetković (2018)** — "Ghost and tachyon free Poincaré gauge theories", [arXiv:1812.02675](https://arxiv.org/abs/1812.02675) | 450 ghost-and-tachyon-free PGT critical cases; systematic catalog |
| **Aoki (2020)** — "Non-linearly ghost-free higher curvature gravity", [arXiv:2009.11739](https://arxiv.org/abs/2009.11739) | Ghost-free PGT via equivalence to ghost-free massive bigravity |
| **Bahamonde et al. (2025)** — "Cosmology of Cubic Poincaré Gauge gravity", [arXiv:2506.17017](https://arxiv.org/abs/2506.17017) | Cubic torsion invariants eliminate ghost pathologies in axial/vector sectors |
| **Bahamonde et al. (2026)** — "Coupling Electromagnetism to Torsion", [arXiv:2507.02362](https://arxiv.org/abs/2507.02362) | Non-minimal torsion-EM coupling; black holes with spin-charge interactions |

## Non-Minimal Torsion-EM Coupling

| Reference | Relevance |
| --- | --- |
| **Hehl & Obukhov (2000)** — "How does the EM field couple to gravity?", [arXiv:gr-qc/0001010](https://arxiv.org/abs/gr-qc/0001010) | Canonical classification of all EM-gravity coupling structures in metric-affine geometry |
| **Rubilar et al. (2003)** — "Torsion nonminimally coupled to the EM field and birefringence", [arXiv:gr-qc/0305049](https://arxiv.org/abs/gr-qc/0305049) | Proves T²F² couplings make light sensitive to torsion; vacuum birefringence |
| **Itin (2003)** — "Maxwell's field coupled nonminimally to quadratic torsion", [arXiv:gr-qc/0307063](https://arxiv.org/abs/gr-qc/0307063) | Two families: F·T² (modifies Maxwell) and F²·T² (modifies constitutive tensor); induced axion field |

## Dark Photon / Kinetic Mixing Portal

| Reference | Relevance |
| --- | --- |
| **Holdom (1986)** — "Two U(1)'s and ε charge shifts", Phys. Lett. B 166, 196 | Original theorem: pure kinetic mixing between two U(1) gauge fields is equivalent to a field redefinition for matter-free sectors. Establishes that kinetic mixing alone is observable only via coupling to charged matter (millicharges). |
| **Pospelov (2008)** — "Secluded U(1) below the weak scale", [arXiv:0811.1030](https://arxiv.org/abs/0811.1030) | Canonical dark-photon-portal Lagrangian `L ⊃ -¼V² + ½κ V_{μν}F^Y_{μν} + ½m_V² V²` with positive-sign kinetic mixing convention. Effective coupling to EM current after diagonalisation. |
| **Redondo (2008)** — "Helioscope bounds on hidden sector photons", [arXiv:0801.1527](https://arxiv.org/abs/0801.1527) | Explicit `P(γ→γ') = 4ε² · sin²((m² − ω_p²)·L/4ω)` formula including plasma frequency; rule of thumb for in-medium vs vacuum conversion. |
| **An, Pospelov, Pradler (2013)** — "New stellar constraints on dark photons", [arXiv:1302.3884](https://arxiv.org/abs/1302.3884) | Stellar dark-photon bounds; eq. (2.4)–(2.6): vacuum and in-medium conversion formulas. |
| **Fabbrichesi, Gabrielli, Lanfranchi (2020)** — "The Dark Photon", [arXiv:2005.01515](https://arxiv.org/abs/2005.01515) | Comprehensive review of dark-photon phenomenology, conventions, and experimental constraints. |
| **Caputo, Millar, O'Hare, Vitagliano (2021)** — "Dark photon limits: a cookbook", [arXiv:2105.04565](https://arxiv.org/abs/2105.04565) | Practical oscillation formulae for dark-photon experiments including haloscope / DM-radio analyses. |

## Cosmology Program

The reference set for the TIDAL → Cosmology program (`docs/COSMOLOGY_PROGRAM.md`,
umbrella #488): spectator perturbations on a CAMB ΛCDM background, packaged as a Cobaya
extension. All TeX sources are local under `literature/` — read them there rather than
re-fetching.

### CMB theory canon

| Reference | Relevance |
| --- | --- |
| **Seljak & Zaldarriaga (1996)** — "A line-of-sight integration approach to CMB anisotropies", ApJ 469, 437, [arXiv:astro-ph/9603033](https://arxiv.org/abs/astro-ph/9603033) | The line-of-sight method every modern Boltzmann code uses: `C_ℓ` as sources × transfer functions instead of a full hierarchy. The route WS4 takes for per-rung observables, and the structure a birefringence rotation must be placed *inside* (see Murai below). |
| **Ma & Bertschinger (1995)** — "Cosmological perturbation theory in the synchronous and conformal Newtonian gauges", ApJ 455, 7, [arXiv:astro-ph/9506072](https://arxiv.org/abs/astro-ph/9506072) | The reference statement of the standard-sector perturbation equations in both gauges, and the dictionary between them. WS2 derives the new sector's equations against this convention. |
| **Moss (2026)** — "nanoCMB: a minimal CMB power spectrum calculator in Python", [arXiv:2602.23466](https://arxiv.org/abs/2602.23466) | A complete CMB solver in ~1400 readable lines — the onboarding reference for the pipeline as a whole, and an independent cross-check for the O0 gate (sub-percent vs CAMB/nanoCMB, `2 ≤ ℓ ≤ 2500`). |
| **Hahn et al. (2023)** — "DISCO-DJ I: a differentiable Einstein-Boltzmann solver for cosmology", [arXiv:2311.03291](https://arxiv.org/abs/2311.03291) | Einstein-Boltzmann in JAX with gradients throughout; the differentiable-host candidate weighed in the H3 integration-target decision. |
| **Blas, Lesgourgues & Tram (2011)** — "CLASS II: Approximation schemes", [arXiv:1104.2933](https://arxiv.org/abs/1104.2933) | Where the tight-coupling and free-streaming approximations come from, and what they buy (1069 s → 19.4 s). Sets WS3's performance budget: 10× slower than CAMB is acceptable, 100× is fatal. |

### Numerical methods for oscillatory ODEs (WS3)

WS3 replaces `expm(M·t)` with a solver for time-dependent `M(η)`. The ladder is
exponential midpoint → Magnus → adiabatic/WKB switching → transfer matrices; these are its
sources, and together they mark the research gap (**every published WKB-switching solver is
scalar-only; no matrix RKWKB solver exists**).

| Reference | Relevance |
| --- | --- |
| **Blanes, Casas, Oteo & Ros (2009)** — "The Magnus expansion and some of its applications", Phys. Rep. 470, 151, [arXiv:0810.5488](https://arxiv.org/abs/0810.5488) | The review for ladder rung 2: 4th-order Gauss-Legendre Magnus, batched over `k`, which reduces exactly to `expm(M t)` for constant `M` — giving a machine-precision regression test against the existing solver. |
| **Agocs, Handley, Lasenby & Hobson (2019)** — "An efficient method for solving highly oscillatory ODEs" (oscode / RKWKB), [arXiv:1906.01421](https://arxiv.org/abs/1906.01421) | Rung 3: adiabatic/WKB regime switching for `kη ≫ 1`, without which Magnus alone must resolve every oscillation (`∫‖M‖ds < π`, `‖M‖ ~ k`). Scalar-only — half of the WS3 research gap. |
| **Agocs & Barnett (2022)** — "An adaptive spectral method for oscillatory second-order linear ODEs with frequency-independent cost" (riccati), [arXiv:2212.06924](https://arxiv.org/abs/2212.06924) | The successor to oscode: cost independent of frequency. Same scalar-only limitation, so the same gap for a coupled system. |
| **Ioannisian & Smirnov (2008)** — "Describing neutrino oscillations in matter with Magnus expansion", [arXiv:0803.1967](https://arxiv.org/abs/0803.1967) | The closest existing analogue to what WS3 needs: Magnus/WKB applied to a *matrix* oscillation problem with a varying background. The working template for a matrix RKWKB scheme. |
| **Haddadin & Handley (2018)** — "Rapid numerical solutions for the Mukhanov-Sasaki equation", [arXiv:1809.11095](https://arxiv.org/abs/1809.11095) | Rung 4: piecewise-analytic transfer matrices for a cosmological mode equation — the fallback if WKB switching cannot be lifted to matrices. |

### Particle spectrum / polology (WS6)

| Reference | Relevance |
| --- | --- |
| **Barker et al. (2026)** — "Numerical polology: towards next-generation model-building for cosmology", [arXiv:2606.30785](https://arxiv.org/abs/2606.30785) | The WS6 algorithm source. States the case this program is built on: symbolic spectrum computation (PSALTer, [arXiv:2406.09500](https://arxiv.org/abs/2406.09500)) scales poorly through expression swell, so sampling *arbitrary* Lagrangians requires the numerical route. Minkowski-only is correct and sufficient — the spectrum screens in vacuo. |

### Cosmic birefringence and polarization (O4)

| Reference | Relevance |
| --- | --- |
| **Komatsu (2022)** — "New physics from the polarised light of the cosmic microwave background", Nat. Rev. Phys. 4, 452, [arXiv:2202.13919](https://arxiv.org/abs/2202.13919) | The review that frames O4: what a parity-odd sector does to `EB`/`TB`, how the measurement is actually made, and which systematics (miscalibrated polarization angles) it must be separated from. |
| **Das, Mohanty & Prasanna (2009)** — "Constraints on background torsion from birefringence of CMB polarization", [arXiv:0908.0629](https://arxiv.org/abs/0908.0629) | The only existing torsion → CMB-rotation paper, and the reason O4 is a *distinctive* channel: torsion-induced rotation scales as `ν²`, unlike the frequency-independent axion case — so the two are separable by the frequency data below. |
| **Minami & Komatsu (2020)** — "New extraction of the cosmic birefringence from the Planck 2018 polarization data", PRL 125, 221301, [arXiv:2011.11254](https://arxiv.org/abs/2011.11254) | The first `β ≠ 0` hint, and the method that made it possible: solving for the instrument's polarization-angle miscalibration simultaneously with `β`. |
| **Eskilt & Komatsu (2022)** — "Improved constraints on cosmic birefringence from the WMAP and Planck CMB polarization data", [arXiv:2205.13962](https://arxiv.org/abs/2205.13962) | Tests the frequency dependence across WMAP + Planck and finds `β` frequency-independent — which is exactly the observable that bounds the `ν²` torsion channel. |
| **Eskilt (2026)** — "Cosmic birefringence from a joint analysis of ACT and Planck", [arXiv:2608.06480](https://arxiv.org/abs/2608.06480) | The current measurement, `β = 0.277° ± 0.057°` (4.8σ) — O4's target number. |
| **Murai et al. (2022)** — "Isotropic cosmic birefringence from early dark energy", [arXiv:2209.07804](https://arxiv.org/abs/2209.07804) | Where the naive `EB` formula fails: when the rotation is time-dependent it cannot be applied as a post-processing rotation of the final spectra but must sit inside the line-of-sight integral. WS4's rule — post-processing is exact only for constant, isotropic, frequency-independent rotation. |
| **Cai & Guan (2021)** — "Computing microwave background polarization power spectra from cosmic birefringence" (class_rot), [arXiv:2111.14199](https://arxiv.org/abs/2111.14199) | A public implementation of rotated spectra; WS4's cross-validation target alongside CAMB and nanoCMB. |

### Spectator methodology

| Reference | Relevance |
| --- | --- |
| **Amendola, Ballesteros & Pettorino (2014)** — "Effects of modified gravity on B-mode polarization", [arXiv:1405.7004](https://arxiv.org/abs/1405.7004) | The O2 precedent: modified tensor propagation (friction, mass, speed) carried through the tensor transfer function to B-modes — the same chain O2 builds from a Lagrangian instead of a parametrization. |
| **Kushwaha & Jain (2025)** — "Constraining circular polarization of high-frequency gravitational waves with CMB", [arXiv:2502.12517](https://arxiv.org/abs/2502.12517) | The V-mode channel: graviton → photon conversion imprinting circular polarization. Supports the H7 argument that distinctive channels (V-modes, `E→B`) can beat the `N_eff` bound where broadband conversion cannot. |

### Cosmological magnetic fields (O3's assumed background)

O3 has no signal without an external coherent magnetic field — the mixing term is linear in
it. These are the sources behind `docs/cosmology/magnetic_field_background.md`, which
records what the expanding-universe Gertsenshtein literature assumes and what O3 adopts.

| Reference | Relevance |
| --- | --- |
| **Subramanian (2016)** — "The origin, evolution and signatures of primordial magnetic fields", Rep. Prog. Phys. 79, 076901, [arXiv:1504.02311](https://arxiv.org/abs/1504.02311) | The energy-budget statement that makes a primordial field a legitimate spectator: `r_B = B₀²/(8πρ_γ0) ≈ 10⁻⁷ B₋₉²`, equipartition with the CMB only at `B₀ ~ 3.2 μG`. The number O3 reports per run. |
| **Durrer & Neronov (2013)** — "Cosmological magnetic fields: their generation, evolution and observation", A&A Rev. 21, 62, [arXiv:1303.7121](https://arxiv.org/abs/1303.7121) | The conventions every conversion paper inherits: flux freezing giving `B ∝ a⁻²` and `ρ_B ∝ a⁻⁴`, the smoothing definition `B_λ² = 8π P_B(2π/λ)/λ³`, the coherence length, the 1 Mpc normalization, and the observational window. |
| **Domcke & Garcia-Cely (2021)** — "Potential of radio telescopes as high-frequency GW detectors", PRL 126, 021104, [arXiv:2006.01161](https://arxiv.org/abs/2006.01161) | The template for O3's line-of-sight rate: patch averaging with `Δz₀ = min[λ_EQ, λ_B⁰]`, the transverse `2/3`, the `𝓕` boost for a realistic small-scale spectrum, and the `ΔN_eff` inequality the program enforces. |
| **Addazi et al. (2024)** — "Resonant graviton-photon conversion with stochastic magnetic field in the expanding universe", [arXiv:2401.15965](https://arxiv.org/abs/2401.15965) | The stochastic-field treatment and the reason the domain model is a choice, not a default: conversion grows linearly in `λ_B` when `λ_B ≲ l_osc`, peaking at `λ_B = π l_osc`. |
| **Dolgov & Ejlli (2013)** — "Conversion of relic gravitational waves into photons in cosmological magnetic fields", [arXiv:1211.0500](https://arxiv.org/abs/1211.0500) | The one paper in this set that engages with a backreaction bound (Caprini–Durrer) — and explicitly sets it aside in favor of the weaker CMB bounds. The reason the program computes `r_B` rather than citing an assertion. |
| **Paoletti et al. (2022)** — "Constraints on primordial magnetic fields from their impact on the ionization history with Planck 2018", [arXiv:2204.06302](https://arxiv.org/abs/2204.06302) | The only source in the surveyed set that writes `ρ_B/ρ_γ` down explicitly, and the current CMB constraint methodology. |

The program's spectator-limit justification — the clean double expansion in which order 0
gives the background equations (discarded in favour of CAMB's solution), order 1 tadpoles
vanish on-shell, and order 2 retains the full quadratic action including all mixing terms —
is stated in Cembranos et al., [arXiv:2302.08186](https://arxiv.org/abs/2302.08186), listed
under Gertsenshtein Effect above and also local.

## Verification & Validation Methodology

| Reference                                                                                                                             | Relevance                                                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| **Roache (1998)** — "Verification and Validation in Computational Science and Engineering", Hermosa Publishers                        | V&V methodology: Richardson extrapolation, Grid Convergence Index (GCI), Method of Manufactured Solutions (MMS) |
| **NASA GRC** — "Examining Spatial (Grid) Convergence". [grc.nasa.gov](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html) | Standard grid convergence tutorial                                                                              |
| **AIAA G-077-1998** — "Guide for Verification and Validation of Computational Fluid Dynamics Simulations"                             | Industry standard for V&V                                                                                       |

## PML / Absorbing Boundaries

| Reference                                                                                                                            | Relevance                                       |
| ------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------- |
| **Bérenger (1994)** — "A perfectly matched layer for the absorption of electromagnetic waves", J. Computational Physics 114, 185–200 | Original PML formulation                        |
| **Johnson (2007)** — "Notes on Perfectly Matched Layers (PMLs)", MIT. [PDF](https://math.mit.edu/~stevenj/18.369/spring09/pml.pdf)   | Clear tutorial on complex coordinate stretching |

## Data Standards

| Reference                                                                                                                                        | Relevance                                                                  |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------- |
| **XDMF** — "XDMF Model and Format". [xdmf.org](https://www.xdmf.org/index.php/XDMF_Model_and_Format)                                             | Standard descriptor for HDF5-backed PDE output (ParaView/VisIt compatible) |
| **Schroeder et al. (2021)** — "Automated dimensional analysis for PDEs within FEniCS/UFL", [arXiv:2601.06535](https://arxiv.org/html/2601.06535) | Automated unit/dimension checking for variational forms                    |

## Topological Field Theory (Chern-Simons)

| Reference                                                                                                         | Relevance                                                          |
| ----------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| **Deser, Jackiw & Templeton (1982)** — "Topologically massive gauge theories", Annals of Physics 140, 372–411     | Original formulation of topologically massive gauge theory in 2+1D |
| **Dunne (1999)** — "Aspects of Chern-Simons Theory", [arXiv:hep-th/9902115](https://arxiv.org/abs/hep-th/9902115) | Pedagogical review of Chern-Simons field theory                    |

## Numerical Methods / ODE-DAE Solvers

| Reference                                                                                                                                                                                                            | Relevance                                                                                |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| **SUNDIALS** — Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and Differential/Algebraic Equation Solvers", ACM TOMS 31(3), 363–396 (2005). [LLNL](https://computing.llnl.gov/projects/sundials)                     | IDA (DAE) and CVODE (BDF adaptive ODE) backends via scikit-sundae                        |
| **Hairer, Lubich & Wanner (2006)** — _Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations_, Springer, 2nd ed.                                                       | Störmer-Verlet (leapfrog) symplectic integrator theory; shadow Hamiltonian error analysis |
| **Courant, Friedrichs & Lewy (1928)** — "Über die partiellen Differenzengleichungen der mathematischen Physik", Math. Ann. 100, 32–74                                                                                | CFL stability condition for explicit time-stepping                                       |
| **Dormand & Prince (1980)** — "A family of embedded Runge-Kutta formulae", J. Comp. Appl. Math. 6, 19–26                                                                                                            | DOP853 embedded RK method (available via scipy `solve_ivp`)                              |
| **Goldstein, Poole & Safko (2002)** — _Classical Mechanics_, 3rd ed., Addison-Wesley                                                                                                                                 | Legendre transform, canonical momenta, Hamiltonian formulation (Ch. 8)                   |

## Reproducibility

| Reference                                                                                       | Relevance                                                           |
| ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| **Snakemake** — Mölder et al., "Sustainable data analysis with Snakemake", F1000Research (2021) | Workflow management for parameter sweeps and reproducible pipelines |
