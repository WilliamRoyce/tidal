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
| **Domcke & Garcia-Cely (2023)** — "A simple derivation of the Gertsenshtein effect", [arXiv:2301.02072](https://arxiv.org/abs/2301.02072)                                                                            | Modern pedagogical derivation; thin-magnet formula P ~ (κ B L / 2)²                                           |
| **Hwang & Noh (2023)** — "On graviton-photon conversions in magnetic environments", [arXiv:2310.04150](https://arxiv.org/abs/2310.04150)                                                                              | Proper EM field definitions in curved spacetime; formulation-dependent graviton mass; critique of ad-hoc derivations |
| **Berlin et al. (2024)** — "Numerical analysis of resonant axion-photon mixing", [arXiv:2405.08865](https://arxiv.org/abs/2405.08865)                                                                                | Numerical methods for resonant graviton-photon / axion-photon mixing — directly analogous to TIDAL's use case |
| **Ejlli et al. (2019)** — "Upper limits on the amplitude of ultra-high-frequency gravitational waves", Eur. Phys. J. C 79, 1032 (2019). [Springer](https://link.springer.com/article/10.1140/epjc/s10052-019-7542-5) | Experimental bounds on Gertsenshtein conversion                                                               |
| **Ejlli (2020)** — "Graviton-photon mixing: exact solution in a constant magnetic field", [arXiv:2004.02714](https://arxiv.org/abs/2004.02714) | First exact (non-perturbative) Gertsenshtein solution |
| **Palessandro (2024)** — "Graviton-Photon Oscillations as a Probe of Quantum Gravity", [arXiv:2405.01407](https://arxiv.org/abs/2405.01407) | Extended Gertsenshtein treatment (**contains normalisation error** corrected by TIDAL) |
| **Hwang & Noh (2024)** — "Graviton-photon conversions in Euler-Heisenberg nonlinear electrodynamics", [arXiv:2405.11786](https://arxiv.org/abs/2405.11786) | Nonlinear QED corrections to Gertsenshtein; chiral GW propagation |
| **Graviton-photon oscillation in general modified gravity (2023)** — [arXiv:2302.08186](https://arxiv.org/abs/2302.08186) | Gertsenshtein in Horndeski, massive gravity, cosmic backgrounds |

## Torsion in Gauge Gravity

| Reference | Relevance |
| --- | --- |
| **Blagojević & Hehl (2013)** — *Gauge Theories of Gravitation*. Imperial College Press | Comprehensive PGT textbook |
| **Shapiro (2002)** — "Physical Aspects of the Space-Time Torsion", [arXiv:hep-th/0103093](https://arxiv.org/abs/hep-th/0103093) | Torsion phenomenology, propagation constraints, ghost-tachyon conditions |
| **Hehl et al. (1976)** — "General Relativity with Spin and Torsion", Rev. Mod. Phys. 48:393 | Foundational Einstein-Cartan review |
| **Sezgin & van Nieuwenhuizen (1980)** — "New ghost-free gravity Lagrangians", Phys. Rev. D 21:3269 | Linearised PGT propagating modes, ghost-free parameter conditions |
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
