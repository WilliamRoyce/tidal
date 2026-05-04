# Genre Conventions: The Physics Software Package Paper

This document synthesises the conventions of the physics-software-package paper genre, drawing primarily on three Barker papers: PSALTer original (2406.09500, Phys. Rev. D), PSALTer v2 (2506.02111, targeting a journal format with a JHEP-style layout), and Hamilcar (2512.25007, JHEP). The BHL HiGGS paper (2206.00658) provides a predecessor that illuminates how the genre evolved. Differences between the three papers are noted where they are instructive.

---

## (a) Title Conventions

All three Barker software papers place the software name first or make it the dominant element: "PSALTer: Particle Spectrum for Any Tensor Lagrangian", "The particle spectra of parity-violating theories: a less radical approach and an upgrade of PSALTer", "Fast Poisson brackets and constraint algebras in canonical gravity". The first title gives the tool an acronym that is immediately expanded; the second leads with the physics problem and identifies the tool as an upgrade; the third describes the capability without naming the tool in the title (Hamilcar appears only in the abstract and body).

The pattern: software papers benefit from a title that either embeds the tool name or encodes the physical capability the tool unlocks. Pure acronym-only titles ("We present PSALTer") are avoided; the title always gestures at the physics. A subtitle or colon structure is common.

HiGGS (2206.00658) uses "Supercomputers against strong coupling in gravity with curvature and torsion" — here the computational capability (supercomputers) comes first, and the physics problem comes second. This is the most evocative title of the set, and it positions the paper as a computational physics paper as much as a mathematical physics paper.

---

## (b) Problem-Centric Opening: Why Does Physics Need This Tool?

Every Barker software paper opens by arguing for necessity, not by announcing capability. The first section of text establishes that a particular physical calculation is important, that it is currently expensive or infeasible without the tool, and that prior attempts either scaled poorly or were not made available to the community. The tool is introduced as the solution to a previously articulated problem, not as an achievement in its own right.

PSALTer (2406.09500) opens with "In theoretical physics it is often necessary to extract the propagating degrees of freedom from a specified model", moves through the difficulty for higher-rank tensors, and reaches the tool's introduction only after establishing that "no computer algebra implementation of the algorithm was made available to the community." The tool fills an explicit gap.

Hamilcar (2512.25007) opens with the Dirac–Bergmann algorithm's importance for constraint counting, then explains that "certain steps in canonical analysis would benefit from automation", principally the computation of Poisson brackets "between arbitrary functionals of phase-space variables." The tool automates what was previously done by hand or not done at all.

This problem-first framing is not universal in software papers from other groups; it is a Barker stylistic choice that makes the papers read more like physics papers than user manuals. The practical effect is that the paper can be understood — and cited — by readers who will never use the software, because the physics motivation stands independently.

---

## (c) Architecture and Capability Section Structure

The main technical section of a Barker software paper is organised by *what the software does*, not by *how it was implemented*. The PSALTer paper has a section titled "Examples with code" and a separate section "Theoretical development"; Hamilcar has "Explicit examples" covering both usage and theory. In each case, the user-facing API comes before the internal theoretical machinery.

Within the architecture section, the sub-structure is:
1. Loading and installation (brief — a named paragraph, not a subsection, unless the paper is a formal manual)
2. Geometric environment / pre-defined objects (a table or enumeration, always present)
3. Key function definitions (API reference style — command syntax, arguments, options, types)
4. Syntax highlighting conventions (for papers with code listings)

The PSALTer original (2406.09500) is more tutorial-like: it walks through examples linearly from scalar to vector to tensor theories, with the code listing for each. Hamilcar is more reference-like: it separates the API description from the worked examples, giving a complete account of each function before showing it in action. The PSALTer v2 (2506.02111) sits between the two — it presents the upgraded algorithm as a theoretical development first, then demonstrates it with worked examples including new parity-violating cases.

**Evolution across versions:** PSALTer v1 is primarily a demonstration paper — the theory section is secondary to the examples. Hamilcar is more evenly weighted toward theory and examples, reflecting Barker's growing confidence that the community needs both. PSALTer v2 explicitly labels an algorithmic advance as a contribution distinct from the software, showing that genre maturity involves articulating conceptual contributions separately from implementation.

---

## (d) Validation Section Conventions

Validation in Barker software papers proceeds through a hierarchy of tests:

1. **Analytic limits**: The simplest possible case (scalar field, massless Klein–Gordon, Maxwell electromagnetism) is run first, and the output is compared to the textbook result. The goal is to demonstrate that the tool reproduces what is already known.

2. **Reference model reproduction**: A well-known theory with a non-trivial spectrum (Proca, linearised GR, Einstein–Cartan theory) is run, and the output is compared to published results from prior literature. Citations to the specific prior computations are included.

3. **Pathological cases**: A "sick" theory — one that should fail unitarity — is run, and the tool is shown to correctly identify the ghost or tachyon. This demonstrates that the tool does not silently pass bad theories.

4. **Scale-up / non-trivial case**: A theory that would be impractical to compute manually (higher-spin, higher-rank, or full PGT) is run, demonstrating the tool's core advantage. Timing or computational cost data may be included.

The PSALTer paper (2406.09500) progresses through exactly this hierarchy: massless scalar, massive scalar, Maxwell, Proca, sick Maxwell, sick Proca, then higher-rank theories including the full PGT spin connection. Hamilcar (2512.25007) validates against pure GR (the canonical example), R² theory (a non-trivial extension), and then the Goroff–Sagnotti two-loop GR action (the "scale-up" case).

**What is absent from validation sections:** convergence tests in the numerical sense, unit tests in the software engineering sense, or timing benchmarks with formal error bars. The validation is physics-completeness testing, not software QA in a computer-science framework.

---

## (e) Example / Case Study Section Conventions

Examples in Barker software papers serve a dual purpose: they validate the tool, and they are themselves physics results that may be novel. The PSALTer paper derives the full particle spectrum of PGT theories as part of its worked examples — the examples section is simultaneously a user tutorial and a research contribution.

Each example follows the pattern:
- Name the theory (bold paragraph heading or subsection)
- State the action / Lagrangian in display form
- Show the code for inputting the Lagrangian (for software papers)
- Show the output (always as a figure for PSALTer; as an equation or formula for Hamilcar)
- Interpret the output in prose, citing prior literature for comparison

The convention is to begin with the simple and progress toward the complex, where "complexity" means both physical content (more degrees of freedom, more sectors) and computational demand (more coupled, less tractable analytically). The reader is expected to follow the simple cases completely, and to trust the complex cases on the basis of the validated simple ones.

In PSALTer v2 (2506.02111), the examples section is divided into two explicit parts: a pedagogical set of toy models (p-forms) followed by the main physics result (most general parity-indefinite Einstein–Cartan theory). This bipartite structure — pedagogy then research — is cleaner than the original PSALTer's more uniform progression, and reflects genre maturation.

---

## (f) Handling Implementation Detail

Implementation detail in Barker software papers is consistently relegated to appendices and footnotes. The main text mentions algorithmic choices at a high level ("the Moore–Penrose pseudoinverse", "parallelised across available CPU cores") but does not expose loop structures, data types, or memory management. This is a physics-community convention: the audience trusts the implementation is correct once validation passes.

Where code is shown in the main text, it serves as user documentation, not implementation documentation. The PSALTer paper shows user-session code (how to call `ParticleSpectrum[]`) but not the internal routines that compute the SPOs. Hamilcar shows how to call `PoissonBracket[]` and `FindAlgebra[]` but not the variational derivative engine underneath.

Long code listings are handled via `\lstinputlisting{}` from external `.tex` files (one file per listing). This keeps the main source clean and allows listings to be updated without touching the manuscript. The Hamilcar paper uses this convention consistently throughout.

Appendices carry: installation instructions (always present), mathematical proofs of the underlying algorithm, and index conventions. The HiGGS/supercomputers paper (2206.00658) has the most extensive appendix structure, with the full PGT constraint algebra tabulated there — recognising that the constraint brackets are primary data that the paper makes available but cannot display inline.

---

## (g) First Release Scope

Each Barker software paper explicitly states what is and is not included in the current release, and what is deferred to future work. This scoping is not apologetic; it is presented as a deliberate choice rather than a limitation. The PSALTer paper closes its abstract with "The initial release allows for parity-preserving operators constructed from fields of up to rank three: this functionality will be extended in future versions." The conclusions section restates this and gives a research roadmap for extensions (parity-odd operators, complex fields, automated identification of special cases).

The convention is: state the scope positively as "what we achieve" first, then state the extension as "what we will do next". Do not lead with "we could not include X" — lead with "we include A, B, C, which covers the most important physical cases; D and E are left for future work."

In Hamilcar (2512.25007), the scope statement is more precisely calibrated: the paper explicitly says that Hamilcar "completely deprecates" HiGGS, the predecessor package, and explains why the new package is strictly more general. This successor-declares-predecessor-obsolete move is unusual but honest — it tells the community where to direct their efforts.

---

## Cross-Paper Evolution Notes

**From PSALTer v1 to v2:** The v1 paper is largely a tool demonstration; the physics results (particle spectra) are illustrative rather than novel per se. By v2, Barker explicitly separates a new algorithmic contribution (avoiding radicals when computing ghost conditions) from the software upgrade. The papers are becoming more confident about articulating intellectual contributions alongside tool contributions.

**From HiGGS to Hamilcar:** HiGGS (2206.00658) is presented as an HPC survey instrument — the emphasis is on scale. Hamilcar is presented as a general-purpose algebraic tool. HiGGS lacks a formal API section; Hamilcar has a complete one. The genre understanding has matured: software papers now need an API reference section, not just an existence proof.

**Persistent genre signature across all four papers:** Every paper in this group ends its examples section with a non-trivial case that would have been impractical without the tool. This "flagship example" structure closes the loop between the motivating problem (stated in the introduction) and the delivered capability (demonstrated in the examples). The flagship example is always a theory from the literature, not a synthetic benchmark.

---

## Appendix Conventions for TIDAL Software Documentation

The five sections above cover the software-paper genre at the level of a full standalone paper (PSALTer, Hamilcar). TIDAL's software contribution appears in Appendices A–E of a physics report rather than as a standalone paper. The conventions differ from a standalone software paper in important ways. This section documents how to write Apps A–E.

**The general rule**: each appendix is a focused technical document that supports a specific claim in the main body. The main text makes a claim ("TIDAL auto-selects the Fourier modal solver") and cites the appendix; the appendix documents how and why. Each appendix opens with a one-sentence purpose statement.

### (h) App A — Pipeline Architecture (~1800 words + diagram)

The audience for App A is a reader who wants to understand the overall system before examining any subsystem. Model: the HiGGS paper (2206.00658) appendix, which documents the full PGT constraint algebra and serves as primary data that the paper makes available.

**Structure**:
1. One-paragraph overview (3–4 sentences): what TIDAL does, what its inputs and outputs are, what it automates. Write this as the appendix's opening — a précis of the whole pipeline.
2. Pipeline diagram (TikZ): the symbolic stage (Wolfram/xAct) → JSON handoff → numerical stage (Python). Keep the diagram schematic; do not show code-level details.
3. One named paragraph per pipeline stage: (i) Symbolic (xAct/xPert, Euler–Lagrange, component decomposition); (ii) JSON handoff (format, what fields are present); (iii) Solver dispatch (auto-selection logic); (iv) Output (simulation data, measurement types). Each named paragraph is 4–6 sentences.
4. One named paragraph on dependencies: cite each package with a reference (xAct, xPert, xTras, SUNDIALS/scikit-sundae, numpy, scipy, Mathematica). "Standard tools" without references is not acceptable.

**What to avoid**: describing the algorithm's internal structure in App A — that is App C. App A documents the architecture from the *user's perspective*: what goes in, what the tool does at each stage, what comes out.

### (i) App B — Symbolic Stage (~1600 words)

The audience for App B is a reader who wants to understand the Wolfram/xAct component: how it receives a TOML Lagrangian, what symbolic operations it performs, what it outputs.

**Structure**:
1. One paragraph on the xAct environment: what packages are loaded, why xPert is used for metric perturbations (not xTens alone), and why xTras is needed for the Einstein–Cartan coupling.
2. One named paragraph per major symbolic step: (i) Lagrangian input from TOML; (ii) Euler–Lagrange derivation (E-L velocity form, not Legendre transform); (iii) component decomposition (DecomposeToComponents, cross-field handling); (iv) JSON export.
3. For each step, show representative Wolfram code *only* if it serves as user documentation (how to call the function, not how it is implemented internally). The PSALTer paper convention: show user-session code, not internal routines.
4. Close with a paragraph on failure modes and diagnostics: what happens when the derivation times out, when indices are inconsistent, when the JSON is malformed.

**Code listing conventions**: use `\lstinputlisting{}` from external `.wls` files rather than inline listings, following the Hamilcar convention. This keeps the appendix source clean and allows listings to be updated independently.

### (j) App C — Numerical Stage (~2000 words)

The audience for App C is a reader who wants to understand the Python solver backends and the selection logic. The Fourier modal solver is the primary backend and receives the longest treatment; the others are documented briefly.

**Structure**:
1. One opening paragraph: the solver selection hierarchy (modal → IDA → CVODE → leapfrog → scipy) and the auto-selection logic.
2. **Fourier modal solver** (the primary backend, ~800 words): describe the eigendecomposition approach (constant-coefficient systems → exact matrix exponential; position-dependent → Krylov expm_multiply). State what "machine precision" means in this context. Describe the constraint elimination via Fourier Schur complement. Cite the mathematical foundations (matrix exponential methods: Higham 2008, Moler–van Loan). One figure showing convergence vs. resolution.
3. **Analytical Jacobian** (~400 words): describe the three-tier structure (dense/sparse/GMRES), with one sentence on when each is active. State the speedup factors.
4. **Other backends** (~200 words total): IDA for DAE systems, CVODE for adaptive ODE, leapfrog for symplectic, scipy as general-purpose fallback. One named paragraph of 3–4 sentences each.
5. Close with a paragraph on the JSON-to-solver interface: how the coefficient matrix is extracted from the JSON spec and how solver selection is dispatched.

### (k) App D — Validation Suite (~1200 words)

Follow the Barker validation hierarchy (from PSALTer 2406.09500 and Hamilcar 2512.25007):

1. **Analytic limits** (named paragraph): massless Klein–Gordon, Maxwell electromagnetism. TIDAL output compared to textbook result. One figure per case.
2. **Reference model reproduction** (named paragraph): the Boccaletti formula for Einstein–Maxwell Gertsenshtein conversion. Show the simulation vs. analytic comparison with the 0.04% agreement stated explicitly. This is the most important validation result.
3. **Pathological cases** (named paragraph): the trace-channel ghost, identified in §4.5. App D verifies that TIDAL correctly diagnoses ghost instability (exponential growth in the time-domain simulation) at the known unstable parameter point. This validates the ghost-diagnosis methodology.
4. **Convergence** (named paragraph): convergence test for the Fourier modal solver. Show error vs. grid resolution for one theory (Einstein–Maxwell). State the convergence rate (spectral/machine-precision for the modal solver).

**Caption convention for App D figures**: each caption must include (i) what the figure shows, (ii) the parameter values used, and (iii) the key numerical result or agreement statement. For Gertsenshtein figures, include the $B_0$ value and $t_\mathrm{end}$ explicitly.

### (l) App E — HPC Infrastructure and Reproducibility (~1000 words)

The audience for App E is a reader who wants to reproduce the campaign results or understand the computational resources required.

**Structure**:
1. One paragraph on the CSD3 setup: partition (sapphire/icelake), node configuration, TIDAL version, Python/scipy/sundials versions. Cite CSD3 (DiRAC; the specific HPC facility).
2. One table: campaign summary. Columns: campaign label, HPC job ID, theory, parameters swept, number of runs, wall time, key result. This table is the primary reproducibility artefact.
3. One paragraph on reproducibility: TIDAL version at campaign time (commit hash or version tag), where to find the input JSONs and output data. Reference the Zenodo DOI (if minted) or the GitHub repository tag.
4. One paragraph on data archiving: where simulation outputs are stored, file format (HDF5, CSV, JSON), and how to re-run a single campaign point with the provided inputs.

**Reproducibility statement template** (from Hamilcar): "All results were produced with [Software] version [X.Y.Z] at commit [hash], available at [URL]. Input files and output data are archived at [DOI/URL]."
