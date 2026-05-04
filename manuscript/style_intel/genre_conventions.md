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

## The BHL Physics-Calculation Paper

The five sections above cover the software-paper genre. §1, §2, §4, and §5 of the TIDAL report are primarily physics-theory content. Their genre conventions are distinct and are documented here, drawn from: 2406.12826, 2101.02645, 2510.08201, 2407.09598, 2303.11094, 2205.13534, and 2510.17094.

### (h) Macro-Structure of a BHL Physics Paper

A BHL physics paper at letter or short-article length (≤ 15 pages) has the following macro-structure:

1. **Abstract** (~180 words): phenomenon → contribution → scope statement or main finding.
2. **Introduction** (one section, ~600–800 words): four-stage scaffold (Frontier → Deep Problem → Prior Art → Hook). Ends with a structural paragraph mapping the remaining sections.
3. **Framework / Theory** (one or two sections): establishes the gauge structure, defines fields and notation, derives or states the key equations of motion or particle spectrum. Subsections are organised by the logical chain of the derivation, not by topic alone.
4. **Results / Analysis** (one or two sections): applies the framework. Opens with validation (known case reproduced), then advances to new results. Null results and constraints are reported before speculative findings.
5. **Conclusions / Discussion** (one section, ~400–600 words): three paragraphs — what we found, what it means for the broader programme, what should be done next. The final paragraph is always forward-looking but specific: not "future work is important" but "the next step is to test whether X holds at nonlinear order."
6. **Appendices**: proofs, extended tables, computational details that are primary data but would interrupt the main narrative.

The TIDAL report departs from this macro-structure only by inserting a "Computational approach" section (§3) between theory and results — otherwise it follows the template exactly.

### (i) How BHL Papers Introduce Gauge Theory Frameworks

The theory section of a BHL gauge-gravity paper follows a strict pedagogical sequence:
(a) Recall the Riemannian baseline (Levi-Civita connection, standard GR) in one or two sentences.
(b) State the generalisation: "In Poincaré gauge theory, the connection is promoted to an independent field." No more than two sentences.
(c) Define the two field strengths (curvature, torsion) in paired adjacent equations with nearly identical structure.
(d) State the most general quadratic action in the field strengths, with named coupling constants.
(e) Introduce the irreducible decomposition only after the action — never before.

Step (a) is present even in papers where the reader is assumed to know GR. Its purpose is to fix notation for the paper, not to teach GR. Step (e) is often in a subsection or named paragraph of its own ("Irreducible torsion", "Torsion decomposition").

### (j) How BHL Papers Handle a Parameter Space Result

When a calculation involves scanning a space of coupling constants or model parameters (as in §4 of the report):

- The **parameter space is always defined before any scan results are presented**. A table or list naming every parameter, its range, and its physical meaning appears before the first result figure.
- **Positive validation is always the first result**, even if the main scientific contribution is a null or a constraint.
- **Constraints are stated as closed inequalities**, not as upper/lower bounds from a scan: "ghost freedom requires $a_1 > 0$", not "our scan found no ghosts above $a_1 = 0$".
- **Figures carry the primary result load in long parameter scans**; prose identifies the key feature visible in each figure (the boundary of the stable region, the torsion-independent line) without repeating the number already readable in the plot.
- When the scan result is a null (no amplification found across the parameter space), the result is quantified as a precision bound: "conversion probability matches the torsion-free prediction to within $\delta < 10^{-5}$ across the entire $(a_1, a_3)$ plane". The prose then explains what physical mechanism enforces the null.

### (k) Discussion Section Conventions in BHL Physics Papers

The discussion section of a BHL physics paper (2303.11094, 2406.12826, 2510.08201) has a recognisable three-move structure:

1. **Restatement of findings** (~1 paragraph): a plain-language summary of what was found, without equations. Uses "we have shown that...", "we found that...", "our main finding is that...". Does not hedge retrospectively — findings that were reported as provisional in the results section are stated confidently in the conclusions.
2. **Implications for the programme** (~1–2 paragraphs): what the findings imply for the broader theoretical programme. This is where the "three constructive paths" framing belongs in the TIDAL report — not in the results section. The standard move is: "If our null result is robust, it has the following implication for X... however, three avenues remain open: [list]."
3. **Specific next steps** (~1 paragraph): one or two concrete, specific research questions that follow directly from the paper's findings. Not vague ("future work should explore...") but actionable: "The natural next step is to apply TIDAL to the full non-minimal PGT Lagrangian (App. B), which contains the torsion–curvature cross-term that the present campaign excluded by construction." This closing paragraph should name a specific theory, a specific observable, or a specific calculational extension.
