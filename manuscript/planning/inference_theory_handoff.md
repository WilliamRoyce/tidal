# Inference-theory handoff for the manuscript appendix

## 1. Orientation

You are writing manuscript prose for the appendix that covers nested sampling
and Bayesian inference theory as used in this project. The target file is
either:

- additional preliminary paragraphs prepended to
  `manuscript/sections/appendices/inference_architecture.tex` (App J), before
  the current v3-specific paragraphs, or
- a sibling appendix dedicated to the theory, with App J then citing it.

Pick (a) if the theory paragraphs fit inside App J's word budget without
crowding the v3 content; pick (b) if length pushes past ~1500 words.

**Style anchors** (see §6):

- `\paragraph*{Title}` headings, not `\subsection`.
- Semantic CamelCase labels (`\label{EvidenceIntegral}`), never `\label{eq:1}`.
- Prose before and after every equation.
- First-person plural ("we") for reported work.
- Cite on first mention of any external tool/paper.
- Cross-references via `\cref{}` with descriptive labels.

**Audience**: a physicist who knows Bayesian inference but may not be familiar
with our non-standard repurposing of it.

**The critical framing the prose MUST establish up front**: there is no
observed dataset. $L(\theta)$ is a deterministic physics score computed by the
PDE solver — the amplification ratio of graviton–photon conversion relative to
the General-Relativistic (Gertsenshtein) baseline — not $P(d_{\mathrm{obs}}\mid
\theta)$. Without making this explicit early, every subsequent sentence about
$\mathcal{Z}$ and $D_{\mathrm{KL}}$ misleads the reader into reading it as
ordinary parameter estimation against data. It is not.

## 2. Existing infrastructure to reuse

- `docs/tex/inference.tex` — substantial (~530-line) project-level documentation
  of the inference machinery: CLI, module layout, prior types (uniform,
  log-uniform, normal, arctan-uniform), sampling backends (PolyChord +
  anesthetic), output dataclasses. **This is implementation reference, not
  appendix prose** — paraphrase sparingly and do not copy. The appendix
  should focus on theory and physical interpretation, not CLI surface.
- `manuscript/sections/appendices/inference_architecture.tex` — current App J
  scaffold with six TODO paragraphs (Motivation, Soft-penalty floor,
  Compactified priors, Cubed-sphere joint prior, Architecture validation,
  Linearised-theory validity). These cover the v3-specific architecture; the
  theory paragraphs you write are upstream of them.
- `tidal/inference/_nested.py` — PolyChord backend, `recommend_nlive` heuristic,
  result loading via anesthetic.
- `tidal/inference/_likelihood.py` — defines the score function modes
  (`maximize`, `minimize`, `extremize`, `gaussian`, `threshold`), the validity
  cap at $P_{\max} > 0.5$, and the score expression
  $\log L = \pm\log(P_{\max}/\mathrm{baseline})$.
- `tidal/inference/_prior.py` — implements the prior transforms from $[0,1]^d$
  to physical parameter space.
- `tidal/inference/_importance.py` — implements $D_{\mathrm{KL}}$, $d_G$, and
  per-parameter marginal $D_{\mathrm{KL}}$ via anesthetic.

## 3. Theory content the appendix must cover

Each subsection below is a self-contained brief for one manuscript paragraph
(or merge nearby ones if the structure benefits). Each brief gives: the key
idea, the equation(s) to typeset, prose pointers, and the citation.

### 3.1 Bayesian setup and the evidence integral

**Key idea**: introduce the quantity nested sampling targets. The Bayesian
evidence is

$$
\mathcal{Z} = \int L(\theta)\,\pi(\theta)\,\mathrm{d}\theta.
$$

In ordinary Bayesian inference $L(\theta) = P(d_{\mathrm{obs}}\mid\theta)$ and
$\mathcal{Z} = P(d_{\mathrm{obs}}\mid M)$ — the probability of the observed
data under model $M$. Sketch the curse of dimensionality (the prior volume
grows exponentially, the high-likelihood region shrinks to a tiny corner) to
motivate why naive Monte Carlo fails.

**Citation**: `skilling2004nested` (Skilling 2004, AIP Conf. Proc. 735) — see
§4 for bibtex.

### 3.2 Skilling's reformulation: prior volume contraction

**Key idea**: reformulate the $d$-dimensional integral as a 1D integral over
prior volume. Define the prior mass enclosed by the likelihood contour at
threshold $L^*$:

$$
X(L^*) = \int_{L(\theta) > L^*}\pi(\theta)\,\mathrm{d}\theta.
$$

Then

$$
\mathcal{Z} = \int_0^1 L(X)\,\mathrm{d}X.
$$

State the live-points algorithm in two or three sentences: $N_{\mathrm{live}}$
points are drawn from the prior; at each step the point with the lowest
$L$-value is replaced by a new draw from the prior conditioned on
$L > L_{\min}$; the prior volume contracts geometrically, with the surviving
volume after $i$ iterations $X_i \approx \exp(-i/N_{\mathrm{live}})$.

Note the termination criterion succinctly: iteration stops when the remaining
live points contribute negligibly to the accumulated $\mathcal{Z}$, i.e. when
$L_{\max}\,X_{\mathrm{rem}} \ll \mathcal{Z}_{\mathrm{accum}}$. The surviving
live points are the highest-scoring points found; they contribute only as a
lump-sum correction at termination because the surviving volume is
exponentially small.

**Citation**: `skilling2004nested`.

### 3.3 PolyChord and slice sampling

**Key idea**: the constrained-prior draw step is the bottleneck. Ellipsoidal
rejection (MultiNest-style) scales poorly in high $d$. PolyChord replaces it
with slice sampling along $\mathcal{O}(d)$ random directions per replacement.
The cost scales as $\mathcal{O}(d^2)$, enabling tractable inference in our
4–8 dimensional coupling spaces.

State two implementation facts the reader needs to interpret our settings:

- $\mathrm{num\_repeats} = 5d$ default — sets the slice-sampling chain length
  per replacement.
- PolyChord operates internally on the unit hypercube $[0,1]^d$. The
  `prior_transform` map carries the shape of $\pi(\theta)$ — uniform,
  log-uniform, arctan-uniform — out into physical coupling space.

**Citation**: `handley2015polychord` (Handley, Hobson & Lasenby 2015, MNRAS 453).

### 3.4 The non-standard repurposing: simulation-based inference

**This is the most important paragraph.** It establishes the framing that
makes the rest of the appendix interpretable.

**Key idea**: in this work there is no observed gravitational-wave or
photon-flux measurement to fit. The PDE solver maps each coupling vector
$\theta$ to a deterministic conversion probability $P_{\max}(\theta)$. We
construct a score function from $P_{\max}$ and use it in place of a
statistical likelihood.

Define the **amplification factor** relative to the Gertsenshtein baseline:

$$
A(\theta) \;=\; \frac{P_{\max}(\theta)}{\sin^2(\kappa B_0 t/2)}.
$$

Define the score function used as the PolyChord likelihood:

$$
\log L(\theta) \;=\; \pm\log A(\theta).
$$

The sign selects between three implemented modes (see
`tidal/inference/_likelihood.py`):

- `maximize` ($+\log A$): rewards $\theta$ producing larger conversion than GR.
- `minimize` ($-\log A$): rewards $\theta$ producing smaller conversion than GR.
- `gaussian:T:σ` ($-\tfrac{1}{2}\{(P_{\max}-T)/\sigma\}^2$): used when a target
  $P_{\max}$ is known, e.g. fitting to a measured signal. Not currently the
  campaign default.

The prose must NOT slide into "fitting data" or "constraining parameters from
observations" language. The right language is *surveying* or *exploring* the
coupling space: which combinations of couplings produce amplified
graviton–photon conversion, and which suppress it.

The closest standard name for this mode of analysis is **prior-predictive
analysis** or **simulation-based inference** (the latter from the SBI
literature). The appendix may use either term; "simulation-based inference"
is more familiar to the gravitational-physics community.

### 3.5 Physical interpretation of $\mathcal{Z}$

**Key idea**: with $L(\theta) = A(\theta)$ (taking the `maximize` sign for
concreteness),

$$
\mathcal{Z} \;=\; \int A(\theta)\,\pi(\theta)\,\mathrm{d}\theta \;=\; \langle
A\rangle_\pi.
$$

$\mathcal{Z}$ is therefore the **prior-averaged amplification** — the expected
amplification factor if one drew a random coupling vector from the prior and
ran the PDE simulation.

Three regimes for the reader:

- $\log\mathcal{Z}\approx 0$ ($\mathcal{Z}\approx 1$): the model is GR-like on
  average across the coupling space. *Null* result.
- $\log\mathcal{Z} \gg 0$: enhancement is **generic** — the typical sample
  from the prior already produces amplification.
- $\log\mathcal{Z} \ll 0$: the model on average suppresses conversion below
  the GR baseline.

The crucial caveat to state explicitly: this is **not** a Bayes factor in the
conventional sense. We are not comparing two models against the same observed
dataset. The number is a prior-predictive summary of the amplification
distribution, and its scale is set by the chosen prior over couplings. A wide
prior dilutes the average; a narrow prior concentrated on amplifying
parameters inflates it. This is desirable — the prior encodes what we
consider physically reasonable couplings, and $\mathcal{Z}$ asks whether
amplification occurs within that space.

### 3.6 Physical interpretation of $D_{\mathrm{KL}}$

**Key idea**: $\mathcal{Z}$ alone does not say whether amplification (if
present) comes from specific coupling values or is broadly distributed. The
Kullback–Leibler divergence from prior to score-weighted posterior fills that
gap:

$$
D_{\mathrm{KL}}\bigl(p(\theta\mid A)\,\|\,\pi(\theta)\bigr) \;=\; \int
p(\theta\mid A)\,\log\!\left[\frac{p(\theta\mid A)}{\pi(\theta)}\right]
\mathrm{d}\theta,
$$

where the "posterior" is the prior reweighted by the score:

$$
p(\theta\mid A) \;\propto\; A(\theta)\,\pi(\theta).
$$

Two regimes:

- $D_{\mathrm{KL}}\approx 0$: the posterior is indistinguishable from the
  prior — amplification (or its absence) is uniform across the parameter
  space. No coupling region drives the result.
- $D_{\mathrm{KL}}$ large: the posterior is concentrated on a specific
  coupling region. The result is *fine-tuned*.

State the diagnostic value clearly: $D_{\mathrm{KL}}$ is the quantity that
discriminates fine-tuned signals from generic behaviour. A high
$\log\mathcal{Z}$ paired with small $D_{\mathrm{KL}}$ means amplification
happens everywhere in the prior; a high $\log\mathcal{Z}$ paired with large
$D_{\mathrm{KL}}$ means amplification requires specific parameter values.

The implementation also provides per-parameter marginal $D_{\mathrm{KL}}$
(`tidal/inference/_importance.py`): the information gain projected onto a
single coupling at a time. This decomposition identifies which individual
couplings carry the signal.

**Citation**: `handley2019anesthetic` (Handley 2019, JOSS 4(37)) — the tool
that computes $D_{\mathrm{KL}}$ from the PolyChord output.

### 3.7 Bayesian model dimensionality $d_G$

**Key idea**: a second derived quantity quantifies *how many couplings* are
effectively constrained:

$$
d_G \;=\; 2\,\mathrm{Var}\bigl[\log L \mid p(\theta\mid A)\bigr].
$$

If $d_G \approx 0$, no coupling is constrained by the score (consistent with
small $D_{\mathrm{KL}}$). If $d_G$ approaches or exceeds the dimensionality
of the coupling vector, multiple couplings are jointly constrained.

In the appendix the simplest description is: $d_G$ counts the effective
number of couplings on which the score function depends within the prior.
Small $d_G$ paired with small $D_{\mathrm{KL}}$ is the appendix's compact
signature of a "no-constraint" outcome — language that §4 Results will appeal
back to.

### 3.8 Validity regime and gates

**Key idea**: the score function $A(\theta) = P_{\max}/\sin^2(\kappa B_0 t/2)$
is only meaningful in the linearised, perturbative regime in which our PDE
solver operates. State the two conditions explicitly:

- $P_{\max} \ll 1$ (perturbative conversion probability),
- $\kappa B_0 t \ll 1$ (small-argument limit of the GR baseline; equivalently
  $P_{\mathrm{GR}} \ll 1$).

Outside this regime, $A$ acquires a spurious $B_0$ dependence and the
linearised conversion measurement breaks down.

Two hard physics gates are implemented in
`tidal/inference/_likelihood.py`:

- **Tachyonic gate**: parameter regions with positive-real-part eigenvalues
  (growing solutions) are excluded. In v2 this was a hard $-\infty$
  rejection at $\gamma_{\mathrm{eff}} > 0.15$; in v3 it is replaced by the
  gradient-aware soft-penalty floor described in App J §"Soft-penalty floor"
  (forward-reference `\cref{SoftPenaltyFloor}`).
- **Hwang–Noh perturbativity gate**: samples with $P_{\max} > 0.5$ are
  excluded as having crossed the perturbative validity boundary. This is
  always applied, independent of v2/v3 architecture.

Forward-reference `\cref{PerturbativeRegime}` (App G) for the full physical
discussion of why these boundaries matter.

**Citation**: Hwang & Noh (Phys. Rev. D 65, 124010, 2002) for the
perturbativity boundary — see §4 for bibtex caveat.

## 4. Required bibliography additions

Before citing, you must add these entries to `manuscript/references.bib`. Use
INSPIRE-HEP for bibtex per the project standard
(`feedback_inspire_for_bibtex.md`): fetch from
`https://inspirehep.net/api/literature?q=arxiv:<id>&format=bibtex`.

| Cite key | Reference |
|----------|-----------|
| `skilling2004nested` | Skilling, J. (2004), "Nested Sampling", AIP Conf. Proc. 735, 395–405 |
| `handley2015polychord` | Handley, Hobson, Lasenby (2015), MNRAS 453(4), 4384–4398 (arXiv:1506.00171) |
| `handley2019anesthetic` | Handley, W. (2019), "anesthetic: nested sampling visualization", JOSS 4(37), 1414 (arXiv:1905.04768) |
| `ashton2022nested` | Ashton et al. (2022), Nature Rev. Methods Primers 2:39 (arXiv:2205.15570) — recommended general-audience citation if a primer reference is wanted |
| `hwangnoh2002` | Hwang & Noh (2002), Phys. Rev. D 65, 124010 — **verify** the correct paper. The existing `references.bib` has Hwang–Noh entries from 2023/2024 (`hwangnoh2023emdef`, `hwangnoh2023graviton`, `hwangnoh2024nonlinear`); these are different works. The 2002 boundary may live in a different paper; cross-check `tidal/inference/_likelihood.py` for the source of the 0.5 constant. |
| `trotta2008bayes` | Trotta, R. (2008), Contemp. Phys. 49, 71 (arXiv:0803.4089) — optional, if the prose invokes the Jeffreys scale for Bayes-factor interpretation. Probably skip given §3.5's "not a Bayes factor" caveat. |

`barker2024poincare` and the Barker software corpus
(`barker2023higgs`, `barker2024psalter`, `barker2025psalter2`,
`barker2025hamilcar`) are already present.

## 5. Result anchors — explicitly deferred

The appendix is **theory only**. Campaign runs are ongoing and no specific
$\log\mathcal{Z}$, $D_{\mathrm{KL}}$, or $d_G$ value should be cited as an
anchor at this stage. The appendix should:

- Describe what each quantity means and how to read its sign and magnitude.
- Forward-reference `\cref{Results}` for actual values once those land in §4.
- Use schematic phrasing — "$\log\mathcal{Z}\approx 0$ would indicate a model
  that is GR-like on average across the prior" — not "we observe
  $\log\mathcal{Z} = X$".

If you are tempted to use real numbers from `hpc_results/` as didactic
examples, resist: campaign results may shift, and the appendix should not
need revision when they do.

## 6. Style guidance distilled from existing appendices

- Use `\paragraph*{Title}` headings, not `\subsection`.
- Semantic CamelCase labels — examples that should exist after writing:
  `EvidenceIntegral`, `PriorVolumeContraction`, `SliceSampling`,
  `SimulationBasedInference`, `EvidencePhysical`, `KLDivergence`,
  `BayesianDimensionality`, `ValidityGates`.
- Prose before and after every equation. Never a bare displayed equation
  followed immediately by another.
- Cite on first mention of every external tool/paper.
- First-person plural ("we") for reported work.
- Cross-references via `\cref{}` with descriptive labels — forward-ref to
  `\cref{PerturbativeRegime}` (App G) for the validity caveat and to
  `\cref{SoftPenaltyFloor}` (App J) for the v3 gate replacement.
- Pattern of existing TODO comments in App J:
  `% Source: <path-or-doc>`. Preserve this convention for any TODOs you
  leave behind.

## 7. What to omit

Resist scope creep — explicitly NOT in scope:

- Algorithm internals beyond what is needed to make $\mathcal{Z}$ and
  $D_{\mathrm{KL}}$ interpretable. The reader knows slice sampling exists;
  the appendix does not implement PolyChord.
- Code architecture beyond cross-references to `docs/tex/inference.tex` and
  `tidal/inference/_*.py`. The appendix is theory; the CLI and class
  hierarchy belong in the project doc.
- Specific campaign physics interpretation. That lives in §4 Results.
- $\mathrm{Var}[\log L]$ as a separate concept from $d_G$ — fold into the
  $d_G$ paragraph.
- The full Ashton et al. (2022) primer content. Cite it as a general-audience
  reference; do not paraphrase its examples (cosmology, gravitational-wave
  astronomy, materials science) — those are distractions for our reader.

## 8. Verification when done

- The handoff is self-contained: the appendix should be writable from this
  document alone, plus the cited source files.
- Every equation in the appendix prose either cites a paper or labels itself
  as a definition introduced here.
- No specific result values appear (theory only).
- The new bibtex entries are added before any `\cite{}` is written.
- `\cref{}` cross-references to `PerturbativeRegime` (App G) and
  `SoftPenaltyFloor` (App J) resolve when the full document compiles.
