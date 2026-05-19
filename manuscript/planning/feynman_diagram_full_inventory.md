# Feynman-diagram perspective for graviton–photon–background couplings

## Full archived inventory and analysis

> **Status:** Archive of the broader investigation that motivated the
> focused App H plan in [feynman_diagram_perspective.md](feynman_diagram_perspective.md).
> Captures the literature survey, EFT-critique discussion, and paper-wide
> diagram inventory so the analysis is not lost, even though only the
> App H subset is being acted on in the current pass.

## 1. Origin and scope of the question

App G's J9 argument
([perturbative_regime.tex:87–116](../sections/appendices/perturbative_regime.tex#L87-L116))
already uses QFT-flavoured prose:

- "Each graviton–photon coupling vertex is linear in the background
  $\bar{\FieldF{}}\propto\Bzero$"
- "An odd number of vertex insertions is needed to interconvert the two
  channels"
- $\Pconv = c_2\Bzero^2 + c_4\Bzero^4 + \dots$

The user flagged 1811.03002 (Jones & Singleton 2018) as a possible
template for adding Feynman-style figures supporting this argument, but
broadened the scope: the real question is whether the diagrammatic
perspective is the right language **for the paper as a whole**, and
where in particular it would genuinely enhance the exposition — not just
for the J9 argument.

## 2. What 1811.03002 actually does

- Title: *Interaction between gravitational radiation and electromagnetic
  radiation*. Review, not vertex-rules paper.
- Diagrams are drawn for $\text{graviton}+\text{graviton}\to\text{photon}+\text{photon}$
  (verbatim from abstract).
- Does **not** draw the $h$–$A$–$\bar{F}$ tadpole vertex App G's J9
  argument depends on. The Gertsenshtein process appears in this review
  through *classical* mode-mixing language, not diagrams.
- Functions as **proof of concept** that diagrammatic language is
  acceptable in this corner of the literature, not as a structural model.

## 3. The wider Gertsenshtein literature is deliberately non-diagrammatic

| Paper | Local path | Diagrams? |
|-------|-----------|-----------|
| Palessandro & Rothman 2023 (`2301.02072`) | `literature/2301.02072/main.tex` | None. Neutrino-mixing analogy; classical mode equations. |
| Hwang & Noh 2023 (`2310.04150`) | `literature/2310.04150/main.tex` | None. Coupled EOMs with explicit metric backreaction. |
| Lella et al. 2024 (`2406.17853`) | `literature/2406.17853/main.tex` | None visible from grep. |
| Domcke et al. 2025 (`2507.16609`) | `literature/2507.16609/Gert.tex` | One *schematic* of emission geometry (line 108). Not a Feynman vertex. Paper deliberately framed as classical scattering on curved spacetime. |
| Burgess 2003 (`gr-qc_0311082`) | `literature/gr-qc_0311082/GRET-jhep.tex` | Heavy: ~25 vertex/diagram references, full Feynman-rule chapter — but for self-energy / vertex corrections in EFT, not the magnetic vertex. |

Pattern: modern Gertsenshtein literature (the immediate field) is
deliberately non-diagrammatic; the EFT-of-gravity literature (the
broader context) is heavily diagrammatic. **Adopting diagrammatic
language is a stylistic choice that aligns the paper with the EFT/QFT
tradition rather than the classical-mode-mixing tradition.**

## 4. Critiques of the quantised-graviton picture

Three threads, with verdicts on whether they bear on this manuscript:

- **EFT validity (Donoghue gr-qc/9512024, Burgess gr-qc/0311082).**
  Mainstream view: gravitons are well-defined effective particles below
  $M_P$. GW frequencies are vastly sub-Planckian. EFT is predictive with
  enormous margin. **Not a critique that bears on us.** Footnote suffices.
- **Metric-backreaction (Hwang–Noh 2310.04150).** Genuine substantive
  caveat, but it concerns *classical* identification of $F_{\mu\nu}$ in
  curved spacetime. A Feynman diagram does not engage this issue;
  classical or diagrammatic, the same field-tensor choice has to be made.
  Manuscript already touches this at [theory.tex:276–279](../sections/theory.tex#L276-L279).
- **Trans-Planckian quantisation (Palessandro 2405.01407).** Contested,
  aimed at frequencies vastly above GW. Out of scope.

**Net**: critiques don't argue against using diagrammatic language; they
argue that the diagrammatic language is a *visualisation*, not a
*derivation*. This is the framing the focused memo recommends.

## 5. Paper-wide inventory of candidate diagram locations

Eight candidate locations were surveyed; ordered by load-bearing-ness.

### 5a. Load-bearing (diagram earns its keep)

**(i) App H "Cubic completions"** —
[pgt_enumeration.tex:269–291](../sections/appendices/pgt_enumeration.tex#L269-L291)

The single highest-value location. Four parity-even cubic operators
generate four distinct vertex topologies after background insertion:

- $\mathcal{R}\times\mathcal{F}^2$ ($\lambda_{1..4}$) — curvature–photon mixing
  (additional Gertsenshtein channel)
- $\mathcal{T}^2\times\mathcal{F}$ ($\kappa_{1..3}$) — torsion-mass shift
  $\propto\Bzero$
- $\mathcal{R}^2\times\mathcal{F}$ ($\rho_{1..8}$) — graviton self-energy in
  $\Bzero$
- $\nabla\mathcal{T}\times\mathcal{F}^2$ ($\mu_{1..5}$) — direct
  torsion–photon conversion channel

**This is the focus of the current implementation pass.**

**(ii) App H "Chern–Simons torsion–EM couplings"** —
[pgt_enumeration.tex:293–318](../sections/appendices/pgt_enumeration.tex#L293-L318)

Three irreducible parity-odd $T$-$A$-$F$ couplings (one per irreducible
torsion sector). **Also acted on in the current pass.**

**(iii) Results §4.3 "Observed channel structure"** —
[results.tex:151–224](../sections/results.tex#L151-L224)

Headline empirical finding: $h_\times \leftrightarrow a_\parallel$ is
torsion-independent in the surveyed propagating-PGT class. The claim is
fundamentally about *which vertex topologies connect which channels* —
the $h$–$A$–$\bar{F}$ minimal vertex connects $h_\times \leftrightarrow
a_\parallel$; the $T$–$A$–$\bar{F}$ non-minimal vertex connects to the
trace sector $\{h_+, h_4, h_9, a_\perp\}$ but **not** to $h_\times
\leftrightarrow a_\parallel$. A channel-vertex coincidence diagram would
crystallize an observation that is currently purely algebraic. **Future
pass.** (Different deliverable: tabular cross-classification, not a
Feynman vertex; deferred.)

**(iv) Discussion §5 "Constructive path: non-minimal couplings"** —
[discussion.tex:51–58](../sections/discussion.tex#L51-L58)

Already in QFT prose: "non-minimal torsion–EM couplings ... break the
kinematic identity by introducing a new vertex outside the admissibility
class". A small diagram showing the new vertex (T-A-B̄ from
$\Del{1}\widetilde{\mathcal{R}}_{[\mu\nu]}\FieldF{}$) reused from §4.3
closes the constructive-path loop visually. **Future pass.**

**(v) App G vertex insertions** —
[perturbative_regime.tex:87–116](../sections/appendices/perturbative_regime.tex#L87-L116)

The original locus. A series-of-three sketch (1, 2, 3 insertions) makes
$\Bzero$-counting visual. Lower priority on its own — the algebra is
short and reads cleanly without a figure — but the diagram comes "for
free" if the App H figures establish the vertex vocabulary. **Future
pass.**

### 5b. Illustrative-only (skipped)

**(vi) Theory §2.1 PGT Lagrangian** —
[theory.tex:122–172](../sections/theory.tex#L122-L172). A vertex
catalogue is what App H provides; duplicating in main body costs main-body
figure budget.

**(vii) Theory §2.2 Linearisation kinetic matrix** —
[theory.tex:269–299](../sections/theory.tex#L269-L299). A 2×2 (or 4×4 for
PGT) kinetic-block diagram is just the matrix drawn as boxes-and-lines.
Information-equivalent to the matrix itself; main body figure budget
better spent on Fig. 1 (sweep results) and Fig. 2 (propagating-torsion
map).

**(viii) App F PGT background** —
[pgt_background.tex](../sections/appendices/pgt_background.tex). Currently
kinematic/geometric (tetrads, irreducible torsion, metric-affine
classification). Coupling-vertex Feynman graphs are a thematic break from
the differential-geometry pedagogy. **Reject** the original "App F"
recommendation explicitly with this reason.

## 6. Deep-dive: the PGT non-minimal coupling vertex inventory

The $\Del{1}\widetilde{\mathcal{R}}_{[\mu\nu]}\FieldF{}^{\mu\nu}$
quadratic non-minimal Lagrangian — the surveyed propagating-torsion
case — generates the following vertices after Palatini decomposition
$\widetilde{R} = R^{LC} + \nabla K + KK$ (with $K$ the contortion built
from $T$) and expansion to linear order on the $g=\eta, \bar T = 0,
\bar A_\mu = -\Bzero z\,\delta_{\mu y}$ background:

1. **$h$–$A$–$\bar F$** (minimal, GR): coefficient $\sim 1$, from
   $\sqrt{-g}F^2$ expansion. Standard Gertsenshtein vertex.
2. **$T$–$A$–$\bar F$** (non-minimal, PGT): coefficient $\sim\Del{1}$,
   from $\nabla K \cdot F$ contraction. The "new vertex" of
   [discussion.tex:54](../sections/discussion.tex#L54).
3. **$h$–$T$–$\bar F$** (mixed, PGT): coefficient $\sim\Del{1}$, from
   $\widetilde R^{(1)}_h \cdot F$ contraction (metric-perturbation
   contribution to $\widetilde R$ at linear order).

Plus the cubic-sector vertices in §5a(i) above and the Chern–Simons
vertices in §5a(ii).

The empirical channel-structure observation is exactly the statement
that vertices (2), (3) and the CS family **vanish** when projected onto
the $\{h_\times, a_\parallel\}$ subspace, while (1) survives. A
projection-structure diagram would take a paragraph of dense
kinetic-matrix algebra and make it visually immediate. This is the place
where the Feynman-diagram perspective delivers the most physical content
per square inch, and it is *specific to the non-minimal-coupling PGT
case* — for pure-GR Gertsenshtein only (1) exists and the diagram has
nothing to compare against.

## 7. Why App H, not App F or main body

The user's original instinct was App F. On reflection:

- App F is purely kinematic/geometric (tetrads, irreducible torsion,
  metric-affine classification). A coupling-vertex Feynman graph is
  thematic break.
- Main-body theory §2.2 / §2.3 figure budget is better spent on Figs. 1
  and 2 of `results.tex` (the headline empirical results).
- App H is exactly where the dense vertex-topology enumeration already
  lives in prose. Replacing prose with diagrams there is highest-value.

## 8. Choice of leg-labelling convention

Dynamical external legs are labelled $h$, $A$, $T$ — bare field
symbols, no "$\delta$" prefix. The "$\delta$" used in the linearised
Lagrangian
([theory.tex:217–265](../sections/theory.tex#L217-L265)) is for
distinguishing dynamical from background in *algebraic manipulation*; in
the diagram, the leg-style itself (curly/snake/dashed) plus the explicit
$\bar F$ convention already encode that distinction, so the prefix is
redundant clutter. Standard Feynman-diagram convention agrees. The
figure caption opens with one sentence stating the convention.

## 9. What the focused memo recommends acting on now

- **Figures**: vertex-topology figure(s) for App H §C (cubic completions
  + Chern–Simons couplings). Figure count driven by structural-
  distinction criteria — see
  [feynman_diagram_perspective.md](feynman_diagram_perspective.md).
- **Prose**: figure insertion + connecting sentences + one caveat
  paragraph in App H §C.
- **Bib**: add `burgess2003quantum` (for EFT-context citation in caveat)
  and use existing `hwangnoh2023graviton` (for classical-caveat
  citation).

## 10. What is deferred to future passes

- §5a(iii): channel-vertex coincidence diagram in results §4.3.
- §5a(iv): non-minimal-coupling vertex diagram in discussion §5.
- §5a(v): vertex-insertion series in App G near J9.
- 1811.03002 download (only needed if formally cited; current pass cites
  Burgess + Hwang–Noh, both already local).
