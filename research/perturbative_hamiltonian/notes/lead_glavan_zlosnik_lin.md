# Lead: Glavan-Zlosnik-Lin 2024 (arXiv:2311.17459)

**Citation:** D. Glavan, T. Zlosnik, C. Lin, "Hamiltonian analysis of metric-affine $R^2$ theory", arXiv:2311.17459 (JHEP-style preprint, dated 2024; revised v2 2024-04-29).
**Local file:** `/workspaces/torsion-gertsenshtein/literature/2311.17459/Hamiltonian_metric-affine_R2.tex` (4201 lines)
**Date of read:** 2026-04-27
**Investigator:** Claude (Opus 4.7), Phase 2.1 deep read

---

## Section 1 — GZL's actual construction (verbatim)

### 1.1 Lagrangian and starting framework

The starting action is metric-affine $f(R)$ with $f = R^2$ in $D=4$:

> "$S\bigl[g_{\mu\nu}, \Xi_{\alpha\mu\nu}\bigr] = \int d^4x\,\sqrt{-g}\,R^2$"
> (eq. (Palatini R2), `Hamiltonian_metric-affine_R2.tex:975-977`)

The dynamical fields are the metric $g_{\mu\nu}$ AND the independent affine
connection $\Gamma^\alpha{}_{\mu\nu}$. The connection is decomposed as
$\Gamma = \mathring{\Gamma} + \Xi$, where $\mathring{\Gamma}$ is Levi-Civita
and $\Xi$ is the **distortion tensor** (eq. (connection split),
`tex:411-415`). $\Xi$ is further decomposed (`tex:454-473`) into
contortion $K$ (torsion-built) and disformation $L$ (non-metricity-built).

GZL prefer the projective-invariant $B$-tensor combination:

> "$B_{\alpha\mu\nu} \equiv \Xi_{\mu\alpha\nu} - \Xi_{\nu\alpha\mu} = 2K_{\mu\alpha\nu} + Q_{\mu\alpha\nu} - Q_{\nu\alpha\mu}$"
> (eq. (B tensor), `tex:620-625`)

The form of $R$ in this basis is (`tex:633-640`):

> "$R = \mathring{R} - \mathring{\nabla}_\mu B_\nu{}^{\nu\mu} - \tfrac14 B^\mu{}_{\mu\alpha}B_\nu{}^{\nu\alpha} + \tfrac14 B_{\mu\nu\alpha}B^{\nu\mu\alpha} + \tfrac14 Q^\mu{}_{\mu\alpha}Q_\nu{}^{\nu\alpha} - \tfrac14 Q_{\mu\nu\alpha}Q^{\nu\mu\alpha}$"

Crucially, GZL exploit **three** local symmetries (`tex:519-607`):

1. **Diffeomorphism invariance** (universal).
2. **Projective invariance** (universal in $f(R)$ metric-affine, `tex:546-562`):
   $\Xi_{\alpha\mu\nu} \to \Xi_{\alpha\mu\nu} + g_{\alpha\nu}\zeta_\mu$
   (eq. (projective transformation), `tex:551-557`).
3. **Weyl invariance**, which holds **specifically** for $f(R) = R^2$ in
   $D=4$ (`tex:565-607`):
   $g_{\mu\nu} \to e^\sigma g_{\mu\nu}$, with the matching shift on $\Xi$
   in eq. (Weyl trans), `tex:571-581`.

> "Therefore, this particular case, that we examine in this work, should be treated with special care since it ought to exhibit more first-class constraints compared to the rest of the theories in (\ref{f(R) action})." (`tex:604-607`)

This is GZL's key methodological statement: the $R^2$ case is **special**
because it has more gauge symmetry, hence more first-class constraints,
hence fewer propagating DOF.

### 1.2 Method: Dirac-Bergmann constraint analysis on an ADM-decomposed action

The plan is laid out in the introduction (`tex:319-340`):
ADM decomposition of $(g, \Xi)$ → canonical action → Dirac analysis →
phase-space reduction by solving second-class constraints → equivalence
to Einstein-Hilbert + $\Lambda$.

The constraint analysis (Sec. 5, `tex:2290-2730`) classifies primary,
secondary and tertiary constraints into three sectors (anti-symmetric,
traceless symmetric, scalar trace). The closing tally is

> "$N_{\rm 2nd}=26$ second-class constraints and $N_{\rm 1st}=6$
> first-class constraints. Given that there are $N_{\rm can}=42$
> canonical variables (not including Lagrange multipliers) we can infer
> the number of physical propagating degrees of freedom,
> $N_{\rm phy} = \tfrac12(N_{\rm can} - N_{\rm 2nd} - 2N_{\rm 1st}) = 2$."
> (`tex:2706-2725`)

i.e. just the graviton.

### 1.3 Higher-derivative-term handling: scalar-auxiliary lift

GZL's treatment of the $R^2$ higher-derivative term is the standard
**scalar-auxiliary (Brans-Dicke) lift**, NOT a Dirac-Ostrogradsky
reduction. Crucially this is performed at the Lagrangian level BEFORE
the Hamiltonian analysis:

> "We introduce an auxiliary scalar $\varphi$ that is equal to the Ricci scalar on-shell, which is enforced by a Lagrange multiplier $\lambda$"
> (`tex:802-806`, eq. (intermediate S))

> "$S[g_{\mu\nu},\Xi_{\alpha\mu\nu},\varphi,\lambda] = \int d^Dx\,\sqrt{-g}\,[f(\varphi) + \lambda(\varphi - R)]$"
> (eq. (intermediate S), `tex:807-812`)

Solving for $\lambda = -f'(\varphi)$ produces:

> "$S[g_{\mu\nu},\Xi_{\alpha\mu\nu},\varphi] = \int d^Dx\,\sqrt{-g}\,[f(\varphi) - \varphi f'(\varphi) + f'(\varphi) R]$"
> (eq. (intermediate S 2), `tex:825-832`)

At this point the **`R` is still the metric-affine $R$, not $\mathring{R}$**, and the
$R^2$ term has been replaced by a **linear-in-$R$** term with an
auxiliary scalar coefficient. There is no quadratic-curvature term left
in the Lagrangian. The procedure assumes $f''(\varphi)\not\approx 0$
(`tex:835`).

### 1.4 Constraint structure via field redefinitions

GZL combine the auxiliary-scalar lift with conformal rescaling of the
metric (eq. (g conformal rescaling), `tex:891-894`) and a compensating
shift of $\Xi$ (eq. (3tensor frame trans), `tex:900-908`). For
$f(R)=R^2$ specifically, the auxiliary scalar **disappears entirely**
(`tex:980-994`):

> "It has a very special property that applying to it the field redefinitions (g conformal rescaling) and (3tensor frame trans) completely removes the auxiliary scalar"

leaving:

> "$S = \int d^4x \sqrt{-g}\,(1/\kappa^2)\bigl[\mathring{R} - 1/(4\kappa^2) + \Xi^\mu{}_{\mu\nu}\Xi^{\nu\alpha}{}_\alpha - \Xi^\mu{}_{\nu\alpha}\Xi^{\nu\alpha}{}_\mu\bigr]$"
> (`tex:984-994`)

**This is the key trick that makes the GZL paper tractable.** The
$R^2$ Lagrangian is converted, by symmetry-respecting field redefinitions,
into Einstein-Hilbert + Λ + algebraic-in-$\Xi$ terms. The remaining
$\Xi$-dependence is purely algebraic (no derivatives of $\Xi$), so when
solved for on-shell $\Xi$ drops out, leaving Einstein-Hilbert with
cosmological constant (`tex:999-1008`):

> "$S[g_{\mu\nu}] = \int d^4x \sqrt{-g}\,(1/\kappa^2)[\mathring{R} - 1/(4\kappa^2)]$"

### 1.5 What follows in the Hamiltonian sections

Sections 4-7 of GZL build the canonical formulation directly from the
**post-redefinition** action and verify that the Dirac analysis recovers
the "just the graviton" spectrum. The constraint structure they extract
(primary/secondary/tertiary constraints, first/second class, gauge
generators) lives entirely in the **2nd-order Lagrangian** that came out
of the field-redefinition step, NOT in the original $R^2$ action.

The introduction explicitly tells us this (`tex:319-340`); the structure
of the constraint analysis sections follows the Brans-Dicke-recast
Lagrangian as the input to ADM decomposition.

---

## Section 2 — Constraint-promotion analogue (Q2)

### Verdict: NO published constraint-promotion analogue exists in GZL.

GZL is built around a Lagrangian that has **fixed** the higher-derivative
problem before the Hamiltonian analysis ever begins. The auxiliary
scalar $\varphi$ trick (`tex:807-832`) replaces $R^2$ with $\varphi R$
plus a non-derivative potential. The conformal+distortion field
redefinitions (`tex:891-908`) then absorb $\varphi$ entirely (specific
to the $f=R^2$ case in $D=4$), leaving only **algebraic** dependence on
$\Xi$. There is no small parameter analogous to TIDAL's $b_5$ that
"turns on" $R^2$; the construction is at fixed coupling and the entire
$R^2$ structure is lifted before phase-space analysis.

Searching the local TeX for "limit", "small parameter", "perturbative",
"vanish", "critical":

- "In the limit $D\to 4$" (`tex:885`) — dimensional limit, not coupling.
- "$f''(\varphi)\not\approx 0$" (`tex:835`) — the only place a "$=0$
  limit" is even mentioned; GZL explicitly EXCLUDE this case as
  outside the scope of their procedure.

The condition $f''(\varphi)=0$ is the parameter limit closest to TIDAL's
$b_5\to 0$: $f''=2$ for $f=R^2$ (so non-vanishing), but for a general
$f(R) = R + \epsilon R^2$, $f''(\varphi) = 2\epsilon$ vanishes at
$\epsilon\to 0$. **At that limit the auxiliary-scalar-lift breaks down**
(per GZL's own caveat at `tex:835`), and the construction fails. So GZL
not only do not handle the constraint-promotion limit; their entire
methodology requires staying away from it.

This is a critical structural finding: **GZL's higher-derivative trick
inherits exactly the singularity that TIDAL's b5 → 0 limit hits.** When the
parameter that turns on $R^2$ goes to zero, $f''\to 0$, and the
auxiliary-scalar elimination is non-invertible. GZL handle the case
$f = R^2$ where the coefficient is finite; they do not handle
$f = R + b_5 R^2$ with $b_5 \to 0$.

This is independently confirmed by Bellorin 2025 in the metric-only
quadratic-gravity setting (see Section 4 below):

> "One cannot obtain the ADM Hamiltonian of general relativity as a
> smooth limit of the case we have considered."
> (Bellorin 2025, arXiv:2506.07305 abstract; cited GZL in fn. 17)

---

## Section 3 — Structural transferability (Q3)

### 3.1 Framework comparison

| Aspect | GZL 2024 | TIDAL b5·R̃² PGT |
|---|---|---|
| Geometry | Metric-affine (Palatini): independent connection with both torsion AND non-metricity | Pure Poincaré gauge (PGT): vielbein + spin connection, torsion only, metric-compatible |
| Curvature object | $R$ = Ricci scalar of full $\Gamma$ (`tex:399-402`) | $\tilde{R}$ = Riemann-Cartan Ricci scalar from the spin connection |
| Lagrangian | $f(R) = R^2$ alone (`tex:975-977`) | $L = (1/\kappa^2)\tilde{R} + \alpha_i I_i + b_5 \tilde{R}^2 - F^2/4$ (PGT + EM) |
| Higher-derivative parameter | None — coefficient of $R^2$ is unity | $b_5$ — the controlled small parameter |
| Symmetries used | diff + projective + Weyl (specific to $f=R^2$, $D=4$) | diff + Lorentz only (standard PGT — no projective, no Weyl) |
| DOF spectrum | "just the graviton" (2 DOF, `tex:2706-2725`) | sector-dependent: graviton + torsion modes (axial, trace, tensor-q) |
| Phase-space dim | $N_{\rm can}=42$ canonical (post-redefinition) | 38-component PGT phase space (per round1 synthesis) |

### 3.2 Symmetries: the make-or-break difference

**Projective symmetry** (`tex:546-564`): GZL exploit the fact that
$R$ (Ricci scalar of the FULL connection) is invariant under
$\Xi_{\alpha\mu\nu} \to \Xi_{\alpha\mu\nu} + g_{\alpha\nu}\zeta_\mu$.
This is a 4-component gauge symmetry that fixes the gauge-redundant part
of $\Xi$ and is essential to making the canonical analysis tractable.
TIDAL's PGT does **not** have this symmetry — the spin connection in PGT
transforms under local Lorentz, not under projective shifts. The
projective gauge freedom is unique to metric-affine where the connection
has 64 independent components; in PGT (metric-compatible, vielbein
formulation) this gauge has already been fixed.

**Weyl invariance** (`tex:565-607`): GZL specifically exploit Weyl
invariance of $R^2$ in $D=4$: "this particular case ... ought to
exhibit more first-class constraints compared to the rest of the
theories in (f(R) action))" (`tex:604-607`). Weyl invariance is
WHY $R^2$ propagates only the graviton in metric-affine: the additional
gauge symmetry kills modes that would otherwise propagate.

TIDAL's $b_5 \tilde{R}^2$ PGT is **not** Weyl-invariant. The Einstein-
Hilbert piece $(1/\kappa^2)\tilde{R}$ explicitly breaks Weyl symmetry
(introduces a scale $\kappa$), and the $\alpha_i I_i$ torsion-squared
terms generically break it too. The "just the graviton" outcome of GZL
relies on the local Weyl symmetry that TIDAL's Lagrangian does NOT have.

### 3.3 Phase-space dimension

GZL's phase space (`tex:2708`): $N_{\rm can}=42$ canonical variables
(this is post-field-redefinition; original metric+distortion has more).
After Dirac reduction: 26 second-class + 6 first-class → 2 propagating
DOF.

TIDAL's PGT b5·R̃² phase space (per `MEMORY.md` and round1 synthesis):
38-component vielbein-PGT phase space (10 metric perturbation components
$h_{\mu\nu}$, 24 spin-connection $\omega^{ab}{}_\mu$, after gauge fixing
to suitable gauge).

Different dimensions, different gauge structures. The 38-component
analysis cannot be directly mapped onto GZL's 42-component analysis.

### 3.4 Higher-derivative handling: scalar-auxiliary vs. b5·R̃²

GZL's auxiliary-scalar trick (Brans-Dicke lift) **works** for
$f(R) = R^2$ in metric-affine because:
- $R$ depends on $\Xi$ algebraically (no derivatives of $\Xi$ appear in
  $R$ via $\Xi^2$ terms; $\nabla\Xi$ terms can be absorbed by IBP and
  field redefinitions because of projective symmetry).
- The auxiliary scalar $\varphi$ (=R on-shell) becomes a Lagrange
  multiplier that does NOT propagate, killed by the Weyl-redefined
  $\Xi$ shift (`tex:900-908`).

TIDAL's $b_5 \tilde{R}^2$ in PGT: $\tilde{R}$ depends on the spin
connection $\omega$ via $\tilde{R} \sim \partial\omega + \omega\omega$.
Squaring gives $(\partial\omega)^2$ terms — genuine higher-derivative
kinetic terms for $\omega$. Performing the same scalar-auxiliary trick
gives $\varphi(\tilde{R}) \sim \varphi \partial\omega + \varphi \omega^2$,
which does not eliminate the $\omega$ derivatives — the lift converts
"4th-order in metric perturbation $h$" into "2nd-order plus an auxiliary
scalar coupled to $\partial\omega$", i.e. it migrates the problem from
the metric block to the connection block but does not solve it.

### 3.5 Verdict on transferability

**The GZL technique is structurally tied to two specific gauge
symmetries (projective + Weyl) that TIDAL's b5·R̃² PGT does not have.**
Without those symmetries:
- The auxiliary-scalar lift does not produce a 2nd-order Lagrangian
  (the $\varphi$ does not eliminate cleanly).
- The "just the graviton" outcome does not follow.
- The "field redefinition removes the auxiliary scalar" trick
  (specific to $f=R^2$ in $D=4$ via Weyl) is not available.

The most that transfers is the **methodological framing**: think of the
$R^2$ term as a Lagrangian-level lift to a 2nd-order theory plus an
auxiliary scalar plus algebraic constraints, then do Dirac-Bergmann on
the lifted system. But the specific implementation details (which
symmetries fix which modes, which auxiliary fields can be eliminated,
which constraints are first-class) do not transfer without major
modification.

---

## Section 4 — Forward citations (Q4)

INSPIRE-HEP forward-citation query (recid 2727946) returned **7 citing
papers** as of 2026-04-27. I retrieved abstracts and key claims for the
relevant ones. Three are directly load-bearing for the question whether
GZL's technique extends to TIDAL's case.

### 4.1 Barker & Glavan 2025, "Spectrum of pure R² gravity: full Hamiltonian analysis" (arXiv:2510.08201)

This is by the same Glavan, with Will Barker. **Directly addresses
constraint-promotion at linearisation.** Verbatim from the abstract:

> "We perform a full Hamiltonian constraint analysis of pure
> Ricci-scalar-squared (R²) gravity to clarify recent controversies
> regarding its particle spectrum. While it is well established that
> the full theory consistently propagates three degrees of freedom,
> we confirm that its linearised spectrum around Minkowski spacetime
> is empty. ... The mechanism for this phenomenon is a change in the
> nature of the constraints upon linearisation: ten second-class
> constraints of the full theory become first-class, while the three
> momentum constraints degenerate into a single constraint. ...
> such backgrounds are surfaces of strong coupling in field space,
> where the dynamics of perturbations becomes nonperturbative."
> (Barker-Glavan 2025, abstract)

**Crucially**: this is *metric-only* R² (Stelle-tradition), not metric-
affine, and not PGT. The "second-class becomes first-class" mechanism
is **the closest published analogue to TIDAL's b5→0 constraint
promotion** — same direction (more first-class = more gauge symmetry =
fewer propagating DOF).

But the analogue is **inverted in direction**:
- Barker-Glavan: full nonlinear → linearised: 10 second-class become
  first-class. The constraint-promoting limit is a perturbation
  expansion around Minkowski, not a parameter limit in the Lagrangian.
- TIDAL: $b_5 \to 0$: constraints that were *primary* at $b_5\neq 0$
  become *broken* (or move to different classes). This is a parameter
  limit in the Lagrangian, not a perturbation expansion.

What Barker-Glavan offer **structurally**: an explicit example where
constraints change class *in a controlled limit*, with the verdict
"strong coupling, perturbation theory breaks down."

What they do **NOT** offer: a recipe to compute physical observables
across the singular surface. Their conclusion is that perturbation
theory fails ("becomes nonperturbative"), not that the singular surface
is computable.

This is consistent with TIDAL's own findings: Round 1 Agent A showed
$\det(M)\propto b_5^N$ (rank-jump), and Round 1 Agents D/G showed three
no-go arguments for local first-order auxiliary lifts. Barker-Glavan
2025 confirms this is the right diagnosis ("strongly coupled, non-
perturbative") and inadvertently confirms that no published recipe
handles the regime — including their own.

### 4.2 Karananas 2024, "Particle content of (scalar curvature)² gravities revisited" (arXiv:2407.09598)

> "Studying the spectrum of (pure) R² gravity on Minkowski background
> inevitably results into a Catch-22: any consistent interpretation of
> its particle dynamics dictates that no accidental gauge symmetries
> emerge, a requirement that cannot be fulfilled when the theory is
> studied on Minkowski."
> (Karananas 2024 abstract)

> "making R² gravity infinitely strongly-coupled."
> (Karananas 2024 abstract)

> "Degrees of freedom are reintroduced via interactions, making R²
> gravity infinitely strongly-coupled."
> (Karananas 2024 elaboration)

Note: this paper provides an EH + R² action:

> "$S' = S^{\rm EH} - (1/(12 f^2)) S^{R^2}$"

with linearised quadratic action

> "$S'_2 = (1/2)\int d^4x[(\partial\rho h_{\mu\nu}\partial^\rho h^{\mu\nu} - 2\partial_\mu h^{\mu\nu}\partial^\rho h_{\nu\rho} + 2\partial_\mu h \partial_\nu h^{\mu\nu} - \partial_\mu h\partial^\mu h) + (2/(3 f^2 M_P^2))(\partial^\mu\partial^\nu h_{\mu\nu} - \Box h)^2]$"
> (Karananas 2024)

This **is** an Einstein-Hilbert plus $R^2$/$f^2$ Lagrangian where the $R^2$
coefficient is $1/(12f^2)$ — but the analysis is performed at finite
$f$, with the discontinuity between pure-$R^2$ and pure-EH highlighted
qualitatively. The vanishing-coefficient limit is **not analysed
explicitly**.

### 4.3 Karananas 2024, "Particle content of (scalar curvature)² metric-affine gravity" (arXiv:2408.16818)

This is the metric-affine extension of 2407.09598. Verbatim:

> "Linearizing metric-affine (scalar curvature)² gravity—an 'umbrella'
> theory that includes as special cases the metrical, Einstein-Cartan,
> and Weyl quadratic models—on top of Minkowski spacetime leads to
> (numerous) accidental gauged symmetries. ... Such undesirable
> symmetries are absent already at the leading nontrivial order in
> perturbations on non-flat backgrounds, e.g. de Sitter spacetime,
> which are the appropriate ones for studying the particle dynamics."
> (Karananas 2408.16818 abstract)

**Important inclusion**: "Einstein-Cartan quadratic gravity when the
connection is taken to be metric compatible" (from the body) — this
covers the parity-even sector of TIDAL's PGT b5·R̃². However, no
isolated PGT/torsion analysis is presented; it is one limit of a larger
parameter space. **Vanishing R² coefficient limit is not analysed.**

### 4.4 Bellorin 2025, "Hamiltonian equations of motion of quadratic gravity" (arXiv:2506.07305)

**This is the most directly load-bearing forward citation.** Cites GZL
in footnote 17.

> "We compute explicitly the equations of motion of the Hamiltonian
> formulation of quadratic gravity. This is the theory with the most
> general Lagrangian with terms of quadratic order in the curvature
> tensor (discarding the cosmological constant). ... We compare the
> linear equations with the covariant field equations, finding that,
> if general-relativity terms are active, the linear Hamiltonian
> formulation is valid only if the perturbative spatial metric is
> traceless, a condition that can be freely imposed by recurring to
> an arbitrary function."
> (Bellorin 2025 abstract)

> "One cannot obtain the ADM Hamiltonian of general relativity as a
> smooth limit of the case we have considered."
> (Bellorin 2025, generic-case statement)

This is **the metric-only quadratic-gravity analogue of TIDAL's
b5→0 constraint barrier**. Bellorin's framework is metric-only Stelle
quadratic gravity ($\kappa^{-2}, \alpha, \beta$ for EH, Ricci², R²
coefficients). The discontinuity between "$\alpha,\beta\to 0$" (i.e.,
recover GR) and "generic case" is **explicitly stated** as a non-smooth
limit. Bellorin cites GZL as "the metric-affine variant ... related but
distinct prior work."

**TIDAL b5→0 sits in exactly this metric-only-quadratic-gravity
discontinuity, but in PGT (vielbein, pure torsion).** Bellorin's
published statement that the limit is non-smooth confirms that the
TIDAL barrier is recognised in the metric-only literature; what does
not exist (still) is the PGT-specific extension.

### 4.5 Other citing papers (lower priority for this question)

- **Capozziello 2026 (arXiv:2604.19310)**, "Extrinsic geometry and
  Hamiltonian analysis of symmetric teleparallel gravity" — STG case;
  no constraint-promotion limit, no PGT extension.
- **Aoki 2023 (arXiv:2310.16007)**, "Cosmological Perturbation Theory
  in Metric-Affine Gravity" — perturbation theory around FLRW; "rich
  perturbation spectrum" but no DOF discontinuity at parameter limits.
- **Barker 2026 (arXiv:2512.25007)**, "Fast Poisson brackets and
  constraint algebras in canonical gravity" — algorithmic tool
  ("Hamilcar"); applied only to GR, pure $R^2$, and Goroff-Sagnotti;
  not yet applied to PGT.

### 4.6 Bottom line on forward citations

No published paper applies GZL's technique to PGT. No paper extends
GZL's technique to handle the parameter-vanishing limit. The closest
constraint-promotion analogues — Barker-Glavan 2510.08201 and
Karananas 2407.09598 — recognise the constraint-promotion phenomenon
(in metric-only $R^2$) but do not produce a Hamiltonian recipe
across the singular surface. Their published verdict is "strongly
coupled, non-perturbative", which is **consistent** with TIDAL's
Round 1 no-go arguments rather than a circumvention of them.

Bellorin 2506.07305's metric-only "no smooth GR limit" statement is
the most direct external confirmation that the TIDAL b5→0 barrier
sits in published-recognised territory — but published in metric-only
form, not in PGT form.

---

## Section 5 — Verdict

**(c) — Not transferable.**

GZL's technique applies to a Lagrangian (pure $f(R)$ metric-affine in
$D=4$ Weyl-symmetric form) whose structural properties are fundamentally
different from TIDAL's b5·R̃² PGT in three independent ways, each
sufficient on its own to break the transfer:

1. **Symmetry mismatch.** GZL exploit projective + Weyl gauge symmetries
   that are explicit gauge symmetries of $f(R)$ metric-affine in $D=4$.
   TIDAL's PGT has neither. The "just the graviton" outcome is a
   consequence of the Weyl gauge fixing first-class constraints; without
   that gauge symmetry, the constraint structure is genuinely different
   and richer (`tex:604-607`, "this particular case ... ought to exhibit
   more first-class constraints compared to the rest").

2. **Higher-derivative-trick singularity.** GZL's auxiliary-scalar
   ($\varphi = R$ on-shell) lift requires $f''(\varphi)\not\approx 0$
   (`tex:835`). For $f = R + b_5 R^2$ this is $f''=2b_5$, vanishing
   at $b_5\to 0$. So GZL's own construction is **non-invertible at the
   constraint-promotion limit**. The auxiliary-scalar lift does not
   solve the b5→0 problem — it inherits it.

3. **Geometry mismatch.** GZL's metric-affine connection has 64
   independent components, fixed partly by projective symmetry. TIDAL's
   PGT has a vielbein + spin connection (40 + 24 components), no
   non-metricity, no projective symmetry. The phase-space dimension
   (42 in GZL post-redefinition, 38 in TIDAL after gauge fixing) and
   the gauge-symmetry structure (diff + Lorentz only in TIDAL) are
   incompatible at the level of the Dirac analysis itself.

The most that transfers is **methodological framing**: think of the
$R^2$ term as a Lagrangian-level lift before doing Dirac-Bergmann.
This is the same idea TIDAL's Round 2 Agent C/F (axial-sector Bopp-
Podolsky) and Round 3 Agent J (Curtright Stückelberg) already explored,
with the verdict that the lift either fails (Round 1 Agent D/G no-gos
for h_4/h_7/h_9 metric subspace) or succeeds only for specific sectors
at linear-flat order (axial scalar-auxiliary, Stückelberg at 1+1D toy).

GZL does not provide a missing sector recipe. It does not provide a
parameter-vanishing recipe. And it does not extend to vielbein PGT.

### Concrete next-step recommendation

**Do not invest theoretical work in adapting GZL to TIDAL.** The three
structural mismatches above are independent and severe.

The forward-citation deep dive surfaced one finding worth its own
attention: **Barker-Glavan 2025 (arXiv:2510.08201) and Karananas 2024
(arXiv:2407.09598, 2408.16818)** all describe the constraint-promotion
phenomenon — second-class becoming first-class at linearisation, "infinite
strong coupling," "accidental gauged symmetries" — in metric-only or
metric-affine R²-type Lagrangians. Their published verdict is that the
Minkowski/perturbative limit is non-perturbative.

This is **published evidence that TIDAL's barrier is real and recognised**,
not a missing solution. The right operational use of these papers in
TIDAL's documentation is:

- Cite Barker-Glavan 2510.08201 in `docs/tex/perturbative_reduction_constraint_barrier.tex`
  as published evidence that constraint-class promotion at perturbative
  expansion around a vacuum is a recognised phenomenon ("strong-coupling
  surfaces in field space"). This **strengthens** TIDAL's "this is a
  real barrier" framing.
- Cite Bellorin 2506.07305 footnote 17 as published recognition of GZL
  as the metric-affine prior work, **and** Bellorin's "no smooth GR
  limit" claim as the metric-only quadratic-gravity analogue of TIDAL's
  b5→0 barrier.
- Cite Karananas 2408.16818 as published demonstration that the
  Einstein-Cartan (metric-compatible, parity-even) limit of metric-affine
  scalar-curvature² is one corner of a larger umbrella, with "accidental
  gauged symmetries" obstructing the Minkowski analysis. The TIDAL b5·R̃²
  in PGT is the parity-odd extension of this corner; the published
  obstruction at parity-even is **prior art** for the parity-odd case
  TIDAL hits.

These three forward-citation findings collectively reinforce the
FINAL_ASSESSMENT.md recommendation: **document the limitation, ship the
existing flagship measurement**, and treat the no-go publication
(Publication C in the original investigation framing) as the
defensible theoretical writeup if supervisor requests it. The 25-year-
unsolved framing should be explicitly tied to the higher-derivative-PGT
case (per Meta-N's qualification), with Barker-Glavan 2025 + Karananas
2407.09598 + Bellorin 2506.07305 as published companions confirming the
phenomenon.

**Total work to absorb this finding into TIDAL docs: 1-2 hours.**
**No theoretical extension of GZL is recommended.**

---

## Section 6 — Local-paper line-number references (cumulative)

All references are to `/workspaces/torsion-gertsenshtein/literature/2311.17459/Hamiltonian_metric-affine_R2.tex`.

| Claim | Lines |
|---|---|
| Title and author block | 71-101 |
| Abstract | 112-124 |
| Introduction: f(R) propagates 3 DOF metric, just graviton metric-affine | 198-205 |
| Plan: ADM → canonical → Dirac → phase-space reduction | 319-340 |
| Connection split: Levi-Civita + distortion | 411-415 |
| R in $\Xi$ form | 432-437 |
| Distortion tensor decomposition (K + L) | 454-473 |
| Torsion as antisymmetric part of $\Gamma$ | 478-484 |
| Non-metricity definition | 486-495 |
| Diffeomorphism transformation | 519-544 |
| Projective transformation rule | 546-562 |
| Weyl transformation rule | 565-607 |
| "Special care" for $R^2$ Weyl-invariance | 604-607 |
| B-tensor definition | 620-625 |
| R in B and Q form | 633-640 |
| List of projective-invariant components | 692-712 |
| List of quadratic invariants (Weyl + projective) | 718-751 |
| Equivalence to Einstein-Hilbert (section heading) | 779-781 |
| Auxiliary-scalar lift (eq. intermediate S) | 802-812 |
| Auxiliary-scalar lift after $\lambda$ elimination (eq. intermediate S 2) | 825-832 |
| $f''(\varphi)\not\approx 0$ caveat | 835 |
| Brans-Dicke equivalence with $\sigma=f'(\varphi)$ | 862-886 |
| Conformal rescaling (eq. g conformal rescaling) | 891-894 |
| Distortion field redefinition (eq. 3tensor frame trans) | 900-908 |
| EH + auxiliary scalar after redefinitions | 911-928 |
| Cosmological constant from on-shell $\varphi_0$ | 943-951 |
| Final EH action (general $f$) | 956-966 |
| Pure $R^2$ action (eq. Palatini R2) | 974-978 |
| $R^2$-specific scalar elimination | 980-994 |
| Final EH+Λ action ($f=R^2$ case) | 999-1008 |
| ADM section heading | 1035-1037 |
| Constraint analysis section heading | 2290-2291 |
| Summary of constraint structure (table + DOF count) | 2624-2731 |
| $N_{\rm phy} = 2$ DOF count | 2706-2725 |
| Discussion of more general Weyl-invariant theories beyond R² | 3403-3477 |

---

## Section 7 — Forward-citation references

| Paper | arXiv | Cites GZL? | Constraint promotion at parameter→0? |
|---|---|---|---|
| Barker-Glavan, "Spectrum of pure R² gravity" | 2510.08201 | (same author) | Yes, but at linearisation, not param→0 |
| Karananas, "Particle content of (scalar curvature)² gravities revisited" | 2407.09598 | Likely (not verified — html missing reflist) | Strong-coupling discussion, but param fixed |
| Karananas, "Particle content of (scalar curvature)² metric-affine gravity" | 2408.16818 | Likely | "Accidental gauged symmetries" at Minkowski |
| Bellorin, "Hamiltonian equations of motion of quadratic gravity" | 2506.07305 | Yes (footnote 17) | Yes, "no smooth GR limit" — closest published analogue |
| Capozziello, "Extrinsic geometry and Hamiltonian analysis of symmetric teleparallel gravity" | 2604.19310 | Yes | No |
| Aoki, "Cosmological Perturbation Theory in Metric-Affine Gravity" | 2310.16007 | Yes | No |
| Barker, "Fast Poisson brackets and constraint algebras in canonical gravity" | 2512.25007 | Yes | "Strongly coupled modes" but tool is generic |

INSPIRE-HEP query: `refersto:recid:2727946`, retrieved 2026-04-27.
Total citing papers: 7.
