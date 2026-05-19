# Correction Log W1 — `docs/tex/perturbative_reduction_constraint_barrier.tex`

**Date:** 2026-04-27
**Agent:** Correction Agent W1
**Source of corrections:** `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`,
§ "Concrete documentation corrections" (9 items).
**Outcome:** All 9 corrections applied or noted as N/A. Build verified
(`make` → 255 pages, 1.95 MB PDF, no LaTeX errors, no undefined refs;
9 BibTeX warnings are pre-existing missing-`journal` fields unrelated
to this audit).

---

## Correction 1 — Krupka–Voicu "Theorem 1" → "Definition 1"

**Audit justification (FINAL_ASSESSMENT.md §"Concrete documentation
corrections" item 1, supported by Meta-K K2):**
> "Replace 'Krupka-Voicu Theorem 1' with 'Krupka-Voicu Definition 1
> (canonical variational completion)' throughout. Justification: Meta-K
> verified via `literature/1406.6646/main.tex` that the paper has
> Definition 1, not Theorem 1."

**Local-paper verification:** `literature/1406.6646/main.tex` lines
365–374:
> `\begin{definition}\nThe canonical variational completion of a source form
> $\varepsilon \in \Omega _{n+1}^{r}(Y),$ is the source form $\tau (\varepsilon )$
> given by the difference between the Euler-Lagrange form of the Vainberg--Tonti
> Lagrangian of $\varepsilon $ and $\varepsilon $ itself: …`

A `grep -n -i "definition\|theorem"` confirms only one labelled
`\begin{definition}` and zero `\begin{theorem}` blocks in the body
(line 33 declares `\newtheorem{definition}…`, line 365 is the sole
use).

**Lines changed:** approximately lines 665–676 (new TeX line numbers
post-edit, from the future-work section on Path d).

**Before:**
```latex
  \item \textbf{Krupka--Voicu canonical
    completion}~\cite{KrupkaVoicu2015}: the Vainberg--Tonti homotopy
    formula provides a canonical Lagrangian for any non-variational EOM
    operator, with the Helmholtz residue measuring the obstruction.
```

**After:**
```latex
  \item \textbf{Krupka--Voicu canonical variational completion
    (Definition~1 of \cite{KrupkaVoicu2015})}: the Vainberg--Tonti homotopy
    formula provides a canonical Lagrangian for any non-variational EOM
    operator, with the Helmholtz residue measuring the obstruction.
    % [audit-2026-04-27] Per Meta-K verification of literature/1406.6646/main.tex:
    %   the paper has Definition 1 (Eq. 12), no Theorem 1.
    %   Cited object is the canonical variational completion as defined
    %   by KV; downstream "verifications" are consistency checks of this
    %   definition, not proofs of a theorem.
```

The TeX writeup did not contain any other "Krupka-Voicu Theorem 1"
mentions (`grep -n "Theorem"` returned only an unrelated `\paragraph{Definition.}`
on line 207, which is TIDAL's own definition of constraint promotion).
The "throughout" scope of the correction was therefore satisfied by
this single edit; if any future expansion of the document re-introduces
the "Theorem 1" misnomer, the comment block flags the correct
attribution.

---

## Correction 2 — Repair the BC 2018 Appendix D quote

**Audit justification (FINAL_ASSESSMENT.md item 2, with Meta-K K4 +
Meta-N §3-4 cross-corroboration):**
> "Repair the BC 2018 Appendix D quote: keep the verbatim quote but ADD
> the constructive follow-up sentence that was truncated: 'However,
> since the matrix F for b̄ = 0 is already diagonal, the critical
> parameters cn can be obtained directly from F. The same conclusion
> also holds for the form of H⊥ᶠ.' Reframe the surrounding paragraph:
> BC's Appendix D is constructive for *its* 2nd-order Dirac-Bergmann
> case but does NOT apply to TIDAL's higher-derivative Ostrogradsky case."

**Meta-N §3-4 quote (verbatim from Appendix D, lines 1678-1684 of the
PDF, lines 1884 of `literature/1804.05556/main.tex`):**
> "Now, we have a comment on kind of 'non-analiticity' of the above
> results. Since the assumption $\bar{b} \neq 0$ ensures the regularity
> of the matrix $P$, the diagonal matrix $D$ in (D.2) has no valid
> limit for $\bar{b} \to 0$. Hence, the expressions for $c_n$ when
> $\bar{b} = 0$ cannot be obtained by taking the limit $\bar{b} \to 0$
> of the generic result. **However, since the matrix $F$ for
> $\bar{b} = 0$ is already diagonal, the critical parameters $c_n$ can
> be obtained directly from $F$. The same conclusion also holds for
> the form of $\mathcal{H}^F_\perp$.**"

**Local-paper verification (`literature/1804.05556/main.tex`):**
- Line 658 (within Sec. 3.5 body): "The generic set of the critical
  parameters $c_\pm(F)$, $F=A,B_0,B_1,B_2$, is defined provided the
  parity odd parameters in $F$ do not vanish, see Appendix \ref{appD}.
  Hence, the limit of the final expressions $c_\pm(F)$ when these
  parameters tend to zero is not well defined. However, since in that
  case $F$ is already diagonal, one can identify $c_\pm$ directly from
  $F$." — confirms BC's prose explicitly emphasises the constructive
  workaround in the body.
- Line 1884 (within Appendix D): the disputed quote, verbatim, with
  the constructive follow-up. Both halves are present.

**Meta-N §5 framing (verbatim):**
> "BC 2018 Appendix D is a constructive recipe, but for a different
> class of theories (2nd-order quadratic PGT with parity-odd 2×2
> mixing). The recipe does not transfer to the b5·R̃²-induced
> Ostrogradsky rank-jump that is the actual TIDAL blocker."

Three structural reasons (Meta-N §5, paraphrased): wrong derivative
order (BC 2nd-order Dirac-ADM, TIDAL 4th-order Ostrogradsky); wrong
constraint topology (BC's $c_n \to 0$ adds a primary constraint;
TIDAL's $b_5 \neq 0$ removes algebraic constraints); BC's $\bar{b}$ is
not TIDAL's $b_5$ (different objects, naming coincidence).

**Lines changed:** lines 117–145 (new) replacing lines 117–129 (old),
plus lines 922–937 in the conclusion section.

**Before (history section, lines 117-129):**
```latex
\begin{quote}
``the diagonal matrix $D$ in (D.2) has no valid limit for $\bar{b} \to 0$.
Hence, the expressions for $c_n$ when $\bar{b} = 0$ cannot be obtained by
taking the limit $\bar{b} \to 0$ of the generic result.''
\end{quote}

This is a published, peer-reviewed acknowledgement that the perturbative
Hamiltonian across the $b_5 = 0$ critical surface \emph{does not exist}
for generic parameter points.  TIDAL's LPS \texttt{Throw} on detection of
constraint promotion is therefore not just architecturally honest but
mathematically necessary: the smooth perturbative-Hamiltonian decomposition
that LPS would deliver does not exist as a smooth function of $b_5$, and
no algorithm can produce one.
```

**After (history section):**
```latex
\begin{quote}
``\dots the diagonal matrix $D$ in (D.2) has no valid limit for
$\bar{b} \to 0$.  Hence, the expressions for $c_n$ when $\bar{b} = 0$
cannot be obtained by taking the limit $\bar{b} \to 0$ of the generic
result.  However, since the matrix $F$ for $\bar{b} = 0$ is already
diagonal, the critical parameters $c_n$ can be obtained directly from
$F$.  The same conclusion also holds for the form of $\mathcal{H}^F_\perp$.''
\end{quote}

The BC quote is constructive within its own scope: although the
parametric diagonalisation~(D.2) is singular at $\bar{b} \to 0$, the
critical parameters and $\mathcal{H}^F_\perp$ are recovered directly from
$F$ at $\bar{b} = 0$, where $F$ is already diagonal.  However, BC
2018's framework is built on a 2nd-order Dirac--ADM Legendre transform,
and their "critical parameters" are eigenvalues of small (2×2) mixing
matrices in a single irreducible spin sector. … [three structural
reasons enumerated in TeX] … BC 2018 Appendix~D is a
\emph{neighbouring-but-distinct} case study, not a published recipe
for the higher-derivative-induced rank-jump that constitutes TIDAL's
blocker. The accurate historical statement is: no published Hamiltonian
recipe is known for TIDAL's higher-derivative $b_5\,\tilde{R}^2$ case
across the $b_5 = 0$ critical surface; BC 2018 addresses the
parity-violating 2×2-mixing singular limit constructively but does not
address the Ostrogradsky rank-jump.
```

The conclusion section was also revised (lines 922-937, conclusion
paragraph 2) to drop the over-stated "BC 2018 says no algorithm can
produce one" framing in favour of the more accurate "BC 2018 is
constructive in its own scope but does not apply to TIDAL's
higher-derivative case".

---

## Correction 3 — Verify or correct the "Cabo Bizet-Bartocci" attribution

**Audit justification (FINAL_ASSESSMENT.md item 3, with Meta-K K1
correction):**
> "Verify or remove the Cabo Bizet-Bartocci citation if it appears
> anywhere. Per Meta-K, the original notes cite arXiv:2601.22007
> (Aashish-Saif) and arXiv:2602.12114 (Chan-López et al., NOT
> Cabo-Bizet-Bartocci). If the TeX writeup has 'Cabo Bizet-Bartocci',
> correct to actual authors; otherwise no change needed."

**Found it.** The TeX cites `\cite{CaboBizetBartocci2026}` at
TeX lines 389 and 578. The `references.bib` entry at line 1215 listed:
> author = {Cabo Bizet, A. and Bartocci, G.}

**Verification via Meta-K K1 (verbatim from `meta_review_K_literature_claims.md`):**
> "WebFetch on `https://arxiv.org/abs/2602.12114` (verified
> 2026-04-27): … Authors: E. Chan-López, A. Martín-Ruiz, Jaime Manuel
> Cabrera, Jorge Mauricio Paulin Fuentes."

So the bib entry's author attribution was wrong (the original
investigation appears to have transcribed the authors incorrectly when
making the bib entry, or the entry was pre-populated by an unverified
LLM step).

**Lines changed:** `references.bib` lines 1215-1224.

**Before:**
```bib
@article{CaboBizetBartocci2026,
  author  = {Cabo Bizet, A. and Bartocci, G.},
  title   = {Matrix bordering structure of the {F}addeev--{J}ackiw algorithm},
  journal = {arXiv preprint},
  year    = {2026},
  eprint  = {2602.12114},
  archiveprefix = {arXiv},
  primaryclass = {hep-th},
  note    = {Modern formalisation of FJ via Schur-complement bordering. …}
}
```

**After:**
```bib
% [audit-2026-04-27] Correct authors for arXiv:2602.12114 per Meta-K
% verification (https://arxiv.org/abs/2602.12114).  Previous entry
% mis-attributed the paper to ``Cabo Bizet \& Bartocci''; the actual
% authors are Chan-López, Martín-Ruiz, Cabrera, and Paulin Fuentes.
% The citation key is preserved (CaboBizetBartocci2026) to avoid
% touching the TeX call sites; the human-readable author field is
% corrected.
@article{CaboBizetBartocci2026,
  author  = {Chan-L\'opez, E. and Mart\'in-Ruiz, A. and Cabrera, J. M. and Paulin Fuentes, J. M.},
  title   = {Matrix bordering structure of the {F}addeev--{J}ackiw algorithm: {S}chur complement regularization and symbolic automation},
  journal = {arXiv preprint},
  year    = {2026},
  eprint  = {2602.12114},
  archiveprefix = {arXiv},
  primaryclass = {hep-th},
  note    = {Modern formalisation of FJ via Schur-complement bordering. Termination of the algorithm requires Schur-complement non-degeneracy --- exactly what fails on the constraint-promotion critical surface. Preprint, not yet peer-reviewed.}
}
```

The citation key `CaboBizetBartocci2026` was preserved to avoid
mass-rewriting TeX call sites; the bib entry now contains the correct
author list. The preprint status (not peer-reviewed as of WebFetch
retrieval, per Review 2 §10) is also flagged in the `note` field.

---

## Correction 4 — Downgrade Path B-trace (Barker)

**Audit justification (FINAL_ASSESSMENT.md item 4):**
> "Downgrade Path B-trace from 'established literature (Barker 2024)'
> to 'open research conditional on parity-odd extension of Barker
> 2024'. Justification: Barker 2024 explicitly EXCLUDES parity-odd
> terms (line 109 footnote in `literature/2406.12826/`). Trace torsion
> is identified as a Yang-Mills GAUGE FIELD, not a Goldstone (line
> 168). Verify these via the local TeX file."

**Local-paper verification (`literature/2406.12826/Manuscript.tex`):**
- Line 109 (footnote, verbatim): "In this letter we omit the
  parity-odd invariants only out of simplicity; there are no
  convincing theoretical grounds for excluding them"
- Line 168 (verbatim): "showing that eWGT is the unique scale-invariant
  embedding of PGT. This will identify $T_\mu/3$ with the vector
  $B_\mu$ when expressed in scale-invariant variables, and thereby
  reveal $\partial_{[\mu} T_{\nu]} \partial^{[\mu} T^{\nu]}$ to be a
  Yang-Mills-type term."
- Line 297 (verbatim): "The compensator is purely gauge, so that the
  embedding theory is completely indistinguishable from PGT after
  gauge-fixing."
- Line 300 (verbatim, closing remarks): "Finally, the Yang-Mills-type
  actions in [PGTAction,WGTAction] are restricted to parity-even terms
  for simplicity: the parity-odd extensions should be considered."

The original TeX file contained no explicit "Path B-trace" claim, but
the FINAL_ASSESSMENT correction directs the documentation to record
that any such future construction is conditional on a parity-odd
extension. I added a new paragraph
("Sectoral Stückelberg for trace torsion (open: parity-odd extension)")
in the methods-survey section after the existing Stückelberg paragraph
(approximately line 489 onwards in the new TeX).

**Lines changed:** Inserted a new paragraph at approximately line 489,
within Section "Survey of perturbative-reduction methods…", after the
Stückelberg-lifting paragraph.

**Inserted text (after):**
```latex
\paragraph{Sectoral Stückelberg for trace torsion (open: parity-odd extension).}
A separate constructive lead --- treating the trace torsion
$T_\mu/3$ as a Stückelberg or Goldstone field for the conformal embedding
of PGT in extended Weyl gauge theory --- has been proposed on the basis
of Barker~\cite{barker2024poincare}.  Two facts in the published paper
restrict the applicability for TIDAL's $b_5\,\tilde{R}^2$ case:
(i) Barker identifies $T_\mu/3$ with the Yang--Mills-type Weyl gauge
field …, not with the Goldstone of broken Weyl symmetry; the Goldstone
(compensator) is the \emph{scalar} $\phi$ (verified:
\texttt{literature/2406.12826/Manuscript.tex} line 168, line 297).
(ii) The construction explicitly omits parity-odd invariants ``out of
simplicity; there are no convincing theoretical grounds for excluding
them''~\cite[footnote at line 109]{barker2024poincare}.
TIDAL's $b_5\,\tilde{R}^2$ is parity-odd by construction, so applying
this trace-sector recipe to TIDAL's case requires a parity-odd
extension of Barker that has not been published. Path B-trace is
therefore \emph{open research}, conditional on such an extension ---
not established literature.
```

---

## Correction 5 — Downgrade Path B-tensor-q (Chatzistavrakidis-Ranjbar-Zekoč)

**Audit justification (FINAL_ASSESSMENT.md item 5):**
> "Downgrade Path B-tensor-q from 'applies with caveats' to
> 'conditional on parity-odd extension of Chatzistavrakidis-Ranjbar-Zekoč
> 2024 that does not exist in published literature'."

**Local-paper verification (`literature/2411.16928/Stueckelberg_v2.tex`):**
- Line 49 (abstract, verbatim): "We develop a Stueckelberg mechanism
  for mixed symmetry tensor fields at the linearized level, focusing
  on the massive graviton, the massive (2,1) Curtright field and the
  massive (2,2) field." — confirms linearised, free fields.
- Line 72 (verbatim): "When the mass approaches zero, the fields
  decouple, and the action for the Stueckelberg field represents the
  action for Goldstone modes." — confirms `m → 0` is the Goldstone
  limit (not a constraint-promotion limit).
- The paper text contains no parity-odd terms; all kinetic structures
  are standard Fierz-Pauli / standard (2,1) Curtright with no Pontryagin
  / Chern-Simons-type couplings.

**Bib entry added** (`references.bib`, after Lyakhovich2021):
```bib
@article{ChatzistavrakidisRanjbarZekoc2024,
  author  = {Chatzistavrakidis, A. and Ranjbar, A. and Zeko\v{c}, T.},
  title   = {Higher tensor global symmetries, mixed symmetry fields and the {S}tueckelberg mechanism},
  …
  eprint  = {2411.16928}, …
}
```

**Lines changed:** A new paragraph
("Sectoral Stückelberg for tensor torsion (open: parity-odd extension)")
inserted in the methods-survey section, immediately following the
Path B-trace paragraph from Correction 4.

**Inserted text:**
```latex
\paragraph{Sectoral Stückelberg for tensor torsion (open: parity-odd extension).}
A second constructive lead --- the Stückelberg construction for massive
mixed-symmetry tensor fields (Curtright (2,1) field) in
Chatzistavrakidis, Ranjbar \& Zekoč~\cite{ChatzistavrakidisRanjbarZekoc2024}
--- has been proposed for the tensor-q torsion sector. The published
paper handles only \emph{parity-even free-field} Lagrangians at
linearised order, and the $m \to 0$ limit it analyses is the Goldstone
(symmetry-breaking) limit rather than a constraint-promotion limit
(verified: \texttt{literature/2411.16928/Stueckelberg\_v2.tex},
abstract + Sections 5.1--5.2). Translating the construction to the
parity-odd, non-free $b_5\,\tilde{R}^2$ case requires a parity-odd
extension that does not exist in the published literature. Path
B-tensor-q is therefore likewise \emph{conditional} on an unpublished
parity-odd extension of Chatzistavrakidis-Ranjbar-Zekoč 2024 --- not a
verified construction for TIDAL's actual Lagrangian.
```

---

## Correction 6 — Voicu linearity-in-highest-derivatives gate as separate preflight

**Audit justification (FINAL_ASSESSMENT.md item 6):**
> "Add Voicu linearity-in-highest-derivatives gate as a separate
> preflight, distinct from VT integral convergence. Justification:
> Voicu 2020 (`literature/2009.05459/varcompl_*.tex`) §4 + Appendix A
> identifies two independent failure modes; the agents' 'VT
> convergence' framing bundled them."

**Local-paper verification
(`literature/2009.05459/varcompl_EPJP_format_Rev1_arxFin.tex`):**
- Line 242 (verbatim, from §3.2): "the resulting improper integral
  diverges. This is, e.g., the case of PDE systems … that are
  homogeneous of negative degree smaller or equal to $-1$."
  → defines the homogeneity-degree-$\leq -1$ failure mode for VT
  integral convergence.
- Line 291 (verbatim, from §4): "If a second order PDE system is
  locally variational, then it must be linear in the second order
  derivatives of the dependent variable. Or, the truncated Gauss-Bonnet
  terms $A_{\mu\nu}$ are not linear in the second order derivatives of
  $g_{\mu\nu}$, see Appendix \ref{App:A}, in any dimension."
  → defines the **separate** linearity-in-highest-derivative failure
  mode.
- Line 442 (verbatim, from same Appendix A): "But, any variational PDE
  system which is of second order must be linear in the second order
  derivatives acting on the fundamental dynamical variable [Krupka-book
  p.~147]. Hence the truncated field equations cannot be variational
  and the variation of $A^\mu{}_\mu$ cannot be of second order only,
  but must contain higher derivatives."
  → confirms the linearity gate is independent of and additional to
  the homogeneity gate.

**Lines changed:** Replaced the single-gate "Falsifiable preflight test"
paragraph with a three-gate (P1)/(P2)/(P3) version in the future-work
section (approximately new lines 696–740).

**Before:**
```latex
\textbf{Falsifiable preflight test.}  Compute the Helmholtz
self-adjointness residue $\delta E$ on the Parker--Simon-reduced
order-1 EOM operator …
```

**After:**
```latex
\textbf{Falsifiable preflight tests (three independent gates).}
Direction~(d) is viable for a given PGT theory only if all three of
the following gates hold for the Parker--Simon-reduced order-1 EOM
operator:
\begin{enumerate}
  \item[(P1)] \textbf{Helmholtz self-adjointness}: $\delta E = 0$. …
  \item[(P2)] \textbf{Vainberg--Tonti integral convergence}: …
        (Voicu~\cite[Sec.~3.2]{Voicu2020}; verified in
        \texttt{literature/2009.05459/…} line 242). …
  \item[(P3)] \textbf{Linearity in highest-order derivatives}
        (Voicu~\cite[Sec.~4 + App.~A]{Voicu2020}; verified at the same
        source, lines 291 and 442). … this is a \emph{distinct}
        failure mode from VT integral convergence. …
\end{enumerate}
```

---

## Correction 7 — Add the M_c² → 0 caveat to Path A's verdict

**Audit justification (FINAL_ASSESSMENT.md item 7):**
> "Add the M_c² → 0 caveat to Path A's verdict — L_VT diverges at
> constraint-mass critical surfaces. Justification: Review 1 C5
> sympy-verified in `reviews/scripts_review/`. PGT critical-mass
> surfaces (Karananas 2014, Blagojević 2018) are exactly where Path A
> breaks."

**Review 1 C5 verbatim:**
> "L_VT denominator: `2·M_1⁴·M_2⁴` … 16 terms with `1/M_1²` poles, 6
> terms with `1/M_1⁴` poles … `limit(L_VT, M_1 → 0) = ∞·sign(...)` —
> L_VT genuinely diverges as a constraint mass goes to zero. … For
> PGT, constraint masses M_c are functions of the Lagrangian coupling
> constants. At certain critical surfaces (notably ghost-free PGT
> critical cases — see Karananas 2014, Blagojević 2018), some
> M_c² → 0. Path A's L_VT is undefined on those surfaces."

**Lines changed:** A new "Constraint-mass critical-surface caveat (new
pathology)" paragraph immediately after the new (P1)/(P2)/(P3)
preflight block (approximately lines 742–758 of the new TeX).

**Inserted text (after the preflight gates):**
```latex
\textbf{Constraint-mass critical-surface caveat (new pathology).}
% [audit-2026-04-27] Review 1 C5 (sympy-verified) …
Even when the three preflights pass, the Vainberg--Tonti Lagrangian
$L_{\rm VT}$ for the constraint-promoted sector contains a Routhian
projector denominator proportional to a product of constraint masses
$\prod_c M_c^2$, so $L_{\rm VT}$ contains explicit $1/M_c^2$ and
$1/M_c^4$ poles (Review~1 C5,
\texttt{research/perturbative\_hamiltonian/reviews/scripts\_review/C5\_routhian\_M\_to\_zero.py}).
PGT critical-mass surfaces (where one or more constraint masses vanish,
see~\cite{BlagojevicCvetkovic2018} for the relevant surfaces in
quadratic PGT) are exactly the parameter points where direction~(d)'s
constructive output is undefined. Path~A's verdict is therefore
conditional on $M_c^2 \neq 0$ for every constraint-promoted field $h_c$
at the parameter point of interest; the most physically interesting
PGT critical cases lie on these surfaces, where direction~(d) fails.
```

I did not add a fresh `Karananas2014` bib entry, because the audit
note in `FINAL_ASSESSMENT.md` is paraphrastic and the in-text
reference to Karananas was already in a comment, not in `\cite{}`.
The actual citation in the body uses the existing
`BlagojevicCvetkovic2018` entry, which directly classifies critical
surfaces.

---

## Correction 8 — Tone down F-J cross-validation

**Audit justification (FINAL_ASSESSMENT.md item 8):**
> "Tone down the F-J cross-validation to 'qualitative consistency at
> linearised order; structurally different mechanisms' (Agent F's
> det(M_aux) = 1/b² DIVERGES at b=0 — rank-jump relocated, not
> removed; Agent J's det(H_kin) = 1 - λ_a² is genuinely b5-independent).
> Justification: Review 1 C8 sympy-verified."

**Status: NOT APPLIED — N/A.** The TeX writeup contains no F-J
cross-validation claim. `grep -n -i "cross-valid\|cross valid"
docs/tex/perturbative_reduction_constraint_barrier.tex` returns no
hits. Review 1 C8 audited the Round 3 *synthesis-note* claim, not a
TeX claim. No edit applied to the TeX file under this correction.

I leave this entry as a record so a future reader does not re-introduce
a "structurally parallel" claim into the TeX without revisiting C8.

---

## Correction 9 — Operational primary section for `tidal/measurement/_conversion.py`

**Audit justification (FINAL_ASSESSMENT.md item 9):**
> "Add an operational primary section documenting
> `tidal/measurement/_conversion.py` as the Hamiltonian-based
> conversion observable. The actual headline channel is `h_5 ↔ a_1`
> which is **standard-kinetic, NOT constraint-promoted**, so the
> constraint-promotion barrier does NOT block this measurement.
> Justification: Meta-L confirmed via direct code reading + JSON
> inspection."

**Meta-L verified file evidence (verbatim from
`meta_review_L_pipeline_claims.md`):**
- `tidal/measurement/_conversion.py` line 4-7 (module docstring):
  `P(t) = E_target(t) / E_source(0)`.
- `tidal/measurement/_conversion.py` line 92-95: "Energy is computed
  via the canonical Hamiltonian (kinetic + gradient + mass), including
  the spatial volume element `sqrt|g_spatial|` for curved coordinates."
- JSON inspection (`examples/data/torsion_gertsenshtein.json`):
  - `h_5`: `time_order=2, kinetic_coefficient_symbolic="-kappa^(-2)"`
    (standard graviton).
  - `h_4/h_7/h_9`: `time_order=4, kinetic_coefficient_symbolic="2*b5"`
    (constraint-promoted).
  - All sweep scripts in `examples/torsion_gertsenshtein/` use
    `--source h_5 --target a_1`.
- 276-run dark-photon HPC sweep matches Boccaletti baseline to
  6.7×10⁻⁶ (memory file
  `dark_photon_amplification_campaign_v0.31.md`).

**Independent verification by W1 (this agent):** Read
`tidal/measurement/_conversion.py` lines 1-100 directly. Confirmed:
docstring states `P(t) = E_target(t)/E_source(0)`;
`compute_conversion_probability` calls
`compute_energy_timeseries` (which evaluates
`canonical.hamiltonian_terms`); the result is energy-ratio, not
amplitude-ratio.

**Lines changed:** A new subsection
`\subsection{Operational primary: Hamiltonian-based conversion observable}`
inserted between "Unaffected examples" and "Future work" (approximately
new lines 614–693).

**Inserted text (excerpted; the full subsection is in the TeX file):**
```latex
\subsection{Operational primary: Hamiltonian-based conversion observable}
\label{sec:pr-cb-operational}

The constraint-promotion barrier described above is a real architectural
limit on the \emph{total system energy budget} for theories carrying
$b_5\,\tilde{R}^2$.  It does \emph{not}, however, block TIDAL's flagship
graviton--photon conversion measurement. …

\paragraph{The conversion observable.}
The function \texttt{compute\_conversion\_probability} in
\texttt{tidal/measurement/\_conversion.py} computes
$P(t) = E_{\rm target}(t)/E_{\rm source}(0)$, where the per-field
energies are evaluated via \texttt{compute\_energy\_timeseries}
(\texttt{tidal/measurement/\_energy.py}). The energy expression is the
canonical Hamiltonian (kinetic plus gradient plus mass) drawn from the
JSON \texttt{canonical.hamiltonian\_terms} block …

\paragraph{Why the constraint-promotion barrier does not block this.}
For \texttt{torsion\_gertsenshtein}, the Gertsenshtein-effect headline
sweep uses \texttt{--source h\_5 --target a\_1} … $h_5$ has
\texttt{time\_order = 2} … a standard-kinetic graviton component, \emph{not}
one of the $b_5$-promoted Pais--Uhlenbeck components $h_4, h_7, h_9$.
…

\paragraph{Quantitative validation in production.}
The 276-run dark-photon amplification campaign … reports agreement
with the analytical Boccaletti baseline $P_{\rm GR} =
\sin^2(\kappa B_0 t / 2)$ to a maximum residual of $6.7 \times 10^{-6}$
across the full sweep. …

\paragraph{Where the barrier still matters.}
For total system energy conservation diagnostics … the $h_4/h_7/h_9$
self-energy terms use \texttt{mixed\_2\_0\_0\_0} operators that
resolve via \texttt{\_compute\_acceleration\_from\_eom} (re-evaluating
the EOM at snapshot time). …
```

The Conclusion section was also revised (lines 894–944 of the new
TeX) to reflect that the constraint-promotion barrier matters for
total-system-energy diagnostics but does not block the operational
primary observable. The original conclusion's "Simulation alone is
not measurement-ready" statement was replaced by a more accurate
asymmetric framing: the Hamiltonian-side LPS is necessary for
phase-space-clean total-system observables (blocked) but the headline
$h_5 \leftrightarrow a_1$ conversion observable is sector-clean and
delivers measurement-grade $P(t)$ in production.

---

## Build verification

Ran `cd docs/tex && make` after all edits.

```
This is BibTeX, Version 0.99d (TeX Live 2022/Debian)
The top-level auxiliary file: main.aux
The style file: unsrt.bst
Database file #1: references.bib
Warning--empty journal in dandoy2024constraining
Warning--empty journal in rubilar2003torsion
…
(There were 9 warnings)
Build complete: main.pdf (1950100 bytes, 255 pages)
```

- LaTeX errors: **0**
- Undefined references: **0** (`grep -E "Reference.*undefined" main.log`
  returned nothing).
- Undefined citations: **0** (`grep -E "Citation.*undefined" main.log`
  returned nothing).
- Multiply-defined labels: **0**.
- The 9 BibTeX warnings are pre-existing (`Warning--empty journal`
  for `dandoy2024constraining`, `rubilar2003torsion`,
  `bahamonde2025cubic`, `obukhov2024electrodynamics`,
  `barker2024poincare`, `domcke2025scattering`, `martingarcia2007xact`,
  `FKY2025`, `Dunne2004`); none are caused by this edit pass.
- The new `ChatzistavrakidisRanjbarZekoc2024` bib entry is correctly
  resolved (`grep -c "ChatzistavrakidisRanjbarZekoc2024" main.bbl` = 1).
- The new `\label{sec:pr-cb-operational}` resolves correctly via the
  `\Cref{sec:pr-cb-operational}` reference in the conclusion.

---

## What was NOT done (and why)

- **Correction 8 (F-J cross-validation tone-down)**: not applicable
  to the TeX writeup. The TeX never claimed F-J as parallel mechanisms;
  the audit target was the Round 3 synthesis note, not this TeX file.
  Preserved as a future caveat in this log so a hypothetical future
  expansion does not re-introduce the overstated framing.

- **No new physics claims introduced.** Every edit is a restatement,
  reframing, or downgrade backed by a quoted line from a local TeX
  literature file or a referenced audit document. No "in our view",
  no "expert opinion", no speculative interpolations.

- **No paper text was paraphrased.** The BC 2018 quote is reproduced
  verbatim with the constructive follow-up sentence restored.

- **No content was deleted that the audit did not flag.** The
  retained sections (vDVZ analogue, Boulware-Deser, DHOST, cubic-PGT,
  Stelle, Lyakhovich path (e), HiGGS validation) are unchanged from
  pre-audit.

- **The Karananas 2014 attribution** in the new constraint-mass-caveat
  paragraph is left as a comment (not a `\cite`) because the audit
  note paraphrased it ("PGT critical-mass surfaces (Karananas 2014,
  Blagojević 2018)") without specifying which Karananas paper. The
  in-text `\cite` is `BlagojevicCvetkovic2018`, which is the existing
  bib entry that classifies critical surfaces. If a future expansion
  needs to cite Karananas 2014 directly (likely arXiv:1406.7456 on
  parity-violating PGT spectrum), add a fresh bib entry then.

---

## Summary

All 9 corrections from `FINAL_ASSESSMENT.md` § "Concrete documentation
corrections" were processed:

| # | Status | Lines touched |
|---|--------|---------------|
| 1 | Applied | 1 site (KV definition clarification) |
| 2 | Applied | 2 sites (history quote + conclusion) |
| 3 | Applied | bib entry only (TeX call sites unchanged) |
| 4 | Applied | new paragraph in methods survey |
| 5 | Applied | new paragraph in methods survey + new bib entry |
| 6 | Applied | replaced single preflight with three-gate version |
| 7 | Applied | new caveat paragraph after preflights |
| 8 | N/A — TeX did not contain the audit-flagged claim |
| 9 | Applied | new subsection + revised conclusion |

Build is clean. Audit closed.
