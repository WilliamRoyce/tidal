> **ARCHIVE — verbatim record of the 2026-09-15…16 R-C (#567) planning session plan file.**
>
> Preserved here because it lived only in `~/.claude/plans/` (not version-controlled; its
> backup directory is gitignored) and carries the full reasoning trail: every research pass,
> a pinned evidence ledger (Appendix B), what could not be settled by reading (Appendix C), and
> four rounds of user review that reshaped the plan (Appendix A: eleven points, six comments,
> seven comments, twenty comments).
>
> **Not authoritative.** The living documents supersede it wherever they differ:
> `docs/cosmology/perturbation_tooling.md` (the R-C memo — requirements, methods, tools,
> recommendation and its run ledger) and `docs/COSMOLOGY_PROGRAM.md` (decisions register and
> wave board — D-C is recorded there by the orchestrator, not here). Everything below the rule
> is the plan file exactly as approved, after one post-approval lint pass on the plan itself
> (spelling words added to its `cspell:words` line, three British-spelled or accented tokens
> rewritten, and one sentence that named the forbidden path patterns literally reworded — no
> content change; logged in Appendix A). No redaction was made in the archive. "You" in the text addresses the user
> who dispatched R-C.

---

# R-C (#567) — FRW perturbation requirements, methods, tools: plan

<!-- cspell:words xPand xPert xTras xCoba xTensor xPerm xCore xMAG xIST COPPER HiGGS Hamilcar xPPN FieldsX bimEX xCPS xTerior TexAct xBrauer TraceFree Pitrou Umeh Brizuela Nutma Bahamonde Gigante Valcarcel Gorji Hohmann Heisenberg Kuhn Golovnev Koivisto Nikiforova Damour Chee Toporensky Tretyakov Obukhov Bertschinger Sletmoen SymBoltz Cadabra Zumalacárregui Bellini Sawicki EFTCAMB Gubitosi Gleyzes Langlois Piazza Vernizzi Noller Helpin PSALTer Barker Cembranos Seljak Zaldarriaga Challinor Lasenby Lewis Nicosia Gakis Kiorpelidi Saridakis Mukhanov WXF wolframscript ToxPand SetSlicing DefMetricFields DefMatterFields DefProjectedTensor DefProjectedTensorProperties SplitPerturbations SplitMetric SplitMatter ExtractComponents ExtractOrder ToxPandFromRules SplitGaugeChange DefConformalMetric ConformalWeight VisualizeTensor NewtonGauge SynchronousGauge AnyGauge FLFlat FLCurved DefTensorPerturbation DefMetricPerturbation ExpandPerturbation TorsionCDT RicciScalarCDT ChangeCurvature ChristoffelCDCDT VarD LieD contortion disformation Weitzenböck teleparallel tadpole eikonal Codazzi Bardeen etak hdot cspell pyright ruff normu InducedFrom SignDetOfMetric DummyS DummyV DummyT ToDistortion BreakDistortion DefConnectionPerturbation StartInducedDecomposition pdftotext sympy Duhamel Einstein–Boltzmann adotoa dcaabbea Dicke Fröb BRST recid hersle gentelepert vierbein apstemplate Pfeifer Randjbar Daemi Rubakov Andrade regenerable ONSHELL dgrho syna epsilong refersto xenos Arroja Bartolo Domènech EPJC Feldman Fidler Inomata Pajer Pandcitations Skordis Złośnik cdpost cdpre cona inds -->

> **Plan file for the delegate session that executes R-C.** Read `docs/cosmology/handoffs/R-C.md`
> (the brief) first. This plan is archived verbatim into the repo at the first commit as
> `docs/cosmology/rc_planning_record.md` (R-1 precedent; decided by the user, §8 Q1). Third
> submission: round 1 returned eleven points (Appendix A), round 2 six comments on what the
> lane records, on the physics framing of the residual, and on readability — all applied.

## TL;DR

**Problem (one line).** M3 (the milestone that builds the FRW derivation) must derive, for a
user theory, the linear perturbation equations of its new sector on an FRW background in a
CAMB-named gauge, plus source functions, the background-consistency gate, and the equations in
the forms the solver candidates need (no solver is chosen here — see "Solver forms" below).
Nobody has surveyed what the literature does, which packages exist, and whether the obvious one
(xPand) can carry a torsionful sector.

**Design (plain terms).** This lane produces **a decision memo and a complete record of how it
was reached**, so the next session inherits evidence, not conclusions:

1. `docs/cosmology/perturbation_tooling.md` (≤ ~600 lines, tables over prose) in four parts —
   **Part 1** the equation set as a sourced checklist, each line naming who consumes it;
   **Part 2** how the literature derives such equations; **Part 3** every tool with scope,
   uses, license, release, compatibility, checklist coverage and a torsion verdict; **Part 4** a
   build-on / borrow / build-own recommendation per requirement with cost, license cost, the
   D-A interface it implies, and a probe → outcome → option decision table. It also carries, in
   R-1's format, a "what was verified here and what was not" section, a run ledger and an
   evidence index.
2. `docs/cosmology/rc_planning_record.md` — this plan verbatim, with the revision log, the
   evidence ledger and the not-verified list (the reasoning trail).
3. `scripts/research/perturbations/` — every script, target, printout and installer with a
   provenance table and a "How to reproduce" block; run artifacts in gitignored `third_party/`.
4. `docs/references.md` rows for every paper the memo relies on; GitHub issues for every
   discovery that outlives the memo (as R-1 filed #569–#573).

The evidence comes from **scripts, not sentences**: xPand 0.4.4 and xMAG (with four
dependencies) installed additively into the real userbase; the FRW tensor equation from the
Einstein–Hilbert action at second order, **with the order-1 (tadpole) coefficient computed
and evaluated on the background**, as the headline (nothing is deleted from any equation: the
ε¹ coefficient of the expanded action is printed and checked to vanish once the background
equations hold — that check is what the pipeline will export as the gate); the same equation
from the field equations as the control; the Ma & Bertschinger
scalar equations in both CAMB-named gauges (conformal-Newtonian and synchronous),
machine-checked against the installed `camb.symbolic`;
a signature probe with a source-grounded hypothesis; a torsion probe designed from the papers
that did it (a homogeneous background mode only as an optional, out-of-scope check); and the
torsion sector through
the action route. The harness verdicts are `identical` / `proved-equal` / `proved-different` /
`could-not-decide` (the last counts as **not reproduced**), with negative controls.

**Next step.** Worktree `/tmp/tidal-rc` on `cosmo/rc-perturbation-tooling` → lane gates
(before) + userbase manifest → first commit (plan archive, README, installers, harness,
targets) → draft PR → install xPand → load probe → A2 control → B (MB + CAMB check) → A1
headline → C → D (torsion, incl. background torsion and the action route) → xMAG install +
E (half-day box) → gates (after) + manifest → memo → `docs/references.md` rows → issues → CI
green → report.

**DONE means.** All seven success criteria of `R-C.md` hold from artifacts: Part 1 has a source
and a consumer per line; Part 2 cites paper *and* code per method; Part 3 covers every tool
named in Part 2 and every xAct-family package; both reproductions reach symbolic equality with
a committed negative control; the torsion probe has a definite verdict with its script (xPand
route, action route and, under Q4, the xMAG route); `verify --require-psalter` exit 0 before and
after, transcripts attached (scrubbed); Part 4 recommends per requirement with cost and a
decision table; draft PR from the first commit, `CI <run-id>: <conclusion>`, lane ≤ 2 days;
the userbase gained exactly the six directories of §5.3 and §5.10 and the manifest is identical
outside them and the allow-list of §8.1; the planning record, the README, the references rows
and the issues exist.

## Terms (one sentence of meaning, one of why it matters here)

- **Wave 1 · R-C · M3 · WS2** — Wave 1 is the current batch of work (research memos first);
  R-C is this lane; **M3** is the later milestone in which the FRW derivation is actually
  built, by workstream **WS2** (symbolic derivation). This lane decides *what M3 builds on*; it
  builds no production code.
- **Part 1 / 2 / 3 / 4** — the four sections of the memo this lane writes (`R-C.md`'s
  structure): requirements, methods, tools, recommendation. "Part 4" always means the
  recommendation section of *this* memo, written *in this lane*.
- **D-C** — the decision the memo feeds: which tooling M3's FRW derivation is built on
  (`COSMOLOGY_PROGRAM.md:1007`). The orchestrator records it after this lane merges, as the
  register row "FRW derivation tooling" (`:948-950`); the memo puts the options and their
  costs to the user, it does not decide.
- **Standard sectors** — the ΛCDM species CAMB already evolves: photons, baryons, cold dark
  matter, neutrinos, dark energy, and the metric. The user's Lagrangian contains the new sector
  **and its couplings to the standard fields it touches** — the metric (new gravitational
  interactions) and the photon (modified photon interactions) — and those touched modes are
  evolved inside our block, replacing CAMB's evolution for that channel
  (`repo_reshape.md:233-237`; `R-C.md:43-48`). The sectors it does not touch (baryons, CDM,
  neutrinos, dark energy) reach our equations only as background functions (`a(η)`, `ℋ(η)`,
  `n_e(η)`, …) and as the source terms CAMB supplies (`S_std`).
- **Spectator** — a new sector whose energy density is negligible in the Friedmann equation,
  so CAMB's ΛCDM background `a(η)` is unchanged and only the perturbations on top of it are
  new (`spectator_route.md:54-62`). This is the program's scope: theories that would alter the
  expansion history are out (`COSMOLOGY_PROGRAM.md:200`, "Admissible theories"). The exact
  solution of the full theory need not coincide with ΛCDM; what is required is that the
  difference is negligible to a **derivable tolerance** — the induced source far below the
  signal being computed (`spectator_route.md:107-112`) — which is what the gate below measures,
  and the admissible region of couplings is where that holds. A small nonvanishing background
  torsion that stays below the tolerance is therefore acceptable; one that does not marks the
  theory, or that coupling point, as out of scope.
- **Tadpole / background-consistency gate (#501)** — when the action is expanded in the
  perturbation size ε, the ε¹ coefficient is the theory's own background field equation
  evaluated on the ΛCDM background. "Computing the tadpole" means printing that coefficient
  and evaluating it with the background equations substituted — not removing a term from any
  equation. For an admissible spectator theory it vanishes (or is negligible to a derivable
  tolerance); the pipeline **exports it as a per-theory gate** (`observable_ladder.md:185-193`,
  "a per-theory gate, computed as a background-EOM residual, not an assumption") and never
  propagates it as a source. The reproductions compute it to show the machinery exists.
- **Homogeneous new-sector mode** — FRW *symmetry* permits a torsion background with two
  functions of time (a vector and an axial piece, §2.6). A rung that would *need* such a mode
  (isotropic cosmic birefringence, O4-iso) is **outside the planned scope**: the register keeps
  it as a user-theory question under the gate, not a planned rung (`COSMOLOGY_PROGRAM.md:211`),
  and this lane does not investigate it. What stays in scope is the gate's job of detecting
  when a theory's FRW solution wants such a mode beyond tolerance — the residual on the
  `T̄ = 0` background does not vanish — and marking that theory inadmissible; computing the
  residual needs no representation of the mode itself. Probe D5 (representing the mode) is
  therefore optional and time-boxed, run last if at all, so that Part 4 can state a cost should
  the scope ever widen.
- **Post-Riemannian rewrite** — writing every curvature and covariant derivative of the
  torsionful connection as Levi-Civita quantities plus contortion, then contortion as torsion
  (`K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c − T_{bc}{}^a)`). **Both branches do it**: the
  spectrum branch because PSALTer's derivative is the flat `∂` (`spectrum_design.md:341-351`),
  the FRW branch because xPand's machinery is Levi-Civita-only — which is exactly where a
  package built for independent connections (xMAG) can assist (option O3), and, because xMAG
  also carries non-metricity, the same route would extend to metric-affine perturbations later
  without new machinery. Legacy TIDAL implemented it by
  hand (`_derive.py::_wls_torsion_curvature_decomposition`, `:2164-2189`:
  `ChangeCurvature[L, CDT, CD]` plus the contortion identity, citing Shapiro 2002 eq. 2.9 and
  Hehl et al. 1976);
  Stage 1 ports that (`stage1_engineering_plan.md:519-525`); xMAG's `ToDistortion` +
  `BreakDistortion` package the same step (plus non-metricity). Our field content in both
  branches is metric + torsion (+ vector), the Einstein–Cartan form
  (`spectrum_design.md:353-356`), so the FRW variables are ordinary tensors xPand can carry.
- **Solver-only block (`interfaces_decision.md` §3.8)** — the user's theory file is shared by
  the two derivations (`spectrum_design.md:372-378`: shared = the input model, field
  declarations, coupling roster, term structure; duplicated in each branch = the
  convention-laden symbolic work such as the post-Riemannian rewrite). What the spectrum branch
  never needs — the FRW background (`a(η)`, background field values, a homogeneous mode), how
  each field splits into background + perturbation *on FRW* (on Minkowski the split is
  trivial), and the gauge choice (PSALTer takes the gauge-unfixed Lagrangian and handles gauge
  itself, `:389-396`) — would change the spectrum's content-address (fingerprint) without
  changing the spectrum if it lived in the shared file, so it sits in a separate block with its
  own fingerprint. Part 1 lists what that block must carry; designing it is M3's job.
- **FRW / FLRW** — the homogeneous isotropic expanding background `ds² = a²(η)(dη² − dx²)` in
  our signature. The memo says **FRW** (every design doc does; only `R-C.md` says FLRW).
- **Conformal time η, scale factor a(η), conformal Hubble rate ℋ = a′/a** — CAMB's time
  variable (`t` in `camb.symbolic`); our coefficients keep `a(η)`, `ℋ(η)` symbolic until CAMB's
  table fills them at solve time.
- **SVT decomposition; helicity** — splitting spatial perturbations into scalars, transverse
  vectors and transverse-traceless (TT) tensors, which decouple at linear order on FRW;
  helicity is the Fourier-space label of the same split (0, ±1, ±2; ±3 only for rank-3
  symmetric pieces). The observable channels (tensor → B-modes) are helicity sectors.
- **Eikonal export (#504)** — the second-order wave system rewritten for slowly varying
  amplitudes after factoring out the fast carrier `e^{−iωη}`; CMB photons oscillate ~10²⁶
  times over cosmic history and cannot be integrated directly.
- **Solver forms (what the derivation must emit, without choosing a solver)** — the solver
  research (H3, `solver_design.md` §5 and §7) settled that no candidate is discounted on
  paper: Magnus and matrix-WKB are built together with `rk-adaptive` as the measured baseline
  and a bake-off decides; the O3 photon channel additionally needs the eikonal amplitude form.
  Every candidate consumes a first-order system `y′ = M(η) y (+ s)` with `M` assembled from
  symbolic coefficient expressions on CAMB's η-grid. So the requirement on the derivation is to
  emit **the representations the candidates need** — the second-order coupled block, its
  first-order state form (e.g. `y = (h, h′/k, δT, δT′/k)`), and for the photon channel the
  eikonal amplitude system — as separate sections of one export (`solver_design.md:466-470`
  allows several derived representations per derivation). The memo records the options and
  what form each demands; it commits to none.
- **Distortion, contortion, non-metricity** — the difference between a general affine
  connection and the Levi-Civita one is the distortion tensor; its metric-compatible,
  torsion-carrying part is the contortion (an algebraic function of the torsion tensor); the
  rest comes from non-metricity `∇g ≠ 0` (absent in PGT). Writing `Γ̃ = Γ(LC) + K` lets a
  Levi-Civita-only tool carry torsion as an extra tensor field.
- **`LI[o]` perturbation labels; Lie-derivative labels** — xPert tags a perturbation's order
  with a label index `LI[o]`; xPand adds a second label counting time derivatives, so
  `X[LI[1],LI[2]]` means `X″` at first order. Varying an action after the split must
  integrate these labels by parts in η.
- **`identical` / `proved-equal` / `proved-different` / `could-not-decide`** — the harness
  verdicts: literally the same expression; a difference that canonicalizes to zero; one that
  does not; a timeout. Only the first two count as reproduced.
- **Negative control** — the same comparison with one term of the target deliberately
  altered; it must come out `proved-different`, proving the check can fail.
- **Sentinel line** — a `RC_<STEP>_<KEY>=<value>` line a script prints; success is read from
  sentinels and artifacts, never from a kernel's exit status (unreliable, #561).
- **Shim** — a few-line local edit that only adjusts a version check or renames a symbol so an
  unmodified package loads; never a change of what the package computes (Q3).
- **Two-layer manifest** — sha256 of every userbase file before and after the lane: layer (i)
  must be identical; layer (ii) is the named cache directories that may change, reported file
  by file; any unlisted difference is a finding (§8.1).
- **"The D-A analogue"** — PyTransport/CppTransport also run symbolics once and emit numeric
  code; D-A (decided in R-1) is our version of that split: Wolfram once at derivation time, a
  numeric contract at sampling time.
- **Lane** — the single-license Wolfram kernel; one `wolframscript` machine-wide.

## 1. Context

`R-C.md` is READY (dispatched 2026-09-15, wave board `COSMOLOGY_PROGRAM.md:726`). D-A is
recorded: one committed Wolfram package run by a fixed `driver.wls` on WXF data, so Part 4 must
judge each tool on *whether a committed package can call its functions with data*
(`interfaces_decision.md:210-212`). D-B (Option A′, the shared theory file) means the FRW
branch needs the solver-only block defined above — Part 1 lists what it must carry.
`conventions.md` is canonical: `(+,−,−,−)`, `ε₀₁₂₃ = +1`, CAMB perturbation variables at the
seam, every sign difference vs Ma & Bertschinger attributed to convention or physics (§6).

Two physics constraints frame the survey (`R-C.md:43-52`): standard sectors never enter the
Lagrangian (the split above); and covariant derivatives do not collapse to partials in PGT, so
torsion representability is a required probe for every tool. Integration target (iii):
unmodified CAMB chained to our own solver, so Part 1 also says who consumes each requirement.

## 2. Evidence settled during planning — why each thing was looked at, and what it changed

Every table below opens with the question it answers and closes with the effect on the plan.
Sources are pinned so the memo can cite them without re-deriving; nothing here is left to
re-research.

### 2.1 xPand 0.4.4 — read from the tarball itself

**Why.** xPand is the package the brief names as the obvious candidate; the first plan judged it
from a 2013 paper and a stale 2015 mirror. Whether it even loads on our certified bundle, what
changed since the mirror, and what its own documentation claims decide whether Part 4 can
recommend building on it at all.

| fact | value | why it matters |
| --- | --- | --- |
| tarball `http://www2.iap.fr/users/pitrou/xPand_0.4.4.tar.gz`, 1,603,346 bytes, **sha256 `26e7abcac7bb655235ec39b73850729cf4465748d0d5ea2ad03c0607aef5ceab`**; flat layout `xPand/{xPand.m, Kernel/init.m, xPand.nb, xPand.History, gpl.txt, installation_xpand.txt, Documentation/English/xPandDoc.nb, Examples/*.nb (0–10), Article/xPand.pdf, BenchMarking/, xpand_acknowledged_papers.html, xPandcitations.py}` plus macOS `._*`/`.DS_Store` junk | tarball listing | the installer pins this hash and strips the junk; no nested-folder branch is needed |
| **no Wayback snapshot exists** (`archive.org/wayback/available` → empty) | archive.org API | the source is a personal page; a snapshot is requested at execution (§5.3), outward-facing, asked for in review round 1 |
| `$Version = {"0.4.4",{2025,04,01}}` (`xPand.m:23`); `$xTensorVersionExpected = {"1.1.3",{2018,2,28}}`, `$xPertVersionExpected = {"1.0.6",{2018,2,28}}` (`:27-28`) | `xPand.m` | the minimums are far below our xTensor 1.3.0 / xPert 1.0.6 |
| load gate `If[Not@OrderedQ@Map[Last,{expected, installed}], Throw@Message[General::versions,…]]` (`:94`, `:96`) compares **dates** and throws only when the installed package is *older* | `xPand.m` | **the version gate cannot fire on our bundle**; the only load risk is API drift, so the contingency (Q3) is about drift, not versions. xTensor's history carries one xPand-related fix (`xTensor.History:1690`) |
| 0.4.3 → 0.4.4 diff: 404 diff lines, 136 functional — new Bianchi-specific gauge variants, FRW gauge rules re-listed with unchanged content, the Bianchi structure-constant tensor re-declared with a "TODO to check" comment (2019), GPL v2 → v3, copyright 2012-2025; `xPand.History` has no entry after 0.4.2 | `diff` against the GitHub 0.4.3 `xPand.m` | the 2013 paper and the 0.4.3 source describe the FRW functionality we use; nothing FRW-relevant changed |
| license GPL-3 (`gpl.txt`, header) | tarball | enters Part 4's license cost line (§4.4) |
| `xPandDoc.pdf` (15.3 MB) is a tutorial dump: loading, space-time options, general splitting, Ricci scalar, fluid velocity, stress-energy, gauge change, Bianchi; "To Be Done: null geodesic equation, Boltzmann equation, null geodesic deviation"; public-function list at its end; **no torsion or rank limitation stated**, the only rank statement being `VisualizeTensor` (rank ≤ 2) | `pdftotext` (4068 lines) | Part 3's "documented scope" column now rests on the 0.4.4 documentation, not on the paper alone |
| paper: Pitrou, Roy, Umeh, arXiv:1302.6174, CQG 30 (2013) 165002; INSPIRE recid 1221019, 129 citations; author-maintained list 109; conclusion: "**the current exception being any gravitational theory with torsion**" | TeX, INSPIRE | the authors themselves say torsion is out of scope — Part 2/3 must say how the torsion papers got around that (§2.6), and finding F14 corrects the H7 near-miss table |
| torsion-sector users: Bahamonde et al. 2009.02168 (f(T,B); tetrad-perturbation matrix `main.tex:272-276`; "the xAct packages [… xPand …] were also used", `:302`); Nicosia–Levi Said 2012.11959 (xAct incl. xPand generated the invariant table, `:331`); Barker 2206.00658 (HiGGS "grounded in" the xAct suite incl. xPand, `apstemplate.tex:168`; Minkowski) | TeX | "where it has been used" for a torsionful sector: as bookkeeping beside hand-parametrized tetrads, never as the torsion engine |
| headline users: second-order induced GWs (Inomata–Terada 1912.00785; Yuan et al. 1912.00885; Domènech–Sasaki 2012.14016), Bloomfield 1304.6712, Dai–Pajer 1504.00351/1502.02011, Umeh–Clarkson 1207.2109, Skordis–Złośnik 1905.09465, Arroja–Bartolo 1512.09374, Fidler–Tram 1708.07769, Spiers–Pound 2305.19332 | INSPIRE most-cited citers | the package is mature and widely used for exactly the second-order FRW work M3 needs on the metric side |

**Effect on the plan.** xPand stays the primary candidate **for the metric side** — the 3+1
split, the SVT decomposition, the gauges and the second-order machinery on FRW — despite its
authors' torsion exception, because that exception is about the *connection*: xPand's
conformal and Gauss–Codazzi rules assume the metric's Levi-Civita connection. After the
post-Riemannian rewrite (Terms) the torsion is an ordinary tensor field on that connection,
which is how every torsion paper in §2.6 proceeded and how legacy already works. Torsion
therefore enters through a *combination of existing tools*: xPand for the metric engine, a
decomposition step from another source (xMAG's `ToDistortion`, or the legacy idiom already
ported for Stage 1), and hand-written SVT rules for the rank-3 field — integrating what exists
is the first approach. Building the 3+1 engine ourselves (O2) is the fallback; forking xPand to
carry torsion natively (O6) is the **last resort**, because a fork drifts behind upstream (the
`slegner/CAMB` precedent, #498) and burdens every user with a non-standard package. Probes D
and E decide. The load probe (5.4) tests API
drift, not versions; the "if 0.4.4's arity differs" branches were removed; Part 3's row is
complete.

### 2.2 xPand internals that fix the probes' hypotheses (0.4.4 source)

**Why.** Review round 1 asked that each probe test a stated expectation grounded in the source,
not start from nothing: does xPand's split depend on the metric sign, and what does it do with
a rank-3 antisymmetric tensor?

| what the source says (item and lines) | consequence for the probes |
| --- | --- |
| `SetSlicing[g_?MetricQ, u_, normu_:-1, h_, cd_, {cdpost_, cdpre_}, SpaceTimeType_]` (`xPand.m:1710`) throws unless `SignDetOfMetric@g === -1` (`:1715`, the determinant sign, true for any Lorentzian 4-metric); defines `h` by `DefMetric[1, h, cd, InducedFrom -> {g,u}]` (`:1754`); sets `u·u = normu` (`:1804-1806`); xTensor's induced-metric machinery computes the normal's norm as `Scalar@ContractMetric[g u u]` (`xTensor.m:8706`) and uses `v v / norm` in its projector rules (`:8738-8739`) | mostly-minus is selected by `normu = +1` with `DefMetric[-1, …]` unchanged; **the split itself is sign-generic** |
| `ExtractComponents[expr, h, proj, inds]` (`:3031-3044`) folds over the free indices; the "Time"/`u` projection multiplies by **`−1` for an up index and `+1` for a down index** (`:3042`), i.e. `n·n = −1` hard-coded; no rank restriction (rank > 2 without a projector list only applies vanishing-background rules, `:3049`) | Probe C hypothesis: with `normu = +1` the split runs but "Time" projections come out with the wrong sign — project by hand with `n[-a]`. Probe D: projecting a rank-3 tensor is supported by design |
| `ToxPand` (`:3010-3024`) calls `SplitMatter[uf, duf, -1, h, …]` (`:3014`) — the fluid norm is hard-coded `−1` | in `(+,−,−,−)` build the fluid rules by hand with `SplitMatter[…, +1, …]` and use `ToxPandFromRules` |
| `ToxPandFromRules[expr, RulesList_List, h, n]` (`:2975-2985`) and `SplitPerturbations[expr, ListPairs_List, h]` (`:2883-2960`): the second argument is a **list of rules**; Gauss–Codazzi + `FullToInducedDerivativeAndCDDown`, then `CheckSTFTensors`; nothing rank-specific | hand-written rules for a rank-3 field are the intended interface |
| `DefProjectedTensorProperties::symmetrictensors` (`:191`): "can only add properties to tensors that are, at least, fully symmetric" | SVT properties go on the irreducible pieces (scalars, vectors, rank-2), never on the rank-3 object |
| internal `DummyS/DummyV/DummyT` (`:2066-2073`), rank 0/1/2 helpers used to build derivative-commutation rules | no rank-3 helper exists, so commutation rules for a rank-3 field's derivatives may be missing — Probe D checks the split of `∇T` explicitly |
| `$ConformalTime = True` → primes are d/dη; `ah = a[h]`, `Hh = H[h]` (`:1741`, `:1832-1834`) | our coordinate is conformal time, as required |

**Effect on the plan.** Probe C (5.8) and Probe D (5.9) now test named expectations; the
expected verdicts are written in the plan and either confirmed or overturned by the kernel.

### 2.3 Targets from TeX lines (never from page numbers)

**Why.** A reproduction is only as good as its target. The first plan cited the paper's page
numbers and the paper's printed outputs; reading the TeX showed two of those outputs disagree
with the textbook, so targets now come from Ma & Bertschinger and from a derivation, with the
paper's outputs kept as "expected shapes".

- 1302.6174 `xPand_-_arXiv2.tex`: minimal example `:1350-1360`; its output Out[47] `:1368-1373`
  (`6ℋ² + 6ℋ′ + 6𝒦 + ε(−12ℋ²φ − 12ℋ′φ − 6ℋφ′ + 12𝒦ψ − 18ℋψ′ − 6ψ″ − 2D²φ + 4D²ψ)`);
  Appendix A `:1733-1823`: `SetSlicing[…,"FLFlat"]` `:1737`, `MyToxPand` `:1742`, fluid
  `IndexSet` `:1757-1758`, `MyGR[μ,ν] := EinsteinCD + g Λ/κ − κT` `:1766`, **mixed indices**
  `MyGR[μ,-ν]` `:1770`, order-0 `{Time,Time}` = `Λa²/κ − 3ℋ² + κa²ρ` `:1777-1778`, order-0
  `{Space,Space}` `:1780-1783`, order-1 `{Time,Time}` `:1787-1789`, `{Time,Space}` `:1791-1796`,
  `{Space,Time}` `:1798-1805`, `{Space,Space}` `:1807-1823`.
- These targets are the **known answers the symbolic route must reproduce**: showing that the
  route (action → split → equations, and field equations → split) returns them is what
  qualifies it to derive the unknown equations of a user theory.
- **Two discrepancies between the paper's printed outputs and the textbook**: `:1789` prints
  `−2D_αD^α φ` where the Newtonian-gauge 00 equation carries the Laplacian of the *curvature*
  potential (xPand's ψ; Mukhanov–Feldman–Brandenberger eq. 4.15); `:1809` prints the TT
  friction as `E′ℋ` (coefficient 1) where `a²δG^i_j` gives `E″ + 2ℋE′ − D²E`. So the paper's
  outputs are *expected shapes*, not targets; the kernel adjudicates, and what 0.4.4 prints is
  recorded (finding F15).
- Ma & Bertschinger `literature/astro-ph_9506072/9506072.tex`: conventions `:273-289`,
  synchronous metric `:294-296`, conformal-Newtonian metric `:359-363`, Friedmann `:573-580`,
  synchronous Einstein `:603-618`, conformal-Newtonian Einstein `:644-661` (`ein-cona..d`),
  θ, σ `:623-628`, `T^μ_ν` `:679-685`, δ `:709-710`, gauge transformation `:552-566`; **no
  tensor equation** (scalars only, `:364-370`).
- `camb.symbolic` (camb 2.0.4 in `.venv`): constraints `symbolic.py:360-364`, evolution
  `:380-388`, Newtonian substitutions `:436-452`, `newtonian_gauge` `:455`, `cdm_gauge` `:486`,
  `synchronous_gauge` `:511`, `camb_fortran` `:727`, `compile_source_function_code` `:892`; the
  declared line element `ds² = a²((1+2Ψ_N)dt² − (1−2Φ_N)δ_ij dx^i dx^j)` `:461`.

**Effect on the plan.** `targets.wl` is typed from these lines before any kernel runs (5.2);
the MB reproduction is checked twice, once against MB and once against CAMB's own symbolic
equations in our convention (5.7).

### 2.4 The xAct family, and xMAG's dependency chain

**Why.** The brief requires a row for every xAct-family package, and the user widened the lane
(Q4) to probe xMAG hands-on because it is the one package whose documented scope is a
torsionful connection. Installing it additively requires knowing its exact dependencies and
that none of the targets exists already.

xAct 1.3.0 (29 Dec 2025, GPL) — core: xCore, xPerm, xTensor, xCoba. Contributed (every one
gets a Part-3 row; one line where irrelevant): xPert, Harmonics, Invar, Spinors, xPrint,
SymManipulator, AVF, xTras, TexAct, **xPand**, xTerior, SpinFrames, **xIST/COPPER**, EFTofPNG,
bimEX, **FieldsX**, **xPPN**, SymSpin, **HiGGS**, xBrauer, **xMAG**, **PSALTer**, TInvar,
SpaceSpinors, xCPS, xIdeal, **Hamilcar**.

- **xPert** 1.0.6 (2018-02-28) installed; its header says "only one perturbation structure can
  be defined; only single-parameter metric perturbation theory" — but legacy already declares
  a `DefTensorPerturbation` beside `DefMetricPerturbation` (`_derive.py:2350-2376`), so a
  second perturbed field is not forbidden. Paper 0807.0824, recid 790000, 264 cites.
- **xTras** 1.4.2 (2014-10-30) installed; 1308.3493, 289 cites.
- **xMAG** (Helpin, `github.com/THelpin/xMAG`, main `88026e47…` 2026-01-10). In plain terms:
  by its documented scope it is the one xAct package built for exactly the torsion side of our
  problem — an independent connection with torsion (and non-metricity), the decomposition into
  Levi-Civita plus contortion, perturbations of the connection, variations with respect to it,
  and a hypersurface decomposition — so it is the strongest candidate for the step xPand does
  not do — the expansion step done first, with xPand and our rules producing the perturbation
  equations in the resulting fields afterwards. Two things it is not, neither a disqualifier:
  it has no FRW functions (no single tool does what we are building; the plan chains
  packages), and it is young (version 0.1.0, eleven commits, no paper, a "TODO" in its source)
  — its author maintains several xAct packages, so defects found are drafted as upstream
  issues and, with the user's agreement, filed or sent to the author. It is probed rather than
  trusted (Q4) because assertions get verified before they land. Facts: `xMAG.m` 0.1.0
  (2023-05-18), **GPL-2+** (header), expects xTensor 1.2.0 and xTras 1.0.6 as minimums
  (installed 1.3.0 / 1.4.2 satisfy), Mathematica 8+; needs xTensor, xPert, xPerm, xCore,
  ExpressionManipulation, **`xAct`TraceFree``**, xTras, **`xAct`xBrauer``**. What it offers:
  `DefCovD` with `Master -> met`, `FromMetric -> Null`, `Torsion -> True` defines an
  independent connection with torsion and non-metricity;
  `Distortion/Contorsion/NonMetricity/TorsionVector[covd]`; **`ToDistortion[exp, CD, g]`**
  rewrites curvature and derivatives of the
  independent connection as Levi-Civita quantities plus distortion terms — exactly the step the
  torsion papers do by hand (§2.6), and the same functions cover non-metricity, so the route
  extends to metric-affine perturbations later; `DefConnectionPerturbation[christoffel, pert, param]`
  perturbs the connection with xPert (it does not call `DefMetricPerturbation` — a "To Do"
  comment); `VarD`/`VarL` with respect to the connection; `StartInducedDecomposition`
  (hypersurface decomposition, possibly an independent 3+1 route). No cosmology functions; no
  paper (Helpin's thesis 2407.18019 describes the Brauer packages). **Dependency chain, all five
  targets verified absent on 2026-09-15:** `Applications/SymmetricFunctions/` and
  `Applications/BrauerAlgebra/` (plain packages, from `THelpin/xBrauer_Bundle` main `48be67e1…`
  2026-01-10), `Applications/xAct/xBrauer/` (1.1.0, 2023-11-29, GPL-2+, expects xTensor 1.2.0;
  same bundle), `Applications/xAct/TraceFree/` (`xAct-contrib/TraceFree` master `4e53ab39…`
  2013-11-12; Leo Stein, GPL-2+, 0.1.0; a `TraceFree` option for `DefTensor` plus
  `TraceFreeRules`/`SetTraceFreeRules`), `Applications/xAct/xMAG/`. No LICENSE files in any of
  the repos; the notices are in the `.m` headers.
- **xIST/COPPER** (Noller; `xAct-contrib/xIST`): builds the most general quadratic action for
  scalar-/vector-tensor perturbations on FRW from symmetry (Noether constraints) rather than
  from a covariant Lagrangian; xIST 0.7.3 (2016-03), COPPER 0.8.3 (2016-07); GPL-3; tested on
  Mathematica 9/10; 1604.01396 (79 cites). **Not a candidate component**: it serves an end user
  probing theory space, not the derivation we are building; Part-3 row only.
- **HiGGS** (Barker): Hamiltonian analysis of PGT in tetrad + spin-connection variables on
  Minkowski; deprecated Jan 2026, superseded by **Hamilcar** (GPL-3, 0.0.0-developer,
  2512.25007, local); 2206.00658, 2205.13534, 2101.02645.
- **xPPN** (Hohmann, 2012.14984, EPJC 81 (2021) 504, 11 cites): the homepage returns 404, so the
  code was located on GitHub, `xenos1984/xPPN` (master `dcaabbea…` 2022-10-01; `xPPN.m`
  62.9 KB, examples for GR, Brans–Dicke, New GR, scalar-torsion; no LICENSE file and no
  license/version line in the header; needs xTensor, xPerm, xCore only). It provides
  `SpaceTimeSplit(s)` (3+1 by expanding dummies), background/physical tetrads, torsion rules
  and the Weitzenböck derivative, and non-metricity tensors — **the one xAct package with tetrad and
  torsion 3+1 machinery**, built for the post-Newtonian expansion (Minkowski + velocity
  orders), not for FRW. Could it serve FRW anyway? No: its background objects are hard-wired
  flat (`BkgMetricM4` is η_μν, `BkgMetricT1` is η₀₀ = −1, `BkgTetradT1` is Δ⁰₀ = 1 in its own
  usage messages) and its split is an expansion in velocity orders about Minkowski; its
  variable set (tetrad + flat Weitzenböck connection) is also not ours (metric + torsion
  tensor). Against xMAG, which defines an independent connection on any metric and assumes
  nothing about the background, it is the more limited choice — documentary row.
- **FieldsX** (Fröb, 2008.12422, 29 cites): fermions, gauge fields, BRST; v3 adds frame fields
  and spin connections; GPL-2; no cosmology. Like xIST, a theory-building aid for users rather
  than a derivation component; Part-3 row only. **PSALTer** 2.0.2 `bb45adb0` installed (the
  Minkowski spectrum branch; 2406.09500, recid 2798424).

**Effect on the plan.** The installer for xMAG (5.10) knows its five directories, sources,
commits and licenses; Part 3 has a row for each package; Part 4 gains the hybrid options O3/O4.

### 2.5 Non-xAct comparables

**Why.** The brief asks what other ecosystems do (equations-in solvers, EFT code generators,
inflation code generators, generic CAS) so Part 2 can say what each would leave us to do, and
Part 3 can record licenses and releases.

| tool | facts | what it leaves us |
| --- | --- | --- |
| SymBoltz.jl | MIT; arXiv:2509.24740, A&A 707 (2026) A128 (local TeX); repo `github.com/hersle/SymBoltz.jl` (R-1's pin `3d1f20a3` = v1.7.0 was recorded without a URL — the memo states it); ModelingToolkit; **Newtonian gauge, scalars only**, `(−,+,+,+)`; a new species is a `System(eqs, τ, vars, pars; initialization_eqs)` whose `ρ, P, δ, θ, σ` couple automatically through the Einstein equations (`extended_models/`) | takes equations, never actions; **scalars-only and Newtonian-only, which excludes our channels** (tensor `h_ij`, rank-3 torsion, photon polarization) — at most a future consumer for scalar pieces, never a producer |
| `camb.symbolic` / CAMB 2.0.4 | sympy scalar ΛCDM equations in covariant variables with `cdm_gauge()`, `newtonian_gauge()`, `synchronous_gauge()`, `make_frame_invariant()`; `get_hierarchies()`, `get_scalar_temperature_sources()`; `camb_fortran()` + `compile_source_function_code()`; **the only custom-source hook is `CAMBparams.set_custom_scalar_sources` (`model.py:1111`, scalar sources); no tensor or eikonal custom-source API** (grep of the installed package); `get_time_evolution(q, eta, vars, frame="CDM")` (`results.py:557`); its license file (British-spelled filename upstream) = LGPL-3 code + GFDL-1.3 docs + an arXiv-submission clause (GitHub: NOASSERTION) | the reference for gauge names and variable definitions; a code generator for scalar sources only; everything tensor/eikonal is ours |
| hi_class | Horndeski via α-functions hand-coded in C; 1605.06102 (245 cites), 1909.01828; no license file and no license sentence in the README | no symbolic derivation; no torsion → Horndeski map except special cases (Barker et al. 2006.03581) |
| CLASS | 1104.2932/2933; README: "You can use CLASS freely, provided that in your publications, you cite at least the paper CLASS II"; no license file detected | hand-derived hierarchies |
| EFTCAMB / H-EFTCAMB | EFT of DE in CAMB, Fortran model classes; 1312.5742 (282 cites), 1405.3590; H-EFTCAMB 2603.01662; CAMB's license file + `eftcamb/LICENSE` | same as hi_class |
| EFT of DE (method) | 1210.0201, 1304.4840, 1404.3713 | the α-function mapping is hand-derived per theory |
| CppTransport / PyTransport | model file → GiNaC / SymPy → generated C++; flat FRW inflation; GPL-2+ / GPL-3+; 2018.1 / 2.0 (2017); 1609.00380, 1609.00381 | the D-A analogue; scalar fields only |
| CAMB v2 / astro-ph/9911177 | 2607.14854 (local), 9911177 (5,326 cites) | the consumer of our tables |
| Cadabra2 | 2.5.14 (31 Jul 2025); GPL-3; hep-th/0701238, 1912.08839 | a generic CAS with no cosmological-perturbation machinery |

**Effect on the plan.** Part 2's method list (§4.2) and Part 3's license/release cells are
filled; Part 1's consumer column (§4.1) rests on the CAMB API facts.

### 2.6 Torsion / non-Riemannian perturbations on FRW — how the papers did it (TeX read)

**Why.** Review round 1 asked that Probe D be designed from what prior work actually did, so
that "representable with work" names their work, not ours. The TeX (acknowledgments, appendices)
of every relevant paper was read for the variables, the method and the computer algebra used.

| paper | variables and method | CAS | TeX lines |
| --- | --- | --- | --- |
| Aoki–Bahamonde–Gigante Valcarcel–Gorji **2310.16007** (PRD 110, 024017; 30 cites) | metric-affine on curved FLRW; `Γ̃ = Γ + K + L` (`:166`); torsion and non-metricity are tensor fields with their own perturbations; **FLRW background torsion has a vector and an axial function of time**; **torsion perturbation = 4 scalars {T,B,φ,A} + 4 pseudoscalars {𝒮,ℬ,ϱ,𝒜} + 3 vectors + 3 pseudovectors + 1 TT tensor + 1 TT pseudotensor = 8 + 12 + 4 = 24 d.o.f.** (Table `Table:Torsion`, `:864-878`); explicit parametrization per helicity (`:768-776` h=2, `:798-807` h=1, `:834-842` h=0); second-order action; `(−,+,+,+)` (`:153`) | **none acknowledged in the TeX** (`:1169-1170`) | `CosmoPertMAG.tex` |
| Heisenberg–Hohmann–Kuhn **2311.05495** (69 cites); Heisenberg–Hohmann **2311.05597** | teleparallel; connection perturbation as a tensor `δΓ = ∇δΛ`; SVT counts per family; second-order action; gauge-invariant variables | **none acknowledged** (`7-Discussion.tex:20-25`; `gentelepert.tex:1087-1089`) | TeX |
| Bahamonde et al. **2009.02168** (f(T,B); 38 cites) | tetrad-perturbation matrix in the Weitzenböck gauge (`:272-276`), induced metric (`:278-283`), fluid `δΘ_μν` (`:292-301`); SVT then Fourier | "the xAct packages [… **xPand** …] were also used" (`:302`) | `main.tex` |
| Nicosia–Levi Said **2012.11959** | invariant enumeration; torsion via tetrad | xAct incl. xPand generated the invariant table (`:331`) | `main.tex` |
| Nikiforova–Damour **1804.09215** (18 cites) | tetrad `e^i_μ = e^φ(δ + ε)` and Lorentz connection `A_ijμ = Ā + a_ijμ` on torsionful de Sitter; **background connection carries two nonzero functions** `Ā_{0ab} = e^φ f δ_ab`, `Ā_{abc} = e^φ g ε_abc` (frame indices) (`:724-730`); gauges: symmetric ε, zero shift + longitudinal (`:765-784`); **connection perturbation decomposed by hand into 8 scalars, 6 transverse vectors, 2 TT tensors** with explicit `k_a`, `δ_ab`, `ε_abc` structures (`:819-852`); helicity projection (`:853-873`) | **xAct "instrumental"** (`:2296-2301`) | `NDfinal.tex` |
| Lu–Chee **1601.03943** (JHEP 05 (2016) 024) | PGT with pseudoscalar torsion; scalar perturbations only; vierbein `e^0_μ = δ(1+φ)`, `e^a_μ = aδ(1−ψ)` (`:721-727`, longitudinal gauge, conformal time `:733-735`); torsion carried by two scalar modes `h, f` with **nonzero background** and perturbations `δT_ij0 = δ_ij a² δh`, `δT_ijk = 2ε_ijk a³ δf` (`:744-747`) | **Maple** (`:716-718`) | `2R-PGT-DE.tex` |
| Lu–Chee 1401.2585 | withdrawn ("expressions are unclear"); no source | — | arXiv |
| Barker **2206.00658** | HiGGS on Minkowski; cites xPand as part of the suite (`:168`) | xAct | `apstemplate.tex` |

Also in the map (abstracts): Hohmann 2011.02491; Golovnev–Koivisto 1808.05565 (152 cites);
Hohmann–Pfeifer 2203.01856; Toporensky–Tretyakov 2110.12332; Chen–Ho–Nester–Wang–Yo
0908.3323; Nikiforova–Randjbar-Daemi–Rubakov 0905.3732 (local); Garcia de Andrade
astro-ph/0005519, gr-qc/0202022; Barker et al. 2003.02690 (local), 2006.03581 (local);
Bahamonde et al. 2506.17017 (local); Obukhov 1702.05185 (local). INSPIRE title searches for
"Poincaré gauge" + cosmological + perturbations return 0 hits — the field is thin, as expected
for genuinely new work: the literature supplies parts of the method (the rewrite, the SVT
recipe, which CAS carried it), never a whole method to copy; the memo says which part comes
from where.

**What the precedent says for Probe D.** Every paper that carried torsion through a
Levi-Civita-based tool did the same three things: (a) wrote the torsion (or connection)
perturbation as an explicit **hand-parametrized SVT set — 8 scalars, 6 transverse vectors,
2 TT tensors**; (b) kept the **background torsion nonzero** where symmetry allows it
(1804.09215's `f, g`; 1601.03943's `h, f`; 2310.16007's two functions) — for a spectator theory
this is a homogeneous mode of the new sector, subject to the #501 gate; (c) derived from the
**second-order action** (2310.16007, 2311.05495) or the field equations (1804.09215,
1601.03943). None used xPand's projected-tensor machinery for the torsion itself; the two xPand
users used tetrad matrices and xPand for the metric/bookkeeping side. Two different things are
taken from these papers: their *background* torsion solutions (`f, g`; `h, f`; `T₁, T₂`) belong
to the background problem and are outside our scope beyond negligibility; their *perturbation*
decompositions (2310.16007's helicity table, 1804.09215's 8 + 6 + 2) are the recipe for turning
the rank-3 torsion that a decomposition step (xMAG or the legacy idiom) produces into scalar,
vector and tensor fields a solver can evolve — that is what D3 transcribes. **The same step exists in
this repository already:** legacy TIDAL's `_wls_torsion_curvature_decomposition`
(`_derive.py:2164-2189`) does the Levi-Civita + contortion rewrite by hand on flat space
(`ChangeCurvature[L, CDT, CD]` then `K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c − T_{bc}{}^a)`),
the spectrum branch ports it for Stage 1 (`spectrum_design.md:341-351`), and xMAG's
`ToDistortion` + `BreakDistortion` package it. **Effect:** "representable with work" in Probe D
names exactly (a)–(c), transcribed from 2310.16007's lines; D5 exists because of (b); the
decomposition tool is a choice that binds both branches (§4.4).

### 2.7 Legacy capability (Part 3 item 5)

**Why.** The brief asks what legacy already does for the same field content, to compare against
and not to build on. Knowing that legacy already declares a torsionful xAct connection (and how
it rewrites it) supplies the idioms Probe D reuses.

- `tidal/wolfram/ComponentDecompose.wl`: plain coordinate components (`:1640-1699`, rank ≥ 3 at
  `:1783-1836`); **no SVT, no Fourier, no conformal time, no symbolic `a(η)`/ℋ**; the metric is
  a literal matrix (`_derive.py:79-120`), curvature by `D[]` on it (`CommonUtilities.wl:440-458`);
  signature hard-wired `(−,+,+,+)` (`_derive.py:62-65`).
- A torsionful connection **has** been declared: `_derive.py:506-524` emits
  `DefCovD[CDT[-a], {"#","DT"}, FromMetric -> g, Torsion -> True]`; `:2350-2376` declares
  `T[a,-b,-c]` antisymmetric + `DefTensorPerturbation` with `T̄ = 0` on flat Minkowski; the
  Levi-Civita + contortion rewrite is `ChangeCurvature[L, CDT, CD]` followed by the
  contortion → torsion rule on `ChristoffelCDCDT` (`_derive.py:2164-2189`) and `K̄ = 0`
  (`:2553-2557`) — reused verbatim in Probe D8. Test `tests/test_cli.py:5009-5053` (dry-run).
- `docs/tex/multi_field_perturbation.tex:100-225, 289-336` documents the xPert pipeline — the
  right legacy reference (finding F11), not `background_fields.tex`.
- Frozen oracle `de_sitter_kg` (`tests_cosmo/data/oracles/specs/de_sitter_kg.json`, README
  `:93-107`): 2+1D, metric `e^{2Ht}diag(−1,1,1)`, so `t` **is** conformal time with
  `a = e^{Hη}`; its equation `φ″ + ℋφ′ − ∇²φ + a²m²φ = 0` is the `(D−2)ℋ` friction family with
  D = 3, against which 5.6's scalar control (D = 4, `2ℋ`) is compared documentarily.
- `docs/tex/background_validity.tex:287-310`: `T̄ = 0` is exact only on flat Minkowski +
  uniform `B₀`; on FRW the symmetry-allowed background torsion is the vector + axial pair.

**Effect on the plan.** Legacy comparison is documentary (5.11); D8 reuses legacy's rewrite.

### 2.8 Local literature and the fetch list (Q2)

**Why.** The repo rule is "read TeX locally"; the memo must cite paper *and* code per method.
Present locally: Ma & Bertschinger, SymBoltz, xTras, xPerm, CLASS II, CAMB v2, TorC,
2003.02690, 2006.03581, 2506.17017, 0905.3732, 1702.05185, 2507.02362, Hamilcar 2512.25007,
the PSALTer set. Fetched into the scratchpad during planning and re-fetched into the
worktree's gitignored `literature/` at execution: 1302.6174, 2310.16007, 2311.05495,
2311.05597, 1601.03943, 1804.09215, 2009.02168, 2012.11959, 2206.00658, 2012.14984; plus,
only where the memo cites them for a claim: 0807.0824, 2011.02491, 1808.05565, 1604.01396,
1210.0201, 1304.4840, 1404.3713, 1605.06102, 1909.01828, 1312.5742, 1405.3590, 1609.00380,
1609.00381, astro-ph/9911177, hep-th/0701238, 2008.12422, 2407.18019, 2212.14496, 2203.01856,
2110.12332. `literature/README.md` is **not** regenerated (Q2); `docs/references.md` gets the
curated rows.

### 2.9 Mechanics (from the R-1 precedent; exact commands)

**Why.** R-1 is the precedent for how a research lane commits, proves and reports; every rule
below was read from its scripts and the delegation protocol so that this lane's artifacts are
in the same shape and pass the same gates.

- Worktree: `git worktree add /tmp/tidal-rc -b cosmo/rc-perturbation-tooling feat/cosmology-program`;
  `uv sync --all-extras` inside it; never merge; remove worktree + remote branch after merge
  (`wave_boundary_check.sh:29-30`).
- Gates (each starts a kernel; confirm the lane is free with `pgrep -x WolframKernel`):
  `bash scripts/psalter/ensure_registered.sh` (~1 min, includes verify) and
  `bash scripts/verify-wolfram-setup.sh --require-psalter` (~2 min; **exit 0 required; exit 2 =
  DEGRADED = refusal**). Expected: "All checks passed!" + the one informational xPerm-ldd
  warning R-1 saw.
- Fingerprint: `verify-wolfram-setup.sh:60,216-232` greps `$Version` of exactly
  xCore/xPerm/xTensor/xCoba — an added directory cannot change it (the memo names the four,
  finding F10). xPert 1.0.6 is in no fingerprint — recorded from `xPert.m:1` before and after.
- Install precedent: `scripts/install-psalter.sh` (kernel file test, not `command -v`; copy
  into `Applications/xAct/`; `INSTALLED_COMMIT` stamp). **Never run
  `scripts/install-xact-xcoba.sh`** (`rm -rf xAct`).
- Every kernel launch: `QT_QPA_PLATFORM=offscreen`, throwaway cwd, judged by sentinel lines and
  artifacts, never exit status; scripts self-guard with `pgrep -x WolframKernel` (the hook
  cannot see a bare `.wls`, and resolves `*.sh` tokens against the main clone's cwd, not the
  worktree's).
- Draft PR at the first commit:
  `gh pr create --draft --base feat/cosmology-program --head cosmo/rc-perturbation-tooling --title "docs(cosmology): R-C perturbation tooling survey"`;
  verdict via `gh run view <id> --json conclusion`, written `CI <run-id>: <conclusion>`.
- CI on the new files: `ruff check` + `ruff format --check .` (the `.py` file), pytest incl.
  `tests/test_repo_hygiene.py` (tracked files under `docs/` and `scripts/`; forbids the container
  clone path, a user home path and the Claude project slug — **scrub every pasted transcript**; `git add` then run it),
  cspell strict over `docs/cosmology/**` and `scripts/**` incl. `.wls`/`.wl` (`xPand` is in no
  dictionary → in-file `cspell:words` directives; British spellings hard-flagged). Pre-commit
  adds mypy on the `.py`, `end-of-file-fixer`, `check-added-large-files --maxkb=1000`.
- `docs/README.md` and `handoffs/README.md` rows are the orchestrator's (`r1_planning_record.md:635`):
  ready-to-paste rows go in the report.
- Evidence: run ledger + evidence index inside the memo (R-1 §1.1, §7); payloads in gitignored
  `third_party/perturbations_runs/`; `docs/cosmology/evidence/` is the Tier-1 convention, unused.
- Issues: file one per discovery that outlives the memo, after `gh issue list -S <keyword>`
  (R-1 filed #569–#573); design-document contradictions are reported to the orchestrator in the
  PR and filed as issues labeled for the orchestrator, never edited here.

## 3. Deliverables — the memo and the record of how it was reached

| path | what | why it exists |
| --- | --- | --- |
| `docs/cosmology/perturbation_tooling.md` | the memo (≤ ~600 lines): status blockquote, `cspell:words`, §0 summary, Parts 1–4, **"what was verified here and what was not"** (run ledger: step / what / kernel wall / result; not-tested table with reasons), **evidence index** (digests, sha256s, snapshot URL, printouts), findings routed | the decision input for D-C |
| `docs/cosmology/rc_planning_record.md` | this plan archived verbatim (Q1) with Appendix A (revision log), B (evidence ledger), C (not verified) | the reasoning trail — why each investigation was done and what it changed |
| `docs/references.md` | curated rows under `## Cosmology Program` for every paper the memo relies on, incl. xPand (Q2) | the curated index says *why* each paper is in the library |
| `scripts/research/perturbations/README.md` | status, contents table with a provenance column, "How to reproduce", "Rules this directory follows" | a future session reruns the evidence instead of trusting it |
| `…/install_xpand.sh`, `…/install_xmag.sh`, `…/manifest.sh`, `…/run_lane.sh` | additive installers (§5.3, §5.10), the two-layer manifest (§8.1), the serial lane runner (§5.0) | reproducible install and proof of additivity |
| `…/wolfram/RCSetup.wl`, `probe_load.wls`, `repro_a2_tensor_eom.wls`, `repro_a1_tensor_action.wls`, `repro_b_mb_scalars.wls`, `probe_c_signature.wls`, `probe_d_torsion.wls`, `probe_e_xmag.wls`, `compare.wls`, `targets.wl` | shared setup, the reproductions and probes, the harness, the hand-typed targets with TeX citations | "tool X can do Y" as scripts |
| `…/mb_camb_symbolic.py` + `…/camb_symbolic_newtonian_2.0.4.txt` | prints `camb.symbolic`'s Newtonian- and synchronous-gauge scalar equations; the committed printout the targets cite (regenerable) | the machine check of the MB → our-convention transcription |
| `…/xpand_upstream_issue.md`, `…/patches/xpand-shim.diff` (only if a shim was needed) | drafted upstream issue (never filed; precedent `docs/cosmology/psalter_543_upstream_issue.md`) and the unified diff, also stamped in `INSTALLED_VERSION` | Q3 evidence |
| GitHub issues | one per discovery that outlives the memo (e.g. F14 the near-miss table, F15 the paper's printed outputs, any xPand/xMAG defect found), plus the design-doc contradictions F1–F13 labeled for the orchestrator | the searchable trail |
| userbase | six new directories, each stamped: `Applications/xAct/xPand/` (`INSTALLED_VERSION`), `Applications/SymmetricFunctions/`, `Applications/BrauerAlgebra/`, `Applications/xAct/xBrauer/`, `Applications/xAct/TraceFree/`, `Applications/xAct/xMAG/` (`INSTALLED_COMMIT`) — nothing else written | the hands-on evidence |
| report-back | branch · PR · `CI <id>: <conclusion>` · per-criterion outputs · Part-4 table and decision table · amendments · routed findings and issue numbers · not-tested list · every userbase directory added · which packages D-C would make install requirements (installer promotion, fingerprint, guide) · fetched arXiv ids · ready-to-paste index rows | what the orchestrator needs to record D-C |

## 4. Memo design

### 4.1 Part 1 — requirement checklist (one line each: requirement · source · consumer)

Consumer values, each with its source: **CAMB API** (only `set_custom_scalar_sources`, scalar
sources, and `get_time_evolution(frame=)` for background/standard-mode tables —
`model.py:1111`, `results.py:557`); **WS3 solver** (η-grid assembly of `M(η)` and the source
vector, `solver_design.md:532-533`, `:344-346`); **WS4 line of sight** (tensor `Δ_ℓ^B(k)` and
`C_ℓ^BB`, `observable_ladder.md:141`); **WS2 gate** (validity code, `:139`, `:143`). Where the
consumer changes the requirement, the line says so:

- **Frame and coordinates** (consumer: WS3) — coupled block of the new sector on FRW
  (`R-C.md:83-84`; `observable_ladder.md:136`); conformal time = `camb.symbolic` `t`
  (`repo_reshape.md:443-446, 827`; `conventions.md:75`); `a(η)`, `ℋ(η)` symbolic to the export
  (`repo_reshape.md:828-829`); coefficients evaluated from CAMB's table at solve time
  (`observable_ladder.md:137`); conformal-weight fast path derived from the spec
  (`birefringence_notes.md:68-87`); limits `a → const` and de Sitter (`COSMOLOGY_PROGRAM.md:1060-1062`, #560 caveat).
- **Gauge** (consumer: CAMB seam + WS3) — explicit named input from CAMB's set; projection in
  Wolfram with `camb.symbolic` as reference; spec metadata; seam asserts via
  `get_time_evolution(frame=)` (`repo_reshape.md:455-469, 830-831`). The reproductions exercise
  both CAMB-named gauges (Newtonian and synchronous; CAMB's CDM frame is synchronous with
  covariant variable names, `conventions.md:65`) so the memo shows the interface for each name
  CAMB uses.
- **Field content / connection** (consumer: WS2) — rank-3 torsion and a torsionful connection
  representable (`R-C.md:49-52`); post-Riemannian decomposition about `(η, T̄)`, with a
  **homogeneous new-sector mode** where FRW symmetry allows it and the theory admits it
  (`spectrum_design.md:344-351`; §2.6; O4-iso row `COSMOLOGY_PROGRAM.md:211`); the map from
  PSALTer's field content (tetrad + spin connection or metric + torsion) to the FRW variables
  belongs to the solver-only block (`interfaces_decision.md:413`); the **post-Riemannian
  rewrite is performed in both branches** (`spectrum_design.md:341-351, 372-378`), so the tool
  chosen for it binds Stage 1 as well as M3; derived fields `F = dA` deferrable
  (`multi_field_perturbation.tex:289-336`).
- **Coupled block vs standard sectors** (consumer: WS3; CAMB supplies tables and `S_std`) —
  the split stated (`R-C.md:43-48`; `repo_reshape.md:233-237`); tensor target
  `h″ + 2ℋh′ + k²h = S_std + S_new` with `S_std` = CAMB's anisotropic-stress source and `S_new`
  = the mixing with the new sector's own perturbations, **consumed by our solver, not CAMB**
  (no tensor custom-source API) (`solver_design.md:330-343`); reductions: Cembranos eq. 20
  (`:331-334`), decoupled-limit CAMB BB at `r = 0.036` (`observable_ladder.md:157`).
- **Exports** — per-channel `S(k,η)` derived symbolically per theory (`repo_reshape.md:329-348`),
  **consumed by WS4's own line-of-sight code for tensors** (CAMB's hook is scalar-only) and
  optionally by `set_custom_scalar_sources` for scalar channels, which the spectator scope
  excludes (`spectator_route.md:122-131`); the tensor source is the standard one for a pure
  `h_ij` (Seljak–Zaldarriaga astro-ph/9603033, local; CAMB's Fortran tensor source,
  `solver_design.md:336-343`) and must be **re-derived symbolically when the photon sector is
  modified** (a rotation belongs inside the line-of-sight integral,
  `birefringence_notes.md:181`) — #514, owed by WS2 to WS4; R-C records the requirement and the method sources and
  does not derive it; **the solver forms (Terms) without choosing a
  solver**: the second-order block, its first-order state form `y′ = M(η) y (+ s)` that every
  bake-off candidate consumes (Magnus, matrix-WKB, `rk-adaptive`; `solver_design.md` §7), and
  for the photon channel the eikonal amplitude form `ψ′ = −iM(η)ψ` with carrier, state vector,
  `M` in `(η, a, ℋ, n_e, B, ω)`, Hermitian split, polarization basis and **dropped terms
  exported** (`:472-486`; consumer: WS3 O3 front-end); several representations per derivation
  (`:466-470`); batched-assembly-ready coefficients (`:532-533`). The memo lists what form each
  candidate demands and leaves the choice to the bake-off.
- **Background-consistency gate (#501)** (consumer: WS2 gate) — the order-1 coefficient
  exported as a residual evaluated on the CAMB background; a per-theory **admissibility gate**
  with a derivable tolerance, never a propagated source; specified so it also covers a
  homogeneous new-sector mode (`spectator_route.md:81-91, 107-120`;
  `observable_ladder.md:185-193`; `COSMOLOGY_PROGRAM.md:200, 211`).
- **Conventions** (consumer: every artifact reader) — native CAMB + PSALTer emission;
  `(+,−,−,−)`, `ε₀₁₂₃ = +1` recorded and asserted; CAMB variables verbatim at the seam; every
  MB sign difference attributed (`conventions.md` §1–§3, §6).
- **Solver-only block must carry (listed, not designed)** — background fields, field →
  background + perturbation map, gauge choice, own fingerprint
  (`interfaces_decision.md:408-416, 421-422`).
- **Auxiliary** — energy/Hamiltonian time-dependence (contested, F6); validity flags derivable
  (`observable_ladder.md:143`; `ρ_new/ρ_γ` vs `ΔN_eff`, the spectator condition itself).

Each line ends with the tools that cover it (from Part 3) so Part 4 reads column-wise.

### 4.2 Part 2 — methods (each: input · output · connection assumption · what it leaves us)

1. xPert/xPand line and successors (§2.1, §2.6 users).
2. EFT-of-DE → hi_class / EFTCAMB (α-functions hand-mapped; metric only; no torsion →
   Horndeski map except Barker et al. 2006.03581's special cases).
3. Symmetry-first quadratic action (xIST/COPPER).
4. SymBoltz.jl equations-in (Newtonian, scalars; a consumer of derived equations, never a
   producer).
5. CLASS/CAMB hand-derived hierarchies + `camb.symbolic` (gauge/variable reference;
   scalar source-function code generation).
6. PyTransport/CppTransport symbolic core → generated numerics (the D-A analogue).
7. Torsion / non-Riemannian on FRW (§2.6): Levi-Civita + distortion, symmetry-allowed background
   torsion, hand-parametrized SVT of rank-3, second-order action or field equations,
   gauge-invariant variables; which CAS each used.
8. Cadabra as a generic CAS row.
9. Line-of-sight source construction — Seljak–Zaldarriaga (astro-ph/9603033, local), the
   Hu–White total-angular-momentum method (astro-ph/9702170), CAMB's tensor source in Fortran
   and `camb.symbolic.get_scalar_temperature_sources` as the scalar code-generation precedent —
   the method WS4/M3 will follow for #514; recorded here so the later derivation starts from
   the literature, not from scratch.

Every method cites paper **and** code (repo URL + version/commit).

### 4.3 Part 3 — tools

Columns: tool · documented scope (manual/paper) · where used · license · last release ·
certified-bundle compatibility (14.3.0 × xAct 1.3.0 × xPert 1.0.6) · Part-1 lines covered ·
**torsion verdict**, each labeled *probed here* or *documented only*. Rows: all of §2.4 and §2.5
plus `ComponentDecompose.wl`. One line per irrelevant xAct package.

### 4.4 Part 4 — recommendation per requirement, decision table, D-C options, license cost

For each Part-1 line: build-on / borrow-pieces / build-own · cost (hours or days, with what the
estimate rests on) · interface under D-A (which function the committed package calls, with what
data, what it returns) · **license cost**. Then the **decision table** — one row per probe:

| probe | possible outcomes | favours | moves cost line |
| --- | --- | --- | --- |
| load (5.4) | loads clean / loads after shim / API break | O1–O3 / O1–O3 with a maintenance line / O2 | "xPand maintenance" |
| A1 headline (5.6) | order-2 split + η-variation works / only the pre-split `VarD` route works / neither | O1 as full engine / O1 for the split + own E-L / O2 | "E-L after the split" |
| A2 control (5.5) | mixed-index equation `proved-equal`; lowered-index residual = background equation | validates the gate mechanism for #501 | "background-consistency gate" |
| B + CAMB check (5.7) | four `proved-equal` in both dictionaries / attributed sign flips / disagreement | conventions carriage rule confirmed / rule amended / stop | "CAMB seam transcription" |
| C signature (5.8) | native / split native but projections and fluid norm by hand / not native | O1 native / O1 + 2-site transcription / transcribe at export | "signature transcription" |
| D torsion, xPand route (5.9) | representable with the precedent work / breaks at derivative split / breaks at background torsion | O1 or O3 / O3 / O2 or O4 | "torsion SVT rules", "homogeneous mode" |
| D action route (5.9 D6) | tadpole + TT equations extracted / not | O1 (action route lives) / O2 | "action route for torsion" |
| E xMAG (5.10) | `ToDistortion` + `DefConnectionPerturbation` work with xPand / xMAG loads but does not compose / does not load | O3 / O1 with own distortion step / O1 | "connection → LC + distortion" |

D-C options costed, in the order they are preferred (integrate existing tools first, modify a
tool only as a last resort): **O3** hybrid — xMAG `ToDistortion` for the connection step +
xPand for 3+1/SVT + our SVT rules for the rank-3 field; **O1** xPand with the legacy idiom
(already ported for Stage 1) doing the rewrite instead of xMAG; **O4** xMAG's own
`StartInducedDecomposition` as the 3+1 engine if E(f) shows it is an independent route; **O2**
xPert-only, own 3+1/SVT (legacy-style port, largest build cost); **O5** legacy
`ComponentDecompose` port (rejected on design grounds, costed for completeness); **O6** fork
xPand to carry a torsionful connection natively — **last resort**: a fork drifts behind
upstream (the `slegner/CAMB` precedent, #498) and every user must install a non-standard
package; costed with that maintenance and the GPL consequence below. **Two consequences to
state for every option:** (i) the post-Riemannian rewrite is performed in both branches, so
whichever tool does it (legacy idiom ported, xMAG, or a fork) binds Stage 1 as well as M3 — the
memo says whether the two branches share it or duplicate it; (ii) **whichever packages D-C
adopts become install requirements of the certified configuration** — promoted installers
under `scripts/install-*.sh`, their versions added to `verify-wolfram-setup.sh`'s fingerprint,
a `WOLFRAM_GUIDE.md` step and a re-certification — recommended in Part 4 and routed to the
orchestrator (those paths are not owned by this lane).
**License as a cost line:** the repository is MIT; xAct itself is GPL and is already loaded at
derivation time (PSALTer likewise), so *calling* xPand (GPL-3) or xMAG/xBrauer/TraceFree
(GPL-2+) at derivation time is the existing footing; *vendoring* any of them into the repo,
*porting* their code into our package, or *forking* xPand (O6) would put the derivative under
the GPL — each option states which it needs.

## 5. Wolfram execution plan (lane ≤ 2 days; one kernel; serial)

### 5.0 Discipline for every kernel step

`run_lane.sh <step>` runs one `.wls`: refuses if
`pgrep -x WolframKernel || pgrep -x wolframscript` finds a kernel; exports `QT_QPA_PLATFORM=offscreen`; runs from `mktemp -d`;
`timeout --signal=INT --kill-after=60 <budget>`; tees a scrubbed transcript (`$HOME` → `~`,
repo root → `<repo>`) to `third_party/perturbations_runs/<step>/`; records wall time. Every
`.wls` loads `RCSetup.wl` (versions banner `RC_VERSIONS`, geometry with a `normu` parameter,
fluid, `RCVerdict`, `RCDigest`, WXF export) and prints `RC_<STEP>_DONE`. Success is read from
sentinels and artifacts. Identifiers ASCII; xPand's Unicode heads (Φ, Ψ, ℋ display) are named
as the load probe prints them.

### 5.1 Gates before (kernel ~4 min)

`pgrep` → `bash scripts/psalter/ensure_registered.sh` →
`bash scripts/verify-wolfram-setup.sh --require-psalter`; exit 0 required. Record `xPert.m:1`. `manifest.sh before` (§8.1).

### 5.2 First commit + draft PR (no kernel)

Plan archive, README, installers, `manifest.sh`, `run_lane.sh`, `RCSetup.wl`, `compare.wls`,
`targets.wl` (typed from the TeX lines of §2.3 and from the derivation in §5.5, **fetched
before writing it**), `mb_camb_symbolic.py` + its printout. Rehearse `install_xpand.sh --dest`
against a mock tree in the scratchpad (refusal, junk stripping, additivity diff, stamp) with no
kernel. `git add` → `uv run pytest tests/test_repo_hygiene.py` → cspell pre-check → commit →
`gh pr create --draft …` → PR number.

### 5.3 `install_xpand.sh` (no kernel)

Resolve the userbase without a kernel: `${WOLFRAM_USERBASE:-$HOME/.local/wolfram/userbase}`
asserted by `Applications/xAct/xTensor/xTensor.m` (fallback
`wolframscript -code '$UserBaseDirectory'`). Refuse if `Applications/xAct/xPand/` exists. `curl -fsSL` the tarball
to `mktemp -d`; assert `sha256 == 26e7abca…` (pinned now); extract with `._*` and `.DS_Store`
stripped; assert `xPand/xPand.m` and `xPand/Kernel/init.m`; `cp -a` to
`Applications/xAct/xPand/`; write `INSTALLED_VERSION` (`0.4.4`, source URL, sha256, the
`$Version` line, UTC, `installed_by`, `patch=none`); print the `Applications/xAct` listing
before/after and assert only `xPand/` was added. **Wayback:** request a snapshot with
`curl -s "https://web.archive.org/save/http://www2.iap.fr/users/pitrou/xPand_0.4.4.tar.gz"`
(outward-facing, asked for in review round 1) and record the returned snapshot URL beside the
sha256 in the memo and the stamp. `--print-removal` prints (never runs) the removal line.

### 5.4 `probe_load.wls` — xPand 0.4.4 on the certified bundle (kernel ≤ 10 min)

`Needs["xAct`xPand`"]` under `Check`; print `RC_XPAND_VERSION`, `RC_XTENSOR_VERSION`,
`RC_XPERT_VERSION`, `$MessageList`; print the usage messages of the functions §2.2 names.
Geometry: `DefManifold[M, 4, {α,β,μ,ν,λ,σ}]`, `DefMetric[-1, g[-α,-β], CD]`,
`DefMetricPerturbation[g, dg, ε]`, `SetSlicing[g, n, -1, h, cd, {"|","D"}, "FLFlat"]` (the
0.4.4 signature, `xPand.m:1710`); read back the label layout (`RC_OUT47_LABELS`) and the
definitional lines `RC_DEF_dg_TT`, `RC_DEF_dg_SS`, `RC_DEF_u_S`, `RC_DEF_rho` from
`SplitMetric[g, dg, h, "NewtonGauge"]` and `ToxPand[u[-a], …]` — these fix the MB dictionary
before any comparison. Then the paper's minimal example (`:1350-1360`) compared with Out[47]
(`:1368-1373`, 𝒦 = 0 on FLFlat). **Success:** no messages, `identical`/`proved-equal`, control
`proved-different`. **Failure branch (Q3):** classify (`$MessageList`, `Stack[]`) as version
check / renamed symbol → shim inside `xPand/` (diff committed, stamp updated, upstream-issue
draft written, not filed) or as API break → stop, documentary row, xPert-only fallback.

### 5.5 A2 `repro_a2_tensor_eom.wls` — control from the field equations (kernel ≤ 15 min)

**Why mixed indices (verified by hand, then in the kernel).** With `G_ij = g_ik G^k_j`,
`δG_ij = δg_ik Ḡ^k_j + ḡ_ik δG^k_j`. On FRW in conformal time `Ḡ^k_j = −(2ℋ′+ℋ²)a^{−2}δ^k_j`
and for a TT mode `δg_ij = a²h_ij`, so `δG_ij = a²δG^i_j − (2ℋ′+ℋ²)h_ij`, while
`δT_ij = a²P̄ h_ij`. The extra terms cancel only with the background equation
`2ℋ′ + ℋ² = −κa²P̄`. In vacuum with `a(η)` unspecified the lowered-index comparison **fails for
a physics reason**; `a²δG^i_j = ½(h″ + 2ℋh′ − ∇²h)^i_j` carries no background term (the paper's
Appendix A uses `MyGR[μ,-ν]`, `:1770`). Same at second order: the EH quadratic TT action has
`h²` terms proportional to the background equations unless matter or Λ is included and the
order-0 equations are used. **Choice:** include Λ and a perfect fluid (as `:1757-1766`), solve
the order-0 split for `{ρ̄, ℋ′}` (`bgRules`), and run three things: (i) mixed indices →
`proved-equal` with the target; (ii) lowered indices with `bgRules` → `proved-equal`;
(iii) lowered indices without `bgRules` → the printed residual `∝ (2ℋ′+ℋ²+κa²P̄) E_ab`, reported
as the first-order picture of the gate mechanism (#501): on a background that does not solve
the theory, this is what fails to cancel.
Script: `rules = SplitMetric[g, dg, h, "AnyGauge"]` with scalars and vectors zeroed
(`Φ, Ψ, Bs, Bv, Es, Ev → 0`, fluid perturbations → 0);
`ToxPandFromRules[MyGR[a,-b], rules, h, 1]`; `ExtractComponents[ExtractOrder[ah[]^2 res, 1], h, {"Space","Space"}]`. **Target**
(typed, `targets.wl`, cited to MFB eq. 4.15 and MB `:644-661` for the scalar sector):
`Et[LI[1],LI[2],-a,-b] + 2 Hh[] Et[LI[1],LI[1],-a,-b] − cd[-c]@cd[c]@Et[LI[1],-a,-b]` up to
a numeric normalization the harness reports. **Negative controls:** friction `2ℋ → 3ℋ`; `D²`
dropped.

### 5.6 A1 `repro_a1_tensor_action.wls` — the headline: action → order 1 → order 2 → equation

The programme derives from Lagrangians, so this is M3's real path; it is **measured, not boxed
away** — every obstruction is timed and becomes a Part-4 cost line.
Theory: Einstein–Hilbert + Λ + a canonical scalar `χ` with potential `V(χ)` (a genuine matter
action, so the background is on-shell without a fluid action):
`L = √−g [R/(2κ) − Λ/κ − ½ ∇χ·∇χ − V(χ)]`, background `χ̄(η)`, `a(η)` unspecified, `V` an unspecified function.
Steps: (1) `L2 = ToxPandFromRules[L, Join[rules, {dχ[LI[o_]] :> X[LI[o]]}], h, 2]` (falls back
to `ExpandPerturbation@Perturbed[…, 2]` + `Conformal` + `SplitPerturbations` if the density
is not handled in one call; `√−g → a⁴√−g_conf` by hand if `Conformal` does not take `Detg`);
(2) **order 1 = the tadpole**: `ExtractOrder[L2, 1]` selects the ε¹ coefficient of the
expanded action (it removes nothing from the ε² equations, which are derived separately in
step 3); that coefficient must be a combination of the Friedmann equations and the `χ̄`
equation times the perturbations — printed as `RC_A1_TADPOLE`, then evaluated with the
order-0 rules substituted and checked to vanish (`RC_A1_TADPOLE_ONSHELL=0`); this is the
quantity the #501 gate exports and evaluates on the CAMB background for a user theory; (3) order 2, TT sector, `bgRules` applied → the quadratic TT action
`∝ a²[(E′)² − (DE)²]`; (4) **Euler–Lagrange after the split**: levels `Et0, Et1, Et2` replace
`Et[LI[1],LI[n],…]`, spatial parts by `VarD[Etk[-a,-b], cd]`, time parts by a linear `dEta`
(Leibniz; `dEta[ah] = ah Hh`, `dEta[Hh] = Hh′`; `dEta[h[..]]` from a read-back of
`ToxPand[LieD[n][h]]`), `EL = P0 − dEta[P1] + dEta[dEta[P2]]`; **scalar control first** (the
massless `χ` TT-free sector must give `X″ + 2ℋX′ − D²X`, the `de_sitter_kg` family with D = 4);
(5) compare with A2's result and the target. Also run the pre-split route
(`VarD[dg[LI[1],-a,-b], CD]` on the second-order action, then split) as the cheap cross-check; the two must be
`proved-equal`. Budget: 3.5 h wall on Day 1 with a 45-min kernel timeout per attempt; if step 4
is not `proved-equal` by then, the exact expression shapes and the time spent are the
measurement, A2 + the pre-split route stand as the action-route result, and the cost line is
written from that.

### 5.7 B `repro_b_mb_scalars.wls` — Ma & Bertschinger, machine-checked twice (kernel ≤ 30 min)

Flat FRW, `DefMatterFields[u, du, h]`, `T = (ρ+P) u u + P g` (+ a traceless projected stress
`Π_ab` to exercise σ, B2), `MyGR[a,-b]`, `ToxPand[…, "NewtonGauge", 1]`; order 0 → both
Friedmann equations (MB `:573-580`); order 1 → `{Time,Time}`, `{Time,Space}`, trace and
traceless `{Space,Space}`, vectors and tensors zeroed. **Dictionary** from the read-back lines
of 5.4: MB ψ ↔ xPand Φ (lapse), MB φ ↔ xPand Ψ (curvature), δρ, δP, v_i ↔ the fluid heads,
`4πG = κ/2`, `k² → −D²`, `θ = ik^jv_j ↔ D^jv_j`, σ ↔ the `k̂_i k̂_j` projection; 23b and 23d are
compared after applying `D^i` and `D^iD^j` to the 0i and traceless-ij equations. Both sides
`(−,+,+,+)` here, so residual signs are variable definitions or Fourier factors, each attributed.
**Second check, in our convention:** `mb_camb_symbolic.py` prints `newtonian_gauge(constraints)`,
`newtonian_gauge(pert_eqs)`, `newtonian_gauge(total_eqs)` and the synchronous forms from the
installed `camb.symbolic` (declared `(+,−,−,−)`, `symbolic.py:461`); the printout is committed;
`targets.wl` transcribes those lines (cited) with `Φ_N ↔ xPand Ψ`, `Ψ_N ↔ xPand Φ`,
`dgrho = a²κδρ`, `adotoa = ℋ`, `k → −D²`; the kernel compares — this is the check `conventions.md`
§5–§6 assigns to the MB → CAMB conversion, and it settles by machine that mixed-index scalar
equations are invariant under `g → −g` (rule 1) with identical variable definitions (rule 2).
**Synchronous gauge too:** the same script with `"SynchronousGauge"`, compared with MB's
synchronous equations (`:603-618`, `ein-syna..d`, variables `h, η`) and with `camb.symbolic`'s
`synchronous_gauge` output and its CDM frame, so both CAMB-named gauges are exercised and the
memo shows the interface for each. **Success:** four `identical`/`proved-equal` per dictionary
and per gauge + one negative control each.

### 5.8 C `probe_c_signature.wls` — hypothesis from source, then the kernel (kernel ≤ 15 min)

Hypothesis (§2.2): `SetSlicing[g, n, +1, h, cd, …]` with `DefMetric[-1, …]` is accepted (only
the det sign is checked) and `SplitPerturbations` is sign-generic (the projector uses
`v v / norm`), but `ExtractComponents`'s "Time" projection (`:3042`) and `ToxPand`'s fluid norm
(`:3014`) hard-code `n·n = −1`. Test: rerun the minimal example with `normu = +1`; compare the
scalar `SplitPerturbations` result to the `−1` result under `R → −R` and the variable flips of
`conventions.md` §6 rule 2; then project a vector once with `ExtractComponents` and once by
hand with `n[-a]` and report the sign. **Verdict:** native / split native + two hand-done
sites / not native.

### 5.9 D `probe_d_torsion.wls` — torsion on FRW, designed from the precedent (kernel ≤ 60 min)

D1 declare `DefTensor[T[a,-b,-c], M, Antisymmetric[{-b,-c}]]`,
`DefTensorPerturbation[dT[LI[o],a,-b,-c], T[a,-b,-c], M]` beside `dg` (legacy makes the same pair, `_derive.py:2350-2376`).
D2 **24-count** by xCoba components (`DefChart`, `ToCanonical` of all `T[{i},{-j},{-k}]`,
`Length@DeleteDuplicates` of the nonzero ones modulo sign = 24).
D3 **SVT rules transcribed from 2310.16007 `:768-776, :798-807, :834-842`** (with `N → 1`,
`γ_ij → h_ij`, `𝒟 → cd`, their `ε_ijk → n^d epsilong[-d,-a,-b,-c]`): 8 scalars {T,B,φ,A} ∪
{𝒮,ℬ,ϱ,𝒜}, 6 transverse vectors, 2 TT tensors, each a `DefProjectedTensor` of rank ≤ 2 (which
is what xPand's property machinery accepts); mechanical checks: antisymmetry of the rule's
right-hand side, 16 symbols, 8 + 2·6 + 2·2 = 24 (cross-checked against 1804.09215 `:819-852`).
D4 **projection**: `ExtractComponents` over all eight `{Time|Space}³` patterns must evaluate
(supported by design, `:3031-3044`).
D5 **homogeneous background torsion mode (optional, ≤ 30 min, run after D6–D8 if time
remains; outside the planned scope — Terms)**: declare the FRW-symmetric
background `T̄^a_{bc} = 2 T1[η] n_[b h_c]^a + 2 T2[η] ε^a_{bcd} n^d` with `T1, T2` background-only
projected scalars (`SpaceTimesOfDefinition -> {"Background"}`, like `a[h]`), as a background
rule in the list `SplitPerturbations` accepts (`ProjectionAndBackgroundRuleQ`, `:2891`); run
`T·T` and `∇_dT` at orders 0 and 1 and assert the split carries `T1, T2` and their primes
through; the order-1 term with `T̄ ≠ 0` is the tadpole the #501 gate would evaluate — printed.
If the background rule is refused, the exact refusal is the verdict and the gap goes to the
not-tested table with its Part-4 consequence.
D6 **action route for torsion**: a torsion-quadratic Lagrangian
`L_T = √−g [c₁ T_{abc}T^{abc} + c₂ T_{abc}T^{bac} + c₃ T_a T^a]` with `T̄ = 0` (then with D5's `T̄`) → `ToxPandFromRules[…, 2]`
→ order-1 term (zero for `T̄ = 0`, nonzero for `T̄ ≠ 0`, printed) → order-2 TT sector → the
5.6 Euler–Lagrange operator → the equations for the two TT torsion pieces; assert the helicity
sectors decouple and the TT equations are diagonal in the two pieces. This is what M3 needs,
not a decoupling check of `T·T` alone.
D7 **derivative split**: `cd`/`CD` of `dT` through `SplitPerturbations` must leave no
residual `CD` and only `D` and primes; one hand-checkable piece (`∇_0 T^0_{0i}` on FRW =
`∂_η T^0_{0i} − ℋ T^0_{0i}`) compared by the harness.
D8 **torsionful connection**: legacy's
`DefCovD[CDT, {"#","DT"}, FromMetric -> g, Torsion -> True]`; `Quiet@Check[ToxPandFromRules[RicciScalarCDT[], …], $Failed]` → record the failure mode
(the authors' stated exception); then legacy's rewrite `ChangeCurvature[…, CDT, CD]` + the
contortion → torsion rule (`_derive.py:2164-2189`) → rerun D7 on the result; identity check:
`δR̃ − δR` equals the split of the linear torsion term.
**Verdicts** (tensor field, background mode, connection — separately): representable /
representable with work (naming (a)–(c) of §2.6 and the sites) / not representable; each
pointing at a sentinel; an unreached check is reported unreached.

### 5.10 `install_xmag.sh` + E `probe_e_xmag.wls` (Q4; box half a day; kernel ≤ 45 min)

Installer: 5.3's discipline; pinned `git fetch --depth 1 <sha>` of
`THelpin/xBrauer_Bundle@48be67e1…`, `xAct-contrib/TraceFree@4e53ab39…`, `THelpin/xMAG@88026e47…`; refuses if any of the
five targets exists; copies `SymmetricFunctions/`, `BrauerAlgebra/` → `Applications/`;
`xBrauer/`, `TraceFree/` (repo root is the package), `xMAG/` → `Applications/xAct/`;
`INSTALLED_COMMIT` in each; `--print-removal`. **No patch to xBrauer**; an xMAG shim follows Q3.
Probe: (a) `Needs["xAct`xMAG`"]` — versions, messages (expects xTensor 1.2.0 / xTras 1.0.6
minimums, satisfied); stop per Q4 if it demands other versions or writes outside its
directories. (b)
`DefCovD[CDt[-a], {"~","∇~"}, Master -> g, FromMetric -> Null, Torsion -> True]`; inspect `Distortion`, `Contorsion`, `TorsionVector`. (c)
`ToDistortion[RicciScalarCDt[], CD, g]` → assert no `CDt`-curvature remains. (d) `DefMetricPerturbation[g, dg, ε]` +
`DefConnectionPerturbation[ChristoffelCDt, dΓ, ε]` — does xPert accept both? (e) end-to-end:
`ToDistortion` of D6's invariants written with `CDt` → perturb to order 2 → xPand split with
D3's rules for the contortion perturbation → assert the split runs and matches D6. (f) one call
of `StartInducedDecomposition` with its usage printed — an independent 3+1 route or not.
**Verdict** as in 5.9, plus license (GPL-2+ headers; no LICENSE files) and last release
(commits 2026-01-10 / 2013-11-12; no tags) in the Part-3 row.

### 5.11 Legacy comparison (Part 3 item 5; no kernel)

Capability table (§2.7) beside xPand/xMAG; legacy's own torsion declaration shown by
`uv run pytest tests/test_cli.py -k torsion_covd_dry_run`; the frozen `de_sitter_kg` equation
against 5.6's scalar control (same action, same volume factor, `(D−2)ℋ` friction). **No
`tidal derive`.**

### 5.12 Gates after + manifest (kernel ~4 min)

`bash scripts/verify-wolfram-setup.sh --require-psalter` exit 0; `manifest.sh after`; layer (i)
identical, layer (ii) per §8.1; both transcripts scrubbed into the memo's evidence index.

### 5.13 `compare.wls` and the verdict rule

Modeled on `compare_waveoperators.wls`: `SameQ` →
`TimeConstrained[ToCanonical@ContractMetric[expr − target], 120]` → `Simplify` fallback
(120 s, one retry at 300 s) →
`identical / proved-equal / proved-different / could-not-decide`; digests = SHA256 of the canonical form
after `ScreenDollarIndices`, printed twice, never the verdict; exit 0 only if every result is
`identical`/`proved-equal` **and** every control is `proved-different`. **`could-not-decide`
counts as not reproduced** and is reported as such.

### 5.14 Budget and abort criteria

| step | kernel (est.) | wall incl. writing | abort / branch |
| --- | --- | --- | --- |
| 5.1 gates before + manifest | 4 min | 15 min | exit ≠ 0 → stop, report |
| 5.2 commit + PR | 0 | 1 h | — |
| 5.3 install xPand | 0 | 15 min | sha/shape mismatch → stop |
| 5.4 load probe | 2–10 min | 30 min | Q3 branch |
| 5.5 A2 control | 5–10 min | 45 min | > 30 min kernel → reduce |
| 5.7 B + CAMB check (both gauges) | 15–25 min | 2.25 h | dictionary read-back missing → do not compare |
| 5.6 A1 headline | 20–60 min | **3.5 h, measured** | see 5.6 |
| 5.8 C | 5–10 min | 30 min | — |
| 5.9 D (D5 optional, last) | 30–60 min | 2.5 h | each sub-check reports reached/unreached |
| 5.10 xMAG install + E | 30–45 min | **half-day box** | Q4 stop conditions |
| 5.12 gates after + manifest | 4 min | 15 min | exit ≠ 0 → investigate first |
| memo, references rows, README, issues, CI | 0 | 6 h (lane released) | — |

Day 1: 5.1 → 5.2 → 5.3 → 5.4 → 5.5 → 5.7 → 5.6 → 5.8 (kernel ≈ 1.5 h). Day 2: 5.9 → 5.10 →
5.12; the lane is released at 5.12. The memo, the `docs/references.md` rows, the issues and CI
need no kernel and may finish after the lane is released. Literature TeX is fetched before
`targets.wl` is written (5.2), not while kernels run.

## 6. Findings routed to the orchestrator (reported, filed as issues, never edited here)

- F1 `R-C.md:30, :57` cite `COSMOLOGY_PROGRAM.md:302-310` for the H7 xPand paragraph; it is at
  `:317-325`.
- F2 `R-C.md:62` calls `solver_design.md:330-352` "§2"; it is under §4, and §5 (`:461-506`)
  is the section that owes WS2 the eikonal contract.
- F3 `primer.md:130, :251` use mostly-plus and call φ a potential; `conventions.md` fixes
  `(+,−,−,−)` and φ = Weyl potential.
- F4 `spectrum_design.md:292-294` still says the solver branch may adopt Ma–Bertschinger
  conventions; the §4.3 heading still asserts the withdrawn "opposite signatures" premise.
- F5 "standard sectors enter only as sources" (`R-C.md:45-47`) vs "coupled standard mode is
  inside the block" (`repo_reshape.md:233-237`) — Part 1 states the split.
- F6 the `t = 0.0` energy bug is "critical path" (`observable_ladder.md:138`), "moot"
  (`repo_reshape.md:800-803`) and "must be fixed" (`:835-837`).
- F7 `gauge_certificate` dropped in `repo_reshape.md:808-810`, live in
  `solver_design.md:312, 547-548`.
- F8 O4-iso needs a homogeneous new-sector mode yet the admissible-theories row scopes to the
  assumed background; the gate must be specified so it covers that mode (Probe D5).
- F9 SymBoltz.jl pinned by commit with no repository URL in R-1's ledger.
- F10 "the fingerprint checks the four existing packages" — the memo names them.
- F11 `R-C.md:65-66` points at `background_fields.tex` for legacy's perturbation pipeline; the
  right file is `multi_field_perturbation.tex`.
- F12 FLRW (R-C.md) vs FRW (everything else).
- F13 `docs/references.md` has no xPand row (added under Q2).
- F14 xPand's authors state torsion theories are the exception — the H7 near-miss table
  (`spectator_route.md:137`, "any theory") overstates.
- F15 the xPand paper's Appendix A prints `−2D²φ` (`:1789`) and `E′ℋ` (`:1809`) where the
  textbook has `D²ψ` and `2ℋE′`; the kernel decides what 0.4.4 prints (memo finding).

## 7. Risks

| risk | mitigation |
| --- | --- |
| xPand API drift on xTensor 1.3.0 (the version gate cannot fire, §2.1) | `probe_load` prints `$MessageList` and the failing symbol; Q3 branches pre-decided |
| second perturbation structure beside `dg` (`dT`, `dΓ`, `dχ`) | legacy precedent (`_derive.py:2350-2376`); B2's `Π_ab` rehearses rank 2 before D needs rank 3 |
| Euler–Lagrange over Lie-derivative labels is awkward (5.6) | scalar control first; time measured, obstruction recorded as the cost line; pre-split `VarD` route and A2 stand |
| `ToCanonical`/`Simplify` slow at order 2 | TT-only rules; `TimeConstrained` (120 s, one 300 s retry); `could-not-decide` = not reproduced |
| front-end hang | `QT_QPA_PLATFORM=offscreen` and `timeout` on every launch |
| hygiene test rejects pasted `/home/...` | print-time scrub; `git add` then `uv run pytest tests/test_repo_hygiene.py` |
| cspell strict | in-file directives; American English; no British quotes from abstracts |
| guard blind spots (bare `.wls`, worktree path resolution) | `run_lane.sh` and both installers self-guard with `pgrep -x WolframKernel` |
| the known answer hides a background residual (lowered indices, off-shell background) | mixed-index runs as the clean target; lowered-index residual reported deliberately (5.5); order-0 rules solved from the split, never typed |
| `SameQ` across kernels on abstract-index expressions | `proved-equal` via `ToCanonical` is the primary tier; digests via `ScreenDollarIndices` |
| lane overrun | budget table (§5.14); abort criteria per step; memo written after the lane is released |
| the record is thinner than the work | Appendix A/B/C archived; memo's verified/not-verified ledger; provenance README; issues per discovery (§3) |

## 8. Decisions made by the user at planning (2026-09-15)

- **Q1 — plan archive: yes.** Archive the approved plan verbatim as
  `docs/cosmology/rc_planning_record.md` in the first commit; owned paths extended by exactly
  that file. R-1 precedent: header "ARCHIVE — verbatim … Not authoritative" naming the living
  documents (`perturbation_tooling.md`, the register in `COSMOLOGY_PROGRAM.md`); body
  byte-identical below the rule (checked with `diff`); spelling fixed in the plan first;
  redactions disclosed; hygiene test run after `git add`.
- **Q2 — literature: fetch + curated rows, no README regeneration.** Fetch only what the memo
  cites for a claim into the worktree's gitignored `literature/`; read TeX, not abstracts. Add
  curated rows to `docs/references.md` under `## Cosmology Program` (line 155), including the
  missing xPand row; `docs/references.md` joins the owned paths. **Do not touch
  `literature/README.md`** (`index_literature.py` rebuilds from disk; a fresh worktree has none
  of the ~163 papers and would silently delete every other row). List the fetched arXiv ids in
  the report; the orchestrator copies the directories into the main checkout and regenerates.
- **Q3 — xPand patch: shim only, three conditions.** A version check or small compatibility
  shim inside `Applications/xAct/xPand/` only; a real API break against xTensor 1.3.0 / xPert
  1.0.6 → stop patching, record the exact failure, documentary row + xPert-only reproductions,
  report the Part-4 cost change. Evidence: unified diff committed under
  `scripts/research/perturbations/`, stamped in `INSTALLED_VERSION`, stated in the memo;
  `verify --require-psalter` exit 0 before and after; and a sha256 manifest of every userbase
  file outside the new directories, before and after, shown identical (verify fingerprints
  only four packages). Draft the upstream issue, do not file it (precedent
  `docs/cosmology/psalter_543_upstream_issue.md`), under `scripts/research/perturbations/`.
- **Q4 — xMAG: extend the fence, install additively and probe.** Rationale (user): xMAG is the
  one package whose documented scope is a torsionful connection; "tool X can represent torsion"
  is a script, not a sentence; a documentary row would leave the most relevant candidate
  untested exactly where the recommendation rests. Conditions: new directories only (all five
  verified absent); nothing existing modified; **no patch to xBrauer**; the one manifest covers
  both installs and must show only the new directories were added; stop at the load probe if
  xMAG or xBrauer need different xTensor/xPert versions or would write outside their own
  directories (row becomes documentary with the exact failure); record license and last
  release; timebox half a day inside the ≤ 2-day lane; report every directory added.

### 8.1 Manifest rules, fixed before the before-manifest is taken

Userbase top level (listed 2026-09-15): `.activation_backup/`,
`ApplicationData/{CloudObject,Credentials}/`, `Applications/xAct/` (939 files), `Autoload/PacletManager/`,
`FrontEnd/{init.m,14.3_Caches}/`, `Kernel/init.m`, `Licensing/mathpass`,
`Paclets/{Cached,Configuration,Repository,Temporary}/` (2639 files), `SearchIndices/`, `SystemFiles/`.
- **Layer (i) — must be byte-identical before/after:** `Applications/**` minus the six new
  directories, `Kernel/**`, `Licensing/**`, `Autoload/**`, `SystemFiles/**`,
  `.activation_backup/**`.
- **Layer (ii) — allowed to change, reported file by file:** `Paclets/**` (paclet manager
  state and caches, modified by every kernel start), `FrontEnd/**` (front-end caches and
  `init.m`), `SearchIndices/**`, `ApplicationData/CloudObject/**` and
  `ApplicationData/Credentials/**` (cloud-session state that the kernel maintains on its own;
  we never touch it — the login is left alone).
- **Any difference outside both lists is a finding**, never something classified afterwards.
  The six new directories are listed by name in the "after" report.

## Appendix A — revision log (planning)

### Round 0 (2026-09-15, exploration and first submission)

- Three read-only surveys (design docs → requirement lines and 14 routed findings; lane
  mechanics from the R-1 precedent; legacy decomposer, local literature, installed userbase)
  plus a targeted read of Ma & Bertschinger's TeX; external facts pinned (xPand 0.4.4 and its
  user list, the xAct family, SymBoltz.jl, `camb.symbolic`, comparables' licenses, INSPIRE
  citation data, the torsion-on-FRW papers).
- User decisions Q1–Q4 (§8); Q4 overrode the recommendation (xMAG probed, not documented).
- Post-decision research: xMAG's dependency chain resolved, commits pinned, targets verified
  absent, xMAG's API read from source.
- A Plan-agent review of the kernel sequence was still running at submission (item 1 of
  round 1); its report arrived with the review and is reconciled below.

### Round 1 (2026-09-15, the plan was returned with eleven points)

| # | proposed (round 0) | feedback | research done | outcome |
| --- | --- | --- | --- | --- |
| 1 | Appendix C deferred the tarball, signature, rank-3, CAS, xPPN reads to the lane | plan mode is read-only, not light on research; settle everything reading can settle | tarball downloaded (sha256, layout, version fields, gate mechanics, `SetSlicing` signature, 0.4.3→0.4.4 diff, docs via `pdftotext`); `ExtractComponents`, `SplitPerturbations`, `ToxPand`, `DefProjectedTensorProperties` read; xTensor's `InducedFrom` read; TeX of the ten papers fetched and grepped for CAS; xPPN found on GitHub; xMAG contradiction fixed; Plan-agent report reconciled | §2.1–§2.3, §2.6 rewritten; Appendix C holds only lane items |
| 2 | Probe D designed from first principles | read how prior work carried torsion through xPand/xAct first | 2009.02168 (tetrad matrix + xPand), 2012.11959, 2206.00658, 2310.16007 (parametrization lines, table), 1804.09215 (8+6+2 by hand, xAct), 1601.03943 (Maple, two torsion modes, nonzero background) | §2.6 precedent paragraph; 5.9 D3 transcribes 2310.16007's rules; "work" = (a)–(c) |
| 3 | A2 compared lowered-index `δG_ij` in vacuum on `a(η)` | lowered indices carry `(2ℋ′+ℋ²)E` terms cancelled only on-shell; use mixed indices or matter + background equations; same at second order | derived by hand in 5.5; the paper's Appendix uses `MyGR[μ,-ν]` | 5.5 runs mixed, lowered-on-shell and lowered-off-shell (residual reported); 5.6 uses EH + Λ + scalar with `bgRules` |
| 4 | action route a 3-h side box; A2 headline | the programme derives from Lagrangians: A1 with the order-1 extraction is the headline, measured not boxed; Probe D on the action route too | — | 5.6 rewritten as the headline with the tadpole step and a scalar control; D6 added |
| 5 | background torsion not probed | `T̄ ≠ 0` is a requirement (F8, 2310.16007) | 2310.16007's FRW torsion background, 1804.09215's `f, g`, 1601.03943's `h, f` | D5 added; Part-1 line amended; not-tested table if refused |
| 6 | MB → our convention as a written table; targets by page; `could-not-decide` unstated | machine-check with `camb.symbolic`; TeX-line citations; `could-not-decide` = not reproduced | camb 2.0.4 found in `.venv`; `constraints`, `pert_eqs`, `Newtonian_subs` read (`symbolic.py:355-452`); 1302.6174 TeX lines pinned; two paper-vs-textbook discrepancies found (F15) | 5.7 second check + committed printout; §2.3; 5.13 rule |
| 7 | Part 1 lines without consumers | say whether CAMB's API, WS3 or WS4 consumes each line | installed CAMB grepped: only `set_custom_scalar_sources` (`model.py:1111`), no tensor/eikonal hook; `get_time_evolution(frame=)` | §4.1 consumer column; tensor `S_new` and the eikonal export marked ours |
| 8 | Part 4 without a decision table; three D-C options; no license cost | probe → outcomes → option → cost line; add xMAG routes; license per option | licenses pinned (xPand GPL-3; xMAG/xBrauer/TraceFree GPL-2+; CAMB LGPL-3+GFDL; Cadabra GPL-3; CLASS/hi_class no license file) | §4.4 decision table, O1–O5, license cost paragraph |
| 9 | terms defined by reference | define each plainly where it first appears | — | Terms section rewritten with meaning + relevance |
| 10 | manifest layer (ii) classified after the run | fix the allow-list by directory name beforehand | userbase listed to depth 2; recently modified files inspected | §8.1 |
| 11 | memo length unstated; risk table cited §5.9; targets while kernels run; no Wayback; review log missing | ≤ ~600 lines; fix reference; fetch TeX before `targets.wl`; Wayback snapshot; log the round | Wayback checked: no snapshot exists → SPN request at execution | §3, §5.2, §5.3, §7, this table |

**Plan-agent reconciliation (its report arrived with round 1).** Adopted: the version-gate
mechanics (verified independently in `xPand.m:94,96`), mixed indices with the lowered-index
residual as a deliberate output, `DefMetric[-1,…]` + `normu` for the signature probe (its
"never `DefMetric[1,…]`" point stands: the first argument is the determinant sign), the
paper-syntax read, the legacy `ChangeCurvature` idioms (`_derive.py:2164-2189`), the
`de_sitter_kg` documentary bridge (`(D−2)ℋ`), the `ToCanonical`-first verdict with
`ScreenDollarIndices` digests, the scalar control before the tensor Euler–Lagrange, the
two-phase sha256 (now pinned directly), the B2 anisotropic-stress step, the guard's worktree
blind spot. Declined: an executable `uninstall_xpand.sh` (fence: no removal is run; a
`--print-removal` line suffices and the orchestrator's post-merge check covers the directories);
its claim of a 1.1.5 `SeparateMetric` restriction and of two 2018 xPand-specific xTensor fixes
(one entry found, `xTensor.History:1690`; the rest not found — not cited); its four questions
(answered: no uninstall script; shim policy is Q3; no 2+1 stretch unless Day 2 has slack after
E; A1 is measured per 5.6, not extended into D/E time).

### Round 2 (2026-09-15, six comments at the second submission)

| # | comment | outcome |
| --- | --- | --- |
| 1 | the output must not be only the memo — record every investigation and finding with evidence, so future work has them and knows why the decision was made | TL;DR and §3 rewritten around "the memo **and** the record": planning archive with Appendices A–C, the memo's verified/not-verified ledger and evidence index, the provenance README, `docs/references.md` rows, GitHub issues per discovery (F1–F15 and any tool defect), a "record is thinner than the work" risk row |
| 2 | the residual "would drive our modes" — but admissible theories satisfy the background (spectators); we do not test theories that alter the expansion history | reframed everywhere: the tadpole is the **background-consistency gate** (#501) that admissible spectator theories pass, exported to check that, never propagated as a source (`observable_ladder.md:185-193`; `spectator_route.md:81-91`; `COSMOLOGY_PROGRAM.md:200`); "spectator" and "homogeneous new-sector mode" defined; D5 reframed as representing a zero mode of the new sector, subject to the same gate |
| 3 | "Part 4" of what, in which wave? | Terms now define Part 1–4 (this memo, this lane), D-C (recorded by the orchestrator after merge), Wave 1 / M3 / WS2 |
| 4 | what is the solver-only block and why separate? | Terms explain: the shared theory file feeds two derivations; FRW-only inputs would change the spectrum's fingerprint without changing the spectrum, so they get their own block and fingerprint (`interfaces_decision.md` §3.8) |
| 5 | are "standard sectors" neutrinos, dark matter, etc.? | Terms define them (photons, baryons, CDM, neutrinos, dark energy, the metric) and how they enter |
| 6 | §2 too terse to review — say why each thing was investigated and its impact | every §2 subsection now opens with "Why" and closes with "Effect on the plan"; tables gained a "why it matters"/"consequence" column |

### Round 3 (2026-09-15, seven comments at the third submission)

| # | comment | outcome |
| --- | --- | --- |
| 1 | the eikonal export: do not commit to a solver approach early; bear the solver research in mind and what form it demands | new Terms entry "Solver forms": the derivation emits the representations every bake-off candidate consumes (second-order block, first-order `y′ = M(η)y`, and the eikonal amplitude form for the photon channel) and the memo records what each demands without choosing (H3: Magnus and matrix-WKB together, `rk-adaptive` baseline, bake-off decides; `solver_design.md` §5, §7, `:466-470`); TL;DR and the Part-1 exports line reworded |
| 2 | is xMAG good for many torsion aspects? | §2.4 now opens the xMAG entry in plain terms: by documented scope, yes (independent torsionful connection, Levi-Civita + contortion decomposition, connection perturbations and variations, an induced decomposition); but unpublished, 0.1.0, no cosmology functions, untested on our bundle — hence probed, not trusted |
| 3 | "extracted" — not deleting a term? | clarified in the TL;DR, Terms and 5.6 step (2): the ε¹ coefficient of the expanded action is computed and evaluated on the background; nothing is removed from any equation |
| 4 | does a homogeneous torsion background contradict the unchanged-background assumption? | Terms entry rewritten: FRW symmetry *permits* the mode; our scope requires it to leave `a(η)` unchanged (spectator condition) and to satisfy the theory's equations on ΛCDM (#501 gate); for most admissible theories it is zero; it is defined because O4-iso needs one and the register treats that as a user-theory question under the gate; D5 checks representability only |
| 5 | the particle spectrum also needs some of the "solver-only" inputs, e.g. decomposing torsion out of the covariant derivatives | corrected: the post-Riemannian rewrite is performed in **both** branches (`spectrum_design.md:341-351`, listed as duplicated symbolic work at `:372-378`); the solver-only block holds only the FRW background, the FRW field → background + perturbation split, and the gauge (PSALTer takes the gauge-unfixed Lagrangian, `:389-396`); new Terms entry "Post-Riemannian rewrite"; Part 1 and §4.4 record that the decomposition tool binds both branches |
| 6 | xPand primary although its authors exclude torsion? modify it, or integrate with torsion/MAG packages? | §2.1 "Effect" rewritten: xPand is primary for the metric side; torsion enters as a tensor field after the rewrite (how every §2.6 paper and legacy proceed) via a combination with xMAG or the legacy idiom (O3/O1); forking xPand to carry torsion natively added as **O6** with its maintenance and GPL cost; Probes D and E decide |
| 7 | is the papers' by-hand step what legacy TIDAL implemented by hand? | yes — `_derive.py::_wls_torsion_curvature_decomposition` (`:2164-2189`: `ChangeCurvature[L, CDT, CD]` + the contortion identity, citing Shapiro 2002 eq. 2.9 and Hehl et al. 1976), ported for Stage 1 (`stage1_engineering_plan.md:519-525`); stated in Terms, §2.6 and §2.7 |

### Round 4 (2026-09-16, twenty comments at the fourth submission)

| # | comment | outcome |
| --- | --- | --- |
| 1 | standard sectors "never enter the Lagrangian" — but couplings to the photon and metric do, and we replace CAMB's evolution of those | Terms "Standard sectors" rewritten: the Lagrangian carries the new sector and its couplings to the touched standard fields (metric, photon), which are evolved in our block; only untouched sectors arrive as tables and `S_std` |
| 2 | negligible in the Friedmann equation — so a not-exactly-ΛCDM solution is fine if the change is negligible, e.g. a small background torsion | Terms "Spectator" extended: negligible to a derivable tolerance (`spectator_route.md:107-112`), measured by the gate; a small sub-tolerance background torsion is acceptable |
| 3 | a rung needing a nonzero background is outside this project's scope | Terms "Homogeneous new-sector mode" rewritten: O4-iso is outside the planned scope; the gate detects theories that want such a mode and marks them inadmissible; D5 demoted to optional, last, ≤ 30 min |
| 4 | xPand Levi-Civita-only — assistance from xMAG? | yes: stated in the "Post-Riemannian rewrite" term (option O3) |
| 5 | xMAG's distortion functions could extend to non-metricity perturbations | noted in Terms and §2.4: the same functions cover non-metricity, so the route extends to metric-affine theories later |
| 6 | forking xPand should not be the first approach; integrate existing tools; forks drift (slegner) | §2.1 "Effect" and §4.4 reordered: O3 (integration) first, O6 (fork) explicitly last resort with the drift and user-burden cost |
| 7 | conventions in xPand must be consistent with the rest of the project | Probe C's purpose restated: the committed package emits `(+,−,−,−)`, `ε₀₁₂₃ = +1` and CAMB variables either natively or by transcription at export, recorded and asserted (`conventions.md`) |
| 8 | targets = things the symbolic route must reproduce? | yes: sentence added to §2.3 |
| 9 | options for different gauges, CAMB convention first | reproduction B now runs Newtonian **and** synchronous, compared with MB and with `camb.symbolic`'s `synchronous_gauge`/CDM frame; Part-1 gauge line amended; budget +10 min |
| 10 | expansion step first (xMAG), then other packages | that is O3, now the preferred option |
| 11 | xMAG maturity acceptable; file issues or contact the author | §2.4 reworded: not a disqualifier; defects drafted as upstream issues and, with the user's agreement, filed or sent to the author |
| 12 | no single tool does it all; chain packages | §2.4 reworded accordingly |
| 13 | xIST/COPPER is for end users probing theories, not for us | marked "not a candidate component; Part-3 row only" |
| 14 | could xPPN be used on FRW, or is it worse than xMAG? | answered in §2.4 from its usage messages: backgrounds hard-wired flat, velocity-order split, tetrad + flat-connection variables — more limited than xMAG; documentary row |
| 15 | FieldsX likewise a user aid | marked as such |
| 16 | SymBoltz scalars-only is a limiting factor | §2.5 row: scalars-only and Newtonian-only excludes our channels; at most a future consumer for scalar pieces |
| 17 | background torsion solutions — background only, or relevant for turning the decomposition into simulable fields? | §2.6 now separates the two: background solutions are out of scope beyond negligibility; the perturbation decompositions are the recipe D3 transcribes |
| 18 | papers show which tools were used — follow their footsteps | kept; §2.6 "precedent" paragraph is the design basis for Probe D |
| 19 | a thin field is expected; learn parts, do not expect a full method | §2.6 sentence rewritten |
| 20 | if xPand/xMAG prove key, make them package requirements | §4.4 consequence (ii): promoted installers, fingerprint, guide, re-certification — recommended and routed; report-back row updated |
| 21 | line-of-sight quantities — later derivation work too? | yes: Part-1 exports line names the standard tensor source and the re-derivation needed when the photon sector is modified (#514, WS2 → WS4); Part 2 gains item 9 (LOS method literature) |

## Appendix B — evidence ledger (pinned during planning)

- xPand: tarball sha256 `26e7abcac7bb655235ec39b73850729cf4465748d0d5ea2ad03c0607aef5ceab`
  (1,603,346 bytes; no Wayback snapshot yet); `xPand.m` lines cited in §2.1–§2.2; GitHub mirror
  0.4.3 diff; `xPandDoc.pdf` via `pdftotext` (4068 lines); INSPIRE recid 1221019 (129 cites),
  `refersto:recid:1221019` for citing papers (the `refersto:arxiv:` form is **not** applied as a
  filter); xPert recid 790000, xTras 1249377, xIST 1444311, PSALTer 2798424.
- xMAG chain: `THelpin/xMAG@88026e47…`, `THelpin/xBrauer_Bundle@48be67e1…`,
  `xAct-contrib/TraceFree@4e53ab39…`; xPPN `xenos1984/xPPN@dcaabbea…`.
- Papers (TeX in the scratchpad, re-fetched to `literature/` at execution): 1302.6174,
  2310.16007, 2311.05495, 2311.05597, 1601.03943, 1804.09215, 2009.02168, 2012.11959,
  2206.00658, 2012.14984; Ma & Bertschinger local.
- CAMB 2.0.4 in `.venv` (`camb/model.py`, `results.py`, `symbolic.py` lines cited in §2.3, §2.5).
- Userbase listing 2026-09-15 (§8.1); all six install targets absent.

## Appendix C — needs the lane or a write (nothing else is deferred)

- Whether xPand 0.4.4 **loads** on the certified bundle (only a kernel can tell; the version
  gate is proven not to fire).
- What 0.4.4 prints for the two Appendix-A lines with textbook discrepancies (F15).
- The outcome of every probe and reproduction (5.4–5.10), including whether `SplitPerturbations`
  accepts a background rule for `T̄` and whether the Euler–Lagrange step after the split works.
- The Wayback snapshot URL (created at execution).
- The userbase manifest before/after (a kernel run is needed to see the caches move).

## Appendix D — the two post-merge passes, and what each was asked to settle

This record is the archive of R-C's own planning. Two passes ran after it, both from plans the
user reviewed; their plans are archived here in summary so the trail is complete.

**Pass 1, the follow-up (2026-09-17, merged with the lane).** Asked for after the user
challenged the xMAG findings on the PSALTer #543 precedent: *be certain the failures were not
our setup*. Outcome: the xMAG verdict was overturned by a Tier-1 replay of the author's own
notebook cells (13 `identical`, 2 `proved-equal`); four of six reported negatives were our
calling forms; the "opposite sign" was xMAG's silent load-time `$RiemannSign = -1`; both
session hazards were separated and located; the project's curvature conventions were adopted
and recorded with kernel evidence (memo §7); four defects in our own harness were found and
fixed, with the third-party protocol written to prevent a repeat.

**Pass 2, the closing pass (2026-09-17, PR #588).** Asked for after the user required that the
deferred items stay visible in documentation and issue comments, and that the five remaining
recommendations be resolved. Five items, and what the runs did to each:

| item | outcome |
| --- | --- |
| the slice epsilon had no stated orientation | it was never ours to choose: PSALTer's own code carries the dictionary (`DefGeometry.m:60-63`). Adopted and tested in both signatures; the kernel proves neither sign, which is stated plainly |
| the family-A import map was asserted, not run | run, with the counting established block by block and a Lagrangian odd in the torsion as the control that must flip |
| the lost simplification was attributed to xMAG | measured to be **xBrauer's**, in one kernel in three states, with xMAG exonerated and two rival mechanisms excluded |
| two documentation gaps | §7.4 (what one `Needs` changes, per package) written; `RCSplitGuarded` moved into the harness and `f2` re-run from the committed bytes |
| one filed claim was broader than its evidence | #582's comment corrected to the exact 20 theory files |

**Two things the closing pass added that were not in its plan**, both because a run said so: a
new xPand defect (#589, the hard-coded slice determinant sign, which only shows up in the
project's own signature), and the finding that a multi-line definition without an outer bracket
had silently truncated the lane's own parametrization — amended at the instruction site in
`CLAUDE.md` and added to the third-party protocol.

**One thing it could not settle**, recorded rather than glossed: the import map with a
non-vanishing background torsion. Supplying a background rule for `Tor` alongside the
perturbation rule makes every order-1 piece vanish for each invariant tried, so the
half-applied map was never exercised. Memo §1.2 carries it with an owner.
*(Retracted 2026-09-18: the vanishing was never xPand's. The section that made those calls
passed four arguments to a three-argument helper; Mathematica left each call unevaluated, and
`ExtractOrder` of an unevaluated call is `0`. See Pass 3.)*

**Pass 3, the #591 follow-up (2026-09-18).** Asked for after #591 was routed to M3 as a tooling
gap that "usually fails and sometimes silently empties the first order": *find out what really
happens when xPand is given a nonzero background value, and how to resolve it*; move the
supervisor agenda to 18 September; build guards so that mistakes of the kinds above stop
recurring; and check every claim in the upstream reports before anything is sent. The user
fixed the scope at planning: **zero background torsion, exactly** — the new sector modifies
the perturbations only, on the unchanged ΛCDM background — to be raised at the supervisor
meeting. Outcome:

| item | outcome |
| --- | --- |
| what xPand does with a background value | a value that is a sum of terms along the normal (the rank-3 torsion shape) is filed as a projected background and used while the split rules are prepared, where it fails, and the messages point elsewhere; zero, and one-term vector and rank-2 backgrounds, work (F9) |
| the "silently empty first order" | ours — a four-argument call to a three-argument helper; retracted in #591, the memo and here |
| background torsion | out of scope by the user's decision (2026-09-17); #591 closed with the condition that reopens it: a theory of interest that does not admit zero background torsion on FRW |
| the xBrauer mechanism, withdrawn in this pass's own plan | confirmed one variable at a time (F10); the withdrawal — in the plan only, never committed or filed — rested on a comparison that changed two things |
| the reports already filed | xBrauer_Bundle #2 correct; xMAG #2 has three errors, with a correction comment drafted; every other report checked claim by claim, its snippets run as a reader would |
| guards | `wl_lint.py` before every kernel and in CI, `RCOrder`, `PACKAGE_FACTS.md`, protocol rules 9–10 and the snippet checker, each seen to fail first |

**Tension with round 4 of this record, left to the orchestrator:** the Terms section above
("Spectator", lines 113-121) accepts "a small nonvanishing background torsion that stays below the
tolerance". The Pass 3 scope is zero, exactly. The archived text is not edited; the decision
and its reopening condition are in the memo §1.2 and on #591.
