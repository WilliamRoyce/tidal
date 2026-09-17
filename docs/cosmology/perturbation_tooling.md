# FRW perturbation tooling memo — requirements, methods, tools, recommendation (R-C, #567)

<!-- cspell:words xPand xPert xTras xCoba xTensor xPerm xCore xMAG xIST COPPER HiGGS Hamilcar xPPN FieldsX bimEX xCPS xTerior TexAct xBrauer TraceFree Pitrou Umeh Brizuela Marugán Bahamonde Gigante Valcarcel Gorji Hohmann Heisenberg Kuhn Golovnev Koivisto Nikiforova Damour Chee Toporensky Tretyakov Odintsov Obukhov Bertschinger Sletmoen SymBoltz Cadabra Zumalacárregui Bellini Sawicki EFTCAMB Gubitosi Gleyzes Langlois Piazza Vernizzi Noller Helpin PSALTer Barker Cembranos Seljak Zaldarriaga Challinor Lasenby Lewis Nicosia Gakis Kiorpelidi Saridakis Mukhanov Feldman Brandenberger Peeters Seery Mulryne Ronayne Fröb WXF wolframscript ToxPand ToxPandFromRules SetSlicing DefMetricFields DefMatterFields DefProjectedTensor SplitPerturbations SplitMetric SplitMatter ExtractComponents ExtractOrder DefCovD DefTensorPerturbation DefMetricPerturbation DefConnectionPerturbation ToDistortion BreakDistortion BreakContorsion StartInducedDecomposition ChangeCurvature ChangeCovD ChristoffelCDCDT RicciScalarCDT TorsionCDT ContorsionCDT VarD CommuteCDSafe InducedFrom DefChart Detg normu tadpole eikonal contortion Weitzenböck teleparallel cspell adotoa dgrho hdot Kleg Kstd NOLAP ONSHELL OFFSHELL onshell offshell RCEL RCSort RCTry RCVerdict RCNormalize RCDigest ldd Einstein–Cartan Riemann–Cartan Shapiro Hehl SymManipulator Invar Spinors xPrint AVF SpinFrames EFTofPNG SymSpin TInvar SpaceSpinors xIdeal Harmonics BRST Ferreira Skordis Złośnik cona dcaabbea hersle nonproj recid regenerable syna tlsv xenos CDCDL KLEGACY RICCISIGN Torsioncd epsilonh kilj unsimplified RIEMANNSIGN TORSIONSIGN MAGChristoffelQ TorsionToDistortion FrozenMetricQ MasterOf EinsteinToRicci UndefChart UndefCovD PDBc Kstd RCSigns RCSetSigns RCNewMessages RCStackOn -->

> **Status: EVIDENCE COMPLETE, RECOMMENDATION FOR D-C — 2026-09-16, corrected 2026-09-17.**
> The follow-up lane of 2026-09-17 re-examined every negative xMAG finding against the
> author's own documented calls and stored outputs: four were this lane's calling forms or
> its connection declaration, one was an undocumented convention, and the package
> reproduces its documentation on this bundle (§1.1 T1, §3.3). The recommendation is
> unchanged (O1′) and now carries xMAG as an oracle test; §7 adds the curvature conventions
> the project had never stated, and §8 the amendments and the harness defects behind them.
> Original status line: **Status: EVIDENCE COMPLETE, RECOMMENDATION FOR D-C — 2026-09-16.** Written by research lane
> R-C for the orchestrator, who records D-C after merge; nothing here is a decision. Every
> "works"/"does not work" below is a kernel run judged by sentinel lines (§1.1) with the
> script committed under `scripts/research/perturbations/`; every citation is a TeX line
> pinned in `rc_planning_record.md` §2 (the verbatim planning archive, with the four review
> rounds). The four lane installs are additive and stamped; the userbase manifest outside
> them is byte-identical before and after (§6). Terms (spectator, tadpole gate, solver forms,
> post-Riemannian rewrite, `identical`/`proved-equal`) are defined in the planning record's
> Terms section and used here without repetition.

## 0. Summary

**What M3 must build** (Part 1): for a user theory of metric + torsion (+ vector) on a CAMB
ΛCDM background, in a CAMB-named gauge and in `(+,−,−,−)`, the coupled linear block of the
new sector and the standard modes it touches, the tadpole residual that gates admissibility,
the tensor and scalar source terms, and three export forms of the same equations (second
order, first-order state form, eikonal amplitude form for the photon channel).

**What the evidence says** (Parts 2–3, decision table in Part 4):

- **xPand 0.4.4 runs unmodified on the certified bundle** (no shim; Q3 never triggered) and
  is the right engine for the metric side: the 3+1 split, the SVT decomposition, both
  CAMB-named gauges, second order. The FRW tensor equation is reproduced from the
  Einstein–Hilbert **action** (headline A1, `proved-equal`, tadpole computed and vanishing on
  shell) and from the field equations (A2, `identical`); Ma & Bertschinger's eight scalar
  equations reproduce in both gauges; the `00` constraints match `camb.symbolic`'s own text in
  CAMB's signature (`proved-equal`, `c = 1`).
- **The project's signature works natively** with `n·n = +1`; two package sites hard-code
  `−1` and are replaced by a one-line projection and by the fluid split with the right norm
  (Probe C); the map to the mostly-plus literature has exactly three named sites.
- **Torsion is representable through xPand** as a rank-3 perturbed tensor with hand SVT rules
  (24 = 8 + 6·2 + 2·2, transcribed from the papers that did it): projections, the
  second-order action, the TT equations of both torsion tensor modes, and the derivative
  split all run; the connection itself must be rewritten first (the authors' stated exception,
  observed: xPand leaves `ChristoffelCDCDT` unsplit). Two session hazards were found and
  bounded (§3.2).
- **xMAG works as documented, in its own curvature convention.** The follow-up lane
  (2026-09-17) replayed the targeted cells of the author's own `DefCovD_xMAG.nb` against his
  stored outputs: 13 `identical`, 2 `proved-equal` (dummy names), 0 other. Its contortion is
  the standard one, its `ToDistortion` route and the by-hand xTensor identity agree under
  either curvature convention, and the vectorial-torsion textbook check passes in general
  dimension and at four. What R-C reported as an "opposite sign" was a comparison across two
  kernels with different values of the global `$RiemannSign`, which **xMAG sets to −1 at
  load, silently** (`xMAG.m:110`; xTensor's default is +1): under the project's convention
  the two routes are `proved-equal`, and the cross-convention comparison reproduces R-C's
  result exactly. Three of R-C's other four xMAG items were our calling forms (§3.3); the
  package-side findings are the undocumented load-time sign, the loss of xPand's
  transverse/traceless simplification in a shared kernel, and the collapse of xMAG's own
  machinery when a connection is declared without `Master` (§3.2). It is an **independent
  oracle for the rewrite now, and the component the moment non-metricity or connection
  variations enter** (Part 4).
- **Legacy's post-Riemannian rewrite is wrong** (`_derive.py:2164-2181`): its "contortion"
  is not a contortion (wrong slot order) and carries the wrong sign relative to xTensor's
  `ChristoffelCDCDT`; every `R̃` torsion term it produced differs from the correct identity
  (#582). Stage 1's port of that text (`stage1_engineering_plan.md:519-525`) must port the
  corrected identity instead.

**Recommendation (Part 4):** option **O1′** — xPand for the metric engine, the pure-xTensor
post-Riemannian rewrite with the kernel-verified identity, the lane's SVT rule set for the
rank-3 torsion, and the lane's η-space Euler–Lagrange operator; build-own only the exports
(solver forms, eikonal, line-of-sight sources), which no tool provides. xPand becomes an
install requirement (installer promotion, fingerprint, guide); the xMAG chain stays a research
install, used as the **oracle test** of our rewrite (an independent implementation of the same
identity, compared in a test, never called by the pipeline). Estimated M3 tooling cost on this
route: **≈ 6 engineer-days** (Part 4 table), unchanged by the follow-up.

## 1. What was verified here, and what was not

### 1.1 Run ledger (one Wolfram kernel at a time; Wolfram 14.3.0, xAct 1.3.0, xPert 1.0.6, xPand 0.4.4)

Every step: `run_lane.sh <step>` (guard, offscreen front end, throwaway cwd, hard timeout,
scrubbed transcript under `third_party/perturbations_runs/<step>/<utc>/`, gitignored).
Verdicts are the harness's (`RCVerdict`): `identical (SameQ)` / `proved-equal` (canonical
difference is 0) / `proved-different` / `could-not-decide`; `c` is the reported normalization.

| step (script) | what | kernel wall | result |
| --- | --- | --- | --- |
| gates before | `ensure_registered.sh` → `verify --require-psalter` | ~4 min | exit 0, "All checks passed", one informational xPerm-ldd warning; `gates/before_ensure_registered.txt` |
| manifest before | two-layer sha256 of the userbase | — | 3587 files; layer 1 (must be identical) 942 files, digest `9e2c81a4…6e63fbd` |
| install xPand | `install_xpand.sh` — since 2026-09-17 `scripts/install-xpand.sh` (sha256-pinned tarball, junk stripped, stamp) | 0 | 26 files added under `Applications/xAct/xPand/`, `INSTALLED_VERSION` (`patch=none`); Wayback snapshot refused by the container's TLS (not tested) |
| `probe_load` | load, usages, geometry, the paper's minimal example (`:1350-1373`) | 8 s | no load messages; `OUT47 identical (SameQ)`, control `proved-different`; label layout and gauge rules read back |
| A2 `repro_a2_tensor_eom` | tensor equation from `G^a_b + Λδ^a_b/κ − κT^a_b`, Λ + fluid | 13 s | mixed indices: `identical (SameQ)`, `c = 1`; controls (`2ℋ→3ℋ`, no `D²`) `proved-different`; unchanged on shell; lowered indices: `proved-equal` on shell (`c = a²`), off-shell residual `= (2Λa⁴/κ − 2a²ℋ² − 4a²ℋ′ − 2κa⁴P̄) E_ab` (the background equation × E: the #501 mechanism) |
| B `repro_b_mb_scalars` | Ma & Bertschinger `ein-cona..d` (`:644-661`) and `ein-syna..d` (`:603-618`) | 20 s | 8/8 `proved-equal`: Newtonian `c = 2, −2, 2, 1`; synchronous `c = 2, −2, −1, 1`; 8/8 controls `proved-different`; velocity sign flipped `proved-different`; Friedmann 1 both gauges (`c = −3`); traceless parts traceless |
| A1 `repro_a1_tensor_action` (**headline**) | `√−g[R/2κ − Λ/κ − ½(∂χ)² − V(χ)]` → split to order 2 → η-space Euler–Lagrange | 37 s (splits 0.6 + 1.2 + 5.6 s) | scalar control `χ″ + 2ℋχ′ − D²χ` `proved-equal` (`c = −a²`); tadpole: variation w.r.t. lapse = Friedmann 1 (`c = a²/κ`), w.r.t. curvature potential `= (a²/κ)(−3F₁ − 6F₂)` (`True`), w.r.t. χ = background χ equation (`c = −a²`); background rules solved from the tadpole equal the typed ones (`True`); all three vanish on shell (`True`); TT: `E″ + 2ℋE′ − D²E` `proved-equal` (`c = −a²/κ`), controls `proved-different`; off-shell residual = background equation × E |
| C `probe_c_signature` (`normu=1`) | xPand in `(+,−,−,−)`; CAMB's equations in CAMB's signature | 12 s | `SetSlicing[…, +1, …]` accepted, no messages, `n·n = 1`; minimal example `proved-equal` under the three-site map (without the Laplacian site: residual `4εD²φ + 8εD²ψ`; without the lapse flip: `proved-different`); `ExtractComponents` Time projection `−V0` (hard-coded `−1`), by hand `+V0`; built-in `ToxPand` fluid does not run, `SplitMatter[…, +1, …]` does; Friedmann by hand `c = −1`; Newtonian `00` vs map `c = −1`; **`camb.symbolic` Newtonian constraint [1] `proved-equal c = 1`, synchronous constraint [1] `proved-equal c = 1`**, controls `proved-different`; the CAMB traceless check was not decided (my projection kept the vector and tensor sectors) |
| D `probe_d_torsion` | rank-3 torsion beside `dg`; SVT rules; projections; action route; connection | 8 s (+ kinetic split 50 s in the first run) | D1 declared, no messages; D3 16 symbols = 24 functions, rule antisymmetric (`0`), TT and transverse properties hold; D4 all eight `Time/Space` projections evaluate; D6 order-2 split of `c₁T² + c₂TT + c₃T_aT^a + c₄(∇T)²`: 275 terms, tadpole `0` for `T̄ = 0`, no residual `CD`/`dTor`, TT sector: algebraic terms diagonal in the two TT modes, the kinetic term couples them through `4c₄a²ℋ ε^{ab}{}_c D^c TT2` (parity-even, `ℋ`-suppressed); Euler–Lagrange for each mode has `TT″`, `ℋ TT′`, `D²TT`; scalar sector carries all eight scalars; D2 xCoba count `48` (antisymmetric pairs not identified by `ToCanonical` on component indices; combinatorial `24`); D7/D8 moved to the clean-kernel run below after the session hazard |
| DL `probe_d_limits` (fresh kernel) | rank map; hand check; connection; the three rewrites; the hazard | 57 s | splits run for rank 0, 1, 2, 3, 4 of the torsion field with a pure-metric sanity split intact between them; `n^d n_a n^b h^c_e ∇_d T^a_bc = a (δT^0_{0i})′` `proved-equal` (`c = a`, no friction); `RicciScalarCDT` splits only to unsplit `ChristoffelCDCDT`/`Perturbation[ChristoffelCDCDT]` objects; xTensor: `CDT_b v^a − CD_b v^a = −ChristoffelCDCDT^a_bs v^s`; standard `K` is a contortion (`True`), legacy's is not (`False`); correct `R̃ = R + ¼T_abc T^abc + ½T_abc T^bac + T^a_a^b T^c_bc − 2∇_b T^a_a^b` (vector torsion: `R − 6v² + 6∇·v`); legacy's rewrite `proved-different`; xMAG's rewrite (imported WXF) `proved-different` from both signs of the standard substitution: torsion part negated; `δR̃ − δR = (18ℋ TS3 + 6TS3′ + 2D²TS1 + 6ℋD²TS4 + 2D²TS4′ + 4D²TS8)/a` (total derivatives: no tadpole at `T̄ = 0`); after `DefChart` the sanity split is `BROKEN` (`InducedFrom::unknown`) |
| install xMAG chain | `install_xmag.sh` (three pinned commits) | 0 | 38 files in five new directories, `INSTALLED_COMMIT` stamps, additivity asserted |
| E `probe_e_xmag`, E2 `probe_e2_xmag_induced` | xMAG beside xPand | 10–16 s, 8 s | loads (2 s, no messages), versions satisfied; xPand's `DefProjectedTensor` and split still work after loading and after xMAG's `DefCovD`; Riemann–Cartan connection defined with `ContorsionCDT`, `PerturbationChristoffelCDT`; `BreakContorsion` = standard `K` (`c = 1`); `ChristoffelCDCDT = −ContorsionCDT`; `ToDistortion` 0.03 s, free of `CDT` objects after `BreakContorsion`; `BreakDistortion` → `Null` on this connection; a second torsion connection → `ToDistortion` `Null`; `DefConnectionPerturbation` coexists with `dg` (`True`) but `Perturbation[TorsionCDT]` stays unexpanded; `ChangeCurvature` on a plain connection returns its input; `StartInducedDecomposition` fails on both connection types, with and without xPand's slicing (`DefCovD::invalid`, `TorsionQ::unknown`); xPand's split stays `ok` after loading xMAG, after xMAG's `DefCovD` and after `ToDistortion` + `BreakContorsion`, and is `BROKEN` (`InducedFrom::unknown`) from the moment a second connection is defined (`E_XPAND_SPLIT_AFTER_SECOND_CONNECTION_CDL=BROKEN`, run `e/20260916T114954Z`); the end-to-end split of xMAG's rewrite was therefore done in the clean kernel (DL, WXF hand-over) |
| **T1** `tier1_xmag` (follow-up) | xMAG loaded first, the author's own cells replayed against his stored outputs; then the sign experiment and the three calls R-C got wrong | 20 s | Tier-1: **13 `identical`, 2 `proved-equal`** (dummy names; canonical difference 0), 0 other; `SIGNS_AFTER_XMAG_LOAD={-1,1,1,1,1}`; A1/A2: xMAG's route and the by-hand identity `proved-equal` under **both** conventions; A3: cross-convention `proved-different` with the torsion part exactly negated (`proved-equal`); A4: vectorial torsion `R̃ − R = s_R(−(d−1)(d−2)v² + 2(d−1)∇·v)` `proved-equal` in general `d` and at `d = 4`, both signs; A5: the scalar commutator equals `−s_T T^c{}_{ab}∇_c f` (`proved-equal`), wrong-sign control `proved-different`; A6: `K_std` `identical` to `BreakContorsion`, both defining properties `0`; A7: recomputed under +1, the author's stored Riemann cell differs exactly by the negated torsion part; A8: torsion antisymmetric in its last two slots, contortion in its outer pair (family B); C1: two-argument `BreakDistortion` `identical` to `BreakContorsion`, three-argument form fails; C2: `StartInducedDecomposition` **runs** on the general connection with one-character symbols and silently sets `$ExtrinsicKSign = $AccelerationSign = −1`; C3: both `VarD` connection variations run |
| **F2** `f2_hazards` (follow-up) | the two "xPand breaks" cases, the twin control, the caller-side wrapper, the shared-kernel control | 21 s | chart case: split `BROKEN`, guarded split `ok` and `identical` to the reference, wrong-guard control `BROKEN`, `UndefChart` restores; xCoba count `24` with the sorted-pair dedupe (48 without); `Master`-less connection: split `BROKEN`, one-argument `EinsteinToRicci` `BROKEN` and the two-argument form `ok`, the `Master -> g` twin leaves everything intact, the strict guard recovers it, `UndefCovD` restores; **shared kernel: loading xMAG turns the transverse/traceless probe from `{0,0,0}` into three unsimplified expressions**, so the scalar sector still agrees `identical` and the residual is confined to vector/tensor terms |
| **F3** `f3_changecurvature` (follow-up) | `ChangeCurvature` on a `Master`-less torsion connection, plain vs xMAG | 8 s | plain xTensor expands it (`{ChristoffelCDCDL, g, RicciCD}`); with xMAG loaded it returns its input for any `Master`-less connection, and the one-argument `TorsionToDistortion` and `MAGChristoffelQ` are inert on it; with such a connection in the session even the healthy `Master -> g` twin fails, and `UndefCovD` restores it |
| **F5/F6** `f5_sign_audit`, `f6_psalter_signs` (follow-up) | the five sign globals per package | 8 s, 7 s | plain xTensor `{1,1,1,1,1}`; unchanged by xPand and by `SetSlicing` (`$ExtrinsicKSign = 1`); unchanged by PSALTer, which defines no curvature at all (`$Metrics = {G}`); xMAG flips `$RiemannSign` to −1 at load; `RCSetSigns` restores |
| gates after | `verify --require-psalter` | ~2 min | exit 0, "All checks passed"; `gates/after_verify.txt` |
| manifest after + diff | | — | layer 1 **identical** (942 files, `9e2c81a4…6e63fbd`); layer 2: 4 paclet-manager files changed (allowed list); new: 64 files in exactly the six directories (`xPand 26, SymmetricFunctions 7, BrauerAlgebra 9, xBrauer 9, TraceFree 5, xMAG 8`) |

Sentinel names quoted in this memo (`RC_…`) are grep-able in the transcripts; the WXF/txt
payloads of every reported expression sit beside them.

### 1.2 Not tested here, and why

| item | why not | consequence |
| --- | --- | --- |
| Wayback snapshot of the xPand tarball | archive.org refuses the container's TLS handshake (`tlsv1 alert access denied`) | the stamp records the failure; request from an unrestricted machine (orchestrator) |
| D5, a homogeneous background torsion mode | optional and outside the planned scope (planning record, Terms; round 4) | Part 4 gives no cost for O4-iso; the gate covers the mode without representing it (#579) |
| curved FRW (`"FLCurved"`), Bianchi | flat ΛCDM is the program's background | none for M3 |
| xMAG's induced decomposition **as an FRW 3+1 engine** (O4) | it runs when called as documented (T1 C2); exercising it as an alternative to xPand's slicing was outside the time box | O4 is no longer withdrawn, only untested as an engine, with a cost line in Part 4 |
| the CAMB traceless-equation transcription | projection kept vectors/tensors; MB's `ein-cond` covers the same equation in mostly plus | none |
| xPand + PSALTer in one kernel | the two branches run in separate kernels by design (D-A) | none |
| scalar and vector torsion Euler–Lagrange equations | the TT sector is the O2 channel; the scalar sector's symbols were checked present (D6) | M3 work, same operator |
| performance at scale | largest split here 50 s (order-2 kinetic torsion Lagrangian, 275 terms) | fine for M3; a full PGT Lagrangian is minutes, not hours |

## Part 1 — the equation set as a checklist

Consumers: **CAMB API** (`set_custom_scalar_sources`, scalar only, `model.py:1111`;
`get_time_evolution(frame=)`, `results.py:557`), **WS3 solver** (`solver_design.md:344-346,
532-533`), **WS4 line of sight** (`observable_ladder.md:141`), **WS2 gate** (`:139, :143`).
"Covered by" names the tool the lane ran for it (Part 3).

| # | requirement | source | consumer | covered by |
| --- | --- | --- | --- | --- |
| 1 | coupled linear block of the new sector on flat FRW, conformal time `η` = CAMB `t`, `a(η)`, `ℋ(η)` symbolic to the export, coefficients from CAMB's table at solve time; limits `a → const`, de Sitter | `R-C.md:83-84`; `repo_reshape.md:443-446, 827-829`; `conventions.md:75`; `COSMOLOGY_PROGRAM.md:1060-1062` | WS3 | xPand (`$ConformalTime`, `ah`, `Hh`; A1, A2, B) |
| 2 | gauge as an explicit CAMB-named input; projection in Wolfram with `camb.symbolic` as reference; gauge in spec metadata; seam assertion via `get_time_evolution(frame=)` | `repo_reshape.md:455-469, 830-831`; `conventions.md:65` | CAMB seam, WS3 | xPand `"NewtonGauge"`, `"SynchronousGauge"` (B); CAMB check (C) |
| 3 | rank-3 torsion and a torsionful connection; post-Riemannian rewrite about `(η, T̄ = 0)`; the field map from PSALTer's content to FRW variables in the solver-only block; the rewrite performed in both branches | `R-C.md:49-52`; `spectrum_design.md:341-351, 372-378`; `interfaces_decision.md:413` | WS2 | xPand + hand rules (D, DL); pure xTensor rewrite (DL); **not** xMAG (E) |
| 4 | the split of standard sectors: touched modes (metric, photon) inside the block; untouched sectors as tables and `S_std`; tensor target `h″ + 2ℋh′ + k²h = S_std + S_new` consumed by our solver (no tensor custom-source API in CAMB) | `R-C.md:43-48`; `repo_reshape.md:233-237`; `solver_design.md:330-343`; #576 | WS3 | xPand (A1/A2 give the operator; sources are ours) |
| 5 | the background-consistency gate (#501): the order-1 coefficient exported as a residual on the CAMB background, admissibility with a derivable tolerance, specified to cover a homogeneous new-sector mode | `spectator_route.md:81-91, 107-120`; `observable_ladder.md:185-193`; `COSMOLOGY_PROGRAM.md:200, 211`; #579 | WS2 gate | A1's tadpole extraction (`ExtractOrder[…,1]` + variation); A2's residual |
| 6 | exports: per-channel `S(k,η)`; the standard tensor source re-derived when the photon sector is modified (#514); the solver forms: second-order block, first-order state form `y′ = M(η)y (+s)`, eikonal amplitude form for the photon channel with dropped terms exported; several representations per derivation; batched-assembly-ready coefficients | `repo_reshape.md:329-348`; `birefringence_notes.md:181`; `solver_design.md:466-486, 532-533` | WS4, WS3 | none — build-own (Part 4) |
| 7 | conventions: native CAMB + PSALTer emission; `(+,−,−,−)`, `ε₀₁₂₃ = +1` asserted; CAMB variables verbatim at the seam; every MB sign difference attributed | `conventions.md` §1–§3, §6 | every reader | xPand with `normu = +1` (C); B's dictionary |
| 8 | the solver-only block carries: background fields, field → background + perturbation map, gauge choice, its own fingerprint | `interfaces_decision.md:408-416, 421-422` | WS2 | listed, not designed |
| 9 | auxiliary: validity flags (`ρ_new/ρ_γ`, `ΔN_eff`); energy/Hamiltonian time dependence (contested, #577) | `observable_ladder.md:143` | WS2 gate | ours |

## Part 2 — how the literature derives such equations

| method | input → output | connection | paper | code | what it leaves us |
| --- | --- | --- | --- | --- | --- |
| xPert/xPand line | covariant action or field equations → 3+1 split, SVT, gauge-fixed equations on FRW to any order | Levi-Civita of the metric | Brizuela–Martín-García–Mena Marugán arXiv:0807.0824; Pitrou–Roy–Umeh 1302.6174 | `xPert` 1.0.6 (xAct 1.3.0); `xPand` 0.4.4 (tarball, sha256 `26e7abca…`) | torsion after the rewrite; exports |
| second-order action for a torsion / connection perturbation on FRW, hand-parametrized SVT (8 + 6 + 2), gauge-invariant variables | action → quadratic action → equations | metric-affine or teleparallel, by hand | Aoki–Bahamonde–Gigante Valcarcel–Gorji 2310.16007 (`:768-878`); Heisenberg–Hohmann–Kuhn 2311.05495; Heisenberg–Hohmann 2311.05597; Nikiforova–Damour 1804.09215 (`:819-852`) | none acknowledged (2310.16007 `:1169-1170`); xAct "instrumental" (1804.09215 `:2296-2301`) | the SVT recipe (transcribed into Probe D3) |
| field equations with a nonzero background torsion, scalar sector | PGT field equations → scalar perturbations | Riemann–Cartan | Lu–Chee 1601.03943 (`:721-747`) | Maple (`:716-718`) | out of scope (background mode) |
| xPand as bookkeeping beside a tetrad matrix | tetrad perturbation → metric → xPand | teleparallel | Bahamonde et al. 2009.02168 (`:272-302`); Nicosia–Levi Said–Gakis 2012.11959 (`:331`) | xAct incl. xPand | how xPand has been used near torsion: never as the torsion engine |
| EFT of dark energy → α-functions hand-coded | action class → α(η) functions → Boltzmann code | metric only | Gubitosi–Piazza–Vernizzi 1210.0201; Gleyzes–Langlois–Piazza–Vernizzi 1304.4840; Bellini–Sawicki 1404.3713 | hi_class (1605.06102, 1909.01828; no license file); EFTCAMB (1312.5742, 1405.3590) | no derivation; no torsion → Horndeski map except special cases (Barker et al. 2006.03581) |
| symmetry-first quadratic action | field content + symmetries → most general quadratic action | metric | Lagos–Baker–Ferreira–Noller 1604.01396 | xIST 0.7.3 / COPPER 0.8.3 (GPL-3) | a user's theory-space tool, not a derivation component |
| equations-in Einstein–Boltzmann | hand-derived species equations → solver | — | Sletmoen 2509.24740 | SymBoltz.jl v1.7.0 (`hersle/SymBoltz.jl@3d1f20a3`, MIT), scalars only, Newtonian only | a possible consumer of scalar pieces; never a producer |
| hand-derived hierarchies + symbolic scalar source generation | — | — | Lewis–Challinor–Lasenby astro-ph/9911177; Ma–Bertschinger astro-ph/9506072 | CAMB 2.0.4 (`camb.symbolic`, `camb_fortran`, `compile_source_function_code`; LGPL-3 + GFDL) | the gauge/variable reference; scalar sources only |
| symbolic core → generated numerics (the D-A analogue) | model file → GiNaC/SymPy → C++ | flat FRW, scalars | Seery 1609.00380; Mulryne–Ronayne 1609.00381 | CppTransport 2018.1 (GPL-2+), PyTransport 2.0 (GPL-3+) | the pattern, not the physics |
| generic field-theory CAS | — | any | Peeters hep-th/0701238 | Cadabra2 2.5.14 (GPL-3) | no cosmological machinery |
| line-of-sight sources | Boltzmann hierarchy → `S(k,η)` → transfer functions | — | Seljak–Zaldarriaga astro-ph/9603033; Hu–White astro-ph/9702170 | CAMB's Fortran tensor source; `get_scalar_temperature_sources` | the method WS4/M3 follows for #514 |

## Part 3 — tools

Columns: documented scope · where used · license · last release · certified-bundle
compatibility (14.3.0 × xAct 1.3.0 × xPert 1.0.6) · Part-1 lines covered · torsion verdict
(**probed here** unless marked *documented only*).

| tool | scope | used where | license | release | compatibility | covers | torsion verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **xPand 0.4.4** (Pitrou–Roy–Umeh) | 3+1 split of perturbed tensors on FRW/Bianchi, SVT, gauges, conformal time, any order; "torsion theories the current exception" | 129 citing papers (INSPIRE 1221019): induced GWs, second-order CMB, Skordis–Złośnik | GPL-3 | 0.4.4, 2025-04-01 (`$Version`); history silent after 0.4.2 | **loads and runs unmodified**; version gate compares dates and cannot fire; xTensor 1.3.0 API fine | 1, 2, 3 (after rewrite), 4 (operator), 5, 7 | **representable with work**: the torsion as a rank-3 perturbed tensor + hand SVT rules; rank 0–4 splits; connection objects must be rewritten first (DL); session must stay free of metric-less `CovD`s (#583) |
| **xPert 1.0.6** | metric perturbation theory to any order | 264 cites (recid 790000) | GPL | 2018-02-28 | installed, in xPand's chain | 1, 3 | a second perturbed tensor (`DefTensorPerturbation` beside `dg`) works (D1, E) |
| xTensor 1.3.0 / xCore / xPerm 1.2.4 / xCoba 0.8.6 | abstract tensors, canonicalization, components | everywhere | GPL | 2025-12-29 | certified (fingerprinted four) | all | `DefCovD[…, Torsion -> True]` is the torsionful connection; `ChangeCurvature` + the standard contortion is the correct rewrite (DL) |
| xTras 1.4.2 | invariants, `VarD` helpers | xMAG dependency | GPL | 2014-10-30 | installed | 5 | — |
| **xMAG 0.1.0** (Helpin) + xBrauer 1.1.0 + TraceFree 0.1.0 + SymmetricFunctions 1.0.0 + BrauerAlgebra 1.1.0 | independent connection with torsion and non-metricity, distortion/contortion decomposition, connection perturbations and variations, induced decomposition | no paper; thesis 2407.18019 for the Brauer packages; two documentation notebooks with stored outputs | GPL-2+ (headers; no LICENSE files) | commits 2026-01-10 / 2013-11-12, no tags | **reproduces its own documentation on this bundle** (T1: 13 `identical`, 2 `proved-equal`); needs `Master -> g`; sets `$RiemannSign = −1` at load and two more signs inside `StartInducedDecomposition`; must not share xPand's kernel | 3 (rewrite) and the connection variation | **representable and correct**: contortion standard, `ToDistortion` equal to the by-hand xTensor identity under either convention; **adopted as the oracle test** of our rewrite, and as the component when non-metricity or connection variations enter (Part 4) |
| xIST 0.7.3 / COPPER 0.8.3 (Noller et al.) | most general quadratic action on FRW from symmetry | 1604.01396, 79 cites | GPL-3 | 2016 | untested (Mathematica 9/10 era) | — | *documented only*: end-user tool, not a component |
| xPPN (Hohmann) | post-Newtonian 3+1 about Minkowski, tetrad + Weitzenböck | 2012.14984, 11 cites | no license line | `xenos1984/xPPN@dcaabbea`, 2022-10-01 | untested | — | *documented only*: backgrounds hard-wired flat; not FRW |
| FieldsX 3 (Fröb) | fermions, gauge fields, BRST, frame fields | 2008.12422, 29 cites | GPL-2 | — | untested | — | *documented only*: theory-building aid |
| HiGGS → Hamilcar (Barker) | Hamiltonian analysis of PGT on Minkowski | 2206.00658; 2512.25007 | GPL-3 | Hamilcar 0.0.0 | untested | — | *documented only*: Minkowski |
| PSALTer 2.0.2 `bb45adb0` | particle spectrum on Minkowski | Stage 1 | GPL | installed | certified | — | the other branch; same rewrite needed there |
| Harmonics, Invar, Spinors, xPrint, SymManipulator, AVF, TexAct, xTerior, SpinFrames, EFTofPNG, bimEX, SymSpin, TInvar, SpaceSpinors, xCPS, xIdeal | (one line each) spherical harmonics; Riemann invariants; spinors; printing; symmetrized tensors; vector fields; TeX; exterior calculus; spin frames; post-Newtonian EFT; bimetric; spin-weighted; tensor invariants; space spinors; conformal; ideals | — | GPL | with xAct 1.3.0 | installed | — | *documented only*: irrelevant to FRW perturbations |
| `camb.symbolic` (CAMB 2.0.4) | sympy scalar ΛCDM equations, three gauges, scalar source code generation | the seam | LGPL-3 + GFDL | 2.0.4 | in `.venv` | 2, 7 | no torsion; **the reference text the lane matched** (C) |
| SymBoltz.jl v1.7.0 | equations-in solver, Newtonian, scalars | R-1 comparable | MIT | v1.7.0 | — | — | *documented only*: excludes our channels |
| hi_class / CLASS / EFTCAMB | hand-coded Boltzmann codes | — | none stated / cite-CLASS-II / CAMB's | — | — | — | *documented only* |
| CppTransport / PyTransport | model → generated numerics | inflation | GPL-2+ / GPL-3+ | 2018.1 / 2.0 | — | — | *documented only* |
| Cadabra2 2.5.14 | generic CAS | — | GPL-3 | 2025-07-31 | — | — | *documented only* |
| legacy `tidal/wolfram/ComponentDecompose.wl` | plain coordinate components, rank ≥ 3, `(−,+,+,+)` hard-wired, literal metric | legacy examples | MIT | v0.53 | — | — | no SVT, no Fourier, no `a(η)`; it does declare a torsionful connection (`uv run pytest tests/test_cli.py -k torsion_covd_dry_run` passes, 1 test, no kernel), but its `R̃` rewrite is wrong (#582), so it is a comparison only and not an oracle for torsion terms |

### 3.1 The xPand facts M3 will rely on (from the source and the runs)

- `SetSlicing[g, n, normu, h, cd, {post, pre}, type]`: `normu = ±1` accepted (only the
  determinant sign is checked, `xPand.m:1715`); `ExtractComponents` uses `−1` for an up index
  (`:3042`) and `ToxPand` passes `−1` as the fluid norm (`:3014`): in `(+,−,−,−)` project
  with `normu · n_a n^b` by hand and use `SplitMatter[u, du, +1, …]` + `ToxPandFromRules` (C).
- The split is a **list of rules** (`ToxPandFromRules[expr, rules, h, order]`); a rank-3
  field enters through a rule `dTor[LI[o_], a_, b_, c_] :> …` whose right-hand side is built
  from projected tensors of rank ≤ 2 (the only rank `DefProjectedTensor` gives properties to);
  pattern names must be manifold indices.
- `Sqrt[-Detg[]]` survives the conformal step as `Sqrt[-(Detg[] a^8)]` with `Detg[]` the flat
  conformal background's determinant: set `Detg[] -> -1` after the split (A1).
- `VarD` on the split density needs auxiliary tensors declared `OrthogonalTo -> {n[a], n[b]}`
  (xPand's induced derivative refuses non-projected arguments, `Validate::nonproj`), and the
  result needs xPand's own derivative commutation (`xAct`xPand`Private`CommuteCDSafe`) before
  transversality rules can act (A1's `RCSort`).
- The Lie-derivative labels `X[LI[order], LI[n]]` make the η-variation a Leibniz operator
  (`dEta`) plus `VarD` for the spatial part: `EL = P₀ − dEta[P₁] + dEta[dEta[P₂]]` (A1, D6).
- `ToCanonical::noident` on the perturbation parameter is benign; `Check` must ignore it.

### 3.2 The session hazards, located

1. **xPand's own** (#583): a covariant derivative with **no metric at all**
   (`MetricOfCovD[cd] === Null` — xCoba's chart derivative `PDBc`, or an affine connection
   declared without `FromMetric`) makes every later split throw `InducedFrom::unknown` and
   return `Null`. The line is `xPand.m:2566` in `ToMetric`, reached from every
   `ToxPandFromRules` through `Conformal` (`:2619`): it selects over `Rest@$CovDs` and asks
   `InducedFrom[MetricOfCovD[#]]`, where xTensor guards the same pattern in its own code
   (`xTensor.m:8951`). Measured workarounds (F2): the caller-side wrapper
   `Block[{$CovDs = Select[$CovDs, # === PD || MetricOfCovD[#] =!= Null &]}, …]` recovers a
   split `identical` to the clean-kernel reference, a wrong-guard control stays broken, and
   `UndefChart` also restores it. **For M3:** count components in another kernel, or use the
   five-line wrapper; the one-line upstream fix is drafted in `xpand_upstream_issue.md`.
2. **Ours** (reported inside #583/#585): a torsionful connection declared **without
   `Master`** breaks the split at a different site — xMAG's replacement of
   `EinsteinToRicci` (`xMAG.m:1154`) tests `FrozenMetricQ[MasterOf[#]]`, and
   `FrozenMetricQ[Null]` (`:1130`) reaches the same `InducedFrom[Null]` fallback
   (`xTensor.m:8302`) from `ToMetric` (`xPand.m:2569`). The twin declared `Master -> g`
   leaves the split, `ToDistortion` and `ChangeCurvature` intact (F2, F3); with the
   `Master`-less one present even the healthy twin fails, and `UndefCovD` restores it.
   xMAG's own `?Master` text makes `Master -> g` the documented declaration. **For M3:**
   every connection gets `Master -> g`.
3. **xMAG in xPand's kernel: no** (the decisive control, F2 §III). Loading xMAG turns the
   elementary transverse/traceless probe from `{0, 0, 0}` into three unsimplified
   expressions, so xPand's simplification of the vector and tensor sectors stops working;
   the scalar sector still comes out `identical` and the residual is confined to those
   terms. Separately, xMAG flips `$RiemannSign` at load, so a split in that kernel is in the
   other curvature convention. **For M3:** one package per kernel with WXF hand-over, which
   is what D-A's committed-package-plus-driver design already implies.

   > **Attribution corrected 2026-09-17 (post-merge exploration, relayed through the
   > orchestrator; #585):** the lost simplification is **xBrauer's**, not xMAG's — xMAG merely
   > pulls xBrauer in as a dependency. The operational rule is unchanged (one package per
   > kernel, WXF hand-over), but the cause names the right package: anything loading the
   > xBrauer chain inherits it, and xMAG alone is not what breaks xPand's transverse and
   > traceless simplifications.

### 3.3 xMAG, used correctly (what this lane had wrong, and what the author documents)

| R-C's call | what it did | the documented call | source |
| --- | --- | --- | --- |
| `BreakDistortion[expr, CDT, g]` | `Null` + `Validate::inhom`: the three-argument body expects non-metricity objects, which a metric-compatible connection never defines | `BreakDistortion[expr, CDT]` — the two-argument dispatcher routes a torsion-only connection to `BreakContorsion`; `identical` to calling it directly (T1 C1) | `xMAG.m:1476` vs `:1664` |
| `StartInducedDecomposition[g, CDT, {{";h", …}, …}, …]` | `DefCovD::invalid: ;h is not a valid Postfix symbol for a derivative`, then a cascade of consequences | one-character postfix symbols, on the **general** connection: `StartInducedDecomposition[g, CD, {{"~","Dh"},{"^","DGh"}}, {nv, hh}]` runs (T1 C2) | `xTensor.m:6591`; the usage says "works only for … GL(dim) independent connection"; the author's `In[175]` uses one-character symbols |
| a second `DefConnectionPerturbation`, then `Perturbation[TorsionCDT]` | the tensor exists, the torsion perturbation never does | it already ran inside `DefCovD`; vary with `VarD[ChristoffelCDT[-a,b,c], cd][L]` as the author does, and build the torsion perturbation as `δΓ − δΓᵀ` | `xMAG.m:425-426, 1079, 1367-1441`; Tutorial `In[300]` |
| a connection without `Master` | §3.2 item 2 | `DefCovD[CDT[-a], {"#","DT"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True]` | `?Master`; the author's `In[10]` |

Two more facts to carry: **load xMAG before defining anything, then re-assert the sign
globals** (`$RiemannSign = 1`); and **`StartInducedDecomposition` silently sets
`$ExtrinsicKSign = $AccelerationSign = −1`** (`xMAG.m:1792`), so print the signs again after
calling it.

## Part 4 — recommendation per requirement, decision table, D-C options, license

### 4.1 Per requirement

Cost is engineer time to production quality inside M3 (tests, docstrings, provenance), with
what the estimate rests on; "interface" is what the committed package calls under D-A.

| # | requirement | build-on / borrow / build-own | cost | rests on | D-A interface |
| --- | --- | --- | --- | --- | --- |
| 1 | FRW block, conformal time, symbolic `a`, `ℋ` | **build on xPand** | 0.5 d (wrap `SetSlicing`, `DefMetricFields`, `DefMatterFields` in the package) | probe_load, A1, A2, B ran in seconds | `ToxPandFromRules[L, rules, h, 2]` on the theory's Lagrangian density (WXF in), split expression out |
| 2 | CAMB-named gauge | **build on xPand** (`"NewtonGauge"`, `"SynchronousGauge"` = CAMB's CDM frame variables) | 0.5 d (name map + `camb.symbolic` assertion at derivation time) | B (8/8), C (CAMB `00` both gauges) | `SplitMetric[g, dg, h, gauge]` chosen by the user's gauge string |
| 3 | torsion field + connection, post-Riemannian rewrite | **build on xTensor** for the rewrite (`ChangeCurvature` + `ChristoffelCDCDT → −K_std`, with the two kernel checks as tests); **borrow** the SVT rule set (Probe D3, from 2310.16007); build-own the rule generator for general rank | 1.5 d (rewrite 0.5, rules 1) | DL (`KSTD_IS_A_CONTORTION`, vector-torsion check), D3–D6 | rule list `Join[SplitMetric[…], torsionRules, {T̄ -> 0}]`; rewrite applied before the split |
| 4 | coupled block with the touched standard modes | **build on xPand** for the operator; sources ours | 1 d | A1/A2 (operator); `S_std` from CAMB tables | the split's order-2 density → Euler–Lagrange per field |
| 5 | #501 gate | **borrow A1's machinery** (order-1 coefficient + variation w.r.t. lapse, curvature potential, each new field) | 0.5 d | A1 tadpole (`ONSHELL_ALL_ZERO=True`), A2 residual | `ExtractOrder[L2, 1]` → `RCEL` → residual expression exported with `a`, `ℋ` symbolic |
| 6 | exports: solver forms, eikonal, LOS sources | **build-own** (no tool provides them; the split output is the input) | 1.5 d for the two solver forms + coefficient batching; eikonal (#504) and LOS (#514) are their own tickets | `solver_design.md` §5, §7 | package function from the Euler–Lagrange system to `y′ = M(η)y (+s)` |
| 7 | conventions native | **build on xPand with `normu = +1`** + two by-hand sites | 0.5 d (projection helper, fluid split with norm, assertions in the artifact header) | C (map at three sites; CAMB `00` identical in form) | `SetSlicing[g, n, +1, …]`; header asserts `n·n = +1`, `ε₀₁₂₃ = +1` |
| 8 | solver-only block contents | listed (Part 1) — design is M3's | — | `interfaces_decision.md` §3.8 | — |
| 9 | validity flags | build-own | 0.5 d | — | — |
| | **total tooling on O1′** | | **≈ 6 d** | | |

### 4.2 Decision table — what each probe returned, what it moves

| probe | outcome (this lane) | favors | cost line moved |
| --- | --- | --- | --- |
| load (5.4) | loads clean, no shim | O1–O3 open | "xPand maintenance": 0 |
| A1 headline (5.6) | order-2 split + η-variation `proved-equal`; three machinery facts needed (§3.1), all cheap | **O1 as full engine** | "E-L after the split": ≈ 1 d, not "unknown" |
| A2 control (5.5) | mixed `identical`; lowered residual = background × E | validates the #501 mechanism | "gate": 0.5 d |
| B + CAMB (5.7) | 8/8 MB; CAMB `00` in both gauges `c = 1` in CAMB's signature | conventions rule confirmed: run natively, no table | "CAMB seam transcription": 0 beyond the assertion |
| C signature (5.8) | split native; two hand-done sites | O1 native | "signature transcription": 0.5 d |
| D torsion, xPand route (5.9) | representable with the precedent work (SVT rules); rank 0–4 split; connection objects need the rewrite first | O1 or O3 | "torsion SVT rules": 1 d |
| D action route (D6) | tadpole 0 at `T̄ = 0`; TT equations for both modes; scalar sector present | action route lives for torsion | "action route for torsion": 0 extra (same operator) |
| DL rewrite | correct identity established in pure xTensor; legacy's wrong (#582) | **O1′** (pure xTensor rewrite) | "connection → LC + torsion": 0.5 d + the Stage-1 correction |
| E xMAG (5.10) | loads; contortion right; four of its six negative items turned out to be our calling forms | superseded by T1 | — |
| **T1 xMAG called as documented** | reproduces the author's stored outputs; agrees with the by-hand identity under either convention; the sign difference was the undocumented load-time global | **O1′ with xMAG as the oracle**; O3 becomes a later adoption rather than a blocked one | "xMAG convention audit": 0, done; "oracle test": 0.25 d |
| **F2/F3 hazards** | one xPand robustness defect with a wrapper that recovers the split identically; one configuration error of ours; xMAG must not share xPand's kernel | O1′ unchanged; the wrapper is five lines in our package | "xPand `$CovDs` wrapper": 0.25 d |

### 4.3 D-C options, in the order the evidence ranks them

1. **O1′ — xPand + pure-xTensor rewrite + lane rule sets (recommended).** xPand for 3+1,
   SVT, gauges, second order, `normu = +1`; the post-Riemannian rewrite by
   `ChangeCurvature[L, CDT, CD]` with `ChristoffelCDCDT → −½(T^a_bc + T_b^a_c + T_c^a_b)`,
   guarded by the two kernel tests; Probe D3's SVT rules; A1's Euler–Lagrange operator.
   Cost ≈ 6 d (§4.1). Consequence (i): the rewrite is shared by both branches as one
   committed function — Stage 1 must **not** port legacy's text (`stage1_engineering_plan.md:519-525`).
   Consequence (ii): xPand becomes an install requirement: promote `install_xpand.sh` to
   `scripts/install-xpand.sh`, add `xPand.m`'s `$Version` to `verify-wolfram-setup.sh`'s
   fingerprint (today it checks xCore, xPerm, xTensor, xCoba only; #580), a `WOLFRAM_GUIDE.md`
   step, re-certification. Routed to the orchestrator (paths not owned here). **Done
   2026-09-17** by the orchestrator: the installer is `scripts/install-xpand.sh`, `verify`
   gained its own xPand check (DEGRADED when absent or different, watched failing), the guide is
   a seven-step path, and the configuration was re-certified.
2. **O3 — xMAG for the connection algebra, xPand for the split: adopt when needed, no
   longer blocked.** The follow-up removed the reason it was deferred — there is no sign
   defect to settle. What it does not buy is anything O1′ lacks for metric-compatible
   torsion, so the recommendation stays O1′ **with xMAG as the oracle test** (0.25 d), and
   O3 is adopted the moment **non-metricity** (`∇g ≠ 0`) or a **Palatini-style variation
   with respect to the independent connection** enters: 1 d plus promoting the five-package
   install, one package per kernel, WXF hand-over.
3. **O2 — xPert only, own 3+1/SVT.** Not needed: xPand's engine passed every probe. Cost
   would be 10–15 d (a rewrite of what xPand does), kept as the fallback if upstream xPand
   ever breaks against a future xTensor (the version gate cannot warn).
4. **O4 — xMAG's induced decomposition as the 3+1 engine.** No longer withdrawn: it runs
   when called as documented (T1 C2). Untested as an FRW engine and not recommended —
   xPand's slicing passed every reproduction — but evaluating it would cost a day, not be
   impossible.
5. **O5 — port legacy's decomposer.** Rejected on design grounds and now on correctness
   (#582).
6. **O6 — fork xPand for torsion.** Not needed: no shim, no patch; torsion enters through
   rules xPand accepts by design. A fork would drift (the `slegner/CAMB` precedent, #498)
   and put the derivative under GPL-3.

**License as a cost line.** The repository is MIT. xAct and PSALTer are GPL and already
loaded at derivation time; calling xPand (GPL-3) the same way is the existing footing and
costs nothing. O1′ vendors and ports nothing from GPL code: the SVT rules and the
Euler–Lagrange operator are the lane's own. O3 would add four GPL-2+ packages as
requirements; O6 would put a fork under GPL-3.

**The D-A interface, concretely.** One committed package function per Part-1 line, called by
the fixed `driver.wls` on WXF data: theory density in → split density out; split density +
gauge string in → Euler–Lagrange system out; the same system in → `y′ = M(η)y (+s)` out with
`a`, `ℋ` symbolic; order-1 coefficient in → gate residual out. Every function returns
expressions in xPand's projected-tensor heads, which the exporter maps to CAMB names at the
seam (the dictionary of B and C).

## 5. Findings routed to the orchestrator (issues filed; nothing edited here)

| id | finding | issue |
| --- | --- | --- |
| F3, F4 | `primer.md` mostly-plus/φ; `spectrum_design.md` still cites MB conventions and the withdrawn opposite-signature premise | #575 |
| F5 | "standard sectors only as sources" vs "coupled standard mode inside the block" | #576 |
| F6 | the `t = 0.0` energy bug is critical / moot / must-fix in three places | #577 |
| F7 | `gauge_certificate` dropped vs live | #578 |
| F8 | O4-iso needs a homogeneous mode; the admissible-theories row and the #501 gate spec must say how it is covered | #579 |
| F1, F2, F9–F12 | `R-C.md` line citations, FLRW/FRW, SymBoltz URL, the four fingerprinted packages | #580 |
| F14 | the H7 near-miss table says xPand handles "any theory" | #581 |
| F13 | `docs/references.md` had no xPand row | added in this lane (Q2) |
| F15 | the xPand paper's Appendix A prints `−2D²φ` and `E′ℋ` where 0.4.4 prints `D²ψ` and `2ℋE′` | #584 |
| L1 | legacy's contortion identity is not a contortion and has the wrong sign; Stage 1 ports it | **#582** (priority high) |
| L2 | xPand: a metric-less `CovD` in the session breaks every later split | #583; draft `xpand_upstream_issue.md` |
| L3 | xMAG on the certified bundle. **Corrected by the follow-up:** four of the six reported items were our calling forms or our connection declaration; the package-side findings are the undocumented load-time `$RiemannSign = −1`, the loss of xPand's transverse/traceless simplification in a shared kernel, the `Master`-less collapse, and the three-argument `BreakDistortion` precondition | #585 (rewritten); draft `xmag_upstream_issue.md` |
| L5 | the project had **no stated curvature convention**; adopted here as xTensor's defaults, with the import map for family-A sources (§7) | #586 (extended); `conventions.md` §2 at merge |
| L6 | four defects in our own harness made earlier probes misreport; the protocol that prevents a repeat is in `scripts/research/perturbations/README.md` (§8) | fixed in this lane; no issue |
| L4 | xPand under `(+,−,−,−)`: the three-site map and the two by-hand sites — the rule `conventions.md` §6 should carry | #586 |

## 6. Evidence index

- **Installs.** xPand: `Applications/xAct/xPand/INSTALLED_VERSION` — `0.4.4`, tarball
  sha256 `26e7abcac7bb655235ec39b73850729cf4465748d0d5ea2ad03c0607aef5ceab`,
  `installed_utc=2026-09-16T10:37:41Z`, `patch=none`, Wayback line records the TLS refusal.
  xMAG chain: `INSTALLED_COMMIT` in each of `Applications/SymmetricFunctions`,
  `Applications/BrauerAlgebra` (bundle `48be67e1c9037650f661c9d8ae538e807e407ba1`),
  `Applications/xAct/xBrauer` (same), `Applications/xAct/TraceFree`
  (`4e53ab3996f3be5da312576f1a24fb7be2ab5ddf`), `Applications/xAct/xMAG`
  (`88026e47af1f8651f999775339e5c080d82d9d04`). Transcripts:
  `third_party/perturbations_runs/install/{xpand_install.txt, xmag_install.txt}`.
- **Gates.** `gates/before_ensure_registered.txt` (exit 0, end `2026-09-16T10:36:21Z`),
  `gates/after_verify.txt` (exit 0). Both show "All checks passed!" and the one
  informational warning R-1 saw.
- **Manifest.** `manifest/{before,after}.{all,layer1,layer2,new}.sha256`, `manifest/diff.txt` and, after the follow-up lane, `manifest/diff_followup.txt`:
  `MANIFEST_LAYER1=identical` (942 files, `9e2c81a485c218355a53f81f424553bf24bbd77b30c0cf18f545f9ecf6e63fbd`
  before and after), `MANIFEST_LAYER2_CHANGED_FILES=4` (`Paclets/Configuration/*` and
  `Paclets/…` manager data, on the allowed list), `MANIFEST_NEW_AFTER_FILES=64` in the six
  listed directories only — identically before and after the follow-up lane, which added
  no file to the userbase.
- **Digests** (SHA256 of the canonical form after `ScreenDollarIndices`; fingerprints, not
  verdicts): Out[47] `e446cd3c…`; A2 mixed `b5ee18a0…`, lowered on shell `13924a59…`; MB
  Newtonian a–d `aeca68da…`, `0083bf07…`, `62431812…`, `4bfe054d…`; synchronous a–d
  `61e494e7…`, `126e8cc0…`, `81545f20…`, `757fc13e…`; A1 scalar control `8b4f0529…`, A1 TT
  `8a4a07ab…`. Full values in the run directories' `*.txt` beside the WXF payloads.
- **Targets.** Ma & Bertschinger `literature/astro-ph_9506072/9506072.tex` `:573-580,
  :603-618, :644-661`; MFB eq. 4.15 for the tensor operator; CAMB's own text in
  `scripts/research/perturbations/camb_symbolic_newtonian_2.0.4.txt` (regenerable by
  `mb_camb_symbolic.py`; camb 2.0.4); the SVT parametrization `literature/2310.16007` `:768-878`.
- **Scripts** (all under `scripts/research/perturbations/`, README has the provenance table
  and the reproduction block): `install_xpand.sh`, `install_xmag.sh`, `manifest.sh`,
  `run_lane.sh`, `wolfram/RCSetup.wl`, `wolfram/probe_load.wls`, `wolfram/repro_a2_tensor_eom.wls`,
  `wolfram/repro_b_mb_scalars.wls`, `wolfram/repro_a1_tensor_action.wls`,
  `wolfram/probe_c_signature.wls`, `wolfram/probe_d_torsion.wls`, `wolfram/probe_d_limits.wls`,
  `wolfram/probe_e_xmag.wls`, `wolfram/probe_e2_xmag_induced.wls`, and from the follow-up
  lane `wolfram/RCSetupCore.wl`, `wolfram/tier1_xmag.wls`, `wolfram/f2_hazards.wls`,
  `wolfram/f3_changecurvature.wls`, `wolfram/f5_sign_audit.wls`,
  `wolfram/f6_psalter_signs.wls`; plus `mb_camb_symbolic.py`,
  `xpand_upstream_issue.md`, `xmag_upstream_issue.md`.
- **Literature fetched into the worktree's `literature/`** (31 ids; the orchestrator copies
  them into the main checkout and regenerates `literature/README.md`, Q2): 1302.6174,
  0807.0824, 2310.16007, 2311.05495, 2311.05597, 1601.03943, 1804.09215, 2009.02168,
  2012.11959, 2206.00658, 2012.14984, 1604.01396, 2008.12422, 2407.18019, 2212.14496,
  1210.0201, 1304.4840, 1404.3713, 1605.06102, 1909.01828, 1312.5742, 1405.3590, 1609.00380,
  1609.00381, astro-ph/9911177, astro-ph/9702170, hep-th/0701238, 2011.02491, 1808.05565,
  2203.01856, 2110.12332. Curated rows: `docs/references.md` § "FRW perturbation tooling".

## 7. Curvature conventions: what each source uses, what the project adopts

Nothing in the project stated a curvature convention before this lane: `conventions.md` fixes
the signature and the ε orientation only, PSALTer defines no curvature at all and sets none of
xTensor's sign globals (F6), and CAMB carries the Einstein equations in a fixed form with no
Riemann tensor. xPand follows xTensor's defaults (its paper `:312, :337`); xMAG changes one of
them at load. **Adopted, and to be quoted into `conventions.md` §2 by the orchestrator:**

> Curvature signs are xTensor's defaults, `$RiemannSign = $RicciSign = $TorsionSign =
> $epsilonSign = +1` (`xTensor.m:287-289`, defaults at `:1837-1843`), with the derivative
> index in the middle slot of the connection and of the contortion. They are asserted at the
> start of every kernel, re-asserted after every `Needs` and after any induced decomposition,
> and printed into every artifact header. xPand needs no adjustment; **xMAG sets
> `$RiemannSign = −1` when it loads** (`xMAG.m:110`, undocumented) and
> `$ExtrinsicKSign = $AccelerationSign = −1` inside `StartInducedDecomposition`
> (`xMAG.m:1792`), so a kernel that uses xMAG resets them.

### 7.1 The statements, as the kernel printed them

Each line is grep-able in the named run directory under `third_party/perturbations_runs/`.

| statement | sentinel and value | run |
| --- | --- | --- |
| xTensor's own definitions | `RC_F5_RIEMANNSIGN_USAGE`: "Riemann[-a,-b,-c,d] = $RiemannSign * ( PD[-b][Christoffel[d,-a,-c] + …)"; `RC_F5_RICCISIGN_USAGE`: "Ricci[-a,-b] = $RicciSign * Riemann[-a,-c,-b,c]"; `RC_F5_TORSIONSIGN_USAGE`: "cd[-a]@cd[-b]@f[] - cd[-b]@cd[-a]@f[] = - $TorsionSign Torsioncd[c,-a,-b] cd[-c]@f[]" | `f5/` |
| the five globals, per package | `RC_SIGNS_PLAIN_XTENSOR={1, 1, 1, 1, 1}`; `RC_SIGNS_AFTER_XPAND={1, 1, 1, 1, 1}`; `RC_SIGNS_AFTER_XPAND_SETSLICING={1, 1, 1, 1, 1}`; `RC_SIGNS_AFTER_PSALTER={1, 1, 1, 1, 1}`; `RC_SIGNS_AFTER_XMAG_LOAD={-1, 1, 1, 1, 1}`; after `StartInducedDecomposition`, `{1, 1, 1, -1, -1}` | `f5/`, `f6/`, `t1/` |
| the connection difference | `RC_DL_CDT_MINUS_CD_ON_A_VECTOR=-(ChristoffelCDCDT[ia, -ib, is]*vv[-is])` with `RC_DL_CHRISTOFFELCDCDT_IS_MINUS_K=True`; with xMAG loaded, `RC_T1_CHRISTOFFELCDCDT_AUTORULE=-ContorsionCDT[i1, -i2, -i3]` | `dl/20260916T114120Z`, `t1/` |
| the contortion identity | `RC_T1_A6_KSTD_IDENTITY_INPUTFORM=HoldForm[Kstd[a, b, c]] -> (TorsionCDT[a, b, c] + TorsionCDT[b, a, c] + TorsionCDT[c, a, b])/2`; `RC_T1_A6_KSTD_EQUALS_XMAG_CONTORSION=identical (SameQ)`; `RC_T1_A6_KSTD_ANTISYM_MINUS_T=0`; `RC_T1_A6_KSTD_METRICITY=0`; legacy's expression `RC_DL_KLEGACY_IS_A_CONTORTION=False` | `t1/`, `dl/20260916T114120Z` |
| the Einstein–Cartan check | `RC_T1_A4_VECTOR_TORSION_VS_TEXTBOOK_PLUS1=proved-equal` and `RC_T1_A4_AT_DIM4_PLUS1=proved-equal` for `R̃ − R = s_R(−(d−1)(d−2) v·v + 2(d−1)∇·v)`, i.e. `R − 6v² + 6∇·v` at `d = 4`; the same under `−1` with the sign carried | `t1/` |
| the torsion sign | `RC_T1_A5_SCALAR_COMMUTATOR=-(TorsionCDT[-is, -i1, -i2]*CDT[is][ff[]])` with `…_VS_MINUS_T_GRAD_F=proved-equal` and the control `…_PLUS_T_GRAD_F=proved-different` | `t1/` |
| the slot family | `RC_T1_A8_TORSION_IS_ANTISYM_IN_LAST_TWO=identical (SameQ)`; `RC_T1_A8_CONTORSION_ANTISYM_PAIR={identical (SameQ), proved-different}` — the outer pair, not the first two | `t1/` |
| the `(+,−,−,−)` map and the two by-hand sites | `RC_C_OUT47_VS_SIGNATURE_MAP=proved-equal`; without the Laplacian site `proved-different` with residual `4εD²φ + 8εD²ψ`; the control without the lapse flip `proved-different`; `RC_C_TIME_PROJECTION_OF_V0_n=-V0` against `RC_C_TIME_PROJECTION_BY_HAND_OK=True`; `RC_C_SPLITMATTER_WITH_NORMU_RAN=True`; `RC_C_CAMB_NEWTON_00_VS_XPAND=proved-equal up to c=1`; `RC_C_CAMB_SYNC_00_VS_XPAND=proved-equal up to c=1` | `c/20260916T114552Z` |

### 7.2 The dictionary, and the import map for family-A sources

Two families of contortion appear in the literature, differing in which slot of the connection
carries the derivative index. Both are correct; no software we use is in family A.

| source | signature | Riemann / Ricci vs xTensor | torsion / contortion | displays `R̃ = R + …`? |
| --- | --- | --- | --- | --- |
| xPand paper 1302.6174 | (−,+,+,+) `:369, :680` | same slot order and sign `:312`; Ricci "second and fourth" `:337` | not stated | no |
| Ma & Bertschinger | (−,+,+,+) `:295, :361` | not stated | not stated | no |
| Aoki et al. 2310.16007 | (−,+,+,+) `:153`; ε₀₁₂₃ = +1 `:201` | other slot order `:181`; `R̃_{μν} = R̃^λ{}_{μλν}` `:185` | `T = 2Γ̃_{[μν]}` `:160`; **family A** `:168` | no |
| Nikiforova–Damour 1804.09215 | (−,+,+,+) `:221`; ε^{0123} = +1 `:245` | frame form `:286`; `R_{ij} = η^{kl}R_{kilj}` `:296` | **family A** `:344`, `K_{ijk} = −K_{jik}` `:341` | in words `:256-260` |
| Heisenberg–Hohmann–Kuhn 2311.05495 | not stated (mostly plus) | other slot order `:23`; `R_{μν} = R^λ{}_{μλν}` `:36` | derivative slot stated via `Q` `:14`; **family B** `:18` | **yes** `:38, :46` |
| Shapiro hep-th/0103093 | not stated | other slot order `:659`; `R̃_{τβ} = R̃^α{}_{ταβ}` `:687` | **family A** `:623`, `K_{αβγ} = −K_{βαγ}` `:633` | **yes** `:695, :736` |
| Barker 2206.00658 | (+,−,−,−) `:265` | other slot order `:130` | derivative slot last `:219`; **family B** `:226` | no |
| PSALTer code and paper | (+,−,−,−) `PSALTer.m:101`, `2406.09500:102` | flat, no curvature; sets no sign global (F6) | the user's field strengths | no |
| xTensor 1.3.0, hence xPand and xMAG | — | the reference (the blockquote above) | **family B**: `T^a{}_{bc} = s_T(Γ^a{}_{bc} − Γ^a{}_{cb})`, `K` antisymmetric in its outer pair | through `ChangeCurvature` |

**Import map, applied once at import with the source cited and a test** (as
`conventions.md:27-29` already requires): a torsion tensor read from a family-A source is
**minus** xTensor's `TorsionCDT`, its Riemann slot order is remapped, and its contortion is
re-expressed as `K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c + T_c{}^a{}_b)`. This applies to
Probe D3's SVT parametrization, transcribed from Aoki et al.: the map flips the sign of its 24
potentials, which leaves representability and the counting untouched but must travel into M3's
rule set.

### 7.3 Every convention this lane touched, and where it now stands

| convention | settled | where |
| --- | --- | --- |
| signature `(+,−,−,−)`, and how it enters xPand (`normu = +1`, determinant sign −1) | yes | `conventions.md` §2; Probe C; §7.1 |
| ε orientation `ε₀₁₂₃ = +1` (`$epsilonSign`) | yes | `conventions.md` §2; F5 |
| xPand's slice tensor `epsilonh` relative to `ε` and to `n` | **not checked** | §1.2; needed before Probe D3's parity-odd rules are used in anger |
| Riemann sign, Ricci contraction, Ricci scalar sign | **adopted here** | the blockquote above; F5, F6 |
| torsion sign, the connection's derivative slot, the contortion family, the family-A import map | **adopted here** | §7.1, §7.2; T1 A5, A6, A8 |
| extrinsic curvature and acceleration signs in the 3+1 split | recorded: xPand `+1`, xMAG's induced decomposition sets both to `−1` | F5, T1 C2 |
| perturbation definition, xPand's projected variables, the MB dictionary | yes, by run | Probe B, `probe_load` |
| conformal time, `a(η)`, `ℋ`, CAMB names and gauges | yes | `conventions.md` §3; Probe C |
| Fourier convention (`k² ↔ −∇²`, `θ`, the velocity sign) | yes, by run | Probe B's velocity-sign control; Probe C |
| units `κ = 8πG`, `c = 1` | yes | Part 1 |

## 8. Amendments to the plan, and the harness defects behind them

- **Targets were typed inline** in each comparing script with the TeX line in a comment; the
  planned separate `targets.wl` was never written, and `mb_camb_symbolic.py` now points at
  the two scripts that cite its printout. The planned `compare.wls` was not written either:
  `RCVerdict` in `RCSetupCore.wl` implements the verdict rule, and the dangling `compare`
  step has been removed from `run_lane.sh`.
- **A1's pre-split `VarD` cross-check was dropped**: the on-shell action route already
  reproduces the field-equation route, which is what the cross-check was for.
- **Four defects in this lane's own harness** made the first xMAG probes misreport, and are
  fixed in `RCSetupCore.wl`: `RCTry` used `Check`, so any message — including the benign
  `DefMetric::old` that the author's own run prints — became a failure; `$MessageList = {}`
  is a no-op on a protected symbol, so every message list printed was cumulative rather than
  per-call (`RCNewMessages` takes the tail instead); the harness file was read before any
  xAct package, so its short names created shadowing `Global` symbols and every verdict in
  one run was computed by inert functions; and a comment containing `*)` closed early,
  truncating the file silently. Two more were in the probes: `CommuteCovDs` was called with
  its index pair reversed (`xTensor.m:6218` matches outer, then inner), and one comparison
  ran against an unexpanded symbol. **The protocol that prevents a repeat is in
  `scripts/research/perturbations/README.md`, "Working with a third-party xAct package".**
- **The xCoba component count** of the rank-3 torsion is 24, not the 48 first reported: the
  dedupe must key on the sorted antisymmetric pair (F2).
