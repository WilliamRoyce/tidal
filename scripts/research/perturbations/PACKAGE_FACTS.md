# Verified facts about the packages this lane uses

<!-- cspell:words xPand xMAG xBrauer xPert xTensor xPerm xCoba xTras TraceFree PSALTer ToxPandFromRules SetSlicing DefProjectedTensor DefTensorPerturbation DefMetricPerturbation ExtractOrder ExtractComponents SplitMatter ProjectionAndBackgroundRuleQ RulesCovDsOfTensor Projectorh InducedFrom MetricOfCovD CommuteCDSafe CommuteCovDs Detg normu epsilonh OrthogonalTo ProjectedWith MetricOfTensor FirstMetricQ MasterOf ContractMetric ContractMetric1 SeparateMetric IndicesDown DefCovD BreakDistortion BreakContorsion ToDistortion ChangeCurvature StartInducedDecomposition EinsteinToRicci Contorsion Distorsi SetConnectionRelations ConnectionRelations StringLength RCOrder RCTry RCNewMessages RCSigns RCSetSigns noident inhom partw Anticommutator RightCosetRepresentative ChristoffelCTCDk ChristoffelCDkCT CDk UpValue TvPW Invar tbar UNCOMPUTED TOCANONICAL EPSH CDCDL nonproj -->

Read the section for a package **before writing code against it**. Every entry is something a
kernel showed, with the run or source line that showed it; nothing here is inferred from a
name or from how another package behaves. When a fact turns out wrong, correct it here and say
so in the memo's §8 — do not delete it silently.

Two rules sit above every section:

- **Print any value the package chooses; never assume it.** The perturbation parameter's name,
  the sign globals, the context a bare name resolves to, the metric a tensor reports — each of
  these has been different from what we assumed at least once.
- **A result proves it was computed before anything is read off it.** Mathematica returns a
  call it cannot evaluate unchanged, and a later step turns that into a plausible answer (see
  the first section). Use `RCOrder`, not `ExtractOrder`, on a split.

Runs are `third_party/perturbations_runs/<step>/<utc>/transcript.txt` (gitignored; re-create
them with `run_lane.sh <step>`). Line numbers are the installed versions: xTensor 1.3.0,
xPand 0.4.4, xMAG `88026e47`, xBrauer `48be67e1`.

## The Wolfram Language itself — where it fails without saying so

| fact | shown by |
| --- | --- |
| A call that matches no definition — for example, a helper called with four arguments when it takes three — is returned **unevaluated**, with no message. `ExtractOrder` of such an unevaluated call is `0`: an "empty first order" from a split that never ran | `f9` part E: `RC_F9_E_FOUR_ARGUMENT_CALL_UNCOMPUTED=True`, `RC_F9_E_EXTRACTORDER_OF_UNCOMPUTED_IS_ZERO=True`. The origin of the retracted #591 claim. `wl_lint.py` (`ARITY`) refuses it |
| A newline ends a statement as soon as the expression is complete. `f[x_] := a + b` followed by a line `+ c` defines `a + b`; the next line is a separate statement that does nothing | `f7`'s first run (8 potentials instead of 16); `CLAUDE.md`. `wl_lint.py` (`MULTILINE`) |
| `*)` inside a comment closes it; the rest of the file is then parsed as code or silently truncated | memo §8. `wl_lint.py` (`BALANCE`) |
| A line is **parsed before it is evaluated**. In ``Needs["xAct`xTensor`"]; $RiemannSign`` the name is read before xTensor loads and becomes an unrelated ``Global`$RiemannSign`` | `upstream/check_snippets.py` on the snippet filed in xMAG #2: ``got 'Global`$RiemannSign'``. `wl_lint.py` (`NEEDS`) |
| The same applies to a harness file read **before** a package: its short names become `Global` symbols that shadow the package's, and every function using them is inert | memo §8; `RCSetupCore.wl` must be read after the package a step is about |
| `$MessageList = {}` does nothing (the symbol is protected); message lists printed after it are cumulative | memo §8. Use `RCNewMessages[tag, expr]` |
| `Check` treats every message as failure, including benign ones such as `DefMetric::old`, which the xMAG author's own successful run prints | memo §8. `RCTry` does not use `Check` |

## Names that change meaning when a package loads

| fact | shown by |
| --- | --- |
| Every xAct package exports its own `$Version`. After xPand loads, a bare `$Version` is **xPand's**; after xMAG, **xMAG's**. Write ``System`$Version`` | `f11/20260918T133113Z`: `RC_F11_BARE_VERSION_CONTEXT_AFTER_XPAND=xAct`xPand``, `…_AFTER_XMAG=xAct`xMAG``. The harness printed xPand's version as `RC_WOLFRAM_VERSION` in every run before 2026-09-18. `wl_lint.py` (`SHADOW`) |
| Deliberate xAct replacements of System names: `Symmetrize`, `Commutator`, `Anticommutator`, `Up`, `Down` (xTensor); `Cycles`, `RightCosetRepresentative` (xPerm) | `f11`: `RC_F11_SHARED_WITH_SYSTEM_AFTER_XPAND` |
| xPand exports `a`, `H`, `P`, `T`, `k`, `Es`, `Et`, `Ev`, `Bs`, `Bv`, `Ls`, `Lv`, `Vs`, `Vv`, `V0`, `Xi`, `CS`, `Connection`, `nt`, `av` — never reuse them | `RCSetup.wl` header |

## xTensor 1.3.0

| fact | shown by |
| --- | --- |
| The sign globals default to `+1`. `RCSigns` prints five, in this order: `$RiemannSign`, `$RicciSign`, `$TorsionSign`, `$ExtrinsicKSign`, `$AccelerationSign` (`$epsilonSign` is also `+1`) | `xTensor.m:1837-1843`; `f5`: `RC_SIGNS_PLAIN_XTENSOR={1, 1, 1, 1, 1}` |
| A derivative's postfix symbol must be exactly one character | `xTensor.m:6591`; `t1` C2 |
| `CommuteCovDs` takes the index pair outer-then-inner | `xTensor.m:6218`; memo §8 |
| A tensor declared with `ProjectedWith -> {h[…]}` gets its own automatic projector rule, which contracts `h` **when the product is built** — before `ContractMetric` is ever called. A comparison that adds `ProjectedWith` changes two things at once | `xTensor.m:4258`; `f10`: `RC_F10_CONTRACTED_WHEN_BUILT={TvPW}` |
| `ToCanonical` alone contracts the induced metric into a plain projected tensor, with or without xBrauer | `f8`: `RC_F8_BASE_TOCANONICAL_ALONE_CONTRACTS_H=True` |
| A product of two slice epsilons is the number `6` at construction; a substitution applied afterwards has nothing to act on | `f7`: `RC_F7A_EPSH_SQUARE_RAW_IS_SIX=True` |
| `ChangeCurvature` on a connection without `Master` expands normally in plain xTensor | `f3`: `{ChristoffelCDCDL, g, RicciCD}` |

## xPand 0.4.4

| fact | shown by |
| --- | --- |
| The perturbation parameter is `\[Epsilon]` (`$PerturbationParameter`) **whatever symbol is passed to `DefMetricPerturbation`**. A "did the split run" test must look for `$PerturbationParameter`, never for the name you chose | `f9` part A: `RC_F9_A_PERTURBATION_PARAMETER=\[Epsilon]`; the #591 diagnosis repeated this mistake once |
| Load xPand first; nothing session-wide changes at load (all five signs stay `+1`, also after `SetSlicing`) | `f5`; every run's `RC_SIGNS_AFTER_XPAND` |
| `SetSlicing` accepts `normu = ±1`; `ExtractComponents` multiplies an upper index by `−1` for "Time" (`:3042`) and `ToxPand` passes `−1` as the fluid norm (`:3014`), so with `n·n = +1` project by hand and use `SplitMatter[u, du, +1, …]` | probe `c`; `upstream/xpand_item4_time_projection.wls` |
| The induced metric is declared with determinant sign `+1` whatever `normu` is (`:1754`), so with `n·n = +1` `epsilonh` squares to `+6` against `−6` for `ε_abcd n^d` | `f7` at `normu=+1`; #589; `upstream/xpand_item2_slice_sign.wls` |
| A covariant derivative with **no metric** anywhere in the session (an xCoba chart's, or a connection without `FromMetric`) makes every later split return `Null` with `InducedFrom::unknown` (`ToMetric`, `:2566`) | `f2`; #583; `upstream/xpand_item1_chart.wls` |
| **Background values.** A rule for a field's background is filed as a "projected background" when the induced projector annihilates its value (`ProjectionAndBackgroundRuleQ`, `:2881`) and is then used early, while the split rules are prepared (`RulesCovDsOfTensor`, `:2867`). A value that is a **sum** of terms along the normal — the shape a rank-3 torsion background takes — leaves an unexpanded sum inside `Projectorh`; `ToCanonical` fails, the error is caught at `:2840`, a `Null` enters the rule list, and the split fails later with `ReplaceAll::reps` and `Validate::inhom` | `f9` parts B, C: `RC_F9_B_PREBUILD_TORSION_SUM_HAS_NULL=True`, `RC_F9_C_TORSION_SUM_BACKGROUND_RAN=False`; #591; `upstream/xpand_item3_background_rule.wls` |
| A zero background works for every shape; so does a one-term background along the normal for a vector or a rank-2 field (the vector's order 0 is `−c5 tbar² a⁴`, as by hand). A two-term torsion background is not filed as a projected background and completes without messages, but its output has **not** been checked against anything | `f9` parts C, D |
| A failed split leaves the session unharmed: the zero-background split afterwards is `identical` to before | `f9` part F: `RC_F9_F_UNCHANGED_AFTER_FAILURE=identical (SameQ)` |
| `ToCanonical::noident` on the perturbation parameter is benign; `RCTry` quiets it, so a script that needs to see it must not wrap the call in `RCTry` | `RCSetupCore.wl` (`RCTry` docstring) |
| After the conformal step `Sqrt[-Detg[]]` survives as `Sqrt[-(Detg[] a^8)]`: set `Detg[] -> -1` after the split | `a1` |
| `VarD` on a split density needs auxiliary tensors declared `OrthogonalTo -> {n[a], n[b]}` (`Validate::nonproj` otherwise), then xPand's own `CommuteCDSafe` | `a1` (`RCSort`) |
| The four-index projection rule at `:1809` is in the source but never fires in our sessions | `f7`: `RC_F7A_XPAND_1809_RULE_FIRES=False` |

## xPert 1.0.6

| fact | shown by |
| --- | --- |
| A second perturbed tensor (`DefTensorPerturbation`) beside the metric perturbation works | probes `d`, `e` |

## xMAG (`88026e47`)

| fact | shown by |
| --- | --- |
| **Sets `$RiemannSign = −1` at load**, silently (`:110`). Re-assert the signs with `RCSetSigns` after the `Needs` | `f5`, `t1`, `f11`: `RC_SIGNS_F11_END={-1, 1, 1, 1, 1}` |
| `StartInducedDecomposition` sets `$ExtrinsicKSign = $AccelerationSign = −1` (`:1792`) | `t1` C2 |
| Load xMAG before defining anything, as its author does | `t1` (the author's cells replay 13 `identical`, 2 `proved-equal`) |
| Declare every connection `Master -> g`. Without it `ChangeCurvature` returns its input, `EinsteinToRicci` (replaced by xMAG, `:1154`) fails, and even a healthy connection in the same session breaks | `f2`, `f3`; xMAG's `?Master` |
| For a metric-compatible torsionful connection use `BreakDistortion[expr, CDT]` (two arguments); the three-argument form expects non-metricity and returns `Null` | `t1` C1; `xMAG.m:1476` vs `:1664` |
| The contortion-trace relation is never installed: stored under `ConnectionRelations[covd, Contorsion]` (`:1038`), looked up under `Distorsi` (`:1044`) and `Distortion` (`:556`) | `upstream/xmag_issue2_followup.md`, checked by `check_snippets.py` |
| Conversion rules are attached to Christoffel names ordered by `StringLength`, so for some name pairs (`CDk` with `CT`) the rule sits on `ChristoffelCTCDk` while `ChangeCurvature` produces `ChristoffelCDkCT` and nothing converts | a kernel check on 2026-09-18, described in `upstream/xmag_issue2_followup.md`; not a lane step |
| Pulls in xBrauer, TraceFree, xTras and Invar; see their sections | memo §7.4; `f8` for xBrauer |

## xBrauer (`48be67e1`, arrives with xMAG)

| fact | shown by |
| --- | --- |
| **`ContractMetric` stops contracting the induced metric `h` into tensors that have no `Master`**, which is every field xPand's `DefProjectedTensor` makes; `g` still contracts. The re-installed rules (`:1947-1974`) require `MetricOfTensor[T] === h`, or `MetricOfTensor[T] === Null` with `h` the first metric — and `h` is second in `$Metrics = {g, h, gah2}` | `f8`; `f10`: guard value predicts every outcome (`RC_F10_GUARD_PREDICTS_CONTRACTMETRIC_AFTER_XBRAUER=True`, with `RC_F10_CONTROL_GUARD_IS_NOT_CONSTANT=True`) |
| So xPand's transverse and traceless simplifications stop working in a kernel that has loaded xBrauer (or xMAG): **one package per kernel** | `f2`, `f8` |
| `SeparateMetric`, and with it xPand's `IndicesDown`, stops separating a derivative index or an `epsilonh` index (`:1866`, `:1874`) | `f8` |
| `MetricOfTensor[T[inds]]` keys on the first slot (`:316`), so a tensor whose first slot is a label (`LI`) gives `Null` | [xBrauer_Bundle #2](https://github.com/THelpin/xBrauer_Bundle/issues/2), snippet re-run as filed |
| `:1958` has its guard inside `ScalarsOfProduct[…]`, unlike its three siblings | source; xBrauer_Bundle #2 |

## xTras, Invar, TraceFree

| fact | shown by |
| --- | --- |
| Invar sets `$CommuteCovDsOnScalars = False` at load (`Invar.m:204`, announced) — it silently blocked a scalar-commutator check until set back | `t1` A5 |
| TraceFree appends a `TraceFree` option to `DefTensor` (`TraceFree.m:202-204`) | source |

## PSALTer 2.0.2

| fact | shown by |
| --- | --- |
| Sets none of the sign globals and defines no curvature (`$Metrics = {G}`) | `f6` |
| Declares its own slice epsilon and the dictionary to the four-index one (`DefGeometry.m:41-46, 60-63`); this is the adopted convention | `f7`; memo §7 |
