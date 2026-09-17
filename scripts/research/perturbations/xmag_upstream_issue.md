# Draft report for the xMAG author — NOT FILED

<!-- cspell:words xMAG xBrauer TraceFree xPand xTras xPert Helpin ToDistortion BreakDistortion BreakContorsion ToContorsion DefConnectionPerturbation StartInducedDecomposition contortion DefCovD ChangeCurvature Riemann-Cartan Accelerationn CDCDL Tetrahh epsilonhh changeRiemann TorsionToDistortion MAGChristoffelQ EinsteinToRicci FrozenMetricQ MasterOf ConnectionRelations chriscdmetcovd Distorsi StringLength xTensor ExtrinsicKSign AccelerationSign RiemannSign Tensorx chriscovdmetcovd -->

Drafted by research lane R-C (#567) on 2026-09-16 and **rewritten on 2026-09-17** after a
Tier-1 replay showed that most of the first draft's items were the lane's own calling forms.
The user decided at planning that anything found in xMAG is drafted here and filed or sent to
the author only with their agreement (planning record §8, Q4). Everything below is from
`THelpin/xMAG@88026e47` (0.1.0, 2023-05-18) with `THelpin/xBrauer_Bundle@48be67e1` and
`xAct-contrib/TraceFree@4e53ab39`, on Wolfram 14.3.0 and xAct 1.3.0 (xTensor 1.3.0,
xTras 1.4.2, xPert 1.0.6), unmodified. Scripts:
`scripts/research/perturbations/wolfram/{tier1_xmag,f2_hazards,f3_changecurvature}.wls`;
transcripts under `third_party/perturbations_runs/{t1,f2,f3}/`.

## What was verified working

The targeted cells of `Documentation/DefCovD_xMAG.nb` were replayed in the author's order,
with his names, and compared with the outputs stored in the notebook: **13 `identical`, 2
`proved-equal`** (dummy-index names only, canonical difference 0), nothing else. That covers
`DefCovD` for all three geometries, `ChristoffelCDT → BreakChristoffel → ToContorsion →
BreakContorsion`, `ToContorsion` of `RiemannCDT`, `ToDistortion` of `ChristoffelCD`, the
two-argument `BreakDistortion`, and `ToDistortion` of `RiemannCD`. In addition:

- `BreakContorsion[ContorsionCDT[a,-b,-c], CDT]` gives `½(T^a_bc + T_b{}^a{}_c + T_c{}^a{}_b)`,
  which is a contortion (its antisymmetric part is the torsion, and it is metric-compatible —
  both checked to `0`).
- `ToDistortion` + `BreakContorsion` of `RicciScalarCDT[]` equals the same rewrite done by
  hand with xTensor's `ChangeCurvature` and that contortion, **under either value of
  `$RiemannSign`** (`proved-equal` in both cases).
- The vectorial-torsion check `R̃ − R = s_R(−(d−1)(d−2)v·v + 2(d−1)∇·v)` passes in general
  dimension and at `d = 4`.
- `StartInducedDecomposition` runs on a general connection with one-character postfix
  derivative symbols, and `VarD[ChristoffelCDT[-a,b,c], cd][L]` varies with respect to the
  connection as the tutorial shows.

So the report below is about documentation and robustness, not about wrong results.

## 1. `$RiemannSign` is set globally at load, and nothing says so

`xMAG.m:110-111` sets `$RiemannSign = -1` and `$RicciSign = 1` at load. `$RiemannSign` is
xTensor's session-wide global, whose default is `+1` (`xTensor.m:1837`), and it is read at
call time by `changeRiemann` (xTensor's at `:6194`, xMAG's replacement at `xMAG.m:1258`),
`CommuteCovDs` and the curvature relations. Nothing in the README, the usage strings or the
two documentation notebooks mentions it; a `grep` for `RiemannSign` across them returns
nothing.

The consequence for a user: every result xMAG produces, and every result **any other xAct
package produces in the same session**, is in the opposite curvature convention to the one
that session had before `Needs["xAct`xMAG`"]`. Comparing a result computed with xMAG loaded
against one computed without it silently differs by the sign of the torsion part of `R̃` —
which is exactly the mistake this lane made and then attributed to xMAG.

Suggested: one line in `ToDistortion::usage` and in the README stating the convention, or a
`Print` at load like the other variable changes xAct packages announce
("** Variable $CovDFormat changed from Prefix to Postfix"), or — best — localizing the
choice so that it is not a session-wide side effect. `StartInducedDecomposition` likewise
sets `$ExtrinsicKSign = $AccelerationSign = -1` (`xMAG.m:1792`) without announcing it.

## 2. A connection declared without `Master` silently disables the package

`?Master` says `Master -> met` is what "triggers the definitions and relations appropriate"
for each geometry, and every notebook example passes it. When it is omitted, however, the
failure is silent and contagious:

- `ChangeCurvature[RicciScalarCDL[], CDL, CD]` returns its input unchanged and prints
  nothing, where plain xTensor (same connection, xMAG not loaded) expands it into
  `{ChristoffelCDCDL, g, RicciCD}`.
- `TorsionToDistortion[TorsionCDL[a,-b,-c]]` and `MAGChristoffelQ[CDL, g]` stay unevaluated:
  the one-argument dispatcher `xMAG.m:1705-1707` folds over `$CovDs` and calls
  `TorsionToDistortion[expr, CDL, MasterOf[CDL]]`, and no definition matches a `Null` master
  (`:1615` requires `metric_?MetricQ`).
- Worse, that inert call is produced **inside** xMAG's `changeRiemann` replacement, so a
  `Master`-less connection anywhere in the session breaks `ToDistortion` on a perfectly
  healthy `Master -> g` connection as well (`Validate::inhom`); `UndefCovD` on the offender
  restores it.

Suggested: a definition of `TorsionToDistortion[expr, covd_, Null]` (and its siblings) that
returns `expr`, or a message from `DefCovD` when a torsionful connection is defined without
`Master`.

## 3. The three-argument `BreakDistortion` has an unstated precondition

For a metric-compatible connection, `BreakDistortion[expr, CDT, g]` matches the body at
`xMAG.m:1476`, which builds a rule out of `Distortion` and `NonMetricity` tensors that such a
connection never defines; the result is `Null` with `Validate::inhom`. The two-argument form
dispatches correctly (`:1664` → `BreakContorsion`) and gives the right answer. Suggested: a
`NonMetricityQ` guard on the three-argument definition, or a note in the usage.

## 4. Loading xMAG disables other packages' automatic simplifications — but the cause is xBrauer, not xMAG

**Moved to its own draft, `xbrauer_upstream_issue.md`, and re-attributed.** The symptom stands:
in a kernel with xPand's FRW slicing the elementary checks `h^{ab}∇_b E^{(V)}_a`,
`h^{ab}∇_b E^{(T)}_{ac}` and `h^{ab}E^{(T)}_{ab}` all canonicalize to `0`, and after
`Needs["xAct`xMAG`"]`, with nothing else changed, none of them simplifies. The first draft
suspected the deletion of xTensor's generic definitions in `deflistablexTensorxMAGCovDs`
(`xMAG.m:1134-1177`) and said so without isolating it. A three-state measurement
(`f8_contractmetric.wls`, run `f8/20260917T114838Z`) has now isolated it, and **xMAG is not
the cause**:

- `Needs["xAct`xBrauer`"]` **alone** reproduces the loss in full
  (`F8_VERDICT_REGRESSION_REPRODUCED_BY_XBRAUER_ALONE=True`). xBrauer arrives through xMAG's
  own `BeginPackage`, which is the only reason the symptom appeared when xMAG was loaded.
- loading xMAG on top changes **nothing** about metric contraction: `SubValues[ContractMetric1]`
  has the identical count and hash in the two states, and every probe agrees expression by
  expression (`F8_XMAG_ADDS_NOTHING_TO_CONTRACTMETRIC1=True`,
  `F8_XMAG_PROBES_MATCH_XBRAUER={True, True, True, True}`).
- the first draft's suspect does fire, but on a different function: `EinsteinToRicci` gains two
  `DownValues` and a new hash only after xMAG loads, and its `DownValues` are byte-identical
  before and after xBrauer. That replacement is not what breaks the simplifications.

The mechanism, the reproduction and the suggested one-line widening of the guard are in the
xBrauer draft. Both drafts are addressed to the same author and can be sent together.

## 5. Small things found by reading

- `xMAG.m:1044` passes `ConnectionRelations[covd, Distorsi]` (a typo for `Distortion`), so
  the contortion-trace relation `K^a{}_{ab} = -T_b` at `:1038` is never installed for a
  metric-compatible torsionful connection.
- `xMAG.m:1349` and `:1356` clear `chriscdmetcovd`, which is undefined; the intended symbol
  is `chriscovdmetcovd`.
- `xMAG.m:1283, :1297, :1310` choose which of `Christoffel[cd1,cd2]`/`Christoffel[cd2,cd1]`
  gets the automatic conversion by comparing `StringLength` of the two derivative names,
  while xTensor stores the pair sorted by symbol order; for some naming choices (a
  three-character metric derivative and a two-character connection) the two disagree and the
  conversion silently never fires.
- `xMAG.m:1764`: the generic `ToContorsion` has an `If` with no else branch, so it can
  return `Null` for a torsionful non-metric connection reached through that definition.
