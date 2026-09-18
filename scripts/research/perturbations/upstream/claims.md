# Claim-by-claim checks for the upstream reports

<!-- cspell:words xPand xMAG xBrauer xPert xTensor xPerm xCoba MetricOfTensor FirstMetricQ MasterOf ContractMetric ContractMetric1 SeparateMetric SeparateMetric1 SeparateMetric2 IndicesDown InducedMetricQ DefProjectedTensor ProjectedWith ProjectionAndBackgroundRuleQ RulesCovDsOfTensor Projectorh ExtractComponents SplitMatter ToxPand ToMetric InducedFrom MetricOfCovD SetSlicing Distorsi Contorsion SetConnectionRelations ConnectionRelations chriscdmetcovd chriscovdmetcovd ToContorsion StringLength ChristoffelCTCDk ChristoffelCDkCT BreakDistortion BreakContorsion ChangeCurvature CDCDL -->

Every sentence in a report that asserts something, how it was checked, and the result. A report
is ready only when every row reads **verified**. How to re-run the snippet checks:
`python3 check_snippets.py` in this directory (one kernel at a time; it refuses if one is live).

Where a row cites a lane step (`f7`, `f8`, `f9`, `f10`, `a2`, `b`), the script is
`../wolfram/<script>.wls`, run with `../run_lane.sh <step>`.

## xMAG #2 — as filed on 2026-09-18, and the follow-up comment

| claim in the filed issue | how checked | result |
| --- | --- | --- |
| Replaying the author's `DefCovD_xMAG.nb` cells reproduces his stored outputs "exactly" | lane step `t1` (R-C follow-up) | 13 identical, 2 equal up to dummy-index names — **nuance added in the follow-up** |
| The rewrite matches a hand rewrite under either sign convention | `t1`, arms A1–A8 | verified |
| Environment line (xTras 1.4.2, the three pinned commits) | each package's qualified `$Version` in one kernel (2026-09-18); `INSTALLED_COMMIT` files; upstream HEAD compared via the GitHub API | verified; our pins are upstream HEAD |
| Item 1: xMAG sets `$RiemannSign = -1` at load; xTensor's default is +1 | source `xMAG.m:110`, `xTensor.m:1837`; kernel | verified |
| Item 1: the snippet ``Needs["xAct`xTensor`"]; $RiemannSign`` prints `1` | run as filed | **wrong** — prints ``Global`$RiemannSign``; corrected in the follow-up |
| `StartInducedDecomposition` sets `$ExtrinsicKSign = $AccelerationSign = -1` at `:1792` | source | verified |
| Item 2: `ChangeCurvature[RicciScalarCDL[], CDL, CD]` returns its input unchanged, no message | kernel, fresh session | verified |
| Item 2: plain xTensor expands the same call | kernel without xMAG | verified: `{ChristoffelCDCDL, g, RicciCD}` |
| Item 2: afterwards `ToDistortion` fails with `Validate::inhom` on a healthy `Master -> g` connection | kernel: works before, fails after | verified |
| Item 2: the snippet is self-contained | read | **no** — assumes a manifold and metric; setup given in the follow-up |
| Item 3: three-argument `BreakDistortion` → `Null` with `Validate::inhom`; two-argument form works | kernel | verified |
| Item 4: `Distorsi` at `:1044` is a typo "for `Distortion`" | source `:1036-1044`, `:556`; kernel | **wrong key** — the relation is stored under `Contorsion`; `SetConnectionRelations` (`:556`) uses the empty `Distortion` key too; corrected in the follow-up |
| Item 4: so the contortion-trace relation is never installed | kernel: trace unreduced after `DefCovD` and after `SetConnectionRelations`; reduces once installed under `Contorsion` | verified |
| Item 4: `:1349`, `:1356` clear the undefined `chriscdmetcovd` | source | verified |
| Item 4: `StringLength` ordering means a conversion never fires for some names | kernel: `CDk`/`CT` — conversion on `ChristoffelCTCDk`, `ChangeCurvature` yields `ChristoffelCDkCT`; control `CD`/`CDT` converts | verified (was an inference; example added in the follow-up) |
| Item 4: `:1764` has an `If` without an else branch | source | verified |
| "Related": the contraction problem is reported on xBrauer_Bundle | — | **premature** until `xbrauer_issue_contraction.md` is filed; the follow-up links it |
| Every snippet in the follow-up comment | `check_snippets.py` | all as expected |

## xBrauer_Bundle #2 — as filed on 2026-09-18

| claim | how checked | result |
| --- | --- | --- |
| The snippet prints `g` for `T[-a]` and `Null` for `dT[LI[1], -a]` | run exactly as filed, line by line | verified |
| The indexed form at `:316` uses `{inds}[[1]]` | source | verified |
| The usage says the function returns the metric used to raise and lower indices | source `:214` | verified |
| No downstream failure traced in the minimal case | kernel: contraction with `g` into `dT` still works | verified as stated |
| `:1958` has the guard inside `ScalarsOfProduct[…]`; the three siblings do not | source; kernel reads the four sites: `{"]", "&", "]", "]"}` | verified |

## xBrauer contraction issue — `xbrauer_issue_contraction.md`, not yet filed

| claim | how checked | result |
| --- | --- | --- |
| Every snippet output (seven commented lines) | `check_snippets.py` | all as expected; the one descriptive line (separated with `q`) read by eye |
| xPand's splits keep terms that should vanish, with no message | `f8` (xPand alone vs xBrauer alone), `f10` | verified |
| `:1938-1940` issue three `Unset`s; `:1947-1974` install twelve guarded rules | source; `f8` counts `SubValues` 20 → 28, twelve mentioning `FirstMetricQ` | verified (the number of rules actually removed was never measured, and the text does not claim it) |
| Without `Master`, `MetricOfTensor[tensor]` is `Null` (`:315`), so the condition reduces to `FirstMetricQ[q]`, which is `False` | source; `f10`: guard `False` for the no-`Master` tensors, `True` with `Master` | verified |
| The guard alone decides the outcome for every tensor that reaches `ContractMetric` | `f10`: `F10_GUARD_PREDICTS_CONTRACTMETRIC_AFTER_XBRAUER=True`, with the control that the guard is not constant | verified |
| `ProjectedWith` hides the problem by contracting when the product is built | `f10`: `TvPW` contracted when built | verified |
| `:1866` removes the catch-all `SeparateMetric1`; `:1874` replaces `SeparateMetric2`; its derivative branch separates only when `MasterOf[covd] === metric` | source | verified |
| xPand's `IndicesDown` changes with it | `f8` | verified |
| Not related to #2 (head form vs indexed form) | source; `f10` head-form values | verified |
| The suggested fix | — | **not tested**, and the text says so |

## xPand email — `xpand_email.txt` with four attachments, not yet sent

| claim | how checked | result |
| --- | --- | --- |
| Ma and Bertschinger's scalar equations reproduced in both gauges, up to each equation's normalization | lane step `b` (R-C): 8/8 `proved-equal`, with controls | verified |
| The tensor-mode equation reproduced exactly | lane step `a2`: `identical` | verified |
| A torsion field can be declared and split | Probe D and `f9` | verified |
| "We have not yet derived torsion perturbation equations with it" | the programme's state | true |
| Environment versions | ``System`$Version`` and each package's qualified `$Version`, printed in one kernel on 2026-09-18; `f11` prints the kernel's | verified: 14.3.0, xTensor 1.3.0, xPerm 1.2.4, xPert 1.0.6, xPand 0.4.4. (The lane's `RC_WOLFRAM_VERSION` line had printed xPand's version, not the kernel's, until this check found it — `PACKAGE_FACTS.md`) |
| Item 1, as reproduced | attachment 1: 2 controls ok, 2 reported reproduced | verified |
| Item 1: `ToMetric` at `:2566`; xTensor's guard at `xTensor.m:8951` | source | verified |
| Item 2, as reproduced | attachment 2: 2 controls ok, 2 reported reproduced | verified |
| Item 2: `DefMetric[1, h, …]` at `:1754` | source | verified |
| Item 3, as reproduced | attachment 3: 2 controls ok, 2 reported reproduced | verified |
| Item 3: the mechanism (`:2881`, `:2867`, `:2840`, `:2907`) | source; `f9`: the classifier's values, the `Null` in the pre-built rules with its messages | verified |
| Item 4: the `-W0` sign, as reproduced | attachment 4: 1 control ok, 1 reported reproduced | verified |
| Item 4: `-1` hard-coded at `:3042` and `:3014` | source | verified |
| Item 5: the paper prints `-2 D D phi`, `E' H` and `4 H' phi'` | the paper's TeX, `literature/1302.6174/xPand_-_arXiv2.tex:1789, 1809-1811` | verified |
| Item 5: the package gives `-2 D D psi`, `2 H E'`, `4 H' phi` | the author's stored outputs in `Examples/8 Implementation Of Paper.nb`; our `b` and `a2` runs of 0.4.4 against the textbook equations | verified (the third was missed by R-C and recorded on #584 now) |
