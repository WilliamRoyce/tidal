# Report for the xBrauer author — partly filed (2026-09-18)

<!-- cspell:words xBrauer xMAG xPand TraceFree xTras xPert xTensor Helpin ContractMetric ContractMetric1 SeparateMetric SeparateMetric1 SeparateMetric2 MetricOfTensor MasterOf FirstMetricQ InducedMetricQ FrozenMetricQ NormalVectorOf DefProjectedTensor SetSlicing IndicesDown ScalarsOfProduct SubValues UpValues DownValues epsilonh normu SetSlicing ToxPandFromRules Riemann-Cartan Unsets upvalue epsilong -->

**Status 2026-09-18.** The user filed the label-slot and bracket findings as
[xBrauer_Bundle #2](https://github.com/THelpin/xBrauer_Bundle/issues/2), correct as filed. The
contraction and separation regressions below are drafted for filing as
`upstream/xbrauer_issue_contraction.md`, reproduced there with xTensor alone. The contraction
mechanism below is **confirmed** one variable at a time by `wolfram/f10_xbrauer_mechanism.wls`;
the plan for the 2026-09-18 follow-up had withdrawn it (in the plan only; never committed or
filed) on a tensor that also carried `ProjectedWith`, which contracts through its own rule. This file stays as the evidence record.

Drafted by research lane R-C (#567) on 2026-09-17. xBrauer arrived in this project only as a
dependency of xMAG, so the companion draft `xmag_upstream_issue.md` is addressed to the same
author and the two can be sent together. The user decided at planning that anything found in
a third-party package is drafted here and filed or sent only with their agreement (planning
record §8, Q4).

Versions, all unmodified: **xBrauer 1.1.0 (2023-11-29)** from `THelpin/xBrauer_Bundle@48be67e1`
with `BrauerAlgebra 1.1.0` and `SymmetricFunctions 1.0.0`, on Wolfram 14.3.0, xAct 1.3.0
(xTensor 1.3.0, xPert 1.0.6, xTras 1.4.2). The affected package is **xPand 0.4.4**
(Pitrou–Roy–Umeh). Evidence: `scripts/research/perturbations/wolfram/f8_contractmetric.wls`,
transcript `third_party/perturbations_runs/f8/20260917T114838Z/transcript.txt`. That script
measures one kernel in three states — xPand alone, then `Needs["xAct`xBrauer`"]`, then
`Needs["xAct`xMAG`"]` on top — so each item below is attributed by measurement rather than by
reading, and §4 records what is *not* xBrauer's.

## 1. `ContractMetric` no longer contracts an induced metric into a tensor that has no `Master`

**Reproduction.** In a kernel with xPand's flat-FLRW slicing (`SetSlicing[g,n,-1,h,cd,{"|","D"},"FLFlat"]`
then `DefMetricFields[g,dg,h,eps]`), `Evh` is xPand's transverse vector perturbation, declared
through `DefProjectedTensor` and therefore carrying no `Master`:

```
ContractMetric[h[ia, ib] Evh[LI[1], LI[0], -ia]]
```

| state | result | sentinel |
| --- | --- | --- |
| xPand alone | `Evh[LI[1], LI[0], ib]` | `F8_BASE_S1_P1={"no-h", ...}` |
| after `Needs["xAct`xBrauer`"]` | `Evh[LI[1], LI[0], -ia] h[ia, ib]`, unchanged | `F8_XBRAUER_S1_P1={"HAS-H", ...}` |

Both slot orders are affected (`F8_XBRAUER_SLOT1_CONTRACTS=False`,
`F8_XBRAUER_SLOT2_CONTRACTS=False`), and contraction with the **first** metric `g` still works
(`F8_VERDICT_G_STILL_CONTRACTS=True`).

**Mechanism.** `xBrauer.m:1938-1940` `Unset`s three of xTensor's `ContractMetric1` rules and
`:1947-1974` re-installs twelve carrying an extra condition, for the tensor slot

```
MetricOfTensor[tensor] === metric || (MetricOfTensor[tensor] === Null && FirstMetricQ[metric])
```

and the same shape with `MasterOf[covd]` for the derivative slot. Measured in the kernel:
`SubValues[ContractMetric1]` goes from **20 to 28**, and of the rules present afterwards
**12 mention `FirstMetricQ`**, 8 `MetricOfTensor` and 4 `MasterOf`
(`F8_SUB_BASE_CONTRACTMETRIC1={20,...}`, `F8_SUB_XBRAUER_CONTRACTMETRIC1={28,...}`,
`F8_GUARDCOUNT_XBRAUER={8, 4, 12, 4}`).

For xPand the condition is False on both legs:

| quantity | value | why |
| --- | --- | --- |
| `$Metrics` | `{g, h, gah2}` | xPand appends the induced metric after the ambient one (`xPand.m:1754`) |
| `FirstMetricQ[h]` | `False` | `h` is not `FirstMetricOfVBundle` |
| `MetricOfTensor[Evh]`, `MetricOfTensor[Eth]` | `Null` | a `DefProjectedTensor` field has no `Master` |
| the whole condition for `h` and `Evh` | **`False`** | `F8_GUARD_XBRAUER_VALUE_FOR_H_AND_EVH=False` |

**The discriminator.** `epsilonh`, xTensor's epsilon of the same induced metric, *does* carry
`Master -> h`. With the same metric, the same slot and the same pipeline it still contracts:

| tensor | `MasterOf` | `MetricOfTensor` | condition | `ContractMetric[h[ia,ib] T[-ia,...]]` |
| --- | --- | --- | --- | --- |
| `epsilonh` | `h` | `h` | `True` | contracts |
| `Evh` | `Null` | `Null` | `False` | unchanged |

(`F8_GUARD_XBRAUER_MASTEROF={h, g, Null, Null, h}` for `{cd, CD, Evh, Eth, epsilonh}`;
`F8_XBRAUER_E1_P1={"no-h",...}` against `F8_XBRAUER_S1_P1={"HAS-H",...}`;
`F8_VERDICT_MASTER_IS_THE_DISCRIMINATOR=True`.)

**Two rival explanations excluded.** `ContractMetric[expr, h]`, with the metric named so the
changed public driver at `:1906` is bypassed, fails identically
(`F8_VERDICT_DRIVER_EXCLUDED=True`); and the residual still contains an **uncontracted** `h`,
so this is the contraction not happening rather than a rule-ordering change downstream
(`F8_VERDICT_RULE_ORDERING_EXCLUDED=True`).

**Consequence for a package loaded earlier.** xPand installs the transversality and
tracelessness of its vector and tensor perturbations as `UpValues` that match an
already-contracted pattern. Those `UpValues` are intact — the counts for `Eth` and `Evh` are
42 and 42 in every state, and an already-contracted argument still goes to `0`
(`F8_VERDICT_UPVALUES_ALIVE=True`) — but they never see a contracted pattern any more, so

```
ToCanonical[ContractMetric[h[ia, ib] cd[-ib]@Evh[LI[1], LI[0], -ia]]]
```

returns `0` before xBrauer and a nonzero expression carrying `h` afterwards
(`F8_BASE_D_P4` against `F8_XBRAUER_D_P4`). Every 3+1 split in that kernel then carries terms
that are zero by transversality. The effect is invisible: no message, no failure, only a
larger answer.

**Suggested fix, minimal.** Let the `Null` fallback accept an induced metric as well as the
first one, i.e.

```
(MetricOfTensor[tensor] === Null && (FirstMetricQ[metric] || InducedMetricQ[metric]))
```

with the corresponding change for `MasterOf[covd]`. That keeps the behavior xBrauer needs for
frozen metrics while restoring xTensor's for a tensor that names no metric of its own.

## 2. `SeparateMetric` and, through it, xPand's `IndicesDown` stop separating a derivative index

`xBrauer.m:1866` `Unset`s xTensor's `SeparateMetric1` catch-all (`xTensor.m:8682`, which
returned `expr` unchanged for a metric that is neither the first nor the one asked for) and
dispatches induced and frozen metrics to `SeparateMetric2` instead; `:1874` then replaces
`SeparateMetric2`, whose derivative branch tests only `MasterOf[covd] === metric` with **no**
`Null` fallback of the kind the tensor branch has.

| expression | xPand alone | after xBrauer |
| --- | --- | --- |
| `SeparateMetric[g][cd[ia]@Esh[LI[1],LI[0]], ia]` | `h[ia,is] cd[-is][Esh[...]]` | unchanged input |
| `SeparateMetric[g][epsilonh[ia,-ib,-ic], ia]` | `epsilonh[-is,-ib,-ic] g[ia,is]` | unchanged input |
| `SeparateMetric[g][Evh[LI[1],LI[0],ia], ia]` | separated | still separated |
| `IndicesDown[cd[ia]@Esh[LI[1],LI[0]]]` | `h[ia,is] cd[-is][Esh[...]]` | unchanged input |

(`F8_SEP_*`, `F8_IDOWN_*`.) The `Master`-less tensor still separates because its fallback
fires with `FirstMetricQ[g] = True`; the derivative and `epsilonh` cases do not, because
`MasterOf[cd] = h` and `MetricOfTensor[epsilonh] = h`, neither of which equals `g` and neither
of which is `Null`. xPand's `IndicesDown` (`xPand.m:2412`) is built on `SeparateMetric` and
inherits the change. `DownValues[SeparateMetric]` also goes from 0 to 1 and
`SubValues[SeparateMetric2]` from 2 to 5.

## 3. A one-character defect at `xBrauer.m:1958`

Four of the re-installed rules use `ScalarsOfProduct`. In three of them the guard is a separate
conjunct; in the fourth it is swallowed into the argument:

```
:1951   && ScalarsOfProduct[prod][metric[b,c]] && (MetricOfTensor[tensor]===metric || ...);
:1958   && ScalarsOfProduct[prod][metric[b,c]  && (MetricOfTensor[tensor]===metric || ...)];
```

Read from the installed source in the kernel, the character following `metric[b,c]` across the
four sites is `{"]", "&", "]", "]"}` (`F8_XBRAUER_SOURCE_SCALARSOFPRODUCT_SHAPES={4, {"]", "&", "]", "]"}}`,
`F8_XBRAUER_1958_PAREN_DEFECT_PRESENT=True`). That rule — the product form, first index of the
metric — therefore cannot evaluate its condition as the other five do.

## 4. What is *not* xBrauer's

Recorded because this lane first attributed all of it to xMAG, and the three-state measurement
separates them:

- the load-time `$RiemannSign = -1` is **xMAG's** (`xMAG.m:110-111`). xBrauer leaves all five
  of xTensor's sign globals at `+1` (`RC_SIGNS_XBRAUER={1, 1, 1, 1, 1}`).
- the replacement of `EinsteinToRicci` is **xMAG's**: its `DownValues` count is 5 with the same
  hash before and after xBrauer, and 7 with a different hash only after xMAG loads.
- xMAG adds nothing of its own to the contraction behavior: `SubValues[ContractMetric1]` has
  the identical count **and hash** in the xBrauer and xMAG states
  (`F8_XMAG_ADDS_NOTHING_TO_CONTRACTMETRIC1=True`,
  `F8_CONTRACTMETRIC1_SUBVALUE_HASHES={"e03874d3b655", "e03874d3b655"}`), and every probe agrees
  expression by expression (`F8_XMAG_PROBES_MATCH_XBRAUER={True, True, True, True}`).
- `AutomaticRules` are undamaged: xPand's `g h -> h` rule still fires
  (`F8_XBRAUER_AUTOMATICRULES_SURVIVE=True`), so this is not a `$Rules` problem.
- xBrauer also replaces `xAct`xTensor`Private`NormalVectorOf` (`:1899-1901`) and
  `xAct`xTensor`Private`PRJ` (`:1928-1930`); both were measured as changed and neither is the
  cause of §1.

## Notes

Nothing under `Applications/xAct/` was modified at any point; every measurement above is from
the installed bytes. Our own workaround is caller-side and does not touch the package.
