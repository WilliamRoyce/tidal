# Draft upstream report for xPand 0.4.4 — NOT FILED

<!-- cspell:words xPand xCoba xMAG SplitPerturbations ToxPandFromRules DefChart ExtractComponents ToxPand SplitMatter normu Pitrou changetoinducedmetric MetricCovDQ FrozenMetricQ MasterOf EinsteinToRicci -->

Drafted by research lane R-C (#567) on 2026-09-16 under the same rule as
`docs/cosmology/psalter_543_upstream_issue.md`: written for the maintainers, kept in the
repository, filed only if the user decides to. Contact: the xPand page (Cyril Pitrou, IAP).
Everything below was observed on Wolfram 14.3.0, xAct 1.3.0 (xTensor 1.3.0, xPert 1.0.6),
xPand 0.4.4 installed from `http://www2.iap.fr/users/pitrou/xPand_0.4.4.tar.gz`
(sha256 `26e7abcac7bb655235ec39b73850729cf4465748d0d5ea2ad03c0607aef5ceab`), unmodified.

## 1. `SplitPerturbations` returns `Null` once a metric-less covariant derivative exists

**Steps.** Standard flat-FRW slicing (`SetSlicing[g, n, -1, h, cd, {"|","D"}, "FLFlat"]`),
`DefMetricPerturbation`, `DefMetricFields`; a split such as
`ToxPandFromRules[RicciScalarCD[], SplitMetric[g, dg, h, "AnyGauge"], h, 1]` works. Then
`Needs["xAct`xCoba`"]; DefChart[Bc, M4, {0,1,2,3}, {tt[],xx[],yy[],zz[]}]`, which defines the
parallel derivative `PDBc` with `MetricOfCovD[PDBc] === Null`. The same split now prints
`InducedFrom::unknown: Unknown metric Null` and returns `Null`.

**Observed.** `$CovDs = {PD, CD, cd, CDah2, CDT, PDBc}`,
`MetricOfCovD /@ $CovDs = {Null, g, h, gah2, g, Null}`; every later split fails the same
way. (A superficially similar symptom from another package is a different site; see the
note below.)

**Expected.** The split should ignore covariant derivatives that do not appear in the
expression, or report which derivative it could not handle.

**Where to look.** `ToMetric` (`xPand.m:2563`) is reached from every `ToxPandFromRules`
through `Conformal` (`:2619`, `:2657`), and at `:2566` it builds

```
$CovDsNotInduced = Select[Rest@$CovDs, InducedFrom[MetricOfCovD[#]] === Null &]
```

`MetricOfCovD` returns `Null` for a derivative defined without `FromMetric`, and
`InducedFrom[Null]` falls through to `xTensor.m:8302`, which throws. xTensor guards the same
pattern in its own code at `xTensor.m:8951`
(`changetoinducedmetric[covd_, metric_, Null] := metric`). The one-line fix would be

```
$CovDsNotInduced = Select[Rest@$CovDs, MetricOfCovD[#] =!= Null && InducedFrom[MetricOfCovD[#]] === Null &]
```

(`MetricCovDQ` at `xTensor.m:6723` is the same test). A caller-side wrapper that filters
`$CovDs` for the duration of the split recovers a result `identical` to the clean-kernel one,
which is how we work around it meanwhile; a wrong-guard control stays broken, so the
attribution is the guard and not the wrapper.

**A different site, same fallback.** A second symptom we first reported here turned out to
come from another package: with xMAG loaded, its replacement of `EinsteinToRicci`
(`xMAG.m:1154`) tests `FrozenMetricQ[MasterOf[#]]`, and for a connection declared without
`Master` that reaches the same `InducedFrom[Null]` fallback from `ToMetric`'s
`preexpression` line (`xPand.m:2569`). That one is xMAG's (and our declaration's), not
xPand's; it is mentioned only so the two are not confused.

**Reproduction scripts.** `scripts/research/perturbations/wolfram/probe_d_limits.wls`
(sentinels `DL_SANITY_METRIC_SPLIT_7_BEFORE_CHART=ok`,
`DL_SANITY_METRIC_SPLIT_8_AFTER_CHART=BROKEN`) and
`scripts/research/perturbations/wolfram/f2_hazards.wls`, which adds the workaround and its
controls (`F2_I_SPLIT_AFTER_CHART=BROKEN`, `F2_I_GUARDED_SPLIT_AFTER_CHART=ok`,
`F2_I_GUARDED_SPLIT_EQUALS_REFERENCE=identical`, `F2_I_CONTROL_WRONG_GUARD=BROKEN`).

## 2. Two places hard-code the normal's norm `n.n = -1`

Not a defect for the package's stated conventions, recorded because users of the
`(+,-,-,-)` signature will meet it:

- `ExtractComponents` multiplies by `-1` for an up index and `+1` for a down index in its
  `"Time"` projection (`xPand.m:3042`), i.e. it assumes `n.n = -1`. With
  `SetSlicing[g, n, +1, ...]` (accepted; only the determinant sign is checked, `:1715`) the
  `"Time"` component of `V0 n^a` comes back as `-V0` while `n_a V^a / (n.n)` gives `+V0`.
- `ToxPand` passes `-1` as the fluid norm to `SplitMatter` (`:3014`); with `n.n = +1` the
  built-in fluid split does not run. `SplitMatter[u, du, +1, h, gauge, order]` followed by
  `ToxPandFromRules` works.

A `normu` read from the slicing (`SetSlicing` already stores it) at both sites would make
the package signature-agnostic; the split itself already is (verified: the minimal example
and the Newtonian-gauge Einstein equations under `n.n = +1` reproduce the `n.n = -1`
results under the map `R -> -R`, `phi -> -phi` for the same textual gauge rule, and
`D_a D^a -> -D_a D^a`).

## 3. Two misprints in the paper's Appendix A (arXiv:1302.6174)

`:1789` prints `-2 D_a D^a phi` where the package (and the textbook) give the Laplacian of
the curvature potential `psi`; `:1809` prints the tensor friction as `E' H` where the
package gives `2 H E'`. Recorded so that readers do not type the printed outputs as targets.
