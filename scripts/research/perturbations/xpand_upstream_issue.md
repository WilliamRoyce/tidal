# Draft upstream report for xPand 0.4.4 — not yet sent

<!-- cspell:words xPand xCoba xMAG SplitPerturbations ToxPandFromRules DefChart ExtractComponents ToxPand SplitMatter normu Pitrou changetoinducedmetric MetricCovDQ FrozenMetricQ MasterOf EinsteinToRicci epsilong cdpost cdpre -->

**Status 2026-09-18.** xPand has no issue tracker; the report goes to the author by email. The
text to send is `upstream/xpand_email.txt`, with one script per item
(`upstream/xpand_item{1,2,3,4}_*.wls`), checked claim by claim in `upstream/claims.md`. It adds
one item not in this draft — a background rule whose value is a sum of terms along the normal
breaks the preparation of split rules (#591, `wolfram/f9_background_rules.wls`) — and a third
misprint in the paper's Appendix A. This file stays as the evidence record.

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
xPand's; it is mentioned only so the two are not confused. A three-state measurement has
since confirmed that attribution: `EinsteinToRicci` gains `DownValues` and a new hash only
when xMAG loads, and is byte-identical before and after xBrauer
(`f8_contractmetric.wls`, run `f8/20260917T114838Z`). A **third**, unrelated symptom in the
same area — `ContractMetric` refusing to contract xPand's induced metric into its own
projected fields — is xBrauer's and is reported separately in `xbrauer_upstream_issue.md`.
Nothing in that one is xPand's either.

**Reproduction scripts.** `scripts/research/perturbations/wolfram/probe_d_limits.wls`
(sentinels `DL_SANITY_METRIC_SPLIT_7_BEFORE_CHART=ok`,
`DL_SANITY_METRIC_SPLIT_8_AFTER_CHART=BROKEN`) and
`scripts/research/perturbations/wolfram/f2_hazards.wls`, which adds the workaround and its
controls (`F2_I_SPLIT_AFTER_CHART=BROKEN`, `F2_I_GUARDED_SPLIT_AFTER_CHART=ok`,
`F2_I_GUARDED_SPLIT_EQUALS_REFERENCE=identical`, `F2_I_CONTROL_WRONG_GUARD=BROKEN`).

## 2. `SetSlicing` hard-codes the slice determinant sign, so `epsilon[h]` is wrong in `(+,-,-,-)`

**Where.** `xPand.m:1754` declares the induced metric as

```
DefMetric[1, h[-ind1,-ind2], cd, {cdpost, cdpre}, InducedFrom -> {g, u}, PrintAs -> ...]
```

The first argument of `DefMetric` is the sign of the metric's determinant, and it is the
literal `1` regardless of `normu`. With `normu = -1` (mostly plus) that is right: the induced
metric is positive definite. With `normu = +1`, which `SetSlicing` accepts and which is the
setting a mostly-minus project uses, the induced metric is `h_ab = g_ab - n_a n_b`, negative
definite, and the true determinant sign is `-1`.

**Consequence.** `SignDetOfMetric[h]` feeds xTensor's product rule for the epsilon of an
induced metric (`xTensor.m:8134-8137`), so `epsilon[h]` is normalized with the wrong sign, and
it disagrees with the four-index epsilon contracted with the normal. Measured in both
settings, in a kernel that does nothing but slice:

| | `normu = -1` | `normu = +1` |
| --- | --- | --- |
| `n[a] n[-a]` | `-1` | `+1` |
| `SignDetOfMetric[h]`, as declared | `1` | `1` |
| the true sign of `det h` | `1` | `-1` |
| `epsilon[h][-a,-b,-c] epsilon[h][a,b,c]` | `6` | `6` |
| `epsilon[g][-a,-b,-c,-d] n[d] epsilon[g][a,b,c,e] n[-e]` | `6` | `-6` |

Everything else in the mostly-minus session behaves: `n.n` comes out `+1`, the slice trace is
`3`, the epsilon is orthogonal to `n` on every slot and totally antisymmetric, and the
four-index epsilon squares to `-24` in both settings. Only the slice epsilon's normalization
moves, and only by the sign the hard-code gets wrong. Nothing else in our work is affected,
because the parity-even sector never uses `epsilon[h]`.

**Suggested fix.** Pass the computed sign instead of the literal, for example
`DefMetric[If[normu === 1, -1, 1], h[-ind1,-ind2], cd, ...]`, or more generally
`-SignDetOfMetric[g] normu`. A one-line change; the alternative, for a caller, is to work with
`epsilon[g][-a,-b,-c,-d] n[d]` throughout, which is what we do.

**Reproduction.** `scripts/research/perturbations/wolfram/f7_epsilon_and_import_map.wls`,
sentinels `F7A_SIGNDET_H_DECLARED_VS_TRUE`, `F7A_RHS_SQUARE`, `F7A_SQUARE_CONSISTENT`,
`F7A_XPAND_1754_DETERMINANT_HARDCODE_BITES`; runs `f7/20260917T115539Z` (`normu = -1`) and
`f7/20260917T115737Z` (`normu = +1`).

## 3. Two places hard-code the normal's norm `n.n = -1`

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

## 4. Two misprints in the paper's Appendix A (arXiv:1302.6174)

`:1789` prints `-2 D_a D^a phi` where the package (and the textbook) give the Laplacian of
the curvature potential `psi`; `:1809` prints the tensor friction as `E' H` where the
package gives `2 H E'`. Recorded so that readers do not type the printed outputs as targets.
