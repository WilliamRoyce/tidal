# Draft report for the xMAG author — NOT FILED

<!-- cspell:words xMAG xBrauer TraceFree xPand xTras xPert Helpin ToDistortion BreakDistortion BreakContorsion ToContorsion DefConnectionPerturbation StartInducedDecomposition contortion DefCovD ChangeCurvature Riemann-Cartan Accelerationn CDCDL Tetrahh epsilonhh -->

Drafted by research lane R-C (#567) on 2026-09-16. The user decided at planning that
defects found in xMAG are drafted as upstream issues and filed, or sent to the author, only
with the user's agreement (planning record §8, Q4 and round 4 comment 11). Everything below
was observed with `THelpin/xMAG@88026e47` (0.1.0, 2023-05-18), `THelpin/xBrauer_Bundle@48be67e1`
(xBrauer 1.1.0), `xAct-contrib/TraceFree@4e53ab39` (0.1.0), on Wolfram 14.3.0 and xAct 1.3.0
(xTensor 1.3.0, xTras 1.4.2, xPert 1.0.6), unmodified; script
`scripts/research/perturbations/wolfram/probe_e_xmag.wls` and the cross-check
`probe_d_limits.wls` (run directories under `third_party/perturbations_runs/e/` and `dl/`).

## What worked

- Loads cleanly beside xPand (2 s, no messages); the version gates on xTensor 1.2.0 and
  xTras 1.0.6 are satisfied by 1.3.0 and 1.4.2.
- `DefCovD[CDT[-a], {"#","DT"}, Torsion -> True, FromMetric -> g, Master -> g,
  ConnectionRelations -> True]` defines the Riemann-Cartan connection with
  `ContorsionCDT`, `TorsionVectorCDT`, `TFContorsionCDT`, `TFTorsionCDT`,
  `PerturbationChristoffelCDT`, `ChristoffelCDCDT` and the automatic relation
  `ChristoffelCDCDT = -ContorsionCDT`.
- `BreakContorsion[ContorsionCDT[a,-b,-c], CDT]` gives
  `½ (T^a_bc + T_b^a_c + T_c^a_b)`, the standard contortion (checked: its antisymmetric part
  is the torsion and it is metric-compatible).
- `ToDistortion[RicciScalarCDT[], CDT, g]` runs in 0.03 s and returns Levi-Civita
  curvature plus contortion terms; `BreakContorsion` of that result is free of every
  `CDT` object.
- `DefConnectionPerturbation[ChristoffelCDT, dGam, eps]` coexists with xPert's
  `DefMetricPerturbation[g, dg, eps]`: both `dGam[LI[1], a, -b, -c]` and `dg[LI[1], -a, -b]`
  are available afterwards.

## 1. The torsion part of `ToDistortion[RicciScalarCDT[]]` has the opposite sign to xTensor's `ChangeCurvature`

xMAG (after `BreakContorsion`, canonical form):

```
R̃ = R − ¼ T_abc T^abc − ½ T_abc T^bac − T^a_a^b T^c_bc + 2 ∇_b T^a_a^b
```

xTensor alone, `ChangeCurvature[RicciScalarCDT[], CDT, CD]` followed by the substitution
`ChristoffelCDCDT -> −K` with `K^a_bc = ½ (T^a_bc + T_b^a_c + T_c^a_b)` — the sign fixed by
xTensor's own `ChangeCovD[CDT[-b]@v[a], CDT, CD] − CD[-b]@v[a] = −ChristoffelCDCDT^a_bs v^s`
(the relation xMAG also asserts):

```
R̃ = R + ¼ T_abc T^abc + ½ T_abc T^bac + T^a_a^b T^c_bc − 2 ∇_b T^a_a^b
```

The two differ by the sign of every torsion term (`RC_DL_XMAG_MINUS_CORRECT`). For the
vectorial torsion `T^a_bc = δ^a_c v_b − δ^a_b v_c` the second form gives
`R̃ = R − 6 v² + 6 ∇·v`, the textbook Einstein-Cartan result (`R − (2/3) T_μ T^μ + 2 ∇_μ T^μ`
with `T_μ = 3 v_μ`); the first gives `R + 6 v² − 6 ∇·v`. Either xMAG uses a Riemann sign
convention for independent connections that differs from xTensor's default for the same
symbol `RicciScalarCDT`, or there is a sign slip in `ToDistortion`'s curvature rule; the
package's documentation does not say which. A one-line statement of the convention in
`ToDistortion::usage`, or a test against `ChangeCurvature` for the metric-compatible case,
would settle it for users.

## 2. `BreakDistortion` returns `Null` for a metric-compatible connection

With `CDT` as above (torsion, no non-metricity), `BreakDistortion[ToDistortion[RicciScalarCDT[], CDT, g], CDT, g]`
returns `Null` (with `Validate::inhom`), while `BreakContorsion[…, CDT]` works. The usage text
of `BreakDistortion` does not exclude that case; either dispatching to `BreakContorsion` when
`NonMetricityQ[CDT]` is `False`, or a message, would help.

## 3. A second torsion connection breaks `ToDistortion` on the first

After a further plain `DefCovD[CDL[-a], {"&","DL"}, FromMetric -> g, Torsion -> True]`,
`ToDistortion[RicciScalarCDT[], CDT, g]` returns `Null` with two `Validate::inhom`
messages (`RC_E_TODISTORTION_AFTER_SECOND_TORSION_CONNECTION`). One torsion connection per
session appears to be an undocumented assumption.

## 4. `StartInducedDecomposition` does not run on this bundle

`StartInducedDecomposition[g, CDG, {{";h","Dh"},{"%h","DGh"}}, {nn, hh}]` on a general
connection (`FromMetric -> Null, Master -> g, Torsion -> True`), in a kernel with or without
xPand's slicing, stops with `DefCovD::invalid`, `TorsionQ::unknown`, `MakeRule::inhom` after
defining `nn`, `hh`, `epsilonhh`, `Tetrahh`, `ExtrinsicKh`, `Accelerationn`
(`probe_e2_xmag_induced.wls`, `RC_E2_STARTINDUCED_ON_GENERAL_CONNECTION_OK=False`); the same on
the Riemann-Cartan connection. It may be an API drift against xTensor 1.3.0 (the package
expects 1.2.0); the tutorial notebook's call form was used.

## 5. `ChangeCurvature` on a plain torsion connection returns its input with xMAG loaded

With xMAG loaded, `ChangeCurvature[RicciScalarCDL[], CDL, CD]` for a connection defined
without `Master` returns `RicciScalarCDL[]` unchanged, whereas the same call in a kernel
without xMAG expands into `RicciCD`, `ChristoffelCDCDL` and `g` terms. xMAG's overrides of
the curvature-relation functions for torsionful connections (`xMAG.m:1138-1162`) seem to
catch that case.
