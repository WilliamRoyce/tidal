<!-- To post as a comment on https://github.com/THelpin/xMAG/issues/2 . Everything below this
comment is the text to paste. Post it after the xBrauer_Bundle contraction issue is filed, and
replace LINK-TO-CONTRACTION-ISSUE with that issue's URL. Every snippet was run as written by
check_snippets.py; the claims table is claims.md. -->
A few corrections after re-running every snippet in a fresh kernel.

**Item 1 — the snippet as written does not show the change.** On a single line, `Needs[...]; $RiemannSign` reads the name `$RiemannSign` before xTensor is loaded, so it creates an unrelated `Global` symbol and prints that name instead of `1`. This version works however it is entered:

```wolfram
Needs["xAct`xTensor`"]
xAct`xTensor`$RiemannSign    (* 1 *)
Needs["xAct`xMAG`"]
xAct`xTensor`$RiemannSign    (* -1 *)
```

The point stands: xTensor's default is `+1` (`xTensor.m:1837`) and loading xMAG sets it to `-1` (`xMAG.m:110`).

**Item 2 — setup.** The snippet assumes a manifold and metric, for example:

```wolfram
DefManifold[M, 4, {a, b, c, d, e}];
DefMetric[-1, g[-a, -b], CD, {";", "\[Del]"}];
```

**Item 4 — the missing key is `Contorsion`, not `Distortion`, and it is missing in two places.** For a metric-compatible torsionful connection the contortion-trace relation is stored under `ConnectionRelations[covd, Contorsion]` (`:1038`). `DefCovD` installs it under `Distorsi` (`:1044`) and `SetConnectionRelations` under `Distortion` (`:556`); both keys are empty for this connection type, so the relation is never installed:

```wolfram
DefCovD[CDT[-a], {"#", "DT"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True];
ToCanonical[ContorsionCDT[a, -a, -b]]    (* -ContorsionCDT[-b, a, -a], unreduced *)
SetConnectionRelations[CDT];
ToCanonical[ContorsionCDT[a, -a, -b]]    (* -ContorsionCDT[-b, a, -a], still unreduced *)
AutomaticRules[ContorsionCDT, ConnectionRelations[CDT, Contorsion], Verbose -> False];
ToCanonical[ContorsionCDT[a, -a, -b]]    (* -TorsionVectorCDT[-b] *)
```

**Item 4 — a concrete case for the `StringLength` ordering.** With a metric derivative named `CDk` and a torsionful connection `CT` (`Master -> k`), the automatic conversion is attached to `ChristoffelCTCDk`, but `ChangeCurvature[RicciScalarCT[], CT, CDk]` produces `ChristoffelCDkCT`, which is not converted. With `CD` and `CDT` it fires as intended.

**The "Related" note was premature when filed.** The contraction problem is now reported at LINK-TO-CONTRACTION-ISSUE.

**A nuance in the opening:** 13 of the replayed cells matched the stored outputs identically, and 2 matched up to the names of dummy indices.
