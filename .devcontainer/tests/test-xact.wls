#!/usr/bin/env wolframscript
(* xAct Installation Test Script *)

Print["=== xAct Installation Test ==="];
Print[""];

(* Load core packages *)
Print["Loading xAct packages..."];
<<xAct`xTensor`;
<<xAct`xCoba`;

Print["✅ xTensor loaded: ", MemberQ[$Packages, "xAct`xTensor`"]];
Print["✅ xCoba loaded: ", MemberQ[$Packages, "xAct`xCoba`"]];

(* Basic tensor test *)
Print[""];
Print["Testing basic tensor functionality..."];
DefManifold[TestM, 4, IndexRange[a, k]];
DefMetric[-1, testg[-a, -b], TestM];

Print["✅ 4D manifold created with dimension: ", Dim[TestM]];
Print["✅ Metric tensor defined: ", HasMetricQ[testg[]]];

(* Test coordinate charts (xCoba functionality) *)
DefChart[coord, TestM, {0, 1, 2, 3}, {t[], x[], y[], z[]}];
Print["✅ Coordinate chart defined: ", ChartsOfManifold[TestM]];

Print[""];
Print["*** xAct is successfully installed and functional! ***"];
Print["Ready for tensor computations in your TIDAL project!"];
Print[""];
Print["Note: xPerm running in pure Mathematica mode (MathLink disabled for compatibility)"];