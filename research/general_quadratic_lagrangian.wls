#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  Systematic enumeration of all independent scalar invariants
  built from torsionful curvature, torsion, and EM field strength in 4D.

  Uses xTras AllContractions to enumerate index contractions.
  Covers: quadratic, derivative, cubic, and parity-odd sectors.

  References:
    Barker et al. (2406.12826), Shapiro (hep-th/0103093),
    Obukhov et al. (2410.01355), Hayashi & Shirafuji (1979)

  Date: 2026-03-19
*)

(* === Load packages === *)
Print["=== Loading xAct + xTras ==="];
Needs["xAct`xTensor`"];
Needs["xAct`xTras`"];
xAct`xTensor`$DefInfoQ = False;
xAct`xTensor`$UndefInfoQ = False;
Print["Loaded."];

(* === Manifold & metric === *)
DefManifold[M4, 4, {a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, q, r, s}];
DefMetric[-1, metric[-a, -b], CD, PrintAs -> "g"];

(* === Tensors (use GenSet for multi-slot antisymmetry) === *)

(* Torsionful Riemann: antisym in {1,2} and {3,4}, NO pair-exchange *)
DefTensor[RCurv[-a, -b, -c, -d], M4, GenSet[-Cycles[{1, 2}], -Cycles[{3, 4}]]];
(* Torsion: antisym in last two indices *)
DefTensor[Tor[-a, -b, -c], M4, GenSet[-Cycles[{2, 3}]]];
(* EM field strength: antisym *)
DefTensor[FF[-a, -b], M4, Antisymmetric[{1, 2}]];
(* Levi-Civita *)
DefTensor[eps[-a, -b, -c, -d], M4, Antisymmetric[{1, 2, 3, 4}]];

(* Derivative-level *)
DefTensor[DT[-d, -a, -b, -c], M4, GenSet[-Cycles[{3, 4}]]];    (* nabla_d T_{abc} *)
DefTensor[DF[-c, -a, -b], M4, GenSet[-Cycles[{2, 3}]]];        (* nabla_c F_{ab} *)
DefTensor[DRCurv[-e, -a, -b, -c, -d], M4,
  GenSet[-Cycles[{2, 3}], -Cycles[{4, 5}]]];                    (* nabla_e Rtilde_{abcd} *)

Print["All tensors defined. xTensorQ[RCurv] = ", xTensorQ[RCurv]];

(* === Helper === *)
EnumerateSector[name_String, expr_] := Module[{result, count},
  Print["\n--- ", name, " ---"];
  result = Catch[Quiet[AllContractions[expr], {AllContractions::nocontr}]];
  If[!ListQ[result],
    Print["  0 contractions"];
    Return[{name, 0, {}}]
  ];
  count = Length[result];
  Print["  ", count, " contractions"];
  Do[Print["    ", k, ". ", result[[k]] // InputForm], {k, Min[count, 10]}];
  If[count > 10, Print["    ... (", count - 10, " more)"]];
  {name, count, result}
];

(* === RESULTS COLLECTOR === *)
allResults = {};
AppendResult[r_] := AppendTo[allResults, r];

(* ============================================================ *)
Print["\n\n== PART A: QUADRATIC NON-DERIVATIVE =="];
(* ============================================================ *)

AppendResult @ EnumerateSector["R^2",
  RCurv[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h]];
AppendResult @ EnumerateSector["T^2",
  Tor[-a,-b,-c] * Tor[-d,-e,-f]];
AppendResult @ EnumerateSector["F^2",
  FF[-a,-b] * FF[-c,-d]];
AppendResult @ EnumerateSector["R x F",
  RCurv[-a,-b,-c,-d] * FF[-e,-f]];
AppendResult @ EnumerateSector["R x T (odd)",
  RCurv[-a,-b,-c,-d] * Tor[-e,-f,-g]];
AppendResult @ EnumerateSector["T x F (odd)",
  Tor[-a,-b,-c] * FF[-d,-e]];

(* ============================================================ *)
Print["\n\n== PART B: QUADRATIC ONE-DERIVATIVE =="];
(* ============================================================ *)

AppendResult @ EnumerateSector["DT x F",
  DT[-d,-a,-b,-c] * FF[-e,-f]];
AppendResult @ EnumerateSector["T x DF (IBP)",
  Tor[-a,-b,-c] * DF[-d,-e,-f]];
AppendResult @ EnumerateSector["R x DT",
  RCurv[-a,-b,-c,-d] * DT[-e,-f,-g,-h]];
AppendResult @ EnumerateSector["DR x T (IBP)",
  DRCurv[-a,-b,-c,-d,-e] * Tor[-f,-g,-h]];
AppendResult @ EnumerateSector["DT x DT",
  DT[-a,-b,-c,-d] * DT[-e,-f,-g,-h]];
AppendResult @ EnumerateSector["DF x DF",
  DF[-a,-b,-c] * DF[-d,-e,-f]];
AppendResult @ EnumerateSector["DT x T (odd)",
  DT[-a,-b,-c,-d] * Tor[-e,-f,-g]];
AppendResult @ EnumerateSector["DF x F (odd)",
  DF[-a,-b,-c] * FF[-d,-e]];

(* ============================================================ *)
Print["\n\n== PART C: CUBIC NON-DERIVATIVE =="];
Print["(relevant when linearizing around B0 background)"];
(* ============================================================ *)

AppendResult @ EnumerateSector["R x F x F",
  RCurv[-a,-b,-c,-d] * FF[-e,-f] * FF[-g,-h]];
AppendResult @ EnumerateSector["T x T x F",
  Tor[-a,-b,-c] * Tor[-d,-e,-f] * FF[-g,-h]];
AppendResult @ EnumerateSector["R x T x T",
  RCurv[-a,-b,-c,-d] * Tor[-e,-f,-g] * Tor[-h,-i,-j]];
AppendResult @ EnumerateSector["R x R x F",
  RCurv[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * FF[-i,-j]];
AppendResult @ EnumerateSector["T x F x F (odd)",
  Tor[-a,-b,-c] * FF[-d,-e] * FF[-f,-g]];
AppendResult @ EnumerateSector["R x T x F (odd)",
  RCurv[-a,-b,-c,-d] * Tor[-e,-f,-g] * FF[-h,-i]];

(* ============================================================ *)
Print["\n\n== PART D: CUBIC WITH ONE DERIVATIVE =="];
Print["(derivative fixes odd index counts for cubic sectors)"];
(* ============================================================ *)

AppendResult @ EnumerateSector["DT x F x F",
  DT[-a,-b,-c,-d] * FF[-e,-f] * FF[-g,-h]];
AppendResult @ EnumerateSector["T x DF x F",
  Tor[-a,-b,-c] * DF[-d,-e,-f] * FF[-g,-h]];
AppendResult @ EnumerateSector["R x DT x F",
  RCurv[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * FF[-i,-j]];

(* ============================================================ *)
Print["\n\n== PART E: PARITY-ODD (with epsilon) =="];
(* ============================================================ *)

(* Quadratic parity-odd *)
AppendResult @ EnumerateSector["eps R^2",
  eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * RCurv[-i,-j,-k,-l]];
AppendResult @ EnumerateSector["eps T^2",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * Tor[-h,-i,-j]];
AppendResult @ EnumerateSector["eps F^2",
  eps[-a,-b,-c,-d] * FF[-e,-f] * FF[-g,-h]];
AppendResult @ EnumerateSector["eps R x F",
  eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * FF[-i,-j]];

(* Derivative parity-odd *)
AppendResult @ EnumerateSector["eps DT x F",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * FF[-i,-j]];
AppendResult @ EnumerateSector["eps T x DF",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * DF[-h,-i,-j]];

(* Cubic parity-odd *)
AppendResult @ EnumerateSector["eps T x T x F",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * Tor[-h,-i,-j] * FF[-k,-l]];
AppendResult @ EnumerateSector["eps DT x F x F",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * FF[-i,-j] * FF[-k,-l]];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  MASTER SUMMARY"];
Print["========================================"];
(* ============================================================ *)

Do[
  Print["  ", r[[1]], ": ", r[[2]]],
  {r, allResults}
];

(* Export summary to JSON for programmatic use *)
summaryData = (#[[1]] -> #[[2]]) & /@ allResults;
Export[
  FileNameJoin[{DirectoryName[$InputFileName], "enumeration_results.json"}],
  Association @@ summaryData
];
Print["\nResults exported to enumeration_results.json"];
Print["=== DONE ==="];
