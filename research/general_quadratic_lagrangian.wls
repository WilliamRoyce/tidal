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

(* ============================================================ *)
Print["\n\n== PART F: MIXING CLASSIFICATION (B0!=0, Tbar=0) =="];
Print["Automated: count tensor types per sector, classify mixing"];
(* ============================================================ *)

(* Classification rules for B0!=0, Tbar=0 linearization:
   - Count R (curvature) factors → h perturbation
   - Count T (torsion) factors → t perturbation
   - Count F (EM) factors → one can be background F_bar=B0
   A cubic X*Y*Z contributes at O(eps^2) when one factor = background.
   For t-f mixing: need >=1 torsion AND >=1 EM, with one EM -> B0.
   For h-f mixing: need >=1 curvature AND >=1 EM.
   For h-t mixing via B0: need curvature AND torsion AND >=1 EM.
*)

ClassifySector[name_String, count_Integer] := Module[
  {hasR, hasT, hasF, nR, nT, nF, mixing},

  (* Count tensor types from sector name *)
  nR = StringCount[name, "R"] - StringCount[name, "DR"];
  nR = nR + StringCount[name, "DR"];  (* DR still contains curvature *)
  hasR = StringContainsQ[name, "R"];
  hasT = StringContainsQ[name, "T" ~~ Except["x"]];  (* T but not "x" *)
  hasT = hasT || StringContainsQ[name, "DT"];
  hasF = StringContainsQ[name, "F"];

  (* More precise: parse the sector name *)
  hasR = StringContainsQ[name, RegularExpression["(^|\\s|x )D?R"]];
  hasT = StringContainsQ[name, RegularExpression["(^|\\s|x )D?T(?!x)"]];
  hasF = StringContainsQ[name, RegularExpression["(^|\\s|x )D?F"]];

  (* Count F factors *)
  nF = StringCount[name, RegularExpression["(?<![D])F"]];

  mixing = Which[
    count == 0, "none (vanishing)",
    hasR && hasT && hasF, "h-t mixing via B0",
    hasR && hasF && !hasT, "h-f mixing (Gertsenshtein-type)",
    hasT && hasF && !hasR && nF >= 2, "t-f mixing via B0 (cubic)",
    hasT && hasF && !hasR, "t-f mixing",
    hasR && hasT && !hasF, "h-t mixing (B0-independent)",
    hasR && !hasT && !hasF, "curvature self-interaction",
    hasT && !hasR && !hasF, "torsion self-interaction",
    hasF && !hasR && !hasT, "EM self-interaction",
    True, "unclassified"
  ];

  (* Special case: T x T x F with Tbar=0 is t^2*B0 = mass correction *)
  If[StringMatchQ[name, "*T x T x F*"] && !hasR,
    mixing = "torsion mass correction (t^2*B0)"
  ];

  {name, count, mixing}
];

classifiedResults = ClassifySector[#[[1]], #[[2]]] & /@ allResults;

Print["\n  Mixing Classification:"];
Do[
  Print["  ", r[[1]], " (", r[[2]], "): ", r[[3]]],
  {r, classifiedResults}
];

(* ============================================================ *)
Print["\n\n== PART G: CHERN-SIMONS (A-POTENTIAL) COUPLINGS =="];
Print["Gauge-invariant up to boundary: eps V A F"];
(* ============================================================ *)

(* Define EM potential as rank-1 tensor *)
DefTensor[AA[-a], M4, PrintAs -> "A"];

(* The Chern-Simons form: eps^{mnrs} V_m A_n F_{rs}
   where V is a torsion-derived vector.
   We enumerate: eps x TorVector x A x F *)

(* Torsion trace vector: T_mu = T^a_{a mu} — represent as rank-1 *)
DefTensor[Tvec[-a], M4, PrintAs -> "Tv"];

(* Enumerate: eps x Tvec x A x F (4+1+1+2=8 indices) *)
AppendResult @ EnumerateSector["CS: eps Tvec A F",
  eps[-a,-b,-c,-d] * Tvec[-e] * AA[-f] * FF[-g,-h]];

(* Also: eps x Tor x A x F (4+3+1+2=10 indices) — full torsion, not just trace *)
AppendResult @ EnumerateSector["CS: eps Tor A F",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * AA[-h] * FF[-i,-j]];

(* And: eps x DT x A x F with derivatives *)
AppendResult @ EnumerateSector["CS: eps Tor A DF",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * AA[-h] * DF[-i,-j,-k]];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  FINAL MASTER SUMMARY"];
Print["========================================"];
(* ============================================================ *)

Do[
  Print["  ", r[[1]], ": ", r[[2]]],
  {r, allResults}
];

(* Export full classified results *)
classifiedData = Association @@ (
  (#[[1]] -> <|"count" -> #[[2]], "mixing" -> #[[3]]|>) & /@ classifiedResults
);

(* Add CS results *)
Do[
  If[Length[allResults] >= k && allResults[[k, 2]] > 0,
    classifiedData[allResults[[k, 1]]] = <|
      "count" -> allResults[[k, 2]],
      "mixing" -> "Chern-Simons (A-potential)"
    |>
  ],
  {k, Length[allResults] - 2, Length[allResults]}  (* last 3 entries are CS *)
];

Export[
  FileNameJoin[{DirectoryName[$InputFileName], "enumeration_classified.json"}],
  classifiedData
];

(* Also export simple counts *)
summaryData = (#[[1]] -> #[[2]]) & /@ allResults;
Export[
  FileNameJoin[{DirectoryName[$InputFileName], "enumeration_results.json"}],
  Association @@ summaryData
];
Print["\nResults exported to enumeration_results.json and enumeration_classified.json"];
Print["=== DONE ==="];
