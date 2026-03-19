#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  SYSTEMATIC enumeration of the most general Lagrangian via MakeContractionAnsatz.

  Unlike the previous manual approach, this script PROGRAMMATICALLY generates
  ALL possible tensor products from the building blocks {R, T, F, DR, DT, DF}
  and calls AllContractions / MakeContractionAnsatz on each.

  This guarantees completeness: no sector can be missed because the loop
  covers all combinations with even index count.
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

(* === Define all field-strength tensors === *)

(* Torsionful Riemann: antisym in {1,2} and {3,4}, NO pair-exchange *)
DefTensor[RCurv[-a, -b, -c, -d], M4, GenSet[-Cycles[{1, 2}], -Cycles[{3, 4}]]];
(* Torsion: antisym in last two indices *)
DefTensor[Tor[-a, -b, -c], M4, GenSet[-Cycles[{2, 3}]]];
(* EM field strength: antisym *)
DefTensor[FF[-a, -b], M4, Antisymmetric[{1, 2}]];
(* Levi-Civita *)
DefTensor[eps[-a, -b, -c, -d], M4, Antisymmetric[{1, 2, 3, 4}]];

(* Derivative-level *)
DefTensor[DT[-d, -a, -b, -c], M4, GenSet[-Cycles[{3, 4}]]];
DefTensor[DF[-c, -a, -b], M4, GenSet[-Cycles[{2, 3}]]];
DefTensor[DRCurv[-e, -a, -b, -c, -d], M4,
  GenSet[-Cycles[{2, 3}], -Cycles[{4, 5}]]];

(* For Chern-Simons: EM potential *)
DefTensor[AA[-a], M4];
(* Torsion trace vector *)
DefTensor[Tvec[-a], M4];

Print["All tensors defined."];

(* === Index pool === *)
(* We have 18 abstract indices. Assign non-overlapping sets to each tensor factor. *)
allIndices = {a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, q, r, s};

(* === Building blocks with metadata === *)
(* Each entry: {name, nIndices, constructorFunction} *)
(* The constructor takes a list of index symbols and returns the tensor expression *)

blocks = {
  {"R",  4, Function[idx, RCurv[-idx[[1]], -idx[[2]], -idx[[3]], -idx[[4]]]]},
  {"T",  3, Function[idx, Tor[-idx[[1]], -idx[[2]], -idx[[3]]]]},
  {"F",  2, Function[idx, FF[-idx[[1]], -idx[[2]]]]},
  {"DT", 4, Function[idx, DT[-idx[[1]], -idx[[2]], -idx[[3]], -idx[[4]]]]},
  {"DF", 3, Function[idx, DF[-idx[[1]], -idx[[2]], -idx[[3]]]]},
  {"DR", 5, Function[idx, DRCurv[-idx[[1]], -idx[[2]], -idx[[3]], -idx[[4]], -idx[[5]]]]}
};

nBlocks = Length[blocks];

(* === Generate tensor product from a list of block indices === *)
(* blockList: list of indices into `blocks`, e.g. {1,1} for R*R *)
MakeProduct[blockList_List] := Module[
  {expr = 1, offset = 0, bl, nIdx, constructor, myIndices},
  Do[
    bl = blocks[[bi]];
    nIdx = bl[[2]];
    constructor = bl[[3]];
    myIndices = allIndices[[offset + 1 ;; offset + nIdx]];
    expr = expr * constructor[myIndices];
    offset += nIdx;
    , {bi, blockList}
  ];
  expr
];

(* Total indices for a combination *)
TotalIndices[blockList_List] := Total[blocks[[#, 2]] & /@ blockList];

(* Sector name *)
SectorName[blockList_List] := StringRiffle[blocks[[#, 1]] & /@ blockList, " x "];

(* === MAIN ENUMERATION LOOP === *)

allResults = {};

EnumerateIfEven[blockList_List, prefix_String: ""] := Module[
  {nIdx, name, expr, result, count},
  nIdx = TotalIndices[blockList];
  If[prefix == "eps", nIdx += 4];

  (* Skip if odd or too many indices for our pool *)
  If[OddQ[nIdx] || nIdx > Length[allIndices], Return[]];

  name = If[prefix == "", SectorName[blockList],
    prefix <> " " <> SectorName[blockList]];

  (* Build expression *)
  If[prefix == "eps",
    expr = eps[-allIndices[[1]], -allIndices[[2]], -allIndices[[3]], -allIndices[[4]]];
    (* Shift the block indices by 4 for eps *)
    Block[{expr2 = 1, offset = 4, bl, n, constructor, myIdx},
      Do[
        bl = blocks[[bi]];
        n = bl[[2]];
        constructor = bl[[3]];
        myIdx = allIndices[[offset + 1 ;; offset + n]];
        expr2 = expr2 * constructor[myIdx];
        offset += n;
        , {bi, blockList}
      ];
      expr = expr * expr2;
    ],
    expr = MakeProduct[blockList]
  ];

  Print["\n--- ", name, " (", nIdx, " indices) ---"];
  result = Catch[Quiet[AllContractions[expr], {AllContractions::nocontr}]];

  If[!ListQ[result],
    count = 0;
    Print["  0 contractions"],
    count = Length[result];
    Print["  ", count, " contractions"];
    Do[Print["    ", k, ". ", result[[k]] // InputForm], {k, Min[count, 8]}];
    If[count > 8, Print["    ... (", count - 8, " more)"]];
  ];

  AppendTo[allResults, {name, count}];
];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  PART 1: QUADRATIC (all pairs of 2 blocks)"];
Print["========================================"];
(* ============================================================ *)

(* Generate all unordered pairs (with repetition) from 1..nBlocks *)
quadraticPairs = Flatten[Table[{i, j}, {i, 1, nBlocks}, {j, i, nBlocks}], 1];

Print["Total quadratic pairs to check: ", Length[quadraticPairs]];

(* Parity-even *)
Do[EnumerateIfEven[pair], {pair, quadraticPairs}];

(* Parity-odd (with epsilon) *)
Print["\n-- Parity-odd quadratic (with epsilon): --"];
Do[EnumerateIfEven[pair, "eps"], {pair, quadraticPairs}];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  PART 2: CUBIC (all triples of 3 blocks)"];
Print["  (relevant for B0 background linearization)"];
Print["========================================"];
(* ============================================================ *)

(* Generate all unordered triples (with repetition) *)
cubicTriples = Flatten[
  Table[{i, j, k}, {i, 1, nBlocks}, {j, i, nBlocks}, {k, j, nBlocks}], 2];

Print["Total cubic triples to check: ", Length[cubicTriples]];

(* Parity-even cubics *)
Do[EnumerateIfEven[triple], {triple, cubicTriples}];

(* Parity-odd cubics (with epsilon) — only if total indices <= 18 *)
Print["\n-- Parity-odd cubic (with epsilon): --"];
Do[EnumerateIfEven[triple, "eps"], {triple, cubicTriples}];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  PART 3: CHERN-SIMONS (A-potential)"];
Print["========================================"];
(* ============================================================ *)

(* eps x TorsionVector x A x F *)
Print["\n--- CS: eps Tvec A F ---"];
Block[{r},
  r = Catch[Quiet[AllContractions[
    eps[-a,-b,-c,-d] * Tvec[-e] * AA[-f] * FF[-g,-h]]]];
  If[ListQ[r],
    Print["  ", Length[r], " contractions"];
    Do[Print["    ", k, ". ", r[[k]] // InputForm], {k, Length[r]}];
    AppendTo[allResults, {"CS: eps Tvec A F", Length[r]}],
    AppendTo[allResults, {"CS: eps Tvec A F", 0}]
  ]
];

(* eps x FullTorsion x A x F *)
Print["\n--- CS: eps Tor A F ---"];
Block[{r},
  r = Catch[Quiet[AllContractions[
    eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * AA[-h] * FF[-i,-j]]]];
  If[ListQ[r],
    Print["  ", Length[r], " contractions"];
    Do[Print["    ", k, ". ", r[[k]] // InputForm], {k, Min[Length[r], 8]}];
    AppendTo[allResults, {"CS: eps Tor A F", Length[r]}],
    AppendTo[allResults, {"CS: eps Tor A F", 0}]
  ]
];

(* ============================================================ *)
Print["\n\n========================================"];
Print["  MASTER SUMMARY"];
Print["========================================"];
(* ============================================================ *)

(* Filter out zero-count results for the summary *)
Print["\n  Non-zero sectors:"];
nonZero = Select[allResults, #[[2]] > 0 &];
Do[Print["  ", r[[1]], ": ", r[[2]]], {r, nonZero}];

Print["\n  Total non-zero sectors: ", Length[nonZero]];
Print["  Total independent contractions: ", Total[nonZero[[All, 2]]]];
Print["\n  Zero sectors (verified vanishing): ", Length[allResults] - Length[nonZero]];

(* === Mixing Classification === *)
Print["\n\n========================================"];
Print["  MIXING CLASSIFICATION (B0!=0, Tbar=0)"];
Print["========================================"];

ClassifySector[{name_String, count_Integer}] := Module[{mixing, hasR, hasT, hasF, nF},
  hasR = StringContainsQ[name, "R"];
  hasT = StringContainsQ[name, "T"] && !StringMatchQ[name, "*eps*"] || StringContainsQ[name, "DT"] || StringContainsQ[name, "Tvec"];
  (* More careful: T but not just in "eps" *)
  hasT = StringContainsQ[name, RegularExpression["\\bD?T\\b"]];
  hasF = StringContainsQ[name, RegularExpression["\\bD?F\\b"]];
  nF = StringCount[name, RegularExpression["\\bD?F\\b"]];

  mixing = Which[
    count == 0, "vanishing",
    hasR && hasT && hasF, "h-t mixing via B0",
    hasR && hasF && !hasT, If[nF >= 2, "h-f mixing (Gertsenshtein)", "h-f mixing"],
    hasT && hasF && !hasR, If[nF >= 2, "t-f mixing via B0", "t-f mixing"],
    hasR && hasT && !hasF, "h-t mixing (B0-independent)",
    hasR && !hasT && !hasF, "curvature self-interaction",
    hasT && !hasR && !hasF, "torsion self/kinetic",
    hasF && !hasR && !hasT, "EM self/kinetic",
    True, "unclassified"
  ];
  {name, count, mixing}
];

classified = ClassifySector /@ allResults;
nonZeroClassified = Select[classified, #[[2]] > 0 &];

Do[
  Print["  ", r[[1]], " (", r[[2]], "): ", r[[3]]],
  {r, nonZeroClassified}
];

(* === Export === *)
exportData = Association @@ ((#[[1]] -> <|"count" -> #[[2]], "mixing" -> #[[3]]|>) & /@ classified);
Export[
  FileNameJoin[{DirectoryName[$InputFileName], "enumeration_classified.json"}],
  exportData
];

countsOnly = Association @@ ((#[[1]] -> #[[2]]) & /@ allResults);
Export[
  FileNameJoin[{DirectoryName[$InputFileName], "enumeration_results.json"}],
  countsOnly
];

Print["\nExported to enumeration_results.json and enumeration_classified.json"];
Print["=== DONE ==="];
