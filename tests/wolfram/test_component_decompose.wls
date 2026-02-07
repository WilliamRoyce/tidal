#!/usr/bin/env wolframscript
(* Test suite for ComponentDecompose.wl *)
(* Tests higher-rank tensor decomposition (rank 3+) *)

Print["=== ComponentDecompose Test Suite ==="];
Print[""];

(* Load test harness *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

(* Get pipeline path for module loading *)
pipelinePath = FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "..", "torsion_gertsenshtein", "wolfram"
}];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Suppress xAct info messages *)
$DefInfoQ = False;

(* Load pipeline modules *)
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];

(* === Setup: Define test manifolds and tensors === *)

(* Use CDT prefix to avoid conflicts with other test files *)
Quiet[
  If[xTensorQ[CDTM2], UndefManifold[CDTM2]];
  If[xTensorQ[CDTM3], UndefManifold[CDTM3]],
  All
];

(* 2D manifold for basic tests *)
DefManifold[CDTM2, 2, {a, b, c, d, e, f}];
DefMetric[-1, cdtEta2[-a, -b], cdtCD2];
DefChart[cdtCart2, CDTM2, {0, 1}, {t[], x[]}];

(* 3D manifold for higher-dim tests *)
DefManifold[CDTM3, 3, {g, h, i, j, k, l}];
DefMetric[-1, cdtEta3[-g, -h], cdtCD3];
DefChart[cdtCart3, CDTM3, {0, 1, 2}, {t[], x[], y[]}];

(* === Test 1: EnumerateComponentTuples for various ranks === *)
Print["Test 1: EnumerateComponentTuples"];

(* Rank-0 scalar *)
DefTensor[cdtPhi[], CDTM2];
tuples0 = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtPhi, 2];
TestCase["Rank-0 scalar: single empty tuple", {{}}, tuples0];

(* Rank-1 vector *)
DefTensor[cdtV[-a], CDTM2];
tuples1 = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtV, 2];
TestCase["Rank-1 vector in 2D: 2 components", {{0}, {1}}, tuples1];

(* Rank-2 symmetric *)
DefTensor[cdtH[-a, -b], CDTM2, Symmetric[{-a, -b}]];
tuples2s = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtH, 2];
TestCase["Rank-2 symmetric in 2D: 3 components (upper triangle)",
  {{0, 0}, {0, 1}, {1, 1}}, tuples2s];

(* Rank-2 antisymmetric *)
DefTensor[cdtF[-a, -b], CDTM2, Antisymmetric[{-a, -b}]];
tuples2as = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtF, 2];
TestTrue["Rank-2 antisymmetric in 2D: has non-trivial symmetry",
  SymmetryGroupOfTensor[cdtF] =!= StrongGenSet[{}, GenSet[]]];

(* Rank-3 non-symmetric in 2D: should be 2^3 = 8 tuples *)
DefTensor[cdtT[-a, -b, -c], CDTM2];
tuples3 = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtT, 2];
TestCase["Rank-3 non-symmetric in 2D: 8 components",
  8, Length[tuples3]];
TestCase["Rank-3 first tuple is {0,0,0}",
  {0, 0, 0}, tuples3[[1]]];
TestCase["Rank-3 last tuple is {1,1,1}",
  {1, 1, 1}, tuples3[[-1]]];

(* Rank-3 in 3D: should be 3^3 = 27 tuples *)
DefTensor[cdtT3[-g, -h, -i], CDTM3];
tuples3d3 = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtT3, 3];
TestCase["Rank-3 non-symmetric in 3D: 27 components",
  27, Length[tuples3d3]];

(* Rank-4 in 2D: should be 2^4 = 16 tuples *)
DefTensor[cdtR[-a, -b, -c, -d], CDTM2];
tuples4 = TorsionGertsenshtein`ComponentDecompose`EnumerateComponentTuples[cdtR, 2];
TestCase["Rank-4 non-symmetric in 2D: 16 components",
  16, Length[tuples4]];

(* === Test 2: ReplaceHigherRankFieldComponents === *)
Print[""];
Print["Test 2: ReplaceHigherRankFieldComponents"];

(* Build a test expression with rank-3 basis-indexed components *)
(* cdtT[{0,-cdtCart2},{0,-cdtCart2},{0,-cdtCart2}] should -> cdtT0[t[],x[]] *)
coordSyms2 = {t, x};

testExpr3 = cdtT[{0, -cdtCart2}, {0, -cdtCart2}, {0, -cdtCart2}] +
             cdtT[{0, -cdtCart2}, {0, -cdtCart2}, {1, -cdtCart2}] +
             cdtT[{1, -cdtCart2}, {1, -cdtCart2}, {1, -cdtCart2}];

replaced3 = TorsionGertsenshtein`ComponentDecompose`ReplaceHigherRankFieldComponents[
  testExpr3, cdtT, cdtCart2, coordSyms2, 2
];

(* After replacement, expression should contain cdtT0, cdtT1, cdtT7 *)
(* Tuple ordering: {0,0,0}->T0, {0,0,1}->T1, ..., {1,1,1}->T7 *)
TestTrue["Rank-3 replacement: contains cdtT0",
  !FreeQ[replaced3, cdtT0]];
TestTrue["Rank-3 replacement: contains cdtT1",
  !FreeQ[replaced3, cdtT1]];
TestTrue["Rank-3 replacement: contains cdtT7",
  !FreeQ[replaced3, cdtT7]];
TestTrue["Rank-3 replacement: no basis-indexed cdtT remaining",
  FreeQ[replaced3, cdtT[__]]];

(* Test that ReplaceTensorFieldComponents dispatches to the new function *)
replaced3b = TorsionGertsenshtein`ComponentDecompose`ReplaceTensorFieldComponents[
  testExpr3, cdtT, cdtCart2, coordSyms2, 2
];
TestTrue["ReplaceTensorFieldComponents dispatches rank-3 correctly",
  replaced3 === replaced3b];

(* === Test 3: Rank-3 via unified ReplaceTensorFieldComponents === *)
Print[""];
Print["Test 3: ReplaceTensorFieldComponents for rank-4"];

(* Rank-4 replacement *)
testExpr4 = cdtR[{0, -cdtCart2}, {0, -cdtCart2}, {0, -cdtCart2}, {0, -cdtCart2}] +
             cdtR[{1, -cdtCart2}, {1, -cdtCart2}, {1, -cdtCart2}, {1, -cdtCart2}];

replaced4 = TorsionGertsenshtein`ComponentDecompose`ReplaceTensorFieldComponents[
  testExpr4, cdtR, cdtCart2, coordSyms2, 2
];

(* {0,0,0,0}->R0, {1,1,1,1}->R15 *)
TestTrue["Rank-4 replacement: contains cdtR0",
  !FreeQ[replaced4, cdtR0]];
TestTrue["Rank-4 replacement: contains cdtR15",
  !FreeQ[replaced4, cdtR15]];
TestTrue["Rank-4 replacement: no basis-indexed cdtR remaining",
  FreeQ[replaced4, cdtR[__]]];

(* === Test 4: Mixed index signs (chart / -chart) === *)
Print[""];
Print["Test 4: Mixed index signs for rank-3"];

(* Upper indices should also be replaced *)
testExprMixed = cdtT[{0, cdtCart2}, {0, -cdtCart2}, {1, -cdtCart2}];
replacedMixed = TorsionGertsenshtein`ComponentDecompose`ReplaceHigherRankFieldComponents[
  testExprMixed, cdtT, cdtCart2, coordSyms2, 2
];
TestTrue["Rank-3 mixed index signs: replaced correctly",
  FreeQ[replacedMixed, cdtT[__]]];
TestTrue["Rank-3 mixed index signs: maps to cdtT1",
  !FreeQ[replacedMixed, cdtT1]];

(* === Test 5: End-to-end DecomposeToComponents for rank-3 === *)
Print[""];
Print["Test 5: End-to-end DecomposeToComponents for rank-3"];

(* Define a simple rank-3 tensor EOM: T_{abc} = 0 with a source *)
(* We use a trivial equation: cdtCD2[-a][cdtT[-b,-c,-d]] = 0 *)
(* This tests the full pipeline: decompose -> extract -> replace *)

(* Simple identity EOM for rank-3: just the tensor itself *)
eomT = cdtT[-a, -b, -c];

result = Catch[
  TorsionGertsenshtein`ComponentDecompose`DecomposeToComponents[
    eomT, cdtT[-a, -b, -c], cdtCart2
  ]
];

(* Should produce 8 component equations for non-symmetric rank-3 in 2D *)
If[ListQ[result],
  TestCase["DecomposeToComponents rank-3: produces 8 components",
    8, Length[result]];
  TestTrue["DecomposeToComponents rank-3: first component index is 0",
    result[[1, 1]] === 0];
  TestTrue["DecomposeToComponents rank-3: last component index is 7",
    result[[-1, 1]] === 7];
  (* Each component equation should be a scalar (no basis indices left) *)
  TestTrue["DecomposeToComponents rank-3: first component is scalar (cdtT0)",
    !FreeQ[result[[1, 2]], cdtT0]],
  (* If Catch returned a string, it means an error was thrown *)
  TestTrue["DecomposeToComponents rank-3: should not throw",
    False]
];

(* === Test 6: Rank-3 component equations are non-trivial === *)
Print[""];
Print["Test 6: Rank-3 component equations are non-trivial"];

(* Verify that different component tuples produce different scalar functions *)
(* The {0,0,0} component should give cdtT0, the {0,0,1} should give cdtT1, etc. *)
TestTrue["Rank-3 components are distinct scalars",
  Length[Union[result[[All, 2]]]] === 8];

(* Verify each component equation contains its corresponding scalar *)
TestTrue["Component 0 contains cdtT0",
  !FreeQ[result[[1, 2]], cdtT0]];
TestTrue["Component 7 contains cdtT7",
  !FreeQ[result[[8, 2]], cdtT7]];

(* === Test 7: Rank-3 in 3D (27 components) === *)
Print[""];
Print["Test 7: Rank-3 in 3D"];

eomT3 = cdtT3[-g, -h, -i];

resultT3 = Catch[
  TorsionGertsenshtein`ComponentDecompose`DecomposeToComponents[
    eomT3, cdtT3[-g, -h, -i], cdtCart3
  ]
];

If[ListQ[resultT3],
  TestCase["Rank-3 in 3D: produces 27 components",
    27, Length[resultT3]],
  TestTrue["Rank-3 in 3D: should not throw", False]
];

(* === Done === *)
ExitWithTestStatus[];
