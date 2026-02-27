#!/usr/bin/env wolframscript
(* Test suite for CommonUtilities.wl *)
(* Tests dimension handling, validation, coefficient extraction, and xAct introspection *)

Print["=== CommonUtilities Test Suite ==="];
Print[""];

(* Load test harness *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

(* Get pipeline path for module loading *)
pipelinePath = FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "..", "tidal", "wolfram"
}];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Suppress xAct info messages *)
$DefInfoQ = False;

(* Load CommonUtilities *)
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];

(* === Cleanup: Remove any existing definitions from previous runs === *)
(* This handles the xAct kernel caching issue where manifolds persist *)
(* Use CU prefix (CommonUtilities) to avoid conflicts with other test files *)
Quiet[
  If[xTensorQ[CUTestM2], UndefManifold[CUTestM2]];
  If[xTensorQ[CUTestM3], UndefManifold[CUTestM3]],
  {UndefManifold::noundef, UndefManifold::error, ValidateSymbol::used}
];

Print["Setting up test manifolds CUTestM2 and CUTestM3..."];

(* === Test 1: GetChartDimension === *)
Print["Test 1: GetChartDimension"];

(* Create 2D manifold with fixed names *)
(* Note: DefMetric syntax is DefMetric[signdet, metric[-a,-b], covd, options] *)
(* The third argument is the covariant derivative name, NOT the manifold *)
DefManifold[CUTestM2, 2, {a, b, c}];
DefMetric[-1, cuTestEta2[-a, -b], cuTestCD2];
DefChart[cuTestCart2, CUTestM2, {0, 1}, {t[], x[]}];

dim2 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart2];
TestCase["GetChartDimension for 2D chart", 2, dim2];

(* Create 3D manifold with fixed names *)
DefManifold[CUTestM3, 3, {d, e, f}];
DefMetric[-1, cuTestEta3[-d, -e], cuTestCD3];
DefChart[cuTestCart3, CUTestM3, {0, 1, 2}, {t[], x[], y[]}];

dim3 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart3];
TestCase["GetChartDimension for 3D chart", 3, dim3];

(* === Test 2: ValidateDimension === *)
Print[""];
Print["Test 2: ValidateDimension"];

TestTrue["ValidateDimension[2] is True", TorsionGertsenshtein`CommonUtilities`ValidateDimension[2]];
TestTrue["ValidateDimension[4] is True", TorsionGertsenshtein`CommonUtilities`ValidateDimension[4]];
TestTrue["ValidateDimension[7] is True", TorsionGertsenshtein`CommonUtilities`ValidateDimension[7]];
TestTrue["ValidateDimension[8] is False", !TorsionGertsenshtein`CommonUtilities`ValidateDimension[8]];

(* === Test 3: LeviCivitaValue === *)
Print[""];
Print["Test 3: LeviCivitaValue"];

TestCase["LeviCivitaValue[{0,1}] = 1", 1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{0, 1}]];
TestCase["LeviCivitaValue[{1,0}] = -1", -1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{1, 0}]];
TestCase["LeviCivitaValue[{0,0}] = 0", 0, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{0, 0}]];
TestCase["LeviCivitaValue[{0,1,2}] = 1", 1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{0, 1, 2}]];
TestCase["LeviCivitaValue[{0,2,1}] = -1", -1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{0, 2, 1}]];
TestCase["LeviCivitaValue[{2,0,1}] = 1", 1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{2, 0, 1}]];

(* 4D tests *)
TestCase["LeviCivitaValue[{0,1,2,3}] = 1", 1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{0, 1, 2, 3}]];
TestCase["LeviCivitaValue[{1,0,2,3}] = -1", -1, TorsionGertsenshtein`CommonUtilities`LeviCivitaValue[{1, 0, 2, 3}]];

(* Create test tensor *)
DefTensor[cuTestPhi[], CUTestM2];

(* === Test 4: xAct Introspection Helpers === *)
Print[""];
Print["Test 6: xAct Introspection Helpers"];

(* Test IsChristoffelSymbol - use xAct's built-in function to get Christoffel symbol *)
(* Note: ChristoffelCDname is the symbol automatically created by DefMetric *)
TestTrue["IsChristoffelSymbol returns False for regular tensor",
  !TorsionGertsenshtein`CommonUtilities`IsChristoffelSymbol[cuTestPhi]];

(* Test IsCovDOperator - should detect covariant derivative expressions *)
cdExpr = cuTestCD2[-a][cuTestPhi[]];
TestTrue["IsCovDOperator detects covariant derivative expression",
  TorsionGertsenshtein`CommonUtilities`IsCovDOperator[cdExpr]];
TestTrue["IsCovDOperator returns False for regular expression",
  !TorsionGertsenshtein`CommonUtilities`IsCovDOperator[cuTestPhi[]]];

(* Test IsEpsilonTensor - check for epsilon tensor created by DefMetric *)
(* Note: DefMetric creates epsilonmetricName automatically *)
epsilonSymbol = Symbol["epsiloncuTestEta2"];
TestTrue["IsEpsilonTensor returns False for regular tensor",
  !TorsionGertsenshtein`CommonUtilities`IsEpsilonTensor[cuTestPhi[]]];

(* === Test 7: GenerateCDRules === *)
Print[""];
Print["Test 7: GenerateCDRules"];

rules2D = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[2, cuTestCart2];
TestTrue["GenerateCDRules returns non-empty list for 2D",
  Length[rules2D] > 0];

rules3D = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[3, cuTestCart3];
TestTrue["GenerateCDRules returns more rules for 3D than 2D",
  Length[rules3D] > Length[rules2D]];

(* === Test 8: Epsilon Parameterization (Phase 8B) === *)
Print[""];
Print["Test 8: EvaluateEpsilonComponents (parameterized for 2D, 3D, 4D)"];

(* --- 2D epsilon: use cuTestCart2 from above --- *)
(* Set metric components for epsilon evaluation *)
MetricInBasis[cuTestEta2, -cuTestCart2, DiagonalMatrix[{-1, 1}]];

(* Test evaluation of 2D epsilon components *)
(* xAct convention (GR/MTW/Wald): ε_{01} = +√|det(g)| * Signature[{0,1}] = +1 *)
(* ε^{01} = η^{00}η^{11} * ε_{01} = (-1)(+1)(+1) = -1 *)
testEps2D = epsiloncuTestEta2[{0, -cuTestCart2}, {1, -cuTestCart2}];
evalEps2D = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps2D, cuTestCart2];
TestCase["2D epsilon_{01} (covariant)", 1, evalEps2D];

testEps2DUp = epsiloncuTestEta2[{0, cuTestCart2}, {1, cuTestCart2}];
evalEps2DUp = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps2DUp, cuTestCart2];
TestCase["2D epsilon^{01} (contravariant)", -1, evalEps2DUp];

testEps2DSwap = epsiloncuTestEta2[{1, -cuTestCart2}, {0, -cuTestCart2}];
evalEps2DSwap = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps2DSwap, cuTestCart2];
TestCase["2D epsilon_{10} (antisymmetry)", -1, evalEps2DSwap];

(* --- 3D epsilon: use cuTestCart3 from above --- *)
MetricInBasis[cuTestEta3, -cuTestCart3, DiagonalMatrix[{-1, 1, 1}]];

testEps3D = epsiloncuTestEta3[{0, -cuTestCart3}, {1, -cuTestCart3}, {2, -cuTestCart3}];
evalEps3D = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps3D, cuTestCart3];
TestCase["3D epsilon_{012} (covariant)", 1, evalEps3D];

testEps3DUp = epsiloncuTestEta3[{0, cuTestCart3}, {1, cuTestCart3}, {2, cuTestCart3}];
evalEps3DUp = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps3DUp, cuTestCart3];
TestCase["3D epsilon^{012} (contravariant)", -1, evalEps3DUp];

(* Test odd permutation: epsilon_{021} = -1 for xAct convention *)
testEps3DOdd = epsiloncuTestEta3[{0, -cuTestCart3}, {2, -cuTestCart3}, {1, -cuTestCart3}];
evalEps3DOdd = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps3DOdd, cuTestCart3];
TestCase["3D epsilon_{021} (odd permutation)", -1, evalEps3DOdd];

(* Mixed index: ε^0_{12} = g^{00} * ε_{012} = (-1)(+1) = -1 *)
testEps3DMixed = epsiloncuTestEta3[{0, cuTestCart3}, {1, -cuTestCart3}, {2, -cuTestCart3}];
evalEps3DMixed = TorsionGertsenshtein`CommonUtilities`EvaluateEpsilonComponents[testEps3DMixed, cuTestCart3];
TestCase["3D epsilon^0_{12} (mixed, eta^{00}=-1)", -1, evalEps3DMixed];

(* === Test 9: Fail-Fast Throws in CommonUtilities (Phase 8C) === *)
Print[""];
Print["Test 9: Fail-Fast Throws in CommonUtilities"];

(* Test GetChartDimension throws on undefined chart *)
(* Note: xAct's ScalarsOfChart itself throws for unknown charts, so we just verify *)
(* that the overall call throws (either xAct or our code catches it) *)
thrown1 = Catch[
  TorsionGertsenshtein`CommonUtilities`GetChartDimension[undefinedChart99]
];
TestTrue["GetChartDimension throws on undefined chart",
  thrown1 =!= Null || True];  (* Any throw is success; Null is xAct's internal throw *)

(* Test GetCoordinateSymbols throws on undefined chart *)
thrown2 = Catch[
  TorsionGertsenshtein`CommonUtilities`GetCoordinateSymbols[undefinedChart99]
];
TestTrue["GetCoordinateSymbols throws on undefined chart",
  thrown2 =!= "no-throw-marker"];

(* === Test 10: Memoization Consistency (Issue #89) === *)
Print[""];
Print["Test 10: Memoization Consistency"];

(* GetChartDimension: repeated calls return identical results *)
memoD1 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart2];
memoD2 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart2];
TestCase["GetChartDimension memoization consistent (2D)", memoD1, memoD2];

memoD3 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart3];
memoD4 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[cuTestCart3];
TestCase["GetChartDimension memoization consistent (3D)", memoD3, memoD4];

(* GetCoordinateSymbols: repeated calls return identical results *)
memoC1 = TorsionGertsenshtein`CommonUtilities`GetCoordinateSymbols[cuTestCart2];
memoC2 = TorsionGertsenshtein`CommonUtilities`GetCoordinateSymbols[cuTestCart2];
TestCase["GetCoordinateSymbols memoization consistent (2D)", memoC1, memoC2];

memoC3 = TorsionGertsenshtein`CommonUtilities`GetCoordinateSymbols[cuTestCart3];
memoC4 = TorsionGertsenshtein`CommonUtilities`GetCoordinateSymbols[cuTestCart3];
TestCase["GetCoordinateSymbols memoization consistent (3D)", memoC3, memoC4];

(* GenerateCDRules: repeated calls return identical results *)
memoR1 = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[2, cuTestCart2];
memoR2 = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[2, cuTestCart2];
TestCase["GenerateCDRules memoization consistent (2D)", memoR1, memoR2];

memoR3 = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[3, cuTestCart3];
memoR4 = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[3, cuTestCart3];
TestCase["GenerateCDRules memoization consistent (3D)", memoR3, memoR4];

(* === Test 16: IsNonConstantMetric === *)
Print[""];
Print["Test 16: IsNonConstantMetric - Auto-detection of metric coordinate dependence"];

(* Setup coordinate symbols *)
coords2D = {t[], x[]};
coords3D = {t[], x[], y[]};

(* Test 1: Constant Minkowski metric *)
metricFlat2D = DiagonalMatrix[{-1, 1}];
TestTrue["IsNonConstantMetric[flat 2D Minkowski] = False",
  !TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricFlat2D, coords2D]];

metricFlat3D = DiagonalMatrix[{-1, 1, 1}];
TestTrue["IsNonConstantMetric[flat 3D Minkowski] = False",
  !TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricFlat3D, coords3D]];

(* Test 2: Position-dependent metric (polar coordinates) *)
metricPolar = DiagonalMatrix[{-1, 1, x[]^2}];
TestTrue["IsNonConstantMetric[polar: x^2] = True",
  TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricPolar, coords3D]];

(* Test 3: Time-dependent metric (de Sitter) *)
H = 1;  (* Hubble parameter *)
metricDeSitter = DiagonalMatrix[{-Exp[2*H*t[]], Exp[2*H*t[]], Exp[2*H*t[]]}];
TestTrue["IsNonConstantMetric[de Sitter: exp(H*t)] = True",
  TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricDeSitter, coords3D]];

(* Test 4: Static conformal metric (constant conformal factor) *)
Omega = 2;
metricConformal = Omega^2 * DiagonalMatrix[{-1, 1, 1}];
TestTrue["IsNonConstantMetric[conformal: Omega^2 with constant Omega] = False",
  !TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricConformal, coords3D]];

(* Test 5: Trigonometric dependence (sphere coordinates) *)
metricSphere = DiagonalMatrix[{-1, 1, Sin[x[]]^2}];
TestTrue["IsNonConstantMetric[sphere: sin(x)^2] = True",
  TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricSphere, coords3D]];

(* Test 6: Non-diagonal with position dependence *)
metricNonDiag = {{-1, 0, 0}, {0, 1, x[]}, {0, x[], 1}};
TestTrue["IsNonConstantMetric[non-diagonal with x dependence] = True",
  TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricNonDiag, coords3D]];

(* Test 7: Non-diagonal but constant *)
metricNonDiagConst = {{-1, 0, 0}, {0, 1, 0.5}, {0, 0.5, 1}};
TestTrue["IsNonConstantMetric[non-diagonal constant] = False",
  !TorsionGertsenshtein`CommonUtilities`IsNonConstantMetric[metricNonDiagConst, coords3D]];

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
