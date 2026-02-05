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
  "..", "torsion_gertsenshtein", "wolfram"
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
TestTrue["ValidateDimension[5] is False", !TorsionGertsenshtein`CommonUtilities`ValidateDimension[5]];

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

(* === Test 4: ExtractNumericCoefficient === *)
Print[""];
Print["Test 4: ExtractNumericCoefficient"];

(* Create test tensor *)
DefTensor[cuTestPhi[], CUTestM2];

(* Test scalar field coefficient extraction *)
testTerm1 = 3.5 * cuTestPhi[];
coeff1 = TorsionGertsenshtein`CommonUtilities`ExtractNumericCoefficient[testTerm1, "cuTestPhi"];
TestApproxEqual["ExtractNumericCoefficient for 3.5*phi[]", 3.5, coeff1];

testTerm2 = -2 * cuTestPhi[];
coeff2 = TorsionGertsenshtein`CommonUtilities`ExtractNumericCoefficient[testTerm2, "cuTestPhi"];
TestApproxEqual["ExtractNumericCoefficient for -2*phi[]", -2.0, coeff2];

(* Test derivative term coefficient extraction *)
testTerm3 = 4.5 * Derivative[1, 0][phi0][t[], x[]];
coeff3 = TorsionGertsenshtein`CommonUtilities`ExtractNumericCoefficient[testTerm3, "phi0"];
TestApproxEqual["ExtractNumericCoefficient for derivative term", 4.5, coeff3];

(* === Test 5: ExtractCoefficientWithSymbolic === *)
Print[""];
Print["Test 5: ExtractCoefficientWithSymbolic"];

testTerm4 = 5 * cuTestPhi[];
coeffResult4 = TorsionGertsenshtein`CommonUtilities`ExtractCoefficientWithSymbolic[testTerm4, "cuTestPhi"];
TestApproxEqual["ExtractCoefficientWithSymbolic numeric part", 5.0, coeffResult4["numeric"]];
TestCase["ExtractCoefficientWithSymbolic symbolic part is None", None, coeffResult4["symbolic"]];

(* === Test 6: xAct Introspection Helpers === *)
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

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
