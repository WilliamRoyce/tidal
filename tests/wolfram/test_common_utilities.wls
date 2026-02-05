#!/usr/bin/env wolframscript
(* Test suite for CommonUtilities.wl *)
(* Tests CD conversion, epsilon evaluation, coefficient extraction, and xAct introspection *)

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

(* Load CommonUtilities *)
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];

(* === Test 1: GetChartDimension === *)
Print["Test 1: GetChartDimension"];

(* Setup 2D manifold *)
If[!xTensorQ[TestM2],
  DefManifold[TestM2, 2, {a, b, c}];
  DefMetric[-1, testEta2[-a, -b], TestM2];
  DefChart[testCart2, TestM2, {0, 1}, {t[], x[]}]
];

dim2 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[testCart2];
TestCase["GetChartDimension for 2D chart", 2, dim2];

(* Setup 3D manifold *)
If[!xTensorQ[TestM3],
  DefManifold[TestM3, 3, {d, e, f}];
  DefMetric[-1, testEta3[-d, -e], TestM3];
  DefChart[testCart3, TestM3, {0, 1, 2}, {t[], x[], y[]}]
];

dim3 = TorsionGertsenshtein`CommonUtilities`GetChartDimension[testCart3];
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

(* Create test term with numeric coefficient *)
DefTensor[testPhi2[], TestM2];
testTerm1 = 3.5 * testPhi2[];
coeff1 = TorsionGertsenshtein`CommonUtilities`ExtractNumericCoefficient[testTerm1, "testPhi2"];
TestApproxEqual["ExtractNumericCoefficient for 3.5*phi", 3.5, coeff1];

testTerm2 = -2 * testPhi2[];
coeff2 = TorsionGertsenshtein`CommonUtilities`ExtractNumericCoefficient[testTerm2, "testPhi2"];
TestApproxEqual["ExtractNumericCoefficient for -2*phi", -2.0, coeff2];

(* === Test 5: ExtractCoefficientWithSymbolic === *)
Print[""];
Print["Test 5: ExtractCoefficientWithSymbolic"];

testTerm3 = 5 * testPhi2[];
coeffResult3 = TorsionGertsenshtein`CommonUtilities`ExtractCoefficientWithSymbolic[testTerm3, "testPhi2"];
TestApproxEqual["ExtractCoefficientWithSymbolic numeric part", 5.0, coeffResult3["numeric"]];
TestCase["ExtractCoefficientWithSymbolic symbolic part is None", None, coeffResult3["symbolic"]];

(* === Test 6: xAct Introspection Helpers === *)
Print[""];
Print["Test 6: xAct Introspection Helpers"];

(* Test IsChristoffelSymbol *)
christoffel = ChristoffelOfMetric[testEta2[-a, -b]];
TestTrue["IsChristoffelSymbol detects Christoffel symbol",
  TorsionGertsenshtein`CommonUtilities`IsChristoffelSymbol[christoffel]];
TestTrue["IsChristoffelSymbol returns False for non-Christoffel",
  !TorsionGertsenshtein`CommonUtilities`IsChristoffelSymbol[testPhi2]];

(* Test IsCovDOperator *)
testCD2 = CovDOfMetric[testEta2[-a, -b]];
TestTrue["IsCovDOperator detects covariant derivative",
  TorsionGertsenshtein`CommonUtilities`IsCovDOperator[testCD2[-a][testPhi2[]]]];

(* Test IsEpsilonTensor *)
epsilonTensor = epsilontestEta2[-a, -b];
TestTrue["IsEpsilonTensor detects epsilon tensor",
  TorsionGertsenshtein`CommonUtilities`IsEpsilonTensor[epsilonTensor]];

(* === Test 7: ConvertCDToDerivatives === *)
Print[""];
Print["Test 7: ConvertCDToDerivatives"];

(* Create a covariant derivative term *)
DefTensor[testScalar[], TestM2];
cdTerm = testCD2[{0, -testCart2}][testScalar[]];

(* Convert to Derivative form *)
converted = TorsionGertsenshtein`CommonUtilities`ConvertCDToDerivatives[cdTerm, testCart2];

(* The result should contain Derivative, not CD *)
TestTrue["ConvertCDToDerivatives produces Derivative form",
  !FreeQ[converted, Derivative]];
TestTrue["ConvertCDToDerivatives removes CD form",
  FreeQ[converted, testCD2]];

(* === Test 8: GenerateCDRules === *)
Print[""];
Print["Test 8: GenerateCDRules"];

rules2D = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[2, testCart2];
TestTrue["GenerateCDRules returns non-empty list for 2D",
  Length[rules2D] > 0];

rules3D = TorsionGertsenshtein`CommonUtilities`GenerateCDRules[3, testCart3];
TestTrue["GenerateCDRules returns more rules for 3D than 2D",
  Length[rules3D] > Length[rules2D]];

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
