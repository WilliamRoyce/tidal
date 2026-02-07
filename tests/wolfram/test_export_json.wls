#!/usr/bin/env wolframscript
(* Test suite for ExportJSON.wl *)
(* Tests operator detection, term identification, JSON structure building *)

Print["=== ExportJSON Test Suite ==="];
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

(* Load pipeline modules *)
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

(* === Test 1: ParseFieldName === *)
Print["Test 1: ParseFieldName"];

(* Standard format: A_0, phi_1 *)
parsed1 = TorsionGertsenshtein`ExportJSON`ParseFieldName["A_0"];
TestCase["ParseFieldName A_0 base", "A", parsed1["base"]];
TestCase["ParseFieldName A_0 index", 0, parsed1["index"]];
TestCase["ParseFieldName A_0 format", "standard", parsed1["format"]];

parsed2 = TorsionGertsenshtein`ExportJSON`ParseFieldName["phi_3"];
TestCase["ParseFieldName phi_3 base", "phi", parsed2["base"]];
TestCase["ParseFieldName phi_3 index", 3, parsed2["index"]];

(* Tensor format: stress_xy_0 *)
parsed3 = TorsionGertsenshtein`ExportJSON`ParseFieldName["stress_xy_0"];
TestCase["ParseFieldName stress_xy_0 base", "stress_xy", parsed3["base"]];
TestCase["ParseFieldName stress_xy_0 format", "tensor", parsed3["format"]];

(* Compact format: phi0 *)
parsed4 = TorsionGertsenshtein`ExportJSON`ParseFieldName["phi0"];
TestCase["ParseFieldName phi0 base", "phi", parsed4["base"]];
TestCase["ParseFieldName phi0 index", 0, parsed4["index"]];
TestCase["ParseFieldName phi0 format", "compact", parsed4["format"]];

(* Simple format: phi *)
parsed5 = TorsionGertsenshtein`ExportJSON`ParseFieldName["phi"];
TestCase["ParseFieldName phi format", "simple", parsed5["format"]];

(* === Test 2: GetTimeDerivativeOrder === *)
Print[""];
Print["Test 2: GetTimeDerivativeOrder"];

(* Create test expressions with various time derivative orders *)
(* 1+1D: Derivative[dt, dx][f][t, x] - first slot is time *)
testDeriv0 = phi0[t, x];
testDeriv1 = Derivative[1, 0][phi0][t, x];
testDeriv2 = Derivative[2, 0][phi0][t, x];
testDeriv3 = Derivative[0, 2][phi0][t, x];  (* Pure spatial *)

order0 = TorsionGertsenshtein`ExportJSON`GetTimeDerivativeOrder[testDeriv0];
order1 = TorsionGertsenshtein`ExportJSON`GetTimeDerivativeOrder[testDeriv1];
order2 = TorsionGertsenshtein`ExportJSON`GetTimeDerivativeOrder[testDeriv2];
order3 = TorsionGertsenshtein`ExportJSON`GetTimeDerivativeOrder[testDeriv3];

TestCase["GetTimeDerivativeOrder for f[t,x]", 0, order0];
TestCase["GetTimeDerivativeOrder for d_t f", 1, order1];
TestCase["GetTimeDerivativeOrder for d2_t f", 2, order2];
TestCase["GetTimeDerivativeOrder for d2_x f (pure spatial)", 0, order3];

(* === Test 3: ContainsTimeDerivative === *)
Print[""];
Print["Test 3: ContainsTimeDerivative"];

TestTrue["ContainsTimeDerivative with minOrder=1",
  TorsionGertsenshtein`ExportJSON`ContainsTimeDerivative[testDeriv1, 1]];
TestTrue["ContainsTimeDerivative with minOrder=2",
  TorsionGertsenshtein`ExportJSON`ContainsTimeDerivative[testDeriv2, 2]];
TestTrue["ContainsTimeDerivative returns False for spatial",
  !TorsionGertsenshtein`ExportJSON`ContainsTimeDerivative[testDeriv3, 1]];

(* 2+1D case *)
testDeriv2D = Derivative[2, 0, 0][phi0][t, x, y];
TestTrue["ContainsTimeDerivative for 2+1D",
  TorsionGertsenshtein`ExportJSON`ContainsTimeDerivative[testDeriv2D, 2]];

(* === Test 4: IsMixedTimeSpaceDerivative === *)
Print[""];
Print["Test 4: IsMixedTimeSpaceDerivative"];

(* Mixed: d_t d_x *)
testMixed = Derivative[1, 1][phi0][t, x];
TestTrue["IsMixedTimeSpaceDerivative for d_t d_x",
  TorsionGertsenshtein`ExportJSON`IsMixedTimeSpaceDerivative[testMixed]];

(* Pure time *)
TestTrue["IsMixedTimeSpaceDerivative returns False for pure time",
  !TorsionGertsenshtein`ExportJSON`IsMixedTimeSpaceDerivative[testDeriv2]];

(* Pure spatial *)
TestTrue["IsMixedTimeSpaceDerivative returns False for pure spatial",
  !TorsionGertsenshtein`ExportJSON`IsMixedTimeSpaceDerivative[testDeriv3]];

(* 2+1D mixed *)
testMixed2D = Derivative[1, 1, 0][phi0][t, x, y];
TestTrue["IsMixedTimeSpaceDerivative for 2+1D d_t d_x",
  TorsionGertsenshtein`ExportJSON`IsMixedTimeSpaceDerivative[testMixed2D]];

(* === Test 5: IdentifyGradientDirection === *)
Print[""];
Print["Test 5: IdentifyGradientDirection"];

(* 1+1D gradient_x *)
testGradX1D = Derivative[0, 1][phi0][t, x];
gradDir1 = TorsionGertsenshtein`ExportJSON`IdentifyGradientDirection[testGradX1D];
TestCase["IdentifyGradientDirection for 1+1D d_x", "gradient_x", gradDir1];

(* 2+1D gradient_x *)
testGradX2D = Derivative[0, 1, 0][phi0][t, x, y];
gradDir2 = TorsionGertsenshtein`ExportJSON`IdentifyGradientDirection[testGradX2D];
TestCase["IdentifyGradientDirection for 2+1D d_x", "gradient_x", gradDir2];

(* 2+1D gradient_y *)
testGradY2D = Derivative[0, 0, 1][phi0][t, x, y];
gradDir3 = TorsionGertsenshtein`ExportJSON`IdentifyGradientDirection[testGradY2D];
TestCase["IdentifyGradientDirection for 2+1D d_y", "gradient_y", gradDir3];

(* === Test 6: IsSpatialCrossDerivative === *)
Print[""];
Print["Test 6: IsSpatialCrossDerivative"];

(* Cross derivative d_x d_y *)
testCross = Derivative[0, 1, 1][phi0][t, x, y];
TestTrue["IsSpatialCrossDerivative for d_x d_y",
  TorsionGertsenshtein`ExportJSON`IsSpatialCrossDerivative[testCross]];

(* Not cross derivative *)
TestTrue["IsSpatialCrossDerivative returns False for d2_x",
  !TorsionGertsenshtein`ExportJSON`IsSpatialCrossDerivative[testGradX2D]];

(* === Test 7: IdentifySpatialCrossDerivative === *)
Print[""];
Print["Test 7: IdentifySpatialCrossDerivative"];

crossResult = TorsionGertsenshtein`ExportJSON`IdentifySpatialCrossDerivative[testCross];
TestCase["IdentifySpatialCrossDerivative returns cross_derivative_xy",
  "cross_derivative_xy", crossResult];

nonCrossResult = TorsionGertsenshtein`ExportJSON`IdentifySpatialCrossDerivative[testGradX2D];
TestCase["IdentifySpatialCrossDerivative returns False for non-cross",
  False, nonCrossResult];

(* === Test 8: IdentifyDirectionalLaplacian === *)
Print[""];
Print["Test 8: IdentifyDirectionalLaplacian"];

(* Pure d2_x *)
testLapX = Derivative[0, 2, 0][phi0][t, x, y];
lapXResult = TorsionGertsenshtein`ExportJSON`IdentifyDirectionalLaplacian[testLapX];
TestCase["IdentifyDirectionalLaplacian for d2_x", "laplacian_x", lapXResult];

(* Pure d2_y *)
testLapY = Derivative[0, 0, 2][phi0][t, x, y];
lapYResult = TorsionGertsenshtein`ExportJSON`IdentifyDirectionalLaplacian[testLapY];
TestCase["IdentifyDirectionalLaplacian for d2_y", "laplacian_y", lapYResult];

(* Not a directional laplacian (cross derivative) *)
crossLapResult = TorsionGertsenshtein`ExportJSON`IdentifyDirectionalLaplacian[testCross];
TestCase["IdentifyDirectionalLaplacian returns False for cross-derivative",
  False, crossLapResult];

(* === Test 9: ExtractFunctionHeads === *)
Print[""];
Print["Test 9: ExtractFunctionHeads"];

testExpr = 3.0 * Derivative[0, 2][phi0][t, x] + chi0[t, x];
heads = TorsionGertsenshtein`ExportJSON`ExtractFunctionHeads[testExpr];
TestTrue["ExtractFunctionHeads finds phi0",
  MemberQ[heads, "phi0"]];
TestTrue["ExtractFunctionHeads finds chi0",
  MemberQ[heads, "chi0"]];

(* === Test 10: ClassifyOperatorType === *)
Print[""];
Print["Test 10: ClassifyOperatorType"];

(* Identity (no derivative) *)
testIdentity = 5.0 * phi0[t, x];
{opIdentity, mixedIdentity} = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[testIdentity];
TestCase["ClassifyOperatorType for identity", "identity", opIdentity];
TestTrue["ClassifyOperatorType identity is not mixed", !mixedIdentity];

(* Laplacian *)
testLap = Derivative[0, 2][phi0][t, x];
{opLap, mixedLap} = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[testLap];
TestCase["ClassifyOperatorType for laplacian", "laplacian_x", opLap];

(* Mixed time-space *)
{opMixed, isMixed} = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[testMixed];
TestTrue["ClassifyOperatorType detects mixed", isMixed];

(* === Test 11: BuildLHSStructure === *)
Print[""];
Print["Test 11: BuildLHSStructure"];

lhs0 = TorsionGertsenshtein`ExportJSON`BuildLHSStructure["phi_0", 0];
TestCase["BuildLHSStructure order=0 expression", "phi_0", lhs0["expression"]];
TestCase["BuildLHSStructure order=0 time", 0, lhs0["order"]["time"]];

lhs1 = TorsionGertsenshtein`ExportJSON`BuildLHSStructure["phi_0", 1];
TestCase["BuildLHSStructure order=1 expression", "d_t(phi_0)", lhs1["expression"]];

lhs2 = TorsionGertsenshtein`ExportJSON`BuildLHSStructure["phi_0", 2];
TestCase["BuildLHSStructure order=2 expression", "d2_t(phi_0)", lhs2["expression"]];
TestCase["BuildLHSStructure order=2 time", 2, lhs2["order"]["time"]];

(* === Test 12: BuildTermResult === *)
Print[""];
Print["Test 12: BuildTermResult"];

(* Without symbolic coefficient *)
term1 = TorsionGertsenshtein`ExportJSON`BuildTermResult[3.5, "laplacian", "phi_0"];
TestApproxEqual["BuildTermResult coefficient", 3.5, term1["coefficient"]];
TestCase["BuildTermResult operator", "laplacian", term1["operator"]];
TestCase["BuildTermResult field", "phi_0", term1["field"]];
TestTrue["BuildTermResult no symbolic key without arg",
  !KeyExistsQ[term1, "coefficient_symbolic"]];

(* With symbolic coefficient *)
term2 = TorsionGertsenshtein`ExportJSON`BuildTermResult[1.0, "identity", "phi_0", "m2"];
TestTrue["BuildTermResult has symbolic key with arg",
  KeyExistsQ[term2, "coefficient_symbolic"]];
TestCase["BuildTermResult symbolic value", "m2", term2["coefficient_symbolic"]];

(* === Test 13: Fail-Fast Throws in ExportJSON (Phase 8C) === *)
Print[""];
Print["Test 13: Fail-Fast Throws in ExportJSON"];

(* Test ParseFieldName Throw on invalid name *)
thrown1 = Catch[
  TorsionGertsenshtein`ExportJSON`ParseFieldName["123invalid!!!"]
];
TestTrue["ParseFieldName throws on invalid name",
  StringQ[thrown1] && StringContainsQ[thrown1, "ParseFieldName"]];

(* Test ExtractLHSCoefficient Throw on non-derivative term *)
thrown2 = Catch[
  TorsionGertsenshtein`ExportJSON`ExtractLHSCoefficient[phi0[t, x]]
];
TestTrue["ExtractLHSCoefficient throws on non-derivative term",
  StringQ[thrown2] && StringContainsQ[thrown2, "ExtractLHSCoefficient"]];

(* Test ClassifyOperatorType handles order-3 derivative via generic fallback *)
testOrder3 = Derivative[0, 3][phi0][t, x];
result3 = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[testOrder3];
TestTrue["ClassifyOperatorType handles order-3 derivative generically",
  result3[[1]] === "derivative_3_x" && result3[[2]] === False];

(* Note: ParseMultiFieldRHS Throw is a guard for future code paths where *)
(* IdentifyMultiFieldTerm might return Nothing for some terms. Currently *)
(* IdentifyMultiFieldTerm always produces output (uses default field), so *)
(* the zero-terms condition is unreachable from tests without mocking. *)

(* === Test Group: ExtractMassCouplingFromEquations === *)
Print[""];
Print["Test Group: ExtractMassCouplingFromEquations"];

(* Test: Single field with mass term *)
(* Simulates: d2_t(phi_0) = laplacian(phi_0) - 1.0 * identity(phi_0) *)
singleFieldEqs = {
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> "phi_0"|>,
        <|"coefficient" -> -1.0, "operator" -> "identity", "field" -> "phi_0"|>
      }
    |>
  |>
};
singleResult = TorsionGertsenshtein`ExportJSON`ExtractMassCouplingFromEquations[
  singleFieldEqs, {"phi_0"}
];
TestCase["Single field mass_matrix",
  {{1.0}}, singleResult["mass_matrix"]];
TestCase["Single field coupling_matrix",
  {{0.0}}, singleResult["coupling_matrix"]];
TestTrue["Single field no symbolic mass",
  singleResult["mass_matrix_symbolic"] === Null];

(* Test: Multi-field with cross-coupling *)
(* Simulates coupled_scalars: phi has -1.0*phi - 0.5*chi, chi has -0.5*phi - 4.0*chi *)
multiFieldEqs = {
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> 1.0, "operator" -> "laplacian_x", "field" -> "phi_0"|>,
        <|"coefficient" -> -1.0, "operator" -> "identity", "field" -> "phi_0"|>,
        <|"coefficient" -> -0.5, "operator" -> "identity", "field" -> "chi_0"|>
      }
    |>
  |>,
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> 1.0, "operator" -> "laplacian_x", "field" -> "chi_0"|>,
        <|"coefficient" -> -4.0, "operator" -> "identity", "field" -> "chi_0"|>,
        <|"coefficient" -> -0.5, "operator" -> "identity", "field" -> "phi_0"|>
      }
    |>
  |>
};
multiResult = TorsionGertsenshtein`ExportJSON`ExtractMassCouplingFromEquations[
  multiFieldEqs, {"phi_0", "chi_0"}
];
TestCase["Multi-field mass_matrix diagonal",
  {{1.0, 0.0}, {0.0, 4.0}}, multiResult["mass_matrix"]];
TestCase["Multi-field coupling_matrix off-diagonal",
  {{0.0, 0.5}, {0.5, 0.0}}, multiResult["coupling_matrix"]];

(* Test: Massless field — no identity terms *)
masslessEqs = {
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> "A_0"|>
      }
    |>
  |>
};
masslessResult = TorsionGertsenshtein`ExportJSON`ExtractMassCouplingFromEquations[
  masslessEqs, {"A_0"}
];
TestCase["Massless field mass_matrix",
  {{0.0}}, masslessResult["mass_matrix"]];
TestCase["Massless field coupling_matrix",
  {{0.0}}, masslessResult["coupling_matrix"]];

(* Test: Symbolic coefficient *)
symbolicEqs = {
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> -1.0, "operator" -> "identity", "field" -> "phi_0",
          "coefficient_symbolic" -> "-procaMassSquared"|>
      }
    |>
  |>
};
symbolicResult = TorsionGertsenshtein`ExportJSON`ExtractMassCouplingFromEquations[
  symbolicEqs, {"phi_0"}
];
TestCase["Symbolic mass_matrix numeric value",
  {{1.0}}, symbolicResult["mass_matrix"]];
TestTrue["Symbolic mass_matrix_symbolic populated",
  symbolicResult["mass_matrix_symbolic"] =!= Null &&
  symbolicResult["mass_matrix_symbolic"][[1, 1]] === "-procaMassSquared"];

(* Test: Identity term referencing momentum field (pi_0) — should be ignored *)
momentumEqs = {
  <|
    "rhs" -> <|
      "terms" -> {
        <|"coefficient" -> 1.0, "operator" -> "gradient_x", "field" -> "pi_0"|>,
        <|"coefficient" -> -2.0, "operator" -> "identity", "field" -> "A_0"|>
      }
    |>
  |>
};
momentumResult = TorsionGertsenshtein`ExportJSON`ExtractMassCouplingFromEquations[
  momentumEqs, {"A_0"}
];
TestCase["Momentum ref ignored, own mass extracted",
  {{2.0}}, momentumResult["mass_matrix"]];

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
