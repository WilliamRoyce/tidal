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
  "..", "tidal", "wolfram"
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

(* === Test 7: IdentifySpatialCrossDerivative === *)
testCross = Derivative[0, 1, 1][phi0][t, x, y];
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

(* === Test Group: ClassifySpatialProfile (Issue #79) === *)
Print[""];
Print["Test Group: ClassifySpatialProfile"];

csp1 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{1, 0}];
TestCase["ClassifySpatialProfile {1,0} -> gradient_x", "gradient_x", csp1];

csp2 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{0, 1}];
TestCase["ClassifySpatialProfile {0,1} -> gradient_y", "gradient_y", csp2];

csp3 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{2, 0}];
TestCase["ClassifySpatialProfile {2,0} -> laplacian_x", "laplacian_x", csp3];

csp4 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{0, 2}];
TestCase["ClassifySpatialProfile {0,2} -> laplacian_y", "laplacian_y", csp4];

csp5 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{1, 1}];
TestCase["ClassifySpatialProfile {1,1} -> cross_derivative_xy",
  "cross_derivative_xy", csp5];

csp6 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{1, 0, 1}];
TestCase["ClassifySpatialProfile {1,0,1} -> cross_derivative_xz",
  "cross_derivative_xz", csp6];

csp7 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{0, 1, 1}];
TestCase["ClassifySpatialProfile {0,1,1} -> cross_derivative_yz",
  "cross_derivative_yz", csp7];

csp8 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{3, 0}];
TestCase["ClassifySpatialProfile {3,0} -> derivative_3_x",
  "derivative_3_x", csp8];

csp9 = TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{2, 1}];
TestCase["ClassifySpatialProfile {2,1} -> derivative_2x_1y",
  "derivative_2x_1y", csp9];

(* All-zeros should throw *)
cspThrown = Catch[TorsionGertsenshtein`ExportJSON`ClassifySpatialProfile[{0, 0}]];
TestTrue["ClassifySpatialProfile {0,0} throws",
  StringQ[cspThrown] && StringContainsQ[cspThrown, "ClassifySpatialProfile"]];

(* === Test Group: ExtractSpatialOperatorFromMixed (Issue #79) === *)
Print[""];
Print["Test Group: ExtractSpatialOperatorFromMixed"];

(* Define a test symbol for derivative expressions *)
If[!ValueQ[testf], testf[t_, x_] := 0; testf[t_, x_, y_] := 0; testf[t_, x_, y_, z_] := 0];

(* Backward compat: single spatial direction *)
esom1 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 1, 0][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d_x -> gradient_x",
  "gradient_x", esom1];

esom2 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 0, 1][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d_y -> gradient_y",
  "gradient_y", esom2];

(* Bug fixes: multi-spatial mixed derivatives *)
esom3 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 1, 1][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d_x d_y -> cross_derivative_xy",
  "cross_derivative_xy", esom3];

esom4 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 2, 0][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d2_x -> laplacian_x",
  "laplacian_x", esom4];

esom5 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 0, 2][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d2_y -> laplacian_y",
  "laplacian_y", esom5];

(* Generic multi-axis *)
esom6 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 2, 1][testf][t, x, y]];
TestCase["ExtractSpatialOperatorFromMixed d_t d2_x d_y -> derivative_2x_1y",
  "derivative_2x_1y", esom6];

(* 1+1D backward compat *)
esom7 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 1][testf][t, x]];
TestCase["ExtractSpatialOperatorFromMixed 1+1D d_t d_x -> gradient_x",
  "gradient_x", esom7];

(* 3+1D *)
esom8 = TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[1, 0, 0, 1][testf][t, x, y, z]];
TestCase["ExtractSpatialOperatorFromMixed 3+1D d_t d_z -> gradient_z",
  "gradient_z", esom8];

(* Time order != 1 should throw *)
esomThrown = Catch[TorsionGertsenshtein`ExportJSON`ExtractSpatialOperatorFromMixed[
  Derivative[2, 1, 0][testf][t, x, y]]];
TestTrue["ExtractSpatialOperatorFromMixed throws on time order 2",
  StringQ[esomThrown] && StringContainsQ[esomThrown, "Time derivative order"]];

(* === End-to-end: ClassifyOperatorType for mixed derivatives (Issue #79) === *)
Print[""];
Print["Test Group: ClassifyOperatorType mixed derivative fixes"];

{opMixedXY, isMixedXY} = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[
  Derivative[1, 1, 1][testf][t, x, y]];
TestCase["ClassifyOperatorType d_t d_x d_y -> cross_derivative_xy",
  "cross_derivative_xy", opMixedXY];
TestTrue["ClassifyOperatorType d_t d_x d_y sets momentum flag", isMixedXY];

{opMixedLapX, isMixedLapX} = TorsionGertsenshtein`ExportJSON`ClassifyOperatorType[
  Derivative[1, 2, 0][testf][t, x, y]];
TestCase["ClassifyOperatorType d_t d2_x -> laplacian_x",
  "laplacian_x", opMixedLapX];
TestTrue["ClassifyOperatorType d_t d2_x sets momentum flag", isMixedLapX];

(* === Test Group: ExtractDerivativeProfile === *)
Print[""];
Print["Test Group: ExtractDerivativeProfile"];

(* 1+1D *)
edp1 = TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
  Derivative[2, 0][phi0][t, x]];
TestCase["ExtractDerivativeProfile 1+1D d2_t", {2, 0}, edp1];

(* 2+1D *)
edp2 = TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
  Derivative[0, 1, 1][phi0][t, x, y]];
TestCase["ExtractDerivativeProfile 2+1D d_x d_y", {0, 1, 1}, edp2];

(* 3+1D *)
edp3 = TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
  Derivative[1, 0, 0, 1][testf][t, x, y, z]];
TestCase["ExtractDerivativeProfile 3+1D d_t d_z", {1, 0, 0, 1}, edp3];

(* No derivative -> {} *)
edp4 = TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
  phi0[t, x]];
TestCase["ExtractDerivativeProfile no derivative returns {}", {}, edp4];

(* With coefficient *)
edp5 = TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
  3.0 * Derivative[1, 2][phi0][t, x]];
TestCase["ExtractDerivativeProfile with coefficient", {1, 2}, edp5];

(* === Test Group: 3+1D wrapper coverage === *)
Print[""];
Print["Test Group: 3+1D wrapper coverage"];

(* 3+1D gradient_z *)
gradZ3D = TorsionGertsenshtein`ExportJSON`IdentifyGradientDirection[
  Derivative[0, 0, 0, 1][testf][t, x, y, z]];
TestCase["IdentifyGradientDirection 3+1D d_z", "gradient_z", gradZ3D];

(* 3+1D laplacian_z *)
lapZ3D = TorsionGertsenshtein`ExportJSON`IdentifyDirectionalLaplacian[
  Derivative[0, 0, 0, 2][testf][t, x, y, z]];
TestCase["IdentifyDirectionalLaplacian 3+1D d2_z", "laplacian_z", lapZ3D];

(* 3+1D cross_derivative_xz *)
crossXZ3D = TorsionGertsenshtein`ExportJSON`IdentifySpatialCrossDerivative[
  Derivative[0, 1, 0, 1][testf][t, x, y, z]];
TestCase["IdentifySpatialCrossDerivative 3+1D d_x d_z", "cross_derivative_xz", crossXZ3D];

(* 3+1D IsMixedTimeSpaceDerivative *)
TestTrue["IsMixedTimeSpaceDerivative 3+1D d_t d_x",
  TorsionGertsenshtein`ExportJSON`IsMixedTimeSpaceDerivative[
    Derivative[1, 1, 0, 0][testf][t, x, y, z]]];

(* 3+1D IsFirstOrderTimeDerivative *)
TestTrue["IsFirstOrderTimeDerivative 3+1D",
  TorsionGertsenshtein`ExportJSON`IsFirstOrderTimeDerivative[
    Derivative[1, 0, 0, 0][testf][t, x, y, z]]];

(* === Test Group: CountDerivativeOrder === *)
Print[""];
Print["Test Group: CountDerivativeOrder"];

TestCase["CountDerivativeOrder: second spatial derivative",
  2, TorsionGertsenshtein`ExportJSON`CountDerivativeOrder[
    Derivative[0, 2][testf][t, x]]];

TestCase["CountDerivativeOrder: mixed 3-arg sum of orders",
  3, TorsionGertsenshtein`ExportJSON`CountDerivativeOrder[
    Derivative[1, 1, 1][testf][t, x, y]]];

TestCase["CountDerivativeOrder: no derivative returns 0",
  0, TorsionGertsenshtein`ExportJSON`CountDerivativeOrder[
    testf[t, x]]];

(* === Test Group: BuildGenericOperatorName === *)
Print[""];
Print["Test Group: BuildGenericOperatorName"];

TestCase["BuildGenericOperatorName: single-axis third order",
  "derivative_3_x",
  TorsionGertsenshtein`ExportJSON`BuildGenericOperatorName[{3, 0}]];

TestCase["BuildGenericOperatorName: multi-axis mixed",
  "derivative_2x_1y",
  TorsionGertsenshtein`ExportJSON`BuildGenericOperatorName[{2, 1}]];

(* === Test Group: ExtractSpatialDerivativeProfile === *)
Print[""];
Print["Test Group: ExtractSpatialDerivativeProfile"];

TestCase["ExtractSpatialDerivativeProfile: strips time slot",
  {2, 0},
  TorsionGertsenshtein`ExportJSON`ExtractSpatialDerivativeProfile[
    Derivative[1, 2, 0][testf][t, x, y]]];

TestTrue["ExtractSpatialDerivativeProfile: throws on no derivative",
  StringQ[Catch[
    TorsionGertsenshtein`ExportJSON`ExtractSpatialDerivativeProfile[
      testf[t, x]]]]];

(* === Test Group: ExtractDerivativeProfile edge cases === *)
Print[""];
Print["Test Group: ExtractDerivativeProfile edge cases"];

TestCase["ExtractDerivativeProfile: unapplied Derivative",
  {1, 2},
  TorsionGertsenshtein`ExportJSON`ExtractDerivativeProfile[
    Derivative[1, 2][testf]]];

(* === Test Group: FieldToVelocityName — v_field_name convention === *)
Print[""];
Print["Test Group: FieldToVelocityName v_field_name"];

(* Single-base system *)
TestCase["FieldToVelocityName: single base A_0",
  "v_A_0",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["A_0", {"A_0", "A_1", "A_2"}]];

TestCase["FieldToVelocityName: single base A_2",
  "v_A_2",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["A_2", {"A_0", "A_1", "A_2"}]];

(* Multi-base system: v_ prefix + field name regardless of global position *)
TestCase["FieldToVelocityName: multi-base A_0",
  "v_A_0",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["A_0", {"phi_0", "A_0", "A_1", "A_2"}]];

TestCase["FieldToVelocityName: multi-base A_1",
  "v_A_1",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["A_1", {"phi_0", "A_0", "A_1", "A_2"}]];

(* Two vector bases *)
TestCase["FieldToVelocityName: two vectors B_0",
  "v_B_0",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["B_0", {"A_0", "A_1", "A_2", "B_0", "B_1", "B_2"}]];

TestCase["FieldToVelocityName: two vectors B_1",
  "v_B_1",
  TorsionGertsenshtein`ExportJSON`FieldToVelocityName["B_1", {"A_0", "A_1", "A_2", "B_0", "B_1", "B_2"}]];

(* Error case: field not in list *)
TestTrue["FieldToMomentumName: throws on missing field",
  StringQ[Catch[
    TorsionGertsenshtein`ExportJSON`FieldToMomentumName["C_0", {"A_0", "A_1"}]]]];

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
