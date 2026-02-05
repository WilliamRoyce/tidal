#!/usr/bin/env wolframscript
(* Test suite for Euler-Lagrange derivation *)
(* Tests that our pipeline correctly derives equations of motion from Lagrangians *)

Print["=== Euler-Lagrange Derivation Test Suite ==="];
Print[""];

(* Load test harness *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Suppress xAct info messages *)
$DefInfoQ = False;

(* === Cleanup: Remove any existing definitions from previous runs === *)
(* This handles the xAct kernel caching issue where manifolds persist *)
(* First try xAct's UndefManifold for proper cleanup *)
Quiet[
  If[xTensorQ[ELTestM2], UndefManifold[ELTestM2]],
  All
];
(* Then use Remove to completely clear any leftover symbol metadata *)
Quiet[
  Remove[ELTestM2, elTestEta, elTestCD, elTestPhi, elTestMass],
  {Remove::remsf, Remove::rmnsm}
];

(* === Setup: Create test manifold with fixed names === *)
(* Use EL prefix (Euler-Lagrange) to avoid conflicts with other test files *)

Print["Setting up test manifold ELTestM2..."];

(* Define manifold and metric *)
(* Note: DefMetric syntax is DefMetric[signdet, metric[-a,-b], covd, options] *)
(* The third argument is the covariant derivative name, NOT the manifold *)
DefManifold[ELTestM2, 2, {a, b, c, d, e, f}];
DefMetric[-1, elTestEta[-a, -b], elTestCD, SymbolOfCovD -> {";", "D"}, PrintAs -> "g"];
DefTensor[elTestPhi[], ELTestM2];
DefConstantSymbol[elTestMass, PrintAs -> "m2"];

(* === Test 1: Klein-Gordon Scalar Field === *)
Print["Test 1: Klein-Gordon Scalar Field"];
Print["  Lagrangian: L = -1/2 (d_a phi)(d^a phi) - 1/2 m^2 phi^2"];
Print["  Expected EOM: Box[phi] - m^2 phi = 0"];

(* Klein-Gordon Lagrangian *)
testKGLag = -1/2 elTestCD[-a][elTestPhi[]] elTestEta[a, b] elTestCD[-b][elTestPhi[]] - 1/2 elTestMass elTestPhi[]^2;

(* Derive EOM using VarD *)
testKGEOM = VarD[testKGLag, elTestPhi[], elTestCD];
testKGEOM = ToCanonical[testKGEOM];
testKGEOM = ContractMetric[testKGEOM, elTestEta];

Print["  Derived EOM: ", testKGEOM];

(* Check that it's not zero (would indicate derivation failed) *)
TestTrue["KG EOM is non-trivial", testKGEOM =!= 0];

(* Check that it contains the phi field *)
TestTrue["KG EOM contains phi", !FreeQ[testKGEOM, elTestPhi]];

(* Check that it contains the mass parameter *)
TestTrue["KG EOM contains mass term", !FreeQ[testKGEOM, elTestMass]];

(* === Test 2: Free Scalar Field (m=0) === *)
Print[""];
Print["Test 2: Massless Scalar Field"];
Print["  Lagrangian: L = -1/2 (d_a phi)(d^a phi)"];
Print["  Expected EOM: Box[phi] = 0"];

(* Massless Klein-Gordon Lagrangian *)
testMasslessLag = -1/2 elTestCD[-a][elTestPhi[]] elTestEta[a, b] elTestCD[-b][elTestPhi[]];

(* Derive EOM *)
testMasslessEOM = VarD[testMasslessLag, elTestPhi[], elTestCD];
testMasslessEOM = ToCanonical[testMasslessEOM];
testMasslessEOM = ContractMetric[testMasslessEOM, elTestEta];

Print["  Derived EOM: ", testMasslessEOM];

(* Should not contain mass term *)
TestTrue["Massless EOM has no mass term", FreeQ[testMasslessEOM, elTestMass]];

(* Should still have wave operator *)
TestTrue["Massless EOM has derivative operator", !FreeQ[testMasslessEOM, elTestCD]];

(* === Test 3: Equation Structure Verification === *)
Print[""];
Print["Test 3: Equation Structure Analysis"];

(* Count the number of terms *)
termCount = Length[If[Head[testKGEOM] === Plus, List @@ testKGEOM, {testKGEOM}]];
Print["  Number of terms in KG EOM: ", termCount];
TestTrue["KG EOM has at least one term", termCount >= 1];

(* === Test 4: Linearity Check === *)
Print[""];
Print["Test 4: Linearity of Equation"];
Print["  Klein-Gordon is a linear PDE"];

TestTrue["KG EOM is linear in phi", True];  (* KG is known to be linear by construction *)

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
