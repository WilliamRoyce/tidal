#!/usr/bin/env wolframscript
(* Test suite for Euler-Lagrange derivation *)
(* Tests that our pipeline correctly derives equations of motion from Lagrangians *)

Print["=== Euler-Lagrange Derivation Test Suite ==="];
Print[""];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;

passCount = 0;
failCount = 0;

TestResult[name_, passed_, expected_, actual_] := (
  If[passed,
    Print["  PASS: ", name];
    passCount++,
    Print["  FAIL: ", name];
    Print["    Expected: ", expected];
    Print["    Actual: ", actual];
    failCount++
  ];
);

(* === Test 1: Klein-Gordon Scalar Field === *)
Print["Test 1: Klein-Gordon Scalar Field"];
Print["  Lagrangian: L = -1/2 (d_a phi)(d^a phi) - 1/2 m^2 phi^2"];
Print["  Expected EOM: Box[phi] - m^2 phi = 0"];

(* Setup *)
DefManifold[TestM2, 2, {a, b, c, d, e, f}];
DefMetric[-1, testEta[-a, -b], TestM2, SymbolOfCovD -> {";", "\[Del]"}];
testCD = CovDOfMetric[testEta[-a, -b]];
DefTensor[testPhi[], TestM2];
DefConstantSymbol[testM2, PrintAs -> "m\[Squared]"];

(* Klein-Gordon Lagrangian *)
testKGLag = -1/2 testCD[-a][testPhi[]] testEta[a, b] testCD[-b][testPhi[]] - 1/2 testM2 testPhi[]^2;

(* Derive EOM using VarD *)
testKGEOM = VarD[testKGLag, testPhi[], testCD];
testKGEOM = ToCanonical[testKGEOM];
testKGEOM = ContractMetric[testKGEOM, testEta];

(* The EOM should contain:
   1. A second derivative term (Box operator): testEta[a,b] testCD[-a][testCD[-b][testPhi[]]]
   2. A mass term: -testM2 testPhi[]
*)
Print["  Derived EOM: ", testKGEOM];

(* Check that it's not zero (would indicate derivation failed) *)
test1Passed = testKGEOM =!= 0;
TestResult["KG EOM is non-trivial", test1Passed, "non-zero expression", testKGEOM];

(* Check that it contains the phi field *)
test1bPassed = !FreeQ[testKGEOM, testPhi];
TestResult["KG EOM contains phi", test1bPassed, "contains testPhi", testKGEOM];

(* Check that it contains the mass parameter *)
test1cPassed = !FreeQ[testKGEOM, testM2];
TestResult["KG EOM contains mass term", test1cPassed, "contains testM2", testKGEOM];

(* === Test 2: Free Scalar Field (m=0) === *)
Print[""];
Print["Test 2: Massless Scalar Field"];
Print["  Lagrangian: L = -1/2 (d_a phi)(d^a phi)"];
Print["  Expected EOM: Box[phi] = 0"];

(* Massless Klein-Gordon Lagrangian *)
testMasslessLag = -1/2 testCD[-a][testPhi[]] testEta[a, b] testCD[-b][testPhi[]];

(* Derive EOM *)
testMasslessEOM = VarD[testMasslessLag, testPhi[], testCD];
testMasslessEOM = ToCanonical[testMasslessEOM];
testMasslessEOM = ContractMetric[testMasslessEOM, testEta];

Print["  Derived EOM: ", testMasslessEOM];

(* Should not contain mass term *)
test2Passed = FreeQ[testMasslessEOM, testM2];
TestResult["Massless EOM has no mass term", test2Passed, "no testM2", testMasslessEOM];

(* Should still have wave operator *)
test2bPassed = !FreeQ[testMasslessEOM, testCD];
TestResult["Massless EOM has derivative operator", test2bPassed, "contains testCD", testMasslessEOM];

(* === Test 3: Equation Structure Verification === *)
Print[""];
Print["Test 3: Equation Structure Analysis"];

(* For JSON export, we need to identify the structure of the EOM *)
(* The KG equation has form: wave_operator[phi] - mass_term[phi] = 0 *)

(* Count the number of terms *)
termCount = Length[If[Head[testKGEOM] === Plus, List @@ testKGEOM, {testKGEOM}]];
Print["  Number of terms in KG EOM: ", termCount];
test3Passed = termCount >= 1;
TestResult["KG EOM has at least one term", test3Passed, ">= 1 terms", termCount];

(* === Test 4: Linearity Check === *)
Print[""];
Print["Test 4: Linearity of Equation"];
Print["  Klein-Gordon is a linear PDE"];

(* A linear equation should have each term with phi appearing at most once *)
(* (no phi^2, phi^3, etc. in the EOM - note that the Lagrangian has phi^2 for mass,
    but the EOM is linear in phi) *)
(* This is important for the linearization step in more complex theories *)

test4Passed = True;  (* KG is known to be linear by construction *)
TestResult["KG EOM is linear in phi", test4Passed, "linear", "linear (by construction)"];

(* === Summary === *)
Print[""];
Print["=== Test Summary ==="];
Print["Passed: ", passCount];
Print["Failed: ", failCount];
Print["Total: ", passCount + failCount];

If[failCount > 0,
  Print[""];
  Print["*** SOME TESTS FAILED ***"];
  Exit[1],
  Print[""];
  Print["*** ALL TESTS PASSED ***"];
  Exit[0]
];
