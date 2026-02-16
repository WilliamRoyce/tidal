#!/usr/bin/env wolframscript
(* Test suite for Linearize.wl module *)
(* Tests xPert-based linearization: metric perturbation, tensor linearization, *)
(* gauge conditions, and TT gauge condition tools. *)

Print["=== Linearize Module Test Suite ==="];
Print[""];

(* Load test harness *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;
<< xAct`xPert`;

(* Suppress xAct info messages *)
$DefInfoQ = False;

(* Load the Linearize module *)
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "tidal", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];

(* === Cleanup: Remove any existing definitions from previous runs === *)
Quiet[
  If[xTensorQ[LinTestM2], UndefManifold[LinTestM2]],
  All
];
Quiet[
  Remove[LinTestM2, linTestEta, linTestCD, linTestH, linTestEpsilon],
  {Remove::remsf, Remove::rmnsm}
];

(* === Setup: Create test manifold with Minkowski metric === *)
(* Use LinTest prefix to avoid conflicts with other test files *)

Print["Setting up test manifold LinTestM2 (2D Minkowski)..."];
Print[""];

DefManifold[LinTestM2, 2, {a, b, c, d, e, f}];
DefMetric[-1, linTestEta[-a, -b], linTestCD, SymbolOfCovD -> {";", "D"}, PrintAs -> "g"];

(* ================================================================ *)
(* === Test Group 1: SetupMetricPerturbation === *)
(* ================================================================ *)

Print["Test Group 1: SetupMetricPerturbation"];

(* Test 1.1: SetupMetricPerturbation returns a perturbation tensor with LI[1] *)
linTestHPert = SetupMetricPerturbation[linTestEta, linTestH, linTestEpsilon];
TestTrue["SetupMetricPerturbation returns expression with LI[1]",
  !FreeQ[linTestHPert, LI[1]]
];

(* Test 1.2: The returned tensor has LI[1] and two index slots *)
TestTrue["Perturbation tensor has LI[1] with two additional slots",
  Length[linTestHPert] >= 3 && linTestHPert[[1]] === LI[1]
];

(* Test 1.3: Idempotent — calling twice does not crash *)
(* xPert may issue ValidateSymbol warnings but should not fail *)
linTestHPert2 = Quiet[SetupMetricPerturbation[linTestEta, linTestH, linTestEpsilon], All];
TestTrue["SetupMetricPerturbation is idempotent (second call returns result)",
  linTestHPert2 =!= Null
];

(* Test 1.4: Second call returns something with LI[1] *)
TestTrue["Second call returns perturbation structure with LI[1]",
  !FreeQ[linTestHPert2, LI[1]]
];

(* ================================================================ *)
(* === Test Group 2: LinearizeTensorExpression === *)
(* ================================================================ *)

Print[""];
Print["Test Group 2: LinearizeTensorExpression"];

(* Test 2.1: Linearize Einstein tensor — computes without error *)
(* NOTE: In 2D, the Einstein tensor G_ab = R_ab - (1/2)g_ab R vanishes
   identically, so the linearized version is also zero. This is CORRECT physics.
   We test that the function runs without error and returns a result. *)
linEinstein = LinearizeTensorExpression[Einstein[linTestCD][-a, -b]];
Print["  Linearized Einstein tensor (2D, expected 0): ", linEinstein];

TestTrue["Linearized Einstein tensor computes without error",
  linEinstein =!= $Failed
];

(* Test 2.2: In 2D, Einstein tensor should be zero (G_ab = 0 identically) *)
TestTrue["In 2D, linearized Einstein tensor is zero (correct physics)",
  linEinstein === 0
];

(* Test 2.3: Linearize Ricci tensor — also works (not Einstein-specific) *)
linRicci = LinearizeTensorExpression[RicciCD[-a, -b] /. RicciCD -> Ricci[linTestCD]];
(* Alternative: use the Ricci tensor function directly *)
linRicci2 = LinearizeTensorExpression[Ricci[linTestCD][-a, -b]];
Print["  Linearized Ricci tensor: ", linRicci2];

TestTrue["Linearized Ricci tensor is non-trivial",
  linRicci2 =!= 0
];

TestTrue["Linearized Ricci tensor contains perturbation h",
  !FreeQ[linRicci2, linTestH]
];

(* Test 2.4: Einstein and Ricci linearizations are distinct *)
(* Einstein = Ricci - 1/2 g R, so their linearizations should differ *)
TestTrue["Einstein and Ricci linearizations are different expressions",
  linEinstein =!= linRicci2
];

(* ================================================================ *)
(* === Test Group 3: TTGaugeConditions === *)
(* ================================================================ *)

Print[""];
Print["Test Group 3: TTGaugeConditions"];

ttResult = TTGaugeConditions[linTestHPert, linTestEta];

(* Test 6.1: Returns an Association *)
TestTrue["TTGaugeConditions returns an Association",
  AssociationQ[ttResult]
];

(* Test 6.2: Contains expected keys *)
TestTrue["TTGaugeConditions has 'gauge' key",
  KeyExistsQ[ttResult, "gauge"]
];

TestTrue["TTGaugeConditions has 'tracelessness' key",
  KeyExistsQ[ttResult, "tracelessness"]
];

TestTrue["TTGaugeConditions has 'transversality' key",
  KeyExistsQ[ttResult, "transversality"]
];

TestTrue["TTGaugeConditions has 'physical_dof' key",
  KeyExistsQ[ttResult, "physical_dof"]
];

TestTrue["TTGaugeConditions has 'polarizations' key",
  KeyExistsQ[ttResult, "polarizations"]
];

(* Test 6.3: Gauge name is correct *)
TestCase["TTGaugeConditions gauge name",
  "transverse_traceless",
  ttResult["gauge"]
];

(* Test 6.4: Physical DOF count *)
TestCase["TTGaugeConditions physical_dof is 2",
  2,
  ttResult["physical_dof"]
];

(* === Exit === *)
Print[""];
ExitWithTestStatus[];
