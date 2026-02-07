#!/usr/bin/env wolframscript
(* Test suite for Linearize.wl module *)
(* Tests xPert-based linearization: metric perturbation, tensor linearization, *)
(* gauge conditions, and utility functions (IsLinear, SelectLinearTerms). *)

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
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "torsion_gertsenshtein", "wolfram"}];
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
(* === Test Group 3: Backward-Compatible Alias === *)
(* ================================================================ *)

Print[""];
Print["Test Group 3: LinearizedEinsteinEquation (backward-compatible alias)"];

(* Test 3.1: Alias produces same result as LinearizeTensorExpression *)
linEinsteinAlias = LinearizedEinsteinEquation[Einstein[linTestCD][-a, -b]];

TestTrue["LinearizedEinsteinEquation gives same result as LinearizeTensorExpression",
  linEinsteinAlias === linEinstein
];

(* Test 3.2: Alias produces a well-formed result (may be zero in 2D) *)
TestTrue["Alias result is well-formed (not $Failed)",
  linEinsteinAlias =!= $Failed
];

(* ================================================================ *)
(* === Test Group 4: IsLinear and SelectLinearTerms === *)
(* ================================================================ *)

Print[""];
Print["Test Group 4: IsLinear and SelectLinearTerms"];

(* Use plain symbols for polynomial tests *)
Block[{x, aa, bb, cc},

  (* Test 4.1: Linear expression *)
  TestTrue["IsLinear[aa*x + bb, x] is True",
    IsLinear[aa*x + bb, x]
  ];

  (* Test 4.2: Quadratic expression is not linear *)
  TestTrue["IsLinear[x^2, x] is False",
    !IsLinear[x^2, x]
  ];

  (* Test 4.3: Pure constant is linear (degree 0 <= 1) *)
  TestTrue["IsLinear[cc, x] is True (constant is degree 0)",
    IsLinear[cc, x]
  ];

  (* Test 4.4: Mixed polynomial with linear and quadratic terms *)
  TestTrue["IsLinear[aa*x + bb*x^2, x] is False",
    !IsLinear[aa*x + bb*x^2, x]
  ];

  (* Test 4.5: SelectLinearTerms extracts linear part *)
  selectedTerms = SelectLinearTerms[aa*x + bb*x^2 + cc, x];
  Print["  SelectLinearTerms[aa*x + bb*x^2 + cc, x] = ", selectedTerms];
  TestTrue["SelectLinearTerms extracts aa*x from aa*x + bb*x^2 + cc",
    Simplify[selectedTerms - aa*x] === 0
  ];

  (* Test 4.6: SelectLinearTerms on purely linear expression *)
  selectedTerms2 = SelectLinearTerms[3*x + 5, x];
  TestTrue["SelectLinearTerms[3*x + 5, x] gives 3*x",
    Simplify[selectedTerms2 - 3*x] === 0
  ];

  (* Test 4.7: SelectLinearTerms on purely quadratic expression *)
  selectedTerms3 = SelectLinearTerms[x^2 + x^3, x];
  TestTrue["SelectLinearTerms[x^2 + x^3, x] gives 0",
    selectedTerms3 === 0
  ];
];

(* ================================================================ *)
(* === Test Group 5: HasExactlyOneFieldFactor === *)
(* ================================================================ *)

Print[""];
Print["Test Group 5: HasExactlyOneFieldFactor"];

Block[{x, aa},

  TestTrue["HasExactlyOneFieldFactor[aa*x, x] is True",
    HasExactlyOneFieldFactor[aa*x, x]
  ];

  (* HasExactlyOneFieldFactor is designed for product-form PDE terms
     where fields appear as f[args]. Use two explicit field applications. *)
  TestTrue["HasExactlyOneFieldFactor detects two field-related factors",
    !HasExactlyOneFieldFactor[x[t] * x[s], x]
  ];

  TestTrue["HasExactlyOneFieldFactor[aa, x] is False (no field)",
    !HasExactlyOneFieldFactor[aa, x]
  ];
];

(* ================================================================ *)
(* === Test Group 6: TTGaugeConditions === *)
(* ================================================================ *)

Print[""];
Print["Test Group 6: TTGaugeConditions"];

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
