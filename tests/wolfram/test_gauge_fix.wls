#!/usr/bin/env wolframscript
(* Test suite for GaugeFix.wl module *)
(* Tests AddGaugeFixingTerm (core primitive) and BuildLorenzGaugeTerm (preset). *)

Print["=== GaugeFix Module Test Suite ==="];
Print[""];

(* Load test harness *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

(* Load xAct packages *)
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Suppress xAct info messages *)
$DefInfoQ = False;

(* Load the GaugeFix module *)
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "tidal", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "GaugeFix.wl"}]];

(* === Cleanup: Remove any existing definitions from previous runs === *)
Quiet[
  If[xTensorQ[GFTestM2], UndefManifold[GFTestM2]],
  All
];
Quiet[
  Remove[GFTestM2, gfTestEta, gfTestCD, gfTestA],
  {Remove::remsf, Remove::rmnsm}
];

(* === Setup: Create test manifold with Minkowski metric === *)
Print["Setting up test manifold GFTestM2 (2D Minkowski)..."];
Print[""];

DefManifold[GFTestM2, 2, {a, b, c, d, e, f, g, h, i}];
DefMetric[-1, gfTestEta[-a, -b], gfTestCD, SymbolOfCovD -> {";", "D"}, PrintAs -> "eta"];
DefTensor[gfTestA[-a], GFTestM2];

(* ================================================================ *)
(* === Test Group 1: AddGaugeFixingTerm (core primitive)         === *)
(* ================================================================ *)

Print["Test Group 1: AddGaugeFixingTerm"];

(* Test 1.1: Adding zero returns original Lagrangian *)
TestCase["AddGaugeFixingTerm with zero term",
  AddGaugeFixingTerm[L, 0],
  L
];

(* Test 1.2: Adding a term returns L + term *)
TestCase["AddGaugeFixingTerm adds term to Lagrangian",
  AddGaugeFixingTerm[L, gTerm],
  L + gTerm
];

(* Test 1.3: Commutativity *)
TestCase["AddGaugeFixingTerm is commutative with Plus",
  AddGaugeFixingTerm[a + b, c],
  a + b + c
];

(* ================================================================ *)
(* === Test Group 2: BuildLorenzGaugeTerm                        === *)
(* ================================================================ *)

Print[""];
Print["Test Group 2: BuildLorenzGaugeTerm"];

(* Test 2.1: Result is non-zero *)
lorenzTerm = BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, 1];
TestTrue["Lorenz gauge term is non-zero",
  lorenzTerm =!= 0
];

(* Test 2.2: Result contains the field *)
TestTrue["Lorenz gauge term contains field A",
  !FreeQ[lorenzTerm, gfTestA]
];

(* Test 2.3: Result contains the covariant derivative *)
TestTrue["Lorenz gauge term contains covariant derivative",
  !FreeQ[lorenzTerm, gfTestCD]
];

(* Test 2.4: xi=1 is the default (both forms agree) *)
lorenzTermDefault = BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD];
TestCase["Default xi=1 matches explicit xi=1",
  lorenzTermDefault,
  lorenzTerm
];

(* Test 2.5: Different xi values produce different terms *)
lorenzTerm2 = BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, 2];
TestTrue["xi=2 differs from xi=1",
  lorenzTerm =!= lorenzTerm2
];

(* Test 2.6: Invalid xi throws *)
TestTrue["BuildLorenzGaugeTerm throws for xi=0",
  Catch[BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, 0]] =!= Null &&
  StringQ[Catch[BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, 0]]]
];

TestTrue["BuildLorenzGaugeTerm throws for xi=-1",
  StringQ[Catch[BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, -1]]]
];

(* ================================================================ *)
(* === Test Group 3: Integration (AddGaugeFixingTerm + Builder)  === *)
(* ================================================================ *)

Print[""];
Print["Test Group 3: Integration"];

(* Test 3.1: AddGaugeFixingTerm[L, BuildLorenzGaugeTerm[...]] produces L + gauge term *)
lagrangian = L;
gaugeTerm = BuildLorenzGaugeTerm[gfTestA, gfTestEta, gfTestCD, 1];
fixedLagrangian = AddGaugeFixingTerm[lagrangian, gaugeTerm];
TestTrue["Gauge-fixed Lagrangian differs from original",
  fixedLagrangian =!= lagrangian
];

TestTrue["Gauge-fixed Lagrangian contains original L",
  !FreeQ[fixedLagrangian, L]
];

(* ================================================================ *)
(* === Test Group 4: BuildDeDonderGaugeTerm                       === *)
(* ================================================================ *)

Print[""];
Print["Test Group 4: BuildDeDonderGaugeTerm"];

(* Setup: Need a symmetric rank-2 tensor for de Donder gauge *)
Quiet[
  DefTensor[gfTestH[-a, -b], GFTestM2, Symmetric[{-a, -b}]],
  {DefTensor::ntensoraliases}
];

(* Test 4.1: Result is non-zero *)
deDonderTerm = BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD, 1];
TestTrue["De Donder gauge term is non-zero",
  deDonderTerm =!= 0
];

(* Test 4.2: Result contains the field *)
TestTrue["De Donder gauge term contains field h",
  !FreeQ[deDonderTerm, gfTestH]
];

(* Test 4.3: Result contains the covariant derivative *)
TestTrue["De Donder gauge term contains covariant derivative",
  !FreeQ[deDonderTerm, gfTestCD]
];

(* Test 4.4: xi=1 is the default (both forms agree) *)
deDonderDefault = BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD];
TestCase["De Donder default xi=1 matches explicit xi=1",
  deDonderDefault,
  deDonderTerm
];

(* Test 4.5: Different xi values produce different terms *)
deDonderTerm2 = BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD, 2];
TestTrue["De Donder xi=2 differs from xi=1",
  deDonderTerm =!= deDonderTerm2
];

(* Test 4.6: Invalid xi throws *)
TestTrue["BuildDeDonderGaugeTerm throws for xi=0",
  StringQ[Catch[BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD, 0]]]
];

TestTrue["BuildDeDonderGaugeTerm throws for xi=-1",
  StringQ[Catch[BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD, -1]]]
];

(* Test 4.7: Integration — AddGaugeFixingTerm with de Donder *)
dDGaugeTerm = BuildDeDonderGaugeTerm[gfTestH, gfTestEta, gfTestCD, 1];
dDFixedLag = AddGaugeFixingTerm[L, dDGaugeTerm];
TestTrue["De Donder gauge-fixed Lagrangian differs from original",
  dDFixedLag =!= L
];

(* === Print summary === *)
Print[""];
PrintTestSummary[];
