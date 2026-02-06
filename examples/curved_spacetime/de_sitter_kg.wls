#!/usr/bin/env wolframscript
(* De Sitter Klein-Gordon Field - 2+1D Time-Dependent Conformal Factor *)
(* Demonstrates computed Christoffel symbols and Hubble friction from metric *)
(*
   Physics: Klein-Gordon field in 2+1D de Sitter (exponentially expanding) spacetime

   Metric: ds^2 = Omega(t)^2 (-dt^2 + dx^2 + dy^2)  where Omega(t) = exp(H*t)

   The conformal d'Alembertian in (n+1) dimensions with Omega = exp(Ht):
     Box_g phi = Omega^{-2}[-d2_t phi + nabla^2 phi - (n-1)*H*d_t phi]

   For 2+1D (n=2 spatial dimensions):
     Box_g phi - m^2 phi = 0  =>
     d2_t phi = nabla^2 phi - H * d_t phi - m^2 * exp(2Ht) * phi

   Key features:
     - Hubble friction: -(n-1)*H * d_t phi = -H * d_t phi for n=2
     - Time-dependent effective mass: m^2 * Omega^2 = m^2 * exp(2Ht)
     - Christoffel symbols computed by xAct from metric definition
     - Energy NOT conserved: waves lose energy to cosmic expansion

   Note: In 1+1D (n=1), friction coefficient = (n-1)*H = 0 (conformal invariance).
   We use 2+1D specifically because Hubble friction genuinely appears.
*)

Print["=== De Sitter Klein-Gordon 2+1D (Hubble Friction) ==="];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Load pipeline modules *)
Print["Loading pipeline modules..."];
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "torsion_gertsenshtein", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

(* === Step 1: Setup 2+1D de Sitter Spacetime === *)
Print[""];
Print["Step 1: Setting up 2+1D de Sitter spacetime..."];

$DefInfoQ = False;  (* Suppress definition info messages *)

(* Define Hubble parameter as symbolic constant *)
If[!ConstantSymbolQ[dSH],
  DefConstantSymbol[dSH, PrintAs -> "H"]
];

(* Define mass parameter as symbolic constant *)
If[!ConstantSymbolQ[dSm2],
  DefConstantSymbol[dSm2, PrintAs -> "m\[Squared]"]
];

(* Define 2+1D manifold *)
If[!xTensorQ[dSM3],
  DefManifold[dSM3, 3, {a, b, c, d, e, f}]
];

(* Define conformal metric *)
If[!MetricQ[dSg],
  DefMetric[-1, dSg[-a, -b], dSCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y} *)
(* Use t[], x[], y[] (not prefixed) so IsTimeDependentCoefficient detects t[] *)
If[!ChartQ[dScart],
  DefChart[dScart, dSM3, {0, 1, 2}, {t[], x[], y[]}]
];

(* Set metric components: g_ab = Omega^2 eta_ab where Omega = exp(H*t) *)
dSMetricMatrix = Exp[2*dSH*t[]] * DiagonalMatrix[{-1, 1, 1}];
MetricInBasis[dSg, -dScart, dSMetricMatrix];

Print["   Manifold: dSM3 (2+1 dimensional de Sitter)"];
Print["   Conformal factor: Omega(t) = exp(H*t)"];
Print["   Metric: g_ab = exp(2Ht) * diag(-1, +1, +1)"];
Print["   Hubble parameter: H = ", dSH, " (symbolic)"];
Print["   Mass parameter: m^2 = ", dSm2, " (symbolic)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[dSphi],
  DefTensor[dSphi[], dSM3, PrintAs -> "\[Phi]"]
];

Print["   Scalar field: dSphi (Klein-Gordon)"];

(* === Step 3: Construct the Wave Equation === *)
Print[""];
Print["Step 3: Constructing de Sitter Klein-Gordon equation..."];

(*
   Minimally coupled KG in curved spacetime:
     Box_g phi - m^2 phi = 0
   where Box_g = g^{ab} nabla_a nabla_b

   For 2+1D with g_ab = exp(2Ht) * eta_ab (n=2 spatial dims):
     Box_g phi = Omega^{-2}[-d2_t phi + d2_x phi + d2_y phi - (n-1)*H*d_t phi]
               = Omega^{-2}[-d2_t phi + nabla^2 phi - H*d_t phi]

   After normalization (multiply by -Omega^2):
     d2_t phi = nabla^2 phi - H*d_t phi - m^2*Omega^2*phi
*)

waveEq = dSg[b, c] dSCD[-b][dSCD[-c][dSphi[]]] - dSm2 dSphi[];

Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];
Print["   In 2+1D: d2_t phi = nabla^2 phi - H d_t phi - m^2 exp(2Ht) phi"];
Print["                                       ^^^ (n-1)H = H for n=2"];

(* === Step 4: Decompose to Components via Pipeline === *)
Print[""];
Print["Step 4: Decomposing to components using pipeline..."];

componentEqs = DecomposeToComponents[waveEq, dSphi[], dScart, {},
  "ComputeChristoffels" -> True,
  "MetricMatrix" -> dSMetricMatrix];

Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 5: Build and Export JSON === *)
Print[""];
Print["Step 5: Building JSON structure via pipeline..."];

metadata = <|
  "name" -> "De Sitter Klein-Gordon 2+1D (Hubble Friction)",
  "field_prefix" -> "dSphi",
  "lagrangian_expr" -> "-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"}
|>;

(* Build JSON via pipeline *)
jsonOutput = BuildJSONStructure[componentEqs, metadata];

(* Add parameters and metric info to JSON metadata *)
jsonOutput["metadata", "parameters"] = <|"dSH" -> 0.1, "dSm2" -> 1.0|>;
jsonOutput["spacetime", "metric_type"] = "conformal_flat";
jsonOutput["spacetime", "conformal_factor"] = "exp(H*t)";

Print["   JSON structure created via pipeline"];

(* === Step 6: Export JSON === *)
Print[""];
Print["Step 6: Exporting to JSON file..."];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "",
  scriptDir = Directory[]
];
outputDir = FileNameJoin[{scriptDir, "..", "data"}];
If[!DirectoryQ[outputDir],
  CreateDirectory[outputDir]
];
jsonFile = FileNameJoin[{outputDir, "de_sitter_kg.json"}];

(* Export *)
jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

(* === Summary === *)
Print[""];
Print["=== Summary ==="];
Print["From the de Sitter conformal Klein-Gordon equation, the pipeline derived:"];
Print[""];
Print["  Wave equation (2+1D): d2_t phi = nabla^2 phi - H d_t phi - m^2 exp(2Ht) phi"];
Print[""];
Print["  Key features:"];
Print["    1. Conformal metric: ds^2 = exp(2Ht)(-dt^2 + dx^2 + dy^2)"];
Print["    2. Hubble friction: -(n-1)*H d_t phi where n = spatial dimensions"];
Print["       For 2+1D (n=2): -H d_t phi"];
Print["       Note: In 1+1D (n=1), friction = 0 (conformal invariance)"];
Print["    3. Christoffel symbols computed from metric by xAct (not hardcoded)"];
Print["    4. Pipeline handles: actual metric -> LHS normalization -> operator detection"];
Print["    5. Time-dependent mass: m_eff^2(t) = m^2 * exp(2Ht)"];
Print[""];
Print["  Physical interpretation:"];
Print["    - Waves in expanding spacetime lose energy to expansion"];
Print["    - Hubble friction (d_t phi term) causes amplitude decay"];
Print["    - Energy is NOT conserved - lost to cosmic expansion"];
Print[""];
Print["  Runtime parameters:"];
Print["    - dSH: Hubble constant (default 0.1)"];
Print["    - dSm2: Mass squared (default 1.0)"];
Print[""];
Print["*** De Sitter Klein-Gordon 2+1D pipeline complete! ***"];
