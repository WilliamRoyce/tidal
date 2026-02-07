#!/usr/bin/env wolframscript
(* Klein-Gordon Field in Spherical Coordinates *)
(* Demonstrates Christoffel corrections from 3+1D spherical coordinate metric *)
(*
   Physics: Scalar field on flat 3D space using spherical coordinates (r, theta, phi).

   The metric in spherical coordinates:
     ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2 + r^2 sin^2(theta) d(phi)^2

   Using generic coordinate names {t, x, y, z} = {t, r, theta, phi}:
     g_ab = diag(-1, 1, x^2, x^2 Sin[y]^2)

   The d'Alembertian in spherical coordinates:
     Box_g phi = -d2_t phi + d2_r phi + (2/r) d_r phi
                 + (1/r^2) d2_theta phi + (cos(theta)/(r^2 sin(theta))) d_theta phi
                 + (1/(r^2 sin^2(theta))) d2_phi phi

   Key features:
     - gradient_x term: (2/r) d_r phi from radial Christoffels
     - laplacian_y term: (1/r^2) d2_theta phi from g^{yy}
     - gradient_y term: cot(theta)/r^2 from angular Christoffels
     - laplacian_z term: 1/(r^2 sin^2(theta)) from g^{zz}
     - Trigonometric coefficients: exercises Sin/Cos in _mathematica_to_python

   Parameters: spm2 (mass squared)
*)

Print["=== Klein-Gordon in Spherical Coordinates ==="];
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

(* === Step 1: Setup 3+1D Spacetime with Spherical Metric === *)
Print[""];
Print["Step 1: Setting up 3+1D spacetime with spherical coordinate metric..."];

$DefInfoQ = False;

(* Define mass parameter *)
If[!ConstantSymbolQ[spm2],
  DefConstantSymbol[spm2, PrintAs -> "m2"]
];

(* Define 3+1D manifold *)
If[!xTensorQ[spM4],
  DefManifold[spM4, 4, {a, b, c, d, e, f, g, h, i, j}]
];

(* Define metric *)
If[!MetricQ[spg],
  DefMetric[-1, spg[-a, -b], spCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y, z} representing {t, r, theta, phi} *)
If[!ChartQ[spcart],
  DefChart[spcart, spM4, {0, 1, 2, 3}, {t[], x[], y[], z[]}]
];

(* Spherical coordinate metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 + x^2 Sin[y]^2 dz^2 *)
spMetricMatrix = DiagonalMatrix[{-1, 1, x[]^2, x[]^2 Sin[y[]]^2}];
MetricInBasis[spg, -spcart, spMetricMatrix];

Print["   Manifold: spM4 (3+1 dimensional)"];
Print["   Metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 + x^2 Sin[y]^2 dz^2"];
Print["   Coordinates: {t, x, y, z} = {t, r, theta, phi}"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[spphi],
  DefTensor[spphi[], spM4, PrintAs -> "\[Phi]"]
];

(* === Step 3: Construct the Wave Equation === *)
Print[""];
Print["Step 3: Constructing KG equation in spherical coordinates..."];

waveEq = spg[b, c] spCD[-b][spCD[-c][spphi[]]] - spm2 spphi[];

Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];

(* === Step 4: Decompose to Components via Pipeline === *)
Print[""];
Print["Step 4: Decomposing to components using pipeline..."];

componentEqs = DecomposeToComponents[waveEq, spphi[], spcart, {},
  "ComputeChristoffels" -> True,
  "MetricMatrix" -> spMetricMatrix];

Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 5: Build and Export JSON === *)
Print[""];
Print["Step 5: Building JSON structure via pipeline..."];

metadata = <|
  "name" -> "Klein-Gordon in Spherical Coordinates",
  "field_prefix" -> "spphi",
  "lagrangian_expr" -> "-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

jsonOutput = BuildJSONStructure[componentEqs, metadata];

(* Add parameters and metric info *)
jsonOutput["metadata", "parameters"] = <|"spm2" -> 0.0|>;
jsonOutput["spacetime", "metric_type"] = "spherical";
jsonOutput["spacetime", "metric_components"] = "diag(-1, 1, x^2, x^2*Sin[y]^2)";

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
jsonFile = FileNameJoin[{outputDir, "spherical_kg.json"}];

jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

(* === Summary === *)
Print[""];
Print["=== Summary ==="];
Print["From the spherical-coordinate Klein-Gordon equation, the pipeline derived:"];
Print[""];
Print["  Key features:"];
Print["    1. Flat 3D space in spherical coordinates"];
Print["    2. Radial Christoffel: (2/x) gradient_x"];
Print["    3. Polar Christoffel: Cos[y]/(x^2 Sin[y]) gradient_y"];
Print["    4. Angular Laplacians with position-dependent coefficients"];
Print["    5. Trigonometric functions in coefficients (Sin, Cos)"];
Print[""];
Print["*** Klein-Gordon in Spherical Coordinates pipeline complete! ***"];
