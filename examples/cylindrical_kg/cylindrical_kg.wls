#!/usr/bin/env wolframscript
(* Klein-Gordon Field in Cylindrical Coordinates *)
(* Demonstrates Christoffel corrections from 3+1D cylindrical coordinate metric *)
(*
   Physics: Scalar field on flat 3D space using cylindrical coordinates (r, theta, z).

   The metric in cylindrical coordinates:
     ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2 + dz^2

   Using generic coordinate names {t, x, y, z} = {t, r, theta, z}:
     g_ab = diag(-1, 1, x^2, 1)

   The d'Alembertian in cylindrical coordinates:
     Box_g phi = -d2_t phi + d2_r phi + (1/r) d_r phi + (1/r^2) d2_theta phi + d2_z phi

   Key features:
     - Same 1/r and 1/r^2 corrections as polar, plus flat z-direction
     - Demonstrates 3+1D with mixed curved/flat spatial directions
     - laplacian_z has constant coefficient 1 (flat direction)

   Parameters: cylm2 (mass squared)
*)

Print["=== Klein-Gordon in Cylindrical Coordinates ==="];
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

(* === Step 1: Setup 3+1D Spacetime with Cylindrical Metric === *)
Print[""];
Print["Step 1: Setting up 3+1D spacetime with cylindrical coordinate metric..."];

$DefInfoQ = False;

(* Define mass parameter *)
If[!ConstantSymbolQ[cylm2],
  DefConstantSymbol[cylm2, PrintAs -> "m2"]
];

(* Define 3+1D manifold *)
If[!xTensorQ[cylM4],
  DefManifold[cylM4, 4, {a, b, c, d, e, f, g, h, i, j}]
];

(* Define metric *)
If[!MetricQ[cylg],
  DefMetric[-1, cylg[-a, -b], cylCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y, z} representing {t, r, theta, z} *)
If[!ChartQ[cylcart],
  DefChart[cylcart, cylM4, {0, 1, 2, 3}, {t[], x[], y[], z[]}]
];

(* Cylindrical coordinate metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 + dz^2 *)
cylMetricMatrix = DiagonalMatrix[{-1, 1, x[]^2, 1}];
MetricInBasis[cylg, -cylcart, cylMetricMatrix];

Print["   Manifold: cylM4 (3+1 dimensional)"];
Print["   Metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 + dz^2"];
Print["   Coordinates: {t, x, y, z} = {t, r, theta, z}"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[cylphi],
  DefTensor[cylphi[], cylM4, PrintAs -> "\[Phi]"]
];

(* === Step 3: Construct the Wave Equation === *)
Print[""];
Print["Step 3: Constructing KG equation in cylindrical coordinates..."];

waveEq = cylg[b, c] cylCD[-b][cylCD[-c][cylphi[]]] - cylm2 cylphi[];

Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];

(* === Step 4: Decompose to Components via Pipeline === *)
Print[""];
Print["Step 4: Decomposing to components using pipeline..."];

componentEqs = DecomposeToComponents[waveEq, cylphi[], cylcart, {},
  "MetricMatrix" -> cylMetricMatrix];

Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 5: Build and Export JSON === *)
Print[""];
Print["Step 5: Building JSON structure via pipeline..."];

metadata = <|
  "name" -> "Klein-Gordon in Cylindrical Coordinates",
  "field_prefix" -> "cylphi",
  "lagrangian_expr" -> "-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

jsonOutput = BuildJSONStructure[componentEqs, metadata];

(* Add parameters and metric info *)
jsonOutput["metadata", "parameters"] = <|"cylm2" -> 0.0|>;
jsonOutput["spacetime", "metric_type"] = "cylindrical";
jsonOutput["spacetime", "metric_components"] = "diag(-1, 1, x^2, 1)";

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
jsonFile = FileNameJoin[{outputDir, "cylindrical_kg.json"}];

jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

(* === Summary === *)
Print[""];
Print["=== Summary ==="];
Print["From the cylindrical-coordinate Klein-Gordon equation, the pipeline derived:"];
Print[""];
Print["  Key features:"];
Print["    1. Flat 3D space in cylindrical coordinates"];
Print["    2. Radial Christoffel: (1/x) gradient_x"];
Print["    3. Angular Laplacian: (1/x^2) laplacian_y"];
Print["    4. Flat z-direction: laplacian_z with constant coefficient"];
Print["    5. Mixed curved (r, theta) and flat (z) directions"];
Print[""];
Print["*** Klein-Gordon in Cylindrical Coordinates pipeline complete! ***"];
