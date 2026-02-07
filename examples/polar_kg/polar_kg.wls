#!/usr/bin/env wolframscript
(* Klein-Gordon Field in Polar Coordinates *)
(* Demonstrates Christoffel corrections from non-Cartesian coordinate metric *)
(*
   Physics: Scalar field on flat 2D space using polar coordinates (r, theta).

   The metric in polar coordinates:
     ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2

   Using generic coordinate names {t, x, y} = {t, r, theta}:
     g_ab = diag(-1, 1, x^2)

   The d'Alembertian in polar coordinates:
     Box_g phi = -d2_t phi + d2_r phi + (1/r) d_r phi + (1/r^2) d2_theta phi

   Key features:
     - gradient_x term: (1/r) d_r phi from Christoffel symbol Gamma^x_{yy} = -x
     - laplacian_y term: (1/r^2) d2_theta phi from metric component g^{yy} = 1/x^2
     - All corrections derived automatically from the metric by xAct

   Parameters: polm2 (mass squared)
*)

Print["=== Klein-Gordon in Polar Coordinates ==="];
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

(* === Step 1: Setup 2+1D Spacetime with Polar Metric === *)
Print[""];
Print["Step 1: Setting up 2+1D spacetime with polar coordinate metric..."];

$DefInfoQ = False;  (* Suppress definition info messages *)

(* Define mass parameter as symbolic constant *)
If[!ConstantSymbolQ[polm2],
  DefConstantSymbol[polm2, PrintAs -> "m\[Squared]"]
];

(* Define 2+1D manifold *)
If[!xTensorQ[polM3],
  DefManifold[polM3, 3, {a, b, c, d, e, f}]
];

(* Define metric *)
If[!MetricQ[polg],
  DefMetric[-1, polg[-a, -b], polCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y} representing {t, r, theta} *)
If[!ChartQ[polcart],
  DefChart[polcart, polM3, {0, 1, 2}, {t[], x[], y[]}]
];

(* Polar coordinate metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 *)
(* where x = r, y = theta *)
polMetricMatrix = DiagonalMatrix[{-1, 1, x[]^2}];
MetricInBasis[polg, -polcart, polMetricMatrix];

Print["   Manifold: polM3 (2+1 dimensional: time x R^2_polar)"];
Print["   Metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2"];
Print["   Coordinates: {t, x, y} = {t, r, theta}"];
Print["   Mass parameter: m^2 = ", polm2, " (symbolic)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[polphi],
  DefTensor[polphi[], polM3, PrintAs -> "\[Phi]"]
];

Print["   Scalar field: polphi (Klein-Gordon in polar coordinates)"];

(* === Step 3: Construct the Wave Equation === *)
Print[""];
Print["Step 3: Constructing KG equation in polar coordinates..."];

(*
   Minimally coupled KG in curved spacetime:
     Box_g phi - m^2 phi = 0
   where Box_g = g^{ab} nabla_a nabla_b
*)

waveEq = polg[b, c] polCD[-b][polCD[-c][polphi[]]] - polm2 polphi[];

Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];
Print["   Expected: d2_t phi = d2_r phi + (1/r) d_r phi + (1/r^2) d2_theta phi - m^2 phi"];

(* === Step 4: Decompose to Components via Pipeline === *)
Print[""];
Print["Step 4: Decomposing to components using pipeline..."];

componentEqs = DecomposeToComponents[waveEq, polphi[], polcart, {},
  "ComputeChristoffels" -> True,
  "MetricMatrix" -> polMetricMatrix];

Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 5: Build and Export JSON === *)
Print[""];
Print["Step 5: Building JSON structure via pipeline..."];

metadata = <|
  "name" -> "Klein-Gordon in Polar Coordinates",
  "field_prefix" -> "polphi",
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
jsonOutput["metadata", "parameters"] = <|"polm2" -> 1.0|>;
jsonOutput["spacetime", "metric_type"] = "polar";
jsonOutput["spacetime", "metric_components"] = "diag(-1, 1, x^2)";

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
jsonFile = FileNameJoin[{outputDir, "polar_kg.json"}];

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
Print["From the polar-coordinate Klein-Gordon equation, the pipeline derived:"];
Print[""];
Print["  Key features:"];
Print["    1. Flat space, non-Cartesian coordinates: ds^2 = -dt^2 + dx^2 + x^2 dy^2"];
Print["    2. Christoffel correction: (1/x) gradient_x from Gamma^x_{yy} = -x"];
Print["    3. Angular Laplacian: (1/x^2) laplacian_y from g^{yy} = 1/x^2"];
Print["    4. All coefficients derived from metric by xAct (not hardcoded)"];
Print[""];
Print["  Physical interpretation:"];
Print["    - Coordinate x represents radial distance r"];
Print["    - Coordinate y represents azimuthal angle theta"];
Print["    - The (1/r) gradient term accounts for cylindrical wave spreading"];
Print["    - The (1/r^2) angular Laplacian scales angular derivatives correctly"];
Print[""];
Print["  Runtime parameters:"];
Print["    - polm2: Mass squared (default 1.0)"];
Print[""];
Print["*** Klein-Gordon in Polar Coordinates pipeline complete! ***"];
