#!/usr/bin/env wolframscript
(* Klein-Gordon Field on Static 2-Sphere via Stereographic Projection *)
(* Demonstrates position-dependent Christoffel symbols from spatial curvature *)
(*
   Physics: Scalar field on S^2 x R_t (static 2-sphere, time is flat)

   Stereographic projection from the south pole maps the sphere to the plane:
     ds^2 = -dt^2 + Omega^2(x,y) (dx^2 + dy^2)
   where
     Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2)

   The south pole (x=y=0) has Omega = 2, the equator (r=R) has Omega = 1,
   and Omega -> 0 as r -> infinity (approaching the north pole).

   The d'Alembertian on this metric:
     Box_g phi = -d2_t phi + (1/Omega^2) nabla^2 phi + Christoffel corrections

   The key test here:
     - Time is flat: no Hubble friction
     - Spatial curvature: position-dependent wave speed 1/Omega^2
     - Non-trivial Christoffel gradient corrections from d_x Omega, d_y Omega
     - Known analytic: spherical harmonic eigenvalues l(l+1)/R^2

   Parameters: sphR (sphere radius), sphm2 (mass squared)
*)

Print["=== Klein-Gordon on 2-Sphere (Stereographic Projection) ==="];
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

(* === Step 1: Setup 2+1D Spacetime with Sphere Spatial Metric === *)
Print[""];
Print["Step 1: Setting up 2+1D spacetime with 2-sphere spatial metric..."];

$DefInfoQ = False;  (* Suppress definition info messages *)

(* Define sphere radius as symbolic constant *)
If[!ConstantSymbolQ[sphR],
  DefConstantSymbol[sphR, PrintAs -> "R"]
];

(* Define mass parameter as symbolic constant *)
If[!ConstantSymbolQ[sphm2],
  DefConstantSymbol[sphm2, PrintAs -> "m\[Squared]"]
];

(* Define 2+1D manifold *)
If[!xTensorQ[sphM3],
  DefManifold[sphM3, 3, {a, b, c, d, e, f}]
];

(* Define metric *)
If[!MetricQ[sphg],
  DefMetric[-1, sphg[-a, -b], sphCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y} *)
If[!ChartQ[sphcart],
  DefChart[sphcart, sphM3, {0, 1, 2}, {t[], x[], y[]}]
];

(* Conformal factor for stereographic projection *)
(* Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2) *)
sphOmega = 2 sphR^2 / (sphR^2 + x[]^2 + y[]^2);

(* Metric: ds^2 = -dt^2 + Omega^2 (dx^2 + dy^2) *)
(* Note: time component is -1 (NOT multiplied by Omega^2) *)
sphMetricMatrix = DiagonalMatrix[{-1, sphOmega^2, sphOmega^2}];
MetricInBasis[sphg, -sphcart, sphMetricMatrix];

Print["   Manifold: sphM3 (2+1 dimensional: time x S^2)"];
Print["   Conformal factor: Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2)"];
Print["   Metric: ds^2 = -dt^2 + Omega^2 (dx^2 + dy^2)"];
Print["   Sphere radius: R = ", sphR, " (symbolic)"];
Print["   Mass parameter: m^2 = ", sphm2, " (symbolic)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[sphphi],
  DefTensor[sphphi[], sphM3, PrintAs -> "\[Phi]"]
];

Print["   Scalar field: sphphi (Klein-Gordon on S^2)"];

(* === Step 3: Construct the Wave Equation === *)
Print[""];
Print["Step 3: Constructing KG equation on 2-sphere..."];

(*
   Minimally coupled KG in curved spacetime:
     Box_g phi - m^2 phi = 0
   where Box_g = g^{ab} nabla_a nabla_b

   For the product metric ds^2 = -dt^2 + Omega^2 (dx^2 + dy^2):
     g^{tt} = -1, g^{xx} = g^{yy} = 1/Omega^2

   After normalization:
     d2_t phi = (1/Omega^2) nabla^2 phi + (Christoffel gradient terms) - m^2 phi
*)

waveEq = sphg[b, c] sphCD[-b][sphCD[-c][sphphi[]]] - sphm2 sphphi[];

Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];
Print["   Expected: d2_t phi = (1/Omega^2) nabla^2 phi + gradient corrections - m^2 phi"];

(* === Step 4: Decompose to Components via Pipeline === *)
Print[""];
Print["Step 4: Decomposing to components using pipeline..."];

componentEqs = DecomposeToComponents[waveEq, sphphi[], sphcart, {},
  "MetricMatrix" -> sphMetricMatrix];

Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 5: Build and Export JSON === *)
Print[""];
Print["Step 5: Building JSON structure via pipeline..."];

metadata = <|
  "name" -> "Klein-Gordon on 2-Sphere (Stereographic Projection)",
  "field_prefix" -> "sphphi",
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
jsonOutput["metadata", "parameters"] = <|"sphR" -> 1.0, "sphm2" -> 0.0|>;
jsonOutput["spacetime", "metric_type"] = "product_conformal";
jsonOutput["spacetime", "conformal_factor"] = "2*R^2 / (R^2 + x^2 + y^2)";

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
jsonFile = FileNameJoin[{outputDir, "sphere_kg.json"}];

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
Print["From the 2-sphere Klein-Gordon equation, the pipeline derived:"];
Print[""];
Print["  Key features:"];
Print["    1. Static spatial curvature: ds^2 = -dt^2 + Omega^2(x,y)(dx^2 + dy^2)"];
Print["    2. Conformal factor: Omega = 2R^2 / (R^2 + x^2 + y^2)"];
Print["    3. Position-dependent wave speed: 1/Omega^2"];
Print["    4. Christoffel gradient corrections from spatial curvature"];
Print["    5. No Hubble friction (time is flat)"];
Print["    6. All coefficients derived from metric by xAct (not hardcoded)"];
Print[""];
Print["  Physical interpretation:"];
Print["    - Waves propagate slower near south pole (Omega = 2)"];
Print["    - Waves propagate faster near equator (Omega = 1)"];
Print["    - Curvature focuses/defocuses waves"];
Print["    - For l(l+1)/R^2 + m^2 eigenfrequency spectrum"];
Print[""];
Print["  Runtime parameters:"];
Print["    - sphR: Sphere radius (default 1.0)"];
Print["    - sphm2: Mass squared (default 0.0, massless)"];
Print[""];
Print["*** Klein-Gordon on 2-Sphere pipeline complete! ***"];
