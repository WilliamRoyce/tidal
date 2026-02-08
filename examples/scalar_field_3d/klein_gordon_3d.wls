#!/usr/bin/env wolframscript
(* Klein-Gordon 3+1D - Massive scalar field in flat Minkowski spacetime *)
(* Derives the Klein-Gordon equation from its Lagrangian in 4D *)
(*
   Physics: L = -1/2 (nabla_a phi)(nabla^a phi) - 1/2 m^2 phi^2

   In 3+1D Minkowski with signature (-,+,+,+), the Euler-Lagrange equation is:
     Box[phi] - m^2 phi = 0
   which expands to:
     -d^2_t phi + d^2_x phi + d^2_y phi + d^2_z phi - m^2 phi = 0

   Rearranging for the PDE pipeline (solve for d^2_t phi):
     d^2_t phi = d^2_x phi + d^2_y phi + d^2_z phi - m^2 phi

   Output: examples/data/klein_gordon_3d.json
*)

Print["=== Klein-Gordon 3+1D (Massive Scalar Field) ==="];
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

$DefInfoQ = False;

(* ============================================================ *)
(* Step 1: Define 3+1D Minkowski Spacetime *)
(* ============================================================ *)
Print[""];
Print["Step 1: Setting up 3+1D Minkowski spacetime..."];

(* Define 4D manifold (KG-prefixed to avoid symbol conflicts) *)
If[!xTensorQ[kgM4],
  DefManifold[kgM4, 4, {a, b, c, d, e, f, g, h}]
];

(* Define flat Minkowski metric *)
If[!MetricQ[kgEta],
  DefMetric[-1, kgEta[-a, -b], kgCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"]
];

(* Define coordinate chart: {t, x, y, z} *)
If[!ChartQ[kgCart],
  DefChart[kgCart, kgM4, {0, 1, 2, 3}, {t[], x[], y[], z[]}]
];

(* Set flat Minkowski metric components: diag(-1, +1, +1, +1) *)
kgMetricMatrix = DiagonalMatrix[{-1, 1, 1, 1}];
MetricInBasis[kgEta, -kgCart, kgMetricMatrix];

Print["   Manifold: kgM4 (3+1 dimensional)"];
Print["   Metric: eta_ab = diag(-1, +1, +1, +1) (flat Minkowski)"];
Print["   Coordinates: {t, x, y, z}"];

(* ============================================================ *)
(* Step 2: Define the Scalar Field *)
(* ============================================================ *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[phi], DefTensor[phi[], kgM4]];

(* Symbolic mass parameter (resolved at simulation time) *)
m2 = Symbol["m2"];
Print["   Scalar field: phi"];
Print["   Mass parameter: m^2 (symbolic, resolved at runtime)"];

(* ============================================================ *)
(* Step 3: Construct the Klein-Gordon Lagrangian *)
(* ============================================================ *)
Print[""];
Print["Step 3: Constructing Klein-Gordon Lagrangian..."];

(* L = -1/2 (nabla_a phi)(nabla^a phi) - 1/2 m^2 phi^2 *)
KGLagrangian = -1/2 kgCD[-a][phi[]] kgEta[a, b] kgCD[-b][phi[]] - 1/2 m2 phi[]^2;

Print["   L = -1/2 (nabla_a phi)(nabla^a phi) - 1/2 m^2 phi^2"];
Print["   Symbolic form: ", KGLagrangian];

(* ============================================================ *)
(* Step 4: Derive Euler-Lagrange Equation *)
(* ============================================================ *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equation..."];

eom = EulerLagrangeEquation[KGLagrangian, phi[], kgCD];
Print["   EOM: ", eom];

(* ============================================================ *)
(* Step 5: Decompose to Components *)
(* ============================================================ *)
Print[""];
Print["Step 5: Decomposing to coordinate components..."];

componentEqs = DecomposeToComponents[eom, phi[], kgCart];
Print["   Component equations: ", Length[componentEqs]];
Print["   Components: ", componentEqs];

(* ============================================================ *)
(* Step 6: Build and Export JSON *)
(* ============================================================ *)
Print[""];
Print["Step 6: Building and exporting JSON..."];

metadata = <|
  "field_prefix" -> "phi",
  "lagrangian_expr" -> "L = -1/2 eta^{ab} d_a phi d_b phi - 1/2 m^2 phi^2",
  "gauge" -> "none",
  "linearized" -> False,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

jsonOutput = BuildJSONStructure[componentEqs, metadata];
Print["   JSON structure created"];

(* Export JSON *)
outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputPath], CreateDirectory[outputPath]];
jsonFile = FileNameJoin[{outputPath, "klein_gordon_3d.json"}];

jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

Print[""];
Print["*** 3+1D Klein-Gordon derivation complete! ***"];
Print["   EOM: Box[phi] - m^2 phi = 0"];
Print["   d^2_t(phi) = laplacian_x(phi) + laplacian_y(phi) + laplacian_z(phi) - m^2 * phi"];
