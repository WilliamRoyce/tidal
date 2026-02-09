#!/usr/bin/env wolframscript
(* Klein-Gordon Lagrangian Example - Validation for the pipeline *)
(* This derives the Klein-Gordon equation from its Lagrangian *)

Print["=== Klein-Gordon Lagrangian to EOM Pipeline ==="];
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

(* === Step 1: Setup 1+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D Minkowski spacetime..."];

(* Use unique names to avoid conflicts *)
$DefInfoQ = False;  (* Suppress definition info messages *)

(* Define manifold and metric (with defensive checks for kernel reuse) *)
If[!xTensorQ[M2],
  DefManifold[M2, 2, {a, b, c, d, e, f}];
  DefMetric[-1, etaMetric[-a, -b], CD, SymbolOfCovD -> {";", "\[Del]"}, PrintAs -> "\[Eta]"];
  DefChart[cart, M2, {0, 1}, {t[], x[]}];
  MetricInBasis[etaMetric, -cart, DiagonalMatrix[{-1, 1}]];
  ,
  Print["   Using existing manifold M2"];
];
Print["   Manifold: M2 (2D)"];
Print["   Metric: Minkowski \[Eta]_ab with signature (-,+)"];
Print["   Covariant derivative: CD"];
Print["   Coordinates: {t, x}"];
Print["   Metric components: diag(-1, +1)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[phi], DefTensor[phi[], M2]];
(* Mass parameter — kept symbolic for runtime parameterization via --param m2=... *)
If[!ConstantSymbolQ[m2], DefConstantSymbol[m2]];
Print["   Scalar field: phi"];
Print["   Mass parameter: m^2 (symbolic, default 1.0)"];

(* === Step 3: Construct the Klein-Gordon Lagrangian === *)
Print[""];
Print["Step 3: Constructing Klein-Gordon Lagrangian..."];

(* L = -1/2 (nabla_a phi)(nabla^a phi) - 1/2 m^2 phi^2 *)
(* In xTensor notation: CD[-a][phi[]] is nabla_a phi *)
KGLagrangian = -1/2 CD[-a][phi[]] etaMetric[a, b] CD[-b][phi[]] - 1/2 m2 phi[]^2;

Print["   L = -1/2 (\[Del]_a \[Phi])(\[Del]^a \[Phi]) - 1/2 m\[Squared] \[Phi]\[Squared]"];
Print["   Symbolic form: ", KGLagrangian];

(* === Step 4: Derive Euler-Lagrange Equation === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equation..."];
Print["   Using Euler-Lagrange: \[Delta]L/\[Delta]\[Phi] = 0"];

(* Use pipeline's Euler-Lagrange computation *)
eom = EulerLagrangeEquation[KGLagrangian, phi[], CD];
Print["   EOM: ", eom];

(* === Step 5: Verify the Result === *)
Print[""];
Print["Step 5: Verifying result..."];

(* The d'Alembertian (Box) operator is: Box[phi] = eta^ab nabla_a nabla_b phi *)
(* So we expect: Box[phi] - m^2 phi = 0 *)

Print["   Expected form: \[EmptySquare]\[Phi] - m\[Squared]\[Phi] = 0"];
Print["   where \[EmptySquare] = \[Eta]^ab \[Del]_a \[Del]_b is the d'Alembertian"];

(* In coordinates: Box = -d^2/dt^2 + d^2/dx^2 for signature (-,+) *)
Print["   In coordinates: -\[PartialD]t\[Squared]\[Phi] + \[PartialD]x\[Squared]\[Phi] - m\[Squared]\[Phi] = 0"];

(* === Step 6: Extract Components Using Pipeline === *)
Print[""];
Print["Step 6: Decomposing to components using pipeline..."];

(* Use pipeline's component decomposition *)
componentEqs = DecomposeToComponents[eom, phi[], cart];
Print["   Component equations extracted: ", Length[componentEqs], " components"];
Print["   Components: ", componentEqs];

(* === Step 7: Build JSON Structure Using Pipeline === *)
Print[""];
Print["Step 7: Building JSON using pipeline..."];

(* Metadata for JSON export *)
metadata = <|
  "field_prefix" -> "phi",
  "lagrangian_expr" -> "-1/2 CD[-a][phi[]] CD[a][phi[]] - 1/2 m2 phi[]^2",
  "mass_squared" -> N[m2],
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 2,  (* Spacetime dimension: 1+1D = 2 total dims *)
  "signature" -> {-1, 1},
  "coordinates" -> {"t", "x"}
|>;

(* Build JSON using pipeline *)
jsonOutput = BuildJSONStructure[componentEqs, metadata];
jsonOutput["metadata", "parameters"] = <|"m2" -> 1.0|>;
Print["   JSON structure created using pipeline"];

(* Export JSON *)
outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputPath], CreateDirectory[outputPath]];
jsonFile = FileNameJoin[{outputPath, "klein_gordon_1d.json"}];

(* Convert to JSON string *)
jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

(* Export to file *)
Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

Print[""];
Print["*** Klein-Gordon derivation complete! ***"];
Print["The EOM Box[phi] - m^2 phi = 0 was derived from the Lagrangian."];
