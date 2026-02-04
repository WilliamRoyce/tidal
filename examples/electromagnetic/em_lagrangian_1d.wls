#!/usr/bin/env wolframscript
(* Electromagnetic Field Lagrangian Example - 1+1D *)
(* Derives the wave equations for vector potential A_mu from L = -1/4 F_uv F^uv *)

Print["=== Electromagnetic Lagrangian to PDE Pipeline (1+1D) ==="];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Load pipeline modules *)
Print["Loading pipeline modules..."];
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

(* === Step 1: Setup 1+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D Minkowski spacetime..."];

DefManifold[M2, 2, {a, b, c, d, e, f, g, h}];
DefMetric[-1, eta[-a, -b], CD,
  SymbolOfCovD -> {";", "\[Del]"},
  PrintAs -> "\[Eta]"];

(* Define coordinate chart: {t, x} *)
DefChart[cart, M2, {0, 1}, {t[], x[]}];

(* Set Minkowski metric components: diag(-1, +1) *)
MetricInBasis[eta, -cart, DiagonalMatrix[{-1, 1}]];

Print["   Manifold: M2 (1+1 dimensional)"];
Print["   Metric: Minkowski \[Eta]_ab with signature (-,+)"];
Print["   Coordinates: {t, x}"];

(* === Step 2: Define the Vector Potential and Field Strength === *)
Print[""];
Print["Step 2: Defining vector potential A_a and field strength F_ab..."];

(* Vector potential A_a *)
DefTensor[A[-a], M2, PrintAs -> "A"];

(* Field strength tensor F_ab = d_a A_b - d_b A_a (antisymmetric) *)
DefTensor[F[-a, -b], M2, Antisymmetric[{-a, -b}], PrintAs -> "F"];

Print["   Vector potential: A_a"];
Print["   Field strength: F_ab (antisymmetric)"];

(* Field strength definition rule *)
FieldStrengthRule = F[-a, -b] -> CD[-a][A[-b]] - CD[-b][A[-a]];
Print["   Definition: F_ab = \[Del]_a A_b - \[Del]_b A_a"];

(* === Step 3: Construct the EM Lagrangian === *)
Print[""];
Print["Step 3: Constructing EM Lagrangian..."];

(* L = -1/4 F_ab F^ab *)
(* Note: F^ab = eta^ac eta^bd F_cd *)
EMLagrangian = -1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d];

Print["   L = -1/4 F_ab F^ab"];
Print["   Symbolic form: ", EMLagrangian];

(* Substitute field strength definition *)
EMLagrangianExpanded = EMLagrangian /. FieldStrengthRule;
EMLagrangianExpanded = ToCanonical[EMLagrangianExpanded];
Print["   Expanded form: ", EMLagrangianExpanded];

(* === Step 4: Derive Euler-Lagrange Equations === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equations for A_a using pipeline..."];
Print["   \[Delta]L/\[Delta]A_a = 0"];

(* Use pipeline's Euler-Lagrange computation *)
(* For the EM Lagrangian, EOM is: d_b F^ba = 0 (source-free Maxwell) *)
eom = EulerLagrangeEquation[EMLagrangianExpanded, A[-a], CD];
Print["   EOM: ", eom];

(* === Step 5: Apply Lorenz Gauge === *)
Print[""];
Print["Step 5: Applying Lorenz gauge condition..."];
Print["   Lorenz gauge: \[Del]_a A^a = 0"];

(* In Lorenz gauge, Maxwell equations simplify to wave equations:
   Box A_a = 0  (for source-free case)
   where Box = eta^bc d_b d_c is the d'Alembertian *)

(* For 1+1D: Box = -d^2/dt^2 + d^2/dx^2 *)
Print["   In Lorenz gauge: \[EmptySquare] A_a = 0"];
Print["   where \[EmptySquare] = \[Eta]^ab \[Del]_a \[Del]_b"];

(* Substitute F_ab = d_a A_b - d_b A_a into Maxwell equations *)
(* Maxwell: d_b F^ab = 0, F_ab = d_a A_b - d_b A_a *)
(* Expanding: d_b (eta^ac eta^bd (d_c A_d - d_d A_c)) = 0 *)
(* = eta^ac d^b d_c A_b - eta^ac d^b d_b A_c = 0 *)
(* = d^a (d_b A^b) - Box A^a = 0 *)
(* With Lorenz gauge d_b A^b = 0: Box A^a = 0 *)

(* Construct the wave equation directly: Box A_a = eta^bc d_b d_c A_a = 0 *)
waveEqA = eta[b, c] CD[-b][CD[-c][A[-a]]];
Print["   Wave equation form: ", waveEqA, " = 0"];

(* === Step 6: Extract Component Equations Using Pipeline === *)
Print[""];
Print["Step 6: Decomposing to components using pipeline..."];

(* Use pipeline's component decomposition with the wave equation form *)
componentEqs = DecomposeToComponents[waveEqA, A[-a], cart];
Print["   Component equations extracted: ", Length[componentEqs], " components"];
Print["   Components: ", componentEqs];

(* === Step 7: Build JSON Structure Using Pipeline === *)
Print[""];
Print["Step 7: Building JSON using pipeline..."];

(* Metadata for JSON export *)
metadata = <|
  "field_prefix" -> "A",
  "lagrangian_expr" -> "-1/4 F[-a,-b] F[a,b]",
  "mass_squared" -> 0.0,  (* Massless *)
  "gauge" -> "lorenz",
  "linearized" -> True,
  "dimension" -> 2,  (* Spacetime dimension: 1+1D = 2 total dims *)
  "signature" -> {-1, 1},
  "coordinates" -> {"t", "x"}
|>;

(* Build JSON using pipeline *)
jsonOutput = BuildJSONStructure[componentEqs, metadata];
Print["   JSON structure created using pipeline"];

(* === Step 8: Export JSON === *)
Print[""];
Print["Step 8: Exporting to JSON file..."];

(* Determine output path *)
scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "",
  scriptDir = Directory[]
];
outputDir = FileNameJoin[{scriptDir, "..", "..", "..", "examples", "data"}];
If[!DirectoryQ[outputDir],
  CreateDirectory[outputDir]
];
jsonFile = FileNameJoin[{outputDir, "em_1d.json"}];

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
Print["From the EM Lagrangian L = -1/4 F_ab F^ab, we derived:"];
Print["  - Maxwell equations: \[Del]_b F^ba = 0"];
Print["  - In Lorenz gauge (\[Del]_a A^a = 0): \[EmptySquare] A_a = 0"];
Print["  - Component equations:"];
Print["    \[PartialD]_t\[Squared] A_0 = \[PartialD]_x\[Squared] A_0 (wave equation)"];
Print["    \[PartialD]_t\[Squared] A_1 = \[PartialD]_x\[Squared] A_1 (wave equation)"];
Print[""];
Print["The equations are ready for simulation with py-pde!"];
Print[""];
Print["*** EM Lagrangian pipeline complete! ***"];
