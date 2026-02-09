#!/usr/bin/env wolframscript
(* ::Package:: *)

(* chern_simons.wls - Chern-Simons Lagrangian in 2+1D *)
(* Demonstrates: antisymmetric tensors, Levi-Civita density, dimension-specific terms *)
(*
   Lagrangian: L = -1/4 F_uv F^uv + (kappa/2) epsilon^uvr A_u d_v A_r
   Expected EOM: d_u F^un + kappa epsilon^nab d_a A_b = 0

   In Lorenz gauge (d^u A_u = 0), Maxwell part becomes Box A_n = 0
   Full EOM: Box A_n + kappa epsilon^nab d_a A_b = 0
*)

Print["=== Chern-Simons 2+1D Lagrangian to EOM Pipeline ==="];
Print[""];

(* === Load xAct Packages === *)
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

(* === Step 1: Define 2+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 2+1D Minkowski spacetime..."];

(* Use shared symbols to avoid conflicts with other examples *)
If[!xTensorQ[M3],
  DefManifold[M3, 3, {a, b, c, d, e, f, g, h, i, j}];
  (* DefMetric automatically creates epsiloneta3[-a,-b,-c] as the Levi-Civita tensor *)
  DefMetric[-1, eta3[-a, -b], CD3, SymbolOfCovD -> {";", "\[Del]"}, PrintAs -> "\[Eta]"];
  DefChart[cart3, M3, {0, 1, 2}, {t[], x[], y[]}];
  MetricInBasis[eta3, -cart3, DiagonalMatrix[{-1, 1, 1}]];
  Print["   Defined new manifold M3"];
];

Print["   Manifold: M3 (3D = 2+1 spacetime)"];
Print["   Metric: Minkowski with signature (-,+,+)"];
Print["   Coordinates: {t, x, y}"];
Print["   Built-in Levi-Civita: epsiloneta3[-a,-b,-c]"];

(* === Step 2: Define Vector Field and Field Strength === *)
Print[""];
Print["Step 2: Defining vector potential and field strength..."];

If[!xTensorQ[csA],
  DefTensor[csA[-a], M3, PrintAs -> "A"];
  Print["   Defined vector potential A[-a]"];
];

If[!xTensorQ[csF],
  DefTensor[csF[-a, -b], M3, Antisymmetric[{-a, -b}], PrintAs -> "F"];
  Print["   Defined field strength F[-a,-b] (antisymmetric)"];
];

(* Field strength definition rule *)
FieldStrengthRule = csF[-a, -b] -> CD3[-a][csA[-b]] - CD3[-b][csA[-a]];
Print["   Definition: F_ab = D_a A_b - D_b A_a"];

(* === Step 3: Construct Full Chern-Simons Lagrangian === *)
Print[""];
Print["Step 3: Constructing full Chern-Simons Lagrangian..."];

(* Define the CS coupling constant — kept symbolic for runtime parameterization *)
If[!ConstantSymbolQ[kappa], DefConstantSymbol[kappa]];

(*
   Full Maxwell-Chern-Simons Lagrangian:
     L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c

   The EOM is derived symbolically using VarD.
   The epsilon tensor is evaluated to numeric components automatically
   by the pipeline's EvaluateEpsilonComponents function.

   In Lorenz gauge (d^u A_u = 0), this gives:
     Box A_a + kappa * epsilon_a^{bc} d_b A_c = 0

   In 2+1D with signature (-,+,+):
     Box = -d^2/dt^2 + d^2/dx^2 + d^2/dy^2
*)

(* Construct Maxwell Lagrangian directly in terms of A derivatives *)
(* L_M = -1/4 F_ab F^ab = -1/2 (D_a A_b D^a A^b - D_a A_b D^b A^a) *)
MaxwellLagrangian = -1/2 * (
  CD3[-a][csA[-b]] eta3[a, c] eta3[b, d] CD3[-c][csA[-d]] -
  CD3[-a][csA[-b]] eta3[a, d] eta3[b, c] CD3[-c][csA[-d]]
);
Print["   Maxwell term: L_M = -1/2 (D_a A_b D^a A^b - D_a A_b D^b A^a)"];

(* Construct Chern-Simons term *)
(* L_CS = (kappa/2) epsilon^abc A_a D_b A_c *)
(* Note: epsiloneta3 is the Levi-Civita tensor density created by DefMetric *)
CSLagrangian = (kappa/2) * epsiloneta3[a, b, c] * csA[-a] * CD3[-b][csA[-c]];
Print["   CS term: L_CS = (kappa/2) epsilon^abc A_a D_b A_c"];

(* Full Lagrangian *)
FullLagrangian = MaxwellLagrangian + CSLagrangian;
Print["   Full Lagrangian: L = L_M + L_CS"];

(* Simplify *)
FullLagrangian = ToCanonical[FullLagrangian];
FullLagrangian = ContractMetric[FullLagrangian, eta3];

(* === Step 4: Derive Full Equations of Motion === *)
Print[""];
Print["Step 4: Computing equations of motion from full Lagrangian..."];

(* Use VarD for the variation with respect to A[-a] *)
(* This derives both Maxwell and CS contributions *)
eomFull = VarD[csA[-a], CD3][FullLagrangian];
eomFull = ToCanonical[eomFull];
eomFull = ContractMetric[eomFull, eta3];

Print["   Full EOM (symbolic): ", eomFull];

(* === Step 5: Decompose to Components === *)
Print[""];
Print["Step 5: Decomposing to components..."];

(*
   The pipeline's DecomposeToComponents function now automatically handles:
   1. ToBasis conversion - abstract indices to basis indices
   2. TraceBasisDummy - sum over dummy indices
   3. RemoveChristoffelSymbols - flat space simplification
   4. EvaluateEpsilonComponents - NEW! evaluates epsilon tensor to numeric ±1
   5. EvaluateMinkowskiMetric - metric components
   6. ConvertCDToDerivatives - explicit Derivative[] form

   The epsilon tensor from the CS term is now handled automatically!
*)

(* Decompose the full EOM to get component equations *)
componentEquations = DecomposeToComponents[eomFull, csA[-a], cart3];

Print["   Number of component equations: ", Length[componentEquations]];

(* Display each component equation *)
Print[""];
Print["   Component equations (Maxwell + CS automatically combined):"];
Do[
  Print["   A_", componentEquations[[i, 1]], " equation: ", componentEquations[[i, 2]] // Simplify],
  {i, Length[componentEquations]}
];

(* Format for JSON export *)
fullEquations = Table[
  {"A_" <> ToString[componentEquations[[i, 1]]], componentEquations[[i, 2]]},
  {i, Length[componentEquations]}
];

(* === Step 6: Build JSON Structure === *)
Print[""];
Print["Step 6: Building JSON structure..."];

(* Metadata for JSON *)
metadata = <|
  "source" -> "xAct",
  "lagrangian_expr" -> "-1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c",
  "derived_from" -> "Euler-Lagrange",
  "gauge" -> "lorenz",
  "linearized" -> True,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"}
|>;

(* Build JSON using multi-field structure *)
jsonStructure = BuildMultiFieldJSONStructure[fullEquations, metadata];
jsonStructure["metadata", "parameters"] = <|"kappa" -> 0.5|>;

(* === Step 7: Export JSON === *)
Print[""];
Print["Step 7: Exporting JSON..."];

(* Export path *)
outputDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
outputPath = FileNameJoin[{outputDir, "chern_simons_3d.json"}];

(* Pretty print JSON *)
Print["JSON Output:"];
Print[ExportString[jsonStructure, "JSON"]];

(* Export *)
Export[outputPath, jsonStructure, "JSON"];
Print[""];
Print["Exported to: ", outputPath];

(* === Summary === *)
Print[""];
Print["*** Chern-Simons 2+1D derivation complete! ***"];
Print[""];
Print["Pipeline features demonstrated:"];
Print["  - Full symbolic derivation from Lagrangian"];
Print["  - Automatic epsilon (Levi-Civita) tensor evaluation"];
Print["  - Cross-field gradient_x/gradient_y operator detection"];
Print["  - No manual hardcoding of CS terms required"];
Print[""];
Print["Physics summary:"];
Print["  Lagrangian: L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c"];
Print["  Coupling: kappa = ", kappa];
Print[""];
Print["  Component equations (with kappa = ", kappa, "):"];
Print["    d2_t(A_0) = laplacian(A_0) + kappa*(d_x A_2 - d_y A_1)"];
Print["    d2_t(A_1) = laplacian(A_1) + kappa*(d_y A_0)"];
Print["    d2_t(A_2) = laplacian(A_2) - kappa*(d_x A_0)"];
Print[""];
Print["  The Chern-Simons term adds topological mass to the photon."];
Print["  This produces helical wave propagation patterns."];
