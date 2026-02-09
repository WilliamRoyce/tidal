#!/usr/bin/env wolframscript
(* ::Package:: *)

(* scalar_vector_coupling.wls - Stress test: scalar + vector cross-field coupling *)
(* Combines: KG scalar, full Maxwell, Proca mass, Chern-Simons, divergence coupling *)
(*
   Lagrangian:
     L = -1/2 (d_a phi)(d^a phi) - phim2/2 phi^2
         -1/4 F_ab F^ab - Am2/2 A_a A^a
         + kCS/2 epsilon^abc A_a d_b A_c
         + gSV phi (d_a A^a)

   Tests pipeline handling of:
     - Mixed-rank cross-field coupling (scalar + vector)
     - Cross-field first_derivative_t and gradient operators
     - Epsilon tensor + cross-field coupling in same system
     - 4 symbolic constants simultaneously
     - 4x4 mass/coupling matrices
*)

Print["=== Scalar-Vector Coupling 2+1D: Stress Test ==="];
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

$DefInfoQ = False;

(* === Step 1: Define 2+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 2+1D Minkowski spacetime..."];

If[!xTensorQ[svcM3],
  DefManifold[svcM3, 3, {a, b, c, d, e, f, g, h}];
  DefMetric[-1, svcEta[-a, -b], svcCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"];
  DefChart[svcCart, svcM3, {0, 1, 2}, {t[], x[], y[]}];
  MetricInBasis[svcEta, -svcCart, DiagonalMatrix[{-1, 1, 1}]];
  Print["   Defined new manifold svcM3"];
];

Print["   Manifold: svcM3 (3D = 2+1 spacetime)"];
Print["   Metric: Minkowski with signature (-,+,+)"];
Print["   Coordinates: {t, x, y}"];

(* === Step 2: Define Fields and Constants === *)
Print[""];
Print["Step 2: Defining fields and constants..."];

(* Scalar field *)
If[!xTensorQ[svcPhi],
  DefTensor[svcPhi[], svcM3];
  Print["   Defined scalar field phi"];
];

(* Vector field *)
If[!xTensorQ[svcA],
  DefTensor[svcA[-a], svcM3];
  Print["   Defined vector field A[-a]"];
];

(* Symbolic constants — kept symbolic for runtime parameterization *)
If[!ConstantSymbolQ[phim2], DefConstantSymbol[phim2]];
If[!ConstantSymbolQ[Am2], DefConstantSymbol[Am2]];
If[!ConstantSymbolQ[kCS], DefConstantSymbol[kCS]];
If[!ConstantSymbolQ[gSV], DefConstantSymbol[gSV]];
Print["   Constants: phim2, Am2, kCS, gSV (all symbolic)"];

(* === Step 3: Construct the Full Lagrangian === *)
Print[""];
Print["Step 3: Constructing Lagrangian..."];

(* Scalar kinetic + mass *)
ScalarLagrangian = (
  -1/2 svcCD[-a][svcPhi[]] svcEta[a, b] svcCD[-b][svcPhi[]]
  - phim2/2 svcPhi[]^2
);
Print["   Scalar: -1/2 (d phi)^2 - phim2/2 phi^2"];

(* Full Maxwell kinetic: -1/4 F_ab F^ab *)
MaxwellLagrangian = -1/2 * (
  svcCD[-a][svcA[-b]] svcEta[a, c] svcEta[b, d] svcCD[-c][svcA[-d]] -
  svcCD[-a][svcA[-b]] svcEta[a, d] svcEta[b, c] svcCD[-c][svcA[-d]]
);
Print["   Maxwell: -1/4 F_ab F^ab (expanded form)"];

(* Proca mass *)
ProcaMass = -Am2/2 svcA[-a] svcEta[a, b] svcA[-b];
Print["   Proca mass: -Am2/2 A_a A^a"];

(* Chern-Simons topological term *)
CSLagrangian = (kCS/2) * epsilonsvcEta[a, b, c] * svcA[-a] * svcCD[-b][svcA[-c]];
Print["   CS: kCS/2 epsilon^abc A_a d_b A_c"];

(* Divergence coupling: phi * (d_a A^a) *)
CouplingLagrangian = gSV * svcPhi[] * svcEta[a, b] * svcCD[-a][svcA[-b]];
Print["   Coupling: gSV phi (d_a A^a)"];

(* Full Lagrangian *)
FullLagrangian = ScalarLagrangian + MaxwellLagrangian + ProcaMass + CSLagrangian + CouplingLagrangian;
FullLagrangian = ToCanonical[FullLagrangian];
FullLagrangian = ContractMetric[FullLagrangian, svcEta];
Print["   Full Lagrangian assembled and simplified"];

(* === Step 4: Derive Euler-Lagrange Equations === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equations..."];

(* EOM for scalar phi *)
eomPhi = VarD[svcPhi[], svcCD][FullLagrangian];
eomPhi = ToCanonical[eomPhi];
eomPhi = ContractMetric[eomPhi, svcEta];
Print["   EOM phi: ", eomPhi];

(* EOM for vector A *)
eomA = VarD[svcA[-a], svcCD][FullLagrangian];
eomA = ToCanonical[eomA];
eomA = ContractMetric[eomA, svcEta];
Print["   EOM A: ", eomA];

(* === Step 5: Decompose to Components === *)
Print[""];
Print["Step 5: Decomposing to components (cross-field)..."];

(* Decompose phi EOM with A as additional field *)
phiComponents = DecomposeToComponents[eomPhi, svcPhi[], svcCart, {svcA[-a]},
  "MetricMatrix" -> DiagonalMatrix[{-1, 1, 1}]];
Print["   phi components: ", Length[phiComponents]];

(* Decompose A EOM with phi as additional field *)
aComponents = DecomposeToComponents[eomA, svcA[-a], svcCart, {svcPhi[]},
  "MetricMatrix" -> DiagonalMatrix[{-1, 1, 1}]];
Print["   A components: ", Length[aComponents]];

(* === Step 6: Build and Export JSON === *)
Print[""];
Print["Step 6: Building JSON structure..."];

(* Combine field equations for multi-field export *)
fieldEquations = Flatten[{
  Table[{"phi_" <> ToString[phiComponents[[k, 1]]], phiComponents[[k, 2]]},
    {k, Length[phiComponents]}],
  Table[{"A_" <> ToString[aComponents[[k, 1]]], aComponents[[k, 2]]},
    {k, Length[aComponents]}]
}, 1];

(* Metadata *)
metadata = <|
  "source" -> "xAct",
  "lagrangian_expr" -> "-1/2 (d phi)^2 - phim2/2 phi^2 - 1/4 F_ab F^ab - Am2/2 A^2 + kCS/2 eps A dA + gSV phi div(A)",
  "derived_from" -> "Euler-Lagrange",
  "gauge" -> "none",
  "linearized" -> False,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"}
|>;

jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];
jsonStructure["metadata", "parameters"] = <|
  "phim2" -> 1.0, "Am2" -> 0.5, "kCS" -> 0.3, "gSV" -> 0.2
|>;

(* Export *)
outputDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
outputPath = FileNameJoin[{outputDir, "scalar_vector_coupling.json"}];

Print["JSON Output:"];
Print[ExportString[jsonStructure, "JSON"]];

Export[outputPath, jsonStructure, "JSON"];
Print[""];
Print["Exported to: ", outputPath];

Print[""];
Print["*** Scalar-Vector Coupling derivation complete! ***"];
Print[""];
Print["Pipeline features stress-tested:"];
Print["  - Mixed-rank cross-field coupling (scalar + vector)"];
Print["  - Cross-field first_derivative_t and gradient operators"];
Print["  - Epsilon tensor + cross-field coupling in same Lagrangian"];
Print["  - 4 symbolic constants (phim2, Am2, kCS, gSV)"];
Print["  - 4x4 mass/coupling matrices"];
