#!/usr/bin/env wolframscript
(* Coupled Scalar Fields: phi, chi with linear coupling g*phi*chi *)
(* L = 1/2(d phi)^2 - 1/2 m_phi^2 phi^2 + 1/2(d chi)^2 - 1/2 m_chi^2 chi^2 - g phi chi *)
(* Tests: linear coupled systems, mode-mixing, matrix-form PDE generation *)

Print["=== Coupled Scalars Lagrangian to EOM Pipeline ==="];
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

$DefInfoQ = False;  (* Suppress definition info messages *)

(* === Step 1: Setup 1+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D Minkowski spacetime..."];

(* Use same symbol names as other examples to leverage kernel caching *)
(* The symbols M2, eta, CD, cart are shared across examples *)

(* Check if M2 is already defined - if so, skip redefinition *)
If[!xTensorQ[M2],
  DefManifold[M2, 2, {a, b, c, d, e, f}];
  DefMetric[-1, eta[-a, -b], CD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"];
  DefChart[cart, M2, {0, 1}, {t[], x[]}];
  MetricInBasis[eta, -cart, DiagonalMatrix[{-1, 1}]];
  Print["   Defined new manifold M2"];
  ,
  Print["   Using existing manifold M2"];
];

Print["   Manifold: M2 (2D)"];
Print["   Metric: Minkowski with signature (-,+)"];
Print["   Coordinates: {t, x}"];

(* === Step 2: Define Two Scalar Fields === *)
Print[""];
Print["Step 2: Defining scalar fields phi and chi..."];

(* Define scalar fields if not already defined *)
If[!xTensorQ[cplPhi],
  DefTensor[cplPhi[], M2];
];
If[!xTensorQ[cplChi],
  DefTensor[cplChi[], M2];
];

(* Physical parameters (natural units) *)
mPhi2 = 1.0;   (* m_phi^2 = 1 *)
mChi2 = 4.0;   (* m_chi^2 = 4 *)
gCpl = 0.5;    (* coupling constant *)

Print["   Scalar field 1: phi with m_phi^2 = ", mPhi2];
Print["   Scalar field 2: chi with m_chi^2 = ", mChi2];
Print["   Coupling: g = ", gCpl];

(* === Step 3: Construct the Coupled Lagrangian === *)
Print[""];
Print["Step 3: Constructing coupled Lagrangian..."];

(* L = -1/2 (nabla phi)^2 - 1/2 m_phi^2 phi^2 - 1/2 (nabla chi)^2 - 1/2 m_chi^2 chi^2 - g phi chi *)
(* Note: Use explicit + for multiline expressions in Mathematica *)
CoupledLagrangian =
  -1/2 CD[-a][cplPhi[]] eta[a, b] CD[-b][cplPhi[]] - 1/2 mPhi2 cplPhi[]^2 +
  (-1/2 CD[-a][cplChi[]] eta[a, b] CD[-b][cplChi[]] - 1/2 mChi2 cplChi[]^2) +
  (- gCpl cplPhi[] cplChi[]);

Print["   L = -1/2 (\[Del]\[Phi])^2 - 1/2 m\[Phi]^2 \[Phi]^2"];
Print["       -1/2 (\[Del]\[Chi])^2 - 1/2 m\[Chi]^2 \[Chi]^2 - g \[Phi]\[Chi]"];

(* === Step 4: Derive Euler-Lagrange Equations for BOTH Fields === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equations..."];

(* Use VarD from xTensor for functional derivative *)
eomPhiRaw = VarD[cplPhi[], CD][CoupledLagrangian];
eomPhiRaw = ToCanonical[eomPhiRaw];
eomPhiRaw = ContractMetric[eomPhiRaw, eta];
Print["   Raw EOM for phi: ", eomPhiRaw];

eomChiRaw = VarD[cplChi[], CD][CoupledLagrangian];
eomChiRaw = ToCanonical[eomChiRaw];
eomChiRaw = ContractMetric[eomChiRaw, eta];
Print["   Raw EOM for chi: ", eomChiRaw];

(* === Step 5: Decompose to Components Using Pipeline === *)
Print[""];
Print["Step 5: Decomposing to components using pipeline..."];

(* Decompose each field's equation, passing the other field for cross-field transformation *)
(* This ensures that coupling terms like g*phi*chi are properly transformed to coordinates *)
phiComponents = DecomposeToComponents[eomPhiRaw, cplPhi[], cart, {cplChi[]}];
chiComponents = DecomposeToComponents[eomChiRaw, cplChi[], cart, {cplPhi[]}];

Print["   Phi component equations: ", Length[phiComponents], " components"];
Print["   Phi equation: ", phiComponents[[1, 2]]];
Print["   Chi component equations: ", Length[chiComponents], " components"];
Print["   Chi equation: ", chiComponents[[1, 2]]];

(* === Step 6: Build JSON Using Pipeline === *)
Print[""];
Print["Step 6: Building JSON using pipeline..."];

(* Create field equations list for multi-field export *)
(* Format: {{"phi_0", eqPhi}, {"chi_0", eqChi}} *)
fieldEquations = {
  {"phi_0", phiComponents[[1, 2]]},
  {"chi_0", chiComponents[[1, 2]]}
};

(* Metadata for JSON export *)
metadata = <|
  "lagrangian_expr" -> "1/2(d phi)^2 - 1/2 m_phi^2 phi^2 + 1/2(d chi)^2 - 1/2 m_chi^2 chi^2 - g phi chi",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 2,
  "signature" -> {-1, 1},
  "coordinates" -> {"t", "x"},
  "mass_matrix" -> {{N[mPhi2], 0.0}, {0.0, N[mChi2]}},
  "coupling_matrix" -> {{0.0, N[gCpl]}, {N[gCpl], 0.0}}
|>;

(* Build JSON using multi-field pipeline *)
jsonOutput = BuildMultiFieldJSONStructure[fieldEquations, metadata];
Print["   JSON structure created using pipeline"];

(* === Step 7: Export JSON === *)
Print[""];
Print["Step 7: Exporting JSON..."];

outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputPath], CreateDirectory[outputPath]];
jsonFile = FileNameJoin[{outputPath, "coupled_scalars.json"}];

(* Display JSON *)
jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

(* Export to file *)
Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

Print[""];
Print["*** Coupled scalars derivation complete! ***"];
Print["EOMs derived from Lagrangian using pipeline:"];
Print["  d2_t(phi) = laplacian(phi) - m_phi^2 * phi - g * chi"];
Print["  d2_t(chi) = laplacian(chi) - m_chi^2 * chi - g * phi"];
Print[""];
Print["Parameters: m_phi^2 = ", mPhi2, ", m_chi^2 = ", mChi2, ", g = ", gCpl];
