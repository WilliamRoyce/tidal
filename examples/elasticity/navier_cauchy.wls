#!/usr/bin/env wolframscript
(* Navier-Cauchy Equations: Linear Elasticity in 2+1D *)
(* L = (1/2) rho |du/dt|^2 - (1/2) lambda (div u)^2 - mu u_ij u_ij *)
(* where u_ij = (1/2)(d_i u_j + d_j u_i) is the strain tensor *)
(* Tests: Vector PDEs, anisotropic coefficients, directional Laplacians *)

Print["=== Navier-Cauchy (Linear Elasticity) 2+1D Pipeline ==="];
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

(* === Step 1: Setup 2+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 2+1D Minkowski spacetime..."];

(* Define 3D manifold if not already defined *)
If[!xTensorQ[M3],
  DefManifold[M3, 3, {a, b, c, d, e, f, i, j, k, l, m, n}];
  DefMetric[-1, eta3[-a, -b], CD3,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"];
  DefChart[cart3, M3, {0, 1, 2}, {t[], x[], y[]}];
  MetricInBasis[eta3, -cart3, DiagonalMatrix[{-1, 1, 1}]];
  Print["   Defined new manifold M3"];
  ,
  Print["   Using existing manifold M3"];
];

Print["   Manifold: M3 (3D)"];
Print["   Metric: Minkowski with signature (-,+,+)"];
Print["   Coordinates: {t, x, y}"];

(* === Step 2: Define Displacement as Two Scalar Fields === *)
Print[""];
Print["Step 2: Defining displacement components as scalar fields..."];

(* For elasticity, we treat each component of the displacement vector as a scalar field *)
(* This simplifies the variational calculation while preserving the physics *)
If[!xTensorQ[ux],
  DefTensor[ux[], M3, PrintAs -> "u\[LeftBracketingBar]x"];
];
If[!xTensorQ[uy],
  DefTensor[uy[], M3, PrintAs -> "u\[LeftBracketingBar]y"];
];

(* Material parameters (natural units) *)
rho = 1;      (* density *)
lambda = 1;   (* first Lame parameter *)
mu = 1;       (* shear modulus *)

Print["   Displacement x-component: ux[]"];
Print["   Displacement y-component: uy[]"];
Print["   Parameters: rho = ", rho, ", lambda = ", lambda, ", mu = ", mu];
Print["   Derived: lambda + 2*mu = ", lambda + 2*mu, ", lambda + mu = ", lambda + mu];

(* === Step 3: Construct the Elasticity Lagrangian === *)
Print[""];
Print["Step 3: Constructing elasticity Lagrangian..."];

(*
   The full Lagrangian is:
   L = (1/2) rho |du/dt|^2 - (1/2) lambda (div u)^2 - mu * u_ij * u^ij

   where:
   - |du/dt|^2 = (d_t ux)^2 + (d_t uy)^2
   - div u = d_x ux + d_y uy
   - u_ij = (1/2)(d_i u_j + d_j u_i) is the strain tensor

   In 2D spatial: u_ij u^ij = u_xx^2 + 2*u_xy^2 + u_yy^2
   where u_xx = d_x ux, u_yy = d_y uy, u_xy = (1/2)(d_x uy + d_y ux)
*)

(* Time derivatives: d_t ux, d_t uy *)
dtux = CD3[{0, -cart3}][ux[]];
dtuy = CD3[{0, -cart3}][uy[]];

(* Spatial derivatives *)
dxux = CD3[{1, -cart3}][ux[]];  (* d_x ux *)
dyux = CD3[{2, -cart3}][ux[]];  (* d_y ux *)
dxuy = CD3[{1, -cart3}][uy[]];  (* d_x uy *)
dyuy = CD3[{2, -cart3}][uy[]];  (* d_y uy *)

(* Kinetic energy: (1/2) rho |du/dt|^2 *)
KineticTerm = (rho/2) * (dtux^2 + dtuy^2);

(* Divergence term: (1/2) lambda (div u)^2 *)
divU = dxux + dyuy;
DivergenceTerm = (lambda/2) * divU^2;

(* Strain tensor components *)
uxx = dxux;
uyy = dyuy;
uxy = (1/2) * (dxuy + dyux);

(* Strain energy: mu * u_ij * u^ij = mu * (uxx^2 + 2*uxy^2 + uyy^2) *)
StrainTerm = mu * (uxx^2 + 2 * uxy^2 + uyy^2);

(* Full Lagrangian: L = T - V *)
ElasticityLagrangian = KineticTerm - DivergenceTerm - StrainTerm;

Print["   Kinetic: (1/2) rho [(d_t ux)^2 + (d_t uy)^2]"];
Print["   Divergence: (1/2) lambda (d_x ux + d_y uy)^2"];
Print["   Strain: mu (uxx^2 + 2*uxy^2 + uyy^2)"];
Print["   Full Lagrangian: ", Short[ElasticityLagrangian, 3]];

(* === Step 4: Derive Euler-Lagrange Equations === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equations using VarD..."];

(* Compute EOM for ux *)
eomUx = VarD[ux[], CD3][ElasticityLagrangian];
eomUx = ToCanonical[eomUx];
eomUx = ContractMetric[eomUx, eta3];

(* Compute EOM for uy *)
eomUy = VarD[uy[], CD3][ElasticityLagrangian];
eomUy = ToCanonical[eomUy];
eomUy = ContractMetric[eomUy, eta3];

Print["   Raw EOM for ux: ", eomUx];
Print["   Raw EOM for uy: ", eomUy];

(* === Step 5: Decompose to Components Using Pipeline === *)
Print[""];
Print["Step 5: Decomposing to component equations..."];

(* Decompose each equation, passing the other field for cross-field detection *)
uxComponents = DecomposeToComponents[eomUx, ux[], cart3, {uy[]}];
uyComponents = DecomposeToComponents[eomUy, uy[], cart3, {ux[]}];

Print["   ux component equations: ", Length[uxComponents]];
Print["   ux equation: ", uxComponents[[1, 2]]];
Print["   uy component equations: ", Length[uyComponents]];
Print["   uy equation: ", uyComponents[[1, 2]]];

(* === Step 6: Build JSON Using Pipeline === *)
Print[""];
Print["Step 6: Building JSON for Python simulation..."];

(* Create field equations list for multi-field export *)
(* Note: Field names must match the function heads generated by DecomposeToComponents *)
(* ux[] -> ux0[t,x,y], so field name should be "ux_0" to match *)
(* uy[] -> uy0[t,x,y], so field name should be "uy_0" to match *)
fieldEquations = {
  {"ux_0", uxComponents[[1, 2]]},  (* ux component *)
  {"uy_0", uyComponents[[1, 2]]}   (* uy component *)
};

(* Metadata for JSON export *)
metadata = <|
  "lagrangian_expr" -> "(1/2)rho|du/dt|^2 - (1/2)lambda(div u)^2 - mu*u_ij*u_ij",
  "equation_type" -> "Navier-Cauchy",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"},
  "parameters" -> {"rho" -> rho, "lambda" -> lambda, "mu" -> mu},
  "mass_matrix" -> {{0.0, 0.0}, {0.0, 0.0}},
  "coupling_matrix" -> {{0.0, 0.0}, {0.0, 0.0}}
|>;

(* Build JSON using multi-field pipeline *)
jsonOutput = BuildMultiFieldJSONStructure[fieldEquations, metadata];
Print["   JSON structure created using pipeline"];

(* === Step 7: Export JSON === *)
Print[""];
Print["Step 7: Exporting JSON..."];

outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputPath], CreateDirectory[outputPath]];
jsonFile = FileNameJoin[{outputPath, "navier_cauchy_2d.json"}];

(* Display JSON *)
jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

(* Export to file *)
Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

Print[""];
Print["*** Navier-Cauchy derivation complete! ***"];
Print["Expected EOMs (Navier-Cauchy equations):"];
Print["  rho * d2_t(ux) = (lambda+2*mu) * d2_x(ux) + mu * d2_y(ux) + (lambda+mu) * d_x_d_y(uy)"];
Print["  rho * d2_t(uy) = mu * d2_x(uy) + (lambda+2*mu) * d2_y(uy) + (lambda+mu) * d_x_d_y(ux)"];
Print[""];
Print["With rho=1, lambda=1, mu=1:"];
Print["  d2_t(ux) = 3*laplacian_x(ux) + 1*laplacian_y(ux) + 2*cross_derivative_xy(uy)"];
Print["  d2_t(uy) = 1*laplacian_x(uy) + 3*laplacian_y(uy) + 2*cross_derivative_xy(ux)"];
