#!/usr/bin/env wolframscript
(* Linearized Massive Gravity - 2+1D from Einstein - mass term *)
(* Derives linearized equations using xPert: G^(1)_ab - m^2 h_ab = 0 *)
(*
   Physics: Linearized massive gravity in 2+1D flat Minkowski spacetime

   Starting point: Field equation G_ab - m^2 g_ab = 0
   Perturbation: g_ab = eta_ab + epsilon * h_ab (symmetric rank-2 tensor)
   Linearization: G^(1)_ab[h] - m^2 h_ab = 0

   Sign convention: The MINUS on m^2 is essential. With "+", the identity
   terms in the evolution equations get positive coefficient
   (d2_t(h) = +m^2 h + ...), causing exponential growth (tachyonic).
   With "-", we get d2_t(h) = -m^2 h + ... (stable Klein-Gordon).

   In pure 2+1D GR there are no local propagating degrees of freedom
   (the Riemann tensor is fully determined by the Ricci tensor).
   The mass term m^2 h_ab creates a propagating massive mode with:
     - Dispersion relation: omega^2 = k^2 + m^2
     - Group velocity: v_g = k/omega < 1 (subluminal)

   The 6 equations split into:
   - h_0 (h_tt): constraint (no d2_t)
   - h_1, h_2 (h_tx, h_ty): first-order in time
   - h_3, h_4, h_5 (h_xx, h_xy, h_yy): second-order evolution

   Output: massive_gravity_3d.json (gauge-unfixed, 6 equations with mass terms)
*)

Print["=== Linearized Massive Gravity 2+1D ==="];
Print[""];
Print["Field equation: G_ab - m^2 g_ab = 0"];
Print["Linearization: G^(1)_ab[h] - m^2 h_ab = 0"];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;
<< xAct`xPert`;

(* Load pipeline modules *)
Print["Loading pipeline modules..."];
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "tidal", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];
Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];

$DefInfoQ = False;

(* ============================================================ *)
(* Step 1: Define 2+1D Minkowski Background *)
(* ============================================================ *)
Print[""];
Print["Step 1: Setting up 2+1D Minkowski background..."];

(* Define 3D manifold *)
If[!xTensorQ[mgM3],
  DefManifold[mgM3, 3, {a, b, c, d, e, f}]
];

(* Define flat Minkowski metric *)
If[!MetricQ[mgEta],
  DefMetric[-1, mgEta[-a, -b], mgCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"]
];

(* Define coordinate chart: {t, x, y} *)
If[!ChartQ[mgCart],
  DefChart[mgCart, mgM3, {0, 1, 2}, {t[], x[], y[]}]
];

(* Set flat Minkowski metric components: diag(-1, +1, +1) *)
mgMetricMatrix = DiagonalMatrix[{-1, 1, 1}];
MetricInBasis[mgEta, -mgCart, mgMetricMatrix];

Print["   Manifold: mgM3 (2+1 dimensional)"];
Print["   Metric: eta_ab = diag(-1, +1, +1) (flat Minkowski)"];
Print["   Coordinates: {t, x, y}"];

(* ============================================================ *)
(* Step 2: Define Fields and Constants *)
(* ============================================================ *)
Print[""];
Print["Step 2: Defining metric perturbation and mass constant..."];

(* Define plain symmetric tensor for decomposition *)
If[!xTensorQ[mgH],
  DefTensor[mgH[-a, -b], mgM3, Symmetric[{-a, -b}], PrintAs -> "h"]
];

(* Define mass constant *)
If[!ConstantSymbolQ[m2], DefConstantSymbol[m2]];

Print["   h_ab: symmetric rank-2 perturbation (6 components in 3D)"];
Print["   m2: graviton mass squared (DefConstantSymbol)"];

(* ============================================================ *)
(* Step 3: Setup xPert Metric Perturbation *)
(* ============================================================ *)
Print[""];
Print["Step 3: Setting up xPert perturbation..."];

mgHpert = SetupMetricPerturbation[mgEta, mghpert, mgEpsilon];

Print["   Perturbation: g_ab = eta_ab + epsilon * h_ab"];
Print["   Using xPert's DefMetricPerturbation for systematic expansion"];

(* ============================================================ *)
(* Step 4: Linearize G_ab - m^2 g_ab *)
(* ============================================================ *)
Print[""];
Print["Step 4: Computing linearized equations..."];
Print["   Linearizing: Einstein[CD][-a,-b] - m2 * eta[-a,-b]"];

(* The key expression: Einstein tensor - mass term *)
(* LinearizeTensorExpression extracts the O(epsilon) part:
   - Einstein[mgCD][-a,-b] -> G^(1)_ab[h]  (linearized Einstein tensor)
   - -m2 * mgEta[-a,-b] -> -m2 * h_ab  (first-order metric perturbation)
   The minus sign ensures stable evolution: d2_t(h) = -m^2 h + spatial *)
linExpr = LinearizeTensorExpression[Einstein[mgCD][-a, -b] - m2 mgEta[-a, -b]];

Print["   Linearized expression: ", Short[linExpr, 3]];

(* Convert xPert notation to plain tensor for pipeline *)
linExprPlain = linExpr /. mghpert[LI[1], idx__] :> mgH[idx];
linExprPlain = Simplify[linExprPlain];

Print["   Converted to plain tensor: ", Short[linExprPlain, 3]];

(* ============================================================ *)
(* Step 5: Decompose to Component Equations *)
(* ============================================================ *)
Print[""];
Print["Step 5: Decomposing to 6 component equations..."];

componentEqs = DecomposeToComponents[linExprPlain, mgH[-a, -b], mgCart, {},
  "MetricMatrix" -> mgMetricMatrix];

Print["   Generated ", Length[componentEqs], " component equations"];

(* Print component structure *)
Print[""];
Print["   Component structure:"];
Do[
  Module[{idx = eq[[1]], comp = eq[[2]]},
    Print["     h_", idx, ": ", Short[comp, 3]]
  ],
  {eq, componentEqs}
];

(* ============================================================ *)
(* Step 6: Export JSON *)
(* ============================================================ *)
Print[""];
Print["Step 6: Exporting to JSON..."];

(* Build field equations *)
fieldEquations = Table[
  {"h_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]},
  {k, Length[componentEqs]}
];

(* Build metadata *)
metadata = <|
  "name" -> "Linearized Massive Gravity 2+1D",
  "field_prefix" -> "h",
  "lagrangian_expr" -> "Einstein[CD][-a,-b] - m2 * eta[-a,-b]",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"}
|>;

(* Build JSON via pipeline *)
jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];

(* Inject parameter default values *)
jsonStructure = jsonStructure /. ("parameters" -> _) -> ("parameters" -> <|"m2" -> 1.0|>);

Print["   JSON structure created"];

(* Export *)
scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = Directory[]];
outputDir = FileNameJoin[{scriptDir, "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];

jsonFile = FileNameJoin[{outputDir, "massive_gravity_3d.json"}];
jsonString = ExportString[jsonStructure, "JSON"];
Print["   JSON:"];
Print[jsonString];

Export[jsonFile, jsonStructure, "JSON"];
Print["   Exported to: ", jsonFile];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print[""];
Print["=== Summary ==="];
Print[""];
Print["  Linearized massive gravity in 2+1D: G^(1)_ab - m^2 h_ab = 0"];
Print[""];
Print["  6 component equations:"];
Print["    h_0 (h_tt): constraint (no d2_t)"];
Print["    h_1, h_2 (h_tx, h_ty): first-order in time"];
Print["    h_3, h_4, h_5 (h_xx, h_xy, h_yy): second-order evolution"];
Print[""];
Print["  Mass term m^2 creates:"];
Print["    - Propagating mode in 2+1D (pure GR has none)"];
Print["    - Dispersion relation: omega^2 = k^2 + m^2"];
Print["    - Non-zero mass matrix in JSON"];
Print[""];
Print["  Parameter m2 stored symbolically for runtime override:"];
Print["    tidal simulate massive_gravity_3d.json --param m2=2.0"];
Print[""];
Print["*** Massive Gravity derivation complete! ***"];
