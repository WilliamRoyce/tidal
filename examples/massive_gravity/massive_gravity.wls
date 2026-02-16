#!/usr/bin/env wolframscript
(* Linearized Massive Gravity - 2+1D with Fierz-Pauli Mass Term *)
(* Derives linearized equations using xPert and a reference metric bg_ab *)
(*
   Physics: Linearized massive gravity in 2+1D flat Minkowski spacetime

   The Fierz-Pauli mass term is the unique ghost-free linear mass term
   for a spin-2 field.  The full nonlinear equation uses a reference
   metric bg_ab (= Minkowski, not perturbed by xPert):

     G_ab - m^2 (g_ab - bg_ab) + m^2 g_ab g^cd (g_cd - bg_cd) = 0

   Perturbation: g_ab = eta_ab + epsilon * h_ab (symmetric rank-2 tensor)
   At zeroth order (g - bg) = 0, so cross terms vanish.
   At first order this linearizes to:

     G^(1)_ab[h] - m^2 (h_ab - eta_ab h) = 0

   where h = eta^cd h_cd is the trace.

   Ghost-freedom: Fierz-Pauli propagates 5 DOF (not 6).  Wrong relative
   coefficient between h^2 and h_ab h^ab activates a Boulware-Deser ghost.

   The trace term creates off-diagonal mass coupling between the diagonal
   metric components (h_0, h_3, h_5).  Off-diagonal components (h_1, h_2,
   h_4) keep simple -m^2 self-coupling.

   In pure 2+1D GR there are no local propagating degrees of freedom
   (the Riemann tensor is fully determined by the Ricci tensor).
   The mass term creates propagating massive modes with:
     - Dispersion relation: omega^2 = k^2 + m^2
     - Group velocity: v_g = k/omega < 1 (subluminal)

   The 6 equations split into:
   - h_0 (h_tt): constraint (no d2_t), driven by h_3, h_5
   - h_1, h_2 (h_tx, h_ty): first-order in time
   - h_3, h_4, h_5 (h_xx, h_xy, h_yy): second-order evolution

   Output: massive_gravity_3d.json (gauge-unfixed, 6 equations with FP mass)
*)

Print["=== Linearized Massive Gravity 2+1D (Fierz-Pauli) ==="];
Print[""];
Print["Field equation: G_ab - m^2 (g_ab - bg_ab) + m^2 g_ab g^cd (g_cd - bg_cd) = 0"];
Print["Linearization:  G^(1)_ab[h] - m^2 (h_ab - eta_ab h) = 0"];
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
Print["Step 2: Defining fields and constants..."];

(* Define plain symmetric tensor for decomposition *)
If[!xTensorQ[mgH],
  DefTensor[mgH[-a, -b], mgM3, Symmetric[{-a, -b}], PrintAs -> "h"]
];

(* Define background/reference metric — not perturbed by xPert *)
If[!xTensorQ[mgBg],
  DefTensor[mgBg[-a, -b], mgM3, Symmetric[{-a, -b}], PrintAs -> "bg"]
];
(* Explicit zero perturbation: bg is non-dynamical *)
Unprotect[Perturbation];
Perturbation[mgBg[__], ___] := 0;
Protect[Perturbation];

(* Define mass constant *)
If[!ConstantSymbolQ[m2], DefConstantSymbol[m2]];

Print["   h_ab: symmetric rank-2 perturbation (6 components in 3D)"];
Print["   bg_ab: background/reference metric (non-dynamical, = Minkowski)"];
Print["   m2: graviton mass squared (DefConstantSymbol)"];

(* ============================================================ *)
(* Step 3: Setup xPert Metric Perturbation *)
(* ============================================================ *)
Print[""];
Print["Step 3: Setting up xPert perturbation..."];

mgHpert = SetupMetricPerturbation[mgEta, mghpert, mgEpsilon];

Print["   Perturbation: g_ab = eta_ab + epsilon * h_ab"];
Print["   Using xPert's DefMetricPerturbation for systematic expansion"];
Print["   bg_ab is NOT registered with xPert (unperturbed by construction)"];

(* ============================================================ *)
(* Step 4: Linearize with Fierz-Pauli mass term *)
(* ============================================================ *)
Print[""];
Print["Step 4: Computing linearized equations..."];
Print["   Linearizing: G_ab - m2 (g_ab - bg_ab) + m2 g_ab g^cd (g_cd - bg_cd)"];

(* The key expression uses a reference metric bg_ab:
   - Einstein[mgCD][-a,-b] -> G^(1)_ab[h]  (linearized Einstein tensor)
   - m2 (mgEta - mgBg) at first order -> m2 * h_ab  (bg unperturbed)
   - m2 mgEta mgEta[c,d] (mgEta - mgBg) at first order -> m2 * eta_ab * h
   Net mass contribution: -m2 h_ab + m2 eta_ab h = -m2 (h_ab - eta_ab h)
   This is exactly the Fierz-Pauli mass term. *)
linExpr = LinearizeTensorExpression[(
  Einstein[mgCD][-a, -b]
  - m2 (mgEta[-a, -b] - mgBg[-a, -b])
  + m2 mgEta[-a, -b] mgEta[c, d] (mgEta[-c, -d] - mgBg[-c, -d])
)];

Print["   Linearized expression: ", Short[linExpr, 3]];

(* Convert xPert notation to plain tensor for pipeline *)
linExprPlain = linExpr /. mghpert[LI[1], idx__] :> mgH[idx];
(* Replace bg with metric — bg is the background (= Minkowski) by construction *)
linExprPlain = linExprPlain /. mgBg -> mgEta;
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
  "name" -> "Linearized Massive Gravity 2+1D (Fierz-Pauli)",
  "field_prefix" -> "h",
  "lagrangian_expr" -> "G_ab - m2 (g_ab - bg_ab) + m2 g_ab g^cd (g_cd - bg_cd)",
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 3,
  "signature" -> {-1, 1, 1},
  "coordinates" -> {"t", "x", "y"},
  "solve_constraints" -> True,
  "constraint_method" -> "auto",
  "constraint_max_iterations" -> 30,
  "constraint_tolerance" -> 1*^-10,
  "constraint_boundary_conditions" -> <|
    "x" -> <|"type" -> "periodic"|>,
    "y" -> <|"type" -> "periodic"|>
  |>
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
Print["  Linearized massive gravity in 2+1D (Fierz-Pauli):"];
Print["    G^(1)_ab - m^2 (h_ab - eta_ab h) = 0"];
Print[""];
Print["  6 component equations:"];
Print["    h_0 (h_tt): constraint, driven by h_3 and h_5 (trace coupling)"];
Print["    h_1, h_2 (h_tx, h_ty): first-order in time (-m^2 self-coupling)"];
Print["    h_3, h_4, h_5 (h_xx, h_xy, h_yy): second-order evolution"];
Print[""];
Print["  Fierz-Pauli mass term -m^2(h_ab - eta_ab h):"];
Print["    - Ghost-free: propagates 5 DOF (not 6)"];
Print["    - Trace h = -h_0 + h_3 + h_5 couples diagonal components"];
Print["    - Off-diagonal h_1, h_2, h_4: simple -m^2 self-coupling"];
Print["    - Dispersion relation: omega^2 = k^2 + m^2"];
Print[""];
Print["  Parameter m2 stored symbolically for runtime override:"];
Print["    tidal simulate massive_gravity_3d.json --param m2=2.0"];
Print[""];
Print["*** Massive Gravity (Fierz-Pauli) derivation complete! ***"];
