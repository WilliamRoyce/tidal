#!/usr/bin/env wolframscript
(* Linearized Gravity - 3+1D Gravitational Waves from Einstein-Hilbert Action *)
(* Derives gauge-unfixed linearized Einstein equations using xPert *)
(*
   Physics: Linearized vacuum Einstein equations in 3+1D flat Minkowski spacetime

   Starting point: Einstein-Hilbert action S = ∫ R √(-g) d⁴x
   Perturbation: g_μν = η_μν + h_μν (symmetric rank-2 tensor)

   The Einstein tensor G_μν = δS_EH/δg^μν, so its first-order perturbation
   gives exactly the linearized EOM from the EH Lagrangian:
     G^(1)_μν[h] = 0  (vacuum)

   This is NOT gauge-fixed. The 10 equations naturally split into:
   - Hamiltonian constraint: G^(1)_{00} = 0 (no ∂²_t h)
   - Momentum constraints:  G^(1)_{0i} = 0 (first-order in time)
   - Evolution equations:   G^(1)_{ij} = 0 (second-order in time, coupled)

   Gauge conditions are provided as OPTIONAL tools:
   - De Donder gauge: simplifies to □h_μν = 0 for each component
   - TT gauge: reduces to 2 physical DOF (h_+, h_×)

   Output:
   - linearized_gravity.json (gauge-unfixed, 10 coupled equations, mixed time orders)
   - linearized_gravity_dedonder.json (de Donder gauge, 10 uncoupled wave equations)
*)

Print["=== Linearized Gravity 3+1D (Gravitational Waves) ==="];
Print[""];
Print["Starting point: Einstein-Hilbert action S = integral R sqrt(-g) d^4x"];
Print["Linearization: g_ab = eta_ab + h_ab via xPert"];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;
<< xAct`xPert`;

(* Load pipeline modules *)
Print["Loading pipeline modules..."];
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "torsion_gertsenshtein", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];
Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];

$DefInfoQ = False;

(* ============================================================ *)
(* Step 1: Define 3+1D Minkowski Background *)
(* ============================================================ *)
Print[""];
Print["Step 1: Setting up 3+1D Minkowski background..."];

(* Define 4D manifold *)
If[!xTensorQ[gwM4],
  DefManifold[gwM4, 4, {a, b, c, d, e, f, g, h}]
];

(* Define flat Minkowski metric *)
If[!MetricQ[gwg],
  DefMetric[-1, gwg[-a, -b], gwCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "g"]
];

(* Define coordinate chart: {t, x, y, z} *)
If[!ChartQ[gwcart],
  DefChart[gwcart, gwM4, {0, 1, 2, 3}, {t[], x[], y[], z[]}]
];

(* Set flat Minkowski metric components: diag(-1, +1, +1, +1) *)
gwMetricMatrix = DiagonalMatrix[{-1, 1, 1, 1}];
MetricInBasis[gwg, -gwcart, gwMetricMatrix];

Print["   Manifold: gwM4 (3+1 dimensional)"];
Print["   Metric: eta_ab = diag(-1, +1, +1, +1) (flat Minkowski)"];
Print["   Coordinates: {t, x, y, z}"];

(* ============================================================ *)
(* Step 2: Define Metric Perturbation via xPert *)
(* ============================================================ *)
Print[""];
Print["Step 2: Setting up metric perturbation..."];

(* Use the modernized Linearize.wl to set up perturbation *)
gwHpert = SetupMetricPerturbation[gwg, gwh, gwEpsilon];

Print["   Perturbation: g_ab = eta_ab + epsilon * h_ab"];
Print["   h_ab is symmetric rank-2 tensor (10 independent components in 4D)"];
Print["   Using xPert's DefMetricPerturbation for systematic expansion"];

(* ============================================================ *)
(* Step 3: Compute Linearized Einstein Equations (GAUGE-UNFIXED) *)
(* ============================================================ *)
Print[""];
Print["Step 3: Computing linearized Einstein equations (gauge-unfixed)..."];
Print["   This is equivalent to varying the second-order EH action:"];
Print["   G^(1)_ab = delta(R_ab - 1/2 g_ab R) / delta(g^cd) = 0"];
Print[""];

(* Compute the linearized Einstein tensor *)
(* Construct with proper abstract indices from the gwM4 manifold *)
(* LinearizeTensorExpression works on any tensor expression — here we use Einstein *)
linEinstein = LinearizeTensorExpression[Einstein[gwCD][-a, -b]];

Print["   Linearized Einstein tensor computed."];
Print["   Expression: ", linEinstein // Short];

(* ============================================================ *)
(* Step 4: Decompose to Component Equations *)
(* ============================================================ *)
Print[""];
Print["Step 4: Decomposing to 10 component equations..."];

(* The linearized Einstein equation is a rank-2 symmetric tensor equation.
   We decompose each (mu,nu) component into a scalar equation.
   The pipeline's DecomposeToComponents handles this for rank-2 tensors. *)

(* First, we need to extract h from the xPert LI[1] notation.
   The perturbation field is gwh[LI[1], -a, -b].
   DecomposeToComponents expects a plain tensor T[-a,-b] without LI indices.
   We define gwH (distinct from xPert's gwh) for the component decomposition.
   The name gwH produces scalar functions gwH0, gwH1, ... which end in "h"
   matching the field_prefix "h" via ExportJSON's StringEndsQ matching. *)

(* Define a plain symmetric tensor for the decomposition step *)
If[!xTensorQ[gwH],
  DefTensor[gwH[-a, -b], gwM4, Symmetric[{-a, -b}], PrintAs -> "h"]
];

(* Replace xPert notation with plain tensor *)
linEinsteinPlain = linEinstein /. gwh[LI[1], idx__] :> gwH[idx];
linEinsteinPlain = Simplify[linEinsteinPlain];

(* NOTE: Background curvature tensors (Riemann, Ricci, RicciScalar) left symbolic
   by xPert are evaluated automatically by the pipeline's EvaluateCurvatureComponents
   function during component decomposition. For a constant metric (flat Minkowski),
   it computes that all curvature vanishes (∂g=0 ⟹ Γ=0 ⟹ R=0).
   We do NOT hardcode this result. *)

Print["   Converted xPert notation to plain tensor for pipeline"];
Print["   Background curvature will be evaluated from metric during decomposition"];
Print["   Simplified expression: ", Short[linEinsteinPlain, 3]];
Print["   Decomposing rank-2 symmetric tensor in 4D..."];

(* Decompose to components *)
componentEqs = DecomposeToComponents[linEinsteinPlain, gwH[-a, -b], gwcart, {}];

Print["   Generated ", Length[componentEqs], " component equations"];

(* Print component structure *)
Print[""];
Print["   Component structure (constraint vs evolution):"];
Do[
  Module[{idx = eq[[1]], comp = eq[[2]]},
    Print["     h_", idx, ": ", Short[comp, 3]]
  ],
  {eq, componentEqs}
];

(* ============================================================ *)
(* Step 5: Export Gauge-Unfixed JSON (PRIMARY output) *)
(* ============================================================ *)
Print[""];
Print["Step 5: Exporting gauge-unfixed equations to JSON..."];

(* Build metadata *)
metadata = <|
  "name" -> "Linearized Gravity 3+1D (Gauge-Unfixed)",
  "field_prefix" -> "h",
  "lagrangian_expr" -> "R * sqrt(-g)  (Einstein-Hilbert action)",
  "gauge" -> "none",
  "construction" -> "derived",  (* Derived from Lagrangian via xPert pipeline *)
  "linearized" -> True,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

(* Build JSON via pipeline *)
jsonOutput = BuildJSONStructure[componentEqs, metadata];

Print["   JSON structure created (gauge: none)"];

(* Export *)
scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = Directory[]];
outputDir = FileNameJoin[{scriptDir, "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];

jsonFile = FileNameJoin[{outputDir, "linearized_gravity.json"}];
jsonString = ExportString[jsonOutput, "JSON"];
Print["   JSON (gauge-unfixed):"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print["   Exported to: ", jsonFile];

(* ============================================================ *)
(* Step 6 (OPTIONAL): De Donder Gauge — Analytical Construction *)
(* ============================================================ *)
Print[""];
Print["Step 6 (OPTIONAL): De Donder gauge (analytical construction)..."];
Print["   De Donder gauge: div(h_ab) = (1/2) grad(trace(h))"];
Print["   This simplifies each component equation to Box[h_ab] = 0"];
Print[""];

(* IMPORTANT: The de Donder gauge result below is constructed from the
   ANALYTICALLY KNOWN simplification of linearized Einstein equations.
   It is NOT derived by the pipeline's symbolic gauge-fixing machinery.

   The physics:
   In de Donder gauge on flat Minkowski vacuum, the linearized Einstein
   equations simplify dramatically:
     G^(1)_ab = -1/2 Box[h_ab] = 0  (for trace-reversed perturbation)

   For the actual perturbation (not trace-reversed), in de Donder gauge:
     Box[h_ab] - 1/2 eta_ab Box[h] = 0

   With the additional condition that h (trace) satisfies Box[h] = 0,
   each component independently satisfies Box[h_ab] = 0.

   This is used for simulation because the gauge-unfixed equations contain
   constraints (not all components have d2_t), making direct time-stepping
   of all 10 components inappropriate without a constraint-aware solver.
   The de Donder gauge produces 10 uncoupled wave equations suitable for
   standard explicit time-stepping (RK4). *)

(* Build de Donder component equations: d2_t(h_i) = laplacian(h_i) *)
dim = 4;
nComponents = dim * (dim + 1) / 2;  (* 10 for symmetric rank-2 in 4D *)
coords = ScalarsOfChart[gwcart];

(* Create 10 scalar functions for the components *)
(* Sign convention: write as -d²_t f + ∇² f = 0 so that the pipeline's
   LHS extraction (lhsCoeff = -1, |lhsCoeff| = 1) produces RHS = +∇²f.
   This matches the VarD convention where time derivatives appear with negative sign. *)
deDonderComponentEqs = Table[
  Module[{fieldSym, waveEq},
    fieldSym = Symbol["gwH" <> ToString[i]];
    (* Wave equation: -d2_t f + d2_x f + d2_y f + d2_z f = 0 *)
    waveEq = -Derivative[2, 0, 0, 0][fieldSym][Sequence @@ coords]
           + Derivative[0, 2, 0, 0][fieldSym][Sequence @@ coords]
           + Derivative[0, 0, 2, 0][fieldSym][Sequence @@ coords]
           + Derivative[0, 0, 0, 2][fieldSym][Sequence @@ coords];
    {i, waveEq}
  ],
  {i, 0, nComponents - 1}
];

Print["   De Donder: each component satisfies Box[h_i] = 0"];
Print["   Generated ", Length[deDonderComponentEqs], " uncoupled wave equations"];

(* Build de Donder JSON *)
deDonderMetadata = <|
  "name" -> "Linearized Gravity 3+1D (De Donder Gauge)",
  "field_prefix" -> "h",
  "lagrangian_expr" -> "R * sqrt(-g) + gauge fixing term",
  "gauge" -> "de_donder",
  "construction" -> "analytical",  (* NOT derived by pipeline — known physics result *)
  "linearized" -> True,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

deDonderJSON = BuildJSONStructure[deDonderComponentEqs, deDonderMetadata];

deDonderFile = FileNameJoin[{outputDir, "linearized_gravity_dedonder.json"}];
deDonderString = ExportString[deDonderJSON, "JSON"];
Print["   De Donder JSON:"];
Print[deDonderString];

Export[deDonderFile, deDonderJSON, "JSON"];
Print["   Exported to: ", deDonderFile];

(* ============================================================ *)
(* Step 7: TT Gauge Discussion *)
(* ============================================================ *)
Print[""];
Print["Step 7: TT (Transverse-Traceless) gauge discussion..."];
Print[""];
Print["   The TT gauge imposes additional conditions on top of de Donder:"];
Print["   1. Tracelessness: eta^{ab} h_ab = h_00 - h_11 - h_22 - h_33 = 0"];
Print["   2. Transversality: h_{0mu} = 0 for all mu"];
Print[""];
Print["   This reduces 10 DOF to 2 physical polarizations:"];
Print["     h_+ : h_{xx} = -h_{yy}  (plus polarization)"];
Print["     h_x : h_{xy} = h_{yx}   (cross polarization)"];
Print[""];
Print["   For a GW propagating along z:"];
Print["     Only h_{xx}, h_{yy}, h_{xy} are non-zero"];
Print["     With h_{xx} = -h_{yy} (traceless)"];
Print[""];
Print["   TT gauge is imposed via initial conditions, not at the Lagrangian level."];
Print["   The pipeline supports this: set TT-compatible initial data in Python."];

ttInfo = TTGaugeConditions[gwh, gwg];
Print["   TT gauge info: ", ttInfo];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print[""];
Print["=== Summary ==="];
Print[""];
Print["From the Einstein-Hilbert action, via xPert linearization:"];
Print[""];
Print["  1. GAUGE-UNFIXED (primary): linearized_gravity.json"];
Print["     - 10 coupled component equations with mixed time orders"];
Print["     - Hamiltonian constraint (h_00): no d2_t (constraint)"];
Print["     - Momentum constraints (h_0i): first-order in time"];
Print["     - Evolution equations (h_ij): second-order, coupled"];
Print["     - Demonstrates full constraint structure of linearized GR"];
Print[""];
Print["  2. DE DONDER GAUGE (optional): linearized_gravity_dedonder.json"];
Print["     - 10 uncoupled wave equations: Box[h_i] = 0"];
Print["     - Each component propagates independently at c = 1"];
Print[""];
Print["  3. TT GAUGE (conceptual): applied via initial conditions"];
Print["     - Reduces to 2 physical DOF: h_+ and h_x polarizations"];
Print["     - Set in Python simulation, not in Lagrangian derivation"];
Print[""];
Print["  Key design choices:"];
Print["    - Primary output (gauge-unfixed) preserves full constraint structure"];
Print["    - Gauge-unfixed equations derived from EH action via xPert (not manual)"];
Print["    - De Donder result is analytically known, constructed for simulation"];
Print["    - TT gauge imposed via initial conditions, not in Lagrangian"];
Print["    - Curvature evaluated from metric definition, never hardcoded"];
Print[""];
Print["*** Linearized Gravity pipeline complete! ***"];
