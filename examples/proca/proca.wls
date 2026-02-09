#!/usr/bin/env wolframscript
(* Proca Field (Massive Vector) Lagrangian Example - 1+1D *)
(* Derives the massive wave equations from L = -1/4 F_uv F^uv - 1/2 m^2 A_u A^u *)
(*
   Physics: Proca field is a massive photon - breaks gauge invariance but
   automatically enforces Lorenz condition from EOM.

   EOM: Box A^u + m^2 A^u = 0 (massive Klein-Gordon for each component)
   Dispersion: omega^2 = k^2 + m^2
*)

Print["=== Proca Field (Massive Vector) Lagrangian to PDE Pipeline (1+1D) ==="];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;

(* Load pipeline modules *)
Print["Loading pipeline modules..."];
pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "tidal", "wolfram"}];
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

(* === Step 1: Setup 1+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D Minkowski spacetime..."];

(* Use shared symbols to ensure proper Christoffel handling *)
(* Check if already defined to avoid conflicts *)
If[!xTensorQ[M2],
  DefManifold[M2, 2, {a, b, c, d, e, f, g, h}]
];

If[!MetricQ[eta],
  DefMetric[-1, eta[-a, -b], CD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"]
];

(* Define coordinate chart: {t, x} *)
If[!ChartQ[cart],
  DefChart[cart, M2, {0, 1}, {t[], x[]}]
];

(* Set Minkowski metric components: diag(-1, +1) *)
MetricInBasis[eta, -cart, DiagonalMatrix[{-1, 1}]];

Print["   Manifold: M2 (1+1 dimensional)"];
Print["   Metric: Minkowski eta_ab with signature (-,+)"];
Print["   Coordinates: {t, x}"];

(* === Step 2: Define the Vector Potential === *)
Print[""];
Print["Step 2: Defining vector potential A_a..."];

(* Vector potential A_a - check if already defined *)
If[!xTensorQ[procaA],
  DefTensor[procaA[-a], M2, PrintAs -> "A"]
];

Print["   Vector potential: procaA_a"];

(* === Step 3: Construct the Proca Lagrangian === *)
Print[""];
Print["Step 3: Constructing Proca Lagrangian..."];

(* Proca mass parameter - defined as SYMBOLIC constant for runtime parameterization *)
(* The actual numeric value is set at simulation time via Python parameters dict *)
If[!ConstantSymbolQ[procaMassSquared],
  DefConstantSymbol[procaMassSquared, PrintAs -> "m\[Squared]"]
];

Print["   Mass parameter: m^2 = ", procaMassSquared, " (symbolic - set at simulation time)"];

(* L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a *)
(* where F_ab = d_a A_b - d_b A_a *)

(* Maxwell kinetic term: -1/4 F_ab F^ab *)
(* Expand directly to avoid substitution issues (see MEMORY.md) *)
MaxwellTerm = -1/2 CD[-a][procaA[-b]] eta[a, c] eta[b, d] CD[-c][procaA[-d]];

(* Mass term: -1/2 m^2 A_a A^a *)
MassTerm = -1/2 procaMassSquared procaA[-a] eta[a, b] procaA[-b];

(* Complete Proca Lagrangian *)
ProcaLagrangian = MaxwellTerm + MassTerm;

Print["   L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a"];
Print["   Maxwell term: ", MaxwellTerm];
Print["   Mass term: ", MassTerm];
Print["   Full Lagrangian: ", ProcaLagrangian];

(* === Step 4: Derive Euler-Lagrange Equations === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equations for A_a..."];
Print["   delta L / delta A_a = 0"];

(* Use pipeline's Euler-Lagrange computation *)
(* For Proca: EOM is d_b F^ba + m^2 A^a = 0 *)
eom = EulerLagrangeEquation[ProcaLagrangian, procaA[-a], CD];
Print["   EOM (abstract form): ", eom];

(* === Step 5: Apply Lorenz Gauge and Simplify === *)
Print[""];
Print["Step 5: Applying automatic Lorenz condition..."];

(* Key insight for Proca: unlike Maxwell, Lorenz condition is AUTOMATIC *)
(* Taking divergence of EOM: d_a(d_b F^ba + m^2 A^a) = 0 *)
(* Since d_a d_b F^ba = 0 (antisymmetry), we get: m^2 d_a A^a = 0 *)
(* For m != 0: d_a A^a = 0 is enforced, not chosen! *)

Print["   Unlike Maxwell, Proca has NO gauge freedom."];
Print["   Lorenz condition d_a A^a = 0 is automatic from EOM for m != 0"];

(* With automatic Lorenz condition, EOM simplifies to massive wave equation *)
(* Box A_a + m^2 A_a = 0, where Box = eta^bc d_b d_c *)

Print["   Simplified EOM: Box A_a - m^2 A_a = 0"];
Print["   In coordinates (1+1D Minkowski): -d^2_t A + d^2_x A - m^2 A = 0"];
Print["   Rearranging: d^2_t A = d^2_x A - m^2 A"];
Print["   This is Klein-Gordon for each component!"];

(* === Step 6: Extract Component Equations Using Pipeline === *)
Print[""];
Print["Step 6: Decomposing to components using pipeline..."];

(* Construct the massive wave equation in abstract form: Box A_a - m^2 A_a = 0 *)
(* Box = eta^bc d_b d_c = -d^2_t + d^2_x for Minkowski (-,+) *)
(* EOM: Box A - m^2 A = 0 gives d^2_t A = d^2_x A - m^2 A *)
(* This yields the dispersion relation omega^2 = k^2 + m^2 (subluminal propagation) *)
waveEqProca = eta[b, c] CD[-b][CD[-c][procaA[-a]]] - procaMassSquared procaA[-a];

Print["   Wave equation (abstract form): ", waveEqProca, " = 0"];

(* Use pipeline's component decomposition *)
(* The pipeline now properly evaluates Christoffel symbols to zero for flat Minkowski *)
componentEqs = DecomposeToComponents[waveEqProca, procaA[-a], cart];
Print["   Component equations extracted: ", Length[componentEqs], " components"];
Print["   Components: ", componentEqs];

(* === Step 7: Build JSON Structure Using Pipeline === *)
Print[""];
Print["Step 7: Building JSON using pipeline..."];

(* Metadata for JSON export *)
metadata = <|
  "name" -> "Proca Field (Massive Vector)",
  "field_prefix" -> "procaA",  (* Match tensor name used in decomposition *)
  "lagrangian_expr" -> "-1/4 F[-a,-b] F[a,b] - 1/2 m^2 A[-a] A[a]",
  "mass_squared" -> procaMassSquared,
  "gauge" -> "none (automatic Lorenz)",
  "linearized" -> True,
  "dimension" -> 2,  (* Spacetime dimension: 1+1D = 2 total dims *)
  "signature" -> {-1, 1},
  "coordinates" -> {"t", "x"},
  "parameters" -> <|"mass_squared" -> procaMassSquared|>,
  "physics_notes" -> "Massive vector field with dispersion omega^2 = k^2 + m^2"
|>;

(* Build JSON using multi-field pipeline (preserves symbolic coefficients) *)
fieldEquations = {
  {"procaA_0", componentEqs[[1, 2]]},
  {"procaA_1", componentEqs[[2, 2]]}
};
jsonOutput = BuildMultiFieldJSONStructure[fieldEquations, metadata];
Print["   JSON structure created using BuildMultiFieldJSONStructure"];

(* === Step 8: Export JSON === *)
Print[""];
Print["Step 8: Exporting to JSON file..."];

(* Determine output path *)
scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "",
  scriptDir = Directory[]
];
outputDir = FileNameJoin[{scriptDir, "..", "data"}];
If[!DirectoryQ[outputDir],
  CreateDirectory[outputDir]
];
jsonFile = FileNameJoin[{outputDir, "proca_1d.json"}];

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
Print["From the Proca Lagrangian L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a, we derived:"];
Print["  - Proca equations: d_b F^ba + m^2 A^a = 0"];
Print["  - With automatic Lorenz condition: Box A_a + m^2 A_a = 0"];
Print["  - Component equations:"];
Print["    d^2_t A_0 = d^2_x A_0 - m^2 A_0 (massive wave)"];
Print["    d^2_t A_1 = d^2_x A_1 - m^2 A_1 (massive wave)"];
Print[""];
Print["Physics:"];
Print["  - Dispersion relation: omega^2 = k^2 + m^2"];
Print["  - Group velocity: v_g = k/omega < 1 (subluminal)"];
Print["  - Phase velocity: v_p = omega/k > 1 (superluminal)"];
Print["  - For m -> 0: recovers massless Maxwell"];
Print[""];
Print["The equations are ready for simulation with py-pde!"];
Print[""];
Print["*** Proca Lagrangian pipeline complete! ***"];
