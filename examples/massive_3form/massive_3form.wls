#!/usr/bin/env wolframscript
(* Massive 3-Form Field in 3+1D Minkowski Spacetime *)
(* Derives Klein-Gordon equations for each independent component *)
(*
   Physics: Totally antisymmetric rank-3 tensor field C_{abc}
   Lagrangian: L = -1/12 d_a C_{bcd} d^a C^{bcd} - m^2/12 C_{bcd} C^{bcd}
   Symmetry: C_{abc} = -C_{bac} = -C_{acb} = ... (totally antisymmetric)

   In 3+1D Minkowski with signature (-,+,+,+):
     - 64 raw components reduce to 4 independent via antisymmetry
     - Independent components: C_{012}, C_{013}, C_{023}, C_{123}
     - Each satisfies the Klein-Gordon equation: Box C_{abc} - m^2 C_{abc} = 0

   This example demonstrates:
     - Rank-3 tensor support in the pipeline
     - Symmetry reduction via EnumerateComponentTuples (64 -> 4 components)
     - Antisymmetric tensor handling end-to-end

   Output: examples/data/massive_3form.json
*)

Print["=== Massive 3-Form Field in 3+1D ==="];
Print[""];
Print["Lagrangian: L = -1/12 d_a C_bcd d^a C^bcd - m^2/12 C_bcd C^bcd"];
Print["Symmetry: Totally antisymmetric rank-3 tensor"];
Print["Expected: 4 independent components (from 64 raw) via symmetry reduction"];
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

$DefInfoQ = False;

(* ============================================================ *)
(* Step 1: Define 3+1D Minkowski Spacetime *)
(* ============================================================ *)
Print[""];
Print["Step 1: Setting up 3+1D Minkowski spacetime..."];

(* Define 4D manifold (tf prefix to avoid symbol conflicts) *)
If[!xTensorQ[tfM4],
  DefManifold[tfM4, 4, {a, b, c, d, e, f, g, h}]
];

(* Define flat Minkowski metric *)
If[!MetricQ[tfEta],
  DefMetric[-1, tfEta[-a, -b], tfCD,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"]
];

(* Define coordinate chart: {t, x, y, z} *)
If[!ChartQ[tfCart],
  DefChart[tfCart, tfM4, {0, 1, 2, 3}, {t[], x[], y[], z[]}]
];

(* Set flat Minkowski metric components: diag(-1, +1, +1, +1) *)
tfMetricMatrix = DiagonalMatrix[{-1, 1, 1, 1}];
MetricInBasis[tfEta, -tfCart, tfMetricMatrix];

Print["   Manifold: tfM4 (3+1 dimensional)"];
Print["   Metric: eta_ab = diag(-1, +1, +1, +1) (flat Minkowski)"];
Print["   Coordinates: {t, x, y, z}"];

(* ============================================================ *)
(* Step 2: Define Rank-3 Antisymmetric Tensor Field *)
(* ============================================================ *)
Print[""];
Print["Step 2: Defining rank-3 antisymmetric tensor field C..."];

If[!xTensorQ[tfC],
  DefTensor[tfC[-a, -b, -c], tfM4, Antisymmetric[{-a, -b, -c}],
    PrintAs -> "C"]
];

(* Verify symmetry *)
Print["   Tensor: C[-a,-b,-c] on tfM4"];
Print["   Symmetry: Totally antisymmetric (C_{abc} = -C_{bac} etc.)"];
Print["   Symmetry group: ", SymmetryGroupOfTensor[tfC]];

(* Symbolic mass parameter (resolved at simulation time) *)
(* Must use DefConstantSymbol so xAct's VarD and ToCanonical
   recognize m2 as a constant (not a field to vary). *)
If[!ConstantSymbolQ[m2], DefConstantSymbol[m2]];
Print["   Mass parameter: m^2 (symbolic constant, resolved at runtime)"];

(* ============================================================ *)
(* Step 3: Construct the Lagrangian *)
(* ============================================================ *)
Print[""];
Print["Step 3: Constructing Lagrangian..."];

(*
   L = -1/12 d_a C_{bcd} eta^{ah} eta^{be} eta^{bf} eta^{cg} d_h C_{efg}
       - m^2/12 C_{abc} eta^{ae} eta^{bf} eta^{cg} C_{efg}

   IMPORTANT: Use explicit metric tensors (like klein_gordon_3d.wls does)
   rather than automatic index raising. This ensures:
   1. Metric components appear in the EOM for correct sign evaluation
   2. No need for ToCanonical/ContractMetric (which can drop terms)

   The 1/12 = 1/(2*3!) normalization ensures each component gets
   a standard Klein-Gordon equation (without extra numerical factors).
*)
tfLagrangian = (
  -1/12 * tfCD[-d][tfC[-a, -b, -c]] * tfEta[d, h] * tfEta[a, e] * tfEta[b, f] * tfEta[c, g] * tfCD[-h][tfC[-e, -f, -g]]
  - m2/12 * tfC[-a, -b, -c] * tfEta[a, e] * tfEta[b, f] * tfEta[c, g] * tfC[-e, -f, -g]
);

Print["   L = -1/12 d_a C_bcd eta^{ah} eta^{be} eta^{bf} eta^{cg} d_h C_efg"];
Print["       - m^2/12 C_abc eta^{ae} eta^{bf} eta^{cg} C_efg"];
Print["   Symbolic form: ", tfLagrangian];

(* ============================================================ *)
(* Step 4: Derive Euler-Lagrange Equations *)
(* ============================================================ *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equation..."];

(*
   VarD computes the functional derivative delta L / delta C^{abc}.
   xAct handles the antisymmetry of C automatically, including the
   combinatorial factors from the summation convention.
*)
eom = EulerLagrangeEquation[tfLagrangian, tfC[-a, -b, -c], tfCD];

Print["   EOM (abstract): ", eom];
Print["   Expected: Box C_{abc} - m^2 C_{abc} = 0"];

(* ============================================================ *)
(* Step 5: Decompose to Components *)
(* ============================================================ *)
Print[""];
Print["Step 5: Decomposing to coordinate components..."];
Print["   This is where symmetry reduction happens:"];
Print["   EnumerateComponentTuples detects Antisymmetric[{-a,-b,-c}]"];
Print["   64 raw index combinations -> 4 independent components"];
Print[""];

componentEqs = DecomposeToComponents[eom, tfC[-a, -b, -c], tfCart, {}];

Print["   Number of component equations: ", Length[componentEqs]];
Print["   (Expected: 4 from antisymmetric reduction of 4^3 = 64)"];
Print[""];

(* Display component equations *)
Print["   Component equations:"];
Do[
  Print["   C_", componentEqs[[i, 1]], " equation: ",
        componentEqs[[i, 2]] // Simplify],
  {i, Length[componentEqs]}
];

(* Format for JSON export *)
fieldEquations = Table[
  {"C_" <> ToString[componentEqs[[i, 1]]], componentEqs[[i, 2]]},
  {i, Length[componentEqs]}
];

(* ============================================================ *)
(* Step 6: Build and Export JSON *)
(* ============================================================ *)
Print[""];
Print["Step 6: Building and exporting JSON..."];

metadata = <|
  "source" -> "xAct",
  "field_prefix" -> "C",
  "lagrangian_expr" -> "L = -1/12 d_a C_bcd d^a C^bcd - m^2/12 C_bcd C^bcd",
  "derived_from" -> "Euler-Lagrange",
  "gauge" -> "none",
  "linearized" -> False,
  "dimension" -> 4,
  "signature" -> {-1, 1, 1, 1},
  "coordinates" -> {"t", "x", "y", "z"}
|>;

jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];

(* Export *)
outputDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
outputPath = FileNameJoin[{outputDir, "massive_3form.json"}];

Print["JSON Output:"];
Print[ExportString[jsonStructure, "JSON"]];

Export[outputPath, jsonStructure, "JSON"];
Print[""];
Print["Exported to: ", outputPath];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print[""];
Print["*** Massive 3-Form derivation complete! ***"];
Print[""];
Print["Pipeline features demonstrated:"];
Print["  - Rank-3 tensor field definition with Antisymmetric[{-a,-b,-c}]"];
Print["  - Tensor symmetry reduction: 64 raw -> 4 independent components"];
Print["  - EnumerateComponentTuples uses SymmetryGroupOfTensor for rank-3"];
Print["  - Euler-Lagrange equation via VarD for rank-3 tensor"];
Print["  - Full symbolic derivation from Lagrangian (no hardcoded equations)"];
Print[""];
Print["Physics summary:"];
Print["  Lagrangian: L = -1/12 d_a C_bcd d^a C^bcd - m^2/12 C_bcd C^bcd"];
Print["  EOM: Box C_abc = m^2 C_abc (Klein-Gordon for each component)"];
Print[""];
Print["  Independent components (antisymmetric, strictly increasing):"];
Print["    C_0 -> C_{012} (t,x,y)"];
Print["    C_1 -> C_{013} (t,x,z)"];
Print["    C_2 -> C_{023} (t,y,z)"];
Print["    C_3 -> C_{123} (x,y,z)"];
Print[""];
Print["  Each component satisfies: d^2_t C_i = laplacian C_i - m^2 C_i"];
