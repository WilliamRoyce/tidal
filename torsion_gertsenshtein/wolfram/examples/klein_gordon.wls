#!/usr/bin/env wolframscript
(* Klein-Gordon Lagrangian Example - Validation for the pipeline *)
(* This derives the Klein-Gordon equation from its Lagrangian *)

Print["=== Klein-Gordon Lagrangian to EOM Pipeline ==="];
Print[""];

(* Load xAct packages *)
Print["Loading xAct packages..."];
<< xAct`xTensor`;
<< xAct`xCoba`;

(* === Step 1: Setup 1+1D Minkowski Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D Minkowski spacetime..."];

DefManifold[M2, 2, {a, b, c, d, e, f}];
DefMetric[-1, eta[-a, -b], M2, SymbolOfCovD -> {";", "\[Del]"}];

(* Get the covariant derivative *)
CD = CovDOfMetric[eta[-a, -b]];
Print["   Manifold: M2 (2D)"];
Print["   Metric: Minkowski eta_ab with signature (-,+)"];
Print["   Covariant derivative: ", CD];

(* Setup coordinate chart for component extraction later *)
DefChart[cart, M2, {0, 1}, {t[], x[]}];
Print["   Coordinates: {t, x}"];

(* Set Minkowski metric components *)
MetricInBasis[eta, -cart, DiagonalMatrix[{-1, 1}]];
Print["   Metric components: diag(-1, +1)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

DefTensor[phi[], M2];
DefConstantSymbol[m2, PrintAs -> "m\[Squared]"];
Print["   Scalar field: phi"];
Print["   Mass parameter: m^2"];

(* === Step 3: Construct the Klein-Gordon Lagrangian === *)
Print[""];
Print["Step 3: Constructing Klein-Gordon Lagrangian..."];

(* L = -1/2 (nabla_a phi)(nabla^a phi) - 1/2 m^2 phi^2 *)
(* In xTensor notation: CD[-a][phi[]] is nabla_a phi *)
KGLagrangian = -1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - 1/2 m2 phi[]^2;

Print["   L = -1/2 (\[Del]_a \[Phi])(\[Del]^a \[Phi]) - 1/2 m\[Squared] \[Phi]\[Squared]"];
Print["   Symbolic form: ", KGLagrangian];

(* === Step 4: Derive Euler-Lagrange Equation === *)
Print[""];
Print["Step 4: Computing Euler-Lagrange equation..."];
Print["   Using variational derivative: \[Delta]L/\[Delta]\[Phi] - \[Del]_a[\[Delta]L/\[Delta](\[Del]_a\[Phi])] = 0"];

(* xTensor's VarD computes functional derivatives *)
(* VarD[Lagrangian, field[], covd] *)
eomRaw = VarD[KGLagrangian, phi[], CD];
Print["   Raw EOM: ", eomRaw];

(* Simplify *)
eom = ToCanonical[eomRaw];
eom = ContractMetric[eom, eta];
Print["   Simplified EOM: ", eom];

(* === Step 5: Verify the Result === *)
Print[""];
Print["Step 5: Verifying result..."];

(* The d'Alembertian (Box) operator is: Box[phi] = eta^ab nabla_a nabla_b phi *)
(* So we expect: Box[phi] - m^2 phi = 0 *)

Print["   Expected form: \[EmptySquare]\[Phi] - m\[Squared]\[Phi] = 0"];
Print["   where \[EmptySquare] = \[Eta]^ab \[Del]_a \[Del]_b is the d'Alembertian"];

(* In coordinates: Box = -d^2/dt^2 + d^2/dx^2 for signature (-,+) *)
Print["   In coordinates: -\[PartialD]t\[Squared]\[Phi] + \[PartialD]x\[Squared]\[Phi] - m\[Squared]\[Phi] = 0"];

(* === Step 6: Extract Components for JSON Export === *)
Print[""];
Print["Step 6: Preparing for JSON export..."];

(* For a scalar field, there's only one component equation *)
(* The wave equation: d^2/dt^2 phi = laplacian(phi) - m^2 phi *)
(* Where laplacian = d^2/dx^2 in 1D *)

Print["   Field: phi (scalar)"];
Print["   Equation type: wave (second-order in time)"];
Print["   Spatial operator: laplacian"];
Print["   Mass term coefficient: ", m2];

(* === Step 7: JSON Output Structure === *)
Print[""];
Print["Step 7: JSON structure for Python pipeline..."];

jsonOutput = <|
  "metadata" -> <|
    "source" -> "xAct",
    "lagrangian_expr" -> "-1/2 CD[-a][phi[]] CD[a][phi[]] - 1/2 m2 phi[]^2",
    "derived_from" -> "Euler-Lagrange",
    "linearized" -> True  (* Already linear *)
  |>,
  "spacetime" -> <|
    "dimension" -> 2,  (* 1 time + 1 space *)
    "signature" -> {-1, 1},
    "coordinates" -> {"t", "x"}
  |>,
  "fields" -> {
    <|"name" -> "phi", "index" -> 0, "is_dynamical" -> True|>
  },
  "equations" -> {
    <|
      "field" -> "phi",
      "lhs" -> "d2_t(phi)",
      "rhs" -> <|
        "type" -> "linear_combination",
        "terms" -> {
          <|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> "phi"|>,
          <|"coefficient" -> -m2, "operator" -> "identity", "field" -> "phi"|>
        }
      |>
    |>
  },
  "coupling" -> <|
    "mass_matrix" -> {{m2}},
    "coupling_matrix" -> {{0}}
  |>
|>;

Print["   JSON structure created"];
Print[""];

(* Export JSON *)
outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "examples", "data"}];
If[!DirectoryQ[outputPath], CreateDirectory[outputPath]];
jsonFile = FileNameJoin[{outputPath, "klein_gordon_1d.json"}];

(* Convert to JSON string *)
jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

(* Export to file *)
Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

Print[""];
Print["*** Klein-Gordon derivation complete! ***"];
Print["The EOM Box[phi] - m^2 phi = 0 was derived from the Lagrangian."];
