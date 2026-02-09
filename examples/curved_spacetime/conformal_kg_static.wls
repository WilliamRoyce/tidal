#!/usr/bin/env wolframscript
(* Static Conformal Klein-Gordon Example - Phase 1 of Curved Spacetime *)
(*
   Physics: Klein-Gordon scalar field on conformally flat spacetime with CONSTANT
   conformal factor. This validates metric handling before adding time-dependence.

   Metric: ds^2 = Omega^2 * (-dt^2 + dx^2)  where Omega = 2 (constant)

   For constant Omega:
   - Christoffel symbols = 0 (derivatives of constant metric vanish)
   - Effective mass: m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4  =>  m_eff = 2
   - Equation: d2_t phi = d2_x phi - m_eff^2 phi

   Verification: Static conformal with Omega=2, m=1 should match flat space with m=2
*)

Print["=== Static Conformal Klein-Gordon (Curved Spacetime Phase 1) ==="];
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

(* === Step 1: Setup 1+1D Conformally Flat Spacetime === *)
Print[""];
Print["Step 1: Setting up 1+1D conformally flat spacetime..."];

$DefInfoQ = False;  (* Suppress definition info messages *)

(* Define conformal factor - CONSTANT for Phase 1 *)
Omega0 = 2;  (* Constant conformal factor *)
Print["   Conformal factor: Omega = ", Omega0, " (constant)"];

(* Define manifold *)
If[!xTensorQ[confM2],
  DefManifold[confM2, 2, {a, b, c, d, e, f}]
];

(* Define the conformal metric g_ab = Omega^2 * eta_ab *)
(* For Lorentzian signature with constant Omega: still signature -1 *)
If[!MetricQ[confG],
  DefMetric[-1, confG[-a, -b], confCD, SymbolOfCovD -> {";", "\[Del]"}, PrintAs -> "g"]
];

Print["   Manifold: confM2 (1+1D)"];
Print["   Metric: g_ab = Omega^2 * eta_ab (conformally flat)"];
Print["   Covariant derivative: confCD"];

(* Setup coordinate chart *)
If[!ChartQ[confCart],
  DefChart[confCart, confM2, {0, 1}, {t[], x[]}]
];
Print["   Coordinates: {t, x}"];

(* Set conformal metric components: g_ab = Omega^2 * diag(-1, +1) *)
confMetricComponents = Omega0^2 * DiagonalMatrix[{-1, 1}];
MetricInBasis[confG, -confCart, confMetricComponents];
Print["   Metric components: ", confMetricComponents];
Print["   Inverse metric: g^ab = ", Omega0^(-2), " * diag(-1, +1)"];

(* === Step 2: Define the Scalar Field === *)
Print[""];
Print["Step 2: Defining scalar field phi..."];

If[!xTensorQ[confPhi],
  DefTensor[confPhi[], confM2, PrintAs -> "\[Phi]"]
];

(* Mass parameter (m^2 = 1) *)
m2 = 1;
Print["   Scalar field: phi"];
Print["   Bare mass: m^2 = ", m2];
Print["   Expected effective mass: m_eff^2 = m^2 * Omega^2 = ", m2 * Omega0^2];

(* === Step 3: Construct the Minimally Coupled Lagrangian === *)
Print[""];
Print["Step 3: Constructing minimally coupled Lagrangian..."];

(*
   For minimally coupled scalar in curved spacetime:
   L = sqrt(-g) * (-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2)

   For conformal metric g_ab = Omega^2 eta_ab:
   - sqrt(-g) = Omega^d (d=2 in 1+1D)
   - g^ab = Omega^{-2} eta^ab

   So the Lagrangian density becomes:
   L = Omega^2 * (-1/2 Omega^{-2} eta^ab d_a phi d_b phi - 1/2 m^2 phi^2)
     = -1/2 eta^ab d_a phi d_b phi - 1/2 m^2 Omega^2 phi^2

   The kinetic term is UNCHANGED, but mass term is rescaled by Omega^2.
   For constant Omega, Christoffels = 0, so nabla = partial.
*)

(* Construct using xAct's metric (it handles the curved spacetime correctly) *)
(* L = -1/2 CD[-a][phi] g^ab CD[-b][phi] - 1/2 m^2 phi^2 *)
confLagrangian = -1/2 confCD[-a][confPhi[]] confG[a, b] confCD[-b][confPhi[]] - 1/2 m2 confPhi[]^2;

Print["   L = -1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2"];
Print["   For constant Omega=", Omega0, ":"];
Print["     - Kinetic term unchanged (g^ab CD_a phi CD_b phi = Omega^{-2} eta^ab partial_a phi partial_b phi)"];
Print["     - Mass term scaled: m^2 -> m_eff^2 = m^2 * Omega^2 = ", m2 * Omega0^2];

(* === Step 4: Construct the Wave Equation === *)
Print[""];
Print["Step 4: Constructing wave equation..."];

(*
   For conformal metric g_ab = Omega^2 eta_ab:
   Box_g phi = g^ab nabla_a nabla_b phi

   With constant Omega, Christoffel symbols = 0, so nabla = partial.
   The inverse metric g^ab = Omega^{-2} eta^ab appears in the equation.

   The pipeline handles this via MetricMatrix option:
   - EvaluateMetricComponents substitutes actual g^{ab} values (including Omega^{-2})
   - ExportJSON normalizes by the LHS time-derivative coefficient
   - Result: d2_t phi = d2_x phi - m_eff^2 phi where m_eff^2 = m^2 * Omega^2

   This is the general approach - no need to pre-compute mEffSquared manually.
*)

(* Construct wave equation with bare mass m^2 *)
(* The pipeline will handle the conformal factor via MetricMatrix *)
waveEq = confG[b, c] confCD[-b][confCD[-c][confPhi[]]] - m2 confPhi[];
Print["   Wave equation (abstract): Box_g phi - m^2 phi = 0"];
Print["   Pipeline will normalize: d2_t phi = d2_x phi - m_eff^2 phi"];
Print["   where m_eff^2 = m^2 * Omega^2 = ", m2 * Omega0^2];

(* === Step 5: Extract Components === *)
Print[""];
Print["Step 5: Decomposing to components..."];

componentEqs = DecomposeToComponents[waveEq, confPhi[], confCart, {},
  "MetricMatrix" -> confMetricComponents];
Print["   Component equations: ", Length[componentEqs]];
Print["   Result: ", componentEqs];

(* === Step 6: Build and Export JSON === *)
Print[""];
Print["Step 6: Building JSON structure..."];

metadata = <|
  "name" -> "Conformal Klein-Gordon (Static)",
  "field_prefix" -> "confPhi",
  "lagrangian_expr" -> "-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2",
  "conformal_factor" -> Omega0,
  "bare_mass_squared" -> m2,
  "effective_mass_squared" -> m2 * Omega0^2,
  "gauge" -> "none",
  "linearized" -> True,
  "dimension" -> 2,
  "signature" -> {-1, 1},
  "coordinates" -> {"t", "x"},
  "physics_notes" -> "Static conformal factor Omega=2 gives m_eff^2=4. Should match flat KG with m^2=4."
|>;

(* Build JSON *)
jsonOutput = BuildJSONStructure[componentEqs, metadata];

(* Export *)
outputDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "data"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];
jsonFile = FileNameJoin[{outputDir, "conformal_kg_static.json"}];

jsonString = ExportString[jsonOutput, "JSON"];
Print["JSON Output:"];
Print[jsonString];

Export[jsonFile, jsonOutput, "JSON"];
Print[""];
Print["Exported to: ", jsonFile];

(* === Summary === *)
Print[""];
Print["=== Summary ==="];
Print["Conformal metric: g_ab = ", Omega0^2, " * eta_ab"];
Print["For CONSTANT Omega, Christoffel symbols = 0"];
Print["Klein-Gordon equation: Box_g phi - m^2 phi = 0"];
Print["Pipeline handles: actual metric g^{ab} -> LHS normalization -> operator detection"];
Print["Result: d2_t phi = d2_x phi - m_eff^2 phi"];
Print["   where m_eff^2 = m^2 * Omega^2 = ", m2 * Omega0^2];
Print[""];
Print["Verification: This should match flat Klein-Gordon with m^2 = 4"];
Print[""];
Print["*** Static Conformal Klein-Gordon complete! ***"];
