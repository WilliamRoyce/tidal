(* ::Package:: *)
(*
   MODULE: ExportJSON.wl
   PURPOSE: Convert component equations to multi-field JSON format for Python pipeline

   DEPENDENCIES:
     - TorsionGertsenshtein`CommonUtilities` (coefficient extraction, field parsing)

   DATA FLOW:
     ComponentEquations (from ComponentDecompose)
       → ParseMultiFieldRHS (identify operators, extract coefficients)
       → EquationToJSONMultiField (build structured equation)
       → BuildMultiFieldJSONStructure (assemble full JSON)
       → ExportEquationSystem (write to file)

   SUPPORTED OPERATORS:
     - identity: field without derivatives (mass/coupling terms)
     - laplacian: full spatial Laplacian ∇²
     - laplacian_x, laplacian_y, laplacian_z: directional second derivatives
     - gradient_x, gradient_y, gradient_z: first-order spatial gradients
     - cross_derivative_xy, cross_derivative_xz, cross_derivative_yz: mixed spatial
     - first_derivative_t: first-order time derivative (Hubble friction in curved spacetime)

   PDE TYPES (via LHS structure):
     - Elliptic (time_order=0): Poisson, Laplace
     - Parabolic (time_order=1): Heat, diffusion
     - Hyperbolic (time_order=2+): Wave equations

   MOMENTUM GRADIENTS:
     Mixed time-space derivatives (∂_t∂_x φ) are converted to spatial gradients
     of momentum fields: gradient_x(pi_i), where pi_i = ∂_t φ_i.

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`ExportJSON`",
  {"TorsionGertsenshtein`CommonUtilities`"}];

(* Public symbols *)
BuildJSONStructure::usage =
  "BuildJSONStructure[componentEqs, metadata] builds the JSON Association \
structure from component equations and metadata.";

EquationToJSON::usage =
  "EquationToJSON[componentEq, fieldName, fieldIndex] converts a single \
component equation to JSON format.";

BuildMultiFieldJSONStructure::usage =
  "BuildMultiFieldJSONStructure[fieldEquations, metadata] builds JSON for systems \
with multiple independent fields. fieldEquations is a list of {fieldName, equation} pairs.";

ParseFieldName::usage =
  "ParseFieldName[name] parses a field name and returns an Association with \
\"base\", \"index\", and \"format\" keys. Supports multiple formats: \
\"standard\" (A_0), \"tensor\" (stress_xy_0), \"compact\" (phi0), \"simple\" (phi).";

$FieldNameFormats::usage =
  "$FieldNameFormats defines the supported field name format patterns: \
\"standard\" (base_index), \"tensor\" (base_component_index), \
\"compact\" (base+digits), \"simple\" (letters only).";

GetTimeDerivativeOrder::usage =
  "GetTimeDerivativeOrder[term] returns the order of time derivative in a term. \
Returns 0 if no time derivative, 1 for d_t, 2 for d2_t, etc. Time index is always \
the first slot in Derivative[dt, dx, ...].";

DetectLHSTimeOrder::usage =
  "DetectLHSTimeOrder[equation] detects the maximum time derivative order in the equation. \
Used to determine if PDE is elliptic (0), parabolic (1), or hyperbolic (2+).";

BuildLHSStructure::usage =
  "BuildLHSStructure[fieldName, timeOrder] builds a structured LHS representation. \
Returns an Association with \"expression\" and \"order\" keys for flexible PDE types.";

(* Time derivative detection helpers *)
ContainsTimeDerivative::usage =
  "ContainsTimeDerivative[term, minOrder:2] returns True if term contains a time \
derivative of at least minOrder. Time is always the first slot in Derivative[dt, dx, ...]. \
Supports 1+1D (2-arg), 2+1D (3-arg), and 3+1D (4-arg) Derivative forms.";

IsMixedTimeSpaceDerivative::usage =
  "IsMixedTimeSpaceDerivative[term] returns True if term has BOTH time AND space \
derivatives (e.g., Derivative[1,1,0] for d_t d_x). Used to detect terms that should \
be converted to momentum gradients. Supports 1+1D, 2+1D, and 3+1D.";

ExtractLHSCoefficient::usage =
  "ExtractLHSCoefficient[term] extracts the scalar coefficient multiplying a \
Derivative[...][f][...] subexpression. Works for any dimension and any derivative \
pattern. Returns coeff such that term = coeff * Derivative[...][f][...]. \
Used by LHS normalization to divide out non-unit time-derivative coefficients.";

IdentifyGradientDirection::usage =
  "IdentifyGradientDirection[term] returns the gradient direction (\"gradient_x\", \
\"gradient_y\", or \"gradient_z\") based on derivative structure. Time slot must be 0 \
for pure spatial gradients. Defaults to \"gradient_x\" if detection fails.";

(* Spatial derivative classification *)
ExtractSpatialGradientFromMixed::usage =
  "ExtractSpatialGradientFromMixed[term] extracts the spatial gradient direction \
from a mixed time-space derivative term. Returns \"gradient_x\", \"gradient_y\", or \
\"gradient_z\". Warns if term has multiple spatial derivatives (e.g., d_t d_x d_y).";

IsSpatialCrossDerivative::usage =
  "IsSpatialCrossDerivative[term] returns True if term contains a spatial cross-derivative \
(d_x d_y, d_x d_z, or d_y d_z with time slot = 0). Legacy boolean wrapper around \
IdentifySpatialCrossDerivative.";

IdentifySpatialCrossDerivative::usage =
  "IdentifySpatialCrossDerivative[term] identifies spatial cross-derivatives and returns \
the operator name: \"cross_derivative_xy\", \"cross_derivative_xz\", \"cross_derivative_yz\", \
or False if not a cross-derivative. Requires time slot = 0.";

IdentifyDirectionalLaplacian::usage =
  "IdentifyDirectionalLaplacian[term] identifies pure second derivatives in a single \
spatial direction. Returns \"laplacian_x\", \"laplacian_y\", \"laplacian_z\", or False. \
Used for anisotropic equations like Navier-Cauchy where d^2_x and d^2_y have different coefficients.";

(* First-order time derivative detection for curved spacetime *)
IsFirstOrderTimeDerivative::usage =
  "IsFirstOrderTimeDerivative[term] returns True if term contains exactly a first-order \
time derivative (no spatial derivatives). Used to detect Hubble friction terms like \
-2H∂_t φ in curved spacetime. Time index is always the first slot in Derivative[dt, dx, ...].";

(* Term identification and building *)
BuildTermResult::usage =
  "BuildTermResult[coeff, operator, field, symbolicCoeff:None, timeDependent:False] builds a term result \
Association with keys \"coefficient\", \"operator\", \"field\". Optionally includes \
\"coefficient_symbolic\" when symbolicCoeff is not None, and \"time_dependent\" when timeDependent is True.";

IdentifyMultiFieldTerm::usage =
  "IdentifyMultiFieldTerm[term, currentFieldName, allFieldNames] analyzes a \
single term and returns an Association with \"coefficient\", \"operator\", and \"field\" keys. \
Detects: laplacian, directional laplacians, gradients, cross-derivatives, identity, and \
momentum gradient terms from mixed time-space derivatives.";

(* Refactored helper functions (Phase 3, Issue 10) *)
ExtractFunctionHeads::usage =
  "ExtractFunctionHeads[term] extracts all function head names from an expression. \
Returns a list of strings representing function names (e.g., {\"phi0\", \"csA1\"}).";

MatchFieldToHeads::usage =
  "MatchFieldToHeads[functionHeads, allFieldNames, defaultField] matches function heads \
to known field names using case-insensitive and prefix matching. Returns an Association \
with \"field\" and \"head\" keys, or uses defaultField if no match found.";

ExtractTermCoefficient::usage =
  "ExtractTermCoefficient[term, fieldHead, targetField] extracts the numeric coefficient \
and optional symbolic coefficient from a term. Returns {numericCoeff, symbolicCoeff, isTimeDependent}.";

ExtractTermCoefficient::symbolic =
  "Symbolic (non-numeric) coefficient `1` found for field `2`. Storing as symbolic coefficient in JSON.";

ClassifyOperatorType::usage =
  "ClassifyOperatorType[term] classifies the differential operator in a term. Returns \
a list {operatorName, isMixedTimeSpace} where operatorName is one of: \"identity\", \
\"laplacian\", \"laplacian_x/y/z\", \"gradient_x/y/z\", \"cross_derivative_xy/xz/yz\".";

Begin["`Private`"];

(* === JSON Structure Building === *)

(* Single-field BuildJSONStructure delegates to multi-field version *)
(* Converts {index, equation} pairs to {fieldName, equation} pairs *)
BuildJSONStructure[componentEqs_, metadata_Association] := Module[
  {fieldEquations, prefix},
  prefix = metadata["field_prefix"];
  fieldEquations = Table[
    {prefix <> "_" <> ToString[eq[[1]]], eq[[2]]},
    {eq, componentEqs}
  ];
  BuildMultiFieldJSONStructure[fieldEquations, metadata]
];

(* === Multi-Field JSON Structure Building === *)

BuildMultiFieldJSONStructure[fieldEquations_List, metadata_Association] := Module[
  {json, fields, equations, allFieldNames, nFields},

  (* fieldEquations format: {{"phi_0", eqPhi}, {"chi_0", eqChi}, ...} *)
  allFieldNames = fieldEquations[[All, 1]];
  nFields = Length[fieldEquations];

  (* Build fields list *)
  fields = Table[
    <|
      "name" -> fieldEquations[[i, 1]],
      "index" -> i - 1,
      "is_dynamical" -> True
    |>,
    {i, nFields}
  ];

  (* Convert equations to JSON format *)
  (* Pass allFieldNames so cross-field references can be detected *)
  equations = Table[
    EquationToJSONMultiField[
      fieldEquations[[i, 2]],
      fieldEquations[[i, 1]],
      i - 1,
      allFieldNames,
      metadata
    ],
    {i, nFields}
  ];

  (* Build full JSON structure *)
  json = <|
    "metadata" -> <|
      "source" -> "xAct",
      "lagrangian_expr" -> Lookup[metadata, "lagrangian_expr", ""],
      "derived_from" -> "Euler-Lagrange",
      "gauge" -> Lookup[metadata, "gauge", "none"],
      "linearized" -> Lookup[metadata, "linearized", True]
    |>,
    "spacetime" -> <|
      "dimension" -> Lookup[metadata, "dimension", 2],
      "signature" -> Lookup[metadata, "signature", {-1, 1}],
      "coordinates" -> Lookup[metadata, "coordinates", {"t", "x"}]
    |>,
    "fields" -> fields,
    "equations" -> equations,
    "coupling" -> <|
      "mass_matrix" -> Lookup[metadata, "mass_matrix",
        ConstantArray[0.0, {nFields, nFields}]],
      "coupling_matrix" -> Lookup[metadata, "coupling_matrix",
        ConstantArray[0.0, {nFields, nFields}]]
    |>
  |>;

  json
];

(* Equation conversion for multi-field systems *)
(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
EquationToJSONMultiField[componentEq_, fieldName_, fieldIndex_, allFieldNames_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm, lhsTimeOrder, lhsStructure},

  (* Same logic as EquationToJSON but with cross-field awareness *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Detect the time derivative order to determine PDE type *)
  (* elliptic (0), parabolic (1), hyperbolic (2+) *)
  lhsTimeOrder = DetectLHSTimeOrder[componentEq];

  (* For backward compatibility, default to order 2 (hyperbolic) if detection fails *)
  If[lhsTimeOrder == 0 && !FreeQ[componentEq, Derivative],
    (* Equation has derivatives but no time derivatives - likely constraint or elliptic *)
    (* Keep as 0 for elliptic PDEs *)
    Null,
    (* Use detected order for parabolic/hyperbolic *)
    Null
  ];

  (* Use dimension-agnostic helpers for time derivative detection *)
  (* Works for both 1+1D (2-arg Derivative) and 2+1D (3-arg Derivative) *)
  timeDerivTerm = Select[terms, ContainsTimeDerivative[#, lhsTimeOrder] &];
  (* RHS = everything EXCEPT pure time derivatives of detected order *)
  (* Mixed time-space derivatives ARE included - they get converted to momentum gradients *)
  (* by IdentifyMultiFieldTerm (e.g., d_t d_x A -> gradient_x(pi)) *)
  rhs = Total[Select[terms, !ContainsTimeDerivative[#, lhsTimeOrder] &]];

  (* LHS normalization: extract time-derivative coefficient and normalize RHS *)
  (* For |lhsCoeff| = 1: rhs = non-time terms as-is (handles both VarD and direct construction) *)
  (* For non-unit lhsCoeff (curved spacetime): rhs = -non_time_terms / lhsCoeff *)
  (*   Example: lhsCoeff = -Omega^{-2} gives rhs = Omega^2 * non_time_terms *)
  Module[{lhsCoeff},
    lhsCoeff = If[Length[timeDerivTerm] > 0,
      ExtractLHSCoefficient[timeDerivTerm[[1]]],
      -1
    ];
    If[Abs[lhsCoeff] =!= 1,
      rhs = Simplify[-rhs / lhsCoeff]
    ]
  ];

  (* Parse RHS with cross-field detection *)
  rhsTerms = ParseMultiFieldRHS[rhs, fieldName, allFieldNames];

  (* Build structured LHS for flexible PDE types *)
  lhsStructure = BuildLHSStructure[fieldName, lhsTimeOrder];

  <|
    "field" -> fieldName,
    "lhs" -> lhsStructure,  (* Now structured: {"expression": "...", "order": {...}} *)
    "rhs" -> <|
      "type" -> "linear_combination",
      "terms" -> rhsTerms
    |>
  |>
];

(* Parse RHS with cross-field reference detection *)
ParseMultiFieldRHS[eq_, currentFieldName_, allFieldNames_] := Module[
  {terms, parsedTerms},

  terms = If[Head[eq] === Plus, List @@ eq, {eq}];
  parsedTerms = Map[IdentifyMultiFieldTerm[#, currentFieldName, allFieldNames] &, terms];
  parsedTerms = DeleteCases[parsedTerms, Nothing];

  (* Fail explicitly if no terms were parsed - do not inject ghost equations *)
  If[Length[parsedTerms] === 0,
    Throw[StringJoin[
      "ParseMultiFieldRHS: No terms parsed from RHS for field '", currentFieldName, "'. ",
      "RHS expression: ", ToString[eq], ". ",
      "Check that the equation is properly expanded and contains recognizable operators."
    ]]
  ];

  parsedTerms
];

(* === LHS Structure Detection (Phase 2, Issue 6) === *)
(* Supports parabolic (heat: d_t φ), elliptic (Poisson: ∇²φ = f), and hyperbolic (wave: d2_t φ) PDEs *)

(* Get time derivative order from a single term *)
(* Time index is always the first slot in Derivative[dt, dx, ...] *)
GetTimeDerivativeOrder[term_] := Module[{orders},
  (* Extract all time derivative orders (first slot of Derivative) *)
  orders = Cases[term, Derivative[n_, ___][_][___] :> n, {0, Infinity}];
  If[Length[orders] == 0, 0, Max[orders]]
];

(* Detect the maximum time derivative order in an equation *)
(* This determines the PDE type: elliptic (0), parabolic (1), hyperbolic (2+) *)
DetectLHSTimeOrder[equation_] := Module[{terms, maxOrder},
  terms = If[Head[equation] === Plus, List @@ equation, {equation}];
  maxOrder = Max[Map[GetTimeDerivativeOrder, terms]];
  (* Return at least 0 *)
  Max[maxOrder, 0]
];

(* Build structured LHS representation *)
(* Supports variable time derivative orders for different PDE types *)
BuildLHSStructure[fieldName_String, timeOrder_Integer] := <|
  "expression" -> Switch[timeOrder,
    0, fieldName,                                    (* Elliptic: no time derivative *)
    1, "d_t(" <> fieldName <> ")",                   (* Parabolic: first order *)
    2, "d2_t(" <> fieldName <> ")",                  (* Hyperbolic: second order *)
    _, "d" <> ToString[timeOrder] <> "_t(" <> fieldName <> ")"  (* Higher order *)
  ],
  "order" -> <|"time" -> timeOrder, "space" -> 0|>
|>;

(* === Time Derivative Detection Helpers (for LHS/RHS separation) === *)

(* Detect if a term contains time derivatives of at least order minOrder *)
(* Time index is always the first slot in Derivative[dt, dx, ...] *)
(* Supports: 1+1D (2-arg), 2+1D (3-arg), 3+1D (4-arg), and higher dimensions *)
ContainsTimeDerivative[term_, minOrder_:2] := Module[{},
  Which[
    (* 3+1D: Derivative[n, m, p, q] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _, _, _][_][___] /; n >= minOrder], True,
    (* 2+1D: Derivative[n, m, p] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _, _][_][___] /; n >= minOrder], True,
    (* 1+1D: Derivative[n, m] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _][_][___] /; n >= minOrder], True,
    (* Generic: any arity with first slot >= minOrder *)
    !FreeQ[term, Derivative[n_, ___][_][___] /; n >= minOrder], True,
    (* Default: no time derivative of required order *)
    True, False
  ]
];

(* Extract the coefficient of the LHS time derivative term *)
(* For flat spacetime: term = -Derivative[2,0][f][t,x], returns -1 *)
(* For curved spacetime: term = -Exp[-2Ht]*Derivative[2,0,0][f][t,x,y], returns -Exp[-2Ht] *)
(* This generalizes the ad-hoc sign-flip logic to handle arbitrary LHS coefficients *)
ExtractLHSCoefficient[term_] := Module[{derivParts, derivPart},
  derivParts = Cases[{term}, Derivative[__][_][__], {0, Infinity}];
  If[Length[derivParts] == 0,
    Throw[StringJoin[
      "ExtractLHSCoefficient: No Derivative found in LHS term '",
      ToString[term, InputForm], "'. Cannot extract time-derivative coefficient."
    ]]
  ];
  derivPart = derivParts[[1]];
  If[term === derivPart, Return[1]];
  Simplify[term / derivPart]
];

(* Check for mixed time-space derivatives that shouldn't be on RHS *)
(* Returns True if term has BOTH time AND space derivatives *)
(* Supports: 1+1D (2-arg), 2+1D (3-arg), 3+1D (4-arg) *)
IsMixedTimeSpaceDerivative[term_] := Module[{},
  Which[
    (* 3+1D: first slot > 0 AND (any spatial slot > 0) *)
    !FreeQ[term, Derivative[n_, m_, p_, q_][_][___] /; n > 0 && (m > 0 || p > 0 || q > 0)], True,
    (* 2+1D: first slot > 0 AND (second OR third slot > 0) *)
    !FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && (m > 0 || p > 0)], True,
    (* 1+1D: first slot > 0 AND second slot > 0 *)
    !FreeQ[term, Derivative[n_, m_][_][___] /; n > 0 && m > 0], True,
    (* Default: not a mixed derivative *)
    True, False
  ]
];

(* Identify gradient direction from derivative structure *)
(* In 1+1D: Derivative[dt, dx] - second slot is x *)
(* In 2+1D: Derivative[dt, dx, dy] - second slot is x, third is y *)
(* In 3+1D: Derivative[dt, dx, dy, dz] - second=x, third=y, fourth=z *)
(* IMPORTANT: First slot (time) must be 0 for pure spatial gradients *)
IdentifyGradientDirection[term_] := Module[{},
  Which[
    (* 3+1D: Check for gradient_z first (fourth slot nonzero, others zero) *)
    (* Pattern: Derivative[0, 0, 0, n_] where n > 0 *)
    !FreeQ[term, Derivative[0, 0, 0, n_][_][___] /; n > 0], "gradient_z",

    (* 3+1D: gradient_y (third slot nonzero, fourth zero) *)
    (* Pattern: Derivative[0, 0, n_, 0] where n > 0 *)
    !FreeQ[term, Derivative[0, 0, n_, 0][_][___] /; n > 0], "gradient_y",

    (* 3+1D: gradient_x (second slot nonzero, third and fourth zero) *)
    (* Pattern: Derivative[0, n_, 0, 0] where n > 0 *)
    !FreeQ[term, Derivative[0, n_, 0, 0][_][___] /; n > 0], "gradient_x",

    (* 2+1D: Check for gradient_y first (third slot nonzero) *)
    (* Pattern: Derivative[0, 0, n_] where n > 0 - MUST have first slot = 0 *)
    (* This ensures we only match pure spatial gradients, not mixed time-space *)
    !FreeQ[term, Derivative[0, 0, n_][_][___] /; n > 0], "gradient_y",

    (* 2+1D: gradient_x (second slot nonzero, third is zero) *)
    (* Pattern: Derivative[0, n_, 0] where n > 0 - MUST have first slot = 0 *)
    !FreeQ[term, Derivative[0, n_, 0][_][___] /; n > 0], "gradient_x",

    (* 1+1D: only gradient_x exists (second slot) *)
    (* Pattern: Derivative[0, n_] where n > 0 - MUST have first slot = 0 *)
    !FreeQ[term, Derivative[0, n_][_][___] /; n > 0], "gradient_x",

    (* Fail explicitly if gradient direction cannot be determined *)
    True,
      Throw[StringJoin[
        "IdentifyGradientDirection: Cannot determine gradient direction for term '",
        ToString[term], "'. ",
        "Expected Derivative[0, n_, ...] pattern for pure spatial gradient."
      ]]
  ]
];

(* === First-Order Time Derivative Detection (Phase 2 Curved Spacetime) === *)

(* Detect pure first-order time derivative: Derivative[1, 0, ...] with no spatial derivatives *)
(* Used for Hubble friction terms in curved spacetime: -2H∂_t φ *)
(* Returns True only if time order = 1 AND all spatial orders = 0 *)
IsFirstOrderTimeDerivative[term_] := Module[{},
  Which[
    (* 3+1D: Derivative[1, 0, 0, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0, 0, 0][_][___]], True,
    (* 2+1D: Derivative[1, 0, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0, 0][_][___]], True,
    (* 1+1D: Derivative[1, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0][_][___]], True,
    (* Default: not a pure first-order time derivative *)
    True, False
  ]
];

(* === Phase 4: Momentum Gradient Helpers === *)

(* === Flexible Field Name Parsing === *)
(* Supports multiple field naming conventions *)

(* Field name format patterns *)
$FieldNameFormats = <|
  "standard" -> RegularExpression["^([a-zA-Z]+)_([0-9]+)$"],       (* A_0, phi_1 *)
  "tensor" -> RegularExpression["^(.+)_([0-9]+)$"],                (* stress_xy_0, u_x_1 *)
  "compact" -> RegularExpression["^([a-zA-Z]+)([0-9]+)$"],         (* phi0, A1 *)
  "simple" -> RegularExpression["^[a-zA-Z]+$"]                     (* phi, psi *)
|>;

(* Parse field name into base, index, and format *)
ParseFieldName[name_String] := Module[{match},
  (* Try standard format: A_0, phi_1 *)
  match = StringCases[name, $FieldNameFormats["standard"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "standard"|>]
  ];

  (* Try tensor format: stress_xy_0, u_x_1 (greedy match for base) *)
  match = StringCases[name, $FieldNameFormats["tensor"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "tensor"|>]
  ];

  (* Try compact format: phi0, A1 *)
  match = StringCases[name, $FieldNameFormats["compact"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "compact"|>]
  ];

  (* Try simple format: phi, psi (single field, no index) *)
  If[StringMatchQ[name, $FieldNameFormats["simple"]],
    Return[<|"base" -> name, "index" -> 0, "format" -> "simple"|>]
  ];

  (* No format matched - fail explicitly *)
  Throw[StringJoin[
    "ParseFieldName: Cannot parse field name '", name, "'. ",
    "Expected one of: standard (A_0), tensor (stress_xy_0), compact (phi0), or simple (phi)."
  ]]
];

(* Map field name to momentum field name *)
(* For py-pde state [field_0, pi_0, field_1, pi_1, ...], momentum is at odd indices *)
(* Uses ParseFieldName for flexible format support *)
FieldToMomentumName[fieldName_String] := Module[{parsed},
  parsed = ParseFieldName[fieldName];
  "pi_" <> ToString[parsed["index"]]
];

(* Extract spatial gradient direction from mixed time-space derivative *)
(* Derivative[1, 1, 0] (d_t d_x) -> "gradient_x" *)
(* Derivative[1, 0, 1] (d_t d_y) -> "gradient_y" *)
(* Derivative[1, 0, 0, 1] (d_t d_z) -> "gradient_z" *)
(* These represent d_x(pi) or d_y(pi) or d_z(pi) since d_t phi = pi *)
(*
   WARNING: Mixed time + cross-spatial derivatives like d_t d_x d_y or d_t d_x d_z
   cannot be simply represented as gradient(pi). They require cross_derivative(pi).
   For now, we default to gradient_x with a warning.
*)
ExtractSpatialGradientFromMixed[term_] := Module[{},
  (* Check for problematic mixed time + multiple spatial derivatives *)
  If[!FreeQ[term, Derivative[n_, m_, p_, q_][_][___] /; n > 0 && Count[{m, p, q}, _?(# > 0 &)] > 1],
    Print["Warning: Mixed time + cross-spatial derivative in 3+1D detected. ",
          "Cannot exactly represent as momentum gradient. Using gradient_x(pi)."];
  ];
  If[!FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && m > 0 && p > 0],
    Print["Warning: Mixed time + cross-spatial derivative detected (e.g., d_t d_x d_y). ",
          "This term cannot be exactly represented as a simple momentum gradient. ",
          "Approximating as gradient_x(pi). Consider gauge-fixing to eliminate such terms."];
  ];

  Which[
    (* 3+1D: z-only gradient (fourth slot > 0, second and third = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, 0, m_][_][___] /; n > 0 && m > 0], "gradient_z",
    (* 3+1D: y-only gradient (third slot > 0, second and fourth = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, m_, 0][_][___] /; n > 0 && m > 0], "gradient_y",
    (* 3+1D: x-only gradient (second slot > 0, third and fourth = 0, with time) *)
    !FreeQ[term, Derivative[n_, m_, 0, 0][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 2+1D: Check y-only gradient first (third slot > 0, second slot = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, m_][_][___] /; n > 0 && m > 0], "gradient_y",
    (* 2+1D: x-only gradient (second slot > 0, third slot = 0, with time) *)
    !FreeQ[term, Derivative[n_, m_, 0][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 2+1D: Both x and y present - default to gradient_x (warning above) *)
    !FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && m > 0 && p > 0], "gradient_x",
    (* Fallback for other 2+1D/3+1D patterns - second slot nonzero *)
    !FreeQ[term, Derivative[n_, m_, _][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 1+1D: Only x exists (second slot > 0 with first slot > 0) *)
    !FreeQ[term, Derivative[n_, m_][_][___] /; n > 0 && m > 0], "gradient_x",
    (* Fail explicitly if spatial gradient direction cannot be determined *)
    True,
      Throw[StringJoin[
        "ExtractSpatialGradientFromMixed: Cannot extract spatial gradient from mixed derivative '",
        ToString[term], "'. ",
        "Expected pattern Derivative[time>0, spatial...] with at least one nonzero spatial slot."
      ]]
  ]
];

(* Check if term contains a spatial cross-derivative (d_x d_y, d_x d_z, d_y d_z) *)
(* Returns False if not a cross-derivative, or the operator name if it is *)
(* Pattern: Derivative[0, ...] where exactly 2 spatial slots are > 0 *)

(* Legacy boolean version for backward compatibility *)
IsSpatialCrossDerivative[term_] := IdentifySpatialCrossDerivative[term] =!= False;

(* New version that returns the specific operator name *)
IdentifySpatialCrossDerivative[term_] := Module[{},
  Which[
    (* 3+1D: cross_derivative_yz (third and fourth slots > 0, second = 0) *)
    !FreeQ[term, Derivative[0, 0, m_, p_][_][___] /; m > 0 && p > 0], "cross_derivative_yz",
    (* 3+1D: cross_derivative_xz (second and fourth slots > 0, third = 0) *)
    !FreeQ[term, Derivative[0, m_, 0, p_][_][___] /; m > 0 && p > 0], "cross_derivative_xz",
    (* 3+1D: cross_derivative_xy (second and third slots > 0, fourth = 0) *)
    !FreeQ[term, Derivative[0, m_, p_, 0][_][___] /; m > 0 && p > 0], "cross_derivative_xy",
    (* 2+1D: cross_derivative_xy (second and third slots > 0) *)
    !FreeQ[term, Derivative[0, m_, p_][_][___] /; m > 0 && p > 0], "cross_derivative_xy",
    (* Not a cross-derivative *)
    True, False
  ]
];

(* === Phase 5 (Elasticity): Directional Laplacian Detection === *)
(* Identifies pure second derivatives in a single spatial direction *)
(* Returns: "laplacian_x", "laplacian_y", "laplacian_z", or False *)
(* Used for anisotropic equations like Navier-Cauchy where ∂²_x and ∂²_y have different coefficients *)
(* Supports 1+1D (2-arg), 2+1D (3-arg), and 3+1D (4-arg) Derivative forms *)

IdentifyDirectionalLaplacian[term_] := Module[{},
  Which[
    (* 3+1D: Derivative[0, 0, 0, 2] = pure ∂²/∂z² *)
    !FreeQ[term, Derivative[0, 0, 0, 2][_][___]], "laplacian_z",

    (* 3+1D: Derivative[0, 0, 2, 0] = pure ∂²/∂y² (no x or z derivative) *)
    !FreeQ[term, Derivative[0, 0, 2, 0][_][___]], "laplacian_y",

    (* 3+1D: Derivative[0, 2, 0, 0] = pure ∂²/∂x² (no y or z derivative) *)
    !FreeQ[term, Derivative[0, 2, 0, 0][_][___]], "laplacian_x",

    (* 2+1D: Derivative[0, 2, 0] = pure ∂²/∂x² (no y derivative) *)
    !FreeQ[term, Derivative[0, 2, 0][_][___]], "laplacian_x",

    (* 2+1D: Derivative[0, 0, 2] = pure ∂²/∂y² (no x derivative) *)
    !FreeQ[term, Derivative[0, 0, 2][_][___]], "laplacian_y",

    (* 1+1D: Derivative[0, 2] = pure ∂²/∂x² (only spatial dimension) *)
    !FreeQ[term, Derivative[0, 2][_][___]], "laplacian_x",

    (* Default: not a pure directional Laplacian *)
    True, False
  ]
];

(* Helper to build term result with optional symbolic coefficient and time dependence *)
(* Only includes "coefficient_symbolic" key when symbolicCoeff is not None *)
(* Only includes "time_dependent" key when timeDependent is True *)
BuildTermResult[coeff_, op_, field_, symbCoeff_:None, timeDependent_:False] := Module[{result},
  result = <|
    "coefficient" -> N[coeff],
    "operator" -> op,
    "field" -> field
  |>;
  If[symbCoeff =!= None,
    result["coefficient_symbolic"] = symbCoeff
  ];
  If[timeDependent === True,
    result["time_dependent"] = True
  ];
  result
];

(* Check if a coefficient expression depends on coordinate symbols *)
(* Coordinate symbols from xCoba appear as f[] (zero-argument function calls) *)
(* After LHS normalization, the only surviving coordinate is typically time *)
(* Used to flag terms that need time-dependent evaluation in Python *)
IsTimeDependentCoefficient[coeffExpr_] := Module[{coordSyms},
  coordSyms = Cases[coeffExpr, f_Symbol[] :> f, {0, Infinity}];
  Length[coordSyms] > 0
];

(* === Refactored Helper Functions for IdentifyMultiFieldTerm === *)
(* Phase 3, Issue 10: Split 144-line function into smaller, testable units *)

(* Extract all function heads from a term *)
(* For f[args], extract f; for Derivative[n,m][f][args], extract f *)
ExtractFunctionHeads[term_] := Union[Join[
  (* Direct function heads: f[args] -> f *)
  Cases[term, f_Symbol[__] :> ToString[f], {0, Infinity}],
  (* Derivative heads: Derivative[n,m][f][args] -> f *)
  Cases[term, Derivative[__][f_][__] :> ToString[f], {0, Infinity}]
]];

(* Match a field name to function heads in the term *)
(* Returns {matchedFieldName, foundHeadString} or {defaultField, Null} if no match *)
MatchFieldToHeads[functionHeads_List, allFieldNames_List, defaultField_String] := Module[
  {matchedField = defaultField, foundFieldHead = Null},

  (* For each field name like "A_0", "phi_0", find matching function head *)
  (* Strategy: Match BOTH base name and index *)
  (* E.g., "A_0" matches "csA0" (base "A/a" + index "0") *)
  (* E.g., "phi_0" matches "cplPhi0" (base "phi" + index "0") *)

  Do[
    Module[{fieldParts, fieldBase, fieldIndex, headDigits, headBase},
      (* Split field name: "phi_0" -> {"phi", "0"}, "A_1" -> {"A", "1"} *)
      fieldParts = StringSplit[fn, "_"];
      fieldBase = ToLowerCase[First[fieldParts]];
      fieldIndex = Last[fieldParts];

      (* Check if any function head matches BOTH base name AND index *)
      Do[
        (* Extract trailing digits: "cplPhi0" -> "0", "csA2" -> "2" *)
        headDigits = StringCases[head, RegularExpression["\\d+$"]];
        (* Extract base: "cplPhi0" -> "cplphi", "csA2" -> "csa" (lowercase, no digits) *)
        headBase = ToLowerCase[StringReplace[head, RegularExpression["\\d+$"] -> ""]];

        (* Match if: head ends with field base name AND trailing digits match field index *)
        (* Using StringEndsQ instead of StringContainsQ to prevent "A" matching "alpha", "eta" *)
        If[Length[headDigits] > 0 &&
           headDigits[[-1]] === fieldIndex &&
           StringEndsQ[headBase, fieldBase],
          matchedField = fn;
          foundFieldHead = head;
          Break[]
        ],
        {head, functionHeads}
      ];
      If[foundFieldHead =!= Null, Break[]]
    ],
    {fn, allFieldNames}
  ];

  {matchedField, foundFieldHead}
];

(* Extract coefficient from term given the field head *)
(* Returns {numericCoeff, symbolicCoeff, isTimeDependent} *)
(* where symbolicCoeff is None or a string, isTimeDependent is True/False *)
ExtractTermCoefficient[term_, fieldHead_String, targetField_String] := Module[
  {rawCoeff, coefficient = 1.0, symbolicCoeff = None, isTimeDependent = False},

  rawCoeff = term /. {
    (* Replace Derivative[...][field][args] with 1 *)
    Derivative[__][f_][__] /; ToString[f] === fieldHead :> 1,
    (* Replace field[args] with 1 *)
    f_Symbol[__] /; ToString[f] === fieldHead :> 1
  };
  rawCoeff = Simplify[rawCoeff];

  (* Check for time dependence in coefficient *)
  isTimeDependent = IsTimeDependentCoefficient[rawCoeff];

  (* Determine numeric coefficient and symbolic representation *)
  (* Use InputForm for ToString to get clean, machine-parseable strings *)
  Which[
    NumericQ[rawCoeff],
      coefficient = N[rawCoeff];
      symbolicCoeff = None,
    (* Negative symbolic: -m2 -> -1.0, store "-m2" *)
    MatchQ[rawCoeff, Times[-1, _Symbol]],
      coefficient = -1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField],
    (* Positive symbolic: m2 -> 1.0, store "m2" *)
    MatchQ[rawCoeff, _Symbol],
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField],
    (* Try numeric evaluation *)
    NumericQ[Quiet[N[rawCoeff]]],
      coefficient = Quiet[N[rawCoeff]];
      symbolicCoeff = None,
    (* Complex symbolic expression (possibly time-dependent) *)
    True,
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      If[!isTimeDependent,
        Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField]
      ]
  ];

  {coefficient, symbolicCoeff, isTimeDependent}
];

(* Count the total derivative order of a term *)
(* Returns the maximum total order of any Derivative expression found *)
CountDerivativeOrder[term_] := Module[
  {maxOrder},
  maxOrder = 0;
  (* Check for Derivative[n][f][args] patterns (applied derivatives) *)
  Cases[term,
    Derivative[orders__][_][__] :> (maxOrder = Max[maxOrder, Total[{orders}]]),
    {0, Infinity}
  ];
  (* Also check for unapplied Derivative[n][f] patterns *)
  If[maxOrder == 0,
    Cases[term,
      Derivative[orders__][_] :> (maxOrder = Max[maxOrder, Total[{orders}]]),
      {0, Infinity}
    ]
  ];
  maxOrder
];

(* Classify operator type based on derivative structure *)
(* Returns {operatorName, shouldConvertToMomentum} *)
ClassifyOperatorType[term_] := Module[
  {orderDerivatives, dirLap, crossDeriv},

  (* Check for field without derivatives (mass/coupling term) *)
  If[FreeQ[term, Derivative],
    Return[{"identity", False}]
  ];

  (* Check for pure first-order time derivative (Hubble friction in curved spacetime) *)
  (* e.g., -2H∂_t φ in de Sitter: Derivative[1, 0][φ][t, x] with no spatial derivatives *)
  (* This must come BEFORE the mixed derivative check *)
  If[IsFirstOrderTimeDerivative[term],
    Return[{"first_derivative_t", False}]
  ];

  (* Check for mixed time-space derivatives *)
  (* Mixed derivatives like d_t d_x A become gradient_x(pi) since d_t A = pi *)
  If[IsMixedTimeSpaceDerivative[term],
    Return[{ExtractSpatialGradientFromMixed[term], True}]
  ];

  (* Check derivative order for pure spatial derivatives *)
  orderDerivatives = CountDerivativeOrder[term];

  (* Second-order spatial derivatives *)
  If[orderDerivatives == 2,
    (* Priority: directional Laplacians > cross-derivative > generic Laplacian *)
    dirLap = IdentifyDirectionalLaplacian[term];
    If[dirLap =!= False,
      Return[{dirLap, False}]
    ];
    crossDeriv = IdentifySpatialCrossDerivative[term];
    If[crossDeriv =!= False,
      Return[{crossDeriv, False}]
    ];
    Return[{"laplacian", False}]
  ];

  (* First-order spatial derivatives (gradients) *)
  If[orderDerivatives == 1,
    Return[{IdentifyGradientDirection[term], False}]
  ];

  (* Fourth-order spatial derivatives (biharmonic) *)
  If[orderDerivatives == 4,
    Return[{"biharmonic", False}]
  ];

  (* Unknown operator type - fail explicitly *)
  Throw[StringJoin[
    "ClassifyOperatorType: Unsupported derivative order ", ToString[orderDerivatives],
    " in term '", ToString[term, InputForm], "'. ",
    "Supported orders: 0 (identity), 1 (gradient), 2 (laplacian/cross), 4 (biharmonic)."
  ]]
];

(* === Main Function: Identify operator and target field === *)
(* Refactored to use helper functions above *)

IdentifyMultiFieldTerm[term_, currentFieldName_, allFieldNames_] := Module[
  {functionHeads, matchResult, targetField, foundFieldHead,
   coeffResult, coefficient, symbolicCoeff, isTimeDependent,
   operatorResult, operator, isMixedTimeSpace},

  (* Step 1: Extract function heads from term *)
  functionHeads = ExtractFunctionHeads[term];

  (* Step 2: Match field names to function heads *)
  matchResult = MatchFieldToHeads[functionHeads, allFieldNames, currentFieldName];
  targetField = matchResult[[1]];
  foundFieldHead = matchResult[[2]];

  (* Step 3: Extract coefficient (now includes time-dependence check) *)
  If[foundFieldHead =!= Null,
    coeffResult = ExtractTermCoefficient[term, foundFieldHead, targetField];
    coefficient = coeffResult[[1]];
    symbolicCoeff = coeffResult[[2]];
    isTimeDependent = coeffResult[[3]],
    (* Fallback if no field head found *)
    coefficient = 1.0;
    symbolicCoeff = None;
    isTimeDependent = False
  ];

  (* Step 4: Classify operator type *)
  operatorResult = ClassifyOperatorType[term];
  operator = operatorResult[[1]];
  isMixedTimeSpace = operatorResult[[2]];

  (* Convert to momentum field if mixed time-space derivative *)
  If[isMixedTimeSpace,
    targetField = FieldToMomentumName[targetField]
  ];

  (* Build and return result with time-dependence flag *)
  BuildTermResult[coefficient, operator, targetField, symbolicCoeff, isTimeDependent]
];

(* === Equation Conversion === *)

(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
(* Single-field EquationToJSON delegates to multi-field version *)
(* This avoids code duplication while maintaining backward compatibility *)
EquationToJSON[componentEq_, fieldName_, fieldIndex_, metadata_] :=
  EquationToJSONMultiField[componentEq, fieldName, fieldIndex, {fieldName}, metadata];

End[];
EndPackage[];
