(* ::Package:: *)
(* ExportJSON.wl - Export field equations to JSON format for Python pipeline *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`ExportJSON`",
  {"TorsionGertsenshtein`CommonUtilities`"}];

(* Public symbols *)
ExportEquationSystem::usage =
  "ExportEquationSystem[componentEqs, metadata, outputPath] exports the \
component equations to a JSON file compatible with the Python pipeline.";

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

GetTimeDerivativeOrder::usage;
DetectLHSTimeOrder::usage;
BuildLHSStructure::usage;

Begin["`Private`"];

(* === JSON Structure Building === *)

BuildJSONStructure[componentEqs_, metadata_Association] := Module[
  {json, fields, equations},

  (* Extract field information from component equations *)
  fields = Table[
    <|
      "name" -> metadata["field_prefix"] <> "_" <> ToString[eq[[1]]],
      "index" -> eq[[1]],
      "is_dynamical" -> True
    |>,
    {eq, componentEqs}
  ];

  (* Convert equations to JSON format *)
  equations = Table[
    EquationToJSON[
      eq[[2]],
      metadata["field_prefix"] <> "_" <> ToString[eq[[1]]],
      eq[[1]],
      metadata
    ],
    {eq, componentEqs}
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
        ConstantArray[0.0, {Length[componentEqs], Length[componentEqs]}]],
      "coupling_matrix" -> Lookup[metadata, "coupling_matrix",
        ConstantArray[0.0, {Length[componentEqs], Length[componentEqs]}]]
    |>
  |>;

  json
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

  If[Length[timeDerivTerm] > 0 &&
     Head[timeDerivTerm[[1]]] === Times &&
     timeDerivTerm[[1]][[1]] == 1,
    rhs = -rhs
  ];

  (* Parse RHS with cross-field detection *)
  rhsTerms = ParseMultiFieldRHS[rhs, fieldName, allFieldNames, metadata];

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
ParseMultiFieldRHS[eq_, currentFieldName_, allFieldNames_, metadata_] := Module[
  {terms, parsedTerms},

  terms = If[Head[eq] === Plus, List @@ eq, {eq}];
  parsedTerms = Map[IdentifyMultiFieldTerm[#, currentFieldName, allFieldNames, metadata] &, terms];
  parsedTerms = DeleteCases[parsedTerms, Nothing];

  If[Length[parsedTerms] === 0,
    parsedTerms = {<|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> currentFieldName|>}
  ];

  parsedTerms
];

(* === LHS Structure Detection (Phase 2, Issue 6) === *)
(* Supports parabolic (heat: d_t φ), elliptic (Poisson: ∇²φ = f), and hyperbolic (wave: d2_t φ) PDEs *)

GetTimeDerivativeOrder::usage =
  "GetTimeDerivativeOrder[term] returns the order of time derivative in a term. \
Returns 0 if no time derivative, 1 for d_t, 2 for d2_t, etc.";

DetectLHSTimeOrder::usage =
  "DetectLHSTimeOrder[equation] detects the maximum time derivative order in the equation. \
Used to determine if PDE is elliptic (0), parabolic (1), or hyperbolic (2+).";

BuildLHSStructure::usage =
  "BuildLHSStructure[fieldName, timeOrder] builds a structured LHS representation. \
Returns an Association with \"expression\" and \"order\" keys for flexible PDE types.";

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

    (* Default to gradient_x if detection fails *)
    True, "gradient_x"
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

  (* Fallback: treat as simple with index 0 *)
  Print["Warning: Could not parse field name '", name, "'. Treating as simple field with index 0."];
  <|"base" -> name, "index" -> 0, "format" -> "unknown"|>
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
    (* Fallback for other patterns *)
    !FreeQ[term, Derivative[n_, m_, _][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 1+1D: Only x exists (second slot > 0 with first slot > 0) *)
    !FreeQ[term, Derivative[n_, m_][_][___] /; n > 0 && m > 0], "gradient_x",
    (* Default *)
    True, "gradient_x"
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

(* Helper to build term result with optional symbolic coefficient *)
(* Only includes "coefficient_symbolic" key when symbolicCoeff is not None *)
BuildTermResult[coeff_, op_, field_, symbCoeff_:None] := Module[{result},
  result = <|
    "coefficient" -> N[coeff],
    "operator" -> op,
    "field" -> field
  |>;
  If[symbCoeff =!= None,
    result["coefficient_symbolic"] = symbCoeff
  ];
  result
];

(* === Refactored Helper Functions for IdentifyMultiFieldTerm === *)
(* Phase 3, Issue 10: Split 144-line function into smaller, testable units *)

(* Extract all function heads from a term *)
(* For f[args], extract f; for Derivative[n,m][f][args], extract f *)
ExtractFunctionHeads::usage =
  "ExtractFunctionHeads[term] extracts all function head names from an expression. \
Returns a list of strings representing function names (e.g., {\"phi0\", \"csA1\"}).";

ExtractFunctionHeads[term_] := Union[Join[
  (* Direct function heads: f[args] -> f *)
  Cases[term, f_Symbol[__] :> ToString[f], {0, Infinity}],
  (* Derivative heads: Derivative[n,m][f][args] -> f *)
  Cases[term, Derivative[__][f_][__] :> ToString[f], {0, Infinity}]
]];

(* Match a field name to function heads in the term *)
(* Returns {matchedFieldName, foundHeadString} or {defaultField, Null} if no match *)
MatchFieldToHeads::usage =
  "MatchFieldToHeads[functionHeads, allFieldNames, defaultField] finds which field \
a term operates on by matching function heads to field names. Returns \
{matchedFieldName, foundHeadString}.";

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

        (* Match if: head contains field base name AND ends with field index *)
        If[Length[headDigits] > 0 &&
           headDigits[[-1]] === fieldIndex &&
           StringContainsQ[headBase, fieldBase],
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
(* Returns {numericCoeff, symbolicCoeff} where symbolicCoeff is None or a string *)
ExtractTermCoefficient::usage =
  "ExtractTermCoefficient[term, fieldHead, targetField] extracts the coefficient \
from a term by replacing the field with 1. Returns {numericCoeff, symbolicCoeff} \
where symbolicCoeff is None or a string representation.";

ExtractTermCoefficient[term_, fieldHead_String, targetField_String] := Module[
  {rawCoeff, coefficient = 1.0, symbolicCoeff = None},

  rawCoeff = term /. {
    (* Replace Derivative[...][field][args] with 1 *)
    Derivative[__][f_][__] /; ToString[f] === fieldHead :> 1,
    (* Replace field[args] with 1 *)
    f_Symbol[__] /; ToString[f] === fieldHead :> 1
  };
  rawCoeff = Simplify[rawCoeff];

  (* Determine numeric coefficient and symbolic representation *)
  Which[
    NumericQ[rawCoeff],
      coefficient = N[rawCoeff];
      symbolicCoeff = None,
    (* Negative symbolic: -m2 -> -1.0, store "-m2" *)
    MatchQ[rawCoeff, Times[-1, _Symbol]],
      coefficient = -1.0;
      symbolicCoeff = ToString[rawCoeff];
      Message[ExtractNumericCoefficient::symbolic, rawCoeff, targetField],
    (* Positive symbolic: m2 -> 1.0, store "m2" *)
    MatchQ[rawCoeff, _Symbol],
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff];
      Message[ExtractNumericCoefficient::symbolic, rawCoeff, targetField],
    (* Try numeric evaluation *)
    NumericQ[Quiet[N[rawCoeff]]],
      coefficient = Quiet[N[rawCoeff]];
      symbolicCoeff = None,
    (* Complex symbolic expression *)
    True,
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff];
      Message[ExtractNumericCoefficient::symbolic, rawCoeff, targetField]
  ];

  {coefficient, symbolicCoeff}
];

(* Classify operator type based on derivative structure *)
(* Returns {operatorName, shouldConvertToMomentum} *)
ClassifyOperatorType::usage =
  "ClassifyOperatorType[term] determines the differential operator type from the \
derivative structure. Returns {operatorName, isMixedTimeSpace} where operatorName \
is a string like \"laplacian\", \"gradient_x\", \"identity\", etc.";

ClassifyOperatorType[term_] := Module[
  {orderDerivatives, dirLap, crossDeriv},

  (* Check for field without derivatives (mass/coupling term) *)
  If[FreeQ[term, Derivative],
    Return[{"identity", False}]
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

  (* Unknown operator type *)
  {"unknown", False}
];

(* === Main Function: Identify operator and target field === *)
(* Refactored to use helper functions above *)

IdentifyMultiFieldTerm[term_, currentFieldName_, allFieldNames_, metadata_] := Module[
  {functionHeads, matchResult, targetField, foundFieldHead,
   coeffResult, coefficient, symbolicCoeff,
   operatorResult, operator, isMixedTimeSpace},

  (* Step 1: Extract function heads from term *)
  functionHeads = ExtractFunctionHeads[term];

  (* Step 2: Match field names to function heads *)
  matchResult = MatchFieldToHeads[functionHeads, allFieldNames, currentFieldName];
  targetField = matchResult[[1]];
  foundFieldHead = matchResult[[2]];

  (* Step 3: Extract coefficient *)
  If[foundFieldHead =!= Null,
    coeffResult = ExtractTermCoefficient[term, foundFieldHead, targetField];
    coefficient = coeffResult[[1]];
    symbolicCoeff = coeffResult[[2]],
    (* Fallback if no field head found *)
    coefficient = 1.0;
    symbolicCoeff = None
  ];

  (* Step 4: Classify operator type *)
  operatorResult = ClassifyOperatorType[term];
  operator = operatorResult[[1]];
  isMixedTimeSpace = operatorResult[[2]];

  (* Convert to momentum field if mixed time-space derivative *)
  If[isMixedTimeSpace,
    targetField = FieldToMomentumName[targetField]
  ];

  (* Build and return result *)
  BuildTermResult[coefficient, operator, targetField, symbolicCoeff]
];

(* === Equation Conversion === *)

(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
EquationToJSON[componentEq_, fieldName_, fieldIndex_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm, lhsTimeOrder, lhsStructure},

  (* The equation from ComponentDecompose is in the form:
     d²/dx² field - d²/dt² field = 0
     We need to rearrange to: d²/dt² field = d²/dx² field (RHS)

     In component form: Derivative[0,2][f][t,x] - Derivative[2,0][f][t,x] = 0
     Time derivative: Derivative[2,0] or Derivative[n,0] where first slot is time
     Spatial derivative: Derivative[0,m] or Derivative[0,n] where second slot is space *)

  (* Split into additive terms *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Detect the time derivative order to determine PDE type *)
  (* elliptic (0), parabolic (1), hyperbolic (2+) *)
  lhsTimeOrder = DetectLHSTimeOrder[componentEq];

  (* Use dimension-agnostic helpers for time derivative detection *)
  (* Works for both 1+1D (2-arg Derivative) and 2+1D (3-arg Derivative) *)
  timeDerivTerm = Select[terms, ContainsTimeDerivative[#, lhsTimeOrder] &];

  (* RHS = everything except pure time derivatives of detected order, with sign flip if time deriv was negative *)
  (* From eq: spatial - time = 0, we get: time = spatial *)
  (* Mixed time-space derivatives ARE included - they get converted to momentum gradients *)
  rhs = Total[Select[terms, !ContainsTimeDerivative[#, lhsTimeOrder] &]];

  (* If time derivative term was negative, RHS keeps its sign *)
  (* If time derivative term was positive, RHS flips sign *)
  (* In wave equation: Derivative[0,2] - Derivative[2,0] = 0 *)
  (* Time term coefficient is -1, so RHS = spatial (positive) *)
  If[Length[timeDerivTerm] > 0 &&
     Head[timeDerivTerm[[1]]] === Times &&
     timeDerivTerm[[1]][[1]] == 1,  (* positive time derivative coefficient *)
    rhs = -rhs
  ];

  (* Parse RHS for operator identification *)
  rhsTerms = ParseEquationRHS[rhs, fieldName, metadata];

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

ParseEquationRHS[eq_, fieldName_, metadata_] := Module[
  {terms, parsedTerms, fieldPattern},

  (* Parse the actual symbolic structure instead of assuming wave equation *)

  parsedTerms = {};

  (* Extract the field name without index suffix for pattern matching *)
  (* fieldName is like "A_0", extract "A" *)
  fieldPattern = StringSplit[fieldName, "_"][[1]];

  (* Split equation into additive terms *)
  terms = If[Head[eq] === Plus, List @@ eq, {eq}];

  (* Analyze each term to identify operator and coefficient *)
  (* Pass both the pattern for matching and the full field name for output *)
  parsedTerms = Map[IdentifyOperatorTerm[#, fieldPattern, fieldName, metadata] &, terms];

  (* Filter out any Nothing entries *)
  parsedTerms = DeleteCases[parsedTerms, Nothing];

  (* If no terms were identified, use fallback *)
  If[Length[parsedTerms] === 0,
    (* Fallback: assume laplacian term *)
    parsedTerms = {<|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> fieldName|>}
  ];

  parsedTerms
];

(* Identify operator type and coefficient for a single term *)
(* fieldPattern is for matching (e.g., "A"), fullFieldName is for output (e.g., "A_0") *)
IdentifyOperatorTerm[term_, fieldPattern_, fullFieldName_, metadata_] := Module[
  {coefficient, operator, targetField, orderDerivatives},

  (* Initialize *)
  operator = "unknown";
  coefficient = 1.0;
  targetField = fullFieldName;  (* Use the full component name in output *)

  (* Check for field without derivatives (mass-like term) *)
  If[FreeQ[term, Derivative],
    (* This is a mass or identity term *)
    operator = "identity";
    (* Extract coefficient using pattern for matching *)
    coefficient = ExtractNumericCoefficient[term, fieldPattern];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  (* Check for second derivatives (Laplacian or d'Alembertian) *)
  (* In component form after ToBasis, these appear as explicit coordinate derivatives *)

  (* Count derivative order *)
  orderDerivatives = CountDerivativeOrder[term];

  If[orderDerivatives == 2,
    (* Second derivatives - check for directional Laplacians first *)
    Module[{dirLap = IdentifyDirectionalLaplacian[term], crossDeriv},
      If[dirLap =!= False,
        operator = dirLap;
        coefficient = ExtractNumericCoefficient[term, fieldPattern];
        Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
      ];
      (* Check for cross-derivative (xy, xz, yz) *)
      crossDeriv = IdentifySpatialCrossDerivative[term];
      If[crossDeriv =!= False,
        operator = crossDeriv;  (* Returns "cross_derivative_xy", "cross_derivative_xz", or "cross_derivative_yz" *)
        coefficient = ExtractNumericCoefficient[term, fieldPattern];
        Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
      ];
      (* Fallback: generic Laplacian *)
      operator = "laplacian";
      coefficient = ExtractNumericCoefficient[term, fieldPattern];
      Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
    ]
  ];

  If[orderDerivatives == 1,
    (* First derivatives indicate gradient - distinguish x vs y direction *)
    operator = IdentifyGradientDirection[term];
    coefficient = ExtractNumericCoefficient[term, fieldPattern];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  (* If we can't identify it, return unknown *)
  Print["Warning: Could not identify operator for term: ", term];
  <|"coefficient" -> 1.0, "operator" -> "unknown", "field" -> targetField|>
];

(* Count the order of derivatives in a term *)
CountDerivativeOrder[term_] := Module[
  {maxOrder},

  (* Find all Derivative expressions and get maximum order *)
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

(* === File Export === *)

ExportEquationSystem[componentEqs_, metadata_Association, outputPath_String] := Module[
  {json, jsonString},

  (* Build JSON structure *)
  json = BuildJSONStructure[componentEqs, metadata];

  (* Export to file *)
  Export[outputPath, json, "JSON"];

  Print["Exported equation system to: ", outputPath];

  json
];

End[];
EndPackage[];
