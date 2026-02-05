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
EquationToJSONMultiField[componentEq_, fieldName_, fieldIndex_, allFieldNames_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm},

  (* Same logic as EquationToJSON but with cross-field awareness *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Use dimension-agnostic helpers for time derivative detection *)
  (* Works for both 1+1D (2-arg Derivative) and 2+1D (3-arg Derivative) *)
  timeDerivTerm = Select[terms, ContainsTimeDerivative[#, 2] &];
  (* RHS = everything EXCEPT pure 2nd-order time derivatives *)
  (* Mixed time-space derivatives ARE included - they get converted to momentum gradients *)
  (* by IdentifyMultiFieldTerm (e.g., d_t d_x A -> gradient_x(pi)) *)
  rhs = Total[Select[terms, !ContainsTimeDerivative[#, 2] &]];

  If[Length[timeDerivTerm] > 0 &&
     Head[timeDerivTerm[[1]]] === Times &&
     timeDerivTerm[[1]][[1]] == 1,
    rhs = -rhs
  ];

  (* Parse RHS with cross-field detection *)
  rhsTerms = ParseMultiFieldRHS[rhs, fieldName, allFieldNames, metadata];

  <|
    "field" -> fieldName,
    "lhs" -> "d2_t(" <> fieldName <> ")",
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

(* === Time Derivative Detection Helpers (for LHS/RHS separation) === *)

(* Detect if a term contains time derivatives of at least order minOrder *)
(* Time index is always the first slot in Derivative[dt, dx, ...] *)
ContainsTimeDerivative[term_, minOrder_:2] := Module[{},
  Which[
    (* 2+1D: Derivative[n, m, p] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _, _][_][___] /; n >= minOrder], True,
    (* 1+1D: Derivative[n, m] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _][_][___] /; n >= minOrder], True,
    (* Default: no time derivative of required order *)
    True, False
  ]
];

(* Check for mixed time-space derivatives that shouldn't be on RHS *)
(* Returns True if term has BOTH time AND space derivatives *)
IsMixedTimeSpaceDerivative[term_] := Module[{},
  Which[
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
(* IMPORTANT: First slot (time) must be 0 for pure spatial gradients *)
IdentifyGradientDirection[term_] := Module[{},
  Which[
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

(* Map field name to momentum field name *)
(* For py-pde state [field_0, pi_0, field_1, pi_1, ...], momentum is at odd indices *)
(* "A_0" -> "pi_0", "A_1" -> "pi_1", "phi_0" -> "pi_0" *)
(* IMPORTANT: Python expects numeric indices only (pi_0, pi_1, etc.) *)
FieldToMomentumName[fieldName_String] := Module[{parts, idx},
  parts = StringSplit[fieldName, "_"];
  If[Length[parts] >= 2,
    idx = Last[parts];
    (* Validate that index is numeric (digits only) *)
    If[StringMatchQ[idx, DigitCharacter..],
      "pi_" <> idx,  (* Valid: A_0 -> pi_0, phi_1 -> pi_1 *)
      (* Non-numeric suffix - use "0" as fallback with warning *)
      Print["Warning: Field name '", fieldName, "' has non-numeric suffix '", idx,
            "'. Using pi_0 as fallback. For proper momentum mapping, ",
            "use field names like A_0, A_1, phi_0."];
      "pi_0"
    ],
    (* No underscore - use "0" as fallback with warning *)
    Print["Warning: Field name '", fieldName, "' has no numeric suffix. ",
          "Using pi_0 as fallback. For proper momentum mapping, ",
          "use field names like A_0, A_1, phi_0."];
    "pi_0"
  ]
];

(* Extract spatial gradient direction from mixed time-space derivative *)
(* Derivative[1, 1, 0] (d_t d_x) -> "gradient_x" *)
(* Derivative[1, 0, 1] (d_t d_y) -> "gradient_y" *)
(* These represent d_x(pi) or d_y(pi) since d_t phi = pi *)
(*
   WARNING: Derivative[1, 1, 1] (d_t d_x d_y) is a mixed time + cross-spatial derivative.
   This cannot be simply represented as gradient_x(pi) or gradient_y(pi).
   The proper representation would require cross_derivative_xy(pi), which would need
   a more complex operator structure. For now, we default to gradient_x with a warning.
*)
ExtractSpatialGradientFromMixed[term_] := Module[{},
  (* First check for the problematic Derivative[n, m, p] where all are > 0 *)
  If[!FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && m > 0 && p > 0],
    Print["Warning: Mixed time + cross-spatial derivative detected (e.g., d_t d_x d_y). ",
          "This term cannot be exactly represented as a simple momentum gradient. ",
          "Approximating as gradient_x(pi). Consider gauge-fixing to eliminate such terms."];
  ];

  Which[
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

(* Check if term contains a spatial cross-derivative (d_x d_y, NOT d^2_x or d^2_y) *)
(* Pattern: Derivative[0, m, p] where BOTH m > 0 AND p > 0 *)
IsSpatialCrossDerivative[term_] := Module[{},
  !FreeQ[term, Derivative[0, m_, p_][_][___] /; m > 0 && p > 0]
];

(* Identify operator and target field for multi-field systems *)
IdentifyMultiFieldTerm[term_, currentFieldName_, allFieldNames_, metadata_] := Module[
  {coefficient, operator, targetField, orderDerivatives, foundFieldHead,
   functionHeads, fieldBaseName, matchedField},

  operator = "unknown";
  coefficient = 1.0;
  targetField = currentFieldName;
  foundFieldHead = Null;

  (* Extract function heads from the term *)
  (* For f[args], extract f; for Derivative[n,m][f][args], extract f *)
  functionHeads = Union[Join[
    (* Direct function heads: f[args] -> f *)
    Cases[term, f_Symbol[__] :> ToString[f], {0, Infinity}],
    (* Derivative heads: Derivative[n,m][f][args] -> f *)
    Cases[term, Derivative[__][f_][__] :> ToString[f], {0, Infinity}]
  ]];

  (* For each field name like "A_0", "phi_0", find matching function head *)
  (* Strategy: Match BOTH base name and index *)
  (* E.g., "A_0" matches "csA0" (base "A/a" + index "0") *)
  (* E.g., "phi_0" matches "cplPhi0" (base "phi" + index "0") *)
  matchedField = currentFieldName;  (* Default to current field *)

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

  targetField = matchedField;

  (* Extract numeric coefficient by replacing the field function with 1 *)
  If[foundFieldHead =!= Null,
    coefficient = term /. {
      (* Replace Derivative[...][field][args] with 1 *)
      Derivative[__][f_][__] /; ToString[f] === foundFieldHead :> 1,
      (* Replace field[args] with 1 *)
      f_Symbol[__] /; ToString[f] === foundFieldHead :> 1
    };
    coefficient = Simplify[coefficient];
    If[!NumericQ[coefficient], coefficient = 1.0],
    (* Fallback if no field head found *)
    coefficient = 1.0
  ];

  (* Check for field without derivatives (mass/coupling term) *)
  If[FreeQ[term, Derivative],
    operator = "identity";
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  (* === Phase 4: Handle mixed time-space derivatives === *)
  (* Mixed derivatives like d_t d_x A become gradient_x(pi) since d_t A = pi *)
  (* This preserves physics that would otherwise be lost by dropping mixed terms *)
  If[IsMixedTimeSpaceDerivative[term],
    operator = ExtractSpatialGradientFromMixed[term];
    (* Convert field name to momentum field name *)
    targetField = FieldToMomentumName[targetField];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  (* Check for derivative order *)
  orderDerivatives = CountDerivativeOrder[term];

  (* === Phase 4: Handle spatial cross-derivatives === *)
  (* Cross-derivatives like d_x d_y A are NOT laplacian - use special operator *)
  If[orderDerivatives == 2,
    If[IsSpatialCrossDerivative[term],
      operator = "cross_derivative_xy";
      Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>],
      (* Regular 2nd order spatial (laplacian) *)
      operator = "laplacian";
      Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
    ]
  ];

  If[orderDerivatives == 1,
    (* Use helper to distinguish gradient_x vs gradient_y *)
    operator = IdentifyGradientDirection[term];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  <|"coefficient" -> N[coefficient], "operator" -> "unknown", "field" -> targetField|>
];

(* === Equation Conversion === *)

EquationToJSON[componentEq_, fieldName_, fieldIndex_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm},

  (* The equation from ComponentDecompose is in the form:
     d²/dx² field - d²/dt² field = 0
     We need to rearrange to: d²/dt² field = d²/dx² field (RHS)

     In component form: Derivative[0,2][f][t,x] - Derivative[2,0][f][t,x] = 0
     Time derivative: Derivative[2,0] or Derivative[n,0] where first slot is time
     Spatial derivative: Derivative[0,m] or Derivative[0,n] where second slot is space *)

  (* Split into additive terms *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Separate time derivative terms (Derivative[n,0] where n >= 2) from spatial terms *)
  (* In 1+1D: first derivative index is time (t), second is space (x) *)
  (* Derivative[2,0] = d²/dt², Derivative[0,2] = d²/dx² *)

  (* Use dimension-agnostic helpers for time derivative detection *)
  (* Works for both 1+1D (2-arg Derivative) and 2+1D (3-arg Derivative) *)
  timeDerivTerm = Select[terms, ContainsTimeDerivative[#, 2] &];

  (* RHS = everything except pure 2nd-order time derivatives, with sign flip if time deriv was negative *)
  (* From eq: spatial - time = 0, we get: time = spatial *)
  (* Mixed time-space derivatives ARE included - they get converted to momentum gradients *)
  rhs = Total[Select[terms, !ContainsTimeDerivative[#, 2] &]];

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

  <|
    "field" -> fieldName,
    "lhs" -> "d2_t(" <> fieldName <> ")",
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
    (* Second derivatives indicate Laplacian or Box operator *)
    (* In flat space after gauge fixing, typically just Laplacian in spatial coordinates *)
    operator = "laplacian";
    coefficient = ExtractNumericCoefficient[term, fieldPattern];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
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
