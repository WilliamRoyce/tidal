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

  timeDerivTerm = Select[terms, !FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &];
  rhs = Total[Select[terms, FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &]];

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

  (* For each field name like "phi_0", "chi_0", find matching function head *)
  matchedField = currentFieldName;  (* Default to current field *)

  Do[
    (* Extract base name: "phi_0" -> "phi", "chi_0" -> "chi" *)
    fieldBaseName = ToLowerCase[StringSplit[fn, "_"][[1]]];

    (* Check if any function head contains this field name (case insensitive) *)
    Do[
      If[StringContainsQ[ToLowerCase[head], fieldBaseName],
        matchedField = fn;
        foundFieldHead = head;
        Break[]
      ],
      {head, functionHeads}
    ];
    If[foundFieldHead =!= Null, Break[]],
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

  (* Check for second derivatives (Laplacian) *)
  orderDerivatives = CountDerivativeOrder[term];

  If[orderDerivatives == 2,
    operator = "laplacian";
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  If[orderDerivatives == 1,
    operator = "gradient_x";
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

  timeDerivTerm = Select[terms,
    !FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &
  ];

  (* RHS = everything except time derivative, with sign flip if time deriv was negative *)
  (* From eq: spatial - time = 0, we get: time = spatial *)
  (* So RHS is the spatial part *)
  rhs = Total[Select[terms, FreeQ[#, Derivative[n_, 0][_] /; n >= 2] &]];

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
    (* First derivatives indicate gradient *)
    operator = "gradient_x";
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
