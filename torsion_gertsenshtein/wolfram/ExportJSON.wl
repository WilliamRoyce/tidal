(* ::Package:: *)
(* ExportJSON.wl - Export field equations to JSON format for Python pipeline *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`ExportJSON`"];

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
    coefficient = ExtractCoefficientFromTerm[term, fieldPattern];
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
    coefficient = ExtractCoefficientFromTerm[term, fieldPattern];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  If[orderDerivatives == 1,
    (* First derivatives indicate gradient *)
    operator = "gradient_x";
    coefficient = ExtractCoefficientFromTerm[term, fieldPattern];
    Return[<|"coefficient" -> N[coefficient], "operator" -> operator, "field" -> targetField|>]
  ];

  (* If we can't identify it, return unknown *)
  Print["Warning: Could not identify operator for term: ", term];
  <|"coefficient" -> 1.0, "operator" -> "unknown", "field" -> targetField|>
];

(* Extract numeric coefficient from a term *)
ExtractCoefficientFromTerm[term_, fieldName_] := Module[
  {coeff, fieldSymbol, derivPattern},

  (* Create pattern for field function with any arguments *)
  fieldSymbol = Symbol[fieldName];

  (* Replace field and its derivatives with 1 to extract coefficient *)
  coeff = term /. {
    (* Match Derivative[...][f][args] form (applied derivatives) *)
    Derivative[__][f_][__] /; StringContainsQ[ToString[f], fieldName] :> 1,
    (* Match f[args] form (bare field) *)
    f_[__] /; StringMatchQ[ToString[f], fieldName ~~ DigitCharacter ...] :> 1
  };

  (* Simplify and try to get numeric value *)
  coeff = Simplify[coeff];

  (* Try to evaluate - N will work for both numeric and some symbolic expressions *)
  If[NumericQ[coeff],
    Return[coeff]
  ];

  (* If coefficient contains known constant symbols, try to extract numeric part *)
  (* For expressions like -m2, m2, etc., try to factor out the sign *)
  If[MatchQ[coeff, Times[-1, _Symbol]],
    (* Negative symbolic coefficient: -m2 -> use the symbol value if available *)
    coeff = -1.0;  (* Return just the sign, assuming unit coefficient for unknown *)
    Return[coeff]
  ];

  If[MatchQ[coeff, _Symbol],
    (* Positive symbolic coefficient: m2 -> assume unit coefficient *)
    coeff = 1.0;
    Return[coeff]
  ];

  (* Try N as last resort *)
  coeff = Quiet[N[coeff]];
  If[NumericQ[coeff],
    coeff,
    (* Default to 1.0 if all else fails *)
    1.0
  ]
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
