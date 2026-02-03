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
  {terms, rhsTerms},

  (* Parse the equation to extract RHS terms *)
  (* Equation form: d^2/dt^2 field = RHS *)
  (* RHS typically contains: laplacian terms, mass terms, coupling terms *)

  rhsTerms = ParseEquationRHS[componentEq, fieldName, metadata];

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
  {terms, parsedTerms},

  (* Default: wave equation with laplacian and mass *)
  (* For more sophisticated parsing, analyze the actual expression *)

  parsedTerms = {};

  (* Add Laplacian term (coefficient 1.0 for wave equation) *)
  AppendTo[parsedTerms,
    <|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> fieldName|>
  ];

  (* Add mass term if present *)
  If[KeyExistsQ[metadata, "mass_squared"] && metadata["mass_squared"] =!= 0,
    AppendTo[parsedTerms,
      <|"coefficient" -> -N[metadata["mass_squared"]], "operator" -> "identity", "field" -> fieldName|>
    ]
  ];

  parsedTerms
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
