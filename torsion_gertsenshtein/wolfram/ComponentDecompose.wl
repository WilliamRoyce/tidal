(* ::Package:: *)
(* ComponentDecompose.wl - Decompose tensor equations into scalar components *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`ComponentDecompose`",
  {"xAct`xTensor`", "xAct`xCoba`"}];

(* Public symbols *)
DecomposeToComponents::usage =
  "DecomposeToComponents[eom, field, chart] decomposes the tensor equation of \
motion into individual scalar component equations. Returns a list of \
{component_index, component_equation} pairs.";

ExtractOperatorStructure::usage =
  "ExtractOperatorStructure[componentEq, field, chart] analyzes a component \
equation and returns its structure: {coefficient, operator_type, target_field}.";

IdentifyOperator::usage =
  "IdentifyOperator[expr, field, chart] identifies if the expression is a \
known operator (laplacian, identity, gradient, etc.) applied to the field.";

Begin["`Private`"];

(* === Component Decomposition === *)

DecomposeToComponents[eom_, field_, chart_] := Module[
  {dim, components, indices, componentEq, result},

  (* Get the dimension from the chart *)
  dim = Length[IndicesOfChart[chart]];

  (* For a scalar field, there's only one component *)
  If[TensorRank[field] === 0,
    (* Scalar field - evaluate in coordinates *)
    componentEq = ToBasis[chart][eom];
    componentEq = TraceBasisDummy[componentEq];
    Return[{{0, componentEq}}]
  ];

  (* For a vector field A_a, extract each component *)
  If[TensorRank[field] === {1},
    result = Table[
      {
        idx,
        ExtractVectorComponent[eom, field, chart, idx]
      },
      {idx, 0, dim - 1}
    ];
    Return[result]
  ];

  (* For higher-rank tensors, iterate over all index combinations *)
  (* This is a simplified version - extend as needed *)
  Print["Warning: Higher-rank tensor decomposition not fully implemented"];
  {{0, eom}}
];

(* === Vector Component Extraction === *)

ExtractVectorComponent[eom_, field_, chart_, componentIndex_] := Module[
  {componentEq, basis},

  (* Get the coordinate basis *)
  basis = BasisOfChart[chart];

  (* For a vector equation EOM^a = 0, extract the componentIndex-th component *)
  (* This involves contracting with the dual basis vector *)

  componentEq = eom;

  (* Convert to coordinate basis *)
  componentEq = ToBasis[chart][componentEq];

  (* Extract the specific component by setting the index *)
  componentEq = componentEq /. {
    (* Replace free vector index with component value *)
    PatternSequence[a_Symbol /; MemberQ[IndicesOfChart[chart], a]] -> componentIndex
  };

  (* Trace over dummy indices *)
  componentEq = TraceBasisDummy[componentEq];

  componentEq
];

(* === Operator Structure Extraction === *)

ExtractOperatorStructure[componentEq_, field_, chart_] := Module[
  {terms, operatorTerms},

  (* Split equation into additive terms *)
  terms = If[Head[componentEq] === Plus,
    List @@ componentEq,
    {componentEq}
  ];

  (* Analyze each term *)
  operatorTerms = Map[
    AnalyzeTerm[#, field, chart] &,
    terms
  ];

  operatorTerms
];

AnalyzeTerm[term_, field_, chart_] := Module[
  {coefficient, operator, targetField},

  (* Default structure *)
  coefficient = 1;
  operator = "unknown";
  targetField = ToString[field];

  (* Check for mass-like term: coefficient * field *)
  If[FreeQ[term, Derivative] && !FreeQ[term, field],
    coefficient = term /. field -> 1;
    operator = "identity";
    Return[<|"coefficient" -> coefficient, "operator" -> operator, "field" -> targetField|>]
  ];

  (* Check for Laplacian-like term *)
  (* In coordinates, Laplacian appears as second derivatives *)
  If[!FreeQ[term, Derivative[2, 0]] || !FreeQ[term, Derivative[0, 2]],
    (* This is a second derivative - likely part of Laplacian or d'Alembertian *)
    coefficient = ExtractCoefficient[term, field];
    operator = "laplacian";  (* or "spatial_derivative" for more precision *)
    Return[<|"coefficient" -> coefficient, "operator" -> operator, "field" -> targetField|>]
  ];

  (* Fallback *)
  <|"coefficient" -> 1, "operator" -> "unknown", "field" -> targetField|>
];

ExtractCoefficient[term_, field_] := Module[
  {coeff},
  (* Extract numeric coefficient from term *)
  coeff = term /. {field -> 1, _Derivative[__][field] -> 1};
  If[NumericQ[coeff], coeff, 1]
];

(* === Operator Identification === *)

IdentifyOperator[expr_, field_, chart_] := Module[
  {},

  (* Check for d'Alembertian (wave operator) *)
  (* In 1+1D with signature (-,+): Box = -d^2/dt^2 + d^2/dx^2 *)

  (* Check for Laplacian (spatial only) *)
  (* In 1D: Laplacian = d^2/dx^2 *)

  (* Check for identity *)
  If[expr === field || expr === field[],
    Return["identity"]
  ];

  (* Check for gradient *)
  If[MatchQ[expr, Derivative[1, 0][field] | Derivative[0, 1][field]],
    Return["gradient"]
  ];

  "unknown"
];

End[];
EndPackage[];
