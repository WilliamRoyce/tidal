(* ::Package:: *)
(* ComponentDecompose.wl - Decompose tensor equations into scalar components *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`ComponentDecompose`",
  {"xAct`xTensor`", "xAct`xCoba`",
   "TorsionGertsenshtein`CommonUtilities`"}];

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

DecomposeToComponents[eom_, field_, chart_] :=
  DecomposeToComponents[eom, field, chart, {}];

DecomposeToComponents[eom_, field_, chart_, additionalFields_List] := Module[
  {dim, components, indices, componentEq, result, fieldHead, fieldRank, allFieldHeads},

  (* Get the dimension dynamically from the chart *)
  (* ScalarsOfChart returns the coordinate symbols, e.g., {t[], x[]} or {t[], x[], y[]} *)
  dim = Length[ScalarsOfChart[chart]];

  (* Determine field rank using SlotsOfTensor *)
  (* Extract field head from applied form like phi[] or A[-a] *)
  fieldHead = ExtractFieldHead[field];
  fieldRank = Length[SlotsOfTensor[fieldHead]];

  (* Collect all field heads (primary + additional) for coordinate transformation *)
  allFieldHeads = Join[{fieldHead}, ExtractFieldHead /@ additionalFields];

  (* For a scalar field, there's only one component *)
  If[fieldRank === 0,
    (* Scalar field - evaluate in coordinates *)
    componentEq = ToBasis[chart][eom];
    componentEq = TraceBasisDummy[componentEq];

    (* In flat Minkowski space, Christoffel symbols vanish *)
    componentEq = RemoveChristoffelSymbols[componentEq];
    componentEq = Expand[componentEq];

    (* Evaluate epsilon tensor components to numeric ±1 values *)
    (* This handles Chern-Simons and other topological terms with Levi-Civita *)
    componentEq = EvaluateEpsilonComponents[componentEq, chart];
    componentEq = Expand[componentEq];

    (* Get coordinate symbols *)
    Module[{coordSyms},
      coordSyms = GetCoordinateSymbols[chart];

      (* Evaluate metric components for Minkowski signature (-1, +1) *)
      componentEq = EvaluateMinkowskiMetric[componentEq, chart];
      componentEq = Expand[componentEq];

      (* Replace ALL scalar fields with functions of coordinates *)
      (* This ensures cross-field terms are properly transformed *)
      Do[
        With[{fh = afh, cs = coordSyms},
          componentEq = componentEq /. {
            fh[] :> Symbol[ToString[fh] <> "0"][Sequence @@ cs]
          }
        ],
        {afh, allFieldHeads}
      ];

      (* Convert coordinate derivatives to explicit Derivative form *)
      componentEq = ConvertCDToDerivatives[componentEq, chart];
    ];

    (* Expand to get explicit Derivative[...] form *)
    componentEq = Expand[componentEq];
    Return[{{0, componentEq}}]
  ];

  (* For a vector field A_a, extract each component *)
  If[fieldRank === 1,
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
  {componentEq, freeIdx, fieldHead, coordSyms, basisIdx},

  (*
     For a vector equation EOM[-a] = 0, extract the componentIndex-th component.

     Strategy:
     1. Identify the free index in the field (e.g., -a in A[-a])
     2. Replace the free abstract index with a specific basis index {componentIndex, -chart}
     3. Apply ToBasis to expand all remaining abstract indices
     4. TraceBasisDummy to sum over dummy indices
     5. Remove Christoffel symbols (0 in flat Minkowski space)
     6. Evaluate metric components numerically
     7. Convert coordinate derivatives CD[{i, -chart}] to explicit Derivative form
  *)

  (* Step 1: Find the free index in the EOM (the uncontracted index) *)
  fieldHead = ExtractFieldHead[field];

  (* Use IndicesOf to get all indices with their up/down information *)
  Module[{allIndices, freeIndices},
    allIndices = List @@ IndicesOf[][eom];
    (* A free index appears without its partner (ChangeIndex) *)
    (* Contracted indices appear as both up and down forms (e.g., a and -a) *)
    freeIndices = Select[allIndices, !MemberQ[allIndices, ChangeIndex[#]] &];
    (* Take the first free index (for vector equation, there should be exactly one) *)
    If[Length[freeIndices] > 0,
      freeIdx = freeIndices[[1]],
      (* Fallback to field template index if no free index found *)
      freeIdx = Cases[field, _Symbol?AbstractIndexQ, {0, Infinity}][[1]]
    ]
  ];

  (* Step 2: Replace free index with specific basis index *)
  (* Use DownIndexQ to determine if it's covariant (down) or contravariant (up) *)
  basisIdx = If[DownIndexQ[freeIdx],
    {componentIndex, -chart},  (* Covariant: use -chart *)
    {componentIndex, chart}    (* Contravariant: use chart *)
  ];
  componentEq = eom /. freeIdx -> basisIdx;

  (* Step 3: Apply ToBasis to convert remaining abstract indices *)
  componentEq = ToBasis[chart][componentEq];

  (* Step 4: Trace over dummy basis indices *)
  componentEq = TraceBasisDummy[componentEq];

  (* Step 5: Remove Christoffel symbols (0 in flat Minkowski space) *)
  componentEq = RemoveChristoffelSymbols[componentEq];
  componentEq = Expand[componentEq];

  (* Step 6: Evaluate epsilon tensor components to numeric ±1 values *)
  (* This handles Chern-Simons and other topological terms with Levi-Civita *)
  componentEq = EvaluateEpsilonComponents[componentEq, chart];
  componentEq = Expand[componentEq];

  (* Step 7: Get coordinate symbols from chart *)
  coordSyms = GetCoordinateSymbols[chart];

  (* Step 8: Evaluate metric components for Minkowski signature (-1, +1) *)
  componentEq = EvaluateMinkowskiMetric[componentEq, chart];
  componentEq = Expand[componentEq];

  (* Step 9: Convert field components to scalar functions *)
  (* Replace A[{i, -chart}] with symbolic component names like A0, A1 *)
  With[{ch = chart, fh = fieldHead, cs = coordSyms},
    componentEq = componentEq /. {
      fh[{i_Integer, -ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs],
      fh[{i_Integer, ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
    }
  ];

  (* Step 10: Convert coordinate derivatives to explicit Derivative form *)
  componentEq = ConvertCDToDerivatives[componentEq, chart];

  (* Expand *)
  componentEq = Expand[componentEq];

  componentEq
];

(* Helper: Find free (non-dummy) indices in an expression *)
(* Renamed to avoid conflict with xPert's FindFreeIndices *)
FindFreeIndicesLocal[expr_] := Module[
  {allIndices, dummyPairs, freeIndices},

  (* Get all indices from the expression *)
  allIndices = Union[Cases[expr, _?AbstractIndexQ, {0, Infinity}]];

  (* Identify dummy pairs (appearing twice with opposite character) *)
  dummyPairs = Select[allIndices,
    Count[expr, #, {0, Infinity}] > 1 &&
    Count[expr, ChangeIndex[#], {0, Infinity}] > 0 &
  ];

  (* Free indices are those not in dummy pairs *)
  freeIndices = Complement[allIndices, dummyPairs, ChangeIndex /@ dummyPairs];

  freeIndices
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
    coefficient = ExtractNumericCoefficient[term, field];
    operator = "laplacian";  (* or "spatial_derivative" for more precision *)
    Return[<|"coefficient" -> coefficient, "operator" -> operator, "field" -> targetField|>]
  ];

  (* Fallback *)
  <|"coefficient" -> 1, "operator" -> "unknown", "field" -> targetField|>
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
