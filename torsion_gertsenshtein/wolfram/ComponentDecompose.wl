(* ::Package:: *)
(*
   MODULE: ComponentDecompose.wl
   PURPOSE: Decompose tensor equations into scalar component equations

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus)
     - xAct`xCoba` (coordinate bases, ToBasis, TraceBasisDummy)
     - TorsionGertsenshtein`CommonUtilities` (CD conversion, epsilon evaluation)

   DATA FLOW:
     Tensor EOM (abstract indices: A_a, F_{ab})
       → ToBasis (convert to chart basis)
       → TraceBasisDummy (sum over dummy indices)
       → EvaluateChristoffelComponents (flat space: Γ = 0, curved: compute from metric)
       → EvaluateEpsilonComponents (ε tensors → ±1)
       → ConvertCDToDerivatives (CD → Derivative)
       → Extract components (A_0[t,x], A_1[t,x], ...)

   KEY FEATURES:
     - Supports scalar, vector, and rank-2 tensor fields
     - additionalFields parameter for cross-field coupling (e.g., pass {chi[]}
       when decomposing phi to ensure chi transforms to coordinate form)
     - Dimension-agnostic: uses GetChartDimension instead of hardcoded values
     - Automatic epsilon tensor evaluation for Chern-Simons and other topological terms

   OUTPUT FORMAT:
     List of {component_index, component_equation} pairs, where equations are
     in Derivative[dt, dx, ...][field][t, x, ...] form ready for ExportJSON.

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`ComponentDecompose`",
  {"xAct`xTensor`", "xAct`xCoba`",
   "TorsionGertsenshtein`CommonUtilities`"}];

(* Public symbols *)
DecomposeToComponents::usage =
  "DecomposeToComponents[eom, field, chart] decomposes the tensor equation of \
motion into individual scalar component equations. Returns a list of \
{component_index, component_equation} pairs. \
Options: \"ComputeChristoffels\" -> True computes Christoffel symbols from the \
metric definition using xAct (for curved spacetimes). Default is False (flat space). \
\"MetricMatrix\" -> matrix provides an explicit metric matrix for curved spacetime \
evaluation (e.g., Omega^2 * DiagonalMatrix[{-1, 1}]). Default is None (flat Minkowski).";

ExtractTensorComponent::usage =
  "ExtractTensorComponent[eom, field, chart, componentIndices, additionalFields, \
computeChristoffels, metricMatrix] extracts a single component equation from a \
tensor EOM. componentIndices is a list of integer indices specifying which component \
to extract: {} for scalar, {i} for vector, {i,j} for rank-2, {i,j,k} for rank-3, etc. \
This is the unified pipeline used by all ranks.";

ReplaceTensorFieldComponents::usage =
  "ReplaceTensorFieldComponents[expr, fieldHead, chart, coordSyms, dim] replaces \
basis-indexed tensor field components with named scalar functions. Works for any \
rank: scalar (rank 0), vector (rank 1), rank-2 tensors. Extensible to higher ranks \
by adding replacement patterns for the desired rank.";

EnumerateComponentTuples::usage =
  "EnumerateComponentTuples[fieldHead, dim] returns a list of independent component \
index tuples for the given tensor field and dimension. Respects symmetries: \
rank-0 -> {{}}, rank-1 -> {{0},{1},...}, symmetric rank-2 -> upper triangle, \
non-symmetric rank-2 -> all pairs. For rank >= 3, returns all tuples (symmetry \
reduction not yet implemented for rank 3+).";


Begin["`Private`"];

(* === Component Decomposition === *)

(* Options for DecomposeToComponents *)
Options[DecomposeToComponents] = {
  "ComputeChristoffels" -> False,  (* Set True for curved spacetimes *)
  "MetricMatrix" -> None  (* Explicit metric matrix for curved spacetime evaluation *)
};

(* 3-arg signature: eom, field, chart (no additional fields, default options) *)
DecomposeToComponents[eom_, field_, chart_] :=
  DecomposeToComponents[eom, field, chart, {}, "ComputeChristoffels" -> False, "MetricMatrix" -> None];

(* 4-arg signature: eom, field, chart, additionalFields (default options) *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List] :=
  DecomposeToComponents[eom, field, chart, additionalFields, "ComputeChristoffels" -> False, "MetricMatrix" -> None];

(* Full signature with options *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List, opts:OptionsPattern[]] := Module[
  {dim, components, indices, componentEq, result, fieldHead, fieldRank, allFieldHeads,
   computeChristoffels, metricMatrix},

  (* Get option values *)
  computeChristoffels = OptionValue["ComputeChristoffels"];
  metricMatrix = OptionValue["MetricMatrix"];

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

    (* For curved spacetime: use xAct's ChangeCovD + ChristoffelToGradMetric
       to replace CovD and Christoffel symbols with PD and metric derivatives.
       This must happen BEFORE ToBasis, on the abstract tensor expression. *)
    If[computeChristoffels === True && metricMatrix =!= None,
      (* Detect the covariant derivative operator from the equation using xAct's CovDQ *)
      Module[{covdOps, covdOp},
        covdOps = Cases[eom, (f_)[_][_] /; CovDQ[f] :> f, {0, Infinity}] // DeleteDuplicates;
        If[Length[covdOps] > 0,
          covdOp = covdOps[[1]];
          componentEq = ExpandChristoffelsToMetricDerivatives[eom, covdOp, chart],
          componentEq = eom
        ]
      ],
      componentEq = eom
    ];

    componentEq = ToBasis[chart][componentEq];
    componentEq = TraceBasisDummy[componentEq];

    (* For flat spacetime: set all Christoffel symbols to 0 *)
    If[computeChristoffels =!= True,
      componentEq = EvaluateChristoffelComponents[componentEq, chart, False]
    ];
    componentEq = Expand[componentEq];

    (* Evaluate epsilon tensor components to numeric ±1 values *)
    (* This handles Chern-Simons and other topological terms with Levi-Civita *)
    componentEq = EvaluateEpsilonComponents[componentEq, chart];
    componentEq = Expand[componentEq];

    (* Get coordinate symbols *)
    Module[{coordSyms},
      coordSyms = GetCoordinateSymbols[chart];

      (* Evaluate metric components *)
      (* If MetricMatrix is provided, use actual metric values (curved spacetime) *)
      (* Otherwise, use flat Minkowski signature (-1, +1, +1, ...) *)
      If[metricMatrix =!= None,
        componentEq = EvaluateMetricComponents[componentEq, chart, metricMatrix],
        componentEq = EvaluateMinkowskiMetric[componentEq, chart]
      ];
      componentEq = Expand[componentEq];

      (* For curved spacetime: evaluate partial derivatives of metric components *)
      If[metricMatrix =!= None,
        componentEq = EvaluatePDMetric[componentEq, chart, metricMatrix];
        componentEq = Expand[componentEq]
      ];

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

  (* For any tensor field of rank >= 1, use the unified pipeline *)
  If[fieldRank >= 1,
    Module[{componentTuples},
      componentTuples = EnumerateComponentTuples[fieldHead, dim];
      result = Table[
        {
          idx - 1,
          ExtractTensorComponent[eom, field, chart,
            componentTuples[[idx]], additionalFields, computeChristoffels, metricMatrix]
        },
        {idx, 1, Length[componentTuples]}
      ];
      Return[result]
    ]
  ]
];

(* === Unified Tensor Component Extraction === *)
(* Single pipeline for any tensor rank. The old rank-specific functions
   (ExtractVectorComponent, ExtractRank2Component) delegate here. *)

ExtractTensorComponent[eom_, field_, chart_, componentIndices_List,
  additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None] := Module[
  {componentEq, fieldHead, coordSyms, dim, allFieldHeads, rank},

  (*
     Unified extraction pipeline for any-rank tensor EOM:
     1. Find N free indices (N = Length[componentIndices]) and fix to basis values
     2. ToBasis: convert remaining abstract indices to basis
     3. TraceBasisDummy: sum over dummy basis indices
     4. EvaluateChristoffelComponents: Γ = 0 for flat, computed for curved
     5. EvaluateCurvatureComponents: R = 0 for constant metric, computed otherwise
     6. EvaluateEpsilonComponents: ε → ±1
     7. Evaluate metric components
     8. Replace all tensor fields with named scalar functions
     9. ConvertCDToDerivatives: CD → Derivative form
  *)

  fieldHead = ExtractFieldHead[field];
  rank = Length[componentIndices];

  (* Step 1: Find free indices and fix them to componentIndices *)
  Module[{allIndices, freeIdxList, replacements},
    allIndices = List @@ IndicesOf[][eom];
    freeIdxList = Select[allIndices, !MemberQ[allIndices, ChangeIndex[#]] &];

    (* Fallback: extract from field template if expression has no free indices *)
    If[Length[freeIdxList] < rank,
      freeIdxList = Cases[field, _Symbol?AbstractIndexQ, {0, Infinity}]
    ];

    If[Length[freeIdxList] < rank,
      Throw[StringJoin[
        "ExtractTensorComponent: Expected at least ", ToString[rank],
        " free indices, found ", ToString[Length[freeIdxList]],
        " in expression: ", ToString[Short[eom]]
      ]]
    ];

    (* Fix each free index to its corresponding basis value *)
    replacements = Table[
      freeIdxList[[i]] -> If[DownIndexQ[freeIdxList[[i]]],
        {componentIndices[[i]], -chart},
        {componentIndices[[i]], chart}
      ],
      {i, rank}
    ];
    componentEq = eom /. replacements;
  ];

  (* Step 2: ToBasis *)
  componentEq = ToBasis[chart][componentEq];

  (* Step 3: TraceBasisDummy *)
  componentEq = TraceBasisDummy[componentEq];

  (* Step 4: Evaluate Christoffel symbols *)
  componentEq = EvaluateChristoffelComponents[componentEq, chart, computeChristoffels];
  componentEq = Expand[componentEq];

  (* Step 5: Evaluate background curvature tensors from the metric *)
  (* For constant metrics this computes R = 0; for non-constant it leaves them *)
  componentEq = EvaluateCurvatureComponents[componentEq, chart,
    If[metricMatrix =!= None, metricMatrix, None]];
  componentEq = Expand[componentEq];

  (* Step 6: Evaluate epsilon tensors *)
  componentEq = EvaluateEpsilonComponents[componentEq, chart];
  componentEq = Expand[componentEq];

  (* Step 7: Get coordinate symbols and evaluate metric *)
  coordSyms = GetCoordinateSymbols[chart];
  dim = Length[coordSyms];

  If[metricMatrix =!= None,
    componentEq = EvaluateMetricComponents[componentEq, chart, metricMatrix],
    componentEq = EvaluateMinkowskiMetric[componentEq, chart]
  ];
  componentEq = Expand[componentEq];

  If[metricMatrix =!= None,
    componentEq = EvaluatePDMetric[componentEq, chart, metricMatrix];
    componentEq = Expand[componentEq]
  ];

  (* Step 8: Replace ALL tensor fields with named scalar functions *)
  allFieldHeads = Join[{fieldHead}, ExtractFieldHead /@ additionalFields];
  Do[
    componentEq = ReplaceTensorFieldComponents[componentEq, afh, chart, coordSyms, dim],
    {afh, allFieldHeads}
  ];

  (* Step 9: Convert coordinate derivatives to Derivative form *)
  componentEq = ConvertCDToDerivatives[componentEq, chart];

  Expand[componentEq]
];

(* === Backward-Compatible Delegates === *)
(* These call the unified ExtractTensorComponent so existing code still works *)

ExtractVectorComponent[eom_, field_, chart_, componentIndex_,
  additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None] :=
  ExtractTensorComponent[eom, field, chart, {componentIndex},
    additionalFields, computeChristoffels, metricMatrix];

ExtractRank2Component[eom_, field_, chart_, idx1_, idx2_,
  additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None] :=
  ExtractTensorComponent[eom, field, chart, {idx1, idx2},
    additionalFields, computeChristoffels, metricMatrix];

(* === Component Enumeration === *)
(* Returns list of independent component index tuples for any tensor rank *)
(* Respects symmetries where implemented *)

EnumerateComponentTuples[fieldHead_, dim_] := Module[
  {rank, symQ},

  rank = Length[SlotsOfTensor[fieldHead]];

  Which[
    rank === 0,
      {{}},

    rank === 1,
      Table[{i}, {i, 0, dim - 1}],

    rank >= 2,
      (* Check for non-trivial symmetry (symmetric or antisymmetric) *)
      symQ = (SymmetryGroupOfTensor[fieldHead] =!= StrongGenSet[{}, GenSet[]]);
      If[rank === 2,
        If[symQ,
          (* Symmetric rank-2: upper triangle *)
          Flatten[Table[{i, j}, {i, 0, dim - 1}, {j, i, dim - 1}], 1],
          (* Non-symmetric rank-2: all pairs *)
          Flatten[Table[{i, j}, {i, 0, dim - 1}, {j, 0, dim - 1}], 1]
        ],
        (* Rank 3+: all tuples (symmetry reduction not yet implemented) *)
        (* To add symmetry support for rank 3+, use SymmetryGroupOfTensor *)
        (* to identify independent components and eliminate redundant ones *)
        Tuples[Range[0, dim - 1], rank]
      ],

    True,
      {{}}
  ]
];

(* === Rank-Generic Field Component Replacement === *)
(* Replaces basis-indexed tensor fields with named scalar functions *)
(* Extensible: add a new rank branch to support higher-rank tensors *)

ReplaceTensorFieldComponents[expr_, fh_, chart_, coordSyms_, dim_] := Module[
  {rank, result = expr},

  rank = Length[SlotsOfTensor[fh]];

  Which[
    rank === 0,
      result = result /. {
        fh[] :> Symbol[ToString[fh] <> "0"][Sequence @@ coordSyms]
      },

    rank === 1,
      With[{ch = chart, cs = coordSyms},
        result = result /. {
          fh[{i_Integer, -ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs],
          fh[{i_Integer, ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
        }
      ],

    rank === 2,
      result = ReplaceRank2FieldComponents[result, fh, chart, coordSyms, dim],

    True,
      (* Clear extension point for higher ranks *)
      Throw[StringJoin[
        "ReplaceTensorFieldComponents: Rank-", ToString[rank],
        " field replacement is not yet implemented for field '", ToString[fh], "'. ",
        "To add support: extend this function with component enumeration and ",
        "replacement rules for rank-", ToString[rank], " tensors, following the ",
        "pattern of ReplaceRank2FieldComponents."
      ]]
  ];

  result
];

(* Helper: Replace rank-2 field basis indices with flat sequential scalar functions *)
(* For symmetric h in dim=2: h[{0,-ch},{0,-ch}] -> h0[t,x], h[{0,-ch},{1,-ch}] -> h1[t,x], h[{1,-ch},{1,-ch}] -> h2[t,x] *)
ReplaceRank2FieldComponents[expr_, fh_, chart_, coordSyms_, dim_] := Module[
  {result = expr, symQ, pairs},

  symQ = (SymmetryGroupOfTensor[fh] =!= StrongGenSet[{}, GenSet[]]);
  pairs = If[symQ,
    Flatten[Table[{ii, jj}, {ii, 0, dim - 1}, {jj, ii, dim - 1}], 1],
    Flatten[Table[{ii, jj}, {ii, 0, dim - 1}, {jj, 0, dim - 1}], 1]
  ];

  With[{ch = chart, cs = coordSyms},
    Do[
      Module[{pair = pairs[[k]], seqIdx = k - 1, sym},
        sym = Symbol[ToString[fh] <> ToString[seqIdx]];

        (* All index configurations: down-down, up-up, mixed *)
        result = result /. {
          fh[{pair[[1]], -ch}, {pair[[2]], -ch}] :> sym[Sequence @@ cs],
          fh[{pair[[1]], ch}, {pair[[2]], ch}] :> sym[Sequence @@ cs],
          fh[{pair[[1]], ch}, {pair[[2]], -ch}] :> sym[Sequence @@ cs],
          fh[{pair[[1]], -ch}, {pair[[2]], ch}] :> sym[Sequence @@ cs]
        };

        (* For symmetric: also map swapped indices to same component *)
        If[symQ && pair[[1]] =!= pair[[2]],
          result = result /. {
            fh[{pair[[2]], -ch}, {pair[[1]], -ch}] :> sym[Sequence @@ cs],
            fh[{pair[[2]], ch}, {pair[[1]], ch}] :> sym[Sequence @@ cs],
            fh[{pair[[2]], ch}, {pair[[1]], -ch}] :> sym[Sequence @@ cs],
            fh[{pair[[2]], -ch}, {pair[[1]], ch}] :> sym[Sequence @@ cs]
          }
        ]
      ],
      {k, 1, Length[pairs]}
    ]
  ];

  result
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

End[];
EndPackage[];
