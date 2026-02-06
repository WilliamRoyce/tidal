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

  (* For a vector field A_a, extract each component *)
  If[fieldRank === 1,
    result = Table[
      {
        idx,
        ExtractVectorComponent[eom, field, chart, idx, additionalFields, computeChristoffels, metricMatrix]
      },
      {idx, 0, dim - 1}
    ];
    Return[result]
  ];

  (* For a rank-2 tensor field h[-a,-b], extract each independent component *)
  If[fieldRank === 2,
    Module[{symQ, componentPairs},
      (* Check for non-trivial symmetry (symmetric or antisymmetric) *)
      symQ = (SymmetryGroupOfTensor[fieldHead] =!= StrongGenSet[{}, GenSet[]]);

      componentPairs = If[symQ,
        (* Symmetric: upper triangle {i,j} with j >= i *)
        Flatten[Table[{i, j}, {i, 0, dim - 1}, {j, i, dim - 1}], 1],
        (* Non-symmetric: all pairs *)
        Flatten[Table[{i, j}, {i, 0, dim - 1}, {j, 0, dim - 1}], 1]
      ];

      result = Table[
        {
          idx - 1,
          ExtractRank2Component[eom, field, chart,
            componentPairs[[idx, 1]], componentPairs[[idx, 2]],
            additionalFields, computeChristoffels, metricMatrix]
        },
        {idx, 1, Length[componentPairs]}
      ];
      Return[result]
    ]
  ];

  (* Rank 3+: fail explicitly *)
  Throw[StringJoin[
    "DecomposeToComponents: Rank-", ToString[fieldRank],
    " tensor decomposition is not supported. ",
    "Only scalar (rank 0), vector (rank 1), and rank-2 tensor fields are implemented. ",
    "Field: ", ToString[fieldHead]
  ]]
];

(* === Vector Component Extraction === *)

ExtractVectorComponent[eom_, field_, chart_, componentIndex_, additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None] := Module[
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

  (* Step 5: Evaluate Christoffel symbols (0 for flat, or computed from metric for curved) *)
  componentEq = EvaluateChristoffelComponents[componentEq, chart, computeChristoffels];
  componentEq = Expand[componentEq];

  (* Step 6: Evaluate epsilon tensor components to numeric ±1 values *)
  (* This handles Chern-Simons and other topological terms with Levi-Civita *)
  componentEq = EvaluateEpsilonComponents[componentEq, chart];
  componentEq = Expand[componentEq];

  (* Step 7: Get coordinate symbols from chart *)
  coordSyms = GetCoordinateSymbols[chart];

  (* Step 8: Evaluate metric components *)
  If[metricMatrix =!= None,
    componentEq = EvaluateMetricComponents[componentEq, chart, metricMatrix],
    componentEq = EvaluateMinkowskiMetric[componentEq, chart]
  ];
  componentEq = Expand[componentEq];

  (* Step 9: Convert field components to scalar functions *)
  (* Replace A[{i, -chart}] with symbolic component names like A0, A1 *)
  With[{ch = chart, fh = fieldHead, cs = coordSyms},
    componentEq = componentEq /. {
      fh[{i_Integer, -ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs],
      fh[{i_Integer, ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
    }
  ];

  (* Step 9b: Transform additional coupled fields to coordinate form *)
  (* Handles both scalar and vector additional fields *)
  Do[
    Module[{afh = ExtractFieldHead[af], afRank},
      afRank = Length[SlotsOfTensor[afh]];
      If[afRank == 0,
        (* Scalar additional field: replace afh[] with afh0[coords...] *)
        componentEq = componentEq /. {
          afh[] :> Symbol[ToString[afh] <> "0"][Sequence @@ coordSyms]
        },
        (* Vector additional field: replace afh[{i, ±chart}] with afhi[coords...] *)
        With[{ch = chart, cs = coordSyms},
          componentEq = componentEq /. {
            afh[{i_Integer, -ch}] :> Symbol[ToString[afh] <> ToString[Abs[i]]][Sequence @@ cs],
            afh[{i_Integer, ch}] :> Symbol[ToString[afh] <> ToString[Abs[i]]][Sequence @@ cs]
          }
        ]
      ]
    ],
    {af, additionalFields}
  ];

  (* Step 10: Convert coordinate derivatives to explicit Derivative form *)
  componentEq = ConvertCDToDerivatives[componentEq, chart];

  (* Expand *)
  componentEq = Expand[componentEq];

  componentEq
];

(* === Rank-2 Tensor Component Extraction === *)

ExtractRank2Component[eom_, field_, chart_, idx1_, idx2_,
  additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None] := Module[
  {componentEq, fieldHead, coordSyms, dim, allFieldHeads},

  (*
     For a rank-2 equation EOM[-a,-b] = 0, extract the (idx1,idx2) component.

     Strategy is the same as ExtractVectorComponent but with TWO free indices:
     1. Find two free indices in the EOM
     2. Replace them with specific basis indices {idx1, -chart} and {idx2, -chart}
     3-8. Same pipeline: ToBasis → TraceBasisDummy → Christoffels → Epsilon → Metric
     9. Replace rank-2 field components with scalar functions (flat sequential naming)
     10. Transform additional fields (scalar, vector, or rank-2)
     11. ConvertCDToDerivatives
  *)

  fieldHead = ExtractFieldHead[field];

  (* Step 1: Find the two free indices *)
  Module[{allIndices, freeIdxList, basisIdx1, basisIdx2},
    allIndices = List @@ IndicesOf[][eom];
    freeIdxList = Select[allIndices, !MemberQ[allIndices, ChangeIndex[#]] &];
    If[Length[freeIdxList] < 2,
      Throw[StringJoin[
        "ExtractRank2Component: Expected at least 2 free indices, found ",
        ToString[Length[freeIdxList]], " in expression: ", ToString[Short[eom]]
      ]]
    ];

    (* Step 2: Replace free indices with specific basis indices *)
    basisIdx1 = If[DownIndexQ[freeIdxList[[1]]],
      {idx1, -chart}, {idx1, chart}];
    basisIdx2 = If[DownIndexQ[freeIdxList[[2]]],
      {idx2, -chart}, {idx2, chart}];
    componentEq = eom /. {
      freeIdxList[[1]] -> basisIdx1,
      freeIdxList[[2]] -> basisIdx2
    };
  ];

  (* Step 3: Apply ToBasis *)
  componentEq = ToBasis[chart][componentEq];

  (* Step 4: TraceBasisDummy *)
  componentEq = TraceBasisDummy[componentEq];

  (* Step 5: Evaluate Christoffel symbols *)
  componentEq = EvaluateChristoffelComponents[componentEq, chart, computeChristoffels];
  componentEq = Expand[componentEq];

  (* Step 6: Evaluate epsilon tensors *)
  componentEq = EvaluateEpsilonComponents[componentEq, chart];
  componentEq = Expand[componentEq];

  (* Step 7: Get coordinate symbols *)
  coordSyms = GetCoordinateSymbols[chart];
  dim = Length[coordSyms];

  (* Step 8: Evaluate metric components *)
  If[metricMatrix =!= None,
    componentEq = EvaluateMetricComponents[componentEq, chart, metricMatrix],
    componentEq = EvaluateMinkowskiMetric[componentEq, chart]
  ];
  componentEq = Expand[componentEq];

  If[metricMatrix =!= None,
    componentEq = EvaluatePDMetric[componentEq, chart, metricMatrix];
    componentEq = Expand[componentEq]
  ];

  (* Step 9: Convert all field components to scalar functions *)
  allFieldHeads = Join[{fieldHead}, ExtractFieldHead /@ additionalFields];

  Do[
    Module[{fh = afh, fhRank = Length[SlotsOfTensor[afh]]},
      Which[
        fhRank === 0,
          componentEq = componentEq /. {
            fh[] :> Symbol[ToString[fh] <> "0"][Sequence @@ coordSyms]
          },

        fhRank === 1,
          With[{ch = chart, cs = coordSyms},
            componentEq = componentEq /. {
              fh[{i_Integer, -ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs],
              fh[{i_Integer, ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
            }
          ],

        fhRank === 2,
          componentEq = ReplaceRank2FieldComponents[componentEq, fh, chart, coordSyms, dim],

        True,
          Throw["ExtractRank2Component: Unsupported additional field rank " <> ToString[fhRank]]
      ]
    ],
    {afh, allFieldHeads}
  ];

  (* Step 10: Convert coordinate derivatives *)
  componentEq = ConvertCDToDerivatives[componentEq, chart];
  Expand[componentEq]
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
