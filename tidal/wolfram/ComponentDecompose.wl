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

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`ComponentDecompose`",
  {"xAct`xTensor`", "xAct`xCoba`",
   "TorsionGertsenshtein`CommonUtilities`"}];

(* Public symbols *)
DecomposeToComponents::usage =
  "DecomposeToComponents[eom, field, chart] decomposes the tensor equation of \
motion into individual scalar component equations. Returns a list of \
{component_index, component_equation} pairs. \
Options: \"ComputeChristoffels\" -> Automatic (default) auto-detects from metric type, \
or True/False for explicit override. \"MetricMatrix\" -> matrix provides an explicit \
metric matrix for curved spacetime evaluation (e.g., Omega^2 * DiagonalMatrix[{-1, 1}]). \
Default is None (flat Minkowski). Auto-detection: constant metrics have Christoffels = 0, \
non-constant metrics trigger explicit Christoffel computation.";

ExtractTensorComponent::usage =
  "ExtractTensorComponent[eom, field, chart, componentIndices, additionalFields, \
computeChristoffels, metricMatrix] extracts a single component equation from a \
tensor EOM. componentIndices is a list of integer indices specifying which component \
to extract: {} for scalar, {i} for vector, {i,j} for rank-2, {i,j,k} for rank-3, etc. \
This is the unified pipeline used by all ranks.";

ReplaceTensorFieldComponents::usage =
  "ReplaceTensorFieldComponents[expr, fieldHead, chart, coordSyms, dim] replaces \
basis-indexed tensor field components with named scalar functions. Works for any \
rank: scalar (rank 0), vector (rank 1), rank-2 tensors (with symmetry), and \
rank 3+ tensors (all components, symmetry reduction not yet implemented).";

ReplaceHigherRankFieldComponents::usage =
  "ReplaceHigherRankFieldComponents[expr, fieldHead, chart, coordSyms, dim] replaces \
basis-indexed rank-3+ tensor components with flat sequential scalar functions. \
Generates replacement rules for all 2^rank index sign configurations per component \
tuple. Example: rank-3 T in dim=2 produces T0..T7 (8 components).";

EnumerateComponentTuples::usage =
  "EnumerateComponentTuples[fieldHead, dim] returns a list of independent component \
index tuples for the given tensor field and dimension. Respects symmetries: \
rank-0 -> {{}}, rank-1 -> {{0},{1},...}, symmetric rank-2 -> upper triangle, \
non-symmetric rank-2 -> all pairs. For rank >= 3, returns all tuples (symmetry \
reduction not yet implemented for rank 3+).";

(* Error messages *)
DecomposeToComponents::badopt =
  "Invalid value for option \"ComputeChristoffels\": `1`. Expected Automatic, True, or False.";


Begin["`Private`"];

(* === Component Decomposition === *)

(* Options for DecomposeToComponents *)
Options[DecomposeToComponents] = {
  "ComputeChristoffels" -> Automatic,  (* Automatic (default), True, or False *)
  "MetricMatrix" -> None  (* Explicit metric matrix for curved spacetime evaluation *)
};

(* 3-arg signature: eom, field, chart (no additional fields, default options) *)
DecomposeToComponents[eom_, field_, chart_] :=
  DecomposeToComponents[eom, field, chart, {}, "ComputeChristoffels" -> Automatic, "MetricMatrix" -> None];

(* 4-arg signature: eom, field, chart, additionalFields (default options) *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List] :=
  DecomposeToComponents[eom, field, chart, additionalFields, "ComputeChristoffels" -> Automatic, "MetricMatrix" -> None];

(* Full signature with options *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List, opts:OptionsPattern[]] := Module[
  {dim, components, indices, componentEq, result, fieldHead, fieldRank, allFieldHeads,
   computeChristoffels, metricMatrix, computeChristoffelsOption, shouldComputeChristoffels,
   coords},

  (* Get option values *)
  computeChristoffelsOption = OptionValue["ComputeChristoffels"];
  metricMatrix = OptionValue["MetricMatrix"];

  (* Auto-detect whether Christoffel computation is needed *)
  coords = GetCoordinateSymbols[chart];
  shouldComputeChristoffels = Which[
    (* Explicit user override: True *)
    computeChristoffelsOption === True, True,

    (* Explicit user override: False *)
    computeChristoffelsOption === False, False,

    (* Automatic detection (default) *)
    computeChristoffelsOption === Automatic,
      If[metricMatrix === None,
        False,  (* No metric → flat Minkowski space *)
        IsNonConstantMetric[metricMatrix, coords]  (* Auto-detect from metric *)
      ],

    (* Fallback for invalid option *)
    True,
      Message[DecomposeToComponents::badopt, "ComputeChristoffels"];
      False
  ];

  (* Use the auto-detected/overridden value for the rest of the pipeline *)
  computeChristoffels = shouldComputeChristoffels;

  (* Get the dimension dynamically from the chart via memoized wrapper *)
  dim = GetChartDimension[chart];

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

(* === Symmetry Reduction Helpers === *)
(* These are package-private helpers for EnumerateComponentTuples.
   Defined at package scope (not inside Module) to avoid Mathematica's
   Module variable localization issues with cross-referencing SetDelayed
   function definitions.

   IMPORTANT: xAct's Cycles is xAct`xPerm`Cycles, NOT System`Cycles.
   All pattern matching and head comparisons must use the fully-qualified
   name to avoid symbol context mismatches. *)

(* Apply a single xAct cycle {a1,a2,...,an} to a tuple:
   maps position a1->a2, a2->a3, ..., an->a1 *)
applyCycleToTuple[tuple_, cycle_List] := Module[{result = tuple, n = Length[cycle]},
  Do[
    result[[ cycle[[Mod[k, n] + 1]] ]] = tuple[[ cycle[[k]] ]],
    {k, 1, n}
  ];
  result
];

(* Apply an xAct Cycles[{i,j,...}, ...] to a tuple.
   Each argument of Cycles is a single cycle list.
   xAct uses Cycles[{i,j}] (single-brace), NOT Mathematica's
   Cycles[{{i,j}}] (double-brace). *)
applyPermToTuple[tuple_, cyc_xAct`xPerm`Cycles] := Module[{result = tuple},
  Do[result = applyCycleToTuple[result, cyc[[k]]], {k, Length[cyc]}];
  result
];

(* Parse an xAct symmetry generator into {sign, xAct`xPerm`Cycles[...]}
   Generators have the form:
     Cycles[{i,j}]    = symmetric swap (sign +1)
     -Cycles[{i,j}]   = antisymmetric swap (sign -1) *)
extractGenSignAndPerm[gen_] := Which[
  Head[gen] === xAct`xPerm`Cycles, {1, gen},
  Head[gen] === Times && gen[[1]] === -1 && Head[gen[[2]]] === xAct`xPerm`Cycles, {-1, gen[[2]]},
  True, {1, xAct`xPerm`Cycles[{}]}  (* identity fallback *)
];

(* Check if a tuple is the canonical representative under given generators *)
isCanonicalTuple[tuple_, gens_List] := AllTrue[gens, Module[
  {sp = extractGenSignAndPerm[#], sign, perm, permuted},
  sign = sp[[1]];
  perm = sp[[2]];
  permuted = applyPermToTuple[tuple, perm];
  Which[
    (* Antisymmetric with equal indices at swapped positions = zero component *)
    sign == -1 && permuted === tuple, False,
    (* Fixed point under this generator *)
    permuted === tuple, True,
    (* Canonical: original <= permuted lexicographically *)
    OrderedQ[{tuple, permuted}], True,
    (* Not canonical: permuted < original *)
    True, False
  ]
] &];

(* === Component Enumeration === *)
(* Returns list of independent component index tuples for any tensor rank *)
(* Respects symmetries: symmetric → upper triangle, antisymmetric → strictly increasing, etc. *)

EnumerateComponentTuples[fieldHead_, dim_] := Module[
  {rank, symGroup, gens, allTuples},

  rank = Length[SlotsOfTensor[fieldHead]];

  Which[
    rank === 0,
      {{}},

    rank === 1,
      Table[{i}, {i, 0, dim - 1}],

    rank >= 2,
      allTuples = Tuples[Range[0, dim - 1], rank];
      symGroup = SymmetryGroupOfTensor[fieldHead];

      If[symGroup === StrongGenSet[{}, GenSet[]],
        (* No symmetry: keep all tuples *)
        allTuples,

        (* Has symmetry: extract generators and filter to canonical reps *)
        gens = List @@ symGroup[[2]];  (* Extract generators from GenSet[...] *)
        Select[allTuples, isCanonicalTuple[#, gens] &]
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

    rank >= 3,
      result = ReplaceHigherRankFieldComponents[result, fh, chart, coordSyms, dim],

    True,
      (* Should not reach here - all ranks >= 0 are handled *)
      Throw[StringJoin[
        "ReplaceTensorFieldComponents: Unexpected rank ", ToString[rank],
        " for field '", ToString[fh], "'."
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

(* Helper: Replace rank-N (N>=3) field basis indices with flat sequential scalar functions *)
(* Follows the same pattern as ReplaceRank2FieldComponents but for arbitrary rank *)
(* For rank-3 non-symmetric T in dim=2: T[{0,-ch},{0,-ch},{0,-ch}] -> T0[t,x], etc. *)
ReplaceHigherRankFieldComponents[expr_, fh_, chart_, coordSyms_, dim_] := Module[
  {result = expr, rank, tuples},

  rank = Length[SlotsOfTensor[fh]];
  tuples = EnumerateComponentTuples[fh, dim];

  With[{ch = chart, cs = coordSyms},
    Do[
      Module[{tuple = tuples[[k]], seqIdx = k - 1, sym, indexConfigs},
        sym = Symbol[ToString[fh] <> ToString[seqIdx]];

        (* Generate all 2^rank index sign configurations (chart/-chart per slot) *)
        indexConfigs = Tuples[{ch, -ch}, rank];
        Do[
          Module[{pattern},
            pattern = Table[{tuple[[n]], config[[n]]}, {n, rank}];
            result = result /. {fh @@ pattern :> sym[Sequence @@ cs]}
          ],
          {config, indexConfigs}
        ]
      ],
      {k, 1, Length[tuples]}
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
