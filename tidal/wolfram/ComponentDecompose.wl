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
       → SeparateFieldMetrics (undo ContractMetric: V^a → g^{ab} V_{-b})
         [hoisted: applied once per EOM before component loop, not per component]
       → ToBasis (convert to chart basis)
       → **Batched** TraceBasisDummy + Expand + metric evaluation
         (fused to prevent O(dim^{2*rank}) intermediate memory blowup)
       → EvaluateChristoffelComponents (flat: Γ = 0, curved: compute from metric)
       → EvaluateCurvatureComponents (constant metric: R = 0, non-constant: compute)
       → EvaluateEpsilonComponents (ε tensors → ±1)
       → Evaluate remaining metric components (idempotent cleanup)
       → Replace tensor fields with scalar functions
       → ConvertCDToDerivatives (CD → Derivative)

   MEMORY OPTIMIZATION (batched TraceBasisDummy + metric fusion):
     TraceBasisDummy expands dummy index sums, producing O(dim^{2*n_dummy})
     terms per original term. Most off-diagonal metric components are zero
     (flat) or sparse (curved). By fusing TraceBasisDummy + Expand + metric
     evaluation into batches of ~50 input terms, peak memory is bounded by a
     single batch (~800 terms) rather than the full expansion (~970K terms).
     For Einstein-Maxwell cross-coupling in 3+1D flat spacetime:
     568 input terms → 970K expanded → ~200 after metric (99.97% reduction).
     Without batching, the a-field decomposition peaks at 7.6+ GB and OOMs.

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

DecomposeScalarExpression::usage =
  "DecomposeScalarExpression[expr, chart, allFieldHeads, opts] converts a scalar \
tensor expression (e.g., a Hamiltonian density) to component form using the same \
pipeline as DecomposeToComponents. allFieldHeads is a list of ALL tensor field \
heads that appear in the expression (primary + additional + background). \
The result is a single scalar expression in Derivative[...][field][coords] form. \
Options: \"ComputeChristoffels\" -> Automatic, \"MetricMatrix\" -> None.";

SeparateFieldMetrics::usage =
  "SeparateFieldMetrics[expr, chart] applies xAct's SeparateMetric to undo metric \
contractions in tensor expressions. Converts V^a (= g^{ab} V_{-b}) back to explicit \
g^{ab} * V_{-b}, ensuring all field tensor indices are in their canonical (covariant) \
form before ToBasis decomposition. This prevents sign errors when mapping tensor \
components to scalar functions (e.g., V^0 = -V_0 for Minkowski temporal component). \
The metric is automatically extracted from the covariant derivative in the expression.";

ValidateNoUnresolvedBackgrounds::usage =
  "ValidateNoUnresolvedBackgrounds[expr, bgHeads] checks that no abstract background \
field symbols remain unresolved in the component expression.";

BatchedTraceBasisDummyWithMetric::usage =
  "BatchedTraceBasisDummyWithMetric[componentEq, chart, metricMatrix, batchSize, \
backgroundFieldRules] performs TraceBasisDummy in batches to prevent O(dim^{2n}) \
term blowup during Expand. Each batch of ~batchSize terms is traced, expanded, \
and metric-evaluated before accumulating, keeping peak memory bounded. Without \
batching, Einstein-Maxwell cross-coupling exceeds 7 GB; with batching, <1 GB. \
Optional backgroundFieldRules evaluates background field DownValues during the \
fused loop for further memory reduction.";

(* Error messages *)
DecomposeToComponents::badopt =
  "Invalid value for option \"ComputeChristoffels\": `1`. Expected Automatic, True, or False.";


Begin["`Private`"];

(* === Metric Separation for Correct Index Positions === *)
(*
   After EulerLagrange applies ContractMetric, field tensors may have raised indices:
     V^a = g^{ab} V_{-b}  (metric absorbed into field)
   This causes sign errors when mapping to scalars: V^0 = -V_0 for Minkowski.

   SeparateFieldMetrics undoes this by applying xAct's SeparateMetric, making the
   metric contraction explicit: V^a → g^{ab} * V_{-b}. After ToBasis, fields then
   have only covariant (canonical/defined) indices, and the metric factors evaluate
   to numerical values in EvaluateMetricComponents.
*)

SeparateFieldMetrics[expr_, chart_] := Module[
  {covdOps, metric},

  (* Extract the covariant derivative operator from the expression *)
  covdOps = Cases[expr, (f_)[_][_] /; CovDQ[f] :> f, {0, Infinity}] // DeleteDuplicates;
  If[Length[covdOps] == 0, Return[expr]];

  (* Get the metric associated with the CovD via xAct introspection *)
  metric = MetricOfCovD[covdOps[[1]]];

  (* SeparateMetric: V^a → g^{ab} V_{-b}, making metric contraction explicit *)
  (* This is the inverse of ContractMetric applied in EulerLagrange.wl *)
  SeparateMetric[metric][expr]
];


(* === Component Decomposition === *)

(* Options for DecomposeToComponents *)
Options[DecomposeToComponents] = {
  "ComputeChristoffels" -> Automatic,  (* Automatic (default), True, or False *)
  "MetricMatrix" -> None,  (* Explicit metric matrix for curved spacetime evaluation *)
  "SkipTuples" -> {},  (* Component index tuples to skip (e.g. TT-zeroed {0,mu}) *)
  "BackgroundFieldRules" -> {},  (* List of {fieldHead, {comp0, comp1, ...}} for background field eval *)
  "KilledAxes" -> {},  (* Basis indices of transverse axes for plane-wave reduction *)
  "PertFieldHeads" -> {}  (* Perturbation field head symbols whose transverse PD should be zeroed *)
};

(* 3-arg signature: eom, field, chart (no additional fields, default options) *)
DecomposeToComponents[eom_, field_, chart_] :=
  DecomposeToComponents[eom, field, chart, {}, "ComputeChristoffels" -> Automatic, "MetricMatrix" -> None, "SkipTuples" -> {}, "BackgroundFieldRules" -> {}];

(* 4-arg signature: eom, field, chart, additionalFields (default options) *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List] :=
  DecomposeToComponents[eom, field, chart, additionalFields, "ComputeChristoffels" -> Automatic, "MetricMatrix" -> None, "SkipTuples" -> {}, "BackgroundFieldRules" -> {}];

(* Full signature with options *)
DecomposeToComponents[eom_, field_, chart_, additionalFields_List, opts:OptionsPattern[]] := Module[
  {dim, components, indices, componentEq, result, fieldHead, fieldRank, allFieldHeads,
   computeChristoffels, metricMatrix, computeChristoffelsOption, shouldComputeChristoffels,
   coords, backgroundFieldRules},

  (* Get option values *)
  computeChristoffelsOption = OptionValue["ComputeChristoffels"];
  metricMatrix = OptionValue["MetricMatrix"];
  backgroundFieldRules = OptionValue["BackgroundFieldRules"];

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

    (* Separate metric contractions before ToBasis *)
    (* Ensures cross-field vector/tensor terms have covariant indices *)
    componentEq = SeparateFieldMetrics[componentEq, chart];

    (* Apply ToBasis term-by-term to avoid xperm segfault on large sums *)
    If[Head[componentEq] === Plus,
      componentEq = Total[ToBasis[chart] /@ List @@ componentEq],
      componentEq = ToBasis[chart][componentEq]
    ];
    (* TraceBasisDummy also term-by-term for the same reason *)
    If[Head[componentEq] === Plus,
      componentEq = Total[TraceBasisDummy /@ List @@ componentEq],
      componentEq = TraceBasisDummy[componentEq]
    ];

    (* For flat spacetime: set all Christoffel symbols to 0 *)
    If[computeChristoffels =!= True,
      componentEq = EvaluateChristoffelComponents[componentEq, chart, False]
    ];
    componentEq = Expand[componentEq];

    (* Evaluate epsilon tensor components to numeric ±1 values *)
    (* This handles Chern-Simons and other topological terms with Levi-Civita *)
    (* Pass metricMatrix for correct volume factor and index raising in curved spacetimes *)
    componentEq = EvaluateEpsilonComponents[componentEq, chart, metricMatrix];
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

      (* Evaluate background field components and their PD derivatives *)
      If[backgroundFieldRules =!= {},
        Do[
          componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
          {bg, backgroundFieldRules}
        ];
        componentEq = Expand[componentEq]
      ];

      (* Replace ALL fields (any rank) with functions of coordinates *)
      (* This ensures cross-field terms are properly transformed *)
      (* Uses ReplaceTensorFieldComponents which dispatches by rank: *)
      (*   rank 0: fh[] -> fh0[t,x,y]  *)
      (*   rank 1: fh[{i,-chart}] -> fhi[t,x,y]  *)
      (*   rank 2+: full component replacement *)
      (* Pass metricMatrix so rank-2 fallback uses correct metric weights for curved metrics *)
      Do[
        componentEq = ReplaceTensorFieldComponents[componentEq, afh, chart, coordSyms, dim, metricMatrix],
        {afh, allFieldHeads}
      ];

      (* Convert coordinate derivatives to explicit Derivative form *)
      componentEq = ConvertCDToDerivatives[componentEq, chart];

      (* Post-ConvertCDToDerivatives: catch any residual Derivative[...][bgField][...] forms *)
      If[backgroundFieldRules =!= {},
        Do[
          componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
          {bg, backgroundFieldRules}
        ];
        componentEq = Expand[componentEq]
      ];
    ];

    (* Expand to get explicit Derivative[...] form *)
    componentEq = Expand[componentEq];
    Return[{{0, componentEq}}]
  ];

  (* For any tensor field of rank >= 1, use the unified pipeline *)
  If[fieldRank >= 1,
    Module[{allTuples, componentTuples, skipTuples, flatIdxMap, eomSep},
      allTuples = EnumerateComponentTuples[fieldHead, dim];
      (* Build flat index map: tuple -> original 0-based position *)
      flatIdxMap = Association[Table[allTuples[[i]] -> i - 1, {i, Length[allTuples]}]];
      (* Filter out tuples known to be zero (e.g. TT gauge h_{0,mu}) *)
      skipTuples = OptionValue["SkipTuples"];
      componentTuples = allTuples;
      If[skipTuples =!= {},
        componentTuples = Complement[allTuples, skipTuples];
        Print["SkipTuples: skipping ", Length[allTuples] - Length[componentTuples],
              " components, ", Length[componentTuples], " remaining"]
      ];
      (* Pre-apply SeparateFieldMetrics once on the full EOM before extracting components.
         This converts V^a → g^{ab} V_{-b} (undoes ContractMetric) at abstract-tensor level.
         Safe to hoist: SeparateFieldMetrics acts on dummy/contracted indices, not free indices,
         so it commutes with the per-component free-index replacements in ExtractTensorComponent
         step 1. Hoisting eliminates O(rank * dim^rank) redundant calls.
         CRITICAL: Must use a local variable (eomSep), NOT reassign the function parameter
         `eom`. Reassigning a pattern variable like `eom = ...` tries to set a DownValue on
         the expression value (Set::write: Tag Times/Plus is Protected), silently failing.
         This caused VarD's upper-index E^{ab} to bypass SeparateMetric, breaking the
         metric contraction in the box operator and producing wrong-sign equations. *)
      eomSep = SeparateFieldMetrics[eom, chart];

      (* Use Do+AppendTo instead of Table so Share[] can reclaim memory
         between component extractions — critical for cross-field coupling
         cases like Einstein-Maxwell where expressions grow large. *)
      result = {};
      Do[
        AppendTo[result,
          {flatIdxMap[componentTuples[[idx]]],
           ExtractTensorComponent[eomSep, field, chart,
             componentTuples[[idx]], additionalFields, computeChristoffels, metricMatrix, backgroundFieldRules]}
        ];
        Share[];
        Print["  [", Round[MemoryInUse[]/1024.^2], " MB] component ",
              idx, "/", Length[componentTuples]];,
        {idx, 1, Length[componentTuples]}
      ];
      Return[result]
    ]
  ]
];

(* === Batched TraceBasisDummy + Metric Evaluation (Adaptive) ===

   Fuses TraceBasisDummy, Expand, and early metric evaluation into a single
   batched loop to prevent intermediate expression blowup.

   TraceBasisDummy expands dummy basis index sums, producing O(dim^{2*n_dummy})
   intermediate terms per input term. For cross-field coupling (e.g., rank-2 h
   coupling to rank-1 a in Einstein-Maxwell), this creates ~970K terms from
   ~568 input terms. Early metric evaluation eliminates 99.97% by substituting
   η_{ij}=0 for off-diagonal i≠j. Fusing these steps keeps peak memory bounded
   by a single batch rather than the full expansion (~970K).

   Adaptive batch sizing (added 2026-03-12):
   The first batch runs at the default size (50 terms). After processing,
   the observed expansion factor (peak intermediate terms / input terms) is
   measured. If the expansion factor is high enough to exceed the target peak
   of ~2000 intermediate terms, subsequent batches are DECREASED in size.
   The batch size is NEVER increased above the default — Expand[] scales
   super-linearly with batch size due to cross-term interactions during
   simplification (empirically: 2x batch → 10-20x time).

     new_batch_size = Clamp(Floor(2000 / expansion_factor), 5, default)

   This adapts automatically to theory complexity:
   - Simple scalar theories (expansion ~1-5x): batch stays at 50 (default)
   - Mixed rank (expansion ~17x, e.g. flat EM): batch stays at 50 (within target)
   - High-rank cross-coupling (expansion ~40x+, e.g. curved EM): batch → 20-30
   - Extreme coupling (expansion ~100x): batch → 5-10 (very conservative)

   Per-batch diagnostics are printed showing TraceBD/Expand/Eval term counts,
   elapsed time, and memory usage. The initial and adapted batch sizes are
   reported in the summary line.

   Parameters:
   - componentEq_: Plus expression to decompose (additive terms of the EOM)
   - chart_: xCoba chart for basis evaluation
   - metricMatrix_: explicit metric matrix (or None for Minkowski fast path)
   - batchSize_: initial batch size (default 50, adapted after first batch)
   - backgroundFieldRules_: list of {head, bgHead, componentValues} for
     early evaluation of background field partial derivatives

   Ref: Gertsenshtein a-field OOM fix (7.6+ GB peak → <1 GB with batching). *)
BatchedTraceBasisDummyWithMetric[componentEq_, chart_, metricMatrix_, batchSize_:50,
  backgroundFieldRules_List:{}] := Module[
  {inputTerms, nTerms, result, batch, traced, currentBatchSize, expansionFactor,
   peakTerms, targetPeak, batchStart, batchEnd, tBatch},

  (* Single-term case: no batching needed *)
  If[Head[componentEq] =!= Plus,
    traced = TraceBasisDummy[componentEq];
    traced = Expand[traced];
    If[metricMatrix =!= None,
      traced = EvaluateMetricComponents[traced, chart, metricMatrix],
      traced = EvaluateMinkowskiMetric[traced, chart]
    ];
    Print["    step3-FusedTrace: ", Round[MemoryInUse[]/1024.^2], " MB, ",
          If[Head[Expand[traced]]===Plus, Length[Expand[traced]], 1], " terms (1 input, no batch)"];
    Return[Expand[traced]]
  ];

  inputTerms = List @@ componentEq;
  nTerms = Length[inputTerms];
  result = 0;

  (* Adaptive batch sizing: start with the requested batch size, then *)
  (* adjust based on observed expansion factor from the first batch.  *)
  (* Target: keep peak intermediate terms below ~2000 per batch.      *)
  (* IMPORTANT: Expand[] scales super-linearly with batch size due to *)
  (* cross-term interactions — never increase above the default.      *)
  (* Only decrease for high-expansion theories (e.g. curved EM).      *)
  (* Ref: Empirically, Einstein-Maxwell h decomposition generates     *)
  (* O(2000) terms per batch=50; batch=142 takes 20x longer due to   *)
  (* super-linear Expand cost. Keeping batch <= default is critical.  *)
  targetPeak = 2000;
  currentBatchSize = batchSize;
  expansionFactor = 0;
  batchStart = 1;

  While[batchStart <= nTerms,
    batchEnd = Min[batchStart + currentBatchSize - 1, nTerms];
    tBatch = AbsoluteTime[];
    batch = inputTerms[[batchStart ;; batchEnd]];
    (* TraceBasisDummy: sum dummy basis indices for this batch *)
    traced = Total[TraceBasisDummy /@ batch];
    tracedLen = If[Head[traced]===Plus, Length[traced], 1];
    (* Expand: propagate ComponentValue zeros (MetricInBasis + TT + background) *)
    traced = Expand[traced];
    expandLen = If[Head[traced]===Plus, Length[traced], 1];
    (* Early metric evaluation: collapse off-diagonal zeros immediately *)
    If[metricMatrix =!= None,
      traced = EvaluateMetricComponents[traced, chart, metricMatrix],
      traced = EvaluateMinkowskiMetric[traced, chart]
    ];
    traced = Expand[traced];
    evalLen = If[Head[traced]===Plus, Length[traced], 1];
    (* Early background field evaluation: collapse PD derivatives of background fields *)
    If[backgroundFieldRules =!= {},
      Do[traced = EvaluatePDBackgroundField[traced, chart, bg[[1]], bg[[2]],
           If[Length[bg] >= 3, bg[[3]], {}]],
        {bg, backgroundFieldRules}];
      traced = Expand[traced]
    ];
    peakTerms = Max[tracedLen, expandLen];
    (* Diagnostic: show batch statistics *)
    Print["    batch[", batchStart, ":", batchEnd, "/", nTerms, "]: ",
          "TraceBD=", tracedLen, " Expand=", expandLen, " Eval=", evalLen,
          " (", Round[AbsoluteTime[] - tBatch, 0.1], "s, ",
          Round[MemoryInUse[]/1024.^2], " MB)"];
    result += traced;
    (* Release batch memory *)
    batch =.; traced =.;

    (* After the first batch, adapt batch size based on observed expansion.  *)
    (* Only DECREASE — never increase above default. Expand[] is super-     *)
    (* linear: doubling batch size can 10-20x the cost due to cross-term    *)
    (* interactions during simplification. The default batchSize=50 is      *)
    (* empirically optimal for most theories.                               *)
    If[batchStart == 1 && peakTerms > 0,
      expansionFactor = N[peakTerms / (batchEnd - batchStart + 1)];
      If[expansionFactor > 0,
        currentBatchSize = Max[5, Floor[targetPeak / expansionFactor]];
        (* Never increase above default — Expand cost is super-linear *)
        currentBatchSize = Min[currentBatchSize, batchSize];
        If[currentBatchSize != batchSize,
          Print["    adaptive batch: expansion=",
                Round[expansionFactor, 0.1], "x, reducing batch to ", currentBatchSize,
                " (target peak=", targetPeak, ")"]
        ]
      ]
    ];

    batchStart = batchEnd + 1;
  ];

  result = Expand[result];
  Print["    step3-FusedTrace: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[result]===Plus, Length[result], 1], " terms",
        " (", nTerms, " input, batch=", currentBatchSize,
        If[currentBatchSize != batchSize,
          " (reduced from " <> ToString[batchSize] <> ")", ""],
        ")"];
  result
];

(* === Unified Tensor Component Extraction === *)
(* Single pipeline for any tensor rank. *)

ExtractTensorComponent[eom_, field_, chart_, componentIndices_List,
  additionalFields_List:{}, computeChristoffels_:False, metricMatrix_:None,
  backgroundFieldRules_List:{}] := Module[
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

  (* Step 1.5: SeparateFieldMetrics is now hoisted to DecomposeToComponents before the
     component loop. This step is a no-op here (already applied to eom). *)

  (* Step 2: ToBasis — term-by-term to avoid xperm segfault on large sums *)
  If[Head[componentEq] === Plus,
    componentEq = Total[ToBasis[chart] /@ List @@ componentEq],
    componentEq = ToBasis[chart][componentEq]
  ];
  Print["    step2-ToBasis: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  (* Steps 3+3.5+3.6 fused: TraceBasisDummy + Expand + early metric evaluation.
     Batched to prevent O(dim^{2*n_dummy}) intermediate memory blowup.
     Each batch: TraceBasisDummy → Expand → metric eval → accumulate.
     See BatchedTraceBasisDummyWithMetric for details. *)
  componentEq = BatchedTraceBasisDummyWithMetric[componentEq, chart, metricMatrix, 50, backgroundFieldRules];

  (* Step 4: Evaluate Christoffel symbols *)
  (* For flat spacetime (computeChristoffels=False): Γ=0, removes Christoffel terms.
     For curved spacetime: evaluates computed Christoffel components.
     After early metric eval (step 3.6), the expression is already much smaller,
     so this step operates on far fewer terms. *)
  componentEq = EvaluateChristoffelComponents[componentEq, chart, computeChristoffels, metricMatrix];
  componentEq = Expand[componentEq];
  Print["    step4-Christoffel: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  (* Step 5: Evaluate curvature tensors *)
  (* For constant metrics (flat or conformally flat): R=0, no-op.
     For non-constant metrics: substitutes computed Riemann/Ricci components. *)
  componentEq = EvaluateCurvatureComponents[componentEq, chart,
    If[metricMatrix =!= None, metricMatrix, None]];
  componentEq = Expand[componentEq];
  Print["    step5-Curvature: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  (* Step 6: Evaluate epsilon tensors *)
  componentEq = EvaluateEpsilonComponents[componentEq, chart, metricMatrix];
  componentEq = Expand[componentEq];
  Print["    step6-Epsilon: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  (* Step 7: Evaluate remaining metric components (idempotent after step 3.6).
     For most cases this is a no-op since metrics were already evaluated in
     step 3.6. However, Christoffel/curvature evaluation (steps 4-5) may
     introduce NEW metric components for curved spacetimes, so we re-evaluate
     to catch those. For flat spacetime this is guaranteed to be a no-op. *)
  If[metricMatrix =!= None,
    componentEq = EvaluateMetricComponents[componentEq, chart, metricMatrix],
    componentEq = EvaluateMinkowskiMetric[componentEq, chart]
  ];
  componentEq = Expand[componentEq];
  Print["    step7-Metric: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  If[metricMatrix =!= None,
    componentEq = EvaluatePDMetric[componentEq, chart, metricMatrix];
    componentEq = Expand[componentEq]
  ];

  (* Step 7b: Evaluate background field components and PD derivatives *)
  If[backgroundFieldRules =!= {},
    Do[
      componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentEq = Expand[componentEq];
    Print["    step7b-BackgroundField: ", Round[MemoryInUse[]/1024.^2], " MB, ",
          If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"]
  ];

  (* Step 8: Get coordinate symbols and replace ALL tensor fields with named scalar functions *)
  coordSyms = GetCoordinateSymbols[chart];
  dim = Length[coordSyms];

  allFieldHeads = Join[{fieldHead}, ExtractFieldHead /@ additionalFields];
  (* Pass metricMatrix so rank-2 fallback uses correct metric weights for curved metrics *)
  Do[
    componentEq = ReplaceTensorFieldComponents[componentEq, afh, chart, coordSyms, dim, metricMatrix],
    {afh, allFieldHeads}
  ];

  (* Step 9: Convert coordinate derivatives to Derivative form *)
  componentEq = ConvertCDToDerivatives[componentEq, chart];

  (* Step 9b: Catch residual background field Derivative[...][bgField][{mu,-chart}] forms
     that may survive ConvertCDToDerivatives *)
  If[backgroundFieldRules =!= {},
    Do[
      componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentEq = Expand[componentEq]
  ];

  (* Final simplification: Simplify then Expand to ensure cancellations
     (e.g., trace terms in linearized Einstein tensor) while preserving
     the Plus structure that EquationToJSONMultiField expects. *)
  Expand[Simplify[componentEq]]
];

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
(* metricMatrix: optional metric matrix for correct curved-metric index weighting in rank-2 fallback *)

ReplaceTensorFieldComponents[expr_, fh_, chart_, coordSyms_, dim_, metricMatrix_:None] := Module[
  {rank, result = expr},

  rank = Length[SlotsOfTensor[fh]];

  Which[
    rank === 0,
      result = result /. {
        fh[] :> Symbol[ToString[fh] <> "0"][Sequence @@ coordSyms]
      },

    rank === 1,
      With[{ch = chart, cs = coordSyms},
        (* After SeparateMetric, all field indices should be covariant (-ch) *)
        result = result /. {
          fh[{i_Integer, -ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
        };
        (* Safety fallback: if any contravariant field indices remain, warn and replace *)
        If[!FreeQ[result, fh[{_Integer, ch}]],
          Print["WARNING: Unexpected contravariant index for field ", fh,
                ". SeparateMetric may not have fully separated metric contractions. ",
                "Falling back to direct replacement (may produce incorrect signs)."];
          result = result /. {
            fh[{i_Integer, ch}] :> Symbol[ToString[fh] <> ToString[Abs[i]]][Sequence @@ cs]
          }
        ]
      ],

    rank === 2,
      result = ReplaceRank2FieldComponents[result, fh, chart, coordSyms, dim, metricMatrix],

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
(* metricMatrix: if provided (curved metric), non-covariant index fallback uses correct metric weights *)
ReplaceRank2FieldComponents[expr_, fh_, chart_, coordSyms_, dim_, metricMatrix_:None] := Module[
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

        (* Primary: covariant-covariant (canonical after SeparateMetric) *)
        result = result /. {
          fh[{pair[[1]], -ch}, {pair[[2]], -ch}] :> sym[Sequence @@ cs]
        };

        (* For symmetric: also map swapped covariant indices to same component *)
        If[symQ && pair[[1]] =!= pair[[2]],
          result = result /. {
            fh[{pair[[2]], -ch}, {pair[[1]], -ch}] :> sym[Sequence @@ cs]
          }
        ];

        (* Safety fallback: if any non-covariant index configs remain for this pair *)
        If[!FreeQ[result, fh[{pair[[1]], _}, {pair[[2]], _}]] ||
           (symQ && pair[[1]] =!= pair[[2]] && !FreeQ[result, fh[{pair[[2]], _}, {pair[[1]], _}]]),
          Print["WARNING: Non-covariant rank-2 indices for field ", fh,
                " at pair ", pair, ". Falling back to metric-weighted replacement."];
          (* For curved metrics (diagonal): apply correct g^{ia}g^{jb} weights.
             h^{ij} = (g^{ii})(g^{jj}) h_{ij}  (diagonal metric, no sum)
             h^i_j  = g^{ii} h_{ij}
             h_i^j  = g^{jj} h_{ij}
             For flat metric (metricMatrix===None): weights are 1 (historical behavior). *)
          Module[{wUU, wUDi, wDUi, wUUji, wUDji, wDUji, invMatrix, diagInv},
            If[metricMatrix =!= None,
              invMatrix = GetCachedInverseMetric[metricMatrix];
              (* Check diagonal: use diagonal weights; warn if non-diagonal *)
              If[DiagonalMatrixQ[invMatrix],
                diagInv = Diagonal[invMatrix];
                wUU   = diagInv[[pair[[1]]+1]] * diagInv[[pair[[2]]+1]];
                wUDi  = diagInv[[pair[[1]]+1]];  (* h^i_j: raise first index *)
                wDUi  = diagInv[[pair[[2]]+1]];  (* h_i^j: raise second index *)
                wUUji = diagInv[[pair[[2]]+1]] * diagInv[[pair[[1]]+1]];
                wUDji = diagInv[[pair[[2]]+1]];
                wDUji = diagInv[[pair[[1]]+1]],
                (* Non-diagonal: correct weighting requires sum over all components — not implemented.
                   Fall back to unweighted (wrong) behavior with a clear error message. *)
                Print["ERROR: Non-diagonal metric in ReplaceRank2FieldComponents fallback. ",
                      "Metric-weighted index raising not implemented for non-diagonal metrics. ",
                      "Equations may be incorrect for field ", fh, " at pair ", pair, "."];
                wUU = wUDi = wDUi = wUUji = wUDji = wDUji = 1
              ],
              (* Flat metric or None: identity weights (historical behavior) *)
              wUU = wUDi = wDUi = wUUji = wUDji = wDUji = 1
            ];
            result = result /. {
              fh[{pair[[1]], ch}, {pair[[2]], ch}]  :> wUU  * sym[Sequence @@ cs],
              fh[{pair[[1]], ch}, {pair[[2]], -ch}] :> wUDi * sym[Sequence @@ cs],
              fh[{pair[[1]], -ch}, {pair[[2]], ch}] :> wDUi * sym[Sequence @@ cs]
            };
            If[symQ && pair[[1]] =!= pair[[2]],
              result = result /. {
                fh[{pair[[2]], ch}, {pair[[1]], ch}]  :> wUUji * sym[Sequence @@ cs],
                fh[{pair[[2]], ch}, {pair[[1]], -ch}] :> wUDji * sym[Sequence @@ cs],
                fh[{pair[[2]], -ch}, {pair[[1]], ch}] :> wDUji * sym[Sequence @@ cs]
              }
            ]
          ]
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
      Module[{tuple = tuples[[k]], seqIdx = k - 1, sym, indexConfigs, pattern},
        sym = Symbol[ToString[fh] <> ToString[seqIdx]];

        (* Primary: all-covariant (canonical after SeparateMetric) *)
        pattern = Table[{tuple[[n]], -ch}, {n, rank}];
        result = result /. {fh @@ pattern :> sym[Sequence @@ cs]};

        (* Safety fallback: if field still appears with this tuple, try all configs *)
        If[!FreeQ[result, fh],
          indexConfigs = Tuples[{ch, -ch}, rank];
          Do[
            Module[{pat},
              pat = Table[{tuple[[n]], config[[n]]}, {n, rank}];
              If[pat =!= pattern,  (* Skip already-applied covariant config *)
                result = result /. {fh @@ pat :> sym[Sequence @@ cs]}
              ]
            ],
            {config, indexConfigs}
          ]
        ]
      ],
      {k, 1, Length[tuples]}
    ]
  ];

  result
];


(* === Scalar Expression Decomposition (Phase K: Canonical Momentum) === *)
(*
  DecomposeScalarExpression: converts a scalar tensor expression (like the
  Hamiltonian density H) to component form. Uses the same pipeline as
  the scalar case in DecomposeToComponents:
    1. ToBasis → TraceBasisDummy → EvaluateChristoffelComponents
    2. EvaluateEpsilonComponents → Evaluate metric
    3. ReplaceTensorFieldComponents → ConvertCDToDerivatives

  Unlike DecomposeToComponents (which takes a tensor EOM + field to determine
  which components to extract), this function takes a scalar expression and
  returns a single component-form expression.

  allFieldHeads: ALL tensor field heads appearing in the expression.
  This includes primary fields, additional (cross-coupled) fields,
  and background fields.  All must be passed so that
  ReplaceTensorFieldComponents converts them to coordinate functions.
*)

DecomposeScalarExpression[expr_, chart_, allFieldHeads_List] :=
  DecomposeScalarExpression[expr, chart, allFieldHeads,
    "ComputeChristoffels" -> Automatic, "MetricMatrix" -> None, "BackgroundFieldRules" -> {}];

(* DecomposeScalarExpression: memory-efficient batched pipeline.                *)
(* Mirrors the step-by-step approach in DecomposeToComponents (steps 2-9):     *)
(*   ToBasis (term-by-term) → BatchedTraceBasisDummy → Christoffel → Curvature *)
(*   → Epsilon → Metric → PDMetric → BackgroundField → FieldReplace → CDConvert *)
(* The batched approach prevents O(dim^{2K}) memory blowup that occurs when    *)
(* TraceBasisDummy + Expand run on the full expression at once.                *)
(* Ref: OOM crash on radial Gertsenshtein canonical pipeline (term 2/34,       *)
(* Einstein-Hilbert sector in spherical coords, exceeded 10 GB unbatched).     *)
DecomposeScalarExpression[expr_, chart_, allFieldHeads_List, opts:OptionsPattern[DecomposeToComponents]] := Module[
  {componentExpr, metricMatrix, computeChristoffelsOption, shouldComputeChristoffels,
   coords, dim, coordSyms, backgroundFieldRules, computeChristoffels},

  metricMatrix = OptionValue[DecomposeToComponents, {opts}, "MetricMatrix"];
  computeChristoffelsOption = OptionValue[DecomposeToComponents, {opts}, "ComputeChristoffels"];
  backgroundFieldRules = OptionValue[DecomposeToComponents, {opts}, "BackgroundFieldRules"];

  coords = GetCoordinateSymbols[chart];
  dim = GetChartDimension[chart];
  coordSyms = coords;

  (* Auto-detect whether Christoffel computation is needed *)
  shouldComputeChristoffels = Which[
    computeChristoffelsOption === True, True,
    computeChristoffelsOption === False, False,
    computeChristoffelsOption === Automatic,
      If[metricMatrix === None, False, IsNonConstantMetric[metricMatrix, coords]],
    True, False
  ];
  computeChristoffels = shouldComputeChristoffels;

  componentExpr = expr;

  (* Expand Scalar[] wrappers.  Gauge-fixing terms produce                   *)
  (* Scalar[CD_a V^a]^2 — the Scalar wrapper is opaque and prevents ToBasis/ *)
  (* TraceBasisDummy from decomposing the inner expression to chart-basis     *)
  (* components.  Expand by applying ToBasis + TraceBasisDummy to the inner   *)
  (* expression, replacing the Scalar with the resulting component sum.       *)
  (* This must happen before the main ToBasis call.                           *)
  componentExpr = componentExpr //. Scalar[x_] :> Module[{inner},
    inner = ToBasis[chart][x];
    inner = TraceBasisDummy[inner];
    inner
  ];

  (* For curved spacetime: expand Christoffel symbols to metric derivatives *)
  If[shouldComputeChristoffels && metricMatrix =!= None,
    Module[{covdOps, covdOp},
      covdOps = Cases[componentExpr, (f_)[_][_] /; CovDQ[f] :> f, {0, Infinity}] // DeleteDuplicates;
      If[Length[covdOps] > 0,
        covdOp = covdOps[[1]];
        componentExpr = ExpandChristoffelsToMetricDerivatives[componentExpr, covdOp, chart]
      ]
    ]
  ];

  (* Separate metric contractions before ToBasis *)
  (* Ensures all field tensor indices are in canonical (covariant) form *)
  componentExpr = SeparateFieldMetrics[componentExpr, chart];

  (* Step 2: ToBasis — term-by-term to avoid xperm segfault on large sums *)
  If[Head[componentExpr] === Plus,
    componentExpr = Total[ToBasis[chart] /@ List @@ componentExpr],
    componentExpr = ToBasis[chart][componentExpr]
  ];
  Print["    [scalar] step2-ToBasis: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];

  (* NOTE: Pre-TraceBasisDummy plane-wave zeroing was attempted here via three    *)
  (* approaches; ALL caused performance regression.  Root cause verified:        *)
  (* after ToBasis, PD operators retain abstract dummy indices (e.g. -h$34637,   *)
  (* not {2,-chart}) — concrete basis integers only appear INSIDE TraceBasisDummy *)
  (* enumeration.  Therefore:                                                    *)
  (* 1. ComponentValue PD zeroing: +9-56% slower (xAct internal rule tables add  *)
  (*    pattern-matching overhead but DON'T auto-evaluate during Mathematica eval)*)
  (* 2. Direct DownValues (specific integers): +66% slower (adds O(1) hash       *)
  (*    lookups to EVERY TraceBasisDummy evaluation step — metric DownValues work *)
  (*    because metrics are multiplicative factors; PD wraps inner expressions    *)
  (*    that are already evaluated)                                               *)
  (* 3. ReplaceAll on PD[{killedAxis,-chart}]: impossible (indices are abstract)  *)
  (* The per-term plane-wave reduction in _derive.py (after ConvertCDToDerivatives*)
  (* resolves indices) correctly zeros transverse Derivative patterns.            *)

  (* Steps 3+3.5+3.6 fused: TraceBasisDummy + Expand + early metric evaluation.
     Batched to prevent O(dim^{2*n_dummy}) intermediate memory blowup.
     Each batch: TraceBasisDummy → Expand → metric eval → accumulate.
     See BatchedTraceBasisDummyWithMetric for details. *)
  componentExpr = BatchedTraceBasisDummyWithMetric[componentExpr, chart, metricMatrix, 50, backgroundFieldRules];

  (* Step 4: Evaluate Christoffel symbols *)
  If[!computeChristoffels,
    componentExpr = EvaluateChristoffelComponents[componentExpr, chart, False]
  ];
  componentExpr = Expand[componentExpr];
  Print["    [scalar] step4-Christoffel: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];

  (* Step 5: Evaluate curvature tensors *)
  componentExpr = EvaluateCurvatureComponents[componentExpr, chart,
    If[metricMatrix =!= None, metricMatrix, None]];
  componentExpr = Expand[componentExpr];
  Print["    [scalar] step5-Curvature: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];

  (* Step 6: Evaluate epsilon tensor components *)
  componentExpr = EvaluateEpsilonComponents[componentExpr, chart, metricMatrix];
  componentExpr = Expand[componentExpr];
  Print["    [scalar] step6-Epsilon: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];

  (* Step 7: Evaluate remaining metric components *)
  If[metricMatrix =!= None,
    componentExpr = EvaluateMetricComponents[componentExpr, chart, metricMatrix],
    componentExpr = EvaluateMinkowskiMetric[componentExpr, chart]
  ];
  componentExpr = Expand[componentExpr];
  Print["    [scalar] step7-Metric: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];

  (* Evaluate partial derivatives of metric for curved spacetime *)
  If[metricMatrix =!= None,
    componentExpr = EvaluatePDMetric[componentExpr, chart, metricMatrix];
    componentExpr = Expand[componentExpr]
  ];

  (* Step 7b: Evaluate background field components and PD derivatives *)
  If[backgroundFieldRules =!= {},
    Do[
      componentExpr = EvaluatePDBackgroundField[componentExpr, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentExpr = Expand[componentExpr];
    Print["    [scalar] step7b-BackgroundField: ", Round[MemoryInUse[]/1024.^2], " MB, ",
          If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"]
  ];

  (* Step 8: Replace ALL tensor fields with named scalar component functions *)
  Do[
    componentExpr = ReplaceTensorFieldComponents[componentExpr, afh, chart, coordSyms, dim, metricMatrix],
    {afh, allFieldHeads}
  ];

  (* Step 9: Convert coordinate derivatives to Derivative form *)
  componentExpr = ConvertCDToDerivatives[componentExpr, chart];

  (* Step 9b: Catch residual background field Derivative forms *)
  If[backgroundFieldRules =!= {},
    Do[
      componentExpr = EvaluatePDBackgroundField[componentExpr, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentExpr = Expand[componentExpr]
  ];

  (* Reclaim memory before returning *)
  Share[];
  Expand[componentExpr]
];


End[];
EndPackage[];
