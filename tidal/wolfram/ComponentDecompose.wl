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

ComponentEulerLagrange::usage =
  "ComponentEulerLagrange[lagrangian, fieldFuncs, coords] computes Euler-Lagrange \
equations from a COMPONENT (scalar) Lagrangian by standard variational calculus. \
Each field function (e.g., tidalH0, tidalT5) is varied independently via D[] with \
integration by parts for higher-order derivatives. Returns a list of EOM expressions \
(one per field function, in the same order as fieldFuncs). \
For quadratic Lagrangians, produces linear equations.";

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


(* === Scalar Wrapper Expansion === *)
(* xPert wraps contracted tensor subexpressions in Scalar[...], marking  *)
(* them as scalar invariants.  ToBasis[chart] treats Scalar as opaque —  *)
(* it will not enter the wrapper to convert abstract indices to basis     *)
(* form.  This leaves abstract dummy indices (e.g. -g$27013) inside,     *)
(* causing ConvertCDToDerivatives to fail (it expects {i, -chart} form). *)
(*                                                                        *)
(* ExpandScalarWrappers resolves this by applying ToBasis + TraceBasis-   *)
(* Dummy to the contents of each Scalar[], converting abstract indices    *)
(* to concrete basis sums.  The Scalar wrapper is then removed (the      *)
(* result is a plain numeric/symbolic expression).                        *)
(*                                                                        *)
(* Uses bounded iteration (max 5 passes) with fixed-point check to       *)
(* prevent $RecursionLimit if ToBasis re-introduces Scalar wrappers.     *)
(* Zero overhead for Scalar-free expressions (FreeQ early exit).         *)
(*                                                                        *)
(* Used by: ExtractTensorComponent, DecomposeToComponents (rank 0),      *)
(*          DecomposeScalarExpression.                                    *)
ExpandScalarWrappers[expr_, chart_, computeChristoffels_:False] := Module[
  {result = expr, prev, iter = 0, maxIter = 5,
   metric, covd, christoffelPD, hasCDScalars},
  (* Early exit if no Scalar wrappers present *)
  If[FreeQ[result, Scalar], Return[result]];

  (* Detect if any Scalar contents have CDs — gate full staggered       *)
  (* pipeline vs fast path.  R̃²-decomposed torsion theories produce   *)
  (* Scalar[η^{ab}(R_{ab} + ∇∇T...)] with nested CDs that require     *)
  (* Christoffel zeroing + SeparateMetric between ToBasis passes.       *)
  hasCDScalars = !FreeQ[result, Scalar[x_ /; !FreeQ[x, _?CovDQ]]];
  If[hasCDScalars,
    Module[{manifold2, tangentBundle2, metrics2},
      manifold2 = ManifoldOfChart[chart];
      tangentBundle2 = Symbol["Tangent" <> ToString[manifold2]];
      metrics2 = MetricsOfVBundle[tangentBundle2];
      If[ListQ[metrics2] && Length[metrics2] > 0,
        metric = First[metrics2];
        covd = CovDOfMetric[metric];
        christoffelPD = Symbol["Christoffel" <> ToString[covd] <> "PD" <> ToString[chart]],
        hasCDScalars = False  (* No metric → fall back to simple path *)
      ]
    ]
  ];

  (* Suppress Validate::repeated during Scalar expansion.                *)
  (* When Scalar[x] contents share abstract index names with the outer   *)
  (* expression (common in R̃²-decomposed torsion theories), ToBasis on *)
  (* the inner expression triggers Validate::repeated which Throws.     *)
  (* This is a false positive — the indices ARE properly contracted     *)
  (* within each Scalar; the collision is temporary and resolves after  *)
  (* TraceBasisDummy sums over basis indices.                           *)
  While[!FreeQ[result, Scalar] && iter < maxIter,
    prev = result;
    result = result /. Scalar[x_] :> Module[{inner = x},
      (* Apply CD shorthand rules BEFORE ToBasis — while indices are     *)
      (* still in abstract form so MakeRule patterns match. Without this *)
      (* CD[-a]@field[...] survives ToBasis (opaque wrapper) and the     *)
      (* post-ToBasis form has DummyIn indices that MakeRule can't match.*)
      If[ListQ[Global`$CDShorthandRules] && Length[Global`$CDShorthandRules] > 0,
        Do[inner = inner /. rule, {rule, Global`$CDShorthandRules}]
      ];

      If[hasCDScalars && !FreeQ[inner, _?CovDQ],
        (* Full staggered pipeline for Scalar contents with nested CDs.  *)
        (* Ref: supervisor's EuclideanSplinter (commit 4a89164).         *)
        (*                                                               *)
        (* Strategy (3-phase):                                           *)
        (* Phase 1: Record ALL abstract dummy pairs via IndicesOf while  *)
        (*   the expression still has proper xAct tensor structure.     *)
        (*   IndicesOf is theory-agnostic: handles CovD, metrics,       *)
        (*   Christoffels, Basis tensors, delta — ALL index positions.  *)
        (* Phase 2: Staggered ToBasis pipeline resolves most pairs.     *)
        (* Phase 3: Sum-trace any surviving pairs from the Phase 1 set. *)
        (*                                                               *)
        (* Why previous approaches failed:                              *)
        (* - ToBasis: only converts FREE indices, not contracted dummies *)
        (* - TraceBasisDummy: only traces basis-form dummies             *)
        (* - ComponentArray: only expands free indices (scalars → no-op) *)
        (* - Cases-based detection: misses Basis/delta/nested CovD pos.  *)
        (* - RenameDummies: returns identical for canonical expressions  *)
        (* IndicesOf before pipeline + Sum after = complete resolution.  *)

        (* Collect ALL abstract indices in the Scalar content BEFORE the   *)
        (* staggered pipeline, while the expression has proper tensor   *)
        (* structure. Use Cases on ALL subexpressions to find abstract  *)
        (* symbols, then identify dummy pairs by matching up/down.      *)
        (* IndicesOf was tried but returns expression-level results,    *)
        (* not individual index symbols (commit history: eb1c216 etc.). *)
        Module[{allDummyPairs = {}, dimScalar = DimOfManifold[ManifoldOfChart[chart]]},
          Module[{allAbstract, upSet, downSet},
            (* Find ALL abstract symbols anywhere in the expression *)
            allAbstract = Cases[inner, s_Symbol?AbstractIndexQ, {0, Infinity}] // DeleteDuplicates;
            (* An abstract symbol is a dummy if it appears BOTH as +s and -s.
               In xAct: up = s (bare symbol), down = -s (Times[-1, s]).
               Also check DownIndexQ on the arguments of tensor/CovD heads. *)
            upSet = allAbstract;
            (* Build down set from multiple sources *)
            downSet = {};
            (* Source 1: negated symbols in the expression tree *)
            downSet = Union[downSet,
              Cases[inner, -s_Symbol?AbstractIndexQ :> s, {0, Infinity}]];
            (* Source 2: CovD down-index positions *)
            downSet = Union[downSet,
              Cases[inner, (f_?CovDQ)[idx_][_] /; DownIndexQ[idx] :> UpIndex[idx],
                {0, Infinity}]];
            (* Source 3: general DownIndexQ in any position *)
            downSet = Union[downSet,
              Cases[inner, idx_ /; DownIndexQ[idx] && AbstractIndexQ[UpIndex[idx]] :> UpIndex[idx],
                {0, Infinity}]];
            allDummyPairs = Intersection[upSet, downSet];
          ];

          (* Phase 2: Staggered ToBasis pipeline (resolves most pairs).  *)
          inner = Catch[ToCanonical[inner]];
          inner = ContractMetric[inner, metric];
          Module[{prev2 = -1, cur2, iter2 = 0},
            While[iter2 < 8,
              If[Head[inner] === Plus,
                inner = Total[ToBasis[chart] /@ List @@ inner],
                inner = ToBasis[chart][inner]
              ];
              If[!computeChristoffels,
                inner = inner /. christoffelPD -> Zero,
                inner = ToValues[inner]
              ];
              cur2 = If[Head[inner] === Plus, Length[inner], 1];
              If[cur2 == prev2, Break[]];
              prev2 = cur2;
              iter2++;
            ];
            inner = SeparateMetric[metric][inner];
            If[Head[inner] === Plus,
              inner = Total[ToBasis[chart] /@ List @@ inner],
              inner = ToBasis[chart][inner]
            ];
            If[!computeChristoffels,
              inner = inner /. christoffelPD -> Zero,
              inner = ToValues[inner]
            ];
            inner = TraceBasisDummy[inner];
            inner = ToValues[inner];
            inner = ToValues[inner];
          ];

          (* Phase 3: Sum-trace surviving abstract pairs.                 *)
          (* Filter to pairs that still appear in the expression.       *)
          Module[{remaining, replaced},
            remaining = Select[allDummyPairs, !FreeQ[inner, #] &];
            If[Length[remaining] > 0,
              replaced = inner;
              Do[
                replaced = Sum[
                  replaced /. {p -> {i, chart}, ChangeIndex[p] -> {i, -chart}},
                  {i, 0, dimScalar - 1}
                ];
                replaced = Expand[replaced];
                If[!computeChristoffels, replaced = replaced /. christoffelPD -> Zero];
                replaced = ToValues[replaced];
                replaced = Expand[replaced],
                {p, remaining}
              ];
              inner = replaced;
            ]
          ];

          (* Fallback: Cases-based detection if IndicesOf missed any.    *)
          (* This handles edge cases where IndicesOf returns empty or    *)
          (* the expression gained new abstract indices during pipeline. *)
          Module[{upIdx, downIdx, pairedIdx, replaced, iterTrace = 0, prevTrace},
            While[iterTrace < 5,
              prevTrace = inner;
              upIdx = Cases[inner, s_Symbol /; AbstractIndexQ[s], {0, Infinity}] // DeleteDuplicates;
              downIdx = Cases[inner, (f_?CovDQ)[idx_][_] /; DownIndexQ[idx] :> UpIndex[idx], {0, Infinity}] // DeleteDuplicates;
              downIdx = Union[downIdx,
                Cases[inner, (f_?xTensorQ)[___, -s_Symbol, ___] /; AbstractIndexQ[s] :> s, {0, Infinity}]];
              downIdx = Union[downIdx,
                Cases[inner, (f_?MetricQ)[___, -s_Symbol, ___] /; AbstractIndexQ[s] :> s, {0, Infinity}]];
              (* Also check Basis tensor and general DownIndexQ positions *)
              downIdx = Union[downIdx,
                Cases[inner, Basis[-s_Symbol, _] /; AbstractIndexQ[s] :> s, {0, Infinity}]];
              downIdx = Union[downIdx,
                Cases[inner, (f_?xTensorQ)[___, idx_, ___] /; DownIndexQ[idx] && AbstractIndexQ[UpIndex[idx]] :> UpIndex[idx], {0, Infinity}]];
              pairedIdx = Intersection[upIdx, downIdx];
              If[Length[pairedIdx] == 0, Break[]];
              replaced = inner;
              Do[
                replaced = Sum[
                  replaced /. {p -> {i, chart}, ChangeIndex[p] -> {i, -chart}},
                  {i, 0, dimScalar - 1}
                ];
                replaced = Expand[replaced];
                If[!computeChristoffels, replaced = replaced /. christoffelPD -> Zero];
                replaced = ToValues[replaced];
                replaced = Expand[replaced],
                {p, pairedIdx}
              ];
              inner = replaced;
              iterTrace++;
              If[inner === prevTrace, Break[]]
            ]
          ];
        ],
        (* Fast path: simple Scalar contents without CDs *)
        inner = ToBasis[chart][inner];
        inner = TraceBasisDummy[inner];
      ];
      inner
    ];
    (* Normalize negative-integer basis indices to non-negative.          *)
    (* xCoba uses {-n, -chart} for covariant component n (the negative   *)
    (* on the integer mirrors the negative on the chart).  For component  *)
    (* extraction and CD→Derivative conversion, {-n, -chart} ≡ {n,-chart}*)
    (* Normalize so downstream patterns (GenerateCDRules, etc.) match.   *)
    (* Ref: CommonUtilities.wl line 511 — same convention, handled w/ Abs*)
    result = result /. {n_Integer?Negative, s:(chart | -chart)} :> {-n, s};
    iter++;
    If[result === prev, Break[]]  (* Fixed point reached *)
  ];
  (* Note: Validate::repeated is managed by the caller (DecomposeToComponents / *)
  (* StaggeredToBasis), not here — do NOT re-enable it.                         *)
  If[!FreeQ[result, Scalar] && iter >= maxIter,
    Print["WARNING: ExpandScalarWrappers did not fully converge after ",
      maxIter, " iterations. Remaining Scalar[]: ",
      Short[Cases[result, Scalar[x_] :> Short[x, 2], {0, Infinity}], 3]]
  ];
  result
];


(* === Staggered ToBasis Pipeline === *)
(* Ref: supervisor's EuclideanSplinter (commit 4a89164).                   *)
(* Iteratively applies ToBasis + Christoffel zeroing/evaluation to peel    *)
(* covariant derivative layers one at a time, preventing the O(dim^{2K})   *)
(* dummy-index explosion that occurs with a single ToBasis + TraceBasis-   *)
(* Dummy pass.  After staggered ToBasis, ToValues substitutes pre-computed *)
(* ComponentValues (from MetricCompute) and TraceBasisDummy traces only     *)
(* the remaining unresolved dummies — a much smaller set.                  *)
(*                                                                          *)
(* Pipeline: ToBasis → Christoffel{Zero/Eval} → ToBasis → Christoffel →   *)
(*           SeparateMetric → ToBasis → ToValues → TraceBasisDummy →       *)
(*           ToValues                                                       *)
(*                                                                          *)
(* Used by: ExtractTensorComponent, DecomposeScalarExpression.             *)

StaggeredToBasis[expr_, chart_, computeChristoffels_:False] := Module[
  {e = expr, prevLen = -1, curLen, iter = 0, maxIter = 4, metric, covd, christoffelPD},
  (* Suppress Validate::repeated during the staggered pipeline.          *)
  (* R̃-decomposed expressions can have abstract index name collisions   *)
  (* between Scalar wrapper contents and the outer expression after      *)
  (* ExpandScalarWrappers resolves them. These are false positives —    *)
  (* the indices are properly contracted, just temporarily sharing names.*)
  Off[Validate::repeated];

  (* Get the metric associated with this chart.
     Strategy: ManifoldOfChart → Tangent{manifold} → MetricsOfVBundle *)
  Module[{manifold, tangentBundle, metrics},
    manifold = ManifoldOfChart[chart];
    tangentBundle = Symbol["Tangent" <> ToString[manifold]];
    metrics = MetricsOfVBundle[tangentBundle];
    If[!ListQ[metrics] || Length[metrics] == 0,
      (* No metric found — fall back to basic ToBasis + TraceBasisDummy *)
      Print["WARNING: StaggeredToBasis: no metric for chart ", chart];
      If[Head[e] === Plus,
        e = Total[ToBasis[chart] /@ List @@ e],
        e = ToBasis[chart][e]
      ];
      e = TraceBasisDummy[e];
      Return[e]
    ];
    metric = First[metrics];
    covd = CovDOfMetric[metric];
  ];
  (* Determine the Christoffel symbol name: Christoffel{CD}PD{chart} *);
  christoffelPD = Symbol["Christoffel" <> ToString[covd] <> "PD" <> ToString[chart]];

  (* Staggered ToBasis: each pass converts the outermost abstract indices *)
  (* to basis form, then Christoffel zeroing/evaluation collapses terms   *)
  (* before the next pass introduces more basis indices from inner CDs.  *)
  While[iter < maxIter,
    (* ToBasis — term-by-term to avoid xperm segfault on large sums *)
    If[Head[e] === Plus,
      e = Total[ToBasis[chart] /@ List @@ e],
      e = ToBasis[chart][e]
    ];
    (* Zero or evaluate Christoffels introduced by this ToBasis pass *)
    If[!computeChristoffels,
      e = e /. christoffelPD -> Zero,
      (* For curved spacetime: ToValues substitutes MetricCompute-cached Christoffels *)
      e = ToValues[e]
    ];
    curLen = If[Head[e] === Plus, Length[e], 1];
    If[curLen == prevLen, Break[]];  (* Fixed point — no new indices to resolve *)
    prevLen = curLen;
    iter++;
  ];

  (* SeparateMetric: resolve any metric contractions introduced by ToBasis *)
  e = SeparateMetric[metric][e];
  (* Final ToBasis pass for metric-separated terms *)
  If[Head[e] === Plus,
    e = Total[ToBasis[chart] /@ List @@ e],
    e = ToBasis[chart][e]
  ];
  (* Zero Christoffels again after SeparateMetric + ToBasis *)
  If[!computeChristoffels,
    e = e /. christoffelPD -> Zero,
    e = ToValues[e]
  ];

  (* NOTE: Supervisor's EuclideanSplinter calls ToValues BEFORE             *)
  (* TraceBasisDummy, but that only helps for rank>0 expressions where      *)
  (* ComponentArray pre-resolves free indices. For scalar expressions       *)
  (* (canonical Lagrangian), ALL indices are contracted — ToValues cannot   *)
  (* resolve contracted pairs (only TraceBasisDummy can enumerate them).    *)
  (* Tested: Expand+ToValues×2 before TraceBasisDummy → no improvement     *)
  (* (42.7s unchanged). TraceBasisDummy's O(dim^{2K}) cost is fundamental  *)
  (* for K~4 contracted pairs per EH Ricci scalar product term.            *)

  (* TraceBasisDummy: per-term to prevent O(dim^{2K}) explosion.         *)
  (* The full expression may have K=45 contracted dummy pairs across all *)
  (* additive terms, but each individual term has only ~3 pairs.         *)
  (* Per-term: O(N × dim^6) instead of O(dim^90). Same pattern as       *)
  (* BatchedTraceBasisDummyWithMetric (line 1195).                       *)
  If[Head[e] === Plus,
    e = Total[TraceBasisDummy /@ List @@ e],
    e = TraceBasisDummy[e]
  ];
  e = Expand[e];

  (* ToValues: substitute remaining pre-computed ComponentValues *)
  (* (fields, background fields, any not resolved by first pass) *)
  e = ToValues[e];
  e = ToValues[e];  (* Double ToValues: some expressions need two passes *)

  (* Force-resolve remaining abstract-index operators.                    *)
  (* After ToValues + TraceBasisDummy, some CD operators or CD shorthand *)
  (* tensors still have abstract DummyIn indices because:                *)
  (*   - CD wrappers are opaque to ToBasis                              *)
  (*   - CD shorthands (CD1*, CD2*, ...) are excluded from isCDlikeQ    *)
  (* Apply ToBasis + TraceBasisDummy + ToValues iteratively to force     *)
  (* basis-index conversion for BOTH raw CDs AND CD shorthands.         *)
  (* Ref: supervisor's EuclideanSplinter applies ToBasis exhaustively.   *)
  If[!FreeQ[e, _?isCDlikeQ] ||
     !FreeQ[e, (f_?xTensorQ)[__] /; StringMatchQ[ToString[f], "CD" ~~ DigitCharacter ~~ __]],
    Module[{prev2 = -1, cur2, iter2 = 0},
      While[iter2 < 8,
        If[Head[e] === Plus,
          e = Total[ToBasis[chart] /@ List @@ e],
          e = ToBasis[chart][e]
        ];
        If[!computeChristoffels,
          e = e /. christoffelPD -> Zero,
          e = ToValues[e]
        ];
        cur2 = If[Head[e] === Plus, Length[e], 1];
        If[cur2 == prev2, Break[]];
        prev2 = cur2;
        iter2++;
      ];
      e = SeparateMetric[metric][e];
      If[Head[e] === Plus,
        e = Total[ToBasis[chart] /@ List @@ e],
        e = ToBasis[chart][e]
      ];
      If[!computeChristoffels,
        e = e /. christoffelPD -> Zero,
        e = ToValues[e]
      ];
      If[Head[e] === Plus,
        e = Total[TraceBasisDummy /@ List @@ e],
        e = TraceBasisDummy[e]
      ];
      e = Expand[e];
      e = ToValues[e];
      e = ToValues[e];
    ];
  ];

  (* Trace remaining abstract DummyIn indices via IndicesOf.              *)
  (* After all ToBasis passes, some CovD operators (raw or shorthand)    *)
  (* still have abstract DummyIn indices that ToBasis cannot convert     *)
  (* because they're contracted in products. Use IndicesOf to find ALL   *)
  (* abstract dummy pairs, then sum each pair over basis values.         *)
  (* Gated on: any remaining abstract indices in the expression.        *)
  (* Trace abstract DummyIn indices via Cases-based search.              *)
  (* IndicesOf[] returns empty for these expressions (DummyIn indices     *)
  (* from ExpandScalarWrappers are not tracked by xAct's index system). *)
  (* Instead: find DummyIn-allocated symbols (g$NNN) via Cases,         *)
  (* matching them in CovD indices and metric indices to find pairs.    *)
  Module[{dimLocal, cdDownSyms, metricSyms, allUpSyms, paired, replaced},
    dimLocal = DimOfManifold[ManifoldOfChart[chart]];
    (* Find abstract symbols in CovD down-index position *)
    (* Note: xAct may store -g$ as Times[-1, g$] or as a direct minus *)
    cdDownSyms = Cases[e, (f_?CovDQ)[idx_][_] /; DownIndexQ[idx] :> UpIndex[idx], {0, Infinity}] // DeleteDuplicates;
    cdDownSyms = Select[cdDownSyms, AbstractIndexQ[#] && StringMatchQ[ToString[#], "g$" ~~ __] &];
    (* Find abstract symbols in metric positions *)
    metricSyms = Cases[e, (f_?MetricQ)[args__] :>
      Sequence @@ Cases[{args}, s_Symbol /; AbstractIndexQ[s]], {0, Infinity}] // DeleteDuplicates;
    (* Find DummyIn-generated abstract symbols (g$NNN format).             *)
    (* These are the indices from ExpandScalarWrappers that xAct's        *)
    (* ToBasis/TraceBasisDummy cannot resolve. Standard abstract indices  *)
    (* (a, b, c, ...) are left for the CD→Derivative fallback.           *)
    allUpSyms = Cases[e, s_Symbol /; AbstractIndexQ[s] &&
      StringMatchQ[ToString[s], "g$" ~~ __], {0, Infinity}] // DeleteDuplicates;
    (* Pairs: symbols that appear as BOTH up (in metric/tensor) and down (in CD) *)
    paired = Intersection[allUpSyms, cdDownSyms];
    (* Also check for non-CD pairs: DummyIn indices in metrics paired with tensor down-indices *)
    Module[{tensorDownSyms},
      tensorDownSyms = Cases[e, (f_?xTensorQ)[-s_Symbol, ___] /; StringMatchQ[ToString[s], "g$" ~~ __] :> s,
        {0, Infinity}] // DeleteDuplicates;
      paired = Union[paired, Intersection[metricSyms, tensorDownSyms]]
    ];
    If[Length[paired] > 0,
      Print["  Tracing ", Length[paired], " abstract dummy pairs (per-term)"];
      (* Per-term tracing: each additive term has only ~3 abstract pairs,  *)
      (* not 45 across all terms. Sum-trace each term independently to     *)
      (* avoid O(dim^{2K}) explosion on the full expression.               *)
      (* Mathematically safe: Einstein summation dummies are always         *)
      (* contracted WITHIN a single product term, never across additive    *)
      (* terms. Ref: BatchedTraceBasisDummyWithMetric (same principle).    *)
      Module[{terms, tracedTerms = {}, tTerm},
        terms = If[Head[e] === Plus, List @@ e, {e}];
        Do[
          Module[{term = terms[[k]], termPaired, upT, downT, replaced},
            (* Find abstract dummy pairs IN THIS TERM *)
            upT = Cases[term, s_Symbol /; AbstractIndexQ[s] &&
              !MemberQ[coordNames, ToString[s]], {0, Infinity}] // DeleteDuplicates;
            downT = {};
            downT = Union[downT,
              Cases[term, (f_?CovDQ)[idx_][_] /; DownIndexQ[idx] :> UpIndex[idx],
                {0, Infinity}]];
            downT = Union[downT,
              Cases[term, (f_?xTensorQ)[___, -s_Symbol, ___] /; AbstractIndexQ[s] :> s,
                {0, Infinity}]];
            downT = Union[downT,
              Cases[term, idx_ /; DownIndexQ[idx] && AbstractIndexQ[UpIndex[idx]] :> UpIndex[idx],
                {0, Infinity}]];
            downT = Select[downT, AbstractIndexQ];
            termPaired = Intersection[upT, downT];
            If[Length[termPaired] > 0,
              replaced = term;
              Do[
                replaced = Sum[
                  replaced /. {p -> {i, chart}, ChangeIndex[p] -> {i, -chart}},
                  {i, 0, dimLocal - 1}
                ];
                replaced = Expand[replaced];
                If[!computeChristoffels, replaced = replaced /. christoffelPD -> Zero];
                replaced = ToValues[replaced];
                replaced = Expand[replaced],
                {p, termPaired}
              ];
              AppendTo[tracedTerms, replaced],
              AppendTo[tracedTerms, term]
            ]
          ],
          {k, Length[terms]}
        ];
        e = Total[tracedTerms];
      ];
      e = Expand[e];
      If[Head[e] === Plus,
        e = Total[TraceBasisDummy /@ List @@ e],
        e = TraceBasisDummy[e]
      ];
      e = ToValues[e];
      e = ToValues[e];
    ]
  ];

  (* CD[scalar] → Derivative: After ToValues, field tensors are replaced
     by named scalar functions (e.g. H1[t,x,y]). Covariant derivatives of
     scalars equal partial derivatives (no Christoffel correction):
       CD[-{idx,chart}][f[t,x,y,...]] → Derivative[0,...,1,...,0][f][t,x,y,...]
     where the 1 is at position idx+1. Goes directly to Derivative form
     instead of intermediate PD (which isCDlikeQ excludes from conversion).
     This applies to both flat and curved spacetimes.
     Ref: Wald (1984), eq. 3.1.15: ∇_a f = ∂_a f for scalar f. *)
  Module[{dim = DimOfManifold[ManifoldOfChart[chart]]},
    e = e /. (f_)[{idx_Integer, chartSign_}][g_[args___]] /;
        CovDQ[f] && FreeQ[g, _?xTensorQ] && (chartSign === chart || chartSign === -chart) :>
      With[{orders = ReplacePart[ConstantArray[0, dim], idx + 1 -> 1]},
        Derivative[Sequence @@ orders][g][args]
      ];
    (* Also handle CD applied to existing Derivative of scalar functions *)
    e = e /. (f_)[{idx_Integer, chartSign_}][Derivative[orders__][g_][args__]] /;
        CovDQ[f] && FreeQ[g, _?xTensorQ] && (chartSign === chart || chartSign === -chart) :>
      With[{paddedOrders = PadRight[{orders}, dim, 0]},
        With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
          Derivative[Sequence @@ newOrders][g][args]
        ]
      ];
    (* Handle Derivative wrapping CD: Derivative[...][CD[{idx,±chart}]][g[args]]  *)
    (* This arises from abstract dummy Sum tracing where CD indices get replaced  *)
    (* before the CD→Derivative conversion, creating inverted nesting.           *)
    e = e /. Derivative[outerOrds__][(f_)[{idx_Integer, chartSign_}]][g_Symbol[args__]] /;
        CovDQ[f] && (chartSign === chart || chartSign === -chart) :>
      With[{paddedOuter = PadRight[{outerOrds}, dim, 0]},
        With[{mergedOrders = ReplacePart[paddedOuter, idx + 1 -> paddedOuter[[idx + 1]] + 1]},
          Derivative[Sequence @@ mergedOrders][g][args]
        ]
      ];
  ];

  On[Validate::repeated];
  e
];


(* === SplinterToArray: term-by-term projection to component array ===     *)
(*                                                                          *)
(* Projects a single abstract-tensor term (not a sum) through the full      *)
(* splinter pipeline INCLUDING ComponentArray, returning a dim^rank array   *)
(* of concrete scalar expressions.  This is the key technique from the      *)
(* supervisor's EuclideanSplinter (commit 4a89164): ComponentArray converts *)
(* ALL abstract indices to concrete basis indices simultaneously, which     *)
(* enables ToValues to resolve shorthand tensor ComponentValues.            *)
(*                                                                          *)
(* StaggeredToBasis uses TraceBasisDummy instead of ComponentArray, which   *)
(* cannot resolve contracted abstract indices inside opaque CD wrappers.    *)
(* SplinterToArray avoids this by atomically expanding all indices before   *)
(* any dummy-index summation.                                               *)
(*                                                                          *)
(* Used by: DecomposeToComponents (term-by-term branch for R̃² torsion).   *)
(* Ref: supervisor's EuclideanSplinter in SphericalEuclidean.m lines 239.  *)

SplinterToArray[term_, chart_, computeChristoffels_:False] := Module[
  {e = term, metric, covd, christoffelPD},

  (* Get metric and Christoffel symbol *)
  Module[{manifold, tangentBundle, metrics},
    manifold = ManifoldOfChart[chart];
    tangentBundle = Symbol["Tangent" <> ToString[manifold]];
    metrics = MetricsOfVBundle[tangentBundle];
    If[!ListQ[metrics] || Length[metrics] == 0,
      Print["WARNING: SplinterToArray: no metric for chart ", chart];
      Return[ToBasis[chart][e]]
    ];
    metric = First[metrics];
    covd = CovDOfMetric[metric];
  ];
  christoffelPD = Symbol["Christoffel" <> ToString[covd] <> "PD" <> ToString[chart]];

  (* Staggered ToBasis (same as StaggeredToBasis but without TraceBasisDummy) *)
  If[Head[e] === Plus,
    e = Total[ToBasis[chart] /@ List @@ e],
    e = ToBasis[chart][e]
  ];
  If[!computeChristoffels, e = e /. christoffelPD -> Zero, e = ToValues[e]];
  If[Head[e] === Plus,
    e = Total[ToBasis[chart] /@ List @@ e],
    e = ToBasis[chart][e]
  ];
  If[!computeChristoffels, e = e /. christoffelPD -> Zero, e = ToValues[e]];
  e = SeparateMetric[metric][e];
  If[Head[e] === Plus,
    e = Total[ToBasis[chart] /@ List @@ e],
    e = ToBasis[chart][e]
  ];
  If[!computeChristoffels, e = e /. christoffelPD -> Zero, e = ToValues[e]];

  (* KEY: ComponentArray atomically expands ALL abstract indices into a     *)
  (* dim^rank array.  This is what makes shorthand ComponentValues resolve  *)
  (* — ToValues can match the concrete basis-index patterns.               *)
  (* Ref: supervisor's EuclideanSplinter lines 246-249.                    *)
  e = ComponentArray[e];
  e = ToValues[e];
  e = TraceBasisDummy[e];
  e = ToValues[e];

  (* CD[scalar] → Derivative: After ToValues, field tensors are replaced  *)
  (* by named scalar functions. CovD of scalars = partial derivatives.    *)
  (* Must be applied to each cell of the ComponentArray result.           *)
  Module[{dim = DimOfManifold[ManifoldOfChart[chart]]},
    e = Map[Function[{cell},
      cell /. (f_)[{idx_Integer, chartSign_}][g_[args___]] /;
          CovDQ[f] && FreeQ[g, _?xTensorQ] && (chartSign === chart || chartSign === -chart) :>
        With[{orders = ReplacePart[ConstantArray[0, dim], idx + 1 -> 1]},
          Derivative[Sequence @@ orders][g][args]
        ]
    ], e, {-1}];
    e = Map[Function[{cell},
      cell /. (f_)[{idx_Integer, chartSign_}][Derivative[orders__][g_][args__]] /;
          CovDQ[f] && FreeQ[g, _?xTensorQ] && (chartSign === chart || chartSign === -chart) :>
        With[{paddedOrders = PadRight[{orders}, dim, 0]},
          With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
            Derivative[Sequence @@ newOrders][g][args]
          ]
        ]
    ], e, {-1}];
    (* Also handle Derivative-wrapping-CD *)
    e = Map[Function[{cell},
      cell /. Derivative[outerOrds__][(f_)[{idx_Integer, chartSign_}]][g_Symbol[args__]] /;
          CovDQ[f] && (chartSign === chart || chartSign === -chart) :>
        With[{paddedOuter = PadRight[{outerOrds}, dim, 0]},
          With[{mergedOrders = ReplacePart[paddedOuter, idx + 1 -> paddedOuter[[idx + 1]] + 1]},
            Derivative[Sequence @@ mergedOrders][g][args]
          ]
        ]
    ], e, {-1}];
  ];

  e
];


(* === Component Decomposition === *)

(* Options for DecomposeToComponents *)
Options[DecomposeToComponents] = {
  "ComputeChristoffels" -> Automatic,  (* Automatic (default), True, or False *)
  "MetricMatrix" -> None,  (* Explicit metric matrix for curved spacetime evaluation *)
  "SkipTuples" -> {},  (* Component index tuples to skip (e.g. TT-zeroed {0,mu}) *)
  "BackgroundFieldRules" -> {},  (* List of {fieldHead, {comp0, comp1, ...}} for background field eval *)
  "KilledAxes" -> {},  (* Basis indices of transverse axes for plane-wave reduction *)
  "PertFieldHeads" -> {},  (* Perturbation field head symbols whose transverse PD should be zeroed *)
  "TermByTerm" -> False  (* Use term-by-term projection via SplinterToArray (for R̃² torsion) *)
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

  (* Suppress Validate::inhom for R̃²-decomposed torsion expressions.     *)
  (* DummyIn-allocated indices in Scalar[x]^n expansion create different   *)
  (* dummy names per Scalar factor; xAct's Validate sees this as           *)
  (* inhomogeneous in a sum. The indices ARE properly contracted inside     *)
  (* each Scalar — the "inhomogeneity" is cosmetic, not a real error.     *)
  (* Off suppresses the Message; wrapping critical calls in Catch absorbs  *)
  (* the Throw[Null] from xAct's UncatchedValidate.                       *)
  Off[Validate::inhom];

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
    componentEq = SeparateFieldMetrics[componentEq, chart];

    (* Expand Scalar[] wrappers before staggered ToBasis *)
    componentEq = ExpandScalarWrappers[componentEq, chart, computeChristoffels =!= False];

    (* Staggered ToBasis pipeline (replaces old ToBasis + TBD + steps 4-7) *)
    componentEq = StaggeredToBasis[componentEq, chart, computeChristoffels =!= False];
    componentEq = Expand[componentEq];

    (* Background field evaluation — not handled by MetricCompute/ToValues *)
    If[backgroundFieldRules =!= {},
      Do[
        componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
        {bg, backgroundFieldRules}
      ];
      componentEq = Expand[componentEq]
    ];

    (* Replace ALL fields with named scalar functions *)
    Module[{coordSyms},
      coordSyms = GetCoordinateSymbols[chart];
      Do[
        componentEq = ReplaceTensorFieldComponents[componentEq, afh, chart, coordSyms, dim, metricMatrix],
        {afh, allFieldHeads}
      ];

      (* Convert coordinate derivatives to explicit Derivative form *)
      (* Catch Throw from ConvertCDToDerivatives::incomplete — the fallback *)
      (* CD→Derivative rules in the caller handle residual operators.       *)
      Module[{cdResult},
        cdResult = Catch[ConvertCDToDerivatives[componentEq, chart]];
        If[!StringQ[cdResult],
          componentEq = cdResult
          (* else: keep componentEq as-is, fallback below handles it *)
        ]
      ];

      (* Fallback CD→Derivative for residual CovD operators *)
      Module[{dim2 = GetChartDimension[chart], prevExpr2, iter2fb = 0},
        prevExpr2 = componentEq;
        While[iter2fb < 10,
          componentEq = componentEq /. (f_)[{idx_Integer, s_}][g_Symbol[args___]] /;
              CovDQ[f] && (s === chart || s === -chart) :>
            With[{orders = ReplacePart[ConstantArray[0, dim2], idx + 1 -> 1]},
              Derivative[Sequence @@ orders][g][args]
            ];
          componentEq = componentEq /. (f_)[{idx_Integer, s_}][Derivative[orders__][g_][args__]] /;
              CovDQ[f] && (s === chart || s === -chart) :>
            With[{paddedOrders = PadRight[{orders}, dim2, 0]},
              With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
                Derivative[Sequence @@ newOrders][g][args]
              ]
            ];
          componentEq = componentEq /. Derivative[outerOrds__][(f_)[{idx_Integer, s_}]][g_Symbol[args__]] /;
              CovDQ[f] && (s === chart || s === -chart) :>
            With[{paddedOuter = PadRight[{outerOrds}, dim2, 0]},
              With[{mergedOrders = ReplacePart[paddedOuter, idx + 1 -> paddedOuter[[idx + 1]] + 1]},
                Derivative[Sequence @@ mergedOrders][g][args]
              ]
            ];
          iter2fb++;
          If[componentEq === prevExpr2, Break[]];
          prevExpr2 = componentEq;
        ];
      ];

      (* Post-ConvertCDToDerivatives: catch residual BG Derivative forms *)
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

      (* Pre-expand Scalar[] wrappers ONCE on the full EOM before the     *)
      (* per-component loop.  Without hoisting, ExpandScalarWrappers runs *)
      (* inside ExtractTensorComponent for EACH component, redundantly    *)
      (* re-expanding the same Scalar contents dim^rank times.            *)
      (* For R̃-decomposed torsion theories: 6x speedup (h) + 9x (t).   *)
      (* ExpandScalarWrappers: resolve Scalar[...] wrappers from xPert.
         For standard theories: hoist before the component loop (6x speedup).
         For R̃-decomposed torsion: do NOT hoist — the resolved Scalar
         contents share abstract index names with the EOM free indices,
         causing index collisions. Let ExpandScalarWrappers run per-component
         inside ExtractTensorComponent (step 1.3), AFTER free indices are
         fixed to specific basis values. *)
      If[!FreeQ[eomSep, Scalar],
        (* Check if expression has Christoffel symbols from R̃ decomposition *)
        Module[{hasTorsionScalars},
          hasTorsionScalars = !FreeQ[eomSep, Scalar[x_ /; !FreeQ[x, _?CovDQ]]];
          If[!hasTorsionScalars,
            (* Safe to hoist: standard theory without CD inside Scalar *)
            Print["  Pre-expanding Scalar[] wrappers (hoisted)..."];
            eomSep = ExpandScalarWrappers[eomSep, chart, computeChristoffels];
            Print["  Scalar expansion complete: ", If[FreeQ[eomSep, Scalar], "all resolved", "some remain"]],
            (* Torsion theory: skip hoisting, let per-component expansion handle it *)
            Print["  Scalar[] wrappers detected with CD operators — expanding per-component"]
          ]
        ]
      ];

      (* Suppress Validate::repeated AND its Throw through the component     *)
      (* loop. R̃-decomposed expressions have abstract index collisions     *)
      (* between resolved Scalar contents and outer EOM indices that are   *)
      (* false positives resolved by ToBasis in StaggeredToBasis.          *)
      (* xAct's Validate does Message[...] + Throw[Null]; Off suppresses   *)
      (* the message but NOT the Throw, so we also need Catch around the   *)
      (* loop (and override ValidateIndices to prevent the Throw).         *)
      Off[Validate::repeated];
      Off[Validate::inhom];
      Unprotect[xAct`xTensor`Private`ValidateIndices];
      Module[{savedValidate = xAct`xTensor`Private`ValidateIndices},
        xAct`xTensor`Private`ValidateIndices = (True &);

        (* ================================================================ *)
        (* Term-by-term projection (supervisor's approach, commit 4a89164). *)
        (* Split the abstract EOM into additive terms and project each      *)
        (* individually through SplinterToArray (which uses ComponentArray  *)
        (* to atomically expand ALL abstract indices). This resolves CD     *)
        (* shorthand ComponentValues that TraceBasisDummy alone cannot.     *)
        (* Ref: SphericalEuclidean.m lines 415-428.                        *)
        (* ================================================================ *)
        If[OptionValue["TermByTerm"] === True,
          Module[{eomTerms, nTerms, termComp, eomComp,
                  fieldRankLocal, dimLocal, tStart, coordSymsLocal,
                  allFieldHeadsLocal},
            (* ================================================================ *)
            (* Supervisor's approach (SphericalEuclidean.m:415-428):             *)
            (* Split the abstract EOM (WITH free indices) into additive terms.  *)
            (* Project each term through SplinterToArray (which uses            *)
            (* ComponentArray to atomically expand ALL indices to basis form).  *)
            (* Accumulate into a dim^rank array. Post-process each component.  *)
            (*                                                                   *)
            (* Why this is fast: ComponentArray is O(1) per tensor, vs          *)
            (* TraceBasisDummy which is O(dim^{2K}) for K contracted pairs.    *)
            (* Scalar wrappers must be pre-resolved (commit 7666cd2) before    *)
            (* this branch runs, to avoid DummyIn index contamination.         *)
            (* ================================================================ *)
            Print["  Term-by-term projection (SplinterToArray + ComponentArray)..."];
            fieldRankLocal = fieldRank;
            dimLocal = dim;
            coordSymsLocal = GetCoordinateSymbols[chart];
            allFieldHeadsLocal = Join[{fieldHead}, ExtractFieldHead /@ additionalFields];

            (* Split abstract EOM into additive terms.                    *)
            (* Peel CollectTensors wrappers (from ToCanonical) and Expand  *)
            (* to expose the Plus structure for term-by-term projection.  *)
            Module[{eomForSplit = eomSep},
              While[Head[eomForSplit] =!= Plus && Length[eomForSplit] == 1,
                eomForSplit = eomForSplit[[1]]
              ];
              If[Head[eomForSplit] =!= Plus,
                eomForSplit = Expand[eomForSplit];
                (* Re-peel after Expand *)
                While[Head[eomForSplit] =!= Plus && Length[eomForSplit] == 1,
                  eomForSplit = eomForSplit[[1]]
                ]
              ];
              eomTerms = If[Head[eomForSplit] === Plus,
                List @@ eomForSplit, {eomForSplit}];
            ];
            nTerms = Length[eomTerms];
            Print["  EOM has ", nTerms, " additive terms"];

            (* Accumulate component arrays term by term *)
            eomComp = ConstantArray[0, Table[dimLocal, {fieldRankLocal}]];
            Do[
              tStart = AbsoluteTime[];
              termComp = Catch[SplinterToArray[eomTerms[[k]], chart, computeChristoffels]];
              If[StringQ[termComp],
                Print["  Term ", k, "/", nTerms, ": FAILED (", termComp, ")"];
                Continue[]
              ];
              eomComp += termComp;
              If[Mod[k, 5] == 0 || k == nTerms,
                Print["  Term ", k, "/", nTerms, " (",
                      Round[AbsoluteTime[] - tStart, 0.1], "s, ",
                      Round[MemoryInUse[]/1024.^2], " MB)"]
              ];
              Share[];
            , {k, nTerms}];

            (* Simplify accumulated component array *)
            Print["  Simplifying component array..."];
            eomComp = Map[Simplify, eomComp, {fieldRankLocal}];

            (* Post-process each component: field replacement + CD→Derivative *)
            Print["  Post-processing components..."];
            result = {};
            Do[
              Module[{tuple, flatIdx, compExpr},
                tuple = componentTuples[[idx]];
                flatIdx = flatIdxMap[tuple];
                compExpr = eomComp[[Sequence @@ (tuple + 1)]];

                (* Replace tensor fields with named scalar functions *)
                Do[
                  compExpr = ReplaceTensorFieldComponents[compExpr, afh, chart,
                    coordSymsLocal, dimLocal, metricMatrix],
                  {afh, allFieldHeadsLocal}
                ];

                (* Convert CD → Derivative *)
                Module[{cdResult},
                  cdResult = Catch[ConvertCDToDerivatives[compExpr, chart]];
                  If[!StringQ[cdResult], compExpr = cdResult]
                ];
                (* Fallback CD→Derivative *)
                Module[{dimFB = dimLocal, prevFB, iterFB = 0},
                  prevFB = compExpr;
                  While[iterFB < 10,
                    compExpr = compExpr /. (f_)[{idxI_Integer, s_}][g_Symbol[args___]] /;
                        CovDQ[f] && (s === chart || s === -chart) :>
                      With[{orders = ReplacePart[ConstantArray[0, dimFB], idxI + 1 -> 1]},
                        Derivative[Sequence @@ orders][g][args]
                      ];
                    compExpr = compExpr /. (f_)[{idxI_Integer, s_}][Derivative[orders__][g_][args__]] /;
                        CovDQ[f] && (s === chart || s === -chart) :>
                      With[{paddedOrders = PadRight[{orders}, dimFB, 0]},
                        With[{newOrders = ReplacePart[paddedOrders, idxI + 1 -> paddedOrders[[idxI + 1]] + 1]},
                          Derivative[Sequence @@ newOrders][g][args]
                        ]
                      ];
                    (* Derivative-wrapping-CD *)
                    compExpr = compExpr /. Derivative[outerOrds__][(f_)[{idxI_Integer, s_}]][g_Symbol[args__]] /;
                        CovDQ[f] && (s === chart || s === -chart) :>
                      With[{paddedOuter = PadRight[{outerOrds}, dimFB, 0]},
                        With[{mergedOrders = ReplacePart[paddedOuter, idxI + 1 -> paddedOuter[[idxI + 1]] + 1]},
                          Derivative[Sequence @@ mergedOrders][g][args]
                        ]
                      ];
                    iterFB++;
                    If[compExpr === prevFB, Break[]];
                    prevFB = compExpr;
                  ];
                ];
                If[backgroundFieldRules =!= {},
                  Do[compExpr = EvaluatePDBackgroundField[compExpr, chart, bg[[1]], bg[[2]],
                      If[Length[bg] >= 3, bg[[3]], {}]],
                    {bg, backgroundFieldRules}];
                  compExpr = Expand[compExpr]
                ];

                AppendTo[result, {flatIdx, Expand[Simplify[compExpr]]}];
                Print["  [", Round[MemoryInUse[]/1024.^2], " MB] component ",
                      idx, "/", Length[componentTuples]];
              ],
              {idx, 1, Length[componentTuples]}
            ];
          ];
          (* Clean up and return *)
          xAct`xTensor`Private`ValidateIndices = savedValidate;
          Protect[xAct`xTensor`Private`ValidateIndices];
          On[Validate::repeated];
          On[Validate::inhom];
          Return[result]
        ];

        (* ================================================================ *)
        (* Standard per-component extraction (ExtractTensorComponent).      *)
        (* Used for non-torsion theories or when TermByTerm is not needed.  *)
        (* ================================================================ *)
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
        xAct`xTensor`Private`ValidateIndices = savedValidate;
      ];
      Protect[xAct`xTensor`Private`ValidateIndices];
      On[Validate::repeated];
      On[Validate::inhom];
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

    (* Continuous adaptive batching: re-measure expansion factor EVERY batch. *)
    (* Previous approach sampled only the first batch, which fails when      *)
    (* later batches have different coupling structure (e.g., cross-field    *)
    (* terms with higher K). Only DECREASE — never increase above default.  *)
    (* Expand[] is super-linear: 2x batch → 10-20x time.                   *)
    (* Ref: issue #198 — first-batch sampling breaks for multi-field theories. *)
    If[peakTerms > 0,
      expansionFactor = N[peakTerms / (batchEnd - batchStart + 1)];
      If[expansionFactor > 0,
        currentBatchSize = Max[5, Floor[targetPeak / expansionFactor]];
        (* Never increase above default — Expand cost is super-linear *)
        currentBatchSize = Min[currentBatchSize, batchSize];
        If[currentBatchSize != (batchEnd - batchStart + 1),
          Print["    adaptive batch: expansion=",
                Round[expansionFactor, 0.1], "x, batch -> ", currentBatchSize,
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

  (* Step 1.3: Expand Scalar[] wrappers before ToBasis.                    *)
  componentEq = ExpandScalarWrappers[componentEq, chart, computeChristoffels];

  (* Step 1.4: Re-apply shorthand CD[field] → CDfield substitutions.       *)
  (* ExpandScalarWrappers may have introduced new CD[field] operators from  *)
  (* resolved Scalar contents. The shorthand MakeRule variables (toCD1...)  *)
  (* are globals set by _wls_shorthand_cd_tensors in _derive.py.           *)
  If[ListQ[Global`$CDShorthandRules] && Length[Global`$CDShorthandRules] > 0,
    Do[componentEq = componentEq /. rule, {rule, Global`$CDShorthandRules}];
    componentEq = ToCanonical[componentEq];
    componentEq = ContractMetric[componentEq];
  ];

  (* Step 2: Staggered ToBasis + ToValues + TraceBasisDummy                *)
  (* Ref: supervisor's EuclideanSplinter (commit 4a89164).                 *)
  (* Replaces the old single-ToBasis + BatchedTraceBasisDummy + steps 4-7  *)
  (* with a unified pipeline that pre-evaluates ComponentValues BEFORE     *)
  (* TraceBasisDummy, dramatically reducing the dummy-index enumeration.   *)
  componentEq = StaggeredToBasis[componentEq, chart, computeChristoffels];
  componentEq = Expand[componentEq];
  Print["    step2-Splinter: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentEq]===Plus, Length[componentEq], 1], " terms"];

  (* Background field evaluation — not handled by MetricCompute/ToValues *)
  If[backgroundFieldRules =!= {},
    Do[
      componentEq = EvaluatePDBackgroundField[componentEq, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentEq = Expand[componentEq];
    Print["    step-BG: ", Round[MemoryInUse[]/1024.^2], " MB, ",
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
  (* Catch Throw from ConvertCDToDerivatives::incomplete — fallback below *)
  Module[{cdResult9},
    cdResult9 = Catch[ConvertCDToDerivatives[componentEq, chart]];
    If[!StringQ[cdResult9], componentEq = cdResult9]
  ];

  (* Step 9a: Fallback CD→Derivative for residual CovD operators *)
  Module[{dim9 = GetChartDimension[chart], prevExpr9, iter9 = 0},
    prevExpr9 = componentEq;
    While[iter9 < 10,
      componentEq = componentEq /. (f_)[{idx_Integer, s_}][g_Symbol[args___]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{orders = ReplacePart[ConstantArray[0, dim9], idx + 1 -> 1]},
          Derivative[Sequence @@ orders][g][args]
        ];
      componentEq = componentEq /. (f_)[{idx_Integer, s_}][Derivative[orders__][g_][args__]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{paddedOrders = PadRight[{orders}, dim9, 0]},
          With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
            Derivative[Sequence @@ newOrders][g][args]
          ]
        ];
      componentEq = componentEq /. Derivative[outerOrds__][(f_)[{idx_Integer, s_}]][g_Symbol[args__]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{paddedOuter = PadRight[{outerOrds}, dim9, 0]},
          With[{mergedOrders = ReplacePart[paddedOuter, idx + 1 -> paddedOuter[[idx + 1]] + 1]},
            Derivative[Sequence @@ mergedOrders][g][args]
          ]
        ];
      iter9++;
      If[componentEq === prevExpr9, Break[]];
      prevExpr9 = componentEq;
    ];
  ];

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
        ];

        (* Derivative-wrapped forms: Derivative[orders][fh][{basis}]         *)
        (* Arise when D[] differentiates through tensor components that      *)
        (* haven't been replaced with scalar functions yet.                  *)
        If[!FreeQ[result, Derivative[__][fh]],
          Do[
            result = result /. {
              Derivative[orders__][fh][{pair[[1]], sc[[1]]}, {pair[[2]], sc[[2]]}] :>
                Derivative[orders][sym][Sequence @@ cs]
            };
            If[symQ && pair[[1]] =!= pair[[2]],
              result = result /. {
                Derivative[orders__][fh][{pair[[2]], sc[[1]]}, {pair[[1]], sc[[2]]}] :>
                  Derivative[orders][sym][Sequence @@ cs]
              }
            ],
            {sc, Tuples[{ch, -ch}, 2]}
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
  {result = expr, rank, canonTuples, allTuples, symGroup, gens,
   canonMap, signMap},

  rank = Length[SlotsOfTensor[fh]];
  canonTuples = EnumerateComponentTuples[fh, dim];
  allTuples = Tuples[Range[0, dim - 1], rank];

  (* Build canonical map: for EVERY tuple, find its canonical representative *)
  (* and the sign factor from the symmetry group.                           *)
  (* Example: T with Antisymmetric[{2,3}]:                                  *)
  (*   {0,0,1} → canonical {0,0,1}, sign +1                                *)
  (*   {0,1,0} → canonical {0,0,1}, sign -1 (antisymmetric swap)           *)
  (*   {0,0,0} → zero (antisymmetric, equal indices)                        *)
  symGroup = SymmetryGroupOfTensor[fh];
  gens = If[symGroup =!= StrongGenSet[{}, GenSet[]],
    List @@ symGroup[[2]], {}];

  (* For each ALL-tuple, determine: is it canonical? If not, what is the    *)
  (* canonical tuple and sign? For simplicity, check each generator.        *)
  canonMap = Association[];
  signMap = Association[];
  Do[
    Module[{tup = allTuples[[j]], isCanon, canonIdx, sign = 1},
      (* Check if this tuple is zero (e.g., antisymmetric with equal indices) *)
      isCanon = isCanonicalTuple[tup, gens];
      If[isCanon,
        (* This IS a canonical tuple — find its index *)
        canonIdx = FirstPosition[canonTuples, tup, {0}][[1]];
        If[canonIdx > 0,
          canonMap[tup] = canonIdx;
          signMap[tup] = 1
        ]
      ,
        (* Non-canonical: find the canonical partner via permutation *)
        Do[
          Module[{sp = extractGenSignAndPerm[gen], s, perm, permTup},
            s = sp[[1]]; perm = sp[[2]];
            permTup = applyPermToTuple[tup, perm];
            If[isCanonicalTuple[permTup, gens],
              Module[{idx = FirstPosition[canonTuples, permTup, {0}][[1]]},
                If[idx > 0,
                  canonMap[tup] = idx;
                  signMap[tup] = s;
                ]
              ]
            ]
          ],
          {gen, gens}
        ]
      ]
    ],
    {j, Length[allTuples]}
  ];

  (* Now replace: for each tuple (canonical or not), try all chart-sign    *)
  (* configurations. Replace with sign * canonicalScalar[coords].          *)
  With[{ch = chart, cs = coordSyms},
    Do[
      Module[{tup = allTuples[[j]]},
        If[KeyExistsQ[canonMap, tup],
          Module[{sym = Symbol[ToString[fh] <> ToString[canonMap[tup] - 1]],
                  sign = signMap[tup]},
            (* Bare tensor: all chart-sign configs *)
            If[!FreeQ[result, fh],
              Do[
                Module[{pat = Table[{tup[[n]], config[[n]]}, {n, rank}]},
                  result = result /. {fh @@ pat :> sign * sym[Sequence @@ cs]}
                ],
                {config, Tuples[{ch, -ch}, rank]}
              ]
            ];
            (* Derivative-wrapped forms *)
            If[!FreeQ[result, Derivative[__][fh]],
              Do[
                Module[{pat = Table[{tup[[n]], config[[n]]}, {n, rank}]},
                  result = result /. {
                    Derivative[orders__][fh][Sequence @@ pat] :>
                      sign * Derivative[orders][sym][Sequence @@ cs]
                  }
                ],
                {config, Tuples[{ch, -ch}, rank]}
              ]
            ]
          ]
        ,
          (* Tuple not in canonMap = zero by symmetry (e.g., antisymmetric *)
          (* with equal indices). Replace with 0.                         *)
          If[!FreeQ[result, fh],
            Do[
              Module[{pat = Table[{tup[[n]], config[[n]]}, {n, rank}]},
                result = result /. {fh @@ pat :> 0}
              ],
              {config, Tuples[{ch, -ch}, rank]}
            ]
          ];
          If[!FreeQ[result, Derivative[__][fh]],
            Do[
              Module[{pat = Table[{tup[[n]], config[[n]]}, {n, rank}]},
                result = result /. {
                  Derivative[orders__][fh][Sequence @@ pat] :> 0
                }
              ],
              {config, Tuples[{ch, -ch}, rank]}
            ]
          ]
        ]
      ],
      {j, Length[allTuples]}
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

  (* Expand Scalar[] wrappers.  Gauge-fixing terms and R̃-decomposed torsion  *)
  (* produce Scalar[CD_a V^a]^2 and Scalar[eta[a,b]*CD[-a][T[-b,-c,-d]]]    *)
  (* that ToBasis cannot penetrate.  Uses shared ExpandScalarWrappers with   *)
  (* bounded iteration for safety.                                           *)
  componentExpr = ExpandScalarWrappers[componentExpr, chart, shouldComputeChristoffels];

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

  (* Step 2: Staggered ToBasis + ToValues + TraceBasisDummy                *)
  (* Ref: supervisor's EuclideanSplinter (commit 4a89164).                 *)
  (* With pre-computed CD shorthand ComponentValues, ToValues resolves     *)
  (* CDN-field tensors in O(1), reducing TraceBasisDummy dummy pairs.      *)
  (*                                                                        *)
  (* Batched processing for large expressions (>50 terms): process in      *)
  (* adaptive batches through StaggeredToBasis to prevent O(dim^{2K})     *)
  (* blowup when cross-field coupling creates high-K terms.               *)
  (* Continuous adaptive: re-measure each batch, not just the first.       *)
  (* Ref: issue #198 — single-pass fails for multi-field theories.        *)
  If[Head[componentExpr] === Plus && Length[componentExpr] > 50,
    Module[{terms, scalarResult = 0, scalarBatch, scalarBatchResult,
            scalarBatchN = 50, tScalarBatch, scalarBatchTime},
      terms = List @@ componentExpr;
      Print["    [scalar] batched Splinter: ", Length[terms], " terms"];
      Do[
        scalarBatch = Total[terms[[i ;; Min[i + scalarBatchN - 1, Length[terms]]]]];
        tScalarBatch = AbsoluteTime[];
        scalarBatchResult = StaggeredToBasis[scalarBatch, chart, computeChristoffels];
        scalarBatchResult = Expand[scalarBatchResult];
        scalarBatchTime = AbsoluteTime[] - tScalarBatch;
        Print["      batch[", i, ":", Min[i + scalarBatchN - 1, Length[terms]],
              "/", Length[terms], "]: ",
              If[Head[scalarBatchResult]===Plus, Length[scalarBatchResult], 1],
              " terms (", Round[scalarBatchTime, 0.1], "s)"];
        (* Continuous adaptive: halve batch if slow, double if fast *)
        If[scalarBatchTime > 10.0, scalarBatchN = Max[5, Floor[scalarBatchN / 2]]];
        If[scalarBatchTime < 1.0 && scalarBatchN < 50,
          scalarBatchN = Min[50, scalarBatchN * 2]];
        scalarResult += scalarBatchResult;
        Share[];,
        {i, 1, Length[terms], scalarBatchN}
      ];
      componentExpr = scalarResult;
    ],
    componentExpr = StaggeredToBasis[componentExpr, chart, computeChristoffels];
    componentExpr = Expand[componentExpr];
  ];
  Print["    [scalar] step2-Splinter: ", Round[MemoryInUse[]/1024.^2], " MB, ",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"];
  (* Diagnostic: check for residual abstract indices *)
  If[!FreeQ[componentExpr, _?AbstractIndexQ],
    Module[{abstractTerms},
      Print["    WARNING: residual abstract indices after Splinter"];
      abstractTerms = If[Head[componentExpr]===Plus,
        Select[List @@ componentExpr, !FreeQ[#, _?AbstractIndexQ]&],
        {componentExpr}
      ];
      Print["    Abstract terms: ", Length[abstractTerms], "/",
        If[Head[componentExpr]===Plus, Length[componentExpr], 1]];
      If[Length[abstractTerms] > 0,
        Print["    First abstract term: ", Short[abstractTerms[[1]], 3]]
      ];
    ]
  ];

  (* Background field evaluation — not handled by MetricCompute/ToValues *)
  If[backgroundFieldRules =!= {},
    Do[
      componentExpr = EvaluatePDBackgroundField[componentExpr, chart, bg[[1]], bg[[2]], If[Length[bg] >= 3, bg[[3]], {}]],
      {bg, backgroundFieldRules}
    ];
    componentExpr = Expand[componentExpr];
    Print["    [scalar] step-BG: ", Round[MemoryInUse[]/1024.^2], " MB, ",
          If[Head[componentExpr]===Plus, Length[componentExpr], 1], " terms"]
  ];

  (* Step 8: Replace ALL tensor fields with named scalar component functions *)
  Do[
    componentExpr = ReplaceTensorFieldComponents[componentExpr, afh, chart, coordSyms, dim, metricMatrix],
    {afh, allFieldHeads}
  ];

  (* Step 9: Convert coordinate derivatives to Derivative form *)
  (* Catch the Throw from ConvertCDToDerivatives::incomplete — the fallback *)
  (* in Step 9a handles any remaining CovD operators directly via CovDQ.   *)
  Module[{cdResult},
    cdResult = Catch[ConvertCDToDerivatives[componentExpr, chart]];
    If[StringQ[cdResult] && StringContainsQ[cdResult, "ConvertCDToDerivatives"],
      Print["  [scalar] ConvertCDToDerivatives incomplete — using fallback"];,
      componentExpr = cdResult
    ]
  ];

  (* Step 9a: Fallback CD→Derivative for any residual CovD operators.      *)
  (* ConvertCDToDerivatives may leave operators in nested structures where  *)
  (* isCDlikeQ-based pattern matching fails. Direct CovDQ check catches    *)
  (* these. This is safe: CovD of a scalar = partial derivative (Wald 3.1.15). *)
  Module[{dim2 = GetChartDimension[chart], prevExpr2, iter2 = 0},
    prevExpr2 = componentExpr;
    While[iter2 < 10,
      (* Rule 1: CD[idx][field[args]] → Derivative *)
      componentExpr = componentExpr /. (f_)[{idx_Integer, s_}][g_Symbol[args___]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{orders = ReplacePart[ConstantArray[0, dim2], idx + 1 -> 1]},
          Derivative[Sequence @@ orders][g][args]
        ];
      (* Rule 2: CD[idx][Derivative[...][field][args]] → merged Derivative *)
      componentExpr = componentExpr /. (f_)[{idx_Integer, s_}][Derivative[orders__][g_][args__]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{paddedOrders = PadRight[{orders}, dim2, 0]},
          With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
            Derivative[Sequence @@ newOrders][g][args]
          ]
        ];
      (* Rule 3: Derivative[orders][CD[idx]][field[args]] → merged Derivative *)
      (* Catches ∂_t(CD_x[field]) = Derivative[1,0,0][tidalCD[{1,-chart}]][field[...]] *)
      componentExpr = componentExpr /. Derivative[outerOrds__][(f_)[{idx_Integer, s_}]][g_Symbol[args__]] /;
          CovDQ[f] && (s === chart || s === -chart) :>
        With[{paddedOuter = PadRight[{outerOrds}, dim2, 0]},
          With[{mergedOrders = ReplacePart[paddedOuter, idx + 1 -> paddedOuter[[idx + 1]] + 1]},
            Derivative[Sequence @@ mergedOrders][g][args]
          ]
        ];
      iter2++;
      If[componentExpr === prevExpr2, Break[]];
      prevExpr2 = componentExpr;
    ];
  ];

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


(* === Component-Level Euler-Lagrange Equations ===                       *)
(*                                                                        *)
(* Compute field equations by differentiating the COMPONENT Lagrangian    *)
(* using standard variational calculus (Euler-Lagrange formula).          *)
(* This bypasses abstract-index VarD + DecomposeToComponents entirely,    *)
(* avoiding the O(dim^{2K}) TraceBasisDummy bottleneck.                   *)
(*                                                                        *)
(* For a quadratic Lagrangian L = (1/2) Σ q_i M_{ij} q_j, the EOM are   *)
(* linear: EOM_i = Σ_j M_{ij} q_j = 0.  The Euler-Lagrange formula      *)
(* with higher-order derivatives is:                                      *)
(*   EOM_i = Σ_{|α|≥0} (-1)^|α| D^α(∂L/∂(D^α q_i))                    *)
(* where D^α = ∂^{α_1}/∂x_1^{α_1} ... ∂^{α_n}/∂x_n^{α_n}.            *)
(*                                                                        *)
(* Ref: supervisor's approach decomposes L → components → project terms.  *)
(* This goes further: compute EOM at component level, avoiding all        *)
(* abstract-index issues.  Mathematically equivalent (E-L commutes with  *)
(* basis decomposition for complete bases).                               *)
(* ====================================================================== *)

ComponentEulerLagrange[lagrangian_, fieldFuncs_List, coords_List] := Module[
  {nFields = Length[fieldFuncs], nCoords = Length[coords], eom,
   q, L = lagrangian, derivTerms, maxOrder, contribution},

  Print["ComponentEulerLagrange: ", nFields, " fields, ", nCoords, " coordinates"];

  (* Collect all derivative orders appearing in the Lagrangian *)
  (* for each field function.  This determines the maximum     *)
  (* order of the E-L formula needed.                          *)
  eom = Table[0, {nFields}];

  Do[
    q = fieldFuncs[[i]];
    Print["  Field ", i, "/", nFields, ": ", q];

    (* Find all Derivative orders of this field in the Lagrangian *)
    derivTerms = Cases[L,
      Derivative[orders__][f_][args__] /; f === q :> {orders},
      {0, Infinity}] // DeleteDuplicates;

    (* Also include the undifferentiated field (identity operator) *)
    derivTerms = Join[{{Sequence @@ Table[0, nCoords]}}, derivTerms] // DeleteDuplicates;

    maxOrder = Max[Total /@ derivTerms];
    Print["    Derivative orders: ", Length[derivTerms],
          " unique, max total order: ", maxOrder];

    (* Apply the Euler-Lagrange formula:                              *)
    (* EOM_i = Σ_{α} (-1)^|α| * D^α[ ∂L/∂(D^α q_i) ]              *)
    (* where the sum is over all multi-indices α that appear in L.   *)
    contribution = 0;
    Do[
      Module[{alpha = orders, totalOrder, dLdDq, ibpResult},
        totalOrder = Total[alpha];

        (* ∂L/∂(D^α q_i): differentiate L w.r.t. the derivative term *)
        If[totalOrder == 0,
          (* Identity: ∂L/∂q_i *)
          dLdDq = D[L, q @@ coords],
          (* Derivative: ∂L/∂(D^α q_i) *)
          dLdDq = D[L, Derivative[Sequence @@ alpha][q][Sequence @@ coords]]
        ];

        If[dLdDq =!= 0,
          (* Apply (-1)^|α| * D^α to the result (integration by parts) *)
          ibpResult = dLdDq;
          Do[
            If[alpha[[coordIdx]] > 0,
              ibpResult = D[ibpResult, {coords[[coordIdx]], alpha[[coordIdx]]}]
            ],
            {coordIdx, nCoords}
          ];
          contribution += (-1)^totalOrder * ibpResult;
        ];
      ],
      {orders, derivTerms}
    ];

    eom[[i]] = Expand[contribution];
    Print["    EOM terms: ", If[Head[eom[[i]]] === Plus, Length[eom[[i]]], 1]];
  , {i, nFields}];

  (* Return as list of {fieldName, equation} pairs *)
  (* Field name is extracted from the function symbol (e.g. tidalH0 → "h_0") *)
  eom
];


End[];
EndPackage[];
