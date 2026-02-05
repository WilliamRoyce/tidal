(* ::Package:: *)
(* CommonUtilities.wl - Shared utilities for the Lagrangian-to-PDE pipeline *)
(* Part of the torsion-gertsenshtein project *)

BeginPackage["TorsionGertsenshtein`CommonUtilities`",
  {"xAct`xTensor`", "xAct`xCoba`"}];

(* Public symbols *)
RemoveChristoffelSymbols::usage =
  "RemoveChristoffelSymbols[expr] sets all Christoffel symbol terms to zero \
(valid for flat Minkowski space).";

EvaluateMinkowskiMetric::usage =
  "EvaluateMinkowskiMetric[expr, chart] evaluates metric components for \
Minkowski signature (-1, +1, ...) in the given chart.";

GetCoordinateSymbols::usage =
  "GetCoordinateSymbols[chart] returns the coordinate scalar symbols from the \
chart, with fallback to {t[], x[]} for 1+1D.";

ConvertCDToDerivatives::usage =
  "ConvertCDToDerivatives[expr, chart] converts covariant derivative operators \
(CD) to explicit Mathematica Derivative form.";

ExtractFieldHead::usage =
  "ExtractFieldHead[field] extracts the tensor head from an applied field \
like phi[] or A[-a].";

ExtractNumericCoefficient::usage =
  "ExtractNumericCoefficient[term, fieldName] extracts the numeric coefficient \
from a term containing the named field.";

EvaluateEpsilonComponents::usage =
  "EvaluateEpsilonComponents[expr, chart] evaluates Levi-Civita (epsilon) tensor \
components to numeric values ±1 for Minkowski signature (-,+,+,...).";

LeviCivitaValue::usage =
  "LeviCivitaValue[indices] returns the Levi-Civita symbol value: +1 for even \
permutations, -1 for odd permutations, 0 for repeated indices.";

Begin["`Private`"];

(* === Christoffel Symbol Removal === *)
(* In flat Minkowski space, Christoffel symbols vanish *)
(* Uses context-independent pattern matching to handle different xAct contexts *)

RemoveChristoffelSymbols[expr_] :=
  expr /. f_[__] /; StringMatchQ[ToString[f], "*Christoffel*"] -> 0;

(* === Minkowski Metric Evaluation === *)
(* Evaluates metric components for signature (-1, +1, +1, ...) *)
(* Handles both covariant indices {i, chart} and contravariant {i, -chart} *)

EvaluateMinkowskiMetric[expr_, chart_] :=
  expr /. {
    (* Diagonal components: -1 for time (index 0), +1 for space *)
    _Symbol[{0, chart}, {0, chart}] -> -1,
    _Symbol[{1, chart}, {1, chart}] -> 1,
    _Symbol[{2, chart}, {2, chart}] -> 1,
    _Symbol[{3, chart}, {3, chart}] -> 1,
    (* Also handle negative chart (covariant basis) *)
    _Symbol[{0, -chart}, {0, -chart}] -> -1,
    _Symbol[{1, -chart}, {1, -chart}] -> 1,
    _Symbol[{2, -chart}, {2, -chart}] -> 1,
    _Symbol[{3, -chart}, {3, -chart}] -> 1,
    (* Off-diagonal components: all zero (use conditional for efficiency) *)
    _Symbol[{i_Integer, chart}, {j_Integer, chart}] /; i != j -> 0,
    _Symbol[{i_Integer, -chart}, {j_Integer, -chart}] /; i != j -> 0
  };

(* === Coordinate Symbol Extraction === *)
(* Gets coordinate scalars from chart, with fallback for unevaluated cases *)

GetCoordinateSymbols[chart_] := Module[{coordSyms},
  coordSyms = ScalarsOfChart[chart];
  (* If ScalarsOfChart returns unevaluated, use default 1+1D coordinates *)
  If[Head[coordSyms] === ScalarsOfChart,
    coordSyms = {t[], x[]}
  ];
  coordSyms
];

(* === Covariant Derivative to Explicit Derivative Conversion === *)
(* Converts CD[{i, -chart}][f[args]] to Derivative[...][f][args] form *)
(* Handles both covariant (-chart) and contravariant (chart) indices *)
(* Supports 1+1D (indices 0,1) and 2+1D (indices 0,1,2) *)
(* Uses FixedPoint to handle nested derivative expressions *)

ConvertCDToDerivatives[expr_, chart_] := Module[
  {replaceCDfunc, isCDlike},

  (* Function to check if something is a CD-like operator *)
  (* Uses string matching to be context-independent *)
  isCDlike[x_] := StringMatchQ[ToString[Head[x]], "*CD*"];

  replaceCDfunc[e_] := e /. {
    (* === 1+1D rules (2-argument Derivative) === *)
    (* Covariant index: {i, -chart} *)
    f_[{0, -chart}][g_Symbol[args__]] /; isCDlike[f[{0, -chart}]] :>
      Derivative[1, 0][g][args],
    f_[{1, -chart}][g_Symbol[args__]] /; isCDlike[f[{1, -chart}]] :>
      Derivative[0, 1][g][args],
    f_[{0, -chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{0, -chart}]] :>
      Derivative[n + 1, m][g][args],
    f_[{1, -chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{1, -chart}]] :>
      Derivative[n, m + 1][g][args],
    (* Contravariant index: {i, chart} *)
    f_[{0, chart}][g_Symbol[args__]] /; isCDlike[f[{0, chart}]] :>
      Derivative[1, 0][g][args],
    f_[{1, chart}][g_Symbol[args__]] /; isCDlike[f[{1, chart}]] :>
      Derivative[0, 1][g][args],
    f_[{0, chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{0, chart}]] :>
      Derivative[n + 1, m][g][args],
    f_[{1, chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{1, chart}]] :>
      Derivative[n, m + 1][g][args],

    (* === 2+1D rules (3-argument Derivative for index 2 = y) === *)
    (* These rules convert 2-arg Derivative to 3-arg when index 2 is encountered *)
    (* Covariant index 2: {2, -chart} *)
    f_[{2, -chart}][g_Symbol[args__]] /; isCDlike[f[{2, -chart}]] :>
      Derivative[0, 0, 1][g][args],
    f_[{2, -chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{2, -chart}]] :>
      Derivative[n, m, 1][g][args],
    f_[{2, -chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{2, -chart}]] :>
      Derivative[n, m, p + 1][g][args],
    (* Contravariant index 2: {2, chart} *)
    f_[{2, chart}][g_Symbol[args__]] /; isCDlike[f[{2, chart}]] :>
      Derivative[0, 0, 1][g][args],
    f_[{2, chart}][Derivative[n_, m_][g_][args__]] /; isCDlike[f[{2, chart}]] :>
      Derivative[n, m, 1][g][args],
    f_[{2, chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{2, chart}]] :>
      Derivative[n, m, p + 1][g][args],
    (* Also handle 3-arg derivatives with indices 0 and 1 *)
    f_[{0, -chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{0, -chart}]] :>
      Derivative[n + 1, m, p][g][args],
    f_[{1, -chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{1, -chart}]] :>
      Derivative[n, m + 1, p][g][args],
    f_[{0, chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{0, chart}]] :>
      Derivative[n + 1, m, p][g][args],
    f_[{1, chart}][Derivative[n_, m_, p_][g_][args__]] /; isCDlike[f[{1, chart}]] :>
      Derivative[n, m + 1, p][g][args]
  };

  (* Apply repeatedly until fixed point (max 10 iterations) *)
  FixedPoint[replaceCDfunc, expr, 10]
];

(* === Field Head Extraction === *)
(* Extracts tensor head from applied form like phi[] or A[-a] *)

ExtractFieldHead[field_] := If[Head[field] === Symbol, field, Head[field]];

(* === Numeric Coefficient Extraction === *)
(* Extracts numeric coefficient from a term containing the named field *)
(* Unified version combining logic from ComponentDecompose and ExportJSON *)

ExtractNumericCoefficient[term_, fieldName_] := Module[
  {coeff},

  (* Replace field and its derivatives with 1 to extract coefficient *)
  coeff = term /. {
    (* Match Derivative[...][f][args] form (applied derivatives) *)
    Derivative[__][f_][__] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,
    (* Match f[args] form (bare field with numeric suffix like A0, A1) *)
    f_[__] /; StringMatchQ[ToString[f], ToString[fieldName] ~~ DigitCharacter ...] :> 1,
    (* Direct field match *)
    fieldName -> 1,
    _Derivative[__][fieldName] -> 1
  };

  (* Simplify *)
  coeff = Simplify[coeff];

  (* Handle various coefficient forms *)
  Which[
    NumericQ[coeff], coeff,
    (* Negative symbolic coefficient: -m2 -> return -1.0 with warning *)
    MatchQ[coeff, Times[-1, _Symbol]],
      Print["Warning: Symbolic coefficient '", coeff, "' converted to -1.0. ",
            "For proper numeric values, substitute symbols before JSON export."];
      -1.0,
    (* Positive symbolic coefficient: m2 -> return 1.0 with warning *)
    MatchQ[coeff, _Symbol],
      Print["Warning: Symbolic coefficient '", coeff, "' converted to 1.0. ",
            "For proper numeric values, substitute symbols before JSON export."];
      1.0,
    (* Try numeric evaluation as last resort *)
    NumericQ[Quiet[N[coeff]]], Quiet[N[coeff]],
    (* Default to 1.0 if all else fails, with warning *)
    True,
      Print["Warning: Could not extract numeric coefficient from '", coeff,
            "'. Defaulting to 1.0."];
      1.0
  ]
];

(* === Levi-Civita (Epsilon) Tensor Evaluation === *)
(* Evaluates epsilon tensor components to numeric ±1 values *)
(* For Minkowski signature (-,+,+,...), the sign conventions are: *)
(*   - Covariant ε_012... = -√|g| = -1 (in flat space, √|g|=1 for Minkowski) *)
(*   - Contravariant ε^012... = +1/√|g| = +1 *)
(*   - Mixed indices: each raised index contributes metric factor η^ii *)

(* Levi-Civita symbol: +1 for even permutations, -1 for odd, 0 for repeated *)
LeviCivitaValue[indices_List] := Module[{n = Length[indices]},
  If[Length[DeleteDuplicates[indices]] != n,
    0,  (* Repeated index → 0 *)
    Signature[indices]  (* +1 even perm, -1 odd perm *)
  ]
];

(* Helper: Compute metric factor for raising/lowering indices in Minkowski space *)
(* For signature (-,+,+,...): η^00 = -1, η^11 = η^22 = ... = +1 *)
MinkowskiMetricFactor[idx_Integer] := If[idx == 0, -1, 1];

(* Evaluate epsilon tensor components in expression *)
(* Handles epsilon tensors created by xAct's DefMetric (e.g., epsiloneta, epsiloneta3) *)
(* General approach: identify all epsilon patterns and compute based on index positions *)
EvaluateEpsilonComponents[expr_, chart_] := Module[
  {rules, isEpsilon, evaluateEpsilon3, evaluateEpsilon2},

  (* Pattern to check if a symbol name contains "epsilon" *)
  isEpsilon[s_] := StringMatchQ[ToString[s], "*epsilon*", IgnoreCase -> True];

  (*
     For 3D epsilon with mixed indices:
     Start from fully covariant ε_ijk = -Signature[{i,j,k}] (Minkowski convention)
     Each raised index multiplies by η^ii (the metric component)

     Example: ε_{i j}^{k} = ε_{i j m} η^{m k}
       For i=0, j=1, k=0: ε_{0 1 m} η^{m 0} = ε_{0 1 0} η^{0 0} = 0 * (-1) = 0
       For i=0, j=1, k=2: ε_{0 1 m} η^{m 2} = ε_{0 1 2} η^{2 2} = (-1) * (+1) = -1

     But when xAct puts indices, the third slot being chart (not -chart) means raised.
     So {i, -chart}, {j, -chart}, {k, chart} means ε_{i j}^k
  *)

  evaluateEpsilon3[i_, j_, k_, isUp1_, isUp2_, isUp3_] := Module[
    {baseValue, metricFactor},
    (* Base: fully covariant Levi-Civita for Minkowski = -Signature *)
    baseValue = -LeviCivitaValue[{i, j, k}];
    (* Metric factors for raised indices *)
    metricFactor = 1;
    If[isUp1, metricFactor *= MinkowskiMetricFactor[i]];
    If[isUp2, metricFactor *= MinkowskiMetricFactor[j]];
    If[isUp3, metricFactor *= MinkowskiMetricFactor[k]];
    baseValue * metricFactor
  ];

  evaluateEpsilon2[i_, j_, isUp1_, isUp2_] := Module[
    {baseValue, metricFactor},
    baseValue = -LeviCivitaValue[{i, j}];
    metricFactor = 1;
    If[isUp1, metricFactor *= MinkowskiMetricFactor[i]];
    If[isUp2, metricFactor *= MinkowskiMetricFactor[j]];
    baseValue * metricFactor
  ];

  rules = {
    (* 3D epsilon with any combination of up/down indices *)
    (* Match pattern: f[{i, ±chart}, {j, ±chart}, {k, ±chart}] where f is epsilon-like *)
    f_Symbol[{i_Integer, s1_}, {j_Integer, s2_}, {k_Integer, s3_}] /;
      isEpsilon[f] && MemberQ[{chart, -chart}, s1] && MemberQ[{chart, -chart}, s2] && MemberQ[{chart, -chart}, s3] :>
      evaluateEpsilon3[i, j, k, s1 === chart, s2 === chart, s3 === chart],

    (* 2D epsilon with any combination of up/down indices *)
    f_Symbol[{i_Integer, s1_}, {j_Integer, s2_}] /;
      isEpsilon[f] && MemberQ[{chart, -chart}, s1] && MemberQ[{chart, -chart}, s2] :>
      evaluateEpsilon2[i, j, s1 === chart, s2 === chart]
  };

  (* Apply rules repeatedly until no more matches *)
  expr //. rules
];

End[];
EndPackage[];
