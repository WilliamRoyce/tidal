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
(* Uses FixedPoint to handle nested derivative expressions *)

ConvertCDToDerivatives[expr_, chart_] := Module[
  {replaceCDfunc, isCDlike},

  (* Function to check if something is a CD-like operator *)
  (* Uses string matching to be context-independent *)
  isCDlike[x_] := StringMatchQ[ToString[Head[x]], "*CD*"];

  replaceCDfunc[e_] := e /. {
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
      Derivative[n, m + 1][g][args]
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
    (* Negative symbolic coefficient: -m2 -> return -1.0 *)
    MatchQ[coeff, Times[-1, _Symbol]], -1.0,
    (* Positive symbolic coefficient: m2 -> return 1.0 *)
    MatchQ[coeff, _Symbol], 1.0,
    (* Try numeric evaluation as last resort *)
    NumericQ[Quiet[N[coeff]]], Quiet[N[coeff]],
    (* Default to 1.0 if all else fails *)
    True, 1.0
  ]
];

End[];
EndPackage[];
