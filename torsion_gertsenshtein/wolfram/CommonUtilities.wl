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

GetCoordinateSymbols::fallback =
  "Chart `1` returned invalid coordinates. Using default 1+1D {t[], x[]}.";

GetChartDimension::usage =
  "GetChartDimension[chart] returns the spacetime dimension (number of coordinates) \
from the chart. Use this to infer dimension dynamically instead of hardcoding.";

GetChartDimension::fallback =
  "Chart `1` dimension could not be determined. Defaulting to 2 (1+1D).";

ValidateDimension::usage =
  "ValidateDimension[dim] checks if dimension is within supported range (≤4). \
Returns True if valid, False otherwise.";

ValidateDimension::unsupported =
  "Dimension `1` exceeds supported maximum of `2`. CD conversion rules may be incomplete.";

$MaxSupportedDimension::usage =
  "$MaxSupportedDimension is the maximum spacetime dimension (currently 4 for 3+1D) \
supported by the CD conversion rules.";

GenerateCDRules::usage =
  "GenerateCDRules[dim, chart] generates CD to Derivative conversion rules for \
dim-dimensional spacetime. Eliminates hardcoded rule duplication and enables easy \
extension to higher dimensions.";

ConvertCDToDerivatives::usage =
  "ConvertCDToDerivatives[expr, chart] converts covariant derivative operators \
(CD) to explicit Mathematica Derivative form.";

ConvertCDToDerivatives::incomplete =
  "CD conversion did not fully converge after `1` iterations. Remaining CD operators may exist.";

ExtractFieldHead::usage =
  "ExtractFieldHead[field] extracts the tensor head from an applied field \
like phi[] or A[-a].";

ExtractNumericCoefficient::usage =
  "ExtractNumericCoefficient[term, fieldName] extracts the numeric coefficient \
from a term containing the named field.";

ExtractNumericCoefficient::symbolic =
  "Symbolic coefficient `1` found for field `2`. Storing as metadata.";

ExtractCoefficientWithSymbolic::usage =
  "ExtractCoefficientWithSymbolic[term, fieldName] extracts coefficient preserving \
both numeric value and symbolic expression. Returns <|\"numeric\" -> N, \"symbolic\" -> S|>.";

EvaluateEpsilonComponents::usage =
  "EvaluateEpsilonComponents[expr, chart] evaluates Levi-Civita (epsilon) tensor \
components to numeric values ±1 for Minkowski signature (-,+,+,...).";

LeviCivitaValue::usage =
  "LeviCivitaValue[indices] returns the Levi-Civita symbol value: +1 for even \
permutations, -1 for odd permutations, 0 for repeated indices.";

Begin["`Private`"];

(* === xAct Introspection Helpers === *)
(* Use xAct's type-checking functions instead of fragile string matching *)
(* Wrapped in Quiet to handle inputs that xAct doesn't expect *)

(* Check if expr is a Christoffel symbol using xAct introspection *)
IsChristoffelSymbol[expr_] := Quiet[
  TrueQ[xAct`xTensor`ChristoffelQ[expr]],
  {xAct`xTensor`ChristoffelQ::argx, General::stop}
];

(* Check if expr is a covariant derivative operator using xAct introspection *)
IsCovDOperator[expr_] := Quiet[
  TrueQ[xAct`xTensor`CovDQ[Head[expr]]],
  {xAct`xTensor`CovDQ::argx, General::stop}
];

(* Check if expr is an epsilon (Levi-Civita) tensor using xAct introspection *)
(* Falls back to string matching for cases xAct doesn't recognize *)
IsEpsilonTensor[expr_] := Quiet[
  TrueQ[xAct`xTensor`EpsilonQ[expr]] ||
  TrueQ[xAct`xTensor`EpsilonQ[Head[expr]]],
  {xAct`xTensor`EpsilonQ::argx, General::stop}
];

(* === Christoffel Symbol Removal === *)
(* In flat Minkowski space, Christoffel symbols vanish *)
(* Uses xAct introspection for reliable detection *)

RemoveChristoffelSymbols[expr_] :=
  expr /. f_[__] /; IsChristoffelSymbol[f] -> 0;

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
  (* If ScalarsOfChart returns unevaluated, warn user and use default 1+1D coordinates *)
  If[Head[coordSyms] === ScalarsOfChart,
    Message[GetCoordinateSymbols::fallback, chart];
    coordSyms = {t[], x[]}
  ];
  coordSyms
];

(* === Dimension Validation === *)
(* Validates that dimension is within supported range for CD conversion rules *)
$MaxSupportedDimension = 4;  (* Currently 3+1D is the maximum *)

ValidateDimension[dim_Integer] := Module[{},
  If[dim > $MaxSupportedDimension,
    Message[ValidateDimension::unsupported, dim, $MaxSupportedDimension];
    False,
    True
  ]
];

(* === Dimension Inference from Chart === *)
(* Returns the number of coordinates (spacetime dimension) from the chart *)
(* This is the single source of truth for dimension throughout the pipeline *)

GetChartDimension[chart_] := Module[{coords, dim},
  coords = ScalarsOfChart[chart];
  If[Head[coords] === ScalarsOfChart,
    Message[GetChartDimension::fallback, chart];
    2,  (* Default to 1+1D if chart not properly defined *)
    dim = Length[coords];
    ValidateDimension[dim];
    dim
  ]
];

(* === Dynamic CD Rule Generation === *)
(* Generates CD→Derivative rules dynamically for any dimension up to $MaxSupportedDimension *)
(* This eliminates the 80+ hardcoded rules and enables easy extension to higher dimensions *)

GenerateCDRules[dim_Integer, chart_] := Module[
  {rules, isCDlike},

  (* Validate dimension *)
  If[!ValidateDimension[dim], Return[{}]];

  (* CD-like operator check using xAct introspection *)
  (* Falls back to string matching for edge cases xAct doesn't handle *)
  isCDlike[x_] := IsCovDOperator[x] || StringMatchQ[ToString[Head[x]], "*CD*"];

  (* Generate rules for all combinations:
     - Each coordinate index: 0 to dim-1
     - Both chart signs: chart and -chart (contravariant and covariant)
     - Each Derivative arity: 2 to dim (handles promotion from lower to higher arity) *)

  rules = Flatten[Table[
    With[{idx = i, chartSign = s, arity = a},
      {
        (* Rule for applying CD to bare symbol: f[{idx, chartSign*chart}][g_Symbol[args]] *)
        f_[{idx, chartSign*chart}][g_Symbol[args__]] /; isCDlike[f[{idx, chartSign*chart}]] :>
          With[{newOrders = ReplacePart[ConstantArray[0, arity], idx + 1 -> 1]},
            Derivative[Sequence @@ newOrders][g][args]
          ],

        (* Rule for applying CD to existing Derivative of MATCHING arity *)
        f_[{idx, chartSign*chart}][Derivative[orders__][g_][args__]] /;
          isCDlike[f[{idx, chartSign*chart}]] && Length[{orders}] == arity :>
          With[{newOrders = ReplacePart[{orders}, idx + 1 -> {orders}[[idx + 1]] + 1]},
            Derivative[Sequence @@ newOrders][g][args]
          ],

        (* Rule for PROMOTING lower-arity Derivative to this arity (for indices beyond current arity) *)
        If[idx >= arity - 1 && arity > 2,
          f_[{idx, chartSign*chart}][Derivative[orders__][g_][args__]] /;
            isCDlike[f[{idx, chartSign*chart}]] && Length[{orders}] < arity :>
            With[{paddedOrders = PadRight[{orders}, arity, 0]},
              With[{newOrders = ReplacePart[paddedOrders, idx + 1 -> paddedOrders[[idx + 1]] + 1]},
                Derivative[Sequence @@ newOrders][g][args]
              ]
            ],
          Nothing
        ]
      }
    ],
    {i, 0, dim - 1},      (* Coordinate indices *)
    {s, {1, -1}},         (* Chart signs: contravariant (1) and covariant (-1) *)
    {a, 2, dim}           (* Derivative arities *)
  ]];

  (* Remove Nothing entries from conditional rules *)
  DeleteCases[rules, Nothing]
];

(* === Covariant Derivative to Explicit Derivative Conversion === *)
(* Converts CD[{i, -chart}][f[args]] to Derivative[...][f][args] form *)
(* Now uses dynamically generated rules via GenerateCDRules *)
(* Supports any dimension up to $MaxSupportedDimension *)

ConvertCDToDerivatives[expr_, chart_] := Module[
  {dim, rules, result, maxIter = 20, isCDlike},

  (* Get dimension from chart *)
  dim = GetChartDimension[chart];

  (* Generate rules dynamically for this dimension *)
  rules = GenerateCDRules[dim, chart];

  (* CD-like operator check for convergence validation using xAct introspection *)
  (* Falls back to string matching for edge cases xAct doesn't handle *)
  isCDlike[x_] := IsCovDOperator[x] || StringMatchQ[ToString[Head[x]], "*CD*"];

  (* Apply rules repeatedly until fixed point *)
  result = FixedPoint[# /. rules &, expr, maxIter];

  (* Warn if CD-like operators remain after conversion *)
  If[!FreeQ[result, f_ /; isCDlike[f]],
    Message[ConvertCDToDerivatives::incomplete, maxIter]
  ];

  result
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

(* === Symbolic Coefficient Extraction === *)
(* Extracts coefficient preserving both numeric value and symbolic expression *)
(* Returns an Association with "numeric" and "symbolic" keys *)

ExtractCoefficientWithSymbolic[term_, fieldName_] := Module[
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

  (* Return Association with numeric value and symbolic expression *)
  Which[
    NumericQ[coeff],
      <|"numeric" -> N[coeff], "symbolic" -> None|>,
    (* Negative symbolic coefficient: -m2 *)
    MatchQ[coeff, Times[-1, _Symbol]],
      Message[ExtractNumericCoefficient::symbolic, coeff, fieldName];
      <|"numeric" -> -1.0, "symbolic" -> ToString[coeff]|>,
    (* Positive symbolic coefficient: m2 *)
    MatchQ[coeff, _Symbol],
      Message[ExtractNumericCoefficient::symbolic, coeff, fieldName];
      <|"numeric" -> 1.0, "symbolic" -> ToString[coeff]|>,
    (* Try numeric evaluation *)
    NumericQ[Quiet[N[coeff]]],
      <|"numeric" -> Quiet[N[coeff]], "symbolic" -> None|>,
    (* Complex symbolic expression *)
    True,
      Message[ExtractNumericCoefficient::symbolic, coeff, fieldName];
      <|"numeric" -> 1.0, "symbolic" -> ToString[coeff]|>
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
  {rules, isEpsilon, evaluateEpsilon4, evaluateEpsilon3, evaluateEpsilon2},

  (* Check if symbol is an epsilon tensor using xAct introspection *)
  (* Falls back to string matching for edge cases xAct doesn't handle *)
  isEpsilon[s_] := IsEpsilonTensor[s] || StringMatchQ[ToString[s], "*epsilon*", IgnoreCase -> True];

  (*
     For epsilon tensors with mixed indices:
     Start from fully covariant ε_{ijk...} = -Signature[{i,j,k,...}] (Minkowski convention)
     Each raised index multiplies by η^ii (the metric component)

     Example: ε_{i j}^{k} = ε_{i j m} η^{m k}
       For i=0, j=1, k=0: ε_{0 1 m} η^{m 0} = ε_{0 1 0} η^{0 0} = 0 * (-1) = 0
       For i=0, j=1, k=2: ε_{0 1 m} η^{m 2} = ε_{0 1 2} η^{2 2} = (-1) * (+1) = -1

     When xAct puts indices, chart (not -chart) means raised index.
     So {i, -chart}, {j, -chart}, {k, chart} means ε_{i j}^k
  *)

  (* 4D epsilon tensor evaluation (for 3+1D spacetime) *)
  evaluateEpsilon4[i_, j_, k_, l_, isUp1_, isUp2_, isUp3_, isUp4_] := Module[
    {baseValue, metricFactor},
    (* Base: fully covariant Levi-Civita for Minkowski = -Signature *)
    baseValue = -LeviCivitaValue[{i, j, k, l}];
    (* Metric factors for raised indices *)
    metricFactor = 1;
    If[isUp1, metricFactor *= MinkowskiMetricFactor[i]];
    If[isUp2, metricFactor *= MinkowskiMetricFactor[j]];
    If[isUp3, metricFactor *= MinkowskiMetricFactor[k]];
    If[isUp4, metricFactor *= MinkowskiMetricFactor[l]];
    baseValue * metricFactor
  ];

  (* 3D epsilon tensor evaluation (for 2+1D spacetime) *)
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

  (* 2D epsilon tensor evaluation (for 1+1D spacetime) *)
  evaluateEpsilon2[i_, j_, isUp1_, isUp2_] := Module[
    {baseValue, metricFactor},
    baseValue = -LeviCivitaValue[{i, j}];
    metricFactor = 1;
    If[isUp1, metricFactor *= MinkowskiMetricFactor[i]];
    If[isUp2, metricFactor *= MinkowskiMetricFactor[j]];
    baseValue * metricFactor
  ];

  rules = {
    (* 4D epsilon with any combination of up/down indices (for 3+1D spacetime) *)
    f_Symbol[{i_Integer, s1_}, {j_Integer, s2_}, {k_Integer, s3_}, {l_Integer, s4_}] /;
      isEpsilon[f] &&
      MemberQ[{chart, -chart}, s1] && MemberQ[{chart, -chart}, s2] &&
      MemberQ[{chart, -chart}, s3] && MemberQ[{chart, -chart}, s4] :>
      evaluateEpsilon4[i, j, k, l, s1 === chart, s2 === chart, s3 === chart, s4 === chart],

    (* 3D epsilon with any combination of up/down indices (for 2+1D spacetime) *)
    f_Symbol[{i_Integer, s1_}, {j_Integer, s2_}, {k_Integer, s3_}] /;
      isEpsilon[f] &&
      MemberQ[{chart, -chart}, s1] && MemberQ[{chart, -chart}, s2] && MemberQ[{chart, -chart}, s3] :>
      evaluateEpsilon3[i, j, k, s1 === chart, s2 === chart, s3 === chart],

    (* 2D epsilon with any combination of up/down indices (for 1+1D spacetime) *)
    f_Symbol[{i_Integer, s1_}, {j_Integer, s2_}] /;
      isEpsilon[f] && MemberQ[{chart, -chart}, s1] && MemberQ[{chart, -chart}, s2] :>
      evaluateEpsilon2[i, j, s1 === chart, s2 === chart]
  };

  (* Apply rules repeatedly until no more matches *)
  expr //. rules
];

End[];
EndPackage[];
