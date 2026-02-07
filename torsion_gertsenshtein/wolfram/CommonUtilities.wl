(* ::Package:: *)
(*
   MODULE: CommonUtilities.wl
   PURPOSE: Shared utilities for the Lagrangian-to-PDE pipeline

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus, covariant derivatives)
     - xAct`xCoba` (coordinate bases, component extraction)

   KEY RESPONSIBILITIES:
     - CD → Derivative conversion: ConvertCDToDerivatives transforms xAct covariant
       derivatives into explicit Mathematica Derivative[n,m,...][f][t,x,...] form
     - Christoffel evaluation: EvaluateChristoffelComponents evaluates connection terms
       to their proper values (0 for flat Minkowski, non-zero for curved spacetimes)
     - Epsilon tensor evaluation: EvaluateEpsilonComponents converts Levi-Civita
       tensors to numeric ±1 values (supports 2D, 3D, 4D)
     - Metric evaluation: EvaluateMinkowskiMetric handles (-,+,+,...) signature
     - Coefficient extraction: ExtractCoefficientWithSymbolic preserves symbolic info
     - xAct introspection: IsChristoffelSymbol, IsCovDOperator, IsEpsilonTensor

   DIMENSION SUPPORT:
     - 1+1D: 2 coordinates (t, x)
     - 2+1D: 3 coordinates (t, x, y)
     - 3+1D: 4 coordinates (t, x, y, z)
     - $MaxSupportedDimension = 4 (extensible with GenerateCDRules)

   SIGN CONVENTIONS (Minkowski):
     - Metric signature: (-1, +1, +1, ...)
     - Covariant epsilon: ε_{012...} = -1 (time index first)
     - Contravariant epsilon: ε^{012...} = +1

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`CommonUtilities`",
  {"xAct`xTensor`", "xAct`xCoba`"}];

(* Public symbols *)
EvaluateChristoffelComponents::usage =
  "EvaluateChristoffelComponents[expr, chart] evaluates Christoffel symbol components \
to their numeric values for the given chart. For flat Minkowski space with constant \
metric components, all Christoffel symbols evaluate to zero. For curved spacetimes, \
this would give the proper non-zero connection coefficients. \
Use EvaluateChristoffelComponents[expr, chart, True] to compute from xAct's metric \
definition (for non-flat metrics).";

ComputeChristoffelFromMetric::usage =
  "ComputeChristoffelFromMetric[chart, covd] computes Christoffel symbol components \
using xAct's TraceBasisDummy. Returns a 3D array where christoffel[[a+1, b+1, c+1]] = Gamma^a_bc.";

ComputeChristoffelFromMetricMatrix::usage =
  "ComputeChristoffelFromMetricMatrix[chart, metricMatrix] computes Christoffel symbols \
directly from the metric matrix using the standard formula: \
Gamma^a_{bc} = (1/2) g^{ad} (d_b g_{dc} + d_c g_{bd} - d_d g_{bc}). \
Returns a 3D array where christoffel[[a+1, b+1, c+1]] = Gamma^a_bc. \
Works for any metric including time-dependent ones.";

RemoveChristoffelSymbols::usage =
  "RemoveChristoffelSymbols[expr] sets all Christoffel symbol terms to zero \
(valid for flat Minkowski space). Legacy function; prefer EvaluateChristoffelComponents.";

EvaluateMinkowskiMetric::usage =
  "EvaluateMinkowskiMetric[expr, chart] evaluates metric components for \
Minkowski signature (-1, +1, ...) in the given chart.";

EvaluateMetricComponents::usage =
  "EvaluateMetricComponents[expr, chart, metricMatrix] evaluates metric \
components from an explicit metric matrix. Computes the inverse metric and \
substitutes both covariant and contravariant components. Used for curved \
spacetime metrics (e.g., conformal g_ab = Omega^2 eta_ab).";

EvaluatePDMetric::usage =
  "EvaluatePDMetric[expr, chart, metricMatrix] evaluates partial derivatives \
of metric components by substituting the known metric values and differentiating. \
For example, PD_0[g_{11}] -> D[metricMatrix[[2,2]], coords[[1]]].";

ExpandChristoffelsToMetricDerivatives::usage =
  "ExpandChristoffelsToMetricDerivatives[eom, covd, chart] uses xAct's ChangeCovD \
and ChristoffelToGradMetric to replace covariant derivatives and Christoffel \
symbols with partial derivatives and metric components. This is the general \
approach for curved spacetimes where Christoffel symbols are non-trivial.";

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

(* xAct introspection helpers *)
IsChristoffelSymbol::usage =
  "IsChristoffelSymbol[expr] returns True if expr is a Christoffel symbol, using \
xAct's ChristoffelQ. Wrapped in Quiet for robustness with non-xAct inputs.";

IsCovDOperator::usage =
  "IsCovDOperator[expr] returns True if expr is a covariant derivative operator, \
using xAct's CovDQ. Falls back to string matching for non-xAct expressions.";

IsEpsilonTensor::usage =
  "IsEpsilonTensor[expr] returns True if expr is a Levi-Civita (epsilon) tensor, \
using xAct's EpsilonQ. Handles both tensor heads and applied tensors.";

IsCurvatureTensor::usage =
  "IsCurvatureTensor[f] returns True if f is a curvature tensor (Riemann, Ricci, \
RicciScalar, Weyl, Schouten, Cotton, or Kretschner). Uses string matching on the \
xAct-generated tensor name. Christoffel symbols are excluded.";

EvaluateCurvatureComponents::usage =
  "EvaluateCurvatureComponents[expr, chart, metricMatrix:None] evaluates background \
curvature tensor components from the metric definition. For a constant metric \
(all components independent of coordinates), all curvature vanishes by the theorem \
that d_i g_{ab} = 0 implies Gamma = 0 implies R_{abcd} = 0. The function COMPUTES \
this from the metric properties rather than assuming it. For non-constant metrics, \
curvature is left unevaluated (future: compute from Christoffel formula). \
If metricMatrix is None, uses flat Minkowski metric diag(-1,+1,...).";

MinkowskiMetricFactor::usage =
  "MinkowskiMetricFactor[idx] returns the metric factor for index raising/lowering in \
Minkowski space: -1 for time (idx=0), +1 for spatial indices (idx>0).";

Begin["`Private`"];

(* === xAct Introspection Helpers === *)
(* Use xAct's type-checking functions instead of fragile string matching *)
(* Wrapped in Quiet to handle inputs that xAct doesn't expect *)

(* Check if expr is a Christoffel symbol using xAct introspection *)
(* Also uses string matching fallback for xCoba-generated symbols like ChristoffelCDPDcart *)
(* NOTE: xAct's ChristoffelQ often returns False for ChristoffelCDPDchart symbols,
   which are created by xCoba to represent the connection between covariant derivatives.
   The string matching fallback ensures these are properly detected. *)
IsChristoffelSymbol[expr_] := Quiet[
  Module[{xActResult, stringCheck},
    (* First try xAct's ChristoffelQ *)
    xActResult = TrueQ[xAct`xTensor`ChristoffelQ[expr]];
    If[xActResult, Return[True]];

    (* Fallback: check if symbol name contains "Christoffel" *)
    (* This catches ChristoffelCDPDcart and similar xCoba-generated symbols *)
    stringCheck = StringContainsQ[ToString[expr], "Christoffel"];
    stringCheck
  ],
  {xAct`xTensor`ChristoffelQ::argx, General::stop}
];

(* Check if expr is a covariant derivative operator using xAct introspection *)
(* For expr like CD[-a][phi[]], Head[expr] = CD[-a], and Head[Head[expr]] = CD *)
(* We need to extract the base symbol, not the applied form with indices *)
IsCovDOperator[expr_] := Quiet[
  Module[{h = Head[expr], baseSymbol},
    (* If head is already a symbol, check directly *)
    If[Head[h] === Symbol,
      TrueQ[xAct`xTensor`CovDQ[h]],
      (* Otherwise extract the base symbol (e.g., CD from CD[-a]) *)
      baseSymbol = Head[h];
      TrueQ[xAct`xTensor`CovDQ[baseSymbol]]
    ]
  ],
  {xAct`xTensor`CovDQ::argx, General::stop}
];

(* Check if expr is an epsilon (Levi-Civita) tensor using xAct introspection *)
(* Falls back to string matching for cases xAct doesn't recognize *)
IsEpsilonTensor[expr_] := Quiet[
  TrueQ[xAct`xTensor`EpsilonQ[expr]] ||
  TrueQ[xAct`xTensor`EpsilonQ[Head[expr]]],
  {xAct`xTensor`EpsilonQ::argx, General::stop}
];

(* === Christoffel Symbol Evaluation === *)
(* For flat Minkowski space with constant metric, Christoffel symbols vanish naturally.
   We evaluate them properly rather than removing them blindly, so the code will work
   correctly for curved spacetimes where Christoffel symbols are non-zero.

   For Minkowski metric (constant components), the Christoffel formula gives:
     Γ^a_bc = (1/2) g^{ad} (∂_b g_{cd} + ∂_c g_{bd} - ∂_d g_{bc}) = 0
   because all metric derivatives vanish for constant metric.

   Implementation: We evaluate Christoffel components in basis indices to their
   proper numeric values. For flat space this is 0; for curved space it would be
   computed from the metric.

   NOTE: The legacy RemoveChristoffelSymbols function is kept for backward compatibility
   but now uses the improved IsChristoffelSymbol detection which catches xCoba-generated
   symbols like ChristoffelCDPDcart.
*)

(* Compute Christoffel symbols from xAct's metric definition.
   Uses xAct's ToBasis and TraceBasisDummy to extract the symbolically computed values.
   NO hardcoding - all values derived from the metric via xAct.

   PREREQUISITE: Must have called MetricInBasis[metric, -chart, metricMatrix] first
   to set the metric components. Otherwise returns unevaluated expressions.

   API: xCoba computes Christoffels when metric is defined. To extract components:
   1. ToBasis[chart][expr] - decompose abstract tensor to basis form
   2. TraceBasisDummy - evaluate dummy basis indices to numeric values
*)
ComputeChristoffelFromMetric[chart_, covd_] := Module[
  {dim, christoffelExpr, basisExpr, components},

  dim = GetChartDimension[chart];

  (* Get the Christoffel symbol name for this covariant derivative.
     xCoba generates it as Christoffel<CovD><PD><chart> *)
  christoffelExpr = Symbol["Christoffel" <> ToString[covd] <> "PD" <> ToString[chart]];

  (* Extract component values using proper xAct decomposition *)
  components = Table[
    Module[{expr},
      (* Build the Christoffel with specific indices *)
      expr = christoffelExpr[{a, chart}, {-b, -chart}, {-c, -chart}];
      (* TraceBasisDummy evaluates basis components to their numeric/symbolic values *)
      Quiet[Simplify[TraceBasisDummy[expr]], {TraceBasisDummy::unknown}]
    ],
    {a, 0, dim-1}, {b, 0, dim-1}, {c, 0, dim-1}
  ];

  components
];

(* Compute Christoffel symbols directly from the metric matrix using the standard formula:
     Gamma^a_{bc} = (1/2) g^{ad} (d_b g_{dc} + d_c g_{bd} - d_d g_{bc})
   This is fully general and works for any metric, including time-dependent ones. *)
ComputeChristoffelFromMetricMatrix[chart_, metricMatrix_] := Module[
  {dim, coords, invMetric, dg, components},

  dim = Length[metricMatrix];
  coords = GetCoordinateSymbols[chart];

  (* Compute inverse metric *)
  invMetric = Simplify[Inverse[metricMatrix]];

  (* Compute partial derivatives of metric: dg[[i, j, k]] = d_k g_{ij} *)
  (* Indices are 1-based here *)
  dg = Table[
    Simplify[D[metricMatrix[[i, j]], coords[[k]]]],
    {i, dim}, {j, dim}, {k, dim}
  ];

  (* Gamma^a_{bc} = (1/2) g^{ad} (d_b g_{dc} + d_c g_{bd} - d_d g_{bc}) *)
  (* Note: Table indices are 1-based, so a,b,c,d range from 1 to dim *)
  components = Table[
    Simplify[Sum[
      1/2 * invMetric[[aa, dd]] * (dg[[dd, cc, bb]] + dg[[bb, dd, cc]] - dg[[bb, cc, dd]]),
      {dd, dim}
    ]],
    {aa, dim}, {bb, dim}, {cc, dim}
  ];

  components
];

(* Evaluate Christoffel components to their numeric/symbolic values *)
(* Third argument controls behavior:
     None or False -> flat space, all Christoffels = 0
     True -> compute from metric (requires metricMatrix in 4th argument)
     covdSymbol -> use that CovD symbol to compute via xAct's TraceBasisDummy *)
(* Fourth argument (optional): explicit metric matrix for direct computation *)
EvaluateChristoffelComponents[expr_, chart_, covd_:None, metricMatrix_:None] := Module[
  {result, christoffelValues},

  If[covd === None || covd === False,
    (* Flat space: all Christoffels = 0 *)
    (* Use broad pattern to catch Christoffels with unresolved dummy indices.
       TraceBasisDummy doesn't always trace dummies inside products of
       Christoffel symbols (they're not proper tensors), leaving abstract
       indices like j$12345 instead of {Integer, ±chart}. Since ALL
       Christoffels vanish in flat space, match any applied Christoffel. *)
    result = expr /. {
      f_[__] /; IsChristoffelSymbol[f] :> 0
    };
    Return[result]
  ];

  (* Compute Christoffel values *)
  If[covd === True,
    (* Compute directly from the metric matrix using the standard formula *)
    If[metricMatrix === None,
      Throw["EvaluateChristoffelComponents: covd=True requires a MetricMatrix. \
Pass the metric matrix as the 4th argument."]
    ];
    christoffelValues = ComputeChristoffelFromMetricMatrix[chart, metricMatrix],
    (* Explicit CovD symbol provided - use xAct's TraceBasisDummy *)
    christoffelValues = ComputeChristoffelFromMetric[chart, covd]
  ];

  (* Substitute computed values into expression *)
  (* Note: xCoba uses negative indices for covariant components, e.g., {-2, -chart} means index 2 covariant.
     Use Abs[] to extract the actual index value for array lookup. *)
  result = expr /. {
    f_[{a_Integer, s1:(chart | -chart)}, {b_Integer, s2:(chart | -chart)}, {c_Integer, s3:(chart | -chart)}] /;
      IsChristoffelSymbol[f] :> christoffelValues[[Abs[a]+1, Abs[b]+1, Abs[c]+1]]
  };

  Simplify[result]
];

(* Legacy function for backward compatibility - now properly detects xCoba symbols *)
RemoveChristoffelSymbols[expr_] :=
  expr /. f_[__] /; IsChristoffelSymbol[f] -> 0;

(* === Minkowski Metric Evaluation === *)
(* Evaluates metric components for signature (-1, +1, +1, ...) *)
(* Handles both covariant indices {i, chart} and contravariant {i, -chart} *)

EvaluateMinkowskiMetric[expr_, chart_] := Module[
  {dim, minkowskiMatrix},
  (* Build Minkowski metric dynamically based on chart dimension *)
  dim = GetChartDimension[chart];
  minkowskiMatrix = DiagonalMatrix[Join[{-1}, ConstantArray[1, dim - 1]]];
  (* Delegate to the general metric component evaluator *)
  EvaluateMetricComponents[expr, chart, minkowskiMatrix]
];

(* === General Metric Component Evaluation === *)
(* Evaluates metric components from an explicit metric matrix. *)
(* Unlike EvaluateMinkowskiMetric which hardcodes flat (-1,+1,+1,...), *)
(* this computes the inverse metric from the given matrix and substitutes *)
(* both covariant and contravariant components. *)
(* Used for curved spacetime metrics (e.g., conformal g_ab = Omega^2 eta_ab). *)

EvaluateMetricComponents[expr_, chart_, metricMatrix_] := Module[
  {dim, invMatrix, rules},
  dim = Length[metricMatrix];
  invMatrix = Simplify[Inverse[metricMatrix]];
  rules = Flatten[{
    (* Contravariant g^{ij}: positive chart basis *)
    Table[
      _Symbol[{i, chart}, {j, chart}] -> Simplify[invMatrix[[i + 1, j + 1]]],
      {i, 0, dim - 1}, {j, 0, dim - 1}
    ],
    (* Covariant g_{ij}: negative chart basis *)
    Table[
      _Symbol[{i, -chart}, {j, -chart}] -> Simplify[metricMatrix[[i + 1, j + 1]]],
      {i, 0, dim - 1}, {j, 0, dim - 1}
    ]
  }];
  expr /. rules
];

(* === Partial Derivative of Metric Evaluation === *)
(* Evaluates PD_i[g_{jk}] by substituting known metric components and differentiating *)
(* This is needed after ChristoffelToGradMetric expands Christoffels into metric derivatives *)

EvaluatePDMetric[expr_, chart_, metricMatrix_] := Module[
  {dim, coords},
  dim = Length[metricMatrix];
  coords = GetCoordinateSymbols[chart];

  (* Evaluate PD of covariant metric: PD_i[g_{jk}] -> D[g_{jk}, coords[[i+1]]] *)
  (* The pattern matches PDchart[{i, -chart}][metric[{j, sign*chart}, {k, sign*chart}]] *)
  expr /. {
    pd_[{i_Integer, -chart}][g_[{j_Integer, s1:(chart | -chart)}, {k_Integer, s2:(chart | -chart)}]] /;
      IsCovDOperator[pd] || StringMatchQ[ToString[pd], "PD*"] :>
      Simplify[D[metricMatrix[[Abs[j] + 1, Abs[k] + 1]], coords[[i + 1]]]]
  }
];

(* === Christoffel Expansion via xAct === *)
(* Uses xAct's ChangeCovD and ChristoffelToGradMetric to replace covariant derivatives
   and Christoffel symbols with partial derivatives and metric components.
   This is the general approach for curved spacetimes. *)

ExpandChristoffelsToMetricDerivatives[eom_, covd_, chart_] := Module[
  {step1, step2, pdChart},
  (* Get the PD operator for this chart: xCoba names it PD<chartName> *)
  pdChart = Symbol["PD" <> ToString[chart]];

  (* ChangeCovD: CovD -> PD + Christoffel *)
  step1 = ChangeCovD[eom, covd, pdChart];

  (* ChristoffelToGradMetric: Christoffel -> metric derivatives *)
  step2 = ChristoffelToGradMetric[step1];

  step2
];

(* === Coordinate Symbol Extraction === *)
(* Gets coordinate scalars from chart, with fallback for unevaluated cases *)

GetCoordinateSymbols[chart_] := GetCoordinateSymbols[chart] = Module[{coordSyms},
  coordSyms = ScalarsOfChart[chart];
  (* If ScalarsOfChart returns unevaluated, fail with clear error - no silent fallback *)
  If[Head[coordSyms] === ScalarsOfChart,
    Throw[StringJoin[
      "GetCoordinateSymbols: Cannot extract coordinates from chart '", ToString[chart], "'. ",
      "Ensure chart is properly defined with DefChart before calling this function."
    ]]
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

GetChartDimension[chart_] := GetChartDimension[chart] = Module[{coords, dim},
  coords = ScalarsOfChart[chart];
  (* If ScalarsOfChart returns unevaluated, fail with clear error - no silent fallback to 2 *)
  If[Head[coords] === ScalarsOfChart,
    Throw[StringJoin[
      "GetChartDimension: Cannot determine dimension from chart '", ToString[chart], "'. ",
      "Ensure chart is properly defined with DefChart before calling this function."
    ]]
  ];
  dim = Length[coords];
  If[!ValidateDimension[dim],
    Throw[StringJoin[
      "GetChartDimension: Dimension ", ToString[dim], " exceeds maximum supported (",
      ToString[$MaxSupportedDimension], ")."
    ]]
  ];
  dim
];

(* === Dynamic CD Rule Generation === *)
(* Generates CD→Derivative rules dynamically for any dimension up to $MaxSupportedDimension *)
(* This eliminates the 80+ hardcoded rules and enables easy extension to higher dimensions *)

GenerateCDRules[dim_Integer, chart_] := GenerateCDRules[dim, chart] = Module[
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
        (* Only generate when idx < arity, otherwise idx+1 exceeds array bounds *)
        If[idx < arity,
          f_[{idx, chartSign*chart}][g_Symbol[args__]] /; isCDlike[f[{idx, chartSign*chart}]] :>
            With[{newOrders = ReplacePart[ConstantArray[0, arity], idx + 1 -> 1]},
              Derivative[Sequence @@ newOrders][g][args]
            ],
          Nothing
        ],

        (* Rule for applying CD to existing Derivative of MATCHING arity *)
        (* Only generate when idx < arity *)
        If[idx < arity,
          f_[{idx, chartSign*chart}][Derivative[orders__][g_][args__]] /;
            isCDlike[f[{idx, chartSign*chart}]] && Length[{orders}] == arity :>
            With[{newOrders = ReplacePart[{orders}, idx + 1 -> {orders}[[idx + 1]] + 1]},
              Derivative[Sequence @@ newOrders][g][args]
            ],
          Nothing
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
  (* Excludes Christoffel symbols (which contain "CD" in their name but are connection *)
  (* coefficients, not covariant derivative operators) *)
  isCDlike[x_] := Module[{h = Head[x], hStr},
    If[IsChristoffelSymbol[h], Return[False]];
    hStr = ToString[h];
    If[StringMatchQ[hStr, "*Christoffel*"], Return[False]];
    IsCovDOperator[x] || StringMatchQ[hStr, "*CD*"]
  ];

  (* Apply rules repeatedly until fixed point *)
  result = FixedPoint[# /. rules &, expr, maxIter];

  (* Fail if CD-like operators remain after conversion *)
  If[!FreeQ[result, f_ /; isCDlike[f]],
    Message[ConvertCDToDerivatives::incomplete, maxIter];
    Throw[StringJoin[
      "ConvertCDToDerivatives: Conversion did not fully converge after ",
      ToString[maxIter], " iterations. Unconverted CD operators remain in expression. ",
      "This may indicate unsupported derivative patterns or missing conversion rules."
    ]]
  ];

  result
];

(* === Field Head Extraction === *)
(* Extracts tensor head from applied form like phi[] or A[-a] *)

ExtractFieldHead[field_] := If[Head[field] === Symbol, field, Head[field]];

(* === Shared Field Replacement Rules === *)
(* Used by both ExtractNumericCoefficient and ExtractCoefficientWithSymbolic *)
(* Replaces all occurrences of the named field (and its derivatives) with 1, *)
(* leaving just the coefficient. *)

fieldReplacementRules[fieldName_] := {
  (* Match Derivative[...][f][args] form (applied derivatives) *)
  Derivative[__][f_][__] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,
  (* Match f[args] form (bare field with numeric suffix like A0, A1) *)
  f_[__] /; StringMatchQ[ToString[f], ToString[fieldName] ~~ DigitCharacter ...] :> 1,
  (* Match f[] form (scalar tensor with no arguments, like phi[]) *)
  f_[] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,
  (* Match f[index] form (vector tensor with one index, like A[-a]) *)
  f_[_] /; StringContainsQ[ToString[f], ToString[fieldName]] :> 1,
  (* Direct field match *)
  fieldName -> 1,
  _Derivative[__][fieldName] -> 1
};

(* === Numeric Coefficient Extraction === *)
(* Extracts numeric coefficient from a term containing the named field *)
(* Unified version combining logic from ComponentDecompose and ExportJSON *)

ExtractNumericCoefficient[term_, fieldName_] := Module[
  {coeff},

  coeff = Simplify[term /. fieldReplacementRules[fieldName]];

  (* Handle various coefficient forms - fail explicitly for unresolved symbolic coefficients *)
  Which[
    NumericQ[coeff], coeff,
    (* Try numeric evaluation *)
    NumericQ[Quiet[N[coeff]]], Quiet[N[coeff]],
    (* Symbolic coefficient that cannot be resolved numerically - FAIL *)
    True,
      Throw[StringJoin[
        "ExtractNumericCoefficient: Cannot extract numeric coefficient from '",
        ToString[coeff], "' for field '", ToString[fieldName], "'. ",
        "Substitute symbolic coefficients before export, e.g.: expr /. {m2 -> 1.0}"
      ]]
  ]
];

(* === Symbolic Coefficient Extraction === *)
(* Extracts coefficient preserving both numeric value and symbolic expression *)
(* Returns an Association with "numeric" and "symbolic" keys *)

ExtractCoefficientWithSymbolic[term_, fieldName_] := Module[
  {coeff},

  coeff = Simplify[term /. fieldReplacementRules[fieldName]];

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
  {rules, isEpsilon, evaluateEpsilonN},

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

  (* Dimension-agnostic epsilon tensor evaluation *)
  (* Works for any dimension: extracts indices and up/down flags from argument list *)
  evaluateEpsilonN[indices_List, isUpFlags_List] := Module[
    {baseValue, metricFactor},
    (* Base: fully covariant Levi-Civita for Minkowski = -Signature *)
    baseValue = -LeviCivitaValue[indices];
    (* Metric factors for raised indices: each raised index i contributes η^{ii} *)
    metricFactor = Product[
      If[isUpFlags[[k]], MinkowskiMetricFactor[indices[[k]]], 1],
      {k, Length[indices]}
    ];
    baseValue * metricFactor
  ];

  (* Single dimension-agnostic rule: matches any epsilon tensor with integer basis indices *)
  rules = {
    f_Symbol[args__] /;
      isEpsilon[f] &&
      MatchQ[{args}, {Repeated[{_Integer, chart | -chart}]}] :>
      evaluateEpsilonN[
        {args}[[All, 1]],  (* Extract integer indices *)
        Map[# === chart &, {args}[[All, 2]]]  (* True if raised (chart), False if lowered (-chart) *)
      ]
  };

  (* Apply rules repeatedly until no more matches *)
  expr //. rules
];

(* === Curvature Tensor Detection and Evaluation === *)
(* Identifies curvature tensors by xAct-generated naming convention *)
(* and evaluates them from the metric definition, NOT hardcoded values *)

(* Check if symbol is a curvature tensor using xAct naming conventions *)
(* xAct generates names like RiemannCD, RicciCD, RicciScalarCD, etc. *)
IsCurvatureTensor[f_] := Module[{name = ToString[f]},
  StringContainsQ[name,
    "Riemann" | "Ricci" | "Weyl" | "Schouten" | "Cotton" | "Kretschner"
  ] && !StringContainsQ[name, "Christoffel"]
];

(* Evaluate background curvature components from the metric definition *)
(* This COMPUTES curvature from metric properties — no hardcoding of values *)
(*
   Mathematical basis:
   - For a metric with constant components (∂_i g_{ab} = 0 for all i,a,b):
     Γ^a_{bc} = (1/2) g^{ad}(∂_b g_{cd} + ∂_c g_{bd} - ∂_d g_{bc}) = 0
     R^a_{bcd} = ∂_c Γ^a_{db} - ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} - Γ^a_{de} Γ^e_{cb} = 0
   - For non-constant metrics, curvature must be explicitly computed.

   This ensures that if the metric is accidentally wrong (non-flat when
   flat was intended), the curvature will NOT be silently zeroed.
*)
EvaluateCurvatureComponents[expr_, chart_, metricMatrix_:None] := Module[
  {result = expr, dim, coords, localMetric, metricIsConstant},

  dim = GetChartDimension[chart];
  coords = GetCoordinateSymbols[chart];

  (* Determine the metric to analyze *)
  localMetric = If[metricMatrix === None,
    DiagonalMatrix[Join[{-1}, ConstantArray[1, dim - 1]]],
    metricMatrix
  ];

  (* Check if all metric components are independent of coordinates *)
  metricIsConstant = And @@ Flatten[Table[
    FreeQ[localMetric[[i, j]], Alternatives @@ coords],
    {i, dim}, {j, dim}
  ]];

  If[metricIsConstant,
    (* Constant metric: ∂g = 0 ⟹ Γ = 0 ⟹ R_{abcd} = 0 (computed, not assumed) *)
    result = result /. {
      f_[__] /; IsCurvatureTensor[f] :> 0,
      f_[] /; IsCurvatureTensor[f] :> 0
    },
    (* Non-constant metric: curvature is non-trivial *)
    (* Leave curvature tensors for the decomposition pipeline to handle *)
    (* They will be evaluated through ToBasis + TraceBasisDummy + Christoffel evaluation *)
    Null
  ];

  result
];

End[];
EndPackage[];
