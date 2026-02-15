(* ::Package:: *)
(*
   MODULE: Linearize.wl
   PURPOSE: Perturbation theory and linearization for gravitational and field theories

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus)
     - xAct`xPert` (systematic perturbation theory)

   DATA FLOW:
     Full nonlinear theory (e.g., Einstein-Hilbert, Maxwell, etc.)
       → SetupMetricPerturbation (define g = g0 + ε h)
       → LinearizeTensorExpression (perturb any tensor expression via xPert)
       → Gauge-unfixed linearized EOM (abstract index form)
       → TTGaugeConditions (informational: describes TT gauge for initial conditions)

   KEY FEATURES:
     - xPert-based metric perturbation: systematic, not manual
     - LinearizeTensorExpression works on ANY tensor expression (Einstein, Ricci, etc.)
     - Gauge-unfixed output as PRIMARY result (preserves constraint structure)
     - TTGaugeConditions: informational tool describing TT gauge conditions
     - Works for any dimension (2+1D, 3+1D, etc.)
     - Works with scalar, vector, and tensor fields

   USAGE PATTERN:
     (* 1. Set up perturbation *)
     hPert = SetupMetricPerturbation[g, h, epsilon];

     (* 2. Linearize any tensor expression *)
     linEinstein = LinearizeTensorExpression[Einstein[CD][-a,-b]];
     linRicci = LinearizeTensorExpression[RicciCD[-a,-b]];

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`Linearize`", {"xAct`xTensor`", "xAct`xPert`"}];

(* === Primary API: Metric Perturbation and Linearization === *)

SetupMetricPerturbation::usage =
  "SetupMetricPerturbation[metric, pertName, paramName] sets up a first-order metric \
perturbation using xPert's DefMetricPerturbation. Returns the perturbation tensor \
pertName[LI[1], -a, -b]. Idempotent: skips definition if already set up. \
Example: hPert = SetupMetricPerturbation[g, h, epsilon]";

LinearizeTensorExpression::usage =
  "LinearizeTensorExpression[tensorExpr] computes the first-order perturbation of \
any tensor expression using xPert. Works with any xAct tensor expression — Einstein \
tensor, Ricci tensor, Riemann tensor, field strength, etc. Returns the linearized \
expression in abstract index form (set = 0 for vacuum EOM). \
Examples: \
  LinearizeTensorExpression[Einstein[CD][-a,-b]]  (* linearized Einstein equations *) \
  LinearizeTensorExpression[RicciCD[-a,-b]]        (* linearized Ricci tensor *)";

(* === Informational Gauge Tools === *)

TTGaugeConditions::usage =
  "TTGaugeConditions[hPert, metric] returns the transverse-traceless gauge conditions \
as an informational Association: trace h = 0 (traceless) and h_{0mu} = 0 (transverse). \
These are constraints on the solution space, not evolution equations. They reduce \
the 10 DOF of h_ab to 2 physical polarizations (+ and x) for gravitational waves. \
TT gauge is imposed via initial conditions in the simulation, not at the Lagrangian level.";

Begin["`Private`"];

(* ================================================================ *)
(* === Metric Perturbation Setup (xPert) === *)
(* ================================================================ *)

SetupMetricPerturbation[metric_, pertName_Symbol, paramName_Symbol] := Module[
  {},

  (* Define perturbation parameter if not already defined *)
  If[!ConstantSymbolQ[paramName],
    DefConstantSymbol[paramName, PrintAs -> "\[Epsilon]"]
  ];

  (* Define metric perturbation using xPert *)
  (* DefMetricPerturbation[metric, perturbation, parameter] *)
  (* This defines pertName[LI[n], -a, -b] for order n perturbations *)
  (* Idempotent: skip if perturbation already defined for this metric *)
  If[!xTensorQ[pertName],
    DefMetricPerturbation[metric, pertName, paramName]
  ];

  (* Return the first-order perturbation tensor *)
  pertName[LI[1], -a, -b]
];

(* ================================================================ *)
(* === Linearize Tensor Expression (xPert) === *)
(* ================================================================ *)

(*
  General-purpose first-order perturbation of any tensor expression.

  For the Einstein tensor specifically:
    G_ab = R_ab - (1/2) g_ab R is the variational derivative
    of the Einstein-Hilbert action: G_ab = (2/sqrt(-g)) * δS_EH/δg^{ab}.
    Perturbation[G[-a,-b], 1] gives the linearized vacuum EOM.

  But this function works on ANY tensor expression — Ricci tensor,
  Riemann tensor, field strength tensors, etc. The perturbation is
  computed systematically via xPert's Perturbation + ExpandPerturbation.

  The result preserves gauge freedom (no gauge has been chosen).
*)

LinearizeTensorExpression[tensorExpr_] := Module[
  {result},

  (* Compute first-order perturbation using xPert *)
  result = Perturbation[tensorExpr, 1];

  (* Expand the perturbation (converts Perturbation[] wrappers to explicit h terms) *)
  result = ExpandPerturbation[result];

  (* Validate that xPert actually expanded something *)
  (* If no metric perturbation was set up, ExpandPerturbation returns the *)
  (* unexpanded Perturbation[] wrapper unchanged *)
  If[!FreeQ[result, Perturbation],
    Throw[StringJoin[
      "LinearizeTensorExpression: Result still contains unexpanded Perturbation[] wrappers. ",
      "Ensure SetupMetricPerturbation was called before linearizing. ",
      "Expression: ", ToString[Short[result, 3]]
    ]]
  ];

  (* Apply canonical simplifications *)
  result = ToCanonical[result];
  result = ContractMetric[result];

  result
];

(* ================================================================ *)
(* === Gauge Conditions (Informational) === *)
(* ================================================================ *)

TTGaugeConditions[hPert_, metric_] := Module[{},
  (* TT gauge imposes two conditions on top of de Donder:
     1. Tracelessness: h = metric^{ab} h_{ab} = 0
     2. Transversality: h_{0μ} = 0 (purely spatial perturbation)

     Together with de Donder gauge, this leaves only 2 physical DOF
     in 3+1D: the + and × gravitational wave polarizations.

     These are constraints on the solution/initial data, not
     dynamical equations. They can be imposed by:
     a) Choosing initial data satisfying these conditions
     b) Projecting the solution onto the TT subspace *)

  <|
    "gauge" -> "transverse_traceless",
    "tracelessness" -> "metric^{ab} h_{ab} = 0",
    "transversality" -> "h_{0mu} = 0 for all mu",
    "physical_dof" -> 2,
    "polarizations" -> {"h_+", "h_x"}
  |>
];

End[];
EndPackage[];
