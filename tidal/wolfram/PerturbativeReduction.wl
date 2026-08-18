(* ::Package:: *)
(*
   MODULE: PerturbativeReduction.wl
   PURPOSE: Lagrangian Perturbative Substitution (LPS) — substitute higher
            time derivatives in `lagComp` using the order-0 EOM, before the
            Hamiltonian is extracted via Legendre transform.

   PROBLEM ADDRESSED:
     For theories with `b5*R~^2` or similar higher-curvature corrections, the
     Lagrangian carries `(d_t^2 q)^2` terms.  Phase B's IBP-on-time pass
     reduces these as far as possible but is fundamentally unable to lower
     `q-double-dot * X` (the irreducible 2nd-order residue) without the EOM.
     The resulting `canonical.hamiltonian_terms` JSON contains `mixed_2_*`
     operators encoding `q-double-dot` directly, which:

       (1) require finite-differencing the velocity slot for snapshot-based
           measurement (unstable for widely-spaced snapshots), and
       (2) violate the conventional shape of a quadratic Hamiltonian for a
           genuinely 2nd-order theory (which should depend only on (q, q_t)).

   ALGORITHM (Lagrangian Perturbative Substitution, LPS):

     1. ExtractBaseEquations: set every declared small_parameter to zero in
        each E-L equation, giving the order-0 EOM `M_0 q-double-dot = K_0(q, q_t)`
        for every dynamical (time_order >= 2) field.

     2. SolveBaseEqsForSecondTimeDerivatives: solve the order-0 EOMs
        simultaneously for `q-double-dot_i` per dynamical field, producing rules
        `D[q_i[t,x,...], {t,2}] -> expr_i(q, q_t, q_x, ...)`.  Constraint
        fields (time_order=0) are skipped automatically.

     3. BuildHigherDerivativeSubstitutions: chain the rules upward.  For
        order-3, differentiate the order-2 RHS w.r.t. `t` and substitute
        lower-order rules; for order-4, differentiate again; etc., up to
        the maxOrder detected in lagComp.

     4. ReduceLagrangian: apply the substitution rules to `lagComp` via
        `//.`, then Series-truncate over each smallParam in deterministic
        (sorted) order to drop higher-order tails introduced by the
        substitution.  Assert the result has time_order <= 1.

   API CONVENTION:
     Inputs use Wolfram function HEADS for field identification (matching
     TIDAL's `compToFunc[name]` lookup).  The caller resolves string field
     names to heads via `compToFunc` before calling LPS.

   RELATIONSHIP TO v6 PARKER-SIMON:
     The v6 method (in tidal/solver/perturbative_driver.py) operates on the
     *equations* — Pass 0 evolves the order-0 EOM, Pass 1 integrates order-1
     corrections via Duhamel.  v6 never touches the Lagrangian.  LPS is the
     complementary Lagrangian-side step that gives Phase D (the Legendre
     transform) a healthy 2nd-order input, eliminating residual `mixed_2_*`
     operators in the resulting Hamiltonian terms.

   NON-PERTURBATIVE THEORIES:
     ReduceLagrangian is a no-op when `smallParams` is empty (returns input
     unchanged) and also when the input lagComp is already 2nd-order or
     lower in time.

   Part of the TIDAL Lagrangian-to-PDE pipeline.
*)

BeginPackage["TorsionGertsenshtein`PerturbativeReduction`"];

(* Public symbols *)
DetectMaxTimeOrder::usage =
  "DetectMaxTimeOrder[expr, fieldHeads] returns the maximum time-derivative \
order appearing in `expr` for any of the given dynamical-field heads.  \
Returns 0 when no time-derivative subexpression is found.  Used by \
ReduceLagrangian to decide whether substitution is necessary and how many \
higher-order rules to chain.";

ExtractBaseEquations::usage =
  "ExtractBaseEquations[fieldEquationsByHead, smallParams] returns the \
order-0 component E-L equations by setting every smallParam (string symbol \
name) to zero.  Input is a list of {fieldHead, expression} pairs; output \
is the same shape with smallParams substituted.";

SolveBaseEqsForSecondTimeDerivatives::usage =
  "SolveBaseEqsForSecondTimeDerivatives[baseEqsByHead, dynamicalFieldHeads] \
solves the order-0 equations for the second time derivative of each \
dynamical field.  Returns a list of replacement rules \
`Derivative[2, 0, ..., 0][fieldHead][args] -> rhs(q, q_t, ...)`.  Skips \
fields not in `dynamicalFieldHeads`.  Raises an error if the linear \
system is singular (typically a sign that the order-0 theory has gauge \
symmetry).";

BuildHigherDerivativeSubstitutions::usage =
  "BuildHigherDerivativeSubstitutions[fieldEquationsByHead, \
dynamicalFieldHeads, smallParams, maxOrder] returns a list of substitution \
rules covering Derivative[n,0,...] for n = 2, 3, ..., maxOrder for every \
dynamical field.  Higher-order rules are obtained by differentiating \
lower-order RHSs w.r.t. t and applying lower-order rules to eliminate \
intermediate q-double-dot, q-triple-dot factors.";

ReduceLagrangian::usage =
  "ReduceLagrangian[lagComp, fieldEquations, compToFunc, smallParams, \
reductionOrder] applies higher-order time-derivative substitutions to \
`lagComp` and Series-truncates at `reductionOrder` over each smallParam.  \
Returns the reduced Lagrangian with no time derivative of order > 1.  \
No-op when `smallParams` is empty or when `lagComp` is already 1st-order \
in time.  `fieldEquations` is the existing TIDAL Wolfram variable (list \
of {nameString, expression} pairs); `compToFunc` is the existing TIDAL \
Wolfram Association mapping string names to field function heads.";

Begin["`Private`"];

(* === DetectMaxTimeOrder ============================================== *)

DetectMaxTimeOrder[expr_, fieldHeads_List] := Module[
  {derivs, orders},
  derivs = Cases[expr,
    Derivative[ord__][h_][__] /; MemberQ[fieldHeads, h] :> {ord},
    {0, Infinity}];
  If[derivs === {}, Return[0]];
  orders = (#[[1]]) & /@ derivs;
  Max[orders]
];

(* === ExtractBaseEquations =========================================== *)

ExtractBaseEquations[fieldEquationsByHead_List, smallParams_List] := Module[
  {zeroRule},
  zeroRule = Thread[ToExpression /@ smallParams -> 0];
  Table[
    {fieldEquationsByHead[[k, 1]],
     Expand[fieldEquationsByHead[[k, 2]] /. zeroRule]},
    {k, Length[fieldEquationsByHead]}
  ]
];

(* === SolveBaseEqsForSecondTimeDerivatives =========================== *)
(* Returns rules of the form Derivative[2,0,...][h][args] -> rhs.        *)
(* Internally, only fields whose base EOM contains a d_t^2 term are     *)
(* used as solve targets — this skips constraint fields automatically.   *)

SolveBaseEqsForSecondTimeDerivatives[
  baseEqsByHead_List, dynamicalFieldHeads_List
] := Module[
  {targetEqs, ddqVars, sol, fieldHead, eqExpr, args, ddqVarPattern, hasDDQ},
  targetEqs = Select[baseEqsByHead, MemberQ[dynamicalFieldHeads, #[[1]]] &];
  If[targetEqs === {},
    Return[{}]
  ];

  (* Build the d_t^2 variable for each dynamical field by inspecting    *)
  (* the equation expression to find the field's argument signature    *)
  (* AND verifying the equation actually contains d_t^2 of that field. *)
  ddqVars = {};
  Do[
    fieldHead = targetEqs[[k, 1]];
    eqExpr = targetEqs[[k, 2]];
    (* Find ANY occurrence of fieldHead[args__] in the EOM — gives the *)
    (* canonical argument tuple, which is shared across all derivative *)
    (* orders of the same field.                                       *)
    args = FirstCase[eqExpr, fieldHead[a__] :> {a}, Null, {0, Infinity}];
    If[args === Null,
      args = FirstCase[eqExpr,
        Derivative[__][fieldHead][a__] :> {a}, Null, {0, Infinity}]
    ];
    If[args === Null, Continue[]];
    ddqVarPattern = Derivative[
      2, Sequence @@ ConstantArray[0, Length[args] - 1]
    ][fieldHead] @@ args;
    (* Verify the equation actually contains this d_t^2 term — skip    *)
    (* fields whose order-0 EOM is purely algebraic (constraint).       *)
    hasDDQ = !FreeQ[eqExpr, ddqVarPattern];
    If[hasDDQ, AppendTo[ddqVars, ddqVarPattern]],
    {k, Length[targetEqs]}
  ];

  If[ddqVars === {},
    Return[{}]
  ];

  (* Solve the system of `targetEqs[[k,2]] == 0` for the ddq variables. *)
  (* Filter targetEqs to those whose ddq is in ddqVars.                  *)
  sol = Quiet[
    Solve[
      Thread[(#[[2]] & /@ Select[targetEqs,
        Module[{h = #[[1]], a, dvar},
          a = FirstCase[#[[2]], h[aa__] :> {aa}, Null, {0, Infinity}];
          If[a === Null,
            a = FirstCase[#[[2]],
              Derivative[__][h][aa__] :> {aa}, Null, {0, Infinity}]
          ];
          If[a === Null, Return[False, Module]];
          dvar = Derivative[
            2, Sequence @@ ConstantArray[0, Length[a] - 1]
          ][h] @@ a;
          MemberQ[ddqVars, dvar]
        ] &
      ]) == 0],
      ddqVars
    ],
    {Solve::svars}
  ];

  If[sol === {} || sol === {{}},
    Throw[StringJoin[
      "SolveBaseEqsForSecondTimeDerivatives: order-0 EOM system is ",
      "singular for fields ", ToString[dynamicalFieldHeads], ". ",
      "This typically indicates gauge symmetry at the leading order; ",
      "add a TOML [gauge] block to fix the gauge before perturbative ",
      "reduction."
    ]]
  ];

  sol[[1]]
];

(* === BuildHigherDerivativeSubstitutions ============================ *)

BuildHigherDerivativeSubstitutions[
  fieldEquationsByHead_List, dynamicalFieldHeads_List,
  smallParams_List, maxOrder_Integer
] := Module[
  {baseEqs, baseRules, allRules, fieldHead, baseRule, args, rhs, n, k,
   newLhs, newRhs, lhsHeadOp},

  If[maxOrder < 2, Return[{}]];

  baseEqs = ExtractBaseEquations[fieldEquationsByHead, smallParams];
  baseRules = SolveBaseEqsForSecondTimeDerivatives[
    baseEqs, dynamicalFieldHeads];
  allRules = baseRules;

  Do[
    fieldHead = dynamicalFieldHeads[[k]];
    baseRule = FirstCase[allRules,
      r:(Derivative[2, ___][h_][___] -> _) /; h === fieldHead :> r,
      Null];
    If[baseRule === Null, Continue[]];
    args = List @@ baseRule[[1]];

    rhs = baseRule[[2]];
    Do[
      newRhs = D[rhs, args[[1]]] //. allRules;
      lhsHeadOp = Derivative[
        n, Sequence @@ ConstantArray[0, Length[args] - 1]
      ][fieldHead];
      newLhs = lhsHeadOp @@ args;
      AppendTo[allRules, newLhs -> newRhs];
      rhs = newRhs,
      {n, 3, maxOrder}
    ],
    {k, Length[dynamicalFieldHeads]}
  ];

  allRules
];

(* === ReduceLagrangian =============================================== *)

(* === ReduceLagrangian =============================================== *)
(*                                                                       *)
(* SCOPE: This implementation handles the standard "perturbative         *)
(* correction to a fully dynamical theory" case (e.g., Pais-Uhlenbeck    *)
(* perturbations).  It does NOT yet handle the constraint-promotion     *)
(* case where the small parameter promotes an algebraic constraint to a *)
(* dynamical (higher-derivative) field — that case requires solving the *)
(* order-0 algebraic constraint and chain-differentiating the result,   *)
(* which is significantly more complex and has performance issues for   *)
(* multi-field gauge theories (the //. + Series combination on a        *)
(* substituted multi-field lagComp can take >5 minutes on 38-field      *)
(* theories like torsion_gertsenshtein).                                 *)
(*                                                                       *)
(* When the constraint-promotion case is detected (any field with       *)
(* d_t^2 in lagComp but no d_t^2 in its base EOM), LPS skips with a     *)
(* clear diagnostic message.  The Hamiltonian for those theories        *)
(* retains its existing mixed_2_* operators — same behavior as          *)
(* pre-LPS.  See follow-up issue for the constraint-promotion           *)
(* extension.                                                            *)

ReduceLagrangian[
  lagComp_, fieldEquations_List, compToFunc_, smallParams_List,
  reductionOrder_Integer
] := Module[
  {fieldEquationsByHead, allFieldHeads, fieldsWithDtSquaredInLag,
   baseEqs, dynamicalAtOrder0, constraintAtOrder0,
   maxOrder, rules, substituted, truncated, sortedParams,
   residualMaxOrder, hasDDQInBase},

  If[Length[smallParams] === 0,
    Return[lagComp]
  ];

  fieldEquationsByHead = Table[
    {compToFunc[fieldEquations[[k, 1]]], fieldEquations[[k, 2]]},
    {k, Length[fieldEquations]}
  ];
  allFieldHeads = DeleteDuplicates[fieldEquationsByHead[[All, 1]]];

  maxOrder = DetectMaxTimeOrder[lagComp, allFieldHeads];
  If[maxOrder <= 1,
    Print["[LPS] lagComp already <= 1st order in time; no-op"];
    Return[lagComp]
  ];

  (* Identify fields that appear with d_t^2 (or higher) in lagComp.     *)
  fieldsWithDtSquaredInLag = DeleteDuplicates[
    Cases[lagComp,
      Derivative[n_, ___][h_][__] /;
        MemberQ[allFieldHeads, h] && n >= 2 :> h,
      {0, Infinity}]
  ];

  baseEqs = ExtractBaseEquations[fieldEquationsByHead, smallParams];
  hasDDQInBase[fh_] := Module[
    {entry, args, ddqVar},
    entry = FirstCase[baseEqs, {fh, e_} :> e, Null];
    If[entry === Null, Return[False]];
    args = FirstCase[entry, fh[a__] :> {a}, Null, {0, Infinity}];
    If[args === Null,
      args = FirstCase[entry,
        Derivative[__][fh][a__] :> {a}, Null, {0, Infinity}]
    ];
    If[args === Null, Return[False]];
    ddqVar = Derivative[
      2, Sequence @@ ConstantArray[0, Length[args] - 1]
    ][fh] @@ args;
    !FreeQ[entry, ddqVar]
  ];
  dynamicalAtOrder0 = Select[fieldsWithDtSquaredInLag, hasDDQInBase];
  constraintAtOrder0 = Complement[
    fieldsWithDtSquaredInLag, dynamicalAtOrder0];

  Print["[LPS] Reducing lagComp from time-order ", maxOrder,
    " using small_parameters = ", smallParams,
    " (reductionOrder = ", reductionOrder, ")"];
  Print["[LPS] Fields with d_t^2 in lagComp: ",
    Length[fieldsWithDtSquaredInLag],
    " (", Length[dynamicalAtOrder0], " dynamical at order 0, ",
    Length[constraintAtOrder0], " constraint-promoted)"];

  (* Constraint-promotion case: fail loud.  The standard substitution-   *)
  (* based LPS does not apply here — a constraint-at-order-0 field f_c   *)
  (* has no d_t^2 in its base EOM, so we cannot build a substitution    *)
  (* rule for d_t^2 f_c directly.  The full architectural fix requires  *)
  (* solving the order-0 algebraic constraint for f_c, chain-           *)
  (* differentiating the result twice in time, and substituting back     *)
  (* into lagComp.  This was attempted but produced significant         *)
  (* performance regression on multi-field gauge theories (>5min on a   *)
  (* 38-field torsion theory) and may also introduce subtle correctness *)
  (* issues from chain-rule expansion + Series truncation interaction.  *)
  (* Raise the issue clearly so it can be investigated, rather than     *)
  (* producing a partial or silently-degraded result.                    *)
  If[Length[constraintAtOrder0] > 0,
    Throw[StringJoin[
      "ReduceLagrangian: constraint-promotion case detected.  ",
      ToString[Length[constraintAtOrder0]],
      " field(s) appear with d_t^2 in lagComp but are algebraic ",
      "constraints at order 0: ", ToString[constraintAtOrder0], ".  ",
      "The current LPS algorithm handles only the case where every ",
      "field with higher-order time derivatives in lagComp is also ",
      "dynamical (has d_t^2 term) at order 0.  For theories where the ",
      "small parameter promotes an algebraic constraint to a dynamical ",
      "field (e.g., b5*R~^2 in PGT torsion theories: h_4, h_7, h_9 are ",
      "constraints at b5=0 and dynamical at b5!=0), a more sophisticated ",
      "algorithm is needed: solve the order-0 algebraic constraint for ",
      "the field, chain-differentiate the result, and substitute back. ",
      "Investigation needed — see PerturbativeReduction.wl for the ",
      "skipped implementation path."
    ]]
  ];

  rules = BuildHigherDerivativeSubstitutions[
    fieldEquationsByHead, dynamicalAtOrder0, smallParams, maxOrder];
  Print["[LPS] Built ", Length[rules], " substitution rules"];

  substituted = lagComp //. rules;

  sortedParams = Sort[ToExpression /@ smallParams];
  truncated = Fold[
    Normal[Series[#1, {#2, 0, reductionOrder}]] &,
    substituted,
    sortedParams
  ];

  residualMaxOrder = DetectMaxTimeOrder[truncated, allFieldHeads];
  If[residualMaxOrder > 1,
    Throw[StringJoin[
      "ReduceLagrangian: post-substitution lagComp still has ",
      "time-order ", ToString[residualMaxOrder], ".  This is a bug ",
      "in BuildHigherDerivativeSubstitutions (a higher-order rule was ",
      "missed) or an indication that the substitution introduced new ",
      "higher-order terms not covered by the rule chain.  Investigate ",
      "before trusting the resulting Hamiltonian."
    ]]
  ];
  Print["[LPS] Reduction complete; lagComp now <= 1st order in time"];
  truncated
];

End[];
EndPackage[];
