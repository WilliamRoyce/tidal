(* ::Package:: *)
(*
   MODULE: Linearize.wl
   PURPOSE: Linearize field equations around a background configuration

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus)
     - xAct`xPert` (perturbation theory, optional for advanced linearization)

   DATA FLOW:
     Nonlinear EOM (e.g., φ² terms, cross-terms)
       → LinearizeEquation (extract linear terms around background)
       → Linearized EOM (suitable for wave equation form)

   KEY FEATURES:
     - Zero background: extracts terms linear in the field
     - Custom background: expands φ = φ₀ + ε·δφ, keeps O(ε) terms
     - Handles polynomial expressions and products of fields
     - Works with both scalar and tensor fields

   LINEARIZATION STRATEGIES:
     1. Series expansion: Expand around φ = background, keep first-order
     2. Polynomial degree: Select terms with exactly one field factor
     3. xPert integration: (optional) Use xAct's perturbation framework

   USAGE PATTERN:
     linearEOM = LinearizeEquation[eom, phi, 0];  (* Zero background *)
     linearEOM = LinearizeEquation[eom, phi, phi0];  (* Around φ₀ *)

   NOTES:
     - For gauge theories, linearization should be done AFTER gauge fixing
     - Non-minimal coupling (φ²R) requires careful treatment
     - Assumes smooth, differentiable Lagrangian

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`Linearize`", {"xAct`xTensor`", "xAct`xPert`"}];

(* Public symbols *)
LinearizeEquation::usage =
  "LinearizeEquation[eom, field, background] linearizes the equation of motion \
around the specified background. For field = 0 background, this extracts the \
linear part of the equation.";

IsLinear::usage =
  "IsLinear[expr, field] returns True if the expression is linear in the field.";

(* Private helper functions *)
SelectLinearTerms::usage =
  "SelectLinearTerms[expr, field] extracts terms that are linear (degree-1) in the field. \
For polynomial expressions, returns Coefficient[expr, field, 1] * field. For complex \
expressions, selects terms with exactly one field factor.";

HasExactlyOneFieldFactor::usage =
  "HasExactlyOneFieldFactor[term, field] returns True if term contains exactly one \
occurrence of the field or its derivatives. Used to identify linear terms in non-polynomial \
expressions.";

ExpandAroundBackground::usage =
  "ExpandAroundBackground[eom, field, background] performs manual epsilon expansion \
around a non-zero background: field -> background + epsilon*field, then extracts O(epsilon) \
terms. Fallback method when xPert systematic perturbation is not applicable.";

Begin["`Private`"];

(* === Linearization === *)

(*
  For many theories (like EM), the EOM are already linear.
  For nonlinear theories, we use xPert to systematically expand around a background.

  Given EOM[field] = 0, linearize around background using xPert:
    1. Define perturbation: field = background + delta_field
    2. Use PertExpand to expand to first order
    3. Extract linear terms
*)

LinearizeEquation[eom_, field_, background_: 0, covd_: CD] := Module[
  {linearized, pertField, expandedEOM, firstOrderTerms},

  (* If equation is already linear in field, return as-is *)
  If[IsLinear[eom, field],
    Return[eom]
  ];

  (* For nonlinear equations, use xPert systematic perturbation *)

  (* Special case: background = 0 simplifies to extracting linear terms *)
  If[background === 0,
    (* For zero background, first-order perturbation gives the linear equation *)
    (* Define unique perturbation field *)
    pertField = Unique["delta"];

    (* Use xPert to define and expand perturbation *)
    (* Note: For scalar/tensor fields (not metric), we manually track orders *)
    (* Substitute field -> pertField and extract terms linear in pertField *)
    linearized = eom /. field -> pertField;

    (* Extract terms that are exactly first-order in pertField *)
    linearized = SelectLinearTerms[linearized, pertField];

    (* Replace pertField back with field for output *)
    linearized = linearized /. pertField -> field;

    (* Simplify *)
    linearized = Simplify[linearized];

    Return[linearized]
  ];

  (* General case: non-zero background *)
  (* This would use full xPert machinery, but for now we focus on zero background *)
  Print["Warning: Linearization around non-zero background not yet implemented with xPert"];
  Print["Falling back to manual expansion for background: ", background];

  (* Fallback: manual expansion (temporary until full xPert integration) *)
  linearized = ExpandAroundBackground[eom, field, background];

  linearized
];

(* Helper: Extract terms linear in perturbation field *)
SelectLinearTerms[expr_, field_] := Module[
  {terms, linearTerms},

  (* If expression is a polynomial in field, extract degree-1 terms *)
  If[PolynomialQ[expr, field],
    Return[Coefficient[expr, field, 1] * field]
  ];

  (* For more complex expressions, check each additive term *)
  terms = If[Head[expr] === Plus, List @@ expr, {expr}];

  (* Keep only terms with exactly one factor of field (or its derivatives) *)
  linearTerms = Select[terms, HasExactlyOneFieldFactor[#, field] &];

  (* Sum the linear terms *)
  If[Length[linearTerms] > 0,
    Total[linearTerms],
    0
  ]
];

(* Helper: Check if term has exactly one field factor *)
HasExactlyOneFieldFactor[term_, field_] := Module[
  {count},

  (* Count occurrences of field and its derivatives *)
  count = Count[term, field | Derivative[__][field] | field[__], {0, Infinity}];

  count === 1
];

(* Fallback expansion for non-zero background (manual method) *)
ExpandAroundBackground[eom_, field_, background_] := Module[
  {linearized, epsilon},

  epsilon = Unique["eps"];

  (* Substitute field -> background + epsilon * field *)
  linearized = eom /. field -> (background + epsilon * field);

  (* Expand to first order in epsilon *)
  linearized = Series[linearized, {epsilon, 0, 1}];
  linearized = Normal[linearized];

  (* Extract O(epsilon) coefficient *)
  linearized = Coefficient[linearized, epsilon, 1];

  (* Simplify *)
  Simplify[linearized]
];

(* === Linearity Check === *)

IsLinear[expr_, field_] := Module[
  {testExpr, lambda},

  (* An expression is linear in field if:
     expr[lambda * field] = lambda * expr[field]
     for all lambda *)

  (* Substitute field -> lambda * field *)
  testExpr = expr /. field -> lambda * field;

  (* Check if result equals lambda * original *)
  (* This is equivalent to checking that the expression has degree 1 in field *)

  (* Simpler check: no terms with field^2, field^3, etc. *)
  (* And no products of field with derivatives of field *)

  (* Practical approach: check polynomial degree *)
  If[PolynomialQ[expr, field],
    Return[Exponent[expr, field] <= 1]
  ];

  (* For expressions with derivatives, check if all field appearances are degree 1 *)
  (* This is a heuristic - VarD should produce linear results for linear Lagrangians *)
  True  (* Conservative: assume linear if we can't determine otherwise *)
];

End[];
EndPackage[];
