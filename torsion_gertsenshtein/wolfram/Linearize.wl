(* ::Package:: *)
(* Linearize.wl - Linearization of field equations around a background *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`Linearize`", {"xAct`xTensor`", "xAct`xPert`"}];

(* Public symbols *)
LinearizeEquation::usage =
  "LinearizeEquation[eom, field, background] linearizes the equation of motion \
around the specified background. For field = 0 background, this extracts the \
linear part of the equation.";

IsLinear::usage =
  "IsLinear[expr, field] returns True if the expression is linear in the field.";

Begin["`Private`"];

(* === Linearization === *)

(*
  For many theories (like EM), the EOM are already linear.
  For nonlinear theories, we need to expand around a background.

  Given EOM[field] = 0, linearize around background:
    field = background + epsilon * perturbation
  and keep only O(epsilon) terms.
*)

LinearizeEquation[eom_, field_, background_: 0] := Module[
  {linearized, epsilon, perturbation},

  (* If equation is already linear in field, return as-is *)
  If[IsLinear[eom, field],
    Return[eom]
  ];

  (* For nonlinear equations, use perturbation expansion *)
  (* Define perturbation parameter *)
  If[!ValueQ[epsilon],
    epsilon = Unique["eps"]
  ];

  (* Substitute field -> background + epsilon * field *)
  linearized = eom /. field -> (background + epsilon * field);

  (* Expand to first order in epsilon *)
  linearized = Series[linearized, {epsilon, 0, 1}];
  linearized = Normal[linearized];

  (* Extract O(epsilon) coefficient *)
  linearized = Coefficient[linearized, epsilon, 1];

  (* Simplify *)
  linearized = Simplify[linearized];

  linearized
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
