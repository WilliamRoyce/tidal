(* ::Package:: *)
(*
   MODULE: GaugeFix.wl
   PURPOSE: Extensible gauge-fixing terms for the TIDAL pipeline

   DEPENDENCIES:
     - xAct`xTensor` (tensor algebra)

   DATA FLOW:
     Lagrangian + gauge specification -> Lagrangian + gauge-fixing term

   KEY FEATURES:
     - Core primitive: AddGaugeFixingTerm (adds arbitrary term to L)
     - Preset builders: each returns an EXPRESSION, not an action
     - Extensible: new gauge = new Build*GaugeTerm function

   USAGE PATTERN:
     gaugeTerm = BuildLorenzGaugeTerm[A, eta, CD, 1];
     L = AddGaugeFixingTerm[L, gaugeTerm];

   ADDING NEW PRESETS:
     1. Write BuildNewGaugeTerm[field_, metric_, covd_, params___] in this file
     2. Add entry to _GAUGE_PRESETS in tidal/cli/_derive.py
     3. Done.

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`GaugeFix`", {"xAct`xTensor`"}];

(* === Core Primitive === *)

AddGaugeFixingTerm::usage =
  "AddGaugeFixingTerm[lagrangian, gaugeTerm] returns lagrangian + gaugeTerm. \
This is the fundamental operation: any gauge-fixing term is simply added to L \
before Euler-Lagrange derivation.";

(* === Preset Builders === *)

BuildLorenzGaugeTerm::usage =
  "BuildLorenzGaugeTerm[field, metric, covd, xi] returns the Lorenz gauge-fixing \
term -(1/(2*xi))*(metric^{ab} covd_a field_b)^2 for a vector field. \
With xi=1 (Feynman gauge), this reduces Maxwell equations to uncoupled \
wave equations.  Example: BuildLorenzGaugeTerm[A, eta, CD, 1]";

Begin["`Private`"];

(* ================================================================ *)
(* === Core Primitive                                            === *)
(* ================================================================ *)

AddGaugeFixingTerm[lagrangian_, gaugeTerm_] := lagrangian + gaugeTerm;

(* ================================================================ *)
(* === Lorenz Gauge: -(1/2xi)(div A)^2                          === *)
(* ================================================================ *)

BuildLorenzGaugeTerm[field_, metric_, covd_, xi_:1] := Module[
  {divA, gaugeTerm},

  (* Validate xi > 0 *)
  If[!NumericQ[xi] || xi <= 0,
    Throw["BuildLorenzGaugeTerm: xi must be a positive number, got " <> ToString[xi]]
  ];

  (* div A = metric^{ab} covd_a field_b *)
  divA = metric[a, b] covd[-a][field[-b]];

  (* Gauge-fixing term: -(1/2xi)(div A)^2 *)
  gaugeTerm = -(1/(2 xi)) divA^2;

  gaugeTerm
];

End[];
EndPackage[];
