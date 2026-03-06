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
(* These functions use bare abstract indices (a, b, c, ...) that must resolve *)
(* to xAct indices defined by DefManifold in the caller's context.            *)
(* They MUST be defined outside Begin["`Private`"] to avoid Wolfram's package *)
(* scoping from contextualizing them to GaugeFix`Private`a, etc.             *)

BuildLorenzGaugeTerm::usage =
  "BuildLorenzGaugeTerm[field, metric, covd, xi] returns the Lorenz gauge-fixing \
term -(1/(2*xi))*(metric^{ab} covd_a field_b)^2 for a vector field. \
With xi=1 (Feynman gauge), this reduces Maxwell equations to uncoupled \
wave equations.  Example: BuildLorenzGaugeTerm[A, eta, CD, 1]";

BuildDeDonderGaugeTerm::usage =
  "BuildDeDonderGaugeTerm[field, metric, covd, xi] returns the de Donder \
gauge-fixing term -(1/(2*xi))*metric^{be}*(D_b)(D_e) where \
D_b = metric^{ac} covd_a field_{cb} - (1/2) covd_b (metric^{cd} field_{cd}). \
For linearized gravity with xi=1, this reduces Einstein equations to \
uncoupled wave equations.  Example: BuildDeDonderGaugeTerm[h, eta, CD, 1]";

(* ================================================================ *)
(* === Lorenz Gauge: -(1/2xi)(div A)^2                          === *)
(* ================================================================ *)

BuildLorenzGaugeTerm[field_, metric_, covd_, xi_:1] := Module[
  {divA, gaugeTerm},

  (* Validate xi > 0 *)
  If[!NumericQ[xi] || xi <= 0,
    Throw["BuildLorenzGaugeTerm: xi must be a positive number, got " <> ToString[xi]]
  ];

  (* div A = metric^{ab} covd_a field_b, contracted before squaring *)
  (* so that (div A)^2 keeps the contraction inside each Scalar[] copy *)
  divA = metric[a, b] covd[-a][field[-b]];
  divA = ContractMetric[divA, metric];

  (* Gauge-fixing term: -(1/2xi)(div A)^2 *)
  gaugeTerm = -(1/(2 xi)) divA^2;

  gaugeTerm
];

(* ================================================================ *)
(* === de Donder Gauge: -(1/2xi)(D_b)(D^b)                      === *)
(* === D_b = nabla_a h^a{}_b - (1/2) nabla_b h                  === *)
(* ================================================================ *)

BuildDeDonderGaugeTerm[field_, metric_, covd_, xi_:1] := Module[
  {traceH, d1, d2, gaugeTerm},

  (* Validate xi > 0 *)
  If[!NumericQ[xi] || xi <= 0,
    Throw["BuildDeDonderGaugeTerm: xi must be a positive number, got " <> ToString[xi]]
  ];

  (* Trace: h = g^{cd} h_{cd} *)
  traceH = metric[c, d] field[-c, -d];

  (* De Donder vector D_b = g^{ac} nabla_a h_{cb} - (1/2) nabla_b h *)
  (* First copy with free index -b, dummy indices a, c *)
  d1 = metric[a, c] covd[-a][field[-c, -b]] - (1/2) covd[-b][traceH];

  (* Second copy with free index -e, fresh dummy indices f, g *)
  d2 = metric[f, g] covd[-f][field[-g, -e]] - (1/2) covd[-e][traceH];

  (* Square: -(1/2xi) g^{be} D_b D_e — contract the two free indices *)
  gaugeTerm = -(1/(2 xi)) metric[b, e] d1 d2;

  gaugeTerm
];

(* ================================================================ *)
(* === Private implementations                                   === *)
(* ================================================================ *)

Begin["`Private`"];

AddGaugeFixingTerm[lagrangian_, gaugeTerm_] := lagrangian + gaugeTerm;

End[];
EndPackage[];
