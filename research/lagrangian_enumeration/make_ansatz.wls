#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  Generate the actual Lagrangian ansatz via MakeContractionAnsatz/MakeAnsatz
  for the CORE quadratic sectors. This is what the supervisor asked for:
  "the most general Lagrangian up to 2nd order."

  Separate from the full enumeration (general_quadratic_lagrangian.wls)
  to keep runtime manageable.
*)

Needs["xAct`xTensor`"];
Needs["xAct`xTras`"];
xAct`xTensor`$DefInfoQ = False;
xAct`xTensor`$UndefInfoQ = False;

DefManifold[M4, 4, {a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, q, r, s}];
DefMetric[-1, metric[-a, -b], CD, PrintAs -> "g"];

DefTensor[RCurv[-a, -b, -c, -d], M4, GenSet[-Cycles[{1, 2}], -Cycles[{3, 4}]]];
DefTensor[Tor[-a, -b, -c], M4, GenSet[-Cycles[{2, 3}]]];
DefTensor[FF[-a, -b], M4, Antisymmetric[{1, 2}]];
DefTensor[eps[-a, -b, -c, -d], M4, Antisymmetric[{1, 2, 3, 4}]];

Print["=== MakeContractionAnsatz: Core Quadratic Lagrangian ===\n"];

(* ============================================================ *)
(* CORE PARITY-EVEN: R^2 + T^2 + F^2 + RxF (11 terms)        *)
(* ============================================================ *)

Print["--- Parity-even sectors ---"];

rr = AllContractions[RCurv[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h]];
Print["  R x R: ", Length[rr], " terms"];

tt = AllContractions[Tor[-a,-b,-c] * Tor[-d,-e,-f]];
Print["  T x T: ", Length[tt], " terms"];

ff = AllContractions[FF[-a,-b] * FF[-c,-d]];
Print["  F x F: ", Length[ff], " terms"];

rf = AllContractions[RCurv[-a,-b,-c,-d] * FF[-e,-f]];
Print["  R x F: ", Length[rf], " terms"];

coreEven = Join[rr, tt, ff, rf];
Print["\nTotal parity-even core: ", Length[coreEven], " terms\n"];

(* Generate the ansatz with MakeAnsatz *)
lagrangianEven = MakeAnsatz[coreEven, ConstantPrefix -> "c"];

Print["Parity-even Lagrangian ansatz (MakeAnsatz output):"];
Print[lagrangianEven // InputForm];

Print["\nTerm-by-term:"];
Do[
  Print["  c", k, " * ", coreEven[[k]] // InputForm],
  {k, Length[coreEven]}
];

(* ============================================================ *)
(* CORE PARITY-ODD: eps*(R^2 + T^2 + F^2 + RxF) (21 terms)   *)
(* ============================================================ *)

Print["\n\n--- Parity-odd sectors ---"];

err = AllContractions[eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * RCurv[-i,-j,-k,-l]];
Print["  eps R x R: ", Length[err], " terms"];

ett = AllContractions[eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * Tor[-h,-i,-j]];
Print["  eps T x T: ", Length[ett], " terms"];

eff = AllContractions[eps[-a,-b,-c,-d] * FF[-e,-f] * FF[-g,-h]];
Print["  eps F x F: ", Length[eff], " terms"];

erf = AllContractions[eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * FF[-i,-j]];
Print["  eps R x F: ", Length[erf], " terms"];

coreOdd = Join[err, ett, eff, erf];
Print["\nTotal parity-odd core: ", Length[coreOdd], " terms\n"];

lagrangianOdd = MakeAnsatz[coreOdd, ConstantPrefix -> "d"];

Print["Parity-odd Lagrangian ansatz:"];
Print[lagrangianOdd // InputForm];

Print["\nTerm-by-term:"];
Do[
  Print["  d", k, " * ", coreOdd[[k]] // InputForm],
  {k, Length[coreOdd]}
];

(* ============================================================ *)
(* COMBINED ANSATZ                                               *)
(* ============================================================ *)

Print["\n\n=== SUMMARY ==="];
Print["Core parity-even: ", Length[coreEven], " independent terms"];
Print["Core parity-odd:  ", Length[coreOdd], " independent terms"];
Print["Total core:       ", Length[coreEven] + Length[coreOdd], " parameters"];
Print["(Plus 2 linear: Ricci scalar + cosmological constant)"];

(* Also run MakeContractionAnsatz directly to show it works *)
Print["\n\n=== MakeContractionAnsatz (direct call) ==="];
Print["R^2 ansatz:"];
rAnsatz = MakeContractionAnsatz[
  RCurv[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h],
  ConstantPrefix -> "alpha"];
Print[rAnsatz // InputForm];

Print["\nT^2 ansatz:"];
tAnsatz = MakeContractionAnsatz[
  Tor[-a,-b,-c] * Tor[-d,-e,-f],
  ConstantPrefix -> "beta"];
Print[tAnsatz // InputForm];

Print["\n=== DONE ==="];
