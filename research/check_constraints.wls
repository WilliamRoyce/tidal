#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  Check for Dimensional-Dependent Identities (DDIs) and projective
  invariance constraints on the PGT+EM Lagrangian in 4D.

  Part A: DDIs via ConstructDDIs (identities from over-antisymmetrization in D=4)
  Part B: Projective invariance (trace torsion T_μ constraints)
*)

(* === Load packages === *)
Print["=== Loading xAct + xTras ==="];
Needs["xAct`xTensor`"];
Needs["xAct`xTras`"];
xAct`xTensor`$DefInfoQ = False;
xAct`xTensor`$UndefInfoQ = False;
Print["Loaded.\n"];

(* === Manifold & metric === *)
DefManifold[M4, 4, {a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, q, r, s, t, u, v, w, x, y, z}];
DefMetric[-1, metric[-a, -b], CD, PrintAs -> "g"];

(* === Tensors === *)
DefTensor[RCurv[-a, -b, -c, -d], M4, GenSet[-Cycles[{1, 2}], -Cycles[{3, 4}]]];
DefTensor[Tor[-a, -b, -c], M4, GenSet[-Cycles[{2, 3}]]];
DefTensor[FF[-a, -b], M4, Antisymmetric[{1, 2}]];
DefTensor[eps[-a, -b, -c, -d], M4, Antisymmetric[{1, 2, 3, 4}]];
(* Derivative-level *)
DefTensor[DT[-d, -a, -b, -c], M4, GenSet[-Cycles[{3, 4}]]];

Print["All tensors defined.\n"];

(* ============================================================ *)
(* PART A: Dimensional-Dependent Identities                      *)
(* ============================================================ *)

Print["========================================"];
Print["PART A: Dimensional-Dependent Identities"];
Print["========================================\n"];

(* Helper to check DDIs for a sector *)
CheckDDIs[name_String, expr_] := Module[{contractions, ddis, nContr, nDDI},
  Print["--- ", name, " ---"];

  (* First get all contractions *)
  contractions = Quiet[AllContractions[expr], {AllContractions::nocontr}];
  If[contractions === {} || contractions === $Failed,
    Print["  No contractions found."];
    Return[{0, 0}]
  ];
  nContr = Length[contractions];
  Print["  AllContractions: ", nContr, " terms"];

  (* Now check for DDIs *)
  ddis = Quiet[ConstructDDIs[expr], {AllContractions::nocontr}];
  If[ddis === {} || ddis === $Failed,
    nDDI = 0,
    nDDI = Length[ddis]
  ];
  Print["  ConstructDDIs: ", nDDI, " identities"];
  Print["  Independent in 4D: ", nContr - nDDI, " (was ", nContr, ")"];
  If[nDDI > 0,
    Print["  ** DDI REDUCTION FOUND **"];
    Do[Print["    DDI ", k, ": ", ddis[[k]] // InputForm], {k, Min[nDDI, 3]}];
    If[nDDI > 3, Print["    ... and ", nDDI - 3, " more"]];
  ];
  Print[];
  {nContr, nDDI}
];

(* Check parity-odd non-derivative sectors *)
Print["Parity-odd non-derivative sectors:\n"];

{n1, d1} = CheckDDIs["eps R x R (currently 13)",
  eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * RCurv[-i,-j,-k,-l]];

{n2, d2} = CheckDDIs["eps T x T (currently 4)",
  eps[-a,-b,-c,-d] * Tor[-e,-f,-g] * Tor[-h,-i,-j]];

{n3, d3} = CheckDDIs["eps R x F (currently 3)",
  eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * FF[-i,-j]];

(* Check parity-even non-derivative sectors (should find 0 DDIs) *)
Print["Parity-even non-derivative sectors (control check):\n"];

{n4, d4} = CheckDDIs["R x R (currently 6)",
  RCurv[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h]];

{n5, d5} = CheckDDIs["T x T (currently 3)",
  Tor[-a,-b,-c] * Tor[-d,-e,-f]];

(* Check parity-odd derivative sectors *)
Print["Parity-odd derivative sectors:\n"];

{n6, d6} = CheckDDIs["eps DT x F (currently 6)",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * FF[-i,-j]];

{n7, d7} = CheckDDIs["eps R x DT (currently 36)",
  eps[-a,-b,-c,-d] * RCurv[-e,-f,-g,-h] * DT[-i,-j,-k,-l]];

(* eps DT x DT has 12 indices — check *)
{n8, d8} = CheckDDIs["eps DT x DT (currently 38)",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * DT[-i,-j,-k,-l]];


(* ============================================================ *)
(* PART B: Projective Invariance                                 *)
(* ============================================================ *)

Print["\n========================================"];
Print["PART B: Projective Invariance Analysis"];
Print["========================================\n"];

Print["Under projective transformation Gamma -> Gamma + delta^l_m xi_n:"];
Print["  T^l_mn -> T^l_mn + 2 delta^l_[m xi_n]"];
Print["  => T_m (trace) -> T_m + (D-1) xi_m = T_m + 3 xi_m"];
Print["  => S_m (axial) -> S_m (invariant)"];
Print["  => q_lmn (tensor) -> q_lmn (invariant)"];
Print[];
Print["Projective-invariant Lagrangians must be independent of T_m."];
Print["In the component basis (beta_1, beta_2, beta_3):"];
Print["  T^2 = beta_1 T_lmn T^lmn + beta_2 T_lmn T^mln + beta_3 T_m T^m"];
Print["  The trace-dependent part is beta_3 T_m T^m + trace pieces inside beta_1, beta_2"];
Print[];
Print["From the irreducible decomposition:"];
Print["  L_T2 = a q^2 + b S^2 + c T_m T^m"];
Print["  where c = beta_3 + (2 beta_1 + beta_2)/3"];
Print["  Projective invariance => c = 0 => beta_3 = -(2 beta_1 + beta_2)/3"];
Print["  This reduces T^2 from 3 -> 2 independent parameters.\n"];

Print["Similarly for parity-odd torsion (eps T x T):"];
Print["  The 4 contractions include trace-dependent ones."];
Print["  Projective invariance would eliminate those depending on T_m,"];
Print["  reducing from 4 to <= 3 terms.\n"];

Print["For derivative sectors (DT x F, etc.):"];
Print["  nabla T also contains trace: (nabla T)_trace = nabla_a T_b"];
Print["  Projective invariance constrains couplings involving nabla T_m.\n"];

(* ============================================================ *)
(* SUMMARY                                                       *)
(* ============================================================ *)

Print["\n========================================"];
Print["SUMMARY"];
Print["========================================\n"];

totalDDIs = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8;
Print["Total DDIs found: ", totalDDIs];
Print[];

Print["Sector-by-sector:"];
Print["  eps R x R: ", n1, " -> ", n1 - d1, " (", d1, " DDIs)"];
Print["  eps T x T: ", n2, " -> ", n2 - d2, " (", d2, " DDIs)"];
Print["  eps R x F: ", n3, " -> ", n3 - d3, " (", d3, " DDIs)"];
Print["  R x R: ", n4, " -> ", n4 - d4, " (", d4, " DDIs, control)"];
Print["  T x T: ", n5, " -> ", n5 - d5, " (", d5, " DDIs, control)"];
Print["  eps DT x F: ", n6, " -> ", n6 - d6, " (", d6, " DDIs)"];
Print["  eps R x DT: ", n7, " -> ", n7 - d7, " (", d7, " DDIs)"];
Print["  eps DT x DT: ", n8, " -> ", n8 - d8, " (", d8, " DDIs)"];
Print[];

If[totalDDIs == 0,
  Print["** NO DDIs found — all counts are exact in 4D. **"];
  Print["The AllContractions results are already the true 4D counts."],
  Print["** DDIs FOUND — counts need updating. **"];
  Print["Some parity-odd sectors have fewer independent terms in 4D."]
];

Print["\n=== DONE ==="];
