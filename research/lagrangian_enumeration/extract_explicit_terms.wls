#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  Extract ALL explicit index contractions from each sector of the
  general PGT+EM Lagrangian, printed in a format ready for TeX.

  Outputs each contraction as a Mathematica InputForm expression
  that can be manually transcribed to TeX notation.
*)

Needs["xAct`xTensor`"];
Needs["xAct`xTras`"];
xAct`xTensor`$DefInfoQ = False;
xAct`xTensor`$UndefInfoQ = False;

DefManifold[M4, 4, {a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, q, r, s, t, u, v, w}];
DefMetric[-1, metric[-a, -b], CD, PrintAs -> "g"];

(* Core tensors *)
DefTensor[R[-a, -b, -c, -d], M4, GenSet[-Cycles[{1, 2}], -Cycles[{3, 4}]], PrintAs -> "\!\(\*OverscriptBox[\(R\), \(~\)]\)"];
DefTensor[T[-a, -b, -c], M4, GenSet[-Cycles[{2, 3}]], PrintAs -> "T"];
DefTensor[F[-a, -b], M4, Antisymmetric[{1, 2}], PrintAs -> "F"];
DefTensor[eps[-a, -b, -c, -d], M4, Antisymmetric[{1, 2, 3, 4}], PrintAs -> "\[Epsilon]"];

(* Derivative-level tensors *)
DefTensor[DT[-d, -a, -b, -c], M4, GenSet[-Cycles[{3, 4}]], PrintAs -> "\[Del]T"];
DefTensor[DR[-e, -a, -b, -c, -d], M4, GenSet[-Cycles[{2, 3}], -Cycles[{4, 5}]], PrintAs -> "\[Del]R"];

Print["=== Extracting ALL explicit contractions ===\n"];

ExtractSector[name_String, expr_] := Module[{result},
  result = Quiet[AllContractions[expr], {AllContractions::nocontr}];
  If[result === {} || !ListQ[result],
    Print["--- ", name, ": 0 terms ---\n"];
    Return[{}]
  ];
  Print["--- ", name, ": ", Length[result], " terms ---"];
  Do[
    Print["  [", k, "] ", result[[k]] // InputForm],
    {k, Length[result]}
  ];
  Print[];
  result
];

(* ============================================================ *)
Print["=== CORE PARITY-ODD QUADRATIC ===\n"];

ExtractSector["eps R x R",
  eps[-a,-b,-c,-d] * R[-e,-f,-g,-h] * R[-i,-j,-k,-l]];

ExtractSector["eps T x T",
  eps[-a,-b,-c,-d] * T[-e,-f,-g] * T[-h,-i,-j]];

ExtractSector["eps F x F",
  eps[-a,-b,-c,-d] * F[-e,-f] * F[-g,-h]];

ExtractSector["eps R x F",
  eps[-a,-b,-c,-d] * R[-e,-f,-g,-h] * F[-i,-j]];

(* ============================================================ *)
Print["=== DERIVATIVE EXTENSIONS (parity-even) ===\n"];

ExtractSector["DT x F",
  DT[-d,-a,-b,-c] * F[-e,-f]];

ExtractSector["R x DT",
  R[-a,-b,-c,-d] * DT[-e,-f,-g,-h]];

ExtractSector["DT x DT",
  DT[-a,-b,-c,-d] * DT[-e,-f,-g,-h]];

(* ============================================================ *)
Print["=== DERIVATIVE EXTENSIONS (parity-odd) ===\n"];

ExtractSector["eps DT x F",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * F[-i,-j]];

ExtractSector["eps R x DT",
  eps[-a,-b,-c,-d] * R[-e,-f,-g,-h] * DT[-i,-j,-k,-l]];

ExtractSector["eps DT x DT",
  eps[-a,-b,-c,-d] * DT[-e,-f,-g,-h] * DT[-i,-j,-k,-l]];

(* ============================================================ *)
Print["=== CUBIC TERMS ===\n"];

ExtractSector["R x F x F",
  R[-a,-b,-c,-d] * F[-e,-f] * F[-g,-h]];

ExtractSector["T x T x F",
  T[-a,-b,-c] * T[-d,-e,-f] * F[-g,-h]];

ExtractSector["R x R x F",
  R[-a,-b,-c,-d] * R[-e,-f,-g,-h] * F[-i,-j]];

ExtractSector["DT x F x F",
  DT[-a,-b,-c,-d] * F[-e,-f] * F[-g,-h]];

(* ============================================================ *)
Print["=== CUBIC PARITY-ODD ===\n"];

ExtractSector["eps R x F x F",
  eps[-a,-b,-c,-d] * R[-e,-f,-g,-h] * F[-i,-j] * F[-k,-l]];

ExtractSector["eps T x T x F",
  eps[-a,-b,-c,-d] * T[-e,-f,-g] * T[-h,-i,-j] * F[-k,-l]];

ExtractSector["eps F x F x DT",
  eps[-a,-b,-c,-d] * F[-e,-f] * F[-g,-h] * DT[-i,-j,-k,-l]];

Print["\n=== DONE ==="];
