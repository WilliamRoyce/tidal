#!/usr/bin/env wolframscript
(* xPand 0.4.4, item 4: ExtractComponents' "Time" projection multiplies an upper index by -1
   (xPand.m line 3042), which assumes n.n = -1. With n.n = +1 the time component of W0 n^a
   comes out as -W0, while n_a (W0 n^a) / (n.n) is +W0.
   Run: wolframscript -file xpand_item4_time_projection.wls *)
Needs["xAct`xPand`"];
(* "control" lines check something that should work and print ok or MISMATCH; "reported"
   lines check the behavior this report describes and print reproduced or NOT REPRODUCED. *)
show[x_] := ToString[x, InputForm];
control[k_, label_, observed_, expected_] := Print["CHECK ", k, " control, ", label, ": ",
   If[observed === expected, "ok", "MISMATCH (observed " <> show[observed] <> ", expected " <> show[expected] <> ")"]];
reported[k_, label_, observed_, value_] := Print["CHECK ", k, " reported, ", label, ": ",
   If[observed === value, "reproduced", "NOT REPRODUCED (observed " <> show[observed] <> ", reported " <> show[value] <> ")"]];
canon[x_] := ToCanonical[ContractMetric[x]];
DefManifold[M4, 4, {ia, ib, ic, id, ie, ig}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
DefMetricPerturbation[g, dg, eps];
SetSlicing[g, n, 1, h, cd, {"|", "D"}, "FLFlat"];
DefProjectedTensor[W0[], h];
w = W0[LI[1], LI[0]];
byHand = canon[n[-ia] (w n[ia]) / (n[ib] n[-ib])];
control[1, "by hand, n_a (W0 n^a) / (n.n) is W0", byHand, w];
reported[2, "ExtractComponents gives the \"Time\" component of W0 n^a as -W0", canon[ExtractComponents[w n[ia], h, {"Time"}]], -w];
Print["DONE"];
