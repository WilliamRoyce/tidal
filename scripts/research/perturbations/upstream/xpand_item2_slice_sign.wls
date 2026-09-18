#!/usr/bin/env wolframscript
(* xPand 0.4.4, item 2: with n.n = +1 the induced metric is negative definite, but SetSlicing
   declares it with DefMetric[1, h, ...] (xPand.m line 1754), so the slice epsilon is normalized
   with the wrong sign. Run: wolframscript -file xpand_item2_slice_sign.wls *)
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
control[1, "n.n", canon[n[ia] n[-ia]], 1];
reported[2, "SignDetOfMetric[h] is declared +1 although h is negative definite", SignDetOfMetric[h], 1];
reported[3, "so epsilon[h] contracted with itself gives +6", canon[epsilonh[-ia, -ib, -ic] epsilonh[ia, ib, ic]], 6];
control[4, "while the four-index epsilon contracted with the normal, contracted with itself, gives -6",
  canon[epsilong[-ia, -ib, -ic, -id] n[id] epsilong[ia, ib, ic, ie] n[-ie]], -6];
Print["DONE"];
