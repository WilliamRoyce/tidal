#!/usr/bin/env wolframscript
(* xPand 0.4.4, item 1: once any covariant derivative without a metric exists (here the chart
   derivative xCoba defines), every later ToxPandFromRules returns Null with
   InducedFrom::unknown. Run: wolframscript -file xpand_item1_chart.wls
 *)
Needs["xAct`xPand`"];
(* "control" lines check something that should work and print ok or MISMATCH; "reported"
   lines check the behavior this report describes and print reproduced or NOT REPRODUCED. *)
show[x_] := ToString[x, InputForm];
control[k_, label_, observed_, expected_] := Print["CHECK ", k, " control, ", label, ": ",
   If[observed === expected, "ok", "MISMATCH (observed " <> show[observed] <> ", expected " <> show[expected] <> ")"]];
reported[k_, label_, observed_, value_] := Print["CHECK ", k, " reported, ", label, ": ",
   If[observed === value, "reproduced", "NOT REPRODUCED (observed " <> show[observed] <> ", reported " <> show[value] <> ")"]];
newMessages[n_] := DeleteDuplicates[Drop[$MessageList, n]];
DefManifold[M4, 4, {ia, ib, ic, id, ie, ig}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
DefMetricPerturbation[g, dg, eps];
SetSlicing[g, n, -1, h, cd, {"|", "D"}, "FLFlat"];
DefMetricFields[g, dg, h, eps];
rules = SplitMetric[g, dg, h, "AnyGauge"];
before = ToxPandFromRules[RicciScalarCD[], rules, h, 1];
control[1, "the split works before any chart", before =!= Null && !FreeQ[before, $PerturbationParameter], True];
Needs["xAct`xCoba`"];
DefChart[B, M4, {0, 1, 2, 3}, {t[], x[], y[], z[]}];
control[2, "the chart derivative has no metric", MetricOfCovD[PDB], Null];
m0 = Length[$MessageList];
after = Catch[ToxPandFromRules[RicciScalarCD[], rules, h, 1]];
reported[3, "the same split after DefChart returns Null", after, Null];
reported[4, "and prints InducedFrom::unknown", MemberQ[newMessages[m0], HoldForm[InducedFrom::unknown]], True];
Print["DONE"];
