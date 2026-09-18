#!/usr/bin/env wolframscript
(* xPand 0.4.4, item 3: a background rule whose value is a sum of terms along the normal breaks
   the preparation of split rules; a zero background, and single-term backgrounds for vector
   and rank-2 fields, work. Run: wolframscript -file xpand_item3_background_rule.wls *)
Needs["xAct`xPand`"];
(* "control" lines check something that should work and print ok or MISMATCH; "reported"
   lines check the behavior this report describes and print reproduced or NOT REPRODUCED. *)
show[x_] := ToString[x, InputForm];
control[k_, label_, observed_, expected_] := Print["CHECK ", k, " control, ", label, ": ",
   If[observed === expected, "ok", "MISMATCH (observed " <> show[observed] <> ", expected " <> show[expected] <> ")"]];
reported[k_, label_, observed_, value_] := Print["CHECK ", k, " reported, ", label, ": ",
   If[observed === value, "reproduced", "NOT REPRODUCED (observed " <> show[observed] <> ", reported " <> show[value] <> ")"]];
ran[x_] := x =!= Null && x =!= $Failed && !FreeQ[x, $PerturbationParameter];
DefManifold[M4, 4, {ia, ib, ic, id, ie, ig}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
DefMetricPerturbation[g, dg, eps];
SetSlicing[g, n, -1, h, cd, {"|", "D"}, "FLFlat"];
DefMetricFields[g, dg, h, eps];
DefConstantSymbol[c5];
DefConstantSymbol[tbar];
DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
DefTensorPerturbation[dTor[LI[o], ia, -ib, -ic], Tor[ia, -ib, -ic], M4];
DefTensor[Av[-ia], M4];
DefTensorPerturbation[dAv[LI[o], -ia], Av[-ia], M4];
DefProjectedTensor[TS1[], h];
DefProjectedTensor[TS3[], h];
DefProjectedTensor[AS1[], h];
pert = dTor[LI[o_], ia_, ib_, ic_] :> (
   n[ia] (n[ib] cd[ic]@TS1[LI[o]] - n[ic] cd[ib]@TS1[LI[o]])
   + (n[ib] h[ia, ic] - n[ic] h[ia, ib]) TS3[LI[o]]);
pertA = dAv[LI[o_], ia_] :> n[ia] AS1[LI[o]];
metricRules = SplitMetric[g, dg, h, "AnyGauge"];
L = Sqrt[-Detg[]] c5 n[ib] Tor[ia, -ib, -ia];
zero = ToxPandFromRules[L, Join[metricRules, {pert, Tor[ia_, ib_, ic_] -> 0}], h, 1];
control[1, "with a zero background the split runs", ran[zero], True];
m0 = Length[$MessageList];
sum = Catch[ToxPandFromRules[L, Join[metricRules,
     {pert, Tor[ia_, ib_, ic_] :> (n[ib] h[ia, ic] - n[ic] h[ia, ib]) tbar}], h, 1]];
reported[2, "with a background that is a sum of terms along the normal it does not", ran[sum], False];
reported[3, "and TranslatePerm::invalid is among the messages",
  MemberQ[DeleteDuplicates[Drop[$MessageList, m0]], HoldForm[TranslatePerm::invalid]], True];
LA = Sqrt[-Detg[]] c5 Av[ia] Av[-ia];
vec = ToxPandFromRules[LA, Join[metricRules, {pertA, Av[ia_] :> n[ia] tbar}], h, 1];
control[4, "a single-term vector background along the normal works", ran[vec], True];
Print["DONE"];
