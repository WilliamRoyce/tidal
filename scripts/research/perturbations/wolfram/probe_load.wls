#!/usr/bin/env wolframscript
(* cspell:words DefMatterFields DefMetricFields SetSlicing SplitMatter SplitMetric ToxPand gah phiB psiB usages xPand *)
(* probe_load.wls -- does xPand 0.4.4 load on the certified bundle, and does the paper's
   minimal example (arXiv:1302.6174 TeX :1350-1373) reproduce? Also reads back the label
   layout and the Newtonian-gauge rules so the Ma & Bertschinger dictionary is fixed by
   the package, not assumed. R-C (#567). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
Do[Print["RC_USAGE ", f, " :: ", RCUsage["xAct`xPand`", f]],
   {f, {"SetSlicing", "DefMetricFields", "DefMatterFields", "DefProjectedTensor", "SplitMetric", "SplitMatter",
        "ToxPand", "ToxPandFromRules", "SplitPerturbations", "ExtractComponents", "ExtractOrder",
        "Conformal", "ConformalWeight", "DefConformalMetric", "$ConformalTime"}}];
RCSay["LIST_OF_GAUGES", $ListOfGauges];
RCSay["LIST_OF_SPACETYPES", $ListOfSpaceTypes];

RCGeometry[-1, "FLFlat"];
RCSay["TENSORS_AFTER_GEOMETRY", $Tensors];
RCSay["AH_LABELED", ah[]];
RCSay["HH_LABELED", Hh[]];
RCSay["CD_AH_ZERO", cd[-ia]@ah[] === 0];
RCSay["AH_PRIME", ah[LI[0], LI[1]]];
RCSay["CONFORMAL_METRIC_NAME", gah2];

(* --- the paper's minimal example, TeX :1350-1360; output Out[47] :1368-1373, K = 0 on FLFlat --- *)
DefProjectedTensor[phiB[], h];
DefProjectedTensor[psiB[], h];
{tR, MyR} = AbsoluteTiming[ExpandPerturbation@Perturbed[Conformal[g, gah2][RicciScalarCD[]], 1]];
RCSayPlain["MINIMAL_PERTURB_WALL_S", Round[tR, 0.01]];
MyGauge = {dg[LI[ord_], ia_, ib_] :> -2 n[ia] n[ib] phiB[LI[ord]] - 2 h[ia, ib] psiB[LI[ord]]};
{tS, out47} = AbsoluteTiming[SplitPerturbations[ah[]^2 MyR, MyGauge, h]];
RCSayPlain["MINIMAL_SPLIT_WALL_S", Round[tS, 0.01]];
RCSay["OUT47", out47];
RCSay["OUT47_LABELS", Union@Cases[out47, (phiB | psiB | Hh | ah)[__], Infinity]];
target47 = 6 Hh[LI[0], LI[0]]^2 + 6 Hh[LI[0], LI[1]] +
  RCEps (-12 Hh[LI[0], LI[0]]^2 phiB[LI[1], LI[0]] - 12 Hh[LI[0], LI[1]] phiB[LI[1], LI[0]]
       - 6 Hh[LI[0], LI[0]] phiB[LI[1], LI[1]] - 18 Hh[LI[0], LI[0]] psiB[LI[1], LI[1]]
       - 6 psiB[LI[1], LI[2]] - 2 cd[-ia]@cd[ia]@phiB[LI[1], LI[0]] + 4 cd[-ia]@cd[ia]@psiB[LI[1], LI[0]]);
control47 = target47 /. (-18) -> -17;
RCReport["out47", out47, target47, control47];
RCReport["out47_order0", ExtractOrder[out47, 0], 6 Hh[LI[0], LI[0]]^2 + 6 Hh[LI[0], LI[1]],
         6 Hh[LI[0], LI[0]]^2 + 5 Hh[LI[0], LI[1]]];

(* --- definitional read-back: the Newtonian-gauge metric rules and the fluid rules --- *)
rulesN = SplitMetric[g, dg, h, "NewtonGauge"];
RCSay["RULES_NEWTON", rulesN];
rulesA = SplitMetric[g, dg, h, "AnyGauge"];
RCSay["RULES_ANY", rulesA];
DefMatterFields[u, du, h, eps];
rulesM = SplitMatter[u, du, -1, h, "NewtonGauge", 1];
RCSay["RULES_MATTER_NEWTON", rulesM];
RCSay["LIST_BG_ONLY", $ListFieldsBackgroundOnly[h]];
RCSay["LIST_PERT_ONLY", $ListFieldsPerturbedOnly[h]];
RCSay["LIST_BG_AND_PERT", $ListFieldsBackgroundAndPerturbed[h, u]];
RCSay["TENSORS_FINAL", $Tensors];
RCDone["probe_load"];
