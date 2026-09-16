#!/usr/bin/env wolframscript
(* cspell:words DefMatterFields Friedmann NOLAP SplitMatter SplitMetric ToxPand Vshu Vvhu offshell onshell xPand *)
(* repro_a2_tensor_eom.wls -- A2, the CONTROL: the FRW tensor equation from the perturbed
   Einstein equations (R-C #567, plan §5.5). Flat FRW, Lambda + perfect fluid as in the xPand
   paper's Appendix A (arXiv:1302.6174 TeX :1737-1770). Three runs:
     (i)   mixed indices  G^a_b     -> a^2 dG^a_b (TT)  should be  E'' + 2 H E' - D^2 E
     (ii)  lowered G_ab with the order-0 (background) rules substituted -> the same equation
     (iii) lowered G_ab WITHOUT them -> the residual (2H' + H^2 + kappa a^2 P) E_ab, printed:
           the first-order picture of the background-consistency mechanism (#501).
   Target from the derivation in the plan (a^2 dG^i_j = 1/2 (h'' + 2Hh' - D^2 h), h = 2E),
   MFB eq. 4.15; negative controls: friction 2H -> 3H, and D^2 dropped. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap]; DefConstantSymbol[Lam];
DefMatterFields[u, du, h, eps];
DefTensor[Tf[-ia, -ib], M4];
IndexSet[Tf[ia_, ib_], (\[Rho]u[] + Pu[]) u[ia] u[ib] + Pu[] g[ia, ib]];
MyGR[ia_, ib_] := EinsteinCD[ia, ib] + g[ia, ib] Lam/kap - kap Tf[ia, ib];

(* --- rules: AnyGauge metric split with scalars and vectors switched off; fluid perturbations off --- *)
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
rulesTT = SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV;
rulesMat = SplitMatter[u, du, -1, h, "AnyGauge", 1] /. {Vvhu[__] -> 0, Vshu[__] -> 0};
rulesAll = Join[rulesTT, rulesMat];
RCSay["RULES_TT", rulesTT];

zeroFluidPert = {\[Rho]u[LI[1], __] -> 0, Pu[LI[1], __] -> 0};

(* --- mixed indices, order 0 and 1 --- *)
{tM, resM} = AbsoluteTiming[ToxPandFromRules[MyGR[ia, -ib], rulesAll, h, 1] /. zeroFluidPert];
RCSayPlain["A2_MIXED_WALL_S", Round[tM, 0.01]];
bg = ExtractOrder[ah[]^2 resM, 0];
bgTT = RCCanon[ExtractComponents[bg, h, {"Time", "Time"}]];
bgSS = RCCanon[ExtractComponents[bg, h, {"Space", "Space"}]];
RCSay["A2_BG_TIME_TIME", bgTT];
RCSay["A2_BG_SPACE_SPACE", bgSS];
(* the paper's :1777-1783: Lambda a^2/kappa - 3H^2 + kappa a^2 rho ; (Lambda a^2/kappa - H^2 - 2H') h^a_b *)
bgSSscalar = bgSS /. {h[_, _] -> 1, delta[_, _] -> 1};
bgRules = First@Solve[{bgTT == 0, bgSSscalar == 0}, {\[Rho]u[LI[0], LI[0]], Hh[LI[0], LI[1]]}];
RCSay["A2_BG_RULES", bgRules];

ttM = RCCanon[ExtractComponents[ExtractOrder[ah[]^2 resM, 1], h, {"Space", "Space"}]];
RCSay["A2_TT_MIXED", ttM];
target = Eth[LI[1], LI[2], ia, -ib] + 2 Hh[LI[0], LI[0]] Eth[LI[1], LI[1], ia, -ib] - cd[-ic]@cd[ic]@Eth[LI[1], LI[0], ia, -ib];
controlFriction = target /. (2 Hh[LI[0], LI[0]]) -> 3 Hh[LI[0], LI[0]];
controlNoLap = Eth[LI[1], LI[2], ia, -ib] + 2 Hh[LI[0], LI[0]] Eth[LI[1], LI[1], ia, -ib];
{cM, vM} = RCNormalize[ttM, target];
RCSayPlain["A2_MIXED_NORMALIZATION", ToString[cM, InputForm] <> " :: " <> vM];
RCReport["a2_mixed", ttM, If[cM === None, target, cM target], If[cM === None, controlFriction, cM controlFriction]];
RCSayPlain["A2_MIXED_CONTROL_NOLAP", RCVerdict[ttM, If[cM === None, controlNoLap, cM controlNoLap]]];
(* also the mixed-index result with the background rules applied must be unchanged (no background term) *)
RCSayPlain["A2_MIXED_ONSHELL_SAME", RCVerdict[ttM, ttM /. bgRules]];

(* --- lowered indices --- *)
{tL, resL} = AbsoluteTiming[ToxPandFromRules[MyGR[-ia, -ib], rulesAll, h, 1] /. zeroFluidPert];
RCSayPlain["A2_LOWERED_WALL_S", Round[tL, 0.01]];
ttL = RCCanon[ExtractComponents[ExtractOrder[ah[]^2 resL, 1], h, {"Space", "Space"}]];
RCSay["A2_TT_LOWERED_OFFSHELL", ttL];
ttLon = RCCanon[ttL /. bgRules];
RCSay["A2_TT_LOWERED_ONSHELL", ttLon];
targetL = Eth[LI[1], LI[2], -ia, -ib] + 2 Hh[LI[0], LI[0]] Eth[LI[1], LI[1], -ia, -ib] - cd[-ic]@cd[ic]@Eth[LI[1], LI[0], -ia, -ib];
controlL = targetL /. (2 Hh[LI[0], LI[0]]) -> 3 Hh[LI[0], LI[0]];
candsL = Flatten[Outer[Times, {1, -1, 2, -2, 1/2, -1/2}, {1, ah[LI[0], LI[0]]^2, 1/ah[LI[0], LI[0]]^2}]];
{cL, vL} = RCNormalize[ttLon, targetL, candsL];
RCSayPlain["A2_LOWERED_ONSHELL_NORMALIZATION", ToString[cL, InputForm] <> " :: " <> vL];
RCReport["a2_lowered_onshell", ttLon, If[cL === None, targetL, cL targetL], If[cL === None, controlL, cL controlL]];
resid = RCCanon[ttL - ttLon];
RCSay["A2_RESIDUAL_OFFSHELL_MINUS_ONSHELL", resid];
RCSayPlain["A2_RESIDUAL_IS_ZERO", resid === 0];
RCSayPlain["A2_LOWERED_OFFSHELL_VS_TARGET", Last@RCNormalize[ttL, targetL, candsL]];
RCDone["a2"];
