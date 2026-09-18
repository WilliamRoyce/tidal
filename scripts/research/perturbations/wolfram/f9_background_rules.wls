#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPandFromRules SplitMetric DefProjectedTensor DefTensorPerturbation ProjectionAndBackgroundRuleQ RulesCovDsOfTensor Projectorh ExtractOrder noident tbar Detg RCOrder RCSort PREBUILD NOSV normu pabq UNCOMPUTED *)
(* f9_background_rules.wls -- GH #591: what xPand does with a background value for a field,
   and whether the zero-background path this programme uses is safe.

   WHAT WAS CLAIMED, AND WHAT IS TRUE. R-C reported that handing xPand a background torsion
   through the rule list makes the first-order piece vanish "for some rule lists". That was our
   own bug: f7's helper took three arguments and was called with four, so the call was never
   computed, and ExtractOrder of an unevaluated expression is 0. This script records what
   actually happens, with each check's expected answer and a control.

   THE MECHANISM (xPand.m). ToxPandFromRules conformally transforms, perturbs, then calls
   SplitPerturbations with the rule list. There, ProjectionAndBackgroundRuleQ (:2881) picks out
   rules whose left side is a background field and whose right side vanishes under the slice
   projector, and RulesCovDsOfTensor (:2840-2873) substitutes those EARLY, inside each field's
   induced decomposition. A right side that is a SUM of terms along the normal -- the shape a
   rank-3 field's background takes -- leaves an unexpanded sum inside Projectorh[...];
   ToCanonical cannot canonicalize it, xPerm throws, the Catch at :2840 turns that into a Null
   in the rule list, and the messages the user then sees (ReplaceAll::reps, Validate::inhom)
   point elsewhere. A zero right side, or a single-term one, goes through the same step
   harmlessly -- which is why vector and rank-2 backgrounds along the normal work.

   WHY IT DOES NOT AFFECT US. The programme expands about zero background torsion, exactly
   (user decision 2026-09-17, provisional pending the 18 September supervisor meeting). Part D
   runs every kind of split we use at zero background against recorded references, Part E
   reproduces our own bug and shows RCOrder refusing it, and Part F shows a failed background
   attempt does not damage later splits in the same session.

   ToCanonical::noident is deliberately visible here: RCTry quiets it, and in this bug it was
   the first sign of the real failure. Diagnostic calls use Catch directly, not RCTry.
   Usage: run_lane.sh f9 *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[c5];
DefConstantSymbol[tbar];
AA = ah[LI[0], LI[0]];
fixVol[x_] := NoScalar[(x /. Detg[] -> -1) /. Sqrt[y_] :> PowerExpand[Sqrt[y]]];
(* did a split actually run? -- keyed on xPand's own parameter name, never the one we passed *)
ran[x_] := !(x === Null || x === $Failed) && !FreeQ[x, $PerturbationParameter];
SetAttributes[stage, HoldRest];
stage[k_String, e_] := Module[{n = Length[$MessageList], r},
   r = Catch[e];
   RCSay["F9_" <> k <> "_MESSAGES", DeleteDuplicates[Drop[$MessageList, n]]];
   RCSayPlain["F9_" <> k <> "_RAN", ran[r]];
   r];
order[x_, k_] := RCCanon[fixVol[RCOrder[x, k]]];

(* ---------------- the fields and their perturbation rules ---------------- *)
DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
DefTensorPerturbation[dTor[LI[o], ia, -ib, -ic], Tor[ia, -ib, -ic], M4];
DefTensor[Av[-ia], M4];
DefTensorPerturbation[dAv[LI[o], -ia], Av[-ia], M4];
DefTensor[Bt[-ia, -ib], M4, Symmetric[{-ia, -ib}]];
DefTensorPerturbation[dBt[LI[o], -ia, -ib], Bt[-ia, -ib], M4];
DefProjectedTensor[TS1[], h];
DefProjectedTensor[TS3[], h];
DefProjectedTensor[AS1[], h];
DefProjectedTensor[AS2[], h];
DefProjectedTensor[BS1[], h];
DefProjectedTensor[AV1[-ia], h, TensorProperties -> {"Transverse"}];
DefProjectedTensor[T2b[], h, SpaceTimesOfDefinition -> {"Background"}];
pertT = dTor[LI[o_], ia_, ib_, ic_] :> (
   n[ia] (n[ib] cd[ic]@TS1[LI[o]] - n[ic] cd[ib]@TS1[LI[o]])
   + (n[ib] h[ia, ic] - n[ic] h[ia, ib]) TS3[LI[o]]);
pertA = dAv[LI[o_], ia_] :> (n[ia] AS1[LI[o]] + cd[ia]@AS2[LI[o]] + AV1[LI[o], ia]);
pertB = dBt[LI[o_], ia_, ib_] :> (n[ia] n[ib] BS1[LI[o]] + h[ia, ib] AS1[LI[o]]);
(* background rules, one per shape *)
bgTZero = Tor[ia_, ib_, ic_] -> 0;
bgTSum = Tor[ia_, ib_, ic_] :> (n[ib] h[ia, ic] - n[ic] h[ia, ib]) tbar;
bgTTwo = Tor[ia_, ib_, ic_] :> ((n[ib] h[ia, ic] - n[ic] h[ia, ib]) tbar + epsilonh[ib, ic, ia] T2b[]);
bgAZero = Av[ia_] -> 0;
bgAOne = Av[ia_] :> n[ia] tbar;
bgBOne = Bt[ia_, ib_] :> n[ia] n[ib] tbar;
rAny = SplitMetric[g, dg, h, "AnyGauge"];
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
rNoSV = rAny /. zeroSV;
LTodd = Sqrt[-Detg[]] c5 n[ib] Tor[ia, -ib, -ia];
LTq = Sqrt[-Detg[]] c5 Tor[ia, -ib, -ic] Tor[-ia, ib, ic];
LA1 = Sqrt[-Detg[]] c5 n[ia] Av[-ia];
LAq = Sqrt[-Detg[]] c5 Av[ia] Av[-ia];
LBq = Sqrt[-Detg[]] c5 Bt[ia, ib] Bt[-ia, -ib];

(* ================= Part A: the global switches every split relies on ================= *)
RCSay["F9_A_PERTURBATION_PARAMETER", $PerturbationParameter];
RCSay["F9_A_FIRST_ORDER_VECTOR_AND_TENSOR_FLAGS",
  {xAct`xPand`$FirstOrderVectorPerturbations, xAct`xPand`$FirstOrderTensorPerturbations}];

(* ================= Part B: how xPand classifies each background rule ================= *)
pabq[r_] := xAct`xPand`Private`ProjectionAndBackgroundRuleQ[r, h, n];
RCSayPlain["F9_B_CLASSIFIER_EXISTS", Names["xAct`xPand`Private`ProjectionAndBackgroundRuleQ"] =!= {}];
RCSay["F9_B_CLASSIFIED_AS_PROJECTED_BACKGROUND",
  {"torsion zero" -> pabq[bgTZero], "torsion sum along n" -> pabq[bgTSum],
   "torsion two-term" -> pabq[bgTTwo], "vector along n" -> pabq[bgAOne],
   "rank-2 along n" -> pabq[bgBOne], "perturbation rule" -> pabq[pertT]}];
RCSay["F9_B_PROJECTOR_OF_TORSION_SUM", RCCanon[Projector[h][Last[bgTSum]]]];
RCSay["F9_B_PROJECTOR_OF_TORSION_TWO_TERM", RCCanon[Projector[h][Last[bgTTwo]]]];
(* the early substitution itself: the pre-built split rules for each shape *)
pe = ExpandPerturbation[Perturbed[Conformal[g, ConformalMetricName[g, a[h]]][LTodd], 1]];
prebuild[rules_] := Flatten[
   xAct`xPand`Private`RulesCovDsOfTensor[pe, #, Select[rules, pabq[#] === True &], h] & /@ rules];
RCSayPlain["F9_B_PREBUILDER_EXISTS", Names["xAct`xPand`Private`RulesCovDsOfTensor"] =!= {}];
SetAttributes[stageRules, HoldRest];
stageRules[k_String, e_] := Module[{n = Length[$MessageList], r},
   r = Catch[e];
   RCSay["F9_" <> k <> "_MESSAGES", DeleteDuplicates[Drop[$MessageList, n]]];
   r];
pZ = stageRules["B_PREBUILD_TORSION_ZERO", prebuild[{pertT, bgTZero}]];
RCSayPlain["F9_B_PREBUILD_TORSION_ZERO_HAS_NULL", MemberQ[Flatten[{pZ}], Null]];
pS = stageRules["B_PREBUILD_TORSION_SUM", prebuild[{pertT, bgTSum}]];
RCSayPlain["F9_B_PREBUILD_TORSION_SUM_HAS_NULL", pS === Null || MemberQ[Flatten[{pS}], Null]];

(* ================= Part C: the full split for each shape ================= *)
sTSum = stage["C_TORSION_SUM_BACKGROUND", ToxPandFromRules[LTodd, Join[rAny, {pertT, bgTSum}], h, 1]];
sTTwo = stage["C_TORSION_TWO_TERM_BACKGROUND", ToxPandFromRules[LTq, Join[rAny, {pertT, bgTTwo}], h, 1]];
RCSayPlain["F9_C_TORSION_TWO_TERM_NOTE", "not classified as a projected background, so it is applied late and the split completes; its output has not been validated against anything"];
sAOne = stage["C_VECTOR_BACKGROUND_ALONG_N", ToxPandFromRules[LAq, Join[rAny, {pertA, bgAOne}], h, 1]];
RCSay["F9_C_VECTOR_BACKGROUND_ORDER0", order[sAOne, 0]];
RCSayPlain["F9_C_VECTOR_BACKGROUND_ORDER0_VS_HAND",
  RCVerdict[order[sAOne, 0], -c5 tbar^2 AA^4]];
RCSay["F9_C_VECTOR_BACKGROUND_ORDER1", order[sAOne, 1]];
RCSayPlain["F9_C_VECTOR_BACKGROUND_HAS_CROSS_TERM", !FreeQ[order[sAOne, 1], tbar AS1[LI[1], LI[0]]]];
sBOne = stage["C_RANK2_BACKGROUND_ALONG_N", ToxPandFromRules[LBq, Join[rAny, {pertB, bgBOne}], h, 1]];
RCSay["F9_C_RANK2_BACKGROUND_ORDER1", order[sBOne, 1]];
RCSayPlain["F9_C_RANK2_BACKGROUND_HAS_CROSS_TERM", !FreeQ[order[sBOne, 1], tbar BS1[LI[1], LI[0]]]];

(* ================= Part D: the zero-background path the programme uses ================= *)
refTodd = -3 c5 AA^4 TS3[LI[1], LI[0]];
z1 = stage["D_TORSION_ODD_ZERO", ToxPandFromRules[LTodd, Join[rAny, {pertT, bgTZero}], h, 1]];
RCSayPlain["F9_D_TORSION_ODD_VS_REFERENCE", RCVerdict[order[z1, 1], refTodd]];
z2 = stage["D_TORSION_ODD_ZERO_NOSV", ToxPandFromRules[LTodd, Join[rNoSV, {pertT, bgTZero}], h, 1]];
RCSayPlain["F9_D_TORSION_ODD_NOSV_VS_REFERENCE", RCVerdict[order[z2, 1], refTodd]];
z3 = stage["D_TORSION_QUADRATIC_ZERO", ToxPandFromRules[LTq, Join[rNoSV, {pertT, bgTZero}], h, 2]];
RCSay["F9_D_TORSION_QUADRATIC_ORDER2", order[z3, 2]];
RCSayPlain["F9_D_TORSION_QUADRATIC_VS_REFERENCE", RCVerdict[order[z3, 2],
   -6 c5 AA^4 TS3[LI[1], LI[0]]^2 + 2 c5 AA^4 cd[-ia][TS1[LI[1], LI[0]]] cd[ia][TS1[LI[1], LI[0]]]]];
z4 = stage["D_VECTOR_ZERO", ToxPandFromRules[LA1, Join[rAny, {pertA, bgAZero}], h, 1]];
RCSayPlain["F9_D_VECTOR_VS_REFERENCE", RCVerdict[order[z4, 1], -c5 AA^4 AS1[LI[1], LI[0]]]];
rNewt = SplitMetric[g, dg, h, "NewtonGauge"];
z5 = stage["D_METRIC_ONLY", ToxPandFromRules[RicciScalarCD[], rNewt, h, 1]];
RCSayPlain["F9_D_METRIC_ONLY_TERMS", Length[List @@ Expand[NoScalar[RCSort[RCOrder[z5, 1]]]]]];
RCSayPlain["F9_D_METRIC_ONLY_MATCHES_F2_COUNT_7", Length[List @@ Expand[NoScalar[RCSort[RCOrder[z5, 1]]]]] === 7];
(* control: a deliberately removed piece is detected, so the comparisons above can fail *)
RCSayPlain["F9_D_CONTROL_DROPPED_PIECE", RCVerdict[order[z1, 1] /. TS3[__] -> 0, refTodd]];

(* ================= Part E: our bug, reproduced, and the guard that now stops it ================= *)
splitTo3[L_, rules_, ord_] := ToxPandFromRules[L, rules, h, ord];
RCSayPlain["F9_E_LINT_WOULD_REFUSE_A_FOUR_ARGUMENT_CALL", "see wl_lint.py; this line calls it through ToExpression so the file itself stays clean"];
u4 = ToExpression["splitTo3[LTodd, Join[rAny, {pertT, bgTZero}], h, 1]"];
RCSayPlain["F9_E_FOUR_ARGUMENT_CALL_UNCOMPUTED", Head[u4] === splitTo3];
RCSayPlain["F9_E_EXTRACTORDER_OF_UNCOMPUTED_IS_ZERO", ExtractOrder[u4, 1] === 0];
RCSayPlain["F9_E_RCORDER_REFUSES_UNCOMPUTED", RCOrder[u4, 1] === $Failed];
RCSayPlain["F9_E_RCORDER_ACCEPTS_REAL_SPLIT", RCOrder[z1, 1] =!= $Failed];

(* ================= Part F: a failed background attempt does not damage later splits ================= *)
stage["F_FAILING_CALL_AGAIN", ToxPandFromRules[LTodd, Join[rAny, {pertT, bgTSum}], h, 1]];
z6 = stage["F_ZERO_BACKGROUND_AFTER_FAILURE", ToxPandFromRules[LTodd, Join[rAny, {pertT, bgTZero}], h, 1]];
RCSayPlain["F9_F_UNCHANGED_AFTER_FAILURE", RCVerdict[order[z6, 1], order[z1, 1]]];

RCSigns["f9_end"];
RCDone["f9"];
