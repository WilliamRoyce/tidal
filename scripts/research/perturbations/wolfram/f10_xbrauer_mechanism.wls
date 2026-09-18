#!/usr/bin/env wolframscript
(* cspell:words xPand xBrauer MetricOfTensor FirstMetricQ MasterOf ContractMetric ContractMetric1 DefProjectedTensor OrthogonalTo ProjectedWith Evh UpValues *)
(* f10_xbrauer_mechanism.wls -- after loading xBrauer, why does ContractMetric stop contracting
   xPand's induced metric h into xPand's own fields, and what exactly decides it?

   WHAT WAS CLAIMED, AND WHY IT WAS WITHDRAWN. The R-C memo (section 3.2 item 3) and the first
   xBrauer draft said the discriminator is the tensor's Master: epsilonh (Master h) contracts,
   Evh (no Master) does not. This follow-up's plan found a hand-declared tensor with no Master that still
   contracts after xBrauer, so the explanation was withdrawn. But that tensor was declared with
   ProjectedWith -> {h[...]}, which in xTensor installs the tensor's OWN automatic projector rule
   (xTensor.m:4258) -- a second way to contract h that never goes through ContractMetric1. The
   comparison changed two things at once.

   THIS SCRIPT changes one property at a time, in one xPand session, before and after xBrauer:
     TvO   DefTensor, OrthogonalTo only             (no Master, no own projector rule)
     TvPW  DefTensor, OrthogonalTo + ProjectedWith  (own projector rule)
     TvM   DefTensor, OrthogonalTo + Master -> h    (a Master)
     Tp    DefProjectedTensor, as xPand's fields     (xPand's route)
   and records, for each: whether the bare product h.T is already contracted when built;
   whether ContractMetric contracts it before and after xBrauer; whether the tensor carries a
   rule of its own that mentions h; and the value of xBrauer's re-installed guard
   (xBrauer.m:1947-1974), MetricOfTensor[T] === h || (MetricOfTensor[T] === Null &&
   FirstMetricQ[h]).
   The hypothesis it can refute: xBrauer's guard blocks exactly the tensors that have neither a
   Master nor a projector rule of their own -- which is every field xPand's DefProjectedTensor
   makes.
   RESULT (run 2026-09-18): confirmed. For every tensor that reaches ContractMetric the guard
   predicts the outcome exactly, and the guard is not constant across them, so the check could
   have failed. TvPW is contracted by its own rule the moment the product is built, which is why
   the plan's counterexample looked like one. The memo's original explanation stands; the
   withdrawal is itself withdrawn.
   Usage: run_lane.sh f10 *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefTensor[TvO[-ia], M4, OrthogonalTo -> {n[ia]}];
DefTensor[TvPW[-ia], M4, OrthogonalTo -> {n[ia]}, ProjectedWith -> {h[ia, -ib]}];
DefTensor[TvM[-ia], M4, OrthogonalTo -> {n[ia]}, Master -> h];
DefProjectedTensor[Tp[-ia], h];
tensors = {TvO, TvPW, TvM, Tp};
prod[TvO] := h[ia, ib] TvO[-ia];
prod[TvPW] := h[ia, ib] TvPW[-ia];
prod[TvM] := h[ia, ib] TvM[-ia];
prod[Tp] := h[ia, ib] Tp[LI[1], LI[0], -ia];
contracted[x_] := FreeQ[x, h];
ownRuleMentionsH[t_] := Count[UpValues[t], r_ /; !FreeQ[r, h]];
report[tag_String] := Do[
   RCSayPlain["F10_" <> tag <> "_" <> ToString[t],
     ToString[{"built" -> contracted[prod[t]], "ContractMetric" -> contracted[ContractMetric[prod[t]]]},
       InputForm]], {t, tensors}];

(* ---------------- properties that do not depend on xBrauer ---------------- *)
RCSay["F10_METRICS", $Metrics];
RCSay["F10_MASTEROF", {#, MasterOf[#]} & /@ tensors];
RCSay["F10_OWN_RULES_MENTIONING_H", {#, ownRuleMentionsH[#]} & /@ tensors];
RCSay["F10_ORTHOGONAL_TO_N", {#, OrthogonalToVectorQ[n][#]} & /@ tensors];

(* ---------------- before xBrauer ---------------- *)
report["BEFORE"];

(* ---------------- after xBrauer ---------------- *)
RCNewMessages["needs_xbrauer", Needs["xAct`xBrauer`"]];
report["AFTER"];
guard[t_] := (xAct`xBrauer`MetricOfTensor[t] === h ||
   (xAct`xBrauer`MetricOfTensor[t] === Null && xAct`xBrauer`Private`FirstMetricQ[h]));
RCSay["F10_GUARD_VALUE", {#, guard[#]} & /@ tensors];
RCSay["F10_METRICOFTENSOR_HEAD_FORM", {#, xAct`xBrauer`MetricOfTensor[#]} & /@ tensors];
RCSayPlain["F10_FIRSTMETRICQ_H", xAct`xBrauer`Private`FirstMetricQ[h]];

(* ---------------- the verdict, computed ----------------
   A tensor whose own projector rule contracts h when the product is BUILT never reaches
   ContractMetric at all, so the guard cannot matter for it. For every other tensor the
   hypothesis is that the guard alone predicts whether ContractMetric contracts after xBrauer.
   (A first version of this check counted "own rules mentioning h" as a bypass; that is too
   crude -- DefProjectedTensor installs nineteen such rules, none of which contracts h.T -- and
   it reported the hypothesis false for that reason alone.) *)
after[t_] := contracted[ContractMetric[prod[t]]];
builtContracted = Select[tensors, contracted[prod[#]] &];
reachesContractMetric = Complement[tensors, builtContracted];
RCSay["F10_CONTRACTED_WHEN_BUILT", builtContracted];
RCSay["F10_AFTER_VS_GUARD", {#, after[#], guard[#]} & /@ reachesContractMetric];
RCSayPlain["F10_GUARD_PREDICTS_CONTRACTMETRIC_AFTER_XBRAUER",
  And @@ (after[#] === guard[#] & /@ reachesContractMetric)];
RCSayPlain["F10_CONTROL_GUARD_IS_NOT_CONSTANT",
  Length[DeleteDuplicates[guard /@ reachesContractMetric]] > 1];
RCSayPlain["F10_ONE_VARIABLE_CONTRASTS",
  ToString[{"TvO vs TvPW (own projector rule)" -> {after[TvO], after[TvPW]},
            "TvO vs TvM (Master)" -> {after[TvO], after[TvM]},
            "TvO vs Tp (xPand's route)" -> {after[TvO], after[Tp]}}, InputForm]];
RCSigns["f10_end"];
RCDone["f10"];
