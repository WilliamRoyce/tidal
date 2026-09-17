#!/usr/bin/env wolframscript
(* cspell:words xMAG ChangeCurvature ChangeCovD ToDistortion BreakContorsion TorsionToDistortion changeRiemann MasterOf RicciScalarCDL RicciScalarCDM ContorsionCDM ChristoffelCDCDL SafeCanonical RCStackOn RCSigns COVDS DEFCOVD MAGCHRISTOFFELQ UNDEFCOVD xtensor TORSIONTODISTORTION MASTEROF *)
(* f3_changecurvature.wls -- S2, the one behavior R-C reported that reading did not explain:
   with xMAG loaded, ChangeCurvature on a torsionful connection declared WITHOUT Master
   returns its input unchanged and prints nothing. One kernel, no xPand, three stages:
     (1) plain xTensor: the same call on the same kind of connection -- the control that
         says what the answer should look like;
     (2) xMAG loaded: repeat on the pre-existing connection, and on one declared after;
     (3) the same with Master -> g (the documented declaration) -- the twin.
   xMAG replaces xTensor's private changeRiemann (xMAG.m:1258) and pushes the torsion term
   through its own 1-argument TorsionToDistortion, whose dispatcher (xMAG.m:1705) calls
   MasterOf and has no definition for a Null master -- the suspected mechanism. *)
Needs["xAct`xTensor`"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]];
RCSigns["plain_xtensor"];
DefManifold[M4, 4, {ia, ib, ic, id, ie}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
DefTensor[vv[ia], M4];
DefCovD[CDL[-ia], {"&", "DL"}, FromMetric -> g, Torsion -> True];      (* no Master *)
RCSay["F3_1_MASTEROF_CDL", RCTry[MasterOf[CDL]]];
r1 = RCNewMessages["plain_changecurvature", RCTry[ChangeCurvature[RicciScalarCDL[], CDL, CD]]];
RCSay["F3_1_PLAIN_CHANGECURVATURE", r1];
RCSayPlain["F3_1_PLAIN_CHANGECURVATURE_EXPANDED", RCTry[r1 =!= RicciScalarCDL[] && !FreeQ[r1, RicciCD | ChristoffelCDCDL]]];
RCSay["F3_1_PLAIN_HEADS", RCTry[Union@Cases[r1, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]]];
r1v = RCNewMessages["plain_changecovd", RCTry[ToCanonical[ChangeCovD[CDL[-ia]@vv[ib], CDL, CD] - CD[-ia]@vv[ib]]]];
RCSay["F3_1_PLAIN_CDL_MINUS_CD_ON_A_VECTOR", r1v];

(* ---------------- stage 2: xMAG loaded, same connection ---------------- *)
Needs["xAct`xMAG`"];
RCSigns["after_xmag"];
RCStackOn[{HoldForm[Validate::inhom], HoldForm[InducedFrom::unknown]}];
r2 = RCNewMessages["xmag_changecurvature_cdl", RCTry[ChangeCurvature[RicciScalarCDL[], CDL, CD]]];
RCSay["F3_2_XMAG_CHANGECURVATURE_CDL", r2];
RCSayPlain["F3_2_RETURNS_ITS_INPUT", RCTry[r2 === RicciScalarCDL[]]];
RCSayPlain["F3_2_SAME_AS_PLAIN_RESULT", RCTry[RCVerdict[r2, r1]]];
(* the suspected mechanism, called directly *)
RCSay["F3_2_TORSIONTODISTORTION_1ARG_ON_CDL_TORSION",
  RCTry[TorsionToDistortion[TorsionCDL[ia, -ib, -ic]]]];
RCSay["F3_2_TORSIONTODISTORTION_3ARG_NULL_MASTER",
  RCTry[TorsionToDistortion[TorsionCDL[ia, -ib, -ic], CDL, MasterOf[CDL]]]];
RCSay["F3_2_MAGCHRISTOFFELQ_CDL", RCTry[MAGChristoffelQ[CDL, g]]];
(* and a connection declared the same way AFTER xMAG is loaded *)
okN = RCNewMessages["defcovd_nomaster_after_xmag",
  RCTry[(DefCovD[CDN[-ia], {"%", "DN"}, FromMetric -> g, Torsion -> True]; True)]];
RCSayPlain["F3_2_DEFCOVD_NOMASTER_AFTER_XMAG_OK", okN === True];
r2b = RCNewMessages["xmag_changecurvature_cdn", RCTry[ChangeCurvature[RicciScalarCDN[], CDN, CD]]];
RCSayPlain["F3_2_CDN_RETURNS_ITS_INPUT", RCTry[r2b === RicciScalarCDN[]]];

(* ---------------- stage 3: the twin, declared as xMAG documents ---------------- *)
okM = RCNewMessages["defcovd_master_after_xmag",
  RCTry[(DefCovD[CDM[-ia], {"#", "DM"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True]; True)]];
RCSayPlain["F3_3_DEFCOVD_WITH_MASTER_OK", okM === True];
RCSay["F3_3_MASTEROF_CDM", RCTry[MasterOf[CDM]]];
r3 = RCNewMessages["xmag_changecurvature_cdm", RCTry[ChangeCurvature[RicciScalarCDM[], CDM, CD]]];
RCSayPlain["F3_3_CDM_EXPANDED", RCTry[r3 =!= RicciScalarCDM[] && !FreeQ[r3, RicciCD | ContorsionCDM]]];
RCSay["F3_3_CDM_RESULT", r3];
RCSay["F3_3_CDM_HEADS", RCTry[Union@Cases[r3, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]]];
RCSayPlain["F3_3_TODISTORTION_ON_CDM_RAN", RCTry[
  If[FreeQ[BreakContorsion[ToDistortion[RicciScalarCDM[], CDM, g], CDM], $Failed], "ran", "failed"]]];
(* the same Master-less connection is still in the session: does it spoil the healthy one? *)
RCSayPlain["F3_3_TODISTORTION_ON_CDM_WITH_CDL_PRESENT", RCTry[
  If[FreeQ[ToDistortion[RicciScalarCDM[], CDM, g], $Failed], "ran", "failed"]]];
RCSayPlain["F3_3_UNDEFCOVD_CDL_AND_CDN", RCTry[(UndefCovD[CDL]; UndefCovD[CDN]; True)] === True];
RCSayPlain["F3_3_TODISTORTION_ON_CDM_AFTER_REMOVAL", RCTry[
  If[FreeQ[BreakContorsion[ToDistortion[RicciScalarCDM[], CDM, g], CDM], $Failed], "ran", "failed"]]];
RCStackOff[];
RCSigns["at_end"];
RCDone["f3"];
