(* ::Package:: *)
(* cspell:words DefMatterFields DefMetricFields SetSlicing ToxPand WXF inds normu xPand xPert xTensor *)
(* ==============================================================================
   RCSetup.wl -- xPand-side setup for the R-C Wolfram probes (#567)
   ------------------------------------------------------------------------------
   PURPOSE  load xPand on the certified bundle FIRST (so every later symbol in
            this file and in the calling script resolves into the xAct contexts),
            then the package-independent harness (RCSetupCore.wl), an FRW slicing
            with a selectable normal norm, and xPand's own derivative sorting.
            A step that must load another package first gets RCSetupCore.wl only.
   USAGE    Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]]
            from a .wls run as  wolframscript -file X.wls outdir=DIR key=value ...
   NAMES    xPand exports a, H, P, T, k, Es, Et, Ev, Bs, Bv, Ls, Lv, Vs, Vv, V0,
            Xi, CS, Connection, nt, av -- never reuse them. Indices are ia..is.
   ============================================================================== *)

$RCLoadT0 = AbsoluteTime[];
$RCLoadFailed = False;
Check[Needs["xAct`xPand`"], $RCLoadFailed = True];
$RCLoadMessages = $MessageList;
$RCLoadWall = AbsoluteTime[] - $RCLoadT0;

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]];

RCLoad[] := (
  RCSayPlain["LOAD_FAILED", $RCLoadFailed];
  RCSayPlain["LOAD_WALL_S", Round[$RCLoadWall, 0.1]];
  RCSay["LOAD_MESSAGES", $RCLoadMessages];
  RCSay["WOLFRAM_VERSION", $Version];
  RCSay["XTENSOR_VERSION", xAct`xTensor`$Version];
  RCSay["XPERT_VERSION", xAct`xPert`$Version];
  RCSay["XPAND_VERSION", xAct`xPand`$Version];
  RCSay["CONFORMAL_TIME", $ConformalTime];
  RCSay["XPAND_SYMBOL_COUNT", Length[Names["xAct`xPand`*"]]];
  RCSigns["after_xpand"];
);

(* FRW slicing. normu = -1 selects mostly-plus (n.n = -1); +1 selects mostly-minus.
   DefMetric's first argument is the sign of the determinant, -1 for either. *)
RCGeometry[normu_: -1, type_String: "FLFlat"] := (
  DefManifold[M4, 4, {ia, ib, ic, id, ie, ig, ip, iq, ir, is}];
  DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
  DefMetricPerturbation[g, dg, eps];
  SetSlicing[g, n, normu, h, cd, {"|", "D"}, type];
  DefMetricFields[g, dg, h, eps];
  RCEps = $PerturbationParameter;
  RCSay["PERTURBATION_PARAMETER", RCEps];
  RCSayPlain["GEOMETRY", "normu=" <> ToString[normu] <> " type=" <> type];
  RCSigns["after_setslicing"];
);

(* RCSort: xPand's own commutation of the flat induced derivatives (private helper used by
   SplitPerturbations, xPand.m:2922) so that transversality and traceless rules can fire;
   needed after VarD, which returns derivatives in an arbitrary order. *)
RCSort[x_] := RCCanon[Block[{Print}, xAct`xPand`Private`CommuteCDSafe[NoScalar[x], cd]]];

