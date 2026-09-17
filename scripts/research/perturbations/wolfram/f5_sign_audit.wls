#!/usr/bin/env wolframscript
(* cspell:words xPand xMAG PSALTer RCSigns RCSetSigns MTW normu SetSlicing DefMetricFields ExtrinsicKSign AccelerationSign epsilonSign RICCISIGN XTENSOR xtensor RIEMANNSIGN TORSIONSIGN EPSILONSIGN *)
(* f5_sign_audit.wls -- S3: the curvature-sign globals as each package leaves them, and the
   statements docs/cosmology/conventions.md will quote. xTensor defines the three signs at
   xTensor.m:287-289 and defaults them to +1 at :1837-1843; nothing in the project states a
   curvature convention, and xMAG changes one silently at load (xMAG.m:110) and two more
   inside StartInducedDecomposition (xMAG.m:1792). This script prints the five globals in
   each state, in one kernel, in load order: plain -> xPand -> after SetSlicing -> xMAG. *)
Needs["xAct`xTensor`"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]];
RCSigns["plain_xtensor"];
RCSay["F5_XTENSOR_VERSION", xAct`xTensor`$Version];
RCSay["F5_RIEMANNSIGN_USAGE", ToString[xAct`xTensor`$RiemannSign::usage]];
RCSay["F5_RICCISIGN_USAGE", ToString[xAct`xTensor`$RicciSign::usage]];
RCSay["F5_TORSIONSIGN_USAGE", ToString[xAct`xTensor`$TorsionSign::usage]];
RCSay["F5_EPSILONSIGN_USAGE", ToString[xAct`xTensor`$epsilonSign::usage]];
Needs["xAct`xPand`"];
RCSigns["after_xpand"];
RCSay["F5_XPAND_VERSION", xAct`xPand`$Version];
DefManifold[M4, 4, {ia, ib, ic, id, ie, ig, ip, iq, ir, is}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
DefMetricPerturbation[g, dg, eps];
SetSlicing[g, n, -1, h, cd, {"|", "D"}, "FLFlat"];
DefMetricFields[g, dg, h, eps];
RCSigns["after_xpand_setslicing"];
(* the sign xPand uses for the extrinsic curvature, for the record *)
RCSay["F5_XPAND_EXTRINSIC_K_SIGN_IN_USE", RCTry[$ExtrinsicKSign]];
Needs["xAct`xMAG`"];
RCSigns["after_xmag_load"];
RCSay["F5_XMAG_VERSION", xAct`xMAG`$Version];
RCSetSigns[{1, 1, 1}];
RCSigns["after_project_convention_asserted"];
RCDone["f5"];
