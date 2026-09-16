#!/usr/bin/env wolframscript
(* cspell:words xPand xMAG StartInducedDecomposition contortion COVDS *)
(* probe_e2_xmag_induced.wls -- E (f) in isolation: xMAG's StartInducedDecomposition on a
   kernel that has xPand loaded but NO xPand slicing (RCGeometry not called), so that a failure
   can be attributed to xMAG itself or to the coexistence with xPand's n, h, cd (R-C #567, §5.10). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
Check[Needs["xAct`xMAG`"], Print["RC_E2_XMAG_LOAD_FAILED=True"]; RCDone["e2"]; Quit[]];
DefManifold[M4, 4, {ia, ib, ic, id, ie, ig, ip, iq, ir, is}];
DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
(* the tutorial's general connection, then its own induced decomposition, its own names *)
f2 = RCTry[(DefCovD[CDG[-ia], {"%", "DG"}, Torsion -> True, FromMetric -> Null, Master -> g, ConnectionRelations -> True]; True)];
RCSayPlain["E2_GENERAL_MAG_CONNECTION_OK", f2 === True];
f3 = RCTry[(StartInducedDecomposition[g, CDG, {{";h", "Dh"}, {"%h", "DGh"}}, {nn, hh}]; True)];
RCSayPlain["E2_STARTINDUCED_ON_GENERAL_CONNECTION_OK", f3 === True];
RCSay["E2_STARTINDUCED_MESSAGES", $MessageList];
RCSay["E2_TENSORS_DEFINED", Select[$Tensors, StringContainsQ[ToString[#], "hh"] || StringContainsQ[ToString[#], "Extrinsic"] || StringContainsQ[ToString[#], "Acceleration"] || StringContainsQ[ToString[#], "nn"] &]];
RCSay["E2_COVDS_NOW", $CovDs];
If[f3 === True,
  RCSay["E2_NN_NORM", ToCanonical[nn[ia] nn[-ia]]];
  RCSay["E2_HH_TRACE", ToCanonical[hh[ia, -ia]]];
  RCSayPlain["E2_HAS_NORMALIZATION_OPTION", MemberQ[First /@ Options[StartInducedDecomposition], Normalization]];
];
(* the Riemann-Cartan (metric-compatible) case *)
f4 = RCTry[(DefCovD[CDT[-ia], {"#", "DT"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True]; True)];
f5 = RCTry[(StartInducedDecomposition[g, CDT, {{";k", "Dk"}, {"#k", "DTk"}}, {nk, hk}]; True)];
RCSayPlain["E2_STARTINDUCED_ON_RIEMANN_CARTAN_OK", f5 === True];
RCSay["E2_STARTINDUCED_RC_MESSAGES", $MessageList];
RCDone["e2"];
