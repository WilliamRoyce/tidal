#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPand SplitMetric xMAG xBrauer TraceFree xTras Helpin ToDistortion BreakDistortion ToContorsion BreakContorsion DefConnectionPerturbation StartInducedDecomposition contortion epsilonh dTor Shapiro Hehl Riemann-Cartan CDCDL CDLCD CDTCD CHANGECOVD CHRISTOFFELCDCDL CHRISTOFFELCDLCD COMMUTECOVDS COVDS DEFCOVD DGAM Kleg NONMETRICITYQ TORSIONQ XTENSOR DEFCONNECTIONPERTURBATION TODISTORTION BREAKCONTORSION *)
(* probe_e_xmag.wls -- E: xMAG (Helpin, 0.1.0) beside xPand on the certified bundle (R-C #567,
   plan §5.10, user decision Q4). Runs AFTER install_xmag.sh added the five directories.
   (a) load: versions, expected minimums, messages; the Q4 stop condition is a load failure.
   (b) the Riemann-Cartan connection of the project's field content (metric + torsion):
       DefCovD[CDT, Torsion -> True, FromMetric -> g, Master -> g] as in xMAG's
       Documentation/DefCovD_xMAG.nb; the tensors it defines; xMAG's contortion vs the
       identity legacy TIDAL types by hand (Shapiro 2002 eq. 2.9; _derive.py:2164-2189).
   (c) ToDistortion + BreakDistortion of the Riemann-Cartan Ricci scalar: the post-Riemannian
       rewrite as a package call, compared with legacy's ChangeCurvature route (harness verdict).
   (d) DefConnectionPerturbation beside DefMetricPerturbation: does xPert accept both?
   (e) end-to-end: the rewritten Ricci scalar, with the torsion identified with Probe D's field,
       through xPand's split with Probe D's SVT rules; compared with the legacy-route split.
   (f) StartInducedDecomposition: one call, its usage, what it defines (an independent 3+1
       route or not); it documents itself as built for FromMetric -> Null connections. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap];
reached[k_] := RCSayPlain["E_REACHED", k];
heads[x_] := Union@Cases[x, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}];

(* ---------------- (a) load ---------------- *)
$EFailed = False; t0 = AbsoluteTime[];
Check[Needs["xAct`xMAG`"], $EFailed = True];
RCSayPlain["E_XMAG_LOAD_FAILED", $EFailed];
RCSayPlain["E_XMAG_LOAD_WALL_S", Round[AbsoluteTime[] - t0, 0.1]];
RCSay["E_XMAG_LOAD_MESSAGES", $MessageList];
RCSay["E_XMAG_VERSION", xAct`xMAG`$Version];
RCSay["E_XMAG_EXPECTS_XTENSOR_XTRAS", {xAct`xMAG`$xTensorVersionExpected, xAct`xMAG`$xTrasVersionExpected}];
RCSay["E_XBRAUER_VERSION", xAct`xBrauer`$Version];
RCSay["E_TRACEFREE_VERSION", xAct`TraceFree`$Version];
RCSay["E_XTRAS_VERSION", xAct`xTras`$Version];
RCSay["E_PACKAGES_LOADED", Select[$Packages, StringStartsQ[#, "xAct`"] || StringStartsQ[#, "SymmetricFunctions"] || StringStartsQ[#, "BrauerAlgebra"] &]];
RCSay["E_DEFCOVD_OPTIONS_AFTER_XMAG", First /@ Options[DefCovD]];
If[$EFailed, RCSayPlain["E_STOP", "xMAG did not load: Q4 stop condition"]; RCDone["e"]; Quit[]];
(* xPand still usable after xMAG/xBrauer changed DefTensor/DefCovD options? *)
p1 = RCTry[(DefProjectedTensor[qq[], h]; qq[LI[1], LI[0]])];
RCSayPlain["E_XPAND_DEFPROJECTEDTENSOR_STILL_WORKS", p1 =!= $Failed];
(* does xPand's SPLIT still work right after loading xMAG, before any xMAG connection exists? *)
zeroSV0 = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
sanity0 = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV0, h, 1], 1]]];
RCSayPlain["E_XPAND_SPLIT_AFTER_LOADING_XMAG_BEFORE_ANY_CONNECTION", If[sanity0 === $Failed || sanity0 === Null, "BROKEN", "ok"]];
RCSay["E_XPAND_SPLIT_AFTER_LOADING_MESSAGES", $MessageList]; $MessageList = {};
reached["a"];

(* ---------------- (b) the Riemann-Cartan connection ---------------- *)
b = RCTry[(DefCovD[CDT[-ia], {"#", "DT"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True]; True)];
RCSayPlain["E_DEFCOVD_RC_OK", b === True];
sanity1 = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV0, h, 1], 1]]];
RCSayPlain["E_XPAND_SPLIT_AFTER_XMAG_DEFCOVD", If[sanity1 === $Failed || sanity1 === Null, "BROKEN", "ok"]];
RCSay["E_XPAND_SPLIT_AFTER_XMAG_DEFCOVD_MESSAGES", $MessageList]; $MessageList = {};
RCSay["E_DEFCOVD_MESSAGES", $MessageList];
RCSay["E_TENSORS_NAMED_CDT", Select[$Tensors, StringEndsQ[ToString[#], "CDT"] &]];
RCSayPlain["E_TORSIONQ_NONMETRICITYQ", ToString[{TorsionQ[CDT], NonMetricityQ[CDT]}]];
kx = RCTry[ContorsionCDT[ia, -ib, -ic]];
RCSay["E_CONTORSION_CDT_RAW", kx];
kb = RCTry[BreakContorsion[ContorsionCDT[ia, -ib, -ic], CDT]];
RCSay["E_CONTORSION_BROKEN_TO_TORSION", kb];
(* legacy's identity K^a_bc = 1/2 (T^a_bc + T_b^a_c - T_cb^a) on xTensor's TorsionCDT *)
Kleg = 1/2 (TorsionCDT[ia, -ib, -ic] + TorsionCDT[-ib, ia, -ic] - TorsionCDT[-ic, -ib, ia]);
If[kb =!= $Failed, RCSayPlain["E_CONTORSION_VS_LEGACY_IDENTITY", Last@RCNormalize[RCCanon[kb], RCCanon[Kleg], {1, -1, 2, -2, 1/2, -1/2}]]];
(* xTensor's own relation between the two connections, for the record *)
(* xTensor defined ChristoffelCDCDT (the difference of the two connections), not ChristoffelCDTCD (first run, 20260916T112331Z) *)
chr = RCTry[ChristoffelCDCDT[ia, -ib, -ic]];
RCSay["E_CHRISTOFFEL_CDCDT_EVALUATES_TO", chr];
If[chr =!= $Failed && kb =!= $Failed, RCSayPlain["E_CHRISTOFFEL_CDCDT_VS_BROKEN_CONTORSION", RCTry[Last@RCNormalize[RCCanon[BreakContorsion[chr, CDT]], RCCanon[kb], {1, -1}]]]];

reached["b"];

(* ---------------- (c) ToDistortion of the Riemann-Cartan Ricci scalar ---------------- *)
{tc, rt} = AbsoluteTiming[RCTry[ToDistortion[RicciScalarCDT[], CDT, g]]];
RCSayPlain["E_TODISTORTION_RAN", rt =!= $Failed];
RCSayPlain["E_TODISTORTION_WALL_S", Round[tc, 0.01]];
RCSay["E_TODISTORTION_HEADS", If[rt === $Failed, "$Failed", heads[rt]]];
RCSay["E_TODISTORTION_RESULT", rt];
(* BreakDistortion returned Null for this metric-compatible connection (run 20260916T112932Z); the
   documented call for torsion-only connections is BreakContorsion *)
rtD = If[rt === $Failed, $Failed, RCTry[BreakDistortion[rt, CDT, g]]];
RCSay["E_BREAKDISTORTION_RESULT_ON_RC_CONNECTION", rtD];
rtB = If[rt === $Failed, $Failed, RCTry[BreakContorsion[rt, CDT]]];
RCSayPlain["E_BREAKCONTORSION_RAN", rtB =!= $Failed && rtB =!= Null];
sanityFn[tag_] := Module[{r}, r = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV0, h, 1], 1]]];
  RCSayPlain[tag, If[r === $Failed || r === Null, "BROKEN", "ok"]]; RCSay[tag <> "_MESSAGES", $MessageList]; $MessageList = {}];
sanityFn["E_XPAND_SPLIT_AFTER_TODISTORTION_BREAKCONTORSION"];
lcPlus = lcMinus = lcMatch = $Failed;
If[rtB =!= $Failed && rtB =!= Null,
  rtB = RCCanon[ContractMetric[NoScalar[rtB]]];
  RCSay["E_RICCI_SCALAR_RC_AS_LC_PLUS_TORSION", rtB];
  RCSay["E_RICCI_SCALAR_RC_HEADS", heads[rtB]];
  RCSayPlain["E_RICCI_SCALAR_RC_FREE_OF_CDT_OBJECTS", FreeQ[rtB, RicciScalarCDT | RicciCDT | RiemannCDT | ChristoffelCDCDT | ChristoffelCDTCD | DistortionCDT | ContorsionCDT | CDT]];
  RCSayPlain["E_RICCI_SCALAR_RC_HAS_TORSION_AND_LC_RICCI", !FreeQ[rtB, TorsionCDT] && !FreeQ[rtB, RicciScalarCD]];
  (* legacy's route on a PLAIN xTensor connection CDL (no Master), as legacy declares it:
     ChangeCurvature, then the typed contortion identity on ChristoffelCDCDL -- with legacy's
     sign (+K) and the opposite one (-K); xMAG's own relation is ChristoffelCDCDT = -ContorsionCDT
     (printed above). CDL is defined only now: with a second torsion connection present,
     ToDistortion on CDT returned Null with Validate::inhom (run 20260916T113714Z) -- retested below. *)
  bl = RCTry[(DefCovD[CDL[-ia], {"&", "DL"}, FromMetric -> g, Torsion -> True]; True)];
  RCSayPlain["E_PLAIN_XTENSOR_TORSION_CONNECTION_OK", bl === True];
  sanityFn["E_XPAND_SPLIT_AFTER_SECOND_CONNECTION_CDL"];
  rtAgain = RCTry[ToDistortion[RicciScalarCDT[], CDT, g]];
  RCSayPlain["E_TODISTORTION_AFTER_SECOND_TORSION_CONNECTION", If[rtAgain === $Failed || rtAgain === Null, "returns Null/fails (one torsion connection per kernel)", "still works"]];
  RCSay["E_TODISTORTION_AFTER_SECOND_CONNECTION_MESSAGES", $MessageList];
  (* xTensor's OWN convention for the connection-difference tensor, tested on CDL without xMAG:
     (i) [DL_b, DL_c] f via ChangeCovD gives the antisymmetric part of ChristoffelCDLCD or ChristoffelCDCDL;
     (ii) the same commutator via CommuteCovDs gives xTensor's torsion sign;
     legacy's typed K must satisfy K^a_bc - K^a_cb = (the sign xTensor uses) T^a_bc for the symbol it replaces. *)
  If[bl === True,
    DefTensor[ff[], M4];
    commA = RCTry[ToCanonical[ChangeCovD[CDL[-ib]@CDL[-ic]@ff[] - CDL[-ic]@CDL[-ib]@ff[], CDL, CD]]];
    RCSay["E_XTENSOR_COMMUTATOR_VIA_CHANGECOVD", commA];
    commB = RCTry[ToCanonical[CommuteCovDs[CDL[-ib]@CDL[-ic]@ff[], CDL, {-ib, -ic}] - CDL[-ic]@CDL[-ib]@ff[]]];
    RCSay["E_XTENSOR_COMMUTATOR_VIA_COMMUTECOVDS", commB];
    DefTensor[ff2[ia], M4];
    chrDef = RCTry[ToCanonical[ChangeCovD[CDL[-ib]@ff2[ia], CDL, CD] - CD[-ib]@ff2[ia]]];
    RCSay["E_XTENSOR_CHANGECOVD_OF_A_VECTOR_MINUS_LC", If[chrDef === $Failed, "$Failed", chrDef]];
    KlegL = 1/2 (TorsionCDL[ia, -ib, -ic] + TorsionCDL[-ib, ia, -ic] - TorsionCDL[-ic, -ib, ia]);
    antiK = RCCanon[KlegL - ReplaceIndex[KlegL, {-ib -> -ic, -ic -> -ib}]];
    RCSay["E_LEGACY_K_ANTISYMMETRIC_PART", antiK];
    RCSayPlain["E_LEGACY_K_ANTISYM_EQUALS_PLUS_T", RCCanon[antiK - TorsionCDL[ia, -ib, -ic]] === 0];
    (* what does xTensor say the antisymmetric part of ChristoffelCDCDL is? *)
    antiChr = RCTry[ToCanonical[ChristoffelCDCDL[ia, -ib, -ic] - ChristoffelCDCDL[ia, -ic, -ib]]];
    RCSay["E_XTENSOR_ANTISYM_CHRISTOFFELCDCDL", antiChr];
    antiChr2 = RCTry[ToCanonical[ChristoffelCDLCD[ia, -ib, -ic] - ChristoffelCDLCD[ia, -ic, -ib]]];
    RCSay["E_XTENSOR_ANTISYM_CHRISTOFFELCDLCD", antiChr2];
  ];
  If[bl === True,
    lc0 = RCTry[ChangeCurvature[RicciScalarCDL[], CDL, CD]];
    RCSay["E_LEGACY_ROUTE_HEADS_BEFORE_RULE", If[lc0 === $Failed, "$Failed", heads[lc0]]];
    If[lc0 =!= $Failed,
      ident[sgn_] := ChristoffelCDCDL[a1_, a2_, a3_] :> With[{d1 = DummyIn[Tangent[M4]], d2 = DummyIn[Tangent[M4]]},
          sgn/2 (TorsionCDL[a1, a2, a3] + g[a2, d1] g[a1, d2] TorsionCDL[-d1, -d2, a3] - g[a3, d1] g[a1, d2] TorsionCDL[a2, -d1, -d2])];
      toT = TorsionCDL[x__] :> TorsionCDT[x];
      lcPlus = RCCanon[ContractMetric[Expand[lc0 /. ident[1]]] /. toT];
      lcMinus = RCCanon[ContractMetric[Expand[lc0 /. ident[-1]]] /. toT];
      RCSay["E_LEGACY_ROUTE_PLUS_K", lcPlus];
      RCSay["E_LEGACY_ROUTE_MINUS_K", lcMinus];
      RCSayPlain["E_XMAG_VS_LEGACY_PLUS_K", RCTry[RCVerdict[rtB, lcPlus]]];
      RCSayPlain["E_XMAG_VS_LEGACY_MINUS_K", RCTry[RCVerdict[rtB, lcMinus]]];
      RCSay["E_XMAG_MINUS_LEGACY_PLUS_K", RCTry[RCCanon[rtB - lcPlus]]];
      RCSay["E_XMAG_MINUS_LEGACY_MINUS_K", RCTry[RCCanon[rtB - lcMinus]]];
      RCSayPlain["E_XMAG_VS_LEGACY_CONTROL", RCTry[RCVerdict[rtB, lcPlus + TorsionCDT[ia, -ia, -ib] TorsionCDT[ic, -ic, ib]]]];
      lcMatch = Which[StringStartsQ[RCVerdict[rtB, lcPlus], "identical" | "proved-equal"], lcPlus,
                      StringStartsQ[RCVerdict[rtB, lcMinus], "identical" | "proved-equal"], lcMinus, True, $Failed];
      RCSayPlain["E_LEGACY_SIGN_THAT_MATCHES_XMAG", Which[lcMatch === lcPlus, "+K (legacy's typed sign)", lcMatch === lcMinus, "-K (opposite to legacy's typed sign)", True, "neither"]];
    ];
  ];
];
reached["c"];

(* ---------------- (d) a connection perturbation beside the metric perturbation ---------------- *)
RCSay["E_CONNECTION_PERTURBATION_TENSORS_ALREADY_DEFINED", Select[$Tensors, StringContainsQ[ToString[#], "Christoffel"] || StringContainsQ[ToString[#], "Perturbation"] &]];
RCSayPlain["E_DEFCONNECTIONPERTURBATION_USAGE", StringTake[ToString[DefConnectionPerturbation::usage], UpTo[400]]];
d = RCTry[(DefConnectionPerturbation[ChristoffelCDT, dGam, eps]; True)];
RCSayPlain["E_DEFCONNECTIONPERTURBATION_OK", d === True];
sanityFn["E_XPAND_SPLIT_AFTER_DEFCONNECTIONPERTURBATION"];
RCSay["E_DEFCONNECTIONPERTURBATION_MESSAGES", $MessageList];
If[d === True,
  RCSay["E_DGAM_EXAMPLE", dGam[LI[1], ia, -ib, -ic]];
  pt = RCTry[ExpandPerturbation@Perturbation[TorsionCDT[ia, -ib, -ic]]];
  RCSay["E_PERTURBATION_OF_TORSIONCDT", pt];
  pg = RCTry[ExpandPerturbation@Perturbation[g[-ia, -ib]]];
  RCSay["E_PERTURBATION_OF_METRIC_STILL_DG", pg];
  RCSayPlain["E_BOTH_PERTURBATION_STRUCTURES_COEXIST", pt =!= $Failed && pg =!= $Failed && !FreeQ[pg, dg]];
];
reached["d"];

(* ---------------- (e) end-to-end with xPand and Probe D's rules ---------------- *)
If[rtB =!= $Failed && rtB =!= Null,
  DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
  DefTensorPerturbation[dTor[LI[o], ia, -ib, -ic], Tor[ia, -ib, -ic], M4];
  Do[DefProjectedTensor[Symbol["TS" <> ToString[k]][], h], {k, 8}];
  Do[DefProjectedTensor[Symbol["TV" <> ToString[k]][-ia], h, TensorProperties -> {"Transverse"}], {k, 6}];
  Do[DefProjectedTensor[Symbol["TT" <> ToString[k]][-ia, -ib], h, TensorProperties -> {"SymmetricTensor", "Traceless", "Transverse"}], {k, 2}];
  Wp[o_, a_, c_, d_, e_] := h[a, c] TS3[LI[o]] + cd[a]@cd[c]@TS4[LI[o]] + (cd[a]@TV3[LI[o], c] + cd[c]@TV3[LI[o], a])/2 + TT1[LI[o], a, c] + epsilonh[a, c, -d] (cd[d]@TS5[LI[o]] + TV4[LI[o], d]);
  Up[o_, a_, c_, d_, e_] := h[a, c] TS6[LI[o]] + cd[a]@cd[c]@TS7[LI[o]] + (cd[a]@TV5[LI[o], c] + cd[c]@TV5[LI[o], a])/2 + TT2[LI[o], a, c] + epsilonh[a, c, -e] (cd[e]@TS8[LI[o]] + TV6[LI[o], e]);
  torsionRule = dTor[LI[o_], ia_, ib_, ic_] :> With[{id = DummyIn[Tangent[M4]], ie = DummyIn[Tangent[M4]], ig = DummyIn[Tangent[M4]]},
      n[ia] (n[ib] (cd[ic]@TS1[LI[o]] + TV1[LI[o], ic]) - n[ic] (cd[ib]@TS1[LI[o]] + TV1[LI[o], ib]))
    + n[ia] epsilonh[ib, ic, -id] (cd[id]@TS2[LI[o]] + TV2[LI[o], id])
    + (n[ib] Wp[o, ia, ic, ie, ig] - n[ic] Wp[o, ia, ib, ie, ig])
    + epsilonh[ib, ic, -id] Up[o, ia, id, ie, ig]];
  zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
  rulesMetric = SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV;
  rulesT = Join[rulesMetric, {torsionRule}, {Tor[ia_, ib_, ic_] -> 0}];
  toTor = TorsionCDT[a1_, a2_, a3_] :> Tor[a1, a2, a3];
  rtT = rtB /. toTor;
  Export[FileNameJoin[{RCOutDir, "xmag_rtilde_in_tor.wxf"}], rtT, "WXF"];
  Export[FileNameJoin[{RCOutDir, "xmag_rtilde_in_tor.txt"}], ToString[rtT, InputForm], "Text"];
  RCSayPlain["E_RTILDE_EXPORTED", "xmag_rtilde_in_tor.wxf (split in a clean kernel by probe_d_limits.wls xmag=<path>)"];
  RCSay["E_COVDS_BEFORE_SPLIT", $CovDs];
  RCSay["E_METRIC_OF_EACH_COVD", MetricOfCovD /@ $CovDs];
  sanityE = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
  RCSayPlain["E_PURE_METRIC_SPLIT_IN_THIS_KERNEL", If[sanityE === $Failed || sanityE === Null, "BROKEN (InducedFrom::unknown)", "ok"]];
  RCSay["E_PURE_METRIC_SPLIT_MESSAGES", $MessageList];
  {te, e1} = AbsoluteTiming[RCTry[RCSort[ExtractOrder[ToxPandFromRules[rtT, rulesT, h, 1], 1]]]];
  RCSayPlain["E_END_TO_END_SPLIT_RAN", e1 =!= $Failed];
  RCSayPlain["E_END_TO_END_WALL_S", Round[te, 0.01]];
  If[e1 =!= $Failed && e1 =!= Null,
    RCSay["E_END_TO_END_ORDER1", e1];
    RCSayPlain["E_END_TO_END_NO_RESIDUAL_CD_OR_DTOR", FreeQ[e1, CD] && FreeQ[e1, dTor]];
    If[lcMatch =!= $Failed,
      e2 = RCTry[RCSort[ExtractOrder[ToxPandFromRules[lcMatch /. toTor, rulesT, h, 1], 1]]];
      RCSayPlain["E_END_TO_END_XMAG_VS_LEGACY_ROUTE", If[e2 === $Failed, "legacy route split failed", RCTry[RCVerdict[e1, e2]]]];
    ];
    rR = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
    If[e1 =!= Null && rR =!= $Failed && rR =!= Null, RCSay["E_END_TO_END_TORSION_PART", RCTry[RCSort[e1 - rR]]]];
  ];
  (* the second-order torsion sector of R-tilde through xPand: the piece M3 would vary *)
  {te2, e3} = AbsoluteTiming[RCTry[RCSort[ExtractOrder[ToxPandFromRules[Sqrt[-Detg[]] rtT/(2 kap), rulesT, h, 2], 2] /. Detg[] -> -1 /. Sqrt[y_] :> PowerExpand[Sqrt[y]]]]];
  RCSayPlain["E_ORDER2_SPLIT_RAN", e3 =!= $Failed];
  RCSayPlain["E_ORDER2_WALL_S", Round[te2, 0.01]];
  If[e3 =!= $Failed,
    RCSayPlain["E_ORDER2_TERM_COUNT", Length[List @@ Expand[NoScalar[e3]]]];
    RCSay["E_ORDER2_SYMBOLS", Union@Cases[e3, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}]];
    ex3 = Expand[NoScalar[e3]];
    RCSayPlain["E_ORDER2_HAS_METRIC_TORSION_MIXING", !FreeQ[ex3, Eth] && !FreeQ[ex3, TS1 | TS2 | TS3 | TS4 | TS5 | TS6 | TS7 | TS8 | TV1 | TV2 | TV3 | TV4 | TV5 | TV6 | TT1 | TT2] && Length[Cases[ex3, t_Times /; !FreeQ[t, Eth] && !FreeQ[t, TS1 | TS3 | TT1 | TT2 | TV1 | TV3], {0, 1}]] > 0];
  ];
];
reached["e"];

(* ---------------- (f) StartInducedDecomposition ---------------- *)
RCSayPlain["E_STARTINDUCED_USAGE", StringTake[ToString[StartInducedDecomposition::usage], UpTo[500]]];
RCSay["E_STARTINDUCED_OPTIONS", Options[StartInducedDecomposition]];
f1 = RCTry[(StartInducedDecomposition[g, CDT, {{"|h", "Dh"}, {"#h", "DTh"}}, {nv, hh}]; True)];
RCSayPlain["E_STARTINDUCED_ON_RIEMANN_CARTAN_OK", f1 === True];
RCSay["E_STARTINDUCED_MESSAGES", $MessageList];
If[f1 =!= True,
  f2 = RCTry[(DefCovD[CDG[-ia], {"%", "DG"}, Torsion -> True, FromMetric -> Null, Master -> g, ConnectionRelations -> True]; True)];
  RCSayPlain["E_GENERAL_MAG_CONNECTION_OK", f2 === True];
  If[f2 === True,
    f3 = RCTry[(StartInducedDecomposition[g, CDG, {{"|h", "Dh"}, {"%h", "DGh"}}, {nv, hh}]; True)];
    RCSayPlain["E_STARTINDUCED_ON_GENERAL_CONNECTION_OK", f3 === True];
    RCSay["E_STARTINDUCED_MESSAGES_GENERAL", $MessageList];
  ];
];
RCSay["E_TENSORS_DEFINED_BY_STARTINDUCED", Select[$Tensors, StringContainsQ[ToString[#], "hh"] || StringContainsQ[ToString[#], "Extrinsic"] || StringContainsQ[ToString[#], "Acceleration"] &]];
RCSay["E_COVDS_NOW", $CovDs];
reached["f"];
RCDone["e"];
