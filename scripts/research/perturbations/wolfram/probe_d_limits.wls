#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPand SplitMetric epsilonh dTor contortion Shapiro CHRISTOFFELCDCDT COVDS DEFCOVD DTDT KLEGACY KSTD Kleg Kstd RICCISCALAR TTTS vpart *)
(* probe_d_limits.wls -- D, fresh kernel: where xPand's split stops for a rank-3 field
   (R-C #567, §5.9 D7/D8). probe_d_torsion.wls found the rank-4 free-index split returning
   Null with InducedFrom::unknown (run 20260916T113714Z); here each rank is tried in a clean
   sequence with a pure-metric sanity split before and after, so a failure is attributed to the
   object split, not to earlier state. Then the torsionful DefCovD connection and legacy's
   ChangeCurvature route, with the same sanity splits around them. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap]; DefConstantSymbol[c1];
HH = Hh[LI[0], LI[0]]; AA = ah[LI[0], LI[0]];
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
heads[x_] := Union@Cases[x, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}];
split1[expr_] := RCTry[RCSort[ExtractOrder[ToxPandFromRules[expr, rulesT, h, 1], 1]]];
report[tag_, r_] := (RCSayPlain[tag <> "_RAN", r =!= $Failed && r =!= Null];
   RCSay[tag <> "_RESULT", r]; RCSay[tag <> "_HEADS", If[r === $Failed || r === Null, "none", heads[r]]];
   RCSay[tag <> "_MESSAGES", $MessageList]; $MessageList = {});
sanity[tag_] := Module[{r}, r = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
   RCSayPlain[tag, If[r === $Failed || r === Null, "BROKEN", "ok: " <> ToString[Length[List @@ Expand[r]]] <> " terms"]]];
sanity["DL_SANITY_METRIC_SPLIT_0_START"];
report["DL_RANK0_TT", split1[c1 Tor[ia, -ib, -ic] Tor[-ia, ib, ic]]];
report["DL_RANK0_DTDT", split1[CD[-id]@Tor[ia, -ib, -ic] CD[id]@Tor[-ia, ib, ic]]];
sanity["DL_SANITY_METRIC_SPLIT_1_AFTER_SCALARS"];
report["DL_RANK1_TRACE", split1[Tor[ia, -ia, -ib]]];
report["DL_RANK2_DIVERGENCE", split1[CD[-ia]@Tor[ia, -ib, -ic]]];
report["DL_RANK2_NN_PROJECTED_DT", split1[n[ib] n[-ia] CD[-id]@Tor[ia, -ib, -ic]]];
sanity["DL_SANITY_METRIC_SPLIT_2_AFTER_RANK2"];
report["DL_RANK3_T", split1[Tor[ia, -ib, -ic]]];
sanity["DL_SANITY_METRIC_SPLIT_3_AFTER_RANK3"];
report["DL_RANK4_DT", split1[CD[-id]@Tor[ia, -ib, -ic]]];
sanity["DL_SANITY_METRIC_SPLIT_4_AFTER_RANK4"];
(* the hand-checkable piece from the rank-2 route: n^d h^c_e of n^b n_a nabla_d T^a_bc *)
r2 = split1[n[ib] n[-ia] CD[-id]@Tor[ia, -ib, -ic]];
If[r2 =!= $Failed && r2 =!= Null,
  piece = RCTry[RCSort[n[id] h[ic, -ie] r2]];
  RCSay["DL_PIECE_TTTS_BY_HAND", piece];
  vpart = cd[-ie]@TS1[LI[1], LI[0]] + TV1[LI[1], LI[0], -ie];
  vpartP = cd[-ie]@TS1[LI[1], LI[1]] + TV1[LI[1], LI[1], -ie];
  Do[RCSayPlain["DL_PIECE_VS_HAND_FRICTION_" <> ToString[k], RCTry[Last@RCNormalize[piece, vpartP + k HH vpart, {1, -1, 2, -2, 1/2, -1/2, AA, -AA, 2 AA, -2 AA}]]], {k, {-2, -1, 0, 1, 2}}];
];
(* ---------------- the torsionful connection ---------------- *)
d8 = RCTry[(DefCovD[CDT[-ia], {"#", "DT"}, FromMetric -> g, Torsion -> True]; True)];
RCSayPlain["DL_DEFCOVD_TORSION_OK", d8 === True];
sanity["DL_SANITY_METRIC_SPLIT_5_AFTER_DEFCOVD"];
report["DL_RICCISCALAR_CDT_DIRECT", split1[RicciScalarCDT[]]];
sanity["DL_SANITY_METRIC_SPLIT_6_AFTER_CDT_SPLIT"];
(* ---- xTensor's own meaning of ChristoffelCDCDT, in a plain kernel ---- *)
DefTensor[vv[ia], M4];
chg = RCTry[ToCanonical[ChangeCovD[CDT[-ib]@vv[ia], CDT, CD] - CD[-ib]@vv[ia]]];
RCSay["DL_CDT_MINUS_CD_ON_A_VECTOR", chg];
RCSayPlain["DL_CHRISTOFFELCDCDT_IS_MINUS_K", chg === -ChristoffelCDCDT[ia, -ib, -is] vv[is] || ToCanonical[chg + ChristoffelCDCDT[ia, -ib, -is] vv[is]] === 0];
(* ---- the two candidate contortions, checked against the definition ---- *)
Kstd[a_, b_, c_] := 1/2 (TorsionCDT[a, b, c] + TorsionCDT[b, a, c] + TorsionCDT[c, a, b]);        (* K^a_bc = 1/2 (T^a_bc + T_b^a_c + T_c^a_b) *)
Kleg[a_, b_, c_] := Module[{d1 = DummyIn[Tangent[M4]], d2 = DummyIn[Tangent[M4]]},                 (* legacy _derive.py:2174-2178, verbatim slot order *)
   1/2 (TorsionCDT[a, b, c] + g[b, d1] g[a, d2] TorsionCDT[-d1, -d2, c] - g[c, d1] g[a, d2] TorsionCDT[b, -d1, -d2])];
antisym[K_] := RCCanon[K[ia, -ib, -ic] - K[ia, -ic, -ib] - TorsionCDT[ia, -ib, -ic]];
metricity[K_] := RCCanon[K[-ia, -ib, -ic] + K[-ic, -ib, -ia]];
RCSay["DL_KSTD_ANTISYM_MINUS_T", antisym[Kstd]];
RCSay["DL_KSTD_METRICITY", metricity[Kstd]];
RCSay["DL_KLEGACY_ANTISYM_MINUS_T", antisym[Kleg]];
RCSay["DL_KLEGACY_METRICITY", metricity[Kleg]];
RCSayPlain["DL_KSTD_IS_A_CONTORTION", antisym[Kstd] === 0 && metricity[Kstd] === 0];
RCSayPlain["DL_KLEGACY_IS_A_CONTORTION", antisym[Kleg] === 0 && metricity[Kleg] === 0];
RCSay["DL_KLEGACY_MINUS_KSTD", RCCanon[Kleg[ia, -ib, -ic] - Kstd[ia, -ib, -ic]]];
(* ---- the rewrite three ways: legacy as coded; the correct substitution (-Kstd); xMAG's (imported) ---- *)
lc = RCTry[ChangeCurvature[RicciScalarCDT[], CDT, CD]];
RCSay["DL_LC_HEADS", If[lc === $Failed, "$Failed", heads[lc]]];
subst[K_, sgn_] := ChristoffelCDCDT[a1_, a2_, a3_] :> sgn K[a1, a2, a3];
rewrite[K_, sgn_] := RCCanon[ContractMetric[Expand[lc /. subst[K, sgn]]]] /. TorsionCDT[x__] :> Tor[x];
lcLegacy = rewrite[Kleg, 1];
lcCorrect = rewrite[Kstd, -1];
lcStdPlus = rewrite[Kstd, 1];
RCSay["DL_REWRITE_LEGACY_AS_CODED", lcLegacy];
RCSay["DL_REWRITE_CORRECT_MINUS_KSTD", lcCorrect];
RCSay["DL_REWRITE_PLUS_KSTD", lcStdPlus];
RCSayPlain["DL_LEGACY_VS_CORRECT", RCVerdict[lcLegacy, lcCorrect]];
RCSay["DL_LEGACY_MINUS_CORRECT", RCCanon[lcLegacy - lcCorrect]];
RCSayPlain["DL_PLUS_KSTD_VS_CORRECT_QUADRATIC_TERMS_AGREE", RCCanon[(lcStdPlus - lcCorrect) /. CD[_][_] -> 0] === 0];
If[KeyExistsQ[RCArgs, "xmag"],
  xm = RCTry[Import[RCArgs["xmag"], "WXF"]];
  RCSayPlain["DL_XMAG_REWRITE_IMPORTED", xm =!= $Failed];
  If[xm =!= $Failed,
    RCSay["DL_XMAG_REWRITE", xm];
    RCSayPlain["DL_XMAG_VS_CORRECT", RCVerdict[RCCanon[xm], lcCorrect]];
    RCSayPlain["DL_XMAG_VS_LEGACY_AS_CODED", RCVerdict[RCCanon[xm], lcLegacy]];
    RCSayPlain["DL_XMAG_VS_PLUS_KSTD", RCVerdict[RCCanon[xm], lcStdPlus]];
    RCSay["DL_XMAG_MINUS_CORRECT", RCCanon[xm - lcCorrect]];
  ];
];
(* ---- the FRW consequence: split the correct rewrite and legacy's ---- *)
report["DL_CORRECT_REWRITE_SPLIT", split1[lcCorrect]];
report["DL_LEGACY_REWRITE_SPLIT", split1[lcLegacy]];
rC = split1[lcCorrect]; rL = split1[lcLegacy];
rR = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
If[rC =!= $Failed && rC =!= Null && rR =!= $Failed,
  diffR = RCTry[RCSort[rC - rR]];
  RCSay["DL_DELTA_RTILDE_MINUS_DELTA_R_CORRECT", diffR];
  RCSayPlain["DL_DIFFERENCE_IS_LINEAR_IN_TORSION_ONLY", FreeQ[diffR, Eth] && !FreeQ[diffR, TS1 | TS3 | TS4 | TS8]];
  RCSay["DL_DIFFERENCE_SYMBOLS", Union@Cases[diffR, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}]];
  If[rL =!= $Failed && rL =!= Null, RCSay["DL_LEGACY_MINUS_CORRECT_AFTER_SPLIT", RCTry[RCSort[rL - rC]]]];
];
sanity["DL_SANITY_METRIC_SPLIT_7_BEFORE_CHART"];
(* ---- the state hazard: a derivative without a metric (xCoba chart) breaks every later split ---- *)
Needs["xAct`xCoba`"];
ch = RCTry[(DefChart[Bc, M4, {0, 1, 2, 3}, {tt[], xx[], yy[], zz[]}]; True)];
RCSayPlain["DL_CHART_DEFINED", ch === True];
RCSay["DL_COVDS_AFTER_CHART", $CovDs];
RCSay["DL_METRIC_OF_EACH_COVD_AFTER_CHART", MetricOfCovD /@ $CovDs];
sanity["DL_SANITY_METRIC_SPLIT_8_AFTER_CHART"];
RCSay["DL_SANITY_AFTER_CHART_MESSAGES", $MessageList];
RCDone["dl"];
