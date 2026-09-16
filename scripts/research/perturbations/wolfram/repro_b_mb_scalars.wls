#!/usr/bin/env wolframscript
(* cspell:words Bflip CONB DefMatterFields Friedmann Phih Psih SplitMatter SplitMetric ToxPand Vshu Vvhu cona conb conc cond drho offshell onshell rhou syna synb sync synd xPand *)
(* repro_b_mb_scalars.wls -- B: Ma & Bertschinger's scalar Einstein equations from xPand, in
   both CAMB-named gauges (R-C #567, plan §5.7).
   Targets are typed from literature/astro-ph_9506072/9506072.tex, conformal Newtonian
   :644-661 (ein-cona..d) and synchronous :603-618 (ein-syna..d), with MB's Fourier equations
   turned back into real-space operator form (k^2 -> -D^2, theta = i k.v -> D.v) and their
   dictionary into xPand's variables read back by probe_load:
     MB psi (lapse)      = xPand Phih          MB phi (curvature) = xPand Psih
     MB delta rho / P    = rhou[LI[1]] / Pu[LI[1]]     MB v (v_i = d_i v) = Vshu[LI[1]]
     MB synchronous eta  = xPand Psih (synchronous)    MB h = h_ii = 2 D^2 Esh - 6 Psih
     4 pi G = kappa/2 ; dots = primes (LI[.,1]) ; sigma = 0 for a perfect fluid.
   Both sides are (-,+,+,+) here (MB :295, xPand DefMetric[-1,...]); the transcription into
   (+,-,-,-)/CAMB variables is checked separately against camb_symbolic_newtonian_2.0.4.txt.
   Each equation is compared up to an overall numeric normalization (reported) and gets a
   negative control (one coefficient altered). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap]; DefConstantSymbol[Lam];
DefMatterFields[u, du, h, eps];
DefTensor[Tf[-ia, -ib], M4];
IndexSet[Tf[ia_, ib_], (\[Rho]u[] + Pu[]) u[ia] u[ib] + Pu[] g[ia, ib]];
MyGR[ia_, ib_] := EinsteinCD[ia, ib] + g[ia, ib] Lam/kap - kap Tf[ia, ib];
zeroVT = {Bvh[__] -> 0, Evh[__] -> 0, Eth[__] -> 0, Vvhu[__] -> 0};
HH = Hh[LI[0], LI[0]]; HHp = Hh[LI[0], LI[1]]; AA = ah[LI[0], LI[0]];
rho0 = \[Rho]u[LI[0], LI[0]]; P0 = Pu[LI[0], LI[0]];
drho = \[Rho]u[LI[1], LI[0]]; dP = Pu[LI[1], LI[0]];
lap[x_] := cd[-ic]@cd[ic]@x;
cands = {1, -1, 2, -2, 1/2, -1/2, 3, -3, 1/3, -1/3, 4, -4, 6, -6};

runGauge[gauge_String] := Module[{res, bg, bgTT, bgSS, bgRules, o1, cTT, cTS, cSS, tr, tl, t, tSSraw},
  {t, res} = AbsoluteTiming[ToxPand[MyGR[ia, -ib], dg, u, du, h, gauge, 1] /. zeroVT];
  RCSayPlain["B_" <> gauge <> "_WALL_S", Round[t, 0.01]];
  bg = ExtractOrder[AA^2 res, 0];
  bgTT = RCCanon[ExtractComponents[bg, h, {"Time", "Time"}]];
  bgSS = RCCanon[ExtractComponents[bg, h, {"Space", "Space"}]] /. {h[_, _] -> 1, delta[_, _] -> 1};
  RCSay["B_" <> gauge <> "_BG_TT", bgTT];
  RCSay["B_" <> gauge <> "_BG_SS", bgSS];
  (* MB :573-580 Friedmann: H^2 = (8 pi G/3) a^2 rho (Lambda absorbed), H' = -(4 pi G/3) a^2 (rho + 3P) *)
  RCSayPlain["B_" <> gauge <> "_FRIEDMANN1", Last@RCNormalize[bgTT, HH^2 - (kap/3) AA^2 rho0 - Lam AA^2/(3 kap), cands]];
  bgRules = First@Solve[{bgTT == 0, bgSS == 0}, {rho0, HHp}];
  RCSay["B_" <> gauge <> "_BG_RULES", bgRules];
  o1 = ExtractOrder[AA^2 res, 1];
  cTT = RCCanon[ExtractComponents[o1, h, {"Time", "Time"}]];
  cTS = RCCanon[ExtractComponents[o1, h, {"Time", "Space"}]];
  cSS = RCCanon[ExtractComponents[o1, h, {"Space", "Space"}]];
  tr = RCCanon[cSS h[-ia, ib]];
  tl = RCCanon[cSS - (1/3) h[ia, -ib] ReplaceDummies[tr]];
  RCSay["B_" <> gauge <> "_SS_RAW", cSS];
  RCSayPlain["B_" <> gauge <> "_TRACELESS_TRACE_IS_ZERO", RCCanon[tl h[-ia, ib]] === 0];
  RCSay["B_" <> gauge <> "_TT", cTT]; RCSay["B_" <> gauge <> "_TS", cTS];
  RCSay["B_" <> gauge <> "_SS_TRACE", tr]; RCSay["B_" <> gauge <> "_SS_TRACELESS", tl];
  {cTT, cTS, tr, tl, bgRules}];

compare[name_String, expr_, target_, control_] := Module[{c, v},
  {c, v} = RCNormalize[expr, target, cands];
  RCSayPlain[ToUpperCase[name] <> "_NORMALIZATION", ToString[c, InputForm] <> " :: " <> v];
  RCReport[name, expr, If[c === None, target, c target], If[c === None, control, c control]]];

(* ================= conformal Newtonian gauge (MB :644-661) ================= *)
{nTT, nTS, nTr, nTl, nBg} = runGauge["NewtonGauge"];
PhiN = \[Phi]h[LI[1], LI[0]]; PsiN = \[Psi]h[LI[1], LI[0]];      (* MB psi = PhiN ; MB phi = PsiN *)
PsiNp = \[Psi]h[LI[1], LI[1]]; PsiNpp = \[Psi]h[LI[1], LI[2]]; PhiNp = \[Phi]h[LI[1], LI[1]];
V = Vshu[LI[1], LI[0]];
(* ein-cona  k^2 phi + 3H(phi' + H psi) = 4 pi G a^2 dT^0_0 = -4 pi G a^2 drho  ->  -D^2 Psi + 3H(Psi' + H Phi) + (kap/2) a^2 drho = 0 *)
tA = -lap[PsiN] + 3 HH (PsiNp + HH PhiN) + (kap/2) AA^2 drho;
tAc = -lap[PsiN] + 2 HH (PsiNp + HH PhiN) + (kap/2) AA^2 drho;
compare["mb_cona", nTT, tA, tAc];
(* ein-conb  (before the divergence)  d_i(phi' + H psi) + 4 pi G a^2 (rho + P) v_i = 0  with v_i = d_i V *)
tB = cd[-ib][PsiNp + HH PhiN] + (kap/2) AA^2 (rho0 + P0) cd[-ib][V];
tBflip = cd[-ib][PsiNp + HH PhiN] - (kap/2) AA^2 (rho0 + P0) cd[-ib][V];
compare["mb_conb", nTS, tB, tBflip];
RCSayPlain["MB_CONB_FLIPPED_VELOCITY_SIGN", Last@RCNormalize[nTS, tBflip, cands]];
(* ein-conc  phi'' + H(psi' + 2 phi') + (2H' + H^2) psi + (k^2/3)(phi - psi) = 4 pi G a^2 dP  (trace/3) *)
tC = PsiNpp + HH (PhiNp + 2 PsiNp) + (2 HHp + HH^2) PhiN - (1/3) lap[PsiN - PhiN] - (kap/2) AA^2 dP;
tCc = PsiNpp + HH (PhiNp + 2 PsiNp) + (HHp + HH^2) PhiN - (1/3) lap[PsiN - PhiN] - (kap/2) AA^2 dP;
compare["mb_conc", nTr, 3 tC, 3 tCc];
RCSayPlain["MB_CONC_ONSHELL", Last@RCNormalize[RCCanon[nTr /. nBg], RCCanon[3 tC /. nBg], cands]];
(* ein-cond  k^2 (phi - psi) = 12 pi G a^2 (rho + P) sigma = 0  ->  traceless (D^i D_j - delta/3 D^2)(Psi - Phi) = 0 *)
tD = cd[ia]@cd[-ib]@(PsiN - PhiN) - (1/3) h[ia, -ib] lap[PsiN - PhiN];
compare["mb_cond", nTl, tD, cd[ia]@cd[-ib]@(PsiN + PhiN) - (1/3) h[ia, -ib] lap[PsiN + PhiN]];

(* ================= synchronous gauge (MB :603-618) ================= *)
{sTT, sTS, sTr, sTl, sBg} = runGauge["SynchronousGauge"];
etaS = \[Psi]h[LI[1], LI[0]]; etaSp = \[Psi]h[LI[1], LI[1]]; etaSpp = \[Psi]h[LI[1], LI[2]];
Es0 = Esh[LI[1], LI[0]]; Es1 = Esh[LI[1], LI[1]]; Es2 = Esh[LI[1], LI[2]];
hMB = 2 lap[Es0] - 6 etaS; hMBp = 2 lap[Es1] - 6 etaSp; hMBpp = 2 lap[Es2] - 6 etaSpp;
RCSay["B_SYNC_DICTIONARY", {"eta" -> etaS, "h" -> hMB}];
(* ein-syna  k^2 eta - (1/2) H h' = 4 pi G a^2 dT^0_0 = -4 pi G a^2 drho  ->  -D^2 eta - (1/2) H h' + (kap/2) a^2 drho = 0 *)
sA = -lap[etaS] - (1/2) HH hMBp + (kap/2) AA^2 drho;
sAc = -lap[etaS] - HH hMBp + (kap/2) AA^2 drho;
compare["mb_syna", sTT, sA, sAc];
(* ein-synb  k^2 eta' = 4 pi G a^2 (rho + P) theta  ->  0i:  d_i eta' + 4 pi G a^2 (rho + P) v_i = 0 *)
sB = cd[-ib][etaSp] + (kap/2) AA^2 (rho0 + P0) cd[-ib][V];
compare["mb_synb", sTS, sB, cd[-ib][etaSp] - (kap/2) AA^2 (rho0 + P0) cd[-ib][V]];
(* ein-sync  h'' + 2H h' - 2 k^2 eta = -8 pi G a^2 dT^i_i = -24 pi G a^2 dP  ->  h'' + 2H h' + 2 D^2 eta + 3 kap a^2 dP = 0 *)
sC = hMBpp + 2 HH hMBp + 2 lap[etaS] + 3 kap AA^2 dP;
sCc = hMBpp + 3 HH hMBp + 2 lap[etaS] + 3 kap AA^2 dP;
compare["mb_sync", sTr, sC, sCc];
RCSayPlain["MB_SYNC_ONSHELL", Last@RCNormalize[RCCanon[sTr /. sBg], RCCanon[sC /. sBg], cands]];
(* ein-synd  h'' + 6 eta'' + 2H(h' + 6 eta') - 2 k^2 eta = 0 (sigma = 0)  ->  2 D^2 (E'' + 2H E' + eta) = 0 ; traceless part *)
sD = cd[ia]@cd[-ib]@(Es2 + 2 HH Es1 + etaS) - (1/3) h[ia, -ib] lap[Es2 + 2 HH Es1 + etaS];
compare["mb_synd", sTl, sD, cd[ia]@cd[-ib]@(Es2 + 3 HH Es1 + etaS) - (1/3) h[ia, -ib] lap[Es2 + 3 HH Es1 + etaS]];
RCDone["b"];
