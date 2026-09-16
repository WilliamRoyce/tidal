#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPand SetSlicing SplitMetric SplitMatter DefMatterFields normu phiB psiB Vshu Ntlc Phih Psih drho hdot tproj *)
(* probe_c_signature.wls -- C: can xPand run natively in the project's (+,-,-,-) signature?
   (R-C #567, plan §5.8). Hypothesis from the 0.4.4 source: SetSlicing[g, n, +1, ...] is
   accepted (only the determinant sign is checked, xPand.m:1715) and the split is sign-generic
   (xTensor's projector uses v v / norm), but ExtractComponents hard-codes n.n = -1 in its
   "Time" projection (xPand.m:3042) and ToxPand passes -1 as the fluid norm (xPand.m:3014).
   Tests: (1) acceptance and messages; (2) the minimal Ricci example under normu = +1 against
   the normu = -1 result under the signature map g -> -g: R -> -R; the SAME textual gauge rule
   dg = -2 n n phi - 2 h psi now carries phi_probe = -phi_phys (the n_a n_b term flips with
   g_ab while h_ab psi does not); every h^{ab}-contracted pair flips sign, so D_a D^a -> -D_a D^a
   (first run 20260916T110544Z showed exactly the Laplacian terms left over by the map without
   this last site); (3) the "Time" projection of a vector V0 n^a: ExtractComponents vs by hand;
   (4) the fluid: ToxPand's built-in split (norm hard-coded) vs SplitMatter with normu, with the
   perfect fluid and the Einstein equation transcribed to the signature
   (T = (rho+P) u u - normu P g ;  G - normu Lam g/kap = kap T). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
normu = ToExpression[Lookup[RCArgs, "normu", "1"]];
RCSayPlain["C_NORMU", normu];
RCGeometry[normu, "FLFlat"];
RCSay["C_SETSLICING_MESSAGES", $MessageList];
RCSay["C_NN", RCCanon[n[ia] n[-ia]]];
RCSay["C_H_TRACE", RCCanon[h[ia, -ia]]];
RCSay["C_PROJECTOR_CHECK_h_n_zero", RCCanon[h[-ia, -ib] n[ib]]];

(* (2) the paper's minimal example under this norm *)
DefProjectedTensor[phiB[], h]; DefProjectedTensor[psiB[], h];
MyR = ExpandPerturbation@Perturbed[Conformal[g, gah2][RicciScalarCD[]], 1];
MyGauge = {dg[LI[ord_], ia_, ib_] :> -2 n[ia] n[ib] phiB[LI[ord]] - 2 h[ia, ib] psiB[LI[ord]]};
out = SplitPerturbations[ah[]^2 MyR, MyGauge, h];
RCSay["C_OUT47_THIS_NORM", out];
out47minus = 6 Hh[LI[0], LI[0]]^2 + 6 Hh[LI[0], LI[1]] +
  RCEps (-12 Hh[LI[0], LI[0]]^2 phiB[LI[1], LI[0]] - 12 Hh[LI[0], LI[1]] phiB[LI[1], LI[0]]
       - 6 Hh[LI[0], LI[0]] phiB[LI[1], LI[1]] - 18 Hh[LI[0], LI[0]] psiB[LI[1], LI[1]]
       - 6 psiB[LI[1], LI[2]] - 2 cd[-ia]@cd[ia]@phiB[LI[1], LI[0]] + 4 cd[-ia]@cd[ia]@psiB[LI[1], LI[0]]);
flipLap[x_] := x /. cd[-i_][cd[i_][y_]] :> -cd[-i][cd[i][y]];
mapped = -(out47minus /. phiB[l__] :> -phiB[l]);   (* R -> -R ; phi_probe = -phi_phys ; psi unchanged *)
mappedLap = flipLap[mapped];                         (* + every h^{ab} contraction flips: D^2 -> -D^2 *)
RCSayPlain["C_OUT47_VS_MOSTLY_PLUS_RESULT", RCVerdict[out, out47minus]];
RCSayPlain["C_OUT47_VS_SIGNATURE_MAP_WITHOUT_LAPLACIAN_SITE", RCVerdict[out, mapped]];
RCSayPlain["C_OUT47_VS_SIGNATURE_MAP", RCVerdict[out, mappedLap]];
RCSayPlain["C_OUT47_CONTROL_MAP_WITHOUT_PHI_FLIP", RCVerdict[out, flipLap[-out47minus]]];
RCSay["C_OUT47_MAP_RESIDUAL_WITHOUT_LAPLACIAN_SITE", RCCanon[out - mapped]];
RCSay["C_LAPLACIAN_OF_PSI_TERM_COEFF", Coefficient[RCCanon[ExtractOrder[out, 1]], cd[-ia]@cd[ia]@psiB[LI[1], LI[0]]]];

(* (3) the "Time" projection sign *)
DefProjectedTensor[V0[], h];
DefTensor[vv[ia], M4];
IndexSet[vv[ia_], V0[] n[ia]];
tproj = ExtractComponents[vv[ia], h, {"Time"}];
RCSay["C_TIME_PROJECTION_OF_V0_n", RCCanon[tproj]];
RCSayPlain["C_TIME_PROJECTION_SIGN_OK", RCCanon[tproj - V0[LI[0], LI[0]]] === 0];
byHand = RCCanon[normu n[-ia] vv[ia]];               (* the projector along n is n n / (n.n) *)
RCSay["C_TIME_PROJECTION_BY_HAND", byHand];
RCSayPlain["C_TIME_PROJECTION_BY_HAND_OK", RCCanon[byHand - V0[LI[0], LI[0]]] === 0];

(* (4) the fluid *)
DefConstantSymbol[kap]; DefConstantSymbol[Lam];
DefMatterFields[u, du, h, eps];
DefTensor[Tf[-ia, -ib], M4];
IndexSet[Tf[ia_, ib_], (\[Rho]u[] + Pu[]) u[ia] u[ib] - normu Pu[] g[ia, ib]];
MyGR[ia_, ib_] := EinsteinCD[ia, ib] - normu g[ia, ib] Lam/kap - kap Tf[ia, ib];
res = RCTry[ToxPand[MyGR[ia, -ib], dg, u, du, h, "NewtonGauge", 1]];
RCSayPlain["C_TOXPAND_BUILTIN_FLUID_RAN", res =!= $Failed];
If[res =!= $Failed,
  RCSay["C_BG_TIME_TIME_BUILTIN_FLUID_EXTRACTCOMPONENTS", RCCanon[ExtractComponents[ExtractOrder[ah[]^2 res, 0], h, {"Time", "Time"}]]];
  RCSay["C_UU_FROM_BUILTIN_SPLIT", RCCanon[u[ia] u[-ia]]];
];
rulesM = RCTry[SplitMatter[u, du, normu, h, "NewtonGauge", 1]];
RCSayPlain["C_SPLITMATTER_WITH_NORMU_RAN", rulesM =!= $Failed];
If[rulesM =!= $Failed,
  res2 = ToxPandFromRules[MyGR[ia, -ib], Join[SplitMetric[g, dg, h, "NewtonGauge"], rulesM], h, 1];
  o0 = ExtractOrder[ah[]^2 res2, 0];
  bgX = RCCanon[ExtractComponents[o0, h, {"Time", "Time"}]];
  bgHand = RCCanon[normu n[-ia] n[ib] o0];             (* correct projector for this norm *)
  RCSay["C_BG_TIME_TIME_HANDMADE_FLUID_EXTRACTCOMPONENTS", bgX];
  RCSay["C_BG_TIME_TIME_HANDMADE_FLUID_BY_HAND", bgHand];
  friedmannMP = Lam ah[]^2/kap - 3 Hh[LI[0], LI[0]]^2 + kap ah[]^2 \[Rho]u[LI[0], LI[0]];   (* paper :1777-1778, mostly plus *)
  (* the mixed-index equation flips sign as a whole under g -> -g; the correct projection gives -friedmannMP, ExtractComponents (hard-coded -1) the opposite *)
  RCSayPlain["C_FRIEDMANN_BY_HAND_VS_MOSTLY_PLUS", Last@RCNormalize[bgHand, friedmannMP, {1, -1, 2, -2}]];
  RCSayPlain["C_FRIEDMANN_EXTRACTCOMPONENTS_VS_MOSTLY_PLUS", Last@RCNormalize[bgX, friedmannMP, {1, -1, 2, -2}]];
  RCSayPlain["C_FRIEDMANN_NO_SPURIOUS_PRESSURE", FreeQ[bgHand, Pu]];
  o1 = ExtractOrder[ah[]^2 res2, 1];
  tt1 = RCCanon[normu n[-ia] n[ib] o1];
  RCSay["C_ORDER1_TIME_TIME_BY_HAND", tt1];
  (* mostly-plus Newtonian 00 (repro_b, run 20260916T105232Z): kap a^2 drho + 6 H^2 Phi + 6 H Psi' - 2 D^2 Psi ; the map: whole equation -1, Phi -> -Phi, D^2 -> -D^2 *)
  mp00 = kap ah[]^2 \[Rho]u[LI[1], LI[0]] + 6 Hh[LI[0], LI[0]]^2 \[Phi]h[LI[1], LI[0]] + 6 Hh[LI[0], LI[0]] \[Psi]h[LI[1], LI[1]] - 2 cd[-ia]@cd[ia]@\[Psi]h[LI[1], LI[0]];
  RCSayPlain["C_ORDER1_TIME_TIME_VS_SIGNATURE_MAP", Last@RCNormalize[tt1, flipLap[mp00 /. \[Phi]h[l__] :> -\[Phi]h[l]], {1, -1}]];
  RCSayPlain["C_ORDER1_TIME_TIME_CONTROL_NO_PHI_FLIP", Last@RCNormalize[tt1, flipLap[mp00], {1, -1}]];
];
(* (5) CAMB's own scalar equations (camb.symbolic 2.0.4, committed printout
   camb_symbolic_newtonian_2.0.4.txt) transcribed into THIS kernel, which is CAMB's signature:
   CAMB's line element a^2((1+2 Psi_N) dt^2 - (1-2 Phi_N) dx^2) (symbolic.py:461) against the probe's
   rule dg = -2 n n phi - 2 h psi with n_0 = +a, h_ij = -a^2 delta_ij gives Psi_N = -Phih, Phi_N = Psih;
   CAMB's k^2 X is -nabla^2 X = h^{ab} D_a D_b X here (h negative definite); CAMB's delta is the
   total density perturbation, H is the conformal Hubble rate, kappa = 8 pi G. Synchronous gauge:
   eta_s = Psih and h_MB = h_ii = -2 nabla^2 E_p - 6 Psih = 2 D^2 Esh - 6 Psih (E_p = -E_phys). *)
HH = Hh[LI[0], LI[0]]; AA = ah[LI[0], LI[0]];
If[rulesM =!= $Failed,
  drho = \[Rho]u[LI[1], LI[0]]; PhiP = \[Phi]h[LI[1], LI[0]]; PsiP = \[Psi]h[LI[1], LI[0]]; PsiPp = \[Psi]h[LI[1], LI[1]];
  lapC[x_] := cd[-ic]@cd[ic]@x;
  candsC = {1, -1, 2, -2, 1/2, -1/2};
  (* Newtonian constraint [1] (K = 0): -2 k^2 Phi_N - kappa a^2 delta - 6 H (H Psi_N + Phi_N') *)
  cambN00 = -2 lapC[PsiP] - kap AA^2 drho - 6 HH (HH (-PhiP) + PsiPp);
  cambN00c = -2 lapC[PsiP] - kap AA^2 drho - 4 HH (HH (-PhiP) + PsiPp);
  RCSayPlain["C_CAMB_NEWTON_00_VS_XPAND", Last@RCNormalize[tt1, cambN00, candsC]];
  RCSayPlain["C_CAMB_NEWTON_00_CONTROL", Last@RCNormalize[tt1, cambN00c, candsC]];
  RCSayPlain["C_CAMB_NEWTON_00_CONTROL_NO_LAPSE_FLIP", Last@RCNormalize[tt1, cambN00 /. PhiP -> -PhiP, candsC]];
  (* Newtonian traceless evolution [1]: k^2 (Phi_N - Psi_N) - kappa a^2 Pi = 0, Pi = 0 for the perfect fluid -> (D_i D_j - h_ij D^2/3)(Psih + Phih) = 0 *)
  ss = RCSort[h[ic, -ia] h[ib, -id] o1];
  trS = RCSort[ss h[-ic, id]];
  tlS = RCSort[ss - (1/3) h[ic, -id] ReplaceDummies[trS]];
  RCSay["C_NEWTON_TRACELESS_BY_HAND", tlS];
  cambNtl = cd[ic]@cd[-id]@(PsiP + PhiP) - (1/3) h[ic, -id] lapC[PsiP + PhiP];
  cambNtlc = cd[ic]@cd[-id]@(PsiP - PhiP) - (1/3) h[ic, -id] lapC[PsiP - PhiP];
  RCSayPlain["C_CAMB_NEWTON_TRACELESS_VS_XPAND", Last@RCNormalize[tlS, cambNtl, candsC]];
  RCSayPlain["C_CAMB_NEWTON_TRACELESS_CONTROL", Last@RCNormalize[tlS, cambNtlc, candsC]];
  (* synchronous constraint [1]: -2 k^2 eta_s - kappa a^2 delta + H hdot_s *)
  rulesMS = RCTry[SplitMatter[u, du, normu, h, "SynchronousGauge", 1]];
  If[rulesMS =!= $Failed,
    resS = ToxPandFromRules[MyGR[ia, -ib], Join[SplitMetric[g, dg, h, "SynchronousGauge"], rulesMS], h, 1];
    o1S = ExtractOrder[ah[]^2 resS, 1];
    ttS = RCSort[normu n[-ia] n[ib] o1S];
    RCSay["C_SYNC_ORDER1_TIME_TIME_BY_HAND", ttS];
    etaS = \[Psi]h[LI[1], LI[0]]; hdotS = 2 lapC[Esh[LI[1], LI[1]]] - 6 \[Psi]h[LI[1], LI[1]];
    cambS00 = -2 lapC[etaS] - kap AA^2 drho + HH hdotS;
    cambS00c = -2 lapC[etaS] - kap AA^2 drho + 2 HH hdotS;
    RCSayPlain["C_CAMB_SYNC_00_VS_XPAND", Last@RCNormalize[ttS, cambS00, candsC]];
    RCSayPlain["C_CAMB_SYNC_00_CONTROL", Last@RCNormalize[ttS, cambS00c, candsC]];
  ];
];
RCDone["c"];
