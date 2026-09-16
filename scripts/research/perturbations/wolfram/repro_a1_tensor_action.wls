#!/usr/bin/env wolframscript
(* cspell:words DefMatterFields Detg Friedmann Lchi Lfull RCEL SplitMatter SplitMetric ToxPand chib dchi inds offshell onshell tadpole xPand Phih Psih NOLAP nonproj *)
(* repro_a1_tensor_action.wls -- A1, the HEADLINE: action -> order 1 (tadpole) -> order 2 ->
   Euler-Lagrange after the 3+1 split -> the FRW tensor equation (R-C #567, plan §5.6).
   Theory: Einstein-Hilbert + Lambda + a canonical scalar chi with potential V(chi), so the
   background is on-shell without a fluid action:
        L = sqrt(-g) [ R/(2 kap) - Lam/kap - 1/2 g^{ab} d_a chi d_b chi - V(chi) ]   (mostly plus)
   Steps: (1) split the density with ToxPandFromRules; (2) ExtractOrder[.,1] = the tadpole:
   its Euler-Lagrange coefficients w.r.t. the lapse potential, the curvature potential and the
   scalar perturbation ARE the two Friedmann equations and the background chi equation; they are
   printed, solved for {V, H', chi''} (bgRules) and checked against the textbook forms;
   (3) ExtractOrder[.,2] TT sector; (4) an eta-space Euler-Lagrange operator (VarD for the
   spatial derivatives, a Leibniz "dEta" for the Lie-derivative labels) gives the equation; a
   massless-scalar CONTROL validates the operator first; (5) compare with the A2 target up to
   normalization, with negative controls.
   Two facts of the machinery, found in the first attempt (run 20260916T105252Z) and handled
   here: xPand's conformal split leaves the flat conformal background's Detg[] as a symbolic
   constant inside Sqrt[-(Detg[] a^8)] (set to -1: Cartesian flat slices, so sqrt(-g) = a^4);
   and xPand's induced derivative cd refuses to act on a tensor that is not declared orthogonal
   to n (Validate::nonproj), so the auxiliary tensors that stand for the varied field are
   declared with OrthogonalTo/ProjectedWith. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap]; DefConstantSymbol[Lam];
DefMatterFields[u, du, h, eps];                      (* provides the scalar-field object CurlyPhi[LI[p],LI[q]] *)
chiB = \[CurlyPhi];                                  (* background + perturbed projected scalar: chiB[LI[p], LI[q]] *)
DefTensor[chi[], M4];                                (* the 4D field the action is written with *)
DefTensorPerturbation[dchi[LI[o]], chi[], M4];
DefScalarFunction[V];
HH = Hh[LI[0], LI[0]]; HHp = Hh[LI[0], LI[1]]; AA = ah[LI[0], LI[0]];
lap[x_] := cd[-ic]@cd[ic]@x;
(* the volume element after the split: Detg[] is the flat conformal background's determinant *)
fixVol[x_] := NoScalar[(x /. Detg[] -> -1) /. Sqrt[y_] :> PowerExpand[Sqrt[y]]];

(* ---------------- the eta-space Euler-Lagrange operator ---------------- *)
(* dEta: the Lie derivative along n on the split expression. On flat FRW in the conformal
   frame K = 0, so it commutes with cd and kills h, delta, n. Projected tensors carry the
   derivative count in their second label. *)
ClearAll[dEta];
dEta[x_Plus] := dEta /@ x;
dEta[x_Times] := Plus @@ MapIndexed[(ReplacePart[x, First[#2] -> dEta[#1]]) &, List @@ x];
dEta[x_^n_] := n x^(n - 1) dEta[x];
dEta[c_?NumericQ] := 0;
dEta[c_Symbol] /; ConstantSymbolQ[c] := 0;
dEta[cd[i_][x_]] := cd[i][dEta[x]];
dEta[(h | delta | n)[__]] := 0;
dEta[V[x_]] := Derivative[1][V][x] dEta[x];
dEta[Derivative[k_][V][x_]] := Derivative[k + 1][V][x] dEta[x];
dEta[Scalar[x_]] := dEta[x];
dEta[T_[LI[p_], LI[q_], inds___]] := T[LI[p], LI[q + 1], inds];
dEta[0] = 0;

(* aux tensors standing for the derivative levels of the varied field; spatial (orthogonal to n) *)
DefTensor[X0[], M4]; DefTensor[X1[], M4]; DefTensor[X2[], M4];
DefTensor[Y0[], M4]; DefTensor[Y1[], M4]; DefTensor[Y2[], M4];
DefTensor[Z0[], M4]; DefTensor[Z1[], M4]; DefTensor[Z2[], M4];
DefTensor[E0[-ia, -ib], M4, Symmetric[{-ia, -ib}], OrthogonalTo -> {n[ia], n[ib]}, ProjectedWith -> {h[ia, -ic], h[ib, -ic]}];
DefTensor[E1[-ia, -ib], M4, Symmetric[{-ia, -ib}], OrthogonalTo -> {n[ia], n[ib]}, ProjectedWith -> {h[ia, -ic], h[ib, -ic]}];
DefTensor[E2[-ia, -ib], M4, Symmetric[{-ia, -ib}], OrthogonalTo -> {n[ia], n[ib]}, ProjectedWith -> {h[ia, -ic], h[ib, -ic]}];
(* xPand's own bypass (xPand.m:1177) in addition, so cd never questions the aux tensors *)
E0 /: OrthogonalToVectorQ[n][E0] = True; E1 /: OrthogonalToVectorQ[n][E1] = True; E2 /: OrthogonalToVectorQ[n][E2] = True;
RCSayPlain["A1_AUX_E0_ORTHOGONAL_TO_N", OrthogonalToVectorQ[n][E0]];
RCSay["A1_AUX_E0_N_CONTRACTION", E0[-ia, -ib] n[ib]];

(* RCEL[L, field, order, auxList, indsPattern]: Euler-Lagrange derivative of the split density
   L with respect to field[LI[order], LI[k], inds], k = 0,1,2, treating spatial derivatives by
   VarD (integration by parts on the slices) and time derivatives by dEta. *)
RCEL[L_, field_, ord_, {A0_, A1_, A2_}, indsPattern_] := Module[{Lx, P0, P1, P2, back},
  Lx = NoScalar[L] /. {field[LI[ord], LI[0], ii___] :> A0[ii], field[LI[ord], LI[1], ii___] :> A1[ii], field[LI[ord], LI[2], ii___] :> A2[ii]};
  P0 = VarD[A0 @@ indsPattern, cd][Lx];
  P1 = VarD[A1 @@ indsPattern, cd][Lx];
  P2 = VarD[A2 @@ indsPattern, cd][Lx];
  back = {A0[ii___] :> field[LI[ord], LI[0], ii], A1[ii___] :> field[LI[ord], LI[1], ii], A2[ii___] :> field[LI[ord], LI[2], ii]};
  RCCanon[(P0 - dEta[RCCanon[P1 /. back]] + dEta[dEta[RCCanon[P2 /. back]]]) /. back]];

(* ---------------- rules: metric split + the scalar field ---------------- *)
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
zeroVec = {Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0, Eth[__] -> 0};
rulesAny = SplitMetric[g, dg, h, "AnyGauge"];
rulesTT = rulesAny /. zeroSV;
rulesScalarOnly = rulesAny /. zeroVec;                 (* lapse and curvature potentials only: the tadpole run *)
rulesChi = {chi[] -> chiB[LI[0], LI[0]], dchi[LI[o_]] :> chiB[LI[o], LI[0]]};

(* ================= CONTROL: massless scalar on FRW, TT-free (only the scalar field) ================= *)
Lchi = Sqrt[-Detg[]] (-(1/2) g[ia, ib] CD[-ia]@chi[] CD[-ib]@chi[]);
{tc, L2chi} = AbsoluteTiming[fixVol@ToxPandFromRules[Lchi, Join[rulesTT /. Eth[__] -> 0, rulesChi], h, 2]];
RCSayPlain["A1_CONTROL_SPLIT_WALL_S", Round[tc, 0.01]];
RCSay["A1_CONTROL_L0", RCCanon[ExtractOrder[L2chi, 0]]];
RCSay["A1_CONTROL_L1", RCCanon[ExtractOrder[L2chi, 1]]];
RCSay["A1_CONTROL_L2", RCCanon[ExtractOrder[L2chi, 2]]];
elChi = RCEL[ExtractOrder[L2chi, 2], chiB, 1, {X0, X1, X2}, {}];
RCSay["A1_CONTROL_EL", elChi];
targetChi = chiB[LI[1], LI[2]] + 2 HH chiB[LI[1], LI[1]] - lap[chiB[LI[1], LI[0]]];
candsS = Flatten[Outer[Times, {1, -1, 2, -2, 1/2, -1/2}, {1, AA^2, AA^-2, AA^4}]];
{cS, vS} = RCNormalize[elChi, targetChi, candsS];
RCSayPlain["A1_CONTROL_NORMALIZATION", ToString[cS, InputForm] <> " :: " <> vS];
RCReport["a1_control_scalar", elChi, If[cS === None, targetChi, cS targetChi],
         If[cS === None, targetChi /. 2 HH -> 3 HH, cS (targetChi /. 2 HH -> 3 HH)]];
(* the order-1 (tadpole) coefficient of the scalar sector = the background scalar equation *)
tadChi = RCEL[ExtractOrder[L2chi, 1], chiB, 1, {X0, X1, X2}, {}];
RCSay["A1_CONTROL_TADPOLE_EL", tadChi];
RCSayPlain["A1_CONTROL_TADPOLE_IS_BACKGROUND_EOM", Last@RCNormalize[tadChi, chiB[LI[0], LI[2]] + 2 HH chiB[LI[0], LI[1]], candsS]];

(* ================= the full theory: EH + Lambda + scalar with V ================= *)
Lfull = Sqrt[-Detg[]] (RicciScalarCD[]/(2 kap) - Lam/kap - (1/2) g[ia, ib] CD[-ia]@chi[] CD[-ib]@chi[] - V[chi[]]);

(* ---- step 2: the tadpole, with the scalar metric potentials switched on (order 1 only) ---- *)
{t1, L1s} = AbsoluteTiming[fixVol@ToxPandFromRules[Lfull, Join[rulesScalarOnly, rulesChi], h, 1]];
RCSayPlain["A1_TADPOLE_SPLIT_WALL_S", Round[t1, 0.01]];
L1sc = RCCanon[ExtractOrder[L1s, 1]];
RCSay["A1_TADPOLE_SCALAR_SECTOR", L1sc];
tadPhi = RCEL[L1sc, \[Phi]h, 1, {X0, X1, X2}, {}];      (* lapse potential -> Hamiltonian constraint (Friedmann 1) *)
tadPsi = RCEL[L1sc, \[Psi]h, 1, {Y0, Y1, Y2}, {}];      (* curvature potential -> Friedmann 2 *)
tadChiF = RCEL[L1sc, chiB, 1, {Z0, Z1, Z2}, {}];        (* scalar perturbation -> background chi equation *)
RCSay["A1_TADPOLE_EL_PHI", tadPhi]; RCSay["A1_TADPOLE_EL_PSI", tadPsi]; RCSay["A1_TADPOLE_EL_CHI", tadChiF];
chi0 = chiB[LI[0], LI[0]]; chi1 = chiB[LI[0], LI[1]]; chi2 = chiB[LI[0], LI[2]];
(* textbook forms (conformal time, mostly plus, kappa = 8 pi G):
   3 H^2 = kap a^2 rho,  rho = chi'^2/(2 a^2) + V + Lam/kap ;  H' = H^2 - kap chi'^2/2 ;  chi'' + 2 H chi' + a^2 V'(chi) = 0 *)
F1 = 3 HH^2 - kap AA^2 (chi1^2/(2 AA^2) + V[chi0] + Lam/kap);
F2 = HHp - HH^2 + kap chi1^2/2;
F3 = chi2 + 2 HH chi1 + AA^2 Derivative[1][V][chi0];
candsB = Flatten[Outer[Times, {1, -1, 2, -2, 1/2, -1/2, 3, -3, 1/3, -1/3, 6, -6, 1/6, -1/6}, {1, AA^2, AA^4, AA^-2}, {1, 1/kap, kap}]];
RCSayPlain["A1_TADPOLE_PHI_IS_FRIEDMANN1", Last@RCNormalize[tadPhi, F1, candsB]];
RCSayPlain["A1_TADPOLE_PSI_IS_F2_ALONE", Last@RCNormalize[tadPsi, F2, candsB]];
(* the psi variation may return a combination of F1 and F2; test that too *)
(* by hand: EL_psi = (a^2/kap) (-3 F1 - 6 F2)  (first run, 20260916T110544Z) *)
RCSayPlain["A1_TADPOLE_PSI_IS_MINUS_3F1_MINUS_6F2_TIMES_A2_OVER_KAP", Simplify[NoScalar[tadPsi] - (AA^2/kap) (-3 F1 - 6 F2)] === 0];
RCSayPlain["A1_TADPOLE_CHI_IS_BACKGROUND_EOM", Last@RCNormalize[tadChiF, F3, candsB]];
(* bgRules: solved from the tadpole itself (F1 -> V, F2 -> H', F3 -> chi''), then the typed forms as a check *)
bgFromTadpole = Quiet@Check[First@Solve[{NoScalar[tadPhi] == 0, NoScalar[tadPsi] == 0, NoScalar[tadChiF] == 0}, {V[chi0], HHp, chi2}], $Failed];
RCSay["A1_BG_RULES_FROM_TADPOLE", bgFromTadpole];
bgTyped = {V[chi0] -> 3 HH^2/(kap AA^2) - chi1^2/(2 AA^2) - Lam/kap, HHp -> HH^2 - kap chi1^2/2, chi2 -> -2 HH chi1 - AA^2 Derivative[1][V][chi0]};
RCSayPlain["A1_BG_RULES_TADPOLE_VS_TYPED_AGREE",
  If[bgFromTadpole === $Failed, "tadpole solve failed", ToString[Simplify[(Last /@ bgFromTadpole) - (Last /@ bgTyped) /. bgTyped] === {0, 0, 0}]]];
bgRules = If[bgFromTadpole === $Failed, bgTyped, bgFromTadpole];
RCSay["A1_BG_RULES_USED", bgRules];
RCSayPlain["A1_TADPOLE_ONSHELL_ALL_ZERO", Simplify[NoScalar[{tadPhi, tadPsi, tadChiF}] /. bgRules] === {0, 0, 0}];

(* ---- step 3: order 2, TT sector ---- *)
{tf, L2} = AbsoluteTiming[fixVol@ToxPandFromRules[Lfull, Join[rulesTT, rulesChi], h, 2]];
RCSayPlain["A1_SPLIT_WALL_S", Round[tf, 0.01]];
L0 = RCCanon[ExtractOrder[L2, 0]]; L1 = RCCanon[ExtractOrder[L2, 1]]; L2o = RCCanon[ExtractOrder[L2, 2]];
RCSay["A1_L0", L0]; RCSay["A1_L1_TADPOLE_TT_RUN", L1]; RCSay["A1_L2_TT", L2o];
RCSayPlain["A1_L1_HAS_TT", !FreeQ[L1, Eth]];
RCSayPlain["A1_L1_HAS_SCALAR", !FreeQ[L1, chiB[LI[1], __]]];
L2TT = RCCanon[L2o /. chiB[LI[1], __] -> 0];            (* the pure TT quadratic action *)
RCSay["A1_L2_PURE_TT", L2TT];
L2TTon = RCCanon[L2TT /. bgRules];
RCSay["A1_L2_PURE_TT_ONSHELL", L2TTon];
(* the expected quadratic TT action: a^2/(2 kap) [ E' E' - D E D E ] up to total derivatives; test the E-L instead *)

(* ---- step 4: Euler-Lagrange after the split ---- *)
elUnsorted = RCEL[L2TT, Eth, 1, {E0, E1, E2}, {-ia, -ib}];
RCSay["A1_EL_TT_RAW_UNSORTED", elUnsorted];
el = RCSort[elUnsorted];                      (* commute the flat cd's so transversality can act *)
RCSay["A1_EL_TT_RAW", el];
RCSayPlain["A1_EL_TT_SORTING_CHANGED_IT", el =!= RCCanon[elUnsorted]];
elOn = RCSort[el /. bgRules];
RCSay["A1_EL_TT_ONSHELL", elOn];
(* VarD w.r.t. E0[-ia,-ib] returns an object with upper free indices *)
targetT = Eth[LI[1], LI[2], ia, ib] + 2 HH Eth[LI[1], LI[1], ia, ib] - lap[Eth[LI[1], LI[0], ia, ib]];
candsT = Flatten[Outer[Times, {1, -1, 2, -2, 1/2, -1/2, 4, -4, 1/4, -1/4}, {1, AA^2, AA^-2, AA^4}, {1, 1/kap, kap}]];
{cT, vT} = RCNormalize[elOn, targetT, candsT];
RCSayPlain["A1_TT_NORMALIZATION", ToString[cT, InputForm] <> " :: " <> vT];
RCReport["a1_tt_onshell", elOn, If[cT === None, targetT, cT targetT],
         If[cT === None, targetT /. 2 HH -> 3 HH, cT (targetT /. 2 HH -> 3 HH)]];
RCSayPlain["A1_TT_CONTROL_NOLAP", RCVerdict[elOn, If[cT === None, 1, cT] (Eth[LI[1], LI[2], ia, ib] + 2 HH Eth[LI[1], LI[1], ia, ib])]];
resid = RCCanon[el - elOn];
RCSay["A1_TT_OFFSHELL_MINUS_ONSHELL", resid];
RCSayPlain["A1_TT_OFFSHELL_RESIDUAL_IS_ZERO", resid === 0];
RCSayPlain["A1_TT_OFFSHELL_VS_TARGET", Last@RCNormalize[el, targetT, candsT]];
RCDone["a1"];
