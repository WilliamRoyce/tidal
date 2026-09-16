#!/usr/bin/env wolframscript
(* cspell:words DefMatterFields Detg Friedmann Lchi Lfull RCEL SplitMatter SplitMetric ToxPand chib dchi inds offshell onshell tadpole xPand *)
(* repro_a1_tensor_action.wls -- A1, the HEADLINE: action -> order 1 (tadpole) -> order 2 ->
   Euler-Lagrange after the 3+1 split -> the FRW tensor equation (R-C #567, plan §5.6).
   Theory: Einstein-Hilbert + Lambda + a canonical scalar chi with potential V(chi), so the
   background is on-shell without a fluid action:
        L = sqrt(-g) [ R/(2 kap) - Lam/kap - 1/2 g^{ab} d_a chi d_b chi - V(chi) ]   (mostly plus)
   Steps: (1) split the density to second order with ToxPandFromRules; (2) ExtractOrder[.,1]
   = the tadpole; its Euler-Lagrange coefficients are the background equations, printed and
   then used as bgRules; (3) ExtractOrder[.,2] TT sector; (4) an eta-space Euler-Lagrange
   operator (VarD for the spatial derivatives, a Leibniz "dEta" for the Lie-derivative labels)
   gives the equation; a massless-scalar CONTROL validates the operator first; (5) compare with
   the A2 target up to normalization, with negative controls. *)
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

(* aux tensors standing for the derivative levels of the varied field *)
DefTensor[X0[], M4]; DefTensor[X1[], M4]; DefTensor[X2[], M4];
DefTensor[E0[-ia, -ib], M4, Symmetric[{-ia, -ib}]];
DefTensor[E1[-ia, -ib], M4, Symmetric[{-ia, -ib}]];
DefTensor[E2[-ia, -ib], M4, Symmetric[{-ia, -ib}]];

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

(* ---------------- rules: metric TT sector + the scalar field ---------------- *)
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
rulesAny = SplitMetric[g, dg, h, "AnyGauge"];
rulesTT = rulesAny /. zeroSV;
rulesChi = {chi[] -> chiB[LI[0], LI[0]], dchi[LI[o_]] :> chiB[LI[o], LI[0]]};

(* ================= CONTROL: massless scalar on FRW, TT-free (only the scalar field) ================= *)
Lchi = Sqrt[-Detg[]] (-(1/2) g[ia, ib] CD[-ia]@chi[] CD[-ib]@chi[]);
{tc, L2chi} = AbsoluteTiming[ToxPandFromRules[Lchi, Join[rulesTT /. Eth[__] -> 0, rulesChi], h, 2]];
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
{tf, L2} = AbsoluteTiming[ToxPandFromRules[Lfull, Join[rulesTT, rulesChi], h, 2]];
RCSayPlain["A1_SPLIT_WALL_S", Round[tf, 0.01]];
L0 = RCCanon[ExtractOrder[L2, 0]]; L1 = RCCanon[ExtractOrder[L2, 1]]; L2o = RCCanon[ExtractOrder[L2, 2]];
RCSay["A1_L0", L0]; RCSay["A1_L1_TADPOLE", L1]; RCSay["A1_L2_TT", L2o];
RCSayPlain["A1_L1_HAS_TT", !FreeQ[L1, Eth]];
RCSayPlain["A1_L1_HAS_SCALAR", !FreeQ[L1, chiB[LI[1], __]]];
(* the tadpole's Euler-Lagrange coefficient for the scalar perturbation = background chi equation;
   the metric-sector tadpole needs the scalar metric fields, which are switched off here (TT run),
   so the Friedmann equations are taken from the A2 field-equation route instead: *)
tad = RCEL[L1, chiB, 1, {X0, X1, X2}, {}];
RCSay["A1_TADPOLE_EL_CHI", tad];
(* background equations in the scalar-field theory (conformal time, mostly plus, kappa = 8 pi G):
   3 H^2 = kap a^2 rho,  rho = chi'^2/(2 a^2) + V + Lam/kap ;  H' = -(kap/6) a^2 (rho + 3 P),  P = chi'^2/(2a^2) - V - Lam/kap ;
   chi'' + 2 H chi' + a^2 V'(chi) = 0 *)
chi0 = chiB[LI[0], LI[0]]; chi1 = chiB[LI[0], LI[1]]; chi2 = chiB[LI[0], LI[2]];
bgRules = {V[chi0] -> 3 HH^2/(kap AA^2) - chi1^2/(2 AA^2) - Lam/kap,
           HHp -> -(kap/6) AA^2 (2 chi1^2/AA^2 - 2 (3 HH^2/(kap AA^2) - chi1^2/(2 AA^2) - Lam/kap) - 2 Lam/kap) ,
           chi2 -> -2 HH chi1 - AA^2 Derivative[1][V][chi0]};
bgRules = bgRules /. HHp -> Simplify[HHp /. bgRules];
RCSay["A1_BG_RULES_ASSUMED", bgRules];
RCSayPlain["A1_TADPOLE_EL_CHI_ONSHELL", RCCanon[tad /. bgRules]];

el = RCEL[L2o, Eth, 1, {E0, E1, E2}, {-ia, -ib}];
RCSay["A1_EL_TT_RAW", el];
elOn = RCCanon[el /. bgRules];
RCSay["A1_EL_TT_ONSHELL", elOn];
targetT = Eth[LI[1], LI[2], -ia, -ib] + 2 HH Eth[LI[1], LI[1], -ia, -ib] - lap[Eth[LI[1], LI[0], -ia, -ib]];
candsT = Flatten[Outer[Times, {1, -1, 2, -2, 1/2, -1/2, 4, -4, 1/4, -1/4}, {1, AA^2, AA^-2, AA^4}, {1, 1/kap, kap}]];
{cT, vT} = RCNormalize[elOn, targetT, candsT];
RCSayPlain["A1_TT_NORMALIZATION", ToString[cT, InputForm] <> " :: " <> vT];
RCReport["a1_tt_onshell", elOn, If[cT === None, targetT, cT targetT],
         If[cT === None, targetT /. 2 HH -> 3 HH, cT (targetT /. 2 HH -> 3 HH)]];
RCSay["A1_TT_OFFSHELL_MINUS_ONSHELL", RCCanon[el - elOn]];
RCDone["a1"];
