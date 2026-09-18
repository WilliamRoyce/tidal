#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPand SplitMetric xCoba DefChart epsilonh contortion Aoki Bahamonde Gigante Valcarcel Gorji Nikiforova Damour onshell offshell dTor RCEL Lkin Shapiro Hehl CDTCD CHRISTOFFELCDCDT DEFCOVD RICCISCALARCDT TTTS inds projs vpart *)
(* probe_d_torsion.wls -- D: torsion as a perturbation on FRW through xPand (R-C #567, plan §5.9).
   D1 declare a rank-3 torsion field and its perturbation beside dg (legacy makes the same pair,
      tidal/cli/_derive.py:2350-2376);
   D2 24-component count by xCoba (and the combinatorial count as a fallback);
   D3 the SVT parametrization transcribed from Aoki-Bahamonde-Gigante Valcarcel-Gorji
      (arXiv:2310.16007 TeX :768-776 helicity 2, :798-807 helicity 1, :834-842 helicity 0;
      table :864-878): 8 scalars, 6 transverse vectors, 2 TT tensors = 24, cross-checked with
      Nikiforova-Damour (1804.09215 TeX :819-852); written covariantly on the slices
      (N -> 1, gamma -> h, D -> cd, epsilon_ijk -> epsilonh), a-powers absorbed by the
      conformal frame;
   D4 projection of the rank-3 field by ExtractComponents (supported by design, xPand.m:3031);
   D6 the ACTION route: (a) an algebraic torsion-quadratic Lagrangian, (b) with a kinetic term,
      split to second order; the tadpole checked to vanish for T-bar = 0; the TT sector's
      Euler-Lagrange equations for the two TT pieces (decoupled? diagonal?);
   D7 the derivative split (no residual CD; primes and D only), one hand-checkable piece;
   D8 the torsionful connection DefCovD[..., Torsion -> True] (legacy's declaration): what
      xPand does with RicciScalarCDT, then legacy's ChangeCurvature + contortion identity;
   Every xPand call goes through RCTry so a Throw is reported, not fatal (the first run,
   20260916T110641Z, died at D6 on an uncaught Throw inside Quiet@Check).
   Verdicts: representable / representable with work (named) / not representable. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap]; DefConstantSymbol[c1]; DefConstantSymbol[c2]; DefConstantSymbol[c3]; DefConstantSymbol[c4];
HH = Hh[LI[0], LI[0]]; AA = ah[LI[0], LI[0]];
lap[x_] := cd[-ic]@cd[ic]@x;
reached[k_] := RCSayPlain["D_REACHED", k];
fixVol[x_] := NoScalar[(x /. Detg[] -> -1) /. Sqrt[y_] :> PowerExpand[Sqrt[y]]];

(* ---------------- D1: the field and its perturbation ---------------- *)
d1 = RCTry[(DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
            DefTensorPerturbation[dTor[LI[o], ia, -ib, -ic], Tor[ia, -ib, -ic], M4]; True)];
RCSayPlain["D1_DECLARED_BESIDE_DG", d1 === True];
RCSay["D1_MESSAGES", $MessageList];
reached["D1"];

(* D2 (the xCoba component count) runs LAST: DefChart defines a basis derivative without a metric,
   and every later xPand split then fails with InducedFrom::unknown (runs 20260916T113714Z / 114500Z). *)


(* ---------------- D3: SVT parametrization (8 S + 6 V + 2 T) ---------------- *)
Do[DefProjectedTensor[Symbol["TS" <> ToString[k]][], h], {k, 8}];
vecOK = RCTry[(Do[DefProjectedTensor[Symbol["TV" <> ToString[k]][-ia], h, TensorProperties -> {"Transverse"}], {k, 6}]; True)];
RCSayPlain["D3_TRANSVERSE_VECTORS_DECLARED", vecOK === True];
If[vecOK =!= True, Do[DefProjectedTensor[Symbol["TV" <> ToString[k]][-ia], h], {k, 6}]];
Do[DefProjectedTensor[Symbol["TT" <> ToString[k]][-ia, -ib], h, TensorProperties -> {"SymmetricTensor", "Traceless", "Transverse"}], {k, 2}];
RCSay["D3_TT1_TRANSVERSE", RCCanon[cd[-ia]@TT1[LI[1], ia, -ib]]];
RCSay["D3_TT1_TRACELESS", RCCanon[TT1[LI[1], ia, -ia]]];
RCSay["D3_TV1_TRANSVERSE", RCCanon[cd[-ia]@TV1[LI[1], ia]]];
(* general spatial rank-2 pieces W (for T^i_{0j}) and U (for the dual of T^i_{jk}) *)
Wp[o_, a_, c_, d_, e_] := h[a, c] TS3[LI[o]] + cd[a]@cd[c]@TS4[LI[o]] + (cd[a]@TV3[LI[o], c] + cd[c]@TV3[LI[o], a])/2 + TT1[LI[o], a, c] + epsilonh[a, c, -d] (cd[d]@TS5[LI[o]] + TV4[LI[o], d]);
Up[o_, a_, c_, d_, e_] := h[a, c] TS6[LI[o]] + cd[a]@cd[c]@TS7[LI[o]] + (cd[a]@TV5[LI[o], c] + cd[c]@TV5[LI[o], a])/2 + TT2[LI[o], a, c] + epsilonh[a, c, -e] (cd[e]@TS8[LI[o]] + TV6[LI[o], e]);
torsionRule = dTor[LI[o_], ia_, ib_, ic_] :> With[{id = DummyIn[Tangent[M4]], ie = DummyIn[Tangent[M4]], ig = DummyIn[Tangent[M4]]},
    n[ia] (n[ib] (cd[ic]@TS1[LI[o]] + TV1[LI[o], ic]) - n[ic] (cd[ib]@TS1[LI[o]] + TV1[LI[o], ib]))
  + n[ia] epsilonh[ib, ic, -id] (cd[id]@TS2[LI[o]] + TV2[LI[o], id])
  + (n[ib] Wp[o, ia, ic, ie, ig] - n[ic] Wp[o, ia, ib, ie, ig])
  + epsilonh[ib, ic, -id] Up[o, ia, id, ie, ig]];
bgTorsionZero = {Tor[ia_, ib_, ic_] -> 0};
rhs = dTor[LI[1], ia, -ib, -ic] /. torsionRule;
RCSay["D3_RHS_SYMBOLS", Union@Cases[rhs, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}]];
RCSayPlain["D3_SYMBOL_COUNT_16", Length@Union@Cases[rhs, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}] === 16];
antis = RCCanon[(dTor[LI[1], ia, -ib, -ic] + dTor[LI[1], ia, -ic, -ib]) /. torsionRule];
RCSayPlain["D3_ANTISYMMETRY_OK", antis === 0];
RCSay["D3_ANTISYMMETRY_RESIDUAL", antis];
reached["D3"];

(* ---------------- D4: projection of the rank-3 field ---------------- *)
projs = Tuples[{"Time", "Space"}, 3];
Do[With[{p = projs[[k]]},
   pc = RCTry[RCCanon[ExtractComponents[rhs, h, p]]];
   RCSay["D4_PROJ_" <> StringJoin[StringTake[#, 1] & /@ p], If[pc === $Failed, "$Failed", pc]]], {k, Length[projs]}];
reached["D4"];

(* ---------------- D6: the action route for the torsion sector ---------------- *)
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
rulesMetric = SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV;
rulesT = Join[rulesMetric, {torsionRule}, bgTorsionZero];
LTalg = Sqrt[-Detg[]] (c1 Tor[ia, -ib, -ic] Tor[-ia, ib, ic] + c2 Tor[ia, -ib, -ic] Tor[ib, -ia, ic] + c3 Tor[ia, -ib, -ia] Tor[ic, ib, -ic]);
LTkin = Sqrt[-Detg[]] c4 CD[-id]@Tor[ia, -ib, -ic] CD[id]@Tor[-ia, ib, ic];
{t6a, L2a} = AbsoluteTiming[RCTry[fixVol@ToxPandFromRules[LTalg, rulesT, h, 2]]];
RCSayPlain["D6A_ALGEBRAIC_SPLIT_WALL_S", Round[t6a, 0.01]];
RCSayPlain["D6A_ALGEBRAIC_SPLIT_RAN", L2a =!= $Failed];
{t6b, L2b} = AbsoluteTiming[RCTry[fixVol@ToxPandFromRules[LTkin, rulesT, h, 2]]];
RCSayPlain["D6B_KINETIC_SPLIT_WALL_S", Round[t6b, 0.01]];
RCSayPlain["D6B_KINETIC_SPLIT_RAN", L2b =!= $Failed];
L2T = Plus @@ Select[{L2a, L2b}, # =!= $Failed &];
If[L2T =!= 0,
  L1T = RCCanon[ExtractOrder[L2T, 1]]; L2Tq = RCCanon[ExtractOrder[L2T, 2]];
  RCSay["D6_TADPOLE_T0", L1T];
  RCSayPlain["D6_TADPOLE_VANISHES_FOR_TBAR_ZERO", L1T === 0];
  RCSayPlain["D6_L2_HAS_RESIDUAL_CD", !FreeQ[L2Tq, CD]];
  RCSayPlain["D6_L2_HAS_RESIDUAL_DTOR", !FreeQ[L2Tq, dTor]];
  RCSayPlain["D6_L2_HAS_METRIC_TT", !FreeQ[L2Tq, Eth]];
  RCSayPlain["D6_L2_TERM_COUNT", Length[List @@ Expand[L2Tq]]];
  zeroSVtor = Join[Table[Symbol["TS" <> ToString[k]][__] -> 0, {k, 8}], Table[Symbol["TV" <> ToString[k]][__] -> 0, {k, 6}], {Eth[__] -> 0}];
  L2TT = RCCanon[L2Tq /. zeroSVtor];
  RCSay["D6_L2_TT_SECTOR", L2TT];
  RCSayPlain["D6_TT_SECTOR_HAS_TT1_TT2_CROSS_TERM", !FreeQ[Expand[L2TT], TT1[__] TT2[__]] || !FreeQ[Expand[L2TT], TT2[__] TT1[__]]];
  (* Euler-Lagrange for the TT pieces (the same eta-space operator as A1) *)
  ClearAll[dEta];
  dEta[x_Plus] := dEta /@ x;
  dEta[x_Times] := Plus @@ MapIndexed[(ReplacePart[x, First[#2] -> dEta[#1]]) &, List @@ x];
  dEta[x_^m_] := m x^(m - 1) dEta[x];
  dEta[c_?NumericQ] := 0;
  dEta[c_Symbol] /; ConstantSymbolQ[c] := 0;
  dEta[cd[i_][x_]] := cd[i][dEta[x]];
  dEta[(h | delta | n | epsilonh)[__]] := 0;
  dEta[Scalar[x_]] := dEta[x];
  dEta[T_[LI[p_], LI[q_], inds___]] := T[LI[p], LI[q + 1], inds];
  dEta[0] = 0;
  Do[With[{Y = Symbol["Y" <> ToString[k]]},
    DefTensor[Y[-ia, -ib], M4, Symmetric[{-ia, -ib}], OrthogonalTo -> {n[ia], n[ib]}, ProjectedWith -> {h[ia, -ic], h[ib, -ic]}];
    Y /: OrthogonalToVectorQ[n][Y] = True], {k, 0, 2}];
  RCEL[L_, field_, ord_, {A0_, A1_, A2_}, indsPattern_] := Module[{Lx, P0, P1, P2, back},
    Lx = NoScalar[L] /. {field[LI[ord], LI[0], ii___] :> A0[ii], field[LI[ord], LI[1], ii___] :> A1[ii], field[LI[ord], LI[2], ii___] :> A2[ii]};
    P0 = VarD[A0 @@ indsPattern, cd][Lx]; P1 = VarD[A1 @@ indsPattern, cd][Lx]; P2 = VarD[A2 @@ indsPattern, cd][Lx];
    back = {A0[ii___] :> field[LI[ord], LI[0], ii], A1[ii___] :> field[LI[ord], LI[1], ii], A2[ii___] :> field[LI[ord], LI[2], ii]};
    RCSort[(P0 - dEta[RCCanon[P1 /. back]] + dEta[dEta[RCCanon[P2 /. back]]]) /. back]];
  elTT1 = RCTry[RCEL[L2TT, TT1, 1, {Y0, Y1, Y2}, {-ia, -ib}]];
  RCSay["D6_EL_TT1", elTT1];
  elTT2 = RCTry[RCEL[L2TT, TT2, 1, {Y0, Y1, Y2}, {-ia, -ib}]];
  RCSay["D6_EL_TT2", elTT2];
  RCSayPlain["D6_EL_TT1_FREE_OF_TT2", FreeQ[elTT1, TT2]];
  RCSayPlain["D6_EL_TT2_FREE_OF_TT1", FreeQ[elTT2, TT1]];
  RCSayPlain["D6_EL_TT1_HAS_SECOND_TIME_DERIVATIVE", !FreeQ[elTT1, TT1[LI[1], LI[2], __]]];
  RCSayPlain["D6_EL_TT1_HAS_LAPLACIAN", !FreeQ[elTT1, cd[_][cd[_][TT1[__]]]]];
  RCSayPlain["D6_EL_TT1_HAS_FRICTION", !FreeQ[elTT1, Hh[__] TT1[LI[1], LI[1], __]]];
  RCSayPlain["D6_EL_TT2_HAS_SECOND_TIME_DERIVATIVE", !FreeQ[elTT2, TT2[LI[1], LI[2], __]]];
  (* the scalar sector of the same action: does the split carry the eight scalars? *)
  zeroVTtor = Join[Table[Symbol["TV" <> ToString[k]][__] -> 0, {k, 6}], {TT1[__] -> 0, TT2[__] -> 0, Eth[__] -> 0}];
  L2S = RCCanon[L2Tq /. zeroVTtor];
  RCSay["D6_L2_SCALAR_SECTOR_SYMBOLS", Union@Cases[L2S, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}]];
];
reached["D6"];

(* ---------------- D7: the derivative split ---------------- *)
dv = RCTry[RCSort[ExtractOrder[ToxPandFromRules[CD[-id]@Tor[ia, -ib, -ic], rulesT, h, 1], 1]]];
RCSayPlain["D7_SPLIT_RAN", dv =!= $Failed];
If[dv =!= $Failed,
  RCSay["D7_ORDER1_SPLIT", dv];
  RCSay["D7_HEADS", Union@Cases[dv, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]];
  RCSayPlain["D7_NO_RESIDUAL_CD", FreeQ[dv, CD]];
  RCSayPlain["D7_NO_RESIDUAL_DTOR", FreeQ[dv, dTor]];
  RCSayPlain["D7_NO_RESIDUAL_LIED", FreeQ[dv, LieD]];
  RCSayPlain["D7_HAS_PRIMES", !FreeQ[dv, _[LI[1], LI[1], ___]]];
  RCSay["D7_MESSAGES", $MessageList];
  pieceX = If[dv === Null, "split returned Null (see probe_d_limits.wls for the rank map)", RCTry[RCCanon[ExtractComponents[dv, h, {"Time", "Time", "Time", "Space"}]]]];
  RCSay["D7_PIECE_TTTS_EXTRACTCOMPONENTS", pieceX];
  (* by hand, unambiguous slot order: n^d n_a n^b h^c_e  nabla_d T^a_bc *)
  piece = If[dv === Null, Null, RCTry[RCSort[n[id] n[-ia] n[ib] h[ic, -ie] dv]]];
  RCSay["D7_PIECE_TTTS_BY_HAND", piece];
  If[piece =!= $Failed && piece =!= Null && piece =!= 0,
    vpart = cd[-ie]@TS1[LI[1], LI[0]] + TV1[LI[1], LI[0], -ie];
    vpartP = cd[-ie]@TS1[LI[1], LI[1]] + TV1[LI[1], LI[1], -ie];
    Do[RCSayPlain["D7_PIECE_VS_HAND_FRICTION_" <> ToString[k], RCTry[Last@RCNormalize[piece, vpartP + k HH vpart, {1, -1, 2, -2, 1/2, -1/2}]]], {k, {-2, -1, 0, 1, 2}}];
  ];
];
reached["D7"];

(* ---------------- D8: the torsionful connection ---------------- *)
d8 = RCTry[(DefCovD[CDT[-ia], {"#", "DT"}, FromMetric -> g, Torsion -> True]; True)];
RCSayPlain["D8_DEFCOVD_TORSION_OK", d8 === True];
If[d8 === True,
  r1 = RCTry[ToxPandFromRules[RicciScalarCDT[], rulesT, h, 1]];
  RCSayPlain["D8_TOXPAND_ON_RICCISCALARCDT_RAN", r1 =!= $Failed];
  RCSay["D8_TOXPAND_ON_RICCISCALARCDT_RESULT_HEADS", If[r1 === $Failed, "$Failed", Union@Cases[r1, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]]];
  RCSayPlain["D8_RESULT_STILL_HAS_CDT_OBJECTS", !FreeQ[r1, RicciScalarCDT | RicciCDT | RiemannCDT | ChristoffelCDCDT | ChristoffelCDTCD | TorsionCDT | CDT]];
  (* legacy's route: ChangeCurvature to the Levi-Civita connection, contortion -> torsion (Shapiro 2002 eq. 2.9; _derive.py:2164-2189) *)
  lc = RCTry[ChangeCurvature[RicciScalarCDT[], CDT, CD]];
  RCSayPlain["D8_CHANGECURVATURE_RAN", lc =!= $Failed];
  If[lc =!= $Failed,
    RCSay["D8_LC_HEADS_BEFORE_CONTORTION_RULE", Union@Cases[lc, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]];
    lc = lc /. ChristoffelCDCDT[a1_, a2_, a3_] :> With[{d1 = DummyIn[Tangent[M4]], d2 = DummyIn[Tangent[M4]]},
          1/2 (TorsionCDT[a1, a2, a3] + g[a2, d1] g[a1, d2] TorsionCDT[-d1, -d2, a3] - g[a3, d1] g[a1, d2] TorsionCDT[a2, -d1, -d2])];
    lc = ContractMetric[Expand[lc]];
    RCSayPlain["D8_LC_REWRITE_FREE_OF_CHRISTOFFELCDCDT", FreeQ[lc, ChristoffelCDCDT]];
    RCSay["D8_LC_REWRITE_HEADS", Union@Cases[lc, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]];
    RCSay["D8_LC_REWRITE", RCCanon[lc]];
    lcT = lc /. TorsionCDT[a1_, a2_, a3_] :> Tor[a1, a2, a3];
    r2 = RCTry[RCSort[ExtractOrder[ToxPandFromRules[lcT, rulesT, h, 1], 1]]];
    RCSayPlain["D8_SPLIT_OF_REWRITTEN_RICCI_RAN", r2 =!= $Failed && r2 =!= Null];
    RCSay["D8_SPLIT_OF_REWRITTEN_RICCI", r2];
    If[r2 =!= $Failed,
      rR = RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
      RCSayPlain["D8_METRIC_ONLY_SPLIT_STILL_WORKS", rR =!= $Failed && rR =!= Null];
      diffR = If[rR === $Failed || rR === Null || r2 === Null, $Failed, RCTry[RCSort[r2 - rR]]];
      RCSay["D8_DELTA_RTILDE_MINUS_DELTA_R", diffR];
      RCSayPlain["D8_DIFFERENCE_IS_LINEAR_IN_TORSION_ONLY", diffR =!= $Failed && FreeQ[diffR, Eth] && !FreeQ[diffR, TS1 | TS3 | TV1]];
      RCSayPlain["D8_DIFFERENCE_SYMBOLS", ToString[Union@Cases[diffR, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}]]];
    ];
  ];
];
reached["D8"];

(* The optional D5 step (key d5=1), an FRW-symmetric background torsion mode, was removed on
   2026-09-18: it was never run, it is out of scope (the programme expands about zero
   background torsion, exactly), and its background rule has the shape xPand cannot prepare
   (a sum of terms along the normal for a rank-3 field; GH #591). *)
(* ---------------- D2: 24-component count (last, see the note above) ---------------- *)
Needs["xAct`xCoba`"];   (* its own statement: DefChart must be parsed after xCoba is on the context path *)
d2 = RCTry[(DefChart[Bc, M4, {0, 1, 2, 3}, {tt[], xx[], yy[], zz[]}];
   comps = Flatten@Table[ToCanonical@Tor[{i, Bc}, {-j, -Bc}, {-k, -Bc}], {i, 0, 3}, {j, 0, 3}, {k, 0, 3}];
   Length@DeleteDuplicates[DeleteCases[comps, 0] /. Times[-1, x_] :> x])];
RCSayPlain["D2_INDEPENDENT_COMPONENTS_XCOBA", d2];
RCSayPlain["D2_COUNT_IS_24", d2 === 24];
RCSayPlain["D2_INDEPENDENT_COMPONENTS_COMBINATORIAL", 4 Binomial[4, 2]];
reached["D2"];
RCDone["d"];
