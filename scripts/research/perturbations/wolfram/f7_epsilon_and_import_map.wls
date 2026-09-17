#!/usr/bin/env wolframscript
(* cspell:words xPand ToxPandFromRules SplitMetric PSALTer DefGeometry ToEps FromEps epsilonh epsilong SignDetOfMetric epsilonOrientation ConformalWeight DefProjectedTensor TensorProperties DefTensorPerturbation Aoki Bahamonde Gigante Valcarcel Gorji Nikiforova Damour contortion normu EPSH EPSG dTor inds Detg NORMU PSALTER expandGdelta MasterOf FirstMetricOfVBundle ExtractComponents ExtractOrder SIGNDET EPSILONH EPSILONG projs FROMEPS HTRACE Lalg Lkin Lodd TOEPS UPVALUE UPVALUES VBUNDLE XPANDS abce dabc indep *)
(* f7_epsilon_and_import_map.wls -- R-C closing pass, step f7.

   PART A -- THE SLICE EPSILON. xTensor auto-creates epsilon[h] when xPand declares the
   induced metric (xTensor.m:8115-8127: dim-1 slots, OrthogonalTo and ProjectedWith on the
   FIRST slot only) and fixes its square through an UpValue on epsilonh itself
   (xTensor.m:8134-8137, which rewrites the product to SignDetOfMetric[h] expandGdelta[h][..],
   a determinant of h's -- so the raw product is NOT the number 6; only its canonical form is).
   What no package fixes is the ORIENTATION relative to the four-index epsilon and the normal,
   so a formula written in 4D cannot be compared with a split expression until the link is
   stated. The link is not ours to invent: PSALTer, the programme's owner of the epsilon
   convention, already declares its own slice epsilon and both directions of the dictionary --
     PSALTer/Sources/ReloadPackage/DefGeometry.m:41-46  V[-a] with V[-a] V[a] == 1, and
       Eps[-a,-b,-c] Antisymmetric, OrthogonalTo -> {V[a],V[b],V[c]}
     :60  ToEps   = MakeRule[{epsilonG[-a,-b,-c,-d] V[d], Eps[-a,-b,-c]}, ...]
     :62  FromEps = MakeRule[{Eps[-a,-b,-c], epsilonG[-a,-b,-c,-d] V[d]}, ...]
     :71-72  components -$epsilonSign all-up and +$epsilonSign all-down
   i.e. eps3_abc = eps4_abcd n^d, coefficient +1, contracted index in the LAST slot.
   A therefore TRANSCRIBES and tests consistency; it does not derive the sign. Squaring is
   even in the epsilon, and every abstract-index identity xTensor offers is of that kind
   (epsilonOrientation is read only by xCoba component conversion, xTensor.m:8173), so
   F7A_SIGN_FIXABLE_BY_SQUARE is printed as the constant False so no later reader mistakes
   A8 for a derivation.
   Two traps are pinned because each flips every parity-odd term silently: the slot order
   (eps4_abcd n^d = -eps4_dabc n^d, so the natural first-slot writing is the opposite sign)
   and, for the literature, Nikiforova-Damour 1804.09215:245 stating eps^0123 = +1, which is
   eps_0123 = -1 in any Lorentzian 4D signature and hence opposite to Aoki et al.
   2310.16007:201, from whom Probe D3's parity-odd rules were transcribed.
   A8 also measures a defect of xPand's own: xPand.m:1754 declares the induced metric with
   DefMetric[1, h[..], cd, .., InducedFrom -> {g,u}], hard-coding the slice determinant sign
   to +1 whatever normu is. In mostly-minus the slice metric is negative definite and the
   true sign is -1, so the run under normu=+1 is PREDICTED to fail the square test, and that
   failure is xPand's, not a clash with PSALTer -- PSALTer's Eps is a plain DefTensor with no
   determinant UpValue and is self-consistent at Eps.Eps = -6 through FromEps.
   A10 has its own sentinel because xPand.m:2528-2529 carries a live author comment that the
   conformal rule for the epsilon of a frozen metric "is not working well ... really works
   only when metric1 is the ambient metric".

   PART B -- THE FAMILY-A IMPORT MAP, RUN. The memo asserts that mapping an imported family-A
   torsion onto ours (T -> -T) flips the sign of Probe D3's 24 potentials and leaves the
   counting and the action route untouched. Here D3's rule set is rebuilt, tied back to
   probe_d's recorded projections, the 3+3+9+9 = 24 accounting is established block by block,
   the parity-odd half is identified, and both order-2 actions are compared -- with a
   Lagrangian ODD in the torsion as the control that must flip, so the invariance verdict is
   not vacuous. B8 then names the two cases where the sign IS observable.

   Usage: run_lane.sh f7 [normu=1|-1] [kin=1] . normu=+1 is the programme's mostly-minus
   convention and PSALTer's (V.V = +1); normu=-1 is the mostly-plus setting Probe D3 and f2
   were run in. No DefChart anywhere in this file: a chart poisons every later xPand split
   in the same kernel (xPand.m:2566, f2 case (i)). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
F7Normu = ToExpression[Lookup[RCArgs, "normu", "-1"]];
RCSayPlain["F7_NORMU", F7Normu];
RCGeometry[F7Normu, "FLFlat"];
RCSigns["f7_start"];
DefConstantSymbol[c1]; DefConstantSymbol[c2]; DefConstantSymbol[c3]; DefConstantSymbol[c4]; DefConstantSymbol[c5];
AA = ah[LI[0], LI[0]]; HH = Hh[LI[0], LI[0]];
heads[x_] := Union@Cases[x, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}];
syms[x_] := Union@Cases[x, (s_Symbol)[LI[1], ___] :> s, {0, Infinity}];
aCands = {1, -1, AA, 1/AA, AA^2, 1/AA^2, AA^3, 1/AA^3, -AA^3, -1/AA^3};
fixVol[x_] := NoScalar[(x /. Detg[] -> -1) /. Sqrt[y_] :> PowerExpand[Sqrt[y]]];

(* ==================== A1: what the package fixes by itself ==================== *)
RCSayPlain["F7A_EPSILON_SPELLING_EQUIV", epsilon[h] === epsilonh && epsilon[g] === epsilong];
RCSay["F7A_EPSH_SLOTS", RCTry[SlotsOfTensor[epsilonh]]];
RCSay["F7A_EPSG_SLOTS", RCTry[SlotsOfTensor[epsilong]]];
RCSay["F7A_EPSH_SYMMETRY", RCTry[SymmetryGroupOfTensor[epsilonh[-ia, -ib, -ic]]]];
RCSay["F7A_EPSH_MASTER", RCTry[MasterOf[epsilonh]]];
RCSay["F7A_SIGNDET_G_H", RCTry[{SignDetOfMetric[g], SignDetOfMetric[h]}]];
RCSay["F7A_EPSILON_ORIENTATION_AND_SIGN",
  RCTry[{epsilonOrientation[g, AIndex], epsilonOrientation[h, AIndex], xAct`xTensor`$epsilonSign}]];
RCSayPlain["F7A_ORIENTATION_READ_BY_XCOBA_ONLY",
  "epsilonOrientation is consumed by component conversion (xTensor.m:8173), never by abstract-index work"];
RCSay["F7A_NN_AND_HTRACE", RCTry[{RCCanon[n[ia] n[-ia]], RCCanon[h[ia, -ia]]}]];
RCSay["F7A_EPSH_UPVALUE_COUNT", Length[UpValues[epsilonh]]];
RCSay["F7A_EPSH_CD_DERIVATIVE", RCTry[RCCanon[cd[-id]@epsilonh[-ia, -ib, -ic]]]];
RCSay["F7A_EPSH_LIE_ALONG_N", RCTry[RCCanon[LieD[n[id]][epsilonh[-ia, -ib, -ic]]]]];
RCSay["F7A_FIRST_METRIC_OF_VBUNDLE", RCTry[FirstMetricOfVBundle[Tangent[M4], True]]];
RCSay["F7A_METRICS_ORDER", $Metrics];

(* ==================== A2: normalization, raw versus canonical ==================== *)
sqRaw = epsilonh[-ia, -ib, -ic] epsilonh[ia, ib, ic];
RCSayPlain["F7A_EPSH_SQUARE_RAW_IS_SIX", sqRaw === 6];
RCSayPlain["F7A_EPSH_SQUARE_RAW_FREE_OF_EPSH", FreeQ[sqRaw, epsilonh]];
RCSay["F7A_EPSH_SQUARE_RAW_TERM_COUNT", RCTry[Length[List @@ Expand[sqRaw]]]];
sqH = RCTry[RCCanon[sqRaw]];
RCSay["F7A_EPSH_SQUARE_CANON", sqH];
RCSayPlain["F7A_EPSH_SQUARE_CANON_IS_SIX", sqH === 6];
sqG = RCTry[RCCanon[epsilong[-ia, -ib, -ic, -id] epsilong[ia, ib, ic, id]]];
RCSay["F7A_EPSG_SQUARE_CANON", sqG];
RCSayPlain["F7A_EPSG_SQUARE_IS_MINUS_24", sqG === -24];

(* ==================== A3: orthogonality to n, slot by slot ==================== *)
o1 = n[ia] epsilonh[-ia, -ib, -ic]; o2 = n[ib] epsilonh[-ia, -ib, -ic]; o3 = n[ic] epsilonh[-ia, -ib, -ic];
Do[RCSayPlain["F7A_ORTHO_SLOT" <> ToString[First@p] <> "_RAW_ZERO", Last[p] === 0],
   {p, {{1, o1}, {2, o2}, {3, o3}}}];
Do[RCSayPlain["F7A_ORTHO_SLOT" <> ToString[First@p] <> "_CANON_ZERO", RCTry[RCCanon[Last@p]] === 0],
   {p, {{1, o1}, {2, o2}, {3, o3}}}];
RCSayPlain["F7A_HTRACE_ON_EPSH_ZERO", RCTry[RCCanon[h[ia, ib] epsilonh[-ia, -ib, -ic]]] === 0];

(* ==================== A4: total antisymmetry ==================== *)
RCSayPlain["F7A_ANTISYM_SWAP12", RCTry[RCCanon[epsilonh[-ia, -ib, -ic] + epsilonh[-ib, -ia, -ic]]] === 0];
RCSayPlain["F7A_ANTISYM_SWAP23", RCTry[RCCanon[epsilonh[-ia, -ib, -ic] + epsilonh[-ia, -ic, -ib]]] === 0];
RCSayPlain["F7A_CYCLIC_INVARIANT", RCTry[RCCanon[epsilonh[-ia, -ib, -ic] - epsilonh[-ib, -ic, -ia]]] === 0];

(* ==================== A5: the absence of the relation, asserted positively ==================== *)
mixed = epsilonh[-ia, -ib, -ic] epsilong[ia, ib, ic, id] n[-id];
RCSayPlain["F7A_MIXED_PRODUCT_INERT", !FreeQ[mixed, epsilonh] && !FreeQ[mixed, epsilong]];
RCSay["F7A_MIXED_PRODUCT_CANON_HEADS", RCTry[heads[RCCanon[mixed]]]];
ruleStore = Join[UpValues[epsilonh], DownValues[epsilonh], UpValues[epsilong], DownValues[epsilong],
                 If[Names["xAct`xTensor`$Rules"] === {}, {}, xAct`xTensor`$Rules]];
RCSay["F7A_RULE_STORE_SIZE", Length[ruleStore]];
RCSayPlain["F7A_NO_RULE_MENTIONS_BOTH_EPSILONS",
  Count[ruleStore, r_ /; (!FreeQ[r, epsilonh] && !FreeQ[r, epsilong])] === 0];
RCSayPlain["F7A_XPAND_EPSG_U_H_RULE_PRESENT",
  Count[ruleStore, r_ /; (!FreeQ[r, epsilong] && !FreeQ[r, h])] >= 1];
(* the two counts above came out unexpectedly on the first run, so print the rules themselves:
   which entries mention both epsilons, and which mention the four-index one at all *)
RCSay["F7A_RULES_MENTIONING_BOTH_EPSILONS",
  StringTake[ToString[#, InputForm], UpTo[220]] & /@
    Select[ruleStore, (!FreeQ[#, epsilonh] && !FreeQ[#, epsilong]) &]];
RCSay["F7A_RULES_MENTIONING_EPSG_COUNT", Count[ruleStore, r_ /; !FreeQ[r, epsilong]]];
RCSay["F7A_UPVALUES_OF_EPSG_LHS",
  StringTake[ToString[First[#], InputForm], UpTo[110]] & /@ UpValues[epsilong]];
(* the projection rule xPand installs at xPand.m:1809 is
   AutomaticRules[epsilon[g], {epsilong[i1,i2,i3,i4] u[-i1] h[-i2,e] -> epsilong[i1,e,i3,i4] u[-i1]}]
   where u is the vector passed to SetSlicing, i.e. our n. Test its exact shape, not a
   paraphrase, and print whether any stored rule mentions h at all. *)
RCSayPlain["F7A_XPAND_1809_RULE_FIRES",
  RCTry[RCCanon[epsilong[ia, ib, ic, id] n[-ia] h[-ib, ie] - epsilong[ia, ie, ic, id] n[-ia]]] === 0];
RCSayPlain["F7A_ANY_EPSG_RULE_MENTIONS_H",
  Count[ruleStore, r_ /; (!FreeQ[r, epsilong] && !FreeQ[r, h])]];
RCSay["F7A_RULE_STORE_PARTS",
  {Length[UpValues[epsilonh]], Length[DownValues[epsilonh]],
   Length[UpValues[epsilong]], Length[DownValues[epsilong]],
   If[Names["xAct`xTensor`$Rules"] === {}, "no-$Rules", Length[xAct`xTensor`$Rules]]}];
RCSay["F7A_EPSG_U_H_IDEMPOTENCE",
  RCTry[RCCanon[epsilong[ia, ib, ic, id] n[-ia] h[-ib, ie] - epsilong[ia, ie, ic, id] n[-ia]]]];
RCSayPlain["F7A_EPSH_VS_PLUS_EPSG_N_VERDICT",
  RCTry[RCVerdict[epsilonh[-ia, -ib, -ic], n[id] epsilong[-id, -ia, -ib, -ic]]]];
RCSayPlain["F7A_EPSH_VS_MINUS_EPSG_N_VERDICT",
  RCTry[RCVerdict[epsilonh[-ia, -ib, -ic], -n[id] epsilong[-id, -ia, -ib, -ic]]]];

(* ==================== A6: PSALTer's dictionary, and the slot-order control ==================== *)
rhsLast = epsilong[-ia, -ib, -ic, -id] n[id];
rhsFirst = n[id] epsilong[-id, -ia, -ib, -ic];
RCSayPlain["F7A_ADOPTED_RULE",
  "epsilonh[-a,-b,-c] == epsilong[-a,-b,-c,-d] n[d]   (PSALTer DefGeometry.m:60-63, last slot, +1)"];
RCSayPlain["F7A_SLOT_ORDER_SIGN", RCTry[Last@RCNormalize[rhsLast, rhsFirst, {1, -1}]]];
RCSayPlain["F7A_SLOT_ORDER_SUM_IS_ZERO", RCTry[RCCanon[rhsLast + rhsFirst]] === 0];
RCSay["F7A_PSALTER_V_NORM_VS_XPAND_NORMU", {1, F7Normu}];
RCSayPlain["F7A_NORMALIZATIONS_AGREE", F7Normu === 1];

(* ==================== A7: the LHS's properties, checked on the RHS ==================== *)
RCSayPlain["F7A_RHS_ORTHO_SLOT1_ZERO", RCTry[RCCanon[n[ia] rhsLast]] === 0];
RCSayPlain["F7A_RHS_ORTHO_SLOT2_ZERO", RCTry[RCCanon[n[ib] rhsLast]] === 0];
RCSayPlain["F7A_RHS_ORTHO_SLOT3_ZERO", RCTry[RCCanon[n[ic] rhsLast]] === 0];
RCSayPlain["F7A_RHS_ANTISYM_SWAP12",
  RCTry[RCCanon[rhsLast + epsilong[-ib, -ia, -ic, -id] n[id]]] === 0];
RCSayPlain["F7A_RHS_HTRACE_ZERO", RCTry[RCCanon[h[ia, ib] rhsLast]] === 0];

(* ==================== A8: the one quantitative test, and its explicit limit ==================== *)
sqRhs = RCTry[RCCanon[rhsLast (epsilong[ia, ib, ic, ie] n[-ie])]];
RCSay["F7A_RHS_SQUARE", sqRhs];
RCSayPlain["F7A_SQUARE_CONSISTENT", sqRhs === sqH];
RCSayPlain["F7A_S_SQUARED_NORMALIZATION", RCTry[Last@RCNormalize[sqRhs, sqH, {1, -1}]]];
RCSayPlain["F7A_SQUARE_IS_SIGN_BLIND",
  RCTry[RCCanon[(-epsilonh[-ia, -ib, -ic]) (-epsilonh[ia, ib, ic])]] === sqH];
RCSayPlain["F7A_SIGN_FIXABLE_BY_SQUARE", False];
RCSay["F7A_SIGNDET_H_DECLARED_VS_TRUE", RCTry[{SignDetOfMetric[h], If[F7Normu === -1, 1, -1]}]];
RCSayPlain["F7A_XPAND_1754_DETERMINANT_HARDCODE_BITES",
  RCTry[SignDetOfMetric[h]] =!= If[F7Normu === -1, 1, -1]];
RCSayPlain["F7A_MOSTLY_MINUS_FAILURE_IS_XPANDS_NOT_A_CONVENTION_CLASH",
  "xPand.m:1754 declares DefMetric[1,h,..] whatever normu is; PSALTer's Eps is a plain DefTensor and is self-consistent at -6"];

(* ==================== A9: the dictionary as a usable rule, and the round trip ==================== *)
toEps = RCTry[MakeRule[{epsilong[-ia, -ib, -ic, -id] n[id], epsilonh[-ia, -ib, -ic]},
                       MetricOn -> All, ContractMetrics -> True]];
fromEps = RCTry[MakeRule[{epsilonh[-ia, -ib, -ic], epsilong[-ia, -ib, -ic, -id] n[id]},
                         MetricOn -> All, ContractMetrics -> True]];
RCSayPlain["F7A_MAKERULE_TOEPS_RAN", toEps =!= $Failed];
RCSayPlain["F7A_MAKERULE_FROMEPS_RAN", fromEps =!= $Failed];
RCSay["F7A_MAKERULE_COUNTS", {Length[Flatten@{toEps}], Length[Flatten@{fromEps}]}];
If[toEps =!= $Failed && fromEps =!= $Failed,
  RCSayPlain["F7A_ROUNDTRIP_VERDICT",
    RCTry[RCVerdict[RCCanon[(epsilonh[-ia, -ib, -ic] /. fromEps) /. toEps], epsilonh[-ia, -ib, -ic]]]];
  RCSay["F7A_CONVERTED_SQUARE", RCTry[RCCanon[sqRaw /. fromEps]]];
  RCSayPlain["F7A_CONVERTED_SQUARE_MATCHES", RCTry[RCCanon[sqRaw /. fromEps]] === sqH];
];

(* ==================== PART B: Probe D3's rule set, both contortion families ==================== *)
DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
DefTensorPerturbation[dTor[LI[o], ia, -ib, -ic], Tor[ia, -ib, -ic], M4];
Do[DefProjectedTensor[Symbol["TS" <> ToString[k]][], h], {k, 8}];
Do[DefProjectedTensor[Symbol["TV" <> ToString[k]][-ia], h, TensorProperties -> {"Transverse"}], {k, 6}];
Do[DefProjectedTensor[Symbol["TT" <> ToString[k]][-ia, -ib], h,
    TensorProperties -> {"SymmetricTensor", "Traceless", "Transverse"}], {k, 2}];
RCSayPlain["F7B_TT_DECLARATIONS_FIRE",
  RCTry[{RCCanon[cd[-ia]@TT1[LI[1], ia, -ib]], RCCanon[TT1[LI[1], ia, -ia]],
         RCCanon[cd[-ia]@TV1[LI[1], ia]]}] === {0, 0, 0}];
(* PARENTHESIZED ON PURPOSE. Without the outer bracket the parser ends the definition at the
   first newline where the expression is already complete, so the first run of this script
   defined Wp as its first two terms only and Part B counted 8 potentials instead of 16.
   That is the CLAUDE.md pitfall "parenthesize multiline Lagrangians"; probe_d_torsion.wls
   avoided it by keeping each of these on one line. *)
Wp[o_, aa_, cc_, dd_, ee_] := (
   h[aa, cc] TS3[LI[o]] + cd[aa]@cd[cc]@TS4[LI[o]]
   + (cd[aa]@TV3[LI[o], cc] + cd[cc]@TV3[LI[o], aa])/2 + TT1[LI[o], aa, cc]
   + epsilonh[aa, cc, -dd] (cd[dd]@TS5[LI[o]] + TV4[LI[o], dd]));
Up[o_, aa_, cc_, dd_, ee_] := (
   h[aa, cc] TS6[LI[o]] + cd[aa]@cd[cc]@TS7[LI[o]]
   + (cd[aa]@TV5[LI[o], cc] + cd[cc]@TV5[LI[o], aa])/2 + TT2[LI[o], aa, cc]
   + epsilonh[aa, cc, -ee] (cd[ee]@TS8[LI[o]] + TV6[LI[o], ee]));
(* mkRule[s]: s = +1 is ours, family B (xTensor's slot order); s = -1 is what an imported
   family-A parametrization becomes under the memo's import map T_A = -T_B. *)
mkRule[s_] := dTor[LI[o_], ia_, ib_, ic_] :> s With[
    {id = DummyIn[Tangent[M4]], ie = DummyIn[Tangent[M4]], ig = DummyIn[Tangent[M4]]},
    n[ia] (n[ib] (cd[ic]@TS1[LI[o]] + TV1[LI[o], ic]) - n[ic] (cd[ib]@TS1[LI[o]] + TV1[LI[o], ib]))
  + n[ia] epsilonh[ib, ic, -id] (cd[id]@TS2[LI[o]] + TV2[LI[o], id])
  + (n[ib] Wp[o, ia, ic, ie, ig] - n[ic] Wp[o, ia, ib, ie, ig])
  + epsilonh[ib, ic, -id] Up[o, ia, id, ie, ig]];
ruleB = mkRule[1]; ruleA = mkRule[-1];
bgTorsionZero = {Tor[ia_, ib_, ic_] -> 0};
rhs = dTor[LI[1], ia, -ib, -ic] /. ruleB;

(* ---- B1: provenance, tied back to probe_d's recorded run (d/20260916T113923Z) ---- *)
RCSayPlain["F7B_SYMBOL_COUNT_16", Length[syms[rhs]] === 16];
RCSay["F7B_SYMBOLS", syms[rhs]];
RCSayPlain["F7B_ANTISYMMETRY_OK",
  RCTry[RCCanon[(dTor[LI[1], ia, -ib, -ic] + dTor[LI[1], ia, -ic, -ib]) /. ruleB]] === 0];
projs = Tuples[{"Time", "Space"}, 3];
blockSyms = {};
Do[With[{p = projs[[k]], tag = StringJoin[StringTake[#, 1] & /@ projs[[k]]]},
   pc = RCTry[RCCanon[ExtractComponents[rhs, h, p]]];
   AppendTo[blockSyms, {tag, If[pc === $Failed, {}, syms[pc]]}];
   RCSay["F7B_PROJ_" <> tag <> "_SYMBOLS", If[pc === $Failed, "$Failed", syms[pc]]]], {k, Length[projs]}];
RCSay["F7B_BLOCK_SYMBOLS", blockSyms];
allBlock = Union @@ (Last /@ blockSyms);
RCSayPlain["F7B_BLOCK_UNION_IS_16", Length[allBlock] === 16];
RCSayPlain["F7B_NO_SYMBOL_CROSSES_BLOCKS", Total[Length /@ (Last /@ blockSyms)] === 16];
dofOf[ss_] := Total[If[StringMatchQ[ToString[#], "TS" ~~ __], 1,
                       If[StringMatchQ[ToString[#], "TV" ~~ __], 2, 2]] & /@ ss];
RCSay["F7B_BLOCK_DOF", {First[#], dofOf[Last[#]]} & /@ blockSyms];
(* The eight projections are not eight independent blocks: the torsion is antisymmetric in
   its last two slots, so TTS and TST are the same components up to sign, as are STS and SST.
   Canonicalizing a tag by sorting its last two characters pairs them, and the independent
   blocks then add up to the papers' 24. Summing all eight gives 36 and is double counting. *)
canonTag[t_String] := StringTake[t, 1] <> StringJoin[Sort[{StringTake[t, {2}], StringTake[t, {3}]}]];
indepBlocks = DeleteDuplicatesBy[blockSyms, canonTag[First[#]] &];
RCSay["F7B_INDEPENDENT_BLOCKS", {canonTag[First[#]], dofOf[Last[#]]} & /@ indepBlocks];
RCSayPlain["F7B_DOF_TOTAL_ALL_EIGHT_PROJECTIONS", Total[dofOf[Last[#]] & /@ blockSyms]];
RCSayPlain["F7B_DOF_TOTAL_INDEPENDENT", Total[dofOf[Last[#]] & /@ indepBlocks]];
RCSayPlain["F7B_DOF_TOTAL_IS_24", Total[dofOf[Last[#]] & /@ indepBlocks] === 24];

(* ---- B3: the parity split -- which half the Part A convention governs ---- *)
all16 = syms[rhs];
keepOnly[s_] := RCCanon[rhs /. Join[
   Cases[all16, x_ /; x =!= s :> (x[y___] -> 0)]]];
epsCount[x_] := Union@Cases[x, t_Times :> Count[t, epsilonh[__]], {0, Infinity}];
oddSet = Select[all16, (RCTry[!FreeQ[keepOnly[#], epsilonh]]) === True &];
RCSay["F7B_PARITY_ODD_SYMBOLS", oddSet];
RCSay["F7B_PARITY_EVEN_SYMBOLS", Complement[all16, oddSet]];
RCSayPlain["F7B_PARITY_ODD_DOF", dofOf[oddSet]];
RCSayPlain["F7B_PARITY_EVEN_DOF", dofOf[Complement[all16, oddSet]]];
RCSayPlain["F7B_TS8_TV6_EVEN_BY_EPS_EPS_COLLAPSE",
  RCTry[FreeQ[keepOnly[TS8], epsilonh] && FreeQ[keepOnly[TV6], epsilonh]]];

(* ---- B4: the two readings of the flip, and that they coincide ---- *)
flipAll = Join[Cases[all16, x_ :> (x[y___] :> -x[y])]];
RCSayPlain["F7B_MINUS_MAP_EQUALS_POTENTIAL_FLIP",
  RCTry[RCVerdict[RCCanon[dTor[LI[1], ia, -ib, -ic] /. ruleA], RCCanon[rhs /. flipAll]]]];
RCSayPlain["F7B_RHS_A_IS_MINUS_RHS_B", RCTry[RCVerdict[RCCanon[dTor[LI[1], ia, -ib, -ic] /. ruleA], -rhs]]];
RCSayPlain["F7B_SYMBOL_SETS_IDENTICAL", syms[dTor[LI[1], ia, -ib, -ic] /. ruleA] === all16];
RCSayPlain["F7B_ANTISYMMETRY_A_OK",
  RCTry[RCCanon[(dTor[LI[1], ia, -ib, -ic] + dTor[LI[1], ia, -ic, -ib]) /. ruleA]] === 0];

(* ---- B5: the action route in both families ---- *)
zeroSV = {\[Phi]h[__] -> 0, \[Psi]h[__] -> 0, Bsh[__] -> 0, Bvh[__] -> 0, Esh[__] -> 0, Evh[__] -> 0};
rulesMetric = SplitMetric[g, dg, h, "AnyGauge"] /. zeroSV;
rulesTB = Join[rulesMetric, {ruleB}, bgTorsionZero];
rulesTA = Join[rulesMetric, {ruleA}, bgTorsionZero];
LTalg = Sqrt[-Detg[]] (c1 Tor[ia, -ib, -ic] Tor[-ia, ib, ic] + c2 Tor[ia, -ib, -ic] Tor[ib, -ia, ic]
        + c3 Tor[ia, -ib, -ia] Tor[ic, ib, -ic]);
LTkin = Sqrt[-Detg[]] c4 CD[-id]@Tor[ia, -ib, -ic] CD[id]@Tor[-ia, ib, ic];
LTodd = Sqrt[-Detg[]] c5 n[ib] Tor[ia, -ib, -ia];
splitTo[L_, rules_, ord_] := RCTry[fixVol@ToxPandFromRules[L, rules, h, ord]];
{tAB, LalgB} = AbsoluteTiming[splitTo[LTalg, rulesTB, 2]];
{tAA, LalgA} = AbsoluteTiming[splitTo[LTalg, rulesTA, 2]];
RCSayPlain["F7B_ALG_SPLIT_WALL_S", {Round[tAB, 0.01], Round[tAA, 0.01]}];
RCSayPlain["F7B_ALG_SPLIT_RAN_BOTH", LalgB =!= $Failed && LalgA =!= $Failed];
If[LalgB =!= $Failed && LalgA =!= $Failed,
  L1B = RCCanon[ExtractOrder[LalgB, 1]]; L1A = RCCanon[ExtractOrder[LalgA, 1]];
  L2B = RCCanon[ExtractOrder[LalgB, 2]]; L2A = RCCanon[ExtractOrder[LalgA, 2]];
  RCSayPlain["F7B_TADPOLES_VANISH_BOTH", L1B === 0 && L1A === 0];
  RCSayPlain["F7B_ALG_ORDER2_TERM_COUNTS",
    {Length[List @@ Expand[NoScalar[L2B]]], Length[List @@ Expand[NoScalar[L2A]]]}];
  RCSayPlain["F7B_ALG_ORDER2_A_VS_B_VERDICT", RCTry[RCVerdict[L2A, L2B]]];
  RCSayPlain["F7B_ALG_ORDER2_NORMALIZATION", RCTry[Last@RCNormalize[L2A, L2B, {1, -1}]]];
  RCSayPlain["F7B_ALG_ORDER2_DIGESTS_EQUAL", RCDigest[L2A] === RCDigest[L2B]];
  zeroSVtor = Join[Table[Symbol["TS" <> ToString[k]][__] -> 0, {k, 8}],
                   Table[Symbol["TV" <> ToString[k]][__] -> 0, {k, 6}], {Eth[__] -> 0}];
  L2TT = RCCanon[L2B /. zeroSVtor]; L2TTm = RCCanon[L2A /. zeroSVtor];
  RCSayPlain["F7B_TT_SECTOR_A_VS_B_VERDICT", RCTry[RCVerdict[L2TTm, L2TT]]];
  RCSayPlain["F7B_TT_CROSS_TERM_PRESENT",
    RCTry[!FreeQ[Expand[NoScalar[L2TT]], TT2] && !FreeQ[Expand[NoScalar[L2TT]], TT1]]];
  RCSayPlain["F7B_PROBE_D_CROSS_TERM_SENTINEL_WAS_A_FALSE_NEGATIVE",
    RCTry[FreeQ[Expand[NoScalar[L2TT]], TT1[__] TT2[__]] && !FreeQ[Expand[NoScalar[L2TT]], TT1] && !FreeQ[Expand[NoScalar[L2TT]], TT2]]];
  (* A10 follow-through: the parity-odd sector's a-powers against the parity-even sector's *)
  ahPowers[x_] := Union@Cases[Expand[NoScalar[x]], AA^p_. :> p, {0, Infinity}];
  onlyTS[k_Integer] := RCCanon[L2B /. Cases[all16, x_ /; x =!= Symbol["TS" <> ToString[k]] :> (x[y___] -> 0)]];
  RCSay["F7A10_AH_POWERS_TS1_PARITY_EVEN", RCTry[ahPowers[onlyTS[1]]]];
  RCSay["F7A10_AH_POWERS_TS2_PARITY_ODD", RCTry[ahPowers[onlyTS[2]]]];
  RCSayPlain["F7A10_PARITY_ODD_AH_POWERS_MATCH_PARITY_EVEN",
    RCTry[ahPowers[onlyTS[1]] === ahPowers[onlyTS[2]]]];
  RCSayPlain["F7A10_EPSH_SURVIVES_THE_SPLIT", RCTry[!FreeQ[L2B, epsilonh]]];
];
If[Lookup[RCArgs, "kin", "0"] === "1",
  {tKB, LkinB} = AbsoluteTiming[splitTo[LTkin, rulesTB, 2]];
  {tKA, LkinA} = AbsoluteTiming[splitTo[LTkin, rulesTA, 2]];
  RCSayPlain["F7B_KIN_SPLIT_WALL_S", {Round[tKB, 0.01], Round[tKA, 0.01]}];
  RCSayPlain["F7B_KIN_SPLIT_RAN_BOTH", LkinB =!= $Failed && LkinA =!= $Failed];
  If[LkinB =!= $Failed && LkinA =!= $Failed,
    RCSayPlain["F7B_KIN_ORDER2_A_VS_B_VERDICT",
      RCTry[RCVerdict[RCCanon[ExtractOrder[LkinA, 2]], RCCanon[ExtractOrder[LkinB, 2]]]]]];
];

(* ---- B7: the control that makes the invariance verdict non-vacuous ---- *)
LoddB = splitTo[LTodd, rulesTB, 1]; LoddA = splitTo[LTodd, rulesTA, 1];
RCSayPlain["F7B_ODD_CONTROL_SPLIT_RAN_BOTH", LoddB =!= $Failed && LoddA =!= $Failed];
If[LoddB =!= $Failed && LoddA =!= $Failed,
  O1B = RCCanon[ExtractOrder[LoddB, 1]]; O1A = RCCanon[ExtractOrder[LoddA, 1]];
  RCSayPlain["F7B_ODD_CONTROL_NONZERO", O1B =!= 0];
  RCSayPlain["F7B_ODD_CONTROL_A_VS_B_VERDICT", RCTry[RCVerdict[O1A, O1B]]];
  RCSayPlain["F7B_ODD_CONTROL_IS_CLEAN_SIGN_FLIP", RCTry[Last@RCNormalize[O1A, O1B, {1, -1}]]];
  RCSay["F7B_ODD_CONTROL_ORDER1", O1B];
];

(* ---- B8: where the sign IS observable -- a background torsion with a half-applied map ---- *)
DefProjectedTensor[T1b[], h, SpaceTimesOfDefinition -> {"Background"}];
DefProjectedTensor[T2b[], h, SpaceTimesOfDefinition -> {"Background"}];
bgTorsion = {Tor[ia_, ib_, ic_] :> (n[ib] h[ia, ic] - n[ic] h[ia, ib]) T1b[] + epsilonh[ib, ic, ia] T2b[]};
bgMinus = bgTorsion /. {T1b[x___] :> -T1b[x], T2b[x___] :> -T2b[x]};
tp = RCTry[RCCanon[ExtractOrder[splitTo[LTalg, Join[rulesMetric, {ruleB}, bgTorsion], h, 1], 1]]];
RCSayPlain["F7B_BG_SPLIT_RAN", tp =!= $Failed && tp =!= Null];
If[tp =!= $Failed && tp =!= Null,
  tpP = RCTry[RCCanon[ExtractOrder[splitTo[LTalg, Join[rulesMetric, {ruleA}, bgTorsion], h, 1], 1]]];
  tpF = RCTry[RCCanon[ExtractOrder[splitTo[LTalg, Join[rulesMetric, {ruleA}, bgMinus], h, 1], 1]]];
  RCSayPlain["F7B_BG_TADPOLE_NONZERO", tp =!= 0];
  RCSay["F7B_BG_ORDER1", tp];
  RCSayPlain["F7B_BG_FLIP_PERTURBATION_ONLY", RCTry[Last@RCNormalize[tpP, tp, {1, -1}]]];
  RCSayPlain["F7B_BG_FLIP_BOTH", RCTry[Last@RCNormalize[tpF, tp, {1, -1}]]];
  (* if the order-1 cross term vanishes for this background and this invariant the two
     verdicts above are vacuous, so try each coupling on its own and say which one carries it *)
  Do[With[{cc = cpl},
     RCSay["F7B_BG_ORDER1_" <> ToString[cc],
       RCTry[RCCanon[ExtractOrder[splitTo[
         Sqrt[-Detg[]] cc (Tor[ia, -ib, -ic] Tor[-ia, ib, ic]),
         Join[rulesMetric, {ruleB}, bgTorsion], h, 1], 1]]]]],
    {cpl, {c1}}];
  RCSay["F7B_BG_ORDER1_TRACE_INVARIANT",
    RCTry[RCCanon[ExtractOrder[splitTo[
      Sqrt[-Detg[]] c3 (Tor[ia, -ib, -ia] Tor[ic, ib, -ic]),
      Join[rulesMetric, {ruleB}, bgTorsion], h, 1], 1]]]];
  RCSay["F7B_BG_ORDER1_ODD_LAGRANGIAN",
    RCTry[RCCanon[ExtractOrder[splitTo[LTodd, Join[rulesMetric, {ruleB}, bgTorsion], h, 1], 1]]]];
];
RCSigns["f7_end"];
RCDone["f7"];
