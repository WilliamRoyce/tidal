#!/usr/bin/env wolframscript
(* cspell:words xMAG xBrauer TraceFree xTras xPand Helpin ToDistortion ToContorsion BreakContorsion BreakDistortion BreakChristoffel StartInducedDecomposition DefConnectionPerturbation ChangeCurvature ChangeCovD CommuteCovDs ChristoffelCDCDT ContorsionCDT TorsionCDT RicciScalarCDT RiemannCDT DistortionCD Kstd VarDt NotebookImport InterpretationBox CellLabel CellContext RCSigns RCNewMessages contortion Riemann-Cartan Detg BYHAND CDNM COVDS Riemanncd XTENSOR Kstd RCSetSigns RCStackOn *)
(* tier1_xmag.wls -- S1 of the R-C follow-up: xMAG called the DOCUMENTED way, judged
   against the AUTHOR'S OWN stored outputs, then the curvature-sign experiment.
   Motivated by PSALTer #543: reproduce the author's example before claiming a defect.
   xMAG is loaded FIRST (as in his In[1]), so this script takes RCSetupCore.wl only and
   never loads xPand.

   Part A -- Tier-1 replay. NotebookImport pulls the Input cells and the stored Output
   cells of Documentation/DefCovD_xMAG.nb with their CellLabels; the targeted labels are
   replayed in the author's order with the author's names, and each stored Out[k] is
   compared with ours by RCVerdict. A mismatch in a SETUP cell means API drift against
   xTensor 1.3.0 (his kernel was 1.2.0) and stops every behavioral claim.

   Part B -- the sign experiment. xMAG sets the GLOBAL $RiemannSign = -1 at load
   (xMAG.m:110, undocumented; xTensor's default is +1, xTensor.m:1837) and both its own
   changeRiemann override (xMAG.m:1258) and xTensor's (xTensor.m:6194) read it at call
   time. R-C's "xMAG's torsion part has the opposite sign" compared a -1 result with a
   +1 reference computed in another kernel. Arms A1-A8 settle that, with controls.

   Part C -- the three calls R-C got wrong: 2-argument BreakDistortion (xMAG.m:1664),
   single-character derivative symbols for StartInducedDecomposition (xTensor.m:6591),
   and VarDt for the connection variation (xMAG.m:1428). Each with the failing form as
   its own control. *)

(* xMAG FIRST, as in the author's In[1] -- and before RCSetupCore.wl, whose short xAct
   names must resolve into the xAct contexts (see the guard in that file). *)
Print["RC_SIGNS_BEFORE_ANY_NEEDS=xTensor not loaded"];
$T1T0 = AbsoluteTime[];
Needs["xAct`xMAG`"];                                      (* the author's In[1] *)
$T1Failed = !MemberQ[$Packages, "xAct`xMAG`"];            (* Check would trip on benign messages *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]];
RCSayPlain["T1_XMAG_LOAD_FAILED", $T1Failed];
RCSayPlain["T1_XMAG_LOAD_WALL_S", Round[AbsoluteTime[] - $T1T0, 0.1]];
RCSay["T1_XMAG_VERSION", xAct`xMAG`$Version];
RCSay["T1_XTENSOR_VERSION", xAct`xTensor`$Version];
RCSay["T1_XTRAS_VERSION", xAct`xTras`$Version];
RCSay["T1_XBRAUER_VERSION", xAct`xBrauer`$Version];
RCSay["T1_TRACEFREE_VERSION", xAct`TraceFree`$Version];
RCSigns["after_xmag_load"];                               (* expect {-1, 1, 1, 1, 1} *)
If[$T1Failed, RCSayPlain["T1_STOP", "xMAG did not load"]; RCDone["t1"]; Quit[]];
$XMAGDoc = FileNameJoin[{$UserBaseDirectory, "Applications", "xAct", "xMAG", "Documentation"}];
RCSayPlain["T1_DOC_DIR_EXISTS", DirectoryQ[$XMAGDoc]];

(* ================= Part A: Tier-1 replay of the author's own cells ================= *)
(* Read the notebook as expressions: Input cells carry the code, Output cells the stored
   result inside InterpretationBox[boxes, expr]; CellLabel gives In[k]/Out[k]. *)
$nbFile = FileNameJoin[{$XMAGDoc, "DefCovD_xMAG.nb"}];
nb = RCTry[Get[$nbFile]];
RCSayPlain["T1_NOTEBOOK_READ", nb =!= $Failed];
labelOf[c_] := With[{l = Cases[c, HoldPattern[CellLabel -> s_String] :> s, {0, Infinity}]},
  If[l === {}, "", First[l]]];
numOf[lab_String] := ToExpression[StringJoin[Select[Characters[lab], DigitQ]]];
cells = Cases[nb, Cell[__], {0, Infinity}];
inCells = Select[cells, MatchQ[#, Cell[_, "Input", ___]] && StringStartsQ[labelOf[#], "In["] &];
outCells = Select[cells, MatchQ[#, Cell[_, "Output", ___]] && StringStartsQ[labelOf[#], "Out["] &];
printCells = Select[cells, MatchQ[#, Cell[_, "Print", ___]] &];
RCSayPlain["T1_INPUT_CELL_COUNT", Length[inCells]];
RCSayPlain["T1_OUTPUT_CELL_COUNT", Length[outCells]];
RCSay["T1_INPUT_LABELS", labelOf /@ inCells];

(* The author's kernel wrote its symbols in $CellContext`; ours are in Global`. Map the
   stored expression's symbol names into this session before comparing anything. *)
toHere[HoldComplete[e_]] := HoldComplete[e] /. s_Symbol :>
  With[{n = SymbolName[Unevaluated[s]]},
    If[Context[Unevaluated[s]] === "$CellContext`", Symbol["Global`" <> n], s]];

(* Both Input and Output cells are read with ToExpression on their boxes: xAct's
   InterpretationBox carries the real expression, and ToExpression honours it (taking the
   first InterpretationBox by hand picks up nested sub-boxes -- first run's Out[2] = 1). *)
outExpr[c_] := RCTry[toHere[ToExpression[First[c] /. BoxData[x_] :> x, StandardForm, HoldComplete]]];
storedOut = Association @@ (numOf[labelOf[#]] -> outExpr[#] & /@ outCells);
RCSay["T1_STORED_OUTPUT_LABELS", Sort@Keys[storedOut]];
RCSayPlain["T1_STORED_OUTPUTS_READABLE", Count[Values[storedOut], _HoldComplete]];

(* An Input cell may hold several lines separated by \[IndentingNewLine]; each is one
   evaluation whose result is the next Out[k]. *)
inLines[c_] := Module[{boxes = First[c] /. BoxData[x_] :> x, frags},
  boxes = If[Head[boxes] === List, boxes, {boxes}];
  frags = DeleteCases[boxes, "\[IndentingNewLine]" | "\n" | " "];
  RCTry[ToExpression[#, StandardForm, HoldComplete] & /@ frags]];

(* the author's targeted cells: setup, the contortion chain, Riemann, the general
   connection, the 2-argument BreakDistortion, and the Ricci scalar *)
targets = {1, 2, 3, 4, 10, 13, 23, 46, 48, 57, 59, 62};
$T1Informational = {2};   (* Options@DefMetric: a list of held rules, not an expression to compare *)
$lastOut = Null;
runLine[lab_String, held_, k_Integer] := Module[{res},
  res = RCNewMessages[lab, TimeConstrained[ReleaseHold[held /. HoldPattern[Out[]] :> $lastOut], 120, $Aborted]];
  $lastOut = res;
  If[KeyExistsQ[storedOut, k] && MatchQ[storedOut[k], _HoldComplete],
    Module[{stored = ReleaseHold[storedOut[k]], v},
      If[MemberQ[$T1Informational, k],
        RCSayPlain["T1_OUT" <> ToString[k] <> "_INFORMATIONAL", "option/usage list, not compared"],
        v = RCTry[RCVerdict[res, stored]];
        RCSayPlain["T1_OUT" <> ToString[k] <> "_VERDICT", v];
        RCSay["T1_OUT" <> ToString[k] <> "_OURS", res];
        RCSay["T1_OUT" <> ToString[k] <> "_STORED", stored];
        (* a verdict other than identical is usually only the dummy-index names: the
           author's declared a1 vs our generated b4$nnn. SameDummies + ScreenDollarIndices
           settles it, and the canonical difference is printed either way. *)
        If[!StringStartsQ[ToString[v], "identical"],
          RCSayPlain["T1_OUT" <> ToString[k] <> "_VERDICT_SAMEDUMMIES",
            RCTry[RCVerdict[SameDummies[res], SameDummies[stored]]]];
          RCSay["T1_OUT" <> ToString[k] <> "_CANONICAL_DIFFERENCE",
            RCTry[ScreenDollarIndices[RCCanon[SameDummies[res - stored]]]]]]]],
    RCSay["T1_OUT" <> ToString[k] <> "_OURS_NO_STORED", res]];
  res];
Do[Module[{c = SelectFirst[inCells, numOf[labelOf[#]] === k &], lines, outK = k},
   If[MissingQ[c], RCSayPlain["T1_IN" <> ToString[k] <> "_MISSING", True],
     lines = inLines[c];
     If[lines === $Failed || MemberQ[lines, $Failed],
       RCSayPlain["T1_IN" <> ToString[k] <> "_UNPARSED", True],
       Do[runLine["in" <> ToString[k] <> "_line" <> ToString[i], lines[[i]], outK];
          outK = outK + 1, {i, Length[lines]}]]]], {k, targets}];
RCSay["T1_TENSORS_AFTER_REPLAY", Select[$Tensors, StringContainsQ[ToString[#], "CDT"] || StringContainsQ[ToString[#], "CDNM"] &]];
RCSay["T1_COVDS_AFTER_REPLAY", $CovDs];
RCSigns["after_replay"];
RCSayPlain["T1_PART_A_DONE", True];

(* ================= Part B: the curvature-sign experiment ================= *)
(* The author's g, cd, CDT (metric-compatible, with torsion) and CD (general) are already
   defined by the replay. Two routes to R-tilde, both in THIS kernel:
     xmag[]  = BreakContorsion @ ToDistortion[RicciScalarCDT[], CDT, g]     (the package)
     byHand[] = ChangeCurvature[RicciScalarCDT[], CDT, cd] with the contortion identity
                substituted by hand (xTensor only).
   NOTE why the first attempt (t1/20260917T100425Z) compared the wrong thing: with xMAG
   loaded, ChristoffelcdCDT is ALREADY rewritten to -ContorsionCDT by its automatic rule
   (xMAG.m:1311-1316), so a substitution keyed on ChristoffelcdCDT never fires and the
   by-hand route kept unbroken contortion. Key the substitution on ContorsionCDT. *)
DefTensor[vv[i1], M];
Kstd[a_, b_, c_] := 1/2 (TorsionCDT[a, b, c] + TorsionCDT[b, a, c] + TorsionCDT[c, a, b]);
xmag[] := RCCanon[ContractMetric[NoScalar[BreakContorsion[ToDistortion[RicciScalarCDT[], CDT, g], CDT]]]];
byHand[] := RCCanon[ContractMetric[NoScalar[Expand[
   ChangeCurvature[RicciScalarCDT[], CDT, cd] /. ContorsionCDT[a1_, a2_, a3_] :> Kstd[a1, a2, a3]]]]];
RCSay["T1_CHRISTOFFEL_NAMES_PRESENT", Select[$Tensors, StringContainsQ[ToString[#], "Christoffel"] &]];
RCSay["T1_CHRISTOFFELCDCDT_AUTORULE", RCTry[ChristoffelcdCDT[i1, -i2, -i3]]];
RCSay["T1_CDT_MINUS_CD_ON_A_VECTOR", RCTry[ToCanonical[ChangeCovD[CDT[-i1]@vv[i2], CDT, cd] - cd[-i1]@vv[i2]]]];
RCSayPlain["T1_CDT_MINUS_CD_IS_PLUS_CONTORSION", RCTry[
  RCVerdict[ToCanonical[ChangeCovD[CDT[-i1]@vv[i2], CDT, cd] - cd[-i1]@vv[i2]],
            ContorsionCDT[i2, -i1, -i3] vv[i3]]]];

(* A1: both routes under xMAG's own convention (the comparison Probe E never made) *)
RCSigns["arm_A1"];
x1 = RCTry[xmag[]]; p1 = RCTry[byHand[]];
RCSay["T1_A1_XMAG_MINUS1", x1];
RCSay["T1_A1_BYHAND_MINUS1", p1];
RCSayPlain["T1_A1_XMAG_VS_BYHAND_MINUS1", RCTry[RCVerdict[x1, p1]]];

(* A2: the project's convention asserted; both routes again *)
RCSetSigns[{1, 1, 1}];
x3 = RCTry[xmag[]]; p3 = RCTry[byHand[]];
RCSay["T1_A2_XMAG_PLUS1", x3];
RCSayPlain["T1_A2_XMAG_VS_BYHAND_PLUS1", RCTry[RCVerdict[x3, p3]]];

(* A3: the cross-convention control -- this is the comparison R-C actually made *)
RCSayPlain["T1_A3_CONTROL_CROSS_CONVENTION", RCTry[RCVerdict[x1, p3]]];
RCSayPlain["T1_A3_TORSION_PART_IS_NEGATED", RCTry[
  RCVerdict[RCCanon[x1 - RicciScalarcd[]], RCCanon[-(p3 - RicciScalarcd[])]]]];
RCSay["T1_A3_TORSION_PART_MINUS1", RCTry[RCCanon[x1 - RicciScalarcd[]]]];
RCSay["T1_A3_TORSION_PART_PLUS1", RCTry[RCCanon[x3 - RicciScalarcd[]]]];

(* A4: the textbook Einstein-Cartan check with vectorial torsion, both conventions.
   T^a_bc = delta^a_c v_b - delta^a_b v_c gives, in dim dimensions and with xTensor's
   signs, R-tilde - R = s_R (-(dim-1)(dim-2) v.v + 2(dim-1) div v); at dim = 4 that is
   the textbook -6 v^2 + 6 div v. *)
(* xMAG's ConnectionRelations rewrite the trace T^a_{ba} into TorsionVectorCDT
   (xMAG.m:1038), so the vectorial-torsion substitution must cover the trace too:
   T^a_{ba} = (dim - 1) v_b for this ansatz. *)
vecRule = {TorsionCDT[a1_, a2_, a3_] :> delta[a1, a3] vv[a2] - delta[a1, a2] vv[a3],
           TorsionVectorCDT[a1_] :> (dim - 1) vv[a1]};
vecTarget[s_] := s (-(dim - 1) (dim - 2) Scalar[vv[i1] vv[-i1]] + 2 (dim - 1) cd[-i1]@vv[i1]);
Do[With[{sg = sgn}, RCSetSigns[{sg, 1, 1}];
   Module[{got = RCTry[RCCanon[ContractMetric[((xmag[] /. vecRule) - RicciScalarcd[])]]]},
     RCSay["T1_A4_VECTOR_TORSION_" <> If[sg === 1, "PLUS1", "MINUS1"], got];
     RCSayPlain["T1_A4_VECTOR_TORSION_VS_TEXTBOOK_" <> If[sg === 1, "PLUS1", "MINUS1"],
       RCTry[RCVerdict[got, RCCanon[vecTarget[sg]]]]];
     RCSayPlain["T1_A4_AT_DIM4_" <> If[sg === 1, "PLUS1", "MINUS1"],
       RCTry[RCVerdict[got /. dim -> 4, RCCanon[sg (-6 Scalar[vv[i1] vv[-i1]] + 6 cd[-i1]@vv[i1])] /. dim -> 4]]]]],
  {sgn, {1, -1}}];

(* A5: the SCALAR commutator, which is exactly how xTensor defines $TorsionSign
   (xTensor.m:289): cd[-a]@cd[-b]@f - cd[-b]@cd[-a]@f = -$TorsionSign T^c_ab cd[-c]@f.
   No contortion is involved, so this is an independent check of the torsion convention. *)
RCSetSigns[{1, 1, 1}];
DefTensor[ff[], M];
(* Invar turns $CommuteCovDsOnScalars off when it loads (visible in the load banner), and
   with it off xTensor never applies its own scalar-commutator rule (xTensor.m:6890), so
   the difference stays unsorted and no verdict is possible. Turn it on for this arm. *)
RCSayPlain["T1_A5_COMMUTE_ON_SCALARS_AS_LOADED", $CommuteCovDsOnScalars];
$CommuteCovDsOnScalars = True;
scalarComm[] := ToCanonical[CDT[-i1]@CDT[-i2]@ff[] - CDT[-i2]@CDT[-i1]@ff[]];
RCSay["T1_A5_SCALAR_COMMUTATOR", RCTry[scalarComm[]]];
RCSayPlain["T1_A5_SCALAR_COMMUTATOR_VS_MINUS_T_GRAD_F", RCTry[
  RCVerdict[scalarComm[], ToCanonical[-TorsionCDT[i3, -i1, -i2] CDT[-i3]@ff[]]]]];
RCSayPlain["T1_A5_CONTROL_PLUS_T_GRAD_F", RCTry[
  RCVerdict[scalarComm[], ToCanonical[TorsionCDT[i3, -i1, -i2] CDT[-i3]@ff[]]]]];
RCSayPlain["T1_A5_CONTROL_WRONG_INDEX_ORDER_IS_INERT", RCTry[
  ToCanonical[CommuteCovDs[CDT[-i2]@CDT[-i3]@vv[i1], CDT, {-i2, -i3}] - CDT[-i2]@CDT[-i3]@vv[i1]] === 0]];

(* A6: the identity text and its two defining properties, for conventions.md *)
RCSay["T1_A6_KSTD_IDENTITY_INPUTFORM", HoldForm[Kstd[a, b, c]] -> Kstd[a, b, c]];
RCSay["T1_A6_CONTORSION_BROKEN", RCTry[BreakContorsion[ContorsionCDT[i1, -i2, -i3], CDT]]];
RCSayPlain["T1_A6_KSTD_EQUALS_XMAG_CONTORSION", RCTry[
  RCVerdict[RCCanon[Kstd[i1, -i2, -i3]], RCCanon[BreakContorsion[ContorsionCDT[i1, -i2, -i3], CDT]]]]];
RCSay["T1_A6_KSTD_ANTISYM_MINUS_T", RCTry[RCCanon[Kstd[i1, -i2, -i3] - Kstd[i1, -i3, -i2] - TorsionCDT[i1, -i2, -i3]]]];
RCSay["T1_A6_KSTD_METRICITY", RCTry[RCCanon[Kstd[-i1, -i2, -i3] + Kstd[-i3, -i2, -i1]]]];

(* A7: with the project's convention asserted, does xMAG still reproduce the author's
   Riemann cell (Out[24], stored under -1)? The torsion part must be negated. *)
RCSetSigns[{1, 1, 1}];
With[{stored = If[MatchQ[storedOut[24], _HoldComplete], ReleaseHold[storedOut[24]], $Failed]},
  If[stored =!= $Failed,
    Module[{ours = RCTry[ToContorsion[RiemannCDT[-i1, -i2, -i3, i4]]]},
      RCSayPlain["T1_A7_OUT24_UNDER_PLUS1_VS_STORED", RCTry[RCVerdict[ours, stored]]];
      RCSayPlain["T1_A7_OUT24_TORSION_PART_NEGATED", RCTry[
        RCVerdict[RCCanon[ours - Riemanncd[-i1, -i2, -i3, i4]],
                  RCCanon[-(stored - Riemanncd[-i1, -i2, -i3, i4])]]]];
      RCSay["T1_A7_OUT24_UNDER_PLUS1", ours]]]];

(* A8: the family-A -> family-B import maps that Probe D3's rules need.
   xTensor's torsion is twice the antisymmetric part of the connection difference in its
   LAST TWO slots (the derivative index is the middle one); a family-A paper (Shapiro,
   Aoki et al., Nikiforova-Damour) antisymmetrizes the other pair, so its T is minus ours. *)
RCSay["T1_A8_TORSION_FROM_CHRISTOFFEL", RCTry[ChangeTorsion[TorsionCDT[i1, -i2, -i3], CDT, cd]]];
RCSayPlain["T1_A8_TORSION_IS_ANTISYM_IN_LAST_TWO", RCTry[
  RCVerdict[RCCanon[TorsionCDT[i1, -i2, -i3] + TorsionCDT[i1, -i3, -i2]], 0]]];
RCSayPlain["T1_A8_CONTORSION_ANTISYM_PAIR", RCTry[
  {RCVerdict[RCCanon[ContorsionCDT[-i1, -i2, -i3] + ContorsionCDT[-i3, -i2, -i1]], 0],
   RCVerdict[RCCanon[ContorsionCDT[-i1, -i2, -i3] + ContorsionCDT[-i2, -i1, -i3]], 0]}]];
RCSayPlain["T1_PART_B_DONE", True];

(* ================= Part C: the three calls R-C got wrong ================= *)
(* C1: BreakDistortion, two arguments (xMAG.m:1664 dispatches to BreakContorsion) *)
rt = RCTry[ToDistortion[RicciScalarCDT[], CDT, g]];
b2 = RCNewMessages["break2", RCTry[BreakDistortion[rt, CDT]]];
RCSayPlain["T1_C1_BREAKDISTORTION_2ARG_RAN", b2 =!= $Failed && b2 =!= Null];
RCSayPlain["T1_C1_BREAKDISTORTION_2ARG_EQUALS_BREAKCONTORSION", RCTry[
  RCVerdict[RCCanon[b2], RCCanon[BreakContorsion[rt, CDT]]]]];
b3 = RCNewMessages["break3", RCTry[BreakDistortion[rt, CDT, g]]];
RCSayPlain["T1_C1_CONTROL_BREAKDISTORTION_3ARG", If[b3 === $Failed || b3 === Null, "Null/failed (R-C's error reproduced)", "ran"]];

(* C2: StartInducedDecomposition. Its own usage says it "works only for a non-Riemannian
   geometry with independent metric met and GL(dim) independent connection", and the
   author's In[175] calls it on his GENERAL connection -- R-C called it on a
   metric-compatible one with two-character derivative symbols (xTensor.m:6591 wants one). *)
RCSay["T1_C2_USAGE", RCUsage["xAct`xMAG`", "StartInducedDecomposition"]];
RCSay["T1_C2_OPTIONS", Options[StartInducedDecomposition]];
c2 = RCNewMessages["startinduced_general",
  RCTry[(StartInducedDecomposition[g, CD, {{"~", "Dh"}, {"^", "DGh"}}, {nv, hh}]; True)]];
RCSayPlain["T1_C2_GENERAL_CONNECTION_OK", c2 === True];
RCSigns["after_startinduced"];                            (* xMAG sets two more to -1 *)
If[c2 === True,
  RCSay["T1_C2_TENSORS_DEFINED", Select[$Tensors, StringContainsQ[ToString[#], "hh"] || StringContainsQ[ToString[#], "nv"] &]];
  RCSay["T1_C2_COVDS", $CovDs];
  RCSay["T1_C2_NORM_OF_NORMAL", RCTry[ToCanonical[nv[i1] nv[-i1]]]];
  RCSay["T1_C2_INDUCED_TRACE", RCTry[ToCanonical[hh[i1, -i1]]]];
  RCSay["T1_C2_INDUCED_DECOMPOSITION_OF_A_VECTOR", RCTry[InducedDecomposition[vv[i1], {hh, nv}]]];
  RCSay["T1_C2_RIEMANN_OF_INDUCED_CONNECTION", RCTry[ToDistortion[RiemannCDhh[-i1, -i2, -i3, i4]]]];
];
c2bad = RCNewMessages["startinduced_twochar",
  RCTry[(StartInducedDecomposition[g, CD, {{"~~", "Dk"}, {"^^", "DGk"}}, {nk, hk}]; True)]];
RCSayPlain["T1_C2_CONTROL_TWO_CHAR_SYMBOL", If[c2bad === True, "ran (unexpected)", "refused (R-C's error reproduced)"]];

(* C3: the connection variation. The author varies with respect to the Christoffel with
   xTensor's VarD (Tutorial_Examples_xMAG.nb In[300]: VarD[ChristoffelCD[-i1,i2,i3],cd][L]);
   VarDt is internal to DefConnectionVariation. R-C called DefConnectionPerturbation a
   second time (it already ran inside DefCovD, xMAG.m:1079) and then read
   Perturbation[TorsionCDT] -- which xMAG never defines (it perturbs the connection only). *)
RCSay["T1_C3_PERTURBATION_TENSORS", Select[$Tensors, StringContainsQ[ToString[#], "Perturbation"] &]];
RCSay["T1_C3_DEFCONNECTIONPERTURBATION_USAGE", RCUsage["xAct`xMAG`", "DefConnectionPerturbation"]];
v1 = RCNewMessages["vard_connection",
  RCTry[VarD[ChristoffelCDT[-i1, i2, i3], cd][Sqrt[-Detg[]] RicciScalarCDT[]]]];
RCSayPlain["T1_C3_VARD_CONNECTION_RAN", v1 =!= $Failed && v1 =!= Null];
RCSay["T1_C3_VARD_CONNECTION_RESULT", v1];
RCSay["T1_C3_VARD_CONNECTION_HEADS", If[v1 === $Failed || v1 === Null, "none",
  Union@Cases[v1, s_Symbol[___] /; xTensorQ[s] :> s, {0, Infinity}]]];
(* and the same for the general connection, the case the author shows *)
v2 = RCNewMessages["vard_general",
  RCTry[VarD[ChristoffelCD[-i1, i2, i3], cd][Sqrt[-Detg[]] RicciScalarCD[]]]];
RCSayPlain["T1_C3_VARD_GENERAL_RAN", v2 =!= $Failed && v2 =!= Null];
RCSay["T1_C3_VARD_GENERAL_RESULT", RCTry[ToDistortion[v2]]];
RCSayPlain["T1_PART_C_DONE", True];
RCSigns["at_end"];
RCDone["t1"];
