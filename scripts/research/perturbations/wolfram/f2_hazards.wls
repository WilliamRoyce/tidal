#!/usr/bin/env wolframscript
(* cspell:words xPand xMAG ToxPandFromRules SplitMetric SplitPerturbations ToMetric ToDistortion BreakContorsion EinsteinToRicci FrozenMetricQ MetricOfCovD InducedFrom UndefCovD UndefChart DefChart PDBc xCoba RCStackOn RCSplitGuarded Kstd RCSigns normu CDCDL CDCDM DEFCOVD COVDS MASTEROF dedupe TRANSVERSALITY *)
(* f2_hazards.wls -- S2 of the R-C follow-up: the two "xPand breaks" cases traced to their
   lines, told apart by a twin control, and the caller-side workarounds measured.
   R-C reported one hazard; the sources say there are two different sites and only one is
   xPand's fault:
     (i)  xPand.m:2566 (ToMetric, reached from every ToxPandFromRules through Conformal at
          :2619) selects over Rest@$CovDs and asks InducedFrom[MetricOfCovD[#]] -- which
          throws for a derivative with no metric at all (xCoba's chart derivative PDBc).
          xTensor guards the same pattern at xTensor.m:8951. THIS IS XPAND'S.
     (ii) with xMAG loaded, its replacement of EinsteinToRicci (xMAG.m:1154) tests
          FrozenMetricQ[MasterOf[#]] and FrozenMetricQ[Null] (xMAG.m:1130) reaches the same
          InducedFrom[Null] fallback (xTensor.m:8302) -- but only for a connection declared
          WITHOUT Master, which xMAG's own ?Master text says is not how to declare one.
          THIS WAS OUR CONFIGURATION: the CDM twin control below has Master -> g and must
          leave the split intact.
   Also here: the decisive control for sharing one kernel with xMAG (xMAG sets the global
   $RiemannSign = -1 at load, so an xPand split in that kernel is in the opposite curvature
   convention), and the xCoba component count with the dedupe R-C got wrong (48, not 24). *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];
DefConstantSymbol[kap];
(* The reference must be sign-sensitive and nonzero: a TT-only first-order Ricci scalar
   vanishes identically (first run of this script printed 0), so use the Newtonian-gauge
   SCALAR sector -- the same expression the paper's minimal example prints. *)
rulesMetric = SplitMetric[g, dg, h, "NewtonGauge"];
split[] := RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]];
splitOK[r_] := If[r === $Failed || r === Null, "BROKEN", "ok"];
R0 = split[];
RCSayPlain["F2_REFERENCE_SPLIT", splitOK[R0]];
RCSay["F2_REFERENCE_SPLIT_TERMS", If[R0 === $Failed || R0 === Null, 0, Length[List @@ Expand[NoScalar[R0]]]]];
RCSay["F2_COVDS_AT_START", $CovDs];
(* elementary probe of the simplifications xPand installs with DefMetricFields: the
   divergence of the transverse vector field and of the TT tensor field must vanish *)
transProbe[] := RCTry[{ToCanonical[ContractMetric[h[ia, ib] cd[-ib]@Evh[LI[1], LI[0], -ia]]],
                       ToCanonical[ContractMetric[h[ia, ib] cd[-ib]@Eth[LI[1], LI[0], -ia, -ic]]],
                       ToCanonical[ContractMetric[h[ia, ib] Eth[LI[1], LI[0], -ia, -ib]]]}];
RCSay["F2_TRANSVERSALITY_PROBE_AT_START", transProbe[]];
RCSay["F2_METRIC_OF_EACH_COVD_AT_START", MetricOfCovD /@ $CovDs];

(* ---------------- case (i): xPand's own hazard, a chart derivative ---------------- *)
RCStackOn[{HoldForm[InducedFrom::unknown]}];
Needs["xAct`xCoba`"];
ch = RCNewMessages["defchart", RCTry[(DefChart[Bc, M4, {0, 1, 2, 3}, {tt[], xx[], yy[], zz[]}]; True)]];
RCSayPlain["F2_I_CHART_DEFINED", ch === True];
RCSay["F2_I_COVDS_AFTER_CHART", $CovDs];
RCSay["F2_I_METRIC_OF_EACH_COVD_AFTER_CHART", MetricOfCovD /@ $CovDs];
rChart = RCNewMessages["split_after_chart", split[]];
RCSayPlain["F2_I_SPLIT_AFTER_CHART", splitOK[rChart]];
(* the guarded caller: drop the metric-less derivatives for the duration of the split.
   $CovDs is an ordinary global (xTensor.m:830), and this is exactly the one-line guard
   xPand.m:2566 is missing -- emulated on our side, with xPand unmodified. *)
RCSplitGuarded[] := Block[{$CovDs = Select[$CovDs, # === PD || MetricOfCovD[#] =!= Null &]},
  RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]]];
(* the stricter variant: xMAG's own predicates need a Master as well as a metric, so a
   wrapper that also drops Master-less derivatives covers both hazards at once *)
RCSplitGuardedStrict[] := Block[{$CovDs = Select[$CovDs,
    # === PD || (MetricOfCovD[#] =!= Null && MasterOf[#] =!= Null) &]},
  RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]]];
rGuard = RCNewMessages["guarded_split_chart", RCSplitGuarded[]];
RCSayPlain["F2_I_GUARDED_SPLIT_AFTER_CHART", splitOK[rGuard]];
RCSayPlain["F2_I_GUARDED_SPLIT_EQUALS_REFERENCE", If[rGuard === $Failed || rGuard === Null, "n/a", RCTry[RCVerdict[rGuard, R0]]]];
(* negative control: dropping an UNRELATED derivative must not fix it *)
RCSplitWrongGuard[] := Block[{$CovDs = DeleteCases[$CovDs, CDah2]},
  RCTry[RCSort[ExtractOrder[ToxPandFromRules[RicciScalarCD[], rulesMetric, h, 1], 1]]]];
RCSayPlain["F2_I_CONTROL_WRONG_GUARD", splitOK[RCNewMessages["wrong_guard", RCSplitWrongGuard[]]]];
(* the xCoba component count of a rank-3 antisymmetric torsion, deduped properly:
   R-C compared ToCanonical'd component expressions and got 48 by counting each
   antisymmetric pair twice; keying on the sorted index pair gives the 24 of the papers. *)
DefTensor[Tor[ia, -ib, -ic], M4, Antisymmetric[{-ib, -ic}]];
comps = Flatten[Table[{{i, j, k}, ToCanonical@Tor[{i, Bc}, {-j, -Bc}, {-k, -Bc}]}, {i, 0, 3}, {j, 0, 3}, {k, 0, 3}], 2];
nonzero = Select[comps, Last[#] =!= 0 &];
RCSayPlain["F2_I_XCOBA_NONZERO_COMPONENTS", Length[nonzero]];
RCSayPlain["F2_I_XCOBA_INDEPENDENT_COMPONENTS",
  Length@DeleteDuplicates[{First[#][[1]], Sort[Rest[First[#]]]} & /@ nonzero]];
RCSayPlain["F2_I_COUNT_IS_24", Length@DeleteDuplicates[{First[#][[1]], Sort[Rest[First[#]]]} & /@ nonzero] === 24];
(* undefine the chart and show the split recovers without any guard *)
RCSayPlain["F2_I_UNDEFCHART_RAN", RCTry[(UndefChart[Bc]; True)] === True];
RCSayPlain["F2_I_SPLIT_AFTER_UNDEFCHART", splitOK[RCNewMessages["split_after_undefchart", split[]]]];
RCStackOff[];

(* ---------------- case (ii): xMAG, with and without Master ---------------- *)
RCSigns["before_xmag"];
Needs["xAct`xMAG`"];
RCSigns["after_xmag_load"];                     (* xMAG.m:110 flips $RiemannSign to -1 *)
RCSay["F2_TRANSVERSALITY_PROBE_AFTER_XMAG", transProbe[]];
RCSayPlain["F2_II_SPLIT_AFTER_XMAG_LOAD_ONLY", splitOK[RCNewMessages["split_after_xmag_load", split[]]]];
RCSayPlain["F2_II_SPLIT_AFTER_XMAG_LOAD_VS_REFERENCE_UNDER_MINUS1",
  RCTry[RCVerdict[RCNewMessages["split_after_xmag_minus1", split[]], R0]]];
(* case (ii) is about Master, not about the sign: put the project's convention back so the
   comparisons below isolate the connection declaration *)
RCSetSigns[{1, 1, 1}];
RCSayPlain["F2_II_SPLIT_AFTER_SIGN_RESET_VS_REFERENCE",
  RCTry[RCVerdict[RCNewMessages["split_after_reset", split[]], R0]]];
(* the twin: the SAME connection declared the documented way (Master -> g) and the way
   R-C declared it (no Master). Each is defined in turn and the split retested. *)
okM = RCNewMessages["defcovd_master",
  RCTry[(DefCovD[CDM[-ia], {"#", "DM"}, Torsion -> True, FromMetric -> g, Master -> g, ConnectionRelations -> True]; True)]];
RCSayPlain["F2_II_DEFCOVD_WITH_MASTER_OK", okM === True];
RCSay["F2_II_MASTEROF_CDM", RCTry[MasterOf[CDM]]];
rM = RCNewMessages["split_with_master", split[]];
RCSayPlain["F2_II_SPLIT_WITH_MASTERFUL_CONNECTION", splitOK[rM]];
RCSayPlain["F2_II_TODISTORTION_WITH_MASTER_RAN", RCTry[FreeQ[BreakContorsion[ToDistortion[RicciScalarCDM[], CDM, g], CDM], $Failed]]];
RCStackOn[{HoldForm[InducedFrom::unknown], HoldForm[Validate::inhom]}];
okL = RCNewMessages["defcovd_nomaster",
  RCTry[(DefCovD[CDL[-ia], {"&", "DL"}, Torsion -> True, FromMetric -> g]; True)]];
RCSayPlain["F2_II_DEFCOVD_WITHOUT_MASTER_OK", okL === True];
RCSay["F2_II_MASTEROF_CDL", RCTry[MasterOf[CDL]]];
RCSay["F2_II_METRIC_OF_EACH_COVD_NOW", MetricOfCovD /@ $CovDs];
RCSay["F2_II_MASTER_OF_EACH_COVD_NOW", MasterOf /@ $CovDs];
(* the single-call witness: xMAG's replacement of EinsteinToRicci folds over $CovDs *)
RCSayPlain["F2_II_EINSTEINTORICCI_ONE_ARG",
  splitOK[RCNewMessages["einsteintoricci_1arg", RCTry[EinsteinToRicci[RicciScalarCD[]]]]]];
RCSayPlain["F2_II_EINSTEINTORICCI_TWO_ARG",
  splitOK[RCNewMessages["einsteintoricci_2arg", RCTry[EinsteinToRicci[RicciScalarCD[], CD]]]]];
rL = RCNewMessages["split_without_master", split[]];
RCSayPlain["F2_II_SPLIT_WITH_MASTERLESS_CONNECTION", splitOK[rL]];
RCSayPlain["F2_II_TODISTORTION_ON_CDM_AFTER_CDL", RCTry[
  If[FreeQ[BreakContorsion[ToDistortion[RicciScalarCDM[], CDM, g], CDM], $Failed], "ran", "failed"]]];
rLG = RCNewMessages["guarded_split_masterless", RCSplitGuarded[]];
RCSayPlain["F2_II_GUARDED_SPLIT_WITH_MASTERLESS", splitOK[rLG]];
rLGS = RCNewMessages["strict_guard_masterless", RCSplitGuardedStrict[]];
RCSayPlain["F2_II_STRICT_GUARDED_SPLIT_WITH_MASTERLESS", splitOK[rLGS]];
RCSayPlain["F2_II_STRICT_GUARD_EQUALS_REFERENCE", If[rLGS === $Failed || rLGS === Null, "n/a", RCTry[RCVerdict[rLGS, R0]]]];
RCSayPlain["F2_II_UNDEFCOVD_CDL_RAN", RCTry[(UndefCovD[CDL]; True)] === True];
RCSayPlain["F2_II_SPLIT_AFTER_UNDEFCOVD", splitOK[RCNewMessages["split_after_undefcovd", split[]]]];
RCStackOff[];

(* ---------------- section III: can xMAG share xPand's kernel? Three separable causes:
     the sign global, the extra connection objects, and xMAG's replacements of xTensor's
     curvature-relation functions. Each is switched one at a time against the reference. *)
RCSigns["section_III_start"];
RCSayPlain["F2_III_COVDS_NOW", ToString[$CovDs]];
(* (a) sign only: xMAG loaded, CDM present, convention as xMAG left it *)
RCSetSigns[{-1, 1, 1}];
rMinus = RCNewMessages["III_minus1", split[]];
RCSayPlain["F2_III_A_SPLIT_UNDER_MINUS1", splitOK[rMinus]];
RCSayPlain["F2_III_A_VS_REFERENCE", If[rMinus === $Failed || rMinus === Null, "n/a", RCTry[RCVerdict[rMinus, R0]]]];
RCSayPlain["F2_III_A_IS_CLEAN_OVERALL_FLIP", If[rMinus === $Failed || rMinus === Null, "n/a",
  RCTry[Last@RCNormalize[rMinus, R0, {1, -1}]]]];
(* (b) convention restored, everything else unchanged *)
RCSetSigns[{1, 1, 1}];
rPlus = RCNewMessages["III_plus1", split[]];
RCSayPlain["F2_III_B_SPLIT_AFTER_SIGN_RESET", splitOK[rPlus]];
RCSayPlain["F2_III_B_VS_REFERENCE", If[rPlus === $Failed || rPlus === Null, "n/a", RCTry[RCVerdict[rPlus, R0]]]];
RCSay["F2_III_B_DIFFERENCE_FROM_REFERENCE", If[rPlus === $Failed || rPlus === Null, "n/a", RCTry[RCCanon[rPlus - R0]]]];
(* is the residual confined to the vector and tensor fields (whose transversality and
   tracelessness rules are what stopped firing), leaving the scalar sector intact? *)
killVT = {Bvh[__] -> 0, Evh[__] -> 0, Eth[__] -> 0};
RCSayPlain["F2_III_B_SCALAR_SECTOR_VS_REFERENCE", If[rPlus === $Failed || rPlus === Null, "n/a",
  RCTry[RCVerdict[RCCanon[rPlus /. killVT], RCCanon[R0 /. killVT]]]]];
RCSayPlain["F2_III_B_RESIDUAL_IS_VECTOR_TENSOR_ONLY", If[rPlus === $Failed || rPlus === Null, "n/a",
  RCTry[FreeQ[RCCanon[rPlus - R0], \[Phi]h | \[Psi]h | Bsh | Esh]]]];
RCSay["F2_III_B_TRANSVERSALITY_PROBE_NOW", transProbe[]];
(* (c) the extra connection removed as well: is anything left of the difference? *)
RCSayPlain["F2_III_C_UNDEFCOVD_CDM_RAN", RCTry[(UndefCovD[CDM]; True)] === True];
rBare = RCNewMessages["III_bare", split[]];
RCSayPlain["F2_III_C_VS_REFERENCE", If[rBare === $Failed || rBare === Null, "n/a", RCTry[RCVerdict[rBare, R0]]]];
RCSay["F2_III_REFERENCE_SPLIT", R0];
RCSay["F2_III_A_SPLIT_UNDER_MINUS1_EXPRESSION", rMinus];
RCSigns["at_end"];
RCDone["f2"];
