#!/usr/bin/env wolframscript
(* cspell:words xPand xMAG xBrauer TraceFree MetricOfTensor FirstMetricQ FrozenMetricQ MasterOf SeparateMetric SeparateMetric1 SeparateMetric2 ContractMetric1 IndicesDown IndicesUp NormalVectorOf epsilonh AutomaticRules DownValues SubValues UpValues normu ScalarsOfProduct EinsteinToRicci FirstMetricOfVBundle CONTRACTMETRIC SEPARATEMETRIC NORMALVECTOROF EINSTEINTORICCI INDICESDOWN TOCANONICAL noWs pkg prod invmetric UPVALUE UPVALUES Unsets EPSH *)
(* f8_contractmetric.wls -- R-C closing pass, step f8: WHICH package stops xPand simplifying
   the transverse and traceless sectors, and by which mechanism.

   R-C reported the loss under the heading "loading xMAG disables other packages' automatic
   simplifications" and suspected xMAG's own replacement of xTensor's curvature relations.
   Post-merge source reading moved the blame to xBrauer, which xMAG pulls in through its
   BeginPackage, and the memo carries that as an interim note. This script replaces the
   reading with a measurement, in ONE kernel and THREE states, so the attribution is computed
   rather than asserted:
     base      xPand + FRW slicing only
     xbrauer   Needs["xAct`xBrauer`"] alone -- the control that exonerates or convicts xMAG
     xmag      Needs["xAct`xMAG`"] on top -- does xMAG add anything of its own?
   The mechanism under test (verified line by line in the installed sources):
     xBrauer.m:1938-1940 Unsets three of xTensor's ContractMetric1 rules (xTensor.m:8584-8586
       second-slot pairing and :8590-8592 first-slot pairing are two triples with IDENTICAL
       left-hand sides, so how many the three Unsets remove is NOT decidable by reading --
       F8_SUB_*_CONTRACTMETRIC1 and the slot-order probes settle it);
     :1947-1974 re-installs twelve with the extra guard
       MetricOfTensor[tensor] === metric || (MetricOfTensor[tensor] === Null && FirstMetricQ[metric])
       for the tensor slot, and the same shape with MasterOf[covd] for the derivative slot;
     xPand appends the induced metric h after g (xPand.m:1754), so FirstMetricQ[h] is False,
       and xPand's projected fields carry no Master, so the guard is False for h and True for
       g -- the g-versus-h asymmetry that makes this diagnosable.
   Three rival mechanisms inside xBrauer are excluded explicitly rather than assumed away:
     (a) the public driver, :1906, which changes $NonFrozenMetrics to $Metrics -- excluded by
         calling ContractMetric[x, h] with the metric named;
     (b) xAct`xTensor`Private`NormalVectorOf, replaced at :1899-1901;
     (c) rule ORDERING rather than the guard -- excluded by checking that the residual still
         contains an uncontracted h, not merely that the answer is nonzero.
   Also measured, because they are independent regressions and one is an improvement:
     xBrauer.m:1866 DELETES xTensor's SeparateMetric1 catch-all (xTensor.m:8682, which
       returned expr unchanged for a non-first metric) and dispatches induced and frozen
       metrics to SeparateMetric2 instead; :1874 replaces SeparateMetric2 itself, whose covd
       branch has NO MasterOf === Null fallback, which is what reaches xPand's IndicesDown
       (xPand.m:2412).
     xBrauer.m:1958 carries a one-character defect: its guard is swallowed into the argument
       of ScalarsOfProduct, reading ScalarsOfProduct[prod][metric[b,c] && (guard)] where its
       five siblings read ScalarsOfProduct[prod][metric[b,c]] && (guard).
   Every private symbol is reached through an existence test first, and a MISSING value is a
   legitimate result: MetricOfTensor is xBrauer's own public function and does NOT EXIST
   before xBrauer loads (zero occurrences in xTensor.m), and there are three private
   FirstMetricQ copies (xTensor, xBrauer, xMAG).
   Usage: run_lane.sh f8 . No arguments; no DefChart; no splits, so no ToxPand cost. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCGeometry[-1, "FLFlat"];

(* ---------------- guarded access to private and package-dependent symbols ---------------- *)
RCHas[full_String] := Names[full] =!= {};
RCNum[full_String, kind_String] := If[RCHas[full],
  RCTry[ToExpression["Length[" <> kind <> "[" <> full <> "]]"]], "MISSING"];
RCHashOf[full_String, kind_String] := If[RCHas[full],
  RCTry[StringTake[Hash[ToString[ToExpression[kind <> "[" <> full <> "]"], InputForm],
                        "SHA256", "HexString"], 12]], "MISSING"];
RCEvalIf[nameStr_String, exprStr_String] := If[RCHas[nameStr],
  RCTry[ToExpression[exprStr]], "MISSING:" <> nameStr];
noWs[s_String] := StringDelete[s, Whitespace];
subTextCount[full_String, pat_String] := If[!RCHas[full], "MISSING",
  RCTry[Count[ToExpression["SubValues[" <> full <> "]"],
    r_ /; StringContainsQ[noWs[ToString[r, InputForm]], pat]]]];
(* xBrauer.m:1958 carries a one-character defect: its guard is swallowed into the argument of
   ScalarsOfProduct, reading ScalarsOfProduct[prod][metric[b,c] && (guard)] where its five
   siblings read ScalarsOfProduct[prod][metric[b,c]] && (guard). A stored SubValue does not
   show this (the rendered rule loses the distinction), so the check is done on the package
   SOURCE, located through FindFile rather than any absolute path: one character after
   "metric[b,c]" decides, "]" for the sibling shape and "&" for the swallowed guard. *)
sopSourceShapes[] := Module[{f, lines, hits},
  (* FindFile resolves the context to Kernel/init.m, which only Gets the real file; step up *)
  f = RCTry[FindFile["xAct`xBrauer`"]];
  If[!StringQ[f], Return["FindFile-failed"]];
  If[FileBaseName[f] === "init",
    f = FileNameJoin[{ParentDirectory[DirectoryName[f]], "xBrauer.m"}]];
  If[!FileExistsQ[f], Return["source-not-found"]];
  lines = RCTry[Import[f, "Lines"]];
  If[!ListQ[lines], Return["Import-failed"]];
  hits = Select[lines, StringContainsQ[#, "ScalarsOfProduct"] &];
  {Length[hits], Flatten[StringCases[StringDelete[#, Whitespace],
     "ScalarsOfProduct[prod][metric[b,c]" ~~ x_ :> x] & /@ hits]}];

watchSub = {"xAct`xTensor`Private`ContractMetric1", "xAct`xTensor`Private`SeparateMetric1",
            "xAct`xTensor`Private`SeparateMetric2"};
watchDown = {"xAct`xTensor`Private`PRJ", "xAct`xTensor`Private`NormalVectorOf",
             "xAct`xTensor`ContractMetric", "xAct`xTensor`ToCanonical",
             "xAct`xTensor`SeparateMetric", "xAct`xTensor`EinsteinToRicci",
             "xAct`xPand`IndicesDown"};
shortName[s_String] := ToUpperCase[Last@StringSplit[s, "`"]];
snap[tag_] := (
  Do[RCSayPlain["F8_SUB_" <> ToUpperCase[tag] <> "_" <> shortName[s],
       ToString[{RCNum[s, "SubValues"], RCHashOf[s, "SubValues"]}, InputForm]], {s, watchSub}];
  Do[RCSayPlain["F8_DOWN_" <> ToUpperCase[tag] <> "_" <> shortName[s],
       ToString[{RCNum[s, "DownValues"], RCHashOf[s, "DownValues"]}, InputForm]], {s, watchDown}];
  RCSayPlain["F8_UP_" <> ToUpperCase[tag],
    ToString[{Length[UpValues[Eth]], Length[UpValues[Evh]], Length[UpValues[Esh]],
              Length[UpValues[epsilonh]], Length[UpValues[h]]}, InputForm]];
  (* the tuple is {rules mentioning MetricOfTensor, MasterOf, FirstMetricQ, ScalarsOfProduct} *)
  RCSayPlain["F8_GUARDCOUNT_" <> ToUpperCase[tag],
    ToString[{subTextCount["xAct`xTensor`Private`ContractMetric1", "MetricOfTensor"],
              subTextCount["xAct`xTensor`Private`ContractMetric1", "MasterOf"],
              subTextCount["xAct`xTensor`Private`ContractMetric1", "FirstMetricQ"],
              subTextCount["xAct`xTensor`Private`ContractMetric1", "ScalarsOfProduct"]},
             InputForm]];
  RCSigns[tag];
);

(* ---------------- the probe set ---------------- *)
(* S1 pairs the metric's FIRST slot  (the xTensor.m:8590 / xBrauer.m:1954 triple)
   S2 pairs the metric's SECOND slot (the xTensor.m:8584 / xBrauer.m:1947 triple)
   G1 is S1 with the FIRST metric g, where FirstMetricQ is True
   E1 contracts h into epsilonh, whose MasterOf IS h: same metric, same slot, Master present
   D  is f2's own expression, which exercises the tensor slot and the derivative slot at once
   U  is already contracted and carries no metric at all: the transversality UpValue itself *)
probeList[] := {
  {"S1", h[ia, ib] Evh[LI[1], LI[0], -ia]},
  {"S2", h[ib, ia] Evh[LI[1], LI[0], -ia]},
  {"G1", g[ia, ib] Evh[LI[1], LI[0], -ia]},
  {"E1", h[ia, ib] epsilonh[-ia, -ic, -id]},
  {"D", h[ia, ib] cd[-ib]@Evh[LI[1], LI[0], -ia]},
  {"U", cd[-ia]@Evh[LI[1], LI[0], ia]},
  {"UT", cd[-ia]@Eth[LI[1], LI[0], ia, -ib]},
  {"UTR", Eth[LI[1], LI[0], ia, -ia]}};
pipe[1][x_] := RCTry[ContractMetric[x]];
pipe[2][x_] := RCTry[ContractMetric[x, h]];
pipe[3][x_] := RCTry[ToCanonical[x]];
pipe[4][x_] := RCTry[ToCanonical[ContractMetric[x]]];
state[tag_, nm_, k_] := Module[{r},
  r = pipe[k][Last@SelectFirst[probeList[], First[#] === nm &]];
  RCSayPlain["F8_" <> ToUpperCase[tag] <> "_" <> nm <> "_P" <> ToString[k],
    ToString[{If[r === $Failed, "FAILED", If[FreeQ[r, h], "no-h", "HAS-H"]],
              If[r === $Failed, "FAILED", If[FreeQ[r, g], "no-g", "has-g"]],
              r === 0}, InputForm]];
  r];
probeAll[tag_] := Do[state[tag, nm, k], {nm, First /@ probeList[]}, {k, {1, 2, 3, 4}}];
showOne[tag_, nm_] := RCSay["F8_" <> ToUpperCase[tag] <> "_" <> nm <> "_P1_EXPR", pipe[1][
  Last@SelectFirst[probeList[], First[#] === nm &]]];

(* the second-regression probes *)
sepProbes[tag_] := (
  RCSay["F8_SEP_" <> ToUpperCase[tag] <> "_G_ON_TENSOR", RCTry[SeparateMetric[g][Evh[LI[1], LI[0], ia], ia]]];
  RCSay["F8_SEP_" <> ToUpperCase[tag] <> "_G_ON_CD", RCTry[SeparateMetric[g][cd[ia]@Esh[LI[1], LI[0]], ia]]];
  RCSay["F8_SEP_" <> ToUpperCase[tag] <> "_G_ON_EPSH", RCTry[SeparateMetric[g][epsilonh[ia, -ib, -ic], ia]]];
  RCSay["F8_SEP_" <> ToUpperCase[tag] <> "_H_ON_TENSOR", RCTry[SeparateMetric[h][Evh[LI[1], LI[0], ia], ia]]];
  RCSay["F8_IDOWN_" <> ToUpperCase[tag] <> "_ON_CD", RCTry[IndicesDown[cd[ia]@Esh[LI[1], LI[0]]]]];
  RCSay["F8_IDOWN_" <> ToUpperCase[tag] <> "_ON_EPSH", RCTry[IndicesDown[epsilonh[ia, -ib, -ic]]]];
  RCSay["F8_IDOWN_" <> ToUpperCase[tag] <> "_ON_TENSOR", RCTry[IndicesDown[Evh[LI[1], LI[0], ia]]]];
);
guardProbes[tag_] := (
  RCSay["F8_GUARD_" <> ToUpperCase[tag] <> "_METRICOFTENSOR",
    {RCEvalIf["xAct`xBrauer`MetricOfTensor", "xAct`xBrauer`MetricOfTensor[Evh]"],
     RCEvalIf["xAct`xBrauer`MetricOfTensor", "xAct`xBrauer`MetricOfTensor[Eth]"],
     RCEvalIf["xAct`xBrauer`MetricOfTensor", "xAct`xBrauer`MetricOfTensor[epsilonh]"],
     RCEvalIf["xAct`xBrauer`MetricOfTensor", "xAct`xBrauer`MetricOfTensor[h]"],
     RCEvalIf["xAct`xBrauer`MetricOfTensor", "xAct`xBrauer`MetricOfTensor[g]"]}];
  RCSay["F8_GUARD_" <> ToUpperCase[tag] <> "_FIRSTMETRICQ",
    {RCEvalIf["xAct`xBrauer`Private`FirstMetricQ", "xAct`xBrauer`Private`FirstMetricQ[h]"],
     RCEvalIf["xAct`xBrauer`Private`FirstMetricQ", "xAct`xBrauer`Private`FirstMetricQ[g]"],
     RCEvalIf["xAct`xTensor`Private`FirstMetricQ", "xAct`xTensor`Private`FirstMetricQ[h]"],
     RCEvalIf["xAct`xTensor`Private`FirstMetricQ", "xAct`xTensor`Private`FirstMetricQ[g]"],
     RCEvalIf["xAct`xMAG`Private`FirstMetricQ", "xAct`xMAG`Private`FirstMetricQ[h]"]}];
  RCSay["F8_GUARD_" <> ToUpperCase[tag] <> "_MASTEROF",
    RCTry[{MasterOf[cd], MasterOf[CD], MasterOf[Evh], MasterOf[Eth], MasterOf[epsilonh]}]];
  RCSay["F8_GUARD_" <> ToUpperCase[tag] <> "_FIRSTMETRICOFVBUNDLE", RCTry[FirstMetricOfVBundle[Tangent[M4], True]]];
  RCSay["F8_GUARD_" <> ToUpperCase[tag] <> "_METRICS_ORDER", $Metrics];
  RCSayPlain["F8_GUARD_" <> ToUpperCase[tag] <> "_VALUE_FOR_H_AND_EVH",
    ToString[RCEvalIf["xAct`xBrauer`MetricOfTensor",
      "xAct`xBrauer`MetricOfTensor[Evh] === h || (xAct`xBrauer`MetricOfTensor[Evh] === Null && xAct`xBrauer`Private`FirstMetricQ[h])"], InputForm]];
  RCSayPlain["F8_GUARD_" <> ToUpperCase[tag] <> "_VALUE_FOR_H_AND_EPSH",
    ToString[RCEvalIf["xAct`xBrauer`MetricOfTensor",
      "xAct`xBrauer`MetricOfTensor[epsilonh] === h || (xAct`xBrauer`MetricOfTensor[epsilonh] === Null && xAct`xBrauer`Private`FirstMetricQ[h])"], InputForm]];
);

(* ================= C1: baseline, xPand and the slicing only ================= *)
RCSay["F8_CONTEXTS_OF_FIELDS", {Context[Eth], Context[Evh], Context[epsilonh], Context[h]}];
RCSay["F8_CONTRACTMETRIC_NAMES", Names["*`ContractMetric"]];
RCSay["F8_TOCANONICAL_NAMES", Names["*`ToCanonical"]];
RCSayPlain["F8_AUTOMATICRULES_G_H_BASE", RCTry[RCCanon[g[ia, ib] h[-ib, -ic] - h[ia, -ic]]] === 0];
snap["base"]; guardProbes["base"]; probeAll["base"]; sepProbes["base"];
showOne["base", "S1"]; showOne["base", "D"];
baseS1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "S1" &]];
baseE1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "E1" &]];
baseG1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "G1" &]];
baseD4 = pipe[4][Last@SelectFirst[probeList[], First[#] === "D" &]];
baseU1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "U" &]];
baseUV = {Length[UpValues[Eth]], Length[UpValues[Evh]]};
baseSubN = RCNum["xAct`xTensor`Private`ContractMetric1", "SubValues"];
RCSayPlain["F8_BASE_TRANSVERSALITY_RECOVERS_F2", baseD4 === 0 && baseU1 === 0];
RCSayPlain["F8_BASE_TOCANONICAL_ALONE_CONTRACTS_H",
  RCTry[FreeQ[pipe[3][Last@SelectFirst[probeList[], First[#] === "S1" &]], h]]];

(* ================= C2: xBrauer alone ================= *)
RCNewMessages["needs_xbrauer", Needs["xAct`xBrauer`"]];
RCSayPlain["F8_XBRAUER_LOADED", MemberQ[$Packages, "xAct`xBrauer`"]];
RCSay["F8_XBRAUER_VERSION", RCEvalIf["xAct`xBrauer`$Version", "xAct`xBrauer`$Version"]];
snap["xbrauer"]; guardProbes["xbrauer"]; probeAll["xbrauer"]; sepProbes["xbrauer"];
showOne["xbrauer", "S1"]; showOne["xbrauer", "D"];
brS1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "S1" &]];
brS2 = pipe[1][Last@SelectFirst[probeList[], First[#] === "S2" &]];
brS1m = pipe[2][Last@SelectFirst[probeList[], First[#] === "S1" &]];
brE1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "E1" &]];
brG1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "G1" &]];
brD4 = pipe[4][Last@SelectFirst[probeList[], First[#] === "D" &]];
brU1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "U" &]];
brUV = {Length[UpValues[Eth]], Length[UpValues[Evh]]};
brSubN = RCNum["xAct`xTensor`Private`ContractMetric1", "SubValues"];
brSubHash = RCHashOf["xAct`xTensor`Private`ContractMetric1", "SubValues"];
RCSayPlain["F8_XBRAUER_SUBVALUE_DELTA", If[IntegerQ[brSubN] && IntegerQ[baseSubN], brSubN - baseSubN, "n/a"]];
RCSayPlain["F8_XBRAUER_SLOT1_CONTRACTS", RCTry[FreeQ[brS1, h]]];
RCSayPlain["F8_XBRAUER_SLOT2_CONTRACTS", RCTry[FreeQ[brS2, h]]];
RCSayPlain["F8_XBRAUER_AUTOMATICRULES_SURVIVE", RCTry[RCCanon[g[ia, ib] h[-ib, -ic] - h[ia, -ic]]] === 0];
RCSay["F8_XBRAUER_SOURCE_SCALARSOFPRODUCT_SHAPES", sopSourceShapes[]];
RCSayPlain["F8_XBRAUER_1958_PAREN_DEFECT_PRESENT",
  MemberQ[Last[sopSourceShapes[]], "&"]];

(* ================= C9: xMAG on top ================= *)
RCNewMessages["needs_xmag", Needs["xAct`xMAG`"]];
RCSayPlain["F8_XMAG_LOADED", MemberQ[$Packages, "xAct`xMAG`"]];
RCSigns["after_xmag_load"];
RCSetSigns[{1, 1, 1}];           (* hold the sign fixed so it cannot confound the comparison *)
snap["xmag"]; guardProbes["xmag"]; probeAll["xmag"]; sepProbes["xmag"];
showOne["xmag", "S1"]; showOne["xmag", "D"];
mgS1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "S1" &]];
mgE1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "E1" &]];
mgG1 = pipe[1][Last@SelectFirst[probeList[], First[#] === "G1" &]];
mgD4 = pipe[4][Last@SelectFirst[probeList[], First[#] === "D" &]];
mgSubN = RCNum["xAct`xTensor`Private`ContractMetric1", "SubValues"];
mgSubHash = RCHashOf["xAct`xTensor`Private`ContractMetric1", "SubValues"];
RCSayPlain["F8_XMAG_ADDS_NOTHING_TO_CONTRACTMETRIC1", mgSubN === brSubN && mgSubHash === brSubHash];
RCSayPlain["F8_CONTRACTMETRIC1_SUBVALUE_HASHES", ToString[{brSubHash, mgSubHash}, InputForm]];
RCSayPlain["F8_XMAG_PROBES_MATCH_XBRAUER",
  RCTry[{mgS1 === brS1, mgE1 === brE1, mgG1 === brG1, mgD4 === brD4}]];

(* ================= C10: the attribution, computed ================= *)
uvAlive = (baseUV === brUV) && (baseU1 === 0) && (brU1 === 0);
guardFalse = (RCEvalIf["xAct`xBrauer`MetricOfTensor",
  "xAct`xBrauer`MetricOfTensor[Evh] === h || (xAct`xBrauer`MetricOfTensor[Evh] === Null && xAct`xBrauer`Private`FirstMetricQ[h])"] === False);
gWorks = RCTry[FreeQ[brG1, g]] === True;
masterDisc = (RCTry[FreeQ[brE1, h]] === True) && (RCTry[FreeQ[brS1, h]] === False);
driverExcl = RCTry[FreeQ[brS1m, h]] === False;
orderExcl = RCTry[!FreeQ[brS1, h]] === True;
regressed = (baseD4 === 0) && (brD4 =!= 0);
RCSayPlain["F8_VERDICT_UPVALUES_ALIVE", uvAlive];
RCSayPlain["F8_VERDICT_REGRESSION_REPRODUCED_BY_XBRAUER_ALONE", regressed];
RCSayPlain["F8_VERDICT_GUARD_IS_FALSE_FOR_H", guardFalse];
RCSayPlain["F8_VERDICT_G_STILL_CONTRACTS", gWorks];
RCSayPlain["F8_VERDICT_MASTER_IS_THE_DISCRIMINATOR", masterDisc];
RCSayPlain["F8_VERDICT_DRIVER_EXCLUDED", driverExcl];
RCSayPlain["F8_VERDICT_RULE_ORDERING_EXCLUDED", orderExcl];
RCSayPlain["F8_ATTRIBUTION",
  If[And[uvAlive, regressed, guardFalse, gWorks, masterDisc, driverExcl, orderExcl],
    "xBrauer.m:1947-1974 guard (MetricOfTensor/MasterOf with FirstMetricQ) blocks h-contraction into Master-less projected tensors; xMAG only pulls xBrauer in",
    "NOT ESTABLISHED -- read the individual verdict sentinels above"]];
RCSigns["f8_end"];
RCDone["f8"];
