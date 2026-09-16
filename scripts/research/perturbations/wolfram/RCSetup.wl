(* ::Package:: *)
(* cspell:words DefMatterFields DefMetricFields SetSlicing ToxPand WXF inds normu xPand xPert xTensor *)
(* ==============================================================================
   RCSetup.wl -- shared definitions for the R-C Wolfram probes (#567)
   ------------------------------------------------------------------------------
   PURPOSE  load xPand on the certified bundle FIRST (so every later symbol in
            this file and in the calling script resolves into the xAct contexts),
            set up an FRW slicing with a selectable normal norm, and provide the
            verdict / digest / export helpers every probe uses. Research code.
   USAGE    Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]]
            from a .wls run as  wolframscript -file X.wls outdir=DIR key=value ...
   NAMES    xPand exports a, H, P, T, k, Es, Et, Ev, Bs, Bv, Ls, Lv, Vs, Vv, V0,
            Xi, CS, Connection, nt, av -- never reuse them. Indices are ia..is.
   VERDICTS RCVerdict[expr, target] -> "identical (SameQ)" | "proved-equal ..." |
            "proved-different" | "could-not-decide ...". Digests (SHA256 of the
            canonical InputForm) are fingerprints only, never the verdict.
   ============================================================================== *)

$RCLoadT0 = AbsoluteTime[];
$RCLoadFailed = False;
Check[Needs["xAct`xPand`"], $RCLoadFailed = True];
$RCLoadMessages = $MessageList;
$RCLoadWall = AbsoluteTime[] - $RCLoadT0;

RCArgs = Association @@ (Rule @@ StringSplit[#, "=", 2] & /@
   Select[Rest[$ScriptCommandLine], StringContainsQ[#, "="] &]);
RCOutDir = Lookup[RCArgs, "outdir", Directory[]];
If[!DirectoryQ[RCOutDir], CreateDirectory[RCOutDir]];

RCSay[key_String, val_] := Print["RC_", key, "=", ToString[val, InputForm]];
RCSayPlain[key_String, val_] := Print["RC_", key, "=", val];

RCLoad[] := (
  RCSayPlain["LOAD_FAILED", $RCLoadFailed];
  RCSayPlain["LOAD_WALL_S", Round[$RCLoadWall, 0.1]];
  RCSay["LOAD_MESSAGES", $RCLoadMessages];
  RCSay["WOLFRAM_VERSION", $Version];
  RCSay["XTENSOR_VERSION", xAct`xTensor`$Version];
  RCSay["XPERT_VERSION", xAct`xPert`$Version];
  RCSay["XPAND_VERSION", xAct`xPand`$Version];
  RCSay["CONFORMAL_TIME", $ConformalTime];
  RCSay["XPAND_SYMBOL_COUNT", Length[Names["xAct`xPand`*"]]];
);

RCUsage[name_String] := StringTake[ToString[ToExpression["xAct`xPand`" <> name <> "::usage"]], UpTo[400]];

(* FRW slicing. normu = -1 selects mostly-plus (n.n = -1); +1 selects mostly-minus.
   DefMetric's first argument is the sign of the determinant, -1 for either. *)
RCGeometry[normu_: -1, type_String: "FLFlat"] := (
  DefManifold[M4, 4, {ia, ib, ic, id, ie, ig, ip, iq, ir, is}];
  DefMetric[-1, g[-ia, -ib], CD, {";", "\[Del]"}];
  DefMetricPerturbation[g, dg, eps];
  SetSlicing[g, n, normu, h, cd, {"|", "D"}, type];
  DefMetricFields[g, dg, h, eps];
  RCEps = $PerturbationParameter;
  RCSay["PERTURBATION_PARAMETER", RCEps];
  RCSayPlain["GEOMETRY", "normu=" <> ToString[normu] <> " type=" <> type];
);

RCCanon[x_] := NoScalar[ToCanonical[ContractMetric[NoScalar[x]]]];
RCDigest[x_] := Hash[ToString[ScreenDollarIndices[RCCanon[x]], InputForm], "SHA256", "HexString"];

RCVerdict[expr_, target_, budget_: 120] := Module[{d, s},
  If[expr === target, Return["identical (SameQ)"]];
  d = TimeConstrained[RCCanon[expr - target], budget, $Aborted];
  Which[
    d === $Aborted, "could-not-decide (ToCanonical timed out)",
    d === 0, "proved-equal (ToCanonical of the difference is 0)",
    True, s = TimeConstrained[Simplify[d], budget, $Aborted];
      Which[s === $Aborted, "could-not-decide (Simplify timed out)",
            s === 0, "proved-equal (Simplify of the canonical difference is 0)",
            True, "proved-different"]]];

(* Find a numeric normalization c with expr == c target; returns {c, verdict}. *)
RCNormalize[expr_, target_, cands_List: {1, -1, 2, -2, 1/2, -1/2, 4, -4, 1/4, -1/4}] :=
  Module[{hits},
    hits = Select[cands, StringStartsQ[RCVerdict[expr, # target, 60], "identical" | "proved-equal"] &];
    If[hits === {}, {None, "proved-different (no candidate normalization)"},
       {First[hits], "proved-equal up to c=" <> ToString[First[hits]]}]];

RCReport[name_String, expr_, target_, control_: None] := Module[{v, vc = "n/a", assoc},
  v = RCVerdict[expr, target];
  If[control =!= None, vc = RCVerdict[expr, control]];
  RCSayPlain[ToUpperCase[name] <> "_VERDICT", v];
  RCSayPlain[ToUpperCase[name] <> "_CONTROL_VERDICT", vc];
  RCSayPlain[ToUpperCase[name] <> "_DIGEST", RCDigest[expr]];
  RCSayPlain[ToUpperCase[name] <> "_DIGEST_STABLE", RCDigest[expr] === RCDigest[expr]];
  assoc = <|"label" -> name, "expr" -> ToString[RCCanon[expr], InputForm],
            "target" -> ToString[target, InputForm], "verdict" -> v, "controlVerdict" -> vc,
            "digest" -> RCDigest[expr], "utc" -> DateString["ISODateTime"]|>;
  Export[FileNameJoin[{RCOutDir, name <> ".wxf"}], assoc, "WXF"];
  Export[FileNameJoin[{RCOutDir, name <> ".txt"}], ToString[RCCanon[expr], InputForm], "Text"];
  v];

RCDone[step_String] := Print["RC_", ToUpperCase[step], "_DONE"];
