(* ::Package:: *)
(* cspell:words WXF inds normu xPand xPert xTensor xMAG PSALTer SetSlicing RCSigns RCStackOn wrsym RCSetSigns RCNewMessages *)
(* ==============================================================================
   RCSetupCore.wl -- package-independent harness for the R-C Wolfram probes (#567)
   ------------------------------------------------------------------------------
   PURPOSE  the verdict / digest / export / message / sign helpers, loading NO
            package, so a step can load xMAG (or PSALTer, or nothing) first and
            still use the harness. RCSetup.wl = this file + xPand + geometry.
   USAGE    Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]]
   VERDICTS RCVerdict[expr, target] -> "identical (SameQ)" | "proved-equal ..." |
            "proved-different" | "could-not-decide ...". Digests (SHA256 of the
            canonical InputForm) are fingerprints only, never the verdict.
   SIGNS    RCSigns[tag] prints xTensor's five sign globals. They are session-wide
            and read at call time, so a mixed session is a silent wrong answer:
            print them at start, after every Needs (xMAG sets $RiemannSign = -1,
            xMAG.m:110) and after any StartInducedDecomposition (it sets
            $ExtrinsicKSign and $AccelerationSign to -1, xMAG.m:1792).
   MESSAGES $MessageList is Protected: "$MessageList = {}" is a no-op that only
            prints Set::wrsym, which is why every message list printed by the
            first R-C runs was cumulative. Use RCNewMessages instead.
   ============================================================================== *)

(* GUARD. This file names xAct functions (ToCanonical, ContractMetric, NoScalar,
   ScreenDollarIndices) by their short names, so it must be read AFTER an xAct package is
   loaded -- otherwise the parser creates Global` symbols that shadow xAct's and every
   verdict is computed by inert functions that return their input unevaluated. That is what
   happened in tier1_xmag.wls runs 095608Z-100425Z: Part A was right, Parts B and C were
   silently meaningless. Load order: Needs[<the package this step is about>] first, then
   this file. *)
If[!MemberQ[$Packages, "xAct`xTensor`"],
  Print["RC_SETUP_ERROR=RCSetupCore.wl was read before any xAct package (load Needs[\"xAct`...\"] first)"];
  Quit[1]];

RCArgs = Association @@ (Rule @@ StringSplit[#, "=", 2] & /@
   Select[Rest[$ScriptCommandLine], StringContainsQ[#, "="] &]);
RCOutDir = Lookup[RCArgs, "outdir", Directory[]];
If[!DirectoryQ[RCOutDir], CreateDirectory[RCOutDir]];

RCSay[key_String, val_] := Print["RC_", key, "=", ToString[val, InputForm]];
RCSayPlain[key_String, val_] := Print["RC_", key, "=", val];

(* The sign globals live in xAct`xTensor`. Referring to them by their short names
   BEFORE xTensor is loaded creates Global` copies that then shadow xAct's -- which is
   what broke the first run of this script (t1/20260917T095608Z). Always go through the
   full context name, and say so when xTensor is not loaded yet. *)
$RCSignNames = {"$RiemannSign", "$RicciSign", "$TorsionSign", "$ExtrinsicKSign", "$AccelerationSign"};
RCSignValues[] := If[!MemberQ[$Packages, "xAct`xTensor`"], "xTensor not loaded",
  ToExpression["xAct`xTensor`" <> #] & /@ $RCSignNames];
RCSigns[tag_String] := RCSayPlain["SIGNS_" <> ToUpperCase[tag], ToString[RCSignValues[], InputForm]];
(* RCSetSigns[{sR, sr, sT}] sets the three curvature/torsion globals in xTensor's own
   context (never through a short name) and prints what it did. *)
RCSetSigns[{sR_, sr_, sT_}] := (
  ToExpression["xAct`xTensor`$RiemannSign = " <> ToString[sR]];
  ToExpression["xAct`xTensor`$RicciSign = " <> ToString[sr]];
  ToExpression["xAct`xTensor`$TorsionSign = " <> ToString[sT]];
  RCSigns["after_set"]);

(* RCNewMessages[expr]: evaluate expr and print only the messages IT generated
   ($MessageList is Protected, so it cannot be reset; take the tail instead). *)
(* HoldRest, not HoldAllComplete: with the latter the tag argument is held too, so
   tag_String never matches a computed name (first run of tier1_xmag.wls). *)
SetAttributes[RCNewMessages, HoldRest];
RCNewMessages[tag_String, expr_] := Module[{n = Length[$MessageList], res},
  res = expr;
  RCSay["MESSAGES_" <> ToUpperCase[tag], Drop[$MessageList, n]];
  res];

(* RCStackOn[msgs]: while active, print the xAct/Global stack whenever one of the
   named messages fires -- how a throw site is located without a debugger. *)
RCStackOn[msgs_List] := (
  $RCStackMsgs = msgs;
  Internal`AddHandler["Message", RCStackHandler]);
RCStackOff[] := Internal`RemoveHandler["Message", RCStackHandler];
RCStackHandler[Hold[Message[name_, ___], True]] := If[MemberQ[$RCStackMsgs, HoldForm[name]] || $RCStackMsgs === All,
  RCSay["STACK_" <> StringReplace[ToString[HoldForm[name]], {"::" -> "_", "$" -> ""}],
    Take[Select[Stack[_], !FreeQ[#, Symbol] &], UpTo[12]]]];
RCStackHandler[___] := Null;

RCVersions[] := (
  RCSay["WOLFRAM_VERSION", $Version];
  RCSay["XTENSOR_VERSION", xAct`xTensor`$Version];
  RCSigns["at_load"];
);

RCUsage[context_String, name_String] := StringTake[ToString[ToExpression[context <> name <> "::usage"]], UpTo[400]];

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
       {First[hits], "proved-equal up to c=" <> ToString[First[hits], InputForm]}]];


(* RCTry: run one call with its messages VISIBLE, turning a Throw or a message into $Failed
   (Check alone does not catch xPand's Throw@Message idiom, xPand.m:1789). *)
SetAttributes[RCTry, HoldAll];
(* A message is NOT a failure: xAct prints benign ones (DefMetric::old appears in the
   author's own successful run of StartInducedDecomposition). Using Check here made a
   working call look broken in the t1 runs of 2026-09-17. Only an uncaught Throw is a failure --
   xPand and xMAG both use the Throw@Message idiom (xPand.m:1789, xTensor.m:8302), which
   throws Null, so a Null result from Catch is reported as $Failed. Messages stay visible
   and are attributed per call by RCNewMessages. *)
RCTry[expr_] := Module[{r},
  r = Catch[Quiet[expr, {ToCanonical::noident, General::stop}]];
  If[r === Null, $Failed, r]];

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
