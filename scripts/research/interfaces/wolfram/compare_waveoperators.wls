#!/usr/bin/env wolframscript
(* ==============================================================================
   R-1 prototype: are the wave operators from different Stage-1 routes identical?

   Usage:
     QT_QPA_PLATFORM=offscreen wolframscript -file compare_waveoperators.wls LABEL=FILE [LABEL=FILE ...]
   FILE is either PSALTer's ParticleSpectrograph<Name>.mx (the association is Global`<Name>)
   or our spectrum.wxf (key "wave_operator").

   Verdict = in-kernel SameQ against the FIRST file (the reference), three-valued when it
   fails: proved-equal (Simplify of every element difference is 0), proved-different, or
   could-not-decide. Digests are printed as fingerprints only: a checksum is never the verdict.
   A digest is computed twice per operator to show BinarySerialize is stable in one kernel.
   ============================================================================== *)

Needs["xAct`PSALTer`"];

specs = Rest[$ScriptCommandLine];
If[Length[specs] < 2, Print["COMPARE_ERROR=need at least two LABEL=FILE arguments"]; Exit[2]];

load[file_String] := Which[
  StringEndsQ[file, ".mx"],
    With[{name = StringReplace[FileBaseName[file], "ParticleSpectrograph" -> ""]},
      Get[file]; Symbol["Global`" <> name][xAct`PSALTer`WaveOperator]],
  StringEndsQ[file, ".wxf"], Import[file, "WXF"]["wave_operator"],
  True, $Failed];

digest[x_] := Hash[BinarySerialize[x], "SHA256", "HexString"];

ops = Association@Map[
  With[{parts = StringSplit[#, "=", 2]}, parts[[1]] -> load[parts[[2]]]] &, specs];

KeyValueMap[
  (Print["COMPARE_DIGEST ", #1, "=", digest[#2]];
   Print["COMPARE_DIGEST_STABLE ", #1, "=", digest[#2] === digest[#2]];
   Print["COMPARE_SECTORS ", #1, "=", If[ListQ[#2], Length[#2], "NOT-A-LIST"]]) &, ops];

refLabel = First@Keys[ops]; ref = ops[refLabel];
verdicts = KeyValueMap[
  Function[{label, op},
    Module[{v},
      v = Which[
        SameQ[op, ref], "identical (SameQ)",
        Dimensions[op] =!= Dimensions[ref], "proved-different (dimensions " <> ToString[Dimensions[op]] <> " vs " <> ToString[Dimensions[ref]] <> ")",
        True, Module[{diffs = TimeConstrained[Simplify /@ Flatten[op - ref], 120, $Aborted]},
          Which[diffs === $Aborted, "could-not-decide (Simplify timed out)",
            AllTrue[diffs, PossibleZeroQ], "proved-equal (not SameQ; every element difference simplifies to 0)",
            True, "proved-different"]]];
      Print["COMPARE_VERDICT ", label, " vs ", refLabel, ": ", v];
      v]],
  KeyDrop[ops, refLabel]];

If[AllTrue[verdicts, StringStartsQ[#, "identical"] || StringStartsQ[#, "proved-equal"] &],
  Print["COMPARE PASSED"]; Exit[0],
  Print["COMPARE FAILED"]; Exit[1]]
