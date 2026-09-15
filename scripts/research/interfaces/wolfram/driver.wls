#!/usr/bin/env wolframscript
(* ==============================================================================
   R-1 prototype: the FIXED Stage-1 driver (never generated, identical for every theory)

   Usage (launched by proto/derive_proto.py, never typed by a user):
     QT_QPA_PLATFORM=offscreen wolframscript -file driver.wls <theory.wxf> <outDir> [--validate-only]

   Exit codes: 0 = spectrum.wxf written ("STAGE1 PASSED"), 1 = the theory was refused or the
   derivation failed ("STAGE1_ERROR=<tag>: <message>"), 2 = usage error.
   The caller must still verify the artifact: Wolfram 14.3 can segfault on shutdown after
   correct output, and PSALTer contains a bare Quit[] (ConjectureInverse.m:37-39) that
   exits 0 without a message. Exit status alone is never the verdict.

   Error signalling lives HERE, once: the whole body is inside Catch, so an uncaught
   Throw from PSALTer can no longer exit 0 (#561). Because this file is fixed, adding that
   no longer changes "every generated script".
   ============================================================================== *)

args = Rest[$ScriptCommandLine];
validateOnly = MemberQ[args, "--validate-only"];
args = DeleteCases[args, "--validate-only"];
If[Length[args] =!= 2,
  Print["STAGE1_ERROR=Usage: driver.wls <theory.wxf> <outDir> [--validate-only]"]; Exit[2]];
{theoryFile, outDir} = args;
packageFile = FileNameJoin[{DirectoryName[$InputFileName], "Stage1Proto.wl"}];

(* PSALTer freezes its output directory when it loads (PSALTer.m:49-56): set it first. *)
SetDirectory[outDir];

result = Catch[
  Needs["xAct`PSALTer`"];
  Get[packageFile];
  R1Proto`RunStage1[Import[theoryFile, "WXF"], outDir, "ValidateOnly" -> validateOnly]
];

Which[
  AssociationQ[result],
    Print["STAGE1 PASSED"]; Exit[0],
  FailureQ[result],
    Print["STAGE1_ERROR=", result[[1]], ": ", result["Message"]]; Exit[1],
  True,
    Print["STAGE1_ERROR=Unexpected: ", ToString[Short[result, 3], InputForm]]; Exit[1]
]
