#!/usr/bin/env wolframscript
(* cspell:words xPand xMAG xPerm xTensor xBrauer *)
(* f11_shadowed_names.wls -- which bare names stop meaning the System symbol once a package loads.

   WHY. Every earlier run of this lane printed RC_WOLFRAM_VERSION={"0.4.4", {2025, 4, 1}}: the
   harness wrote a bare $Version, and RCSetup.wl reads it after xPand is loaded, so the name was
   xPand's own $Version, not the kernel's. No message; a plausible value under a wrong label. The
   harness now writes System`$Version, and wl_lint.py refuses a bare $Version (defect SHADOW).
   This step is the evidence for both: it lists, after xPand and again after xMAG, every short
   name an xAct context shares with System`, and the context that the bare name resolves to.
   Expected: $Version resolves to the last-loaded package (xPand, then xMAG); the others are
   xAct's deliberate replacements (Symmetrize, Commutator, Up, Down, ...), recorded, not traps.
   Usage: run_lane.sh f11 *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetup.wl"}]];
RCLoad[];
RCSayPlain["F11_HARNESS_VERSION_IS_THE_KERNELS", StringQ[System`$Version] && StringStartsQ[System`$Version, ToString[System`$VersionNumber]]];
short[s_String] := Last[StringSplit[s, "`"]];
systemNames = short /@ Names["System`*"];
xActContexts[] := Select[Contexts[],
   StringMatchQ[#, "xAct`" ~~ Except["`"] .. ~~ "`"] || MemberQ[{"BrauerAlgebra`", "SymmetricFunctions`"}, #] &];
sharedWithSystem[] := Union @@ (Intersection[short /@ Names[# <> "*"], systemNames] & /@ xActContexts[]);
report[tag_String] := (
  RCSay["F11_SHARED_WITH_SYSTEM_" <> tag, {#, Context[#]} & /@ sharedWithSystem[]];
  RCSayPlain["F11_BARE_VERSION_CONTEXT_" <> tag, Context["$Version"]]);

report["AFTER_XPAND"];
RCNewMessages["needs_xmag", Needs["xAct`xMAG`"]];
report["AFTER_XMAG"];
RCSayPlain["F11_BARE_VERSION_IS_NEVER_THE_KERNELS",
  Context["$Version"] =!= "System`"];
RCSigns["f11_end"];
RCDone["f11"];
