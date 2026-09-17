#!/usr/bin/env wolframscript
(* cspell:words PSALTer RCSigns epsilonSign DefGeometry COVDS xtensor RCSignValues *)
(* f6_psalter_signs.wls -- S3, second half: PSALTer's kernel. PSALTer is flat (it declares
   no curvature) and its sources set none of xTensor's sign globals; it only READS
   $epsilonSign when it gives epsilonG its components (Sources/ReloadPackage/DefGeometry.m).
   This confirms that nothing on the spectrum branch constrains the curvature convention.
   Own kernel: PSALTer starts subkernels and defines its own M4/G/CD. *)
Needs["xAct`xTensor`"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "RCSetupCore.wl"}]];
RCSigns["plain_xtensor"];
Needs["xAct`PSALTer`"];
RCSigns["after_psalter"];
RCSay["F6_PSALTER_PACKAGES", Select[$Packages, StringStartsQ[#, "xAct`"] &]];
RCSay["F6_PSALTER_VERSION_SYMBOLS", Select[Names["xAct`PSALTer`*Version*"], True &]];
RCSay["F6_SIGNATURE_USAGE", RCTry[StringTake[ToString[xAct`PSALTer`G::usage], UpTo[200]]]];
RCSay["F6_PSALTER_METRICS", RCTry[$Metrics]];
RCSay["F6_PSALTER_COVDS", RCTry[$CovDs]];
RCSayPlain["F6_PSALTER_RIEMANN_IS_VANISHING", RCTry[ToString[
  VanishingQ /@ Select[{RiemannCD, RicciCD, RicciScalarCD}, xTensorQ]]]];
RCSayPlain["F6_PSALTER_SETS_NO_SIGN_GLOBAL", RCSignValues[] === {1, 1, 1, 1, 1}];
RCDone["f6"];
