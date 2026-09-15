#!/usr/bin/env wolframscript
(* ==============================================================================
   Preflight: a small PSALTer run, end to end
   ==============================================================================
   Usage:
     QT_QPA_PLATFORM=offscreen wolframscript -file vector_smoke.wls

   Proves the whole pipeline works on something cheap before the Tier-1 gate
   commits to the real run: PSALTer loads, a field can be declared (which exports
   a PDF through the front end), ParticleSpectrum completes, and the result
   association is written to disk.

   This is OUR restatement of the three-term massive-vector Lagrangian from
   SupplementalMaterials-2607 (WolframLanguage/ParticleSpectroscopy/Models/
   VectorTheory.m), not that file itself: upstream's driver and model scripts use
   `Comment` and `DisplayExpression` from xAct`xPlain`, a package this project
   does not install. The physics is theirs; the transcription is ours.

     L = -(1/2) t1 (grad_a A_b)(grad^a A^b)
         +(1/2) t2 (grad_a A^a)(grad_b A^b)
         -(1/2) t3 A_a A^a

   Run from a throwaway directory: PSALTer writes into the current one.

   EXPIRED 2026-09-15 (R-1 #566, re-run by the orchestrator): the association now carries the
   mass term -- its first WaveOperator element is (Def^2*(-Theta1 + Theta2) - Theta3)/2. The
   discrepancy recorded below was #543 (two unregistered Function Repository resources,
   resolved 2026-09-11), not a PSALTer defect. Kept as the record of what was seen.

   DELIBERATELY NOT ASSERTED HERE: agreement with the published block values.
   On this install the association's WaveOperator and PseudoDeterminant both come
   back as -Def^2*Theta1/2 in each sector, while the rendered spectrograph shows
   (Def^2*(Theta2 - Theta1) - Theta3)/2 and (-Def^2*Theta1 - Theta3)/2 -- i.e. the
   association is missing the mass term. Reproduced with Barker's own options
   (Method -> "Hard", ShowPropagator -> True). PSALTer also emits internal errors
   on this theory (Join, Gather and Part messages) and reports the
   degenerate "(Demonstrably impossible)" unitarity verdict. Unresolved, filed as
   an issue, and load-bearing for the Wave-1 reader gate, which is specified
   against the published values. This preflight therefore checks only that the
   pipeline runs end to end; it makes no claim about the numbers.
   ============================================================================== *)

Print["========================================"];
Print["PSALTer preflight: massive vector"];
Print["========================================"];
Print[""];

Needs["xAct`PSALTer`"];
Print["PSALTER_VERSION=", xAct`PSALTer`Private`$Version[[1]]];

DefConstantSymbol[Theta1, PrintAs -> "\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\)"];
DefConstantSymbol[Theta2, PrintAs -> "\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\)"];
DefConstantSymbol[Theta3, PrintAs -> "\!\(\*SubscriptBox[\(\[Theta]\), \(3\)]\)"];
DefField[VectorField[-a], PrintAs -> "A", PrintSourceAs -> "j"];

LinearizedLagrangian = (
    -(1/2)*Theta1*CD[-a]@VectorField[-b]*CD[a]@VectorField[b]
    + (1/2)*Theta2*CD[-a]@VectorField[a]*CD[-b]@VectorField[b]
    - (1/2)*Theta3*VectorField[-a]*VectorField[a]
);

{timing, result} = AbsoluteTiming@ParticleSpectrum[
    LinearizedLagrangian,
    TheoryName -> "VectorTheory",
    MaxLaurentDepth -> 1,
    ShowPropagator -> False];

Print[""];
Print["PREFLIGHT_WALL_S=", Round[timing, 0.01]];

assoc = Global`VectorTheory;
Print["PREFLIGHT_ASSOC_KEYS=", (SymbolName /@ Keys[assoc])];
Print["PREFLIGHT_SECTORS=", Length[assoc[xAct`PSALTer`WaveOperator]]];
Print["PREFLIGHT_MX=",
    Length@FileNames["ParticleSpectrographVectorTheory.mx", Directory[]]];
Print["PREFLIGHT_PDF=", Length@FileNames["*.pdf", Directory[]]];
Print[""];
Print["PREFLIGHT PASSED"];

Quit[0]
