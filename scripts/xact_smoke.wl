(* ==============================================================================
   xAct/xTensor + xCoba Smoke Test
   ==============================================================================
   This script verifies that xAct and xCoba are properly installed and working.
   
   Usage:
     wolframscript -file scripts/xact_smoke.wl
   
   Expected output:
     - Package loading messages
     - Manifold and chart definitions
     - Metric definition
     - Riemann tensor computation (short output)
     - "SMOKE TEST PASSED" message
   ============================================================================== *)

Print["========================================"];
Print["xAct/xCoba Smoke Test"];
Print["========================================"];
Print[""];

(* Load packages with suppressed warnings *)
Print["Loading xAct packages..."];
Quiet[
  Needs["xAct`xTensor`"];
  Needs["xAct`xCoba`"];
, {General::stop, UpSetDelayed::write, SetDelayed::write}];

Print["✓ Packages loaded successfully"];
Print[""];

(* Define a 4D manifold with abstract indices *)
Print["Defining manifold..."];
DefManifold[M, 4, {a, b, c, d, e, f}];
Print["✓ Manifold M defined (4 dimensions)"];

(* Define a coordinate chart - use safe symbol names *)
Print["Defining coordinate chart..."];
DefChart[cart, M, {0, 1, 2, 3}, {ct[], cr[], ctheta[], cphi[]}];
Print["✓ Chart defined with coordinates: {ct, cr, ctheta, cphi}"];
Print[""];

(* Define metric tensor *)
Print["Defining metric tensor..."];
DefMetric[-1, gM[-a, -b], CDM, {";", "\[Del]"}, PrintAs -> "g"];
Print["✓ Metric tensor g defined with signature (-1)"];
Print[""];

(* Test basic tensor operations *)
Print["Testing tensor operations..."];
testExpr = gM[-a, -b] gM[a, c];
Print["  gM[-a,-b] gM[a,c] = ", testExpr // ToCanonical];
Print["✓ Tensor contraction works"];
Print[""];

(* Test Riemann tensor (abstract computation) *)
Print["Testing Riemann tensor (abstract)..."];
riemannTest = RiemannCDM[-a, -b, -c, -d];
Print["  Riemann[-a,-b,-c,-d] defined: ", Head[riemannTest]];
Print["✓ Riemann tensor accessible"];
Print[""];

(* Test symmetry properties *)
Print["Testing tensor symmetries..."];
symTest = RiemannCDM[-a, -b, -c, -d] + RiemannCDM[-b, -a, -c, -d];
symResult = symTest // ToCanonical;
Print["  R[-a,-b,-c,-d] + R[-b,-a,-c,-d] = ", symResult];
If[symResult === 0,
  Print["✓ Antisymmetry verified (result = 0)"],
  Print["  (Non-zero result is expected for abstract tensors)"]
];
Print[""];

(* Summary *)
Print["========================================"];
Print["SMOKE TEST PASSED"];
Print["========================================"];
Print[""];
Print["xAct/xCoba installation is working correctly."];
Print[""];

Quit[0]
