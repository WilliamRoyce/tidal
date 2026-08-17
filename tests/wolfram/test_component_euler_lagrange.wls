#!/usr/bin/env wolframscript
(* Tests for ComponentEulerLagrange's "Measure" option (GH #394).

   Standard Euler-Lagrange is only valid on a Lagrangian DENSITY sqrt|g| * L.
   Before #394 the pipeline varied the scalar component Lagrangian directly, so
   any first-derivative term arising from sqrt|g| was lost (spatially
   curvilinear metrics) or sign-flipped (time-dependent metrics).

   These cases are the A/B comparisons from the #394 investigation, checked
   against analytic ground truth. *)

Print["=== ComponentEulerLagrange Measure Test Suite (GH #394) ==="];
Print[""];

Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

pipelinePath = FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "..", "tidal", "wolfram"
}];

<< xAct`xTensor`;
<< xAct`xCoba`;
$DefInfoQ = False;

Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];

(* Coefficient of a given derivative in an EOM, after Expand. *)
coeffOf[eom_, deriv_] := Simplify[Coefficient[Expand[Simplify[eom]], deriv]];

(* ------------------------------------------------------------------ *)
(* 1. Polar: g = diag(-1, 1, x^2), sqrt|g| = x                         *)
(*    Correct radial Laplacian carries (1/x) d_x.                      *)
(* ------------------------------------------------------------------ *)

polarCoords = {t, x, y};
polarL = -1/2 (
    -D[pf[t, x, y], t]^2
     + D[pf[t, x, y], x]^2
     + (1/x^2) D[pf[t, x, y], y]^2
  ) - pm2/2 pf[t, x, y]^2;

polarNoMeasure = ComponentEulerLagrange[polarL, {pf}, polarCoords][[1]];
polarWithMeasure = ComponentEulerLagrange[polarL, {pf}, polarCoords,
                     "Measure" -> x][[1]];

TestTrue["polar: WITHOUT measure the (1/x) d_x term is absent (the #394 bug)",
  coeffOf[polarNoMeasure, Derivative[0, 1, 0][pf][t, x, y]] === 0];

TestCase["polar: WITH measure x, d_x coefficient is 1/x",
  1/x,
  coeffOf[polarWithMeasure, Derivative[0, 1, 0][pf][t, x, y]]];

TestCase["polar: second-derivative terms are unaffected by the measure",
  coeffOf[polarNoMeasure, Derivative[0, 2, 0][pf][t, x, y]],
  coeffOf[polarWithMeasure, Derivative[0, 2, 0][pf][t, x, y]]];

(* ------------------------------------------------------------------ *)
(* 2. Reduced radial (spherical_kg_1d): sqrt|g| factors to r^2.        *)
(*    Correct EOM is d2_t = d2_r + (2/r) d_r - m^2.                     *)
(* ------------------------------------------------------------------ *)

radialCoords = {t, x};
radialL = -1/2 (-D[rf[t, x], t]^2 + D[rf[t, x], x]^2) - rm2/2 rf[t, x]^2;

radialWithMeasure = ComponentEulerLagrange[radialL, {rf}, radialCoords,
                      "Measure" -> x^2][[1]];

TestCase["radial: WITH measure r^2, d_r coefficient is 2/r",
  2/x,
  coeffOf[radialWithMeasure, Derivative[0, 1][rf][t, x]]];

(* Killed-coordinate factors are constants w.r.t. the surviving coordinates,
   so an unfactored measure must give the identical result. *)
radialUnfactored = ComponentEulerLagrange[radialL, {rf}, radialCoords,
                     "Measure" -> x^2 Sin[yy]][[1]];

TestCase["radial: unfactored measure r^2 Sin[y] agrees with factored r^2",
  0,
  Simplify[radialWithMeasure - radialUnfactored]];

(* ------------------------------------------------------------------ *)
(* 3. de Sitter: g = a^2 eta with a = Exp[H t], sqrt|g| = Exp[3 H t].   *)
(*    Ground truth (covariant Box phi, and the standard FRW result      *)
(*    (D-2)H in D=3) gives friction coefficient -H.  WITHOUT the        *)
(*    measure the pipeline produced +2H -- anti-damping.                *)
(* ------------------------------------------------------------------ *)

dsCoords = {t, x, y};
dsL = -1/2 (
    -(1/Exp[2 dH t]) D[df[t, x, y], t]^2
     + (1/Exp[2 dH t]) (D[df[t, x, y], x]^2 + D[df[t, x, y], y]^2)
  );

dsNoMeasure = ComponentEulerLagrange[dsL, {df}, dsCoords][[1]];
dsWithMeasure = ComponentEulerLagrange[dsL, {df}, dsCoords,
                  "Measure" -> Exp[3 dH t]][[1]];

TestCase["de Sitter: WITH measure, friction coefficient is -H (ground truth)",
  -dH,
  Simplify[
    coeffOf[dsWithMeasure, Derivative[1, 0, 0][df][t, x, y]] /
    coeffOf[dsWithMeasure, Derivative[2, 0, 0][df][t, x, y]] * (-1)]];

TestTrue["de Sitter: WITHOUT measure the friction sign is wrong (the #394 bug)",
  Module[{ratio},
    ratio = Simplify[
      coeffOf[dsNoMeasure, Derivative[1, 0, 0][df][t, x, y]] /
      coeffOf[dsNoMeasure, Derivative[2, 0, 0][df][t, x, y]] * (-1)];
    Simplify[ratio - (-dH)] =!= 0]];

(* ------------------------------------------------------------------ *)
(* 4. Control: a constant measure must be a no-op.                     *)
(*    This is what keeps flat and constant-conformal metrics on a       *)
(*    bit-identical path -- the guard is "free of coordinates", not     *)
(*    "equal to 1".                                                     *)
(* ------------------------------------------------------------------ *)

flatCoords = {t, x};
flatL = -1/2 (-D[ff[t, x], t]^2 + D[ff[t, x], x]^2) - fm2/2 ff[t, x]^2;

TestCase["control: constant measure 4 gives the same EOM as no measure",
  0,
  Simplify[
    ComponentEulerLagrange[flatL, {ff}, flatCoords][[1]]
    - ComponentEulerLagrange[flatL, {ff}, flatCoords, "Measure" -> 4][[1]]]];

TestCase["control: measure 1 gives the same EOM as no measure",
  0,
  Simplify[
    ComponentEulerLagrange[flatL, {ff}, flatCoords][[1]]
    - ComponentEulerLagrange[flatL, {ff}, flatCoords, "Measure" -> 1][[1]]]];

Print[""];
ExitWithTestStatus[];
