#!/usr/bin/env wolframscript
(* Tests for CanonicalizeExpFraction (GH #378).
   Separated from test_export_json.wls because that suite hits a pre-existing
   uncaught Throw partway through and would prevent these tests from running. *)

Print["=== CanonicalizeExpFraction Test Suite ==="];
Print[""];

Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

pipelinePath = FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "..", "tidal", "wolfram"
}];

<< xAct`xTensor`;
<< xAct`xCoba`;

Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

cef = TorsionGertsenshtein`ExportJSON`CanonicalizeExpFraction;

(* --- No-op cases (no Power[E, _]) --- *)
TestCase["no E^ is no-op (polynomial)", (1 + a)/(1 - a), cef[(1 + a)/(1 - a)]];
TestCase["no E^ is no-op (constant)",   42,              cef[42]];

(* --- Single E^ in denominator: passthrough; this case is handled
       downstream by Python's _invert_exp_denominator. --- *)
TestCase["single E^ in denom is passthrough (handled by Python)",
  1/E^(x^2),
  cef[1/E^(x^2)]];

(* --- The EH dual-Gaussian pattern (equations[1].rhs.terms[7], GH #378) ---
   Use x[] (xCoba coordinate function), which is how the real pipeline
   serializes coordinates — bare `x` masked an earlier bug where Cancel
   silently failed for x[] but worked for x. *)
ehCoeff = E^((zc1^2 + zc2^2 + x[]^2)/sigB^2) /
  (E^((zc1^2 + zc2^2 + x[]^2)/sigB^2)
   - 2*Bpeak^2*E^((zc2^2 + 2*zc1*x[])/sigB^2)*(rho - 8*sigma)
   - 2*Bpeak^2*E^((zc1^2 + 2*zc2*x[])/sigB^2)*(rho - 8*sigma)
   + 4*Bpeak^2*E^((zc1^2 + zc2^2 + 2*zc1*x[] + 2*zc2*x[])/(2*sigB^2))*(rho - 8*sigma));

ehCoeffCanon = cef[ehCoeff];

(* The numerator E^((zc1^2+zc2^2+x[]^2)/sigB^2) should be gone — that's the
   factor that overflows at box edges. *)
TestTrue["EH: numerator E^A factored out (numerator becomes 1)",
  Numerator[ehCoeffCanon] === 1];

ehParams = {x[] -> Null, Bpeak -> 0.01, sigB -> 5.0, zc1 -> 50.0,
            zc2 -> 150.0, rho -> 0.01, sigma -> 0.01};
ehParams = DeleteCases[ehParams, x[] -> _];  (* x[] substituted per-test below *)

(* At x = zc1, original form still evaluates (exp(2500) is huge but finite ratio). *)
TestTrue["EH: numeric agreement at x = 50 (Gaussian centre)",
  Module[{a, b},
    a = N[ehCoeff /. Join[{x[] -> 50.0}, ehParams], 30];
    b = N[ehCoeffCanon /. Join[{x[] -> 50.0}, ehParams], 30];
    Abs[(a - b)/a] < 1.0*^-10
  ]];

(* At x = 100 (box edge for L=200 centred at 0), un-canonicalised form
   overflows machine doubles; canonical form must remain finite and close
   to 1 (no Gaussian support there). *)
TestTrue["EH: canonical form finite at box edge x=100",
  Module[{b},
    b = N[ehCoeffCanon /. Join[{x[] -> 100.0}, ehParams], 30];
    NumericQ[b] && Abs[b - 1.0] < 1.0*^-6
  ]];

TestTrue["EH: canonical form finite at box edge x=-100",
  Module[{b},
    b = N[ehCoeffCanon /. Join[{x[] -> -100.0}, ehParams], 30];
    NumericQ[b] && Abs[b - 1.0] < 1.0*^-6
  ]];

(* Machine-double round-trip: canonical form must not produce NaN/Inf at edges
   under MachinePrecision (this is what trips the Python evaluator). *)
(* --- Wrapped denominator: outer Times around the Plus (terms 0-6 in
       the EH JSON have shape (C * E^B) / (sigB^2 * (E^A - ... )). --- *)
ehWrapped = (4*Bpeak^2*E^((zc1^2 + zc2^2 + 2*zc1*x[] + 2*zc2*x[])/(2*sigB^2))*rho*zc1) /
  (sigB^2 * (E^((zc1^2 + zc2^2 + x[]^2)/sigB^2)
   - 2*Bpeak^2*E^((zc2^2 + 2*zc1*x[])/sigB^2)*(rho - 8*sigma)
   - 2*Bpeak^2*E^((zc1^2 + 2*zc2*x[])/sigB^2)*(rho - 8*sigma)
   + 4*Bpeak^2*E^((zc1^2 + zc2^2 + 2*zc1*x[] + 2*zc2*x[])/(2*sigB^2))*(rho - 8*sigma)));

ehWrappedCanon = cef[ehWrapped];

TestTrue["EH wrapped: MachinePrecision at x=100 is finite",
  Module[{b},
    b = ehWrappedCanon /. Join[{x[] -> 100.0}, ehParams];
    NumericQ[b] && b < Infinity && b > -Infinity && !(b === Indeterminate)
  ]];

(* --- equations[2,3] pattern: dominant exp is (zc1-x[])²/sigB² (same
       degree-2 in x[] as the negative-leading -(zc1-x[])²/sigB²; the
       sign-aware ranking must pick the positive one). --- *)
ehEq2Term7 = -(E^((zc1 - x[])^2/sigB^2) /
  (-E^((zc1 - x[])^2/sigB^2)
   + 2*Bpeak^2*E^((zc1 - x[])^2/sigB^2)*
     (E^(-((zc1 - x[])^2/sigB^2))
      + E^(-((zc2 - x[])^2/sigB^2))
      - 2/E^((zc1^2 + zc2^2 - 2*zc1*x[] - 2*zc2*x[] + 2*x[]^2)/(2*sigB^2)))
     *rho));

ehEq2Canon = cef[ehEq2Term7];

TestTrue["EH eq[2,3]: MachinePrecision at x=100 is finite",
  Module[{b},
    b = ehEq2Canon /. Join[{x[] -> 100.0}, ehParams];
    NumericQ[b] && b < Infinity && b > -Infinity && !(b === Indeterminate)
  ]];

TestTrue["EH eq[2,3]: MachinePrecision at x=-100 is finite",
  Module[{b},
    b = ehEq2Canon /. Join[{x[] -> -100.0}, ehParams];
    NumericQ[b] && b < Infinity && b > -Infinity && !(b === Indeterminate)
  ]];

TestTrue["EH wrapped: MachinePrecision at x=-100 is finite",
  Module[{b},
    b = ehWrappedCanon /. Join[{x[] -> -100.0}, ehParams];
    NumericQ[b] && b < Infinity && b > -Infinity && !(b === Indeterminate)
  ]];

TestTrue["EH: MachinePrecision evaluation at x=100 is finite",
  Module[{b},
    b = ehCoeffCanon /. Join[{x[] -> 100.0}, ehParams];
    NumericQ[b] && b < Infinity && b > -Infinity && !(b === Indeterminate)
  ]];

Print[""];
ExitWithTestStatus[];
