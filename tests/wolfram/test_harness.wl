(* ::Package:: *)
(*
   MODULE: test_harness.wl
   PURPOSE: Test framework for Wolfram pipeline modules

   USAGE:
     << "tests/wolfram/test_harness.wl";

     TestCase["Test name", expected, actual];
     TestTrue["Test name", booleanExpr];
     TestNoError["Test name", expr];

     PrintTestSummary[];
     ExitWithTestStatus[];  (* Exit[0] if all pass, Exit[1] otherwise *)

   FEATURES:
     - Track pass/fail counts
     - Pretty-print test results
     - Support for multiple assertion types
     - Exit codes for CI integration
*)

BeginPackage["TestHarness`"];

(* Public symbols *)
TestCase::usage =
  "TestCase[name, expected, actual] tests that actual equals expected.";

TestTrue::usage =
  "TestTrue[name, expr] tests that expr evaluates to True.";

TestNoError::usage =
  "TestNoError[name, expr] tests that expr evaluates without errors.";

TestApproxEqual::usage =
  "TestApproxEqual[name, expected, actual, tolerance:0.0001] tests numeric equality within tolerance.";

PrintTestSummary::usage =
  "PrintTestSummary[] prints the final test summary with pass/fail counts.";

ExitWithTestStatus::usage =
  "ExitWithTestStatus[] exits with code 0 if all tests passed, 1 otherwise.";

ResetTestCounts::usage =
  "ResetTestCounts[] resets pass/fail counters to 0.";

$PassCount::usage = "Count of passed tests.";
$FailCount::usage = "Count of failed tests.";

Begin["`Private`"];

(* Global counters *)
$PassCount = 0;
$FailCount = 0;

ResetTestCounts[] := (
  $PassCount = 0;
  $FailCount = 0;
);

(* Helper to record result *)
RecordResult[name_String, passed_?BooleanQ, msg_String:""] := (
  If[passed,
    Print["  PASS: ", name];
    $PassCount++,
    Print["  FAIL: ", name];
    If[StringLength[msg] > 0, Print["    ", msg]];
    $FailCount++
  ];
  passed
);

(* Basic equality test *)
TestCase[name_String, expected_, actual_] := Module[{passed},
  passed = (expected === actual);
  RecordResult[name, passed,
    If[!passed,
      "Expected: " <> ToString[expected, InputForm] <> "\n    Actual: " <> ToString[actual, InputForm],
      ""
    ]
  ]
];

(* Boolean test *)
TestTrue[name_String, expr_] := Module[{passed},
  passed = TrueQ[expr];
  RecordResult[name, passed,
    If[!passed, "Expression did not evaluate to True: " <> ToString[expr, InputForm], ""]
  ]
];

(* No-error test (catches $Failed and messages) *)
TestNoError[name_String, expr_] := Module[{result, passed},
  result = Quiet[Check[expr, $Failed]];
  passed = (result =!= $Failed);
  RecordResult[name, passed,
    If[!passed, "Expression raised an error or returned $Failed", ""]
  ]
];

(* Approximate equality test for numerics *)
TestApproxEqual[name_String, expected_?NumericQ, actual_?NumericQ, tolerance_:0.0001] := Module[{passed},
  passed = Abs[expected - actual] < tolerance;
  RecordResult[name, passed,
    If[!passed,
      "Expected: " <> ToString[N[expected]] <> " (tolerance: " <> ToString[tolerance] <> ")\n    Actual: " <> ToString[N[actual]],
      ""
    ]
  ]
];

(* Non-numeric case *)
TestApproxEqual[name_String, expected_, actual_, tolerance_:0.0001] := Module[{},
  RecordResult[name, False, "TestApproxEqual requires numeric arguments"]
];

(* Summary printing *)
PrintTestSummary[] := (
  Print[""];
  Print["=== Test Summary ==="];
  Print["Passed: ", $PassCount];
  Print["Failed: ", $FailCount];
  Print["Total: ", $PassCount + $FailCount];
  If[$FailCount > 0,
    Print[""];
    Print["*** SOME TESTS FAILED ***"],
    Print[""];
    Print["*** ALL TESTS PASSED ***"]
  ];
);

(* Exit with appropriate status code *)
ExitWithTestStatus[] := (
  PrintTestSummary[];
  If[$FailCount > 0, Exit[1], Exit[0]]
);

End[];
EndPackage[];
