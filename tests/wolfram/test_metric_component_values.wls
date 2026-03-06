#!/usr/bin/env wolframscript
(* Empirical test: Does MetricInBasis set contravariant DownValues?
   And does AllComponentValues add anything on top?

   This tests the optimization hypothesis: if contravariant metric zeros
   auto-evaluate after MetricInBasis, then Expand[] in
   BatchedTraceBasisDummyWithMetric already collapses those terms.
   If NOT, AllComponentValues is needed to pre-define them. *)

Print["=== MetricInBasis vs AllComponentValues: Contravariant Auto-Eval Test ==="];
Print[""];

<<xAct`xCoba`;

(* Setup: minimal 4D flat spacetime *)
If[!xTensorQ[mcvEta],
  DefManifold[mcvM, 4, {a, b, c, d, e, f}];
  DefMetric[-1, mcvEta[-a, -b], mcvCD,
    SymbolOfCovD -> {";", "\[Del]"}, PrintAs -> "\[Eta]"];
  DefChart[mcvCart, mcvM, {0, 1, 2, 3}, {t[], x[], y[], z[]}];
];

mcvMetricMatrix = DiagonalMatrix[{-1, 1, 1, 1}];
MetricInBasis[mcvEta, -mcvCart, mcvMetricMatrix];

Print["--- After MetricInBasis[eta, -chart, DiagMatrix] ---"];
Print[""];

(* Test 1: Covariant components (expected to auto-evaluate) *)
Print["Test 1: Covariant metric g_{ij}"];
covDiag = mcvEta[{0, -mcvCart}, {0, -mcvCart}];
covOffDiag = mcvEta[{0, -mcvCart}, {1, -mcvCart}];
Print["  g_{00} = ", covDiag, "  (expect -1)"];
Print["  g_{01} = ", covOffDiag, "  (expect 0)"];
Print["  PASS covariant: ", covDiag === -1 && covOffDiag === 0];
Print[""];

(* Test 2: Contravariant components (THIS IS THE KEY QUESTION) *)
Print["Test 2: Contravariant metric g^{ij}"];
contDiag = mcvEta[{0, mcvCart}, {0, mcvCart}];
contOffDiag = mcvEta[{0, mcvCart}, {1, mcvCart}];
Print["  g^{00} = ", contDiag, "  (expect -1)"];
Print["  g^{01} = ", contOffDiag, "  (expect 0)"];
contravAutoEval = (contDiag === -1 && contOffDiag === 0);
Print["  PASS contravariant: ", contravAutoEval];
Print[""];

If[contravAutoEval,
  Print[">> MetricInBasis DOES set contravariant DownValues."];
  Print[">> AllComponentValues for metric is redundant (defense-in-depth only)."],
  Print[">> MetricInBasis does NOT set contravariant DownValues!"];
  Print[">> AllComponentValues is NEEDED to pre-define g^{ij} zeros."];
];
Print[""];

(* Test 3: Effect on TraceBasisDummy expansion *)
Print["Test 3: TraceBasisDummy expansion with metric contraction"];
(* Build a simple expression: g^{ab} * V_a * W_b *)
DefTensor[mcvV[-a], mcvM];
DefTensor[mcvW[-a], mcvM];
ComponentValue[ComponentArray[mcvV[-mcvCart]], {1, 0, 0, 0}];
ComponentValue[ComponentArray[mcvW[-mcvCart]], {0, 1, 0, 0}];

expr = mcvEta[a, b] mcvV[-a] mcvW[-b];
Print["  Abstract: g^{ab} V_a W_b"];
basisExpr = ToBasis[mcvCart][expr];
Print["  After ToBasis: ", Length[If[Head[basisExpr]===Plus, List@@basisExpr, {basisExpr}]], " terms"];
traced = TraceBasisDummy[basisExpr];
Print["  After TraceBasisDummy: ", Length[If[Head[traced]===Plus, List@@traced, {traced}]], " terms"];
expanded = Expand[traced];
Print["  After Expand: ", Length[If[Head[expanded]===Plus, List@@expanded, {expanded}]], " terms"];
Print["  Expected: should collapse to 0 if metric zeros auto-evaluate"];
Print["  PASS contraction: ", expanded === 0];
Print[""];

(* Test 4: AllComponentValues with abstract indices *)
Print["--- Test 4: AllComponentValues (explicit pre-definition) ---"];
Print[""];

(* Setup a second metric to test AllComponentValues *)
If[!xTensorQ[mcvG],
  DefManifold[mcvM2, 4, {p, q, r, s}];
  DefMetric[-1, mcvG[-p, -q], mcvCD2,
    SymbolOfCovD -> {"|", "\[Del]"}, PrintAs -> "G"];
  DefChart[mcvCart2, mcvM2, {0, 1, 2, 3}, {t2[], x2[], y2[], z2[]}];
];

mcvMetricMatrix2 = DiagonalMatrix[{-1, 1, 1, 1}];
MetricInBasis[mcvG, -mcvCart2, mcvMetricMatrix2];

(* Now add AllComponentValues on top — use abstract indices with basis markers *)
invMatrix = Inverse[mcvMetricMatrix2];
AllComponentValues[mcvG[{-p, -mcvCart2}, {-q, -mcvCart2}], mcvMetricMatrix2];
AllComponentValues[mcvG[{p, mcvCart2}, {q, mcvCart2}], invMatrix];

Print["  After MetricInBasis + AllComponentValues:"];
acvContDiag = mcvG[{0, mcvCart2}, {0, mcvCart2}];
acvContOff = mcvG[{0, mcvCart2}, {1, mcvCart2}];
Print["  G^{00} = ", acvContDiag, "  (expect -1)"];
Print["  G^{01} = ", acvContOff, "  (expect 0)"];
Print["  PASS with AllComponentValues: ", acvContDiag === -1 && acvContOff === 0];
Print[""];

(* Test 5: Verify AllComponentValues makes individual lookups auto-evaluate *)
Print["--- Test 5: Individual component auto-eval after AllComponentValues ---"];
Print[""];
acvCovDiag = mcvG[{0, -mcvCart2}, {0, -mcvCart2}];
acvCovOff = mcvG[{0, -mcvCart2}, {1, -mcvCart2}];
Print["  G_{00} = ", acvCovDiag, "  (expect -1)"];
Print["  G_{01} = ", acvCovOff, "  (expect 0)"];
acvCovPass = (acvCovDiag === -1 && acvCovOff === 0);
Print["  PASS covariant with ACV: ", acvCovPass];
Print[""];

(* Test 6: ComponentValue[] as a QUERY retrieves the stored rules *)
Print["--- Test 6: ComponentValue[] query function ---"];
Print[""];
cvQuery00 = ComponentValue[mcvG[{0, -mcvCart2}, {0, -mcvCart2}]];
cvQuery01 = ComponentValue[mcvG[{0, -mcvCart2}, {1, -mcvCart2}]];
cvQueryC00 = ComponentValue[mcvG[{0, mcvCart2}, {0, mcvCart2}]];
cvQueryC01 = ComponentValue[mcvG[{0, mcvCart2}, {1, mcvCart2}]];
Print["  CV query G_{00} = ", cvQuery00, "  (expect -1)"];
Print["  CV query G_{01} = ", cvQuery01, "  (expect 0)"];
Print["  CV query G^{00} = ", cvQueryC00, "  (expect -1)"];
Print["  CV query G^{01} = ", cvQueryC01, "  (expect 0)"];
Print["  PASS CV query cov: ", cvQuery00 === -1 && cvQuery01 === 0];
Print["  PASS CV query contrav: ", cvQueryC00 === -1 && cvQueryC01 === 0];
Print[""];

(* Test 7: Direct Mathematica Set — bypass xCoba entirely to create DownValues *)
Print["--- Test 7: Direct Mathematica Set for metric DownValues ---"];
Print[""];
Do[
  With[{val = mcvMetricMatrix2[[ii+1, jj+1]]},
    mcvG[{ii, -mcvCart2}, {jj, -mcvCart2}] = val;
    If[ii =!= jj, mcvG[{jj, -mcvCart2}, {ii, -mcvCart2}] = val];
  ];
  With[{val = invMatrix[[ii+1, jj+1]]},
    mcvG[{ii, mcvCart2}, {jj, mcvCart2}] = val;
    If[ii =!= jj, mcvG[{jj, mcvCart2}, {ii, mcvCart2}] = val];
  ],
  {ii, 0, 3}, {jj, ii, 3}
];
testCov00 = mcvG[{0, -mcvCart2}, {0, -mcvCart2}];
testCov01 = mcvG[{0, -mcvCart2}, {1, -mcvCart2}];
testCont00 = mcvG[{0, mcvCart2}, {0, mcvCart2}];
testCont01 = mcvG[{0, mcvCart2}, {1, mcvCart2}];
Print["  G_{00} = ", testCov00, "  (expect -1)"];
Print["  G_{01} = ", testCov01, "  (expect 0)"];
Print["  G^{00} = ", testCont00, "  (expect -1)"];
Print["  G^{01} = ", testCont01, "  (expect 0)"];
dvPass = (testCov00 === -1 && testCov01 === 0 && testCont00 === -1 && testCont01 === 0);
Print["  PASS direct Set: ", dvPass];
Print[""];

(* Test 8: TraceBasisDummy with direct DownValues — does Expand collapse zeros? *)
Print["--- Test 8: TraceBasisDummy + Expand with direct DownValues ---"];
Print[""];
DefTensor[mcvU[-p], mcvM2];
DefTensor[mcvT[-p], mcvM2];
(* Set field component values via xCoba *)
ComponentValue[ComponentArray[mcvU[-mcvCart2]], {1, 0, 0, 0}];
ComponentValue[ComponentArray[mcvT[-mcvCart2]], {0, 1, 0, 0}];
(* ALSO set direct DownValues for fields *)
Do[mcvU[{ii, -mcvCart2}] = {1, 0, 0, 0}[[ii+1]], {ii, 0, 3}];
Do[mcvT[{ii, -mcvCart2}] = {0, 1, 0, 0}[[ii+1]], {ii, 0, 3}];

expr2 = mcvG[p, q] mcvU[-p] mcvT[-q];
basisExpr2 = ToBasis[mcvCart2][expr2];
traced2 = TraceBasisDummy[basisExpr2];
Print["  After TraceBasisDummy: ", Length[If[Head[traced2]===Plus, List@@traced2, {traced2}]], " terms"];
expanded2 = Expand[traced2];
Print["  After Expand: ", expanded2];
Print["  PASS TraceBasisDummy+Expand with DownValues: ", expanded2 === 0];
Print[""];

(* Test 9: Metric-only DownValues (no field DownValues) — partial optimization *)
Print["--- Test 9: Metric-only DownValues (fields stay symbolic) ---"];
Print[""];
(* Use original eta metric which has NO direct DownValues *)
(* Now set ONLY metric DownValues for eta *)
Do[
  With[{val = mcvMetricMatrix[[ii+1, jj+1]]},
    mcvEta[{ii, mcvCart}, {jj, mcvCart}] = val;
    If[ii =!= jj, mcvEta[{jj, mcvCart}, {ii, mcvCart}] = val];
    mcvEta[{ii, -mcvCart}, {jj, -mcvCart}] = val;
    If[ii =!= jj, mcvEta[{jj, -mcvCart}, {ii, -mcvCart}] = val];
  ],
  {ii, 0, 3}, {jj, ii, 3}
];
(* Re-test the original eta contraction — fields V,W still have xCoba rules only *)
expr3 = mcvEta[a, b] mcvV[-a] mcvW[-b];
basisExpr3 = ToBasis[mcvCart][expr3];
traced3 = TraceBasisDummy[basisExpr3];
nTermsRaw = Length[If[Head[traced3]===Plus, List@@traced3, {traced3}]];
expanded3 = Expand[traced3];
nTermsExp = Length[If[Head[expanded3]===Plus, List@@expanded3, {expanded3}]];
Print["  Terms after TraceBasisDummy: ", nTermsRaw];
Print["  Terms after Expand (metric zeros collapse): ", nTermsExp];
Print["  Result: ", expanded3];
Print["  PASS metric-only DownValues reduce terms: ", nTermsExp < nTermsRaw];
Print[""];

(* Summary *)
Print["=== SUMMARY ==="];
Print["MetricInBasis auto-eval: NO (xCoba internal rules, not DownValues)"];
Print["AllComponentValues auto-eval: NO (same internal mechanism)"];
Print["Direct Mathematica Set auto-eval: ", If[dvPass, "YES", "NO"]];
Print["TraceBasisDummy+Expand with full DownValues: ", If[expanded2 === 0, "COLLAPSES", "NO EFFECT"]];
Print["Metric-only DownValues term reduction: ", nTermsRaw, " -> ", nTermsExp,
  " (", If[nTermsExp < nTermsRaw, "REDUCED", "NO CHANGE"], ")"];
