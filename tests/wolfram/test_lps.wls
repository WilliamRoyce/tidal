#!/usr/bin/env wolframscript
(* Test suite for PerturbativeReduction.wl (Lagrangian Perturbative Substitution) *)

Print["=== PerturbativeReduction (LPS) Test Suite ==="];
Print[""];

Get[FileNameJoin[{DirectoryName[$InputFileName], "test_harness.wl"}]];

pipelinePath = FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "..", "tidal", "wolfram"
}];

Get[FileNameJoin[{pipelinePath, "PerturbativeReduction.wl"}]];

ReduceLagrangian = TorsionGertsenshtein`PerturbativeReduction`ReduceLagrangian;
DetectMaxTimeOrder = TorsionGertsenshtein`PerturbativeReduction`DetectMaxTimeOrder;
ExtractBaseEquations = TorsionGertsenshtein`PerturbativeReduction`ExtractBaseEquations;

(* === Test 1: DetectMaxTimeOrder === *)
Print["Test 1: DetectMaxTimeOrder"];

(* Lagrangian: 1/2 q'^2 + b5 (q'')^2  →  max time-order is 2 *)
phi[t_, x_] := phi[t, x];
expr1 = (1/2)*Derivative[1, 0][phi][t, x]^2 + b5*Derivative[2, 0][phi][t, x]^2;
TestCase[
  "DetectMaxTimeOrder finds order 2 in (q'')^2",
  2,
  DetectMaxTimeOrder[expr1, {phi}]
];

(* Lagrangian with 1st-order only *)
expr2 = (1/2)*Derivative[1, 0][phi][t, x]^2 - phi[t, x]^2;
TestCase[
  "DetectMaxTimeOrder finds order 1 for 1st-order L",
  1,
  DetectMaxTimeOrder[expr2, {phi}]
];

(* No time derivatives *)
expr3 = phi[t, x]^2 + Derivative[0, 1][phi][t, x]^2;
TestCase[
  "DetectMaxTimeOrder finds order 0 when no time derivatives",
  0,
  DetectMaxTimeOrder[expr3, {phi}]
];

(* === Test 2: ExtractBaseEquations === *)
Print[""];
Print["Test 2: ExtractBaseEquations"];

(* Single equation: q'' + 2*b5*q'''' - V'(q) = 0 *)
eq1 = {phi, Derivative[2, 0][phi][t, x] + 2*b5*Derivative[4, 0][phi][t, x]
  - phi[t, x]^3};
baseEqs = ExtractBaseEquations[{eq1}, {"b5"}];
TestCase[
  "ExtractBaseEquations sets small_param to 0",
  Derivative[2, 0][phi][t, x] - phi[t, x]^3,
  baseEqs[[1, 2]]
];

(* === Test 3: No-op for non-perturbative theory === *)
Print[""];
Print["Test 3: ReduceLagrangian no-op for empty smallParams"];

lagSimple = (1/2)*Derivative[1, 0][phi][t, x]^2 - (1/2)*phi[t, x]^2;
fieldEqs = {{"phi", Derivative[2, 0][phi][t, x] + phi[t, x]}};
compToFunc = <|"phi" -> phi|>;
result = ReduceLagrangian[lagSimple, fieldEqs, compToFunc, {}, 1];
TestCase[
  "Empty smallParams: input returned unchanged",
  lagSimple,
  result
];

(* === Test 4: No-op for already-1st-order Lagrangian === *)
Print[""];
Print["Test 4: ReduceLagrangian no-op for 1st-order L"];

lag1stOrder = (1/2)*Derivative[1, 0][phi][t, x]^2 - (1/2)*phi[t, x]^2;
result = ReduceLagrangian[lag1stOrder, fieldEqs, compToFunc, {"b5"}, 1];
TestCase[
  "Already 1st-order: input returned unchanged",
  lag1stOrder,
  result
];

(* === Test 5: Pais-Uhlenbeck single-field substitution === *)
Print[""];
Print["Test 5: Pais-Uhlenbeck-style L with dynamical higher-derivative field"];

(* L = 1/2 q'^2 - 1/2 V0(q) - b5 (q'')^2                                   *)
(* The b5 term is the perturbative correction.  At b5=0, the EOM is        *)
(*   q'' = -V0'(q)                                                          *)
(* (assuming V0 is the order-0 potential).  LPS substitutes q'' in the b5  *)
(* term using this rule, yielding L_red = 1/2 q'^2 - 1/2 V0(q) - b5 V0'(q)^2. *)
(* For testing, use V0(q) = q (linear), so V0'(q) = 1 and (q'')^2 -> 1.    *)

V0[q_] := q;
lagPU = (1/2)*Derivative[1, 0][phi][t, x]^2 - (1/2)*V0[phi[t, x]]
  - b5*Derivative[2, 0][phi][t, x]^2;

(* The order-0 EOM (Euler-Lagrange of L at b5=0) for V0(q)=q is:           *)
(*   d/dt(dL/dq') - dL/dq = 0  ->  q'' + 1/2 = 0  ->  q'' = -1/2           *)
(* Wait — for V0(q)=q, dL/dq = -1/2, so q'' = -1/2.                        *)
fieldEqsPU = {{"phi",
  Derivative[2, 0][phi][t, x] + (1/2)
  + 2*b5*Derivative[4, 0][phi][t, x]
}};
compToFuncPU = <|"phi" -> phi|>;

result = ReduceLagrangian[lagPU, fieldEqsPU, compToFuncPU, {"b5"}, 1];
(* After substitution: q''^2 -> (-1/2)^2 = 1/4; the b5*(q'')^2 term        *)
(* becomes b5*1/4 (a constant — passes Series-truncation at order 1).      *)
expected = (1/2)*Derivative[1, 0][phi][t, x]^2 - (1/2)*phi[t, x] - b5*(1/4);
TestCase[
  "Pais-Uhlenbeck reduction matches analytical result",
  Simplify[expected - result],
  0
];

(* Verify resulting L has time-order <= 1 *)
TestCase[
  "Reduced L has time-order <= 1",
  True,
  DetectMaxTimeOrder[result, {phi}] <= 1
];

(* === Test 6: Constraint-promotion case throws clearly === *)
Print[""];
Print["Test 6: Constraint-promotion case (LPS Throws clear error)"];

(* L has q'' but the order-0 EOM for phi is purely algebraic (no d_t^2  *)
(* of phi at b5=0).  This is the constraint-promotion case: phi's       *)
(* kinetic coefficient is proportional to b5, so at b5=0 phi satisfies  *)
(* an algebraic constraint rather than a wave equation.                  *)
(*                                                                       *)
(* The standard substitution-based LPS does not apply here.  LPS Throws *)
(* with a clear message naming the affected field(s) so the issue can   *)
(* be investigated.                                                      *)

p[t_, x_] := p[t, x];
lagConstraint = -1/2*phi[t, x]*p[t, x] - (1/2)*b5*Derivative[2, 0][phi][t, x]^2;
fieldEqsConstraint = {
  {"phi", -p[t, x] - b5*Derivative[4, 0][phi][t, x]},
  {"p", -phi[t, x]}
};
compToFuncConstraint = <|"phi" -> phi, "p" -> p|>;

(* Catch the Throw and verify its message names the constraint-promotion *)
caught = Catch[
  ReduceLagrangian[lagConstraint, fieldEqsConstraint,
    compToFuncConstraint, {"b5"}, 1]
];
TestCase[
  "Constraint-promotion case: throws with descriptive message",
  True,
  StringQ[caught] && StringContainsQ[caught, "constraint-promotion"]
];

(* === Summary === *)
Print[""];
PrintTestSummary[];
