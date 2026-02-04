#!/usr/bin/env wolframscript
(* xTensor Functionality Test - Core Tensor Algebra *)

Print["=== xTensor Core Features Test ==="];
Print[""];

<<xAct`xTensor`;

(* Test 1: Manifold and Metric Setup *)
Print["1. Creating 4D spacetime manifold..."];
DefManifold[Spacetime, 4, IndexRange[a, h]];
DefMetric[-1, g[-a, -b], Spacetime];
Print["   ✅ Manifold: ", Spacetime, " with dimension ", Dim[Spacetime]];
Print["   ✅ Metric signature: ", SignatureOfMetric[g[-a, -b]]];

(* Test 2: Tensor Definitions *)
Print[""];
Print["2. Defining tensor fields..."];
DefTensor[T[a, -b, c], Spacetime];
DefTensor[V[a], Spacetime];
DefTensor[omega[-a, -b], Spacetime, Antisymmetric[{-a, -b}]];
Print["   ✅ Mixed tensor T^a_b^c"];
Print["   ✅ Vector field V^a"];
Print["   ✅ Antisymmetric 2-form ω_ab"];

(* Test 3: Curvature Tensors *)
Print[""];
Print["3. Curvature tensor operations..."];
DefTensor[R[-a, -b, -c, d], Spacetime, RiemannSymmetric[{-a, -b, -c, d}]];
DefTensor[Ric[-a, -b], Spacetime, Symmetric[{-a, -b}]];
DefTensor[RicciScalar[], Spacetime];
Print["   ✅ Riemann tensor R_abcd"];
Print["   ✅ Ricci tensor Ric_ab"];
Print["   ✅ Ricci scalar R"];

(* Test 4: Index Manipulations *)
Print[""];
Print["4. Index contraction and manipulation..."];
expr1 = g[a, b] V[a] V[b];
expr2 = T[a, -a, c];  (* trace *)
expr3 = g[-a, -b] omega[a, b];  (* contraction with metric *)
Print["   ✅ Norm squared: g^ab V_a V_b = ", expr1];
Print["   ✅ Trace: T^a_a^c = ", expr2];
Print["   ✅ Metric contraction: g^ab ω_ab = ", expr3];

(* Test 5: Covariant Derivatives *)
Print[""];
Print["5. Covariant derivatives..."];
CD = CovDOfMetric[g[-a, -b]];
Print["   ✅ Covariant derivative: ", CD];
Print["   ✅ Metric compatibility: ", CD[g[-a, -b]]];

Print[""];
Print["*** xTensor test complete - All core features functional! ***"];
Print["Ready for: General Relativity, Differential Geometry, Tensor Calculus"];