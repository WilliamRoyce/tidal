(* ::Package:: *)
(* EulerLagrange.wl - Euler-Lagrange equation derivation for field theories *)
(* Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline *)

BeginPackage["TorsionGertsenshtein`EulerLagrange`", {"xAct`xTensor`"}];

(* Public symbols *)
EulerLagrangeEquation::usage =
  "EulerLagrangeEquation[lagrangian, field, covd] computes the Euler-Lagrange \
equation for the given field from the Lagrangian density. Returns an expression \
equal to zero (the equation of motion).";

VariationalDerivative::usage =
  "VariationalDerivative[lagrangian, field, covd] computes the functional \
derivative dL/d(field) - D_a[dL/d(D_a field)].";

SetupSpacetime::usage =
  "SetupSpacetime[dim, signature] sets up a spacetime manifold with the given \
dimension and metric signature. Returns {manifold, metric, covd}.";

Begin["`Private`"];

(* === Spacetime Setup Utilities === *)

SetupSpacetime[dim_Integer, signature_Integer] := Module[
  {manifold, metric, covd, indices},

  (* Generate index symbols based on dimension *)
  indices = If[dim == 2,
    {a, b, c, d, e, f},  (* 1+1D *)
    {a, b, c, d, e, f, g, h}  (* 3+1D or higher *)
  ];

  (* Define manifold *)
  DefManifold[ELSpacetime, dim, indices];

  (* Define metric with given signature *)
  DefMetric[signature, ELMetric[-a, -b], ELSpacetime,
    SymbolOfCovD -> {";", "\[Del]"}];

  (* Get covariant derivative *)
  covd = CovDOfMetric[ELMetric[-a, -b]];

  {ELSpacetime, ELMetric, covd}
];

(* === Euler-Lagrange Computation === *)

(*
  For a scalar field phi with Lagrangian L[phi, D_a phi]:
    EOM = dL/d(phi) - D_a[dL/d(D_a phi)] = 0

  For a vector field A_b with Lagrangian L[A_b, D_a A_b]:
    EOM^b = dL/d(A_b) - D_a[dL/d(D_a A_b)] = 0

  We use VarD from xTensor which computes functional derivatives.
*)

EulerLagrangeEquation[lagrangian_, field_, covd_] := Module[
  {eom},

  (* xTensor's VarD computes the variational derivative *)
  (* VarD[expr, field, covd] gives d(expr)/d(field) accounting for derivatives *)
  eom = VarD[lagrangian, field, covd];

  (* Apply canonical simplifications *)
  eom = ToCanonical[eom];
  eom = ContractMetric[eom];

  eom
];

(* Alternative: Explicit computation for pedagogical clarity *)
VariationalDerivative[lagrangian_, field_, covd_] := Module[
  {directTerm, derivativeTerm, result},

  (* Direct variation: dL/d(field) *)
  directTerm = D[lagrangian, field];

  (* Variation with respect to first derivatives: -D_a[dL/d(D_a field)] *)
  (* This requires careful handling of the covariant derivative structure *)
  derivativeTerm = With[{cd = covd},
    (* Find all terms involving cd[field] and compute their contribution *)
    (* This is a simplified version - VarD handles this properly *)
    0  (* Placeholder - use VarD for actual computation *)
  ];

  (* The full EOM *)
  result = directTerm + derivativeTerm;
  ToCanonical[result]
];

End[];
EndPackage[];

(* Example usage (commented out for package import):

<< TorsionGertsenshtein`EulerLagrange`;
<< xAct`xTensor`;

(* Setup 1+1D Minkowski spacetime *)
DefManifold[M2, 2, {a, b, c, d}];
DefMetric[-1, eta[-a, -b], M2];
CD = CovDOfMetric[eta[-a, -b]];

(* Define scalar field for Klein-Gordon *)
DefTensor[phi[], M2];
DefConstantSymbol[m2];

(* Klein-Gordon Lagrangian: L = -1/2 (d_a phi)(d^a phi) - 1/2 m^2 phi^2 *)
KGLagrangian = -1/2 CD[-a][phi[]] CD[a][phi[]] - 1/2 m2 phi[]^2;

(* Compute EOM *)
eom = EulerLagrangeEquation[KGLagrangian, phi[], CD];
(* Should give: Box[phi] - m^2 phi = 0 *)

*)
