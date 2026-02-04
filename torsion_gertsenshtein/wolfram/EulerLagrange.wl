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

  (* Try xTensor's VarD first *)
  (* VarD uses curried syntax: VarD[field, covd][expr] *)
  eom = VarD[field, covd][lagrangian];

  (* Check if VarD returned unevaluated - if so, use manual computation *)
  If[Head[eom] === VarD,
    (* VarD failed - fall back to manual Euler-Lagrange *)
    Print["VarD returned unevaluated - using manual Euler-Lagrange computation"];
    eom = VariationalDerivative[lagrangian, field, covd],

    (* VarD worked - apply canonical simplifications *)
    eom = ToCanonical[eom];
    eom = ContractMetric[eom]
  ];

  eom
];

(* Manual Euler-Lagrange Computation *)
(* Implements: δL/δfield = ∂L/∂field - ∇_a(∂L/∂(∇_a field)) = 0 *)
(* This provides a fallback if VarD is not working *)

VariationalDerivative[lagrangian_, field_, covd_] := Module[
  {directTerm, derivativeTerm, result, isScalar, derivLag, a, b},

  (* Direct term: ∂L/∂field *)
  directTerm = D[lagrangian, field];

  (* Determine if scalar or vector field *)
  (* Scalar field: phi[] - has no indices *)
  (* Vector field: A[-a] - has one index *)
  (* Check for indices by looking for AbstractIndexQ patterns *)
  isScalar = FreeQ[field, _?AbstractIndexQ];

  (* Derivative term: -∇_a (∂L/∂(∇_a field)) *)
  derivativeTerm = If[isScalar,
    (* Scalar field: φ[] *)
    (* Compute -∇_a (∂L/∂(∇_a φ)) *)
    derivLag = D[lagrangian, covd[-a][field]];
    -covd[-a][derivLag],

    (* Vector or tensor field *)
    (* For vector A[-a]: Compute -∇_b (∂L/∂(∇_b A_a)) *)
    derivLag = D[lagrangian, covd[-b][field]];
    -covd[-b][derivLag]
  ];

  (* Euler-Lagrange equation: δL/δfield = 0 *)
  result = directTerm + derivativeTerm;

  (* Simplify using xTensor *)
  result = ToCanonical[result];
  result = ContractMetric[result];

  result
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
