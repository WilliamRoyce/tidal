(* ::Package:: *)
(*
   MODULE: EulerLagrange.wl
   PURPOSE: Derive equations of motion from Lagrangian densities via Euler-Lagrange equations

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus, variational derivatives)

   DATA FLOW:
     Lagrangian density L[φ, ∂φ, ∂²φ, ...]
       → EulerLagrangeEquation (uses xAct VarD for functional derivative δL/δφ)
       → Equation of motion (tensor form)

   KEY FEATURES:
     - Handles scalar and vector fields
     - Supports first and second-order derivative terms in Lagrangian
     - Works with abstract indices (converts to components via ComponentDecompose)
     - Automatic simplification via xAct's NoScalar and SortCovDs

   USAGE PATTERN:
     {manifold, metric, cd} = SetupSpacetime[dim, signature];
     DefTensor[phi[], manifold];  (* Define scalar field *)
     L = -1/2 cd[-a][phi[]] cd[a][phi[]];  (* Lagrangian *)
     eom = EulerLagrangeEquation[L, phi, cd];  (* → ∂_a∂^a φ = 0 *)

   IMPLICIT ASSUMPTIONS:
     - Lagrangian is a scalar density (√-g factor assumed or handled separately)
     - Fields are dynamical (not background)
     - Boundary terms vanish (integration by parts)

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`EulerLagrange`", {"xAct`xTensor`"}];

(* Public symbols *)
EulerLagrangeEquation::usage =
  "EulerLagrangeEquation[lagrangian, field, covd] computes the Euler-Lagrange \
equation for the given field from the Lagrangian density. Returns an expression \
equal to zero (the equation of motion).";

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

  (* Use xTensor's VarD for the variational derivative *)
  (* VarD uses curried syntax: VarD[field, covd][expr] *)
  eom = VarD[field, covd][lagrangian];

  (* Fail if VarD returned unevaluated — do NOT silently fall back *)
  If[Head[eom] === VarD,
    Throw[StringJoin[
      "EulerLagrangeEquation: VarD returned unevaluated for field '",
      ToString[field], "'. ",
      "Ensure xAct is loaded and the field is properly defined on the manifold. ",
      "Check that the field has the correct rank and the covariant derivative matches the metric."
    ]]
  ];

  (* VarD succeeded - apply canonical simplifications *)
  eom = ToCanonical[eom];
  eom = ContractMetric[eom];

  eom
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
