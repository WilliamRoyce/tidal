(* ::Package:: *)
(*
   MODULE: EulerLagrange.wl
   PURPOSE: Derive equations of motion from Lagrangian densities

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus, variational derivatives)

   DATA FLOW:
     Lagrangian density L[φ, ∂φ, ∂²φ, ...]
       → EulerLagrangeEquation (uses xAct VarD for functional derivative δL/δφ)
       → Equation of motion (tensor form)

   NOTE: Canonical momentum and Hamiltonian computation is now performed at
   the component level in _derive.py (_wls_canonical_phase_a / _wls_canonical_phase_b). The old
   abstract-index CanonicalMomentum/LegendreTransformH functions have been
   removed — they were superseded by the component-level Legendre transform
   which correctly handles gauge theories (Proca, Chern-Simons, etc.).

   KEY FEATURES:
     - Handles scalar and vector fields
     - Supports first and second-order derivative terms in Lagrangian
     - Works with abstract indices (converts to components via ComponentDecompose)
     - Automatic simplification via xAct's NoScalar and SortCovDs

   REFERENCES:
     - ADM formalism with xAct: arXiv:2210.10103
     - xCPS package: github.com/juanmargalef/xCPS
     - VarD with auxiliary tensors: xAct examples (Lagrangian-variation-xPert-VarD.nb)

   USAGE PATTERN:
     DefManifold[M2, 2, {a, b, c, d}];
     DefMetric[-1, eta[-a, -b], CD];
     DefTensor[phi[], M2];  (* Define scalar field *)
     L = -1/2 CD[-a][phi[]] CD[a][phi[]];  (* Lagrangian *)
     eom = EulerLagrangeEquation[L, phi, CD];  (* → ∂_a∂^a φ = 0 *)

   IMPLICIT ASSUMPTIONS:
     - Lagrangian is a scalar density (√-g factor assumed or handled separately)
     - Fields are dynamical (not background)
     - Boundary terms vanish (integration by parts)

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`EulerLagrange`", {"xAct`xTensor`"}];

(* Public symbols *)
EulerLagrangeEquation::usage =
  "EulerLagrangeEquation[lagrangian, field, covd] computes the Euler-Lagrange \
equation for the given field from the Lagrangian density. Returns an expression \
equal to zero (the equation of motion).";

Begin["`Private`"];

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
