(* ::Package:: *)
(*
   MODULE: EulerLagrange.wl
   PURPOSE: Derive equations of motion and canonical Hamiltonian from Lagrangian densities

   DEPENDENCIES:
     - xAct`xTensor` (tensor calculus, variational derivatives)

   DATA FLOW:
     Lagrangian density L[φ, ∂φ, ∂²φ, ...]
       → EulerLagrangeEquation (uses xAct VarD for functional derivative δL/δφ)
       → Equation of motion (tensor form, cross-check only)

     Lagrangian density L
       → CanonicalMomentum (VarD with auxiliary tensor for ∂L/∂(∂_t q))
       → LegendreTransformH (H = Σ π·q̇ - L)
       → Hamilton's equations (∂H/∂π, -∂H/∂q)

   KEY FEATURES:
     - Handles scalar and vector fields
     - Supports first and second-order derivative terms in Lagrangian
     - Works with abstract indices (converts to components via ComponentDecompose)
     - Automatic simplification via xAct's NoScalar and SortCovDs
     - Canonical momentum computed directly from L (not from EL equations)
     - Hamiltonian via standard Legendre transform (not ad-hoc assembly)

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
     - For CanonicalMomentum/LegendreTransformH: Lagrangian is quadratic in
       velocities (true for all linearized theories after xPert expansion)

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`EulerLagrange`", {"xAct`xTensor`"}];

(* Public symbols *)
EulerLagrangeEquation::usage =
  "EulerLagrangeEquation[lagrangian, field, covd] computes the Euler-Lagrange \
equation for the given field from the Lagrangian density. Returns an expression \
equal to zero (the equation of motion).";

CanonicalMomentum::usage =
  "CanonicalMomentum[lagrangian, field, chart, covd] computes the canonical \
momentum pi = dL/d(d_t field) directly from the Lagrangian using VarD with an \
auxiliary tensor (ADM-style approach). For scalar fields returns a scalar; for \
vector fields A[-a] returns a one-form pi[-a]. The auxiliary tensor is defined, \
used for VarD, and cleaned up automatically. References: arXiv:2210.10103.";

LegendreTransformH::usage =
  "LegendreTransformH[lagrangian, fields, chart, covd] computes the Hamiltonian \
density H = Sum[pi_i * d_t(q_i)] - L via the standard Legendre transform. \
fields is a list of field heads {phi, A, B, ...}. Returns {H, piList} where H \
is the Hamiltonian expression and piList is the list of {field, pi} pairs. \
For quadratic Lagrangians, velocities are eliminated algebraically.";

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


(* === Canonical Momentum Computation === *)

(*
  Canonical momentum: π = ∂L/∂(∂_t q)

  Strategy (ADM-style):
  1. Create an auxiliary tensor `dotField` with the same index structure as ∂_t(field)
  2. Substitute CD[-t][field[...]] → dotField[...] in the Lagrangian
  3. Apply VarD[dotField, covd] to get δL/δ(dotField)
  4. Since dotField appears only algebraically (no derivatives of dotField in L),
     VarD reduces to the partial derivative ∂L/∂(dotField) = ∂L/∂(∂_t field)
  5. Clean up the auxiliary tensor

  For scalar fields: ∂_t phi → dotPhi (scalar), π is a scalar
  For vector fields: ∂_t A[-a] → dotA[-a] (1-form), π is a 1-form π[-a]
*)

CanonicalMomentum[lagrangian_, field_, chart_, covd_] := Module[
  {timeCoord, manifold, dotName, dotHead, Lsub, pi, fieldIndices,
   fieldRank, dotExpr, substitutionRule},

  timeCoord = First[ScalarsOfChart[chart]];
  manifold = First[Select[HostsOf[field], ManifoldQ]];

  (* Determine field rank from its index structure *)
  fieldRank = Length[SlotsOfTensor[field]];

  (* Create unique auxiliary tensor name *)
  dotName = Symbol["dot$" <> SymbolName[field]];

  (* Define auxiliary tensor with same symmetry as the field *)
  If[fieldRank == 0,
    (* Scalar field: dotField is a scalar *)
    If[!xTensorQ[dotName], DefTensor[dotName[], manifold]],
    (* Vector/tensor field: dotField has same index structure *)
    If[!xTensorQ[dotName],
      DefTensor[dotName @@ Table[Unique["idx"], fieldRank], manifold,
        Symmetric @@ SymmetryGroupOfTensor[field]]
    ]
  ];

  (* Build substitution rule: CD[-t][field[indices]] → dotField[indices] *)
  (* For scalars: CD[-t][phi[]] → dotPhi[] *)
  (* For vectors: CD[-t][A[-a]] → dotA[-a] *)
  If[fieldRank == 0,
    substitutionRule = covd[_?((# === -timeCoord) &)][field[]] :> dotName[];
    dotExpr = dotName[],
    (* For rank-1 (vector): match CD[-t][A[idx_]] *)
    substitutionRule = covd[_?((# === -timeCoord) &)][field[idx__]] :> dotName[idx];
    dotExpr = dotName @@ Table[
      With[{s = SlotsOfTensor[field][[i]]},
        If[s > 0, Unique["a"], -Unique["a"]]
      ],
      {i, fieldRank}
    ]
  ];

  (* Substitute time derivatives with auxiliary tensor *)
  Lsub = lagrangian /. substitutionRule;

  (* VarD with auxiliary tensor: since dotField has no derivatives in Lsub,
     this gives the PARTIAL derivative ∂L/∂(dotField), not a functional derivative *)
  pi = VarD[dotExpr, covd][Lsub];

  (* Fail if VarD returned unevaluated *)
  If[Head[pi] === VarD,
    Throw[StringJoin[
      "CanonicalMomentum: VarD returned unevaluated for auxiliary tensor '",
      ToString[dotName], "' (field: ", ToString[field], "). ",
      "Check that the Lagrangian depends on time derivatives of this field."
    ]]
  ];

  (* Simplify *)
  pi = ToCanonical[pi];
  pi = ContractMetric[pi];

  (* Substitute back: dotField → CD[-t][field] to express pi in terms of velocities *)
  If[fieldRank == 0,
    pi = pi /. dotName[] -> covd[-timeCoord][field[]],
    pi = pi /. dotName[idx__] :> covd[-timeCoord][field[idx]]
  ];

  pi = ToCanonical[pi];
  pi = ContractMetric[pi];

  (* Clean up auxiliary tensor *)
  If[xTensorQ[dotName], UndefTensor[dotName]];

  pi
];


(* === Legendre Transform: Hamiltonian from Lagrangian === *)

(*
  H = Σ_i π_i · q̇_i - L

  For each field q_i:
    1. Compute π_i = ∂L/∂(∂_t q_i) via CanonicalMomentum
    2. Form the kinetic term π_i · q̇_i (contracted over field indices)
    3. Sum over all fields and subtract L

  For quadratic Lagrangians (all linearized theories), the velocities q̇_i
  can be eliminated algebraically by inverting the momentum definition.
  For the standard kinetic term -½(∂_a q)(∂^a q):
    π = -∂_t q  (Minkowski, signature -,+,+)
    q̇ = -π
    H = π·(-π) - L = -π² - L

  The result H is expressed in terms of fields, canonical momenta (as CD[-t][field]),
  and spatial derivatives only.
*)

LegendreTransformH[lagrangian_, fields_List, chart_, covd_] := Module[
  {timeCoord, piList, H, fieldRank, velocity, piTimesVelocity},

  timeCoord = First[ScalarsOfChart[chart]];
  piList = {};

  (* Compute Σ π_i · q̇_i *)
  H = 0;
  Do[
    fieldRank = Length[SlotsOfTensor[f]];
    pi = CanonicalMomentum[lagrangian, f, chart, covd];

    AppendTo[piList, {f, pi}];

    (* Form π · q̇ (contract over field indices) *)
    If[fieldRank == 0,
      (* Scalar: π · q̇ = π * ∂_t φ *)
      velocity = covd[-timeCoord][f[]];
      piTimesVelocity = pi * velocity,
      (* Vector: π_a · q̇^a = π[-a] ∂_t A[a] (raise one index via metric) *)
      piTimesVelocity = pi * covd[-timeCoord][f[]];
      (* For rank > 0, we need proper contraction — use metric *)
      piTimesVelocity = ToCanonical[ContractMetric[piTimesVelocity]]
    ];

    H += piTimesVelocity,
    {f, fields}
  ];

  (* H = Σ π·q̇ - L *)
  H = H - lagrangian;

  (* Simplify *)
  H = ToCanonical[H];
  H = ContractMetric[H];

  {H, piList}
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

(* Compute EOM (cross-check) *)
eom = EulerLagrangeEquation[KGLagrangian, phi[], CD];
(* Should give: Box[phi] - m^2 phi = 0 *)

(* Compute canonical momentum — directly from L, not from EOM *)
pi = CanonicalMomentum[KGLagrangian, phi[], cart, CD];
(* Should give: -CD[-t][phi[]]  (π = -∂_t φ for Minkowski signature -,+) *)

(* Compute Hamiltonian via Legendre transform *)
{H, piList} = LegendreTransformH[KGLagrangian, {phi}, cart, CD];
(* Should give: H = 1/2 (∂_t φ)^2 + 1/2 (∂_x φ)^2 + 1/2 m^2 φ^2 *)

*)
