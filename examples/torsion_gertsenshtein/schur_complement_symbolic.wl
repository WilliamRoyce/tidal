(* ============================================================ *)
(* Symbolic Schur Complement Analysis of Torsion-Mediated       *)
(* Gertsenshtein Amplification                                  *)
(*                                                              *)
(* Reads the derived JSON spec and computes closed-form         *)
(* expressions for the effective coupling mu_eff and mass       *)
(* m2_eff in the h5<->a1 Gertsenshtein channel after            *)
(* eliminating torsion constraint fields via the Fourier-space  *)
(* Schur complement.                                            *)
(*                                                              *)
(* Two independent constraint blocks mediate the coupling:      *)
(*   Block A (3x3): {t_0, t_15, t_22}                          *)
(*   Block B (5x5): {t_3, t_4, t_7, t_8, t_18}                *)
(*                                                              *)
(* Usage:                                                       *)
(*   wolframscript -file schur_complement_symbolic.wl            *)
(* ============================================================ *)

(* --- Configuration --- *)
jsonPath = FileNameJoin[{
  DirectoryName[$InputFileName],
  "..", "data", "torsion_gertsenshtein_nonminimal.json"
}];

Print["Loading JSON spec from: ", jsonPath];

(* --- 1. Import and Parse JSON --- *)
json = Import[jsonPath, "RawJSON"];

equations = json["equations"];

(* Build lookup: fieldName -> equation data *)
eqMap = Association @@ (
  (#["field"] -> #) & /@ equations
);

(* --- 2. Define Symbols --- *)
(* Parameters (already used in coefficient_symbolic strings) *)
(* kappa, B0, alpha1, alpha2, alpha3, delta1 are implicit *)
(* Fourier wavenumber *)
k;
(* Eigenvalue parameter (for time derivatives) *)
\[Lambda];

(* --- 3. Operator Fourier Multipliers --- *)
opMultiplier["identity"] = 1;
opMultiplier["laplacian_x"] = -k^2;
opMultiplier["gradient_x"] = I k;
opMultiplier["first_derivative_t"] = \[Lambda];
opMultiplier["d2_t"] = \[Lambda]^2;

(* --- 4. Parse equation terms into symbolic matrix entries --- *)
(* Each term contributes: coefficient * multiplier -> A[eqField, termField] *)

parseCoefficient[coeff_String] := ToExpression[coeff, InputForm];
parseCoefficient[Null] := 1;
parseCoefficient[None] := 1;
parseCoefficient[_Missing] := 1;

(* Extract matrix contribution from a single term *)
parseTerm[term_Association] := Module[{coeffSym, coeff, op, field, mult},
  coeffSym = Lookup[term, "coefficient_symbolic", Missing["KeyAbsent"]];
  coeff = parseCoefficient[coeffSym];
  op = term["operator"];
  field = term["field"];
  mult = opMultiplier[op];
  If[Head[mult] === opMultiplier,
    Print["WARNING: unknown operator ", op];
    mult = 1;
  ];
  {field, coeff * mult}
];

(* Build symbolic coupling map for an equation:
   eqField -> Association[targetField -> total symbolic coefficient] *)
buildEquationCouplings[eqField_String] := Module[{eq, terms, result},
  eq = eqMap[eqField];
  terms = eq["rhs"]["terms"];
  result = Association[];
  Do[
    Module[{parsed, field, contrib},
      parsed = parseTerm[term];
      field = parsed[[1]];
      contrib = parsed[[2]];
      If[KeyExistsQ[result, field],
        result[field] += contrib,
        result[field] = contrib
      ];
    ],
    {term, terms}
  ];
  result
];

Print["\n=== PARSING EQUATIONS ==="];

(* --- 5. Build Block A: {t_0, t_15, t_22} --- *)
(* These couple to h_5 and a_1 via spatial operators *)
Print["\n--- Block A: {t_0, t_15, t_22} ---"];

blockAFields = {"t_0", "t_15", "t_22"};

(* S_cc for Block A: 3x3 constraint self-coupling *)
SccA = Table[
  Module[{couplings, field},
    couplings = buildEquationCouplings[blockAFields[[i]]];
    field = blockAFields[[j]];
    If[KeyExistsQ[couplings, field],
      couplings[field],
      0
    ]
  ],
  {i, 3}, {j, 3}
];

Print["S_cc (Block A):"];
Print[MatrixForm[SccA // Simplify]];

(* S_cd for Block A: constraint <- dynamical coupling *)
(* Sources from h_5 and a_1 *)
ScdA = Table[
  Module[{couplings, field},
    couplings = buildEquationCouplings[blockAFields[[i]]];
    field = {"h_5", "a_1"}[[j]];
    If[KeyExistsQ[couplings, field],
      couplings[field],
      0
    ]
  ],
  {i, 3}, {j, 2}
];

Print["\nS_cd (Block A, sources from {h_5, a_1}):"];
Print[MatrixForm[ScdA // Simplify]];

(* A_dc for Block A: dynamical <- constraint coupling *)
(* How t_0, t_15, t_22 appear in h_5 and a_1 equations *)
(* For h_5 (2nd order): torsion terms appear in the v_h5 equation *)
(* For a_1 (2nd order): torsion terms appear in the v_a1 equation *)

h5Couplings = buildEquationCouplings["h_5"];
a1Couplings = buildEquationCouplings["a_1"];

(* h_5 equation coupling TO torsion (these go in the v_h5 row) *)
AdcA = Table[
  Module[{couplings, field},
    couplings = {h5Couplings, a1Couplings}[[i]];
    field = blockAFields[[j]];
    If[KeyExistsQ[couplings, field],
      couplings[field],
      0
    ]
  ],
  {i, 2}, {j, 3}
];

(* CRITICAL: Normalize h_5 row by kinetic coefficient.
   The h_5 equation in JSON has kinetic_coefficient_symbolic = -1/kappa^2,
   meaning the actual equation is (-1/kappa^2) d2_t h5 = RHS.
   To get d2_t h5 = normalized_RHS, divide RHS by -1/kappa^2 = multiply by -kappa^2.
   Row 1 of AdcA is the h_5 row; row 2 is a_1 (no kinetic coeff, stays unchanged). *)
kineticH5 = -1/kappa^2;
AdcA[[1]] = AdcA[[1]] / kineticH5;

Print["\nA_dc (Block A, {h_5,a_1} <- {t_0,t_15,t_22}), kinetic-normalized:"];
Print[MatrixForm[AdcA // Simplify]];

(* --- 6. Compute Block A Schur Complement --- *)
Print["\n=== BLOCK A SCHUR COMPLEMENT ==="];

SccAInv = Inverse[SccA] // Simplify;
Print["\nInverse(S_cc_A):"];
Print[MatrixForm[SccAInv]];

(* Schur complement correction to {h_5, a_1} block *)
(* correction[i,j] = sum_c1,c2 A_dc[i,c1] * SccInv[c1,c2] * S_cd[c2,j] *)
correctionA = AdcA . SccAInv . ScdA // Simplify;

Print["\nSchur correction from Block A (2x2, rows={h_5,a_1}, cols={h_5,a_1}):"];
Print[MatrixForm[correctionA]];

(* --- 7. Build Block B: {t_3, t_4, t_7, t_8, t_18} --- *)
Print["\n--- Block B: {t_3, t_4, t_7, t_8, t_18} ---"];

blockBFields = {"t_3", "t_4", "t_7", "t_8", "t_18"};

(* S_cc for Block B: 5x5 *)
SccB = Table[
  Module[{couplings, field},
    couplings = buildEquationCouplings[blockBFields[[i]]];
    field = blockBFields[[j]];
    If[KeyExistsQ[couplings, field],
      couplings[field],
      0
    ]
  ],
  {i, 5}, {j, 5}
];

Print["S_cc (Block B):"];
Print[MatrixForm[SccB // Simplify]];

(* Block B sources: h_5 via first_derivative_t, v_a_1 via gradient_x *)
(* Note: v_a_1 is a velocity, not a field position. In the first-order
   system representation, this appears differently. For now, extract
   whatever dynamical fields appear as sources. *)

(* Sources from h_5 and v_a_1 *)
ScdB = Table[
  Module[{couplings, field},
    couplings = buildEquationCouplings[blockBFields[[i]]];
    field = {"h_5", "v_a_1"}[[j]];
    If[KeyExistsQ[couplings, field],
      couplings[field],
      0
    ]
  ],
  {i, 5}, {j, 2}
];

Print["\nS_cd (Block B, sources from {h_5, v_a_1}):"];
Print[MatrixForm[ScdB // Simplify]];

(* Block B: how these torsion fields appear in h_5 and a_1 equations *)
(* h_5 references t_4, t_8 via first_derivative_t *)
(* a_1 references v_t_18, v_t_4 via gradient_x *)

(* For the first-order system, torsion velocities v_t appear in a_1 eq.
   The velocity of constraint field c is: v_c = recovery @ d' = recovery @ A @ d
   This creates an implicit coupling that the modal solver handles via
   the (I - vel_coupling)^-1 factor. *)

(* Extract coupling from h_5 to t_4, t_8 (via first_derivative_t = lambda) *)
(* and from a_1 to v_t_4, v_t_18 (via gradient_x) *)

(* h5 equation: coupling to Block B torsion fields *)
adcBfromH5 = Table[
  Module[{field},
    field = blockBFields[[j]];
    If[KeyExistsQ[h5Couplings, field],
      h5Couplings[field],
      0
    ]
  ],
  {j, 5}
];

(* Normalize h5 coupling by kinetic coefficient *)
adcBfromH5 = adcBfromH5 / kineticH5;

(* a1 equation: coupling to Block B torsion VELOCITIES *)
(* The a_1 equation references v_t_4 and v_t_18 *)
velCoupA1 = Association[];
Do[
  Module[{parsed},
    parsed = parseTerm[term];
    If[StringMatchQ[parsed[[1]], "v_t_*"],
      velCoupA1[parsed[[1]]] = parsed[[2]];
    ];
  ],
  {term, eqMap["a_1"]["rhs"]["terms"]}
];

(* Map velocity names to constraint field indices *)
adcBvelFromA1 = Table[
  Module[{velName},
    velName = "v_" <> blockBFields[[j]];
    If[KeyExistsQ[velCoupA1, velName],
      velCoupA1[velName],
      0
    ]
  ],
  {j, 5}
];

Print["\nA_dc (h_5 <- Block B torsion fields):"];
Print[adcBfromH5 // Simplify];
Print["\nA_dc (a_1 <- Block B torsion VELOCITIES):"];
Print[adcBvelFromA1 // Simplify];

(* --- 8. Compute Block B Schur Complement --- *)
Print["\n=== BLOCK B SCHUR COMPLEMENT ==="];

SccBInv = Inverse[SccB] // Simplify;
Print["\nInverse(S_cc_B):"];
Print[MatrixForm[SccBInv]];

(* Block B constraint recovery: t_B = -SccBInv . ScdB . {h_5, v_a_1} *)
recoveryB = -SccBInv . ScdB // Simplify;
Print["\nRecovery matrix (5x2, rows=Block B fields, cols={h_5, v_a_1}):"];
Print[MatrixForm[recoveryB]];

(* Block B correction to h_5 equation (via field coupling):
   h_5 gets terms from t_4, t_8 (via first_derivative_t=lambda and gradient_x)
   These are FIELD-POSITION couplings, so correction = adcBfromH5 . (-SccBInv) . ScdB *)
corrBh5 = adcBfromH5 . (-SccBInv) . ScdB // Simplify;

Print["\nBlock B field correction to h_5 eq (row, cols={h_5, v_a_1}):"];
Print[corrBh5];

(* Block B correction to a_1 equation (via VELOCITY coupling):
   a_1 gets terms from v_t_4, v_t_18 (via gradient_x)
   Constraint velocity: v_c = d/dt[c] = d/dt[-SccBInv . ScdB . d]
   At eigenvalue lambda: v_c = lambda * recovery . d
   So the velocity contribution is: adcBvelFromA1 . (lambda * recovery) . d
   = lambda * adcBvelFromA1 . (-SccBInv . ScdB) . d *)
corrBa1vel = \[Lambda] * adcBvelFromA1 . (-SccBInv) . ScdB // Simplify;

Print["\nBlock B velocity correction to a_1 eq (row, cols={h_5, v_a_1}):"];
Print[corrBa1vel];

(* Combined Block B correction to the 2x2 {h_5, a_1} block:
   Row 1 (h_5 / v_h_5 equation): gets field correction from Block B
   Row 2 (a_1 / v_a_1 equation): gets velocity correction from Block B
   But note: ScdB has cols {h_5, v_a_1}, not {h_5, a_1}.
   The v_a_1 column couples to the velocity slot, not the field slot.
   This creates a non-standard matrix structure.

   For the 4x4 block [h_5, v_h_5, a_1, v_a_1]:
   - corrBh5 contributes to the v_h_5 row (acceleration equation)
   - corrBa1vel contributes to the v_a_1 row (acceleration equation)
   But the column indices are {h_5, v_a_1} not {h_5, a_1}.

   So Block B creates coupling: v_h_5 <- h_5 and v_h_5 <- v_a_1
   And: v_a_1 <- h_5 and v_a_1 <- v_a_1 (via lambda)

   The v_a_1 <- v_a_1 term is an IMPLICIT coupling (modifies the LHS). *)

correctionB = {{corrBh5[[1]], corrBh5[[2]]}, {corrBa1vel[[1]], corrBa1vel[[2]]}};
Print["\nBlock B combined correction (2x2, rows={v_h_5, v_a_1}, cols={h_5, v_a_1}):"];
Print[MatrixForm[correctionB // Simplify]];

(* --- 9. Combined Effective Coupling --- *)
Print["\n=== COMBINED EFFECTIVE COUPLING ==="];

(* Block A correction (rows={v_h5, v_a1}, cols={h_5, a_1}):
   These modify the acceleration equations for h_5 and a_1 *)
Print["\n--- From Block A (cols={h_5, a_1}) ---"];
Print["v_h5 <- h_5 (mass shift):  ", correctionA[[1,1]] // Simplify];
Print["v_h5 <- a_1 (coupling):    ", correctionA[[1,2]] // Simplify];
Print["v_a1 <- h_5 (coupling):    ", correctionA[[2,1]] // Simplify];
Print["v_a1 <- a_1 (mass shift):  ", correctionA[[2,2]] // Simplify];

(* Block B correction (rows={v_h5, v_a1}, cols={h_5, v_a_1}):
   NOTE: Block B couples to v_a_1, not a_1! This creates implicit LHS. *)
Print["\n--- From Block B (cols={h_5, v_a_1}) ---"];
Print["v_h5 <- h_5:    ", corrBh5[[1]] // Simplify];
Print["v_h5 <- v_a_1:  ", corrBh5[[2]] // Simplify];
Print["v_a1 <- h_5:    ", corrBa1vel[[1]] // Simplify];
Print["v_a1 <- v_a_1:  ", corrBa1vel[[2]] // Simplify];

(* The baseline Gertsenshtein h_5<->a_1 coupling:
   BEFORE normalization: v_h5 eq has B0 * I*k * a_1
   AFTER normalization: v_h5 eq has (-kappa^2) * B0 * I*k * a_1 = -B0*kappa^2 * I*k * a_1
   For v_a1: no kinetic coeff on a_1, so stays B0 * I*k * h_5 *)
Print["\n--- Combined effective coupling ---"];
muGRh5 = -kappa^2 * B0 * I * k; (* h5<-a1 coupling AFTER normalization *)
muGRa1 = B0 * I * k;             (* a1<-h5 coupling (no normalization needed) *)
Print["mu_GR(h5<-a1) = ", muGRh5, " (after kinetic normalization)"];
Print["mu_GR(a1<-h5) = ", muGRa1];

(* Total coupling shift h5<-a1: Block A (cols h_5,a_1) + Block B (no a_1 col) *)
dmuH5fromA1 = correctionA[[1,2]];
Print["Delta_mu(h5<-a1) from Block A = ", dmuH5fromA1 // Simplify];
(* Block B contributes to h5 only via h_5 column, not a_1 *)

(* Total a1 mass shift: Block A + Block B implicit *)
dm2A1 = correctionA[[2,2]];
dm2A1implB = corrBa1vel[[2]]; (* implicit: modifies LHS, lambda-dependent *)
Print["\nDelta_m2(a1) from Block A = ", dm2A1 // Simplify];
Print["Delta_implicit(v_a1<-v_a1) from Block B = ", dm2A1implB // Simplify];
Print["  (This is a LHS modification: (1 - D)*v_a1' = RHS)"];

(* Total h5 mass shift: Block A + Block B *)
dm2H5 = correctionA[[1,1]] + corrBh5[[1]];
Print["\nDelta_m2(h5) from Block A + B = ", dm2H5 // Simplify];

(* Total coupling h5->a1 shift (in v_a1 equation): Block A + Block B *)
dmuA1fromH5 = correctionA[[2,1]] + corrBa1vel[[1]];
Print["Delta_mu(a1<-h5) from Block A + B = ", dmuA1fromH5 // Simplify];

(* Effective coupling = GR baseline + Schur correction.
   The Schur complement SUBTRACTS the correction from the system matrix: A_eff = A - A_dc.Scc^-1.S_cd
   So the effective coupling element becomes: mu_GR - correction *)
muEffH5 = muGRh5 - dmuH5fromA1;
muEffA1 = muGRa1 - dmuA1fromH5;
Print["\nmu_eff(h5<-a1) = ", muEffH5 // Simplify];
Print["mu_eff(a1<-h5) = ", muEffA1 // Simplify];
Print["mu_ratio(h5<-a1) = mu_eff/mu_GR = ", Simplify[muEffH5/muGRh5]];
Print["mu_ratio(a1<-h5) = mu_eff/mu_GR = ", Simplify[muEffA1/muGRa1]];

(* Effective masses.
   After kinetic normalization (dividing h5 RHS by -1/kappa^2 = multiplying by -kappa^2):
   h5 baseline: d2_t h5 = 1*laplacian(h5) + (-kappa^2*B0^2/2)*h5 + coupling*a1
   In Fourier: m2(h5)_baseline = -k^2 - kappa^2*B0^2/2  (NEGATIVE = oscillatory)
   The Schur correction dm2H5 has already been normalized (AdcA row 1 was divided by kineticH5). *)
m2H5eff = -k^2 - kappa^2*B0^2/2 - dm2H5;
m2A1eff = -k^2 - dm2A1;
Print["\nm2_eff(h5) = -k^2 - kappa^2*B0^2/2 - Delta = ", m2H5eff // Simplify];
Print["m2_eff(a1) = -k^2 - Delta_A = ", m2A1eff // Simplify];
Print["Implicit LHS factor for a1: 1 - ", dm2A1implB // Simplify];

(* --- 10. TeX Export --- *)
Print["\n=== TEX EXPORT ==="];

Print["\n% Block A mass matrix S_cc:"];
Print["S_{cc,A} = ", TeXForm[SccA // Simplify]];

Print["\n% Block A inverse:"];
Print["S_{cc,A}^{-1} = ", TeXForm[SccAInv // Simplify]];

Print["\n% Block A correction matrix (rows={v_h5,v_a1}, cols={h_5,a_1}):"];
Print["\\Delta_A = ", TeXForm[correctionA // Simplify]];

Print["\n% Block B mass matrix S_cc:"];
Print["S_{cc,B} = ", TeXForm[SccB // Simplify]];

Print["\n% Block B field correction to v_h5 (cols={h_5, v_a_1}):"];
Print["\\Delta_{B,h5} = ", TeXForm[corrBh5 // Simplify]];

Print["\n% Block B velocity correction to v_a1 (cols={h_5, v_a_1}):"];
Print["\\Delta_{B,a1} = ", TeXForm[corrBa1vel // Simplify]];

Print["\n% Effective coupling h5<-a1:"];
Print["\\mu_{\\text{eff}}^{h_5 \\leftarrow a_1} = ", TeXForm[muEffH5 // Simplify]];

Print["\n% Effective coupling a1<-h5:"];
Print["\\mu_{\\text{eff}}^{a_1 \\leftarrow h_5} = ", TeXForm[muEffA1 // Simplify]];

Print["\n% Effective h5 mass^2:"];
Print["m^2_{\\text{eff}}(h_5) = ", TeXForm[m2H5eff // Simplify]];

Print["\n% Effective a1 mass^2:"];
Print["m^2_{\\text{eff}}(a_1) = ", TeXForm[m2A1eff // Simplify]];

(* --- 11. Numerical Verification --- *)
Print["\n=== NUMERICAL VERIFICATION ==="];
Print["Evaluating at alpha1=0, alpha2=-0.6, alpha3=1.0, kappa=1, B0=0.0001, delta1=0.5, k=0.503:"];

testRules = {
  alpha1 -> 0, alpha2 -> -0.6, alpha3 -> 1.0,
  kappa -> 1, B0 -> 0.0001, delta1 -> 0.5, k -> 0.5027,
  \[Lambda] -> 0.5027  (* eigenvalue ~ k for the dominant mode *)
};

Print["\n--- Block A ---"];
Print["  SccA = ", SccA /. testRules // N];
Print["  correctionA = ", correctionA /. testRules // N];

Print["\n--- Block B ---"];
Print["  SccB = ", SccB /. testRules // N];
Print["  corrBh5 = ", corrBh5 /. testRules // N];
Print["  corrBa1vel = ", corrBa1vel /. testRules // N];

Print["\n--- Combined ---"];
Print["  mu_eff(h5<-a1) = ", muEffH5 /. testRules // N];
Print["  mu_eff(a1<-h5) = ", muEffA1 /. testRules // N];
Print["  m2_eff(h5) = ", m2H5eff /. testRules // N];
Print["  m2_eff(a1) = ", m2A1eff /. testRules // N];

(* Expected from Python Schur complement at (d1=0.5, a2=-0.6, k=0.503):
   coupling shift (Im) ~ -5.3e-5  (this includes BOTH blocks)
   mass shift ~ +0.81 *)
Print["\n  Expected from Python:"];
Print["    coupling shift (Im of v_h5<-a1) ~ -5.3e-5"];
Print["    a1 mass shift ~ +0.81"];

(* Also test at amplification peak *)
Print["\n--- Amplification peak: alpha2=-0.82, delta1=1.0 ---"];
peakRules = {
  alpha1 -> 0, alpha2 -> -0.82, alpha3 -> 1.0,
  kappa -> 1, B0 -> 0.0001, delta1 -> 1.0, k -> 0.5027,
  \[Lambda] -> 0.5027
};
Print["  mu_eff(h5<-a1) = ", muEffH5 /. peakRules // N];
Print["  mu_eff(a1<-h5) = ", muEffA1 /. peakRules // N];
Print["  Python: mu_eff/mu_GR ~ -1.92"];

Print["\nDone."];

Quit[];
