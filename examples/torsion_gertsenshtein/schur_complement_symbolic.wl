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

Print["\nA_dc (Block A, {h_5,a_1} <- {t_0,t_15,t_22}):"];
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
AdcB_h5 = Table[
  Module[{field},
    field = blockBFields[[j]];
    If[KeyExistsQ[h5Couplings, field],
      h5Couplings[field],
      0
    ]
  ],
  {j, 5}
];

(* a1 equation: coupling to Block B torsion VELOCITIES *)
(* The a_1 equation references v_t_4 and v_t_18 *)
a1VelCouplings = Association[];
Do[
  Module[{parsed},
    parsed = parseTerm[term];
    If[StringMatchQ[parsed[[1]], "v_t_*"],
      a1VelCouplings[parsed[[1]]] = parsed[[2]];
    ];
  ],
  {term, eqMap["a_1"]["rhs"]["terms"]}
];

(* Map velocity names to constraint field indices *)
AdcB_a1_vel = Table[
  Module[{velName},
    velName = "v_" <> blockBFields[[j]];
    If[KeyExistsQ[a1VelCouplings, velName],
      a1VelCouplings[velName],
      0
    ]
  ],
  {j, 5}
];

Print["\nA_dc (h_5 <- Block B torsion fields):"];
Print[AdcB_h5 // Simplify];
Print["\nA_dc (a_1 <- Block B torsion VELOCITIES):"];
Print[AdcB_a1_vel // Simplify];

(* --- 8. Compute Block B Schur Complement --- *)
Print["\n=== BLOCK B SCHUR COMPLEMENT ==="];

SccBInv = Inverse[SccB] // Simplify;

(* Block B correction to h_5 equation (via field coupling) *)
(* h_5 correction from Block B: AdcB_h5 . SccBInv . ScdB *)
correctionB_h5 = AdcB_h5 . SccBInv . ScdB // Simplify;

Print["\nBlock B correction to h_5 (row vector, cols={h_5, v_a_1}):"];
Print[correctionB_h5];

(* Block B correction to a_1 equation (via velocity coupling) *)
(* This is more complex: a_1 couples to torsion VELOCITIES, which require
   the recovery matrix and A_reduced to evaluate. For now, report the
   raw velocity coupling. *)
correctionB_a1_vel = AdcB_a1_vel . SccBInv . ScdB // Simplify;

Print["\nBlock B correction to a_1 via velocities (row, cols={h_5, v_a_1}):"];
Print[correctionB_a1_vel];

(* --- 9. Combined Effective Coupling --- *)
Print["\n=== COMBINED EFFECTIVE COUPLING ==="];

(* Block A contributes to the v_h5 equation (h_5 row of correction)
   and v_a1 equation (a_1 row of correction).

   The 2x2 correction from Block A:
     correctionA[[1,1]] : correction to (v_h5 <- h_5) [h_5 mass shift]
     correctionA[[1,2]] : correction to (v_h5 <- a_1) [h_5<-a_1 coupling shift]
     correctionA[[2,1]] : correction to (v_a1 <- h_5) [a_1<-h_5 coupling shift]
     correctionA[[2,2]] : correction to (v_a1 <- a_1) [a_1 mass shift]
*)

Print["\n--- From Block A ---"];
Print["h_5 mass shift (v_h5 <- h_5):    ", correctionA[[1,1]] // Simplify];
Print["h_5<-a_1 coupling shift:          ", correctionA[[1,2]] // Simplify];
Print["a_1<-h_5 coupling shift:          ", correctionA[[2,1]] // Simplify];
Print["a_1 mass shift (v_a1 <- a_1):    ", correctionA[[2,2]] // Simplify];

(* The baseline Gertsenshtein coupling is: B0 * I * k (from gradient_x) *)
Print["\n--- Effective coupling formula ---"];
muGR = B0 * I * k;
Print["mu_GR (baseline) = ", muGR];

(* The coupling shift from Block A *)
dmuA = correctionA[[1,2]];
Print["Delta_mu from Block A = ", dmuA // Simplify];

(* Effective coupling *)
muEff = muGR - dmuA;  (* Schur complement SUBTRACTS the correction *)
Print["mu_eff = mu_GR - Delta_mu_A = ", muEff // Simplify];

(* a_1 mass shift from Block A *)
dm2A = correctionA[[2,2]];
Print["\nDelta_m2 from Block A = ", dm2A // Simplify];

(* Effective a_1 mass (baseline is -k^2 from laplacian) *)
m2GR = -k^2;
m2Eff = m2GR - dm2A;  (* Schur complement subtracts *)
Print["m2_eff(a_1) = -k^2 - Delta_m2_A = ", m2Eff // Simplify];

(* --- 10. TeX Export --- *)
Print["\n=== TEX EXPORT ==="];

Print["\n% Effective coupling (Block A contribution):"];
Print["\\mu_{\\text{eff}} = ", TeXForm[muEff // Simplify]];

Print["\n% Effective a_1 mass^2 (Block A contribution):"];
Print["m^2_{\\text{eff}}(a_1) = ", TeXForm[m2Eff // Simplify]];

Print["\n% Block A Schur complement correction matrix:"];
Print["\\Delta_{\\text{A}} = ", TeXForm[correctionA // Simplify]];

Print["\n% Inverse of Block A mass matrix:"];
Print["S_{cc,A}^{-1} = ", TeXForm[SccAInv // Simplify]];

(* --- 11. Numerical Verification --- *)
Print["\n=== NUMERICAL VERIFICATION ==="];
Print["Evaluating at alpha1=0, alpha2=-0.6, alpha3=1.0, kappa=1, B0=0.0001, delta1=0.5, k=0.503:"];

testRules = {
  alpha1 -> 0, alpha2 -> -0.6, alpha3 -> 1.0,
  kappa -> 1, B0 -> 0.0001, delta1 -> 0.5, k -> 0.5027
};

Print["  SccA = ", SccA /. testRules // N];
Print["  SccAInv = ", SccAInv /. testRules // N];
Print["  correctionA = ", correctionA /. testRules // N];
Print["  mu_GR = ", muGR /. testRules // N];
Print["  mu_eff = ", muEff /. testRules // N];
Print["  m2_eff = ", m2Eff /. testRules // N];

(* Compare with Python: at alpha2=-0.6, delta1=0.5, k=0.503:
   coupling shift = -5.32e-5 (imaginary part)
   mass shift = +0.81 *)
Print["\n  Expected from Python:"];
Print["    coupling shift (Im) ~ -5.3e-5"];
Print["    mass shift ~ +0.81"];

Print["\nDone."];

Quit[];
