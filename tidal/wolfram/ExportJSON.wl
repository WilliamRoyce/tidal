(* ::Package:: *)
(*
   MODULE: ExportJSON.wl
   PURPOSE: Convert component equations to multi-field JSON format for Python pipeline

   DEPENDENCIES:
     - None (standalone — uses only built-in Mathematica)

   DATA FLOW:
     ComponentEquations (from ComponentDecompose)
       → ParseMultiFieldRHS (identify operators, extract coefficients)
       → EquationToJSONMultiField (build structured equation)
       → BuildMultiFieldJSONStructure (assemble full JSON)
       → ExportEquationSystem (write to file)

   SUPPORTED OPERATORS:
     - identity: field without derivatives (mass/coupling terms)
     - laplacian: full spatial Laplacian ∇²
     - laplacian_x, laplacian_y, laplacian_z: directional second derivatives
     - gradient_x, gradient_y, gradient_z: first-order spatial gradients
     - cross_derivative_xy, cross_derivative_xz, cross_derivative_yz: mixed spatial
     - first_derivative_t: first-order time derivative (Hubble friction in curved spacetime)

   PDE TYPES (via LHS structure):
     - Elliptic (time_order=0): Poisson, Laplace
     - Parabolic (time_order=1): Heat, diffusion
     - Hyperbolic (time_order=2+): Wave equations

   MOMENTUM GRADIENTS:
     Mixed time-space derivatives (∂_t∂_x φ) are converted to spatial gradients
     of momentum fields: gradient_x(pi_i), where pi_i = ∂_t φ_i.

   Part of the TIDAL Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`ExportJSON`"];

(* Public symbols *)
BuildJSONStructure::usage =
  "BuildJSONStructure[componentEqs, metadata] builds the JSON Association \
structure from component equations and metadata.";

BuildMultiFieldJSONStructure::usage =
  "BuildMultiFieldJSONStructure[fieldEquations, metadata] builds JSON for systems \
with multiple independent fields. fieldEquations is a list of {fieldName, equation} pairs.";

ParseFieldName::usage =
  "ParseFieldName[name] parses a field name and returns an Association with \
\"base\", \"index\", and \"format\" keys. Supports multiple formats: \
\"standard\" (A_0), \"tensor\" (stress_xy_0), \"compact\" (phi0), \"simple\" (phi).";

$FieldNameFormats::usage =
  "$FieldNameFormats defines the supported field name format patterns: \
\"standard\" (base_index), \"tensor\" (base_component_index), \
\"compact\" (base+digits), \"simple\" (letters only).";

GetTimeDerivativeOrder::usage =
  "GetTimeDerivativeOrder[term] returns the order of time derivative in a term. \
Returns 0 if no time derivative, 1 for d_t, 2 for d2_t, etc. Time index is always \
the first slot in Derivative[dt, dx, ...].";

DetectLHSTimeOrder::usage =
  "DetectLHSTimeOrder[equation] detects the maximum time derivative order in the equation. \
Used to determine if PDE is elliptic (0), parabolic (1), or hyperbolic (2+).";

BuildLHSStructure::usage =
  "BuildLHSStructure[fieldName, timeOrder] builds a structured LHS representation. \
Returns an Association with \"expression\" and \"order\" keys for flexible PDE types.";

(* Time derivative detection helpers *)
ContainsTimeDerivative::usage =
  "ContainsTimeDerivative[term, minOrder:2] returns True if term contains a time \
derivative of at least minOrder. Time is always the first slot in Derivative[dt, dx, ...]. \
Supports 1+1D (2-arg), 2+1D (3-arg), and 3+1D (4-arg) Derivative forms.";

IsMixedTimeSpaceDerivative::usage =
  "IsMixedTimeSpaceDerivative[term] returns True if term has BOTH time AND space \
derivatives (e.g., Derivative[1,1,0] for d_t d_x). Used to detect terms that should \
be converted to momentum gradients. Supports 1+1D, 2+1D, and 3+1D.";

ExtractLHSCoefficient::usage =
  "ExtractLHSCoefficient[term] extracts the scalar coefficient multiplying a \
Derivative[...][f][...] subexpression. Works for any dimension and any derivative \
pattern. Returns coeff such that term = coeff * Derivative[...][f][...]. \
Used by LHS normalization to divide out non-unit time-derivative coefficients.";

IdentifyGradientDirection::usage =
  "IdentifyGradientDirection[term] returns the gradient direction (\"gradient_x\", \
\"gradient_y\", or \"gradient_z\") based on derivative structure. Time slot must be 0 \
for pure spatial gradients. Defaults to \"gradient_x\" if detection fails.";

(* Spatial derivative classification *)
ClassifySpatialProfile::usage =
  "ClassifySpatialProfile[spatialOrders] maps a spatial derivative order list to a \
canonical operator name. {1,0} -> \"gradient_x\", {2,0} -> \"laplacian_x\", \
{0,2} -> \"laplacian_y\", {1,1} -> \"cross_derivative_xy\", \
{3,0} -> \"derivative_3_x\", {2,1} -> \"derivative_2x_1y\". \
Delegates to BuildGenericOperatorName for non-canonical cases.";

ExtractSpatialOperatorFromMixed::usage =
  "ExtractSpatialOperatorFromMixed[term] classifies the spatial operator in a \
mixed time-space derivative term. Strips the time order (must be exactly 1) and \
classifies the spatial part using ClassifySpatialProfile. Returns canonical operator \
names like \"gradient_x\", \"laplacian_x\", \"cross_derivative_xy\", etc. \
Throws if time order is not 1.";

IdentifySpatialCrossDerivative::usage =
  "IdentifySpatialCrossDerivative[term] identifies spatial cross-derivatives and returns \
the operator name: \"cross_derivative_xy\", \"cross_derivative_xz\", \"cross_derivative_yz\", \
or False if not a cross-derivative. Requires time slot = 0.";

CountDerivativeOrder::usage =
  "CountDerivativeOrder[term] counts the total derivative order of a term. \
Delegates to ExtractDerivativeProfile and returns Total of the profile, \
or 0 if no Derivative pattern is found.";

IdentifyDirectionalLaplacian::usage =
  "IdentifyDirectionalLaplacian[term] identifies pure second derivatives in a single \
spatial direction. Returns \"laplacian_x\", \"laplacian_y\", \"laplacian_z\", or False. \
Used for anisotropic equations like Navier-Cauchy where d^2_x and d^2_y have different coefficients.";

(* Generic derivative order support (Phase 12) *)
ExtractDerivativeProfile::usage =
  "ExtractDerivativeProfile[term] extracts the full derivative order profile from a \
Derivative[dt, dx, dy, ...] pattern. Returns a list {dt, dx, dy, ...} including \
the time slot. E.g., Derivative[1,2,0] -> {1,2,0}. Returns {} if no Derivative found.";

ExtractSpatialDerivativeProfile::usage =
  "ExtractSpatialDerivativeProfile[term] extracts the spatial derivative orders from a \
Derivative[dt, dx, dy, ...] pattern. Returns a list of non-negative integers representing \
the derivative order in each spatial axis (time slot excluded). E.g., Derivative[0,3,0] -> {3,0}. \
Delegates to ExtractDerivativeProfile and takes Rest.";

BuildGenericOperatorName::usage =
  "BuildGenericOperatorName[spatialOrders] builds an operator name string from a list of \
spatial derivative orders. Single-axis: {3,0} -> \"derivative_3_x\", {0,5} -> \"derivative_5_y\". \
Multi-axis: {2,1} -> \"derivative_2x_1y\". Axis names: 1->x, 2->y, 3->z.";

(* First-order time derivative detection for curved spacetime *)
IsFirstOrderTimeDerivative::usage =
  "IsFirstOrderTimeDerivative[term] returns True if term contains exactly a first-order \
time derivative (no spatial derivatives). Used to detect Hubble friction terms like \
-2H∂_t φ in curved spacetime. Time index is always the first slot in Derivative[dt, dx, ...].";

(* Term identification and building *)
BuildTermResult::usage =
  "BuildTermResult[coeff, operator, field, symbolicCoeff:None, timeDependent:False, coordDeps:{}, orderInEps:0] builds a term result \
Association with keys \"coefficient\", \"operator\", \"field\". Optionally includes \
\"coefficient_symbolic\" when symbolicCoeff is not None, \"time_dependent\" when timeDependent is True, \
\"coordinate_dependent\" when coordDeps is non-empty, and \"order_in_eps\" when orderInEps > 0.";

ComputeOrderInEps::usage =
  "ComputeOrderInEps[expr, paramsList] returns the maximum total exponent of the small-parameter \
symbols across monomials of Expand[expr]. Returns 0 when paramsList is empty or expr has no \
dependence on any small parameter. Used by the Phase C perturbative reduction (v6) to tag each \
emitted operator term with its order in the small parameters of L = L0 + eps*L1 + eps^2*L2 + ...";

IdentifyMultiFieldTerm::usage =
  "IdentifyMultiFieldTerm[term, currentFieldName, allFieldNames, smallParams:{}] analyzes a \
single term and returns an Association with \"coefficient\", \"operator\", and \"field\" keys. \
When smallParams is a non-empty list of Wolfram symbols, also tags the term with \"order_in_eps\" \
computed from the symbolic coefficient. Detects: laplacian, directional laplacians, gradients, \
cross-derivatives, identity, and momentum gradient terms from mixed time-space derivatives.";

ConstraintSolverHints::usage =
  "ConstraintSolverHints[fieldName, timeOrder, metadata] builds a constraint_solver \
Association for constraint equations (timeOrder=0) when metadata contains \
\"solve_constraints\" -> True. Returns Nothing for non-constraint equations or when \
solver is not requested. Boundary conditions come from metadata key \
\"constraint_boundary_conditions\".";

ExtractMassCouplingFromEquations::usage =
  "ExtractMassCouplingFromEquations[equations, allFieldNames] scans the already-built \
equation JSON Associations for identity operator terms and builds mass_matrix and \
coupling_matrix. Returns an Association with keys \"mass_matrix\", \"coupling_matrix\", \
and optionally \"mass_matrix_symbolic\" and \"coupling_matrix_symbolic\" when any \
symbolic coefficients are present. Convention: matrix[i][j] = -(coefficient of \
identity(field_j) in equation_i).";


(* Refactored helper functions (Phase 3, Issue 10) *)
ExtractFunctionHeads::usage =
  "ExtractFunctionHeads[term] extracts all function head names from an expression. \
Returns a list of strings representing function names (e.g., {\"phi0\", \"csA1\"}).";

MatchFieldToHeads::usage =
  "MatchFieldToHeads[functionHeads, allFieldNames, defaultField] matches function heads \
to known field names using case-insensitive and prefix matching. Returns an Association \
with \"field\" and \"head\" keys, or uses defaultField if no match found.";

ExtractTermCoefficient::usage =
  "ExtractTermCoefficient[term, fieldHead, targetField, smallParams:{}] extracts the numeric coefficient \
and optional symbolic coefficient from a term. Returns {numericCoeff, symbolicCoeff, isTimeDependent, coordDeps, orderInEps}. \
When smallParams is a non-empty list of Wolfram symbols, orderInEps is the total exponent of those symbols in the raw \
(pre-ToString) coefficient expression; otherwise 0.";

ExtractTermCoefficient::symbolic =
  "Symbolic (non-numeric) coefficient `1` found for field `2`. Storing as symbolic coefficient in JSON.";

ClassifyOperatorType::usage =
  "ClassifyOperatorType[term] classifies the differential operator in a term. Returns \
a list {operatorName, isMixedTimeSpace} where operatorName is one of: \"identity\", \
\"laplacian\", \"laplacian_x/y/z\", \"gradient_x/y/z\", \"cross_derivative_xy/xz/yz\".";

(* === Phase K: Canonical Momentum / Hamiltonian Export === *)

ClassifyHamiltonianFactor::usage =
  "ClassifyHamiltonianFactor[factor, allFieldNames] classifies a single field-dependent \
factor in a Hamiltonian term. Returns an Association with \"field\" (component name) \
and \"operator\" keys. operator is one of: \"identity\" (bare field), \
\"time_derivative\" (maps to canonical momentum at evaluation time), \
\"gradient_x\"/\"gradient_y\"/\"gradient_z\" (spatial gradients), \
\"laplacian_x\"/\"laplacian_y\"/\"laplacian_z\" (second spatial derivatives), \
or \"cross_derivative_xy\"/etc. for mixed spatial derivatives.";

ParseHamiltonianExpression::usage =
  "ParseHamiltonianExpression[componentExpr, allFieldNames] parses a component-form \
Hamiltonian density (scalar expression in Derivative[...][f][...] form) into a list \
of structured quadratic term Associations. Each term has keys \"coefficient\", \
\"factor_a\" and \"factor_b\" (each with \"field\" and \"operator\" subkeys). \
Only quadratic terms (exactly 2 field factors) are parsed; non-quadratic terms are \
skipped with a warning. Used by the canonical momentum pipeline (Phase K) to export \
the Hamiltonian for Python-side energy evaluation.";

ParseMultiFieldRHS::usage =
  "ParseMultiFieldRHS[eq, currentFieldName, allFieldNames, smallParams:{}] parses a linear \
combination of operators applied to fields into a list of term Associations. \
Each term has keys \"coefficient\", \"operator\", and \"field\" (plus optional \
\"order_in_eps\" when smallParams is supplied). Used by both the equation export \
pipeline and the canonical field rate computation (Phase K).";


Begin["`Private`"];

(* === JSON Structure Building === *)

(* Single-field BuildJSONStructure delegates to multi-field version *)
(* Converts {index, equation} pairs to {fieldName, equation} pairs *)
BuildJSONStructure[componentEqs_, metadata_Association] := Module[
  {fieldEquations, prefix},
  prefix = metadata["field_prefix"];
  fieldEquations = Table[
    {prefix <> "_" <> ToString[eq[[1]]], eq[[2]]},
    {eq, componentEqs}
  ];
  BuildMultiFieldJSONStructure[fieldEquations, metadata]
];

(* === Auto-compute mass/coupling matrices from equation terms === *)
(* Scans the already-built equation Associations for identity operator terms *)
(* and extracts mass_matrix (diagonal: own-field mass) and coupling_matrix *)
(* (off-diagonal: cross-field coupling). *)
(* Convention: matrix[i][j] = -(coefficient of identity(field_j) in equation_i) *)
(* This ensures mass^2 is positive for standard Lagrangians where RHS has -m^2*phi. *)

ExtractMassCouplingFromEquations[equations_List, allFieldNames_List] := Module[
  {nFields, mass, coupling, massSymbolic, couplingSymbolic,
   terms, fieldRef, j, coeff, symb, hasSymbolicMass, hasSymbolicCoupling},

  nFields = Length[allFieldNames];
  mass = ConstantArray[0.0, {nFields, nFields}];
  coupling = ConstantArray[0.0, {nFields, nFields}];
  massSymbolic = ConstantArray[Null, {nFields, nFields}];
  couplingSymbolic = ConstantArray[Null, {nFields, nFields}];

  Do[
    terms = equations[[i]]["rhs"]["terms"];
    Do[
      If[term["operator"] === "identity",
        fieldRef = term["field"];
        j = FirstPosition[allFieldNames, fieldRef, Missing["NotFound"]];
        If[!MissingQ[j],
          j = j[[1]];
          coeff = -N[term["coefficient"]];
          symb = Lookup[term, "coefficient_symbolic", Null];
          If[i === j,
            mass[[i, j]] += coeff;
            If[symb =!= Null, massSymbolic[[i, j]] = symb],
            coupling[[i, j]] += coeff;
            If[symb =!= Null, couplingSymbolic[[i, j]] = symb]
          ]
        ]
      ],
      {term, terms}
    ],
    {i, nFields}
  ];

  hasSymbolicMass = AnyTrue[Flatten[massSymbolic], # =!= Null &];
  hasSymbolicCoupling = AnyTrue[Flatten[couplingSymbolic], # =!= Null &];

  <|
    "mass_matrix" -> mass,
    "coupling_matrix" -> coupling,
    "mass_matrix_symbolic" -> If[hasSymbolicMass, massSymbolic, Null],
    "coupling_matrix_symbolic" -> If[hasSymbolicCoupling, couplingSymbolic, Null]
  |>
];

(* === Multi-Field JSON Structure Building === *)

BuildMultiFieldJSONStructure[fieldEquations_List, metadata_Association] := Module[
  {json, fields, equations, allFieldNames, nFields, autoMatrices, couplingSection,
   workingEqs},

  (* fieldEquations format: {{"phi_0", eqPhi}, {"chi_0", eqChi}, ...} *)
  workingEqs = fieldEquations;
  allFieldNames = workingEqs[[All, 1]];
  nFields = Length[workingEqs];

  (* === Cross-field d²_t equation reassignment ===
     In multi-field tensor systems, a component equation may contain d²_t of a
     DIFFERENT field than the one it's nominally assigned to. For example, in
     Fierz-Pauli massive gravity the xx-component equation has d²_t(h_yy) and
     the yy-component equation has d²_t(h_xx). Detect this and reassign so each
     field gets the equation that actually evolves it.

     Component E-L guarantees correct field→equation assignment by
     construction (∂L/∂q_i → EOM_i), so no swap logic is needed. *)

  (* Build fields list — enrich with tensor metadata when available *)
  fields = Table[
    Module[{entry, fname, tensorMeta},
      fname = workingEqs[[i, 1]];
      entry = <|
        "name" -> fname,
        "index" -> i - 1,
        "is_dynamical" -> True
      |>;
      (* Inject tensor component metadata if supplied by derive pipeline *)
      tensorMeta = Lookup[Lookup[metadata, "component_metadata", <||>], fname, Null];
      If[tensorMeta =!= Null,
        entry["tensor_head"] = tensorMeta["head"];
        entry["tensor_rank"] = tensorMeta["rank"];
        entry["tensor_indices"] = tensorMeta["indices"];
      ];
      entry
    ],
    {i, nFields}
  ];

  (* Convert equations to JSON format *)
  (* Pass allFieldNames so cross-field references can be detected *)
  equations = {};
  Do[
    Print["  JSON eq ", i, "/", nFields, ": ", workingEqs[[i, 1]],
      " (LeafCount=", LeafCount[workingEqs[[i, 2]]], ")"];
    AppendTo[equations,
      EquationToJSONMultiField[
        workingEqs[[i, 2]],
        workingEqs[[i, 1]],
        i - 1,
        allFieldNames,
        metadata
      ]
    ];
    Share[];
    ClearSystemCache[];,
    {i, nFields}
  ];

  (* Auto-compute mass/coupling matrices from equation identity terms *)
  autoMatrices = ExtractMassCouplingFromEquations[equations, allFieldNames];

  (* Warn if user metadata provides matrices — they will be ignored *)
  If[KeyExistsQ[metadata, "mass_matrix"],
    Print["WARNING: Ignoring user-provided mass_matrix in metadata. ",
          "Using auto-computed values from equation identity terms."]
  ];
  If[KeyExistsQ[metadata, "coupling_matrix"],
    Print["WARNING: Ignoring user-provided coupling_matrix in metadata. ",
          "Using auto-computed values from equation identity terms."]
  ];

  (* Only emit symbolic matrices — numeric matrices are redundant because
     the Wolfram pipeline cannot evaluate DefConstantSymbol values at export
     time (they produce ±1.0 shape factors, not actual parameter values).
     Python recomputes correct numeric matrices from terms + parameters. *)
  couplingSection = <||>;
  If[autoMatrices["mass_matrix_symbolic"] =!= Null,
    couplingSection["mass_matrix_symbolic"] = autoMatrices["mass_matrix_symbolic"]
  ];
  If[autoMatrices["coupling_matrix_symbolic"] =!= Null,
    couplingSection["coupling_matrix_symbolic"] = autoMatrices["coupling_matrix_symbolic"]
  ];

  (* Build full JSON structure *)
  Module[{jsonMeta, smallParamsList},
    jsonMeta = <|
      "source" -> "xAct",
      "lagrangian_expr" -> Lookup[metadata, "lagrangian_expr", ""],
      "derived_from" -> "Euler-Lagrange",
      "gauge" -> Lookup[metadata, "gauge", "none"],
      "linearized" -> Lookup[metadata, "linearized", True]
    |>;
    (* When [perturbation] is configured, emit the small-parameter names +    *)
    (* truncation order so Python's PerturbativeSolver knows which symbols to *)
    (* treat as small and the default --perturbative-order from the theory.   *)
    smallParamsList = Lookup[metadata, "small_parameters", {}];
    If[Length[smallParamsList] > 0,
      jsonMeta["perturbation"] = <|
        "small_parameters" -> Map[ToString[#, InputForm] &, smallParamsList],
        "order" -> Lookup[metadata, "perturbation_order", 1]
      |>
    ];
    json = <|
      "metadata" -> jsonMeta,
      "spacetime" -> <|
        "dimension" -> Lookup[metadata, "dimension", 2],
        "signature" -> Lookup[metadata, "signature", {-1, 1}],
        "coordinates" -> Lookup[metadata, "coordinates", {"t", "x"}]
      |>,
      "fields" -> fields,
      "equations" -> equations,
      "coupling" -> couplingSection
    |>;
    json
  ]
];

(* Equation conversion for multi-field systems *)
(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
EquationToJSONMultiField[componentEq_, fieldName_, fieldIndex_, allFieldNames_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm, lhsTimeOrder, lhsStructure, kineticCoeffStr, smallParams},
  kineticCoeffStr = None;
  (* Small-parameter symbols for order-in-eps tagging. The metadata key       *)
  (* "small_parameters" is a list of Wolfram symbols (e.g. {b5, rho})         *)
  (* supplied by _derive.py from the [perturbation] TOML section. Empty by   *)
  (* default for non-perturbative theories.                                   *)
  smallParams = Lookup[metadata, "small_parameters", {}];

  (* Same logic as EquationToJSON but with cross-field awareness *)
  terms = If[Head[componentEq] === Plus, List @@ componentEq, {componentEq}];

  (* Detect the time derivative order of the CURRENT field to determine PDE type *)
  (* Field-aware: only considers time derivatives of fieldName, not cross-field terms *)
  (* elliptic (0), parabolic (1), hyperbolic (2+) *)
  lhsTimeOrder = DetectLHSTimeOrder[componentEq, fieldName];

  (* Separate LHS (own-field time derivatives) from RHS *)
  If[lhsTimeOrder == 0,
    (* Elliptic/constraint: no time derivatives of this field — entire equation is RHS *)
    timeDerivTerm = {};
    rhs = componentEq;
    ,
    (* Parabolic/Hyperbolic: separate own-field time derivatives from RHS *)
    (* Only time derivatives of the CURRENT field are LHS candidates *)
    (* Cross-field time derivatives (e.g., d2_t(h_4) in h_0's equation) stay on RHS *)
    timeDerivTerm = Select[terms, ContainsOwnTimeDerivative[#, fieldName, lhsTimeOrder] &];
    (* RHS = everything EXCEPT own-field time derivatives of detected order *)
    (* Mixed time-space derivatives ARE included - they get converted to momentum gradients *)
    (* by IdentifyMultiFieldTerm (e.g., d_t d_x A -> gradient_x(pi)) *)
    rhs = Total[Select[terms, !ContainsOwnTimeDerivative[#, fieldName, lhsTimeOrder] &]];
  ];

  (* Safety net: if own-field time-derivative COEFFICIENTS sum to zero,
     reclassify as constraint. Catches edge cases where Simplify didn't
     fully cancel phantom d²_t terms in ComponentDecompose. *)
  If[lhsTimeOrder > 0 && Length[timeDerivTerm] > 0,
    Module[{coeffSum, coeffList, isPhantom},
      coeffList = ExtractLHSCoefficient /@ timeDerivTerm;
      coeffSum = If[LeafCount[Total[coeffList]] < 500, Simplify[Total[coeffList]], Total[coeffList]];
      isPhantom = (coeffSum === 0) || PossibleZeroQ[coeffSum];
      If[!isPhantom,
        Module[{syms, numVal},
          syms = DeleteDuplicates[Cases[coeffSum, _Symbol, {0, Infinity}]];
          numVal = Quiet[N[coeffSum /. Thread[syms -> Table[Prime[k]*E/7, {k, Length[syms]}]]]];
          If[NumericQ[numVal] && Abs[numVal] < 1e-10, isPhantom = True]
        ]
      ];
      If[isPhantom,
        Print["  NOTE: d" <> ToString[lhsTimeOrder] <> "_t(" <> fieldName <>
              ") coefficients cancel (sum=0). Reclassifying as constraint."];
        lhsTimeOrder = 0;
        timeDerivTerm = {};
        rhs = componentEq;
      ]
    ]
  ];

  (* LHS normalization: extract time-derivative coefficient and normalize RHS        *)
  (* For |lhsCoeff| = 1: rhs = non-time terms as-is                                *)
  (* For non-unit lhsCoeff with coord dependence (e.g. 1/Omega[r[]]): divide here  *)
  (* For non-unit lhsCoeff that is a pure parameter (e.g. xi, kappa^2):            *)
  (*   DO NOT divide — keep RHS unnormalized and record kinetic_coefficient_symbolic *)
  (*   so Python can handle the xi=0 degenerate case (field becomes constraint).    *)
  (*   Discrimination: FreeQ[lhsCoeff, _[]] is True for pure parameters,           *)
  (*   False for metric/coordinate-dependent expressions like Omega[r[]].           *)
  (*                                                                                *)
  (* #301 fix: lhsCoeff is the SUM of coefficients across every own-d²_t term, not  *)
  (* just the first. Without this, asymmetric kinetic corrections (e.g.             *)
  (* Euler-Heisenberg `(F·F̃)²` on a pure magnetic background) are silently dropped  *)
  (* — ρ happens to survive because `(F·F)²` is Lorentz-covariant around B̄ and the  *)
  (* correction rides on the spatial side too, but σ has no spatial partner. See    *)
  (* GitHub issue #301 and docs/tex/perturbative_reduction.tex §Validation.         *)
  Module[{lhsCoeff},
    lhsCoeff = If[Length[timeDerivTerm] > 0,
      Module[{rawSum},
        rawSum = Total[ExtractLHSCoefficient /@ timeDerivTerm];
        (* Expand rather than Simplify so the emitted string is a sum of
           monomials (no parenthesised sub-sums). The Python-side
           split_small_parameter_kinetic helper (#301 Phase 3) walks the
           AST term-by-term and cannot handle (a + b)*c forms — Expand
           distributes such structures up front. *)
        If[LeafCount[rawSum] < 500, Expand[rawSum], rawSum]
      ],
      -1
    ];
    If[Abs[lhsCoeff] =!= 1,
      If[!NumericQ[lhsCoeff] && FreeQ[lhsCoeff, _[]],
        (* Pure-parameter kinetic coeff (e.g., xi): keep RHS unnormalized.           *)
        (* Python normalizes at runtime so xi=0 can be handled as a constraint.     *)
        rhs = -rhs;
        kineticCoeffStr = ToString[lhsCoeff, InputForm],
        (* Coordinate-dependent (e.g., 1/Omega^2) or numeric: divide here as before *)
        rhs = -rhs / lhsCoeff
      ]
    ];
    (* ALWAYS Expand the RHS to ensure Plus structure for ParseMultiFieldRHS.
       Total[] (line 395) may trigger Mathematica's auto-factoring, collapsing
       separate linear terms into a single Times with multiple field heads.
       Without Expand, ParseMultiFieldRHS misidentifies these as bilinear
       (field-dependent) coefficients. This is critical for any multi-component
       system where the LHS coefficient is ±1 (standard wave equations).
       NOTE: Skip Expand for large expressions (LeafCount > 1000) — R̃²
       4th-order products in 4D exceed $RecursionLimit and crash the kernel
       via TerminatedEvaluation (uncatchable). For such expressions the
       component E-L already Expand'd each EOM term, so Plus structure holds. *)
    If[LeafCount[rhs] < 1000,
      rhs = Expand[rhs]
    ]
  ];

  (* Collect like field-derivative terms before parsing.
     After constraint substitution + Expand, the same operator-field pair
     may appear in multiple additive terms (e.g., 9 x laplacian(a_1) with
     separate coefficients).  Collect groups them into single terms with
     summed coefficients, reducing JSON bloat and preserving symbolic forms.
     Ref: #231 — uncombined duplicates from constraint elimination. *)
  (* Group like field-derivative terms by extracting each term's
     "field shape" (which field head + which derivative orders) and
     summing terms with the same shape.  Uses GatherBy instead of
     Collect because Mathematica's Collect uses structural template
     matching that fails for Derivative[__][f][__] patterns. *)
  Module[{fHeads, termList, getShape, grouped},
    fHeads = Union[
      Cases[rhs, f_Symbol[__] :> f, {0, Infinity}, Heads -> False],
      Cases[rhs, Derivative[__][f_][__] :> f, {0, Infinity}, Heads -> False]
    ];
    fHeads = Select[fHeads, Function[h,
      Module[{match},
        match = MatchFieldToHeads[{ToString[h]}, allFieldNames, ""];
        match[[1]] =!= ""
      ]
    ]];
    If[Length[fHeads] > 0 && LeafCount[rhs] < 50000 && Head[rhs] === Plus,
      termList = List @@ rhs;
      (* Extract field-derivative shape: {fieldHead, derivativeOrders} *)
      getShape[term_] := Module[{shapes},
        shapes = Join[
          Cases[term, f_Symbol[args__] /; MemberQ[fHeads, f] :> {f, {}}, {0, Infinity}],
          Cases[term, Derivative[orders__][f_][args__] /; MemberQ[fHeads, f] :> {f, {orders}}, {0, Infinity}]
        ];
        Sort[shapes]
      ];
      grouped = GatherBy[termList, getShape];
      rhs = Total[Total /@ grouped];
    ]
  ];

  (* Parse RHS with cross-field detection *)
  rhsTerms = ParseMultiFieldRHS[rhs, fieldName, allFieldNames, smallParams];

  (* Validate that no emitted term exceeds the declared perturbation_order.
     Under TIDAL's linearized-theory constraint, order_in_eps > 1 should
     never occur: a quadratic Lagrangian with a single linear coupling ε
     produces EOM coefficients with at most one factor of ε. A higher tag
     means the Lagrangian contains products of small parameters (e.g. ε·δ)
     or higher-power couplings (e.g. ε²) which are outside TIDAL's design
     scope and will be silently dropped by PerturbativeSolver when solving
     at order=1. *)
  Module[{declaredOrder, highOrderTerms},
    declaredOrder = Lookup[metadata, "perturbation_order", 1];
    highOrderTerms = Select[rhsTerms, Lookup[#, "order_in_eps", 0] > declaredOrder &];
    If[Length[highOrderTerms] > 0,
      Print["  WARNING [" <> fieldName <> "]: " <> ToString[Length[highOrderTerms]] <>
            " term(s) have order_in_eps > " <> ToString[declaredOrder] <>
            " (declared perturbation_order). TIDAL's linearized-theory pipeline " <>
            "only supports order_in_eps <= 1. These terms will be silently ignored " <>
            "by PerturbativeSolver. Check your Lagrangian for products of small " <>
            "parameters or higher-power couplings, which are not supported."]
    ]
  ];

  (* NOTE: The pi_ sign hack that was previously here (negating coefficients of
     momentum terms in wave equations) has been removed. The root cause — metric
     contractions absorbing into field tensors (V^a = g^{ab} V_b) and then mapping
     both covariant and contravariant to the same scalar — is now fixed in
     ComponentDecompose.wl via SeparateMetric applied before ToBasis. *)

  (* Build structured LHS for flexible PDE types *)
  lhsStructure = BuildLHSStructure[fieldName, lhsTimeOrder];

  (* Annotate LHS with kinetic coefficient when RHS was left unnormalized       *)
  (* (param-based kinetic coeff). Python uses this to normalize at runtime.    *)
  If[kineticCoeffStr =!= None,
    lhsStructure["kinetic_coefficient_symbolic"] = kineticCoeffStr
  ];

  Module[{result, constraintHints},
    result = <|
      "field" -> fieldName,
      "lhs" -> lhsStructure,  (* Now structured: {"expression": "...", "order": {...}} *)
      "rhs" -> <|
        "type" -> "linear_combination",
        "terms" -> rhsTerms
      |>
    |>;

    (* Add constraint_solver hints for elliptic equations when enabled *)
    constraintHints = ConstraintSolverHints[fieldName, lhsTimeOrder, metadata];
    If[constraintHints =!= Nothing,
      result["constraint_solver"] = constraintHints
    ];

    result
  ]
];

(* === RHS Term Parsing === *)

(* Parse RHS with cross-field reference detection *)
ParseMultiFieldRHS[eq_, currentFieldName_, allFieldNames_, smallParams_List:{}] := Module[
  {terms, parsedTerms},

  terms = If[Head[eq] === Plus, List @@ eq, {eq}];
  parsedTerms = Map[IdentifyMultiFieldTerm[#, currentFieldName, allFieldNames, smallParams] &, terms];
  parsedTerms = DeleteCases[parsedTerms, Nothing];

  (* Fail explicitly if no terms were parsed - do not inject ghost equations *)
  If[Length[parsedTerms] === 0,
    Throw[StringJoin[
      "ParseMultiFieldRHS: No terms parsed from RHS for field '", currentFieldName, "'. ",
      "RHS expression: ", ToString[eq], ". ",
      "Check that the equation is properly expanded and contains recognizable operators."
    ]]
  ];

  parsedTerms
];

(* === LHS Structure Detection (Phase 2, Issue 6) === *)
(* Supports parabolic (heat: d_t φ), elliptic (Poisson: ∇²φ = f), and hyperbolic (wave: d2_t φ) PDEs *)

(* Get time derivative order from a single term *)
(* Time index is always the first slot in Derivative[dt, dx, ...] *)
GetTimeDerivativeOrder[term_] := Module[{orders},
  (* Extract all time derivative orders (first slot of Derivative) *)
  orders = Cases[term, Derivative[n_, ___][_][___] :> n, {0, Infinity}];
  If[Length[orders] == 0, 0, Max[orders]]
];

(* Detect the maximum time derivative order in an equation *)
(* This determines the PDE type: elliptic (0), parabolic (1), hyperbolic (2+) *)
DetectLHSTimeOrder[equation_] := Module[{terms, maxOrder},
  terms = If[Head[equation] === Plus, List @@ equation, {equation}];
  maxOrder = Max[Map[GetTimeDerivativeOrder, terms]];
  (* Return at least 0 *)
  Max[maxOrder, 0]
];

(* === Field-Aware LHS Detection (for multi-field cross-coupled equations) === *)

(* Check if a function head string matches a specific field name *)
(* Uses same StringEndsQ+digit logic as MatchFieldToHeads *)
(* Example: FunctionMatchesField["gwH0", "h_0"] -> True *)
(* Example: FunctionMatchesField["gwH4", "h_0"] -> False *)
FunctionMatchesField[headStr_String, fieldName_String] := Module[
  {fieldParts, fieldBase, fieldIndex, headDigits, headBase},
  fieldParts = StringSplit[fieldName, "_"];
  If[Length[fieldParts] < 2, Return[False]];
  fieldBase = ToLowerCase[First[fieldParts]];
  fieldIndex = Last[fieldParts];
  headDigits = StringCases[headStr, RegularExpression["\\d+$"]];
  headBase = ToLowerCase[StringReplace[headStr, RegularExpression["\\d+$"] -> ""]];
  Length[headDigits] > 0 && headDigits[[-1]] === fieldIndex && StringEndsQ[headBase, fieldBase]
];

(* Field-aware overload: only considers time derivatives of the specified field *)
(* This is critical for multi-field systems where cross-field time derivatives *)
(* (e.g., d2_t(h_4) appearing in h_0's equation) must NOT be classified as LHS *)
DetectLHSTimeOrder[equation_, fieldName_String] := Module[
  {terms, fieldTermOrders, maxOrder},
  terms = If[Head[equation] === Plus, List @@ equation, {equation}];
  (* For each term, get time derivative order ONLY if it applies to the current field *)
  fieldTermOrders = Map[
    Function[term, Module[{derivs},
      (* Only match PURE time derivatives (all spatial orders = 0).         *)
      (* Mixed time-space like Derivative[2,2,0] must stay on RHS —        *)
      (* they represent ∂²_t ∂²_x terms from R̃², not pure ∂²_t.         *)
      derivs = Cases[term,
        Derivative[n_, rest___][f_][___] /;
          n > 0 && AllTrue[{rest}, # === 0 &] &&
          FunctionMatchesField[ToString[f], fieldName] :> n,
        {0, Infinity}];
      If[Length[derivs] == 0, 0, Max[derivs]]
    ]],
    terms
  ];
  maxOrder = If[Length[fieldTermOrders] == 0, 0, Max[fieldTermOrders]];
  Max[maxOrder, 0]
];

(* Build structured LHS representation *)
(* Supports variable time derivative orders for different PDE types *)
BuildLHSStructure[fieldName_String, timeOrder_Integer] := <|
  "expression" -> Switch[timeOrder,
    0, fieldName,                                    (* Elliptic: no time derivative *)
    1, "d_t(" <> fieldName <> ")",                   (* Parabolic: first order *)
    2, "d2_t(" <> fieldName <> ")",                  (* Hyperbolic: second order *)
    _, "d" <> ToString[timeOrder] <> "_t(" <> fieldName <> ")"  (* Higher order *)
  ],
  "order" -> <|"time" -> timeOrder, "space" -> 0|>
|>;

(* === Time Derivative Detection Helpers (for LHS/RHS separation) === *)

(* Detect if a term contains time derivatives of at least order minOrder *)
(* Time index is always the first slot in Derivative[dt, dx, ...] *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile *)
ContainsTimeDerivative[term_, minOrder_:2] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  Length[profile] >= 1 && First[profile] >= minOrder
];

(* Field-aware time derivative check: only matches derivatives of the CURRENT field *)
(* Unlike ContainsTimeDerivative, this ignores cross-field time derivatives *)
(* Example: ContainsOwnTimeDerivative[-Derivative[2,0,0,0][gwH4][t,x,y,z], "h_0", 2] -> False *)
(* Example: ContainsOwnTimeDerivative[-Derivative[2,0,0,0][gwH0][t,x,y,z], "h_0", 2] -> True *)
ContainsOwnTimeDerivative[term_, fieldName_String, minOrder_Integer] := Module[
  {matchingDerivs},
  (* Only match PURE time derivatives (all spatial orders = 0).           *)
  (* Mixed time-space like Derivative[2,2,0] are spatial operators, not   *)
  (* LHS time derivatives. They stay on RHS for ParseMultiFieldRHS.      *)
  matchingDerivs = Cases[term,
    Derivative[n_, rest___][f_][___] /;
      n >= minOrder && AllTrue[{rest}, # === 0 &] &&
      FunctionMatchesField[ToString[f], fieldName],
    {0, Infinity}];
  Length[matchingDerivs] > 0
];

(* Extract the coefficient of the LHS time derivative term *)
(* For flat spacetime: term = -Derivative[2,0][f][t,x], returns -1 *)
(* For curved spacetime: term = -Exp[-2Ht]*Derivative[2,0,0][f][t,x,y], returns -Exp[-2Ht] *)
(* This generalizes the ad-hoc sign-flip logic to handle arbitrary LHS coefficients *)
ExtractLHSCoefficient[term_] := Module[{derivParts, derivPart},
  derivParts = Cases[{term}, Derivative[__][_][__], {0, Infinity}];
  If[Length[derivParts] == 0,
    Throw[StringJoin[
      "ExtractLHSCoefficient: No Derivative found in LHS term '",
      ToString[term, InputForm], "'. Cannot extract time-derivative coefficient."
    ]]
  ];
  derivPart = derivParts[[1]];
  If[term === derivPart, Return[1]];
  Module[{ratio = term / derivPart},
    If[LeafCount[ratio] < 500, Simplify[ratio], ratio]
  ]
];

(* Check for mixed time-space derivatives that shouldn't be on RHS *)
(* Returns True if term has BOTH time AND space derivatives *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile *)
IsMixedTimeSpaceDerivative[term_] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  (* Only first-order time + spatial qualifies as mixed (d_t f = v_f).    *)
  (* Higher-order time (d²_t d_x, d²_t d²_x) are genuine higher-order   *)
  (* mixed operators from R̃² Lagrangians, not velocity conversions.     *)
  Length[profile] >= 2 && First[profile] === 1 && Max[Rest[profile]] > 0
];

(* Identify gradient direction from derivative structure *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile + ClassifySpatialProfile *)
IdentifyGradientDirection[term_] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  If[Length[profile] < 2,
    Throw[StringJoin[
      "IdentifyGradientDirection: Cannot determine gradient direction for term '",
      ToString[term], "'. ",
      "Expected Derivative[0, n_, ...] pattern for pure spatial gradient."
    ]]
  ];
  ClassifySpatialProfile[Rest[profile]]
];

(* === First-Order Time Derivative Detection (Phase 2 Curved Spacetime) === *)

(* Detect pure first-order time derivative: Derivative[1, 0, ...] with no spatial derivatives *)
(* Used for Hubble friction terms in curved spacetime: -2H∂_t φ *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile *)
IsFirstOrderTimeDerivative[term_] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  Length[profile] >= 2 && First[profile] == 1 && Max[Rest[profile]] == 0
];

(* === Phase 4: Momentum Gradient Helpers === *)

(* === Flexible Field Name Parsing === *)
(* Supports multiple field naming conventions *)

(* Field name format patterns *)
$FieldNameFormats = <|
  "standard" -> RegularExpression["^([a-zA-Z]+)_([0-9]+)$"],       (* A_0, phi_1 *)
  "tensor" -> RegularExpression["^(.+)_([0-9]+)$"],                (* stress_xy_0, u_x_1 *)
  "compact" -> RegularExpression["^([a-zA-Z]+)([0-9]+)$"],         (* phi0, A1 *)
  "simple" -> RegularExpression["^[a-zA-Z]+$"]                     (* phi, psi *)
|>;

(* Parse field name into base, index, and format *)
ParseFieldName[name_String] := Module[{match},
  (* Try standard format: A_0, phi_1 *)
  match = StringCases[name, $FieldNameFormats["standard"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "standard"|>]
  ];

  (* Try tensor format: stress_xy_0, u_x_1 (greedy match for base) *)
  match = StringCases[name, $FieldNameFormats["tensor"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "tensor"|>]
  ];

  (* Try compact format: phi0, A1 *)
  match = StringCases[name, $FieldNameFormats["compact"] -> {"$1", "$2"}];
  If[Length[match] > 0,
    Return[<|"base" -> match[[1, 1]], "index" -> ToExpression[match[[1, 2]]], "format" -> "compact"|>]
  ];

  (* Try simple format: phi, psi (single field, no index) *)
  If[StringMatchQ[name, $FieldNameFormats["simple"]],
    Return[<|"base" -> name, "index" -> 0, "format" -> "simple"|>]
  ];

  (* No format matched - fail explicitly *)
  Throw[StringJoin[
    "ParseFieldName: Cannot parse field name '", name, "'. ",
    "Expected one of: standard (A_0), tensor (stress_xy_0), compact (phi0), or simple (phi)."
  ]]
];

(* Map field name to velocity slot name: v_{fieldName}.                  *)
(* This matches the StateLayout naming convention in the Python solver.  *)
(* For example: FieldToVelocityName["A_0", {"phi_0","A_0","A_1","A_2"}] *)
(*   -> "v_A_0"                                                          *)
FieldToVelocityName[fieldName_String, allFieldNames_List] := Module[{pos},
  pos = FirstPosition[allFieldNames, fieldName];
  If[MissingQ[pos],
    Throw[StringJoin[
      "FieldToVelocityName: field '", fieldName,
      "' not found in field list: ", ToString[allFieldNames]
    ]]
  ];
  "v_" <> fieldName
];

(* Classify spatial operator in a mixed time-space derivative *)
(* Strips the time order (must be exactly 1) and classifies the spatial part *)
(* using ClassifySpatialProfile. Since d_t phi = v, mixed derivatives become *)
(* spatial operators on the velocity field: d_t d_x phi = d_x(v), *)
(* d_t d_x d_y phi = d_x d_y(pi) = cross_derivative_xy(pi), etc. *)
ExtractSpatialOperatorFromMixed[term_] := Module[
  {profile, timeOrder, spatialOrders},
  profile = ExtractDerivativeProfile[term];
  If[Length[profile] == 0,
    Throw[StringJoin[
      "ExtractSpatialOperatorFromMixed: No Derivative pattern found in term '",
      ToString[term, InputForm], "'."
    ]]
  ];
  timeOrder = First[profile];
  spatialOrders = Rest[profile];
  (* Validate time order is exactly 1 *)
  (* d_t phi = pi is well-defined, but d2_t phi is the LHS (acceleration) *)
  If[timeOrder != 1,
    Throw[StringJoin[
      "ExtractSpatialOperatorFromMixed: Time derivative order is ",
      ToString[timeOrder], " (expected 1) in term '",
      ToString[term, InputForm], "'. ",
      "Only first-order time mixed derivatives can be converted to ",
      "momentum spatial operators (d_t phi = pi)."
    ]]
  ];
  ClassifySpatialProfile[spatialOrders]
];

(* Check if term contains a spatial cross-derivative (d_x d_y, d_x d_z, d_y d_z) *)
(* Returns False if not a cross-derivative, or the operator name if it is *)
(* Pattern: Derivative[0, ...] where exactly 2 spatial slots are > 0 *)

(* Returns the specific operator name *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile + ClassifySpatialProfile *)
IdentifySpatialCrossDerivative[term_] := Module[{profile, result},
  profile = ExtractDerivativeProfile[term];
  If[Length[profile] < 2 || First[profile] =!= 0, Return[False]];
  result = ClassifySpatialProfile[Rest[profile]];
  If[StringMatchQ[result, "cross_derivative_" ~~ __], result, False]
];

(* === Phase 5 (Elasticity): Directional Laplacian Detection === *)
(* Identifies pure second derivatives in a single spatial direction *)
(* Returns: "laplacian_x", "laplacian_y", "laplacian_z", or False *)
(* Dimension-agnostic: delegates to ExtractDerivativeProfile + ClassifySpatialProfile *)
IdentifyDirectionalLaplacian[term_] := Module[{profile, result},
  profile = ExtractDerivativeProfile[term];
  If[Length[profile] < 2 || First[profile] =!= 0, Return[False]];
  result = ClassifySpatialProfile[Rest[profile]];
  If[StringMatchQ[result, "laplacian_" ~~ _], result, False]
];

(* Helper to build term result with optional symbolic coefficient and coordinate dependence *)
(* Only includes "coefficient_symbolic" key when symbolicCoeff is not None *)
(* Only includes "time_dependent" key when timeDependent is True *)
(* Only includes "coordinate_dependent" key when coordDeps is non-empty *)
(* Only includes "order_in_eps" key when orderInEps > 0 *)
BuildTermResult[coeff_, op_, field_, symbCoeff_:None, timeDependent_:False, coordDeps_:{}, orderInEps_:0] := Module[{result},
  result = <|
    "coefficient" -> N[coeff],
    "operator" -> op,
    "field" -> field
  |>;
  If[symbCoeff =!= None,
    result["coefficient_symbolic"] = symbCoeff
  ];
  If[timeDependent === True,
    result["time_dependent"] = True
  ];
  If[Length[coordDeps] > 0,
    result["coordinate_dependent"] = coordDeps
  ];
  If[IntegerQ[orderInEps] && orderInEps > 0,
    result["order_in_eps"] = orderInEps
  ];
  result
];

(* Compute total order of an expression in a list of small-parameter symbols. *)
(* Returns the maximum total exponent across monomials of the Expand'd form.  *)
(* Example: ComputeOrderInEps[b5^2 + eps*b5, {b5, eps}] returns 2.            *)
(* Returns 0 when the params list is empty or the expression is independent   *)
(* of all small parameters. Used by Phase C (perturbative reduction v6) to    *)
(* tag each emitted OperatorTerm with its order in the small parameters.     *)
ComputeOrderInEps[expr_, paramsList_List] := Module[{expanded, monomials, monomialOrder},
  If[Length[paramsList] === 0, Return[0]];
  expanded = Expand[expr];
  monomials = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  monomialOrder[m_] := Total[Max[Exponent[m, #], 0] & /@ paramsList];
  Max[Append[monomialOrder /@ monomials, 0]]
];

(* Return the list of coordinate names that a coefficient depends on *)
(* Coordinate symbols from xCoba appear as f[] (zero-argument function calls) *)
(* Returns e.g. {"t"}, {"x", "y"}, or {} for constant coefficients *)
IsCoordinateDependentCoefficient[coeffExpr_] := Module[{coordSyms},
  coordSyms = Cases[coeffExpr, f_Symbol[] :> f, {0, Infinity}];
  Union[ToString /@ coordSyms]
];

(* Backward-compatible boolean: True if any coordinate dependence *)
IsTimeDependentCoefficient[coeffExpr_] :=
  Length[IsCoordinateDependentCoefficient[coeffExpr]] > 0;

(* === Refactored Helper Functions for IdentifyMultiFieldTerm === *)
(* Phase 3, Issue 10: Split 144-line function into smaller, testable units *)

(* Extract all function heads from a term *)
(* For f[args], extract f; for Derivative[n,m][f][args], extract f *)
ExtractFunctionHeads[term_] := Union[Join[
  (* Direct function heads: f[args] -> f *)
  Cases[term, f_Symbol[__] :> ToString[f], {0, Infinity}],
  (* Derivative heads: Derivative[n,m][f][args] -> f *)
  Cases[term, Derivative[__][f_][__] :> ToString[f], {0, Infinity}]
]];

(* Match a field name to function heads in the term *)
(* Returns {matchedFieldName, foundHeadString} or {defaultField, Null} if no match *)
MatchFieldToHeads[functionHeads_List, allFieldNames_List, defaultField_String] := Module[
  {matchedField = defaultField, foundFieldHead = Null},

  (* Priority 1: exact head -> field name mapping via Global`$FieldHeadMapping. *)
  (* This handles velocity function symbols whose names use absolute field    *)
  (* indices that don't match the field name's own subscript index.           *)
  (* Uses Global` context since this package is in Private`.                  *)
  (* Safe no-op when $FieldHeadMapping is undefined.                          *)
  If[AssociationQ[Global`$FieldHeadMapping],
    Do[
      If[KeyExistsQ[Global`$FieldHeadMapping, head] &&
         MemberQ[allFieldNames, Global`$FieldHeadMapping[head]],
        matchedField = Global`$FieldHeadMapping[head];
        foundFieldHead = head;
        Break[]
      ],
      {head, functionHeads}
    ];
    If[foundFieldHead =!= Null, Return[{matchedField, foundFieldHead}]]
  ];

  (* Priority 2: heuristic base+index matching *)
  (* Strategy: Match BOTH base name and index *)
  (* E.g., "A_0" matches "csA0" (base "A/a" + index "0") *)
  (* E.g., "phi_0" matches "cplPhi0" (base "phi" + index "0") *)

  Do[
    Module[{fieldParts, fieldBase, fieldIndex, headDigits, headBase},
      (* Split field name: "phi_0" -> {"phi", "0"}, "A_1" -> {"A", "1"} *)
      fieldParts = StringSplit[fn, "_"];
      fieldBase = ToLowerCase[First[fieldParts]];
      fieldIndex = Last[fieldParts];

      (* Check if any function head matches BOTH base name AND index *)
      Do[
        (* Extract trailing digits: "cplPhi0" -> "0", "csA2" -> "2" *)
        headDigits = StringCases[head, RegularExpression["\\d+$"]];
        (* Extract base: "cplPhi0" -> "cplphi", "csA2" -> "csa" (lowercase, no digits) *)
        headBase = ToLowerCase[StringReplace[head, RegularExpression["\\d+$"] -> ""]];

        (* Match if: head ends with field base name AND trailing digits match field index *)
        (* Using StringEndsQ instead of StringContainsQ to prevent "A" matching "alpha", "eta" *)
        If[Length[headDigits] > 0 &&
           headDigits[[-1]] === fieldIndex &&
           StringEndsQ[headBase, fieldBase],
          matchedField = fn;
          foundFieldHead = head;
          Break[]
        ],
        {head, functionHeads}
      ];
      If[foundFieldHead =!= Null, Break[]]
    ],
    {fn, allFieldNames}
  ];

  {matchedField, foundFieldHead}
];

(* Extract coefficient from term given the field head *)
(* Returns {numericCoeff, symbolicCoeff, isTimeDependent, coordDeps} *)
(* where symbolicCoeff is None or a string, isTimeDependent is True/False, *)
(* coordDeps is a list of coordinate name strings e.g. {"t"}, {"x", "y"} *)
ExtractTermCoefficient[term_, fieldHead_String, targetField_String, smallParams_List:{}] := Module[
  {rawCoeff, coefficient = 1.0, symbolicCoeff = None,
   isTimeDependent = False, coordDeps = {}, orderInEps = 0},

  rawCoeff = term /. {
    (* Replace Derivative[...][field][args] with 1 *)
    Derivative[__][f_][__] /; ToString[f] === fieldHead :> 1,
    (* Replace field[args] with 1 *)
    f_Symbol[__] /; ToString[f] === fieldHead :> 1
  };
  If[!NumericQ[rawCoeff] && LeafCount[rawCoeff] < 500, rawCoeff = Simplify[rawCoeff]];

  (* Compute total order in small parameters before string conversion         *)
  (* (needs the raw symbolic expression, not its ToString representation).    *)
  orderInEps = ComputeOrderInEps[rawCoeff, smallParams];

  (* Check for coordinate dependence in coefficient *)
  coordDeps = IsCoordinateDependentCoefficient[rawCoeff];
  isTimeDependent = MemberQ[coordDeps, "t"];

  (* Determine numeric coefficient and symbolic representation *)
  (* Use InputForm for ToString to get clean, machine-parseable strings *)
  Which[
    NumericQ[rawCoeff],
      coefficient = N[rawCoeff];
      symbolicCoeff = None,
    (* Negative symbolic: -m2 -> -1.0, store "-m2" *)
    MatchQ[rawCoeff, Times[-1, _Symbol]],
      coefficient = -1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField],
    (* Numeric * symbolic product: -2*m2 -> -2.0, store "-2*m2" *)
    Head[rawCoeff] === Times && NumericQ[First[rawCoeff]],
      coefficient = N[First[rawCoeff]];
      symbolicCoeff = ToString[rawCoeff, InputForm];
      Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField],
    (* Positive symbolic: m2 -> 1.0, store "m2" *)
    MatchQ[rawCoeff, _Symbol],
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField],
    (* Try numeric evaluation *)
    NumericQ[Quiet[N[rawCoeff]]],
      coefficient = Quiet[N[rawCoeff]];
      symbolicCoeff = None,
    (* Complex symbolic expression (possibly coordinate-dependent) *)
    True,
      coefficient = 1.0;
      symbolicCoeff = ToString[rawCoeff, InputForm];
      If[Length[coordDeps] == 0,
        Message[ExtractTermCoefficient::symbolic, rawCoeff, targetField]
      ]
  ];

  {coefficient, symbolicCoeff, isTimeDependent, coordDeps, orderInEps}
];

(* Count the total derivative order of a term *)
(* Returns the total order of the Derivative expression found, or 0 if none *)
(* Delegates to ExtractDerivativeProfile to avoid duplicating extraction logic *)
CountDerivativeOrder[term_] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  (* Count SPATIAL derivative order only (exclude time slot = first).     *)
  (* Mixed time-space terms (time ≥ 2) are handled by BuildMixed-       *)
  (* OperatorName BEFORE this function is called. The terms reaching     *)
  (* here have time order 0 (pure spatial) or 1 (handled by IsMixed).   *)
  If[Length[profile] <= 1, 0, Total[Rest[profile]]]
];

(* === Generic Derivative Order Support === *)

(* Extract full derivative order profile from a term *)
(* Returns {dt, dx, dy, ...} including the time slot *)
(* Returns {} if no Derivative pattern found *)
ExtractDerivativeProfile[term_] := Module[{profiles = {}},
  Cases[term,
    Derivative[orders__][_][__] :> AppendTo[profiles, {orders}],
    {0, Infinity}
  ];
  If[Length[profiles] == 0,
    Cases[term,
      Derivative[orders__][_] :> AppendTo[profiles, {orders}],
      {0, Infinity}
    ]
  ];
  If[Length[profiles] == 0, Return[{}]];
  profiles[[1]]
];

(* Extract spatial derivative orders from a term *)
(* Delegates to ExtractDerivativeProfile and takes Rest (strips time slot) *)
(* Throws if no Derivative pattern found *)
ExtractSpatialDerivativeProfile[term_] := Module[{profile},
  profile = ExtractDerivativeProfile[term];
  If[Length[profile] == 0,
    Throw[StringJoin[
      "ExtractSpatialDerivativeProfile: No Derivative pattern found in term '",
      ToString[term, InputForm], "'."
    ]]
  ];
  Rest[profile]
];

(* Build operator name from spatial derivative orders *)
(* Single-axis: {3,0} -> "derivative_3_x", {0,5} -> "derivative_5_y" *)
(* Multi-axis: {2,1} -> "derivative_2x_1y" *)
BuildGenericOperatorName[spatialOrders_List] := Module[
  {nonzero, axisNames = {"x", "y", "z", "w", "v", "u"}, parts},

  (* Collect {axisIndex, order} pairs for nonzero spatial orders *)
  nonzero = {};
  Do[
    If[i <= Length[spatialOrders] && spatialOrders[[i]] > 0,
      If[i > Length[axisNames],
        Throw["BuildGenericOperatorName: Spatial dimension " <> ToString[i] <>
              " exceeds maximum supported (" <> ToString[Length[axisNames]] <> ")."]
      ];
      AppendTo[nonzero, {i, spatialOrders[[i]]}]
    ],
    {i, 1, Length[spatialOrders]}
  ];

  If[Length[nonzero] == 0,
    Throw[StringJoin[
      "BuildGenericOperatorName: No nonzero spatial derivatives in profile {",
      StringJoin[Riffle[ToString /@ spatialOrders, ", "]], "}."
    ]]
  ];

  If[Length[nonzero] == 1,
    (* Single axis: derivative_N_x *)
    StringJoin["derivative_", ToString[nonzero[[1, 2]]], "_", axisNames[[nonzero[[1, 1]]]]],
    (* Multi-axis: derivative_2x_1y *)
    parts = Table[
      StringJoin[ToString[pair[[2]]], axisNames[[pair[[1]]]]],
      {pair, nonzero}
    ];
    StringJoin["derivative_", Riffle[parts, "_"]]
  ]
];

(* Map spatial derivative order list to canonical operator name *)
(* Uses canonical names for common cases (gradient, laplacian, cross_derivative) *)
(* and delegates to BuildGenericOperatorName for higher-order/mixed cases *)
ClassifySpatialProfile[spatialOrders_List] := Module[
  {totalOrder, nonzeroPositions, axisNames = {"x", "y", "z", "w", "v", "u"}},

  totalOrder = Total[spatialOrders];
  nonzeroPositions = Flatten[Position[spatialOrders, _?(# > 0 &)]];

  If[totalOrder == 0,
    Throw["ClassifySpatialProfile: All spatial orders are zero. " <>
          "No spatial derivative to classify."]
  ];

  Which[
    (* Single axis with order 1 -> gradient *)
    Length[nonzeroPositions] == 1 && spatialOrders[[nonzeroPositions[[1]]]] == 1,
      "gradient_" <> axisNames[[nonzeroPositions[[1]]]],

    (* Single axis with order 2 -> directional laplacian *)
    Length[nonzeroPositions] == 1 && spatialOrders[[nonzeroPositions[[1]]]] == 2,
      "laplacian_" <> axisNames[[nonzeroPositions[[1]]]],

    (* Two axes, each with order 1 -> cross_derivative *)
    Length[nonzeroPositions] == 2 &&
      spatialOrders[[nonzeroPositions[[1]]]] == 1 &&
      spatialOrders[[nonzeroPositions[[2]]]] == 1,
      "cross_derivative_" <> axisNames[[nonzeroPositions[[1]]]] <>
                             axisNames[[nonzeroPositions[[2]]]],

    (* General case -> delegate to BuildGenericOperatorName *)
    True,
      BuildGenericOperatorName[spatialOrders]
  ]
];

(* Build operator name for mixed time-space derivatives.                 *)
(* Encodes the full profile: {2,2,0} → "mixed_T_S2_S0" (T=time²,      *)
(* S2=spatial order 2 in x, S0=spatial order 0 in y).                    *)
(* The solver's RHSEvaluator resolves "mixed_T_S..." operator names by  *)
(* applying the indicated spatial derivatives to the field and its       *)
(* time derivatives. For now, these terms contribute to the JSON but    *)
(* may need solver-side support for proper time integration.            *)
BuildMixedOperatorName[fullProfile_List] := Module[
  {timeOrder, spatialOrders, axisNames = {"x", "y", "z", "w"},
   parts},
  timeOrder = First[fullProfile];
  spatialOrders = Rest[fullProfile];
  parts = {"mixed_T" <> If[timeOrder > 1, ToString[timeOrder], ""]};
  Do[
    If[i <= Length[spatialOrders] && spatialOrders[[i]] > 0,
      AppendTo[parts, "S" <> ToString[spatialOrders[[i]]] <>
        If[i <= Length[axisNames], axisNames[[i]], "d" <> ToString[i]]]
    ],
    {i, Length[spatialOrders]}
  ];
  (* Pure time derivative on RHS (cross-field d²_t): name "d2_t" etc. *)
  If[Total[spatialOrders] === 0,
    "d" <> ToString[timeOrder] <> "_t",
    StringJoin[Riffle[parts, "_"]]
  ]
];

(* Classify operator type based on derivative structure *)
(* Returns {operatorName, shouldConvertToMomentum} *)
ClassifyOperatorType[term_] := Module[
  {orderDerivatives, dirLap, crossDeriv},

  (* Check for field without derivatives (mass/coupling term) *)
  If[FreeQ[term, Derivative],
    Return[{"identity", False}]
  ];

  (* Check for pure first-order time derivative (Hubble friction in curved spacetime) *)
  (* e.g., -2H∂_t φ in de Sitter: Derivative[1, 0][φ][t, x] with no spatial derivatives *)
  (* This must come BEFORE the mixed derivative check *)
  If[IsFirstOrderTimeDerivative[term],
    Return[{"first_derivative_t", False}]
  ];

  (* Check for mixed time-space derivatives (first-order time only).      *)
  (* d_t d_x A becomes gradient_x(pi) since d_t A = pi (velocity form). *)
  If[IsMixedTimeSpaceDerivative[term],
    Return[{ExtractSpatialOperatorFromMixed[term], True}]
  ];

  (* Higher-order time derivatives on the RHS (cross-field or mixed).    *)
  (* These arise from R̃² Lagrangians: d²_t d²_x h, d²_t h (cross-field).*)
  (* Build a generic operator name encoding the full derivative profile. *)
  Module[{fullProfile = ExtractDerivativeProfile[term]},
    If[Length[fullProfile] >= 1 && First[fullProfile] >= 2,
      Return[{BuildMixedOperatorName[fullProfile], False}]
    ]
  ];

  (* Check derivative order for pure spatial derivatives *)
  orderDerivatives = CountDerivativeOrder[term];

  (* Second-order spatial derivatives *)
  If[orderDerivatives == 2,
    (* Priority: directional Laplacians > cross-derivative > generic Laplacian *)
    dirLap = IdentifyDirectionalLaplacian[term];
    If[dirLap =!= False,
      Return[{dirLap, False}]
    ];
    crossDeriv = IdentifySpatialCrossDerivative[term];
    If[crossDeriv =!= False,
      Return[{crossDeriv, False}]
    ];
    Return[{"laplacian", False}]
  ];

  (* First-order spatial derivatives (gradients) *)
  If[orderDerivatives == 1,
    Return[{IdentifyGradientDirection[term], False}]
  ];

  (* Fourth-order spatial derivatives (biharmonic) *)
  If[orderDerivatives == 4,
    Return[{"biharmonic", False}]
  ];

  (* Generic higher-order derivatives (Phase 12) *)
  (* Handles orders 3, 5, 6+ by extracting per-axis derivative profile *)
  Module[{spatialProfile, genericName},
    spatialProfile = ExtractSpatialDerivativeProfile[term];
    genericName = BuildGenericOperatorName[spatialProfile];
    Return[{genericName, False}]
  ]
];

(* === Main Function: Identify operator and target field === *)
(* Refactored to use helper functions above *)

IdentifyMultiFieldTerm[term_, currentFieldName_, allFieldNames_, smallParams_List:{}] := Module[
  {functionHeads, matchResult, targetField, foundFieldHead,
   coeffResult, coefficient, symbolicCoeff, isTimeDependent, coordDeps, orderInEps,
   operatorResult, operator, isMixedTimeSpace},

  (* Skip zero terms — they arise from plane-wave reduction or vanishing *)
  (* components. Without this check, ClassifyOperatorType[0] returns     *)
  (* {"identity", False} and ExtractTermCoefficient returns 1.0,         *)
  (* producing the trivial equation field = 1.0 * identity(field).       *)
  If[term === 0 || term === 0., Return[Nothing]];

  (* Step 1: Extract function heads from term *)
  functionHeads = ExtractFunctionHeads[term];

  (* Step 2: Match field names to function heads *)
  matchResult = MatchFieldToHeads[functionHeads, allFieldNames, currentFieldName];
  targetField = matchResult[[1]];
  foundFieldHead = matchResult[[2]];

  (* Step 3: Extract coefficient (includes coordinate-dependence check) *)
  If[foundFieldHead =!= Null,
    coeffResult = ExtractTermCoefficient[term, foundFieldHead, targetField, smallParams];
    coefficient = coeffResult[[1]];
    symbolicCoeff = coeffResult[[2]];
    isTimeDependent = coeffResult[[3]];
    coordDeps = coeffResult[[4]];
    orderInEps = coeffResult[[5]],
    (* Fallback if no field head found *)
    coefficient = 1.0;
    symbolicCoeff = None;
    isTimeDependent = False;
    coordDeps = {};
    orderInEps = 0
  ];

  (* Step 4: Classify operator type *)
  operatorResult = ClassifyOperatorType[term];
  operator = operatorResult[[1]];
  isMixedTimeSpace = operatorResult[[2]];

  (* Convert to velocity field if mixed time-space derivative *)
  If[isMixedTimeSpace,
    targetField = FieldToVelocityName[targetField, allFieldNames]
  ];

  (* Step 5: Detect nonlinear (field-dependent) coefficients *)
  (* In a linear term, the field head should appear exactly once (bare or in Derivative). *)
  (* More than one occurrence means the coefficient depends on the field itself. *)
  Module[{fieldOccurrences, termResult},
    fieldOccurrences = Length[Cases[term,
      (f_Symbol[__] /; MemberQ[allFieldNames, ToString[f]]) |
      (Derivative[__][f_][__] /; MemberQ[allFieldNames, ToString[f]]),
      {0, Infinity}
    ]];

    (* Build term result *)
    termResult = BuildTermResult[coefficient, operator, targetField, symbolicCoeff, isTimeDependent, coordDeps, orderInEps];

    (* Attach warning if nonlinear *)
    If[fieldOccurrences > 1,
      Print["WARNING: Nonlinear (field-dependent) coefficient detected in term for field '",
        targetField, "'. The linear PDE solver will treat this as a constant coefficient."];
      termResult["warning"] = "nonlinear_coefficient"
    ];

    termResult
  ]
];

(* === Constraint Solver Hints (Issue #91) === *)
(* Builds constraint_solver JSON section for elliptic constraint equations *)
(* Only produces output when: timeOrder == 0 AND metadata has "solve_constraints" -> True *)
(* Returns Nothing when not applicable, so it integrates cleanly with Association building *)

ConstraintSolverHints[fieldName_String, timeOrder_Integer, metadata_Association] := Module[
  {enableSolver, bcHints, bcAssoc, method, maxIter, tol},

  (* Only applicable to constraint equations (time_order = 0) *)
  If[timeOrder =!= 0, Return[Nothing]];

  enableSolver = TrueQ[Lookup[metadata, "solve_constraints", False]];
  If[!enableSolver, Return[Nothing]];

  (* Read configurable solver parameters from metadata with defaults *)
  method = Lookup[metadata, "constraint_method", "auto"];
  maxIter = Lookup[metadata, "constraint_max_iterations", 30];
  tol = Lookup[metadata, "constraint_tolerance", 1*^-10];

  (* Build boundary conditions from metadata *)
  bcHints = Lookup[metadata, "constraint_boundary_conditions", <||>];

  (* Convert to JSON-compatible format *)
  (* Input format: <|"x" -> <|"type" -> "periodic"|>, "y" -> <|"type" -> "dirichlet", "value" -> 0.0|>|> *)
  bcAssoc = If[AssociationQ[bcHints],
    bcHints,
    <||>
  ];

  <|
    "enabled" -> True,
    "method" -> method,
    "max_iterations" -> maxIter,
    "tolerance" -> tol,
    "boundary_conditions" -> bcAssoc
  |>
];


(* === Phase K: Hamiltonian Term Parsing === *)
(*
  After computing H via LegendreTransformH and decomposing to component form
  via DecomposeScalarExpression, the Hamiltonian is a sum of quadratic terms
  in component fields and their derivatives.  Each term has the structure:

    coefficient * factor_a * factor_b

  where factor_a, factor_b are each:
    - fieldHead[t, x, y] (identity operator)
    - Derivative[dt, dx, dy][fieldHead][t, x, y] (derivative operator)

  For squared terms: coefficient * factor^2 (factor_a = factor_b).

  ClassifyHamiltonianFactor identifies the field name and operator type.
  ParseHamiltonianExpression parses the full expanded expression.
*)

(* Classify a single field-dependent factor in a Hamiltonian term *)
ClassifyHamiltonianFactor[factor_, allFieldNames_List] := Module[
  {functionHeads, matchResult, fieldName, profile, timeOrder, spatialOrders, operator},

  (* Extract function head and match to field name *)
  functionHeads = ExtractFunctionHeads[factor];
  matchResult = MatchFieldToHeads[functionHeads, allFieldNames, "unknown"];
  fieldName = matchResult[[1]];

  If[fieldName === "unknown",
    Throw[StringJoin[
      "ClassifyHamiltonianFactor: Cannot match factor '",
      ToString[factor, InputForm], "' to any known field name: ",
      ToString[allFieldNames]
    ]]
  ];

  (* Classify the operator based on derivative structure *)
  If[FreeQ[factor, Derivative],
    (* No derivative: identity (bare field value) *)
    operator = "identity",
    (* Has derivative: classify *)
    profile = ExtractDerivativeProfile[factor];
    If[Length[profile] == 0,
      Throw[StringJoin[
        "ClassifyHamiltonianFactor: Derivative pattern found but profile empty for '",
        ToString[factor, InputForm], "'."
      ]]
    ];
    timeOrder = First[profile];
    spatialOrders = Rest[profile];
    Which[
      (* Pure first-order time derivative: maps to canonical momentum at evaluation *)
      timeOrder == 1 && Max[spatialOrders] == 0,
        operator = "time_derivative",
      (* Pure spatial derivative *)
      timeOrder == 0,
        operator = ClassifySpatialProfile[spatialOrders],
      (* Mixed time-space: should not appear in properly simplified H *)
      True,
        Print["WARNING: Mixed time-space derivative in Hamiltonian factor: ",
          ToString[factor, InputForm]];
        operator = "mixed_" <> ToString[timeOrder] <> "_" <>
          StringJoin[Riffle[ToString /@ spatialOrders, "_"]]
    ]
  ];

  <|"field" -> fieldName, "operator" -> operator|>
];


(* Parse a single expanded term of the Hamiltonian into structured form *)
(* Returns an Association with coefficient, factor_a, factor_b, or Nothing if not quadratic *)
ParseSingleHamiltonianTerm[term_, fieldHeads_List, allFieldNames_List] := Module[
  {factors, fieldFactors = {}, coeffFactors = {},
   factorA, factorB, coefficient, numCoeff, symbolicCoeff = Null},

  (* Split Times expression into individual factors *)
  factors = If[Head[term] === Times, List @@ term, {term}];

  (* Classify each factor as field-dependent or coefficient *)
  Do[
    Which[
      (* Squared derivative: Power[Derivative[...][f][...], 2] *)
      MatchQ[factor, Power[Derivative[__][f_][__], n_Integer] /; MemberQ[fieldHeads, f] && n == 2],
        AppendTo[fieldFactors, factor[[1]]];
        AppendTo[fieldFactors, factor[[1]]],

      (* Squared bare field: Power[f[...], 2] *)
      MatchQ[factor, Power[f_Symbol[__], n_Integer] /; MemberQ[fieldHeads, f] && n == 2],
        AppendTo[fieldFactors, factor[[1]]];
        AppendTo[fieldFactors, factor[[1]]],

      (* Single derivative field factor: Derivative[...][f][...] *)
      MatchQ[factor, Derivative[__][f_][__] /; MemberQ[fieldHeads, f]],
        AppendTo[fieldFactors, factor],

      (* Single bare field factor: f[...] *)
      MatchQ[factor, f_Symbol[__] /; MemberQ[fieldHeads, f]],
        AppendTo[fieldFactors, factor],

      (* Everything else is part of the coefficient *)
      True,
        AppendTo[coeffFactors, factor]
    ],
    {factor, factors}
  ];

  (* Must have exactly 2 field factors for a quadratic term *)
  If[Length[fieldFactors] != 2,
    If[Length[fieldFactors] > 0,
      Print["WARNING: Hamiltonian term has ", Length[fieldFactors],
        " field factors (expected 2): ", ToString[Short[term], InputForm]]
    ];
    Return[Nothing]
  ];

  (* Compute coefficient (product of non-field factors) *)
  coefficient = If[Length[coeffFactors] > 0, Times @@ coeffFactors, 1];

  (* Separate numeric and symbolic parts *)
  If[NumericQ[coefficient],
    numCoeff = N[coefficient],
    (* Try N[] for simple symbolic expressions *)
    If[NumericQ[Quiet[N[coefficient]]],
      numCoeff = Quiet[N[coefficient]],
      numCoeff = 1.0;
      symbolicCoeff = ToString[coefficient, InputForm]
    ]
  ];

  (* Classify both field factors *)
  factorA = ClassifyHamiltonianFactor[fieldFactors[[1]], allFieldNames];
  factorB = ClassifyHamiltonianFactor[fieldFactors[[2]], allFieldNames];

  (* Build result *)
  Module[{result, coordDeps, termClass},
    (* Classify as self-energy or interaction based on base field *)
    termClass = If[factorA["field"] === factorB["field"], "self", "interaction"];
    result = <|
      "coefficient" -> numCoeff,
      "factor_a" -> factorA,
      "factor_b" -> factorB,
      "term_class" -> termClass
    |>;
    If[symbolicCoeff =!= Null,
      result["coefficient_symbolic"] = symbolicCoeff;
      (* Emit coordinate_dependent when symbolic coefficient contains spatial coords *)
      coordDeps = Select[IsCoordinateDependentCoefficient[coefficient], # =!= "t" &];
      If[Length[coordDeps] > 0,
        result["coordinate_dependent"] = coordDeps
      ]
    ];
    result
  ]
];


(* Parse the full expanded Hamiltonian expression into structured quadratic terms *)
ParseHamiltonianExpression[componentExpr_, allFieldNames_List, torsionPertName_String:""] := Module[
  {terms, fieldHeads, result},

  (* Discover all function heads that correspond to known fields *)
  fieldHeads = Union[Join[
    Cases[componentExpr, f_Symbol[__] :> f, {0, Infinity}],
    Cases[componentExpr, Derivative[__][f_][__] :> f, {0, Infinity}]
  ]];

  (* Filter to only heads that match known field names *)
  fieldHeads = Select[fieldHeads, Function[h,
    Module[{match},
      match = MatchFieldToHeads[{ToString[h]}, allFieldNames, ""];
      match[[1]] =!= ""
    ]
  ]];

  If[Length[fieldHeads] == 0,
    Throw[StringJoin[
      "ParseHamiltonianExpression: No known field heads found in expression. ",
      "Expected fields: ", ToString[allFieldNames], ". ",
      "Expression (short): ", ToString[Short[componentExpr, 3]]
    ]]
  ];

  (* Ensure expression is fully expanded before splitting into Plus terms.
     The Legendre transform H = Sum pi*vel - L should produce an expanded
     expression, but defensive Expand prevents silent failures from
     auto-factoring (same rationale as the Expand in EquationToJSONMultiField).
     NOTE: Skip Expand for very large expressions — R̃² 4th-order products
     exceed $RecursionLimit and crash via TerminatedEvaluation (uncatchable).
     Threshold increased from 1000 to 100000: the self-energy Hamiltonian
     for non-minimal torsion theories can have LeafCount ~5000-50000 and
     MUST be expanded for correct quadratic-term parsing. Only R̃²
     theories with 4th-order products (LeafCount >>100K) need to skip. *)
  If[LeafCount[componentExpr] < 100000,
    componentExpr = Expand[componentExpr]
  ];

  (* Split into additive terms *)
  terms = If[Head[componentExpr] === Plus, List @@ componentExpr, {componentExpr}];

  (* Parse each term *)
  result = Map[ParseSingleHamiltonianTerm[#, fieldHeads, allFieldNames] &, terms];
  result = DeleteCases[result, Nothing];

  (* Hamiltonian sector filter: exclude terms referencing torsion fields.
     $tidalHamiltonianFilter and $tidalTorsionHead are set in the WLS by
     _derive.py for torsion theories. This keeps only GW+EM self-energy
     and interaction terms — torsion self-energy and torsion cross-terms
     are dropped since they don't contribute to the conversion measurement
     C₀ = P/B₀² (which uses per-field self-energy only).
     See: hamiltonian_filter_design.md *)
  (* Torsion filter: passed as explicit argument to avoid Wolfram package
     scoping issues with global $-prefixed variables. *)
  If[StringLength[torsionPertName] > 0,
    Module[{nBefore = Length[result], tPert, nAfter},
      tPert = torsionPertName;
      (* Filter torsion terms from parsed Hamiltonian *)
      result = Select[result, Function[term,
        Module[{fieldA, fieldB, keep},
          fieldA = term[["factor_a", "field"]];
          fieldB = term[["factor_b", "field"]];
          keep = !StringMatchQ[fieldA, tPert ~~ "_" ~~ __] &&
                 !StringMatchQ[fieldB, tPert ~~ "_" ~~ __];
          keep
        ]
      ]];
      nAfter = Length[result];
      Print["  Hamiltonian torsion filter: ", nBefore, " -> ",
        nAfter, " terms (torsion fields excluded)"];
    ]
  ];

  If[Length[result] == 0,
    Throw[StringJoin[
      "ParseHamiltonianExpression: No quadratic terms found. ",
      "The Hamiltonian should be quadratic in fields/momenta. ",
      "Expression (short): ", ToString[Short[componentExpr, 3]]
    ]]
  ];

  result
];


End[];
EndPackage[];
