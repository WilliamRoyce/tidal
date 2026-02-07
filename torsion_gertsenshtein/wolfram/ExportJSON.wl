(* ::Package:: *)
(*
   MODULE: ExportJSON.wl
   PURPOSE: Convert component equations to multi-field JSON format for Python pipeline

   DEPENDENCIES:
     - TorsionGertsenshtein`CommonUtilities` (coefficient extraction, field parsing)

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

   Part of the torsion-gertsenshtein Lagrangian-to-PDE pipeline
*)

BeginPackage["TorsionGertsenshtein`ExportJSON`",
  {"TorsionGertsenshtein`CommonUtilities`"}];

(* Public symbols *)
BuildJSONStructure::usage =
  "BuildJSONStructure[componentEqs, metadata] builds the JSON Association \
structure from component equations and metadata.";

EquationToJSON::usage =
  "EquationToJSON[componentEq, fieldName, fieldIndex] converts a single \
component equation to JSON format.";

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
ExtractSpatialGradientFromMixed::usage =
  "ExtractSpatialGradientFromMixed[term] extracts the spatial gradient direction \
from a mixed time-space derivative term. Returns \"gradient_x\", \"gradient_y\", or \
\"gradient_z\". Warns if term has multiple spatial derivatives (e.g., d_t d_x d_y).";

IsSpatialCrossDerivative::usage =
  "IsSpatialCrossDerivative[term] returns True if term contains a spatial cross-derivative \
(d_x d_y, d_x d_z, or d_y d_z with time slot = 0). Legacy boolean wrapper around \
IdentifySpatialCrossDerivative.";

IdentifySpatialCrossDerivative::usage =
  "IdentifySpatialCrossDerivative[term] identifies spatial cross-derivatives and returns \
the operator name: \"cross_derivative_xy\", \"cross_derivative_xz\", \"cross_derivative_yz\", \
or False if not a cross-derivative. Requires time slot = 0.";

IdentifyDirectionalLaplacian::usage =
  "IdentifyDirectionalLaplacian[term] identifies pure second derivatives in a single \
spatial direction. Returns \"laplacian_x\", \"laplacian_y\", \"laplacian_z\", or False. \
Used for anisotropic equations like Navier-Cauchy where d^2_x and d^2_y have different coefficients.";

(* Generic derivative order support (Phase 12) *)
ExtractSpatialDerivativeProfile::usage =
  "ExtractSpatialDerivativeProfile[term] extracts the spatial derivative orders from a \
Derivative[dt, dx, dy, ...] pattern. Returns a list of non-negative integers representing \
the derivative order in each spatial axis (time slot excluded). E.g., Derivative[0,3,0] -> {3,0}.";

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
  "BuildTermResult[coeff, operator, field, symbolicCoeff:None, timeDependent:False, coordDeps:{}] builds a term result \
Association with keys \"coefficient\", \"operator\", \"field\". Optionally includes \
\"coefficient_symbolic\" when symbolicCoeff is not None, \"time_dependent\" when timeDependent is True, \
and \"coordinate_dependent\" when coordDeps is non-empty.";

IdentifyMultiFieldTerm::usage =
  "IdentifyMultiFieldTerm[term, currentFieldName, allFieldNames] analyzes a \
single term and returns an Association with \"coefficient\", \"operator\", and \"field\" keys. \
Detects: laplacian, directional laplacians, gradients, cross-derivatives, identity, and \
momentum gradient terms from mixed time-space derivatives.";

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
  "ExtractTermCoefficient[term, fieldHead, targetField] extracts the numeric coefficient \
and optional symbolic coefficient from a term. Returns {numericCoeff, symbolicCoeff, isTimeDependent, coordDeps}.";

ExtractTermCoefficient::symbolic =
  "Symbolic (non-numeric) coefficient `1` found for field `2`. Storing as symbolic coefficient in JSON.";

ClassifyOperatorType::usage =
  "ClassifyOperatorType[term] classifies the differential operator in a term. Returns \
a list {operatorName, isMixedTimeSpace} where operatorName is one of: \"identity\", \
\"laplacian\", \"laplacian_x/y/z\", \"gradient_x/y/z\", \"cross_derivative_xy/xz/yz\".";

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
  {json, fields, equations, allFieldNames, nFields, autoMatrices, couplingSection},

  (* fieldEquations format: {{"phi_0", eqPhi}, {"chi_0", eqChi}, ...} *)
  allFieldNames = fieldEquations[[All, 1]];
  nFields = Length[fieldEquations];

  (* Build fields list *)
  fields = Table[
    <|
      "name" -> fieldEquations[[i, 1]],
      "index" -> i - 1,
      "is_dynamical" -> True
    |>,
    {i, nFields}
  ];

  (* Convert equations to JSON format *)
  (* Pass allFieldNames so cross-field references can be detected *)
  equations = Table[
    EquationToJSONMultiField[
      fieldEquations[[i, 2]],
      fieldEquations[[i, 1]],
      i - 1,
      allFieldNames,
      metadata
    ],
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

  couplingSection = <|
    "mass_matrix" -> autoMatrices["mass_matrix"],
    "coupling_matrix" -> autoMatrices["coupling_matrix"]
  |>;
  If[autoMatrices["mass_matrix_symbolic"] =!= Null,
    couplingSection["mass_matrix_symbolic"] = autoMatrices["mass_matrix_symbolic"]
  ];
  If[autoMatrices["coupling_matrix_symbolic"] =!= Null,
    couplingSection["coupling_matrix_symbolic"] = autoMatrices["coupling_matrix_symbolic"]
  ];

  (* Build full JSON structure *)
  json = <|
    "metadata" -> <|
      "source" -> "xAct",
      "lagrangian_expr" -> Lookup[metadata, "lagrangian_expr", ""],
      "derived_from" -> "Euler-Lagrange",
      "gauge" -> Lookup[metadata, "gauge", "none"],
      "linearized" -> Lookup[metadata, "linearized", True]
    |>,
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
];

(* Equation conversion for multi-field systems *)
(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
EquationToJSONMultiField[componentEq_, fieldName_, fieldIndex_, allFieldNames_, metadata_] := Module[
  {terms, rhsTerms, rhs, timeDerivTerm, lhsTimeOrder, lhsStructure},

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

  (* LHS normalization: extract time-derivative coefficient and normalize RHS *)
  (* For |lhsCoeff| = 1: rhs = non-time terms as-is (handles both VarD and direct construction) *)
  (* For non-unit lhsCoeff (curved spacetime): rhs = -non_time_terms / lhsCoeff *)
  (*   Example: lhsCoeff = -Omega^{-2} gives rhs = Omega^2 * non_time_terms *)
  Module[{lhsCoeff},
    lhsCoeff = If[Length[timeDerivTerm] > 0,
      ExtractLHSCoefficient[timeDerivTerm[[1]]],
      -1
    ];
    If[Abs[lhsCoeff] =!= 1,
      rhs = Simplify[-rhs / lhsCoeff]
    ]
  ];

  (* Parse RHS with cross-field detection *)
  rhsTerms = ParseMultiFieldRHS[rhs, fieldName, allFieldNames];

  (* Build structured LHS for flexible PDE types *)
  lhsStructure = BuildLHSStructure[fieldName, lhsTimeOrder];

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

(* Parse RHS with cross-field reference detection *)
ParseMultiFieldRHS[eq_, currentFieldName_, allFieldNames_] := Module[
  {terms, parsedTerms},

  terms = If[Head[eq] === Plus, List @@ eq, {eq}];
  parsedTerms = Map[IdentifyMultiFieldTerm[#, currentFieldName, allFieldNames] &, terms];
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
      derivs = Cases[term,
        Derivative[n_, ___][f_][___] /; FunctionMatchesField[ToString[f], fieldName] :> n,
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
(* Supports: 1+1D (2-arg), 2+1D (3-arg), 3+1D (4-arg), and higher dimensions *)
ContainsTimeDerivative[term_, minOrder_:2] := Module[{},
  Which[
    (* 3+1D: Derivative[n, m, p, q] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _, _, _][_][___] /; n >= minOrder], True,
    (* 2+1D: Derivative[n, m, p] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _, _][_][___] /; n >= minOrder], True,
    (* 1+1D: Derivative[n, m] - first slot n is time *)
    !FreeQ[term, Derivative[n_, _][_][___] /; n >= minOrder], True,
    (* Generic: any arity with first slot >= minOrder *)
    !FreeQ[term, Derivative[n_, ___][_][___] /; n >= minOrder], True,
    (* Default: no time derivative of required order *)
    True, False
  ]
];

(* Field-aware time derivative check: only matches derivatives of the CURRENT field *)
(* Unlike ContainsTimeDerivative, this ignores cross-field time derivatives *)
(* Example: ContainsOwnTimeDerivative[-Derivative[2,0,0,0][gwH4][t,x,y,z], "h_0", 2] -> False *)
(* Example: ContainsOwnTimeDerivative[-Derivative[2,0,0,0][gwH0][t,x,y,z], "h_0", 2] -> True *)
ContainsOwnTimeDerivative[term_, fieldName_String, minOrder_Integer] := Module[
  {matchingDerivs},
  matchingDerivs = Cases[term,
    Derivative[n_, ___][f_][___] /; n >= minOrder && FunctionMatchesField[ToString[f], fieldName],
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
  Simplify[term / derivPart]
];

(* Check for mixed time-space derivatives that shouldn't be on RHS *)
(* Returns True if term has BOTH time AND space derivatives *)
(* Supports: 1+1D (2-arg), 2+1D (3-arg), 3+1D (4-arg) *)
IsMixedTimeSpaceDerivative[term_] := Module[{},
  Which[
    (* 3+1D: first slot > 0 AND (any spatial slot > 0) *)
    !FreeQ[term, Derivative[n_, m_, p_, q_][_][___] /; n > 0 && (m > 0 || p > 0 || q > 0)], True,
    (* 2+1D: first slot > 0 AND (second OR third slot > 0) *)
    !FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && (m > 0 || p > 0)], True,
    (* 1+1D: first slot > 0 AND second slot > 0 *)
    !FreeQ[term, Derivative[n_, m_][_][___] /; n > 0 && m > 0], True,
    (* Default: not a mixed derivative *)
    True, False
  ]
];

(* Identify gradient direction from derivative structure *)
(* In 1+1D: Derivative[dt, dx] - second slot is x *)
(* In 2+1D: Derivative[dt, dx, dy] - second slot is x, third is y *)
(* In 3+1D: Derivative[dt, dx, dy, dz] - second=x, third=y, fourth=z *)
(* IMPORTANT: First slot (time) must be 0 for pure spatial gradients *)
IdentifyGradientDirection[term_] := Module[{},
  Which[
    (* 3+1D: Check for gradient_z first (fourth slot nonzero, others zero) *)
    (* Pattern: Derivative[0, 0, 0, n_] where n > 0 *)
    !FreeQ[term, Derivative[0, 0, 0, n_][_][___] /; n > 0], "gradient_z",

    (* 3+1D: gradient_y (third slot nonzero, fourth zero) *)
    (* Pattern: Derivative[0, 0, n_, 0] where n > 0 *)
    !FreeQ[term, Derivative[0, 0, n_, 0][_][___] /; n > 0], "gradient_y",

    (* 3+1D: gradient_x (second slot nonzero, third and fourth zero) *)
    (* Pattern: Derivative[0, n_, 0, 0] where n > 0 *)
    !FreeQ[term, Derivative[0, n_, 0, 0][_][___] /; n > 0], "gradient_x",

    (* 2+1D: Check for gradient_y first (third slot nonzero) *)
    (* Pattern: Derivative[0, 0, n_] where n > 0 - MUST have first slot = 0 *)
    (* This ensures we only match pure spatial gradients, not mixed time-space *)
    !FreeQ[term, Derivative[0, 0, n_][_][___] /; n > 0], "gradient_y",

    (* 2+1D: gradient_x (second slot nonzero, third is zero) *)
    (* Pattern: Derivative[0, n_, 0] where n > 0 - MUST have first slot = 0 *)
    !FreeQ[term, Derivative[0, n_, 0][_][___] /; n > 0], "gradient_x",

    (* 1+1D: only gradient_x exists (second slot) *)
    (* Pattern: Derivative[0, n_] where n > 0 - MUST have first slot = 0 *)
    !FreeQ[term, Derivative[0, n_][_][___] /; n > 0], "gradient_x",

    (* Fail explicitly if gradient direction cannot be determined *)
    True,
      Throw[StringJoin[
        "IdentifyGradientDirection: Cannot determine gradient direction for term '",
        ToString[term], "'. ",
        "Expected Derivative[0, n_, ...] pattern for pure spatial gradient."
      ]]
  ]
];

(* === First-Order Time Derivative Detection (Phase 2 Curved Spacetime) === *)

(* Detect pure first-order time derivative: Derivative[1, 0, ...] with no spatial derivatives *)
(* Used for Hubble friction terms in curved spacetime: -2H∂_t φ *)
(* Returns True only if time order = 1 AND all spatial orders = 0 *)
IsFirstOrderTimeDerivative[term_] := Module[{},
  Which[
    (* 3+1D: Derivative[1, 0, 0, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0, 0, 0][_][___]], True,
    (* 2+1D: Derivative[1, 0, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0, 0][_][___]], True,
    (* 1+1D: Derivative[1, 0] - first-order time, no space *)
    !FreeQ[term, Derivative[1, 0][_][___]], True,
    (* Default: not a pure first-order time derivative *)
    True, False
  ]
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

(* Map field name to momentum field name *)
(* For py-pde state [field_0, pi_0, field_1, pi_1, ...], momentum is at odd indices *)
(* Uses ParseFieldName for flexible format support *)
FieldToMomentumName[fieldName_String] := Module[{parsed},
  parsed = ParseFieldName[fieldName];
  "pi_" <> ToString[parsed["index"]]
];

(* Extract spatial gradient direction from mixed time-space derivative *)
(* Derivative[1, 1, 0] (d_t d_x) -> "gradient_x" *)
(* Derivative[1, 0, 1] (d_t d_y) -> "gradient_y" *)
(* Derivative[1, 0, 0, 1] (d_t d_z) -> "gradient_z" *)
(* These represent d_x(pi) or d_y(pi) or d_z(pi) since d_t phi = pi *)
(*
   WARNING: Mixed time + cross-spatial derivatives like d_t d_x d_y or d_t d_x d_z
   cannot be simply represented as gradient(pi). They require cross_derivative(pi).
   For now, we default to gradient_x with a warning.
*)
ExtractSpatialGradientFromMixed[term_] := Module[{},
  (* Check for problematic mixed time + multiple spatial derivatives *)
  If[!FreeQ[term, Derivative[n_, m_, p_, q_][_][___] /; n > 0 && Count[{m, p, q}, _?(# > 0 &)] > 1],
    Print["Warning: Mixed time + cross-spatial derivative in 3+1D detected. ",
          "Cannot exactly represent as momentum gradient. Using gradient_x(pi)."];
  ];
  If[!FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && m > 0 && p > 0],
    Print["Warning: Mixed time + cross-spatial derivative detected (e.g., d_t d_x d_y). ",
          "This term cannot be exactly represented as a simple momentum gradient. ",
          "Approximating as gradient_x(pi). Consider gauge-fixing to eliminate such terms."];
  ];

  Which[
    (* 3+1D: z-only gradient (fourth slot > 0, second and third = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, 0, m_][_][___] /; n > 0 && m > 0], "gradient_z",
    (* 3+1D: y-only gradient (third slot > 0, second and fourth = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, m_, 0][_][___] /; n > 0 && m > 0], "gradient_y",
    (* 3+1D: x-only gradient (second slot > 0, third and fourth = 0, with time) *)
    !FreeQ[term, Derivative[n_, m_, 0, 0][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 2+1D: Check y-only gradient first (third slot > 0, second slot = 0, with time) *)
    !FreeQ[term, Derivative[n_, 0, m_][_][___] /; n > 0 && m > 0], "gradient_y",
    (* 2+1D: x-only gradient (second slot > 0, third slot = 0, with time) *)
    !FreeQ[term, Derivative[n_, m_, 0][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 2+1D: Both x and y present - default to gradient_x (warning above) *)
    !FreeQ[term, Derivative[n_, m_, p_][_][___] /; n > 0 && m > 0 && p > 0], "gradient_x",
    (* Fallback for other 2+1D/3+1D patterns - second slot nonzero *)
    !FreeQ[term, Derivative[n_, m_, _][_][___] /; n > 0 && m > 0], "gradient_x",
    (* 1+1D: Only x exists (second slot > 0 with first slot > 0) *)
    !FreeQ[term, Derivative[n_, m_][_][___] /; n > 0 && m > 0], "gradient_x",
    (* Fail explicitly if spatial gradient direction cannot be determined *)
    True,
      Throw[StringJoin[
        "ExtractSpatialGradientFromMixed: Cannot extract spatial gradient from mixed derivative '",
        ToString[term], "'. ",
        "Expected pattern Derivative[time>0, spatial...] with at least one nonzero spatial slot."
      ]]
  ]
];

(* Check if term contains a spatial cross-derivative (d_x d_y, d_x d_z, d_y d_z) *)
(* Returns False if not a cross-derivative, or the operator name if it is *)
(* Pattern: Derivative[0, ...] where exactly 2 spatial slots are > 0 *)

(* Legacy boolean version for backward compatibility *)
IsSpatialCrossDerivative[term_] := IdentifySpatialCrossDerivative[term] =!= False;

(* New version that returns the specific operator name *)
IdentifySpatialCrossDerivative[term_] := Module[{},
  Which[
    (* 3+1D: cross_derivative_yz (third and fourth slots > 0, second = 0) *)
    !FreeQ[term, Derivative[0, 0, m_, p_][_][___] /; m > 0 && p > 0], "cross_derivative_yz",
    (* 3+1D: cross_derivative_xz (second and fourth slots > 0, third = 0) *)
    !FreeQ[term, Derivative[0, m_, 0, p_][_][___] /; m > 0 && p > 0], "cross_derivative_xz",
    (* 3+1D: cross_derivative_xy (second and third slots > 0, fourth = 0) *)
    !FreeQ[term, Derivative[0, m_, p_, 0][_][___] /; m > 0 && p > 0], "cross_derivative_xy",
    (* 2+1D: cross_derivative_xy (second and third slots > 0) *)
    !FreeQ[term, Derivative[0, m_, p_][_][___] /; m > 0 && p > 0], "cross_derivative_xy",
    (* Not a cross-derivative *)
    True, False
  ]
];

(* === Phase 5 (Elasticity): Directional Laplacian Detection === *)
(* Identifies pure second derivatives in a single spatial direction *)
(* Returns: "laplacian_x", "laplacian_y", "laplacian_z", or False *)
(* Used for anisotropic equations like Navier-Cauchy where ∂²_x and ∂²_y have different coefficients *)
(* Supports 1+1D (2-arg), 2+1D (3-arg), and 3+1D (4-arg) Derivative forms *)

IdentifyDirectionalLaplacian[term_] := Module[{},
  Which[
    (* 3+1D: Derivative[0, 0, 0, 2] = pure ∂²/∂z² *)
    !FreeQ[term, Derivative[0, 0, 0, 2][_][___]], "laplacian_z",

    (* 3+1D: Derivative[0, 0, 2, 0] = pure ∂²/∂y² (no x or z derivative) *)
    !FreeQ[term, Derivative[0, 0, 2, 0][_][___]], "laplacian_y",

    (* 3+1D: Derivative[0, 2, 0, 0] = pure ∂²/∂x² (no y or z derivative) *)
    !FreeQ[term, Derivative[0, 2, 0, 0][_][___]], "laplacian_x",

    (* 2+1D: Derivative[0, 2, 0] = pure ∂²/∂x² (no y derivative) *)
    !FreeQ[term, Derivative[0, 2, 0][_][___]], "laplacian_x",

    (* 2+1D: Derivative[0, 0, 2] = pure ∂²/∂y² (no x derivative) *)
    !FreeQ[term, Derivative[0, 0, 2][_][___]], "laplacian_y",

    (* 1+1D: Derivative[0, 2] = pure ∂²/∂x² (only spatial dimension) *)
    !FreeQ[term, Derivative[0, 2][_][___]], "laplacian_x",

    (* Default: not a pure directional Laplacian *)
    True, False
  ]
];

(* Helper to build term result with optional symbolic coefficient and coordinate dependence *)
(* Only includes "coefficient_symbolic" key when symbolicCoeff is not None *)
(* Only includes "time_dependent" key when timeDependent is True *)
(* Only includes "coordinate_dependent" key when coordDeps is non-empty *)
BuildTermResult[coeff_, op_, field_, symbCoeff_:None, timeDependent_:False, coordDeps_:{}] := Module[{result},
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
  result
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

  (* For each field name like "A_0", "phi_0", find matching function head *)
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
ExtractTermCoefficient[term_, fieldHead_String, targetField_String] := Module[
  {rawCoeff, coefficient = 1.0, symbolicCoeff = None,
   isTimeDependent = False, coordDeps = {}},

  rawCoeff = term /. {
    (* Replace Derivative[...][field][args] with 1 *)
    Derivative[__][f_][__] /; ToString[f] === fieldHead :> 1,
    (* Replace field[args] with 1 *)
    f_Symbol[__] /; ToString[f] === fieldHead :> 1
  };
  rawCoeff = Simplify[rawCoeff];

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

  {coefficient, symbolicCoeff, isTimeDependent, coordDeps}
];

(* Count the total derivative order of a term *)
(* Returns the maximum total order of any Derivative expression found *)
CountDerivativeOrder[term_] := Module[
  {maxOrder},
  maxOrder = 0;
  (* Check for Derivative[n][f][args] patterns (applied derivatives) *)
  Cases[term,
    Derivative[orders__][_][__] :> (maxOrder = Max[maxOrder, Total[{orders}]]),
    {0, Infinity}
  ];
  (* Also check for unapplied Derivative[n][f] patterns *)
  If[maxOrder == 0,
    Cases[term,
      Derivative[orders__][_] :> (maxOrder = Max[maxOrder, Total[{orders}]]),
      {0, Infinity}
    ]
  ];
  maxOrder
];

(* === Phase 12: Generic Derivative Order Support === *)
(* Extract spatial derivative orders from a term *)
(* Returns a list of non-negative integers: {dx_order, dy_order, dz_order, ...} *)
(* Time slot (first in Derivative[dt, dx, ...]) is excluded *)
ExtractSpatialDerivativeProfile[term_] := Module[{profiles = {}},
  Cases[term,
    Derivative[orders__][_][__] :> Module[{orderList = Rest[{orders}]},
      AppendTo[profiles, orderList]
    ],
    {0, Infinity}
  ];
  If[Length[profiles] == 0,
    (* Also check unapplied Derivative[orders][f] *)
    Cases[term,
      Derivative[orders__][_] :> Module[{orderList = Rest[{orders}]},
        AppendTo[profiles, orderList]
      ],
      {0, Infinity}
    ]
  ];
  If[Length[profiles] == 0,
    Throw[StringJoin[
      "ExtractSpatialDerivativeProfile: No Derivative pattern found in term '",
      ToString[term, InputForm], "'."
    ]]
  ];
  profiles[[1]]
];

(* Build operator name from spatial derivative orders *)
(* Single-axis: {3,0} -> "derivative_3_x", {0,5} -> "derivative_5_y" *)
(* Multi-axis: {2,1} -> "derivative_2x_1y" *)
BuildGenericOperatorName[spatialOrders_List] := Module[
  {nonzero, axisNames = {"x", "y", "z"}, parts},

  (* Collect {axisIndex, order} pairs for nonzero spatial orders *)
  nonzero = {};
  Do[
    If[i <= Length[spatialOrders] && spatialOrders[[i]] > 0,
      AppendTo[nonzero, {i, spatialOrders[[i]]}]
    ],
    {i, 1, Min[Length[spatialOrders], 3]}
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

  (* Check for mixed time-space derivatives *)
  (* Mixed derivatives like d_t d_x A become gradient_x(pi) since d_t A = pi *)
  If[IsMixedTimeSpaceDerivative[term],
    Return[{ExtractSpatialGradientFromMixed[term], True}]
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

IdentifyMultiFieldTerm[term_, currentFieldName_, allFieldNames_] := Module[
  {functionHeads, matchResult, targetField, foundFieldHead,
   coeffResult, coefficient, symbolicCoeff, isTimeDependent, coordDeps,
   operatorResult, operator, isMixedTimeSpace},

  (* Step 1: Extract function heads from term *)
  functionHeads = ExtractFunctionHeads[term];

  (* Step 2: Match field names to function heads *)
  matchResult = MatchFieldToHeads[functionHeads, allFieldNames, currentFieldName];
  targetField = matchResult[[1]];
  foundFieldHead = matchResult[[2]];

  (* Step 3: Extract coefficient (includes coordinate-dependence check) *)
  If[foundFieldHead =!= Null,
    coeffResult = ExtractTermCoefficient[term, foundFieldHead, targetField];
    coefficient = coeffResult[[1]];
    symbolicCoeff = coeffResult[[2]];
    isTimeDependent = coeffResult[[3]];
    coordDeps = coeffResult[[4]],
    (* Fallback if no field head found *)
    coefficient = 1.0;
    symbolicCoeff = None;
    isTimeDependent = False;
    coordDeps = {}
  ];

  (* Step 4: Classify operator type *)
  operatorResult = ClassifyOperatorType[term];
  operator = operatorResult[[1]];
  isMixedTimeSpace = operatorResult[[2]];

  (* Convert to momentum field if mixed time-space derivative *)
  If[isMixedTimeSpace,
    targetField = FieldToMomentumName[targetField]
  ];

  (* Build and return result with coordinate-dependence info *)
  BuildTermResult[coefficient, operator, targetField, symbolicCoeff, isTimeDependent, coordDeps]
];

(* === Constraint Solver Hints (Issue #91) === *)
(* Builds constraint_solver JSON section for elliptic constraint equations *)
(* Only produces output when: timeOrder == 0 AND metadata has "solve_constraints" -> True *)
(* Returns Nothing when not applicable, so it integrates cleanly with Association building *)

ConstraintSolverHints[fieldName_String, timeOrder_Integer, metadata_Association] := Module[
  {enableSolver, bcHints, bcAssoc},

  (* Only applicable to constraint equations (time_order = 0) *)
  If[timeOrder =!= 0, Return[Nothing]];

  enableSolver = TrueQ[Lookup[metadata, "solve_constraints", False]];
  If[!enableSolver, Return[Nothing]];

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
    "method" -> "poisson",
    "boundary_conditions" -> bcAssoc
  |>
];

(* === Equation Conversion === *)

(* Phase 2, Issue 6: Now supports parabolic (d_t), elliptic (no time), and hyperbolic (d2_t) PDEs *)
(* Single-field EquationToJSON delegates to multi-field version *)
(* This avoids code duplication while maintaining backward compatibility *)
EquationToJSON[componentEq_, fieldName_, fieldIndex_, metadata_] :=
  EquationToJSONMultiField[componentEq, fieldName, fieldIndex, {fieldName}, metadata];

End[];
EndPackage[];
