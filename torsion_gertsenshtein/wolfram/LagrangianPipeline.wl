(* ::Package:: *)
(* LagrangianPipeline.wl - Main entry point for Lagrangian to PDE pipeline *)
(* Part of the torsion-gertsenshtein project *)

BeginPackage["TorsionGertsenshtein`LagrangianPipeline`",
  {"xAct`xTensor`", "xAct`xCoba`", "xAct`xPert`",
   "TorsionGertsenshtein`EulerLagrange`",
   "TorsionGertsenshtein`Linearize`",
   "TorsionGertsenshtein`ComponentDecompose`",
   "TorsionGertsenshtein`ExportJSON`"}];

(* Public symbols *)
ProcessLagrangian::usage =
  "ProcessLagrangian[lagrangian, field, covd, chart, opts] derives field \
equations from a Lagrangian and exports them to JSON format. Returns the \
JSON Association.";

SetupMinkowskiND::usage =
  "SetupMinkowskiND[dim] sets up dim-dimensional Minkowski spacetime. \
dim=2 gives 1+1D with {t,x}, dim=4 gives 3+1D with {t,x,y,z}. \
Returns {manifold, metric, covd, chart}.";

SetupMinkowski1D::usage =
  "SetupMinkowski1D[] sets up 1+1D Minkowski spacetime with coordinates {t, x}. \
Equivalent to SetupMinkowskiND[2]. Returns {manifold, metric, covd, chart}.";

SetupMinkowski3D::usage =
  "SetupMinkowski3D[] sets up 3+1D Minkowski spacetime with coordinates {t, x, y, z}. \
Equivalent to SetupMinkowskiND[4]. Returns {manifold, metric, covd, chart}.";

Begin["`Private`"];

(* === Spacetime Setup Helpers === *)

(* Parametric Minkowski spacetime setup *)
SetupMinkowskiND[dim_Integer] := Module[
  {manifoldSym, metricSym, chartSym, covd, indices, coords, signature},

  (* Generate unique symbols for this setup based on spatial dimension *)
  manifoldSym = Symbol["$Minkowski" <> ToString[dim - 1] <> "D"];
  metricSym = Symbol["$eta" <> ToString[dim - 1] <> "D"];
  chartSym = Symbol["$cart" <> ToString[dim - 1] <> "D"];

  (* Generate index symbols based on dimension *)
  indices = Take[{a, b, c, d, e, f, g, h, i, j, k, l}, Min[dim + 4, 12]];

  (* Generate coordinate symbols *)
  coords = Switch[dim,
    2, {t[], x[]},
    4, {t[], x[], y[], z[]},
    _, Table[Symbol["x" <> ToString[i]][], {i, 0, dim - 1}]
  ];

  (* Metric signature: (-1, +1, +1, ...) *)
  signature = DiagonalMatrix[Prepend[ConstantArray[1, dim - 1], -1]];

  (* Define manifold *)
  DefManifold[manifoldSym, dim, indices];

  (* Define metric with Minkowski signature *)
  DefMetric[-1, metricSym[-a, -b], manifoldSym,
    SymbolOfCovD -> {";", "\[Del]"},
    PrintAs -> "\[Eta]"];

  (* Get covariant derivative *)
  covd = CovDOfMetric[metricSym[-a, -b]];

  (* Define coordinate chart *)
  DefChart[chartSym, manifoldSym, Range[0, dim - 1], coords];

  (* Set metric components *)
  MetricInBasis[metricSym, -chartSym, signature];

  {manifoldSym, metricSym, covd, chartSym}
];

(* Backward compatibility wrappers *)
SetupMinkowski1D[] := SetupMinkowskiND[2];
SetupMinkowski3D[] := SetupMinkowskiND[4];

(* === Main Pipeline === *)

Options[ProcessLagrangian] = {
  "OutputPath" -> None,
  "Metadata" -> <||>,
  "Verbose" -> True
};

ProcessLagrangian[lagrangian_, field_, covd_, chart_, opts:OptionsPattern[]] := Module[
  {eom, linearizedEOM, componentEqs, metadata, json, outputPath, verbose},

  verbose = OptionValue["Verbose"];
  outputPath = OptionValue["OutputPath"];
  metadata = OptionValue["Metadata"];

  (* Step 1: Derive Euler-Lagrange equations *)
  If[verbose, Print["Step 1: Deriving Euler-Lagrange equations..."]];
  eom = EulerLagrangeEquation[lagrangian, field, covd];
  If[verbose, Print["  EOM: ", eom]];

  (* Step 2: Linearize (if needed) *)
  If[verbose, Print["Step 2: Linearizing equations..."]];
  linearizedEOM = LinearizeEquation[eom, field, 0];
  If[verbose, Print["  Linearized EOM: ", linearizedEOM]];

  (* Step 3: Decompose into components *)
  If[verbose, Print["Step 3: Decomposing into components..."]];
  componentEqs = DecomposeToComponents[linearizedEOM, field, chart];
  If[verbose,
    Print["  Found ", Length[componentEqs], " component equation(s)"];
    Do[Print["    Component ", eq[[1]], ": ", eq[[2]]], {eq, componentEqs}]
  ];

  (* Step 4: Build JSON structure *)
  If[verbose, Print["Step 4: Building JSON structure..."]];

  (* Merge user metadata with defaults *)
  metadata = Join[
    <|
      "field_prefix" -> If[TensorRank[field] === 0, "phi", "A"],
      "dimension" -> Length[IndicesOfChart[chart]],
      "signature" -> DiagonalMatrix[{-1}~Join~ConstantArray[1, Length[IndicesOfChart[chart]] - 1]] // Diagonal,
      "coordinates" -> Map[ToString, IndicesOfChart[chart]],
      "lagrangian_expr" -> ToString[lagrangian],
      "linearized" -> True
    |>,
    metadata
  ];

  json = BuildJSONStructure[componentEqs, metadata];

  (* Step 5: Export to file (if path provided) *)
  If[outputPath =!= None,
    If[verbose, Print["Step 5: Exporting to ", outputPath, "..."]];
    Export[outputPath, json, "JSON"];
    If[verbose, Print["  Export complete."]];
  ];

  If[verbose, Print["Pipeline complete!"]];

  json
];

End[];
EndPackage[];

(* Usage example (commented out for package import):

<< TorsionGertsenshtein`LagrangianPipeline`;

(* Setup 1+1D Minkowski spacetime *)
{M, eta, CD, cart} = SetupMinkowski1D[];

(* Define scalar field *)
DefTensor[phi[], M];
DefConstantSymbol[m2];

(* Klein-Gordon Lagrangian *)
KGLag = -1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - 1/2 m2 phi[]^2;

(* Process and export *)
json = ProcessLagrangian[KGLag, phi[], CD, cart,
  "OutputPath" -> "klein_gordon_1d.json",
  "Metadata" -> <|"mass_squared" -> m2|>
];

*)
