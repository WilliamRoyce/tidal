#!/usr/bin/env wolframscript
(* Electrostatics: Poisson Equation from Lagrangian *)
(* Demonstrates constraint solver integration (Issue #91) *)
(*
   Physics: Electrostatic potential from a charge distribution.

   Lagrangian: L = -1/2 eta^{ij} d_i phi d_j phi - rho * phi
   (spatial-only Lagrangian, i,j run over spatial indices only)

   Euler-Lagrange equation:
     delta L / delta phi = 0
     => nabla^2 phi + rho = 0
     => nabla^2 phi = -rho  (Poisson equation)

   This is an elliptic constraint equation (time_order = 0).
   The constraint solver uses py-pde's solve_poisson_equation to solve it.

   Boundary conditions: Dirichlet phi=0 (grounded conducting box).

   NOTE: Since this is a purely spatial equation (no time dynamics for phi),
   the Wolfram pipeline derives it as a constraint. The 'solve_constraints'
   metadata flag enables the Python-side elliptic solver.
*)

Print["=== Electrostatics: Poisson Equation ==="];
Print[""];

(* For this simple example, the JSON is constructed analytically *)
(* since the Poisson equation nabla^2 phi = -rho is well-known. *)
(* A full xAct derivation would use: *)
(*   L = -1/2 CD[-i][phi[]] eta[i,j] CD[-j][phi[]] - rho[] phi[] *)
(* which gives EOM: nabla^2 phi + rho = 0 *)

Print["The electrostatic Lagrangian:"];
Print["  L = -1/2 (nabla phi)^2 - rho * phi"];
Print[""];
Print["Euler-Lagrange equation:"];
Print["  nabla^2 phi + rho = 0"];
Print["  => nabla^2 phi = -rho  (Poisson equation)"];
Print[""];

(* Build JSON structure manually with constraint solver hints *)
Print["Building JSON with constraint_solver hints..."];

json = <|
  "metadata" -> <|
    "source" -> "analytical",
    "lagrangian_expr" -> "1/2 (nabla phi)^2 + rho * phi",
    "derived_from" -> "Euler-Lagrange",
    "gauge" -> "none",
    "linearized" -> True,
    "construction" -> "analytical"
  |>,
  "spacetime" -> <|
    "dimension" -> 3,
    "signature" -> {-1, 1, 1},
    "coordinates" -> {"t", "x", "y"}
  |>,
  "fields" -> {
    <|"name" -> "phi", "index" -> 0, "is_dynamical" -> True|>,
    <|"name" -> "rho", "index" -> 1, "is_dynamical" -> True|>
  },
  "equations" -> {
    <|
      "field" -> "phi",
      "lhs" -> <|"expression" -> "phi", "order" -> <|"time" -> 0, "space" -> 0|>|>,
      "rhs" -> <|
        "type" -> "linear_combination",
        "terms" -> {
          <|"coefficient" -> 1.0, "operator" -> "laplacian", "field" -> "phi"|>,
          <|"coefficient" -> 1.0, "operator" -> "identity", "field" -> "rho"|>
        }
      |>,
      "constraint_solver" -> <|
        "enabled" -> True,
        "method" -> "poisson",
        "boundary_conditions" -> <|
          "x" -> <|"type" -> "dirichlet", "value" -> 0.0|>,
          "y" -> <|"type" -> "dirichlet", "value" -> 0.0|>
        |>
      |>
    |>,
    <|
      "field" -> "rho",
      "lhs" -> <|"expression" -> "d2_t(rho)", "order" -> <|"time" -> 2, "space" -> 0|>|>,
      "rhs" -> <|
        "type" -> "linear_combination",
        "terms" -> {
          <|"coefficient" -> 0.0, "operator" -> "identity", "field" -> "rho"|>
        }
      |>
    |>
  },
  "coupling" -> <|
    "mass_matrix" -> {{0.0, 0.0}, {0.0, 0.0}},
    "coupling_matrix" -> {{0.0, 0.0}, {0.0, 0.0}}
  |>
|>;

(* Export JSON *)
outputPath = FileNameJoin[{DirectoryName[$InputFileName], "..", "data", "electrostatics_2d.json"}];
Export[outputPath, json, "JSON"];
Print["JSON exported to: ", outputPath];

Print[""];
Print["=== Done ==="];
Print["Run the simulation with:"];
Print["  python examples/electrostatics/electrostatics_simulation.py"];
