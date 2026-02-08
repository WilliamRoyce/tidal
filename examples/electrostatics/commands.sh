#!/usr/bin/env bash
# CLI equivalents for the Electrostatics 2D example
# See also: electrostatics.wls (manual derivation), electrostatics_simulation.py (Python simulation)
#
# NOTE: Electrostatics solves the Poisson equation (no time evolution).
# The `tg simulate` command is designed for time-dependent PDEs and cannot
# directly solve elliptic/static problems. Use the Python script instead.

set -euo pipefail

# Derive equations (manual .wls only — Poisson equation)
# tg derive electrostatics.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/electrostatics_2d.json

# NOTE: `tg simulate` is not applicable — electrostatics is a boundary-value
# problem (Poisson equation), not an initial-value problem. See the Python
# script electrostatics_simulation.py for the iterative relaxation solver.
