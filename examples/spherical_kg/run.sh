#!/usr/bin/env bash
# CLI equivalents for the Spherical Klein-Gordon 3+1D example
# See also: spherical_kg.wls (manual derivation), spherical_kg_simulation.py (Python simulation)
#
# NOTE: The derive step uses spherical coordinates (r, theta, phi) with a
# coordinate-dependent metric and trigonometric coefficient functions.
# Both TOML and manual .wls derivation are supported.
# The simulation works via CLI with --bc and --ic formula flags.
#
# To run manually:  cd examples/spherical_kg

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/spherical_kg.json

# Run simulation (Gaussian shell at r=3, Neumann in r and theta, periodic in phi)
# Coordinates: x=r, y=theta, z=phi
tidal simulate ../data/spherical_kg.json \
  --param spm2=0.5 \
  --grid-shape 64 \
  --bounds 0.5:8,0.05:3.09,0:6.283185 \
  --bc neumann,neumann,periodic \
  --ic formula \
  --ic-formula "np.exp(-(x - 3.0)**2 / 0.72)" \
  --t-end 5.0 \
  --dt 0.01
