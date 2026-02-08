#!/usr/bin/env bash
# CLI equivalents for the Curved Spacetime examples
# See also: de_sitter_kg.wls, conformal_kg_static.wls (manual derivations)
#           de_sitter_simulation.py, conformal_kg_simulation.py (Python simulations)
#
# To run manually:  cd examples/curved_spacetime && tg derive de_sitter.toml

set -euo pipefail
cd "$(dirname "$0")"

### De Sitter Klein-Gordon (2+1D, expanding universe) ###

# Derive equations from Lagrangian (requires wolframscript)
tg derive de_sitter.toml

# Inspect the equation system
tg inspect ../data/de_sitter_kg.json

# Run simulation (2D Gaussian pulse with Hubble friction)
tg simulate ../data/de_sitter_kg.json \
  --param dSH=0.1 --param dSm2=1.0 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-width 3.0 \
  --t-end 20.0 \
  --dt 0.01

### Conformal Static Klein-Gordon (1+1D, constant conformal factor) ###

# Derive equations from Lagrangian (requires wolframscript)
tg derive conformal_static.toml

# Inspect the equation system
tg inspect ../data/conformal_kg_static.json

# Run simulation (Gaussian pulse, effective mass from conformal factor)
tg simulate ../data/conformal_kg_static.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 20.0 \
  --dt 0.002
