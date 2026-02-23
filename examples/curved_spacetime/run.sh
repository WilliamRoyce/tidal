#!/usr/bin/env bash
# Curved Spacetime — Full derive → inspect → simulate → plot pipeline
#
# Physics: Two variants of Klein-Gordon in curved spacetime:
#   1. De Sitter (2+1D): Expanding universe with Hubble friction H.
#      Time-dependent coefficients from FRW metric: a(t) = exp(Ht).
#   2. Conformal static (1+1D): Constant conformal factor Ω².
#
# Running this script:
#   cd examples/curved_spacetime && bash run.sh
#
# Or run each step manually:
#   tidal derive de_sitter.toml
#   tidal inspect ../data/de_sitter_kg.json
#   tidal simulate ../data/de_sitter_kg.json --param dSH=0.1 --param dSm2=1.0 \
#     --grid-shape 64 --bounds 0:50 --periodic --ic gaussian --ic-width 3.0 \
#     --t-end 20.0 --output ../data/de_sitter_output
#   tidal plot ../data/de_sitter_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

### De Sitter Klein-Gordon (2+1D, expanding universe) ###

# Derive equations from Lagrangian (requires wolframscript)
tidal derive de_sitter.toml

# Inspect the equation system
tidal inspect ../data/de_sitter_kg.json

# Run simulation (2D Gaussian pulse with Hubble friction)
tidal simulate ../data/de_sitter_kg.json \
  --param dSH=0.1 --param dSm2=1.0 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-width 3.0 \
  --t-end 20.0 \
  --output ../data/de_sitter_output

# Visualize de Sitter results (plots saved into the simulation output directory)
tidal plot ../data/de_sitter_output --type snapshot --time-index 0 --quiet
tidal plot ../data/de_sitter_output --type snapshot --time-index -1 --quiet
tidal plot ../data/de_sitter_output --type amplitude --overlay 'exp(-0.1*t)' --quiet
tidal plot ../data/de_sitter_output --type profile --cross-section y=25.0 --quiet

### Conformal Static Klein-Gordon (1+1D, constant conformal factor) ###

# Derive equations from Lagrangian (requires wolframscript)
tidal derive conformal_static.toml

# Inspect the equation system
tidal inspect ../data/conformal_kg_static.json

# Run simulation (Gaussian pulse, effective mass from conformal factor)
tidal simulate ../data/conformal_kg_static.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 20.0 \
  --output ../data/conformal_kg_output

# Visualize conformal KG results (plots saved into the simulation output directory)
tidal plot ../data/conformal_kg_output --type heatmap --quiet
tidal plot ../data/conformal_kg_output --type profile --quiet
tidal plot ../data/conformal_kg_output --type amplitude --quiet
