#!/usr/bin/env bash
# CLI equivalents for the Elasticity (Navier-Cauchy) 2+1D example
# See also: navier_cauchy.wls (manual derivation), elasticity_from_lagrangian.py (Python simulation)
#
# NOTE: The elasticity Lagrangian uses component-level expressions (not covariant
# tensor notation) because the Lame parameters break Lorentz invariance. The derive
# step requires the manual .wls script; no TOML config is provided.

set -euo pipefail

# Derive equations (manual .wls only — no TOML equivalent)
# tg derive navier_cauchy.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/navier_cauchy_2d.json

# Run simulation (Gaussian pulse in ux displacement)
# Parameters rho, lambda, mu are baked into the JSON as numeric coefficients
tg simulate ../data/navier_cauchy_2d.json \
  --grid-shape 64 \
  --bounds 0:10 \
  --periodic \
  --ic gaussian \
  --ic-component ux_0 \
  --ic-width 1.0 \
  --t-end 3.0 \
  --dt 0.005
