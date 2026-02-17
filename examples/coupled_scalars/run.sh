#!/usr/bin/env bash
# Coupled Scalars 1+1D — Full derive → inspect → simulate → measure pipeline
#
# Physics: Two coupled Klein-Gordon scalar fields (phi, chi) in 1+1D.
# The coupling term gCpl transfers energy from phi to chi over time,
# producing Rabi-like oscillations in the conversion probability P(t).
# The mixing length L_mix = pi/omega_dom characterizes the oscillation period.
#
# See also:
#   coupled_from_lagrangian.py   — Python simulation script
#   measure_conversion.py        — detailed measurement analysis
#
# Running this script:
#   cd examples/coupled_scalars && uv run bash run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Step 1: Derive equations from the Lagrangian (requires wolframscript)
#   uv run tidal derive theory.toml
#
#   # Step 2: Inspect the generated equation system
#   uv run tidal inspect ../data/coupled_scalars.json
#
#   # Step 3: Simulate with Gaussian IC in phi only (chi starts at zero)
#   uv run tidal simulate ../data/coupled_scalars.json \
#     --param mPhi2=1.0 --param mChi2=4.0 --param gCpl=0.5 \
#     --grid-shape 256 --bounds 0:100 --periodic \
#     --ic gaussian --ic-component phi_0 --ic-center 30.0 --ic-width 5.0 \
#     --t-end 20.0 --dt 0.01 --output ../data/coupled_scalars_output
#
#   # Step 4: Measure conversion probability and characteristic mixing length
#   uv run tidal measure ../data/coupled_scalars_output \
#     --what conversion,mixing --source phi_0 --target chi_0
#
#   # Step 5: Spectral conversion P(k,t) — which modes participate in mixing?
#   uv run tidal measure ../data/coupled_scalars_output \
#     --what spectral_conversion --source phi_0 --target chi_0
#
#   # Step 6: Dispersion relation omega(k) — wave frequency vs wavenumber
#   uv run tidal measure ../data/coupled_scalars_output \
#     --what dispersion --source phi_0
#
#   # Step 7: Combined plot with all measurements
#   uv run tidal measure ../data/coupled_scalars_output \
#     --what conversion,mixing,spectral_conversion,dispersion \
#     --source phi_0 --target chi_0 \
#     --output ../data/coupled_scalars_measurement.png

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (2 fields: phi_0, chi_0)
tidal inspect ../data/coupled_scalars.json

# Step 3: Run simulation
# Off-center Gaussian in phi only; coupling transfers energy to chi over time
tidal simulate ../data/coupled_scalars.json \
  --param mPhi2=1.0 --param mChi2=4.0 --param gCpl=0.5 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center 30.0 \
  --ic-width 5.0 \
  --t-end 20.0 \
  --dt 0.01 \
  --output ../data/coupled_scalars_output

# Step 4: Measure conversion probability and characteristic mixing length
# P(t) = E_chi(t) / E_phi(0) tracks what fraction of phi's energy went to chi
# L_mix = pi / omega_dom is the half-period of the dominant oscillation in P(t)
tidal measure ../data/coupled_scalars_output \
  --what conversion,mixing \
  --source phi_0 --target chi_0

# Step 5: Spectral conversion P(k,t)
# Shows which Fourier modes participate in the energy conversion over time
tidal measure ../data/coupled_scalars_output \
  --what spectral_conversion \
  --source phi_0 --target chi_0

# Step 6: Dispersion relation omega(k)
# Extracts wave frequency vs wavenumber via spacetime FFT
tidal measure ../data/coupled_scalars_output \
  --what dispersion \
  --source phi_0

# Step 7: Combined measurement plot (all panels in one figure)
tidal measure ../data/coupled_scalars_output \
  --what conversion,mixing,spectral_conversion,dispersion \
  --source phi_0 --target chi_0 \
  --output ../data/coupled_scalars_measurement.png
