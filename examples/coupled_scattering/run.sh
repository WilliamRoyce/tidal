#!/usr/bin/env bash
# Coupled Scalar Scattering 2+1D — Full derive → inspect → simulate → measure pipeline
#
# Physics: Two coupled scalar fields (phi, chi) in 2+1D Minkowski spacetime,
# coupled through a spatially localized Gaussian background field:
#   G(x,y) = g0 * exp(-(x^2 + y^2) / (2 R^2))
#
# The coupling is non-zero only near the origin. An incident phi wave packet
# scatters off the coupling region, partially converting to chi radiation.
# This demonstrates parametric mode conversion — directly analogous to the
# Gertsenshtein effect (phi=photon, chi=graviton, G=background B-field).
#
# See also:
#   coupled_scattering_simulation.py  — Python simulation with custom ICs + plots
#
# Running this script:
#   cd examples/coupled_scattering && uv run bash run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Step 1: Derive equations from the Lagrangian (requires wolframscript)
#   uv run tidal derive theory.toml
#
#   # Step 2: Inspect the generated equation system
#   uv run tidal inspect ../data/coupled_scattering.json
#
#   # Step 3: Simulate — phi wave packet hits coupling region, converts to chi
#   uv run tidal simulate ../data/coupled_scattering.json \
#     --param mPhi2=1.0 --param mChi2=4.0 --param g0=1.0 --param R=8.0 \
#     --grid-shape 128 --bounds=-50:50,-50:50 --periodic \
#     --ic gaussian --ic-component phi_0 --ic-center=-25.0,0.0 --ic-width 4.0 \
#     --t-end 40.0 --dt 0.02 --output ../data/coupled_scattering_output
#
#   # Step 4: Measure conversion probability and mixing length
#   uv run tidal measure ../data/coupled_scattering_output \
#     --what conversion,mixing --source phi_0 --target chi_0
#
#   # Step 5: Spectral conversion P(k,t)
#   uv run tidal measure ../data/coupled_scattering_output \
#     --what spectral_conversion --source phi_0 --target chi_0
#
#   # Step 6: Dispersion relation omega(k)
#   uv run tidal measure ../data/coupled_scattering_output \
#     --what dispersion --source phi_0
#
#   # Step 7: Combined measurement plot
#   uv run tidal measure ../data/coupled_scattering_output \
#     --what conversion,mixing,spectral_conversion,dispersion \
#     --source phi_0 --target chi_0 \
#     --output ../data/coupled_scattering_measurement.png

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (2 fields: phi_0, chi_0)
tidal inspect ../data/coupled_scattering.json

# Step 3: Run simulation
# Gaussian wave packet in phi at x=-25, propagating rightward (k0=3).
# Chi starts as a blob at origin. Coupling G(x,y) is Gaussian, radius R=8.
tidal simulate ../data/coupled_scattering.json \
  --param mPhi2=1.0 --param mChi2=4.0 --param g0=1.0 --param R=8.0 \
  --grid-shape 128 \
  --bounds=-50:50,-50:50 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center=-25.0,0.0 \
  --ic-width 4.0 \
  --t-end 40.0 \
  --dt 0.02 \
  --output ../data/coupled_scattering_output

# Step 4: Measure conversion probability and characteristic mixing length
# P(t) = E_chi(t) / E_phi(0) tracks what fraction of phi's energy went to chi
tidal measure ../data/coupled_scattering_output \
  --what conversion,mixing \
  --source phi_0 --target chi_0

# Step 5: Spectral conversion P(k,t)
# Shows which Fourier modes participate in the energy conversion over time
tidal measure ../data/coupled_scattering_output \
  --what spectral_conversion \
  --source phi_0 --target chi_0

# Step 6: Dispersion relation omega(k)
# Extracts wave frequency vs wavenumber via spacetime FFT
tidal measure ../data/coupled_scattering_output \
  --what dispersion \
  --source phi_0

# Step 7: Combined measurement plot (all panels in one figure)
tidal measure ../data/coupled_scattering_output \
  --what conversion,mixing,spectral_conversion,dispersion \
  --source phi_0 --target chi_0 \
  --output ../data/coupled_scattering_measurement.png

# NOTE: Energy diagnostics (virial) will raise ValueError for this example
# because the coupling G(x,y) is position-dependent. The simulation and
# conversion measurements work correctly — only virial energy is unsupported.
echo "Pipeline complete. See coupled_scattering_simulation.py for custom IC + plots."
