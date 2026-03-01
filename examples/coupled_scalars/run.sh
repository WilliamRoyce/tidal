#!/usr/bin/env bash
# Coupled Scalars 1+1D — Full derive → inspect → simulate → measure pipeline
#
# Physics: Two coupled Klein-Gordon scalar fields (phi, chi) in 1+1D.
# The coupling term gCpl transfers energy from phi to chi over time,
# producing Rabi-like oscillations in the conversion probability P(t).
# The mixing length L_mix = pi/omega_dom characterizes the oscillation period.
#
# Running this script:
#   cd examples/coupled_scalars && bash run.sh
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
#     --ic-amplitude 1.0 \
#     --t-end 20.0 --output ../data/coupled_scalars_output
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
#     --what energy,conservation,conversion,mixing,spectral_conversion,dispersion \
#     --source phi_0 --target chi_0 \
#     --output ../data/coupled_scalars_output/measurement.png

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (2 fields: phi_0, chi_0)
tidal inspect ../data/coupled_scalars.json

# ============================================================================
# Run 1: Asymmetric — phi only (chi starts at zero)
# ============================================================================
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
  --ic-amplitude 1.0 \
  --t-end 20.0 \
  --output ../data/coupled_scalars_output

# --- Run 1 measurements ---

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
  --what energy,conservation,conversion,mixing,spectral_conversion,dispersion \
  --source phi_0 --target chi_0 \
  --output ../data/coupled_scalars_output/measurement.png

# --- Run 1 plots ---
tidal plot ../data/coupled_scalars_output --type heatmap --field phi_0 --quiet
tidal plot ../data/coupled_scalars_output --type heatmap --field chi_0 --quiet
tidal plot ../data/coupled_scalars_output --type amplitude --quiet
tidal plot ../data/coupled_scalars_output --type energy --quiet
tidal plot ../data/coupled_scalars_output --type compare --quiet
tidal plot ../data/coupled_scalars_output --type hamiltonian --quiet
tidal plot ../data/coupled_scalars_output --type conservation --quiet

echo ""
echo "=== Run 1 (asymmetric) complete ==="

# ============================================================================
# Run 2: Symmetric — both fields excited (using --ic-field)
# ============================================================================
# What happens when both fields start with energy? Using --ic-field, we set
# chi_0 to a Gaussian at x=70 (far from phi's Gaussian at x=30). The
# conversion dynamics are qualitatively different: both fields exchange
# energy bidirectionally from the start, rather than one-way transfer.

tidal simulate ../data/coupled_scalars.json \
  --param mPhi2=1.0 --param mChi2=4.0 --param gCpl=0.5 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center 30.0 \
  --ic-width 5.0 \
  --ic-amplitude 1.0 \
  --ic-field 'chi_0:exp(-(x - 70)**2 / 50)' \
  --t-end 20.0 \
  --output ../data/coupled_scalars_symmetric

# --- Run 2 measurements ---
tidal measure ../data/coupled_scalars_symmetric \
  --what conversion,mixing \
  --source phi_0 --target chi_0

tidal measure ../data/coupled_scalars_symmetric \
  --what energy,conservation,conversion,mixing \
  --source phi_0 --target chi_0 \
  --output ../data/coupled_scalars_symmetric/measurement.png

# --- Run 2 plots ---
tidal plot ../data/coupled_scalars_symmetric --type heatmap --field phi_0 \
  --title "phi (symmetric excitation)" --quiet
tidal plot ../data/coupled_scalars_symmetric --type heatmap --field chi_0 \
  --title "chi (symmetric excitation)" --quiet
tidal plot ../data/coupled_scalars_symmetric --type amplitude --quiet
tidal plot ../data/coupled_scalars_symmetric --type hamiltonian --quiet
tidal plot ../data/coupled_scalars_symmetric --type conservation --quiet

echo ""
echo "=== Run 2 (symmetric) complete ==="
echo ""
echo "Compare Run 1 (phi only) vs Run 2 (both fields excited)."
echo "  Asymmetric:  ../data/coupled_scalars_output/"
echo "  Symmetric:   ../data/coupled_scalars_symmetric/"
