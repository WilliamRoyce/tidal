#!/usr/bin/env bash
# Coupled Scalar Scattering — Plane-wave Reduced 1+1D Pipeline
#
# Physics: Two coupled scalar fields in 3+1D Minkowski spacetime, reduced
# to 1+1D via plane-wave ansatz (d/dx = d/dy = 0, propagation along z).
# Coupling is a localized Gaussian background: G(z) = g0 * exp(-z^2/(2R^2))
#
#   d2t(phi) = d2z(phi) - mPhi2 phi - G(z) chi
#   d2t(chi) = d2z(chi) - mChi2 chi - G(z) phi
#
# Analogy: phi=photon, chi=graviton, G(z)=background magnetic field.
# A phi wave packet (Gaussian envelope with carrier k=3) starts at x=-20,
# well outside the coupling region at the origin. The wavevector gives it
# rightward velocity, so it propagates cleanly into the coupling region
# and excites chi radiation (Gertsenshtein-effect analogue).
#
# The 1D reduction eliminates geometric spreading and 2D interference,
# giving clean conversion signals and unambiguous heatmaps.
#
# Stability condition: mPhi2 * mChi2 > g0^2 (positive-definite mass matrix).
#
# Running this script:
#   cd examples/coupled_scattering && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
# tidal derive theory.toml

# Inspect the equation system (1+1D reduced, 2 dynamical fields)
tidal inspect ../data/coupled_scattering.json

# ============================================================================
# Run 1: Resonant conversion (mPhi2 = mChi2 = 1.0)
# ============================================================================
# Equal masses give maximum mode conversion. A phi wave packet (Gaussian
# envelope, carrier k=3) starts at x=-20, well outside the coupling region
# (origin, R=5). It propagates rightward into the coupling and excites
# chi radiation via parametric conversion.

OUT1=../data/coupled_scattering_resonant

tidal simulate ../data/coupled_scattering.json \
  --param mPhi2=1.0 --param mChi2=1.0 --param g0=0.5 --param R=5.0 \
  --grid-shape 512 \
  --bounds=-40:40 \
  --periodic \
  --scheme leapfrog \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 80.0 \
  --output "$OUT1"

# --- Resonant plots ---

# Spacetime heatmaps: phi_0 and chi_0 side by side
# phi_0 shows the incident plane wave; chi_0 shows the generated wave
tidal plot "$OUT1" --type heatmap --field phi_0 \
  --title "phi (incident wave, resonant)" \
  --output "$OUT1/heatmap_phi_0.png"

tidal plot "$OUT1" --type heatmap --field chi_0 \
  --title "chi (generated wave, resonant)" \
  --output "$OUT1/heatmap_chi_0.png"

# Peak amplitude vs time: shows energy flowing from phi to chi
tidal plot "$OUT1" --type amplitude --fields phi_0,chi_0 \
  --title "Amplitude evolution (resonant, mPhi2=mChi2=1)" \
  --output "$OUT1/amplitude.png"

# Profile snapshots showing wave structure at key times
tidal plot "$OUT1" --type profile --field phi_0 \
  --time-indices 0,25,50,75,100 \
  --title "phi profile evolution (resonant)" \
  --output "$OUT1/profile_phi_0.png"

tidal plot "$OUT1" --type profile --field chi_0 \
  --time-indices 0,25,50,75,100 \
  --title "chi profile evolution (resonant)" \
  --output "$OUT1/profile_chi_0.png"

# Hamiltonian energy decomposition (kinetic + gradient + mass + coupling)
tidal plot "$OUT1" --type hamiltonian --fields phi_0,chi_0 \
  --title "Hamiltonian energy (resonant)" \
  --output "$OUT1/hamiltonian.png"

# Energy conservation (should be tight — tests position-dependent Hamiltonian)
tidal plot "$OUT1" --type conservation \
  --title "Energy conservation (resonant)" \
  --output "$OUT1/conservation.png"

# Initial vs final profile comparison
tidal plot "$OUT1" --type compare --fields phi_0,chi_0 \
  --title "Initial vs final (resonant)" \
  --output "$OUT1/compare.png"

# Conversion probability: P(t) = E_chi(t) / E_phi(0)
tidal measure "$OUT1" \
  --what conversion,mixing \
  --source phi_0 --target chi_0

echo ""
echo "=== Run 1 (resonant) complete. Plots: $OUT1/ ==="

# ============================================================================
# Run 2: Non-degenerate (mChi2 = 4.0, mass mismatch)
# ============================================================================
# Unequal masses reduce conversion efficiency. Same Gaussian pulse from
# x=-20 enters the coupling region, but the heavier chi field responds
# more weakly and disperses differently.

OUT2=../data/coupled_scattering_nondegenerate

tidal simulate ../data/coupled_scattering.json \
  --param mPhi2=1.0 --param mChi2=4.0 --param g0=0.5 --param R=5.0 \
  --grid-shape 512 \
  --bounds=-40:40 \
  --periodic \
  --scheme leapfrog \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 80.0 \
  --output "$OUT2"

# --- Non-degenerate plots ---

tidal plot "$OUT2" --type heatmap --field phi_0 \
  --title "phi (incident, non-degenerate)" \
  --output "$OUT2/heatmap_phi_0.png"

tidal plot "$OUT2" --type heatmap --field chi_0 \
  --title "chi (generated, non-degenerate mChi2=4)" \
  --output "$OUT2/heatmap_chi_0.png"

tidal plot "$OUT2" --type amplitude --fields phi_0,chi_0 \
  --title "Amplitude evolution (non-degenerate, mChi2=4)" \
  --output "$OUT2/amplitude.png"

tidal plot "$OUT2" --type conservation \
  --title "Energy conservation (non-degenerate)" \
  --output "$OUT2/conservation.png"

tidal measure "$OUT2" \
  --what conversion,mixing \
  --source phi_0 --target chi_0

echo ""
echo "=== Run 2 (non-degenerate) complete. Plots: $OUT2/ ==="

# ============================================================================
# Run 3: Strong coupling (g0 = 0.9, near stability threshold)
# ============================================================================
# Stronger coupling (g0^2=0.81 vs threshold at mPhi2*mChi2=1.0) produces
# more dramatic chi excitation. Same pulse from x=-20 enters the coupling
# region with much larger energy transfer.

OUT3=../data/coupled_scattering_strong

tidal simulate ../data/coupled_scattering.json \
  --param mPhi2=1.0 --param mChi2=1.0 --param g0=0.9 --param R=5.0 \
  --grid-shape 512 \
  --bounds=-40:40 \
  --periodic \
  --scheme leapfrog \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 80.0 \
  --output "$OUT3"

# --- Strong coupling plots ---

tidal plot "$OUT3" --type heatmap --field phi_0 \
  --title "phi (strong coupling g0=0.9)" \
  --output "$OUT3/heatmap_phi_0.png"

tidal plot "$OUT3" --type heatmap --field chi_0 \
  --title "chi (strong coupling g0=0.9)" \
  --output "$OUT3/heatmap_chi_0.png"

tidal plot "$OUT3" --type amplitude --fields phi_0,chi_0 \
  --title "Amplitude evolution (strong coupling)" \
  --output "$OUT3/amplitude.png"

tidal plot "$OUT3" --type profile --field phi_0 \
  --time-indices 0,25,50,75,100 \
  --title "phi profile (strong coupling)" \
  --output "$OUT3/profile_phi_0.png"

tidal plot "$OUT3" --type profile --field chi_0 \
  --time-indices 0,25,50,75,100 \
  --title "chi excitation (strong coupling)" \
  --output "$OUT3/profile_chi_0.png"

tidal plot "$OUT3" --type conservation \
  --title "Energy conservation (strong coupling)" \
  --output "$OUT3/conservation.png"

tidal measure "$OUT3" \
  --what conversion,mixing \
  --source phi_0 --target chi_0

echo ""
echo "=== Run 3 (strong coupling) complete. Plots: $OUT3/ ==="

echo ""
echo "All runs complete."
echo "  Resonant:        $OUT1/"
echo "  Non-degenerate:  $OUT2/"
echo "  Strong coupling: $OUT3/"
