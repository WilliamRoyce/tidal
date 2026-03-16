#!/usr/bin/env bash
# Gertsenshtein Scattering Sweep — Wavepacket through Localized B-field
#
# A graviton wavepacket (h_7) starts at x=-80, propagates through a localized
# Gaussian B-field B(x) = Bpeak*exp(-x²/(2R²)) centered at x=0, and partially
# converts to photons (a_2).  Sweeps Bpeak to trace the Boccaletti sin² curve.
#
# Uses periodic BCs with large domain [-200, 200] to enable the modal solver
# (auto-selected).  The B-field decays to float64 zero at boundaries:
#   exp(-200²/50) ≈ 10⁻³⁴⁷ → passes periodic coefficient continuity check.
#
# The modal solver uses expm_multiply (Al-Mohy & Higham 2011) for the full
# convolution matrix with position-dependent coefficients.  This is 5.7x faster
# than CVODE and handles the non-normal gradient coupling correctly.
#
# NOTE: IC noise replicates (--ic-perturbation) are NOT used here because
# the modal solver's expm_multiply path has genuinely growing modes for
# the a_0 constraint field in Fourier space.  Noise on a_0 projects onto
# these modes and causes overflow.  The clean IC (only h_7 excited) avoids
# this by having zero projection onto unstable subspace.
# For ensemble UQ, use --scheme cvode instead of modal.
#
# Analytical target (Boccaletti 1970, massless vacuum Gertsenshtein):
#   P = sin²(kappa * Bpeak * R * sqrt(pi/2))
#     = sin²(kappa/2 * integral(B(x) dx, -inf, +inf))
# This formula is EXACT for massless fields (Delta_k = 0) regardless of
# B-field profile shape.  Ref: Boccaletti et al. (1970, Nuovo Cimento 70B:129).
#
# Parameter range:
#   Bpeak in [0.01, 0.30] (20 points)
#   For kappa=1, R=5: argument ranges from 0.063 (P≈0.004) to 1.88 (P≈0.95)
#   This covers the full weak-field → near-maximum conversion regime.
#
# NOTE: constant "Bpeak" (not "B0_peak") — underscores are Mathematica patterns.
# NOTE: negative --bounds and --ic-center require = syntax: --bounds=-200:200
#
# Ref: Boccaletti et al. (1970, Nuovo Cimento 70B:129)
#      Al-Mohy & Higham (2011, SIAM J. Sci. Comput. 33(2):488)
#      Trefethen & Embree (2005, Spectra and Pseudospectra, Ch. 14)
#
# Running:
#   cd examples/gertsenshtein && bash sweep_scattering.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein Scattering Sweep ==="
echo "Wavepacket h_7 through Gaussian B-field, sweeping Bpeak"
echo ""

tidal sweep ../data/gertsenshtein_localized.json \
  --sweep "Bpeak=0.01:0.30:20" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  "--bounds=-200:200" \
  --periodic \
  --ic gaussian \
  --ic-wavevector 2.0 \
  --ic-amplitude 0.1 \
  --ic-width 5.0 \
  "--ic-center=-80.0" \
  --ic-component h_7 \
  --t-end 120.0 \
  --param kappa=1.0 --param R=5.0 \
  --output ../data/gertsenshtein_sweep_scattering

echo ""

# Plot with error bands + Boccaletti analytical overlay
# 3-panel figure: TIDAL (mean ± σ) | analytical sin² | |error|
# Ref: Boccaletti et al. (1970, Nuovo Cimento 70B:129)
tidal plot ../data/gertsenshtein_sweep_scattering --type sweep \
  --metric P_final \
  --title "Gertsenshtein scattering: P(Bpeak) — modal solver (expm_multiply)" \
  --overlay 'sin(kappa * Bpeak * R * sqrt(pi/2))**2' \
  --output ../data/gertsenshtein_sweep_scattering/sweep.png --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: examples/data/gertsenshtein_sweep_scattering/"
echo "Plot:    examples/data/gertsenshtein_sweep_scattering/sweep.png"
