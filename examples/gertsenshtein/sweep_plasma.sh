#!/usr/bin/env bash
# Gertsenshtein Plasma Detuning Sweep — Phase F1 Validation
#
# Sweeps the effective photon mass (plasma frequency squared, omegaP2)
# at fixed B0 and measures conversion probability P(graviton -> photon).
#
# Analytical target (Raffelt & Stodolsky 1988, PRD 37:1237):
#   P = (mu/omega_m)^2 * sin^2(omega_m * D)
#   where:
#     mu = kappa * B0 / 2 = 0.05        (coupling strength, kappa=1, B0=0.1)
#     Delta = omegaP2 / (2 * omega)      (detuning, omega ~ k ~ 2.011)
#     omega_m = sqrt(Delta^2 + mu^2)     (mixing frequency)
#
# At omegaP2 = 0 (vacuum): P = sin^2(0.05 * 50) = sin^2(2.5)
# As omegaP2 increases: P suppressed as (mu/Delta)^2
# MSW-like resonance: not applicable for massless graviton (m_g = 0),
#   but the detuning curve maps the full transition from resonant to
#   off-resonant regimes.
#
# IC: monochromatic plane wave cos(k*z) on h_7 (graviton h_+)
# k = 2*pi*32/100 ~ 2.011 (32 wavelengths in domain [0,100])
#
# Running:
#   cd examples/gertsenshtein && bash sweep_plasma.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein Plasma Detuning Sweep (Phase F1) ==="
echo ""

tidal sweep ../data/gertsenshtein_plasma.json \
  --sweep "omegaP2=0.0:1.0:30" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-wavevector 2.0106 \
  --ic-amplitude 0.1 \
  --ic-component h_7 \
  --t-end 50.0 \
  --param kappa=1.0 \
  --param B0=0.1 \
  --output ../data/gertsenshtein_sweep_plasma

echo ""

# Plot sweep results: P(omegaP2) vs analytical detuning curve
tidal plot ../data/gertsenshtein_sweep_plasma --type sweep \
  --metric P_final \
  --title "Gertsenshtein P(omegaP2) — plasma detuning" \
  --output ../data/gertsenshtein_sweep_plasma/sweep.png --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: examples/data/gertsenshtein_sweep_plasma/"
echo "Plot:   examples/data/gertsenshtein_sweep_plasma/sweep.png"
