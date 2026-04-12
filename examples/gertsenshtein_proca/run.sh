#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — Validation Run
#
# Tests the basic functionality of the Proca/plasma Gertsenshtein theory:
# 1. mA2=0: must recover standard Gertsenshtein P = sin^2(kappa*B0*D/2)
# 2. mA2>0: detuned conversion (lower P than vacuum)
#
# Running:
#   cd examples/gertsenshtein_proca && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json

echo "=== Gertsenshtein with Massive Photon: Validation ==="
echo ""

# --- Test 1: mA2=0 (massless vacuum — standard Gertsenshtein) ---
echo "--- Test 1: mA2=0 (should recover vacuum Gertsenshtein) ---"
OUT1=/tmp/tidal_tests/gert_proca_massless

tidal simulate "$JSON" \
  --param kappa=1.0 --param B0=0.1 --param mA2=0.0 \
  --grid-shape 512 --bounds 0:100 --periodic \
  --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT1"

echo ""
tidal measure "$OUT1" --what conversion --source h_7 --target a_2 \
  --param kappa=1.0 --param B0=0.1 --param mA2=0.0

echo ""
echo "Expected: P ≈ sin²(1.0 * 0.1 * 50.0 / 2) = sin²(2.5) ≈ 0.358"
echo ""

# --- Test 2: mA2=0.5 (detuned — suppressed conversion) ---
echo "--- Test 2: mA2=0.5 (detuned — should be suppressed) ---"
OUT2=/tmp/tidal_tests/gert_proca_detuned

tidal simulate "$JSON" \
  --param kappa=1.0 --param B0=0.1 --param mA2=0.5 \
  --grid-shape 512 --bounds 0:100 --periodic \
  --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT2"

echo ""
tidal measure "$OUT2" --what conversion --source h_7 --target a_2 \
  --param kappa=1.0 --param B0=0.1 --param mA2=0.5

echo ""
# Raffelt-Stodolsky: Delta = -0.5/(2*2.0106) ≈ -0.124
# mu = kappa*B0/2 = 0.05
# Delta_osc = sqrt(0.124^2 + 0.05^2) ≈ 0.134
# sin^2(2*theta) = (0.05/0.134)^2 ≈ 0.139
# P = 0.139 * sin^2(0.134*50/2) = 0.139 * sin^2(3.35) ≈ 0.139 * 0.035 ≈ 0.005
echo "Expected: P ≈ 0.005 (strongly suppressed by mass detuning)"
echo ""

# --- Test 3: Conservation check ---
echo "--- Conservation check (mA2=0) ---"
tidal measure "$OUT1" --what conservation \
  --param kappa=1.0 --param B0=0.1 --param mA2=0.0

echo ""
echo "Done. Results in: /tmp/tidal_tests/gert_proca_*"
