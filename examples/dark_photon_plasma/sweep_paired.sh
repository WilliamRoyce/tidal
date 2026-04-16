#!/usr/bin/env bash
# Paired Plasma-Only Baseline: Dark Photon Plasma Model with δ_m = 0
#
# Runs the same (mA², mT²) parameter space as sweep_mc.sh but with
# kinetic mixing turned off (deltam=0). This gives the plasma-only
# Gertsenshtein conversion P(h→a) at each point, which serves as the
# baseline for computing the dark-photon amplification factor:
#
#   A_dark = P_full(δ≠0) / P_plasma(δ=0)
#
# A_dark > 1 means the dark photon enhances total conversion.
# A_dark < 1 means the dark photon drains energy from the photon channel.
# A_dark = 1 means the dark photon has no effect (either off-resonance
# or in the Holdom-trivial vacuum limit).
#
# Use with: analyze_sweep.py <full_sweep_dir> --paired-dir <this_dir>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/dark_photon_plasma.json"
OUTPUT="${SCRIPT_DIR}/../data/dark_photon_plasma_paired"

# Number of samples (should match sweep_mc.sh)
N_SAMPLES="${1:-1000}"

# Fixed physics parameters (same as sweep_mc.sh)
KAPPA=1.0
B0=0.01
K0=2.0
XI=1.0

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

# Same mass bounds as sweep_mc.sh, but deltam FIXED at 0
MA2_BOUNDS="0.01:5.0:2"
MT2_BOUNDS="0.01:5.0:2"

echo "=== Plasma Dark Photon: Paired Baseline (δ_m = 0) ==="
echo "Samples: ${N_SAMPLES}"
echo "Fixed: kappa=${KAPPA}, B0=${B0}, xi=${XI}, deltam=0 (no mixing)"
echo "Sweep: mA2 ∈ [0.01, 5], mT2 ∈ [0.01, 5]"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "mA2=${MA2_BOUNDS}" \
  --sweep "mT2=${MT2_BOUNDS}" \
  --sweep-strategy latin_hypercube \
  --n-samples "${N_SAMPLES}" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" --param "xi=${XI}" \
  --param "deltam=0" --param "mT2=1.0" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "=== Paired baseline complete ==="
echo "Results: ${OUTPUT}"
echo ""
echo "Use with: python analyze_sweep.py <full_sweep> --paired-dir ${OUTPUT}"
