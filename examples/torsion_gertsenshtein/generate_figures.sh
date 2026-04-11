#!/usr/bin/env bash
# Generate all 7 publication figures for the torsion-Gertsenshtein campaign.
#
# Figures 1-6 use the tidal CLI plotting infrastructure.
# Figure 7 uses a standalone script (raw NPY time-series, no CLI equivalent).
#
# Usage:
#   bash generate_figures.sh [DATA_DIR] [OUTPUT_DIR]
#
# Defaults:
#   DATA_DIR   = /tmp/tidal_sweeps_v2
#   OUTPUT_DIR = DATA_DIR/figures

set -euo pipefail

DATA_DIR="${1:-/tmp/tidal_sweeps_v2}"
OUT_DIR="${2:-${DATA_DIR}/figures}"
mkdir -p "$OUT_DIR"

echo "Generating publication figures..."
echo "  Data:   $DATA_DIR"
echo "  Output: $OUT_DIR"
echo

# --- Fig 1: Phase 1 overview (3x3 grid, which parameters matter) ---
echo "Figure 1: Phase 1 Overview"
uv run tidal plot "$DATA_DIR/phase1" \
  --type campaign --metric log10_A --layout 3x3 \
  --scatter --log-y --classify "active:log10_A>0.01" \
  --title "Phase 1: Which Parameters Affect Conversion?" \
  --dpi 200 --output "$OUT_DIR/fig1_phase1_overview.png"

# --- Fig 2: Active parameters detail (delta1, zeta2, zeta3) ---
echo "Figure 2: Active Parameters Detail"
uv run tidal plot "$DATA_DIR/phase1" \
  --type campaign --metric log10_A --layout 1x3 \
  --scatter --annotate-extremes --log-y \
  --title "Active Coupling Parameters: Conversion Modification" \
  --dpi 200 --output "$OUT_DIR/fig2_active_detail.png"

# --- Fig 3: Phase 2 pairwise interactions ---
echo "Figure 3: Pairwise Interactions"
uv run tidal plot "$DATA_DIR/phase2" \
  --type campaign --metric log10_A --layout 2x3 \
  --title "Phase 2: Pairwise Parameter Interactions" \
  --dpi 200 --output "$OUT_DIR/fig3_pairwise.png"

# --- Fig 4: Chi toxicity (divergence rates) ---
echo "Figure 4: Chi Toxicity"
# Use the Phase 2 LHS data for per-pair divergence analysis
uv run tidal plot "$DATA_DIR/phase2" \
  --type campaign --metric log10_A --layout 3x4 \
  --title "Phase 2: Divergence Rates by Parameter Pair" \
  --dpi 200 --output "$OUT_DIR/fig4_chi_toxicity.png"

# --- Fig 5a: Phase 3 tornado (sensitivity) ---
echo "Figure 5a: Tornado (Sensitivity)"
uv run tidal plot "$DATA_DIR/phase3/lhs1000" \
  --type sweep-tornado --metric log10_A \
  --title "Phase 3: Parameter Sensitivity (1000 LHS Samples)" \
  --dpi 200 --output "$OUT_DIR/fig5a_tornado.png"

# --- Fig 5b: Phase 3 histogram (distribution) ---
echo "Figure 5b: Histogram (Distribution)"
uv run tidal plot "$DATA_DIR/phase3/lhs1000" \
  --type sweep-histogram --metric log10_A \
  --title "Phase 3: Amplification Distribution" \
  --dpi 200 --output "$OUT_DIR/fig5b_distribution.png"

# --- Fig 6: Phase 3 scatter matrix (active params) ---
echo "Figure 6: Scatter Matrix"
uv run tidal plot "$DATA_DIR/phase3/lhs1000" \
  --type sweep-scatter --metric log10_A \
  --params delta1,zeta2,zeta3 --log-diagonal \
  --divergent-center 0 --cmap RdYlBu_r \
  --title "Phase 3: Active Parameter Space" \
  --dpi 200 --output "$OUT_DIR/fig6_scatter_matrix.png"

# --- Fig 7: Frequency modification (standalone, needs raw NPY data) ---
echo "Figure 7: Frequency Modification"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$SCRIPT_DIR/plot_frequency_modification.py" \
  "/tmp/tidal_tests" "$OUT_DIR/fig7_frequency_modification.png" \
  || echo "  (Skipped: simulation time-series data not found)"

echo
echo "All figures saved to: $OUT_DIR"
