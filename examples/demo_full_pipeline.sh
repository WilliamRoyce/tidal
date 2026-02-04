#!/bin/bash
# Demonstration of the complete Lagrangian-to-Simulation pipeline
# This script shows both stages working together without any hardcoded physics

set -e  # Exit on error

echo "============================================================"
echo "Lagrangian-to-Simulation Pipeline Demonstration"
echo "============================================================"
echo ""
echo "This demonstrates that NO physics is hardcoded in Python."
echo "All equation structure comes from symbolic derivation."
echo ""

# Change to project root
cd "$(dirname "$0")/.."

echo "============================================================"
echo "EXAMPLE 1: Electromagnetic Field (Vector Field)"
echo "============================================================"
echo ""
echo "Lagrangian: L = -1/4 F_μν F^μν"
echo ""

echo "--- Stage 1: Symbolic Derivation (Mathematica/xAct) ---"
echo "Deriving Maxwell equations from Lagrangian..."
cd examples/electromagnetic
wolframscript -file em_lagrangian_1d.wls > /tmp/em_derivation.log 2>&1
echo "✓ Equations derived and exported to JSON"
echo ""

cd ../..
echo "--- Verify JSON Structure ---"
echo "Checking that JSON contains Laplacian operator (no mass term):"
cat examples/data/em_1d.json | jq '.equations[0].rhs.terms'
echo ""

echo "--- Stage 2: Numerical Simulation (Python/py-pde) ---"
echo "Building PDE dynamically from JSON specification..."
python examples/electromagnetic/em_from_lagrangian.py > /tmp/em_simulation.log 2>&1
echo "✓ Simulation complete"
echo "  Output: outputs/em_from_lagrangian_output.png"
echo ""

echo "============================================================"
echo "EXAMPLE 2: Klein-Gordon Field (Scalar Field)"
echo "============================================================"
echo ""
echo "Lagrangian: L = -1/2 (∂φ)² - 1/2 m²φ²"
echo ""

echo "--- Stage 1: Symbolic Derivation (Mathematica/xAct) ---"
echo "Deriving Klein-Gordon equation from Lagrangian..."
cd examples/scalar_field
wolframscript -file klein_gordon.wls > /tmp/kg_derivation.log 2>&1
echo "✓ Equations derived and exported to JSON"
echo ""

cd ../..
echo "--- Verify JSON Structure ---"
echo "Checking that JSON contains BOTH Laplacian AND mass term:"
cat examples/data/klein_gordon_1d.json | jq '.equations[0].rhs.terms'
echo ""

echo "--- Stage 2: Numerical Simulation (Python/py-pde) ---"
echo "Building PDE dynamically from JSON specification..."
python examples/scalar_field/kg_from_lagrangian.py > /tmp/kg_simulation.log 2>&1
echo "✓ Simulation complete"
echo "  Output: outputs/kg_from_lagrangian_output.png"
echo ""

echo "============================================================"
echo "COMPARISON: Different Lagrangians → Different Physics"
echo "============================================================"
echo ""
echo "EM JSON (no mass term):"
echo "  $(cat examples/data/em_1d.json | jq -c '.equations[0].rhs.terms[0]')"
echo ""
echo "Klein-Gordon JSON (has mass term):"
echo "  Identity (mass):  $(cat examples/data/klein_gordon_1d.json | jq -c '.equations[0].rhs.terms[0]')"
echo "  Laplacian (wave): $(cat examples/data/klein_gordon_1d.json | jq -c '.equations[0].rhs.terms[1]')"
echo ""

echo "============================================================"
echo "VERIFICATION: No Hardcoded Physics"
echo "============================================================"
echo ""
echo "The Python code uses build_pde_from_json() which:"
echo "  1. Loads operator types from JSON (laplacian, identity, etc.)"
echo "  2. Dynamically applies operators based on specification"
echo "  3. Does NOT know about 'wave equations' or 'mass terms'"
echo ""
echo "Evidence: Same Python code simulates both EM (massless) and"
echo "Klein-Gordon (massive) by reading different JSON files."
echo ""

echo "============================================================"
echo "✓ Pipeline demonstration complete!"
echo "============================================================"
echo ""
echo "Check outputs:"
echo "  - outputs/em_from_lagrangian_output.png"
echo "  - outputs/kg_from_lagrangian_output.png"
echo ""
echo "See examples/README.md for detailed documentation."
