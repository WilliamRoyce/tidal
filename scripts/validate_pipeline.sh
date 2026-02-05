#!/bin/bash
# End-to-end pipeline validation: Lagrangian -> JSON -> PDE simulation
#
# Usage: ./scripts/validate_pipeline.sh
#
# Runs a complete pipeline validation:
# 1. Derives equations from Lagrangian (Wolfram)
# 2. Validates JSON was created
# 3. Runs Python pipeline tests

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

echo "=== Pipeline Validation ==="
echo ""

# Step 1: Run a simple derivation
echo "Step 1: Deriving Klein-Gordon equations..."
echo ""
wolframscript -file "$PROJECT_ROOT/examples/scalar_field/klein_gordon.wls"

# Step 2: Validate JSON was created
JSON_FILE="$PROJECT_ROOT/examples/data/klein_gordon_1d.json"
echo ""
echo "Step 2: Validating JSON output..."
if [ ! -f "$JSON_FILE" ]; then
    echo "ERROR: JSON file not created: $JSON_FILE"
    exit 1
fi
echo "  JSON created: $JSON_FILE"

# Validate JSON is valid
if ! python3 -c "import json; json.load(open('$JSON_FILE'))" 2>/dev/null; then
    echo "ERROR: JSON file is not valid JSON"
    exit 1
fi
echo "  JSON is valid"

# Step 3: Run Python pipeline tests that use this JSON
echo ""
echo "Step 3: Running Python pipeline tests..."
uv run pytest tests/test_em_pipeline.py -v

echo ""
echo "=== Pipeline validation passed ==="
exit 0
