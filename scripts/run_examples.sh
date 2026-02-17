#!/bin/bash
# Run all Wolfram example derivations to regenerate JSON files
#
# Usage: ./scripts/run_examples.sh
#
# Runs each example script in examples/ to regenerate equations JSON.
# Useful after modifying Wolfram pipeline modules.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Running Example Derivations ==="
echo ""

EXAMPLES=(
    "examples/scalar_field/klein_gordon.wls"
    "examples/electromagnetic/em_lagrangian_1d.wls"
    "examples/coupled_scalars/coupled_scalars.wls"
    "examples/chern_simons/chern_simons.wls"
    "examples/elasticity/navier_cauchy.wls"
    "examples/proca/proca.wls"
    "examples/curved_spacetime/de_sitter_kg.wls"
    "examples/curved_spacetime/conformal_kg_static.wls"
    "examples/sphere_kg/sphere_kg.wls"
    "examples/polar_kg/polar_kg.wls"
    "examples/spherical_kg/spherical_kg.wls"
    "examples/cylindrical_kg/cylindrical_kg.wls"
    "examples/scalar_field_3d/klein_gordon_3d.wls"
    "examples/scalar_vector_coupling/scalar_vector_coupling.wls"
    "examples/massive_3form/massive_3form.wls"
    "examples/massive_gravity/massive_gravity.wls"
    "examples/coupled_proca/coupled_proca.wls"
)

PASSED=0
FAILED=0
FAILED_EXAMPLES=()

for example in "${EXAMPLES[@]}"; do
    EXAMPLE_PATH="$PROJECT_ROOT/$example"

    if [ ! -f "$EXAMPLE_PATH" ]; then
        echo "SKIP: $example (file not found)"
        continue
    fi

    echo "Running: $example"
    if wolframscript -file "$EXAMPLE_PATH"; then
        PASSED=$((PASSED + 1))
        echo "  OK"
    else
        FAILED=$((FAILED + 1))
        FAILED_EXAMPLES+=("$example")
        echo "  FAILED"
    fi
    echo ""
done

echo "=== Example Derivation Summary ==="
echo "Passed: $PASSED"
echo "Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Failed examples:"
    for e in "${FAILED_EXAMPLES[@]}"; do
        echo "  - $e"
    done
    exit 1
fi

echo ""
echo "=== All examples completed successfully ==="
exit 0
