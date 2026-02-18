#!/bin/bash
# Run all example derivations via `tidal derive` to regenerate JSON files
#
# Usage: ./scripts/run_examples.sh
#
# Derives equations from every example TOML config.
# Useful after modifying Wolfram pipeline modules.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Running Example Derivations ==="
echo ""

EXAMPLES=(
    "examples/scalar_field/theory.toml"
    "examples/electromagnetic/theory.toml"
    "examples/coupled_scalars/theory.toml"
    "examples/chern_simons/theory.toml"
    "examples/elasticity/theory.toml"
    "examples/proca/theory.toml"
    "examples/curved_spacetime/de_sitter.toml"
    "examples/curved_spacetime/conformal_static.toml"
    "examples/sphere_kg/theory.toml"
    "examples/polar_kg/theory.toml"
    "examples/spherical_kg/theory.toml"
    "examples/cylindrical_kg/theory.toml"
    "examples/scalar_field_3d/theory.toml"
    "examples/scalar_vector_coupling/theory.toml"
    "examples/massive_3form/theory.toml"
    "examples/massive_gravity/theory.toml"
    "examples/coupled_proca/theory.toml"
    "examples/electrostatics/theory.toml"
    "examples/gravitational_waves/theory.toml"
    "examples/coupled_scattering/theory.toml"
    "examples/scalar_potential_well/theory.toml"
    "examples/vector_background/theory.toml"
    "examples/proca_background/theory.toml"
)

PASSED=0
FAILED=0
FAILED_EXAMPLES=()

for toml in "${EXAMPLES[@]}"; do
    TOML_PATH="$PROJECT_ROOT/$toml"

    if [ ! -f "$TOML_PATH" ]; then
        echo "SKIP: $toml (file not found)"
        continue
    fi

    echo "Running: tidal derive $toml"
    if (cd "$PROJECT_ROOT" && tidal derive "$toml"); then
        PASSED=$((PASSED + 1))
        echo "  OK"
    else
        FAILED=$((FAILED + 1))
        FAILED_EXAMPLES+=("$toml")
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
