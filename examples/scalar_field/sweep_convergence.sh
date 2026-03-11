#!/usr/bin/env bash
# Grid Convergence Study — Klein-Gordon Scalar Field
#
# Runs convergence studies on the simplest TIDAL example (free scalar field)
# to verify spatial convergence order at different FD accuracy levels.
#
# Three sweeps demonstrate the --fd-order feature:
#   1. fd_order=2 (3-point stencil): error decreases as O(dx^2)
#   2. fd_order=4 (5-point stencil, default): error decreases as O(dx^4)
#   3. fd_order=6 (7-point stencil): error decreases as O(dx^6)
#
# The log-log slope of the convergence plot confirms the theoretical rate.
# At N=128, order 4 achieves the same accuracy as order 2 at N=512.
#
# Ref: Fornberg (1988), Mathematics of Computation 51(184), pp. 699-706.
#
# Running this script:
#   cd examples/scalar_field && bash sweep_convergence.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Klein-Gordon Convergence Study ==="
echo ""

# Study 1: 2nd-order FD — O(dx^2) convergence
echo "--- FD order 2 (3-point stencil) ---"
tidal sweep ../data/klein_gordon_1d.json \
  --converge "32,64,128,256" \
  --param m2=1.0 \
  --measure conservation \
  --bounds 0:50 \
  --periodic \
  --ic gaussian --ic-center 25.0 --ic-width 3.0 \
  --t-end 20.0 \
  --fd-order 2 \
  --output ../data/sweep_convergence_kg_o2

# Study 2: 4th-order FD (default) — O(dx^4) convergence
echo "--- FD order 4 (5-point stencil, default) ---"
tidal sweep ../data/klein_gordon_1d.json \
  --converge "32,64,128,256" \
  --param m2=1.0 \
  --measure conservation \
  --bounds 0:50 \
  --periodic \
  --ic gaussian --ic-center 25.0 --ic-width 3.0 \
  --t-end 20.0 \
  --output ../data/sweep_convergence_kg_o4

# Study 3: 6th-order FD — O(dx^6) convergence
echo "--- FD order 6 (7-point stencil) ---"
tidal sweep ../data/klein_gordon_1d.json \
  --converge "32,64,128,256" \
  --param m2=1.0 \
  --measure conservation \
  --bounds 0:50 \
  --periodic \
  --ic gaussian --ic-center 25.0 --ic-width 3.0 \
  --t-end 20.0 \
  --fd-order 6 \
  --output ../data/sweep_convergence_kg_o6

echo ""
echo "Results written to:"
echo "  ../data/sweep_convergence_kg_o2/  (FD order 2)"
echo "  ../data/sweep_convergence_kg_o4/  (FD order 4, default)"
echo "  ../data/sweep_convergence_kg_o6/  (FD order 6)"
echo ""
echo "Plot the convergence:"
echo "  tidal plot ../data/sweep_convergence_kg_o2/ --type convergence --metric max_energy_error"
echo "  tidal plot ../data/sweep_convergence_kg_o4/ --type convergence --metric max_energy_error"
echo "  tidal plot ../data/sweep_convergence_kg_o6/ --type convergence --metric max_energy_error"
