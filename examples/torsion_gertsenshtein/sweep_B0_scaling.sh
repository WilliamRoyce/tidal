#!/usr/bin/env bash
# B₀ scaling verification for torsion-mediated Gertsenshtein amplification
#
# Goal: Verify C₀ = P/B₀² is constant at the amplification peak,
# confirming that the effect is a coupling modification (linear in B₀²),
# not a nonlinear backreaction effect.
#
# Fixed parameters: delta1=1.0, alpha2=-0.82 (amplification peak),
#   alpha1=0, alpha3=1.0, kappa=1
# Swept: B0 from 1e-6 to 1e-2 (log scale, 9 points)
#
# Expected: C₀ = P/B₀² should be approximately constant across all B₀.
# Deviations at large B₀ indicate nonlinear regime (P > 0.1).

set -euo pipefail

SPEC="examples/data/torsion_gertsenshtein_nonminimal.json"
OUTPUT="examples/data/nonminimal_B0_scaling"

uv run tidal sweep "$SPEC" \
  --sweep "B0=1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2" \
  --measure peak_conversion --source h_5 --target a_1 \
  --grid-shape 256 --bounds "0:100" --periodic --t-end 50 \
  --param delta1=1.0 --param alpha2=-0.82 --param alpha1=0 --param alpha3=1.0 --param kappa=1 \
  --ic plane-wave --ic-amplitude 0.1 --ic-component h_5 \
  --output "$OUTPUT" --parallel 4 --resume

echo ""
echo "Postprocessing..."
echo "B0,P_max,C0" > "${OUTPUT}/B0_scaling.csv"

# Extract P_max for each B0 and compute C0
python3 -c "
import csv, math
with open('${OUTPUT}/results.csv') as f:
    rows = list(csv.DictReader(f))
print(f'B0 scaling results ({len(rows)} points):')
print(f'{\"B0\":>10} {\"P_max\":>12} {\"C0=P/B0^2\":>14} {\"status\":>8}')
for r in rows:
    b0 = float(r['B0'])
    status = r.get('run_status', 'success')
    try:
        pmax = float(r.get('P_max', ''))
        c0 = pmax / b0**2
        print(f'{b0:10.1e} {pmax:12.4e} {c0:14.2f} {status:>8}')
    except (ValueError, TypeError):
        print(f'{b0:10.1e} {\"N/A\":>12} {\"N/A\":>14} {status:>8}')
"
