#!/usr/bin/env bash
# v3 D2.0 Bahamonde amp — sub-theory of T5 with β,δ free, ξ,χ,ζ pruned
#
# v3 architecture: see docs/V3_ARCHITECTURE.md.
# Sub-theory partition retained per docs/lagrangian_depruning_audit.md:
# Bahamonde et al.'s paper studies propagating-torsion-mass-only (no kinetic
# trace, no derivative-torsion-EM cross-terms).  Sub-theory = literature
# comparison choice, not analytical-inertness pruning, so partition stays.
#
# Reference v2: 28598736 (log Z = +0.616 ± 0.001, ESS=877).

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d20_bahamonde_amp_v3 --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --param xi=0 --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 400 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d20_bahamonde_amp_v3'
