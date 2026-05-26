#!/usr/bin/env bash
# Phase E overnight LANDSCAPE-quality run for ALL Wave 1+2+3 theories in parallel.
#
# Runs inside one standard-queue sapphire allocation. Each theory's amp+sup
# chains run in parallel as independent mpirun calls, sharing the node pool.
# Designed for ~12-14h walltime on 7 sapphire nodes (~784 cores).
#
# Theory roster (8 theories × amp+sup = 16 chains):
#   Small (ndim≤6, 32 ranks each):
#     - E.T4 (Ricci-EM nonminimal, 4D)
#     - E.T1 (Dark-Photon-Plasma, 4D)
#     - E.T5.0 (Bahamonde 5D — restricted from general_nonminimal)
#     - E.T5.1 (Barker 6D — restricted from general_nonminimal)
#   Medium (ndim=8-9, 32 ranks each):
#     - E.NP (non-prop torsion, 8D — ξ pinned to 0 on general_nonminimal)
#     - E.T5.2 (Shapiro 8D — restricted)
#     - E.T5.3 (Complete 9D)
#   Large (ndim=18, 128 ranks):
#     - E.T7s (complete-even 18D — full chi sector)
#
# Empirical rate (post-#379, pre-#384 caching): ~4 ndead/min at 32 ranks.
# Convergence at LANDSCAPE (nlive=25·ndim): ~10·nlive ndead, so:
#   ndim=4: ~4 hr, ndim=8: ~8 hr, ndim=9: ~10 hr, ndim=18 @ 128 ranks: ~12 hr
#
# Resource budget per chain:
#   14 chains × 32 ranks = 448 cores
#   2 E.T7s chains × 128 ranks = 256 cores
#   Total: 704 cores ≈ 7 sapphire nodes (784 cores)

set -euo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"
GEOMETRY_ENV="${TIDAL_ROOT}/scripts/hpc_submit_drafts/v3e_localised/_geometry.env"

# shellcheck disable=SC1090
source "${GEOMETRY_ENV}"

cd "${TIDAL_ROOT}"

# RESULTS_DIR is set by the sapphire sbatch template wrapper.
: "${RESULTS_DIR:?RESULTS_DIR must be set by sbatch template}"

# Shared physics flags (LHS gauge etc. baked into the JSONs).
COMMON_PHYS="\
  --baseline-formula sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2 \
  --soft-floor-noise 1.0 \
  --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
  --grid-shape ${N} --bounds=0:${L} --periodic \
  --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
  --t-end ${T_END} --snapshots ${SNAPSHOTS} \
  --measure conversion,peak_conversion --source h_5 --target a_1 \
  --method nested --sampler polychord --precision-criterion 0.01"

declare -A PIDS=()

launch_chain() {
  local tag="$1"             # e.g. "t4_amp"
  local data="$2"            # JSON path
  local ranks="$3"
  local nlive="$4"
  local nrep="$5"
  local likelihood="$6"      # "P_max:maximize" or "P_max:minimize"
  local priors="$7"          # all --prior flags as one string
  local params_extra="${8:-}"  # extra --param flags (for pinning, e.g. xi=0)

  local outdir="${RESULTS_DIR}/e_${tag}"
  local logfile="${RESULTS_DIR}/e_${tag}.log"

  echo "=== [$(date +%H:%M:%S)] Launching e_${tag}: ${ranks} ranks, nlive=${nlive}, num_repeats=${nrep} ==="
  # Use srun (not mpirun) for SLURM-aware rank/core carving across multi-node
  # allocation — multiple parallel chains within one sbatch require srun's
  # --exclusive resource carving to avoid core/memory oversubscription.
  # --mpi=pmi2 chosen for intel-oneapi-mpi compatibility on CSD3.
  # shellcheck disable=SC2086
  srun --exclusive --ntasks="${ranks}" --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
    tidal sample "${data}" \
    ${priors} \
    --likelihood "${likelihood}" \
    ${COMMON_PHYS} \
    ${params_extra} \
    --nlive "${nlive}" --num-repeats "${nrep}" \
    --output "${outdir}" \
    > "${logfile}" 2>&1 &
  PIDS["${tag}"]=$!
}

# ─── Theory priors (BSM couplings only — geometry pinned via --param) ───────
# Wave 1
PRIORS_T4="--prior alpha1=arctan_uniform:-89:89 --prior alpha2=arctan_uniform:-89:89 \
  --prior alpha3=arctan_uniform:-89:89 --prior delta1=arctan_uniform:-89:89"

PRIORS_T1="--prior alpha3=log_uniform:1e-3:1e3 --prior xi=log_uniform:1e-3:1e3 \
  --prior deltam=arctan_uniform:-89:89 --prior mA2=log_uniform:1e-3:1e3"

# Wave 2 (Bahamonde family — all share general_nonminimal JSON)
PRIORS_T50="--prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
  --prior beta3=arctan_uniform:-89:89 --prior xi=log_uniform:1e-3:1e3 \
  --prior delta1=arctan_uniform:-89:89"
PARAMS_T50="--param chi=0.0 --param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0"

PRIORS_T51="${PRIORS_T50} --prior chi=arctan_uniform:-89:89"
PARAMS_T51="--param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0"

PRIORS_T52="${PRIORS_T50} --prior zeta1=arctan_uniform:-89:89 \
  --prior zeta2=arctan_uniform:-89:89 --prior zeta3=arctan_uniform:-89:89"
PARAMS_T52="--param chi=0.0"

PRIORS_T53="${PRIORS_T50} --prior chi=arctan_uniform:-89:89 \
  --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 \
  --prior zeta3=arctan_uniform:-89:89"
PARAMS_T53=""

PRIORS_NP="--prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
  --prior beta3=arctan_uniform:-89:89 --prior delta1=arctan_uniform:-89:89 \
  --prior chi=arctan_uniform:-89:89 \
  --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 \
  --prior zeta3=arctan_uniform:-89:89"
PARAMS_NP="--param xi=0.0"

# Wave 3 — E.T7s priors (chi1..chi10 fan-out)
PRIORS_T7S="${PRIORS_T50} \
  --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 --prior zeta3=arctan_uniform:-89:89 \
  --prior chi1=arctan_uniform:-89:89 --prior chi2=arctan_uniform:-89:89 --prior chi3=arctan_uniform:-89:89 \
  --prior chi4=arctan_uniform:-89:89 --prior chi5=arctan_uniform:-89:89 --prior chi6=arctan_uniform:-89:89 \
  --prior chi7=arctan_uniform:-89:89 --prior chi8=arctan_uniform:-89:89 --prior chi9=arctan_uniform:-89:89 \
  --prior chi10=arctan_uniform:-89:89"
PARAMS_T7S=""

# ─── Launch all 16 chains in parallel ───────────────────────────────────────
echo "=== [$(date +%H:%M:%S)] OVERNIGHT BATCH START — 16 chains in parallel ==="

JSON_T4="examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"
JSON_T1="examples/data/dark_photon_plasma_e_dual_gaussian.json"
JSON_T5="examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json"
JSON_T7S="examples/data/torsion_gertsenshtein_complete_even_e_dual_gaussian.json"

launch_chain t4_amp     "${JSON_T4}"  32  100  8  "P_max:maximize" "${PRIORS_T4}"  ""
launch_chain t4_sup     "${JSON_T4}"  32  100  8  "P_max:minimize" "${PRIORS_T4}"  ""
launch_chain t1_amp     "${JSON_T1}"  32  100  8  "P_max:maximize" "${PRIORS_T1}"  ""
launch_chain t1_sup     "${JSON_T1}"  32  100  8  "P_max:minimize" "${PRIORS_T1}"  ""
launch_chain t5_amp     "${JSON_T5}"  32  125 10  "P_max:maximize" "${PRIORS_T50}" "${PARAMS_T50}"
launch_chain t5_sup     "${JSON_T5}"  32  125 10  "P_max:minimize" "${PRIORS_T50}" "${PARAMS_T50}"
launch_chain t51_amp    "${JSON_T5}"  32  150 12  "P_max:maximize" "${PRIORS_T51}" "${PARAMS_T51}"
launch_chain t51_sup    "${JSON_T5}"  32  150 12  "P_max:minimize" "${PRIORS_T51}" "${PARAMS_T51}"
launch_chain np_amp     "${JSON_T5}"  32  200 16  "P_max:maximize" "${PRIORS_NP}"  "${PARAMS_NP}"
launch_chain np_sup     "${JSON_T5}"  32  200 16  "P_max:minimize" "${PRIORS_NP}"  "${PARAMS_NP}"
launch_chain t52_amp    "${JSON_T5}"  32  200 16  "P_max:maximize" "${PRIORS_T52}" "${PARAMS_T52}"
launch_chain t52_sup    "${JSON_T5}"  32  200 16  "P_max:minimize" "${PRIORS_T52}" "${PARAMS_T52}"
launch_chain t53_amp    "${JSON_T5}"  32  225 18  "P_max:maximize" "${PRIORS_T53}" "${PARAMS_T53}"
launch_chain t53_sup    "${JSON_T5}"  32  225 18  "P_max:minimize" "${PRIORS_T53}" "${PARAMS_T53}"
# E.T7s SKIPPED in this 4-node overnight (QOS cpu=448 limit). Run separately
# via a dedicated job: budget 128 ranks × 2 chains × 18h would need cpul QOS
# (cpu=512, 7-day) — see scripts/hpc_submit_drafts/v3e_localised/run_overnight_t7s.sh
# launch_chain t7s_amp    "${JSON_T7S}" 128 450 36  "P_max:maximize" "${PRIORS_T7S}" "${PARAMS_T7S}"
# launch_chain t7s_sup    "${JSON_T7S}" 128 450 36  "P_max:minimize" "${PRIORS_T7S}" "${PARAMS_T7S}"

echo "=== [$(date +%H:%M:%S)] All 16 chains launched — waiting ==="
echo "PIDs: ${PIDS[*]}"

# Wait for all; record per-chain exit status.
for tag in "${!PIDS[@]}"; do
  if wait "${PIDS[$tag]}"; then
    echo "=== [$(date +%H:%M:%S)] e_${tag} OK ==="
  else
    echo "=== [$(date +%H:%M:%S)] e_${tag} FAILED (exit $?) ==="
  fi
done

echo
echo "=== [$(date +%H:%M:%S)] OVERNIGHT BATCH DONE ==="
date
