#!/usr/bin/env bash
# Phase E interactive amp+sup parallel runner
#
# Runs the amp + sup PolyChord chains for one or more Phase E theories in
# PARALLEL inside a single INTR sapphire allocation, maximizing HPC throughput.
#
# Usage (from local devcontainer):
#   1. Book a sapphire INTR slot (1 hour):
#        bash scripts/hpc_shuttle.sh interactive
#      → returns <jobid>
#   2. SSH into the compute node:
#        bash scripts/hpc_shuttle.sh attach <jobid>
#   3. On the compute node, source the Phase E geometry and run this script:
#        bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3e_localised/run_interactive_amp_sup.sh <theory_tag>
#      where <theory_tag> ∈ {t1, t4, t5, np, t7s, t8s}
#
# Within one slot, amp+sup chains share the 112 sapphire cores: 32 ranks each
# (64 total), leaving ~48 cores headroom for OS/IO. Memory: 192GB/node ÷ 64
# ranks ≈ 3GB/rank — well above the E.T2 OOM-at-1.5GB observation.
#
# To chain multiple theories in one slot, run them sequentially:
#   bash ...run_interactive_amp_sup.sh t4    # ~30 min
#   bash ...run_interactive_amp_sup.sh t5    # ~30 min
# (each call runs amp+sup PARALLEL but DIFFERENT THEORIES SEQUENTIAL)

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <theory_tag>"
  echo "  theory_tag: t1 | t4 | t5 | np | t7s | t8s"
  exit 1
fi

THEORY_TAG="$1"
TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"
PHASE_E_DIR="${TIDAL_ROOT}/scripts/hpc_submit_drafts/v3e_localised"

cd "${TIDAL_ROOT}"

# Source the canonical Phase E geometry.
# shellcheck disable=SC1091
source "${PHASE_E_DIR}/_geometry.env"

# Map tag → theory metadata.
case "${THEORY_TAG}" in
  t1)
    DATA="examples/data/dark_photon_plasma_e_dual_gaussian.json"
    PRIORS="--prior alpha3=log_uniform:1e-3:1e3 --prior xi=log_uniform:1e-3:1e3 \
            --prior deltam=arctan_uniform:-89:89 --prior mA2=log_uniform:1e-3:1e3"
    PARAMS=""
    NLIVE=100 ; NREP=8
    ;;
  t4)
    DATA="examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"
    PRIORS="--prior alpha1=arctan_uniform:-89:89 --prior alpha2=arctan_uniform:-89:89 \
            --prior alpha3=arctan_uniform:-89:89 --prior delta1=arctan_uniform:-89:89"
    PARAMS=""
    NLIVE=100 ; NREP=8
    ;;
  t5)
    DATA="examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json"
    PRIORS="--prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
            --prior beta3=arctan_uniform:-89:89 --prior xi=log_uniform:1e-3:1e3 \
            --prior delta1=arctan_uniform:-89:89"
    PARAMS="--param chi=0.0 --param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0"
    NLIVE=125 ; NREP=10
    ;;
  np)
    DATA="examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json"
    PRIORS="--prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
            --prior beta3=arctan_uniform:-89:89 --prior delta1=arctan_uniform:-89:89 \
            --prior chi=arctan_uniform:-89:89 \
            --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 \
            --prior zeta3=arctan_uniform:-89:89"
    PARAMS="--param xi=0.0"
    NLIVE=200 ; NREP=16
    ;;
  t7s)
    DATA="examples/data/torsion_gertsenshtein_complete_even_e_dual_gaussian.json"
    PRIORS="--prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
            --prior beta3=arctan_uniform:-89:89 --prior xi=log_uniform:1e-3:1e3 \
            --prior delta1=arctan_uniform:-89:89 \
            --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 \
            --prior zeta3=arctan_uniform:-89:89"
    PARAMS="--param chi1=0.0 --param chi2=0.0 --param chi3=0.0 --param chi4=0.0 \
            --param chi5=0.0 --param chi6=0.0 --param chi7=0.0 --param chi8=0.0 \
            --param chi9=0.0 --param chi10=0.0"
    NLIVE=200 ; NREP=16
    ;;
  *)
    echo "Unknown theory tag: ${THEORY_TAG}"
    exit 1
    ;;
esac

# ── Environment setup (compute node) ────────────────────────────────────────
echo "=== [$(date +%H:%M:%S)] Phase E interactive amp+sup runner — theory=${THEORY_TAG} ==="
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
if [ ! -d "${LOCAL_SP}/site-packages" ]; then
  echo "Extracting venv tarball..."
  mkdir -p "${LOCAL_SP}"
  tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"
fi

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# shellcheck disable=SC1091
. /etc/profile.d/modules.sh
module load rhel8/default-icl   2>/dev/null || true
module load intel/2019.3.199    2>/dev/null || true
module load intel-oneapi-mpi    2>/dev/null || true

python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py; print('venv OK')"

export RESULTS_DIR="${RESULTS_DIR:-${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}}"
mkdir -p "${RESULTS_DIR}/e_${THEORY_TAG}_amp" "${RESULTS_DIR}/e_${THEORY_TAG}_sup"

# Shared physics flags (sourced from _geometry.env).
PHYS=" \
  --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
  ${PARAMS} \
  --baseline-formula sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2 \
  --soft-floor-noise 1.0 \
  --grid-shape ${N} --bounds=0:${L} --periodic \
  --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
  --t-end ${T_END} --snapshots ${SNAPSHOTS} \
  --measure conversion,peak_conversion --source h_5 --target a_1 \
  --method nested --sampler polychord --nlive ${NLIVE} --num-repeats ${NREP} \
  --precision-criterion 0.01"

# ── Launch amp + sup in PARALLEL ────────────────────────────────────────────
# 32 ranks per chain × 2 chains = 64 cores total (out of 112 sapphire).
# Memory: 192GB/node ÷ 64 ranks ≈ 3GB/rank, well above the E.T2 OOM-at-1.5GB.
echo "=== [$(date +%H:%M:%S)] Launching e_${THEORY_TAG}_amp (32 ranks) in background ==="
# shellcheck disable=SC2086
mpirun -n 32 tidal sample ${DATA} ${PRIORS} ${PHYS} \
  --likelihood "P_max:maximize" \
  --output "${RESULTS_DIR}/e_${THEORY_TAG}_amp" \
  > "${RESULTS_DIR}/e_${THEORY_TAG}_amp.log" 2>&1 &
PID_AMP=$!

echo "=== [$(date +%H:%M:%S)] Launching e_${THEORY_TAG}_sup (32 ranks) in background ==="
# shellcheck disable=SC2086
mpirun -n 32 tidal sample ${DATA} ${PRIORS} ${PHYS} \
  --likelihood "P_max:minimize" \
  --output "${RESULTS_DIR}/e_${THEORY_TAG}_sup" \
  > "${RESULTS_DIR}/e_${THEORY_TAG}_sup.log" 2>&1 &
PID_SUP=$!

echo "=== [$(date +%H:%M:%S)] Both chains launched: amp=${PID_AMP} sup=${PID_SUP} — waiting for both ==="
wait ${PID_AMP} ${PID_SUP}
echo "=== [$(date +%H:%M:%S)] DONE — both chains finished ==="
date
