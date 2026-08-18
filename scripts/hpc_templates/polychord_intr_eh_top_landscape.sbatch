#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-10,cpu-r-15

# EH+Gert + Topological control term landscape.
#
# Theory: gertsenshtein_eh_top.json (derived 2026-05-25).
# Lagrangian: (1/κ²)R - (1/4)F² + (ρ/8)(F²)² + (σ/8)(FF̃)² + λ_top * FF̃
#
# lambdatop * F*F~ is a Pontryagin topological density (total derivative).
# It MUST NOT appear in the EOM → posterior in lambdatop must be exactly flat
# (prior recovery, D_KL = 0). This is a methodological null control:
#   - If flat: pipeline correctly identifies topological terms as inert ✓
#   - If not flat: bug in symbolic derivation or integration pipeline ✗
#
# ndim=3 (rho, sigma, lambdatop). nlive=75 (25×ndim), num_repeats=6 (2×ndim).
# Est. wall time: ~5-10 min.
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_eh_top_landscape.sbatch \
#     --ntasks 8 --nodes 1 --name eh-top-landscape \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=3
NLIVE=75   # 25 * NDIM
NREP=6     # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== EH+Gert+Top landscape (ndim=3, nlive=75, nrep=6, t_end=1) ==="
echo "=== Control: lambdatop (Pontryagin density) posterior must be flat ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

module load rhel8/default-icl 2>/dev/null || true
module load intel/2019.3.199 2>/dev/null || true
module load intel-oneapi-mpi 2>/dev/null || true

echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from the extracted venv.
Refresh: bash scripts/hpc_install_polychord.sh && bash scripts/hpc_refresh_venv_tar.sh
EOF
    exit 3
fi

mkdir -p "${RESULTS_DIR}"
cd "${TIDAL_ROOT}"

DATA="examples/data/gertsenshtein_eh_top.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 3 priors: rho, sigma (EH physical couplings), lambdatop (topological null control)
PRIOR_ARGS="--prior rho=arctan_uniform:-89:89
  --prior sigma=arctan_uniform:-89:89
  --prior lambdatop=arctan_uniform:-89:89"

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching EH+Top amp (${NTASKS_PER_CHAIN}r) + sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    ndim=3, nlive=75, num_repeats=6, t_end=1"
echo "    Expected: rho/sigma constrained; lambdatop posterior = prior (flat)"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/eh_top_amp_v1" \
    > "${RESULTS_DIR}/eh_top_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/eh_top_sup_v1" \
    > "${RESULTS_DIR}/eh_top_sup.log" 2>&1 &
PID_B=$!

echo "    EH+Top amp pid=${PID_A}  sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — amp rc=${RC_A}  sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
