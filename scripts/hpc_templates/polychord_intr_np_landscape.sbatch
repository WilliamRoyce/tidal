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

# Non-propagating torsion landscape — D2.3 with xi=0 pinned.
#
# Theory: torsion_gertsenshtein_general_nonminimal.json (same as D2.0-D2.3).
#   xi is fixed at 0.0 → torsion becomes purely algebraic (no kinetic (∇T)² term).
#   All torsion-EM couplings remain free: beta1-3, delta1, chi, zeta1-3. ndim=8.
#
# Physics question: does torsion-EM coupling alone (without propagating torsion)
#   produce Gertsenshtein enhancement? If NP is non-null but D2.3 is null, the
#   coupling structure is sufficient. If NP is null but D2.3 is non-null, torsion
#   propagation (xi≠0) is the necessary ingredient.
#
# nlive=200 (25×ndim), num_repeats=16 (2×ndim), t_end=1.
# NTASKS_PER_CHAIN=38 → 76 total tasks, 1 icelake node.
# Est. wall time: ~15-25 min.
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_np_landscape.sbatch \
#     --ntasks 76 --nodes 1 --name np-landscape \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=8
NLIVE=200  # 25 * NDIM
NREP=16    # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== Non-propagating torsion landscape (xi=0, ndim=8, nlive=200, nrep=16, t_end=1) ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

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

DATA="examples/data/torsion_gertsenshtein_general_nonminimal.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01 --param xi=0.0
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 8 priors — xi pinned to 0 via --param above, not a free prior
PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior delta1=arctan_uniform:-89:89
  --prior chi=arctan_uniform:-89:89
  --prior zeta1=arctan_uniform:-89:89
  --prior zeta2=arctan_uniform:-89:89
  --prior zeta3=arctan_uniform:-89:89"

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching NP amp (${NTASKS_PER_CHAIN}r) + sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    xi=0 (pinned), ndim=8, nlive=200, num_repeats=16, t_end=1, soft_floor_logl=-15"
echo "    Free params: beta1-3, delta1, chi, zeta1-3 (all torsion-EM couplings)"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/np_amp_v1" \
    > "${RESULTS_DIR}/np_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/np_sup_v1" \
    > "${RESULTS_DIR}/np_sup.log" 2>&1 &
PID_B=$!

echo "    NP amp pid=${PID_A}  sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — NP amp rc=${RC_A}  NP sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
