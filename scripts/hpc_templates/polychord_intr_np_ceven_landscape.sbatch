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

# NP Complete-Even: non-propagating T7 torsion (xi=0 pinned).
#
# Theory: torsion_gertsenshtein_complete_even.json with --param xi=0.0
# Lagrangian: T7 basis (beta1-3, delta1, zeta1-3, chi1-10) but xi=0 pinned →
#   algebraic/Einstein-Cartan-like torsion; no independent torsion wave equation.
#   Full parity-even interaction sector (chi1-10, R̃×∂T basis) remains free.
#
# Science question: Does T7's amp logZ=+10.30 (v2, 29682868) require torsion
#   propagation (xi), or is the chi1-10 coupling structure alone sufficient?
#   Direct comparison: NP.T7 (ndim=17, xi=0) vs T7v2 (ndim=18, xi free).
#   If NP.T7 matches T7v2: enhancement is driven by chi couplings, not kinetics.
#   If NP.T7 < T7v2: xi (propagating torsion DOF) is the necessary ingredient.
#
# Also extends the baseline NP run (29700462, ndim=8, general_nonminimal) which
#   only had T5-level couplings (beta1-3, delta1, chi, zeta1-3). NP.T7 adds the
#   chi1-10 (R̃×∂T) sector absent from NP.T5 — richer non-propagating theory.
#
# ndim=17, nlive=425 (25×ndim), num_repeats=34 (2×ndim).
# Est. wall time: ~25-35 min at 76 ranks (38r amp + 38r sup in parallel).
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_np_ceven_landscape.sbatch \
#     --ntasks 76 --nodes 1 --name np-ceven \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=17
NLIVE=425   # 25 * NDIM
NREP=34     # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== NP Complete-Even (T7 + xi=0, ndim=17, nlive=425, nrep=34, t_end=1) ==="
echo "=== Science: isolate chi1-10 coupling contribution vs T7v2 (xi free) ==="
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

DATA="examples/data/torsion_gertsenshtein_complete_even.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01 --param xi=0.0
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 17 priors: T7 complete-even minus xi (pinned to 0 above)
PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior delta1=arctan_uniform:-89:89
  --prior zeta1=arctan_uniform:-89:89
  --prior zeta2=arctan_uniform:-89:89
  --prior zeta3=arctan_uniform:-89:89
  --prior chi1=arctan_uniform:-89:89
  --prior chi2=arctan_uniform:-89:89
  --prior chi3=arctan_uniform:-89:89
  --prior chi4=arctan_uniform:-89:89
  --prior chi5=arctan_uniform:-89:89
  --prior chi6=arctan_uniform:-89:89
  --prior chi7=arctan_uniform:-89:89
  --prior chi8=arctan_uniform:-89:89
  --prior chi9=arctan_uniform:-89:89
  --prior chi10=arctan_uniform:-89:89"

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching NP-CEven amp (${NTASKS_PER_CHAIN}r) + sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    ndim=17, nlive=425, num_repeats=34, t_end=1, xi=0 pinned"
echo "    Comparison: T7v2 amp logZ=+10.30 (xi free, ndim=18)"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/np_ceven_amp_v1" \
    > "${RESULTS_DIR}/np_ceven_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/np_ceven_sup_v1" \
    > "${RESULTS_DIR}/np_ceven_sup.log" 2>&1 &
PID_B=$!

echo "    NP-CEven amp pid=${PID_A}  sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — amp rc=${RC_A}  sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
