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

# T8 Complete-Odd PGT landscape — INTR queue, t_end=1, first attempt.
#
# ndim=53 (d16/xt24/xt32/xt34/zt4/kappa vanish in plane-wave reduction).
# nlive=1325 (25×ndim), num_repeats=106 (2×ndim).
# t_end=1: modal solver expm safe to γ>709. Amp+sup run in parallel (2×38 ranks).
#
# If the 1hr INTR wall is hit before convergence, resubmit a follow-up
# pointing --output at this job's results dir and add --read-resume.
#
# Theory: torsion_gertsenshtein_complete_odd.json
#   T6 base (beta1-3, xi, delta1, chi, zeta1-3, d14-15, d17, d19-21, zt1-3, zt5-6)
#   + xt1-36 (ε·R̃×∂T sector; xt24/xt32/xt34 vanish).
#   Total active params: 53.
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_t8_landscape.sbatch \
#     --ntasks 76 --nodes 1 --name t8-landscape \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=53
NLIVE=1325  # 25 * NDIM
NREP=106    # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== T8 Complete-Odd PGT landscape (t_end=1, ndim=53, nlive=1325, nrep=106) ==="
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

DATA="examples/data/torsion_gertsenshtein_complete_odd.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 53 priors — d16/xt24/xt32/xt34/zt4/kappa excluded (vanish in plane-wave reduction)
PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior xi=log_uniform:1e-3:1e3
  --prior delta1=arctan_uniform:-89:89
  --prior chi=arctan_uniform:-89:89
  --prior zeta1=arctan_uniform:-89:89
  --prior zeta2=arctan_uniform:-89:89
  --prior zeta3=arctan_uniform:-89:89
  --prior d14=arctan_uniform:-89:89
  --prior d15=arctan_uniform:-89:89
  --prior d17=arctan_uniform:-89:89
  --prior d19=arctan_uniform:-89:89
  --prior d20=arctan_uniform:-89:89
  --prior d21=arctan_uniform:-89:89
  --prior zt1=arctan_uniform:-89:89
  --prior zt2=arctan_uniform:-89:89
  --prior zt3=arctan_uniform:-89:89
  --prior zt5=arctan_uniform:-89:89
  --prior zt6=arctan_uniform:-89:89
  --prior xt1=arctan_uniform:-89:89
  --prior xt2=arctan_uniform:-89:89
  --prior xt3=arctan_uniform:-89:89
  --prior xt4=arctan_uniform:-89:89
  --prior xt5=arctan_uniform:-89:89
  --prior xt6=arctan_uniform:-89:89
  --prior xt7=arctan_uniform:-89:89
  --prior xt8=arctan_uniform:-89:89
  --prior xt9=arctan_uniform:-89:89
  --prior xt10=arctan_uniform:-89:89
  --prior xt11=arctan_uniform:-89:89
  --prior xt12=arctan_uniform:-89:89
  --prior xt13=arctan_uniform:-89:89
  --prior xt14=arctan_uniform:-89:89
  --prior xt15=arctan_uniform:-89:89
  --prior xt16=arctan_uniform:-89:89
  --prior xt17=arctan_uniform:-89:89
  --prior xt18=arctan_uniform:-89:89
  --prior xt19=arctan_uniform:-89:89
  --prior xt20=arctan_uniform:-89:89
  --prior xt21=arctan_uniform:-89:89
  --prior xt22=arctan_uniform:-89:89
  --prior xt23=arctan_uniform:-89:89
  --prior xt25=arctan_uniform:-89:89
  --prior xt26=arctan_uniform:-89:89
  --prior xt27=arctan_uniform:-89:89
  --prior xt28=arctan_uniform:-89:89
  --prior xt29=arctan_uniform:-89:89
  --prior xt30=arctan_uniform:-89:89
  --prior xt31=arctan_uniform:-89:89
  --prior xt33=arctan_uniform:-89:89
  --prior xt35=arctan_uniform:-89:89
  --prior xt36=arctan_uniform:-89:89"

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching T8 amp (${NTASKS_PER_CHAIN}r) + T8 sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    ndim=53, nlive=1325, num_repeats=106, t_end=1, soft_floor_logl=-15"
echo "    Est. convergence: ~3300 dead pts"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/t8_amp_v1" \
    > "${RESULTS_DIR}/t8_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/t8_sup_v1" \
    > "${RESULTS_DIR}/t8_sup.log" 2>&1 &
PID_B=$!

echo "    T8 amp pid=${PID_A}  T8 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — T8 amp rc=${RC_A}  T8 sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
