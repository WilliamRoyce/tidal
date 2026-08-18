#!/usr/bin/env bash
# Landscape pass for Stage A (DPP) and D1 (Ricci-EM) — nlive=25×ndim.
#
# Priority: complete ALL landscape surveys before any publication-quality runs.
# D2.x landscape already done (jobs 29468539/29468763/29471255).
# This covers the two remaining theories.
#
# Run INSIDE an interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_landscape_batch_stagea_d1.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_landscape_stagea_d1.log 2>&1 &
#
# Strategy:
#   Phase 1: Stage A amp (76r, ~2 min) + Stage A sup (76r, ~2 min) — sequential, fast
#   Phase 2: D1 amp (38r) + D1 sup (38r) in parallel — est ~60 min each at 38r
#            (grid=512 makes each sim 4× slower than D2.x; expect SIGTERM at ~90%)
#   Resume:  run_landscape_batch_d1_resume.sh with PREV_JOB_ID=<this jobid>
#
# After the session: bash scripts/hpc_shuttle.sh pull <jobid>

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

echo "=== [$(date +%H:%M:%S)] Setting up PolyChord environment ==="
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"
export TIDAL_MODAL_BACKEND=jax

. /etc/profile.d/modules.sh
module load rhel8/default-icl   2>/dev/null || true
module load intel/2019.3.199    2>/dev/null || true
module load intel-oneapi-mpi    2>/dev/null || true

python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py; print('venv OK')"

cd "${TIDAL_ROOT}"
export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
echo "RESULTS_DIR=${RESULTS_DIR}"
echo ""

# ── Phase 1a: Stage A amp landscape (76 ranks, ~2 min) ───────────────────────
echo "=== [$(date +%H:%M:%S)] Phase 1a: Stage A amp landscape (76 ranks) ==="
mpirun -n 76 tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:1e-3:1e3" \
    --prior "deltam=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output "${RESULTS_DIR}/stage_a_amp_v3"
echo "=== [$(date +%H:%M:%S)] Phase 1a done ==="
echo ""

# ── Phase 1b: Stage A sup landscape (76 ranks, ~2 min) ───────────────────────
echo "=== [$(date +%H:%M:%S)] Phase 1b: Stage A sup landscape (76 ranks) ==="
mpirun -n 76 tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:1e-3:1e3" \
    --prior "deltam=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output "${RESULTS_DIR}/stage_a_sup_v3"
echo "=== [$(date +%H:%M:%S)] Phase 1b done ==="
echo ""

# ── Phase 2: D1 amp + D1 sup landscape (38+38 ranks, parallel) ───────────────
echo "=== [$(date +%H:%M:%S)] Phase 2: D1 amp + D1 sup landscape (38+38 ranks) ==="
echo "    est ~60 min each at 38r (grid=512); expect SIGTERM at ~90%"
echo "    resume with run_landscape_batch_d1_resume.sh PREV_JOB_ID=${SLURM_JOB_ID}"

mpirun -n 38 tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output "${RESULTS_DIR}/d1_amp_v3" \
    > "${RESULTS_DIR}/d1_amp.log" 2>&1 &
PID_A=$!

mpirun -n 38 tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output "${RESULTS_DIR}/d1_sup_v3" \
    > "${RESULTS_DIR}/d1_sup.log" 2>&1 &
PID_B=$!

echo "    D1 amp pid=${PID_A}  D1 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Phase 2 done — D1 amp rc=${RC_A}  D1 sup rc=${RC_B} ==="
echo ""
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
echo "=== If D1 SIGTERM'd: export PREV_JOB_ID=${SLURM_JOB_ID}"
echo "===                  bash run_landscape_batch_d1_resume.sh ==="
