#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks=1
#SBATCH --time={{TIME}}
#SBATCH --output={{REMOTE_ROOT}}/logs/{{JOB_NAME}}_%j.log
#SBATCH --exclude=cpu-r-10,cpu-r-15

# Reusable template for PolyChord nested-sampling jobs on CSD3.  Expects
# the caller to have pushed tidal with `hpc_shuttle.sh push` and to have
# PolyChordLite cloned at ${REMOTE_ROOT}/PolyChordLite.
#
# The heavy lifting in this template is the "build PolyChord + inject
# into site-packages" dance we need because (a) pypolychord is yanked
# from PyPI so we compile from source, and (b) tidal runs against a
# /tmp-extracted site-packages copy (avoids /rds Lustre hangs — see
# issue #269 and session log 2026-04-19).

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
POLYDIR="${TIDAL_ROOT}/PolyChordLite"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname) ==="
date

# 1. Extract the pre-built venv tarball from /home (NFS, avoids /rds Lustre
#    hangs observed on some compute nodes — cpu-r-10/cpu-r-15 excluded above).
echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# 2. Build PolyChord.  `intel` module provides ifort for libchord.so;
#    `setup.py build_ext` needs CC=gcc because the intel module pollutes
#    CC with icx, which distutils can't introspect (fails with "icx -v"
#    non-zero exit).
echo "--- Building PolyChord ---"
cd "${POLYDIR}"
module load intel/2019.3.199 2>/dev/null || true
make all MPI=0 2>&1 | tail -2
CC=gcc python setup.py build_ext --inplace 2>&1 | tail -3

# 3. Inject compiled extension into local site-packages.  Layout matters:
#    - site-packages/pypolychord/                  — Python package
#    - site-packages/_pypolychord.*.so             — top-level (the
#      package does `import _pypolychord`, NOT relative)
#    - site-packages/pypolychord/lib/libchord.so   — matches the compiled
#      rpath `$ORIGIN/pypolychord/lib`
echo "--- Installing pypolychord into local site-packages ---"
SP_DIR="${LOCAL_SP}/site-packages"
rm -rf "${SP_DIR}/pypolychord" "${SP_DIR}"/_pypolychord*.so
cp -r "${POLYDIR}/pypolychord" "${SP_DIR}/"
cp "${POLYDIR}"/_pypolychord*.so "${SP_DIR}/"
mkdir -p "${SP_DIR}/pypolychord/lib"
cp "${POLYDIR}/lib/libchord.so" "${SP_DIR}/pypolychord/lib/"
python -c "from pypolychord import run_polychord; import anesthetic; print('imports OK')"

# 4. Run the caller's tidal command.  Use --likelihood-backend memory
#    (default since v0.33.x) for ~10x per-eval speedup and
#    --no-clustering / --num-repeats=N*ndim / --precision-criterion=0.01
#    for Handley-recommended exploration settings on unimodal posteriors.
echo "--- Running tidal command ---"
cd "${TIDAL_ROOT}"
{{COMMAND}}

echo "=== DONE ==="
date
