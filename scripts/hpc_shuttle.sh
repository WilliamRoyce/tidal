#!/usr/bin/env bash
# hpc_shuttle.sh — local <-> CSD3 shuttle for the tidal project.
#
# Never edit files on the remote directly. Never rebuild the dev environment
# there. Never poll squeue/sinfo in a loop. Never retry on SSH auth failure
# (Fail2Ban 20-min block). See docs: CLAUDE.md § HPC Workflow.

set -euo pipefail

# --- Config ---------------------------------------------------------------
readonly HOST="${HPC_HOST:-csd3}"                      # from ~/.ssh/config
readonly REMOTE_ROOT="${HPC_ROOT:-/rds/user/wr286/hpc-work/tidal}"
readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly JOBS_FILE="${REPO_ROOT}/scripts/.hpc_jobs"
readonly TEMPLATE_DIR="${REPO_ROOT}/scripts/hpc_templates"

# --- Helpers --------------------------------------------------------------
die() { echo "error: $*" >&2; exit 1; }
note() { echo "==> $*" >&2; }

# Wrap any remote ssh/rsync to fail fast on auth errors (no retry).
remote_exec() {
  ssh -o BatchMode=no -o NumberOfPasswordPrompts=1 "$HOST" "$@"
}

check_master() {
  # Verify the ControlMaster is live so commands don't trigger a fresh MFA.
  if ! ssh -O check "$HOST" 2>/dev/null; then
    cat >&2 <<EOF
error: no live SSH ControlMaster for '$HOST'.

Open one interactively first (you'll be asked for CRSid password + TOTP):

    ssh $HOST

Then re-run this command. Don't retry on auth failure — Fail2Ban will
block this IP for 20 minutes after repeated failures.
EOF
    exit 2
  fi
}

# --- Subcommands ----------------------------------------------------------

cmd_push() {
  check_master
  note "rsync local -> ${HOST}:${REMOTE_ROOT}"
  rsync -az --delete \
    --exclude='.git/' \
    --exclude='__pycache__/' \
    --exclude='.pytest_cache/' \
    --exclude='.ruff_cache/' \
    --exclude='.venv/' \
    --exclude='.claude*/' \
    --exclude='/tmp/' \
    --exclude='/build/' \
    --exclude='*.egg-info/' \
    --exclude='tidal/wolfram/' \
    --exclude='*.wls' \
    --exclude='/PolyChordLite/' \
    --exclude='/literature/' \
    --exclude='/third_party/' \
    --exclude='/xCobaInspirationForClaude/' \
    --exclude='/manuscript/' \
    --exclude='/docs/build/' \
    --exclude='/research/' \
    --exclude='*.png' \
    -e "ssh" \
    "${REPO_ROOT}/" "${HOST}:${REMOTE_ROOT}/"
}

cmd_setup() {
  check_master
  note "harvesting existing sbatch templates from remote"
  remote_exec "ls -1 \$HOME/*.sbatch 2>/dev/null; ls -1 /usr/local/Cluster-Docs/SLURM/ 2>/dev/null" || true

  note "mybalance (for account selection):"
  remote_exec "mybalance" || echo "  (mybalance not on PATH yet?)"

  note "creating ${REMOTE_ROOT}/{slurm_logs,.venv}"
  remote_exec "mkdir -p ${REMOTE_ROOT}/slurm_logs"

  cmd_push

  note "creating remote venv + installing tidal (python-only, no wolfram extras)"
  # python/3.11.0-icl is available under BOTH default-ccl (login) AND
  # default-sar (sapphire runtime), so build and runtime use the same python.
  remote_exec "source /etc/profile.d/modules.sh && \
    module purge && module load rhel8/default-ccl && \
    module load python/3.11.0-icl && \
    cd ${REMOTE_ROOT} && \
    rm -rf .venv && \
    python -m venv .venv && \
    .venv/bin/pip install --upgrade pip wheel && \
    .venv/bin/pip install -e . && \
    .venv/bin/tidal --version"
}

cmd_resolve_account() {
  # Prefer DiRAC > SL2 > SL3. Echo the chosen account to stdout.
  check_master
  local balance
  balance="$(remote_exec "mybalance" 2>/dev/null || true)"
  # Extract project names from typical mybalance columns.
  # Heuristic: pick first DiRAC project, else first SL2, else first SL3.
  local dirac sl2 sl3
  dirac="$(awk '/[Dd]i[Rr][Aa][Cc]/ {print $1; exit}' <<<"$balance")"
  sl2="$(awk '/SL2-CPU/ {print $1; exit}' <<<"$balance")"
  sl3="$(awk '/SL3-CPU/ {print $1; exit}' <<<"$balance")"
  local chosen="${dirac:-${sl2:-${sl3:-}}}"
  [[ -n "$chosen" ]] || die "could not parse any project from mybalance — pass --account explicitly"
  note "account resolved: $chosen (DiRAC=${dirac:-none} SL2=${sl2:-none} SL3=${sl3:-none})"
  echo "$chosen"
}

cmd_submit() {
  local template="" cmd="" nodes=1 time="01:00:00" name="tidal" account="" ntasks=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --template) template="$2"; shift 2 ;;
      --cmd)      cmd="$2"; shift 2 ;;
      --nodes)    nodes="$2"; shift 2 ;;
      --ntasks)   ntasks="$2"; shift 2 ;;
      --time)     time="$2"; shift 2 ;;
      --name)     name="$2"; shift 2 ;;
      --account)  account="$2"; shift 2 ;;
      *) die "unknown arg: $1" ;;
    esac
  done
  [[ -f "$template" ]] || die "template not found: $template"
  [[ -n "$cmd" ]] || die "--cmd required"
  [[ -n "$account" ]] || account="$(cmd_resolve_account)"

  # Default ntasks to 112 * nodes (sapphire cores/node) unless overridden.
  # MPI templates (e.g. scripts/hpc_templates/polychord_intr.sbatch) need
  # the caller to pass --ntasks N explicitly so ${MPIRUN_PREFIX} launches
  # exactly N ranks.  tidal sweep uses multiprocessing.Pool
  # (single-node shared-memory), so --nodes=1 is always correct for sweeps.
  if [[ -z "$ntasks" ]]; then
    ntasks=$(( nodes * 112 ))
  fi
  local rendered
  rendered="$(sed \
    -e "s|{{JOB_NAME}}|${name}|g" \
    -e "s|{{ACCOUNT}}|${account}|g" \
    -e "s|{{NODES}}|${nodes}|g" \
    -e "s|{{NTASKS}}|${ntasks}|g" \
    -e "s|{{TIME}}|${time}|g" \
    -e "s|{{COMMAND}}|${cmd}|g" \
    -e "s|{{REMOTE_ROOT}}|${REMOTE_ROOT}|g" \
    "$template")"

  note "submitting to sbatch on ${HOST}"
  local jobid
  jobid="$(remote_exec "cd ${REMOTE_ROOT} && sbatch --parsable" <<<"$rendered")"
  jobid="${jobid%%;*}"  # strip ;cluster suffix if any
  [[ -n "$jobid" ]] || die "sbatch returned empty job id"
  mkdir -p "$(dirname "$JOBS_FILE")"
  echo "$(date -Iseconds) $jobid $name $template $cmd" >> "$JOBS_FILE"
  echo "$jobid"
}

cmd_status() {
  check_master
  local jobid="${1:-}"
  if [[ -n "$jobid" ]]; then
    remote_exec "squeue -j $jobid -o '%.10i %.12j %.8T %.10M %.6D %R'"
  else
    remote_exec "squeue -u \$USER -o '%.10i %.12j %.8T %.10M %.6D %R'"
  fi
}

cmd_tail() {
  check_master
  local jobid="${1:?jobid required}"
  local follow="${2:-}"
  local logfile="${REMOTE_ROOT}/slurm_logs/slurm-${jobid}.out"
  if [[ "$follow" == "--follow" || "$follow" == "-f" ]]; then
    remote_exec "tail -f ${logfile}"
  else
    remote_exec "tail -n 200 ${logfile}"
  fi
}

cmd_htop() {
  check_master
  local jobid="${1:?jobid required}"
  local node
  node="$(remote_exec "squeue -j $jobid -h -o %N")"
  [[ -n "$node" ]] || die "job $jobid not running or has no node assigned"
  note "attaching htop on compute node: $node"
  ssh -t "$HOST" "ssh -t $node htop"
}

cmd_pull() {
  check_master
  local jobid="" all="" src=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --all)  all="--all"; shift ;;
      --src)  src="$2"; shift 2 ;;
      *)      jobid="$1"; shift ;;
    esac
  done
  [[ -n "$jobid" ]] || die "jobid required"

  if [[ -z "$src" ]]; then
    # Parse the --output PATH recorded for this jobid in .hpc_jobs.
    # Records are space-separated: <iso-date> <jobid> <name> <template> <full-cmd>
    # The full-cmd contains "--output <path>" which we extract.
    if [[ -f "$JOBS_FILE" ]]; then
      local line
      line="$(awk -v j="$jobid" '$2==j' "$JOBS_FILE" | tail -1)"
      if [[ -n "$line" ]]; then
        src="$(sed -n 's/.*--output \([^ ]*\).*/\1/p' <<<"$line")"
      fi
    fi
    [[ -n "$src" ]] || die "could not resolve remote output path for job $jobid; pass --src PATH explicitly"
    note "resolved remote src: $src"
  fi
  src="${src%/}/"  # ensure trailing slash for rsync

  local dst="${REPO_ROOT}/hpc_results/${jobid}/"
  mkdir -p "$dst"
  if [[ "$all" == "--all" ]]; then
    note "WARNING: --all pulls raw sim data; may be large"
    rsync -az -e ssh "${HOST}:${src}" "$dst"
  else
    note "pulling lightweight artefacts only (figures, csvs, summary/result/sweep jsons, nested-sampling chains)"
    rsync -az -e ssh \
      --include='*/' \
      --include='figures/***' \
      --include='*.csv' \
      --include='results.json' \
      --include='sweep.json' \
      --include='*_summary.json' \
      --include='SWEEP_RESULTS.md' \
      --include='inference.json' \
      --include='importance.json' \
      --include='*.png' \
      --include='_chains/*.txt' \
      --include='_chains/*.stats' \
      --include='_chains/*.paramnames' \
      --exclude='_chains/*.resume' \
      --exclude='*' \
      "${HOST}:${src}" "$dst"
  fi
  note "pulled to: $dst"
}

cmd_cancel() {
  check_master
  local jobid="${1:?jobid required}"
  remote_exec "scancel $jobid"
}

cmd_shell() {
  check_master
  ssh -t "$HOST" "cd ${REMOTE_ROOT} && exec bash -l"
}

usage() {
  cat <<EOF
Usage: $0 <subcommand> [args...]

Subcommands:
  push                              rsync local tidal tree -> remote (Python only, no Wolfram)
  setup                             one-time: harvest templates, create remote venv, install tidal
  submit --template T --cmd C [opts]  render sbatch template and submit via ssh
         [--nodes N] [--ntasks N] [--time HH:MM:SS] [--name X] [--account P]
  status [jobid]                    one-shot squeue
  tail <jobid> [--follow|-f]        tail remote slurm log
  htop <jobid>                      attach htop on the compute node
  pull <jobid> [--all] [--src PATH] rsync lightweight artefacts back (--all for raw data).
                                    Source path auto-parsed from --output in .hpc_jobs;
                                    override with --src for ad-hoc jobs.
  cancel <jobid>                    scancel
  shell                             interactive ssh into remote tidal dir

Env:
  HPC_HOST (default: csd3)
  HPC_ROOT (default: ~/rds/hpc-work/tidal)
EOF
}

main() {
  local sub="${1:-}"; shift || true
  case "$sub" in
    push)    cmd_push "$@" ;;
    setup)   cmd_setup "$@" ;;
    submit)  cmd_submit "$@" ;;
    status)  cmd_status "$@" ;;
    tail)    cmd_tail "$@" ;;
    htop)    cmd_htop "$@" ;;
    pull)    cmd_pull "$@" ;;
    cancel)  cmd_cancel "$@" ;;
    shell)   cmd_shell "$@" ;;
    resolve-account) cmd_resolve_account ;;
    ""|-h|--help) usage ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
