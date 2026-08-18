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
    --exclude='/results/' \
    --exclude='/logs/' \
    --exclude='/slurm_logs/' \
    --exclude='/hpc_results/' \
    -e "ssh" \
    "${REPO_ROOT}/" "${HOST}:${REMOTE_ROOT}/"
  _assert_version_sync || die \
    "HPC tidal version drifted. Re-install + refresh tarball before submitting jobs:
  ssh $HOST 'cd ${REMOTE_ROOT} && source .venv/bin/activate && pip install -e . --quiet'
  bash scripts/hpc_refresh_venv_tar.sh"
}

_assert_version_sync() {
  # Verify HPC package metadata matches local source version.  rsync
  # of source alone leaves .dist-info/ stale (editable install records
  # the version at install time).  A stale metadata on its own is
  # usually cosmetic, but the tarball built from that metadata is
  # consumed by compute-node jobs — if the tarball carries old
  # site-packages/tidal bytecode, jobs quietly run pre-fix code
  # against new configs (see #307).  Fail fast when versions diverge.
  local local_version hpc_version tarball_version
  local_version="$(awk -F'"' '/^version/{print $2; exit}' "${REPO_ROOT}/pyproject.toml" 2>/dev/null || true)"
  [[ -n "$local_version" ]] || { note "WARN: could not read local tidal version"; return 0; }
  hpc_version="$(remote_exec "cd ${REMOTE_ROOT} && source .venv/bin/activate 2>/dev/null && pip show tidal 2>/dev/null | awk '/^Version:/{print \$2}'" 2>/dev/null | tr -d '[:space:]')"
  tarball_version="$(remote_exec "tar tf /home/wr286/venv_site.tar 2>/dev/null | grep -oE 'tidal-[0-9.]+\\.dist-info' | head -1 | sed 's/tidal-//;s/\\.dist-info//'" 2>/dev/null | tr -d '[:space:]')"
  note "version check: local=${local_version} hpc_venv=${hpc_version:-UNKNOWN} hpc_tarball=${tarball_version:-UNKNOWN}"
  if [[ -n "$hpc_version" && "$hpc_version" != "$local_version" ]]; then
    note "WARN: HPC venv tidal is ${hpc_version}, local is ${local_version}"
    return 1
  fi
  if [[ -n "$tarball_version" && "$tarball_version" != "$local_version" ]]; then
    note "WARN: HPC tarball tidal is ${tarball_version}, local is ${local_version}"
    return 1
  fi
  return 0
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
  # mybalance columns: user usage | ACCOUNT usage | limit available
  # The account name is in column 4 (after the first '|' separator).
  dirac="$(awk '/[Dd]i[Rr][Aa][Cc]/ {print $4; exit}' <<<"$balance")"
  sl2="$(awk '/SL2-CPU/ {print $4; exit}' <<<"$balance")"
  sl3="$(awk '/SL3-CPU/ {print $4; exit}' <<<"$balance")"
  local chosen="${dirac:-${sl2:-${sl3:-}}}"
  [[ -n "$chosen" ]] || die "could not parse any project from mybalance — pass --account explicitly"
  note "account resolved: $chosen (DiRAC=${dirac:-none} SL2=${sl2:-none} SL3=${sl3:-none})"
  echo "$chosen"
}

cmd_submit() {
  local template="" cmd="" nodes=1 time="01:00:00" name="tidal" account="" ntasks="" campaign=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --template) template="$2"; shift 2 ;;
      --cmd)      cmd="$2"; shift 2 ;;
      --nodes)    nodes="$2"; shift 2 ;;
      --ntasks)   ntasks="$2"; shift 2 ;;
      --time)     time="$2"; shift 2 ;;
      --name)     name="$2"; shift 2 ;;
      --account)  account="$2"; shift 2 ;;
      --campaign) campaign="$2"; shift 2 ;;
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
  # In sed replacement strings, '&' expands to the matched text and
  # the delimiter (here '|') terminates the replacement. Escape both
  # in user-supplied values so '&&', 'tail|less', etc. in --cmd survive.
  # Also collapse shell line-continuation sequences (\ + newline + spaces)
  # to a single space: literal newlines in the sed replacement text would
  # terminate the s command and cause "unterminated s command" errors.
  local safe_cmd safe_root safe_account
  safe_cmd="$(printf '%s' "$cmd" | sed ':a;N;$!ba;s/\\\n[[:space:]]*/ /g')"
  safe_cmd="${safe_cmd//\\/\\\\}"      # escape backslashes first
  safe_cmd="${safe_cmd//&/\\&}"        # then escape ampersands
  safe_cmd="${safe_cmd//|/\\|}"        # then escape the sed delimiter
  safe_root="${REMOTE_ROOT//&/\\&}"
  safe_root="${safe_root//|/\\|}"
  safe_account="${account//&/\\&}"
  safe_account="${safe_account//|/\\|}"

  # Campaign dir: stable path across INTR resume rounds.  When --campaign NAME
  # is given, {{CAMPAIGN_DIR_SETUP}} in the template becomes an export + mkdir
  # that sets CAMPAIGN_DIR to hpc_results/campaigns/NAME on the HPC.  Use
  # ${CAMPAIGN_DIR} in --cmd --output to share chains across rounds.
  # PolyChord writes .resume checkpoints by default; the second round adds
  # --read-resume to tidal sample to continue from the checkpoint.
  local camp_setup
  if [[ -n "$campaign" ]]; then
    local camp_path="${REMOTE_ROOT}/hpc_results/campaigns/${campaign}"
    camp_setup="export CAMPAIGN_DIR=\"${camp_path}\"; mkdir -p \"${camp_path}\""
    camp_setup="${camp_setup//&/\\&}"
    camp_setup="${camp_setup//|/\\|}"
  else
    camp_setup="# (no campaign dir — output goes to \${RESULTS_DIR})"
  fi

  local rendered
  rendered="$(sed \
    -e "s|{{JOB_NAME}}|${name}|g" \
    -e "s|{{ACCOUNT}}|${safe_account}|g" \
    -e "s|{{NODES}}|${nodes}|g" \
    -e "s|{{NTASKS}}|${ntasks}|g" \
    -e "s|{{TIME}}|${time}|g" \
    -e "s|{{COMMAND}}|${safe_cmd}|g" \
    -e "s|{{REMOTE_ROOT}}|${safe_root}|g" \
    -e "s|{{CAMPAIGN_DIR_SETUP}}|${camp_setup}|g" \
    "$template")"

  note "submitting to sbatch on ${HOST}"
  local jobid
  jobid="$(remote_exec "cd ${REMOTE_ROOT} && sbatch --parsable" <<<"$rendered")"
  jobid="${jobid%%;*}"  # strip ;cluster suffix if any
  [[ -n "$jobid" ]] || die "sbatch returned empty job id"
  mkdir -p "$(dirname "$JOBS_FILE")"
  # Expand ${CAMPAIGN_DIR} in the recorded cmd so pull <jobid> can resolve the
  # --output path without needing pull-campaign NAME.
  local recorded_cmd="$cmd"
  if [[ -n "$campaign" ]]; then
    recorded_cmd="${cmd//\$\{CAMPAIGN_DIR\}/${camp_path}}"
  fi
  echo "$(date -Iseconds) $jobid $name $template $recorded_cmd" >> "$JOBS_FILE"
  echo "$jobid"
  # User-visible reminder: HPC jobs that crash on startup typically die
  # within 5-30 s (template bugs, MPI init failures, missing modules).
  # `wait` would block for the full INTR queue time before discovering
  # the job is already dead.  `check` does a single sacct + log-tail in
  # one shot.  See feedback_hpc_check_fast_fail memory for context.
  note "tip: ~30 s after submit, run \`bash $0 check $jobid\` to confirm the job didn't fast-fail before \`wait\`-blocking"
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

cmd_check() {
  # One-shot post-submit fast-fail check.  Catches jobs that crashed on
  # startup (template bugs, unbound variables, MPI init failures, missing
  # modules, CLI flag typos) without burning the INTR queue's wait time.
  # Single sacct + single log read; respects the no-polling rule.
  check_master
  local jobid="${1:?jobid required}"
  remote_exec "sacct -j $jobid --format=JobID,State,Reason,Start,End,ExitCode -P 2>/dev/null | head -8"
  echo "---"
  local logfile="${REMOTE_ROOT}/slurm_logs/slurm-${jobid}.out"
  remote_exec "if [[ -f ${logfile} ]]; then echo '== last 30 lines of slurm log =='; tail -30 ${logfile}; else echo '(slurm log not yet written: job either still pending or has not started)'; fi"
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
    note "pulling lightweight artifacts only (figures, csvs, summary/result/sweep jsons, nested-sampling chains)"
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

cmd_pull_campaign() {
  check_master
  local name="${1:-}"; [[ -n "$name" ]] || die "usage: pull-campaign NAME"
  local src="${REMOTE_ROOT}/hpc_results/campaigns/${name}"
  local dst="${REPO_ROOT}/hpc_results/campaigns/${name}"
  mkdir -p "${dst}"
  note "pulling campaign artifacts: ${name}"
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
    "${HOST}:${src}/" "${dst}/"
  note "pulled to: ${dst}/"
}

cmd_wait() {
  # Wait for a job to start, then tail its log.  Uses file-existence polling
  # (NOT squeue) so no Slurm controller load is generated.  The wait loop runs
  # locally — no persistent background process is left on the cluster.
  # Per CSD3 admin guidance: if squeue must be called in a loop, minimum 120 s
  # between calls; prefer file-existence checks instead (man squeue PERFORMANCE).
  check_master
  local jobid="${1:?jobid required}"
  local interval="${2:-120}"
  local logfile="${REMOTE_ROOT}/slurm_logs/slurm-${jobid}.out"
  note "waiting for slurm-${jobid}.out to appear (checking every ${interval}s via stat — no squeue)"
  note "Loop runs locally; zero persistent processes left on the cluster."
  while ! remote_exec "test -f ${logfile}" 2>/dev/null; do
    note "not started yet — sleeping ${interval}s locally"
    sleep "${interval}"
  done
  note "log file appeared — job has started. Tailing (Ctrl-C to detach):"
  remote_exec "tail -30f ${logfile}"
}

cmd_cancel() {
  check_master
  local jobid="${1:?jobid required}"
  remote_exec "scancel $jobid"
}

cmd_interactive() {
  local time="01:00:00" name="tidal-interactive" account="" nodes=1 ntasks=112
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --time)    time="$2";    shift 2 ;;
      --name)    name="$2";    shift 2 ;;
      --account) account="$2"; shift 2 ;;
      --nodes)   nodes="$2";   shift 2 ;;
      --ntasks)  ntasks="$2";  shift 2 ;;
      *) die "unknown arg: $1" ;;
    esac
  done
  check_master
  [[ -n "$account" ]] || account="$(cmd_resolve_account)"

  local template="${TEMPLATE_DIR}/interactive.sbatch"
  [[ -f "$template" ]] || die "template not found: $template"

  local safe_root safe_account
  safe_root="${REMOTE_ROOT//&/\\&}"; safe_root="${safe_root//|/\\|}"
  safe_account="${account//&/\\&}"; safe_account="${safe_account//|/\\|}"

  local rendered
  rendered="$(sed \
    -e "s|{{JOB_NAME}}|${name}|g" \
    -e "s|{{ACCOUNT}}|${safe_account}|g" \
    -e "s|{{NODES}}|${nodes}|g" \
    -e "s|{{NTASKS}}|${ntasks}|g" \
    -e "s|{{TIME}}|${time}|g" \
    -e "s|{{REMOTE_ROOT}}|${safe_root}|g" \
    "$template")"

  note "submitting interactive allocation (nodes=${nodes} ntasks=${ntasks} time=${time})"
  local jobid
  jobid="$(remote_exec "cd ${REMOTE_ROOT} && sbatch --parsable" <<<"$rendered")"
  jobid="${jobid%%;*}"
  [[ -n "$jobid" ]] || die "sbatch returned empty job id"

  local results_dir="${REMOTE_ROOT}/hpc_results/${jobid}"
  mkdir -p "$(dirname "$JOBS_FILE")"
  # Record without brackets so pull's sed parser extracts --output PATH cleanly.
  echo "$(date -Iseconds) $jobid $name ${template} --output ${results_dir}" >> "$JOBS_FILE"
  echo "$jobid"
  note "interactive allocation submitted: job $jobid"
  note "tip: ~30 s after submit, run \`bash $0 check $jobid\` to confirm startup"
  note "wait for start: bash $0 wait $jobid"
  note "then attach:    bash $0 attach $jobid"
  note "cancel early:   bash $0 cancel $jobid"
}

cmd_attach() {
  check_master
  local jobid="${1:?jobid required}"
  # Get state and nodelist in one squeue call.  squeue -o %N returns reason strings
  # like '(Resources)' for pending jobs (not empty), so we must check state explicitly.
  local sq_out state node
  sq_out="$(remote_exec "squeue -j $jobid -h -o '%T,%N'" 2>/dev/null || true)"
  state="$(cut -d, -f1 <<<"$sq_out")"
  node="$(cut -d, -f2 <<<"$sq_out")"
  if [[ "$state" != "RUNNING" || -z "$node" ]]; then
    die "job $jobid is not RUNNING (state: ${state:-unknown}); use: bash $0 wait $jobid"
  fi
  local results_dir="${REMOTE_ROOT}/hpc_results/${jobid}"
  note "attaching to compute node: $node (job $jobid)"
  note "RESULTS_DIR: $results_dir"
  # SSH jump: devcontainer -> login node -> compute node.
  # Local vars ($node, ${REMOTE_ROOT}, ${results_dir}) are expanded here before the
  # string reaches the login node; the single-quoted inner command then runs verbatim
  # on the compute node.  exec bash (not bash -l) inherits the env without re-sourcing
  # profile scripts that could reset the module environment.
  ssh -t "$HOST" "ssh -t $node \
    'source /etc/profile.d/modules.sh 2>/dev/null || true; \
     module purge 2>/dev/null || true; \
     module load rhel8/default-sar 2>/dev/null || true; \
     module load python/3.11.0-icl 2>/dev/null || true; \
     cd ${REMOTE_ROOT}; \
     source .venv/bin/activate; \
     export RESULTS_DIR=${results_dir}; \
     exec bash'"
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
         [--campaign NAME]            stable campaign dir for PolyChord resume runs
                                      sets \${CAMPAIGN_DIR} in the sbatch; use in --cmd --output
  status [jobid]                    one-shot squeue (never in a loop — 120s min between calls)
  check <jobid>                     one-shot post-submit fast-fail check (sacct + log tail);
                                    use ~30s after submit before any wait-blocking
  wait <jobid> [interval_s]         wait for job to start via file-existence (NOT squeue), then tail -f
  tail <jobid> [--follow|-f]        tail remote slurm log
  htop <jobid>                      attach htop on the compute node
  pull <jobid> [--all] [--src PATH] rsync lightweight artifacts back (--all for raw data).
                                    Source path auto-parsed from --output in .hpc_jobs;
                                    override with --src for ad-hoc jobs.
  pull-campaign NAME                pull campaign artifacts from hpc_results/campaigns/NAME
                                    (use after resume runs complete with === DONE ===)
  interactive [--time HH:MM:SS] [--nodes N] [--ntasks N] [--name X] [--account P]
                                    book a full sapphire node via sbatch sleep (INTR QOS);
                                    allocation persists across SSH disconnects; follow with
                                    'wait <jobid>' then 'attach <jobid>'
  attach <jobid>                    SSH to the compute node of a running interactive job;
                                    modules loaded, venv activated, \$RESULTS_DIR set
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
    check)   cmd_check "$@" ;;
    tail)    cmd_tail "$@" ;;
    htop)    cmd_htop "$@" ;;
    pull)    cmd_pull "$@" ;;
    pull-campaign) cmd_pull_campaign "$@" ;;
    wait)    cmd_wait "$@" ;;
    interactive) cmd_interactive "$@" ;;
    attach)      cmd_attach "$@" ;;
    cancel)  cmd_cancel "$@" ;;
    shell)   cmd_shell "$@" ;;
    resolve-account) cmd_resolve_account ;;
    ""|-h|--help) usage ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
