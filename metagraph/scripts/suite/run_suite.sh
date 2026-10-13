#!/usr/bin/env bash
# One command for a whole suite run (DESIGN §5, §7, §8.5): the arms against the staging host, never two at
# once, then merge, compare and report (md + pdf).
#
#   run_suite.sh [--manifest M | --smoke] [--paired | --sequential]
#                [--mcp-cmd TEMPLATE | --mcp-cmd default | --no-mcp | --mcp-run DIR]
#                [--areas P,M] [--limit N] [--skip-local] [--skip-remote] [--resume RUN_DIR]
#                [--compare EARLIER_RUN_DIR] [--by-design PATH] [--pause S] [--max-requests N]
#                [--warm-every N] [--allow-stub] [--no-pdf]
#
# --paired (the default): per scenario in manifest order, MCP cold (with --mcp-cmd) -> direct-local (one
# run_direct.py step on mex over one persistent ssh connection) -> direct-remote (a step from here through the
# public URL) -> MCP warm; one resource sampler on mex for the whole loop. --sequential: one arm after the
# other (direct-local on mex with its own sampler, then direct-remote, then the MCP arm), for runs that skip
# arms or read a finished MCP run directory (--mcp-run DIR, copied to arms/mcp).
#
# --mcp-cmd TEMPLATE: the async session's MCP runner (scripts/bench/suite_mcp_runner.py), called per scenario
# and repeat with {ids} (comma-separated scenario ids), {repeat} (the repeat index to run: 0 cold, >= 1 warm),
# {repeats} ({repeat} + 1), {out} (RUN_DIR/arms/mcp) and {manifest} (the run's manifest) replaced by
# shell-quoted values; run with `bash -c` and this environment unchanged. Only the template is recorded,
# never the expanded command or any environment value. A full run (not --smoke) uses the default template
# DEFAULT_MCP_CMD below (SUITE_MCP_CMD replaces it, '' for none) unless --mcp-cmd, --no-mcp or --mcp-run is
# given; --mcp-cmd default asks for it in a smoke run; --resume keeps the run's template. The default runner
# needs METAGRAPH_ADMIN_KEY_FILE (a file holding the staging admin key) in the environment, e.g.
#   METAGRAPH_ADMIN_KEY_FILE=/Users/raetsch/tmp/mgss/secrets/staging_admin_key run_suite.sh
# and judges with checks.py of this directory (--checks-dir).
#
# A run directory is runs/<UTC>_<release-short>/ under this directory (SUITE_RUNS overrides). Every
# invocation starts a new run directory, except --resume RUN_DIR, which continues that run (its mode unless
# --paired / --sequential is given; its manifest copy unless --manifest is given): each arm skips what it
# recorded. Ctrl-C stops the step in flight cleanly (a mex runner gets SIGINT), pulls what was recorded and
# leaves the run resumable.
#
# Every invocation first takes the lock ~/suite-dev/LOCK on mex (pid, host, start time; a lock whose process
# is gone is removed) and refuses to start while any run_direct.py runs under ~/suite-dev. --max-requests N
# caps the direct arms' requests of this invocation in total (closing capabilities GETs included);
# --smoke defaults it to 40 and the manifest to manifest/suite_smoke_v1.json. checks.py next to this script is
# used when present (--checks checks); without it the stub runs only for manifests that expect nothing but
# the HTTP status, unless --allow-stub.
#
# On mex only ~/suite-dev/ is used (scripts, manifest, runs, LOCK) and docker is only read (docker stats).
# Environment: SUITE_MEX (srv-metagraph@mex.inf.ethz.ch; local:DIR runs the mex side in DIR, for tests),
# SUITE_MEX_DIR (suite-dev), SUITE_MEX_PY (python3), SUITE_LOCAL_BASE (http://localhost:42010),
# SUITE_REMOTE_BASE (the public refseq33m-experimental URL), SUITE_RUNS, PYTHON (python3), SUITE_POLL_S (15),
# SUITE_SAMPLER (1; 0: no mex sampler besides direct-local's own in sequential mode), SUITE_WAIT_RUNNER_S
# (150), SUITE_SSH_CONTROL_DIR (/tmp/suite-ssh-<uid>).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3}"
MEX="${SUITE_MEX:-srv-metagraph@mex.inf.ethz.ch}"
MEX_DIR="${SUITE_MEX_DIR:-suite-dev}"
MEX_PY="${SUITE_MEX_PY:-python3}"
LOCAL_BASE="${SUITE_LOCAL_BASE:-http://localhost:42010}"
REMOTE_BASE="${SUITE_REMOTE_BASE:-https://metagraph.ethz.ch:8080/api/refseq33m-experimental}"
RUNS="${SUITE_RUNS:-$HERE/runs}"
POLL_S="${SUITE_POLL_S:-15}"
SAMPLER="${SUITE_SAMPLER:-${SUITE_REMOTE_SAMPLER:-1}}"
WAIT_RUNNER_S="${SUITE_WAIT_RUNNER_S:-150}"
CTL_DIR="${SUITE_SSH_CONTROL_DIR:-/tmp/suite-ssh-$(id -u)}"
MEX_LOCAL=""
if [[ "$MEX" == local:* ]]; then MEX_LOCAL="${MEX#local:}"; fi
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=30 -o LogLevel=ERROR)
if [[ -z "$MEX_LOCAL" ]]; then
  # One persistent connection for every ssh / rsync of this invocation (closed on exit).
  mkdir -p "$CTL_DIR" && chmod 700 "$CTL_DIR"
  SSH_OPTS+=(-o ControlMaster=auto -o "ControlPath=$CTL_DIR/%C" -o ControlPersist=300)
fi
# The MCP arm of a full run: the async session's runner (search-service scripts/bench/suite_mcp_runner.py, branch
# feat/suite-mcp-runner), judging with this directory's checks.py; --check-queue pauses while the service queue
# holds other work (the Logan chunks aside) and needs the admin key from the environment.
ASYNC_MCP_CMD="/Users/raetsch/git/services/metagraph-search-service/.venv/bin/python /Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_runner.py --manifest {manifest} --ids {ids} --repeat {repeat} --out {out} --check-queue --queue-ignore sra-logan-chunks --checks-dir $(printf '%q' "$HERE")"
DEFAULT_MCP_CMD="${SUITE_MCP_CMD-$ASYNC_MCP_CMD}"
# The mex runners of the suite, as pgrep / pkill patterns ([s] keeps the pattern from matching its own shell).
RUNNER_PATTERN='[s]cripts/run_direct\.py'

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }
log() { printf '[run_suite %s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }
die() { log "error: $*"; exit 2; }
quote() { local out="" a; for a in "$@"; do out+="$(printf '%q' "$a") "; done; printf '%s' "$out"; }
# A shell command on mex, in its home directory (stderr dropped, no stdin).
mex() {
  if [[ -n "$MEX_LOCAL" ]]; then (cd "$MEX_LOCAL" && bash -c "$*") 2>/dev/null < /dev/null
  else ssh "${SSH_OPTS[@]}" "$MEX" "$*" 2>/dev/null < /dev/null; fi
}
# The same with stderr kept (a runner step: its own stderr goes to runner.log on mex).
mex_step() {
  if [[ -n "$MEX_LOCAL" ]]; then (cd "$MEX_LOCAL" && bash -c "$*") < /dev/null
  else ssh "${SSH_OPTS[@]}" "$MEX" "$*" < /dev/null; fi
}
mex_path() { if [[ -n "$MEX_LOCAL" ]]; then printf '%s/%s' "$MEX_LOCAL" "$1"; else printf '%s:%s' "$MEX" "$1"; fi; }
mex_rsync() {
  if [[ -n "$MEX_LOCAL" ]]; then rsync -a "$@"; else rsync -a -e "ssh ${SSH_OPTS[*]}" "$@"; fi
}

MANIFEST=""; SMOKE=0; AREAS=""; LIMIT=""; SKIP_LOCAL=0; SKIP_REMOTE=0; MCP_RUN=""; MCP_CMD=""; COMPARE=""
RESUME=""; BY_DESIGN=""; NO_PDF=0; PAUSE="1"; MAX_REQUESTS=""; MODE=""; ALLOW_STUB=0; WARM_EVERY=10; NO_MCP=0
MCP_GIVEN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --smoke) SMOKE=1; shift ;;
    --paired) MODE=paired; shift ;;
    --sequential) MODE=sequential; shift ;;
    --areas) AREAS="$2"; shift 2 ;;
    --limit) LIMIT="$2"; shift 2 ;;
    --skip-local) SKIP_LOCAL=1; shift ;;
    --skip-remote) SKIP_REMOTE=1; shift ;;
    --mcp-run) MCP_RUN="$2"; shift 2 ;;
    --mcp-cmd) MCP_CMD="$2"; MCP_GIVEN=1; shift 2 ;;
    --no-mcp) NO_MCP=1; MCP_CMD=""; shift ;;
    --compare) COMPARE="$2"; shift 2 ;;
    --resume) RESUME="$2"; shift 2 ;;
    --by-design) BY_DESIGN="$2"; shift 2 ;;
    --pause) PAUSE="$2"; shift 2 ;;
    --max-requests) MAX_REQUESTS="$2"; shift 2 ;;
    --warm-every) WARM_EVERY="$2"; shift 2 ;;
    --allow-stub) ALLOW_STUB=1; shift ;;
    --no-pdf) NO_PDF=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown option $1" ;;
  esac
done
# The smoke load rule: at most ~40 requests per smoke run.
if [[ $SMOKE == 1 && -z "$MAX_REQUESTS" ]]; then MAX_REQUESTS=40; fi
[[ -z "$MCP_RUN" || -d "$MCP_RUN" ]] || die "--mcp-run $MCP_RUN is not a directory"
[[ -z "$MCP_RUN" || -z "$MCP_CMD" ]] || die "--mcp-run and --mcp-cmd exclude each other"
[[ $NO_MCP == 0 || $MCP_GIVEN == 0 ]] || die "--no-mcp and --mcp-cmd exclude each other"
[[ "$MCP_CMD" != default ]] || MCP_CMD="$DEFAULT_MCP_CMD"
[[ -z "$COMPARE" || -d "$COMPARE" ]] || die "--compare $COMPARE is not a directory"

# --- the run directory ------------------------------------------------------------------------------------
if [[ -n "$RESUME" ]]; then
  [[ -d "$RESUME" ]] || die "--resume $RESUME is not a directory"
  RUN_DIR="$(cd "$RESUME" && pwd)"
  RUN_ID="$(basename "$RUN_DIR")"
  if [[ -z "$MANIFEST" ]]; then
    for cand in "$RUN_DIR/manifest.json" "$RUN_DIR/arms/direct-local/manifest.json" \
                "$RUN_DIR/arms/direct-remote/manifest.json"; do
      if [[ -f "$cand" ]]; then MANIFEST="$cand"; break; fi
    done
  fi
  if [[ -z "$MODE" ]]; then
    MODE="$("$PY" -c 'import json, sys
try:
    print(json.load(open(sys.argv[1]))["sessions"][-1].get("mode") or "")
except (OSError, ValueError, KeyError, IndexError):
    print("")' "$RUN_DIR/orchestrator.json" 2>/dev/null || true)"
  fi
  if [[ $MCP_GIVEN == 0 && $NO_MCP == 0 && -z "$MCP_RUN" ]]; then
    # The run's own MCP template (the last invocation's), so a resume pairs the same arms.
    MCP_CMD="$("$PY" -c 'import json, sys
try:
    print(json.load(open(sys.argv[1]))["sessions"][-1].get("mcp_cmd") or "")
except (OSError, ValueError, KeyError, IndexError):
    print("")' "$RUN_DIR/orchestrator.json" 2>/dev/null || true)"
  fi
  log "resuming $RUN_ID"
else
  if [[ -z "$MANIFEST" ]]; then
    if [[ $SMOKE == 1 ]]; then MANIFEST="$HERE/manifest/suite_smoke_v1.json"; else MANIFEST="$HERE/manifest/suite_v1.json"; fi
  fi
  [[ -n "$MANIFEST" && -f "$MANIFEST" ]] || die "no manifest (${MANIFEST:-none}); pass --manifest"
  log "release: one GET /capabilities on $REMOTE_BASE (names the run directory)"
  RELEASE="$("$PY" - "$REMOTE_BASE" <<'PYEOF' 2>/dev/null || true
import json, sys, urllib.request
with urllib.request.urlopen(sys.argv[1].rstrip("/") + "/capabilities", timeout=30) as r:
    print(json.loads(r.read().decode("utf-8")).get("release") or "")
PYEOF
)"
  RUN_ID="$("$PY" "$HERE/suite_common.py" stamp)_$("$PY" "$HERE/suite_common.py" release-short "${RELEASE:-}")"
  RUN_DIR="$RUNS/$RUN_ID"
  [[ ! -e "$RUN_DIR" ]] || die "$RUN_DIR exists already"
fi
MODE="${MODE:-paired}"
[[ "$MODE" == paired || "$MODE" == sequential ]] || die "unknown mode $MODE"
if [[ -z "$RESUME" && $SMOKE == 0 && $MCP_GIVEN == 0 && $NO_MCP == 0 && -z "$MCP_RUN" ]]; then
  MCP_CMD="$DEFAULT_MCP_CMD"   # a full run includes the MCP arm
fi
if [[ -n "$MCP_CMD" ]]; then
  [[ "$MCP_CMD" == *"{ids}"* && "$MCP_CMD" == *"{out}"* ]] || die "--mcp-cmd needs the placeholders {ids} and {out}"
  MCP_EXE="${MCP_CMD%% *}"
  if [[ "$MCP_EXE" == /* && ! -x "$MCP_EXE" ]]; then
    die "the MCP runner's interpreter $MCP_EXE is missing (--no-mcp runs without the MCP arm, --mcp-cmd names another)"
  fi
  if [[ "$MCP_CMD" == *--check-queue* && -z "${METAGRAPH_ADMIN_KEY_FILE:-}" && -z "${METAGRAPH_ADMIN_KEY:-}" ]]; then
    die "the MCP runner checks the service queue and needs METAGRAPH_ADMIN_KEY_FILE (or METAGRAPH_ADMIN_KEY) in the environment"
  fi
fi
[[ "$MODE" == sequential || -z "$MCP_RUN" ]] || die "--mcp-run reads a finished MCP run: use it with --sequential (or --mcp-cmd)"
[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || die "no manifest (${MANIFEST:-none}); pass --manifest"
MANIFEST="$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")"

# --- the checks (review M3) -----------------------------------------------------------------------------------
if [[ -f "$HERE/checks.py" ]]; then CHECKS=checks; else CHECKS=stub; fi
if [[ $CHECKS == stub && $ALLOW_STUB == 0 ]]; then
  BEYOND="$("$PY" "$HERE/suite_common.py" beyond-status "$MANIFEST" "$AREAS" "$LIMIT")"
  [[ -z "$BEYOND" ]] || die "no checks.py next to this script, and the selected scenarios expect checks the stub cannot evaluate (${BEYOND:0:200}); add lane S1's checks.py, or pass --allow-stub for a run that checks the HTTP status only"
fi

mkdir -p "$RUN_DIR/arms"
COMMIT="$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
DIRTY="$(git -C "$HERE" status --porcelain -- . 2>/dev/null | grep -v '^?? runs/' | head -1 || true)"
[[ -z "$DIRTY" ]] || COMMIT="$COMMIT+dirty"
log "run directory $RUN_DIR ($MODE)"
log "manifest $MANIFEST; scripts $COMMIT; checks $CHECKS"

ARGS=(--manifest "$MANIFEST" --pause "$PAUSE" --scripts-commit "$COMMIT" --run-id "$RUN_ID" --checks "$CHECKS"
      --warm-every "$WARM_EVERY")
[[ -z "$AREAS" ]] || ARGS+=(--areas "$AREAS")
[[ -z "$LIMIT" ]] || ARGS+=(--limit "$LIMIT")
[[ $ALLOW_STUB == 0 ]] || ARGS+=(--allow-stub)
ARM_REL="runs/$RUN_ID/arms/direct-local"
LOCAL_ARGS=("${ARGS[@]}")
LOCAL_ARGS[1]="manifest/$RUN_ID.json"

# --- state for the traps --------------------------------------------------------------------------------------
LOCK_HELD=0; LOCK_TOKEN=""; SAMPLER_PID=""; HB_REL=""; REMOTE_PID=""; LOOP_STARTED=0; CLEANED=0
NEW_RUN=1; [[ -z "$RESUME" ]] || NEW_RUN=0
ARMS_STARTED=0

stop_sampler() {
  [[ -n "$SAMPLER_PID" ]] || return 0
  mex "kill -TERM $SAMPLER_PID" || true
  local i
  for i in 1 2 3 4 5 6; do
    mex "kill -0 $SAMPLER_PID" || break
    sleep 1
  done
  SAMPLER_PID=""
}
# The samplers' files of this run on mex (every invocation's; merge_run deduplicates): review M5.
pull_samples() {
  if mex "test -d $MEX_DIR/runs/$RUN_ID/samples"; then
    mkdir -p "$RUN_DIR/samples"
    mex_rsync --exclude '.heartbeat' "$(mex_path "$MEX_DIR/runs/$RUN_ID/samples/")" "$RUN_DIR/samples/" 2>/dev/null \
      || log "pulling the samples from $MEX failed (they stay there; the next invocation pulls them)"
  fi
}
pull_local_arm() {
  if mex "test -d $MEX_DIR/$ARM_REL"; then
    mkdir -p "$RUN_DIR/arms/direct-local"
    mex_rsync "$(mex_path "$MEX_DIR/$ARM_REL/")" "$RUN_DIR/arms/direct-local/" 2>/dev/null \
      || log "pulling direct-local from $MEX failed"
  fi
}
release_lock() {
  [[ $LOCK_HELD == 1 ]] || return 0
  mex "grep -q $LOCK_TOKEN $MEX_DIR/LOCK && rm -f $MEX_DIR/LOCK" || true
  LOCK_HELD=0
}
close_master() {
  [[ -z "$MEX_LOCAL" ]] || return 0
  ssh -o "ControlPath=$CTL_DIR/%C" -O exit "$MEX" >/dev/null 2>&1 || true
}
cleanup() {
  [[ $CLEANED == 0 ]] || return 0
  CLEANED=1
  stop_sampler
  pull_samples
  if [[ $NEW_RUN == 1 && $ARMS_STARTED == 0 ]]; then
    # Refused before any arm ran (lock, runner, options): leave no empty run directory behind.
    mex "rmdir $MEX_DIR/runs/$RUN_ID/arms $MEX_DIR/runs/$RUN_ID" || true
    rmdir "$RUN_DIR/arms" "$RUN_DIR" 2>/dev/null || true
  fi
  release_lock
  close_master
}
orch_set() {  # orch_set KEY JSON [KEY JSON ...]: fields of this invocation's session in orchestrator.json
  "$PY" - "$RUN_DIR/orchestrator.json" "$@" <<'PYEOF' || true
import json, os, sys
path, args = sys.argv[1], sys.argv[2:]
with open(path) as fh:
    data = json.load(fh)
session = data["sessions"][-1]
for key, value in zip(args[::2], args[1::2]):
    session[key] = json.loads(value)
with open(path + ".tmp", "w") as fh:
    json.dump(data, fh, indent=1, sort_keys=True)
os.replace(path + ".tmp", path)
PYEOF
}
# Waits until no run_direct.py of this suite runs on mex; 1 when one still does after `$1` seconds.
wait_no_runner() {
  local waited=0 pids
  while true; do
    pids="$(mex "pgrep -u \$(id -u) -f '$RUNNER_PATTERN' | tr '\n' ' '" || true)"
    [[ -n "${pids// /}" ]] || return 0
    (( waited < $1 )) || { log "run_direct.py still running on $MEX (pid(s) $pids)"; return 1; }
    sleep 5
    waited=$((waited + 5))
  done
}
interrupted() {
  trap - INT TERM
  log "interrupted"
  if [[ $LOCK_HELD == 1 ]]; then
    # The runner of the step in flight (paired) or of the arm (sequential) stops after its request.
    mex "pkill -INT -u \$(id -u) -f '$RUNNER_PATTERN .*--run-id $RUN_ID'" || true
    if ! wait_no_runner 120; then
      log "the mex runner still runs (a long request may take up to 30 min); what it recorded so far is pulled, the rest stays on mex; --resume waits for it or attaches to it"
    fi
    pull_local_arm
  fi
  [[ $LOOP_STARTED == 0 ]] || orch_set loop "{\"mode\": \"$MODE\", \"state\": \"interrupted\", \"scenarios\": ${N_PLAN:-0}, \"completed\": ${COMPLETED:-0}}"
  cleanup
  log "resume with: $0 --resume $RUN_DIR"
  exit 130
}
trap interrupted INT TERM
trap cleanup EXIT

# --- mex: reachable, clock, lock, no other runner, scripts -----------------------------------------------------
mex "mkdir -p $MEX_DIR/scripts $MEX_DIR/manifest $MEX_DIR/runs/$RUN_ID/arms" || die "ssh $MEX failed"
T_A="$("$PY" -c 'import time; print("%.3f" % time.time())')"
T_MEX="$(mex "$MEX_PY -c 'import time; print(\"%.3f\" % time.time())'" || true)"
T_B="$("$PY" -c 'import time; print("%.3f" % time.time())')"
CLOCK_OFFSET="$("$PY" "$HERE/suite_common.py" clock-offset "$T_A" "$T_MEX" "$T_B")"
log "mex clock offset ${CLOCK_OFFSET:-unknown} s"

take_lock() {
  local host now_utc now_epoch json held verdict attempt
  host="$(hostname)"
  LOCK_TOKEN="$("$PY" -c 'import uuid; print(uuid.uuid4().hex)')"
  now_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  now_epoch="$("$PY" -c 'import time; print("%.3f" % time.time())')"
  json="$(printf '{"pid": %d, "host": "%s", "started_utc": "%s", "started_epoch": %s, "run_id": "%s", "mode": "%s", "token": "%s"}' \
    "$$" "$host" "$now_utc" "$now_epoch" "$RUN_ID" "$MODE" "$LOCK_TOKEN")"
  for attempt in 1 2 3; do
    if mex "(set -o noclobber; printf '%s\n' $(quote "$json") > $MEX_DIR/LOCK)"; then
      LOCK_HELD=1
      log "lock ~/$MEX_DIR/LOCK taken on $MEX"
      return 0
    fi
    held="$(mex "cat $MEX_DIR/LOCK" || true)"
    [[ -n "$held" ]] || continue
    verdict="$("$PY" "$HERE/suite_common.py" lock-check "$held" "$host")"
    case "$verdict" in
      stale*) log "removing a stale lock (${verdict#stale })"; mex "rm -f $MEX_DIR/LOCK" || true ;;
      *) die "the suite lock ~/$MEX_DIR/LOCK on $MEX is ${verdict#* }; another invocation drives the host (remove the file only if that process is gone)" ;;
    esac
  done
  die "could not take the lock ~/$MEX_DIR/LOCK on $MEX"
}
take_lock

# Any run_direct.py under ~/suite-dev (an earlier invocation's, cut off): never start another arm next to it.
# A sequential --resume attaches to its own direct-local runner instead.
ATTACH_PID=""
if [[ $MODE == sequential && -n "$RESUME" && $SKIP_LOCAL == 0 ]]; then
  LIVE="$(mex "cat $MEX_DIR/$ARM_REL/state.json 2>/dev/null" | "$PY" -c '
import json, sys
try:
    s = json.load(sys.stdin)
except ValueError:
    s = {}
print(s.get("pid", "") if s.get("state") == "running" else "")' 2>/dev/null || true)"
  PIDS="$(mex "pgrep -u \$(id -u) -f '$RUNNER_PATTERN' | tr '\n' ' '" || true)"
  if [[ "$LIVE" =~ ^[0-9]+$ && "${PIDS// /}" == "$LIVE" ]]; then ATTACH_PID="$LIVE"; fi
fi
if [[ -z "$ATTACH_PID" ]] && ! wait_no_runner "$WAIT_RUNNER_S"; then
  die "a run_direct.py runs under ~/$MEX_DIR on $MEX (an earlier invocation's); wait for it, stop it with kill -INT, or --resume its run"
fi

mex_rsync --delete --exclude 'runs/' --exclude 'tests/' --exclude '__pycache__/' --exclude '*.pyc' \
  "$HERE/" "$(mex_path "$MEX_DIR/scripts/")" || die "rsync of the scripts to $MEX failed"
mex_rsync "$MANIFEST" "$(mex_path "$MEX_DIR/manifest/$RUN_ID.json")" || die "rsync of the manifest to $MEX failed"

"$PY" - "$RUN_DIR/orchestrator.json" "$MANIFEST" "$COMMIT" "$SMOKE" "$AREAS" "$LIMIT" "$SKIP_LOCAL" \
    "$SKIP_REMOTE" "$MCP_RUN" "$MCP_CMD" "$COMPARE" "$PAUSE" "$MAX_REQUESTS" "$LOCAL_BASE" "$REMOTE_BASE" "$MEX" \
    "$MODE" "$CHECKS" "$ALLOW_STUB" "$WARM_EVERY" "$CLOCK_OFFSET" "$(hostname)" <<'PYEOF'
import json, os, sys, time
path = sys.argv[1]
keys = ["manifest", "scripts_commit", "smoke", "areas", "limit", "skip_local", "skip_remote", "mcp_run",
        "mcp_cmd", "compare", "pause_s", "max_requests", "local_base", "remote_base", "mex", "mode", "checks",
        "allow_stub", "warm_every", "mex_clock_offset_s", "host"]
session = dict(zip(keys, sys.argv[2:]))
try:
    session["mex_clock_offset_s"] = float(session["mex_clock_offset_s"])
except ValueError:
    session["mex_clock_offset_s"] = None
session["started_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
data = {"sessions": []}
if os.path.exists(path):
    with open(path) as fh:
        data = json.load(fh)
data["sessions"].append(session)
with open(path + ".tmp", "w") as fh:
    json.dump(data, fh, indent=1, sort_keys=True)
os.replace(path + ".tmp", path)
PYEOF
pull_samples

# The direct-local arm's directory on mex: pushed back when only the local copy survived, else pulled (a resume
# plans from it).
prepare_local_arm() {
  if [[ -d "$RUN_DIR/arms/direct-local" ]] && ! mex "test -f $MEX_DIR/$ARM_REL/run.json"; then
    log "the mex copy of direct-local is gone: pushing the local one back"
    mex "mkdir -p $MEX_DIR/$ARM_REL"
    mex_rsync "$RUN_DIR/arms/direct-local/" "$(mex_path "$MEX_DIR/$ARM_REL/")"
  else
    mex "mkdir -p $MEX_DIR/$ARM_REL"
    pull_local_arm
  fi
}
start_sampler() {  # start_sampler LABEL MAX_SECONDS [heartbeat]
  [[ $SAMPLER == 1 ]] || return 0
  local rel cmd hb_args=()
  rel="runs/$RUN_ID/samples/$1-$("$PY" "$HERE/suite_common.py" stamp).jsonl"
  if [[ "${3:-}" == heartbeat ]]; then
    HB_REL="runs/$RUN_ID/samples/.heartbeat"
    hb_args=(--heartbeat "$HB_REL" --heartbeat-max-age 3600)
  fi
  cmd="$(quote "$MEX_PY" scripts/sampler.py --out "$rel" --label "$1" --max-seconds "$2" ${hb_args[@]+"${hb_args[@]}"})"
  SAMPLER_PID="$(mex "cd $MEX_DIR && mkdir -p runs/$RUN_ID/samples && touch runs/$RUN_ID/samples/.heartbeat && (nohup $cmd > /dev/null 2>&1 < /dev/null & echo \$!)" || true)"
  if [[ "$SAMPLER_PID" =~ ^[0-9]+$ ]]; then
    log "$1: sampler on $MEX (pid $SAMPLER_PID)"
  else
    log "$1: no sampler on $MEX (got '$SAMPLER_PID'); this phase has no resource samples"
    SAMPLER_PID=""
  fi
}
heartbeat() { [[ -z "$SAMPLER_PID" || -z "$HB_REL" ]] || mex "touch $MEX_DIR/$HB_REL" || true; }

# The MCP runner (async's), per scenario ids and repeat; the environment passes through unchanged.
mcp_step() {  # mcp_step IDS REPEAT
  local cmd rc=0
  cmd="$("$PY" - "$MCP_CMD" "$1" "$2" "$RUN_DIR/arms/mcp" "$MANIFEST" <<'PYEOF'
import shlex, sys
template, ids, repeat, out, manifest = sys.argv[1:6]
values = {"ids": ids, "repeat": repeat, "repeats": str(int(repeat) + 1), "out": out, "manifest": manifest}
for key, value in values.items():
    template = template.replace("{%s}" % key, shlex.quote(value))
print(template)
PYEOF
)"
  mkdir -p "$RUN_DIR/arms/mcp"
  log "mcp: $( [[ $2 == 0 ]] && echo cold || echo warm ) repeat $2 of ${1:0:120}"
  bash -c "$cmd" < /dev/null >&2 || rc=$?
  sleep "$PAUSE"
  return $rc
}

USED=0          # the direct arms' requests in this invocation
RESERVE=0       # requests kept back for the closing GETs (paired mode: one per direct arm)
# One run_direct.py invocation; its progress goes to the log, its request count into USED.
direct_call() {  # direct_call local|remote EXTRA_ARGS...
  local arm="$1" out="" rc=0 status requests cap_args=()
  shift
  if [[ -n "$MAX_REQUESTS" ]]; then
    local left=$((MAX_REQUESTS - USED - RESERVE))
    (( left > 0 )) || return 3
    cap_args=(--max-requests "$left")
  fi
  if [[ $arm == local ]]; then
    out="$(mex_step "cd $MEX_DIR && $(quote "$MEX_PY" scripts/run_direct.py --arm direct-local --base "$LOCAL_BASE" \
      --out "$ARM_REL" "${LOCAL_ARGS[@]}" "$@" --status-line ${cap_args[@]+"${cap_args[@]}"}) 2>> $ARM_REL/runner.log")" || rc=$?
  else
    out="$("$PY" "$HERE/run_direct.py" --arm direct-remote --base "$REMOTE_BASE" --out "$RUN_DIR/arms/direct-remote" \
      "${ARGS[@]}" "$@" --status-line ${cap_args[@]+"${cap_args[@]}"} < /dev/null)" || rc=$?
  fi
  printf '%s\n' "$out" | grep -v '^RUN_DIRECT_STATUS ' | sed '/^$/d; s/^/    /' >&2 || true
  status="$(printf '%s\n' "$out" | sed -n 's/^RUN_DIRECT_STATUS //p' | tail -1)"
  requests="$(printf '%s' "$status" | sed -n 's/.*"requests": \([0-9][0-9]*\).*/\1/p')"
  USED=$((USED + ${requests:-0}))
  if [[ -z "$status" && $rc == 0 ]]; then rc=1; fi
  return $rc
}

FINAL_RC=0
# --- paired mode ----------------------------------------------------------------------------------------------
run_paired() {
  local arms="" entry sid l r mc mw rep rc fails_local=0 fails_remote=0 fails_mcp=0 since_pull=0
  [[ $SKIP_LOCAL == 1 ]] || { arms+="direct-local,"; RESERVE=$((RESERVE + 1)); }
  [[ $SKIP_REMOTE == 1 ]] || { arms+="direct-remote,"; RESERVE=$((RESERVE + 1)); }
  [[ -z "$MCP_CMD" ]] || arms+="mcp,"
  [[ -n "$arms" ]] || die "nothing to run (both direct arms skipped and no --mcp-cmd)"
  [[ $SKIP_LOCAL == 1 ]] || prepare_local_arm
  local plan=()
  while IFS= read -r entry; do plan+=("$entry"); done < <("$PY" "$HERE/suite_common.py" paired-plan "$MANIFEST" \
    "$RUN_DIR" "$AREAS" "$LIMIT" "$WARM_EVERY" "$arms")
  N_PLAN=${#plan[@]}
  COMPLETED=0
  local max_s=86400; [[ $SMOKE == 0 ]] || max_s=1800
  start_sampler paired "$max_s" heartbeat
  LOOP_STARTED=1
  LOOP_STATE=running
  orch_set loop "{\"mode\": \"paired\", \"state\": \"running\", \"scenarios\": $N_PLAN, \"completed\": 0}"
  log "paired: $N_PLAN scenario(s); arms ${arms%,}"
  for entry in ${plan[@]+"${plan[@]}"}; do
    read -r sid l r mc mw <<< "$entry"
    heartbeat
    if [[ $mc == 1 ]]; then
      rc=0; mcp_step "$sid" 0 || rc=$?
      if [[ $rc != 0 ]]; then fails_mcp=$((fails_mcp + 1)); log "mcp: exit $rc on $sid"; else fails_mcp=0; fi
    fi
    if [[ $l == 1 ]]; then
      log "direct-local: $sid"
      rc=0; direct_call local --only "$sid" || rc=$?
      if [[ $rc == 3 ]]; then LOOP_STATE=capped; break; fi
      if [[ $rc != 0 ]]; then
        fails_local=$((fails_local + 1)); log "direct-local: exit $rc on $sid (see arms/direct-local/runner.log)"
        # An ssh that failed may have left the step running: never start the next arm next to it.
        wait_no_runner "$WAIT_RUNNER_S" || { LOOP_STATE=failed; break; }
      else
        fails_local=0
      fi
    fi
    if [[ $r == 1 ]]; then
      log "direct-remote: $sid"
      rc=0; direct_call remote --only "$sid" || rc=$?
      if [[ $rc == 3 ]]; then LOOP_STATE=capped; break; fi
      if [[ $rc != 0 ]]; then fails_remote=$((fails_remote + 1)); log "direct-remote: exit $rc on $sid"; else fails_remote=0; fi
    fi
    if [[ "$mw" != "-" ]]; then
      for rep in ${mw//,/ }; do
        rc=0; mcp_step "$sid" "$rep" || rc=$?
        if [[ $rc != 0 ]]; then fails_mcp=$((fails_mcp + 1)); log "mcp: exit $rc on $sid"; else fails_mcp=0; fi
      done
    fi
    if (( fails_local >= 3 || fails_remote >= 3 || fails_mcp >= 3 )); then
      log "three failing steps of one arm in a row: stopping the loop"
      LOOP_STATE=failed
      break
    fi
    COMPLETED=$((COMPLETED + 1))
    since_pull=$((since_pull + 1))
    if [[ $SKIP_LOCAL == 0 && $since_pull -ge 25 ]]; then pull_local_arm; since_pull=0; fi
  done
  [[ $LOOP_STATE != running ]] || LOOP_STATE=done
  [[ $LOOP_STATE != capped ]] || log "paired: stopped at --max-requests $MAX_REQUESTS ($USED direct requests); --resume continues"
  # The closing GET /capabilities of each direct arm (one request kept back for each) and its state over the
  # whole selection.
  if [[ $SKIP_LOCAL == 0 ]]; then
    RESERVE=$((RESERVE - 1))
    direct_call local --finish || log "direct-local: the closing step failed"
  fi
  if [[ $SKIP_REMOTE == 0 ]]; then
    RESERVE=$((RESERVE - 1))
    direct_call remote --finish || log "direct-remote: the closing step failed"
  fi
  [[ $SKIP_LOCAL == 1 ]] || pull_local_arm
  stop_sampler
  orch_set loop "{\"mode\": \"paired\", \"state\": \"$LOOP_STATE\", \"scenarios\": $N_PLAN, \"completed\": $COMPLETED, \"direct_requests\": $USED}"
  log "paired: $LOOP_STATE, $COMPLETED/$N_PLAN scenario(s), $USED direct request(s)"
  [[ $LOOP_STATE != failed ]] || FINAL_RC=1
}

# --- sequential mode ------------------------------------------------------------------------------------------
run_sequential() {
  local rc state requests
  if [[ $SKIP_LOCAL == 0 ]]; then
    prepare_local_arm
    local cap_args=()
    [[ -z "$MAX_REQUESTS" ]] || cap_args=(--max-requests "$MAX_REQUESTS")
    if [[ -n "$ATTACH_PID" ]]; then
      REMOTE_PID="$ATTACH_PID"
      log "direct-local: attaching to the runner still running on $MEX (pid $REMOTE_PID)"
    else
      local cmd
      cmd="$(quote "$MEX_PY" scripts/run_direct.py --arm direct-local --base "$LOCAL_BASE" --out "$ARM_REL" \
        --sample-docker "${LOCAL_ARGS[@]}" ${cap_args[@]+"${cap_args[@]}"})"
      log "direct-local: starting on $MEX against $LOCAL_BASE"
      REMOTE_PID="$(mex "cd $MEX_DIR && mkdir -p $ARM_REL && (nohup $cmd >> $ARM_REL/runner.log 2>&1 < /dev/null & echo \$!)" || true)"
      if [[ ! "$REMOTE_PID" =~ ^[0-9]+$ ]]; then
        # The runner may have started even so (a stray line in the ssh output): find it before giving up.
        REMOTE_PID="$(mex "pgrep -u \$(id -u) -f '$RUNNER_PATTERN .*--run-id $RUN_ID' | head -1" || true)"
        [[ "$REMOTE_PID" =~ ^[0-9]+$ ]] || die "could not start the runner on $MEX"
      fi
    fi
    log "direct-local: mex pid $REMOTE_PID; polling every ${POLL_S}s"
    local silent=0 st progress
    while true; do
      # An explicit alive / dead answer: a failed ssh says neither, and must never read as "finished" (the
      # remote arm would then run while direct-local still does).
      st="$(mex "if kill -0 $REMOTE_PID 2>/dev/null; then echo alive; else echo dead; fi; cat $MEX_DIR/$ARM_REL/state.json 2>/dev/null" || true)"
      case "$st" in
        alive*|dead*) silent=0 ;;
        *)
          silent=$((silent + 1))
          [[ $silent -lt 40 ]] || die "no answer from $MEX for $((silent * POLL_S)) s; the runner may still run there (pid $REMOTE_PID): --resume $RUN_DIR once it is reachable"
          log "direct-local: no answer from $MEX ($silent); retrying"
          sleep "$POLL_S"
          continue ;;
      esac
      progress="$(printf '%s' "$st" | "$PY" -c '
import json, sys
text = sys.stdin.read()
i = text.find("{")
try:
    s = json.loads(text[i:]) if i >= 0 else {}
except ValueError:
    s = {}
print("%s %s/%s %s" % (s.get("state", "?"), s.get("done", "?"), s.get("planned", "?"), s.get("current") or ""))
' 2>/dev/null || true)"
      log "direct-local: $progress"
      case "$st" in alive*) sleep "$POLL_S" ;; *) break ;; esac
    done
    REMOTE_PID=""
    pull_local_arm
    state="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("state"))' \
      "$RUN_DIR/arms/direct-local/state.json" 2>/dev/null || echo unknown)"
    requests="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("requests_session") or 0)' \
      "$RUN_DIR/arms/direct-local/state.json" 2>/dev/null || echo 0)"
    USED=$((USED + requests))
    log "direct-local: $state ($requests request(s))"
    [[ "$state" == "done" ]] || log "direct-local did not finish ($state): see arms/direct-local/runner.log; --resume continues it"
  fi
  if [[ $SKIP_REMOTE == 0 ]]; then
    local max_s=43200; [[ $SMOKE == 0 ]] || max_s=1800
    start_sampler direct-remote "$max_s"
    log "direct-remote: from $(hostname) against $REMOTE_BASE"
    rc=0; direct_call remote || rc=$?
    stop_sampler
    case "$rc" in
      0) log "direct-remote: done" ;;
      3) log "direct-remote: stopped at --max-requests $MAX_REQUESTS (incomplete; --resume continues it)" ;;
      *) log "direct-remote: exit $rc (see above)" ;;
    esac
  fi
  if [[ -n "$MCP_CMD" ]]; then
    # The MCP arm last: repeat 0 of every scenario, then each warm repeat.
    local arms="mcp" entry sid l r mc mw rep cold="" warm1=""
    while IFS= read -r entry; do
      read -r sid l r mc mw <<< "$entry"
      [[ $mc != 1 ]] || cold+="$sid,"
      [[ ",$mw," != *",1,"* ]] || warm1+="$sid,"
    done < <("$PY" "$HERE/suite_common.py" paired-plan "$MANIFEST" "$RUN_DIR" "$AREAS" "$LIMIT" "$WARM_EVERY" "$arms")
    [[ -z "$cold" ]] || mcp_step "${cold%,}" 0 || log "mcp: the cold pass failed"
    [[ -z "$warm1" ]] || mcp_step "${warm1%,}" 1 || log "mcp: the warm pass failed"
  fi
  if [[ -n "$MCP_RUN" ]]; then
    log "mcp: copying $MCP_RUN"
    mkdir -p "$RUN_DIR/arms/mcp"
    rsync -a "$MCP_RUN/" "$RUN_DIR/arms/mcp/"
  fi
}

ARMS_STARTED=1
if [[ $MODE == paired ]]; then run_paired; else run_sequential; fi
cleanup

# --- merge, compare, report -----------------------------------------------------------------------------------
"$PY" "$HERE/merge_run.py" "$RUN_DIR"
CMP_ARGS=("$RUN_DIR" --no-merge)
[[ -z "$BY_DESIGN" ]] || CMP_ARGS+=(--by-design "$BY_DESIGN")
"$PY" "$HERE/compare.py" "${CMP_ARGS[@]}"
REP_ARGS=("$RUN_DIR")
[[ -z "$COMPARE" ]] || REP_ARGS+=(--compare "$COMPARE")
[[ $NO_PDF == 0 ]] || REP_ARGS+=(--no-pdf)
"$PY" "$HERE/report.py" "${REP_ARGS[@]}" || log "report: the PDF was not produced (see above); report.md is there"
log "done: $RUN_DIR"
printf '%s\n' "$RUN_DIR"
exit "$FINAL_RC"
