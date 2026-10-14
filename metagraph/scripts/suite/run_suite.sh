#!/usr/bin/env bash
# One command for a whole suite run (DESIGN §5, §7, §8.5): the arms against the staging host, never two at
# once, then merge, compare and report (md + pdf).
#
#   run_suite.sh [--manifest M | --smoke] [--paired | --sequential]
#                [--mcp-cmd TEMPLATE | --mcp-cmd default | --no-mcp | --mcp-run DIR]
#                [--areas P,M] [--ids FILE] [--limit N] [--skip-local | --local] [--skip-remote | --remote]
#                [--resume RUN_DIR] [--compare EARLIER_RUN_DIR] [--by-design PATH] [--pause S] [--max-requests N]
#                [--warm-every N] [--checks auto|checks|stub] [--allow-stub | --no-allow-stub] [--no-pdf]
#                [--gate POLICY | --no-gate]
#   run_suite.sh --reevaluate RUN_DIR [--reeval-checks DIR] [--reeval-mcp-auth anonymous|api_key]
#                [--compare EARLIER_RUN_DIR] [--by-design PATH] [--no-pdf] [--gate POLICY]
#
# --ids FILE: only these scenarios (ids separated by commas, blanks or new lines, '#' comments; an id the manifest
# lacks is refused), in manifest order, with --areas / --limit applied on top. --resume without --ids, --areas,
# --limit or --gate runs the scenarios of the run's whole plan (every invocation's selection: a plain resume after a
# narrower targeted one finishes the run; a resume that names a narrower selection says so). A resume keeps the
# last invocation's mode, --checks, MCP template, --warm-every and --pause unless given again, and the run's FIRST
# invocation's flags: --skip-local, --skip-remote, --allow-stub, --smoke (its request cap of 40 per invocation and
# its sampler bound) and --gate, so a plain `--resume RUN_DIR` continues the run as it was started and is gated as
# it was; a resume that names one (or its negation --local, --remote, --no-allow-stub, --no-gate) applies it to that
# invocation only. Not --max-requests, which caps one invocation. Every
# invocation writes the run's plan, RUN_DIR/plan.json, before any arm runs: each (scenario, arm, repeat) its
# selection plans (gate.py plan); one that adds executions is added, the run's plan is the union.
#
# --gate POLICY (gate.py, a policy of gate_policy.py): the strict deploy gate. Without --ids / --areas / --limit the
# selection is the policy's (its scenarios less those excluded on every arm this invocation runs); the merge treats
# differing manifest identities as an error (merge_run.py --strict-identity), gate.py judges the run after the
# comparison (RUN_DIR/gate.json, which the report shows) and the exit status is the gate's: 0 PASS, 1 FAIL (2: the
# policy unreadable). With --reevaluate it judges the re-evaluated set (gate.reevaluated.json). Without --gate the
# exit status is the orchestration's, as before: a FAIL in the report still exits 0.
#
# --reevaluate RUN_DIR: offline, no lock and no request to any host: judges the run's recorded answers again with
# checks.py of this directory (or of --reeval-checks DIR) (reevaluate.py), then merge_run / compare / report
# --reevaluated: executions.reevaluated.jsonl, compare.reevaluated.json, report.reevaluated.md / .pdf beside the
# recorded files, which stay. --reeval-mcp-auth: reevaluate.py --mcp-auth, the MCP auth rule 2 assumes for lines
# whose record states none (a run before async's 2a8acbaf); recorded in the re-evaluation's mark.
#
# --paired (the default): per scenario in manifest order, MCP cold (with --mcp-cmd) -> direct-local (one
# run_direct.py step on mex over one persistent ssh connection) -> direct-remote (a step from here through the
# public URL) -> MCP warm; one resource sampler on mex for the whole loop. --sequential: one arm after the
# other (direct-local on mex with its own sampler, then direct-remote, then the MCP arm), for runs that skip
# arms or read a finished MCP run directory (--mcp-run DIR, copied to arms/mcp).
#
# --mcp-cmd TEMPLATE: the async session's MCP runner (scripts/bench/suite_mcp_runner.py), called per scenario
# and repeat with {ids} (comma-separated scenario ids), {repeat} (the repeat index to run: 0 cold, >= 1 warm),
# {repeats} ({repeat} + 1), {out} (RUN_DIR/arms/mcp), {manifest} (the run's manifest) and {checks_dir} (the run's
# snapshot RUN_DIR/scripts, below) replaced by
# shell-quoted values; run with `bash -c` and this environment unchanged. Only the template is recorded,
# never the expanded command or any environment value. A full run (not --smoke) uses the default template
# DEFAULT_MCP_CMD below (SUITE_MCP_CMD replaces it, '' for none) unless --mcp-cmd, --no-mcp or --mcp-run is
# given; --mcp-cmd default asks for it in a smoke run; --resume keeps the run's template. The default runner
# needs METAGRAPH_ADMIN_KEY_FILE (a file holding the staging admin key) in the environment, e.g.
#   METAGRAPH_ADMIN_KEY_FILE=/Users/raetsch/tmp/mgss/secrets/staging_admin_key run_suite.sh
# and judges with the run's snapshot (--checks-dir {checks_dir}).
#
# The run's snapshot: the run's first invocation copies what judges and runs the arms (suite_common.SNAPSHOT_FILES:
# run_direct.py, sampler.py, suite_common.py, checks.py, checks_stub.py, service_by_design.json) into
# RUN_DIR/scripts with their sha256 (SNAPSHOT.json). direct-remote runs from it, it is what goes to mex for
# direct-local, and the MCP runner gets it as {checks_dir}; a resume keeps it. A change of this directory during
# the run, or before a resume, changes nothing for the run (a new run, or --reevaluate, judges with it); the
# runners record the sha256 of the checks they imported, and merge_run.py states every arm judged with others.
#
# A run directory is runs/<UTC>_<release-short>/ under this directory (SUITE_RUNS overrides). Every
# invocation starts a new run directory, except --resume RUN_DIR, which continues that run (its mode unless
# --paired / --sequential is given; its manifest copy unless --manifest is given): each arm skips what it
# recorded. Ctrl-C stops the step in flight cleanly (a mex runner gets SIGINT), pulls what was recorded and
# leaves the run resumable.
#
# Every invocation first takes the lock ~/suite-dev/LOCK on mex (pid, host, start time; a lock whose process
# is gone is removed) and refuses to start while any run_direct.py runs under ~/suite-dev, both before any
# request to the host (a new run's release GET comes after them). --max-requests N
# caps the direct arms' requests of this invocation in total (closing capabilities GETs included);
# --smoke defaults it to 40 and the manifest to manifest/suite_smoke_v4.json. --checks picks the direct arms'
# checks: auto (the default; --resume keeps the run's choice) runs checks.py when it is next to this script,
# else the stub; checks requires checks.py; stub runs checks_stub.py even beside checks.py. The stub runs only
# for manifests that expect nothing but the HTTP status, unless --allow-stub.
#
# On mex only ~/suite-dev/ is used (scripts, manifest, runs, LOCK) and docker is only read (docker stats).
# Environment: SUITE_MEX (srv-metagraph@mex.inf.ethz.ch; local:DIR runs the mex side in DIR, for tests),
# SUITE_MEX_DIR (suite-dev), SUITE_MEX_PY (python3), SUITE_LOCAL_BASE (http://localhost:42010),
# SUITE_REMOTE_BASE (the public refseq33m-experimental URL), SUITE_RUNS, PYTHON (python3), SUITE_POLL_S (15),
# SUITE_SAMPLER (1; 0: no mex sampler besides direct-local's own in sequential mode), SUITE_WAIT_RUNNER_S
# (150), SUITE_WAIT_OWN_RUNNER_S (1860: how long a --resume waits for its own run's runner, the long class's client
# timeout and a minute), SUITE_SSH_CONTROL_DIR (/tmp/suite-ssh-<uid>), SUITE_HOST_ID (the lock's id of this machine; default
# macOS's LocalHostName, else the node name), SUITE_HEARTBEAT_S (300; the paired sampler's heartbeat interval),
# SUITE_RUNNER_PATTERN (the pgrep pattern of the mex runners; the offline tests narrow it to their own).
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
WAIT_OWN_RUNNER_S="${SUITE_WAIT_OWN_RUNNER_S:-1860}"
CTL_DIR="${SUITE_SSH_CONTROL_DIR:-/tmp/suite-ssh-$(id -u)}"
# This invocation's own control socket (review N1): a shared one would let another invocation's exit (`-O exit`
# in its cleanup) end this one's master and with it the ssh step in flight.
CTL_PATH="$CTL_DIR/%C-$$"
MEX_LOCAL=""
if [[ "$MEX" == local:* ]]; then MEX_LOCAL="${MEX#local:}"; fi
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=30 -o LogLevel=ERROR)
if [[ -z "$MEX_LOCAL" ]]; then
  # One persistent connection for every ssh / rsync of this invocation, and only of this one (closed on exit).
  mkdir -p "$CTL_DIR" && chmod 700 "$CTL_DIR"
  SSH_OPTS+=(-o ControlMaster=auto -o "ControlPath=$CTL_PATH" -o ControlPersist=300)
fi
# The MCP arm of a full run: the async session's runner (search-service scripts/bench/suite_mcp_runner.py, branch
# feat/suite-mcp-runner), judging with this directory's checks.py; --check-queue pauses while the service queue
# holds other work (the Logan chunks aside) and needs the admin key from the environment.
ASYNC_MCP_CMD="/Users/raetsch/git/services/metagraph-search-service/.venv/bin/python /Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_runner.py --manifest {manifest} --ids {ids} --repeat {repeat} --out {out} --check-queue --queue-ignore sra-logan-chunks --checks-dir {checks_dir}"
DEFAULT_MCP_CMD="${SUITE_MCP_CMD-$ASYNC_MCP_CMD}"
# The mex runners of the suite, as pgrep / pkill patterns ([s] keeps the pattern from matching its own shell).
# SUITE_RUNNER_PATTERN narrows it for the offline tests, whose "mex" is this machine: there a real run's ssh step
# (its command line holds scripts/run_direct.py) would otherwise read as a runner of the test's own.
RUNNER_PATTERN="${SUITE_RUNNER_PATTERN:-[s]cripts/run_direct\.py}"

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }
log() { printf '[run_suite %s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }
die() { log "error: $*"; exit 2; }
quote() { local out="" a; for a in "$@"; do out+="$(printf '%q' "$a") "; done; printf '%s' "$out"; }
# A shell command on mex, in its home directory (stderr dropped, no stdin).
mex() {
  if [[ -n "$MEX_LOCAL" ]]; then (cd "$MEX_LOCAL" && bash -c "$*") 2>/dev/null < /dev/null
  else ssh "${SSH_OPTS[@]}" "$MEX" "$*" 2>/dev/null < /dev/null; fi
}
# The same for a runner step (its own stderr goes to runner.log on mex): what else the session prints on stderr
# (mex's login banner, an ssh error) is shown only when the step fails.
mex_step() {
  local rc=0 errf="$RUN_DIR/.mex-step.err"
  if [[ -n "$MEX_LOCAL" ]]; then (cd "$MEX_LOCAL" && bash -c "$*") < /dev/null 2> "$errf" || rc=$?
  else ssh "${SSH_OPTS[@]}" "$MEX" "$*" < /dev/null 2> "$errf" || rc=$?; fi
  if [[ $rc != 0 ]]; then tail -5 "$errf" | sed 's/^/    ssh: /' >&2; fi
  rm -f "$errf"
  return $rc
}
mex_path() { if [[ -n "$MEX_LOCAL" ]]; then printf '%s/%s' "$MEX_LOCAL" "$1"; else printf '%s:%s' "$MEX" "$1"; fi; }
mex_rsync() {
  if [[ -n "$MEX_LOCAL" ]]; then rsync -a "$@"; else rsync -a -e "ssh ${SSH_OPTS[*]}" "$@"; fi
}

MANIFEST=""; SMOKE=0; AREAS=""; LIMIT=""; SKIP_LOCAL=0; SKIP_REMOTE=0; MCP_RUN=""; MCP_CMD=""; COMPARE=""
RESUME=""; BY_DESIGN=""; NO_PDF=0; PAUSE="1"; MAX_REQUESTS=""; MODE=""; ALLOW_STUB=0; WARM_EVERY=10; NO_MCP=0
MCP_GIVEN=0; CHECKS_MODE=""; REEVALUATE=""; REEVAL_CHECKS=""; REEVAL_MCP_AUTH=""; IDS_FILE=""; IDS=""; GATE=""
WARM_GIVEN=0; RESUME_SELECTION=0; PAUSE_GIVEN=0; GIVEN=""
# the run's flags this invocation names ('' none; else the option, a negation included): they apply to it only
LOCAL_SET=""; REMOTE_SET=""; STUB_SET=""; GATE_SET=""
# set_once VAR OPTION: an option and its negation exclude each other
set_once() {
  local cur="${!1}"
  [[ -z "$cur" || "$cur" == "$2" ]] || die "--$cur and --$2 exclude each other"
  printf -v "$1" '%s' "$2"
}
while [[ $# -gt 0 ]]; do
  # the options this invocation gives (recorded; a resume restores the run's own for the others)
  case "$1" in --*) GIVEN+="${1#--},";; esac
  case "$1" in
    --reevaluate) REEVALUATE="$2"; shift 2 ;;
    --reeval-checks) REEVAL_CHECKS="$2"; shift 2 ;;
    --reeval-mcp-auth) REEVAL_MCP_AUTH="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --smoke) SMOKE=1; shift ;;
    --paired) MODE=paired; shift ;;
    --sequential) MODE=sequential; shift ;;
    --areas) AREAS="$2"; shift 2 ;;
    --ids) IDS_FILE="$2"; shift 2 ;;
    --gate) set_once GATE_SET gate; GATE="$2"; shift 2 ;;
    --no-gate) set_once GATE_SET no-gate; GATE=""; shift ;;
    --limit) LIMIT="$2"; shift 2 ;;
    --skip-local) set_once LOCAL_SET skip-local; SKIP_LOCAL=1; shift ;;
    --local) set_once LOCAL_SET local; SKIP_LOCAL=0; shift ;;
    --skip-remote) set_once REMOTE_SET skip-remote; SKIP_REMOTE=1; shift ;;
    --remote) set_once REMOTE_SET remote; SKIP_REMOTE=0; shift ;;
    --mcp-run) MCP_RUN="$2"; shift 2 ;;
    --mcp-cmd) MCP_CMD="$2"; MCP_GIVEN=1; shift 2 ;;
    --no-mcp) NO_MCP=1; MCP_CMD=""; shift ;;
    --compare) COMPARE="$2"; shift 2 ;;
    --resume) RESUME="$2"; shift 2 ;;
    --by-design) BY_DESIGN="$2"; shift 2 ;;
    --pause) PAUSE="$2"; PAUSE_GIVEN=1; shift 2 ;;
    --max-requests) MAX_REQUESTS="$2"; shift 2 ;;
    --warm-every) WARM_EVERY="$2"; WARM_GIVEN=1; shift 2 ;;
    --allow-stub) set_once STUB_SET allow-stub; ALLOW_STUB=1; shift ;;
    --no-allow-stub) set_once STUB_SET no-allow-stub; ALLOW_STUB=0; shift ;;
    --checks) CHECKS_MODE="$2"; shift 2 ;;
    --no-pdf) NO_PDF=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown option $1" ;;
  esac
done
[[ -z "$MCP_RUN" || -d "$MCP_RUN" ]] || die "--mcp-run $MCP_RUN is not a directory"
[[ -z "$MCP_RUN" || -z "$MCP_CMD" ]] || die "--mcp-run and --mcp-cmd exclude each other"
[[ $NO_MCP == 0 || $MCP_GIVEN == 0 ]] || die "--no-mcp and --mcp-cmd exclude each other"
[[ "$MCP_CMD" != default ]] || MCP_CMD="$DEFAULT_MCP_CMD"
[[ -z "$COMPARE" || -d "$COMPARE" ]] || die "--compare $COMPARE is not a directory"
[[ -z "$CHECKS_MODE" || "$CHECKS_MODE" =~ ^(auto|checks|stub)$ ]] || die "--checks takes auto, checks or stub"
[[ -z "$IDS_FILE" || -f "$IDS_FILE" ]] || die "--ids $IDS_FILE is not a file"
check_gate() {
  [[ -f "$GATE" ]] || die "--gate $GATE is not a file"
  GATE="$(cd "$(dirname "$GATE")" && pwd)/$(basename "$GATE")"
  "$PY" "$HERE/gate.py" policy-ids "$GATE" > /dev/null || die "--gate $GATE is not a gate policy"
}
[[ -z "$GATE" ]] || check_gate

# --- the offline re-evaluation of a recorded run (lane S4): no lock, no ssh, no request to any host ----------
if [[ -n "$REEVALUATE" ]]; then
  [[ -d "$REEVALUATE/arms" ]] || die "--reevaluate $REEVALUATE has no arms/ directory"
  [[ -z "$REEVAL_CHECKS" || -f "$REEVAL_CHECKS/checks.py" ]] || die "--reeval-checks $REEVAL_CHECKS has no checks.py"
  REEVALUATE="$(cd "$REEVALUATE" && pwd)"
  RE_ARGS=("$REEVALUATE")
  [[ -z "$REEVAL_CHECKS" ]] || RE_ARGS+=(--checks "$REEVAL_CHECKS")
  [[ -z "$REEVAL_MCP_AUTH" ]] || RE_ARGS+=(--mcp-auth "$REEVAL_MCP_AUTH")
  [[ -z "$BY_DESIGN" ]] || RE_ARGS+=(--by-design "$BY_DESIGN")
  "$PY" "$HERE/reevaluate.py" "${RE_ARGS[@]}"
  MERGE_RC=0
  if [[ -n "$GATE" ]]; then
    "$PY" "$HERE/merge_run.py" "$REEVALUATE" --reevaluated --strict-identity || MERGE_RC=$?
    [[ $MERGE_RC == 0 || $MERGE_RC == 3 ]] || exit "$MERGE_RC"
  else
    "$PY" "$HERE/merge_run.py" "$REEVALUATE" --reevaluated
  fi
  CMP_ARGS=("$REEVALUATE" --no-merge --reevaluated)
  [[ -z "$BY_DESIGN" ]] || CMP_ARGS+=(--by-design "$BY_DESIGN")
  "$PY" "$HERE/compare.py" "${CMP_ARGS[@]}"
  GATE_RC=0
  if [[ -n "$GATE" ]]; then
    "$PY" "$HERE/gate.py" "$REEVALUATE" "$GATE" --reevaluated >&2 || GATE_RC=$?
    [[ $MERGE_RC == 0 || $GATE_RC != 0 ]] || GATE_RC=1
  fi
  REP_ARGS=("$REEVALUATE" --reevaluated)
  [[ -z "$COMPARE" ]] || REP_ARGS+=(--compare "$COMPARE")
  [[ $NO_PDF == 0 ]] || REP_ARGS+=(--no-pdf)
  "$PY" "$HERE/report.py" "${REP_ARGS[@]}" || log "report: the PDF was not produced (see above); report.reevaluated.md is there"
  log "re-evaluated: $REEVALUATE"
  printf '%s\n' "$REEVALUATE"
  if [[ -n "$GATE" ]]; then
    log "gate: $( [[ $GATE_RC == 0 ]] && echo PASS || echo "FAIL (exit $GATE_RC)" ) (gate.reevaluated.json)"
    exit "$GATE_RC"
  fi
  exit 0
fi

# --- the run directory ------------------------------------------------------------------------------------
if [[ -n "$RESUME" ]]; then
  [[ -d "$RESUME" ]] || die "--resume $RESUME is not a directory"
  RUN_DIR="$(cd "$RESUME" && pwd)"
  RUN_ID="$(basename "$RUN_DIR")"
  ARM_REL="runs/$RUN_ID/arms/direct-local"
  if [[ -z "$MANIFEST" ]]; then
    for cand in "$RUN_DIR/manifest.json" "$RUN_DIR/arms/direct-local/manifest.json" \
                "$RUN_DIR/arms/direct-remote/manifest.json"; do
      if [[ -f "$cand" ]]; then MANIFEST="$cand"; break; fi
    done
  fi
  # The run's own options (its last invocation's, orchestrator.json sessions[-1]) for those this one does not give,
  # so a plain resume continues the run as it was started: the same arms, checks, repeats, pace, load cap and gate.
  RESTORED=""
  session_value() {  # session_value INDEX KEY: what the run's invocation INDEX (0 the first, -1 the last) recorded
    "$PY" -c 'import json, sys
try:
    v = json.load(open(sys.argv[1]))["sessions"][int(sys.argv[2])].get(sys.argv[3])
    print("" if v is None else v)
except (OSError, ValueError, KeyError, IndexError, AttributeError):
    print("")' "$RUN_DIR/orchestrator.json" "$1" "$2" 2>/dev/null || true
  }
  last_session() { session_value -1 "$1"; }
  first_session() { session_value 0 "$1"; }
  restore() {  # restore NAME VALUE: log what was restored
    RESTORED+="$1 $2, "
  }
  if [[ -z "$MODE" ]]; then MODE="$(last_session mode)"; [[ -z "$MODE" ]] || restore mode "$MODE"; fi
  if [[ -z "$CHECKS_MODE" ]]; then
    # the run's own checks, so a resume judges its executions alike
    CHECKS_MODE="$(last_session checks)"; [[ -z "$CHECKS_MODE" ]] || restore checks "$CHECKS_MODE"
  fi
  # The run's own selection (the union of every invocation's plan, below) unless this one names its own.
  [[ -n "$IDS_FILE" || -n "$AREAS" || -n "$LIMIT" || -n "$GATE" ]] || RESUME_SELECTION=1
  if [[ $WARM_GIVEN == 0 ]]; then
    # its warm_every likewise (the repeats it planned)
    WARM_EVERY="$(last_session warm_every)"; WARM_EVERY="${WARM_EVERY:-10}"; restore warm_every "$WARM_EVERY"
  fi
  if [[ $MCP_GIVEN == 0 && $NO_MCP == 0 && -z "$MCP_RUN" ]]; then
    # The run's own MCP template (the last invocation's), so a resume pairs the same arms; a template of before
    # {checks_dir} that judged with this directory judges with the run's snapshot instead.
    MCP_CMD="$(last_session mcp_cmd)"; restore mcp "$( [[ -n "$MCP_CMD" ]] && echo "its template" || echo none)"
    OLD_CHECKS_DIR="--checks-dir $(printf '%q' "$HERE")"
    NEW_CHECKS_DIR='--checks-dir {checks_dir}'
    if [[ "$MCP_CMD" == *"$OLD_CHECKS_DIR"* ]]; then MCP_CMD="${MCP_CMD//"$OLD_CHECKS_DIR"/$NEW_CHECKS_DIR}"; fi
  fi
  if [[ $PAUSE_GIVEN == 0 ]]; then
    v="$(last_session pause_s)"; if [[ -n "$v" ]]; then PAUSE="$v"; restore pause "$v"; fi
  fi
  # The run's flags are its FIRST invocation's (a resume naming one, or its negation, changes that invocation only,
  # so one `--resume --skip-local` never drops direct-local from the later resumes).
  for opt in skip_local skip_remote allow_stub smoke; do
    if [[ "$(first_session $opt)" == 1 ]]; then
      case $opt in
        skip_local) [[ -n "$LOCAL_SET" ]] || { SKIP_LOCAL=1; restore skip-local ""; } ;;
        skip_remote) [[ -n "$REMOTE_SET" ]] || { SKIP_REMOTE=1; restore skip-remote ""; } ;;
        allow_stub) [[ -n "$STUB_SET" ]] || { ALLOW_STUB=1; restore allow-stub ""; } ;;
        smoke) [[ $SMOKE == 1 ]] || { SMOKE=1; restore smoke ""; } ;;
      esac
    fi
  done
  if [[ -z "$GATE_SET" ]]; then
    v="$(first_session gate)"
    if [[ -n "$v" ]]; then
      [[ -f "$v" ]] || die "the run was gated by $v, which is gone: pass --gate POLICY"
      GATE="$v"; check_gate; restore gate "$GATE"
    fi
  fi
  log "resuming $RUN_ID$( [[ -z "$RESTORED" ]] || printf '; the run'"'"'s own options: %s' "${RESTORED%, }")"
else
  if [[ -z "$MANIFEST" ]]; then
    if [[ $SMOKE == 1 ]]; then MANIFEST="$HERE/manifest/suite_smoke_v4.json"; else MANIFEST="$HERE/manifest/suite_v4.json"; fi
  fi
  [[ -n "$MANIFEST" && -f "$MANIFEST" ]] || die "no manifest (${MANIFEST:-none}); pass --manifest"
  # A new run is named after the lock is taken (the release GET below): review N4.
  RUN_ID=""; RUN_DIR=""
fi
MODE="${MODE:-paired}"
[[ "$MODE" == paired || "$MODE" == sequential ]] || die "unknown mode $MODE"
# The smoke load rule: at most ~40 requests per smoke invocation (a resumed smoke run's included).
if [[ $SMOKE == 1 && -z "$MAX_REQUESTS" ]]; then MAX_REQUESTS=40; fi
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
  [[ "$MCP_CMD" == *"{checks_dir}"* ]] || log "warning: the MCP template names no {checks_dir}: its runner judges with other checks than the run's snapshot (the merge states it, the gate fails on it)"
fi
[[ "$MODE" == sequential || -z "$MCP_RUN" ]] || die "--mcp-run reads a finished MCP run: use it with --sequential (or --mcp-cmd)"
[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || die "no manifest (${MANIFEST:-none}); pass --manifest"
MANIFEST="$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")"
if [[ -n "$IDS_FILE" ]]; then
  IDS="$("$PY" "$HERE/gate.py" ids "$IDS_FILE" --manifest "$MANIFEST")" || die "--ids $IDS_FILE: see above"
elif [[ $RESUME_SELECTION == 1 ]]; then
  # every scenario any invocation of the run planned (plan.json, else orchestrator.json's sessions): a plain resume
  # after a narrower targeted one finishes the run; '' (no plan) is the whole manifest, as before
  IDS="$("$PY" "$HERE/gate.py" resume-ids "$RUN_DIR" --manifest "$MANIFEST")" || die "the run's plan could not be read"
elif [[ -n "$GATE" && -z "$AREAS" && -z "$LIMIT" ]]; then
  # the policy's selection, for the arms this invocation runs
  GATE_ARMS=""
  [[ $SKIP_LOCAL == 1 ]] || GATE_ARMS+="direct-local,"
  [[ $SKIP_REMOTE == 1 ]] || GATE_ARMS+="direct-remote,"
  [[ -z "$MCP_CMD" && -z "$MCP_RUN" ]] || GATE_ARMS+="mcp,"
  IDS="$("$PY" "$HERE/gate.py" policy-ids "$GATE" --arms "${GATE_ARMS%,}")" || die "--gate $GATE: see above"
  [[ -n "$IDS" ]] || die "--gate $GATE selects no scenario for the arms ${GATE_ARMS%,}"
fi

# --- the checks (review M3) -----------------------------------------------------------------------------------
case "${CHECKS_MODE:-auto}" in
  auto) if [[ -f "$HERE/checks.py" ]]; then CHECKS=checks; else CHECKS=stub; fi ;;
  checks) [[ -f "$HERE/checks.py" ]] || die "--checks checks: no checks.py next to this script"; CHECKS=checks ;;
  stub) CHECKS=stub ;;
esac
if [[ $CHECKS == stub && $ALLOW_STUB == 0 ]]; then
  BEYOND="$("$PY" "$HERE/suite_common.py" beyond-status "$MANIFEST" "$AREAS" "$LIMIT" "$IDS")"
  if [[ -f "$HERE/checks.py" ]]; then WHY="--checks stub"; else WHY="no checks.py next to this script"; fi
  [[ -z "$BEYOND" ]] || die "the checks stub runs ($WHY), and the selected scenarios expect checks it cannot evaluate (${BEYOND:0:200}); run lane S1's checks.py, or pass --allow-stub for a run that checks the HTTP status only"
fi

# --- state for the traps --------------------------------------------------------------------------------------
LOCK_HELD=0; LOCK_TOKEN=""; LOCK_HOST=""; LOCK_STARTED_UTC=""; LOCK_STARTED_EPOCH=""
SAMPLER_PID=""; HB_REL=""; REMOTE_PID=""; LOOP_STARTED=0; CLEANED=0; ARMS_STARTED=0
CREATED_RUN_DIR=0; CREATED_MEX_RUN=0; ARM_REL="${ARM_REL:-}"

# Touches the paired sampler's heartbeat every SUITE_HEARTBEAT_S (300) seconds while this invocation lives (review
# N7): one step longer than the sampler's heartbeat max age (an MCP step waiting on the service queue, a long-class
# request) no longer ends the sampling, while a sampler whose orchestrator died still stops after that age.
HB_KEEPER=""
start_heartbeat_keeper() {
  [[ -n "$SAMPLER_PID" && -n "$HB_REL" ]] || return 0
  local every="${SUITE_HEARTBEAT_S:-300}" me=$$
  (
    trap - INT TERM EXIT
    n=0
    while kill -0 "$me" 2>/dev/null; do
      sleep 1
      n=$((n + 1))
      if (( n >= every )); then mex "touch $MEX_DIR/$HB_REL" || true; n=0; fi
    done
  ) > /dev/null 2>&1 < /dev/null &
  HB_KEEPER=$!
}
stop_heartbeat_keeper() {
  [[ -n "$HB_KEEPER" ]] || return 0
  kill "$HB_KEEPER" 2>/dev/null || true
  wait "$HB_KEEPER" 2>/dev/null || true
  HB_KEEPER=""
}
stop_sampler() {
  stop_heartbeat_keeper
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
  local out
  out="$(mex "if grep -q $LOCK_TOKEN $MEX_DIR/LOCK 2>/dev/null; then rm -f $MEX_DIR/LOCK; fi; echo released" || true)"
  if [[ "$out" != *released* ]]; then
    log "the lock ~/$MEX_DIR/LOCK on $MEX could not be released (no answer); the next invocation removes it as stale once this process is gone"
  fi
  LOCK_HELD=0
}
close_master() {  # this invocation's master only (its own ControlPath)
  [[ -z "$MEX_LOCAL" ]] || return 0
  ssh -o "ControlPath=$CTL_PATH" -O exit "$MEX" >/dev/null 2>&1 || true
}
cleanup() {
  [[ $CLEANED == 0 ]] || return 0
  CLEANED=1
  stop_sampler
  [[ -z "$RUN_ID" ]] || pull_samples
  if [[ $ARMS_STARTED == 0 ]]; then
    # Refused before any arm ran: leave no empty run directory behind that this invocation created.
    [[ $CREATED_MEX_RUN == 0 ]] || mex "rmdir $MEX_DIR/runs/$RUN_ID/arms $MEX_DIR/runs/$RUN_ID" || true
    if [[ $CREATED_RUN_DIR == 1 ]]; then
      rm -f "$RUN_DIR/manifest.json"
      for f in "$RUN_DIR/scripts/"*; do [[ ! -f "$f" ]] || rm -f "$f"; done   # this invocation's snapshot
      rmdir "$RUN_DIR/scripts" "$RUN_DIR/arms" "$RUN_DIR" 2>/dev/null || true
    fi
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
# The pids of this suite's run_direct.py processes on mex in RUNNER_PIDS ('' for none); returns 1 when mex gave no
# answer (the marker is missing), which says nothing about the runner.
runner_pids() {
  local out
  out="$(mex "printf 'runners:'; pgrep -u \$(id -u) -f '$RUNNER_PATTERN' | tr '\n' ' '; echo" || true)"
  [[ "$out" == *runners:* ]] || return 1
  RUNNER_PIDS="$(printf '%s' "${out##*runners:}" | tr -s ' \n' ' ' | sed 's/^ //; s/ $//')"
}
# Waits until no run_direct.py of this suite runs on mex: 0 then; 1 when one still runs after `$1` seconds; 2 when
# mex did not answer the last check after `$1` seconds. A failing ssh is retried every 5 s within that time and
# never reads as "no runner" (review N3: the next arm would start next to a runner cut off from its ssh).
wait_no_runner() {
  local waited=0 last="" silent=0
  RUNNER_PIDS=""
  while true; do
    if runner_pids; then
      [[ -n "$RUNNER_PIDS" ]] || return 0
      last=running
    else
      last=silent
      silent=$((silent + 1))
      (( silent % 6 != 1 )) || log "no answer from $MEX while checking for a running run_direct.py; retrying for up to $1 s"
    fi
    if (( waited >= $1 )); then
      if [[ $last == running ]]; then
        log "run_direct.py still running on $MEX (pid(s) $RUNNER_PIDS)"
        return 1
      fi
      log "no answer from $MEX for $waited s: whether a run_direct.py still runs there is unknown"
      return 2
    fi
    sleep 5
    waited=$((waited + 5))
  done
}
# This run's runners on mex in OWN_PIDS and all of the suite's in RUNNER_PIDS, from one answer ('' for none);
# returns 1 when mex gave no answer.
own_runner_pids() {
  local out
  out="$(mex "printf 'runners:'; pgrep -u \$(id -u) -f '$RUNNER_PATTERN' | tr '\n' ' '; printf 'own:'; pgrep -u \$(id -u) -f '$RUNNER_PATTERN .*--run-id $RUN_ID' | tr '\n' ' '; echo" || true)"
  [[ "$out" == *runners:*own:* ]] || return 1
  OWN_PIDS="$(printf '%s' "${out##*own:}" | tr -s ' \n' ' ' | sed 's/^ //; s/ $//')"
  out="${out%%own:*}"
  RUNNER_PIDS="$(printf '%s' "${out##*runners:}" | tr -s ' \n' ' ' | sed 's/^ //; s/ $//')"
}
# A resume next to its own run's runner (a paired step cut off by Ctrl-C ends after its request in flight, which
# may take up to the long class's client timeout): waits for it, with a log line every minute, instead of
# refusing. 0: no runner of this run, or another's runs too (wait_no_runner judges that); 1: still running after
# `$1` seconds; 2: mex did not answer.
wait_own_runner() {
  local waited=0 next=0
  while true; do
    own_runner_pids || return 2
    [[ -n "$OWN_PIDS" && "$OWN_PIDS" == "$RUNNER_PIDS" ]] || return 0
    (( waited < $1 )) || return 1
    if (( waited >= next )); then
      log "this run's runner still runs on $MEX (pid(s) $OWN_PIDS): it ends after its request in flight; waiting for it ($waited of at most $1 s)"
      next=$((next + 60))
    fi
    sleep 5
    waited=$((waited + 5))
  done
}
interrupted() {
  trap - INT TERM
  log "interrupted"
  if [[ $LOCK_HELD == 1 && -n "$RUN_ID" ]]; then
    # The runner of the step in flight (paired) or of the arm (sequential) stops after its request.
    mex "pkill -INT -u \$(id -u) -f '$RUNNER_PATTERN .*--run-id $RUN_ID'" || true
    local wrc=0
    wait_no_runner 120 || wrc=$?
    case $wrc in
      1) if [[ $MODE == sequential ]]; then
           log "the mex runner still runs (it stops after its request in flight, which may take up to 30 min); what it recorded so far is pulled, the rest stays on mex; --resume attaches to it"
         else
           log "the mex runner still runs (it stops after its request in flight, which may take up to 30 min); what it recorded so far is pulled, the rest stays on mex; --resume waits for it (up to $WAIT_OWN_RUNNER_S s), then continues the run"
         fi ;;
      2) log "no answer from $MEX: the mex runner may still run there; --resume looks for it once mex answers and waits for it" ;;
    esac
    pull_local_arm
  fi
  [[ $LOOP_STARTED == 0 ]] || loop_record interrupted
  cleanup
  [[ -z "$RUN_DIR" || ! -d "$RUN_DIR" ]] || log "resume with: $0 --resume $RUN_DIR"
  exit 130
}
trap interrupted INT TERM
trap cleanup EXIT

# --- mex: reachable, the lock, no other runner (before any request to the host: review N4) ----------------------
mex "mkdir -p $MEX_DIR" || die "ssh $MEX failed"

# The lock document of this invocation (its run id once the run is named).
lock_json() {
  printf '{"pid": %d, "host": "%s", "hostname": "%s", "started_utc": "%s", "started_epoch": %s, "run_id": "%s", "mode": "%s", "token": "%s"}' \
    "$$" "$LOCK_HOST" "$(hostname)" "$LOCK_STARTED_UTC" "$LOCK_STARTED_EPOCH" "$RUN_ID" "$MODE" "$LOCK_TOKEN"
}
take_lock() {
  local json held verdict attempt brc
  # A stable id of this machine, not its network name (a DHCP name changes on a VPN reconnect): review N6.
  LOCK_HOST="$("$PY" "$HERE/suite_common.py" host-id)"
  LOCK_TOKEN="$("$PY" -c 'import uuid; print(uuid.uuid4().hex)')"
  LOCK_STARTED_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  LOCK_STARTED_EPOCH="$("$PY" -c 'import time; print("%.3f" % time.time())')"
  json="$(lock_json)"
  for attempt in 1 2 3; do
    if mex "(set -o noclobber; printf '%s\n' $(quote "$json") > $MEX_DIR/LOCK)"; then
      LOCK_HELD=1
      log "lock ~/$MEX_DIR/LOCK taken on $MEX"
      return 0
    fi
    held="$(mex "cat $MEX_DIR/LOCK" || true)"
    [[ -n "$held" ]] || continue
    verdict="$("$PY" "$HERE/suite_common.py" lock-check "$held" "$LOCK_HOST")"
    case "$verdict" in
      stale*)
        # Removed only if it still holds what was judged stale (review N6).
        log "removing a stale lock (${verdict#stale })"
        brc=0
        mex "$("$PY" "$HERE/suite_common.py" lock-break-cmd "$MEX_DIR" "$held")" || brc=$?
        [[ $brc == 0 ]] || { log "the stale lock was not removed (it changed, or another removal is in progress); retrying"; sleep 1; } ;;
      *) die "the suite lock ~/$MEX_DIR/LOCK on $MEX is ${verdict#* }; another invocation drives the host (remove the file only if that process is gone)" ;;
    esac
  done
  die "could not take the lock ~/$MEX_DIR/LOCK on $MEX (a directory ~/$MEX_DIR/LOCK.break there is a stale-lock removal in progress, or one left by a cut-off invocation: remove it when no invocation runs)"
}
# Puts the run id into this invocation's lock once the run is named (only while the lock is still its own).
lock_update() {
  local json
  json="$(lock_json)"
  mex "grep -q $LOCK_TOKEN $MEX_DIR/LOCK && printf '%s\n' $(quote "$json") > $MEX_DIR/LOCK.$LOCK_TOKEN && mv -f $MEX_DIR/LOCK.$LOCK_TOKEN $MEX_DIR/LOCK" \
    || log "the run id could not be written into the lock (it stays held)"
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
if [[ -z "$ATTACH_PID" ]]; then
  WRC=0
  if [[ -n "$RESUME" ]]; then
    wait_own_runner "$WAIT_OWN_RUNNER_S" || WRC=$?
    case $WRC in
      1) die "this run's runner (pid(s) $OWN_PIDS) still runs on $MEX after $WAIT_OWN_RUNNER_S s: stop it with kill -INT on $MEX (it records its request in flight), then --resume $RUN_DIR again" ;;
      2) die "no answer from $MEX while checking for this run's runner: refusing to start (an ssh failure never reads as no runner)" ;;
    esac
  fi
  wait_no_runner "$WAIT_RUNNER_S" || WRC=$?
  case $WRC in
    0) ;;
    1) die "a run_direct.py runs under ~/$MEX_DIR on $MEX (another invocation's, pid(s) $RUNNER_PIDS); wait for it, or stop it with kill -INT on $MEX (it records its request in flight) and --resume that run later" ;;
    *) die "no answer from $MEX while checking for a running run_direct.py: refusing to start (an ssh failure never reads as no runner)" ;;
  esac
fi

# --- a new run: named from the host's release, now that the lock is held (review N4) -------------------------
if [[ -z "$RESUME" ]]; then
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
  ARM_REL="runs/$RUN_ID/arms/direct-local"
  lock_update
fi

if [[ ! -d "$RUN_DIR" ]]; then mkdir -p "$RUN_DIR/arms"; CREATED_RUN_DIR=1; else mkdir -p "$RUN_DIR/arms"; fi
# The run's manifest from the start, so a --resume finds it even when no arm's copy reached this side (merge_run
# keeps an existing one).
[[ -f "$RUN_DIR/manifest.json" ]] || cp "$MANIFEST" "$RUN_DIR/manifest.json"
COMMIT="$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
DIRTY="$(git -C "$HERE" status --porcelain -- . 2>/dev/null | grep -v '^?? runs/' | head -1 || true)"
[[ -z "$DIRTY" ]] || COMMIT="$COMMIT+dirty"
log "run directory $RUN_DIR ($MODE)"
log "manifest $MANIFEST; scripts $COMMIT; checks $CHECKS"
# The run's snapshot (header): taken by its first invocation (or by the first resume of a run older than it), kept.
SNAP="$RUN_DIR/scripts"
SESSIONS_BEFORE="$("$PY" -c 'import json, sys
try:
    print(len(json.load(open(sys.argv[1]))["sessions"]))
except (OSError, ValueError, KeyError, TypeError):
    print(0)' "$RUN_DIR/orchestrator.json" 2>/dev/null || echo 0)"
SNAP_STATE="$("$PY" "$HERE/suite_common.py" snapshot "$RUN_DIR" "$HERE" "$COMMIT" "$SESSIONS_BEFORE")" \
  || die "the run's snapshot of its scripts could not be taken in $SNAP"
case "$SNAP_STATE" in
  created*)
    if [[ $SESSIONS_BEFORE == 0 ]]; then log "scripts: the run's snapshot in $SNAP"
    else log "scripts: the run's snapshot in $SNAP, taken by this resume (its $SESSIONS_BEFORE earlier invocation(s) judged with the checks of their time)"; fi ;;
  *)
    SNAP_CHANGED="$(printf '%s\n' "$SNAP_STATE" | sed -n 's/^changed //p' | tr '\n' ' ')"
    SNAP_LIVE="$(printf '%s\n' "$SNAP_STATE" | sed -n 's/^live-differs //p' | tr '\n' ' ')"
    [[ -z "$SNAP_CHANGED" ]] || log "warning: the run's snapshot changed after it was taken (${SNAP_CHANGED% }): the merge states which arms judged with what"
    [[ -z "$SNAP_LIVE" ]] || log "scripts: ${SNAP_LIVE% } of $HERE differ(s) from the run's snapshot: the run keeps judging with its snapshot (a new run, or --reevaluate, judges with these)" ;;
esac

ARGS=(--manifest "$MANIFEST" --pause "$PAUSE" --scripts-commit "$COMMIT" --run-id "$RUN_ID" --checks "$CHECKS"
      --warm-every "$WARM_EVERY")
[[ -z "$AREAS" ]] || ARGS+=(--areas "$AREAS")
[[ -z "$IDS" ]] || ARGS+=(--ids "$IDS")
[[ -z "$LIMIT" ]] || ARGS+=(--limit "$LIMIT")
[[ $ALLOW_STUB == 0 ]] || ARGS+=(--allow-stub)
LOCAL_ARGS=("${ARGS[@]}")
LOCAL_ARGS[1]="manifest/$RUN_ID.json"

# --- mex: the run's directories, the clock ---------------------------------------------------------------------
[[ -n "$RESUME" ]] || mex "test -e $MEX_DIR/runs/$RUN_ID" || CREATED_MEX_RUN=1
mex "mkdir -p $MEX_DIR/scripts $MEX_DIR/manifest $MEX_DIR/runs/$RUN_ID/arms" || die "ssh $MEX failed"
T_A="$("$PY" -c 'import time; print("%.3f" % time.time())')"
# mex's clock from `date` (no interpreter start-up inside the measured interval: review N8); python where date
# has no %N (BSD date, the local test mode).
T_MEX="$(mex "t=\$(date +%s.%N 2>/dev/null); case \"\$t\" in ''|*[!0-9.]*) $MEX_PY -c 'import time; print(\"%.3f\" % time.time())' ;; *) echo \"\$t\" ;; esac" || true)"
T_B="$("$PY" -c 'import time; print("%.3f" % time.time())')"
CLOCK_OFFSET="$("$PY" "$HERE/suite_common.py" clock-offset "$T_A" "$T_MEX" "$T_B")"
log "mex clock offset ${CLOCK_OFFSET:-unknown} s"

mex_rsync --delete --exclude '__pycache__/' --exclude '*.pyc' \
  "$SNAP/" "$(mex_path "$MEX_DIR/scripts/")" || die "rsync of the run's snapshot to $MEX failed"
mex_rsync "$MANIFEST" "$(mex_path "$MEX_DIR/manifest/$RUN_ID.json")" || die "rsync of the manifest to $MEX failed"

"$PY" - "$RUN_DIR/orchestrator.json" "$MANIFEST" "$COMMIT" "$SMOKE" "$AREAS" "$LIMIT" "$SKIP_LOCAL" \
    "$SKIP_REMOTE" "$MCP_RUN" "$MCP_CMD" "$COMPARE" "$PAUSE" "$MAX_REQUESTS" "$LOCAL_BASE" "$REMOTE_BASE" "$MEX" \
    "$MODE" "$CHECKS" "$ALLOW_STUB" "$WARM_EVERY" "$CLOCK_OFFSET" "$(hostname)" "$IDS" "$IDS_FILE" "$GATE" \
    "$GIVEN" "${RESTORED:-}" "$SNAP_STATE" <<'PYEOF'
import json, os, sys, time
path = sys.argv[1]
keys = ["manifest", "scripts_commit", "smoke", "areas", "limit", "skip_local", "skip_remote", "mcp_run",
        "mcp_cmd", "compare", "pause_s", "max_requests", "local_base", "remote_base", "mex", "mode", "checks",
        "allow_stub", "warm_every", "mex_clock_offset_s", "host", "ids", "ids_file", "gate", "given", "restored",
        "snapshot"]
session = dict(zip(keys, sys.argv[2:]))
# the run's snapshot: taken by this invocation or kept, and what differed (suite_common.py snapshot)
lines = session["snapshot"].splitlines()
session["snapshot"] = {"state": lines[0] if lines else None,
                       "changed": [x.split(" ", 1)[1] for x in lines if x.startswith("changed ")],
                       "live_differs": [x.split(" ", 1)[1] for x in lines if x.startswith("live-differs ")],
                       "mcp_checks_dir": ("snapshot" if "{checks_dir}" in session["mcp_cmd"] else "template")
                       if session["mcp_cmd"] else None}
# the options given on the command line, and those a resume took from the run's last invocation
session["given"] = [x for x in session["given"].split(",") if x]
session["restored"] = [x.strip() for x in session["restored"].split(",") if x.strip()]
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
# The run's plan before any arm runs (gate.py plan): what this invocation's selection must record.
PLAN_ARMS=""
[[ $SKIP_LOCAL == 1 ]] || PLAN_ARMS+="direct-local,"
[[ $SKIP_REMOTE == 1 ]] || PLAN_ARMS+="direct-remote,"
[[ -z "$MCP_CMD" && -z "$MCP_RUN" ]] || PLAN_ARMS+="mcp,"
SESSION_N="$("$PY" -c 'import json, sys; print(len(json.load(open(sys.argv[1]))["sessions"]) - 1)' "$RUN_DIR/orchestrator.json")"
"$PY" "$HERE/gate.py" plan "$RUN_DIR" --manifest "$MANIFEST" --arms "${PLAN_ARMS%,}" --areas "$AREAS" --limit "$LIMIT" \
  --ids "$IDS" --warm-every "$WARM_EVERY" --mode "$MODE" --session "$SESSION_N" >&2 \
  || die "the run's plan could not be written"
if [[ -n "$RESUME" && $RESUME_SELECTION == 0 ]]; then
  # A resume that names its own selection: say when it leaves out scenarios of the run's plan (the run then stays
  # INCOMPLETE until a resume runs them; the report checks the whole plan).
  NARROW="$("$PY" - "$HERE" "$RUN_DIR" "$MANIFEST" "$AREAS" "$LIMIT" "$IDS" <<'PYEOF' || true
import json, sys
sys.path.insert(0, sys.argv[1])
import gate, suite_common as sc
run_dir, manifest = sys.argv[2], sc.read_json(sys.argv[3])
areas = [a for a in sys.argv[4].split(",") if a] or None
limit = int(sys.argv[5]) if sys.argv[5] else None
ids = [i for i in sys.argv[6].split(",") if i] or None
selected = set(s.get("id") for _, s in sc.select_scenarios(manifest, areas, limit, ids))
planned = gate.resume_ids(run_dir, manifest)
out = [i for i in planned if i not in selected]
if out:
    print(json.dumps({"selected": len(selected), "planned": len(planned), "left_out": len(out)}))
PYEOF
)"
  if [[ -n "$NARROW" ]]; then
    orch_set narrower_than_plan "$NARROW"
    log "this invocation's selection leaves out scenarios of the run's plan ($NARROW): they are not run now, and the run reads INCOMPLETE while any of them is not recorded; a plain --resume $RUN_DIR runs the run's whole plan"
  fi
fi

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
  cmd="$("$PY" - "$MCP_CMD" "$1" "$2" "$RUN_DIR/arms/mcp" "$MANIFEST" "$SNAP" <<'PYEOF'
import shlex, sys
template, ids, repeat, out, manifest, checks_dir = sys.argv[1:7]
values = {"ids": ids, "repeat": repeat, "repeats": str(int(repeat) + 1), "out": out, "manifest": manifest,
          "checks_dir": checks_dir}
for key, value in values.items():
    template = template.replace("{%s}" % key, shlex.quote(value))
print(template)
PYEOF
)"
  mkdir -p "$RUN_DIR/arms/mcp"
  log "mcp: $( [[ $2 == 0 ]] && echo cold || echo warm ) repeat $2 of ${1:0:120}"
  # The runner's output on the terminal and in arms/mcp/runner.log, so a runner that crashes leaves its traceback
  # in the run directory (the expanded command is never written: the template is recorded instead). tee ignores
  # SIGINT (-i): on Ctrl-C the runner stops after its execution in flight and writes its last lines, which a tee
  # killed by the same Ctrl-C would turn into a BrokenPipeError and exit 1 instead of 130.
  printf '[run_suite %s] mcp: repeat %s of %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$2" "$1" >> "$RUN_DIR/arms/mcp/runner.log"
  bash -c "$cmd" < /dev/null 2>&1 | tee -i -a "$RUN_DIR/arms/mcp/runner.log" >&2 || rc=$?
  printf '[run_suite %s] mcp: exit %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" >> "$RUN_DIR/arms/mcp/runner.log"
  sleep "$PAUSE"
  return $rc
}

MCP_ERRORS=0    # MCP executions in a row that ended in error (or that the runner never recorded)
MCP_MISSING=()  # the MCP executions a runner step did not record, as SCENARIO:rREPEAT
# After an MCP step of the paired loop (async's rules): stop on the runner's exit 2 / 3 (a start refusal: key
# file, the queue wait ran out) or 130 (interrupted), or on three executions in a row that ended `error` (a
# runner-level failure: transport, poll deadline, busy past the retries), were not recorded, or recorded an outcome
# outside the shared vocabulary (suite_common.outcome_known: `unknown`, never a success); never on `fail`
# (exit 1 also means a failing scenario, which is a result). An execution that ended `error` with its job
# cancelled may still run on the host briefly: wait until jobs.jsonl shows it terminal or its time budget has
# passed since the cancellation, so one request is in flight at a time. Returns 1 to stop the loop.
mcp_after() {  # mcp_after IDS REPEAT RC
  local rc="$3" res outcomes wait_s o said=0 i=0 ids
  case "$rc" in
    2|3) LOOP_REASON="mcp_exit_$rc"; log "error: the MCP runner refused to run (exit $rc): stopping the loop"; return 1 ;;
    130) LOOP_REASON=mcp_interrupted; log "the MCP runner was interrupted (exit 130): stopping the loop"; return 1 ;;
  esac
  while true; do
    res="$("$PY" "$HERE/suite_common.py" mcp-step "$RUN_DIR/arms/mcp" "$MANIFEST" "$1" "$2" 2>/dev/null || echo "missing 0")"
    outcomes="${res% *}"; wait_s="${res##* }"
    [[ "$wait_s" =~ ^[0-9]+$ && "$wait_s" != 0 ]] || break
    if [[ $said == 0 ]]; then
      log "mcp: a cancelled job of ${1:0:80} may still run on the host: waiting up to $wait_s s before the next step"
      said=1
    fi
    sleep $(( wait_s < 10 ? wait_s : 10 ))
  done
  ids=(${1//,/ })
  for o in ${outcomes//,/ }; do
    if [[ $o == error || $o == missing || $o == unknown ]]; then MCP_ERRORS=$((MCP_ERRORS + 1)); else MCP_ERRORS=0; fi
    [[ $o != missing ]] || MCP_MISSING+=("${ids[$i]:-?}:r$2")
    [[ $o != unknown ]] || log "mcp: ${ids[$i]:-?} r$2 recorded an outcome outside the shared vocabulary (merge_run reads it as error)"
    i=$((i + 1))
  done
  if (( MCP_ERRORS >= 3 )); then
    LOOP_REASON=mcp_errors
    log "error: $MCP_ERRORS MCP executions in a row ended in error, were not recorded or recorded an unknown outcome: stopping the loop"
    return 1
  fi
  return 0
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
    # its stderr into arms/direct-remote/runner.log, as direct-local's on mex (shown when the step fails)
    local errf="$RUN_DIR/.remote-step.err"
    mkdir -p "$RUN_DIR/arms/direct-remote"
    out="$("$PY" "$SNAP/run_direct.py" --arm direct-remote --base "$REMOTE_BASE" --out "$RUN_DIR/arms/direct-remote" \
      "${ARGS[@]}" "$@" --status-line ${cap_args[@]+"${cap_args[@]}"} < /dev/null 2> "$errf")" || rc=$?
    cat "$errf" >> "$RUN_DIR/arms/direct-remote/runner.log" 2>/dev/null || true
    if [[ $rc != 0 ]]; then tail -5 "$errf" | sed 's/^/    runner: /' >&2; fi
    rm -f "$errf"
  fi
  printf '%s\n' "$out" | grep -v '^RUN_DIRECT_STATUS ' | sed '/^$/d; s/^/    /' >&2 || true
  status="$(printf '%s\n' "$out" | sed -n 's/^RUN_DIRECT_STATUS //p' | tail -1)"
  requests="$(printf '%s' "$status" | sed -n 's/.*"requests": \([0-9][0-9]*\).*/\1/p')"
  USED=$((USED + ${requests:-0}))
  if [[ -z "$status" && $rc == 0 ]]; then rc=1; fi
  return $rc
}

FINAL_RC=0
FAILED_STEPS=()   # the steps that failed: ARM:SCENARIO:EXIT (paired), ARM:-:EXIT or STATE (sequential)
json_list() { local out="" x; for x in "$@"; do out+="\"$x\","; done; printf '[%s]' "${out%,}"; }
# The loop record of this invocation (orchestrator.json sessions[-1].loop: the report reads its state and reason).
# Paired: a scenario is completed when every step of it succeeded. Both modes: failed_steps names each step that
# did not, mcp_missing (the first 50, mcp_missing_count all) the MCP executions a runner step did not record.
loop_record() {  # loop_record STATE
  local counts="" missing=()
  [[ $MODE != paired ]] || counts="\"scenarios\": ${N_PLAN:-0}, \"completed\": ${COMPLETED:-0}, "
  missing=(${MCP_MISSING[@]+"${MCP_MISSING[@]:0:50}"})
  orch_set loop "{\"mode\": \"$MODE\", \"state\": \"$1\", \"reason\": \"${LOOP_REASON:-}\", $counts\"direct_requests\": ${USED:-0}, \"failed_steps\": $(json_list ${FAILED_STEPS[@]+"${FAILED_STEPS[@]}"}), \"mcp_missing\": $(json_list ${missing[@]+"${missing[@]}"}), \"mcp_missing_count\": ${#MCP_MISSING[@]}}"
}
# A loop that ran through without recording what it planned (a failed step, an MCP execution the runner did not
# record) ends incomplete, never done.
settle_loop_state() {
  [[ $LOOP_STATE != running ]] || LOOP_STATE=done
  if [[ $LOOP_STATE == done && ${#FAILED_STEPS[@]} -gt 0 ]]; then LOOP_STATE=incomplete; LOOP_REASON=failed_steps
  elif [[ $LOOP_STATE == done && ${#MCP_MISSING[@]} -gt 0 ]]; then LOOP_STATE=incomplete; LOOP_REASON=mcp_missing; fi
  [[ $LOOP_STATE != failed && $LOOP_STATE != incomplete ]] || FINAL_RC=1
}
loop_summary() {
  (( ${#FAILED_STEPS[@]} == 0 )) || printf '; %d failed step(s): %s' "${#FAILED_STEPS[@]}" "${FAILED_STEPS[*]}"
  (( ${#MCP_MISSING[@]} == 0 )) || printf '; %d MCP execution(s) not recorded by the runner: %s' \
    "${#MCP_MISSING[@]}" "${MCP_MISSING[*]:0:10}"
}
# --- paired mode ----------------------------------------------------------------------------------------------
run_paired() {
  local arms="" entry sid l r mc mw rep rc wrc fails_local=0 fails_remote=0 since_pull=0 runner_live=0 stop=0
  local step_failed nmissing
  [[ $SKIP_LOCAL == 1 ]] || { arms+="direct-local,"; RESERVE=$((RESERVE + 1)); }
  [[ $SKIP_REMOTE == 1 ]] || { arms+="direct-remote,"; RESERVE=$((RESERVE + 1)); }
  [[ -z "$MCP_CMD" ]] || arms+="mcp,"
  [[ -n "$arms" ]] || die "nothing to run (both direct arms skipped and no --mcp-cmd)"
  [[ $SKIP_LOCAL == 1 ]] || prepare_local_arm
  local plan=()
  while IFS= read -r entry; do plan+=("$entry"); done < <("$PY" "$HERE/suite_common.py" paired-plan "$MANIFEST" \
    "$RUN_DIR" "$AREAS" "$LIMIT" "$WARM_EVERY" "$arms" "$IDS")
  N_PLAN=${#plan[@]}
  COMPLETED=0
  local max_s=86400; [[ $SMOKE == 0 ]] || max_s=1800
  start_sampler paired "$max_s" heartbeat
  start_heartbeat_keeper
  LOOP_STARTED=1
  LOOP_STATE=running
  LOOP_REASON=""
  loop_record running
  log "paired: $N_PLAN scenario(s); arms ${arms%,}"
  for entry in ${plan[@]+"${plan[@]}"}; do
    read -r sid l r mc mw <<< "$entry"
    heartbeat
    step_failed=0
    nmissing=${#MCP_MISSING[@]}
    if [[ $mc == 1 ]]; then
      rc=0; mcp_step "$sid" 0 || rc=$?
      [[ $rc == 0 ]] || log "mcp: exit $rc on $sid"
      mcp_after "$sid" 0 "$rc" || { LOOP_STATE=failed; break; }
    fi
    if [[ $l == 1 ]]; then
      log "direct-local: $sid"
      rc=0; direct_call local --only "$sid" || rc=$?
      if [[ $rc == 3 ]]; then LOOP_STATE=capped; break; fi
      if [[ $rc != 0 ]]; then
        fails_local=$((fails_local + 1)); step_failed=1; FAILED_STEPS+=("direct-local:$sid:$rc")
        log "direct-local: exit $rc on $sid (see arms/direct-local/runner.log)"
        # An ssh that failed may have left the step running: never start the next arm next to it, and never read
        # a failing ssh as "the runner is gone" (review N3).
        wrc=0
        wait_no_runner "$WAIT_RUNNER_S" || wrc=$?
        if [[ $wrc == 1 ]]; then
          runner_live=1; LOOP_STATE=failed; LOOP_REASON=runner_still_running
          log "error: the direct-local step of $sid still runs on $MEX after $WAIT_RUNNER_S s: stopping the loop (--resume $RUN_DIR once it ended)"
          break
        elif [[ $wrc != 0 ]]; then
          runner_live=1; LOOP_STATE=failed; LOOP_REASON=mex_unreachable
          log "error: $MEX does not answer after the direct-local step of $sid failed, so its runner may still run there: stopping the loop, no further request (--resume $RUN_DIR once mex answers)"
          break
        fi
      else
        fails_local=0
      fi
    fi
    if [[ $r == 1 ]]; then
      log "direct-remote: $sid"
      rc=0; direct_call remote --only "$sid" || rc=$?
      if [[ $rc == 3 ]]; then LOOP_STATE=capped; break; fi
      if [[ $rc != 0 ]]; then
        fails_remote=$((fails_remote + 1)); step_failed=1; FAILED_STEPS+=("direct-remote:$sid:$rc")
        log "direct-remote: exit $rc on $sid"
      else
        fails_remote=0
      fi
    fi
    if [[ "$mw" != "-" ]]; then
      for rep in ${mw//,/ }; do
        rc=0; mcp_step "$sid" "$rep" || rc=$?
        [[ $rc == 0 ]] || log "mcp: exit $rc on $sid"
        mcp_after "$sid" "$rep" "$rc" || { stop=1; break; }
      done
      if [[ $stop == 1 ]]; then LOOP_STATE=failed; break; fi
    fi
    if (( fails_local >= 3 || fails_remote >= 3 )); then
      log "three failing steps of one arm in a row: stopping the loop"
      LOOP_STATE=failed; LOOP_REASON=failing_steps
      break
    fi
    (( ${#MCP_MISSING[@]} == nmissing )) || step_failed=1
    [[ $step_failed == 1 ]] || COMPLETED=$((COMPLETED + 1))
    since_pull=$((since_pull + 1))
    if [[ $SKIP_LOCAL == 0 && $since_pull -ge 25 ]]; then pull_local_arm; since_pull=0; fi
  done
  settle_loop_state
  [[ $LOOP_STATE != capped ]] || log "paired: stopped at --max-requests $MAX_REQUESTS ($USED direct requests); --resume continues"
  # The closing GET /capabilities of each direct arm (one request kept back for each) and its state over the
  # whole selection; none while a direct-local runner may still be on the host (review N3).
  if [[ $runner_live == 1 ]]; then
    log "paired: no closing steps (a direct-local runner may still run on $MEX); --resume finishes the run"
  else
    if [[ $SKIP_LOCAL == 0 ]]; then
      RESERVE=$((RESERVE - 1))
      direct_call local --finish || log "direct-local: the closing step failed"
    fi
    if [[ $SKIP_REMOTE == 0 ]]; then
      RESERVE=$((RESERVE - 1))
      direct_call remote --finish || log "direct-remote: the closing step failed"
    fi
  fi
  [[ $SKIP_LOCAL == 1 ]] || pull_local_arm
  stop_sampler
  loop_record "$LOOP_STATE"
  log "paired: $LOOP_STATE, $COMPLETED/$N_PLAN scenario(s) completed, $USED direct request(s)$(loop_summary)"
}

# --- sequential mode ------------------------------------------------------------------------------------------
run_sequential() {
  local rc state reason requests capped=0
  LOOP_STARTED=1
  LOOP_STATE=running
  LOOP_REASON=""
  loop_record running
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
    reason="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("reason") or "")' \
      "$RUN_DIR/arms/direct-local/state.json" 2>/dev/null || true)"
    requests="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("requests_session") or 0)' \
      "$RUN_DIR/arms/direct-local/state.json" 2>/dev/null || echo 0)"
    USED=$((USED + requests))
    log "direct-local: $state ($requests request(s))"
    if [[ "$state" == interrupted && "$reason" == max_requests* ]]; then
      capped=1
    elif [[ "$state" != done ]]; then
      FAILED_STEPS+=("direct-local:-:$state")
      log "direct-local did not finish ($state): see arms/direct-local/runner.log; --resume continues it"
    fi
  fi
  if [[ $SKIP_REMOTE == 0 ]]; then
    local max_s=43200; [[ $SMOKE == 0 ]] || max_s=1800
    start_sampler direct-remote "$max_s"
    log "direct-remote: from $(hostname) against $REMOTE_BASE"
    rc=0; direct_call remote || rc=$?
    stop_sampler
    case "$rc" in
      0) log "direct-remote: done" ;;
      3) capped=1; log "direct-remote: stopped at --max-requests $MAX_REQUESTS (incomplete; --resume continues it)" ;;
      *) FAILED_STEPS+=("direct-remote:-:$rc"); log "direct-remote: exit $rc (see above)" ;;
    esac
  fi
  if [[ -n "$MCP_CMD" ]]; then
    # The MCP arm last: repeat 0 of every scenario, then one pass per warm repeat index (1, then 2 for the
    # scenarios the paired plan gives three MCP repeats: suite_common.mcp_repeats), so a sequential run records the
    # same MCP executions as a paired one.
    local arms="mcp" entry sid l r mc mw rep cold="" maxrep=0 pass
    local warm=()
    while IFS= read -r entry; do
      read -r sid l r mc mw <<< "$entry"
      [[ $mc != 1 ]] || cold+="$sid,"
      if [[ "$mw" != "-" ]]; then
        for rep in ${mw//,/ }; do
          warm[$rep]="${warm[$rep]:-}$sid,"
          (( rep <= maxrep )) || maxrep=$rep
        done
      fi
    done < <("$PY" "$HERE/suite_common.py" paired-plan "$MANIFEST" "$RUN_DIR" "$AREAS" "$LIMIT" "$WARM_EVERY" "$arms" "$IDS")
    local go=1
    if [[ -n "$cold" ]]; then
      rc=0; mcp_step "${cold%,}" 0 || rc=$?
      [[ $rc == 0 ]] || log "mcp: the cold pass exited $rc"
      # the paired loop's rules between the passes: a cancelled job's wait, the stop on a refusal or errors
      mcp_after "${cold%,}" 0 "$rc" || { go=0; LOOP_STATE=failed; log "mcp: no warm pass"; }
    fi
    for (( rep = 1; rep <= maxrep; rep++ )); do
      pass="${warm[$rep]:-}"
      [[ -n "$pass" && $go == 1 ]] || continue
      rc=0; mcp_step "${pass%,}" "$rep" || rc=$?
      [[ $rc == 0 ]] || log "mcp: the warm pass of repeat $rep exited $rc"
      mcp_after "${pass%,}" "$rep" "$rc" || { go=0; LOOP_STATE=failed; log "mcp: no further warm pass"; }
    done
  fi
  if [[ -n "$MCP_RUN" ]]; then
    log "mcp: copying $MCP_RUN"
    mkdir -p "$RUN_DIR/arms/mcp"
    rsync -a "$MCP_RUN/" "$RUN_DIR/arms/mcp/"
  fi
  [[ $LOOP_STATE != running || $capped == 0 ]] || LOOP_STATE=capped
  settle_loop_state
  loop_record "$LOOP_STATE"
  log "sequential: $LOOP_STATE$( [[ -z "$LOOP_REASON" ]] || printf ' (%s)' "$LOOP_REASON"), $USED direct request(s)$(loop_summary)"
}

ARMS_STARTED=1
if [[ $MODE == paired ]]; then run_paired; else run_sequential; fi
cleanup

# --- merge, compare, gate, report -----------------------------------------------------------------------------
MERGE_RC=0
if [[ -n "$GATE" ]]; then
  "$PY" "$HERE/merge_run.py" "$RUN_DIR" --strict-identity || MERGE_RC=$?
  [[ $MERGE_RC == 0 || $MERGE_RC == 3 ]] || exit "$MERGE_RC"
else
  "$PY" "$HERE/merge_run.py" "$RUN_DIR"
fi
CMP_ARGS=("$RUN_DIR" --no-merge)
[[ -z "$BY_DESIGN" ]] || CMP_ARGS+=(--by-design "$BY_DESIGN")
"$PY" "$HERE/compare.py" "${CMP_ARGS[@]}"
GATE_RC=0
if [[ -n "$GATE" ]]; then
  "$PY" "$HERE/gate.py" "$RUN_DIR" "$GATE" >&2 || GATE_RC=$?
  [[ $MERGE_RC == 0 || $GATE_RC != 0 ]] || GATE_RC=1
fi
REP_ARGS=("$RUN_DIR")
[[ -z "$COMPARE" ]] || REP_ARGS+=(--compare "$COMPARE")
[[ $NO_PDF == 0 ]] || REP_ARGS+=(--no-pdf)
"$PY" "$HERE/report.py" "${REP_ARGS[@]}" || log "report: the PDF was not produced (see above); report.md is there"
log "done: $RUN_DIR"
printf '%s\n' "$RUN_DIR"
if [[ -n "$GATE" ]]; then
  log "gate: $( [[ $GATE_RC == 0 ]] && echo PASS || echo "FAIL (exit $GATE_RC)" ) (gate.json)"
  exit "$GATE_RC"
fi
exit "$FINAL_RC"
