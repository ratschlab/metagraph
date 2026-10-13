#!/usr/bin/env bash
# One command for a whole suite run (DESIGN §8.5): direct-local on mex with the resource sampler, then
# direct-remote from here through the public URL (a sampler on mex records that phase too), then the MCP arm's
# run directory if given, then merge, compare and report (md + pdf). The arms run one after the other, never
# two at once against the host.
#
#   run_suite.sh [--manifest M | --smoke] [--areas P,M] [--limit N] [--skip-local] [--skip-remote]
#                [--mcp-run DIR] [--compare EARLIER_RUN_DIR] [--resume RUN_DIR] [--by-design PATH]
#                [--pause S] [--max-requests N] [--no-pdf]
#
# A run directory is runs/<UTC>_<release-short>/ under this directory (SUITE_RUNS overrides). Every
# invocation starts a new run directory, except --resume RUN_DIR, which continues that run: each arm skips
# the executions it already recorded (the manifest is the run's own copy unless --manifest is given; another
# manifest is refused). Ctrl-C stops the arm in flight cleanly (the mex runner gets SIGINT), pulls what it
# recorded, and leaves the run resumable.
#
# On mex only ~/suite-dev/ is used (scripts, manifest, runs) and docker is only read (docker stats).
# Environment: SUITE_MEX (srv-metagraph@mex.inf.ethz.ch), SUITE_MEX_DIR (suite-dev), SUITE_MEX_PY (python3),
# SUITE_LOCAL_BASE (http://localhost:42010), SUITE_REMOTE_BASE (the public refseq33m-experimental URL),
# SUITE_RUNS, PYTHON (python3), SUITE_POLL_S (15), SUITE_REMOTE_SAMPLER (1; 0: no sampler on mex during the
# remote arm, so a --skip-local run never contacts mex).
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
REMOTE_SAMPLER="${SUITE_REMOTE_SAMPLER:-1}"
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=30 -o LogLevel=ERROR)

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }
log() { printf '[run_suite %s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }
die() { log "error: $*"; exit 2; }
mex() { ssh "${SSH_OPTS[@]}" "$MEX" "$@" 2>/dev/null; }
quote() { local out="" a; for a in "$@"; do out+="$(printf '%q' "$a") "; done; printf '%s' "$out"; }

MANIFEST=""; SMOKE=0; AREAS=""; LIMIT=""; SKIP_LOCAL=0; SKIP_REMOTE=0; MCP_RUN=""; COMPARE=""; RESUME=""
BY_DESIGN=""; NO_PDF=0; PAUSE="1"; MAX_REQUESTS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --smoke) SMOKE=1; shift ;;
    --areas) AREAS="$2"; shift 2 ;;
    --limit) LIMIT="$2"; shift 2 ;;
    --skip-local) SKIP_LOCAL=1; shift ;;
    --skip-remote) SKIP_REMOTE=1; shift ;;
    --mcp-run) MCP_RUN="$2"; shift 2 ;;
    --compare) COMPARE="$2"; shift 2 ;;
    --resume) RESUME="$2"; shift 2 ;;
    --by-design) BY_DESIGN="$2"; shift 2 ;;
    --pause) PAUSE="$2"; shift 2 ;;
    --max-requests) MAX_REQUESTS="$2"; shift 2 ;;
    --no-pdf) NO_PDF=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown option $1" ;;
  esac
done
# The smoke load rule: at most ~40 requests per arm and smoke run.
if [[ $SMOKE == 1 && -z "$MAX_REQUESTS" ]]; then MAX_REQUESTS=40; fi
[[ -z "$MCP_RUN" || -d "$MCP_RUN" ]] || die "--mcp-run $MCP_RUN is not a directory"
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
  log "resuming $RUN_ID"
else
  if [[ -z "$MANIFEST" ]]; then
    if [[ $SMOKE == 1 ]]; then MANIFEST="$HERE/manifest/suite_smoke_v1.json"; else MANIFEST="$HERE/manifest/suite_v1.json"; fi
  fi
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
[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || die "no manifest (${MANIFEST:-none}); pass --manifest"
MANIFEST="$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")"
mkdir -p "$RUN_DIR/arms"
COMMIT="$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
DIRTY="$(git -C "$HERE" status --porcelain -- . 2>/dev/null | grep -v '^?? runs/' | head -1 || true)"
[[ -z "$DIRTY" ]] || COMMIT="$COMMIT+dirty"
log "run directory $RUN_DIR"
log "manifest $MANIFEST; scripts $COMMIT"

"$PY" - "$RUN_DIR/orchestrator.json" "$MANIFEST" "$COMMIT" "$SMOKE" "$AREAS" "$LIMIT" "$SKIP_LOCAL" \
    "$SKIP_REMOTE" "$MCP_RUN" "$COMPARE" "$PAUSE" "$MAX_REQUESTS" "$LOCAL_BASE" "$REMOTE_BASE" "$MEX" <<'PYEOF'
import json, os, sys, time
path = sys.argv[1]
keys = ["manifest", "scripts_commit", "smoke", "areas", "limit", "skip_local", "skip_remote", "mcp_run",
        "compare", "pause_s", "max_requests", "local_base", "remote_base", "mex"]
session = dict(zip(keys, sys.argv[2:]))
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

ARGS=(--manifest "$MANIFEST" --pause "$PAUSE" --scripts-commit "$COMMIT" --run-id "$RUN_ID")
[[ -z "$AREAS" ]] || ARGS+=(--areas "$AREAS")
[[ -z "$LIMIT" ]] || ARGS+=(--limit "$LIMIT")
[[ -z "$MAX_REQUESTS" ]] || ARGS+=(--max-requests "$MAX_REQUESTS")

REMOTE_PID=""
SAMPLER_PID=""
cleanup() {
  if [[ -n "$SAMPLER_PID" ]]; then
    mex "kill -TERM $SAMPLER_PID" || true
    SAMPLER_PID=""
  fi
}
interrupted() {
  log "interrupted"
  if [[ -n "$REMOTE_PID" ]]; then
    log "stopping the mex runner (pid $REMOTE_PID) cleanly"
    mex "kill -INT $REMOTE_PID" || true
    for _ in $(seq 1 60); do
      mex "kill -0 $REMOTE_PID" || break
      sleep 2
    done
    rsync -a -e "ssh ${SSH_OPTS[*]}" "$MEX:$MEX_DIR/runs/$RUN_ID/arms/direct-local/" \
      "$RUN_DIR/arms/direct-local/" 2>/dev/null || true
  fi
  cleanup
  log "resume with: $0 --resume $RUN_DIR"
  exit 130
}
trap interrupted INT TERM
trap cleanup EXIT

sync_scripts() {
  mex "mkdir -p $MEX_DIR/scripts $MEX_DIR/manifest $MEX_DIR/runs/$RUN_ID/arms" || die "ssh $MEX failed"
  rsync -a -e "ssh ${SSH_OPTS[*]}" --include='*.py' --include='*.json' --exclude='*' "$HERE/" \
    "$MEX:$MEX_DIR/scripts/" || die "rsync of the scripts to $MEX failed"
  rsync -a -e "ssh ${SSH_OPTS[*]}" "$MANIFEST" "$MEX:$MEX_DIR/manifest/$RUN_ID.json" \
    || die "rsync of the manifest to $MEX failed"
}

# --- direct-local on mex ----------------------------------------------------------------------------------
if [[ $SKIP_LOCAL == 0 || ( $SKIP_REMOTE == 0 && $REMOTE_SAMPLER == 1 ) ]]; then
  sync_scripts
fi
if [[ $SKIP_LOCAL == 0 ]]; then
  ARM_REL="runs/$RUN_ID/arms/direct-local"
  if [[ -d "$RUN_DIR/arms/direct-local" ]] && ! mex "test -d $MEX_DIR/$ARM_REL"; then
    log "the mex copy of direct-local is gone: pushing the local one back before resuming"
    mex "mkdir -p $MEX_DIR/$ARM_REL"
    rsync -a -e "ssh ${SSH_OPTS[*]}" "$RUN_DIR/arms/direct-local/" "$MEX:$MEX_DIR/$ARM_REL/"
  fi
  LOCAL_ARGS=("${ARGS[@]}")
  LOCAL_ARGS[1]="manifest/$RUN_ID.json"
  CMD="$(quote "$MEX_PY" scripts/run_direct.py --arm direct-local --base "$LOCAL_BASE" --out "$ARM_REL" \
    --sample-docker "${LOCAL_ARGS[@]}")"
  # A runner of an earlier invocation still alive on mex (this script was cut off): attach to it.
  LIVE_PID="$(mex "cat $MEX_DIR/$ARM_REL/state.json 2>/dev/null" | "$PY" -c '
import json, sys
try:
    s = json.load(sys.stdin)
except ValueError:
    s = {}
print(s.get("pid", "") if s.get("state") == "running" else "")
' 2>/dev/null || true)"
  if [[ "$LIVE_PID" =~ ^[0-9]+$ ]] && mex "kill -0 $LIVE_PID"; then
    REMOTE_PID="$LIVE_PID"
    log "direct-local: attaching to the runner still running on $MEX (pid $REMOTE_PID)"
  else
    log "direct-local: starting on $MEX against $LOCAL_BASE"
    REMOTE_PID="$(mex "cd $MEX_DIR && mkdir -p $ARM_REL && (nohup $CMD >> $ARM_REL/runner.log 2>&1 < /dev/null & echo \$!)")"
    [[ "$REMOTE_PID" =~ ^[0-9]+$ ]] || die "could not start the runner on $MEX (got '$REMOTE_PID')"
  fi
  log "direct-local: mex pid $REMOTE_PID; polling every ${POLL_S}s"
  SILENT=0
  while true; do
    # An explicit alive / dead answer: a failed ssh says neither, and must never read as "finished" (the
    # remote arm would then run while direct-local still does).
    STATUS="$(mex "if kill -0 $REMOTE_PID 2>/dev/null; then echo alive; else echo dead; fi; cat $MEX_DIR/$ARM_REL/state.json 2>/dev/null" || true)"
    case "$STATUS" in
      alive*|dead*) SILENT=0 ;;
      *)
        SILENT=$((SILENT + 1))
        [[ $SILENT -lt 40 ]] || die "no answer from $MEX for $((SILENT * POLL_S)) s; the runner may still run there (pid $REMOTE_PID): --resume $RUN_DIR once it is reachable"
        log "direct-local: no answer from $MEX ($SILENT); retrying"
        sleep "$POLL_S"
        continue ;;
    esac
    PROGRESS="$(printf '%s' "$STATUS" | "$PY" -c '
import json, sys
text = sys.stdin.read()
i = text.find("{")
try:
    s = json.loads(text[i:]) if i >= 0 else {}
except ValueError:
    s = {}
print("%s %s/%s %s" % (s.get("state", "?"), s.get("done", "?"), s.get("planned", "?"), s.get("current") or ""))
' 2>/dev/null || true)"
    log "direct-local: $PROGRESS"
    case "$STATUS" in alive*) sleep "$POLL_S" ;; *) break ;; esac
  done
  REMOTE_PID=""
  rsync -a -e "ssh ${SSH_OPTS[*]}" "$MEX:$MEX_DIR/$ARM_REL/" "$RUN_DIR/arms/direct-local/" \
    || die "pulling direct-local from $MEX failed"
  STATE="$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("state"))' \
    "$RUN_DIR/arms/direct-local/state.json" 2>/dev/null || echo unknown)"
  log "direct-local: $STATE"
  [[ "$STATE" == "done" ]] || log "direct-local did not finish ($STATE): see arms/direct-local/runner.log; --resume continues it"
fi

# --- direct-remote from here --------------------------------------------------------------------------------
if [[ $SKIP_REMOTE == 0 ]]; then
  SAMPLES_REL="runs/$RUN_ID/remote-samples/samples-$("$PY" "$HERE/suite_common.py" stamp).jsonl"
  MAX_S=43200; [[ $SMOKE == 0 ]] || MAX_S=1800
  SCMD="$(quote "$MEX_PY" scripts/sampler.py --out "$SAMPLES_REL" --label direct-remote --max-seconds "$MAX_S")"
  SAMPLER_PID=""
  if [[ $REMOTE_SAMPLER == 1 ]]; then
    SAMPLER_PID="$(mex "cd $MEX_DIR && mkdir -p runs/$RUN_ID/remote-samples && (nohup $SCMD > /dev/null 2>&1 < /dev/null & echo \$!)" || true)"
  fi
  if [[ "$SAMPLER_PID" =~ ^[0-9]+$ ]]; then
    log "direct-remote: sampler on $MEX (pid $SAMPLER_PID)"
  else
    log "direct-remote: no sampler on $MEX (got '$SAMPLER_PID'); the remote arm has no resource samples"
    SAMPLER_PID=""
  fi
  log "direct-remote: from $(hostname) against $REMOTE_BASE"
  RC=0
  "$PY" "$HERE/run_direct.py" --arm direct-remote --base "$REMOTE_BASE" --out "$RUN_DIR/arms/direct-remote" \
    "${ARGS[@]}" || RC=$?
  if [[ -n "$SAMPLER_PID" ]]; then
    mex "kill -TERM $SAMPLER_PID" || true
    SAMPLER_PID=""
    sleep 3
    mkdir -p "$RUN_DIR/arms/direct-remote"
    if rsync -a -e "ssh ${SSH_OPTS[*]}" "$MEX:$MEX_DIR/$SAMPLES_REL" "$RUN_DIR/arms/direct-remote/.samples.part" \
        2>/dev/null; then
      cat "$RUN_DIR/arms/direct-remote/.samples.part" >> "$RUN_DIR/arms/direct-remote/samples.jsonl"
      rm -f "$RUN_DIR/arms/direct-remote/.samples.part"
    fi
  fi
  [[ $RC == 130 ]] && { log "direct-remote interrupted"; log "resume with: $0 --resume $RUN_DIR"; exit 130; }
  case "$RC" in
    0) log "direct-remote: done" ;;
    3) log "direct-remote: stopped at --max-requests $MAX_REQUESTS (incomplete; --resume continues it)" ;;
    *) log "direct-remote: exit $RC (see above)" ;;
  esac
fi

# --- the MCP arm (the async session's runner writes the same schema) -----------------------------------------
if [[ -n "$MCP_RUN" ]]; then
  log "mcp: copying $MCP_RUN"
  mkdir -p "$RUN_DIR/arms/mcp"
  rsync -a "$MCP_RUN/" "$RUN_DIR/arms/mcp/"
fi

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
