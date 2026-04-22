#!/usr/bin/env bash
# Parallel runner for metagraph unit_tests using Google Test sharding.
#
# Each shard runs in its own CWD ($ROOT/shard$i) so that temp files written
# to "." by third-party libs (sshash in particular dumps sshash.tmp.run_*
# into tmp_dirname, which defaults to ".") don't collide across shards.
# A shared $ROOT/tests symlink points at the real tests/ dir so the tests'
# hardcoded "../tests/data" still resolves. Tests' own scratch output goes
# to a per-process temp dir via test_dump_dir().
#
# Usage:
#   unit_tests_parallel.sh [-j N] [--help] [unit_tests args...]
#
# Env vars:
#   UNIT_TESTS       path to the unit_tests binary (default: ./unit_tests)

set -u -o pipefail

GREEN=$'\033[0;32m'
RED=$'\033[0;31m'
RST=$'\033[0m'

print_help() {
    sed -n '2,/^$/ { s/^# \{0,1\}//; p; }' "${BASH_SOURCE[0]}"
}

WORKERS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
UNIT_TESTS="${UNIT_TESTS:-}"

# Parse our flags; everything else is forwarded to unit_tests.
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)   WORKERS="$2"; shift 2 ;;
        -j*)         WORKERS="${1#-j}"; shift ;;
        --jobs=*)    WORKERS="${1#--jobs=}"; shift ;;
        -h|--help)   print_help; exit 0 ;;
        --)          shift; break ;;
        *)           break ;;
    esac
done

if ! [[ "$WORKERS" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: invalid -j value: $WORKERS" >&2
    exit 2
fi

# Resolve UNIT_TESTS
if [[ -z "$UNIT_TESTS" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -x "$SCRIPT_DIR/unit_tests" ]]; then
        UNIT_TESTS="$SCRIPT_DIR/unit_tests"
    elif [[ -x "./unit_tests" ]]; then
        UNIT_TESTS="$(pwd)/unit_tests"
    else
        echo "error: unit_tests binary not found; set UNIT_TESTS=<path>" >&2
        exit 2
    fi
fi

if [[ ! -x "$UNIT_TESTS" ]]; then
    echo "error: not executable: $UNIT_TESTS" >&2
    exit 2
fi
UNIT_TESTS="$(cd "$(dirname "$UNIT_TESTS")" && pwd)/$(basename "$UNIT_TESTS")"

# Locate the tests/ dir relative to the binary (metagraph/build/unit_tests
# -> metagraph/tests). Shards reference it via a shared symlink so hardcoded
# "../tests/data" lookups still resolve from each per-shard CWD.
TESTS_DIR="$(dirname "$UNIT_TESTS")/../tests"
if [[ ! -d "$TESTS_DIR" ]]; then
    echo "error: tests dir not found: $TESTS_DIR" >&2
    exit 2
fi
TESTS_DIR="$(cd "$TESTS_DIR" && pwd)"

# Per-shard workdirs and logs. On normal exit, clean up; on Ctrl-C/TERM,
# also kill any still-running shards so they don't outlive the wrapper.
ROOT="$(mktemp -d "${TMPDIR:-/tmp}/unit_tests_parallel.XXXXXX")" \
    || { echo "error: mktemp -d failed" >&2; exit 2; }
# Belt-and-suspenders guard for `rm -rf "$ROOT"`. Only ever delete paths that
# still match the pattern we created. Catches a hypothetical scenario where
# $ROOT somehow gets clobbered between creation and cleanup.
safe_cleanup() {
    case "$ROOT" in
        */unit_tests_parallel.??????)
            [[ -d "$ROOT" ]] && rm -rf "$ROOT"
            ;;
    esac
}
ln -s "$TESTS_DIR" "$ROOT/tests"
for ((i = 0; i < WORKERS; i++)); do
    mkdir -p "$ROOT/shard$i"
done
declare -a PIDS=()
on_exit() { safe_cleanup; }
on_interrupt() {
    trap - EXIT INT TERM
    # SIGTERM first, give them a moment, then SIGKILL anything still alive.
    # Must fully reap before rm -rf — otherwise a slow-dying child can write
    # new files into $ROOT after we walk it, leaving orphans behind.
    for pid in "${PIDS[@]}"; do
        kill -TERM "$pid" 2>/dev/null || true
    done
    for _ in 1 2 3 4 5 6 7 8 9 10; do
        alive=0
        for pid in "${PIDS[@]}"; do
            kill -0 "$pid" 2>/dev/null && alive=1 && break
        done
        [[ $alive -eq 0 ]] && break
        sleep 0.2
    done
    for pid in "${PIDS[@]}"; do
        kill -0 "$pid" 2>/dev/null && kill -KILL "$pid" 2>/dev/null || true
    done
    for pid in "${PIDS[@]}"; do
        wait "$pid" 2>/dev/null || true
    done
    safe_cleanup
    exit 130
}
trap on_exit EXIT
trap on_interrupt INT TERM

echo "unit_tests: $UNIT_TESTS"
echo "tests dir:  $TESTS_DIR"
echo "shards:     $WORKERS"
echo "workdir:    $ROOT"
echo

# Launch shards.
echo "Running $WORKERS shards in parallel..."
START=$(date +%s)

declare -a LOGS
for ((i = 0; i < WORKERS; i++)); do
    LOG="$ROOT/shard$i.log"
    LOGS[$i]="$LOG"
    # TEST_DUMP_DIR pins test_dump_dir() inside $ROOT so transient dump
    # files get swept along with the rest of the shard workdir, including
    # on Ctrl+C when C++ atexit handlers wouldn't have run.
    # `exec` replaces the subshell with unit_tests so $! points at the real
    # binary — otherwise `kill $!` just kills the subshell and orphans the
    # child, which then reparents to init and keeps running.
    (
        cd "$ROOT/shard$i" \
            && GTEST_TOTAL_SHARDS="$WORKERS" GTEST_SHARD_INDEX="$i" \
               TEST_DUMP_DIR="$ROOT/shard$i/dump" \
               exec "$UNIT_TESTS" --gtest_color=yes "$@"
    ) > "$LOG" 2>&1 &
    PIDS[$i]=$!
done

SHARDS_FAILED=0
declare -a STATUS_COLOR
declare -a DONE
for ((i = 0; i < WORKERS; i++)); do DONE[$i]=0; done

# Drain shards in completion order so each shard's full log lands in the
# output as soon as that shard finishes. Using `kill -0` + `sleep` polling
# instead of `wait -n -p` so we stay compatible with stock macOS bash 3.2.
# Relies on bash's internal job table: once SIGCHLD marks a background child
# dead, `kill -0 <pid>` returns non-zero even before we `wait` to reap it, so
# we're not at risk of busy-looping on a zombie.
remaining=$WORKERS
while (( remaining > 0 )); do
    progressed=0
    for ((i = 0; i < WORKERS; i++)); do
        [[ "${DONE[$i]}" == "1" ]] && continue
        kill -0 "${PIDS[$i]}" 2>/dev/null && continue   # still running
        rc=0
        wait "${PIDS[$i]}" || rc=$?
        DONE[$i]=1
        remaining=$(( remaining - 1 ))
        progressed=1
        secs=$(( $(date +%s) - START ))
        if [[ $rc -eq 0 ]]; then
            STATUS_COLOR[$i]="${GREEN}PASSED${RST}"
        else
            STATUS_COLOR[$i]="${RED}FAILED${RST}"
            SHARDS_FAILED=$((SHARDS_FAILED + 1))
        fi
        echo
        echo "=============== Shard ${i} (${STATUS_COLOR[$i]}, ${secs}s, $((WORKERS - remaining))/${WORKERS}) ==============="
        cat "${LOGS[$i]}"
    done
    (( progressed == 0 && remaining > 0 )) && sleep 1
done

echo
echo "=== Per-shard summary ==="
for i in "${!PIDS[@]}"; do
    SUMMARY=$(grep -aE '\[==========\].*ran' "${LOGS[$i]}" | tail -1 \
              | sed -E $'s/\x1b\\[[0-9;]*m//g; s/^\\[==========\\] //' || true)
    echo "Shard $i: ${STATUS_COLOR[$i]}  $SUMMARY"
done

# Aggregate gtest-style summary across all shards.
TOTAL_TESTS=0
TOTAL_PASSED=0
TOTAL_FAILED=0
FAILED_LIST=()
CRASHED_SHARDS=()
for i in "${!PIDS[@]}"; do
    PLAIN=$(sed -E $'s/\x1b\\[[0-9;]*m//g' "${LOGS[$i]}")

    # "[==========] N tests from M test suites ran. (X ms total)"
    RAN_LINE=$(echo "$PLAIN" | grep -E '^\[==========\].*ran\.' | tail -1 || true)
    if [[ -n "$RAN_LINE" ]]; then
        N=$(echo "$RAN_LINE" | sed -E 's/^\[==========\] *([0-9]+).*/\1/')
        TOTAL_TESTS=$((TOTAL_TESTS + N))
    elif [[ "${STATUS_COLOR[$i]}" == *FAILED* ]]; then
        # Shard exited non-zero without producing a gtest summary line —
        # a crash, an abort in global setup, or a filter that matched nothing.
        # Flag it explicitly so the aggregate doesn't look clean.
        CRASHED_SHARDS+=("$i")
    fi

    # "[  PASSED  ] N tests."
    P_LINE=$(echo "$PLAIN" | grep -E '^\[  PASSED  \] [0-9]+ tests?\.' | tail -1 || true)
    if [[ -n "$P_LINE" ]]; then
        P_N=$(echo "$P_LINE" | sed -E 's/^\[  PASSED  \] *([0-9]+).*/\1/')
        TOTAL_PASSED=$((TOTAL_PASSED + P_N))
    fi

    # "[  FAILED  ] N tests, listed below:" plus its "[  FAILED  ] Name..." lines
    F_HDR=$(echo "$PLAIN" | grep -E '^\[  FAILED  \] [0-9]+ tests?, listed below:' || true)
    if [[ -n "$F_HDR" ]]; then
        N_FAIL=$(echo "$F_HDR" | sed -E 's/^\[  FAILED  \] *([0-9]+).*/\1/' | head -1)
        TOTAL_FAILED=$((TOTAL_FAILED + N_FAIL))
        while IFS= read -r line; do
            FAILED_LIST+=("$line")
        done < <(echo "$PLAIN" | awk '
            /^\[  FAILED  \] [0-9]+ tests?, listed below:/ { cap=1; next }
            /^[0-9]+ FAILED TESTS?/ { cap=0 }
            cap && /^\[  FAILED  \] / { print }')
    fi
done

plural() { [[ "$1" -eq 1 ]] && echo test || echo tests; }

ELAPSED=$(( $(date +%s) - START ))
echo
echo "${GREEN}[==========]${RST} $TOTAL_TESTS $(plural $TOTAL_TESTS) ran. (${ELAPSED} s wall, across $WORKERS shards)"
echo "${GREEN}[  PASSED  ]${RST} $TOTAL_PASSED $(plural $TOTAL_PASSED)."
if [[ $TOTAL_FAILED -gt 0 ]]; then
    echo "${RED}[  FAILED  ]${RST} $TOTAL_FAILED $(plural $TOTAL_FAILED), listed below:"
    for line in "${FAILED_LIST[@]}"; do
        echo "${RED}[  FAILED  ]${RST}${line#\[  FAILED  \]}"
    done
    echo
    [[ $TOTAL_FAILED -eq 1 ]] \
        && echo "${RED}1 FAILED TEST${RST}" \
        || echo "${RED}$TOTAL_FAILED FAILED TESTS${RST}"
fi
if (( ${#CRASHED_SHARDS[@]} > 0 )); then
    echo
    echo "${RED}shard(s) exited abnormally without a test summary: ${CRASHED_SHARDS[*]}${RST}"
    echo "  (the test count above omits whatever these shards would have run)"
fi

[[ $SHARDS_FAILED -gt 0 ]] && exit 1
exit 0
