#!/usr/bin/env bash
# Parallel runner for metagraph unit_tests using Google Test sharding.
#
# Each shard runs in its own working directory with a private copy of
# tests/data, so shards don't race on transient dump files.
#
# Usage:
#   unit_tests_parallel.sh [-j N] [unit_tests args...]
#
# Env vars:
#   UNIT_TESTS       path to the unit_tests binary (default: ./unit_tests)
#   TESTS_DATA       path to tests/data (default: ../tests/data relative to
#                    the unit_tests binary's dir)

set -u

WORKERS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
UNIT_TESTS="${UNIT_TESTS:-}"
TESTS_DATA="${TESTS_DATA:-}"

# Parse -j / --jobs
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            WORKERS="$2"
            shift 2
            ;;
        -j*)
            WORKERS="${1#-j}"
            shift
            ;;
        --jobs=*)
            WORKERS="${1#--jobs=}"
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            break
            ;;
    esac
done

if ! [[ "$WORKERS" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: invalid -j value: $WORKERS" >&2
    exit 2
fi

# Resolve UNIT_TESTS
if [[ -z "$UNIT_TESTS" ]]; then
    # Look next to this script, then in CWD
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

# Resolve TESTS_DATA
if [[ -z "$TESTS_DATA" ]]; then
    TESTS_DATA="$(dirname "$UNIT_TESTS")/../tests/data"
fi
if [[ ! -d "$TESTS_DATA" ]]; then
    echo "error: tests/data dir not found: $TESTS_DATA" >&2
    exit 2
fi
TESTS_DATA="$(cd "$TESTS_DATA" && pwd)"

# Per-shard workdir root (cleaned on exit)
ROOT="$(mktemp -d -t unit_tests_parallel.XXXXXX)"
trap 'rm -rf "$ROOT"' EXIT

echo "unit_tests: $UNIT_TESTS"
echo "tests/data: $TESTS_DATA"
echo "shards:     $WORKERS"
echo "workdir:    $ROOT"
echo

# Provision shard dirs in parallel (cp is the slow part).
echo "Preparing shard working directories..."
for ((i = 0; i < WORKERS; i++)); do
    (
        SHARD_BUILD="$ROOT/shard$i/build"
        SHARD_DATA="$ROOT/shard$i/tests/data"
        mkdir -p "$SHARD_BUILD" "$(dirname "$SHARD_DATA")"
        cp -r "$TESTS_DATA" "$SHARD_DATA"
    ) &
done
wait

# Launch shards.
echo "Running $WORKERS shards in parallel..."
START=$(date +%s)

declare -a PIDS
declare -a LOGS
for ((i = 0; i < WORKERS; i++)); do
    SHARD_BUILD="$ROOT/shard$i/build"
    LOG="$ROOT/shard$i.log"
    LOGS[$i]="$LOG"
    (
        cd "$SHARD_BUILD" \
            && GTEST_TOTAL_SHARDS="$WORKERS" GTEST_SHARD_INDEX="$i" \
               "$UNIT_TESTS" "$@" > "$LOG" 2>&1
    ) &
    PIDS[$i]=$!
done

FAILED=0
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        SUMMARY=$(grep -E '^\[==========\]' "${LOGS[$i]}" | tail -1)
        echo "Shard $i: PASSED  $SUMMARY"
    else
        FAILED=$((FAILED + 1))
        echo "Shard $i: FAILED  (see ${LOGS[$i]}):"
        tail -40 "${LOGS[$i]}" | sed 's/^/    /'
    fi
done

ELAPSED=$(( $(date +%s) - START ))
echo
echo "Elapsed: ${ELAPSED}s across $WORKERS shards"
if [[ $FAILED -gt 0 ]]; then
    echo "$FAILED shard(s) failed"
    exit 1
fi
echo "All shards passed"
exit 0
