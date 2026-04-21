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
               "$UNIT_TESTS" --gtest_color=yes "$@" > "$LOG" 2>&1
    ) &
    PIDS[$i]=$!
done

# ANSI colors for the wrapper's own status lines (gtest handles its own).
GREEN=$'\033[0;32m'
RED=$'\033[0;31m'
RST=$'\033[0m'

FAILED=0
declare -a STATUS
declare -a STATUS_COLOR
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        STATUS[$i]="PASSED"
        STATUS_COLOR[$i]="${GREEN}PASSED${RST}"
    else
        STATUS[$i]="FAILED"
        STATUS_COLOR[$i]="${RED}FAILED${RST}"
        FAILED=$((FAILED + 1))
    fi
done

# Dump each shard's full output so test-level logs (timings, [ OK ] / [ FAILED ]
# lines, any test stdout) are visible. Sequential dump keeps lines from
# interleaving.
for i in "${!PIDS[@]}"; do
    echo
    echo "=============== Shard $i (${STATUS_COLOR[$i]}) ==============="
    cat "${LOGS[$i]}"
done

echo
echo "=== Per-shard summary ==="
for i in "${!PIDS[@]}"; do
    SUMMARY=$(grep -aE '\[==========\].*ran' "${LOGS[$i]}" | tail -1 \
              | sed -E $'s/\x1b\\[[0-9;]*m//g; s/^\\[==========\\] //')
    echo "Shard $i: ${STATUS_COLOR[$i]}  $SUMMARY"
done

# Aggregate gtest-style summary across all shards.
TOTAL_TESTS=0
TOTAL_PASSED=0
TOTAL_FAILED=0
declare -a FAILED_LIST
for i in "${!PIDS[@]}"; do
    PLAIN=$(sed -E $'s/\x1b\\[[0-9;]*m//g' "${LOGS[$i]}")

    # "[==========] N tests from M test suites ran. (X ms total)"
    RAN_LINE=$(echo "$PLAIN" | grep -E '^\[==========\].*ran\.' | tail -1)
    if [[ -n "$RAN_LINE" ]]; then
        N=$(echo "$RAN_LINE" | sed -E 's/^\[==========\] *([0-9]+).*/\1/')
        TOTAL_TESTS=$((TOTAL_TESTS + N))
    fi

    # "[  PASSED  ] N tests."
    P_LINE=$(echo "$PLAIN" | grep -E '^\[  PASSED  \] [0-9]+ tests?\.' | tail -1)
    if [[ -n "$P_LINE" ]]; then
        TOTAL_PASSED=$((TOTAL_PASSED + $(echo "$P_LINE" | sed -E 's/^\[  PASSED  \] *([0-9]+).*/\1/')))
    fi

    # "[  FAILED  ] N tests, listed below:" plus its "[  FAILED  ] Name..." lines
    F_HDR=$(echo "$PLAIN" | grep -E '^\[  FAILED  \] [0-9]+ tests?, listed below:')
    if [[ -n "$F_HDR" ]]; then
        N_FAIL=$(echo "$F_HDR" | sed -E 's/^\[  FAILED  \] *([0-9]+).*/\1/' | head -1)
        TOTAL_FAILED=$((TOTAL_FAILED + N_FAIL))
        # Lines between the header and the "N FAILED TEST(S)" trailer.
        while IFS= read -r line; do
            FAILED_LIST+=("$line")
        done < <(echo "$PLAIN" | awk '
            /^\[  FAILED  \] [0-9]+ tests?, listed below:/ { cap=1; next }
            /^[0-9]+ FAILED TESTS?/ { cap=0 }
            cap && /^\[  FAILED  \] / { print }')
    fi
done

# Match gtest's singular/plural.
plural() { [[ "$1" -eq 1 ]] && echo test || echo tests; }

ELAPSED=$(( $(date +%s) - START ))
echo
echo "${GREEN}[==========]${RST} $TOTAL_TESTS $(plural $TOTAL_TESTS) ran. (${ELAPSED} s wall, across $WORKERS shards)"
echo "${GREEN}[  PASSED  ]${RST} $TOTAL_PASSED $(plural $TOTAL_PASSED)."
if [[ $TOTAL_FAILED -gt 0 ]]; then
    echo "${RED}[  FAILED  ]${RST} $TOTAL_FAILED $(plural $TOTAL_FAILED), listed below:"
    for line in "${FAILED_LIST[@]}"; do
        # Re-colour the leading "[  FAILED  ]" on each listed name.
        echo "${RED}[  FAILED  ]${RST}${line#\[  FAILED  \]}"
    done
    echo
    [[ $TOTAL_FAILED -eq 1 ]] \
        && echo "${RED}1 FAILED TEST${RST}" \
        || echo "${RED}$TOTAL_FAILED FAILED TESTS${RST}"
fi

[[ $FAILED -gt 0 ]] && exit 1
exit 0
