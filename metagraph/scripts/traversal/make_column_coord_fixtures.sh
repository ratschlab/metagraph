#!/usr/bin/env bash
#
# make_column_coord_fixtures.sh -- the column-label fixtures of the depth gate (a column-kind
# fixture with at least 16 chains per run; DESIGN-traverse-graphlet.md §26.5). mini_refseq's
# taxid columns never reach that regime (at most 12 chains in a run on the recorded seeds),
# while refseq33m's taxid columns hold many genomes, so a run there carries a full list at the
# cap of 16 routinely. These fixtures give every column 16 records of its sequence, so that
# every run of a trace walk carries 16 chains:
#
#   lockstep  NCOLS_LOCKSTEP (200) columns c000.fa ..., each 16 records of ONE shared LENGTH-bp
#             sequence: the columns walk in lockstep, one path, every label a run of 16 chains
#   divcol    NCOLS_DIVCOL (100) columns, each 16 records of its own sequence: a shared CORE-bp
#             core between column-specific flanks of (LENGTH - CORE) / 2 bp, so a walk from
#             the core splits into one branch per column at each end of it
#
#   graph:      DBGSuccinct, k=31, BASIC mode
#   annotation: RowDiff<BRWT> with k-mer coordinates, one column per file (labels c000.fa, ...:
#               annotated inside fa/ on relative names, so the labels and every text and
#               account measured on them do not depend on where the fixture is built), row-diff
#               anchors sparse (--max-path-length MAX_PATH)
#               -> <out>/cols.row_diff_brwt_coord.annodbg
#
# usage: make_column_coord_fixtures.sh [lockstep|divcol|both] [BUILD_DIR] [METAGRAPH_BIN]
#   BUILD_DIR      default: <repo>/metagraph/build (writes BUILD_DIR/coord_lockstep and
#                  BUILD_DIR/coord_divcol, where MiniRefSeqWide.DISABLED_CoordinatesDepthByRegime
#                  looks for them when run from BUILD_DIR)
#   METAGRAPH_BIN  default: <repo>/metagraph/build/metagraph
# env:
#   NCOLS_LOCKSTEP (200), NCOLS_DIVCOL (100), RECORDS (16), LENGTH (3000), CORE (600),
#   MAX_PATH (1000), THREADS (8)
#
# The sequences are pseudo-random from fixed seeds, so the fixtures are the same on every build.
# NOTE: there is no `timeout` command on macOS; this script does not use one.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_MG_DIR="$(cd "$HERE/../.." && pwd)"
WHICH="${1:-both}"
BUILD_DIR="${2:-$REPO_MG_DIR/build}"
MG="${3:-$REPO_MG_DIR/build/metagraph}"
NCOLS_LOCKSTEP="${NCOLS_LOCKSTEP:-200}"
NCOLS_DIVCOL="${NCOLS_DIVCOL:-100}"
RECORDS="${RECORDS:-16}"
LENGTH="${LENGTH:-3000}"
CORE="${CORE:-600}"
MAX_PATH="${MAX_PATH:-1000}"
THREADS="${THREADS:-8}"
K=31

case "$WHICH" in
    lockstep|divcol|both) ;;
    *) echo "usage: $0 [lockstep|divcol|both] [BUILD_DIR] [METAGRAPH_BIN]" >&2; exit 2 ;;
esac
[ -x "$MG" ] || { echo "metagraph binary not executable: $MG" >&2; exit 1; }
MG="$(cd "$(dirname "$MG")" && pwd)/$(basename "$MG")"
mkdir -p "$BUILD_DIR"
BUILD_DIR="$(cd "$BUILD_DIR" && pwd)"

build_fixture() {
    local kind="$1" ncols="$2"
    local out="$BUILD_DIR/coord_$kind"
    rm -rf "$out"
    mkdir -p "$out/fa" "$out/columns" "$out/rd_columns"
    local log="$out/build.log"

    python3 - "$out" "$kind" "$ncols" "$RECORDS" "$LENGTH" "$CORE" <<'PYEOF'
import random, sys
out, kind, ncols, records, length, core = sys.argv[1], sys.argv[2], *map(int, sys.argv[3:])
rng = random.Random(20261006)
def rand(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))
if kind == 'lockstep':
    shared = rand(length)
    seqs = [shared] * ncols
else:
    flank = (length - core) // 2
    middle = rand(core)
    seqs = [rand(flank) + middle + rand(length - core - flank) for _ in range(ncols)]
with open(out + '/all.fa', 'w') as every:
    for c, seq in enumerate(seqs):
        with open('%s/fa/c%03d.fa' % (out, c), 'w') as f:
            for r in range(records):
                f.write('>c%03d_r%d\n%s\n' % (c, r, seq))
        every.write('>c%03d\n%s\n' % (c, seq))
PYEOF

    {
        echo "== graph"
        "$MG" build -v -p "$THREADS" --mode basic --graph succinct --state stat -k "$K" \
            --mem-cap-gb 2 -o "$out/graph_k31" "$out/all.fa"
        echo "== annotation with coordinates (a column per file, labelled by its relative name)"
        (cd "$out/fa" && "$MG" annotate -v -p "$THREADS" -i ../graph_k31.dbg --separately \
            --anno-filename --coordinates --mem-cap-gb 2 -o ../columns c*.fa)
        echo "== row_diff stages 0/1/2, --max-path-length $MAX_PATH"
        for stage in 0 1 2; do
            "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff --row-diff-stage "$stage" \
                --coordinates --max-path-length "$MAX_PATH" -i "$out/graph_k31.dbg" \
                --mem-cap-gb 4 -o "$out/rd_columns/out" "$out"/columns/*.column.annodbg
        done
        echo "== row_diff_brwt_coord"
        "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff_brwt_coord --greedy \
            -i "$out/graph_k31.dbg" -o "$out/cols" "$out"/rd_columns/*.column.annodbg
    } > "$log" 2>&1 || { tail -40 "$log" >&2; exit 1; }

    for f in graph_k31.dbg graph_k31.dbg.anchors cols.row_diff_brwt_coord.annodbg; do
        [ -f "$out/$f" ] || { echo "missing $out/$f (see $log)" >&2; exit 1; }
    done
    echo "$out: $ncols columns x $RECORDS records of $LENGTH bp ($kind)"
}

if [ "$WHICH" = lockstep ] || [ "$WHICH" = both ]; then
    build_fixture lockstep "$NCOLS_LOCKSTEP"
fi
if [ "$WHICH" = divcol ] || [ "$WHICH" = both ]; then
    build_fixture divcol "$NCOLS_DIVCOL"
fi
