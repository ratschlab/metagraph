#!/usr/bin/env bash
#
# make_wide_coord_fixture.sh -- the wide record-coordinates fixture of measurement M2
# (DESIGN-traverse-graphlet.md §18.8): the same
# few-kb sequence in many records under ONE column, so that every run of a trace walk carries
# thousands of coordinate chains (column labels), and the same records as header labels
# through a CoordToHeader (many labels with one chain each). mini_refseq has at most 6 chains
# a run with header labels and 14 with column labels, and the many-label fixture one chain
# per label: neither reaches the regime where the coordinates dominate a seed's account and
# text, which is what the delivery reserve's coordinate share (coordinate_account_per_text_byte)
# has to bound.
#
#   graph:      DBGSuccinct, k=31, BASIC mode (one path of LENGTH - 30 k-mers)
#   annotation: RowDiff<BRWT> with k-mer coordinates, one column (wide.fa), anchors sparse
#               (--max-path-length MAX_PATH: a row is reconstructed along long row-diff paths)
#               -> <out>/wide.row_diff_brwt_coord.annodbg
#   sidecar:    CoordToHeader (.seqs): record headers r0 .. r<N-1> as header labels
#
# usage: make_wide_coord_fixture.sh [OUT_DIR] [METAGRAPH_BIN]
#   OUT_DIR        default: <repo>/metagraph/build/wide_coord
#   METAGRAPH_BIN  default: <repo>/metagraph/build/metagraph
# env:
#   RECORDS   number of records (default 5000)
#   LENGTH    length of the shared sequence in bp (default 3000)
#   MAX_PATH  row-diff --max-path-length (default 1000: sparse anchors)
#   THREADS   (default 8)
#
# The sequence is pseudo-random from a fixed seed, so the fixture is the same on every build.
# NOTE: there is no `timeout` command on macOS; this script does not use one.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_MG_DIR="$(cd "$HERE/../.." && pwd)"
OUT="${1:-$REPO_MG_DIR/build/wide_coord}"
MG="${2:-$REPO_MG_DIR/build/metagraph}"
RECORDS="${RECORDS:-5000}"
LENGTH="${LENGTH:-3000}"
MAX_PATH="${MAX_PATH:-1000}"
THREADS="${THREADS:-8}"
K=31

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
[ -x "$MG" ] || { echo "metagraph binary not executable: $MG" >&2; exit 1; }
MG="$(cd "$(dirname "$MG")" && pwd)/$(basename "$MG")"
rm -rf "$OUT"/wide* "$OUT"/graph_k31.dbg* "$OUT"/rd_columns "$OUT"/build.log
mkdir -p "$OUT/rd_columns"
LOG="$OUT/build.log"

python3 - "$OUT/wide.fa" "$RECORDS" "$LENGTH" <<'PYEOF'
import random, sys
path, records, length = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(20261005)
seq = ''.join(rng.choice('ACGT') for _ in range(length))
with open(path, 'w') as f:
    for i in range(records):
        f.write(f'>r{i}\n{seq}\n')
PYEOF

{
    echo "== graph"
    "$MG" build -v -p "$THREADS" --mode basic --graph succinct --state stat -k "$K" \
        --mem-cap-gb 2 -o "$OUT/graph_k31" "$OUT/wide.fa"
    # Both annotate calls run inside OUT on the relative name wide.fa: --anno-filename makes
    # the file name the column's label, and an absolute path would make the label, and with it
    # the text and account of M2's column rows, depend on where the fixture was built (the
    # column row's text would be 14,916 B in one directory and 14,738 B in another)
    echo "== annotation with coordinates (one column, labelled wide.fa)"
    (cd "$OUT" && "$MG" annotate -v -p 1 -i graph_k31.dbg --anno-filename --coordinates \
        --mem-cap-gb 2 -o wide wide.fa)
    echo "== row_diff stages 0/1/2, --max-path-length $MAX_PATH"
    for stage in 0 1 2; do
        "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff --row-diff-stage "$stage" \
            --coordinates --max-path-length "$MAX_PATH" -i "$OUT/graph_k31.dbg" --mem-cap-gb 4 \
            -o "$OUT/rd_columns/out" "$OUT/wide.column.annodbg"
    done
    echo "== row_diff_brwt_coord"
    "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff_brwt_coord --greedy \
        -i "$OUT/graph_k31.dbg" -o "$OUT/wide" "$OUT/rd_columns/wide.column.annodbg"
    echo "== CoordToHeader"
    (cd "$OUT" && "$MG" annotate -v -p "$THREADS" -i graph_k31.dbg --anno-filename \
        --index-header-coords -o wide wide.fa)
} > "$LOG" 2>&1 || { tail -40 "$LOG" >&2; exit 1; }

for f in graph_k31.dbg graph_k31.dbg.anchors graph_k31.dbg.rd_succ wide.row_diff_brwt_coord.annodbg wide.seqs; do
    [ -f "$OUT/$f" ] || { echo "missing $OUT/$f (see $LOG)" >&2; exit 1; }
done
echo "$OUT: $RECORDS records of $LENGTH bp under one column (graph_k31.dbg, wide.row_diff_brwt_coord.annodbg, wide.seqs)"
