#!/usr/bin/env bash
#
# build_mini_refseq.sh -- build a small MetaGraph index that mirrors the public
# "refseq33m" index format, from real NCBI RefSeq records, for testing the
# labeled-traversal feature.
#
#   graph:      DBGSuccinct, k=31, BASIC mode, forward strand only, suffix ranges
#               indexed (--index-ranges 12), dummy k-mers NOT masked
#               -> <out>/graph_k31.dbg (+ .anchors / .rd_succ for RowDiff)
#   annotation: RowDiff<BRWT> with k-mer COORDINATES, one column per taxid FASTA,
#               BRWT relaxed, columns relabeled to the bare taxid
#               -> <out>/annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg
#   sidecar:    CoordToHeader mapping (file coords -> bare accession headers)
#               -> <out>/annotation.relaxed.relabeled.seqs
#
# usage: build_mini_refseq.sh [OUT_DIR] [METAGRAPH_BIN] [ACCESSIONS_TSV]
#   OUT_DIR        default: <repo>/metagraph/build/mini_refseq
#   METAGRAPH_BIN  default: <repo>/metagraph/build/metagraph
#   ACCESSIONS_TSV default: <this dir>/mini_refseq_accessions.tsv
#                  columns: accession taxid length group refseq33m_fw_mask refseq33m_rc_mask title
# env:
#   THREADS        parallelism for build/transform steps (default 8)
#   REUSE_RAW=1    reuse <out>/raw/<accession>.fa downloads if present (default: re-download)
#
#   MANIFEST_ONLY=1  only (re)write the index manifest of an existing <out> (no download,
#                  no build): for a fixture built before the manifest existed
#
# Steps: download (NCBI E-utilities, <=3 req/s, retries) -> headers rewritten to bare
# accession.version, grouped into <out>/fasta/<taxid>.fa -> build graph -> annotate with
# coordinates -> row_diff stages 0/1/2 (--coordinates) -> row_diff_brwt_coord -> relax_brwt
# -> rename columns to taxids -> CoordToHeader sidecar in final column order -> sanity
# queries with blaNDM-1 (fw + rc) -> <out>/manifest.json (the build record) -> the index
# manifest <out>/annotation.relaxed.relabeled.manifest.json (DESIGN-traverse-graphlet.md
# §3.1: every file of the served bundle with size and sha256; its digest is the index_fp
# that `--index-manifest` makes the server report).
#
# NOTE: there is no `timeout` command on macOS; this script does not use one.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_MG_DIR="$(cd "$HERE/../.." && pwd)"           # <repo>/metagraph

OUT="${1:-$REPO_MG_DIR/build/mini_refseq}"
MG="${2:-$REPO_MG_DIR/build/metagraph}"
TSV="${3:-$HERE/mini_refseq_accessions.tsv}"
THREADS="${THREADS:-8}"
REUSE_RAW="${REUSE_RAW:-0}"
MANIFEST_ONLY="${MANIFEST_ONLY:-0}"
K=31

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
MG="$(cd "$(dirname "$MG")" && pwd)/$(basename "$MG")"
TSV="$(cd "$(dirname "$TSV")" && pwd)/$(basename "$TSV")"

case "$OUT" in *[[:space:]]*) echo "OUT_DIR must not contain whitespace (rename rules are whitespace separated)" >&2; exit 1;; esac
[ -x "$MG" ] || { echo "metagraph binary not executable: $MG" >&2; exit 1; }
[ -f "$TSV" ] || { echo "accession table not found: $TSV" >&2; exit 1; }

ANNO_BASE="annotation"
FINAL_BASE="annotation.relaxed.relabeled"
FINAL_ANNO="$OUT/$FINAL_BASE.row_diff_brwt_coord.annodbg"
GRAPH="$OUT/graph_k31.dbg"
CMDLOG="$OUT/commands.log"

log()  { printf '\n[%s] == %s\n' "$(date +%H:%M:%S)" "$*"; }
run()  { printf '%q ' "$@" >> "$CMDLOG"; printf '\n' >> "$CMDLOG"; "$@"; }
note() { printf '# %s\n' "$*" >> "$CMDLOG"; }

# The index manifest (DESIGN-traverse-graphlet.md §3.1). The identity of an index is the
# content of the files a server loads -- graph, annotation and every sidecar -- not its
# metadata: two annotations over one graph can agree on k, columns and counts and still
# disagree on which column carries which sequence. Hashing is done here, at build time;
# the server reads the digests instead of hashing the bundle at every start. index_fp =
# sha256 over the lines "<path>\t<size>\t<sha256>\n" sorted by path (bytes), exactly as
# `metagraph traverse|server_query --index-manifest` recomputes and cross-checks it.
write_index_manifest() {
    log "index manifest (sha256 of every file of the bundle)"
    python3 - "$OUT" "$MG" "$K" "$TSV" "$FINAL_BASE" <<'PYEOF'
import hashlib, json, os, subprocess, sys
out, mg, k, tsv, final_base = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def entry(rel):
    path = os.path.join(out, rel)
    if not os.path.isfile(path):
        sys.exit(f'index manifest: missing bundle file {path}')
    return {'path': rel, 'size': os.path.getsize(path), 'sha256': digest(path)}

# what `-i graph_k31.dbg -a <annotation>` loads: the graph with its RowDiff sidecars, the
# annotation and the CoordToHeader sidecar
bundle = ['graph_k31.dbg', 'graph_k31.dbg.anchors', 'graph_k31.dbg.rd_succ',
          final_base + '.row_diff_brwt_coord.annodbg', final_base + '.seqs']
files = sorted((entry(rel) for rel in bundle), key=lambda e: e['path'].encode())
canonical = ''.join(f"{e['path']}\t{e['size']}\t{e['sha256']}\n" for e in files)
index_fp = hashlib.sha256(canonical.encode()).hexdigest()
fasta = sorted(f for f in os.listdir(os.path.join(out, 'fasta')) if f.endswith('.fa'))
version = subprocess.run([mg, '--version'], capture_output=True, text=True).stdout.strip().splitlines()
manifest = {
    'format': 'metagraph-index-manifest',
    'version': 1,
    'index_fp': index_fp,
    'files': files,
    # metadata, not identity: rebuilding the same bytes with another binary is the same index
    'builder': {'script': 'scripts/traversal/build_mini_refseq.sh',
                'metagraph_version': version[0] if version else '', 'k': k, 'mode': 'basic'},
    'inputs': {'accession_table': {'path': os.path.basename(tsv), 'sha256': digest(tsv)},
               'fasta': [entry(os.path.join('fasta', f)) for f in fasta]},
}
path = os.path.join(out, final_base + '.manifest.json')
with open(path, 'w') as f:
    json.dump(manifest, f, indent=1)
    f.write('\n')
print(f'{path}: index_fp {index_fp}')
PYEOF
}

if [ "$MANIFEST_ONLY" = "1" ]; then
    write_index_manifest
    exit 0
fi

# ---------------------------------------------------------------------------
# 0. clean previous index artifacts (keep work/ and, with REUSE_RAW=1, raw/)
# ---------------------------------------------------------------------------
log "output dir $OUT  (metagraph: $MG, threads: $THREADS)"
rm -rf "$OUT"/graph_k31.dbg* "$OUT"/annotation.* "$OUT"/rd_columns "$OUT"/intermediate \
       "$OUT"/fasta "$OUT"/sanity "$OUT"/rename_cols.txt "$OUT"/manifest.json "$OUT"/temp_*
[ "$REUSE_RAW" = "1" ] || rm -rf "$OUT/raw"
mkdir -p "$OUT/fasta" "$OUT/rd_columns" "$OUT/intermediate" "$OUT/sanity"
: > "$CMDLOG"
note "metagraph: $MG ($("$MG" --version 2>&1 | head -1))"
note "accession table: $TSV"

# ---------------------------------------------------------------------------
# 1+2. download records, rewrite headers to bare accession.version, group per taxid
# ---------------------------------------------------------------------------
log "download RefSeq records from NCBI E-utilities"
note "python3 fetch (embedded): esummary verification + efetch rettype=fasta, grouped into fasta/<taxid>.fa"
python3 - "$TSV" "$OUT" <<'PYEOF'
import csv, json, os, sys, time, urllib.request, urllib.parse, urllib.error

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
MIN_INTERVAL = 0.4          # <= 3 requests / s without an API key
_last = [0.0]


def eutils(endpoint, params, retries=6):
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        wait = MIN_INTERVAL - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        try:
            req = urllib.request.Request(EUTILS + endpoint, data=data,
                                         headers={'User-Agent': 'metagraph-mini-refseq/1.0'})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            _last[0] = time.time()
            return body
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as e:
            _last[0] = time.time()
            print(f'  [{endpoint}] attempt {attempt + 1} failed: {e}', file=sys.stderr)
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f'{endpoint} failed after {retries} attempts')


def parse_fasta(text):
    recs, hdr, seq = [], None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if hdr is not None:
                recs.append((hdr, ''.join(seq)))
            hdr, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, ''.join(seq)))
    return recs


def main(tsv, out_dir):
    rows = list(csv.DictReader(open(tsv), delimiter='\t'))
    accs = [r['accession'] for r in rows]
    expected = {r['accession']: r for r in rows}
    raw_dir = os.path.join(out_dir, 'raw')
    fasta_dir = os.path.join(out_dir, 'fasta')
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(fasta_dir, exist_ok=True)

    # esummary: verify length + taxid of every accession against the table
    summary = {}
    for i in range(0, len(accs), 100):
        d = json.loads(eutils('esummary.fcgi', {'db': 'nuccore', 'retmode': 'json',
                                                 'id': ','.join(accs[i:i + 100])}))
        for uid in d['result']['uids']:
            x = d['result'][uid]
            summary[x['accessionversion']] = x
    bad = []
    for a, r in expected.items():
        s = summary.get(a)
        if s is None:
            bad.append(f'{a}: not found in esummary')
        elif int(s['slen']) != int(r['length']) or str(s['taxid']) != str(r['taxid']):
            bad.append(f"{a}: esummary slen={s['slen']} taxid={s['taxid']} vs table {r['length']}/{r['taxid']}")
    if bad:
        print('\n'.join(bad), file=sys.stderr)
        sys.exit('esummary does not match the accession table')
    print(f'esummary OK for {len(accs)} accessions, total {sum(int(r["length"]) for r in rows)} bp')

    # efetch FASTA in batches of 10 (cached per record under raw/)
    todo = [a for a in accs if not os.path.exists(os.path.join(raw_dir, a + '.fa'))]
    for i in range(0, len(todo), 10):
        batch = todo[i:i + 10]
        text = eutils('efetch.fcgi', {'db': 'nuccore', 'rettype': 'fasta', 'retmode': 'text',
                                      'id': ','.join(batch)}).decode()
        got = {hdr.split()[0]: seq for hdr, seq in parse_fasta(text)}
        for a in batch:
            if a not in got:
                sys.exit(f'{a}: missing from efetch response')
            seq = got[a]
            if len(seq) != int(expected[a]['length']):
                sys.exit(f'{a}: fetched {len(seq)} bp, expected {expected[a]["length"]}')
            with open(os.path.join(raw_dir, a + '.fa'), 'w') as f:
                f.write('>' + a + '\n')              # bare accession.version header
                for j in range(0, len(seq), 70):
                    f.write(seq[j:j + 70] + '\n')
        print(f'  fetched {min(i + 10, len(todo))}/{len(todo)}', flush=True)

    # group per taxid (table order), forward strand as downloaded
    by_tax = {}
    for r in rows:
        by_tax.setdefault(r['taxid'], []).append(r['accession'])
    for taxid, members in by_tax.items():
        path = os.path.join(fasta_dir, taxid + '.fa')
        with open(path, 'w') as out:
            for a in members:
                with open(os.path.join(raw_dir, a + '.fa')) as f:
                    out.write(f.read())
        print(f'  {path}: {len(members)} records')
    print(f'wrote {len(by_tax)} taxid FASTA files to {fasta_dir}')


main(sys.argv[1], sys.argv[2])
PYEOF

FASTAS=()
while IFS= read -r f; do FASTAS+=("$f"); done < <(ls "$OUT"/fasta/*.fa | sort)
echo "input FASTA files (${#FASTAS[@]}):"; printf '  %s\n' "${FASTAS[@]}"

# ---------------------------------------------------------------------------
# 3. graph: k=31, BASIC, succinct, forward strand only, indexed suffix ranges,
#    dummy k-mers kept (no .edgemask, like refseq33m)
# ---------------------------------------------------------------------------
log "build graph"
run "$MG" build -v -p "$THREADS" --mode basic --graph succinct --state stat -k "$K" \
    --index-ranges 12 --mem-cap-gb 2 -o "$OUT/graph_k31" "${FASTAS[@]}" 2>&1 | grep -v 'dummy prefix of length'
[ -f "$GRAPH" ] || { echo "graph not built" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 4. annotation: one column per taxid file, with k-mer coordinates
#    (-p 1 so the column order is the command-line order and the fixture is deterministic)
# ---------------------------------------------------------------------------
log "annotate with coordinates"
run "$MG" annotate -v -p 1 -i "$GRAPH" --anno-filename --coordinates --mem-cap-gb 2 \
    -o "$OUT/$ANNO_BASE" "${FASTAS[@]}" 2>&1 | grep -v -E '^\s*metagraph annotate|To enable per-sequence'

log "row_diff stages 0/1/2 (with coordinates)"
for stage in 0 1 2; do
    run "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff --row-diff-stage "$stage" --coordinates \
        -i "$GRAPH" --mem-cap-gb 4 -o "$OUT/rd_columns/out" "$OUT/$ANNO_BASE.column.annodbg" 2>&1 | tail -2
done
rm -f "$OUT/rd_columns/$ANNO_BASE.row_count" "$OUT/rd_columns/$ANNO_BASE.row_reduction"
for f in anchors rd_succ; do [ -f "$GRAPH.$f" ] || { echo "missing $GRAPH.$f" >&2; exit 1; }; done

log "row_diff_brwt_coord"
run "$MG" transform_anno -v -p "$THREADS" --anno-type row_diff_brwt_coord --greedy -i "$GRAPH" \
    -o "$OUT/$ANNO_BASE" "$OUT/rd_columns/$ANNO_BASE.column.annodbg" 2>&1 | tail -2

log "relax BRWT"
run "$MG" relax_brwt -v -p "$THREADS" --relax-arity 32 -o "$OUT/$ANNO_BASE.relaxed" \
    "$OUT/$ANNO_BASE.row_diff_brwt_coord.annodbg" 2>&1 | tail -1

log "relabel columns to bare taxid"
"$MG" stats -a "$OUT/$ANNO_BASE.relaxed.row_diff_brwt_coord.annodbg" --print-col-names "$GRAPH" 2>/dev/null \
    | sed -n '/^Number of columns/,$p' | grep '\.fa$' \
    | awk '{n=$1; sub(/.*\//,"",n); sub(/\.fa$/,"",n); print $1, n}' > "$OUT/rename_cols.txt"
note "rename rules written to $OUT/rename_cols.txt (<fasta path> <taxid>)"
cat "$OUT/rename_cols.txt"
run "$MG" transform_anno -v -p "$THREADS" --rename-cols "$OUT/rename_cols.txt" -o "$OUT/$FINAL_BASE" \
    "$OUT/$ANNO_BASE.relaxed.row_diff_brwt_coord.annodbg" 2>&1 | tail -1
[ -f "$FINAL_ANNO" ] || { echo "final annotation not built" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 5. CoordToHeader sidecar, FASTA files passed in the FINAL column order
# ---------------------------------------------------------------------------
log "CoordToHeader sidecar (.seqs) in final column order"
"$MG" stats -a "$FINAL_ANNO" --print-col-names "$GRAPH" 2>/dev/null \
    | sed -n '/^Number of columns/,$p' | tail -n +2 > "$OUT/sanity/column_order.txt"
ORDERED=()
while IFS= read -r t; do [ -n "$t" ] && ORDERED+=("$OUT/fasta/$t.fa"); done < "$OUT/sanity/column_order.txt"
[ "${#ORDERED[@]}" -eq "${#FASTAS[@]}" ] || { echo "column count mismatch" >&2; exit 1; }
echo "final column order:"; printf '  %s\n' "${ORDERED[@]}"
run "$MG" annotate -v -p "$THREADS" -i "$GRAPH" --anno-filename --index-header-coords \
    -o "$OUT/$FINAL_BASE" "${ORDERED[@]}" 2>&1 | tail -1
[ -f "$OUT/$FINAL_BASE.seqs" ] || { echo "sidecar not built" >&2; exit 1; }

# move construction-only intermediates out of the way (refseq33m ships only the final files)
mv "$OUT/$ANNO_BASE.column.annodbg" "$OUT/$ANNO_BASE.column.annodbg.coords" "$OUT/rd_columns" \
   "$OUT/$ANNO_BASE.row_diff_brwt_coord.annodbg" "$OUT/$ANNO_BASE.relaxed.row_diff_brwt_coord.annodbg" \
   "$OUT/$ANNO_BASE.linkage" "$OUT/rename_cols.txt" "$OUT/intermediate/" 2>/dev/null || true
mv "$GRAPH.succ" "$GRAPH.pred" "$GRAPH.succ_boundary" "$GRAPH.pred_boundary" "$OUT/intermediate/" 2>/dev/null || true

# ---------------------------------------------------------------------------
# 6. sanity checks with blaNDM-1 (forward + reverse complement)
# ---------------------------------------------------------------------------
log "sanity queries"
cat > "$OUT/sanity/ndm1.fa" <<'EOF'
>ndm1
ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGCATTAGCCGCTGCATTGATGCTGAGCGGGTGCATGCCCGGTGAAATCCGCCCGACGATTGGCCAGCAAATGGAAACTGGCGACCAACGGTTTGGCGATCTGGTTTTCCGCCAGCTCGCACCGAATGTCTGGCAGCACACTTCCTATCTCGACATGCCGGGTTTCGGGGCAGTCGCTTCCAACGGTTTGATCGTCAGGGATGGCGGCCGCGTGCTGGTGGTCGATACCGCCTGGACCGATGACCAGACCGCCCAGATCCTCAACTGGATCAAGCAGGAGATCAACCTGCCGGTCGCGCTGGCGGTGGTGACTCACGCGCATCAGGACAAGATGGGCGGTATGGACGCGCTGCATGCGGCGGGGATTGCGACTTATGCCAATGCGTTGTCGAACCAGCTTGCCCCGCAAGAGGGGATGGTTGCGGCGCAACACAGCCTGACTTTCGCCGCCAATGGCTGGGTCGAACCAGCAACCGCGCCCAACTTTGGCCCGCTCAAGGTATTTTACCCCGGCCCCGGCCACACCAGTGACAATATCACCGTTGGGATCGACGGCACCGACATCGCTTTTGGTGGCTGCCTGATCAAGGACAGCAAGGCCAAGTCGCTCGGCAATCTCGGTGATGCCGACACTGAGCACTACGCCGCGTCAGCGCGCGCGTTTGGTGCGGCGTTCCCCAAGGCCAGCATGATCGTGATGAGCCATTCCGCCCCCGATAGCCGCGCCGCAATCACTCATACGGCCCGCATGGCCGACAAGCTGCGCTGA
EOF
python3 - "$OUT/sanity/ndm1.fa" "$OUT/sanity/ndm1_rc.fa" <<'PYEOF'
import sys
seq = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>'))
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
open(sys.argv[2], 'w').write('>ndm1_rc\n' + ''.join(comp[c] for c in reversed(seq)) + '\n')
PYEOF

S="$OUT/sanity"
run "$MG" query -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1.fa"    > "$S/q_fw_default.tsv" 2>/dev/null
run "$MG" query -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1_rc.fa" > "$S/q_rc_default.tsv" 2>/dev/null
run "$MG" query --query-mode coords -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1.fa"    > "$S/q_fw_coords.tsv" 2>/dev/null
run "$MG" query --query-mode coords -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1_rc.fa" > "$S/q_rc_coords.tsv" 2>/dev/null
run "$MG" query --json --query-mode coords -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1.fa" > "$S/q_fw_coords.json" 2>/dev/null
run "$MG" query --query-mode signature --min-kmers-fraction-label 0.001 -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1.fa"    > "$S/q_fw_signature.tsv" 2>/dev/null
run "$MG" query --query-mode signature --min-kmers-fraction-label 0.001 -i "$GRAPH" -a "$FINAL_ANNO" "$S/ndm1_rc.fa" > "$S/q_rc_signature.tsv" 2>/dev/null
run "$MG" stats -a "$FINAL_ANNO" --print-col-names "$GRAPH" > "$S/stats.txt" 2>&1
run "$MG" stats "$OUT/$FINAL_BASE.seqs" >> "$S/stats.txt" 2>&1
run "$MG" stats --count-dummy "$GRAPH" 2>&1 | grep -i dummy >> "$S/stats.txt"

echo "-- default mode, forward:";  tr '\t' '\n' < "$S/q_fw_default.tsv" | tail -n +3 | head -5; echo "   ..."
echo "-- coords mode, forward:";   tr '\t' '\n' < "$S/q_fw_coords.tsv"  | tail -n +3 | head -3; echo "   ..."

# ---------------------------------------------------------------------------
# 7. verify against refseq33m masks + write manifest.json
# ---------------------------------------------------------------------------
log "verify masks vs refseq33m and write manifest"
python3 - "$TSV" "$OUT" "$MG" "$K" <<'PYEOF'
import csv, json, os, re, subprocess, sys, time
tsv, out, mg, k = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
S = os.path.join(out, 'sanity')
rows = list(csv.DictReader(open(tsv), delimiter='\t'))

def parse_signature(path):
    res = {}
    for line in open(path):
        for p in line.rstrip('\n').split('\t')[2:]:
            m = re.match(r'<(.+?)>:(\d+):([xo0-9]+):(\d+)$', p)
            if m:
                res[m.group(1)] = {'kmer_count': int(m.group(2)), 'mask': m.group(3), 'score': int(m.group(4))}
    return res

def parse_default(path):
    return {m.group(1): int(m.group(2))
            for line in open(path) for p in line.rstrip('\n').split('\t')[2:]
            for m in [re.match(r'<(.+?)>:(\d+)$', p)] if m}

fw, rc = parse_signature(f'{S}/q_fw_signature.tsv'), parse_signature(f'{S}/q_rc_signature.tsv')
fw_def, rc_def = parse_default(f'{S}/q_fw_default.tsv'), parse_default(f'{S}/q_rc_default.tsv')
problems, per_acc = [], []
for r in rows:
    a = r['accession']
    exp_fw, exp_rc = r['refseq33m_fw_mask'], r['refseq33m_rc_mask']
    got_fw = fw.get(a, {}).get('mask', '-'); got_rc = rc.get(a, {}).get('mask', '-')
    ok = (got_fw == exp_fw) and (got_rc == exp_rc)
    if not ok:
        problems.append(f'{a} ({r["group"]}): fw {got_fw} vs refseq33m {exp_fw}; rc {got_rc} vs {exp_rc}')
    if r['group'] == 'negative' and (a in fw or a in rc):
        problems.append(f'{a}: negative control carries NDM k-mers')
    per_acc.append({**{k_: r[k_] for k_ in ('accession', 'taxid', 'length', 'group', 'title')},
                    'refseq33m_fw_mask': exp_fw, 'refseq33m_rc_mask': exp_rc,
                    'mini_fw_mask': got_fw, 'mini_rc_mask': got_rc,
                    'mini_fw_kmers': fw_def.get(a, 0), 'mini_rc_kmers': rc_def.get(a, 0), 'mask_match': ok})
labels_are_accessions = all(re.match(r'^(NZ_|NC_|NG_)', l) for l in list(fw) + list(rc))
if not labels_are_accessions:
    problems.append('query labels are not accession headers (is the .seqs sidecar in place?)')
n_fw_full = sum(1 for v in fw.values() if v['mask'] == 'x783')
n_rc_full = sum(1 for v in rc.values() if v['mask'] == 'x783')
stats_txt = open(f'{S}/stats.txt').read()
cols = [l.strip() for l in open(f'{S}/column_order.txt') if l.strip()]
sanity = {
    'query': 'blaNDM-1 (813 bp, 783 k-mers)',
    'labels_are_accession_headers': labels_are_accessions,
    'forward_hits': len(fw), 'forward_full_x783': n_fw_full,
    'rc_hits': len(rc), 'rc_full_x783': n_rc_full,
    'masks_identical_to_refseq33m': len(rows) - len([p for p in problems if 'refseq33m' in p]),
    'negative_controls_clean': not any('negative control' in p for p in problems),
    'graph_k': re.search(r'^k: (\d+)', stats_txt, re.M).group(1),
    'graph_mode': re.search(r'^mode: (\w+)', stats_txt, re.M).group(1),
    'graph_nodes': int(re.search(r'^nodes \(k\): (\d+)', stats_txt, re.M).group(1)),
    'indexed_suffix_length': int(re.search(r'indexed suffix length: (\d+)', stats_txt).group(1)),
    'num_columns': len(cols), 'column_order': cols,
    'seqs_total_sequences': int(re.search(r'total sequences: (\d+)', stats_txt).group(1)),
    'seqs_total_kmers': int(re.search(r'total k-mers: (\d+)', stats_txt).group(1)),
    'problems': problems,
}
json.dump(sanity, open(f'{S}/summary.json', 'w'), indent=1)
cmds = [l.rstrip('\n') for l in open(f'{out}/commands.log')]
version = subprocess.run([mg, '--version'], capture_output=True, text=True).stdout.strip().splitlines()[0]
manifest = {
    'name': 'mini_refseq', 'mirrors': 'refseq33m (DBGSuccinct k=31 BASIC, RowDiff<BRWT> coordinates, taxid columns, CoordToHeader .seqs)',
    'built_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'metagraph_version': version, 'k': k, 'mode': 'basic',
    'graph': 'graph_k31.dbg', 'graph_sidecars': ['graph_k31.dbg.anchors', 'graph_k31.dbg.rd_succ'],
    'annotation': 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg', 'coord_to_header': 'annotation.relaxed.relabeled.seqs',
    'fasta_dir': 'fasta', 'num_records': len(rows), 'total_bp': sum(int(r['length']) for r in rows),
    'groups': {g: sum(1 for r in rows if r['group'] == g) for g in ('ndm_fw', 'ndm_rc', 'ndm_snp', 'ndm_both', 'negative')},
    'columns': cols, 'records': per_acc, 'commands': cmds, 'sanity': sanity,
}
json.dump(manifest, open(f'{out}/manifest.json', 'w'), indent=1)
print(json.dumps({k_: v for k_, v in sanity.items() if k_ not in ('column_order',)}, indent=1))
if problems:
    print('\n'.join(problems), file=sys.stderr)
    sys.exit('SANITY CHECK FAILED')
print(f'\nOK: {out}/manifest.json written')
PYEOF

write_index_manifest

log "done"
ls -la "$OUT"
