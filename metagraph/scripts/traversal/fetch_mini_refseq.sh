#!/usr/bin/env bash
#
# fetch_mini_refseq.sh -- fetch the published copy of the mini index instead of building it.
#
# build_mini_refseq.sh builds the mini index from NCBI's records; a rebuild reproduces the
# graph but not the annotation byte for byte (its row-diff anchors are picked anew), so the
# fixtures answered on one build differ from what another build answers in the decode work.
# The copy those fixtures were answered on (index_fp bea44d60...) is published as the
# pre-release test-data-mini-refseq-v1 of ratschlab/metagraph; this script downloads it,
# checks it and unpacks it where the tests look for the mini index, without NCBI.
#
# usage: fetch_mini_refseq.sh [DEST]
#        fetch_mini_refseq.sh --asset      print the archive's name (a cache key) and exit
#   DEST  default: <repo>/metagraph/build/mini_refseq, where every test and script looks by
#         default (the unit and integration tests: mini_refseq/ in their working directory,
#         build/; the integration and Python tests also take $METAGRAPH_MINI_REFSEQ)
# env:
#   METAGRAPH_TEST_DATA  where the archive is kept between runs
#                        (default: ${XDG_CACHE_HOME:-~/.cache}/metagraph-test-data)
#
# Idempotent: an archive in METAGRAPH_TEST_DATA with the pinned sha256 is not downloaded
# again, and a DEST that already holds this index is left as it is. A DEST holding anything
# else (a local build, say) is refused, never overwritten.
#
# The archive holds the index bundle, its index manifest, the build record and fasta/. What
# else of build_mini_refseq.sh's output the repository reads is written here as that script
# writes it: the blaNDM-1 query of the MiniRefSeq unit tests and of the real-index suite
# (sanity/ndm1.fa, ndm1_rc.fa; the sequence taken from the script) and the records one file
# each (raw/<accession>.fa, read by the real-index suite; fasta/<taxid>.fa is their
# concatenation).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_MG_DIR="$(cd "$HERE/../.." && pwd)"           # <repo>/metagraph

# the published copy: the archive's sha256 is pinned here, so a changed release fails
RELEASE_URL="https://github.com/ratschlab/metagraph/releases/download/test-data-mini-refseq-v1"
ASSET="mini_refseq-bea44d60.tar.gz"
ASSET_SHA256="1b9e07dd6d1767466205a40c5739ce47c86255abab65df71cf6a46b0e4f6fd8b"
INDEX_FP="bea44d604ffe08d60b348b630b4a35bd0ae0b428615727345ba191228c7ef5e7"
INDEX_MANIFEST="annotation.relaxed.relabeled.manifest.json"

if [ "${1:-}" = "--asset" ]; then
    echo "$ASSET"
    exit 0
fi
case "${1:-}" in -*) echo "usage: $0 [DEST] | --asset" >&2; exit 2;; esac

DEST="${1:-$REPO_MG_DIR/build/mini_refseq}"
CACHE="${METAGRAPH_TEST_DATA:-${XDG_CACHE_HOME:-$HOME/.cache}/metagraph-test-data}"

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1"; else shasum -a 256 "$1"; fi \
        | cut -d' ' -f1
}

download() {   # URL FILE
    curl -fsSL --retry 5 --retry-delay 3 --retry-all-errors --connect-timeout 30 -o "$2" "$1"
}

# Does directory $1 hold the published index? Its manifest is the pinned one, every file the
# manifest lists (the bundle and the source FASTA) has the size and sha256 listed, and the
# files derived from build_mini_refseq.sh and fasta/ (sanity/, raw/) are as derived. With
# --write-derived, those are written (into a fresh extraction) instead of compared. Prints
# what differs; python3 -I because the files it reads come from a download.
check_index() {   # DIR [--write-derived]
    python3 -I - "$1" "$INDEX_MANIFEST" "$INDEX_FP" "$HERE/build_mini_refseq.sh" "${2:-}" <<'PYEOF'
import hashlib, json, os, re, sys

root, manifest_name, index_fp, build_script, mode = sys.argv[1:6]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


problems = []
try:
    with open(os.path.join(root, manifest_name)) as f:
        manifest = json.load(f)
except (OSError, ValueError) as e:
    sys.exit(f'{root}: no index manifest ({e})')
# index_fp as build_mini_refseq.sh and `metagraph --index-manifest` compute it
canonical = ''.join(f"{e['path']}\t{e['size']}\t{e['sha256']}\n"
                    for e in sorted(manifest['files'], key=lambda e: e['path'].encode()))
if manifest.get('index_fp') != index_fp or hashlib.sha256(canonical.encode()).hexdigest() != index_fp:
    problems.append(f'index_fp {manifest.get("index_fp")}, not the published {index_fp}')
fasta = manifest.get('inputs', {}).get('fasta', [])
for e in manifest['files'] + fasta:
    path = os.path.join(root, e['path'])
    if not os.path.isfile(path):
        problems.append(f'missing: {e["path"]}')
    elif os.path.getsize(path) != e['size'] or sha256(path) != e['sha256']:
        problems.append(f'differs from the manifest: {e["path"]}')
if problems:
    print('\n'.join(f'  {p}' for p in problems), file=sys.stderr)
    sys.exit(1)

# the query of build_mini_refseq.sh's sanity step, as it writes it
m = re.search(r"/sanity/ndm1\.fa\" <<'EOF'\n>ndm1\n([ACGT]+)\nEOF\n", open(build_script).read())
if not m:
    sys.exit(f'{build_script}: the blaNDM-1 query (sanity/ndm1.fa) is not where it was')
seq = m.group(1)
rc = seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))
derived = {'sanity/ndm1.fa': f'>ndm1\n{seq}\n', 'sanity/ndm1_rc.fa': f'>ndm1_rc\n{rc}\n'}
# the records as the script downloads them, one file each under raw/ (">accession" and the
# sequence in lines of 70); each fasta/<taxid>.fa is its records' files concatenated
for e in fasta:
    with open(os.path.join(root, e['path'])) as f:
        for record in re.split(r'(?m)^(?=>)', f.read()):
            if record:
                accession = record[1:record.index('\n')]
                derived[f'raw/{accession}.fa'] = record
for rel, text in derived.items():
    path = os.path.join(root, rel)
    if mode == '--write-derived':
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write(text)
    elif not os.path.isfile(path) or open(path).read() != text:
        problems.append(f'missing or different: {rel}')

if problems:
    print('\n'.join(f'  {p}' for p in problems), file=sys.stderr)
    sys.exit(1)
PYEOF
}

note_where() {
    echo "mini index (index_fp $INDEX_FP) in $DEST"
    if [ "$(cd "$DEST" && pwd -P)" != "$(cd "$REPO_MG_DIR/build/mini_refseq" 2>/dev/null && pwd -P)" ]; then
        echo "  not the default place: the integration and Python tests read it with"
        echo "  METAGRAPH_MINI_REFSEQ=$DEST, the unit tests as mini_refseq/ in their working directory"
    fi
}

# --------------------------------------------------------------------------- DEST present
if [ -d "$DEST" ] && [ -z "$(ls -A "$DEST")" ]; then
    rmdir "$DEST"                       # an empty directory holds nothing to keep
fi
if [ -e "$DEST" ]; then
    if check_index "$DEST" 2>/dev/null; then
        echo "up to date: $DEST"
        note_where
        exit 0
    fi
    echo "$DEST exists and does not hold the published mini index:" >&2
    check_index "$DEST" >&2 || true
    echo "remove it, or give another DEST" >&2
    exit 1
fi

# --------------------------------------------------------------------------- the archive
mkdir -p "$CACHE"
ARCHIVE="$CACHE/$ASSET"
if [ -f "$ARCHIVE" ] && [ "$(sha256_of "$ARCHIVE")" = "$ASSET_SHA256" ]; then
    echo "archive from the cache: $ARCHIVE"
else
    echo "download $RELEASE_URL/$ASSET"
    part="$ARCHIVE.part.$$"
    trap 'rm -f "$part" "$part.sha256"' EXIT
    download "$RELEASE_URL/$ASSET.sha256" "$part.sha256"
    # the release's own checksum file must name the pinned digest
    read -r published name < "$part.sha256" || true
    if [ "${published:-}" != "$ASSET_SHA256" ] || [ "${name:-}" != "$ASSET" ]; then
        echo "$ASSET.sha256 of the release states '${published:-} ${name:-}'," \
             "not the pinned $ASSET_SHA256" >&2
        exit 1
    fi
    download "$RELEASE_URL/$ASSET" "$part"
    got="$(sha256_of "$part")"
    if [ "$got" != "$ASSET_SHA256" ]; then
        echo "$ASSET: sha256 $got, expected $ASSET_SHA256" >&2
        exit 1
    fi
    mv "$part" "$ARCHIVE"
    mv "$part.sha256" "$ARCHIVE.sha256"
    trap - EXIT
fi

# --------------------------------------------------------------------------- unpack
# into a directory beside DEST, checked, then renamed: DEST is never left half written
mkdir -p "$(dirname "$DEST")"
WORK="$(mktemp -d "$(dirname "$DEST")/.fetch_mini_refseq.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$ARCHIVE" -C "$WORK"
[ -d "$WORK/mini_refseq" ] || { echo "$ASSET holds no mini_refseq/" >&2; exit 1; }
check_index "$WORK/mini_refseq" --write-derived || { echo "$ASSET: not the published index" >&2; exit 1; }
[ ! -e "$DEST" ] || { echo "$DEST appeared while unpacking" >&2; exit 1; }
mv "$WORK/mini_refseq" "$DEST"
echo "unpacked $ASSET"
note_where
