#!/usr/bin/env bash
#
# fetch_mini_refseq.sh -- fetch the mini RefSeq index the traversal end-to-end tests run on.
#
# The index (42 blaNDM-carrying RefSeq records, ~9.3 Mbp; k = 31 basic DBGSuccinct,
# RowDiff<BRWT> with k-mer coordinates, a .seqs header index; the format of the public
# refseq33m index) is published as the pre-release test-data-mini-refseq-v1 of
# ratschlab/metagraph. This script downloads it, checks its pinned sha256, unpacks it, and
# writes the blaNDM-1 query of the tests (sanity/ndm1.fa).
#
# usage: fetch_mini_refseq.sh [DEST]
#   DEST  default: <repo>/metagraph/build/mini_refseq, where the unit tests look for it
#         (mini_refseq/ in their working directory)
# env:
#   METAGRAPH_TEST_DATA  where the archive is kept between runs
#                        (default: ${XDG_CACHE_HOME:-~/.cache}/metagraph-test-data)
#
# Idempotent: a cached archive with the pinned sha256 is not downloaded again, and an existing
# DEST is left as it is.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_MG_DIR="$(cd "$HERE/../.." && pwd)"

RELEASE_URL="https://github.com/ratschlab/metagraph/releases/download/test-data-mini-refseq-v1"
ASSET="mini_refseq-bea44d60.tar.gz"
ASSET_SHA256="1b9e07dd6d1767466205a40c5739ce47c86255abab65df71cf6a46b0e4f6fd8b"

DEST="${1:-$REPO_MG_DIR/build/mini_refseq}"
CACHE="${METAGRAPH_TEST_DATA:-${XDG_CACHE_HOME:-$HOME/.cache}/metagraph-test-data}"

if [ -e "$DEST" ]; then
    echo "exists, left as it is: $DEST"
    exit 0
fi

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1"; else shasum -a 256 "$1"; fi \
        | cut -d' ' -f1
}

mkdir -p "$CACHE"
ARCHIVE="$CACHE/$ASSET"
if [ ! -f "$ARCHIVE" ] || [ "$(sha256_of "$ARCHIVE")" != "$ASSET_SHA256" ]; then
    echo "download $RELEASE_URL/$ASSET"
    part="$ARCHIVE.part.$$"
    trap 'rm -f "$part"' EXIT
    curl -fsSL --retry 5 --retry-delay 3 --retry-all-errors --connect-timeout 30 \
        -o "$part" "$RELEASE_URL/$ASSET"
    got="$(sha256_of "$part")"
    if [ "$got" != "$ASSET_SHA256" ]; then
        echo "$ASSET: sha256 $got, expected $ASSET_SHA256" >&2
        exit 1
    fi
    mv "$part" "$ARCHIVE"
    trap - EXIT
fi

# unpacked beside DEST, then renamed: DEST is never left half written
mkdir -p "$(dirname "$DEST")"
WORK="$(mktemp -d "$(dirname "$DEST")/.fetch_mini_refseq.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$ARCHIVE" -C "$WORK"
[ -d "$WORK/mini_refseq" ] || { echo "$ASSET holds no mini_refseq/" >&2; exit 1; }

# the blaNDM-1 gene (813 bp), the query of the tests
mkdir -p "$WORK/mini_refseq/sanity"
cat > "$WORK/mini_refseq/sanity/ndm1.fa" <<'EOF'
>ndm1
ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGCATTAGCCGCTGCATTGATGCTGAGCGGGTGCATGCCCGGTGAAATCCGCCCGACGATTGGCCAGCAAATGGAAACTGGCGACCAACGGTTTGGCGATCTGGTTTTCCGCCAGCTCGCACCGAATGTCTGGCAGCACACTTCCTATCTCGACATGCCGGGTTTCGGGGCAGTCGCTTCCAACGGTTTGATCGTCAGGGATGGCGGCCGCGTGCTGGTGGTCGATACCGCCTGGACCGATGACCAGACCGCCCAGATCCTCAACTGGATCAAGCAGGAGATCAACCTGCCGGTCGCGCTGGCGGTGGTGACTCACGCGCATCAGGACAAGATGGGCGGTATGGACGCGCTGCATGCGGCGGGGATTGCGACTTATGCCAATGCGTTGTCGAACCAGCTTGCCCCGCAAGAGGGGATGGTTGCGGCGCAACACAGCCTGACTTTCGCCGCCAATGGCTGGGTCGAACCAGCAACCGCGCCCAACTTTGGCCCGCTCAAGGTATTTTACCCCGGCCCCGGCCACACCAGTGACAATATCACCGTTGGGATCGACGGCACCGACATCGCTTTTGGTGGCTGCCTGATCAAGGACAGCAAGGCCAAGTCGCTCGGCAATCTCGGTGATGCCGACACTGAGCACTACGCCGCGTCAGCGCGCGCGTTTGGTGCGGCGTTCCCCAAGGCCAGCATGATCGTGATGAGCCATTCCGCCCCCGATAGCCGCGCCGCAATCACTCATACGGCCCGCATGGCCGACAAGCTGCGCTGA
EOF

mv "$WORK/mini_refseq" "$DEST"
echo "unpacked $ASSET into $DEST"
