#!/usr/bin/env bash

MCCORTEX="/cluster/home/mikhaika/stuff/mccortex/bin/mccortex31"


if [ $# -ne 2 ]; then
  echo "Usage: $0 <INPUT_MCCORTEX_DUMP> <OUT_FASTA>"
  echo "Example: $0 /cluster/work/grlab/projects/metagenome/data/BIGSI/dumps/DRR000001.ctx \
    /cluster/work/grlab/projects/metagenome/data/BIGSI/dumps/DRR000001.unitigs.fasta.gz"
  exit 1
fi

GRAPH_DUMP="$1"
UNITIGS="${2%.gz}"

# check if the file already exists
if [ -f "${UNITIGS}.gz" ]; then
  echo "Already done: $UNITIGS"
  exit 0
fi

if [[ $GRAPH_DUMP == *.bz2 ]]; then
  bunzip2 $GRAPH_DUMP
  GRAPH_DUMP="${GRAPH_DUMP%.bz2}"
fi

SIZE="$(wc -c $GRAPH_DUMP | awk '{print int($1 / 10**6 * 2.5 + 150)}')"

if ! $MCCORTEX check -m "${SIZE}MB" "$GRAPH_DUMP"; then
  echo "ERROR: Dump broken $GRAPH_DUMP" >&2
  exit 1
fi

if ! $MCCORTEX unitigs -m "${SIZE}MB" -o "$UNITIGS" "$GRAPH_DUMP"; then
  echo "ERROR: Cannot pull unitigs" >&2
  exit 1
fi

gzip "$UNITIGS"
rm "$GRAPH_DUMP"

echo "$UNITIGS"
