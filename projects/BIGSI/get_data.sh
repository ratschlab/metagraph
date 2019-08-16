#!/usr/bin/env bash

FINAL_DEST="/cluster/work/grlab/projects/metagenome/data/BIGSI"
BUFFER="/cluster/work/grlab/projects/metagenome/data/BIGSI/dumps"
URL_BASE="ftp://ftp.ebi.ac.uk/pub/software/bigsi/nat_biotech_2018/ctx"

if [ $# -ne 1 ]; then
  echo "Usage: $0 <DUMP_URL>" >&2
  echo "Example: $0 DRR000001" >&2
  exit 1
fi

ID="$1"
DUMP="${ID}.ctx"
# .../DRR000/DRR000001/cleaned/DRR000001.ctx
URL="$URL_BASE/${ID:0:6}/$ID/cleaned/$DUMP"
UNITIGS="$ID.unitigs.fasta.gz"

# check if the file already exists
if [ -f $FINAL_DEST/$UNITIGS ]; then
  echo "Already done: $UNITIGS"
  exit 0
fi

if ! ./download_mccortex_graph.sh "$URL" "$BUFFER" 2>/dev/null; then
  DUMP="${DUMP}.bz2"
  URL="${URL}.bz2"
  if ! ./download_mccortex_graph.sh "$URL" "$BUFFER"; then
    echo "ERROR: Cannot download $ID" >&2
    exit 1
  fi
fi

while [ $(bjobs | wc -l) -gt 15000 ]; do
  sleep 5
done

SIZE="$(wc -c $BUFFER/$DUMP | awk '{print int($1 / 10**6 * 3 + 200)}')"

bsub -J "pull_$DUMP" \
     -oo "${BUFFER}/${DUMP}.lsf" \
     -W 3:00 -n 1 -R "rusage[mem=${SIZE}] span[hosts=1]" \
    "./pull_unitigs.sh $BUFFER/$DUMP $FINAL_DEST/$UNITIGS"
