#!/usr/bin/env bash


if [ $# -ne 2 ]; then
  echo "Usage: $0 <DUMP_URL> <DEST>" >&2
  echo "Example: $0 \
    \"ftp://ftp.ebi.ac.uk/pub/software/bigsi/nat_biotech_2018/ctx/DRR000/DRR000001/cleaned/DRR000001.ctx.bz2\" \
    /cluster/work/grlab/projects/metagenome/data/BIGSI/dumps" >&2
  exit 1
fi

URL="$1"
FILE="$(basename $URL)"
DEST="$2"

# check if the file already exists
if [ -f $DEST/$FILE ]; then
  echo "Already downloaded: $FILE"
  exit 0
fi

# download file to buffer zone
if ! wget -q -O "$DEST/$FILE" "$URL"; then
  echo "ERROR: Download failed." >&2
  echo "Command: wget -q -O $DEST/$FILE $URL" >&2
  rm "$DEST/$FILE"
  exit 1
fi

echo "Downloaded: $FILE"
