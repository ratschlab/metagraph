#!/usr/bin/env bash


if [ $# -ne 1 ]; then
  echo "Usage: $0 <FILE.sra>"
  echo "Example: $0 path/to/file/DRR042288.sra"
  exit 1
fi


FILE="$1"
OUTDIR="$(dirname $FILE)"


if [ ! -f $FILE ] || ! vdb-validate -x $FILE 2> /dev/null; then
  echo "ERROR: $FILE file corrupt"
  exit 1
fi

mkdir -p "$OUTDIR/statistics"
sra-stat --xml --statistics "$FILE" > "$OUTDIR/statistics/$(basename ${FILE%.sra}).xml"

fastq-dump --split-spot --skip-technical --gzip -O "$OUTDIR" "$FILE" && rm "$FILE"
#fastq-dump --split-spot --skip-technical -Z "$FILE" | gzip -9 > "${FILE%.sra}_.fastq.gz"
