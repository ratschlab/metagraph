#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <ID>"
  exit 1
fi

OUT="$HOME/metagenome/benchmark_TCGA/nobackup"

if [ -d "$OUT/${1}" ]; then
  exit 0
fi

./gdc-client download -t ~/gdc-user-token.2018-07-03T15_27_25.227Z.txt -n 10 -d ${OUT}/ ${1} > /dev/null

for x in $OUT/${1}/*.bam; do
  if [ -f $x ] && [ ! -f ${x%.bam}.fa.gz ]; then
    bsub -J "process_bam_${dir}_bam" -W 72:00 -n 1 -R "rusage[mem=1000]" "./process_bam.sh $x"
  fi
done
