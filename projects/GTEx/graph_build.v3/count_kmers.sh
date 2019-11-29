#!/usr/bin/bash


KMC="$(dirname ${BASH_SOURCE[0]})/../../../metagraph/build/KMC/kmc"


if [ $# -lt 4 ]; then
    echo -e "Usage:\n$0 <k> <sampleID> <outdir> <samples>" >&2
    exit 1
fi

K="$1"
shift
SAMPLEID="$1"
shift
OUT="$1"
shift
SAMPLES="$@"
cutoff=1
q=0.05
num_threads=2

if [ $cutoff -eq 0 ]; then
  echo "Error: cutoff is too small" >&2
  exit 0
fi

### trim sample
module load seqtk
rm -f ${TMPDIR}/${SAMPLEID}.files
outfiles=""
for sample in $SAMPLES
do
    samplebase=$(basename $sample)
    out=${TMPDIR}/${samplebase%.gz}.trim_q${q}.gz
    echo $out >> ${TMPDIR}/${SAMPLEID}.files
    echo trimming $sample into $out
    outfiles="$outfiles $out"
    (seqtk trimfq -q $q $sample | gzip > ${out}) &
done
wait

mkdir -p "${TMPDIR}/${SAMPLEID}.cache"
#filename="$(basename $FILE)"
$KMC -k$K -m10 -ci$cutoff -fq -t$num_threads "@"${TMPDIR}/${SAMPLEID}.files $OUT/${SAMPLEID}.k${K} ${TMPDIR}/${SAMPLEID}.cache
rm -r "${TMPDIR}/${SAMPLEID}.cache" "${TMPDIR}/${SAMPLEID}.files"
rm $outfiles
