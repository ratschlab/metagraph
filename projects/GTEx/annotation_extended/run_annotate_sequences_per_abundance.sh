#!/bin/bash

set -e

### source paths 
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K="$1"

mem=100000
threads=8
pmem=$(($mem / $threads))

graph=${basedir}/gtex/output_k${K}_trimmed_clean_graph_extended_chunked/graph_k${K}.dbg
seqdir=${basedir}/gtex/output_k${K}_trimmed_clean
seqfile=${basedir}/gtex/output_k${K}_trimmed_clean_graph_extended_chunked/all_input_files.txt

outdir=${basedir}/gtex/output_k${K}_trimmed_extended_annotation
mkdir -p $outdir

TMP=temp_k${K}_annoseq_$(date +%Y%m%d_%H:%M:%S)
mkdir -p $TMP
cd $TMP
cp $seqfile files_to_annotate.txt
for fname in $(find ${seqdir} -name \*.fasta.gz)
do
    fbase=$(basename $fname)
    echo $fname >> files_to_annotate.txt
done
split -l $((${threads} * 8)) files_to_annotate.txt

for fname in x*
do
    echo "cat $(pwd)/${fname} | /usr/bin/time -v $metagraph annotate -v -i ${graph} --outfile-base ${outdir}/ --parallel $threads --anno-filename --separately 2>&1 | tee $(pwd)/annotate_gtex.k${K}.${fname}.log" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $(pwd)/annotate_gtex.k${K}.${fname}.lsf.log
done
