#!/bin/bash

set -e

### source paths
# set basedir
# set metagraph
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K="$1"

mem=160000
threads=16
pmem=$(($mem / $threads))

graph=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
seqdir=${basedir}/gtex/graphs/output_k${K}_trimmed_clean
outdir=${basedir}/gtex/graphs/output_k${K}_trimmed_annotation_per_abundance_quantiles
mkdir -p $outdir
TMP=temp_k${K}_annoseq_$(date +%Y%m%d_%H:%M:%S)
mkdir -p $TMP
cd $TMP
for fname in $(find ${seqdir} -name \*.fasta.gz)
do
    fbase=$(basename $fname)
    if [ ! -f ${outdir}/${fbase}.column.annodbg ]
    then
        echo $fname >> files_to_annotate.txt
    fi
done
split -l $((${threads} * 8)) files_to_annotate.txt

for fname in x*
do
    echo "cat $(pwd)/${fname} | /usr/bin/time -v $metagraph annotate -v -i ${graph} --outfile-base ${outdir}/ --parallel $threads --anno-filename --separately 2>&1 | tee $(pwd)/annotate_gtex.k${K}.${fname}.log" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $(pwd)/annotate_gtex.k${K}.${fname}.lsf.log
done
