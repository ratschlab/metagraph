#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

mem=80000
threads=4
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation
seqdir=${basedir}/metasub/graphs/output_k${K}_cleaned
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph
TMP=temp_k${K}_$(date +%Y%m%d_%H:%M:%S)
mkdir -p $TMP
cd $TMP
touch files_to_annotate.txt
for fname in $(find ${seqdir} -name \*.fasta.gz)
do
    fbase=$(basename $fname)
    if [ ! -f ${outdir}/${fbase}.column.annodbg ]
    then
        echo $fname >> files_to_annotate.txt
    fi
done
split -l $(($threads * 4)) files_to_annotate.txt

for fname in x*
do
    echo "cat $(pwd)/${fname} | /usr/bin/time -v $metagraph annotate -v -i ${graphdir}/graph_merged_k${K}.dbg --outfile-base ${outdir}/ --parallel $threads --anno-filename --separately 2>&1 | tee $(pwd)/annotate_metasub.k${K}.${fname}.log" | bsub -J ms_anno_k${K} -We 24:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $(pwd)/annotate_metasub.k${K}.${fname}.lsf.log
done
