#!/bin/bash

set -e

K=31

mem=320000
threads=18
mem=400000
threads=12
#mem=700000
#threads=24
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metagut/graphs/output_k${K}_graph_chunked
seqdir=${basedir}/metagut/graphs/output_k${K}
outdir=${basedir}/metagut/graphs/output_k${K}_annotation_clean
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/metagraph
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
    echo "cat $(pwd)/${fname} | /usr/bin/time -v $metagraph annotate -v -i ${graphdir}/graph_merged_k${K}.dbg --outfile-base ${outdir}/ --parallel $threads --anno-filename --separately 2>&1 | tee $(pwd)/annotate_metagut.k${K}.${fname}.log" | bsub -J mg_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $(pwd)/annotate_metagut.k${K}.${fname}.lsf.log
done
