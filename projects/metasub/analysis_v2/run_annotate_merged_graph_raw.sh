#!/bin/bash

set -e

K=19

#mem=240000
#threads=18
mem=500000
threads=24
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_raw
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
TMP=temp_raw_k${K}_${RANDOM}
mkdir -p $TMP
cd $TMP
for fname in $(tail -n+2 ../../complete_metadata_extended.clean.v2.csv | cut -f 41 -d ',' | tr ':' '\n')
do
    fbase=$(basename $fname)
    if [ ! -f ${outdir}/${fbase}.column.annodbg ]
    then
        echo $fname >> files_to_annotate.txt
    fi
done
split -l 100 files_to_annotate.txt

for fname in x*
do
    echo "cat $(pwd)/${fname} | /usr/bin/time -v $metagraph annotate -v --fwd-and-reverse -i ${graphdir}/graph_merged_k${K}.dbg --outfile-base ${outdir}/ --parallel $threads --anno-filename --separately 2>&1 | tee $(pwd)/annotate_metasub.k${K}.${fname}.log" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $(pwd)/annotate_metasub.k${K}.${fname}.lsf.log
done
