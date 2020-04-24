#!/bin/bash

set -e

K=41

mem=150000
threads=8
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome
seqdir1=${basedir}/data/gnomad/release_3.0
seqdir2=${basedir}/data/gtex/graphs/output_k41_trimmed_clean_graph_chunked
seqdir3=${basedir}/raw_data/ref_genomes/hg38
seqdir4=${basedir}/raw_data/ref_genomes/gencode_v32
outdir=${basedir}/data/gtex/graphs/output_k${K}_trimmed_clean_graph_extended_chunked
mkdir -p $outdir

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

input=${outdir}/all_input_files.txt
if [ ! -f ${input} ]
then
    find ${seqdir1} -name \*.fasta.gz > ${input}
    ls -1 ${seqdir2}/*.fasta.gz >> ${input}
    ls -1 ${seqdir3}/GRCh38.primary_assembly.genome.fa.gz >> ${input}
    ls -1 ${seqdir4}/gencode.v32.transcripts.fa.gz >> ${input}
fi

for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}
do
    ### get rid of all combinations that are not allowed
    if [ -z $(echo $F | grep -v -E '[ACGT]\\\$(\\\$)+$' | grep -v -E '(\\\$)+[ACGT]+(\\\$)+$' | grep -v -E '[ACGT]+\\\$(\\\$)*[ACGT]+') ]
    then
        continue
    fi

    FF=$(eval "echo $F")
    if [ -f "${outdir}/graph_merged_k${K}.${FF}.dbg.chunk" ]
    then
        echo chunk $FF exists
    else
        echo "submitting $FF"
        echo "cat ${input} | /usr/bin/time -v ${metagraph} build -v --parallel ${threads} -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_k${K} --suffix $F 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J gtxEx_${FF} -oo ${outdir}/build_k${K}_$FF.lsf.log -We 22:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done

