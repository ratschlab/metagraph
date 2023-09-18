#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K=$1

basedir=/cluster/work/grlab/projects/metagenome/data
annotation=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_clean.collect.relabeled.brwt.annodbg
graph=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked/graph_merged_k${K}.dbg
outdir=${basedir}/metasub/queries/
mkdir -p $outdir
query=/cluster/work/grlab/projects/metagenome/raw_data/AMR/CARD_2.0.2/nucleotide_fasta_protein_homolog_model.fasta

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

threads=8
mem=500000
pmem=$((${mem} / ${threads}))

log=${outdir}/amr_CARD_nucleotide_fasta_protein_homolog_model_k${K}.kmer_based.log
out=${outdir}/amr_CARD_nucleotide_fasta_protein_homolog_model_k${K}.kmer_based.tsv
echo "/usr/bin/time -v $metagraph query -i $graph -a ${annotation} --query-mode matches --min-kmers-fraction-label 0.1 -p ${threads} ${query} > $out" | bsub -J ms_q_amr${K} -oo ${log} -We 24:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
