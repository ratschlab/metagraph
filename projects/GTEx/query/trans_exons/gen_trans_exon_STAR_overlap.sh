#!/bin/bash

set -e

hitfile=/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all_hits.result_samples.txt
query=/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all.result.hits.fa
alignments=/cluster/work/grlab/projects/metagenome/data/gtex/alignments_star_trans_exons/results/alignments

while IFS=',' read -r -a line
do
    if [ "${#line[@]}" -lt "3" ]
    then
        continue
    fi
    ### get query seq
    qstring=$(echo "${line[1]}" | cut -f 1 -d ';')
    qseq=$(grep -A1 ${qstring} ${query} | tail -n 1)
    qseq_rc=$(echo $qseq | tr ACGTacgt TGCAtgca | rev)
    qseq_part=$(echo $qseq | cut -c 16-65) ### middle 50 bases
    qseq_part_rc=$(echo $qseq_rc | cut -c 16-65) ### middle 50 bases
    for i in $(seq 2 $((${#line[@]} - 1)) )
    do
        sample=$(echo ${line[${i}]} | cut -f 1 -d '>' | cut -f 2 -d '<')
        bamfile=${alignments}/${sample}.all.bam
        outfile=${alignments}/${sample}.hits/${sample}.all.${qstring}.${qseq_part}.${qseq}.sam
        if [ -f ${bamfile}.bai -a ! -f ${outfile} ]
        then
            mkdir -p $(dirname $outfile)
            echo "module load samtools; samtools view ${bamfile} | grep ${qseq_part} > $outfile" | bsub -M 2000 -n 1 -R "rusage[mem=2000]" -We 4:00 -n 1 -J trans_hits -o /dev/null
        fi
        outfile=${alignments}/${sample}.hits/${sample}.all.${qstring}.${qseq_part_rc}.${qseq_rc}.sam
        if [ -f ${bamfile}.bai -a ! -f ${outfile} ]
        then
            mkdir -p $(dirname $outfile)
            echo "module load samtools; samtools view ${bamfile} | grep ${qseq_part_rc} > $outfile" | bsub -M 2000 -n 1 -R "rusage[mem=2000]" -We 4:00 -n 1 -J trans_hits -o /dev/null
        fi
    done
done < <(cat $hitfile | tr $'\t' ',')

