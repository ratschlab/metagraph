#!/bin/bash

set -e

conda activate  /cluster/work/grlab/projects/GTEx/conda/761e4791

infile=/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all.fa
outfile=/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all.mapped_genome

STAR --genomeDir /cluster/work/grlab/projects/GTEx/genome/GRCh38.p12.star_index --readFilesIn $infile --runThreadN 2 --sjdbOverhang 100 --outFileNamePrefix $outfile --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand cat --outFilterMatchNmin 21 --outFilterScoreMinOverLread 0.33 --outSAMstrandField intronMotif --outSAMmode Full --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --limitSjdbInsertNsj 2000000 --outSAMtype BAM Unsorted --outSAMheaderHD @HD VN:1.4 --outSAMmultNmax 1 --twopassMode Basic
