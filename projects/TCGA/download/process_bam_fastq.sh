#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <FILE.bam>"
  exit 1
fi

FILE=$1

if [ ! -f $FILE ]; then
  echo file "$FILE" not found
  exit 1
fi

SFILE=${FILE%.bam}.sorted.bam
FBASE=$(basename $FILE)
STMP=/cluster/work/grlab/projects/metagenome/raw_data/tcga/tmp/${RANDOM}_${FBASE%.bam}
mkdir -p $STMP

module load samtools
### decide whether file is single or paired end
samtools view -H $FILE > ${FILE%.bam}.txt
pair_status=$(python $(pwd)/get_pair_status.py ${FILE%.bam}.txt)

### process file
if [ "$pair_status" == "paired" ]
then
    bamscript=$(pwd)/bam2fastq.py
    rm -f ${SFILE%.bam}.r1.fq.gz
    rm -f ${SFILE%.bam}.r2.fq.gz
    md5sum $FILE > ${FILE}.md5 && samtools view -F256 -h ${FILE} | samtools sort -n /dev/stdin -T ${STMP}/sorted > ${SFILE} && python ${bamscript} ${SFILE} && gzip ${SFILE%.bam}.r1.fq && gzip ${SFILE%.bam}.r2.fq && samtools view -H $FILE > ${FILE%.bam}.txt && rm $SFILE $FILE && rm -r $STMP && rm ${FILE%.bam}.bai
else
    bamscript=$(pwd)/bam2fastq_single.py 
    rm -f ${SFILE%.bam}.r1.fq.gz
    md5sum $FILE > ${FILE}.md5 && samtools view -F256 -h ${FILE} | samtools sort -n /dev/stdin -T ${STMP}/sorted > ${SFILE} && python ${bamscript} ${SFILE} && gzip ${SFILE%.bam}.r1.fq && samtools view -H $FILE > ${FILE%.bam}.txt && rm $SFILE $FILE && rm -r $STMP && rm ${FILE%.bam}.bai
fi

