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

module load samtools
md5sum $FILE > ${FILE}.md5 && samtools view -F256 $FILE | awk '{ print "@"$1"_"$2"\n"$10"\n+\n"$11 }' | gzip > ${FILE%.bam}.fastq.gz && samtools view -H $FILE > ${FILE%.bam}.txt && rm $FILE

