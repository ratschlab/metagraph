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

samtools view -F256 $FILE | awk '{ print ">"$1"_"$2"\n"$10}' | gzip > ${FILE%.bam}.fa.gz && samtools view -H $FILE > ${FILE%.bam}.txt && rm $FILE

