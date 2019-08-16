#!/bin/env bash


if [ $# -ne 2 ]; then
    echo -e "Usage:\n$0 <fastq.gz> <trim_param>" >&2
    exit 1
fi

FILE="$1"

if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

trim_param=$2

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

trimmed_filename="${FILE%.fastq.gz}.trimfq_${trim_param}.fastq.gz"

if [ ! -f $trimmed_filename ]; then
    seqtk trimfq -q $trim_param "$FILE" | gzip > $trimmed_filename
    ret_code=$?
    if [ $ret_code -eq 0 ]; then
        rm $FILE
        echo OK
    fi
fi
