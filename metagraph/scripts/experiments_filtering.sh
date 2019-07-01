#!/bin/env bash


KMC="$(dirname ${BASH_SOURCE[0]})/../build/KMC/kmc"
Metagraph="$(dirname ${BASH_SOURCE[0]})/../build/metagraph"

if [ $# -ne 4 ]; then
    echo -e "Usage:\n$0 <fastq.gz> <trim_param> <abundance> <threshold>" >&2
    exit 1
fi

FILE="$1"

if [ ! -f $FILE ]; then
    echo "File does not exist"
    exit 1
fi

trim_param=$2
cutoff=$3
threshold=$4

if [ ${FILE: -4} == '.bz2' ]; then
    echo "Error: gzip compressed file expected"
    exit 1
fi

filename="$(basename $FILE)"
trimmed_filename="${filename%.fastq.gz}.trimfq_${trim_param}.fastq.gz"

if [ ! -f $trimmed_filename ]; then
    seqtk trimfq -q $trim_param "$FILE" | gzip > $trimmed_filename
fi


tmp_dir="./kmc_${filename}_${trim_param}_${cutoff}_${threshold}"
mkdir -p $tmp_dir

time $KMC -k21 -m5 -fq -t10 -ci${cutoff} $trimmed_filename \
                                         "$tmp_dir/${trimmed_filename}.kmc" \
                                         "$tmp_dir" \
                                         > "${trimmed_filename%.fastq.gz}_${cutoff}_${threshold}.stats"

filtered_filename="${trimmed_filename%.fastq.gz}.filter_k20_s$(($cutoff-1))_${threshold}.fastq.gz"
if [ ! -f $filtered_filename ]; then
    cp $trimmed_filename $tmp_dir
    $Metagraph filter --parallel 6 -k 20 --noise-freq $(($cutoff-1)) --noisy-kmers $threshold --kmc "$tmp_dir/$trimmed_filename"
    $Metagraph filter --generate-fastq -k 20 --noise-freq $(($cutoff-1)) --noisy-kmers $threshold "$tmp_dir/$trimmed_filename"
    mv $tmp_dir/$filtered_filename ./$filtered_filename
    rm $tmp_dir/$trimmed_filename
fi

time $KMC -k21 -m5 -fq -t10 -ci${cutoff} $filtered_filename \
                                         "$tmp_dir/database" \
                                         "$tmp_dir" \
                                         > "${trimmed_filename%.fastq.gz}_${cutoff}_${threshold}_filtered.stats"

rm -r $tmp_dir
