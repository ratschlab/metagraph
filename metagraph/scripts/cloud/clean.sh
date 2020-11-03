#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client

. common.sh

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - location of the dbg file to clean
#  - the location where the cleaned file bill be placed

# check the command-line arguments
if [ "$#" -ne 6 ]; then
      cmd=("$@")
	    echo_err "Usage: clean.sh <sra_id> <input_file> <output_file> <num_singletons> <fallback> <cores>, called with ${cmd[*]}"
	    exit 1
fi

set -e # exit on error

sra_number=$1
input_file=$2
output_file=$3
num_singletons=$4
fallback=$5
cores=$6
execute metagraph clean -v -p "$cores" --min-count "${fallback}" --prune-unitigs 0 --fallback "${fallback}" --prune-tips 62 --to-fasta -o "${output_file}" "${input_file}"
rm -rf $(dirname "${input_file}")
