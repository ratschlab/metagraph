#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client

. common.sh

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - location of the dbg file to clean
#  - the location where the cleaned file bill be placed

# check the command-line arguments
if [ "$#" -ne 3 ]; then
      cmd=("$@")
	    echo_err "Usage: clean.sh <sra_id> <input_file> <output_file>, called with ${cmd[*]}"
	    exit 1
fi

set -e # exit on error

sra_number=$1
input_file=$2
output_file=$3
# --prune-unitigs 0 --fallback 3
execute metagraph clean -v -p 1 --min-count 2  --prune-unitigs 0 --fallback 5 --prune-tips 62 --to-fasta -o "${output_file}" "${input_file}"
rm -rf $(dirname "${input_file}")
