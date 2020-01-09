#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client

# displays a message to stderr in red
function echo_err() {
	RED='\033[0;31m'
	NC='\033[0m'
	echo -e "${RED}Error:${NC} $*" 1>&2;
}

# executes the given command
function execute {
    cmd=("$@")
    echo "Executing ${cmd[*]}"

    # execute the command:
    "${cmd[@]}" || exit 1
}

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - location of the dbg file to clean
#  - the location where the cleaned file bill be placed

# check the command-line arguments
if [ "$#" -ne 3 ]; then
	    echo_err "Usage: clean.sh <sra_id> <input_file> <output_dir>"
	    exit 1
fi

sra_number=$1
input_file=$2
output_dir=$3

# for reasons I don't understand glob expansion doesn't work, so doing it manually by concatenating all downloaded fastq.gz files
input_filenames=""
for i in $(ls -p "${output_dir}/${sra_number}/"); do
  input_filenames="$input_filenames ${output_dir}/${sra_number}/$i"
done
execute metagraph clean -v -p 4 --prune-unitigs 0 --fallback 3 --prune-tips 20 --unitigs --to-fasta -o "${output_dir}" "${input_file}"