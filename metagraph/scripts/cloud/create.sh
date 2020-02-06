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
#  - the location of the downloaded files
#  - the location where the graph will be placed

# check the command-line arguments
if [ "$#" -ne 3 ]; then
	    echo_err "Usage: build.sh <sra_id> <input_dir> <output_dir>"
	    exit 1
fi

# exit on error
set -e

sra_number=$1
input_dir=$2
output_dir=$3

# for reasons I don't understand glob expansion doesn't work, so doing it manually by concatenating all downloaded fastq.gz files
input_filenames=""
for i in $(ls -p "${input_dir}"); do
  input_filenames="$input_filenames ${input_dir}/$i"
done

if [[ -z "${input_filenames// }" ]]; then
  echo_err "No input files given. Good-bye."
  exit 1
fi
mkdir -p "${output_dir}"
tmp_dir="${output_dir}/tmp"
mkdir -p "${tmp_dir}"
execute metagraph build -v -p 4 -k 31 --container vector_disk --canonical --count-kmers -o "${output_dir}/${sra_number}" --mem-cap-gb 3 --tmp_dir ${tmp_dir}  ${input_filenames}
execute rm -rf "${input_dir}"
execute rm -rf ${tmp_dir}
echo "Build script finished."