#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client

. common.sh

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - the location of the downloaded files
#  - the location where the graph will be placed

# check the command-line arguments
if [ "$#" -ne 6 ]; then
	    echo_err "Usage: build.sh <sra_id> <input_dir> <output_dir> <memory_in_GB> <cores> <container_type>"
	    exit 1
fi

# exit on error
set -e

sra_id=$1
download_dir=$2
output_dir=$3
mem_cap_gb=$4
cores=$5
container_type=$6
input_dir="${download_dir}/kmc"

mkdir -p "${output_dir}"
tmp_dir="${output_dir}/tmp"
mkdir -p "${tmp_dir}"


if [ -z "$(ls -A ${input_dir})" ]; then
  echo_err "[$sra_id] No input files given. Good-bye."
  exit 100
fi

exit_code=0
set +e
if execute metagraph build -v -p "$cores" -k 23 --canonical --state small --count-kmers -o "${output_dir}/${sra_id}" --mem-cap-gb ${mem_cap_gb} --disk-swap ${tmp_dir} --disk-cap-gb 400 --count-width 16 ${input_dir}/${sra_id}.kmc.kmc_pre; then
  exit_code=0
else
  exit_code=1
fi
space=$(df -P "$output_dir" | tail -1 | awk '{print $4}')
if (( space < 2*mem_cap_gb )); then # we ran out of disk
  echo_err "[$sra_id] Not enough disk space: $space bytes left. Giving up"
  exit_code=2
fi
execute rm -rf ${tmp_dir}
mkdir -p "${tmp_dir}"

execute metagraph stats -p "$cores" --count-dummy "${output_dir}/${sra_id}.dbg" | grep "edges" > "${output_dir}/${sra_id}.stats"

execute rm -rf "${input_dir}"
if (( exit_code == 0 )); then
  echo "[$sra_id] Build script finished successfully."
else
  echo "[$sra_id] Build script finished with failure."
fi
exit $exit_code
