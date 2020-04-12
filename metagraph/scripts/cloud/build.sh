#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client

. common.sh

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - the location of the downloaded files
#  - the location where the graph will be placed

# check the command-line arguments
if [ "$#" -ne 5 ]; then
	    echo_err "Usage: build.sh <sra_id> <input_dir> <output_dir> <memory_in_GB> <container_type"
	    exit 1
fi

# exit on error
set -e

sra_number=$1
download_dir=$2
output_dir=$3
mem_cap_gb=$4
container_type=$5
input_dir="${download_dir}/kmc"

mkdir -p "${output_dir}"
tmp_dir="${output_dir}/tmp"
mkdir -p "${tmp_dir}"


if [ -z "$(ls -A ${input_dir})" ]; then
  echo_err "No input files given. Good-bye."
  exit 1
fi

exit_code=0
set +e
for i in {1..1}; do  # we now know how much memory we need, so retrying is not necessary
  if execute metagraph build -v -p 4 -k 31 --container "${container_type}" --canonical --state small --count-kmers --no-shrink -o "${output_dir}/${sra_number}" --mem-cap-gb ${mem_cap_gb} --tmp-dir ${tmp_dir} --disk-cap-gb 200 --count-width 16 ${input_dir}/${sra_number}.kmc.kmc_pre; then
    exit_code=0
    break
  else
    exit_code=1
  fi
  space=$(df -P "$output_dir" | tail -1 | awk '{print $4}')
  if (( space < 2*mem_cap_gb )); then # we ran out of disk
    echo_err "Not enough disk space to process $sra_number. $space bytes left. Giving up"
    exit_code=2
    break
  fi
  if (( mem_cap_gb > 2)); then
    mem_cap_gb=$(( mem_cap_gb - 2 ))
  else
    echo_err "Build process for $sra_number died and memory size can't be reduced. Giving up."
    exit_code=1
    break
  fi
  execute rm -rf ${tmp_dir}
  mkdir -p "${tmp_dir}"
done

execute rm -rf "${input_dir}"
if (( exit_code == 0 )); then
  echo "Build script finished successfully."
else
  echo "Build script finished with failure."
fi
exit $exit_code
