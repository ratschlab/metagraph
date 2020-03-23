#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from NCBI, using the data stored in Google Cloud

. common.sh

############ Main Script ###############
# Arguments:
# -  the bucket index (1 to 9) where the SRA resides
#  - the SRA id to be downloaded
#  - the location where the download will be placed

# check the command-line arguments
if [ "$#" -ne 3 ]; then
	    echo_err "Usage: download.sh <sra_bucket_idx> <sra_id> <download_dir>"
	    exit 1
fi

# exit on error
set -e

trap "exit" INT
sra_bucket=$1
sra_number=$2
download_dir=$3

download_dir=${download_dir%/}  # remove trailing slash, if present

output_dir="${download_dir}/${sra_number}"
sra_dir="${output_dir}/sra"
fastq_dir="${output_dir}/fastq"
tmp_dir="${output_dir}/tmp"
kmc_dir="${output_dir}/kmc"

mkdir -p "${output_dir}"
mkdir -p "${sra_dir}"
mkdir -p "${tmp_dir}"
mkdir -p "${fastq_dir}"
mkdir -p "${kmc_dir}"

if ! [ -f "${sra_dir}" ]; then
  execute gsutil -q -u metagraph  cp "gs://sra-pub-run-${sra_bucket}/${sra_number}/*" "${sra_dir}/"
else
  echo "${sra_number} already downloaded"
fi

exit_code=0
for sra_file in $(ls -p "${sra_dir}"); do
  # fasterq-dump is flakey, so trying the dump 3 times before giving up
  if ! (execute_retry fasterq-dump "${sra_dir}/${sra_file}" -f -e 4  -O "${fastq_dir}" -t "${tmp_dir}"); then
    exit_code=1
  fi
  rm -rf "${tmp_dir}/*"
  if (( exit_code != 0 )); then
    echo_err "Download of $sra_number failed"
    exit $exit_code
  fi
done

if [ -z "$(ls -A ${sra_dir})" ]; then
  echo_err "No SRA files available. Good-bye."
  exit 1
fi

kmc_input=${kmc_dir}/sra_file_list
for i in $(ls -p "${fastq_dir}"); do
  echo -e "${fastq_dir}/$i\n" >> "$kmc_input"
done
bin_count=$(( $(ulimit -n) - 10))
kmc_output="${output_dir}/stats"
execute kmc -k31 -ci1 -m2 -fq -cs65535 -n$bin_count -j$kmc_output @${kmc_input} ${kmc_dir}/${sra_number}.kmc $tmp_dir
#TOOD: set ci back to 2 above
rm -rf "${tmp_dir}" "${fastq_dir}"

exit $exit_code
# Note: sra_dir is deleted later in the python client, because we want to measure its size
