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

sra_dir="${download_dir}/${sra_number}/sra"

if ! [ -d $download_dir ]; then
	mkdir -p $download_dir
fi
mkdir -p "${sra_dir}"
mkdir -p "${download_dir}/${sra_number}"


if ! [ -f "${sra_dir}" ]; then
  execute gsutil -q -u metagraph  cp "gs://sra-pub-run-${sra_bucket}/${sra_number}/*" "${sra_dir}/"
else
  echo "${sra_number} already downloaded"
fi
tmp_dir="${download_dir}/tmp/${sra_number}"
mkdir -p "${tmp_dir}"
for sra_file in $(ls -p "${sra_dir}"); do
  execute fasterq-dump "${sra_dir}/${sra_file}" -f -e 4  -O "${download_dir}/${sra_number}" -t "${tmp_dir}"
  rm -rf "${tmp_dir}"
done
# Note: sra_dir is deleted later in the python client, because we want to measure its size