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

trap "exit" INT
sra_bucket=$1
sra_number=$2
download_dir=$3

sra_dir="${download_dir}/sra"

if ! [ -d $download_dir ]; then
	mkdir -p $download_dir
	mkdir -p ${sra_dir}
fi
gsutil -u metagraph  cp -r "gs://sra-pub-run-${bucket}/${sra_number}" "${sra_dir}/"
for sra_file in $(ls -p ${sra_dir}); do
  fasterq-dump ${sra_file} -f -e 10  -O "${download_dir}"
done
rm -rf ${sra_dir}
