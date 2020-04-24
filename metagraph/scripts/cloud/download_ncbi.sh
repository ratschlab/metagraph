#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from NCBI, using the data stored in Google Cloud

. common.sh

# exit with a given exit code; since Python is unable to actually get the exit code, write the exit code into a file
function exit_with() {
  code=$1
  echo $code > "$output_dir/code"
  exit $code
}

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
sra_id=$2
download_dir=$3

download_dir=${download_dir%/}  # remove trailing slash, if present

output_dir="${download_dir}/${sra_id}"
sra_dir="${output_dir}/sra"
fastq_dir="${output_dir}/fastq"
tmp_dir="${output_dir}/tmp"
kmc_dir="${output_dir}/kmc"

mkdir -p "${output_dir}"
mkdir -p "${sra_dir}"
mkdir -p "${tmp_dir}"
mkdir -p "${fastq_dir}"
mkdir -p "${kmc_dir}"

vdb-config --report-cloud-identity no  # the next command should work whether we use the cloud or not
size_str=$(vdb-dump --info "${sra_id}" | grep size | cut -f 2 -d:)
size=${size_str//,}
size="$(echo -e "${size}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

re='^[0-9]+$'
echo "Sra size is $size bytes"
if ! [[ $size =~ $re ]] ; then
   echo_err "Result of vdb-info is not a number"
   exit_with 7
fi
echo $size > "${output_dir}/size"

if (( sra_bucket > 0 )); then  # Get data from GCS
  vdb-config --report-cloud-identity yes  # may not be needed, as we explicitly fetch the file ourselves
  if ! [ -f "${sra_dir}" ]; then
    if ! (execute gsutil -q -u metagraph  cp "gs://sra-pub-run-${sra_bucket}/${sra_id}/*" "${sra_dir}/"); then
      echo_err 'gsutil command failed. Exiting with exit_code 2'
      exit_with 2
    fi
    if [ -z "$(ls -A ${sra_dir})" ]; then
      echo_err "[$sra_id] No SRA files available. Good-bye."
      exit_with 5
    fi
  else
    echo "${sra_id} already downloaded"
  fi

  exit_code=0
  for sra_file in $(ls -p "${sra_dir}"); do
    # fasterq-dump is flakey, so trying the dump 3 times before giving up
    if (( use_gcs )); then
      source="${sra_dir}/${sra_file}"
    else
      source="${sra_id}"
    fi
    if ! (execute_retry fasterq-dump "${source}" -f -e 8  -O "${fastq_dir}" -t "${tmp_dir}"); then
      exit_code=4  # TODO: check if return code 3 is also acceptable as success
    fi
    rm -rf "${tmp_dir}/*"
    if (( exit_code != 0 )); then
      echo_err "[$sra_id] Download failed while running fasterq-dump"
      exit_with $exit_code
    fi
  done
else  # Get data via HTTP (machine must have external ip or internet access via NAT)
    source="${sra_id}"
    vdb-config --report-cloud-identity no  # otherwise NCBI will try to use the cloud and fail
    if ! (execute_retry fasterq-dump "${source}" -f -e 4  -O "${fastq_dir}" -t "${tmp_dir}"); then
      rm -rf "${tmp_dir}/*"
      echo_err "[$sra_id] Download failed while running fasterq-dump"
      exit_with 4
    fi
fi

kmc_input=${kmc_dir}/sra_file_list
for i in $(ls -p "${fastq_dir}"); do
  echo -e "${fastq_dir}/$i\n" >> "$kmc_input"
done
bin_count=$(( $(ulimit -n) - 10))
kmc_output="${output_dir}/stats"
if ! (execute kmc -k31 -ci1 -m2 -fq -cs65535 -t8 -n$bin_count -j"$kmc_output" "@${kmc_input}" "${kmc_dir}/${sra_id}.kmc" "${tmp_dir}"); then
  echo_err "kmc command failed. Exiting with code 6"
  exit_with 6
fi
if ! [ -f $kmc_output ]; then
  echo_err "kmc output '$kmc_output' missing. Exiting with code 6"
  exit_with 6
fi

unique_kmers=$(jq -r ' .Stats | ."#Unique_k-mers"' $kmc_output)
total_kmers=$(jq -r ' .Stats | ."#Total no. of k-mers"' $kmc_output)
coverage=$((total_kmers/unique_kmers))
if ((coverage >= 5)); then
  echo "[$sra_id] Coverage is $coverage, eliminating singletons"
  if ! (execute kmc -k31 -ci2 -m2 -fq -cs65535 -t8 -n$bin_count -j"$kmc_output" "@${kmc_input}" "${kmc_dir}/${sra_id}.kmc" "$tmp_dir"); then
    echo_err "kmc command run #2 failed. Exiting with code 7"
    exit_with 7
  fi
else
  echo "[$sra_id] Coverage is $coverage, keeping singletons"
fi
# edit the KMC output and add the coverage property (this will also eliminate all stuff except for 'Stats')
jq --arg cov $coverage  '.Stats | ."#k-mers_coverage" = $cov' $kmc_output > $tmp_dir/stats
mv $tmp_dir/stats $kmc_output
# just to be on the safe side, set the singleton count to 0 (i.e. "not set") if we don't remove singletons
if ((coverage < 5)); then
  jq '."#k-mers_below_min_threshold" = 0' $kmc_output > $tmp_dir/stats
  mv $tmp_dir/stats $kmc_output
fi

echo singleton_kmers > "${kmc_dir}/${sra_id}.stats"
rm -rf "${tmp_dir}" "${fastq_dir}"
exit_with 0

# Note: sra_dir is deleted later in the python client, because we want to measure its size
