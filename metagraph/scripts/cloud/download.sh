#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from NCBI, using the data stored in Google Cloud

. common.sh

# exit with a given exit code; since Python is unable to actually get the exit code, write the exit code into a file
function exit_with() {
  code=$1
  echo $code > "$output_dir/code"
  exit $code
}

#download the given SRA from ENA
# Vars:
#  sra_id - the sra to download
#  fastq_dir - where to place the downloaded fastq files
function download_ena() {
  len=${#sra_id}
  SUBFOLDER=""
  if (( len > 9 )); then
    SUBFOLDER="/$(printf "%03d" ${sra_id:9:$len})"
  fi
  FOLDER="/vol1/fastq/${sra_id:0:6}$SUBFOLDER"
  execute ascp -QT -k1 -d -P 33001 -i "$ASPERA_SSH" era-fasp@fasp.sra.ebi.ac.uk:"${FOLDER}/${sra_id}/" "${fastq_dir}"
}

# Downloads the given SRA from NCBI. Uses GCS if a bucket is provided, HTTP otherwise
# Vars:
#  bucket - the sra bucket (0 if HTTP)
#  sra_id - the sra to download
#  fastq_dir - where to place the downloaded fastq files
function download_ncbi() {
  sra_dir="${output_dir}/sra"
  mkdir -p "${sra_dir}"

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
    echo "[$sra_id] gsutil finished successfully"
    exit_code=0
    for sra_file in $(ls -p "${sra_dir}"); do
      if ! (execute fasterq-dump "${sra_dir}/${sra_file}" -f -e 4  -O "${fastq_dir}" -t "${tmp_dir}"); then
        echo_err "[$sra_id] fasterq-dump failed, trying fastq-dump. Go get a coffeee."
        if ! (execute fastq-dump "${sra_dir}/${sra_file}" --fasta default  -O "${fastq_dir}" --skip-technical); then
          echo_err "[$sra_id] fastq-dump also failed, this SRA cannot be processed."
          exit_code=4  # TODO: check if return code 3 is also acceptable as success
        fi
      fi
      rm -rf "${tmp_dir}/*"
      if (( exit_code != 0 )); then
        echo_err "[$sra_id] Download failed while running fast(er)q-dump"
        exit_with $exit_code
      fi
    done
  else  # Get data via HTTP (machine must have external ip or internet access via NAT)
      source="${sra_id}"
      vdb-config --report-cloud-identity no  # otherwise NCBI will try to use the cloud and fail
      if ! (execute fasterq-dump "${source}" -f -e 4  -O "${fastq_dir}" -t "${tmp_dir}"); then
        rm -rf "${tmp_dir}/*"
        echo_err "[$sra_id] Download failed while running fasterq-dump"
        exit_with 4
      fi
  fi
  echo "[$sra_id] fasterq-dump command finished successfully"
}

function set_size() {
  size=0
  set +e
  vdb-config --report-cloud-identity no  # the next command should work whether we use the cloud or not
  size_str=$(vdb-dump --info "${sra_id}" | grep size | cut -f 2 -d:)
  size=${size_str//,}
  size="$(echo -e "${size}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
  re='^[0-9]+$'
  echo "Sra size is $size bytes"
  if ! [[ $size =~ $re ]] ; then
     echo_err "Result of vdb-info is not a number"
     size=0
  fi
  set -e
  echo $size > "${output_dir}/size"
}

############ Main Script ###############
# Arguments:
# -  the bucket index (1 to 9) where the SRA resides
#  - the SRA id to be downloaded
#  - the location where the download will be placed

# check the command-line arguments
if [ "$#" -ne 4 ]; then
	    echo_err "Usage: download.sh <source> <sra_id> <download_dir> <sra_bucket_idx>"
	    exit 1
fi

# exit on error
set -e
check_key  # this also sets $ASPERA_SSH

trap "exit" INT
source=$1
sra_id=$2
download_dir=$3
sra_bucket=$4

download_dir=${download_dir%/}  # remove trailing slash, if present
output_dir="${download_dir}/${sra_id}"
fastq_dir="${output_dir}/fastq"
tmp_dir="${output_dir}/tmp"
kmc_dir="${output_dir}/kmc"

mkdir -p "${output_dir}"
mkdir -p "${tmp_dir}"
mkdir -p "${fastq_dir}"
mkdir -p "${kmc_dir}"

if [ "$source" == "ncbi" ]; then
  download_ncbi
else
  download_ena
fi

set_size
# echo "0" > "${output_dir}/size"  # and remove this



kmc_input="${kmc_dir}/sra_file_list"
find "${fastq_dir}" -type f > "$kmc_input"
bytes=$(du -sb "${fastq_dir}" | cut -f1)
MB=$((1000000))
if ((bytes <= 300 * MB)); then
  count=1
elif ((bytes <= 500 * MB)); then
  count=3
elif ((bytes <= 1000 * MB)); then
  count=10
elif ((bytes <= 3000 * MB)); then
  count=20
else
  count=50
fi

echo "Size is $bytes, setting KMC count to $count"


bin_count=$(( $(ulimit -n) - 10))
bin_count=$(( bin_count > 2000 ? 2000 : bin_count))
kmc_output="${output_dir}/stats"
if ! (execute kmc -k23 -ci"${count}" -m2 -fq -cs65535 -t4 -n$bin_count -j"$kmc_output" "@${kmc_input}" "${kmc_dir}/${sra_id}.kmc" "${tmp_dir}"); then
  echo_err "[$sra_id] kmc command failed. Exiting with code 6"
  exit_with 6
fi
if ! [ -f $kmc_output ]; then
  echo_err "[$sra_id] kmc output '$kmc_output' missing. Exiting with code 6"
  exit_with 6
fi
echo "[$sra_id] kmc command finished successfully"

coverage=1
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
