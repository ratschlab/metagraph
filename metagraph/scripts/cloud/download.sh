#!/bin/bash

# Script to download a specific SRA (Sequence Read Archive) directory from ENA (European Nucleutide Archive) using IBM's aspera client
. common.sh
############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - the location where the download will be placed

# check the command-line arguments
if [ "$#" -ne 2 ]; then
	    echo_err "Usage: download.sh <sra_id> <download_dir>"
	    exit 1
fi

# exit on error
set -e

trap "exit" INT
sra_number=$1
download_dir=$2
if ! [ -d "$download_dir" ]; then
	mkdir -p "$download_dir"
fi
check_key  # this also sets $ASPERA_SSH
echo
echo
if (( ${#sra_number} > 9 )); then
  SUBFOLDER="/$(printf "%03d" ${sra_number:9:${#sra_number}})"
else
  SUBFOLDER=""
fi
FOLDER="/vol1/fastq/${sra_number:0:6}$SUBFOLDER"
execute ascp -QT -k1 -d -P 33001 -i "$ASPERA_SSH" era-fasp@fasp.sra.ebi.ac.uk:"${FOLDER}/${sra_number}/" "${download_dir}"

