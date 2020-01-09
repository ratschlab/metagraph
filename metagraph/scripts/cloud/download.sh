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

# check that we can find the aspera ssh key
function check_key {
	declare -a KEY_DIRS=("/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh" "$HOME/.ssh/asperaweb_id_dsa.openssh")

	for KEY_DIR in "${KEY_DIRS[@]}"; do
		if [ -f "$KEY_DIR" ]; then
			echo "Found aspera key in: $KEY_DIR"
			ASPERA_SSH=$KEY_DIR
			return
		fi
	done
	echo_err "Aspera key not found in any of: " "${KEY_DIRS[@]}" "."
	echo_err "Please copy the key from your installation directory to one of the above locations"
	exit 1
}

############ Main Script ###############
# Arguments:
#  - the SRA id to be downloaded using aspera
#  - the location where the download will be placed

# check the command-line arguments
if [ "$#" -ne 2 ]; then
	    echo_err "Usage: download.sh <sra_id> <download_dir>"
	    exit 1
fi

trap "exit" INT
PORT=33001
sra_number=$1
download_dir=$2
if ! [ -d $download_dir ]; then
	mkdir -p $download_dir
fi
check_key
echo
echo
if (( ${#sra_number} > 9 )); then
  SUBFOLDER="/$(printf "%03d" ${sra_number:9:${#sra_number}})"
else
  SUBFOLDER=""
fi
FOLDER="/vol1/fastq/${sra_number:0:6}$SUBFOLDER"
execute ascp -QT -k1 -d -P "$PORT" -i "$ASPERA_SSH" era-fasp@fasp.sra.ebi.ac.uk:"${FOLDER}/${sra_number}/" "${download_dir}"

