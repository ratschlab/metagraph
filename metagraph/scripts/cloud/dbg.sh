#!/bin/bash

# displays a message to stderr in red
function echo_err() {
	RED='\033[0;31m'
	NC='\033[0m'
	echo -e "${RED}Error:${NC} $@" 1>&2;
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

	for KEY_DIR in ${KEY_DIRS[@]}; do 
		if [ -f $KEY_DIR ]; then 
			echo "Found aspera key in: $KEY_DIR"
			ASPERA_SSH=$KEY_DIR
			return 
		fi
	done
	echo_err "Aspera key not found in any of: ${KEY_DIRS[@]}."
	echo_err "Please copy the key from your installation directory to one of the above locations"
	exit 1
}

# make sure the aspera client is available
function check_ascp {
	if ! [ -x "$(command -v ascp)" ]; then
		echo_err "ascp executable not found. Please install by downloading it from" "https://download.asperasoft.com/download/sw/connect/3.9.7/ibm-aspera-connect-3.9.7.175481-linux-g2.12-64.tar.gz"
    	exit 1
  	fi
}

function check_metagraph {
	if ! [ -x "$(command -v metagraph)" ]; then
		echo_err "metagraph executable not found in PATH. Would it be too much to ask to add it?"
    	exit 1
  	fi
}


# check the preprequisites for running the script
function prereq {
	check_ascp
	check_key
	check_metagraph
}

# returns the n-th chunk from the given file
# Arguments:
#  - the file to split into chunks
#  - the number of chunks to split into
#  - the index of the desired chunk (zero-based)
# Returns: 
#  - the name of the file where the n-th chunk was placed, in $INPUT_CHUNK_FILE
function get_chunk {
	file_name=$1
	chunk_count=$2
	chunk_index=$3
	if ! [ -f $file_name ]; then
		echo_err "Could not find file: $file_name"
		exit 1
	fi
	total_lines=$(wc -l < $file_name)
	(( lines_per_chunk = (total_lines + chunk_count - 1) / chunk_count ))
	if (( lines_per_chunk > total_lines )); then
		lines_per_chunk=$total_lines
	fi
	if (( lines_per_chunk == 0 )); then
		lines_per_chunk=1
	fi
	(( first=(lines_per_chunk*chunk_index + 1) ))
	(( last=(lines_per_chunk*(chunk_index + 1)) ))
    INPUT_CHUNK_FILE="/tmp/chunk$3_$2"
    echo "Writing lines $first-$last into $INPUT_CHUNK_FILE"
	sed -n "${first},${last}p" $file_name > $INPUT_CHUNK_FILE
}

############ Main Script ###############
# Arguments:
#. - the name of the file to read ENA IDs from
#. - the total number of chunks
#  - the index of the current chunk (zero-based)

# check the command-line arguments
if [ "$#" -ne 3 ]; then
	    echo_err "Usage: dbg.sh <srr_ids_file> <chunk_count> <current_chunk>"
	    exit 1
fi

trap "exit" INT
prereq # make sure all prerequsites are fulfilled
get_chunk $1 $2 $3 # get the $3-th chunk out of $2
PORT=33001
download_dir="${HOME}/data1"
if ! [ -d $download_dir ]; then
	mkdir -p $download_dir
fi
# read from $INPUT_CHUNK_FILE line by line and download the corresponding ENA submission
while read -r sra_number; do
	echo
	echo
	if (( ${#sra_number} > 9 )); then
		SUBFOLDER="/$(printf "%03d" ${sra_number:9:${#sra_number}})"
	else
		SUBFOLDER="/$sra_number"
	fi
	FOLDER="/vol1/fastq/${sra_number:0:6}$SUBFOLDER"
	if [ -d "$download_dir/$sra_number" ]; then
		echo "Submission $sra_number already downloaded. Skipping"
	else 
		ascp -QT -P "$PORT" -i "$ASPERA_SSH" era-fasp@fasp.sra.ebi.ac.uk:"${FOLDER}/${sra_number}/" "${download_dir}"
    fi
	if [ -f "$HOME/dbg/${sra_number}.dbg" ]; then
		echo "deBrujin graph already generated for $sra_number. Skipping."
	else
		metagraph build -v -p 1 -k 10 --canonical --count-kmers -o "${HOME}/dbg/${sra_number}"  "${download_dir}/${sra_number}/*"
	fi
	if [ -f "${HOME}/dbg/${sra_number}pruned.fasta.gz" ]; then
		echo "deBrujin graph already pruned for $sra_number. Skipping."
	else
		execute metagraph clean -v -p 4 --prune-unitigs 0 --fallback 3 --prune-tips 20 --unitigs --to-fasta -o "${HOME}/dbg/${sra_number}pruned" "${HOME}/dbg/${sra_number}.dbg"
	fi
done < $INPUT_CHUNK_FILE
