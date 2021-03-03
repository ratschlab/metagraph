#!/bin/bash

# displays a message to stderr in red
function echo_err() {
    RED='\033[0;31m'
    NC='\033[0m'
    datetime=$(date "+%Y-%m-%d %H:%M:%S:%3N")
    echo -e "${RED}[$datetime: Error]${NC} $@" 1>&2;
}

function echo_log() {
    GRAY='\033[0;37m'
    NC='\033[0m'
    datetime=$(date "+%Y-%m-%d %H:%M:%S:%3N")
    echo -e "${GRAY}[$datetime: Info]${NC} $@" 1>&2;
}

# executes the given command with 3 retries
function execute {
    cmd=("$@")
    echo_log "Executing ${cmd[*]}"
    n=0
    while true; do
        "${cmd[@]}" && break 
        if ((n == 3)); then
            echo_err "Command ${cmd[*]} failed $n times. Bailing out"
            exit 1
        fi
        n=$((n+1))
       echo_log "Command ${cmd[*]} failed. Retrying..."
    done
}

# check that we can find the aspera ssh key
function check_key {
    declare -a KEY_DIRS=("/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh" "$HOME/.ssh/asperaweb_id_dsa.openssh")

    for KEY_DIR in ${KEY_DIRS[@]}; do 
        if [ -f $KEY_DIR ]; then 
            echo_log "Found aspera key in: $KEY_DIR"
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
        echo_err "ascp executable not found. Please install by downloading it from" "https://downloads.asperasoft.com/en/downloads/62"
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
    local file_name=$1
    local chunk_count=$2
    local chunk_index=$3
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
    echo_log "Writing lines $first-$last into $INPUT_CHUNK_FILE"
    sed -n "${first},${last}p" $file_name > $INPUT_CHUNK_FILE
}

# builds and cleans a deBrujin graph using 'metagraph build' and 'metagraph clean'
# Arguments:
#  - the download directory
#  - the name of the downloaded submission (e.g. SRA134234)
function build_and_clean_graph {
    local download_dir=$1
    local sra_number=$2
    if ! [ -d "$download_dir/$sra_number" ]; then
        echo_err "Directory $download_dir/sra_number not found"
        exit 1
    fi
    if ! [ -d "$HOME/dbg" ]; then
        mkdir -p "$HOME/dbg"
    fi
    if [ -f "$HOME/dbg/${sra_number}.dbg" ]; then
        echo_log "deBrujin graph already generated for $sra_number. Skipping."
    else
        # for reasons I don't understand glob expansion doesn't work, so doing it manually by concatenating all downloaded fastq.gz files
        input_filenames="" 
        for i in $(ls -p "${download_dir}/${sra_number}/"); do
          input_filenames="$input_filenames ${download_dir}/${sra_number}/$i"
        done
        execute metagraph build -v -p 2 -k 31 --mem-cap-gb 5 --mode canonical --count-kmers -o "${HOME}/dbg/${sra_number}"  $input_filenames
    fi

    if [ -f "${HOME}/dbg/${sra_number}pruned.fasta.gz" ]; then
        echo_log "deBrujin graph already pruned for $sra_number. Skipping."
    else
        execute metagraph clean -v -p 2 --prune-unitigs 0 --fallback 3 --prune-tips 62 --unitigs --to-fasta -o "${HOME}/dbg/${sra_number}pruned" "${HOME}/dbg/${sra_number}.dbg"
    fi
    # clean the downloaded files
    # rm -rf ${download_dir}/${sra_number}
}

############ Main Script ###############
# Arguments:
#. - the name of the file to read ENA IDs from
#. - the total number of chunks
#  - the index of the current chunk (zero-based)

# make sure subprocesses are killed on exit
trap "exit" INT TERM
trap "kill 0" EXIT

metagraph_path="$HOME/projects2014-metagenome"
if ! [[ $- == *i* ]]; then # not in interactive mode, need to source bashrc manually
    echo_log "Non-interactive shell detected - setting paths manually"
    PATH=$PATH:~/.aspera/connect/bin/
    PATH="$PATH:${metagraph_path}/metagraph/build"
    # and to update git
    cd $metagraph_path
    git checkout dd/cloud
    git pull origin dd/cloud
fi

prereq # make sure all prerequsites are fulfilled
  # check the command-line arguments
if (( "$#" == 3 )); then
  echo_log "Using command line parameters."
  get_chunk $1 $2 $3 # get the $3-th chunk out of $2
else
  if (( "$#" != 0 )); then
          echo_err "Usage: dbg.sh [<srr_ids_file> <chunk_count> <current_chunk>]"
          exit 1
  fi
  echo_log "No parameters specified, assuming Google Cloud Instance. Retrieving params from the cloud."
  sra_id_file="${metagraph_path}/metagraph/scripts/cloud/SraRunInfo_viridiplantae_500.csv.ids"
  curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" > /tmp/instance_id
  curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/num_instances -H "Metadata-Flavor: Google" > /tmp/num_instances
  num_instances=$(cat /tmp/num_instances)
  current_instance=$(cat /tmp/instance_id)
  get_chunk $sra_id_file $num_instances $current_instance
fi
PORT=33001
download_dir="${HOME}/data"
if ! [ -d $download_dir ]; then
    mkdir -p $download_dir
fi
STARTTIME=$(date +%s)
# read from $INPUT_CHUNK_FILE line by line and download the corresponding ENA submission
i=0
while read -r sra_number; do
    echo -e "\n\n\n==== Round number $i ==========="
    if (( ${#sra_number} > 9 )); then
        last_digits=${sra_number:9:${#sra_number}}
        SUBFOLDER="/$( printf "%03d" $((10#$last_digits)) )"
    else
        SUBFOLDER=""
    fi
    FOLDER="/vol1/fastq/${sra_number:0:6}$SUBFOLDER"
    start=$(date +%s)
    execute ascp -QTd -k 1 -P "$PORT" -i "$ASPERA_SSH" era-fasp@fasp.sra.ebi.ac.uk:"${FOLDER}/${sra_number}/" "${download_dir}"
    end=$(date +%s)
    echo_log "Downloaded data in $(($end-$start)) seconds"
    wait # for the previous background graph build process
    build_and_clean_graph $download_dir $sra_number &
    i=$((i+1))
done < $INPUT_CHUNK_FILE
wait # for the last graph build+clean
ENDTIME=$(date +%s)
echo_log "Graph download+generation for all submissions completed in $(($ENDTIME - $STARTTIME)) seconds"
