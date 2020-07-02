#!/bin/bash

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
        echo "Could not find file: $file_name"
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
#. - the total number of chunks
#  - the index of the current chunk (zero-based)

# make sure subprocesses are killed on exit
function start() {
  if (( "$#" == 2 )); then
    echo "Using command line parameters."
    num_instances=$1
    current_instance=$2
  else
    if (( "$#" != 0 )); then
            echo "Usage: dbg.sh [<chunk_count> <current_chunk>]"
            exit 1
    fi
    echo "No parameters specified, assuming Google Cloud Instance. Retrieving params from the cloud."
    curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" > /tmp/instance_id
    curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/num_instances -H "Metadata-Flavor: Google" > /tmp/num_instances
    num_instances=$(cat /tmp/num_instances)
    current_instance=$(cat /tmp/instance_id)
  fi
  PATH=$PATH:/snap/bin:/usr/local/bin
  gsutil cp gs://mg38/files /tmp
  echo "Getting chunk $current_instance out of $num_instances"
  get_chunk /tmp/files $num_instances $current_instance
  cd /tmp
  rm -rf /tmp/sync
  mkdir /tmp/sync
  # copy  data from Google Storage to a local directory
  cat $INPUT_CHUNK_FILE | gsutil -m cp -r -I  /tmp/sync/
  # leomed:/cluster/work/grlab/projects/metagenome/data/cloudcompute
  rsync -avzm --stats --safe-links --ignore-existing --dry-run --human-readable /tmp/sync/ hex:/cluster/project/raetsch/lab/01/projects/metagraph > /tmp/transfer.log
  sed -i 1d /tmp/transfer.log
  head -n -17 /tmp/transfer.log > /tmp/transf.log

  lines=$(( wc -l ))
  echo "Transferring $lines files"
  rm -rf /tmp/transfers*
  split -l 100 /tmp/transf.log /tmp/transfers
  ls /tmp/transfers* | parallel --line-buffer --verbose -j 5 rsync --progress -av -r --files-from {} /tmp/sync/ hex:/cluster/project/raetsch/lab/01/projects/metagraph
}

export -f start
export -f get_chunk
su ddanciu -c "bash -c start"