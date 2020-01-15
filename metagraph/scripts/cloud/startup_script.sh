#!/bin/bash
# script that gets excecuted when a cloud instance is started

# get the id of the instance from the metadata server and place it into a file
curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" > /tmp/instance_id
curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/server_host -H "Metadata-Flavor: Google" > /tmp/server_host

instance_id=$(cat /tmp/instance_id)
current_instance=$(cat /tmp/server_host)

metagraph_path="$HOME/projects2014-metagenome/metagraph"
cd "$metagraph_path"
git checkout dd/cloud  # TODO: remove this at some point
git pull origin dd/cloud

mkdir -p build
cd build
cmake ..
make -j

cd "$metagraph_path/scripts/cloud"
python3 ./client.py --server_host="${current_instance}" --client_id="${instance_id}"