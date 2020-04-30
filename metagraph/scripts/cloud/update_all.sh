#!/bin/bash

for instance in $(gcloud compute instances list --format='value[separator=","](name)'); do
  gcloud compute ssh "ddanciu@$instance" --zone=us-east1-b --command="cd /home/ddanciu/projects2014-metagenome/metagraph/scripts/cloud; git fetch; git checkout origin/dd/cloud"
done
