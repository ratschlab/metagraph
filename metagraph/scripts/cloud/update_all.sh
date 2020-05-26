#!/bin/bash

for instance in $(gcloud compute instances list --format='value[separator=","](name)'); do
  gcloud compute ssh "ddanciu@$instance" --zone=us-east1-b --command="sudo apt-get purge -y stackdriver-agent"
done
