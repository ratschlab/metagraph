#!/bin/bash

# get the id of the instance from the metadata server and place it into a file
curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/instance_id -H "Metadata-Flavor: Google" > /tmp/instance_id