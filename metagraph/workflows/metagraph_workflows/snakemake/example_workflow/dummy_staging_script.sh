#!/usr/bin/env bash

SAMPLE_ID="$1"
OUTPUT_FILE="$2"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ln -s ${SCRIPT_DIR}/example_data/${SAMPLE_ID}.fa ${OUTPUT_FILE}