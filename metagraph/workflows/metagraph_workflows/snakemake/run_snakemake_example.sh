#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${SCRIPT_DIR}

CORES=2
snakemake --configfile default.yml example_workflow/example.yml -p --cores ${CORES} "$@"
