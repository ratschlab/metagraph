#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${SCRIPT_DIR}

export PYTHONPATH=$(realpath ${SCRIPT_DIR}/../):${PYTHONPATH}

CORES=2
snakemake --configfile default.yml example.yml -p --cores ${CORES} "$@"
