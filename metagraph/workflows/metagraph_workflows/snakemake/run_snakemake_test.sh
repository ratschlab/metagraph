#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${SCRIPT_DIR}

CORES=2
snakemake --configfile default.yml test_workflow/test.yml -p --cores ${CORES} "$@"
