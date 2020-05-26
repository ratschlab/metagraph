#!/usr/bin/env bash

set -e

if [ $# -ne 2 ]; then
    echo -e "Usage:\n$0 <API_DIR> <VENV_DIR>" >&2
    exit 1
fi

API_DIR="$1"
VENV_DIR="$2"

echo "Setting up virtual environment in ${VENV_DIR}"
python3 -m venv ${VENV_DIR}
source ${VENV_DIR}/bin/activate

echo "Installing required packages"
pip install -e ${API_DIR}
pip install parameterized==0.7.1

# signalling, that the environment set up was successful
touch ${VENV_DIR}/DONE
