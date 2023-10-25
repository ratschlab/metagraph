#!/usr/bin/env bash

set -e

if [ $# -ne 2 ]; then
    echo -e "Usage:\n$0 <API_DIR> <VENV_DIR>" >&2
    exit 1
fi

API_DIR="$1"
VENV_DIR="$2"

if [ -f ${VENV_DIR}/DONE ]; then
  echo "Found a previously set up virtual environment in ${VENV_DIR}"
  exit 0
fi
echo "Setting up virtual environment in ${VENV_DIR}"
python3 -m venv ${VENV_DIR}
if [ -f /etc/pip.conf ]; then
    mkdir -p ${VENV_DIR}/pip.conf.d/
    cp /etc/pip.conf ${VENV_DIR}/pip.conf.d/
    echo "/etc/pip.conf exists and has been copied to the virtual environment"
else
    echo "/etc/pip.conf doesn't exist. Continuing without copying pip configuration"
fi
source ${VENV_DIR}/bin/activate

echo "Installing required packages"
pip install -e ${API_DIR}
pip install parameterized==0.7.4

# signalling, that the environment set up was successful
touch ${VENV_DIR}/DONE
