#!/usr/bin/env bash

set -e

BASEDIR=$(dirname "$0")

VENV_DIR=${BASEDIR}/integration_test_env

echo "Setting up virtual environment"
python3 -m venv ${VENV_DIR}
source ${VENV_DIR}/bin/activate

echo "Installing required packages"
pip install -e ${BASEDIR}/../api/python
pip install parameterized==0.7.1

# signalling, that the environment set up was successful
touch ${VENV_DIR}/DONE
