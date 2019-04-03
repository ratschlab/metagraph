#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <ID>"
  exit 1
fi
token="$(cat ~/gdc-user-token.2018-07-03T15_27_25.227Z.txt)"
curl -sS -O -J -H "X-Auth-Token: $token" "https://api.gdc.cancer.gov/data/${1}"
