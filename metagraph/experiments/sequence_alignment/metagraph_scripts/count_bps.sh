#!/bin/bash

FILE_NAME=$1

awk '{ $0; getline line; print length(line) }' $FILE_NAME | awk '{sum += $0} END {print sum}'
