#!/usr/bin/bash

if [ "$#" -ne 1 ] || ! [ -f "$1" ]; then
  echo -e "Usage:\n$0 <logfile>" >&2
  exit 1
fi

logfile=$1

num_files_processed=$(cat $logfile | tr -d '\000' | grep -o "Finished extracting sequences .*" | cut -d' ' -f6 | sort | uniq | wc -l)
echo Files finished: $num_files_processed
size_gb=$(for file in $(cat $logfile | tr -d '\000' | grep -o "Finished extracting sequences .*" | cut -d' ' -f6 | sort | uniq); do ls -l $file; done | cut -d' ' -f 5 | awk '{ sum += $1/2**30; } END { print sum; }')
echo $size_gb Gb

time_sec=$(cat $logfile | tr -d '\000' | grep -o "Finished extracting sequences .*" | tail -1 | cut -d' ' -f8 | tr -d 'sec[,]')
echo $time_sec sec
echo "$(echo "scale=2; $time_sec/$size_gb" | bc) sec/Gb"
