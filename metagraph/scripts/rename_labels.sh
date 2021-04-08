#!/bin/bash

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <METAGRAPH> <ANNOTATION>"
    exit 1
fi

METAGRAPH="$1"
anno_in="$2"

mem=$(ls -s --block-size=1048576 $anno_in | cut -d' ' -f1)
threads=1
pmem=$(($mem / $threads * 5 / 4 + 1000))

extension="${anno_in##*.}"
anno_in_prefix="${anno_in%.*}"
anno_extension="${anno_in_prefix##*.}"
anno_out=${anno_in%.*.*}.relabeled.${anno_extension}.${extension}

### generate label map
$METAGRAPH stats -a ${anno_in} --print-col-names | grep -v info | grep -v Number > ${anno_in}.labels
rev ${anno_in}.labels | cut -f 1 -d '/' | rev | cut -f 1 -d '.' > ${anno_in}.labels_new
paste ${anno_in}.labels ${anno_in}.labels_new > ${anno_in}.labels_map
rm ${anno_in}.labels ${anno_in}.labels_new

echo "5 first renaming rules:"
head -n 5 ${anno_in}.labels_map
echo "Ready to rename the annotation labels for $anno_in"
echo "using the generated mapping ${anno_in}.labels_map"
echo "The result will be written into ${anno_out}"
read -p "Press Y/y to continue..." -n 1 -r
echo # move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi

### rename annotation labels
bsub -J relabel_columns_$(basename $anno_in_prefix) \
     -We 24:00 \
     -n $threads \
     -R "rusage[mem=${pmem}] span[hosts=1]" \
     -o /dev/null \
    "$METAGRAPH transform_anno -v \
        --rename-cols ${anno_in}.labels_map \
        -o ${anno_out} \
        -p $threads \
        $anno_in"

