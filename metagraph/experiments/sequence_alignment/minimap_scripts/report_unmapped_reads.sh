#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <sam file name> <fastq file name>"
    echo "To extract raw reads corresponding to a subset of those reads in the sam file"
    echo "Sam file should either only include mapped reads or only include unmapped reads"
    echo "An example to  filter reads in a sam file"
    echo "samtools view -f 4 short_reads_all/alignment.sam > unmapped_short_reads_all.sam"
    echo "or"
    echo "samtools view -F 4 short_reads_all/alignment.sam > mapped_short_reads_all.sam"
    
    exit
fi

SAM_FILE=$1
FASTQ_FILE_PATH=$2
OUTPUT_FILENAME="raw_${SAM_FILE}"

rm indices_in_original_file.txt; rm indices.txt; rm $OUTPUT_FILENAME;

#awk '{line=$0; split(line,array,"\t"); print(array[1]);}' $SAM_FILE\
#| awk '{line=$0; n=split(line,array,"./"); }' > indices.txt &&
awk '{line=$0; split(line,array,"\t"); print(array[1]);}' $SAM_FILE > indices.txt &&

prev_index=""
while read index; do
    if [ "$prev_index" != "$index" ]; then
#        grep -n "@SRR[0-9]*\.${index} " $FASTQ_FILE_PATH | awk '{line=$0; split(line, array, ":"); print(array[1]);}' >> indices_in_original_file.txt
        grep -n "${index}" $FASTQ_FILE_PATH | awk '{line=$0; split(line, array, ":"); print(array[1]);}' >> indices_in_original_file.txt
    fi
    prev_index=$index
done < indices.txt && 

while read orig_index; do
#        head -n $((${orig_index}+3)) $FASTQ_FILE_PATH | tail -n 4 >> $OUTPUT_FILENAME
        head -n $((${orig_index}+1)) $FASTQ_FILE_PATH | tail -n 2 >> $OUTPUT_FILENAME
done < indices_in_original_file.txt
