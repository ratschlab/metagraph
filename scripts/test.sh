#/usr/bin/env bash

file="$(dirname ${BASH_SOURCE[0]})/../metagraph/tests/data/transcripts_1000.fa"
exe="$(dirname ${BASH_SOURCE[0]})/../metagraph/build/metagengraph"

$exe build -k 12 -o test_graph $file
$exe stats test_graph
$exe annotate -i test_graph -o test_annotation --anno-header $file
$exe stats -a test_annotation
$exe classify -i test_graph -a test_annotation <(cat $file | head -n 1000) | tee test_annotation_out.tsv
[[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
        == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
    && echo Passed: OK
