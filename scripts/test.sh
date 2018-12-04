#/usr/bin/env bash

file="$(dirname ${BASH_SOURCE[0]})/../metagraph/tests/data/transcripts_1000.fa"
file_small="$(dirname ${BASH_SOURCE[0]})/../metagraph/tests/data/transcripts_100.fa"

if [ -f metagengraph ]; then
  exe="./metagengraph"
else
  exe="$(dirname ${BASH_SOURCE[0]})/../metagraph/build/metagengraph"
fi

$exe build -k 12 -o test_graph $file_small
$exe stats test_graph --count-dummy
echo ""
$exe build -k 12 --graph hash -o test_graph $file_small
$exe stats test_graph.orhashdbg
echo ""

$exe annotate -i test_graph.dbg -o test_annotation --anno-header --header-delimiter '|' $file_small
$exe stats -a test_annotation
$exe classify -i test_graph -a test_annotation <(cat $file | tail -n 200) | tee test_annotation_out.tsv
[[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
        == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
    && echo Passed: OK
echo ""

$exe extend -i test_graph.dbg -a test_annotation -o test_graph $file
$exe stats test_graph.dbg --count-dummy
$exe annotate -i test_graph.dbg -a test_graph -o test_annotation --anno-header --header-delimiter '|' $file
$exe stats -a test_annotation
$exe classify -i test_graph -a test_annotation <(cat $file | tail -n 200) | tee test_annotation_out.tsv
[[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
        == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
    && echo Passed: OK
echo ""

$exe annotate -i test_graph.orhashdbg -o test_annotation --anno-header --header-delimiter '|' $file_small
$exe stats -a test_annotation
$exe classify -i test_graph.orhashdbg -a test_annotation <(cat $file | tail -n 200) | tee test_annotation_out.tsv
[[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
        == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
    && echo Passed: OK
echo ""

$exe extend -i test_graph.orhashdbg -a test_annotation -o test_graph $file
$exe stats test_graph.orhashdbg
$exe annotate -i test_graph.orhashdbg -a test_graph -o test_annotation --anno-header --header-delimiter '|' $file
$exe stats -a test_annotation
$exe classify -i test_graph.orhashdbg -a test_annotation <(cat $file | tail -n 200) | tee test_annotation_out.tsv
[[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
        == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
    && echo Passed: OK
echo ""

$exe build -k 12 -o test_graph $file
$exe stats test_graph --count-dummy
$exe annotate -i test_graph -o test_annotation --anno-header --header-delimiter '|' $file
$exe stats -a test_annotation
echo ""

time $exe transform_anno --anno-type column -o test_annotation test_annotation
$exe stats --anno-type column -a test_annotation
echo ""

time $exe transform_anno --anno-type row --sparse -o test_annotation test_annotation
$exe stats --anno-type row -a test_annotation
echo ""

time $exe transform_anno --anno-type brwt -o test_annotation test_annotation
$exe stats --anno-type brwt -a test_annotation
echo ""

time $exe transform_anno --anno-type bin_rel_wt_sdsl -o test_annotation test_annotation
$exe stats --anno-type bin_rel_wt_sdsl -a test_annotation
echo ""

time $exe transform_anno --anno-type bin_rel_wt -o test_annotation test_annotation
$exe stats --anno-type bin_rel_wt -a test_annotation
echo ""

time $exe transform_anno --anno-type flat -o test_annotation test_annotation
$exe stats --anno-type flat -a test_annotation
echo ""

time $exe transform_anno --anno-type rbfish -o test_annotation test_annotation
$exe stats --anno-type rbfish -a test_annotation
echo ""

