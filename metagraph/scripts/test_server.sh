#/usr/bin/env bash

file_small="$(dirname ${BASH_SOURCE[0]})/../tests/data/transcripts_100.fa"

if [ -f metagengraph ]; then
  exe="./metagengraph"
else
  exe="$(dirname ${BASH_SOURCE[0]})/../build/metagengraph"
fi

$exe build -k 12 -o test_server_graph $file_small
$exe annotate -i test_server_graph.dbg -o test_server_annotation --anno-header --header-delimiter '|' $file_small


$exe server_classify -i test_server_graph -a test_server_annotation.column.annodbg --port 5555 \
      > server.log &
SERVER_PID=$!

sleep 1
echo -n "GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTC" | nc -w 1 "127.0.0.1" 5555 \
      | tee client.log

kill $SERVER_PID
