#!/usr/bin/bash

./metagraph query -v -i ~/big_graph/graph_SRR_k20_2/graph_SRR_2.dbg -a ~/big_graph/graph_SRR_k20_2/merged_annotation.column.annodbg --min-kmers-fraction-label 0.8 --labels-delimiter ", " ~/transcripts_1000.fa
