#!/usr/bin/bash

./metagengraph classify -v -i ~/big_graph/graph_SRR_k20_2/graph_SRR_2.dbg -a ~/big_graph/graph_SRR_k20_2/merged_annotation.color.annodbg --discovery-fraction 0.8 --labels-delimiter ", " ~/transcripts_1000.fa
