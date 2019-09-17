  $ #!/usr/bin/env bash
  $ file="$TESTDIR/../tests/data/transcripts_1000.fa"
  $ file_small="$TESTDIR/../tests/data/transcripts_100.fa"
  $ exe="$TESTDIR/../build/metagraph"
  $ 
  $ $exe build -k 12 -o test_graph $file_small
  Graph chunk with 45271 k-mers was built in *sec (glob)

  $ $exe stats test_graph.dbg --count-dummy
  Statistics for graph test_graph.dbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 44919
  canonical mode: no
  ========================================================
  ====================== BOSS STATS ======================
  k: 12
  nodes (k-1): 44425
  edges ( k ): 45271
  state: fast
  W stats: {'$': 35, 'A': 11695, 'C': 11698, 'G': 9911, 'T': 11932}
  F stats: {'$': 5, 'A': 11701, 'C': 11677, 'G': 9975, 'T': 11913}
  dummy source edges: 318
  dummy sink edges: 34
  ========================================================

  $ $exe annotate -i test_graph.dbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 44919
  density: 2.531743e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph.dbg -a test_annotation.column.annodbg <(cat $file_small | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000636676.1|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000489043.1|AL669831.3-213|AL669831.3|183|transcribed_processed_pseudogene|\ttranscribed_processed_pseudogene:183:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000636676.1:OTTHUMT00000489043.1:AL669831.3-213 (esc)
  1\tENST00000447954.2|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000006878.3|AL669831.3-214|AL669831.3|355|processed_transcript|\tprocessed_transcript:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000447954.2:OTTHUMT00000006878.3:AL669831.3-214:355 (esc)
  2\tENST00000423796.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006707.1|AC114498.1-201|AC114498.1|607|lincRNA|\tlincRNA:607:ENST00000423796.1:ENSG00000235146.2:OTTHUMG00000002329.1:OTTHUMT00000006707.1:AC114498.1-201:AC114498.1 (esc)
  3\tENST00000450696.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006706.1|AC114498.1-202|AC114498.1|402|lincRNA|\tlincRNA:ENSG00000235146.2:OTTHUMG00000002329.1:AC114498.1:ENST00000450696.1:OTTHUMT00000006706.1:AC114498.1-202:402 (esc)
  4\tENST00000416931.1|ENSG00000225972.1|OTTHUMG00000002338.1|OTTHUMT00000006720.1|MTND1P23-201|MTND1P23|372|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000416931.1:ENSG00000225972.1:OTTHUMG00000002338.1:OTTHUMT00000006720.1:MTND1P23-201:MTND1P23:372 (esc)
  5\tENST00000457540.1|ENSG00000225630.1|OTTHUMG00000002336.1|OTTHUMT00000006718.1|MTND2P28-201|MTND2P28|1044|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000457540.1:ENSG00000225630.1:OTTHUMG00000002336.1:OTTHUMT00000006718.1:MTND2P28-201:MTND2P28:1044 (esc)
  6\tENST00000414273.1|ENSG00000237973.1|OTTHUMG00000002333.2|OTTHUMT00000006715.2|MTCO1P12-201|MTCO1P12|1543|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543 (esc)
  7\tENST00000621981.1|ENSG00000278791.1|-|-|MIR6723-201|MIR6723|89|miRNA|\t-:miRNA:ENST00000621981.1:ENSG00000278791.1:MIR6723-201:MIR6723:89 (esc)
  8\tENST00000427426.1|ENSG00000229344.1|OTTHUMG00000002334.1|OTTHUMT00000006716.1|MTCO2P12-201|MTCO2P12|682|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000427426.1:ENSG00000229344.1:OTTHUMG00000002334.1:OTTHUMT00000006716.1:MTCO2P12-201:MTCO2P12:682 (esc)
  9\tENST00000467115.1|ENSG00000240409.1|OTTHUMG00000002473.1|OTTHUMT00000007027.1|MTATP8P1-201|MTATP8P1|207|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000467115.1:ENSG00000240409.1:OTTHUMG00000002473.1:OTTHUMT00000007027.1:MTATP8P1-201:MTATP8P1:207 (esc)
  10\tENST00000514057.1|ENSG00000248527.1|OTTHUMG00000002335.2|OTTHUMT00000006717.2|MTATP6P1-201|MTATP6P1|681|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000514057.1:ENSG00000248527.1:OTTHUMG00000002335.2:OTTHUMT00000006717.2:MTATP6P1-201:MTATP6P1:681 (esc)
  11\tENST00000416718.2|ENSG00000198744.5|OTTHUMG00000002337.2|OTTHUMT00000006719.2|MTCO3P12-201|MTCO3P12|547|unprocessed_pseudogene|\tunprocessed_pseudogene:547:ENST00000416718.2:ENSG00000198744.5:OTTHUMG00000002337.2:OTTHUMT00000006719.2:MTCO3P12-201:MTCO3P12 (esc)
  12\tENST00000438434.2|ENSG00000268663.1|OTTHUMG00000002340.3|OTTHUMT00000006722.3|WBP1LP6-201|WBP1LP6|424|processed_pseudogene|\tprocessed_pseudogene:ENST00000438434.2:ENSG00000268663.1:OTTHUMG00000002340.3:OTTHUMT00000006722.3:WBP1LP6-201:WBP1LP6:424 (esc)
  13\tENST00000332831.4|ENSG00000284662.1|OTTHUMG00000002581.3|OTTHUMT00000007334.3|OR4F16-201|OR4F16|995|protein_coding|\tprotein_coding:ENST00000426406.3:ENSG00000284733.1:OTTHUMG00000002860.3:OTTHUMT00000007999.3:OR4F29-201:OR4F29:995:ENST00000332831.4:ENSG00000284662.1:OTTHUMG00000002581.3:OTTHUMT00000007334.3:OR4F16-201:OR4F16 (esc)
  14\tENST00000440782.3|ENSG00000229376.3|OTTHUMG00000057431.3|OTTHUMT00000127611.3|CICP3-201|CICP3|2455|processed_pseudogene|\tprocessed_pseudogene:ENST00000440782.3:ENSG00000229376.3:OTTHUMG00000057431.3:OTTHUMT00000127611.3:CICP3-201:CICP3:2455 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe build -k 12 --canonical -o test_graph_canonical $file_small
  Graph chunk with 87426 k-mers was built in *sec (glob)
  $ $exe stats test_graph_canonical.dbg --count-dummy
  Statistics for graph test_graph_canonical.dbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 86840
  canonical mode: yes
  ========================================================
  ====================== BOSS STATS ======================
  k: 12
  nodes (k-1): 85060
  edges ( k ): 87426
  state: fast
  W stats: {'$': 61, 'A': 22854, 'C': 20874, 'G': 20720, 'T': 22917}
  F stats: {'$': 5, 'A': 22861, 'C': 20847, 'G': 20818, 'T': 22895}
  dummy source edges: 526
  dummy sink edges: 60
  ========================================================
  $ $exe annotate -i test_graph_canonical.dbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 86840
  density: 1.307378e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph_canonical.dbg -a test_annotation.column.annodbg <(cat $file_small | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000636676.1|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000489043.1|AL669831.3-213|AL669831.3|183|transcribed_processed_pseudogene|\ttranscribed_processed_pseudogene:183:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000636676.1:OTTHUMT00000489043.1:AL669831.3-213 (esc)
  1\tENST00000447954.2|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000006878.3|AL669831.3-214|AL669831.3|355|processed_transcript|\tprocessed_transcript:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000447954.2:OTTHUMT00000006878.3:AL669831.3-214:355 (esc)
  2\tENST00000423796.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006707.1|AC114498.1-201|AC114498.1|607|lincRNA|\tlincRNA:607:ENST00000423796.1:ENSG00000235146.2:OTTHUMG00000002329.1:OTTHUMT00000006707.1:AC114498.1-201:AC114498.1 (esc)
  3\tENST00000450696.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006706.1|AC114498.1-202|AC114498.1|402|lincRNA|\tlincRNA:ENSG00000235146.2:OTTHUMG00000002329.1:AC114498.1:ENST00000450696.1:OTTHUMT00000006706.1:AC114498.1-202:402 (esc)
  4\tENST00000416931.1|ENSG00000225972.1|OTTHUMG00000002338.1|OTTHUMT00000006720.1|MTND1P23-201|MTND1P23|372|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000416931.1:ENSG00000225972.1:OTTHUMG00000002338.1:OTTHUMT00000006720.1:MTND1P23-201:MTND1P23:372 (esc)
  5\tENST00000457540.1|ENSG00000225630.1|OTTHUMG00000002336.1|OTTHUMT00000006718.1|MTND2P28-201|MTND2P28|1044|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000457540.1:ENSG00000225630.1:OTTHUMG00000002336.1:OTTHUMT00000006718.1:MTND2P28-201:MTND2P28:1044 (esc)
  6\tENST00000414273.1|ENSG00000237973.1|OTTHUMG00000002333.2|OTTHUMT00000006715.2|MTCO1P12-201|MTCO1P12|1543|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543 (esc)
  7\tENST00000621981.1|ENSG00000278791.1|-|-|MIR6723-201|MIR6723|89|miRNA|\tunprocessed_pseudogene:-:miRNA:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543:ENST00000621981.1:ENSG00000278791.1:MIR6723-201:MIR6723:89 (esc)
  8\tENST00000427426.1|ENSG00000229344.1|OTTHUMG00000002334.1|OTTHUMT00000006716.1|MTCO2P12-201|MTCO2P12|682|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000427426.1:ENSG00000229344.1:OTTHUMG00000002334.1:OTTHUMT00000006716.1:MTCO2P12-201:MTCO2P12:682 (esc)
  9\tENST00000467115.1|ENSG00000240409.1|OTTHUMG00000002473.1|OTTHUMT00000007027.1|MTATP8P1-201|MTATP8P1|207|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000467115.1:ENSG00000240409.1:OTTHUMG00000002473.1:OTTHUMT00000007027.1:MTATP8P1-201:MTATP8P1:207 (esc)
  10\tENST00000514057.1|ENSG00000248527.1|OTTHUMG00000002335.2|OTTHUMT00000006717.2|MTATP6P1-201|MTATP6P1|681|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000514057.1:ENSG00000248527.1:OTTHUMG00000002335.2:OTTHUMT00000006717.2:MTATP6P1-201:MTATP6P1:681 (esc)
  11\tENST00000416718.2|ENSG00000198744.5|OTTHUMG00000002337.2|OTTHUMT00000006719.2|MTCO3P12-201|MTCO3P12|547|unprocessed_pseudogene|\tunprocessed_pseudogene:547:ENST00000416718.2:ENSG00000198744.5:OTTHUMG00000002337.2:OTTHUMT00000006719.2:MTCO3P12-201:MTCO3P12 (esc)
  12\tENST00000438434.2|ENSG00000268663.1|OTTHUMG00000002340.3|OTTHUMT00000006722.3|WBP1LP6-201|WBP1LP6|424|processed_pseudogene|\tprocessed_pseudogene:ENST00000438434.2:ENSG00000268663.1:OTTHUMG00000002340.3:OTTHUMT00000006722.3:WBP1LP6-201:WBP1LP6:424 (esc)
  13\tENST00000332831.4|ENSG00000284662.1|OTTHUMG00000002581.3|OTTHUMT00000007334.3|OR4F16-201|OR4F16|995|protein_coding|\tprotein_coding:ENST00000426406.3:ENSG00000284733.1:OTTHUMG00000002860.3:OTTHUMT00000007999.3:OR4F29-201:OR4F29:995:ENST00000332831.4:ENSG00000284662.1:OTTHUMG00000002581.3:OTTHUMT00000007334.3:OR4F16-201:OR4F16 (esc)
  14\tENST00000440782.3|ENSG00000229376.3|OTTHUMG00000057431.3|OTTHUMT00000127611.3|CICP3-201|CICP3|2455|processed_pseudogene|\tprocessed_pseudogene:ENST00000440782.3:ENSG00000229376.3:OTTHUMG00000057431.3:OTTHUMT00000127611.3:CICP3-201:CICP3:2455 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe build -k 12 --graph hash -o test_graph $file_small
  $ $exe stats test_graph.orhashdbg
  Statistics for graph test_graph.orhashdbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 44919
  canonical mode: no
  ========================================================

  $ $exe annotate -i test_graph.orhashdbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 44919
  density: 2.531743e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph.orhashdbg -a test_annotation.column.annodbg <(cat $file_small | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000636676.1|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000489043.1|AL669831.3-213|AL669831.3|183|transcribed_processed_pseudogene|\ttranscribed_processed_pseudogene:183:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000636676.1:OTTHUMT00000489043.1:AL669831.3-213 (esc)
  1\tENST00000447954.2|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000006878.3|AL669831.3-214|AL669831.3|355|processed_transcript|\tprocessed_transcript:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000447954.2:OTTHUMT00000006878.3:AL669831.3-214:355 (esc)
  2\tENST00000423796.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006707.1|AC114498.1-201|AC114498.1|607|lincRNA|\tlincRNA:607:ENST00000423796.1:ENSG00000235146.2:OTTHUMG00000002329.1:OTTHUMT00000006707.1:AC114498.1-201:AC114498.1 (esc)
  3\tENST00000450696.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006706.1|AC114498.1-202|AC114498.1|402|lincRNA|\tlincRNA:ENSG00000235146.2:OTTHUMG00000002329.1:AC114498.1:ENST00000450696.1:OTTHUMT00000006706.1:AC114498.1-202:402 (esc)
  4\tENST00000416931.1|ENSG00000225972.1|OTTHUMG00000002338.1|OTTHUMT00000006720.1|MTND1P23-201|MTND1P23|372|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000416931.1:ENSG00000225972.1:OTTHUMG00000002338.1:OTTHUMT00000006720.1:MTND1P23-201:MTND1P23:372 (esc)
  5\tENST00000457540.1|ENSG00000225630.1|OTTHUMG00000002336.1|OTTHUMT00000006718.1|MTND2P28-201|MTND2P28|1044|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000457540.1:ENSG00000225630.1:OTTHUMG00000002336.1:OTTHUMT00000006718.1:MTND2P28-201:MTND2P28:1044 (esc)
  6\tENST00000414273.1|ENSG00000237973.1|OTTHUMG00000002333.2|OTTHUMT00000006715.2|MTCO1P12-201|MTCO1P12|1543|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543 (esc)
  7\tENST00000621981.1|ENSG00000278791.1|-|-|MIR6723-201|MIR6723|89|miRNA|\t-:miRNA:ENST00000621981.1:ENSG00000278791.1:MIR6723-201:MIR6723:89 (esc)
  8\tENST00000427426.1|ENSG00000229344.1|OTTHUMG00000002334.1|OTTHUMT00000006716.1|MTCO2P12-201|MTCO2P12|682|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000427426.1:ENSG00000229344.1:OTTHUMG00000002334.1:OTTHUMT00000006716.1:MTCO2P12-201:MTCO2P12:682 (esc)
  9\tENST00000467115.1|ENSG00000240409.1|OTTHUMG00000002473.1|OTTHUMT00000007027.1|MTATP8P1-201|MTATP8P1|207|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000467115.1:ENSG00000240409.1:OTTHUMG00000002473.1:OTTHUMT00000007027.1:MTATP8P1-201:MTATP8P1:207 (esc)
  10\tENST00000514057.1|ENSG00000248527.1|OTTHUMG00000002335.2|OTTHUMT00000006717.2|MTATP6P1-201|MTATP6P1|681|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000514057.1:ENSG00000248527.1:OTTHUMG00000002335.2:OTTHUMT00000006717.2:MTATP6P1-201:MTATP6P1:681 (esc)
  11\tENST00000416718.2|ENSG00000198744.5|OTTHUMG00000002337.2|OTTHUMT00000006719.2|MTCO3P12-201|MTCO3P12|547|unprocessed_pseudogene|\tunprocessed_pseudogene:547:ENST00000416718.2:ENSG00000198744.5:OTTHUMG00000002337.2:OTTHUMT00000006719.2:MTCO3P12-201:MTCO3P12 (esc)
  12\tENST00000438434.2|ENSG00000268663.1|OTTHUMG00000002340.3|OTTHUMT00000006722.3|WBP1LP6-201|WBP1LP6|424|processed_pseudogene|\tprocessed_pseudogene:ENST00000438434.2:ENSG00000268663.1:OTTHUMG00000002340.3:OTTHUMT00000006722.3:WBP1LP6-201:WBP1LP6:424 (esc)
  13\tENST00000332831.4|ENSG00000284662.1|OTTHUMG00000002581.3|OTTHUMT00000007334.3|OR4F16-201|OR4F16|995|protein_coding|\tprotein_coding:ENST00000426406.3:ENSG00000284733.1:OTTHUMG00000002860.3:OTTHUMT00000007999.3:OR4F29-201:OR4F29:995:ENST00000332831.4:ENSG00000284662.1:OTTHUMG00000002581.3:OTTHUMT00000007334.3:OR4F16-201:OR4F16 (esc)
  14\tENST00000440782.3|ENSG00000229376.3|OTTHUMG00000057431.3|OTTHUMT00000127611.3|CICP3-201|CICP3|2455|processed_pseudogene|\tprocessed_pseudogene:ENST00000440782.3:ENSG00000229376.3:OTTHUMG00000057431.3:OTTHUMT00000127611.3:CICP3-201:CICP3:2455 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe build -k 12 --canonical --graph hash -o test_graph_canonical $file_small
  $ $exe stats test_graph_canonical.orhashdbg
  Statistics for graph test_graph_canonical.orhashdbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 86840
  canonical mode: yes
  ========================================================

  $ $exe annotate -i test_graph_canonical.orhashdbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 86840
  density: 1.307378e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph_canonical.orhashdbg -a test_annotation.column.annodbg <(cat $file_small | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000636676.1|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000489043.1|AL669831.3-213|AL669831.3|183|transcribed_processed_pseudogene|\ttranscribed_processed_pseudogene:183:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000636676.1:OTTHUMT00000489043.1:AL669831.3-213 (esc)
  1\tENST00000447954.2|ENSG00000230021.9|OTTHUMG00000191652.2|OTTHUMT00000006878.3|AL669831.3-214|AL669831.3|355|processed_transcript|\tprocessed_transcript:ENSG00000230021.9:OTTHUMG00000191652.2:AL669831.3:ENST00000447954.2:OTTHUMT00000006878.3:AL669831.3-214:355 (esc)
  2\tENST00000423796.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006707.1|AC114498.1-201|AC114498.1|607|lincRNA|\tlincRNA:607:ENST00000423796.1:ENSG00000235146.2:OTTHUMG00000002329.1:OTTHUMT00000006707.1:AC114498.1-201:AC114498.1 (esc)
  3\tENST00000450696.1|ENSG00000235146.2|OTTHUMG00000002329.1|OTTHUMT00000006706.1|AC114498.1-202|AC114498.1|402|lincRNA|\tlincRNA:ENSG00000235146.2:OTTHUMG00000002329.1:AC114498.1:ENST00000450696.1:OTTHUMT00000006706.1:AC114498.1-202:402 (esc)
  4\tENST00000416931.1|ENSG00000225972.1|OTTHUMG00000002338.1|OTTHUMT00000006720.1|MTND1P23-201|MTND1P23|372|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000416931.1:ENSG00000225972.1:OTTHUMG00000002338.1:OTTHUMT00000006720.1:MTND1P23-201:MTND1P23:372 (esc)
  5\tENST00000457540.1|ENSG00000225630.1|OTTHUMG00000002336.1|OTTHUMT00000006718.1|MTND2P28-201|MTND2P28|1044|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000457540.1:ENSG00000225630.1:OTTHUMG00000002336.1:OTTHUMT00000006718.1:MTND2P28-201:MTND2P28:1044 (esc)
  6\tENST00000414273.1|ENSG00000237973.1|OTTHUMG00000002333.2|OTTHUMT00000006715.2|MTCO1P12-201|MTCO1P12|1543|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543 (esc)
  7\tENST00000621981.1|ENSG00000278791.1|-|-|MIR6723-201|MIR6723|89|miRNA|\tunprocessed_pseudogene:-:miRNA:ENST00000414273.1:ENSG00000237973.1:OTTHUMG00000002333.2:OTTHUMT00000006715.2:MTCO1P12-201:MTCO1P12:1543:ENST00000621981.1:ENSG00000278791.1:MIR6723-201:MIR6723:89 (esc)
  8\tENST00000427426.1|ENSG00000229344.1|OTTHUMG00000002334.1|OTTHUMT00000006716.1|MTCO2P12-201|MTCO2P12|682|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000427426.1:ENSG00000229344.1:OTTHUMG00000002334.1:OTTHUMT00000006716.1:MTCO2P12-201:MTCO2P12:682 (esc)
  9\tENST00000467115.1|ENSG00000240409.1|OTTHUMG00000002473.1|OTTHUMT00000007027.1|MTATP8P1-201|MTATP8P1|207|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000467115.1:ENSG00000240409.1:OTTHUMG00000002473.1:OTTHUMT00000007027.1:MTATP8P1-201:MTATP8P1:207 (esc)
  10\tENST00000514057.1|ENSG00000248527.1|OTTHUMG00000002335.2|OTTHUMT00000006717.2|MTATP6P1-201|MTATP6P1|681|unprocessed_pseudogene|\tunprocessed_pseudogene:ENST00000514057.1:ENSG00000248527.1:OTTHUMG00000002335.2:OTTHUMT00000006717.2:MTATP6P1-201:MTATP6P1:681 (esc)
  11\tENST00000416718.2|ENSG00000198744.5|OTTHUMG00000002337.2|OTTHUMT00000006719.2|MTCO3P12-201|MTCO3P12|547|unprocessed_pseudogene|\tunprocessed_pseudogene:547:ENST00000416718.2:ENSG00000198744.5:OTTHUMG00000002337.2:OTTHUMT00000006719.2:MTCO3P12-201:MTCO3P12 (esc)
  12\tENST00000438434.2|ENSG00000268663.1|OTTHUMG00000002340.3|OTTHUMT00000006722.3|WBP1LP6-201|WBP1LP6|424|processed_pseudogene|\tprocessed_pseudogene:ENST00000438434.2:ENSG00000268663.1:OTTHUMG00000002340.3:OTTHUMT00000006722.3:WBP1LP6-201:WBP1LP6:424 (esc)
  13\tENST00000332831.4|ENSG00000284662.1|OTTHUMG00000002581.3|OTTHUMT00000007334.3|OR4F16-201|OR4F16|995|protein_coding|\tprotein_coding:ENST00000426406.3:ENSG00000284733.1:OTTHUMG00000002860.3:OTTHUMT00000007999.3:OR4F29-201:OR4F29:995:ENST00000332831.4:ENSG00000284662.1:OTTHUMG00000002581.3:OTTHUMT00000007334.3:OR4F16-201:OR4F16 (esc)
  14\tENST00000440782.3|ENSG00000229376.3|OTTHUMG00000057431.3|OTTHUMT00000127611.3|CICP3-201|CICP3|2455|processed_pseudogene|\tprocessed_pseudogene:ENST00000440782.3:ENSG00000229376.3:OTTHUMG00000057431.3:OTTHUMT00000127611.3:CICP3-201:CICP3:2455 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe annotate -i test_graph.dbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 44919
  density: 2.531743e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph.dbg -a test_annotation.column.annodbg <(cat $file | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000491020.1|ENSG00000116254.17|OTTHUMG00000000952.6|OTTHUMT00000002828.2|CHD5-205|CHD5|646|nonsense_mediated_decay|\t (esc)
  1\tENST00000484532.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003213.2|RPL22-207|RPL22|449|protein_coding|\t (esc)
  2\tENST00000234875.8|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002830.1|RPL22-201|RPL22|2078|protein_coding|\t (esc)
  3\tENST00000480661.1|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002832.2|RPL22-206|RPL22|2279|retained_intron|\t (esc)
  4\tENST00000497965.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002831.2|RPL22-208|RPL22|786|protein_coding|\t (esc)
  5\tENST00000471204.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003678.3|RPL22-205|RPL22|524|protein_coding|\t (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe extend -i test_graph.dbg -a test_annotation.column.annodbg -o test_graph $file
  $ $exe stats test_graph.dbg --count-dummy
  Statistics for graph test_graph.dbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 540664
  canonical mode: no
  ========================================================
  ====================== BOSS STATS ======================
  k: 12
  nodes (k-1): 472251
  edges ( k ): 541016
  state: dynamic
  W stats: {'$': 215, 'A': 114512, 'C': 155526, 'G': 154184, 'T': 116579}
  F stats: {'$': 5, 'A': 114184, 'C': 155966, 'G': 155009, 'T': 115852}
  dummy source edges: 2620
  dummy sink edges: 214
  ========================================================
  $ $exe annotate -i test_graph.dbg -a test_graph.column.annodbg -o test_annotation --anno-header --header-delimiter '|' $file
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 540664
  density: 3.550326e-03
  representation: column
  ========================================================
  $ $exe query -i test_graph.dbg -a test_annotation.column.annodbg <(cat $file | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000491020.1|ENSG00000116254.17|OTTHUMG00000000952.6|OTTHUMT00000002828.2|CHD5-205|CHD5|646|nonsense_mediated_decay|\tnonsense_mediated_decay:ENSG00000116254.17:OTTHUMG00000000952.6:CHD5:ENST00000491020.1:OTTHUMT00000002828.2:CHD5-205:646 (esc)
  1\tENST00000484532.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003213.2|RPL22-207|RPL22|449|protein_coding|\tprotein_coding:449:ENST00000484532.5:ENSG00000116251.9:OTTHUMG00000000953.5:OTTHUMT00000003213.2:RPL22-207:RPL22 (esc)
  2\tENST00000234875.8|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002830.1|RPL22-201|RPL22|2078|protein_coding|\tprotein_coding:2078:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000234875.8:OTTHUMT00000002830.1:RPL22-201 (esc)
  3\tENST00000480661.1|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002832.2|RPL22-206|RPL22|2279|retained_intron|\tretained_intron:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000480661.1:OTTHUMT00000002832.2:RPL22-206:2279 (esc)
  4\tENST00000497965.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002831.2|RPL22-208|RPL22|786|protein_coding|\tprotein_coding:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000497965.5:OTTHUMT00000002831.2:RPL22-208:786 (esc)
  5\tENST00000471204.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003678.3|RPL22-205|RPL22|524|protein_coding|\tprotein_coding:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000471204.5:OTTHUMT00000003678.3:RPL22-205:524 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe annotate -i test_graph.orhashdbg -o test_annotation --anno-header --header-delimiter '|' $file_small
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  514
  objects: 44919
  density: 2.531743e-02
  representation: column
  ========================================================
  $ $exe query -i test_graph.orhashdbg -a test_annotation.column.annodbg <(cat $file | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000491020.1|ENSG00000116254.17|OTTHUMG00000000952.6|OTTHUMT00000002828.2|CHD5-205|CHD5|646|nonsense_mediated_decay|\t (esc)
  1\tENST00000484532.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003213.2|RPL22-207|RPL22|449|protein_coding|\t (esc)
  2\tENST00000234875.8|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002830.1|RPL22-201|RPL22|2078|protein_coding|\t (esc)
  3\tENST00000480661.1|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002832.2|RPL22-206|RPL22|2279|retained_intron|\t (esc)
  4\tENST00000497965.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002831.2|RPL22-208|RPL22|786|protein_coding|\t (esc)
  5\tENST00000471204.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003678.3|RPL22-205|RPL22|524|protein_coding|\t (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe extend -i test_graph.orhashdbg -a test_annotation.column.annodbg -o test_graph $file
  $ $exe stats test_graph.orhashdbg
  Statistics for graph test_graph.orhashdbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 538182
  canonical mode: no
  ========================================================
  $ $exe annotate -i test_graph.orhashdbg -a test_graph.column.annodbg -o test_annotation --anno-header --header-delimiter '|' $file
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566745e-03
  representation: column
  ========================================================
  $ $exe query -i test_graph.orhashdbg -a test_annotation.column.annodbg <(cat $file | tail -n 200) | tee test_annotation_out.tsv
  0\tENST00000491020.1|ENSG00000116254.17|OTTHUMG00000000952.6|OTTHUMT00000002828.2|CHD5-205|CHD5|646|nonsense_mediated_decay|\tnonsense_mediated_decay:ENSG00000116254.17:OTTHUMG00000000952.6:CHD5:ENST00000491020.1:OTTHUMT00000002828.2:CHD5-205:646 (esc)
  1\tENST00000484532.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003213.2|RPL22-207|RPL22|449|protein_coding|\tprotein_coding:449:ENST00000484532.5:ENSG00000116251.9:OTTHUMG00000000953.5:OTTHUMT00000003213.2:RPL22-207:RPL22 (esc)
  2\tENST00000234875.8|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002830.1|RPL22-201|RPL22|2078|protein_coding|\tprotein_coding:2078:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000234875.8:OTTHUMT00000002830.1:RPL22-201 (esc)
  3\tENST00000480661.1|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002832.2|RPL22-206|RPL22|2279|retained_intron|\tretained_intron:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000480661.1:OTTHUMT00000002832.2:RPL22-206:2279 (esc)
  4\tENST00000497965.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000002831.2|RPL22-208|RPL22|786|protein_coding|\tprotein_coding:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000497965.5:OTTHUMT00000002831.2:RPL22-208:786 (esc)
  5\tENST00000471204.5|ENSG00000116251.9|OTTHUMG00000000953.5|OTTHUMT00000003678.3|RPL22-205|RPL22|524|protein_coding|\tprotein_coding:ENSG00000116251.9:OTTHUMG00000000953.5:RPL22:ENST00000471204.5:OTTHUMT00000003678.3:RPL22-205:524 (esc)
  $ [[ $(awk '{print $1}' test_annotation_out.tsv | wc -w) \
  >         == $(awk '{print $3}' test_annotation_out.tsv | wc -w) ]] \
  >     && echo Passed: OK
  Passed: OK

  $ $exe build -k 12 -o test_graph $file
  Graph chunk with 540348 k-mers was built in *sec (glob)
  $ $exe stats test_graph.dbg --count-dummy
  Statistics for graph test_graph.dbg
  ====================== GRAPH STATS =====================
  k: 12
  nodes (k): 538182
  canonical mode: no
  ========================================================
  ====================== BOSS STATS ======================
  k: 12
  nodes (k-1): 471684
  edges ( k ): 540348
  state: fast
  W stats: {'$': 215, 'A': 114392, 'C': 155313, 'G': 153984, 'T': 116444}
  F stats: {'$': 5, 'A': 114061, 'C': 155752, 'G': 154807, 'T': 115723}
  dummy source edges: 1952
  dummy sink edges: 214
  ========================================================
  $ $exe annotate -i test_graph.dbg -o test_annotation --anno-header --header-delimiter '|' $file
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: column
  ========================================================

  $ $exe transform_anno --anno-type column -o test_annotation test_annotation.column.annodbg
  Skipping conversion: same input and target type: column
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: column
  ========================================================

  $ $exe transform_anno --anno-type row --sparse -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.row.annodbg
  Statistics for annotation test_annotation.row.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: row
  ========================================================

  $ $exe transform_anno -p 4 --anno-type brwt -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.brwt.annodbg
  Statistics for annotation test_annotation.brwt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: brwt
  ========================================================

  $ $exe relax_brwt -p 4 --arity 20 -o test_annotation_relax test_annotation.brwt.annodbg
  $ $exe stats -a test_annotation_relax.brwt.annodbg
  Statistics for annotation test_annotation_relax.brwt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: brwt
  ========================================================

  $ $exe transform_anno -p 4 --anno-type brwt --greedy -o test_annotation_pm test_annotation.column.annodbg
  $ $exe stats -a test_annotation_pm.brwt.annodbg
  Statistics for annotation test_annotation_pm.brwt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: brwt
  ========================================================

  $ $exe relax_brwt -p 4 --arity 20 -o test_annotation_pm_relax test_annotation_pm.brwt.annodbg
  $ $exe stats -a test_annotation_pm_relax.brwt.annodbg
  Statistics for annotation test_annotation_pm_relax.brwt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: brwt
  ========================================================

  $ $exe transform_anno --anno-type bin_rel_wt_sdsl -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.bin_rel_wt_sdsl.annodbg
  Statistics for annotation test_annotation.bin_rel_wt_sdsl.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: bin_rel_wt_sdsl
  ========================================================

  $ $exe transform_anno --anno-type bin_rel_wt -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.bin_rel_wt.annodbg
  Statistics for annotation test_annotation.bin_rel_wt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: bin_rel_wt
  ========================================================

  $ $exe transform_anno --anno-type flat -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.flat.annodbg
  Statistics for annotation test_annotation.flat.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: flat
  ========================================================

  $ $exe transform_anno --anno-type rbfish -o test_annotation test_annotation.column.annodbg
  $ $exe stats -a test_annotation.rbfish.annodbg
  Statistics for annotation test_annotation.rbfish.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: rbfish
  ========================================================

  $ $exe annotate -i test_graph.dbg -o test_annotation --anno-type row --anno-header --header-delimiter '|' $file
  $ $exe stats -a test_annotation.row.annodbg
  Statistics for annotation test_annotation.row.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: row
  ========================================================

  $ $exe transform_anno --anno-type column -o test_annotation test_annotation.row.annodbg
  TODO
  $ $exe stats -a test_annotation.column.annodbg
  Statistics for annotation test_annotation.column.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: column
  ========================================================

  $ $exe transform_anno --anno-type row --sparse -o test_annotation test_annotation.row.annodbg
  Skipping conversion: same input and target type: row
  $ $exe stats -a test_annotation.row.annodbg
  Statistics for annotation test_annotation.row.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: row
  ========================================================

  $ $exe transform_anno --anno-type brwt -o test_annotation test_annotation.row.annodbg
  TODO
  $ $exe stats -a test_annotation.brwt.annodbg
  Statistics for annotation test_annotation.brwt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: brwt
  ========================================================

  $ $exe transform_anno --anno-type bin_rel_wt_sdsl -o test_annotation test_annotation.row.annodbg
  $ $exe stats -a test_annotation.bin_rel_wt_sdsl.annodbg
  Statistics for annotation test_annotation.bin_rel_wt_sdsl.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: bin_rel_wt_sdsl
  ========================================================

  $ $exe transform_anno --anno-type bin_rel_wt -o test_annotation test_annotation.row.annodbg
  $ $exe stats -a test_annotation.bin_rel_wt.annodbg
  Statistics for annotation test_annotation.bin_rel_wt.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: bin_rel_wt
  ========================================================

  $ $exe transform_anno --anno-type flat -o test_annotation test_annotation.row.annodbg
  $ $exe stats -a test_annotation.flat.annodbg
  Statistics for annotation test_annotation.flat.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: flat
  ========================================================

  $ $exe transform_anno --anno-type rbfish -o test_annotation test_annotation.row.annodbg
  $ $exe stats -a test_annotation.rbfish.annodbg
  Statistics for annotation test_annotation.rbfish.annodbg
  =================== ANNOTATION STATS ===================
  labels:  4401
  objects: 538182
  density: 3.566741e-03
  representation: rbfish
  ========================================================
