## Extended graph build
These scripts take the GTEx graph as a basis and build a joint graph also containing variation from
the gnomAD project as well as the human reference genome sequence (hg38). The build process contains
the following steps:

1. Take gnomAD vcfs and build a graph using the hg38 reference sequence: `run_build_gnomad_graph.sh`
(builds one graph per chromosome)
2. Transform the gnomeAD graphs into fasta sequences: `run_gnomad_graph_to_fasta.sh`
3. Transform the input GTEx graph into fasta sequences: `run_gtex_graph_to_fasta.sh`
4. Build the joint extended GTEx graph in individual chunks: `run_build_graph_chunked.sh`
5. Collect the chunks of the joint graph build: `run_build_graph_chunked_collect.sh`
