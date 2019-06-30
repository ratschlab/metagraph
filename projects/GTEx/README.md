## Construction of graphs

#### Generation of individual kmer count files
```bash
graph_build/run_count_kmers.sh
```
The above script relies on calling `graph_build/count_kmers.sh`.

#### Join KMC count files into merged files of 100
```bash
graph_build/run_merge_kmc_dbs_all.sh
```
The above script relies on calling `graph_build/merge_kmc_dbs.sh`.

#### Transform the joint KMC count files into individual graphs
```bash
graph_build/run_build_single_graphs_from_merged_kmc_files.sh
```

#### Extract FASTA sequences from the joint KMC graphs
```bash
graph_build/run_build_single_graphs_from_merged_kmc_files.sh
```

#### Build a joint graph over all samples based on the extracted joint sequences
```bash
run_build_joint_graph_from_sequences_extracted_from_merged_kmc_graphs.sh
```

#### In addition, generation of a single graph per unmerged KMC file as an input for annotation
```bash
graph_build/run_build_single_graphs_from_single_kmc_files_array.sh
```

## Annotation

### Annotation per sample
#### From graphs built per sample (from the per-sample KMC files), extract FASTA sequences
```bash
annotation/run_extract_sequences_array.sh
```

#### Annotate the joint graph with one column per sample based on the extracted FASTA sequences
```bash
annotation/run_annotate_sequences.sh
```

#### Collect individual annotation columns into joint matrix
```bash
annotation/run_collect_annotation_columns.sh
```

#### Transform collected column matrix into BRWT matrix
```bash
annotation/convert_annotation_brwt.sh
```

### Annotation of tissues
#### Collapse/rename individual columns by tissue
```bash
### variant 1 (collapse from samples):
annotation/run_collapse_annotation_columns_by_tissue.sh

### variant 2 (directly construct with tissue label):
annotation/run_annotate_by_tissue.sh
```
(The above scripts rely on preprocessing the metadata with `annotation/split_files_by_tissue.sh`.)

#### Transform column matrix of tissue annotations into BRWT format
```bash
annotation/run_brwt_conversion_by_tissue.sh
```

### Abundance-aware annotation
#### From graphs built per sample (from the per-sample KMC files), extract FASTA sequences filtered by abundace
(currently in 5 different abundance levels)
```bash
annotation/run_extract_abundance_filtered_sequences_array_no_kmc.sh
```
The above script relies on calling `annotation/extract_abundance_filtered_sequence_no_kmc.sh`.

#### Annotate the joint merged graph with sequences of different abundance levels
```bash
annotation/run_annotate_sequences_per_abundance_array.sh
```
(Also exists as non-array script `annotation/run_annotate_sequences_per_abundance.sh`)

#### Collect individual abundance annotation columns into joint matrix (per abundance level)
```bash
annotation/run_collect_annotation_columns_per_abundance.sh
```

#### Collect the matrices of the different expression levels into a joint matrix and convert into BRWT
```bash
annotation/convert_annotation_brwt_abundances.sh
```
