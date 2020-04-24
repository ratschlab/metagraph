## Annotation of GTEx graphs
This directory contains the relevant scripts to generate the annotations of the GTEx graphs.

The annotation can be generated on a per-sample level and on a per abundance level.

### Annotation per abundance
1. Generate a column annotation file for each abundance-level input: `run_annotate_sequences_per_abundance.sh`
2. Combine abundance-level column annotations into joint matrix: `run_annotate_sequences_per_abundance_collect.sh`
3. Convert column annotation into BRWT representation: `run_annotate_sequences_per_abundance_convert_to_brwt.sh`

### Annotation per sample
There are multiple ways a per sample annotation can be generated:

1a. Based on the abundance-level annotation: `run_relabel_abundances_to_samples.sh`
1b. Directly based on sample-level sequence files:
    - first extract sequences: `run_extract_sample_sequences_array.sh`
    - annotate sequence per sample `run_annotate_sequences_per_sample.sh`
2. Convert sample-level annotation to BRWT: `run_annotate_sequences_per_sample_convert_to_brwt.sh`
