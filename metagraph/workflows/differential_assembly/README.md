# MetaGraph count-aware differential assembly

This workflow runs the differential assembly implementation on MetaGraph's
`count_diff_assem` branch. It indexes per-sample k-mer abundances, constructs a
joint primary de Bruijn graph, stores one count-aware annotation column per
sample, tests k-mers (or graph unitigs) for differential abundance, and
assembles the significant subgraphs for both sides of every comparison.

## Pipeline

1. Validate the sample manifest and comparison definitions.
2. Count canonical k-mers in each sample with KMC. Counts of one are retained
   at this stage so cleaning remains an explicit statistical-analysis choice.
3. Build one counted canonical MetaGraph per sample.
4. Export each sample graph as primary FASTA plus the matching
   `.kmer_counts.gz` sidecar.
5. Build a cohort-wide canonical graph, export its primary k-mers, and build
   the final primary joint graph.
6. Annotate that graph separately for every sample with `--count-kmers`. The
   FASTA filename is the sample label used by the differential-assembly JSON.
7. For per-k-mer tests, stream the count-bearing column annotations directly.
   For `test_by_unitig: true`, first merge the columns into the branch's
   count-aware `int_row_flat` representation; unitig aggregation is not
   implemented for multiple column annotations in this branch.
8. Run `metagraph assemble` once per comparison. The resulting FASTA contains
   sequences from both the `in` and `out` differential subgraphs, identified
   by their headers. Optional raw p-values are written to `.all.pvals`.

## Requirements

- MetaGraph built from `count_diff_assem` (the DNA alphabet build)
- KMC
- Snakemake 7 or newer
- GNU `time` (optional; set `gnu_time_cmd: ""` to disable it)

The workflow-only Python dependencies can be installed with
`pip install -r requirements.txt`.

The Biomed installation used by the example configuration is:

```text
/cluster/customapps/biomed/grlab/users/ttanna/metagraph/metagraph/build/metagraph
```

## Configure and run

Work from this directory:

```bash
cd metagraph/workflows/differential_assembly
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --cores 16 --printshellcmds
```

Edit `config/config.yaml` and `config/samples.tsv` first. The manifest has at
least these columns:

```text
sample  reads  phenotype
case_1  /data/case_1.fastq.gz  case
ctrl_1  /data/ctrl_1.fastq.gz  control
```

It is tab-separated. Repeat a sample on multiple rows to supply paired or split
read files; its metadata values must be identical on every row. Each comparison
selects `in_values` and `out_values` from one metadata column. Values are
treated as strings and the two selected sample sets must be non-empty and
disjoint.

Comparisons may also be generated automatically from categorical metadata:

```yaml
auto_comparisons:
  enabled: true
  metadata_columns: [phenotype, study_site]  # empty means all metadata columns
  missing_values: ["", "NA", "NaN", "null"]
  minimum_samples_per_group: 2
```

For every eligible value, this creates a `<column>_<value>_vs_rest`
comparison. The `in` side contains that value and the `out` side contains all
other non-missing values in the same column. Names are made filesystem-safe,
categories below the configured sample minimum are skipped, and explicit and
automatic comparisons may be used together. An optional `assembly` mapping
inside `auto_comparisons` applies the same statistical overrides to every
generated comparison.

The defaults follow the manuscript workflow: count-aware testing, singleton
filtering (`min_count: 2`), a 0.05 significance threshold, and the negative
binomial exact test. Set `test_by_unitig: true` to aggregate per-k-mer p-values
over unitigs with the implementation's Cauchy combination step. All statistical
fields may be overridden inside an individual comparison's `assembly` mapping.

Resource settings use the standard Snakemake `threads`, `mem_mb`, and `runtime`
resources. A cluster profile may map them to scheduler requests. Intermediate
KMC databases and per-sample graph binaries are marked `temp`; pass
`--notemp` when they should be retained.

## Important implementation constraints

- `count_width` must be wide enough for `kmc_max_count`; the workflow validates
  this before constructing the DAG.
- The `.fasta.gz` and `.kmer_counts.gz` files are a pair. Renaming or moving
  only one of them makes MetaGraph silently fall back to counts of one.
- Annotation labels are the exact FASTA paths passed to `metagraph annotate`.
  Differential rule JSON files are generated from those same paths to prevent
  label mismatches.
- Automatically generated comparisons use true one-vs-rest groups. Missing
  values are excluded from both sides rather than being treated as controls.
- The current branch emits one FASTA per comparison, not separate `in` and
  `out` files. With enumeration enabled, headers are prefixed with
  `<comparison>.in.` or `<comparison>.out.` (and optionally `.shared.`).
