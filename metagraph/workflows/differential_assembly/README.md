# Metadiff

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

This pipeline implements a workflow to build and annotate [metagraphs](https://metagraph.ethz.ch/static/docs/index.html) from metagenomic data. 

## Configuration

A configuration file (`config.yaml`) is to be set up in the `config` directory. This file requires the path to the input files, path to the metadata folder, and also the output paths to save the graphs and log files (these output directories are created if they do not exist). It also requires the extension of the input files (can be `.fasta`/`.fa`, `.fastq`/`.fq`, `.fasta.gz`/`.fasta.gz`, `.fastq.gz`/`.fastq.gz`). Further, you can setup the differential assembly parameters (`in_min` and `out_max`) as described [here](https://metagraph.ethz.ch/static/docs/sequence_assembly.html#differential-assembly). Finally, the number of available threads and memory per thread can be specified.

## Metadata

A tab separated metadata file (`.txt` or `.tsv`) with the first column containing sample identifiers has to be supplied at the metadata path. The sample identifiers should be identical to the names of the input files (without extensions). The workflow ssumes the metadata to be categorical, parses each column and creates JSON files for each metadata category in the metadata file. It then uses these JSON files to generate differential assemblies for all metadata categories.

## To Do

For functionally annotating the differential assemblies, a [`DIAMOND blastx`](https://github.com/bbuchfink/diamond) rule will be added.

