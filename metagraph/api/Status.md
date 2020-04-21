# Server API Capabilities

(WIP)


## General Capabilities

* [x] HTTP
* [x] JSON strings as response payload, including error messages
* [x] parallel processing of HTTP queries
* [x] support compression
* [ ] (high) robust against malformed input (currently it still may crash under certain circumstances)
* [ ] (low-medium) implement mechanism to tell client to come back after a while to check whether computation has been performed done
      (useful for longer-running queries)


## Design Principles

* don't return redundant information

## Column Labels

List of column labels as simple JSON list

### Python API Support

Supported.

## Query

Annotate with exact matches.

Takes fasta files

Supported parameters:
* [x] `num_labels` (return top `num` labels)
* [x] `discovery_fraction` [0, 1.0] (minimum fraction of kmers of sequence found in annotation)
* [x] `fast`: fast queries, i.e. first generate query graph (currently only exposed for test purposes)
* [x] `align`: align query first, include aligned sequence and score in result
* [ ] (low) `fwd-and-reverse` (query its reverse complement as well)
* [ ] (low-medium) `print_signature`:  mask of k-mer: mask of present kmers. one mask per sequence


Returns:
Returns list of objects, one for each query sequence. Such an object contains following fields:
* `seq_description`: description from fasta sequence
* `sequence`: aligned sequence if align=True
* `score`: alignment score if align=True
* `results`: result for that sequence, which is a list of objects with the following fields:
    * `sampleCount`: Count
    * `sampleName`: Column name

More capabilities:

* [ ] (low) don't write temporary file for query (or only do so for fast queries)

### Python API Support

Supported. Both "raw" json query as well as result as data frame with the following columns
* `sampleCount`
* `sampleName`
* if multiple queries
    * `seq_description`: values 0..(n-1) where n is the number of sequences queried for
* if align:
    * `score`
    * `sequence`

## Align

Find closest sequence in graph

Takes fasta file.

No parameters currently supported. Currently, processing multiple sequences sequentially.

* [ ] return graph
* [ ] (low) return more than one alignment per sequence?


Returns list of result object for each query sequence. A result object contains following fields:
* `seq_description`: description from fasta sequence
* `score`: alignment score [int]
* `sequence`: aligned sequence


### Python API Support

Supported.


## Python API General Capabilities

* [ ] simplify API for single graph use.
* [ ] in case of many query sequences: generate query graph on client, send graph/contigs to server
