# Server API

## TODO
* [ ] (high) robust against malformed input (currently it still may crash under certain circumstances)
* [ ] (low-medium) implement mechanism to tell client to come back after a while to check whether computation has been done (useful for longer-running queries)
* [ ] (low) `fwd-and-reverse` (query its reverse complement as well)
* [ ] (low-medium) `print_signature`:  mask of k-mer: mask of present kmers. one mask per sequence
* [ ] (low) don't write temporary file for query (or only do so for fast queries)
* [ ] in case of many query sequences: generate query graph on client, send graph/contigs to server
