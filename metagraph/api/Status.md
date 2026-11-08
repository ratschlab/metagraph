# Server API development notes

This checklist records development ideas. See the
[API documentation](../docs/source/api.rst) for usage.

## Checklist
* [ ] (high) investigate the previously reported crashes on malformed input
* [ ] (low-medium) implement mechanism to tell client to come back after a while to check whether computation has been done (useful for longer-running queries)
* [ ] (low) `fwd-and-reverse` (query its reverse complement as well)
* [x] return k-mer match signatures: supported by the `with_signature` search parameter
* [ ] (low) don't write temporary file for query (or only do so for fast queries)
* [ ] in case of many query sequences: generate query graph on client, send graph/contigs to server
