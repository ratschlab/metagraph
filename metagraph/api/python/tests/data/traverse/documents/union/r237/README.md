# One seed across RefSeq r237's 11 graphs, as the service exported it

The round trip of docs/DESIGN-graphlet-merge.md OQ21 (`api/python/tests/test_traverse_union.py`,
`TestR237RoundTrip`). Public RefSeq data: no private sequence or record.

- Source: the staging service (refseq-experimental, release refseq237), 2026-11-16, traversal job
  2cb81308-3733-4512-8f54-6c5455f309a5 (preset walk, no labels, all 11 graphs), the e2e seed
  TRAVERSAL_SEED_16S_200BP (E. coli K-12 MG1655 rrsA, NC_000913.3:4035531-4035730, 200 bp: `seed.fa`).
- `bact_other.mgt.gz`, `bact_pseudomonadota.mgt.gz`: the two answered parts, byte for byte the `.mgt.gz` the
  service's mgt export delivered (`2cb81308-...-96740db9-8e87-4da6-9183-4432d713152d-0.mgt.gz` and
  `2cb81308-...-7ff8cad9-c43e-4144-8ceb-a2100b1951f8-0.mgt.gz`; sha256 96a48164... and 4d436db7...). Each holds
  the standalone document (H, then the J envelope with the pair's paths and identity, then the body); the body
  sha256 (`g.dump(envelope=False)`) is cda902d0... and 6d48f657....
- `absent.json`: the nine pairs whose traversal record says `not_in_graph` (no graphlet), as `{pair, result,
  error}` inputs of `union_pairs()`, built from the job's records (pair paths, `index_fp`, `kmers_present` of
  170, the error text with the graph runs).

The job is deleted; these files are the only copies.
