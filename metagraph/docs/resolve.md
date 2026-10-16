# Resolve query support and select traversal seeds

`metagraph traverse --resolve` and `POST /resolve` expose the resolve/select
stage from `gr/labeled-traversal` (revision `07dfa4875fcf`). They turn partial
query matches into supported intervals and explicit seed sequences with label
sets. Traversal beyond those seeds is a separate increment.

This first extraction profiles **annotation columns** using **k-mer membership**.
A column supporting every k-mer of an interval need not contain that interval
consecutively in one record. Coordinates and `.seqs` are not required; loading
`.seqs` does not change the label level. Header labels and coordinate traces are
rejected, rather than answered using weaker evidence.

## CLI

Use the existing executable, with one JSON request per file:

```sh
metagraph traverse --resolve --mmap -i graph.dbg -a annotation.brwt.annodbg request.json
```

```json
{
  "sequence": "ACGTTGCAACGT",
  "discover": {"kind": "column", "max_labels": 100},
  "select": {"policy": "max_support", "max_seeds": 3, "min_block_bp": 8}
}
```

Instead of `discover`, provide `"labels": ["sample_A", "sample_B"]` to profile
explicit columns. Exactly one is required. Several request files produce one
compact JSON response per file; an erroneous request produces an `error` object,
later files are processed, and the process exits nonzero if any request failed.
There is no extra installed executable.

## HTTP and Python

A single-index `metagraph server_query` accepts the same request as JSON at
`POST /resolve`. It uses the already loaded graph and annotation. Multi-index
server routing is not part of this increment and is explicitly rejected.

```python
from metagraph.client import GraphClient

client = GraphClient("localhost", 5555)
result = client.resolve(
    "ACGTTGCAACGT",
    discover={"kind": "column", "max_labels": 100},
    select={"policy": "max_support", "max_seeds": 3, "min_block_bp": 8},
)
```

`GraphClientJson.resolve` returns the same hierarchical JSON result. This is a
profile/seed response, not a search DataFrame. No web or MCP integration is added.

## Contract

Intervals are zero-based and half-open. `kmer_interval: [a,b]` covers k-mer starts
`a` through `b-1`; its sequence is `query[a:b+k-1]`.

The response reports `k`, `regime`, `num_kmers`, `support: "kmer"`, `graph_runs`,
per-column `labels` with matched-k-mer counts and maximal `runs`, and `candidates`
grouping identical maximal runs. Discovery ranks columns by matched-k-mer count,
breaking ties by column ID. `labels_truncated` reports discarded columns; it is
null when none were discarded. Explicit columns retain request order.

Optional `select` retains the prototype's three policies:

- `max_support` (default): choose intervals supported by the most labels, with
  longer intervals then earlier starts breaking ties; later rounds exclude
  already chosen intervals.
- `longest_first`: choose from the ranked maximal-run candidates.
- `explicit`: use `seeds` containing `kmer_interval` and `labels`; report supplied
  labels that do not cover an interval separately.

Selection returns the sequence, interval, supporting labels, deterministic seed
ID, label-population truncation and overlaps for each seed. Its defaults are
`max_seeds: 10`, `min_block_bp: 0` (at least k), `max_labels_per_seed: 1000`,
`label_order: "hash"`, `sample_seed: 0`, and `merge_overlapping: true`.
`label_order` also accepts `column_id` and `kmers_supported`. Overlapping picks
may be reduced to their common interval, following the prototype's policy.

`min_block_kmers` (default 1) filters candidates. `run_format` accepts `intervals`
(default) or the existing `x<n>o<n>` RLE notation. `bounds.max_query_bp` defaults
to 1,000,000. Counts and limits must be nonnegative integers; limits requiring a
nonempty result budget must be positive. Unknown fields and unsupported modes
are rejected. These are input/result limits, **not** hard memory or execution-time
budgets: discovering top labels still reads all relevant rows and accumulates
support for the columns encountered.

The selection algorithm, seed IDs and response fields come from the prototype.
The annotation adapter uses the current matrix API with batches of 64 distinct
rows; it avoids importing the walker's label oracle, decoding-budget extensions
and service lifecycle. The PR comparison evidence covers this adapter separately.
No topology sidecar, graph format change or dependency update is required.
