# Resolve support and extend annotated graph paths

`metagraph traverse` extends exact seed sequences using the graph topology and its existing annotation. It uses the same binary, graph formats, annotation formats and `--mmap` loading as regular queries. No new index format is required.

The workflow has three parts:

1. **Resolve** a query into intervals present in the graph and intervals supported by each label.
2. **Select** exact seeds from those intervals, choosing the longest intervals, the strongest shared support, or explicit intervals.
3. **Traverse** left, right or both ways from each seed under a strategy. The result records extension paths, supporting labels, branching, label changes and why exploration stopped.

This complements alignment: an aligner can supply an exact graph sequence to resolve and extend. Traversal does not perform approximate alignment or select an alignment scoring model.

## Support semantics

- `support: "kmer"` requires the label on each traversed k-mer. This works with ordinary binary annotations. Repeated k-mers can connect contexts from different occurrences within a label; the result does not certify a contiguous source sequence.
- `support: "trace"` follows consecutive annotation coordinates. It requires coordinate annotations and a BASIC graph. Use `branching.on_reconverge: "keep"`; merging evidence from different parents cannot justify a coordinate trace.
- Header labels, selected with `labels.seed_label_kind: "header"`, also restrict traversal to one record in the existing `.seqs` mapping. Column coordinates alone cannot distinguish adjacent records in a column. Build the mapping with `annotate --anno-filename --index-header-coords` using the final annotation's input files in column order. Record headers must be unique within each file.

The two directions are reported separately. If one side reaches an end, the other can continue. To combine left and right results as source contexts, retain their label and coordinate evidence; arbitrary combinations of arms are not certified source sequences.

## CLI

Save a resolve request as `resolve.json`:

```json
{
  "sequence": "ACGTA",
  "discover": {},
  "select": {"policy": "longest_first"}
}
```

Run it against an existing index:

```sh
metagraph traverse --resolve -i graph.dbg -a annotation.column.annodbg --mmap resolve.json
```

Then supply the selected sequence and labels as a traversal request:

```json
{
  "seeds": [{"sequence": "ACGTA", "labels": ["A", "B"]}],
  "strategy": {
    "direction": "both",
    "support": "kmer",
    "bounds": {"max_extension_bp": 100},
    "output": {"detail": "full"}
  }
}
```

```sh
metagraph traverse -i graph.dbg -a annotation.column.annodbg --mmap traverse.json
```

Omit a seed's `labels` to derive its permitted labels from the annotation. `direction` accepts `left`, `right` and `both`. Strategies can constrain exploration by label support or annotate graph paths, price label switches, limit branching, and bound extension length, time, memory and logical work. Unknown fields and incompatible settings are rejected. The [core specification](SPEC-labeled-traversal-core.md) describes the strategy and result fields.

Multiple request files produce one JSON result per file. `--json` selects compact output. A failed request sets a nonzero exit status; subsequent valid request files are still processed.

## HTTP and Python

The existing `server_query` command exposes:

| Route | Purpose |
| --- | --- |
| `POST /resolve` | Support profiles and optional seed selection |
| `POST /traverse` | Seed extension |
| `GET /traverse/capabilities` | Support modes and limits of an index |
| `GET /capabilities` | Server mode, routes and attempt contract |
| `POST /traverse/cancel` | Cancel or prevent an attempt by ID |
| `GET /traverse/attempt/{id}` | Query a running or retained attempt |

POST bodies use the CLI request format. On a server started from a CSV, select one index with `graph: "name"` or `graphs: ["name"]`. If that name spans several physical graphs, also supply `graph_path`. Different annotations of one graph must have separate names. A missing seed in a selected shard is reported per seed.

Both existing Python clients expose the same operation:

```python
from metagraph.client import GraphClient

client = GraphClient("localhost", 5555)
profile = client.resolve("ACGTA", discover={}, select={"policy": "longest_first"})
seeds = [{key: seed[key] for key in ("seed_id", "sequence", "labels")}
         for seed in profile["selection"]["seeds"]]
result = client.traverse({
    "seeds": seeds,
    "strategy": {"direction": "right", "bounds": {"max_extension_bp": 100}}
})
```

`traversal_capabilities`, `cancel_traversal` and `traversal_attempt` expose the corresponding routes. `timeout` sets the HTTP client's timeout; request budgets control server work. Returned traversal objects remain JSON dictionaries. HTTP errors raise `GraphApiError`, a `RuntimeError` subclass whose `status_code` and `response` preserve cancellation tombstones, attempt usage and other structured error details.

## Budgets, mmap and results

Server defaults cap each seed at 30 seconds, each request at 64 seeds, each seed at 100,000 bases and its derived set at 10,000 labels. The `--traverse-max-*` options configure these caps; the CLI applies the request's own limits. The server reports clamps explicitly.

Walks stopped by a budget return their completed prefix and stop diagnostics. Check `outcome`, arm completeness and `limitations` before interpreting the result as exhaustive. Memory limits bound a model of the working state and output and are soft; deadlines are checked between work chunks, not asynchronously inside every backend operation. A resolve response that cannot be serialized within its requested deadline returns HTTP 503.

Traversal reads the server's existing resident or memory-mapped index. It does not reload an index into RAM per request: `in_ram: true` is rejected. Start the server without `--mmap` when that index should be resident in RAM.

The reverse header dictionary is allocated on the first traversal request for an index with a `.seqs` mapping. Ordinary `/search`, `/align`, `/stats` and capability probes do not build it. That first request's preparation time is reported as `timing.prepare_ms`, before the per-seed work budget. Subsequent requests reuse the dictionary. Its memory scales with the number of distinct headers. A request-local row-diff cache defaults to at most 128 MiB and is further limited by an explicit memory budget; `--traverse-path-cache-mb 0` disables it.

`output.detail` accepts `summary`, `tree`, `full` and `graphlet`. Graphlets contain Gunnar's MGT v1 representation; JSON and graphlets preserve the same traversal evidence. Optional `output.coordinates: true` reports record-coordinate evidence where supported. The [graphlet design reference](DESIGN-traverse-graphlet.md) documents the format.

Index manifests are optional. Without one, responses still identify the graph's regime, size and label metadata, but do not certify that two independently loaded bundles have identical contents. `--index-manifest` and `--index-name` describe a single index; a server CSV can add optional manifest and namespace columns. `traverse --index-inventory -i GRAPH -a ANNOTATION` lists the files a manifest must cover.

## Integration provenance

The resolver, label oracle, walker, budget-aware annotation decoding and wire codec are extracted from `gr/labeled-traversal` at `d5d26baa24b7976627a715d912e9926be38fd127`. The existing server loader and search-cache handler are retained. HTTP traversal is implemented as a separate route module using that loader; the Python methods extend the existing clients.

Alignment and traversal share `AnnotationKeyMapper` for oriented-node mapping and the existing batched matrix interfaces for owned annotation rows. Their label dictionaries, coordinate chains and search budgets retain their separate policies.

The extraction retains the selection corrections from #670: longest-first seeds include every containing support run, explicit seeds respect `min_block_bp`, and explicit seed lists require the explicit selection policy. Graph selection also rejects a mismatched `graph_path` when a name contains only one pair.

Header lookup preparation is deferred to traversal use to avoid imposing its allocation on ordinary query servers. The original branch's per-request RAM-loading scheduler, pattern-search engine, separate Python graphlet workspace and search service are separate integration work. The two imported design references describe that broader branch as well as the traversal contract implemented here.
