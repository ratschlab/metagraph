# Labeled traversal (experimental, schema version 1)

`/resolve` finds where a query sequence is supported by which labels and proposes seeds.
`/traverse` extends seeds in the graph, following only the labels that carry them.

Both are experimental: the schema may still change, and the routes are off unless the server
is started with `--experimental-traverse`. They serve a single index (`-i` / `-a`), not a
server started with a list of graphs.

```
metagraph server_query -i graph.dbg -a annotation.annodbg --experimental-traverse

metagraph traverse -i graph.dbg -a annotation.annodbg request.json            # a /traverse request
metagraph traverse --resolve -i graph.dbg -a annotation.annodbg request.json  # a /resolve request
```

| Route | Method | |
|---|---|---|
| `/resolve` | POST | a support profile of a query, and optionally seeds selected from it |
| `/traverse` | POST | the extension of one or more seeds |
| `/traverse/capabilities` | GET | the index (`k`, regime, alphabet, header labels) and the server's limits |

A request is a JSON object. A field the schema does not list is an error (HTTP 400; the CLI
exits with status 1), so that a misspelled field is never ignored.

## Labels

A label is either an annotation **column** (for example a taxid or a file) or, on an index with
k-mer coordinates and a sequence header index (`annotate --index-header-coords`), one indexed
**sequence** of a column, named by its FASTA **header** (for example an accession). A name is
resolved as a column first, then as a header.

## `/resolve`

```json
{
  "sequence": "ATGGAATTGCCC...",
  "labels": ["NZ_CM008882.1", "562"],
  "discover": {"max_labels": 1000, "kind": "column"},
  "min_block_kmers": 1,
  "select": {"policy": "max_support", "max_seeds": 10, "min_block_bp": 0,
             "max_labels_per_seed": 1000, "merge_overlapping": true}
}
```

| Field | | Default |
|---|---|---|
| `sequence` | the query, at most 1,000,000 bp | required |
| `labels` | the labels to profile, at most 100,000 | exactly one of `labels` and `discover` |
| `discover.max_labels` | keep the labels supporting the most k-mers | 1000 |
| `discover.kind` | `column` or `header` | `column` |
| `min_block_kmers` | the shortest support run that becomes a seed candidate | 1 |
| `select` | also select seeds | none |
| `select.policy` | `max_support`: the interval contained in the most labels' support runs, repeatedly, outside the intervals picked before; `longest_first`: the longest candidates | `max_support` |
| `select.max_seeds`, `select.min_block_bp` | at most that many seeds (at most 1000), each at least that long (at least k) | 10, 0 |
| `select.max_labels_per_seed` | the labels kept per seed, by the FNV-1a hash of their names; the others are counted with a digest | 1000 |
| `select.merge_overlapping` | merge overlapping picks onto their common interval | true |

The response:

```json
{
  "k": 31, "regime": "basic", "num_kmers": 783,
  "graph_runs": [[0, 783]],
  "labels": [{"label": "NZ_CM008882.1", "kind": "header", "kmers_supported": 783,
              "runs": [[0, 783]]}],
  "labels_truncated": null,
  "candidates": [{"kmer_interval": [0, 783], "bp_interval": [0, 813],
                  "labels": ["NZ_CM008882.1"]}],
  "selection": {
    "policy": {...}, "counts": {"candidates": 3, "eligible": 3, "selected": 1},
    "seeds": [{"seed_id": "4c1f...", "sequence": "ATGG...", "kmer_interval": [0, 783],
               "labels": ["NZ_CM008882.1"], "overlaps_with": [],
               "label_population": {"supporting_total": 1, "included": 1, "dropped_count": 0}}]
  }
}
```

Intervals are half-open k-mer intervals `[begin, end)`, except `bp_interval`, which is in
bases. `graph_runs` are the runs of the query's k-mers present in the graph; a label's `runs`
are the maximal runs of k-mers it carries. With discovery, `labels_truncated` states how many
labels there were and where the cut fell. A selected seed is valid `/traverse` input as it
stands, with its `seed_id`.

## `/traverse`

```json
{
  "seeds": [
    {"seed_id": "4c1f...", "sequence": "ATGGAATTGCCC...", "labels": ["NZ_CM008882.1"]}
  ],
  "strategy": {
    "schema_version": 1,
    "direction": "both",
    "labels": {"max_seed_labels": 1000, "seed_label_kind": "auto"},
    "branching": {"max_label_branches": 0},
    "bounds": {"max_extension_bp": 5000, "time_budget_ms": 30000},
    "output": {"sequences": true}
  }
}
```

| Field | | Default |
|---|---|---|
| `seeds` | at most `--traverse-max-seeds` seeds, each at most `--traverse-max-seed-bp` long, every k-mer in the graph | required |
| `seeds[].seed_id` | the caller's id, compared with the id the seed validates to | none |
| `seeds[].labels` | the permitted labels; omitted: derived from the seed | derived |
| `direction` | `both`, `left` or `right` | `both` |
| `labels.max_seed_labels` | the most labels a derived set keeps, in ascending (column, sequence) order (at most 100,000) | 1000 |
| `labels.seed_label_kind` | what a derived label is: `column`, `header`, or `auto` (header on an index with a header index) | `auto` |
| `branching.max_label_branches` | how many times a label may continue on more than one successor: an integer, or `"unlimited"` | 0 |
| `bounds.max_extension_bp` | the extension length per arm | 5000 |
| `bounds.time_budget_ms` | per seed, so a request may take up to its number of seeds times this; clamped to `--traverse-max-time-ms` | `--traverse-max-time-ms` |
| `output.sequences` | include the segments' sequences | true |

### The walk

A seed's permitted labels are the labels named, or the labels carrying every k-mer of the
seed. A named label that does not carry every seed k-mer is dropped and reported with the
k-mer runs it carries.

From the seed's last k-mer (the right arm) and its first (the left arm) the walk extends
level by level, one base per level. Each path keeps the labels that carry every k-mer walked
so far. A label leaves a path with a reason:

| Reason | |
|---|---|
| `dead_end` | the node has no successor (`hairpin`: its only continuation was a hairpin) |
| `label_lost` | successors exist, none carries the label |
| `branch` | the label is on several successors and its branch allowance is spent (`split_limit`: the path's 64 splits are spent) |
| `edge_reuse` | the label's continuation takes an edge this path took before |
| `edge_reuse_rc` | ... in the opposite orientation (an inverted repeat) |
| `rejoined_seed` | the continuation steps back into the seed (for example a circular record) |
| `max_extension_bp` | the extension length was reached |
| `time_budget`, `max_steps`, `max_live_paths`, `max_paths`, `max_output_bp` | a bound stopped the arm (see below) |

Where the labels of a path continue on different successors, the path splits. Paths that
meet again are kept apart, so each arm is a tree of segments.

The walk rule makes the walks finite on a graph with cycles: a path never takes an edge
((k+1)-mer) twice, in canonical graphs in either orientation, never steps back into the seed,
and never takes a self-reverse-complementary step.

**Bounds.** Besides the extension length and the time budget, fixed caps bound a seed's work:
1,000,000 steps, 1000 live paths and 10,000 paths per arm, and 2,000,000 output bases per arm.
The time budget and the steps stop the whole seed; the other caps stop one arm. A cap trips
between two heads of one level: `complete_to_bp` is the extension up to which every admissible
walk is present.

### The response

```json
{
  "schema_version": 1,
  "strategy": {...},
  "results": [
    {
      "seed": {
        "seed_id": "4c1f...", "validated_seed_id": "4c1f...", "seed_id_mismatch": false,
        "length_bp": 813, "num_kmers": 783,
        "labels": ["NZ_CM008882.1"],
        "labels_from_seed": false, "labels_supporting_total": 1,
        "labels_dropped": 0, "labels_dropped_digest": "",
        "dropped_labels": [{"label": "562", "runs": [[0, 231]]}]
      },
      "label_dict": [{"name": "NZ_CM008882.1", "kind": "header", "column": 6, "seq_id": 3}],
      "arms": {
        "left": {...},
        "right": {
          "status": "complete",
          "complete_to_bp": 5000,
          "frontier_remaining": {"live_paths": 0},
          "segments": [
            {"id": 0, "parents": [], "from_bp": 0, "length_bp": 1200, "sequence": "GCAT...",
             "labels": [0], "labels_at_end": [0],
             "events": [{"at_bp": 1200, "type": "label_end", "label": 0, "reason": "dead_end"}]}
          ],
          "splits": [],
          "paths": [{"id": 0, "segments": [0], "length_bp": 1200,
                     "end_labels": [{"label": 0, "reason": "dead_end", "branches": 0}]}],
          "steps": 1200
        }
      },
      "elapsed_ms": 12.5
    }
  ]
}
```

- `labels` and the `label` fields of an arm are indices into `label_dict`.
- `validated_seed_id` is the hex FNV-1a-64 of the case-mapped seed (in canonical graphs the
  smaller of it and its reverse complement) and its sorted labels: the id `/resolve` gives.
- A derived set over `max_seed_labels` states how many labels it dropped
  (`labels_dropped`) and the FNV-1a-64 digest of their names.
- A segment's `sequence` is in natural orientation on both arms; on the left arm a
  segment's bases precede its children's. `from_bp` is the distance of its first base from
  the seed boundary.
- A path is a leaf of the tree: `segments` from the root, `length_bp` from the seed, and
  the labels reaching its end with their reasons (`branches`: how many times the label took
  several successors). `path_reason` is set when a bound ended the path.
- A split is `divergence` when its labels part, `ambiguous` when a label follows several
  children.
- `cap_trigger` (`reason`, `at_bp`, `demand`) is set when a bound stopped the arm
  (`status: "truncated"`).

A seed whose labels cannot be derived — no label carries every k-mer (`no_carrier`), its
cheapest first k-mer carries too many labels (`too_wide`), the time budget ran out
(`time_budget`), or a derived header name is ambiguous (`ambiguous_header`) — is reported in
place, and the other seeds are still walked:

```json
{"seed": {"seed_id": "...", "length_bp": 813, "labels_from_seed": true},
 "error": {"cause": "no_carrier", "message": "...", "limit": 783, "observed": 412}}
```

Every other error (an unknown label, a seed not in the graph, too many seeds) fails the
request.

## Python client

```python
from metagraph.client import GraphClient

client = GraphClient('localhost', 5555)
profile = client.resolve(query, discover={'kind': 'header'},
                         select={'policy': 'longest_first', 'max_seeds': 1})
seed = profile['selection']['seeds'][0]
result = client.traverse([{'seed_id': seed['seed_id'], 'sequence': seed['sequence'],
                           'labels': seed['labels']}],
                         strategy={'bounds': {'max_extension_bp': 2000}})
```

## Origin

This is version 1 of the labeled traversal of `gr/labeled-traversal` (Gunnar Rätsch), with
the same walk and labels for these features. Not included: label switching and loss budgets,
`support: trace` and record coordinates, the annotate mode, memory and work budgets, the
attempt and cancel protocol, reconvergence merging, beams and quorums, continuations, and the
graphlet (MGT) output.
