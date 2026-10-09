# Real-data tests for `metagraph.traverse`

This directory tests the graphlet library (`api/python/metagraph/traverse`,
`docs/DESIGN-traverse-graphlet.md`) end to end against real MetaGraph servers. It checks
retrievals from three real indexes against the server's own `detail: full` JSON and
against `/resolve`. `/resolve` is the search path: it maps k-mers and reads annotation
rows, and shares no code with the walker.

There are two ways to run it:

- **Live** (`test_real_*.py`): needs the three servers, or at least a filled retrieval
  cache in a scratch directory. It covers about 9,400 tests over 581 cached cells.
- **Offline** (`../test_traverse_real_offline.py`): runs in CI with no server and no
  scratch directory. It replays 18 small real retrievals committed under
  `../data/traverse/real/`, using the same conformance and oracle code. See
  [Offline fixtures](#offline-fixtures-ci).

The unit-test discovery of `api/python/tests` (`-p 'test_traverse_*.py'`) does not
recurse into this directory (it has no `__init__.py`). Only the offline module runs
there. pytest (used by `tox.ini`) would collect `test_real_*.py` too; they skip without a
cache. This was not tried with pytest, which is not installed on the machine that built
the fixtures.

## Files

| file | role |
|---|---|
| `realdata.py` | The harness: endpoints and skip helpers, the seed catalog (built from real `/resolve` calls), the strategy matrix (21 named strategies), the retrieval cache, the oracle client (`/resolve`, cached on disk), the record oracle (`oracle_positions`: the bases at reported record coordinates and string counts, from the source FASTA of mini_refseq, cached like the `/resolve` answers), and the independent spellers. CLI: `status`, `catalog`, `fill`, `report`. |
| `test_harness_selfcheck.py` | The harness itself: cache consistency; the two spellers agree with each other and with the library. |
| `test_real_conformance.py` | Conformance of every cached retrieval (9 checks per cell). |
| `test_real_oracle.py` | The search oracle: every cached retrieval re-checked by `/resolve` (10 areas). |
| `test_real_ops.py` | Library operations against the full JSON, the spellers and the oracle. |
| `test_real_continuations.py` | `continuation()`, `next_request()` and `deepen()`, offline and live. |
| `test_real_agent.py` | The MCP tool layer (`mcp_tools`, `store`, `client`) driven like an agent. |
| `test_real_fuzz.py` | Mutation fuzzing of real bodies (parser and tools), in a watched child process. |
| `test_real_scale.py` | The largest retrievals: time, memory, super-linearity. Slow, opt-in. |
| `offline.py` | Builds the offline fixtures (`snapshot`) and loads the suite pointed at them (`Harness`). |

## Running the live suite

### 1. Servers

The suite expects three servers on 127.0.0.1. They were started with these flags:
`-p 8 --address 127.0.0.1 --traverse-max-time-ms 120000 --index-name <name>`, plus
`--index-manifest` for mini_refseq only.

```sh
MG=metagraph/build/metagraph     # any build with /traverse and /resolve
IDX=~/metagraph-indexes
MINI=metagraph/build/mini_refseq # built by scripts/traversal/build_mini_refseq.sh

$MG server_query -i $IDX/sra/graph_primary_small.dbg \
    -a $IDX/sra/annotation.relaxed.row_diff_brwt.annodbg \
    --port 5701 --address 127.0.0.1 -p 8 --traverse-max-time-ms 120000 \
    --index-name sra_primary_small
$MG server_query -i $IDX/uhgg/graph_complete_k31.dbg \
    -a $IDX/uhgg/annotation.relaxed.relabeled.row_diff_brwt.annodbg \
    --port 5702 --address 127.0.0.1 -p 8 --traverse-max-time-ms 120000 \
    --index-name uhgg_k31
$MG server_query -i $MINI/graph_k31.dbg \
    -a $MINI/annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg \
    --port 5703 --address 127.0.0.1 -p 8 --traverse-max-time-ms 120000 \
    --index-name mini_refseq \
    --index-manifest $MINI/annotation.relaxed.relabeled.manifest.json
```

| index | port | regime | labels | identity | load |
|---|---|---|---|---|---|
| `sra` | 5701 | primary, 42 G nodes, RowDiff&lt;BRWT&gt;, no coordinates | 5,184 SRA runs (columns) | `index_fp` null: unverifiable | about 20 s, about 12 GB RSS (22 GB after a day of requests) |
| `uhgg` | 5702 | basic, 9.7 G nodes, RowDiff&lt;BRWT&gt; | 4,644 genomes (columns) | unverifiable | about 9 s, about 7 GB |
| `mini_refseq` | 5703 | basic, RowDiff&lt;BRWT&gt; with k-mer coordinates (the refseq33m format) | 9 columns, accession headers via `.seqs` | `index_fp` from the manifest: verifiable | under 1 s |

A server answers 503 while it loads; the harness treats that as "down" and skips.

### 2. Environment

| variable | default | effect |
|---|---|---|
| `METAGRAPH_REAL_SERVERS` | `sra=127.0.0.1:5701,uhgg=…:5702,mini_refseq=…:5703` | Override addresses (`name=host:port`, or a JSON object). `name=off` disables one index. |
| `METAGRAPH_REAL_SCRATCH` | `~/.cache/metagraph-traverse-real` | Root of `catalog.json`, `cache/`, `oracle_cache/`, `servers/capabilities.json` and `scale_report.json`. Never inside the repo. |
| `METAGRAPH_REAL_MAX_FULL_MB` | 40 | Conformance compares full JSON up to this size; larger cells run the body-only checks. |
| `METAGRAPH_REAL_ORACLE_LEAVES`, `_CLAIMS`, `_CHUNK_BP` | 8, 40, 60000 | Oracle sample sizes, and the bases per `/resolve` call. |
| `METAGRAPH_REAL_OPS_MAX_BODY`, `_MAX_FULL`, `_LEAF_CAP`, `_ORACLE_WALKS`, `_LABEL_CAP`, `_ANNOTATE_LEAF_CAP` | 1 MB, 32 MB, 250, 4, 16, 40 | Scope of `test_real_ops`. |
| `METAGRAPH_REAL_CONT_LEAF_CAP`, `_LIVE_PER_ARM`, `_DEEP_PER_INDEX`, `_BP` | 200, 2, 8, 150 | Scope of `test_real_continuations`. |
| `METAGRAPH_REAL_SLOW` | unset | `1` runs `test_real_scale` (a few minutes, a few GB of RAM). |
| `METAGRAPH_REAL_HEAVY` | 1 | `0` keeps the 31.75 MB body (403 MB of full JSON) off the scale ladder. |
| `METAGRAPH_REAL_NO_XFAIL` | unset | `1` runs scale product-bug tests as plain tests. All are fixed now. |

### 3. Catalog and cache

```sh
cd api/python/tests/real
python3 realdata.py status                 # servers up?, catalog and cache sizes
python3 realdata.py catalog                # seeds from real /resolve calls (--refresh to redo)
python3 realdata.py fill --workers 2       # every (seed, strategy) cell, detail graphlet AND full
python3 realdata.py report
```

The catalog has 38 seeds, plus one 3-seed batch per index (41 entries):

- **SRA:** 16S windows, random gut-genome 50-mers, and a hub 50-mer carried by 508 runs.
- **UHGG:** a 16S window, random 100 bp windows, a contig end with its switch target,
  and a window carried by 7 genomes.
- **mini_refseq:** blaNDM-1, its reverse complement, record windows, a record end, and three
  150 bp windows of an insertion sequence repeated in several records (kind `repeat`,
  `mini_rep_00..02`: only the trace strategies run on them, `NARROW_KINDS`). A catalog made
  before them gains them on the next `catalog` run, the other seeds kept.

Record coordinates (feature level 6): `trace_coords` (trace, `output.coordinates` true,
`max_coordinate_occurrences` "unlimited", radius 2 kb) and `trace_coords_column` (the same on
column labels: kind `column`) run on every mini_refseq seed; the `trace` strategy pins
`output.coordinates: false`, so that a library that asks for coordinates by itself
(`traverse_fetch` on a level-6 server) fetches what the cache holds.

`fill` fetched all 581 cells (SRA 214, UHGG 211, mini_refseq 156) in about 6 minutes.
`time_tight` cells are not deterministic: their two fetches stop at different points.
Cells whose full JSON exceeds 64 MB, or whose bodies exceed 8 MB, are marked heavy.

### 4. Running

```sh
cd api/python/tests/real
python3 -m unittest -v test_real_conformance          # cache only
python3 -m unittest -v test_real_oracle               # cache + oracle cache (or servers)
python3 -m unittest -v test_real_ops test_real_continuations
python3 -m unittest -v test_real_agent test_real_fuzz
METAGRAPH_REAL_SLOW=1 python3 -m unittest -v test_real_scale
python3 -m unittest discover -p 'test_*.py'           # everything but the slow suite
```

Every test skips cleanly when its server is down, its cache cell is missing, or (for
the oracle) the live server serves a different index than the cache was filled from.
Oracle answers are cached on disk, keyed by index identity and payload, so a second
oracle run needs no server. Without servers and cache, everything skips.

## What each area checks

**Conformance** (per cell and seed result):

- `roundtrip`: `parse(text).dump() == text`, `is_canonical`, and the body-only and
  enveloped dumps.
- `counts`: the six A counts recounted from the raw records equal A, the JSON summary
  and the full JSON.
- `rules`: `check_rules(reference=full)` reports nothing.
- `to_json`: equals the server's detail-full result after the §2.5 normalisation, which
  is written independently of the library.
- `invariants`: full-JSON invariants (label_end events, switch provenance, `end_labels`).
- `summary`: the per-seed JSON summary agrees with the body and with the full result,
  including the conservative outcome rule (DESIGN §14).
- `identity`: H identity equals the server's capabilities. mini_refseq is verifiable;
  SRA and UHGG are not.
- `seed`: S equals the request; the derived labels are the catalog's full carriers.
- `save_load`: the `.mgt` file round-trips, and a cut body never parses.

**Oracle** (`/resolve`, on a deterministic sample per arm; coordinates as in the
module docstring):

| area | what `/resolve` confirms |
|---|---|
| `a_present` | Sampled leaf walks, spelled with the seed by the harness's own speller, are fully in the graph. Without bases, the continuations are. |
| `b_claims` | Every sampled claim's displayed stretch `[evidence_from, to_bp)` carries its label on every k-mer; trace claims also under `support: trace`. |
| `b_lost` | At a `label_lost` or silent-switch end inside a segment, the label is absent on the next node and present on the last one. |
| `b_walks` | The top walks: `Walk.sequence` equals the independent spelling; every `labels_full` label covers the walk; every displayed support run is real. |
| `b_stretch` | Annotate: the maximal displayed stretches, recomputed here from the recorded sets, are real. |
| `b_maximal` | Library only: annotate `claims()` hold each label's maximal route under the §5.1 union rule, and nothing beyond it. |
| `c_routes` | A run whose displayed support starts at a merge is `route_only` at that cut, and `routes(spell=True)` spells its own real route carrying the label. |
| `d_annotate` | Recorded label sets equal `/resolve`'s labels at sampled nodes; a cut list is a subset that carries the true total. |
| `e_direct` | For sampled labels, a route of length `direct_bp` exists on which the label covers every k-mer. |
| `e_library` | Library only: `routes()` spells a route of length `direct_bp` for every label. |
| `f_positions` | Record coordinates (cells tagged `coordinates`), against the index's SOURCE records, not the server: the bases at every seed occurrence are the seed, at every run occurrence of every label walk the run's own route bases, at every claim cut at half the arm's depth the claim's own bases (the library's clipping), at every FASTA header interval (1-based closed) the record's; a seed list's total is the seed's count in the label's record(s), an unmarked run's total the count of its chain string, a `lower_bound` run's at most that. A column label's interval is checked in the record(s) whose numbering holds it. |

**The other modules:**

- **ops:** walks, claims, label_walks, routes, support_profile, support_changes,
  splits, FASTA/GFA, subgraph views and compare. Each is checked against expectations
  re-derived from the full JSON, the spellers or the oracle.
- **continuations:** spelling, coordinates and the length rule offline. Live: the
  generated request is accepted, the child continues the parent walk, and its labels
  are checked against `/resolve`.
- **agent:**
  - scripted MCP sessions per index;
  - a sweep of the list tools over every cell (ceilings, paging, totals, evidence
    block, rows equal to the full JSON);
  - cursors and structured errors.
- **fuzz:** truncation, dropped or duplicated records, count flips, huge ranges and
  integers, bad UTF-8, lone surrogates, CRLF and field corruption. Each must raise
  `GraphletFormatError`, or round-trip canonically. Runs in a child process watched for
  time and RSS; the tools must return `format_error` and store nothing.
- **scale:** the largest retrievals measured against generous bounds:
  - parse time, retained heap and `memory_bytes()`;
  - dump, to_json, walks, claims, exports and compare;
  - the store's spool and LRU, and a cold MCP session;
  - exponents fitted against the work W, to catch super-linear growth.

## Offline fixtures (CI)

`../test_traverse_real_offline.py` loads **private copies** of `realdata.py`,
`test_real_conformance.py` and `test_real_oracle.py`. Their harness points at
`../data/traverse/real/`, with every server off (`http` raises; nothing opens a
socket). It then runs the suite's own checks:

- one test per (cell, conformance check);
- one per (cell, oracle area);
- `full_invariants` on the time-budgeted cell, whose graphlet and full fetches stopped
  at different points;
- fixture-level tests: file sizes, coverage of the situations below, and two real-data
  facts. In the hairpin cells the skipped step is the self-reverse-complementary 32-mer
  A16T16; with `skip` the event's 152 labels end `Dh`, and with `follow` 151 of them
  end `edge_reuse_rc` one step later. At the contig end the lineages switch into the
  genome that continues past it.

Two rules keep CI honest:

- A check that **ran** when the fixtures were made (manifest `expect.ran`) **fails** if
  it skips now.
- A `/resolve` payload with no recorded answer **fails** with `OfflineMiss`. It means
  the oracle plan changed, so the library's output changed.

Layout of `../data/traverse/real/`:

- `manifest.json`: rows, why each cell is here, and the `expect` record.
- `capabilities.json`: per index, without process facts.
- `catalog.json.gz`: the seeds.
- per cell `<index>/<cell>.{request,graphlet,full}.json.gz`: the server's exact
  bytes, gzipped.
- per cell `<index>/<cell>.oracle.json.gz`: the recorded `/resolve` answers (and, for a cell with record coordinates, the record oracle's: per interval the digests of the bases, per string its count). They drop
  the unused `candidates` field, and keep the queried sequence only as its length and
  sha256 (the answer's key binds the full payload).

Every gzipped file is under 60 KB; the whole set is 632 KB. A rebuild with unchanged
inputs is byte-identical.

| cell | graphlet / full (KB gz) | checks run | what it covers |
|---|---:|---:|---|
| `sra/sra_hub__strict` | 18.0 / 41.5 | 16 | The default strategy on a 50-mer carried by 508 SRA runs: 508 derived seed labels, merges, runs whose displayed support starts at a merge. |
| `sra/sra_rand50_00__limit1` | 8.5 / 13.8 | 16 | Branch limit 1: merges, merge-closed (m) runs, and 8 split clones carrying their source's `route_bp` (the 2026-10-03 server fix). |
| `sra/sra_rand50_07__switch1` | 2.8 / 3.1 | 15 | One switch (E s): a run entered by switch that ends at the loss budget (B), the silent Lw end of its source, a superseded (Ls) end. |
| `sra/sra_rand50_00__annotate_merge` | 4.9 / 10.0 | 17 | Annotate with merge: P runs, merges with empty partitions, the §5.1 union rule. |
| `sra/sra_rand50_00__cut_lists` | 5.2 / 12.1 | 17 | Cut lists (`max_labels_per_node` 2), `label_evidence: lower_bound`; 94 of 98 sampled nodes are cut, and each total is confirmed. |
| `sra/sra_rand50_01__cut_events` | 3.4 / 6.5 | 15 | Cut evidence: `max_branch_events` 1, `branch_diagnostics: cut`. |
| `sra/sra_16s_OK510341__time_tight` | 3.1 / 4.6 | 14 | A partial, time-budgeted result: `time_budget_ms` 200, both arms stopped (M), walks partial. |
| `sra/sra_batch3__batch_refuse` | 3.1 / 3.9 | 15 | A 3-seed batch with `max_seed_labels` 3: the middle seed is refused per seed, the others are walked. |
| `sra/sra_hairpin__left_only` (derived) | 16.1 / 37.1 | 15 | A skipped hairpin and its 152 `dead_end/hairpin` ends, on a 50-mer 8 bp before the first hairpin of the 2.6 MB body `sra_16s_PZ326290__beam20`. The region is a poly-A/poly-T run. |
| `sra/sra_hairpin__left_only_follow` (derived) | 16.2 / 36.2 | 15 | The same seed with `hairpins: follow`: the followed hairpin, and the reverse-complement retrace ending `edge_reuse_rc` (V). |
| `uhgg/uhgg_conserved__limit2` | 16.3 / 23.1 | 16 | A window carried by 7 genomes, branch limit 2 with merge: merges, merge-closed runs, routes through non-first parents, 21 split clones carrying `route_bp`. |
| `uhgg/uhgg_contig_end__annotate_merge` | 3.4 / 5.3 | 17 | Annotate with merge on the last 100 bp of a contig. |
| `uhgg/uhgg_contig_end__switch1_limit1_keep` (derived) | 11.6 / 15.2 | 14 | The contig end with one switch allowed (target MGYG-HGUT-01052 as `labels.extra`), branch limit 1, keep: 13 switches, into that genome past the end. |
| `mini_refseq/mini_ndm1__trace` | 22.3 / 28.5 | 14 | Accession header labels under `support: trace`, keep; verifiable identity; edge-reuse (U) ends. |
| `mini_refseq/mini_ndm1_rc__switch1` | 8.7 / 12.0 | 16 | blaNDM-1 reverse complement on a basic graph (the other strand): switches between header labels, routes through merges. |
| `mini_refseq/mini_ndm1__annotate_exh` | 4.4 / 7.0 | 17 | Annotate exhaustive on header labels, checked against `/resolve` discover of kind header. |
| `mini_refseq/mini_rep_00__trace_coords` | 8.3 / 9.7 | 16 | Record coordinates of header labels (kind `record`) on a repeat window, "unlimited"; `f_positions` from the recorded source-record answers (digests of the bases, string counts). |
| `mini_refseq/mini_rep_00__trace_coords_column` | 3.2 / 3.4 | 16 | Record coordinates of column labels (kind `column`): the column's k-mer index space, `trace_record_boundaries`. |

"Derived" cells are not in the cache matrix. The cache had no small hairpin (the only
one is in a 2.6 MB body), no `edge_reuse_rc` end anywhere, and no UHGG cell where a
switch actually happens. `snapshot` fetches these three live from seeds it derives, and
stages them under `<scratch>/offline_derived/`.

The 13 recorded skips are data-driven:

- 10 cells have no merge-derived support, so `c_routes` skips;
- the contig-end keep cell has no `label_lost` end inside a segment, so `b_lost` skips;
- the time-budgeted cell skips `to_json` and `invariants`, because its two fetches
  stopped apart; `full_invariants` runs instead.

```sh
python3 -m unittest discover -s api/python/tests -p 'test_traverse_real_offline.py'
cd api/python/tests/real
python3 offline.py check        # every recorded check, offline, as a table
python3 offline.py situations   # what the bodies cover
python3 offline.py snapshot     # rebuild: needs the cache; servers only for derived cells
                                # not staged yet, and for /resolve answers not cached yet
```

`snapshot` copies the cells and records the answers every oracle plan asks for. It then
runs every check offline on the new set and records what ran. If any check fails, the
old fixtures are kept. When an intended library change alters an oracle plan, rebuild
while the servers run and review the diff of `manifest.json`.
