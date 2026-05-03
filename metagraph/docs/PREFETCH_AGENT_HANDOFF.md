# Agent handoff: mmap prefetch work (`mk/prefetch`)

This note is for whoever continues **MADV_WILLNEED prefetch** work and runs benchmarks on a server. Read this first, then open the PR linked below for the full checklist and **Results** table.

## TL;DR

| Item | Value |
|------|--------|
| **Branch** | `mk/prefetch` (starts from `origin/master`) |
| **PR** | [#628](https://github.com/ratschlab/metagraph/pull/628) — implementation plan + phase checkboxes + **Results** table |
| **Code landed (local)** | Phases **0–5** on `mk/prefetch` (see commits below). Phases **6–7** still open per PR. |
| **Next code step** | **Phase 6** (optional `valid_edges_`, data-driven) **or** **Phase 7** (virtual `prefetch()` cleanup on `DeBruijnGraph` / `BinaryMatrix` after more phases merge). |

### Commits on `mk/prefetch` (prefetch stack)

| Phase | Commit | Summary |
|-------|--------|---------|
| **0** | `532ea6de` | `scripts/bench_query_prefetch.py` |
| **0** (doc) | `701d67c7` | This handoff + pointers |
| **1** | `f4cdb48c` | `suffix_ranges` + `utils::madvise_willneed` / `get_mmap_data` + `query_fasta` (CanonicalDBG unwrap) |
| **2** | `9443f9de` | RowDiff `.anchors` / `.rd_succ` mmap + `IRowDiff::prefetch()` + `query_fasta` |
| **3–5** | `e2a8c656` | Bloom `.bloom` mmap + `DBGSuccinct::prefetch_bloom_filter`; BRWT `prefetch_if_dense`; RowDisk / IntRowDisk / CoordRowDisk boundary `MADV_WILLNEED` at batch entry |

## What you should do first (server)

1. **Checkout and build**
   ```bash
   git fetch origin && git checkout mk/prefetch && git pull
   # build metagraph as usual (e.g. cmake + make in your build dir)
   ```
   Use the DNA binary for benchmarks, e.g. `.../build/metagraph_DNA`.

2. **Baseline benchmark (Phase 0 row in PR table)** — compare against **`701d67c7`** (script + docs, **no** C++ prefetch) **or** rebuild that commit for a fair binary match.
   ```bash
   metagraph/scripts/bench_query_prefetch.py \
     --metagraph /path/to/your/build/metagraph_DNA \
     --graphs graphs.csv \
     --query reads.fa \
     --mmap --madv-random \
     -p 4 --threads-each 8 \
     --warmup 1 --repeats 3 \
     --json bench/phase0_baseline.json
   ```
   Optional **cold OS page cache** between repeats (needs root):
   ```bash
   sudo metagraph/scripts/bench_query_prefetch.py ... --drop-cache --json bench/phase0_cold.json
   ```
   **Note:** `sudo` is **only** for `--drop-cache`. Normal runs do not need root.

3. **Record results** in PR #628 → **Results** table (one row per phase you benchmark).

## Benchmark script

- **Path:** `metagraph/scripts/bench_query_prefetch.py`
- **Behavior:** Spawns `metagraph server_query`, POSTs `/search` with FASTA JSON (same shape as `api/python/metagraph/client.py`). Parses `-v` trace lines:
  - **K-mer mapping:** `[Query graph construction] Contigs mapped to the full graph ... in X sec`
  - **Row-diff annotation:** `RD query [...] traversal: X sec, call_rd_rows: Y sec, ...`
- **Output:** Human-readable summary + `--json` for diffing runs.
- **Docs:** Module docstring at top of the script.

### Multi-graph CSV (`server_query`)

One row per index:

```
<name>,<graph_base_path>,<annotation_base_path>
```

Spaces around commas break parsing — keep rows tight. See `server_query` help in `src/cli/config/config.cpp`.

**`graphs` JSON filter:** If the CSV has **more than 10 distinct names**, the server rejects requests **without** a `"graphs"` list. The script’s repeated **`--graph-name`** flags supply that list (one flag per distinct name you want to query).

### Flags that matter for prefetch A/B

- **`--mmap`** — graph/annotation loaded via mmap (prefetch targets exist).
- **`--madv-random`** — enables madvise hints (`utils::with_madvise()`); all `MADV_WILLNEED` paths are gated the same way.

Forward extra server flags with repeated `--server-arg`, e.g.:

```bash
--server-arg --query-batch-size --server-arg 100000
```

## Implementation roadmap (sync with PR #628)

The **authoritative** checklist is in **PR #628**. Status on this branch:

| Phase | Content | Status on `mk/prefetch` |
|-------|---------|-------------------------|
| **0** | Benchmark script | Done (`532ea6de`) |
| **1** | `suffix_ranges` + `madvise_willneed` + `get_mmap_data` + CanonicalDBG unwrap | Done (`f4cdb48c`; originally cherry-picked from `94d7caa60` / `mk/madvise-suffix-ranges`) |
| **2** | RowDiff `anchor_` / `fork_succ_` prefetch | Done (`9443f9de`) |
| **3** | Bloom filter prefetch | Done (`e2a8c656`) |
| **4** | BRWT adaptive `prefetch_if_dense` on `nonzero_rows_` | Done (`e2a8c656`) |
| **5** | RowDisk / IntRowDisk / CoordRowDisk `boundary_` prefetch | Done (`e2a8c656`) |
| **6** | Optional `valid_edges_` | **Not started** (data-driven) |
| **7** | Virtual `prefetch()` cleanup | **Not started** (after several phases land) |

After **each** phase you benchmark: rebuild at that commit, re-run the script with the **same** arguments, save a new `--json`, update the PR **Results** table and tick the checkbox in PR #628.

## Gotchas

1. **Trace level:** The script always passes **`-v`** to `server_query` so trace lines appear.
2. **Per-run attribution:** The script slices server log by buffer offsets between timed requests. **Overlapping** concurrent clients can blur per-run stats; single-client `--parallel 1` on the script side is safest.
3. **`--drop-cache`:** Without root, the script warns and skips cache drop; benchmark still runs. For A/B on warm servers, skipping cache drop is often fine (see PR discussion).
4. **macOS vs Linux build dirs:** Adjust `--metagraph` to your `metagraph_DNA` path.
5. **Single-graph `server_query` + scripts:** The server may log **“Will listen”** before the async graph load finishes; readiness is **HTTP `GET /stats` returning 200** (not only a log grep). Multi-graph mode logs **“Ready to serve queries”** after load.
6. **`strace -e madvise`:** Shows `MADV_WILLNEED` when prefetches run; most wall time is still normal I/O/CPU outside `madvise`.

## Related code pointers

- Server multi-graph load: `src/cli/server.cpp` (CSV parsing, `graphs_cache`)
- Search endpoint: `POST /search` in same file
- **Phase 1:** `src/common/utils/file_utils.{hpp,cpp}`, `src/graph/representation/succinct/dbg_succinct.{hpp,cpp}`, `src/cli/query.cpp`
- **Phase 2:** `src/annotation/binary_matrix/row_diff/row_diff.{hpp,cpp}`, `src/cli/query.cpp`
- **Phases 3–5:** `src/kmer/kmer_bloom_filter.{hpp,cpp}`, `src/graph/representation/succinct/dbg_succinct.{hpp,cpp}`, `src/annotation/binary_matrix/multi_brwt/brwt.{hpp,cpp}`, `src/annotation/binary_matrix/row_disk/row_disk.{hpp,cpp}`, `src/annotation/int_matrix/row_disk/{int_row_disk,coord_row_disk}.{hpp,cpp}`, `src/cli/query.cpp`

## Questions?

If PR #628 body and this file disagree, **prefer PR #628** for checklist/results — update this doc when the workflow changes.
