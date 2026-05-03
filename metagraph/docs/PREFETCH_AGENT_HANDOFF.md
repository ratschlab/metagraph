# Agent handoff: mmap prefetch work (`mk/prefetch`)

This note is for whoever continues **MADV_WILLNEED prefetch** work and runs benchmarks on a server. Read this first, then open the PR linked below for the full checklist and results table.

## TL;DR

| Item | Value |
|------|--------|
| **Branch** | `mk/prefetch` (starts from `origin/master`) |
| **PR** | [#628](https://github.com/ratschlab/metagraph/pull/628) — implementation plan + phase checkboxes + **Results** table |
| **Current branch contents** | Benchmark script only (`scripts/bench_query_prefetch.py`) |
| **Next code step** | Cherry-pick **`94d7caa60`** from `mk/madvise-suffix-ranges` → Phase 1 (suffix_ranges `MADV_WILLNEED` + helpers) |

## What you should do first (server)

1. **Checkout and build**
   ```bash
   git fetch origin && git checkout mk/prefetch && git pull
   # build metagraph as usual (e.g. cmake + make in your build dir)
   ```

2. **Baseline benchmark (Phase 0)** — codebase is effectively **master + bench script** only; no prefetch yet.
   ```bash
   metagraph/scripts/bench_query_prefetch.py \
     --metagraph /path/to/your/build/metagraph \
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

3. **Record results** in PR #628 → **Results** table, row “0 (master baseline)”.

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

### Flags that matter for prefetch A/B

- **`--mmap`** — graph/annotation loaded via mmap (prefetch targets exist).
- **`--madv-random`** — enables madvise hints (`utils::with_madvise()`); Phase 1’s `MADV_WILLNEED` is gated the same way.

Forward extra server flags with repeated `--server-arg`, e.g.:

```bash
--server-arg --query-batch-size --server-arg 100000
```

## Implementation roadmap (after baseline)

The **authoritative** breakdown is in **PR #628 description** (checkboxes + expected deltas per phase). Short map:

| Phase | Content | Source |
|-------|-----------|--------|
| **0** | Benchmark script | Already on `mk/prefetch` |
| **1** | `suffix_ranges` + `madvise_willneed` + `get_mmap_data` + CanonicalDBG unwrap | Cherry-pick **`94d7caa60`** from `mk/madvise-suffix-ranges` |
| **2** | RowDiff `anchor_` / `fork_succ_` prefetch | Implement per PR |
| **3** | Bloom filter prefetch | Implement per PR |
| **4** | BRWT adaptive node prefetch | Implement per PR |
| **5** | RowDisk boundary prefetch | Implement per PR |
| **6** | Optional `valid_edges_` | Data-driven |
| **7** | Virtual `prefetch()` cleanup | After several phases land |

After **each** phase: rebuild, re-run the benchmark with the **same** arguments, save a new `--json`, update the PR **Results** table and tick the checkbox.

### Cherry-pick Phase 1

```bash
git checkout mk/prefetch
git fetch origin mk/madvise-suffix-ranges   # if needed
git cherry-pick 94d7caa60
# resolve conflicts if any, build, benchmark → fill Results row “1”
```

If `mk/madvise-suffix-ranges` is deleted later, use `git cherry-pick 94d7caa60` by SHA from `git log`.

## Gotchas

1. **Trace level:** The script always passes **`-v`** to `server_query` so trace lines appear.
2. **Per-run attribution:** The script slices server log by buffer offsets between timed requests. **Overlapping** concurrent clients can blur per-run stats; single-client `--parallel 1` on the script side is safest.
3. **`--drop-cache`:** Without root, the script warns and skips cache drop; benchmark still runs.
4. **macOS vs Linux build dirs:** User previously used `metagraph/build`; adjust `--metagraph` to your binary path.

## Related code pointers

- Server multi-graph load: `src/cli/server.cpp` (CSV parsing, `graphs_cache`)
- Search endpoint: `POST /search` in same file
- Phase 1 touchpoints (after cherry-pick): `src/common/utils/file_utils.{hpp,cpp}`, `src/graph/representation/succinct/dbg_succinct.{hpp,cpp}`, `src/cli/query.cpp`

## Questions?

If PR #628 body and this file disagree, **prefer PR #628** for checklist/results — update this doc when the workflow changes.
