# Review request: the traversal graphlet design (retrieve once, process locally)

> **Round 3 (design v3).** Your review of v2 is addressed in v3; §13 of the design lists the nine changes, each marked
> *(v3)* in place: views keep original ids and are saved as the unchanged backing body plus a view selector in `J`;
> an index identity (`index_ns` + content fingerprint `index_fp`, §3.1) gates every cross-retrieval join; one float
> algorithm (shortest scientific digits expanded positionally — 0 mismatches between C++ and Python on 20,000 raw-bit
> doubles plus edge cases) with raw-bit and differential golden vectors; evidence derived at the run's anchored
> endpoint and then clipped (`route_only` when the intersection is empty); tagged label selectors; a defined canonical
> form with whole-document fixtures; annotate `label_summary` pseudocode; cursors for walk chains and comparisons.
> Please check the v3 changes and say whether the freeze gate (§2.7) is now sufficient.

> **Round 2 (design v2).** Your review of v1 (`9fc93893`) is addressed in `docs/DESIGN-traverse-graphlet.md` v2; §12
> there is the list (17 items), each change is made in place and marked *(v2)*. In particular: `G.split` stores
> `Split::ambiguous`; `LabelRun` gains terminal `branches`/`loss` (`R` fields; unknown at an earlier cut); the codec
> is pinned (positional shortest round-trip floats, normative SETEXPR bases, UTF-8 byte prefixes at character
> boundaries) with shared golden vectors as the freeze gate (§2.7); `evidence_from` is derived per claim through merge
> partitions and `prev_run` (§5.1), never from `R.route_bp`; both continuation formulas are stated; labels are
> addressed by `{name, ref}` and ambiguous bare names are rejected; cursors, per-tool ceilings, `replay` and the
> oversize policy are defined (§6); entry handles are opaque and separate from body dedup; comparability uses the
> oriented seed sequence and index identity; `label_summary` is given as normative pseudocode matching
> `Walker::summarize()`. Please verify the v2 changes, look for what the corrections themselves broke, and say
> whether MGT v1 can be frozen once the golden vectors pass. The original v1 request follows unchanged.

**What I want from you:** an adversarial review of a *design*, before its implementation hardens. The document is
`docs/DESIGN-traverse-graphlet.md` on branch `gr/labeled-traversal` (head `e5273acc` or later; the design is
committed next to it). You know the traversal feature from rounds 1–3 (`docs/REVIEW-REQUEST-labeled-traversal.md`,
`docs/SPEC-labeled-traversal-core.md` §6–§7); this round is about the delivery format and the local library, not the
walker. Implementation of both tracks is starting in parallel; findings will be folded in as they arrive, so please
rank them by how expensive they become once code exists.

**The architecture under review.** An agent asks many questions about one retrieval. Every question is combinatorics
on a structure the backend has already materialised (`SeedResult`); only the first retrieval needs the index. So the
C++ server returns, once, a lossless compact text per seed ("MGT v1", design §2) plus a small JSON summary (§3), and a
Python library (`api/python/metagraph/traverse/`, §5) parses it into an in-memory graphlet and answers everything
else locally; an MCP server exposes those answers over handles (§6). The library also builds the *next* backend
request (a continuation is a new traversal), never answers it.

**Priorities, in order:**

1. **Losslessness and the `*` convention (§2.2–§2.5).** The format claims `dump(parse(x)) == x` and
   `Graphlet.to_json() == today's detail: full results[i]` after normalisation. Attack both. Is anything in
   `SeedResult` (walker.hpp) not representable? Is any field marked derivable actually ambiguous (two different
   results with the same dump)? The design found one such case already — clone run ids at splits — and paid a walker
   field for it (§4); find the next one. Are the `*` rules (§2.3) correct as stated, in particular `G.end *` in
   constrain mode ("entry ∪ switch targets − runs ending inside") and the annotate root's cut boundary list?
2. **Evidence semantics surviving the round trip (§2.4, §5 "Evidence semantics").** `route_bp` is stored on the run
   (earliest stamp) and per leaf label (the entry's, latest stamp) — the design says they differ on a run that passed
   two merges via non-first parents. Is the local `claims()` (claims not leaves; `evidence_from = max(from_bp,
   route_bp)`) sound under merging? Can any local view name a label against bases it does not support? Is the
   orientation rule (walking order on both arms; natural = reverse of the left concatenation; seed coordinate
   `−(i+1)`) stated so that a reader cannot get it wrong — it was got wrong once in round 2?
3. **The silent run end (§0.4).** A switch source's run ends `LABEL_LOST` with no `label_end` event; the format
   carries `Lw`. Are there other walker outputs that exist only in `growth[]` counters or in `Event::text` and would
   be lost, or double-counted, by deriving events from runs?
4. **Local operations (§5).** `walks()` ranking, `support_profile()`/`support_changes()`, `compare()` restricted to
   `min(complete_to_bp)` with `completeness_scope`, `subgraph()` inheriting `complete_to_bp` — which of these would
   return something a careful reader of the full JSON would not have concluded? Which need the backend after all?
5. **Size and memory (§7).** The estimates are per-record byte costs, not measurements. Where is the design's
   arithmetic weakest? The label-free case has 2,836 label names that are long file paths and up to 827 labels at a
   node — is front coding + delta-coded set expressions enough, or should the format carry a per-arm label alias
   table? What is the realistic Python footprint of a 2,800-segment arm with interned `array('I')` sets?
6. **The MCP protocol (§6).** Handles with TTL, LRU and a disk spool; every answer ≤ 2 KB with `more: N`; names not
   ids; `evidence` block on every local answer. What does an iterating agent still lack? Is anything proposed as
   local that silently needs the index (continuation spelling? label resolution? trace coordinates?)?
7. **Compatibility and cost (§8, §10, §11).** `detail: graphlet` is opt-in; additive JSON fields; one walker field.
   Is the two-track plan realistic, and is the merge point (the T37 conformance test on CLI-generated fixtures) the
   right gate?

**Specific questions.**

- §2.1 chooses a single space as the field separator because the text is embedded in a JSON string (a tab would
  escape to two bytes). Names are the last field and percent-encode only `%`, LF, CR. Is there any name an index can
  contain that breaks this (a leading/trailing space, a name equal to `*` or `.`)?
- §2.2 `R`: run id = ordinal in `ArmResult::runs` order; `end_labels[].run`, `label_summary.runs`, `prev_run` refer to
  it. Is `ArmResult::runs` order deterministic under every `frontier.order`? (Spec §6.8 claims byte-identical
  results across runs; the design relies on it.)
- §2.3 `G.entry *` for merges = union of the partition. `labels_via_parent` lists labels whose *kept entry* came
  through each parent; is the union of those exactly `labels_start` of the merged segment, including labels present
  on several parents?
- §2.5 derives `splits[]` from children grouped by parent and claims a split always has ≥ 2 children; it derives
  `kind = ambiguous` from "a label in ≥ 2 children's entry sets". Is that the walker's definition (`Split::ambiguous`
  is set from the lineage, walker.cpp `process_item`), or can they differ (e.g. under `hairpins: follow`, whose split
  "never counts toward ambiguity")?
- §3 drops `label_summary` from the JSON summary and recomputes it locally. The constrain-mode computation includes
  the rule "a merge does not clamp `direct_bp`"; the annotate-mode one is the parents-first union pass. Are both fully
  specified by the design, or does the library have to read the C++ to reproduce them?
- §4 `LabelRun::segment`: set at `end_run`, at the silent switch-source end, and in `merge_level` for runs closed by a
  merge. Is that every place a run ends or is closed? (Check `commit_entries` clones, `finish_path`, `censor_item`,
  `stop_arm`/`stop_frontier`.)
- §5 `compare(mode='claims')` keys labels by name plus `(kind, column, seq_id)`. Label ids differ between a constrain
  and an annotate retrieval (seed order vs first-seen). Is the join well defined when a header name occurs in two
  columns (the case rounds 2–3 fought over)?
- §6 `traverse_fetch(..., max_graphlet_mb=8)`: what should happen when a dump exceeds it — refuse, truncate with a
  marker, or spool to disk and return a summary only?
- §9 T37 normalises "events sorted within equal `at_bp`" and "`needed_budgets` sorted". Is anything else in today's
  JSON order-dependent but not information (so the conformance test would be flaky without normalising it)?

**Known limitations, do not re-report:** the walker itself and its claims are rounds 1–3; `switch_on: exhausted` is
not part of this design; the MCP server code lives in another repository and is not in scope beyond the tool list;
size numbers in §7 are estimates until `graphlet_measure.py` runs on the SRA server; GFA export uses the k−1 overlap
convention and loses label runs (documented).

**How to report.** For each finding: severity (blocker / major / minor), the design section, a concrete failure
(an input or a `SeedResult` shape → the wrong dump, the wrong local answer, or an ambiguity), and the change to the
design you propose. Separate *confirmed* (you can exhibit the case from the walker's code), *suspected* (needs a
run), and *design disagreement*. If the honest verdict is "the format is sound, the plan is fine, measure before
tuning", say so.
