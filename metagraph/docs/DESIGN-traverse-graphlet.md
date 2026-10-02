# Design: the traversal graphlet — retrieve once, process locally

**Status:** draft v1 (2026-10-02), for review before the implementation hardens. Companion documents:
`SPEC-labeled-traversal-core.md` (the traversal contract this builds on; §7.5 will carry §2 of this document once
accepted), `DESIGN-labeled-traversal-endpoint.md` (the owner's design note), `REVIEW-REQUEST-graphlet-design.md`
(the review prompt for this design).

## Purpose and architecture

An agent exploring a graph locus asks many questions about one retrieval: list the walks, show one sample's route,
where does the support change, compare constrained with label-free, give me FASTA. Every one of them is combinatorics
on a structure of tens of kilobases and a few thousand segments — nothing in them needs the index. The expensive
part, the annotation I/O, happens once. So:

```
agent  <-- small JSON views -->  MCP server (Python)  --- imports --->  metagraph.traverse (this repo, api/python)
                                       |   handles, LRU, disk spool; every follow-up answered locally
                                       v
                              POST /traverse  {strategy.output.detail: "graphlet"}      (one request, one body)
                                       v
                        MetaGraph query server (C++): a lossless compact text per seed + a small JSON summary
```

The backend owns graph work, the label oracle, caps, determinism, the completeness certificate (`complete_to_bp`,
`completeness_scope`, `walk_rule`, `cap_trigger`), the lossless dump and the small diagnostics. The local library owns
parsing, every derivation, ranking and pagination, spelling and orientation, FASTA/GFA, comparisons, sub-tries, handle
lifecycle, size caps on what the agent sees, and building the *next* backend request (a continuation is a new
traversal — the library prepares it, the backend runs it). Nothing in the local layer reads the graph.

Why a custom text and not more JSON views: the measured 3.72 MB compact JSON (0.18 MB gzipped) of a 2,500-leaf
exhaustive trie is ~95 % key names, repeated id lists and structure derivable from a far smaller core; a line-based
text that stores only the primary data is ~0.3–0.5 MB, parses in tens of milliseconds, is self-contained (seed
sequence, provenance, orientation rule in the file), can be `grep`ped, and is written by a streaming writer rather
than a `Json::Value` tree (the tree's deep copies were the 35–40 s tail on 50 MB responses). GFA cannot carry
per-node label runs, label ends with reasons or `route_bp`; it is a local export instead.

Decisions taken in this draft (open to the review): the text is embedded as a JSON string under `detail: graphlet`
(keeps transport, routing, CLI and error path untouched; costs one `\n` escape per line); runs are primary and label
ends derived (halves the biggest block) at the price of one walker field (`LabelRun::segment`, §4); orientation is
stored in walking order on both arms; every derivable field may be written as `*` only when the writer has verified
the rule reproduces the walker's value, so each production dump tests the derivations.

The remainder is the synthesized design (survey → three independent proposals → judged synthesis), lightly edited.

---

All paths relative to `metagraph/`. Line numbers refer to HEAD `373df0a2` and later; treat them as pointers, not anchors.

# 0. Verdict, and four facts that shape the synthesis

Yes to the architecture. `SeedResult` (walker.hpp:505-534) is fully materialised, so the retrieval format is a serializer in `src/cli/traverse.cpp`, not a walker change, with one tiny exception (§4). What I took from each proposal and what I corrected:

1. **Transport is already done.** `process_request(..., compact=true)` (server_utils.cpp:101-127, server.cpp:579/598/626) writes compact JSON and compresses gzip-first. All size targets below are against the *compact* 3.72 MB / 4.55 MB baselines, not the historical 7.6 MB.
2. **Delivery: embed the text in the JSON (`output.detail: "graphlet"`), not a new content type.** Taken from proposals 2/3. It keeps `server_utils.cpp`, `server.cpp` routing, the CLI (`json.loads(stdout)` in `TestTraverseBase._traverse`, test_traverse.py:90-97) and the error path untouched, and gives multi-seed responses for free. Price: field separator must be a space (a tab would JSON-escape to two bytes), and one `\n` escape per line (~15 KB on the SRA case).
3. **Runs are primary, label-end events are derived** (proposal 1's trick, halving the biggest block), **but the anchor needs the walker field** (proposals 1 and 3). Proposal 2's "runs are fully derivable by a sweep" fails on run *ids*: clones created at splits (`commit_entries`, walker.cpp:1727-1735) are identical rows whose ids depend on frontier order, and `end_labels[].run`, `label_summary[].runs`, `prev_run` reference those ids. Reconstructing them means replaying the walker's frontier order under `lowest_loss_first`/`most_supported_first`. Not worth it; four lines in the walker are.
4. **New finding, missed by all three proposals: a run can end with no `label_end` event.** `end_run(arm, src, at, LABEL_LOST)` at walker.cpp:2451 ends a switch *source* silently ("the lineage continues only under other names"); only the `SWITCH` event on the target documents it, and `growth[].label_ends` counts it. So "every ended run ↔ one label_end event" is false in one direction; the `R` record carries a distinct end token (`Lw`) for it, and the anchor must be set there too (§4). Proposal 3's invariant `#label_end <= #runs` still holds.

Additional corrections: `LabelEnd.route_bp` is the *entry's* route_bp (walker.cpp:1496, set at the latest merge the entry came through), while `LabelRun.route_bp` is the *earliest* stamp (`stamp_route`, :1946-1950); they differ on a run that passed two merges via non-first parents, so proposal 1's "route_bp from the run" is wrong for `end_labels` and the leaf record carries it explicitly. Splits always have ≥ 2 children (walker.cpp:2562-2577: `nf == 1` continues the segment), so grouping children by parent reproduces `splits[]` exactly. `revisit` carries `length_bp` (walker.cpp:1882-1884) and `same_distance` (:1923), `hairpin` carries "followed" (:2420, :2790): all dropped by today's JSON, all kept here.

# 1. Requesting it

- `strategy.output.detail`: `summary | tree | full | graphlet` (enumeration at traverse.cpp:329-330). Nothing else is new. `output.sequences` keeps its meaning (`false` → `G` lines without bases, first base kept, continuation sequence explicit); `continuation_bp`, `max_branch_events`, `profile_bin_bp` apply as today.
- The normalized echo gains `output.detail` and `output.timing` (closes the "not echoed" gap; both keys parse, so the echo stays resubmittable verbatim — test_traverse.py:485-493 keeps passing).
- `capabilities` gains `"graphlet_format": 1` and `"detail_levels": ["summary","tree","full","graphlet"]` (in every response and `/traverse/capabilities`).
- CLI unchanged: `metagraph traverse` prints the JSON with the embedded strings. `Graphlet.save()` writes the standalone `.mgt` file (§2, the `J` line).
- Spec §5 placeholders `output.label_lists` / `output.label_runs` (rejected by `Strict` today) are superseded by the graphlet's delta coding: delete from §5.

# 2. The graphlet text format, MGT v1

## 2.1 Lexical rules

- One record per line, `\n`-terminated, fields separated by one space, first field a capital letter. The only free-text field is the **last** field of `L` and `X` (label names) and runs to end of line; in it `%`, LF and CR are percent-encoded (`%25 %0A %0D`), nothing else. Every other field is ASCII without spaces. Seed ids live in the JSON summary, not in the body.
- Integers decimal. Floats shortest round-trip (`std::to_chars` / Python `repr`), integral values without a fraction, `inf` allowed. Booleans `0|1`. `*` = absent / "derive by the rule of this field". `.` = empty list. Lists comma-separated; a list of lists uses `|` (never contains spaces).
- **RANGES**: ascending label ids with consecutive runs collapsed: `0-5,7,9-12`; `.` = empty.
- **SETEXPR** (against a *base set* defined per field): `RANGES` (explicit) | `!` (equal to the base) | `!REMOVED` | `!+ADDED` | `!REMOVED+ADDED` (base minus REMOVED plus ADDED, both RANGES). The encoder emits whichever of explicit/delta is shorter, tie → explicit. Canonical in both encoders, so `dump(parse(x)) == x` byte-exact.
- **Codes.** End reasons (traversal_types.cpp:37-55): `D` dead_end, `L` label_lost, `B` loss_budget, `R` branch, `U` edge_reuse, `V` edge_reuse_rc, `J` rejoined_seed, `T` trace_break, `X` max_extension_bp, `S` max_steps, `P` max_live_paths, `N` max_paths, `O` max_output_bp, `M` time_budget, `W` beam_pruned; a second letter keeps the walker's text qualifier (walker.cpp:2473-2480, quorum texts at the BRANCH ends): `Rm` minority, `Rb` below_min_labels, `Rs` split_limit, `Dh` hairpin, `Ls` superseded, `Lx` switch_sources. Enum and qualifier both survive (today's JSON replaces the enum by the text). Arm `l|r`; status `c|t|p`; scope `p|u` (per_path | united_history); mode `c|a`; support `k|t`; reconverge `m|k`; label kind `c|h`; hairpin `f|s`.

## 2.2 Records

Document order: `H [J] S X* L*`, then per requested arm, left before right: `A B* V*`, then per segment in id order `G P* E* T? C?`, then `R*`; finally `Z`.

```
H mgt 1 <k> <regime> <alphabet> <mode c|a> <support k|t> <reconverge m|k> <cap> <continuation_bp> <seed_index> walk
J <compact JSON>                       # standalone FILES only: the response envelope with results[] reduced to this
                                       # seed's summary (graphlet string removed). Written by the library, never by
                                       # the server. At most once, directly after H.
S <validated_seed_id> <length_bp> <num_kmers> <num_seed_labels> <SEQUENCE>
                                       # the seed as validated (upper case, request orientation): today it is echoed nowhere
X <reason> <a-b,c-d> <name>            # dropped seed label with its k-mer runs on the seed (SeedResult::dropped_labels)
L <c|h> <column> <seq_id|*> <prefix_len> <suffix>
                                       # label id = ordinal (constrain: seed labels then extra; annotate: first seen).
                                       # name = previous L name[:prefix_len] + suffix (front coding). column/seq_id =
                                       # LabelRef (the stable cross-retrieval join key). seq_id * for column labels.
A <l|r> <status c|t|p> <complete_to_bp> <scope p|u> <live_paths> <live_labels> <exact> <max_seen> <nodes_truncated>
  <steps,successor_enumerations,output_bp,pair_evaluations,edge_reuse_probes,reminimisation_rounds,max_reminimisation_rounds>
  <cap_trigger reason,at_bp,segment,live_paths,live_labels,exact | *> <branch_events_total>
  <n_segments> <n_runs> <n_leaves> <n_splits> <n_merges> <bases>          # the last six validate the body
B <from_bp> <max_live_paths> <distinct_live_labels> <live_pairs> <exact> <steps> <divergences> <ambiguous_branches>
  <splits> <reconvergences> <bubbles> <tips> <blocked_repeat> <label_ends code:n,... | .>     # one per growth bin
V <at_bp> <segment> <chars> <labels_per_successor csv> <ambiguous RANGES> <dropped RANGES> <refused char:cause:RANGES;... | .>
                                       # the capped branch events, as stored (the only record of a not-followed successor)
G <parents csv|*> <from_bp> <length_bp> <entry SETEXPR|*> <entry_total|*> <end SETEXPR|*> <partition RANGES|RANGES|...|*> <first_base|*> <bases|*>
                                       # segment id = ordinal within the arm (walker creates parents first, walker.cpp:1351);
                                       # bases in WALKING order; first_base only when bases are absent (sequences:false)
P <from_bp> <to_bp> <total> <SETEXPR>  # annotate: LabelSetRun; base = previous P of this segment, for the first P the entry set
E <at_bp> s <from> <to> <cost>         # switch
E <at_bp> b <char> <reason code> <total> <RANGES>      # blocked successor
E <at_bp> h <char> <total> <RANGES> <f|s>              # hairpin followed | skipped
E <at_bp> v <segment> <delta | =>      # revisit: distance delta, "=" = same_distance (keep mode)
E <at_bp> t <char> <length_bp>   /   E <at_bp> u <length_bp> <alleles>     # tip / bubble: reserved, never emitted today
T <path_reason code|*> <extras label:loss:branches:route_bp,... | .>
                                       # this segment is a LEAF; extras only for labels where any of the three is non-zero
C <n> <loss_used> <branches_used> <labels SETEXPR> [<sequence>]
                                       # continuation (leaves ended by X or a resource reason): the walker's label logic is
                                       # PRIMARY (it is subtle), the spelling is derived: the last n bases of the natural
                                       # seed+flank. <sequence> present only when G carries no bases.
R <segment> <label> <from_bp> <to_bp> <end> <route_bp> <from_label:cost|*> <prev_run|*> <structural_successors|*> [<needed_budget>]
                                       # run id = ordinal (ArmResult::runs order, so end_labels[].run / label_summary.runs /
                                       # prev_run keep their meaning). <segment> = anchor (§4). <end> = two-letter code (ended
                                       # with a label_end event on <segment> at to_bp) | Lw (ended LABEL_LOST silently: the
                                       # lineage switched away; the SWITCH event(s) at to_bp document it) | m (closed by a
                                       # merge at to_bp on <segment>, ended=false). needed_budget only for code B.
Z <line count of the document including this line>
```

## 2.3 The `*` convention (lossless by construction, rules exercised on every dump)

A field marked `|*` has a *rule*. The C++ writer computes the rule's value from the walker's own data and emits `*` **only if it equals the walker's value**, else the explicit value; the reader applies the rule on `*`. So the format never depends on a derivation being right, and every production dump is a test of it. Rules:

- `G.parents *` = root. `G.entry *` (constrain only): root → ids `0..num_seed_labels-1` (root.state = seed labels, walker.cpp:1196-1199); merge → union of the partition. In annotate mode `G.entry` is never `*`; its SETEXPR base is parent[0]'s end set (root: explicit RANGES). Split children in constrain mode: SETEXPR against parent[0]'s end set (which labels followed which branch leaves no event, walker.cpp:2581-2593, so this is primary).
- `G.entry_total *` = |entry| (the only non-`*` case is the annotate root's cut boundary list, walker.cpp:1196, and cut children).
- `G.end *`: constrain → `(entry ∪ {to : E s on this segment}) − {label(r) : R anchored here with to_bp < from_bp+length_bp}` (spec §7.1 "label sets change only through events"; runs ending exactly at the segment end are still in `labels_end`, walker.cpp:2580, 1487); annotate → last `P` set, or entry when the segment has no `P`.
- `G.partition *` = no merge (|parents| ≤ 1) or annotate mode (`labels_via_parent` empty).
- `T.path_reason *` = none (semantic end). `R.from_label:cost *` = entered by seed. `R.prev_run *` = none. `R.structural_successors *` for `Lw`/`m`.

## 2.4 Orientation (the one rule, answering REVIEW-REQUEST §5.1.7)

`walk` in `H`: every `G` stores its bases in walking order on **both** arms, so outward index `i ∈ [from_bp, from_bp+length_bp)` is `bases[i − from_bp]`, and every position field (`from_bp`, `to_bp`, `at_bp`, `route_bp`, `complete_to_bp`) indexes bases directly. Natural orientation: right flank = concat root→leaf; left flank = `reverse(concat root→leaf)` = concat leaf→root of per-segment reversals = `spell_path` (walker.cpp:3220-3232). Whole molecule, natural: `natural(left) + seed + natural(right)`. Seed coordinate of outward `i`: `|seed| + i` (right), `−(i+1)` (left). The serializer reverses `Segment::sequence` of the left arm back to walking order (finalize reversed it at walker.cpp:2877-2878).

## 2.5 Stored vs derived

Stored (primary): bases; DAG (`G`); entry sets where not implied; merge partitions (`labels_via_parent`, not in today's JSON); annotate presence runs (`P`); runs with anchor, end token, route_bp, from/cost, prev_run, structural_successors, needed_budget (`R`); leaf path_reason and per-label `(loss, branches, route_bp)` (`T`); continuation labels/loss/branches/n (`C`); switch/blocked/hairpin/revisit events (`E`); growth bins (`B`); capped branch events incl. refusals (`V`); seed sequence (`S`); dropped labels (`X`); dictionary with column/seq_id (`L`).

Derived by the library (and by `to_json()` to reproduce today's `results[i]`): `children`; leaves (= segments with `T`); `paths[]` (leaf ordinal in segment-id order, chain via `parents[0]`, `length_bp = from_bp + length_bp` — exactly finalize, walker.cpp:2882-2902, so path ids match the server's); `end_labels` = `R` anchored at the leaf with `to_bp == length_bp` (every alive label at the leaf has an ended run, walker.cpp:1488-1491; every run ending at the leaf's last node is a label alive there) + `T` extras; `end_reasons`, `n_labels`; `label_end` events from `R` (code ≠ `Lw`, `m`; at = to_bp); `reconverge` events from `G` with > 1 parent (at = from_bp, segments = parents, walker.cpp:2013-2018); `splits[]` (children with one parent grouped by parent, ordered by (at_bp, first child id); `kind` ambiguous ⇔ a label is in ≥ 2 children's entry sets, constrain only; `labels_before` = |parent end| (constrain) / parent's last `P` total or entry_total (annotate); branches: char = first base, `labels_distinct` = child entry_total, `labels` = entry cut to `cap`); `labels_at_end` = `G.end`; `label_sets` = `P`; `runs[]` = `R`; `label_summary` (constrain: walker.cpp:2906-2946 incl. "a merge does not clamp direct_bp"; annotate: the parents-first union/intersection pass :2949-3010); `needed_budgets` (B-runs; a histogram, order not information); `branch_events_truncated`; continuation sequence; `seed.labels`; every `*_truncated` flag (`total > |list|`).

Not information, normalised by the conformance test: event order among equal `at_bp` (stable sort by `at_bp` only, walker.cpp:2879), `needed_budgets` order, `growth[].label_ends` `null` vs `{}` (traverse.cpp:752 default-constructs; recommend `Json::objectValue`).

## 2.6 Worked excerpt (integration fixture, quorum-pruned right arm, `min_successor_labels: 2`, radius 130; acc2 is the minority at the seed boundary)

```
H mgt 1 31 basic $ACGT c k m 64 1000 0 walk
S 3e6218a89e928eba 120 90 3 TAGCTCAATTTATGTACTGATTGCTTAAAGG…CGGT
L h 0 1 0 acc1
L h 0 2 3 2
L h 0 3 3 3
A r c 130 u 0 0 1 3 0 120,121,120,0,0,0,0 * 1 1 3 1 0 0 120
B 0 1 3 3 1 100 0 0 0 0 0 0 0 R:1
B 100 1 2 2 1 20 0 0 0 0 0 0 0 D:2
V 0 0 GT 2,1 . 1 T:minority:1
G * 0 120 * * * * * GCTGATATATGTGGAATTGGTCCTACGCCAACGG…TGATA
T * .
R 0 0 0 120 D 0 * * 0
R 0 1 0 0 Rm 0 * * 2
R 0 2 0 120 D 0 * * 0
Z 15
```
Reading: root entry `*` = {0,1,2}; `G.end *` = {0,1,2} − {label 1, ends at 0 < 120} = {0,2} = `labels_at_end`; label 1's `Rm` run becomes the `label_end` event (at 0, minority, 2 structural successors); labels 0 and 2 are the leaf's `end_labels` (loss 0, branches 0, route 0); `V` keeps the refused `T` successor and its label. Left arm, radius 10 (both proposals' example): children `G 0 0 10 2 * * * * GTCAGGTGGT` and `G 0 0 10 !2 * * * * TCGTGGATTA`; natural spelling of the first is `TGGTGGACTG` = `LEFT2[-10:]`, which is what today's JSON `sequence` holds.

# 3. The JSON summary (per seed, `detail: graphlet`)

Envelope unchanged (`release`, `capabilities`(+2 keys), `strategy`(+`output.detail/timing`, `clamped`), `walk_rule`, `algorithm_version` stays `traverse-0.2`, `results[]`, `timing?`). Per seed:

```
{"seed": {"seed_id","validated_seed_id","seed_id_mismatch","length_bp","num_kmers","labels_from_seed",
          "labels_supporting_total","labels_dropped","labels_dropped_digest","num_labels","num_seed_labels"},
 "label_mode": "constrain|annotate", "duplicate"?: true,
 "arms": {"left"?: {"status","complete_to_bp","completeness_scope","frontier_remaining":{…},"labels_per_node":{…},
                    "cap_trigger"?:{…},"counters":{7},"branch_events_total",
                    "counts": {"segments","leaves","splits","merges","runs","bases","max_bp",
                               "label_ends":{reason:n},"leaves_by_reason":{path_reason|"semantic":n}}},
          "right"?: {…}},
 "annotation": {"access_path","keys_mapped","rows_requested","direct_reads"}, "timing"?: {…},
 "graphlet": "<MGT text>", "graphlet_bytes": N, "graphlet_lines": N}
```
Derivation failures keep today's `{seed, error}` shape, no `graphlet`. No names, no per-leaf rows, no `label_summary`, no `growth` (in `B`), no `label_dict` (in `L`). Size ≈ 2.5 KB envelope + ~0.7 KB per arm; a 64-seed batch stays ~100 KB of summary.

# 4. The one walker change (marked separately)

`src/graph/traversal/walker.hpp` `LabelRun`: `uint32_t segment = UINT32_MAX;` — the segment on which the run ended or was closed. Set at three sites, no behavioural effect:
1. `end_run()` (walker.cpp:1469) gains a `size_t segment` parameter; `end_label()` (:1479) passes `item.segment`; the silent switch-source end (:2451) passes `item.segment`.
2. `merge_level()`, in the `labels_via_parent` loop (:2000-2008): `else arm.result.runs[e.run].segment = item.segment;` for an entry whose run is not the kept one (that run was closed by this merge).

Justification: the run↔segment association exists only in `Item`/`Entry` state and is dropped from `SeedResult`. From the output alone it is ambiguous whenever two segments end the same label at the same depth (clones at splits are identical rows), it is absent for the silent ends of §0.4 (no event at all), and for ≥ 3-parent merges (proposal 3). Without it, `R` would need the events stored beside it (~+100 KB, the biggest block twice), or a sweep that reconstructs run ids — not reproducible under non-breadth-first frontier orders. Also exposed as `runs[].segment` in the JSON (additive) and asserted in `tests/graph/traversal/test_trie.cpp`: every ended run's anchor holds its `label_end` event at `to_bp` (or a `SWITCH` with `from == label` for silent ends), and every `m` run's anchor is a parent of the segment created at `to_bp`.

# 5. The local library: `api/python/metagraph/traverse/` (stdlib only; pandas lazy)

`setup.py`: `find_packages(include=['metagraph', 'metagraph.*'])` (today's `include=['metagraph']` excludes subpackages from a non-editable install). `metagraph/__init__.py` unchanged (does not import `client`, so `import metagraph.traverse` is pandas-free).

| module | contents |
|---|---|
| `_codec.py` | `encode_ranges/decode_ranges`, `encode_setexpr(ids, base)/decode_setexpr(tok, base)` (shorter-wins, tie explicit), `fmt_num/parse_num`, `pct_escape/unescape`, `REASON/QUAL` code tables, `LabelSetInterner` (content-hashed `array('I')`), `GraphletFormatError(line_no, msg)` |
| `model.py` | slotted dataclasses: `Label(id, kind, column, seq_id, name)`, `SeedInfo(…)`, `Segment(id, parents, from_bp, length_bp, walk, entry, entry_total, end, partition, first_base, presence: list[PresenceRun], events: list[Event], leaf: Leaf|None; derived: children, depth)`, `Event`, `PresenceRun(from_bp, to_bp, labels, total)`, `Leaf(path_reason, extras: dict[label,(loss,branches,route_bp)], continuation: Continuation|None)`, `Run(id, segment, label, from_bp, to_bp, end_code, reason, qualifier, silent, merged, route_bp, from_label, cost, prev_run, structural_successors, needed_budget)`, `Arm(side, status, complete_to_bp, scope, frontier, labels_per_node, cap_trigger, counters, growth, branch_events, segments, runs)`, `Graphlet(format, k, regime, alphabet, mode, support, reconverge, cap, continuation_bp, seed_index, seed: SeedInfo, labels, dropped, arms, summary: dict|None, derived_from)` |
| `parser.py` | `parse(text) -> Graphlet` (one pass; `A` counts and `Z` validated → raises on truncation), `dump(g) -> str` (canonical; `== text`), `Graphlet.from_response(result, response)`, `Graphlet.load(path)/save(path)` (the `J` line) |
| `derive.py` | cached per arm: `leaves, paths, splits, label_end_events, end_labels(leaf), labels_at_end, label_summary, needed_budgets, continuation_sequence(leaf), reconverge_events`; `check_rules(g)` = the §2.3 rules recomputed and compared with explicit fields (the oracle that any `*`-vs-explicit divergence is reported) |
| `ops.py` | the queries below |
| `export.py` | `to_json(g) -> dict` (today's `results[i]`, `detail: full`), `to_fasta`, `to_gfa` |
| `client.py` | `TraverseClient` |
| `store.py` | `GraphletStore` (handles, LRU, disk spool) for the MCP layer |
| `mcp_tools.py` | framework-agnostic tool functions returning dicts ≤ `max_bytes` (the MCP server, wherever it lives, registers them) |
| `frames.py` | `frames(g)` → dict of DataFrames, lazy `import pandas` |

Public API (one line each):

```python
parse(text) -> Graphlet                       # raises GraphletFormatError(line_no, msg); validates A counts and Z
Graphlet.dump() -> str                        # canonical text; equals the input for an unmodified graphlet
Graphlet.from_response(result: dict, response: dict) -> Graphlet   # results[i] + envelope → summary attached
Graphlet.load(path) / .save(path)             # the .mgt file: H, J (envelope with this seed only), body
Graphlet.to_json() -> dict                    # today's detail: full results[i] (natural orientation); the conformance oracle
Graphlet.summary() -> dict                    # agent-facing ≤ 2 KB: per-arm status/complete_to/counts, top labels, caveats
Graphlet.label(name_or_id) -> Label; .label_id(name) -> int
Graphlet.spell(arm, leaf, orientation='natural'|'walk', with_seed=False) -> str
Graphlet.walks(arm, *, top=None, by='support'|'length'|'loss', labels=None, route_consistent=True, min_bp=0) -> list[Walk]
    # Walk(path_id, leaf, segments, length_bp, sequence, path_reason, end_reasons, claims, labels_full, n_alive,
    #      complete, beyond_certified_bp); 'support' ranks by (-len(labels_full), -n_alive, -length_bp, path_id)
Graphlet.claims(arm=None, labels=None, at_most_bp=None, strict=True) -> list[Claim]
    # the §6.9 unit: one per RUN END (claims, not leaves). Claim(label, arm, path_id, segment, from_bp, to_bp, route_bp,
    # evidence_from=max(from_bp, route_bp), kind∈{end, alive, diverged, merged, stretch}, reason, qualifier, end_class∈
    # {lost, dead_end, radius, capped, blocked, pruned, open}, entered_by, from_label, cost, loss, branches, support, exact)
    # strict=True raises IncompleteRecording when annotate nodes_truncated > 0 (a cut list makes the oracle a lower bound)
Graphlet.label_walks(label, arm=None) -> list[LabelWalk]        # per run: (from, to, evidence_from, sequence, end, leaves_below, merged_into)
Graphlet.routes(label, arm) -> list[list[int]]                   # per-label routes through merges via G partitions
Graphlet.support_profile(arm, leaf, kind='displayed'|'route') -> list[SupportRun(from_bp, to_bp, labels, exact)]
Graphlet.support_changes(arm, leaf) -> list[Change(at_bp, added, removed, reasons)]   # where support changes and why
Graphlet.label_summary() -> dict                                 # direct_bp / reach_bp / reentries / runs, both modes
Graphlet.continuation(arm, leaf) -> Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed())
Graphlet.next_request(arm, leaves, bp=None, reduce_budget=True, **overrides) -> dict   # resubmittable /traverse request
Graphlet.subgraph(labels, arm=None, mode='any'|'all') -> Graphlet   # sub-trie, ids kept, inherits complete_to_bp
Graphlet.to_fasta(arm=None, leaves=None, with_seed=True, orientation='natural', width=None) -> str
Graphlet.to_gfa(with_seed=True) -> str          # S per segment on the seed strand, L with k-1 overlap, P per walk, LB/ER tags
Graphlet.compare(other, *, arm=None, labels=None, mode='claims'|'walks'|'labels'|'prefix_subset') -> Comparison
    # keyed by label NAME (+kind/column/seq_id); restricted to min(complete_to_bp); claims cut at D get end_class 'open';
    # Comparison(comparable, reason, depth_used, equal, only_in_a, only_in_b, notes) — notes carry the merge qualification
Graphlet.memory_bytes() -> int
TraverseClient(host, port, api_path=None, *, session=None, timeout=900, release=None)
    .capabilities() -> dict;  .resolve(sequence, *, labels=None, discover=None, select=None, **opts) -> dict
    .traverse(seeds, strategy=None, *, detail='graphlet', timing=True, graph=None) -> TraverseResponse(envelope, graphlets, errors)
    .traverse_raw(request) -> dict;  .deepen(graphlet, arm, leaves, bp=None, **overrides) -> TraverseResponse
    # requests only; explicit Accept-Encoding: gzip, deflate; non-2xx → TraverseError(status, message); 503 → ServerInitializing
GraphletStore(spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800, ttl_disk_s=7*86400)
    .put(response, request) -> handle ('g_'+sha256(text)[:12]); .get(handle) -> Entry (lazy parse, LRU); .free/.list/.save/.load/.sweep
```

Evidence semantics in the round trip: `route_bp` lives on the run (`R`) and per leaf label (`T`); `claims()` reports route support `[from_bp, to_bp)` and displayed-path support `[evidence_from, to_bp)` separately; "contiguous occurrence" is reported only under `support == 't'` (from `H`); per-label end reason = enum + qualifier; left-arm positions index `G` bases directly; completeness (`status`, `complete_to_bp`, `scope`, `cap_trigger`) gates every comparison; label-list cuts carry totals everywhere and `A.nodes_truncated`. Iterative deepening stays a backend call: the library derives the continuation and builds the request; it never answers deepening locally.

# 6. MCP tool surface (functions in `mcp_tools.py`; every return ≤ `max_bytes`, default 2048, `more: N` markers; ids never leave the server, agents see names; every local answer carries `evidence: {complete_to_bp, exact, support, reconverge}`)

Backend-calling: `traverse_capabilities(index)`; `traverse_resolve(index, sequence, labels?|discover?, select?, limit=10)`; `traverse_fetch(index, seed:{sequence, labels?, seed_id?}, strategy, keep=True, max_graphlet_mb=8) -> summary + handle` (the probe is this with tight bounds); `traverse_continue(handle, arm, walk, overrides?, execute=True) -> new summary + parent:{handle, arm, walk, overlap_bp}` (or the request when `execute=False`).

Local over a handle: `graphlet_summary(handle, arm?)`; `graphlet_walks(handle, arm, rank=support|length|loss, n=10, min_bp=0, label?, spell=none|tail|full, tail_bp=60)`; `graphlet_walk(handle, arm, walk)` (segment chain, labels in/out per segment, continuation); `graphlet_support(handle, arm, walk, step?)` (support runs with who left/joined and the reason: a `K`-coded end or "split: N labels took C"); `graphlet_labels(handle, arm?, rank=direct_bp|reach_bp, n=20, min_direct_bp=0, at_bp?, name?)` (table; `name` → that label's runs/routes incl. merge routes); `graphlet_splits(handle, arm, n=10, min_labels_before=2)`; `graphlet_claims(handle, arm, min_bp=0, route_consistent=True, labels?)`; `graphlet_sequence(handle, arm, walk|segment, from, to, orientation, with_seed)` (≤ 16 KB slice); `graphlet_export(handle, format=fasta|gfa|json|mgt, arm?, walks?, path) -> {path, bytes, records}`; `graphlet_compare(a, b, arm, mode) -> Comparison dict` (preconditions: same release, k, validated seed, else `comparable: false`); `graphlet_subtrie(handle, labels, arm?) -> derived handle`; `graphlet_list/free/save/load`. Unknown handle → "replay with traverse_fetch(replay=<handle>)" (the store keeps the request).

# 7. Size and memory (SRA case; estimates from per-record byte costs, to be re-measured by `scripts/traversal/graphlet_measure.py`)

Constrained left arm (2839 segments / 1437 leaves / 1402 splits / 4160 ends / 69 kb bases): `G` 2839 × ~30 B + 69 KB bases ≈ 155 KB; `R` ~4200 × ~24 B ≈ 100 KB; `T` 1437 × ~8 B ≈ 12 KB; `E`/`V`/`B`/`A`/`S`/`L` ≈ 8 KB → **≈ 275 KB ≈ 7.5 % of 3.72 MB** (target < 15 %), ~60-90 KB gzipped vs 180 KB. Rebuilt locally: paths 0.71 MB, splits 0.41 MB, labels_at_end, label_summary, 1437 continuation strings. Label-free (annotate): `G` ≈ 155 KB; `P` ~5000 runs delta-coded ≈ 150 KB; `L` 2836 file-path names front-coded ≈ 100-150 KB (0.42 MB raw: the intrinsic part); `T` 12 KB → **≈ 430-480 KB ≈ 10 % of 4.55 MB**, ~120-160 KB gzipped vs 330 KB. Pessimistic (names share no prefix, fragmented sets) ≈ 1 MB ≈ 22 %; `labels.max_labels_per_node` bounds it by construction, every cut keeps its true count. Default tuned walk: ~2.5 KB body vs 10 KB. Summary < 10 KB per seed (§3). Server: the body is one `std::string` appended straight from `SeedResult` (no `Json::Value` tree for it); peak transient = body + JSON copy + gzip ≈ 3 × body < 1.5 MB. Client: parse ~15k lines in 50-100 ms; model ~5 MB per constrained arm, ~20 MB label-free (interned `array('I')` sets).

# 8. Files to change

C++: `src/cli/traverse.hpp` (`detail` doc; `strategy_to_json(st, cost, detail, timing)`; `std::string graphlet_text(const SeedResult&, const Seed&, const Strategy&, Arm order…)`; `Json::Value seed_result_summary_json(…)`; `constexpr int kGraphletFormatVersion = 1`), `src/cli/traverse.cpp` (enumeration :329-330; echo :530-535; `GraphletWriter` in the anonymous namespace ~350 lines: codec helpers `ranges/encode_setexpr/fmt_num/pct_escape/reason_code`, `header/seed/dropped/labels/arm/growth/branch_events/segment/presence/events/leaf/continuation/runs/trailer`; rule computations for the `*` fields; `arm_to_json`/`seed_result_to_json` branches for `graphlet`; additive JSON fields `runs[].segment`, `segments[].labels_via_parent`, hairpin `followed`, revisit `length_bp`/`same_distance`, `label_dict[].column/seq_id`; `label_ends` → `Json::objectValue`; `capabilities_to_json` +2 keys; `process_traverse_request` passes `seed` and computes `duplicate` before serialising), `src/graph/traversal/walker.hpp/.cpp` (§4, ≤ 6 lines). Not changed: `server.cpp`, `server_utils.cpp` (an optional, separate 400/500 split per spec §10.3 is recommended but independent).

Python: `api/python/setup.py` (find_packages), new `api/python/metagraph/traverse/{__init__,_codec,model,parser,derive,ops,export,client,store,mcp_tools,frames}.py`, `api/python/tests/test_traverse_{codec,parse,derive,ops,client}.py`, `api/python/tests/data/traverse/*.{request.json,full.json,graphlet.json,mgt}`, `scripts/traversal/graphlet_fixtures.py` (regenerate + `--check`), `scripts/traversal/graphlet_measure.py`.

Tests/docs: `integration_tests/test_traverse.py` (new classes, §9), `tests/graph/traversal/test_trie.cpp` (anchor assert), `tests/graph/traversal/test_trie_cases.cpp` (a bubble+switch+hairpin case whose dump feeds the Python fixtures), `docs/SPEC-labeled-traversal-core.md` (§5 detail graphlet + complete echo, delete label_lists/label_runs; §7.1 additive fields, derivable fields noted, same-position events unordered, `label_ends: {}`; new §7.5 "The graphlet (MGT v1)" = §2 verbatim incl. the `*` rules and orientation; §10.3 capabilities keys, transport note; §11 T37-T40; §13 `LabelRun::segment`, silent switch-source ends, orientation decision), `docs/source/api.rst` ("Traversal: retrieve once, query locally"), REVIEW-REQUEST §5.1.7 resolved.

# 9. Test plan

C++ / integration (`integration_tests/test_traverse.py`, CLI and HTTP; `skipUnless metagraph.traverse importable`, the venv installs `api/python` editable):
- T37 round trip: for constrain tuned (merge), constrain exhaustive (keep), annotate exhaustive, annotate beam, quorum, left-only, both arms, caps-tripped `max_steps`, `max_labels_per_node: 1`, 3-seed batch with a derivation failure: run `detail: full` and `detail: graphlet` on the same request; `Graphlet.from_response(r, out).to_json()` equals the full `results[i]` after normalisation (events sorted within equal `at_bp`, `needed_budgets` sorted, `timing` removed); `parse(text).dump() == text`; `A` counts match; `graphlet_lines == Z`; `check_rules()` reports no divergence; `#label_end events ≤ #runs`; every `Lw` run has a `SWITCH` with `from == label` at `to_bp` on its anchor's chain.
- T38 orientation: `spell(left) + seed + spell(right)` reproduces `acc1/acc2/acc3`; left `continuation(leaf).sequence == paths[].continuation.sequence`; resubmitting a locally derived continuation returns 200 and extends.
- T39 oracle through the library: `constrain.compare(annotate, mode='claims').equal` reproduces test_traverse.py:954-978; `prefix_subset` on tuned-vs-exhaustive with `minority` as the omission reason (:621-672 semantics).
- T40 protocol: summary key set; `strategy.output.detail/timing` echoed and resubmittable; `capabilities.graphlet_format == 1`; `sequences: false` → `G` without bases, `first_base` set, `C` with sequence, `splits[].char` reproduced; determinism (two runs byte-equal); gzip transport unchanged (`test_api_traversal_routes_are_compact_and_compressible` stays); body < ½ of compact full JSON on the fixture.
- `test_trie.cpp`: the §4 anchor assert on every trie case (covers merges/switches the integration fixture lacks); the new `test_trie_cases.cpp` fixture (bubble under `merge`, constant-cost switch chain, hairpin) dumped via the CLI into the Python fixture set.

Python unit (`api/python/tests`): codec identities (random sets, shorter-wins, tie→explicit, floats, percent escapes, front coding); fixture parse → `to_json()` ≡ `full.json` (normalised) and `dump()` ≡ text, byte-exact; hand-built DAGs for `derive` (split, merge with partitions and `route_bp`, switch chain incl. a silent `Lw` end, censored leaf, zero-length merged leaf); `claims()` kinds and `evidence_from`; `compare()` equal / prefix_subset / incomparable; `subgraph(['acc2'])` keeps exactly RIGHT2's chain; FASTA/GFA parse back (S/L/P counts, k-1 overlaps, paths spell the walks); `GraphletStore` LRU/TTL with a fake clock; truncated body raises; `python -c "import sys; sys.modules['pandas']=None; import metagraph.traverse"` succeeds.

# 10. Implementation plan: two independent tracks against the agreed grammar

Step 0 (shared, half a day): freeze §2 as spec §7.5; adapt the scratchpad prototype `mgt_proto.py` to this grammar and write the five fixture dumps by hand-conversion from the existing `.json` dumps → the Python track's first fixtures exist before any C++ lands; agree the normalisation rules of §2.5.

Track A — C++ writer (one engineer, ~4-5 days):
1. §4 walker field + `test_trie.cpp` assert (half day; the only change outside `src/cli`).
2. Additive JSON fields, `label_ends: {}`, echo of `detail/timing`, capabilities keys; `detail: graphlet` enumeration (half day).
3. `GraphletWriter`: codec helpers with a gtest in `tests/cli/`; records in document order; the `*` rule computations with the equality check; summary JSON (2 days).
4. Integration tests T37-T40 using the Python parser from Track A's own branch of step 0's prototype if Track B is not merged yet (1 day).
5. Regenerate the Python fixtures from the CLI (`graphlet_fixtures.py`), run `graphlet_measure.py` on the SRA server, put the measured numbers in the spec (half day).

Track B — Python library (one engineer, ~5-6 days, needs only the grammar and the step-0 fixtures):
1. `_codec.py`, `model.py`, `parser.py` with `dump()` identity on the hand-made fixtures (1 day).
2. `derive.py` until `to_json()` matches the `.json` dumps (1.5 days; the place any derivation surprise surfaces).
3. `ops.py`, `export.py` (walks, claims, compare, subgraph, FASTA/GFA, continuation, next_request) (1.5 days).
4. `client.py`, `store.py`, `mcp_tools.py`, `frames.py`; unit tests (1 day).
5. When Track A's fixtures land: swap the hand-made fixtures for the CLI-generated ones, run the full suite, write `docs/source/api.rst` (half day).

Merge point: the T37 conformance test green on the CLI-generated fixtures; after that the MCP server (outside this repo) registers `mcp_tools.py`.

# 11. What changes for existing fields and tests

Nothing breaks: `detail: graphlet` is opt-in; `full/tree/summary` keep their shape (additive fields only; no test asserts exact key sets of `segments`, `events`, `runs`, `label_dict`; the exact-set assertions at test_traverse.py:455 `labels_per_node` and :889-890 `annotation` are untouched); the echo gains two accepted keys (:307-311, :485-493 keep passing); `capabilities` checks are `assertIn` (:711-723); `algorithm_version` stays `traverse-0.2` (:922) because the walk did not change — the format has its own version. `_walks`/`_structural_walks` (:106-153) remain right-arm JSON helpers; the new tests spell through the library. Pre-existing warts left alone: `setup.py` reads the version from `../../../package.json`; `tests/test_helpers.py:41` already fails on `df_from_align_result`.