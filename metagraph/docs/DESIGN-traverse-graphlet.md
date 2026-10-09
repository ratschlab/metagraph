# Design: the traversal graphlet — retrieve once, process locally

**Status:** v5.11 (2026-10-08; v5.11 = the owner's decision #17: the graph's mask and Bloom filter are derived
data outside the index identity, §29; v5.10 = the P3 server items of the review of 2026-10-06, folded into feature level 6
as corrections, §28; v5.9 = the P2 fixes of the review of 2026-10-06, folded into feature level 6, §27;
v5.8 = record coordinates as built, feature level 6, §26, after the review of levels 4–5, §25; v5.7 = the review of pass 5 and the path cache's first reads, §24; v5.6 = pass 5 as implemented, §23; v5.5 = the backend half of stage 4 as implemented, §22; v5.4 = stage 3 as implemented and the answers of its review, §20; v5.3 = the resource contract and library decisions as implemented after implementation reviews 3 and 4, §19; v5.2 = the owner's conservative outcome rule in §14) — **stages 1–3 and the backend half of stage 4 implemented. MGT v1 is FROZEN (2026-10-02, owner's decision after the implementation review: no finding required a format change). Every later change — fixes, stages 2–4 — stays within the v1 records, fields and tokens; a format change requires MGT v2. Spec §7.5 is the normative text where a design excerpt differs** (`3ecbfc47`…`5fbd9057`); the freeze criteria of the fifth review are met (golden vectors and round-trip fixtures pass, size measured on SRA, §7) — MGT v1 freezes on the owner's confirmation; **approved for implementation** by the fifth external review (no further architecture
review needed; MGT v1 freezes once the codec corrections and the round-trip fixtures pass; hard resource guarantees
are advertised only after the corresponding exhaustion and concurrency tests pass). Draft history: v5 (2026-10-02), revised after four external design reviews (of v1 `9fc93893`, v2 `23d109fc`,
v3 `91bda3e9`, v4 `e92f72cf`) and the owner's guarantee requirement; changes are listed in §12 (v2), §13 (v3), §15
(v4) and §16 (v5) and marked where they are made. §14 is the resource and guarantee contract; §14.1 separates what
is frozen with MGT v1 from what is enforced in stages. Not frozen: the fixtures of §2.7 and §14 are the freeze gate. Companion documents:
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
3. **Runs are primary, label-end events are derived** (proposal 1's trick, halving the biggest block), **but the anchor needs the walker field** (proposals 1 and 3). *(v2)* Two more run fields are needed for the same reason — the terminal `loss` and `branches` of a lineage that ends inside a walk exist only in the walker's discarded state (§4) — and `Split::ambiguous` must be stored, not derived (§2.5). Proposal 2's "runs are fully derivable by a sweep" fails on run *ids*: clones created at splits (`commit_entries`, walker.cpp:1727-1735) are identical rows whose ids depend on frontier order, and `end_labels[].run`, `label_summary[].runs`, `prev_run` reference those ids. Reconstructing them means replaying the walker's frontier order under `lowest_loss_first`/`most_supported_first`. Not worth it; four lines in the walker are.
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
- Integers decimal. *(v3)* **Floats, one normative algorithm:** take the *shortest round-trip significant digits* `d₁…dₙ` and decimal exponent of the double — the digits of C++ `std::to_chars(first, last, x, std::chars_format::scientific)` without a precision, equal to those of Python `repr(x)` (both produce the unique shortest digit string, nearest on ties) — and **expand them positionally**: zeros padded up to the decimal point, a `.` only when there is a fractional part, no trailing zeros, `-0` written `0`, `inf` for +∞ (NaN never occurs; all floats in the format are costs, budgets and losses ≥ 0). Python: `format(decimal.Decimal(repr(x)), 'f')` with `-0`/trailing-`.` normalisation. So `1e23` → `100000000000000000000000` and `1.2345678901234568e20` → `123456789012345680000`. (v1's `to_chars`/`repr` disagreed on `0.0001`; v2's `to_chars(..., fixed)` writes the exact binary value — `99999999999999991611392` for `1e23` — which round-trips but is not the same bytes.; a differential probe of this rule — 20,000 non-negative finite doubles drawn over raw bits plus the edge cases above, `DBL_MAX` and subnormals — gave 0 mismatches between the two languages.) Booleans `0|1`. `*` = absent / "derive by the rule of this field". `.` = empty list. Lists comma-separated; a list of lists uses `|` (never contains spaces).
- **RANGES**: ascending label ids with consecutive runs collapsed: `0-5,7,9-12`; `.` = empty.
- *(v2)* **Strings** are UTF-8. `L.prefix_len` counts **bytes** of the previous name's UTF-8 encoding, and the encoder shortens the shared prefix to the nearest UTF-8 character boundary, so the suffix is always valid UTF-8 on its own; a reader reconstructs `prev_bytes[:prefix_len] + suffix_bytes` and decodes. A name may contain spaces and may equal `*` or `.`: it is the remainder of the line after the fixed fields, preserved exactly (after percent-decoding `%25 %0A %0D`).
- **SETEXPR** (against a *base set* defined per field; *(v2)* the bases are normative: `G.entry` → parent[0]'s end set (merges: union of the parents' end sets; root: none — explicit RANGES only); `G.end` → this segment's entry set; `P` → the previous `P` of the segment, the first `P` → the entry set; `C.labels` → the leaf segment's end set): `RANGES` (explicit) | `!` (equal to the base) | `!REMOVED` | `!+ADDED` | `!REMOVED+ADDED` (base minus REMOVED plus ADDED, both RANGES). The encoder emits whichever of explicit/delta is shorter, tie → explicit. Canonical in both encoders, so `dump(parse(x)) == x` byte-exact.
- **Codes.** End reasons (traversal_types.cpp:37-55): `D` dead_end, `L` label_lost, `B` loss_budget, `R` branch, `U` edge_reuse, `V` edge_reuse_rc, `J` rejoined_seed, `T` trace_break, `X` max_extension_bp, `S` max_steps, `P` max_live_paths, `N` max_paths, `O` max_output_bp, `M` time_budget, `W` beam_pruned, *(v4)* `Y` resource_limit (memory or work budget, §14); a second letter keeps the walker's text qualifier (walker.cpp:2473-2480, quorum texts at the BRANCH ends): `Rm` minority, `Rb` below_min_labels, `Rs` split_limit, `Dh` hairpin, `Ls` superseded, `Lx` switch_sources. Enum and qualifier both survive (today's JSON replaces the enum by the text). Arm `l|r`; status `c|t|p`; scope `p|u` (per_path | united_history); mode `c|a`; support `k|t`; reconverge `m|k`; label kind `c|h`; hairpin `f|s`.

## 2.2 Records

Document order: `H [J] S X* L* O Q? K*` *(v5: O, Q, K)*, then per requested arm, left before right: `A B* V*`, then per segment in id order `G P* E* T? C?`, then `R*`; finally `Z`.

```
H mgt 1 <k> <regime> <alphabet> <mode c|a> <support k|t> <reconverge m|k> <cap> <continuation_bp> <seed_index> walk <index_ns|*> <index_fp|*> <index_meta_fp>
                                       # (v3) index identity: <index_ns> = the deployment's name for the index
                                       # ([A-Za-z0-9._-]+, new server flag --index-name; * = unset), <index_fp> = the
                                       # digest of the index bundle's manifest (v4, §3.1; * = no manifest: joins are
                                       # unverifiable), <index_meta_fp> = the metadata hash, a negative check only.
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
  <counters name=value,... >   # (v5) named and extensible, in the walker's order: steps, successor_enumerations,
                               # output_bp, pair_evaluations, edge_reuse_probes, reminimisation_rounds,
                               # max_reminimisation_rounds, refusal_scans, switch_sources_cut, ... — an unknown name
                               # is preserved by the reader, so a newer writer never breaks an older reader
  <cap_trigger reason,at_bp,segment,live_paths,live_labels,exact,demand | *> <branch_events_total>
  <evidence_complete_to_bp|*>   # (v5) branch-event boundary (* = no event dropped)
  <n_segments> <n_runs> <n_leaves> <n_splits> <n_merges> <bases>          # the last six validate the body
B <from_bp> <max_live_paths> <distinct_live_labels> <live_pairs> <exact> <steps> <divergences> <ambiguous_branches>
  <splits> <reconvergences> <bubbles> <tips> <blocked_repeat> <label_ends code:n,... | .>     # one per growth bin
V <at_bp> <segment> <chars> <labels_per_successor csv> <ambiguous RANGES> <dropped RANGES> <refused char:cause:RANGES;... | .>
                                       # the capped branch events, as stored (the only record of a not-followed successor)
G <parents csv|*> <from_bp> <length_bp> <entry SETEXPR|*> <entry_total|*> <end SETEXPR|*> <partition RANGES|RANGES|...|*> <split 0|1|*> <first_base|*> <bases|*>
                                       # segment id = ordinal within the arm (walker creates parents first, walker.cpp:1351);
                                       # bases in WALKING order; first_base only when bases are absent (sequences:false)
                                       # (v2) <split>: Split::ambiguous of the split at this segment's end (1|0), * = the
                                       # segment does not end in a split. PRIMARY: the walker counts predecessor LINEAGES
                                       # before quorum filtering and excludes followed hairpins, so it cannot be derived
                                       # from the children's label sets (disjoint child sets can be ambiguous after a
                                       # switch; overlapping ones can be a divergence).
                                       # (v2) annotate merges: <partition> is one '.' per parent (labels_via_parent empty)
P <from_bp> <to_bp> <total> <SETEXPR>  # annotate: LabelSetRun; base = previous P of this segment, for the first P the entry set
E <at_bp> s <from> <to> <cost>         # switch
E <at_bp> b <char> <reason code> <total> <RANGES>      # blocked successor
E <at_bp> h <char> <total> <RANGES> <f|s>              # hairpin followed | skipped
E <at_bp> v <segment> <delta | =>      # revisit: distance delta, "=" = same_distance (keep mode)
E <at_bp> t <char> <length_bp>   /   E <at_bp> u <length_bp> <alleles>     # tip / bubble: reserved, never emitted today
T <path_reason code|*> <extras label:loss:branches:route_bp,... | .>
                                       # this segment is a LEAF; extras only for labels where any of the three is non-zero
C <n> <loss_used> <branches_used> <labels SETEXPR against the leaf's end set> [<sequence>]
                                       # continuation (leaves ended by X or a resource reason): the walker's label logic is
                                       # PRIMARY (it is subtle), the spelling is derived (v2, both arms stated):
                                       #   right arm: the LAST n bases of  seed + natural(right flank)
                                       #   left arm:  the FIRST n bases of natural(left flank) + seed
                                       # (n <= |seed| + flank length; a continuation may cross into the seed).
                                       # <sequence> present only when G carries no bases.
R <segment> <label> <from_bp> <to_bp> <end> <route_bp> <from_label:cost|*> <prev_run|*> <structural_successors|*> <branches> <loss> [<needed_budget>]
                                       # run id = ordinal (ArmResult::runs order, so end_labels[].run / label_summary.runs /
                                       # prev_run keep their meaning). <segment> = anchor (§4). <end> = two-letter code (ended
                                       # with a label_end event on <segment> at to_bp) | Lw (ended LABEL_LOST silently: the
                                       # lineage switched away; the SWITCH event(s) at to_bp document it) | m (closed by a
                                       # merge at to_bp on <segment>, ended=false). needed_budget only for code B.
                                       # (v2) <branches> <loss>: the lineage's TERMINAL values when the run ended or was
                                       # closed (LabelRun::branches/loss, §4). Valid at to_bp only: a claim cut at an
                                       # earlier depth reports them as unknown, never as the value at the cut.
O <walks c|p|f> <branch_diagnostics c|x> <label_evidence c|l|q> <delivery i|s|p>   # (v5.2) q = qualified
                                       # (v5) the per-seed guarantee dimensions (§14), one per document, before the arms
Q <scope> <resource> <phase> <requested> <effective> <used> <remaining> <actions csv> <message>
                                       # (v5) resource_stop (§14), at most one per document, after O; message is the
                                       # free-text last field (percent-encoded like names)
K <arm l|r|*> <kind> <knob> <limit VALUE> <observed VALUE> <complete_to_bp|*> <extra name=VALUE,...|.> <effect>
                                       # (v5) one stated limitation; arm * = seed-level (seed_labels, server_clamp,
                                       # derivation). (v5.1) VALUE is typed by a one-letter prefix: i:<integer>,
                                       # f:<float, codec rule>, s:<string: every UTF-8 byte outside 0x21..0x7E, '%' and ','
                                       # as uppercase %XX; '=' stays raw — as implemented and in the golden vectors>, or
                                       # u (unlimited) — e.g. the scope limitation's limit is s:merge; <extra> holds
                                       # further fields (server_limit, demand). (v5.1) <effect> is STORED as the
                                       # free-text last field (percent-encoded like names): the C++ builds it from
                                       # values, and templates would have to substitute identically in two languages
Z <line count of the document including this line>
```

## 2.3 The `*` convention (lossless by construction, rules exercised on every dump)

A field marked `|*` has a *rule*. The C++ writer computes the rule's value from the walker's own data and emits `*` **only if it equals the walker's value**, else the explicit value; the reader applies the rule on `*`. So the format never depends on a derivation being right, and every production dump is a test of it. Rules:

- `G.parents *` = root. `G.entry *` (constrain only): root → ids `0..num_seed_labels-1` (root.state = seed labels, walker.cpp:1196-1199); merge → union of the partition. In annotate mode `G.entry` is never `*`; its SETEXPR base is parent[0]'s end set (root: explicit RANGES). Split children in constrain mode: SETEXPR against parent[0]'s end set (which labels followed which branch leaves no event, walker.cpp:2581-2593, so this is primary).
- `G.entry_total *` = |entry| (the only non-`*` case is the annotate root's cut boundary list, walker.cpp:1196, and cut children).
- `G.end *` *(v2, chronological)*: constrain → start from the entry set and apply, in increasing position (ties: ends before switch-ins), every run end anchored on this segment with `to_bp < from_bp+length_bp` (remove its label) and every switch event on it (add its `to`); the result is the end set. v1's set formula removed a label that re-entered after ending (A→B→A inside one segment). Any field may legitimately be written explicitly: the conformance test (T37) checks that every `*` the writer emitted reproduces the walker's value — not that `*` is always emitted. The v1 text read:
  `(entry ∪ {to : E s on this segment}) − {label(r) : R anchored here with to_bp < from_bp+length_bp}` (spec §7.1 "label sets change only through events"; runs ending exactly at the segment end are still in `labels_end`, walker.cpp:2580, 1487); annotate → last `P` set, or entry when the segment has no `P`.
- `G.partition *` = no merge (|parents| ≤ 1), or annotate mode, where `*` decodes to **one empty list per parent** (`labels_via_parent` is empty in annotate mode; *(v3)* stated so the reader reconstructs the right arity).
- *(v3)* **Canonical form.** The writer MUST emit `*` whenever the field's rule reproduces the walker's value, and the explicit value otherwise; every SETEXPR in its shortest form (tie → explicit); RANGES collapsed. `dump()` always writes this canonical form, so `dump(parse(x)) == x` holds for every canonical document — all writer output — while a valid non-canonical document (e.g. an explicit set where `*` would do) parses and dumps to its canonical form. `is_canonical(x)` is part of the library and of T37.
- `T.path_reason *` = none (semantic end). `R.from_label:cost *` = entered by seed. `R.prev_run *` = none. `R.structural_successors *` for `Lw`/`m`.

## 2.4 Orientation (the one rule, answering REVIEW-REQUEST §5.1.7)

`walk` in `H`: every `G` stores its bases in walking order on **both** arms, so outward index `i ∈ [from_bp, from_bp+length_bp)` is `bases[i − from_bp]`, and every position field (`from_bp`, `to_bp`, `at_bp`, `route_bp`, `complete_to_bp`) indexes bases directly. Natural orientation: right flank = concat root→leaf; left flank = `reverse(concat root→leaf)` = concat leaf→root of per-segment reversals = `spell_path` (tests/graph/traversal/walker_paths_for_tests.hpp). Whole molecule, natural: `natural(left) + seed + natural(right)`. Seed coordinate of outward `i`: `|seed| + i` (right), `−(i+1)` (left). The serializer reverses `Segment::sequence` of the left arm back to walking order (finalize reversed it at walker.cpp:2877-2878).

## 2.5 Stored vs derived

Stored (primary): bases; DAG (`G`); entry sets where not implied; merge partitions (`labels_via_parent`, not in today's JSON); annotate presence runs (`P`); runs with anchor, end token, route_bp, from/cost, prev_run, structural_successors, needed_budget (`R`); leaf path_reason and per-label `(loss, branches, route_bp)` (`T`); continuation labels/loss/branches/n (`C`); switch/blocked/hairpin/revisit events (`E`); growth bins (`B`); capped branch events incl. refusals (`V`); seed sequence (`S`); dropped labels (`X`); dictionary with column/seq_id (`L`).

Derived by the library (and by `to_json()` to reproduce today's `results[i]`): `children`; leaves (= segments with `T`); `paths[]` (leaf ordinal in segment-id order, chain via `parents[0]`, `length_bp = from_bp + length_bp` — exactly finalize, walker.cpp:2882-2902, so path ids match the server's); `end_labels` (*(v2)* emitted and derived in ascending label id) = `R` anchored at the leaf with `to_bp == length_bp` (every alive label at the leaf has an ended run, walker.cpp:1488-1491; every run ending at the leaf's last node is a label alive there) + `T` extras; `end_reasons`, `n_labels`; `label_end` events from `R` (code ≠ `Lw`, `m`; at = to_bp); `reconverge` events from `G` with > 1 parent (at = from_bp, segments = parents, walker.cpp:2013-2018); `splits[]` (children with one parent grouped by parent, ordered by (at_bp, first child id); `kind` ambiguous ⇔ the parent's `G.split` is 1 *(v2: stored, see §2.2; v1's "a label in ≥ 2 children's entry sets" was wrong in both directions)*; `labels_before` = |parent end| (constrain) / parent's last `P` total or entry_total (annotate); branches: char = first base, `labels_distinct` = child entry_total, `labels` = entry cut to `cap`); `labels_at_end` = `G.end`; `label_sets` = `P`; `runs[]` = `R`; `label_summary` (constrain: walker.cpp:2906-2946 incl. "a merge does not clamp direct_bp"; annotate: the parents-first union/intersection pass :2949-3010); `needed_budgets` (B-runs; a histogram, order not information); `branch_events_truncated`; continuation sequence; `seed.labels`; every `*_truncated` flag (`total > |list|`).

Not information, normalised by the conformance test: event order among equal `at_bp` (stable sort by `at_bp` only, walker.cpp:2879), `needed_budgets` order, `growth[].label_ends` `null` vs `{}` (traverse.cpp:752 default-constructs; recommend `Json::objectValue`).

## 2.6 Worked excerpt (integration fixture, quorum-pruned right arm, `min_successor_labels: 2`, radius 130; acc2 is the minority at the seed boundary)

```
H mgt 1 31 basic $ACGT c k m 64 1000 0 walk * 9f3c2a71d04be8e5
S 3e6218a89e928eba 120 90 3 TAGCTCAATTTATGTACTGATTGCTTAAAGG…CGGT
L h 0 1 0 acc1
L h 0 2 3 2
L h 0 3 3 3
A r c 130 u 0 0 1 3 0 120,121,120,0,0,0,0 * 1 1 3 1 0 0 120
B 0 1 3 3 1 100 0 0 0 0 0 0 0 R:1
B 100 1 2 2 1 20 0 0 0 0 0 0 0 D:2
V 0 0 GT 2,1 . 1 T:minority:1
G * 0 120 * * * * * * GCTGATATATGTGGAATTGGTCCTACGCCAACGG…TGATA
T * .
R 0 0 0 120 D 0 * * 0 0 0
R 0 1 0 0 Rm 0 * * 2 0 0
R 0 2 0 120 D 0 * * 0 0 0
Z 15
```
Reading: root entry `*` = {0,1,2}; `G.end *` = {0,1,2} − {label 1, ends at 0 < 120} = {0,2} = `labels_at_end`; label 1's `Rm` run becomes the `label_end` event (at 0, minority, 2 structural successors); labels 0 and 2 are the leaf's `end_labels` (loss 0, branches 0, route 0); `V` keeps the refused `T` successor and its label. Left arm, radius 10 (both proposals' example): children `G 0 0 10 2 * * * * * GTCAGGTGGT` and `G 0 0 10 !2 * * * * * TCGTGGATTA` (the left root then carries `split 0`, a divergence); natural spelling of the first is `TGGTGGACTG` = `LEFT2[-10:]`, which is what today's JSON `sequence` holds.

## 2.7 Shared golden vectors — the freeze gate *(v2)*

Before either implementation is treated as authoritative, both must pass the same files byte-exactly:

- `api/python/tests/data/traverse/codec_vectors.tsv` (a gtest and a Python test read it): floats given as **raw IEEE bits** (hex) with the expected text — 0, -0, 1, 0.5, 0.1, 0.0001, 1e-05, 2.5e-07, 1/3, inf, *(v3)* 1e16, 1e22, 1e23, 1.2345678901234568e20, 2⁵³, 2⁵³+2, every 10ᵏ for k = −20…25 and its two neighbouring doubles, `DBL_MAX`, the smallest normal and the smallest subnormal — plus a **differential section** of 10,000 non-negative finite doubles from a seeded generator over raw bits, expected strings produced once by the Python reference and committed; RANGES and SETEXPR against each normative base (shorter-wins, tie → explicit); percent escapes; front-coded names including multi-byte UTF-8 (`éfoo`→`ébar`), names with spaces, and the names `*`, `.`, `c:0`.
- *(v3)* **Whole-document fixtures** `api/python/tests/data/traverse/documents/*.{request.json,full.json,mgt}`: the C++ writer must produce each `.mgt` byte-exactly from its request (gtest + CLI), the library must parse it, `dump()` it back byte-exactly and reproduce `full.json` with `to_json()`. They include the review counterexamples of §9 and the **clipped-merge evidence case** of §5.1 with its expected claims.

MGT v1 is frozen when both sides pass all of them.

# 3. The JSON summary (per seed, `detail: graphlet`)

Envelope unchanged (`release`, `capabilities`(+2 keys), `strategy`(+`output.detail/timing`, `clamped`), `walk_rule`, `algorithm_version` stays `traverse-0.2`, `results[]`, `timing?`). Per seed:

```
{"seed": {"seed_id","validated_seed_id","seed_id_mismatch","length_bp","num_kmers","labels_from_seed",
          "labels_supporting_total","labels_dropped","labels_dropped_digest","num_labels","num_seed_labels"},
 "label_mode": "constrain|annotate", "duplicate"?: true,
 "arms": {"left"?: {"status","complete_to_bp","completeness_scope","frontier_remaining":{…},"labels_per_node":{…},
                    "cap_trigger"?:{…},"counters":{…all counters, by name},"branch_events_total","evidence":{…},"limitations":[…],
                    "counts": {"segments","leaves","splits","merges","runs","bases","max_bp",
                               "label_ends":{reason:n},"leaves_by_reason":{path_reason|"semantic":n}}},
          "right"?: {…}},
 "annotation": {"access_path","keys_mapped","rows_requested","direct_reads"}, "timing"?: {…},
 "outcome": {"walks": "complete|partial|failed", "branch_diagnostics": "complete|cut",
             "label_evidence": "complete|lower_bound", "delivery": "inline|spooled|paged"},   # (v5.1, §14)
 "resource_stop"?: {…},  "limitations": [ {…} ],                     # (v4, §14; per arm too: "evidence", "limitations")
 "graphlet": "<MGT text>", "graphlet_bytes": N, "graphlet_lines": N}
```
Derivation failures keep today's `{seed, error}` shape, no `graphlet`. No names, no per-leaf rows, no `label_summary`, no `growth` (in `B`), no `label_dict` (in `L`). Size ≈ 2.5 KB envelope + ~0.7 KB per arm; a 64-seed batch stays ~100 KB of summary.

## 3.1 Index identity *(v3, corrected in v4)*

`release`, `k`, regime and alphabet do not identify an index: two different annotations over one graph report the
same four values, and their `c:0` would be joined. v3 fixed that with a hash of *metadata* (counts, ordered column
names) — which the v3 review showed is still not identity: two indexes with identical counts and names that swap
which sequences columns A and B annotate get the same fingerprint, and the per-ref name check agrees too.

*(v4)* **`index_fp` is the digest of the immutable index bundle**, not of its metadata:

- The index build writes a **manifest** next to the files (`<prefix>.manifest.json`: every file of the bundle —
  graph, annotation, sidecars (`.seqs`, `.anchors`, `.rd_succ`, `.coords`; never the graph's derived `.edgemask` and
  `.bloom`, decision #17, §29) — with size and sha256, plus
  the builder's version and inputs), and `index_fp` = sha256 over the canonical manifest. Hashing hundreds of GB at
  server start is not an option; hashing at build time is free. For public indexes whose files were not built with
  a manifest, the deployment supplies one (computed once, or from immutable object-store digests such as S3
  checksums) and the server loads it with `--index-manifest`.
- `index_ns` (`--index-name`) names the index for humans and routing; it is not identity.
- Without a manifest the server reports `index_fp: *`, and everything that joins labels across retrievals reports
  *unverifiable* — never equal. The v3 metadata hash survives as `index_meta_fp`: a cheap **negative** check (a
  mismatch proves different indexes), never a positive one.
- Test (freeze gate): two indexes over the same records with identical counts and column names but swapped
  memberships must get different `index_fp` and must compare `unverifiable`/`different`, never equal;
  `scripts/traversal/build_mini_refseq.sh` writes the manifest so the fixture exercises the positive path.

# 4. The walker changes (marked separately) — *(v2: three run fields, was one)*

`src/graph/traversal/walker.hpp` `LabelRun`: `uint32_t segment = UINT32_MAX;` — the segment on which the run ended or was closed. Set at three sites, no behavioural effect:
1. `end_run()` (walker.cpp:1469) gains a `size_t segment` parameter; `end_label()` (:1479) passes `item.segment`; the silent switch-source end (:2451) passes `item.segment`.
2. `merge_level()`, in the `labels_via_parent` loop (:2000-2008): `else arm.result.runs[e.run].segment = item.segment;` for an entry whose run is not the kept one (that run was closed by this merge).
3. *(v2)* `LabelRun` also gains `uint32_t branches` and `double loss`, set at the same sites from the `Entry` that ends or is closed (`end_run()` already receives it). Reason (confirmed by the review with two annotations whose full JSON is identical): branch counts increment *before* quorum filtering, so a label ending inside a walk can carry `branches = 1` or `0` depending on a successor that left no trace in the output; `T` records these values only for leaf labels. Without the fields, `Claim.branches`/`loss` for interior ends would have to be reported as unavailable.

*(v2)* `Split::ambiguous` is already in `SeedResult`; it only has to be serialized (`G.split`, §2.2).

Justification: the run↔segment association exists only in `Item`/`Entry` state and is dropped from `SeedResult`. From the output alone it is ambiguous whenever two segments end the same label at the same depth (clones at splits are identical rows), it is absent for the silent ends of §0.4 (no event at all), and for ≥ 3-parent merges (proposal 3). Without it, `R` would need the events stored beside it (~+100 KB, the biggest block twice), or a sweep that reconstructs run ids — not reproducible under non-breadth-first frontier orders. Also exposed as `runs[].segment` in the JSON (additive) and asserted in `tests/graph/traversal/test_trie.cpp`: every ended run's anchor holds its `label_end` event at `to_bp`; *(v2)* for a silent `Lw` end, at least one run has `prev_run` = this run, starts at `to_bp` and was entered by switch (`from_label` = this run's label) — the switch events themselves may sit on the *children* of the anchor when the source ends at a split, so v1's "a SWITCH on the anchor's chain" was the wrong place to look; and every `m` run's anchor is a parent of the segment created at `to_bp`.

# 5. The local library: `api/python/metagraph/traverse/` (stdlib only; pandas lazy)

`setup.py`: `find_packages(include=['metagraph', 'metagraph.*'])` (today's `include=['metagraph']` excludes subpackages from a non-editable install). `metagraph/__init__.py` unchanged (does not import `client`, so `import metagraph.traverse` is pandas-free).

| module | contents |
|---|---|
| `_codec.py` | `encode_ranges/decode_ranges`, `encode_setexpr(ids, base)/decode_setexpr(tok, base)` (shorter-wins, tie explicit), `fmt_num/parse_num`, `pct_escape/unescape`, `REASON/QUAL` code tables, `LabelSetInterner` (content-hashed `array('I')`), `GraphletFormatError(line_no, msg)` |
| `model.py` | slotted dataclasses: `Label(id, kind, column, seq_id, name)`, `SeedInfo(…)`, `Segment(id, parents, from_bp, length_bp, walk, entry, entry_total, end, partition, split, first_base, presence: list[PresenceRun], events: list[Event], leaf: Leaf|None; derived: children, depth)`, `Event`, `PresenceRun(from_bp, to_bp, labels, total)`, `Leaf(path_reason, extras: dict[label,(loss,branches,route_bp)], continuation: Continuation|None)`, `Run(id, segment, label, from_bp, to_bp, end_code, reason, qualifier, silent, merged, route_bp, from_label, cost, prev_run, structural_successors, branches, loss, needed_budget)`, `Arm(side, status, complete_to_bp, scope, frontier, labels_per_node, cap_trigger, counters, growth, branch_events, segments, runs)`, `Graphlet(format, k, regime, alphabet, index_ns, index_fp, mode, support, reconverge, cap, continuation_bp, seed_index, seed: SeedInfo, labels, dropped, arms, summary: dict|None, derived_from)` |
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
parse(text) -> Graphlet                       # BODY ONLY (v2): no seed_id, derivation metadata, annotation counters or
                                              # replay strategy; to_json(), next_request(), deepen() and summary() raise
                                              # MissingEnvelope on such a graphlet instead of inventing defaults
                                              # raises GraphletFormatError(line_no, msg); validates A counts and Z
Graphlet.dump() -> str                        # canonical text; equals the input for an unmodified graphlet
Graphlet.from_response(result: dict, response: dict) -> Graphlet   # results[i] + envelope → summary attached
Graphlet.load(path) / .save(path)             # the .mgt file: H, J (envelope with this seed only), body
Graphlet.to_json() -> dict                    # today's detail: full results[i] (natural orientation); the conformance oracle
Graphlet.summary() -> dict                    # agent-facing ≤ 2 KB: per-arm status/complete_to/counts, top labels, caveats
Graphlet.label(selector) -> Label                # (v3) TAGGED selectors: {'ref': 'h:<column>:<seq_id>' | 'c:<column>'} |
                                              # {'name': str} | {'id': int}. A bare string is a convenience only: it
                                              # resolves when exactly one label matches it as a ref OR as a name, and
                                              # raises AmbiguousLabel otherwise (a column may be NAMED 'c:0'; annotate mode
                                              # can record two headers 'ACC1' from different columns). Every API result
                                              # carries {name, ref}, never a bare name alone
Graphlet.spell(arm, leaf, orientation='natural'|'walk', with_seed=False) -> str
Graphlet.walks(arm, *, top=None, by='support'|'length'|'loss', labels=None, route_consistent=True, min_bp=0) -> list[Walk]
    # Walk(path_id, leaf, segments, length_bp, sequence, path_reason, end_reasons, claims, labels_full, n_alive,
    #      complete, beyond_certified_bp); 'support' ranks by (-len(labels_full), -n_alive, -length_bp, path_id)
Graphlet.claims(arm=None, labels=None, at_most_bp=None, strict=True) -> list[Claim]
    # the §6.9 unit: one per RUN END (claims, not leaves). Claim(label, arm, path_id, segment, from_bp, to_bp, route_bp,
    # evidence_from (§5.1), kind∈{end, alive, diverged, merged, stretch, route_only}, reason, qualifier, end_class∈
    # {lost, dead_end, radius, capped, blocked, pruned, open}, entered_by, from_label, cost, loss, branches, support, exact)
    # strict=True raises IncompleteRecording when annotate nodes_truncated > 0 (a cut list makes the oracle a lower bound)
Graphlet.label_walks(label, arm=None) -> list[LabelWalk]        # per run: (from, to, evidence_from, sequence, end, leaves_below, merged_into)
Graphlet.routes(label, arm) -> list[list[int]]                   # per-label routes through merges via G partitions
Graphlet.support_profile(arm, leaf, kind='displayed'|'route') -> list[SupportRun(from_bp, to_bp, labels, exact)]
Graphlet.support_changes(arm, leaf) -> list[Change(at_bp, added, removed, reasons)]   # where support changes and why
Graphlet.label_summary() -> dict                                 # direct_bp / reach_bp / reentries / runs, both modes
Graphlet.continuation(arm, leaf) -> Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed())
Graphlet.next_request(arm, leaves, bp=None, reduce_budget=True, **overrides) -> dict   # resubmittable /traverse request
Graphlet.subgraph(selectors, arm=None, mode='any'|'all') -> GraphletView   # (v3) a VIEW over the backing graphlet:
                                              # every id is the backing graphlet's original id (no remapping, no new
                                              # record type); save() writes the UNCHANGED backing body plus, in J,
                                              # {view: {selectors, mode, arm, of: <entry handle or body digest>}}; load()
                                              # restores the view; its completeness is the backing complete_to_bp
                                              # qualified 'for the selected labels' (stated in the view's summary)
Graphlet.to_fasta(arm=None, leaves=None, with_seed=True, orientation='natural', width=None) -> str
Graphlet.to_gfa(with_seed=True) -> str          # S per segment on the seed strand, L with k-1 overlap, P per walk, LB/ER tags
Graphlet.compare(other, *, arm=None, labels=None, mode='claims'|'walks'|'labels'|'prefix_subset') -> Comparison
    # keyed by LabelRef (kind, column, seq_id) — names only for display (v2); restricted to min(complete_to_bp); claims
    # cut at D get end_class 'open'; comparable only for the same index identity (§3.1: equal index_fp, and equal
    # index_ns and release when set; names agreeing per ref) and
    # the same ORIENTED SEED SEQUENCE — not validated_seed_id, which includes the permitted label names and so
    # differs between a constrain and an annotate retrieval of the same seed (v2); results carry the support kind,
    # both strategies and both completeness scopes, and a per_path vs united_history pair is 'qualified', never
    # 'equal' — cutting both to min(complete_to_bp) does not make their termination comparable;
    # Comparison(comparable, reason, depth_used, equal, only_in_a, only_in_b, notes) — notes carry the merge qualification
Graphlet.memory_bytes() -> int
TraverseClient(host, port, api_path=None, *, session=None, timeout=900, release=None)
    .capabilities() -> dict;  .resolve(sequence, *, labels=None, discover=None, select=None, **opts) -> dict
    .traverse(seeds, strategy=None, *, detail='graphlet', timing=True, graph=None) -> TraverseResponse(envelope, graphlets, errors)
    .traverse_raw(request) -> dict;  .deepen(graphlet, arm, leaves, bp=None, **overrides) -> TraverseResponse
    # requests only; explicit Accept-Encoding: gzip, deflate; non-2xx → TraverseError(status, message); 503 → ServerInitializing
GraphletStore(spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800, ttl_disk_s=7*86400)
    # (v2) ENTRY identity and BODY identity are separate: an entry handle is opaque ('g_' + 12 random hex) and
    # owns (normalized request, index identity, envelope, body digest); bodies are deduplicated by sha256 in the
    # spool. Two requests with different loss budgets can produce byte-identical bodies (no switch happened) but
    # continue differently, so a body hash must never be the handle.
    .put(response, request) -> handle (opaque, v2); .get(handle) -> Entry (lazy parse, LRU); .free/.list/.save/.load/.sweep
```

Evidence semantics in the round trip: `route_bp` lives on the run (`R`, the *earliest* stamp) and per leaf label (`T`, the entry's *latest* stamp); `claims()` reports route support `[from_bp, to_bp)` and displayed-path support `[evidence_from, to_bp)` separately, with `evidence_from` derived **per claim and displayed path, never from `R.route_bp`** (v2, §5.1); "contiguous occurrence" is reported only under `support == 't'` (from `H`); per-label end reason = enum + qualifier; left-arm positions index `G` bases directly; completeness (`status`, `complete_to_bp`, `scope`, `cap_trigger`) gates every comparison; label-list cuts carry totals everywhere and `A.nodes_truncated`. Iterative deepening stays a backend call: the library derives the continuation and builds the request; it never answers deepening locally.

## 5.1 Normative derivations *(v2)*

The library must not need the C++ to reproduce these; they are part of the format contract.

**Displayed-path evidence** (`Claim.evidence_from`, `Walk.labels_full`, `walks(by='support')` ranking,
`support_profile(kind='displayed')`). v1 used `max(from_bp, R.route_bp)`, which is wrong after two non-first-parent
merges at `d1 < d2`: the run keeps `d1`, the entry and the leaf record `d2`, and `[d1, d2)` would be attributed to the
label on the displayed path. *(v3)* Provenance is derived **at the run's anchored endpoint** and then clipped — never by evaluating the walk at a cut depth, which skips a later merge and attributes the displayed prefix to a label whose own route spells something else (reproduced: B joins through a non-first parent at depth 14, the displayed path begins `AACC`, B's own route `TCTA`; v2 evaluated at `D = 4` gave `evidence_from = 0`). For the run of a claim, with `p` = the first-parent chain from the root to the run's **anchor** segment (also for merge-closed `m` runs) and `t = to_bp`:

```
evidence_from(lineage, p, t) = 0
cur = the lineage's label at t                       # follow prev_run backwards across switches
for each segment s on p with |parents(s)| > 1, in decreasing from_bp, from_bp(s) <= t:
    label_at = the lineage's INCOMING label at from_bp(s), i.e. before any switch committed at that depth
                                                     # (switches change the name: walk R.prev_run / from_label)
    if label_at not in partition(s)[0]:              # the lineage's kept entry came through a non-first parent
        evidence_from = from_bp(s); break            # the displayed bases before s are not this lineage's
evidence_from = max(evidence_from, from_bp of the claimed run)   # a label never claims bases before its own run
```
For leaf labels the merge-derived part (before the `max` with the run start) equals `T.route_bp` (asserted by T37); interior claims use the same walk. *(v3)* A claim cut at `D < to_bp` reports the displayed-support interval
`[evidence_from, to_bp) ∩ [0, D)`; when that is empty the claim is `route_only` at `D`: the label carries a route
of that length, but not the displayed bases; `routes()` reconstructs and spells the label's own ancestral route
through the partitions when the caller asks for it.
*(v4)* **Run-start guard, applied first:** intersect the run's own interval `[from_bp, to_bp)` with the cut
`[0, D)`. If it is empty — the run starts at or after `D`, e.g. B entered by a switch at depth 10 and the cut is at
5 — the run makes **no claim** at `D` (it does not exist yet); only a run with surviving route support can be
`route_only`. A run of zero length (`from_bp == to_bp`: a label present at the boundary node only, e.g. ended
`minority` at the seed boundary) is a **boundary claim**: reported at its position with `kind: boundary`, no
bases, and only when `from_bp ≤ D`.

**`label_summary`** (constrain), per arm, exactly `Walker::summarize()`:

```
for r in runs (in run order):                  # prev_run(r) < r always, so one forward pass suffices
    root[r] = root[prev_run(r)] if r entered by switch and prev_run(r) set else label(r)
for r in runs:
    summary[label(r)].runs.append(r)
    if not entered_by_switch(r) and from_bp(r) == 0:
        summary[label(r)].direct_bp = max(direct_bp, to_bp(r))     # a merge does not clamp it
    summary[root[r]].reach_bp = max(reach_bp, to_bp(r))           # credited to the lineage ROOT only
for every switch EVENT e:  summary[to(e)].reentries += 1          # events, not switched runs: split clones
                                                                  # copy from_label/cost without a 2nd switch
```

**`label_summary`** (annotate), per arm, exactly `Walker::summarize_annotate()` *(v3, normative)*:

```
for s in segments, for each P run of s, for l in labels(P):      # recorded lists (cut lists: lower bounds)
    summary[l].reach_bp = max(reach_bp, to_bp(P))
alive_end = {}
for s in segments in id order:                                   # parents precede children
    alive = entry(s) if s is the root else ∪ alive_end[p] for p in parents(s)
    for each P run of s in order:
        if alive is empty: break
        still = alive ∩ labels(P)
        for l in alive − still: summary[l].direct_bp = max(direct_bp, from_bp(P))
        alive = still
    for l in alive: summary[l].direct_bp = max(direct_bp, from_bp(s) + length_bp(s))
    alive_end[s] = alive
# runs = []  and  reentries = 0  for every label (no lineages are tracked in annotate mode)
```


# 6. MCP tool surface (functions in `mcp_tools.py`; every return ≤ `max_bytes`; ids never leave the server, agents see `{name, ref}`; every local answer carries `evidence: {complete_to_bp, exact, support, reconverge, scope}`)

*(v2) Paging and sizes.* List tools (`graphlet_walks`, `_claims`, `_splits`, `_labels`, `_support`) take `cursor` and return `next_cursor` (opaque, bound to the entry handle, the tool and its normalized arguments; a cursor presented with other arguments is rejected), plus `total`. `max_bytes` defaults to 2048 for list tools; `graphlet_sequence` is the one tool with its own ceiling (16 KB, stated in its result) and the way to get a long sequence. A single row larger than `max_bytes` is returned alone with `row_truncated: true` and the fields that were cut named (e.g. a walk's spelled tail), never silently shortened. `traverse_fetch(replay=<handle>)` re-runs the entry's stored normalized request against the same index (same release required; a different release is an error, not a silent re-run). *Oversize:* `max_graphlet_mb` is a RAM threshold — a larger body is spooled to disk complete and only the summary and handle are returned; a hard transport/storage limit rejects the fetch explicitly. A truncated MGT document is never stored or returned as a graphlet.

Backend-calling: `traverse_capabilities(index)`; `traverse_resolve(index, sequence, labels?|discover?, select?, limit=10)`; `traverse_fetch(index, seed:{sequence, labels?, seed_id?}, strategy, keep=True, max_graphlet_mb=8) -> summary + handle` (the probe is this with tight bounds); `traverse_continue(handle, arm, walk, overrides?, execute=True) -> new summary + parent:{handle, arm, walk, overlap_bp}` (or the request when `execute=False`).

Local over a handle: `graphlet_summary(handle, arm?)`; `graphlet_walks(handle, arm, rank=support|length|loss, n=10, min_bp=0, label?, spell=none|tail|full, tail_bp=60)`; `graphlet_walk(handle, arm, walk, cursor?)` (segment chain, labels in/out per segment, continuation; *(v3)* the chain is paged by `cursor` like a list); `graphlet_support(handle, arm, walk, step?)` (support runs with who left/joined and the reason: a `K`-coded end or "split: N labels took C"); `graphlet_labels(handle, arm?, rank=direct_bp|reach_bp, n=20, min_direct_bp=0, at_bp?, name?)` (table; `name` → that label's runs/routes incl. merge routes); `graphlet_splits(handle, arm, n=10, min_labels_before=2)`; `graphlet_claims(handle, arm, min_bp=0, route_consistent=True, labels?)`; `graphlet_sequence(handle, arm, walk|segment, from, to, orientation, with_seed)` (≤ 16 KB slice); `graphlet_export(handle, format=fasta|gfa|json|mgt, arm?, walks?, path) -> {path, bytes, records}`; `graphlet_compare(a, b, arm, mode, cursor?) -> Comparison dict` (preconditions *(v2/v3)*: same index identity (§3.1) and the same oriented seed sequence, else `comparable: false`; differing completeness scopes → `comparable: qualified`; *(v3)* returns the counts and the first page of `only_in_a` / `only_in_b` / `differ`, with `next_cursor`; `graphlet_export(format=json, what=compare|walk, …)` writes the complete lists to a file — a truncated field is always recoverable by cursor or export); `graphlet_subtrie(handle, labels, arm?) -> derived handle`; `graphlet_list/free/save/load`. Unknown handle → "replay with traverse_fetch(replay=<handle>)" (the store keeps the request).

**§6 as implemented in stage 1 (the tool layer, `mcp_tools.py`)** — where the code refined the text above:
the evidence block on every local answer is `{complete_to_bp, exact, support, reconverge, scope, outcome,
limitations}` (plus `view`/`qualified` on a derived handle; `exact` requires complete label evidence and no cut
lists); `graphlet_list` is paged (`rows`, `total`, `next_cursor`); `graphlet_walks` takes `route_consistent`, and
walks/claims report `filtered` counts (e.g. `filtered.route_only`) instead of dropping rows silently;
`graphlet_summary(detail='limitations')`; rows carry `labels_in_total`/`labels_out_total`/`n_labels_total`; errors
are structured results with codes `unknown_label`, `backend_error`, `backend_unreachable`, `io_error`,
`path_not_allowed`, `not_replayable`, `view_unsupported`, `not_in_view`, `result_too_large`; exports and saves take
file names confined to the export directory (default `<spool>/exports`); the cursor secret is created in the spool
when none is passed. There is no streaming validator: a parse builds the model once and releases it, so the peak is
the full model (§14's "streaming parse" is not delivered).

# 7. Size and memory (SRA case; estimates from per-record byte costs, to be re-measured by `scripts/traversal/graphlet_measure.py`)

**Measured on the SRA primary index (stage 1, `graphlet_measure.py --server`, 2026-10-02)** — these replace the
estimates below as the reference numbers:

| retrieval | body | gzip | full JSON (compact) | body / JSON | gzip / gzip | summary | parse | model |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 16S (1,396 bp), exhaustive constrain, 300 bp | 460 KB | 53 KB | 3.99 MB | 11.5 % | 27 % | 3.6 KB | 46 ms | 4.4 MB |
| 16S, label-free trie (annotate), 300 bp | 434 KB | 93 KB | 4.59 MB | 9.4 % | 28 % | 2.9 KB | 46 ms | 3.4 MB |
| 16S, default walk, 3 kb | 3 KB | 1 KB | 10 KB | 33 % (fixed overhead: seed, header) | 73 % | 2.4 KB | < 1 ms | < 0.1 MB |
| random 50-mer, label-free beam 20 (most supported), 10 kb | 3.5 MB | 196 KB | 75 MB | 4.7 % | 3 % | 3.6 KB | 0.9 s | 30 MB |

Targets of §7/§10 met on the case they were set for (< 15 % of compact JSON; summary < 10 KB). The interner hit
rate was 0.997 (constrain) and 0.83 (annotate, 2,264 distinct sets) — no per-arm alias table is needed. The small
mini_refseq retrievals (18–29 %) are dominated by fixed per-document overhead.


*(v2) These are estimates under an interning assumption, not bounds.* The model's memory depends on how many label sets are distinct: 5,000 distinct 827-label arrays are ~17 MB for the arrays alone before model objects, interner keys, names and caches, and a synthetic fragmented-set case produced ~17 MB of `SETEXPR` text — so v1's "1 MB pessimistic" is not a bound. `graphlet_measure.py` must report set cardinalities, membership changes per step and the interner hit rate on the SRA retrievals; a per-arm alias table is considered only after those numbers.


Constrained left arm (2839 segments / 1437 leaves / 1402 splits / 4160 ends / 69 kb bases): `G` 2839 × ~30 B + 69 KB bases ≈ 155 KB; `R` ~4200 × ~24 B ≈ 100 KB; `T` 1437 × ~8 B ≈ 12 KB; `E`/`V`/`B`/`A`/`S`/`L` ≈ 8 KB → **≈ 275 KB ≈ 7.5 % of 3.72 MB** (target < 15 %), ~60-90 KB gzipped vs 180 KB. Rebuilt locally: paths 0.71 MB, splits 0.41 MB, labels_at_end, label_summary, 1437 continuation strings. Label-free (annotate): `G` ≈ 155 KB; `P` ~5000 runs delta-coded ≈ 150 KB; `L` 2836 file-path names front-coded ≈ 100-150 KB (0.42 MB raw: the intrinsic part); `T` 12 KB → **≈ 430-480 KB ≈ 10 % of 4.55 MB**, ~120-160 KB gzipped vs 330 KB. Pessimistic (names share no prefix, fragmented sets) ≈ 1 MB ≈ 22 %; `labels.max_labels_per_node` bounds it by construction, every cut keeps its true count. Default tuned walk: ~2.5 KB body vs 10 KB. Summary < 10 KB per seed (§3). Server: the body is one `std::string` appended straight from `SeedResult` (no `Json::Value` tree for it); peak transient = body + JSON copy + gzip ≈ 3 × body < 1.5 MB. Client: parse ~15k lines in 50-100 ms; model ~5 MB per constrained arm, ~20 MB label-free (interned `array('I')` sets).

# 8. Files to change

C++: `src/cli/traverse.hpp` (`detail` doc; `strategy_to_json(st, cost, detail, timing)`; `std::string graphlet_text(const SeedResult&, const Seed&, const Strategy&, Arm order…)`; `Json::Value seed_result_summary_json(…)`; `constexpr int kGraphletFormatVersion = 1`), `src/cli/traverse.cpp` (enumeration :329-330; echo :530-535; `GraphletWriter` in the anonymous namespace ~350 lines: codec helpers `ranges/encode_setexpr/fmt_num/pct_escape/reason_code`, `header/seed/dropped/labels/arm/growth/branch_events/segment/presence/events/leaf/continuation/runs/trailer`; rule computations for the `*` fields; `arm_to_json`/`seed_result_to_json` branches for `graphlet`; additive JSON fields `runs[].segment`, `segments[].labels_via_parent`, hairpin `followed`, revisit `length_bp`/`same_distance`, `label_dict[].column/seq_id`; `label_ends` → `Json::objectValue`; `capabilities_to_json` +2 keys; `process_traverse_request` passes `seed` and computes `duplicate` before serialising), `src/graph/traversal/walker.hpp/.cpp` (§4, ≤ 6 lines). Not changed: `server.cpp`, `server_utils.cpp` (an optional, separate 400/500 split per spec §10.3 is recommended but independent).

Python: `api/python/setup.py` (find_packages), new `api/python/metagraph/traverse/{__init__,_codec,model,parser,derive,ops,export,client,store,mcp_tools,frames}.py`, `api/python/tests/test_traverse_{codec,parse,derive,ops,client}.py`, `api/python/tests/data/traverse/*.{request.json,full.json,graphlet.json,mgt}`, `scripts/traversal/graphlet_fixtures.py` (regenerate + `--check`), `scripts/traversal/graphlet_measure.py`.

Tests/docs: `integration_tests/test_traverse.py` (new classes, §9), `tests/graph/traversal/test_trie.cpp` (anchor assert), `tests/graph/traversal/test_trie_cases.cpp` (a bubble+switch+hairpin case whose dump feeds the Python fixtures), `docs/SPEC-labeled-traversal-core.md` (§5 detail graphlet + complete echo, delete label_lists/label_runs; §7.1 additive fields, derivable fields noted, same-position events unordered, `label_ends: {}`; new §7.5 "The graphlet (MGT v1)" = §2 verbatim incl. the `*` rules and orientation; §10.3 capabilities keys, transport note; §11 T37-T40; §13 `LabelRun::segment`, silent switch-source ends, orientation decision), `docs/source/api.rst` ("Traversal: retrieve once, query locally"), REVIEW-REQUEST §5.1.7 resolved.

# 9. Test plan

C++ / integration (`integration_tests/test_traverse.py`, CLI and HTTP; `skipUnless metagraph.traverse importable`, the venv installs `api/python` editable):
- T37 round trip: for constrain tuned (merge), constrain exhaustive (keep), annotate exhaustive, annotate beam, quorum, left-only, both arms, caps-tripped `max_steps`, `max_labels_per_node: 1`, 3-seed batch with a derivation failure: run `detail: full` and `detail: graphlet` on the same request; `Graphlet.from_response(r, out).to_json()` equals the full `results[i]` after normalisation (events sorted within equal `at_bp`, `needed_budgets` sorted, `timing` removed); `parse(text).dump() == text`; `A` counts match; `graphlet_lines == Z`; `check_rules()` reports no divergence; `#label_end events ≤ #runs`; every `Lw` run is the `prev_run` of at least one run entered by switch at its `to_bp` (§4; the switch events may sit on the anchor's children); `is_canonical(text)` holds for every writer output.
- *(v2)* T37 must also cover the review's counterexamples: an ambiguous split with disjoint child sets after a switch, a divergence with overlapping child sets, a followed-hairpin split, a quorum-filtered split, an A→B→A re-entry inside one segment, a silent switch-source end at a split (switches on the children), two interior ends whose branch counts differ only through a quorum-rejected successor, two runs through two non-first-parent merges (`R.route_bp` ≠ `T.route_bp`), and two same-name headers from different columns in annotate mode. The CI job installs the Python package so T37 cannot silently skip.
- T38 orientation: `spell(left) + seed + spell(right)` reproduces `acc1/acc2/acc3`; left `continuation(leaf).sequence == paths[].continuation.sequence` *(v2: for continuations contained in the flank and for ones crossing into the seed, on both arms)*; resubmitting a locally derived continuation returns 200 and extends.
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

# 12. Changes in v2 (after the external design review of v1)

Confirmed by the reviewer with CLI reproductions or source citations; each change is made in place and marked *(v2)*.

| # | v1 problem | v2 change | where |
|---|---|---|---|
| 1 | split `ambiguous` derived from overlapping child label sets — wrong both ways (lineages before quorum; followed hairpins excluded) | stored: `G.split` | §2.2, §2.5, §4 |
| 2 | interior claims' `branches` (and `loss`) not recoverable; identical full JSON for different values | `LabelRun::branches/loss`, `R` fields; terminal values only, unknown at an earlier cut | §2.2, §4 |
| 3 | codec not interoperable: `to_chars` vs `repr`, unspecified SETEXPR bases for `G.end`/`C.labels`, prefix units | positional shortest round-trip floats; normative bases; UTF-8 byte prefixes at character boundaries; shared golden vectors as freeze gate | §2.1, §2.7 |
| 4 | `evidence_from = max(from_bp, R.route_bp)` attributes `[d1, d2)` after two non-first-parent merges | per-claim derivation through merge partitions and `prev_run`; equals `T.route_bp` at leaves | §5, §5.1 |
| 5 | continuation spelling "last n bases" wrong for the left arm | right: last n of seed+right; left: first n of left+seed | §2.2, T38 |
| 6 | names cannot identify all labels (two `ACC1` headers in annotate mode) | `{name, ref}` everywhere; ambiguous bare names rejected | §5, §6 |
| 7 | `more: N` without a cursor is not paging; 2 KB vs 16 KB contradiction; undefined `replay` | cursors bound to handle+query; per-tool ceilings; oversized rows flagged; `replay` defined | §6 |
| 8 | body hash as entry handle aliases different requests | opaque entry handles; body dedup separate | §5 |
| 9 | `G.end` set formula fails on A→B→A re-entry | chronological rule; explicit fallback legitimate, T37 checks emitted `*` only | §2.3 |
| 10 | silent switch-source assertion looked for the switch on the wrong segment | validated through targets' `prev_run`/start position | §4 |
| 11 | `label_summary` only specified by reference to C++ | normative pseudocode; `reentries` = switch events | §5.1 |
| 12 | body-only parse vs envelope-backed operations conflated | `MissingEnvelope` for operations that need the envelope | §5 |
| 13 | comparability via `validated_seed_id` (differs between modes) | oriented seed sequence + index identity; scope mismatch is `qualified` | §5, §6 |
| 14 | subgraph id semantics; inherited completeness | a view with original ids; completeness qualified by the predicate | §5 |
| 15 | memory "pessimistic 1 MB" not a bound | stated as an estimate under interning; measure first | §7 |
| 16 | oversized retrievals undefined | spool complete + summary, or reject; never truncate | §6 |
| 17 | T37 could skip; missing counterexamples; `end_labels` order; annotate partitions | CI installs the package; counterexamples listed; label order; one empty list per parent | §2.2, §2.5, §9 |

Confirmed by the reviewer as correct in v1: constrain merge partitions union to the merged entry set; the anchor
assignment sites cover every run closure; leaf-end membership is recoverable; annotate root totals preserve the
cut-boundary information; frontier orders have deterministic path-id tie breakers; spaces and literal `*`/`.` in names
are representable when the parser keeps the final-field remainder exactly. The two-track plan stands, with the
grammar and evidence corrections above made before step 0's freeze.

# 13. Changes in v3 (after the external review of v2)

| # | v2 problem | v3 change | where |
|---|---|---|---|
| 1 | `GraphletView.dump()` promised remapped ordinals + an original-id column that no record defines | views keep original ids; `save()` writes the unchanged backing body + the view selector in `J`; no new record | §5 |
| 2 | index identity `(release, k, regime, alphabet)` joins unrelated labels (two annotations over one graph, empty release) | `index_ns` + content fingerprint `index_fp` in `capabilities` and `H`; missing → unverifiable; names checked per ref | §2.2, §3.1, §5, §6 |
| 3 | `to_chars(fixed)` and `Decimal(repr)` differ on `1e23`, `1.2345678901234568e20` | one algorithm: shortest scientific digits, expanded positionally; raw-bit vectors incl. powers of ten, `DBL_MAX`, subnormals; a 10,000-value differential section | §2.1, §2.7 |
| 4 | evidence evaluated at a cut depth skips a later merge (false displayed support) | provenance at the run's anchored endpoint, then intersected with the cut; `route_only` when empty; chain ends at the anchor (incl. `m` runs); incoming label before switches at `d` | §5.1 |
| 5 | ref strings collide with legal names (a column named `c:0`) | tagged selectors; bare strings resolve only when unique as ref or name | §5 |
| 6 | canonical `dump()` vs preserved input encodings | canonical form defined (writer must emit it); identity guaranteed for canonical documents; `is_canonical`; whole-document fixtures | §2.3, §2.7, §9 |
| 7 | stale text: worked records, model table, claims comment, T37 silent-switch assertion; annotate `*` partition arity | all updated; `*` decodes to one empty list per parent | §2.2, §2.3, §2.6, §5, §9 |
| 8 | annotate `label_summary` not normative | pseudocode matching `Walker::summarize_annotate()` | §5.1 |
| 9 | `graphlet_walk` chains and `graphlet_compare` lists not recoverable when large | cursors on both; complete lists by export | §6 |

Confirmed by the review as correct in v2: stored split ambiguity, the chronological `G.end`, terminal run metadata,
both continuation formulas, UTF-8 prefix boundaries, silent-end validation through `prev_run`, opaque entry handles,
`MissingEnvelope`; the full-depth evidence scan is sound with the anchor and switch-boundary qualifications above.

# 14. Resource and guarantee contract *(v4)*

**The principle (the owner's requirement).** Every response certifies what it covers and states every limitation:
either an answer is complete, or a valid shallower answer is returned together with exactly where it stops, what was
cut, why, and which knob or action would go further. Nothing is cut silently, and no resource measure may change the
question being answered without saying so.

**Outcome: four independent dimensions *(v5)*.** v4's single value could not say "partial traversal, spooled",
and "all evidence present" was stronger than any one boolean (an annotate run can reach its radius with every branch
event kept while capped label lists make its summaries lower bounds). Per seed result (JSON `outcome`, MGT `O`):

| dimension | values | guarantee |
|---|---|---|
| `walks` | `complete` · `partial` · `failed` | `complete` only when no walk-class limitation applies: every requested arm is complete to the radius **per path** (keep, or merge where no history was united) and no carrier of the seed was dropped; `partial`: any of `walk_domain`, `scope` with at least one merge, `seed_labels` applies — a valid certified prefix whose limits the `limitations` state; `failed`: no valid traversal exists (structured reason in `limitations`/`resource_stop`) |
| `branch_diagnostics` | `complete` · `cut` | every branch decision and refusal is reported; `cut`: only before each arm's `evidence.complete_to_bp` (`branch_events`) |
| `label_evidence` | `complete` · `lower_bound` · `qualified` | `complete` only when no label-class limitation applies; `lower_bound`: evidence may be *missing* or understated — cut lists (`label_lists`, `inexact_counts`), dropped carriers (`seed_labels`), a cut switch-source list (`switch_sources`), greedy losses (`greedy_losses`); `qualified`: something reported may be *overstated* — a column-label trace across a record boundary (`trace_record_boundaries`), or a permitted set derived from part of the seed (a `derivation` on a walked result: D3, feature level 6, §26; *review of 2026-10-06, X1*); `qualified` wins when both apply |
| `delivery` | `inline` · `spooled` · `paged` | the body is in the response; `spooled`: complete in the spool behind the handle, the response holds the summary; `paged`: delivered in pages. Independent of `walks`, so a partial traversal can be spooled |

*(v5.2, the owner's rule)* **Conservative by construction**: each dimension is `complete` only when no limitation of
its class applies, so a reader can never find a limitation in `limitations` whose dimension still reads `complete`.


**Stated limitations** (already being implemented on the JSON side, spec §7.0): per arm `evidence: {complete,
complete_to_bp}` and `limitations: [{kind, knob, limit, observed, effect, complete_to_bp?}]` with kinds `walk_domain`,
`branch_events` (the stored branch events stop at `max_branch_events`; every branch decision and refusal before the
first dropped event's depth is present — a run cut there says so and names `output.max_branch_events`, which also
accepts `"unlimited"`), `label_lists`, `switch_sources`, `inexact_counts`, `scope`; per seed `seed_labels`,
`server_clamp`. MGT carries them: `A` gains `<evidence_complete_to_bp|*>`, and a new `K` record per limitation
(`K <kind> <knob> <limit> <observed> <complete_to_bp|*>`), so a saved graphlet keeps its own caveats.

**Budgets enforce the bound; depth is the graceful fallback.** A depth limit alone does not protect against a very
wide first level or one enormous annotation row (annotate mode materialises a full row, or every coordinate tuple,
before cutting the recorded list — label_oracle.cpp). So:

- *What counts* *(v5)*. **Memory**: every allocation the request retains at once — the walker's state (heads'
  label states, segments, runs, events, recorded sets), the request's `LabelQuery`/`LabelRecorder` caches, the
  decode scratch of the annotation reads it issues (row-diff dependency rows, coordinate tuples), and the output
  buffers. The index itself (mmapped graph and annotation) and process-wide shared caches are not charged to a
  request; a shared cache that grows because of a request is charged to it. **Work**: *charged work units*, a
  weighted sum of counters the walker already keeps — successor enumerations, rows decoded (weighted by their
  row-diff dependency path length), coordinates mapped, pair evaluations, refusal scans, edge-reuse probes,
  re-minimisation rounds — not CPU time and not elapsed time; `max_steps` alone undercounts (one step can perform
  many source–target comparisons). Elapsed time stays a separate deadline. The deadline and the budgets are checked
  at least every *W* work units (a bounded interval, stated in `capabilities`).
- *Row decoding is budget-aware inside the decoder* *(v5)*. Estimating the *final* row does not bound peak memory:
  row-diff reconstruction loads its dependency rows, whose differences can cancel into a tiny result
  (`row_diff.hpp`). The annotation reads used by traversal go through a decode path that charges each dependency row
  and tuple as it is materialised and aborts the batch, without side effects, when the request's remaining budget
  would be exceeded; the step then stops with `resource_stop.phase = annotation_decode` and the suggestion to use a
  more selective seed or a label-constrained query (narrowing the radius would not help). Until that decode path
  exists (§14.1) the memory bound is *soft* and every response says so (`limitations: memory_bound_soft`, observed
  = the measured overshoot).
- *Atomic commit per head* *(v5)*. "Check before each allocation" is not enough when earlier mutations already
  happened: the walker ends a switch source's run before it allocates the split's children, so a denial in between
  would leave inconsistent runs, events or split records, and lowering `complete_to_bp` does not repair them.
  Processing a head becomes two-phase: **plan** (successors, derived states, label ends, splits, merges — no
  mutation of the result) → **admit** (reserve the plan's allocations *and* the cost of later stopping and
  delivering what it creates, see the next bullet) → **commit** (apply; cannot fail). A head that is not admitted is
  censored exactly like a head beyond a cap today (its walk ends at the last committed level with `resource_limit`),
  and the level it belongs to does not count toward `complete_to_bp`.
- *Delivery is accounted per expansion, not reserved as a fraction* *(v5)*. A fixed reserve fails on a comb-shaped
  trie (one terminating branch per split): segments grow linearly, but materialising every leaf's ancestor chain is
  quadratic. So (a) paths are **lazy** — `PathResult` keeps its leaf; ancestor chains are produced while serialising,
  from parent pointers, through bounded buffers (MGT stores no chains at all; JSON `detail: full` streams them) — and
  (b) admission charges each committed head with the incremental delivery cost of what it adds in the requested
  output format. The guarantee: **every committed prefix remains deliverable within the resources already reserved
  for it.**

**Scopes.** *Locus*: peak working memory and cumulative work across both arms, retries and continuations of one
locus. *Analysis*: total work across loci, an overall deadline, retained graphlet memory and spool storage.
*Service*: concurrent workers and aggregate memory reservations. Requests carry a stable `budget_id` (analysis), a
`locus_id`, and *(v5)* a fresh `attempt_id` per dispatch; splitting a request or retrying never resets an allowance —
memory limits what is held at once, work accumulates across attempts.

*(v5)* **An authoritative ledger with reservations**, not reported usage: two requests that each receive the same
remaining allowance before either reports would both spend it. The MCP layer / search service keeps a ledger per
`budget_id` and `locus_id` (atomic operations, e.g. Redis) and **reserves** an allowance for an attempt before
dispatching it; the reservation is the request's bounds; the backend returns its usage on every response, successful
or not; the ledger **reconciles** reserved against used on the response, and **releases** a reservation on
cancellation, or on lease expiry when a response is lost (the attempt is then charged its full reservation, since
its spending is unknown). The backend enforces the locus scope per request; the ledger enforces the rest.

*(v5.1)* **Charging work and releasing capacity are separate operations.** Consumed work is charged when the response (or the lease expiry) settles the attempt; occupied capacity — the worker slot and the memory reservation — is released only once the backend has actually stopped: on its response, or after the lease, which is only a valid release point because the backend enforces the same deadline itself (the request's `time_budget_ms` plus the server's hard request timeout, so an attempt cannot outlive its lease). A cancellation the backend has not acknowledged releases nothing. *(Review of 2026-10-06, X2: as built, the bound is compared only at the delivery checks, so an attempt outlives it by its run up to the next one — the rest of the piece the bound fell into and, when its walk had not stopped, the walk up to its next poll that reads the clock, the stopped seed's finalisation and the building up to that check —, whose length has no stated bound (§27; "one uninterruptible step" understated it, the review of the fixes found); the lease is a valid release point under that stated assumption, and an answer of `running` or `stopping` past it shows that run still going.)*

**No silent semantic changes.** Enabling a beam, dropping labels, raising a cut, switching `on_reconverge` or the
support kind changes the question or the evidence; none of them is ever applied as a resource measure. When the
agent chooses one, the response says what it changed (`strategy` echo + `limitations`/`scope`).

**The local library has budgets too.** Parsing (a compact MGT expands into large label sets), route enumeration,
comparisons and exports run under allocation and work budgets; an interrupted comparison returns
`comparable: unknown` / `incomplete`, never equality, and an interrupted parse is a **local failure**: the original body (or the handle to it) is kept unchanged and the failure is reported — a shallower graphlet cannot be manufactured from a cut-off parse.

**Continuation is a new traversal.** It does not resume the original exhaustive search (the backend keeps no
frontier or edge history between requests); its result is certified on its own, with the overlap stated.

**Freeze-gate fixtures for this section.** Exhaustion during annotation decoding (including a row-diff row whose
dependencies are dense but whose result is tiny), mid-level traversal, finalisation and serialisation: a valid
shallower result is delivered whenever one exists (`walks: partial` with the right `complete_to_bp` and
`resource_stop`), `failed` only when none does. *(v5)* Allocation denial injected inside a head's commit around a
switch, a split and a merge (the result stays consistent: runs, events, splits and paths agree); a comb-shaped trie
(finalisation stays within the admitted reservation); concurrent attempts against one ledger entry (no double
spending; a lost response charged its reservation); a partial traversal delivered spooled (`walks: partial`,
`delivery: spooled`); a non-zero `refusal_scans` and an unknown future counter round-tripping through `A`; an
interrupted Python parse keeping the original body and reporting a local failure.

## 14.1 What is frozen with MGT v1, and what is enforced in stages *(v5)*

Frozen with the format (the wire contract): the `O`, `Q`, `K` records and the `A` counters/evidence fields; the
`resource_limit` end reason; the four outcome dimensions; the `limitations` kinds (*(v5.4)* their `effect` sentences are free text, which may be reworded within MGT v1; the kinds, knobs and value types are what is frozen; *(v5.8, decision X-C12)* frozen means that no kind a reader knows changes its meaning, knob or value types — the kinds themselves are an **open `[a-z_]` token set**, as are a `K` record's extra field names and a `Q` record's scope, resource, phase and action tokens: a later feature level may add one, and a reader keeps a value it does not know, uninterpreted; the precedent is `memory_bound_soft` (28c51c31), level 6 adds the kind `coordinates` with the extra field `lists_cut` and the action `drop_coordinates` (C12); no record, field or grammar rule changes, which would be MGT v2); the
`resource_stop` structure; `attempt_id`/`budget_id`/`locus_id` in requests and usage in every response. Enforced in
stages, each stated in every response until it lands (`limitations` entries): (1) **now** — the existing caps, the
evidence boundary, stated limitations, the outcome dimensions (being implemented in the JSON); (2) the two-phase
head commit and lazy paths in the walker; (3) the budget-aware row-diff decode path in the annotation library
(`memory_bound_soft` until then); (4) the ledger in the search service; and, separately, (L) work and allocation budgets for the local library's
operations (compare, routes, export; *(v5.3)* not enforced yet, and said so in the library docs — `max_bytes` bounds
the bytes a tool returns, not its computation or peak allocation; when they land, an interrupted comparison reports
`unknown` and an interrupted route or export says it is incomplete). Stage 2 is implemented; how its bounds are
stated is in §19.1. The promise at every stage is the same:
deliver the largest certified prefix that can be safely committed and delivered, and state the limiting resource and
the useful next action — the stages only move hard bounds from *stated* to *enforced*.

# 15. Changes in v4 (after the external review of v3 and the owner's guarantee requirement)

| # | v3 problem | v4 change | where |
|---|---|---|---|
| 1 | `index_fp` hashed metadata: swapped memberships with identical counts/names collide | digest of the immutable index bundle via a build manifest; no manifest → unverifiable; metadata hash kept as a negative check; swapped-membership test | §2.2, §3.1 |
| 2 | `route_only` created claims for runs that start after the cut | run-start guard first; future runs make no claim; zero-length boundary claims defined | §5.1 |
| 3 | no representation of resource stops; no guarantee statement | `outcome` (complete/partial/failed/deferred), `resource_stop`, `resource_limit` end reason (`Y`), per-arm `evidence` + `limitations` (also in MGT: `A` field, `K` records), budgets per locus/analysis/service, reserve for delivery, row-level guard, local-library budgets, exhaustion fixtures in the freeze gate | §2.1, §3, §14 |
| 4 | `max_branch_events` cut refusal evidence silently | boundary `evidence.complete_to_bp` + a `branch_events` limitation naming the knob; `"unlimited"` accepted (implemented now, spec §7.0/§7.2) | §14 |

# 16. Changes in v5 (after the external review of v4)

| # | v4 problem | v5 change | where |
|---|---|---|---|
| 1 | a row guard on the final row does not bound decoding (row-diff dependencies); `max_steps` undercounts work | decode path that charges dependency rows and tuples and aborts without side effects; charged work units; bounded interval between checks; what counts as memory defined; `memory_bound_soft` stated until it lands | §14 |
| 2 | "check before each allocation" can leave inconsistent records (a switch source ended before split children are allocated) | two-phase plan → admit → commit per head; denied heads censored like a cap; allocation-denial fixtures around switches, splits, merges | §14 |
| 3 | a fixed delivery reserve fails on comb-shaped tries (quadratic path materialisation) | lazy paths serialised from parent pointers; delivery cost charged per admitted expansion; every committed prefix deliverable | §14 |
| 4 | wire contract incomplete: no grammar for evidence / `K`; `refusal_scans` cannot round-trip; `outcome`/`resource_stop` unassigned | `A` counters named and extensible + evidence field; `O`, `Q`, `K` records with scope and typed values; effect reconstructed from normative templates | §2.2, §3 |
| 5 | one outcome value conflates walk completeness, branch diagnostics, label evidence and delivery | four independent dimensions (also being implemented in the JSON now) | §14 |
| 6 | aggregate budgets by reported usage allow double spending | ledger with atomic reservations, reconciliation, lease expiry; `attempt_id`; usage on every response; work = charged units | §14 |
| 7 | the parse rule implied a shallower graphlet could be made from a cut-off parse | local failure, original body/handle kept | §14 |
| 8 | the freeze gate missed the above | fixtures added; §14.1 separates the frozen contract from staged enforcement | §14, §14.1 |

# 17. Changes in v5.1 (after the fifth review: approved, three bounded items)

| # | item | change | where |
|---|---|---|---|
| 1 | `K` values could not hold strings (the scope limitation's limit is `merge`); effect templates were promised but not supplied | typed values `i:`/`f:`/`s:`/`u`; `effect` stored as the free-text last field | §2.2 |
| 2 | the JSON summary still showed the scalar outcome | the four dimensions | §3 |
| 3 | cancellation / lease expiry could release capacity while the backend still runs | charging work and releasing capacity separated; capacity released only on the backend's response or after a lease the backend's own deadline guarantees | §14 |

Implementation order (from the review): the C++ writer, the Python parser and the shared conformance fixtures first;
then the staged enforcement of §14.1 (atomic commits and lazy paths, the bounded decode path, the ledger), with every
not-yet-enforced guarantee marked soft in the responses.

# 18. Record coordinates — an additive, v1-compatible extension *(owner's decision, 2026-10-02)*

MGT v1 is frozen, so coordinates do not go into the body. They are an **additive JSON field** of the response, carried
with the envelope the library already keeps (the `J` line of a saved `.mgt`), ignored by v1 readers that do not know
it. Rationale and the v2 alternative are in the discussion of 2026-10-02: the size cost is small (about 2–3× the bytes
of a body record for this part before compression, ~1.2× after gzip, and only for trace retrievals), and the Python
side parses JSON faster than MGT text. Reasons that would justify an MGT v2 are collected in §18.4.

## 18.1 What is reported, and when

Coordinates are reported **only where they are well defined**: under `support: trace` on an index with k-mer
coordinates, where every claim is one coordinate-consecutive occurrence in one record. Under `support: kmer` a stretch
can be stitched from several occurrences, so no single position is honest; nothing is reported, and the response says
so (`coordinates: null` with `coordinates_reason: "support kmer"` / `"index has no coordinates"`).

For **header labels** (an accession via `CoordToHeader`), coordinates are positions **within the record**, 0-based,
in the record's own orientation; for **column labels**, positions in the column's global coordinate space (and the
existing `trace_record_boundaries` limitation applies). Trace runs only on basic (single-strand) graphs, so the
strand is always the record's forward strand on the right arm; on the left arm the walk runs backwards along the
record, which the interval below already accounts for.

## 18.2 The field

Per seed result (JSON, every detail level — `full`, `tree`, `summary` and `graphlet`):

```
"coordinates": {
  "kind": "record" | "column",            // header labels: within the record; column labels: global
  "k": 31,
  "seed": [ {"label": <label id>, "occurrences": [[start, end], ...]} , ... ],
  "arms": {
    "left"|"right": [
      {"run": <run id (the R ordinal / runs[] index)>, "label": <label id>,
       "occurrences": [[start, end], ...]}    // one per live coordinate chain of the run
    ]
  }
}
```

- An occurrence `[start, end)` is a **half-open base interval in record (or column) coordinates**. Normative rule:
  the run's first node (the k-mer entered by step `from_bp + 1`) has k-mer coordinate `c₀`, its last node `c₁`. The
  run's **own bases** — the base each of its steps adds, i.e. the walk's bases `[from_bp, to_bp)` — are record bases
  `[c₀ + k − 1, c₁ + k)` on the right arm, where coordinates increase outward, and `[c₁, c₀ + 1)` on the left arm,
  where they decrease outward (a left-arm step adds the first base of the k-mer it enters). A seed occurrence is
  `[c_first, c_last + k)` over the seed's first and last k-mer.
- Several occurrences per run occur when the label's record repeats the walked sequence (two live coordinate chains);
  they are listed in ascending `start`.
- Runs closed by a merge do not occur (trace and merging are mutually exclusive). Runs entered by a switch carry the
  new label's own coordinates from their start — *(v5.8, amended by decision X-18.2)* except a run entered by a
  switch into a label whose own lineage is still live there: the walker's entry keeps only that label's chains
  that continue from before the switch, so chains starting at the switch node are missing; such a run carries the
  inherited chains only and is marked `lower_bound: true` (decision C-N2 a: the block's `complete` false,
  `runs_lower_bound`), the walk unchanged. It arises only with a switch cost (the default is `forbid`); M1 found
  21 of 3,780 runs, all switch-entered, in 10 of 56 trace cells (§26).
- Bounded like everything else: at most `output.max_coordinate_occurrences` (default 16, accepts `"unlimited"`) per
  run; a cut is stated as a limitation `coordinates` (knob `output.max_coordinate_occurrences`, observed = the true
  count) — never silent. *(v5.8, as built, §26)*: observed is the largest true count among the cut lists, the
  extra field `lists_cut` their number, and each cut list states `occurrences_total`; the block also has `kind`
  `mixed`, `max_occurrences`, `complete`, `runs_lower_bound` and per run `from_bp`, `to_bp`, `chains_ended`
  (C-N1) and `lower_bound`, and the intervals are written from the last node's coordinate c (`[c + k − L, c + k)`
  right, `[c, c + L)` left, L = `to_bp − from_bp`), equivalent to the c₀ form above.

## 18.3 The library

- `Graphlet.coordinates` (from the envelope; `None` with the reason when absent). Body-only graphlets have none and say
  so (`MissingEnvelope` for coordinate queries).
- `Claim.coordinates`: the occurrences of the claimed run(s), **clipped to the claim's displayed interval** (and to a
  cut depth `D`) with the same arithmetic as §18.2; `route_only` claims carry the coordinates of the label's own route.
- `walks(...)`: per walk and label, the record intervals covering the walk's bases; `label_walks()` rows likewise.
- `to_fasta()`: headers gain `acc:start-end` (natural orientation, 1-based closed in the header for readability, stated
  in the docs) when coordinates are present.
- MCP tools: claims/walks/labels rows carry `coordinates` (counted against the byte ceilings; cut with `more`).
- **Oracle**: on mini_refseq the real-index suite fetches each record's sequence from the fixture FASTA and checks that
  the bases at every reported interval equal the spelled claim (natural orientation, both arms) — positional
  verification on top of the search-based one.

## 18.4 Reasons collected for an eventual MGT v2

1. Coordinates in the body (self-contained raw bodies, one conformance regime) — today in the envelope (§18).
2. An extension mechanism (readers ignore unknown records of a reserved form), so later additions need no version
   bump — v1's strict reader rejects unknown records.
Cut a v2 only when more than one strong reason has accumulated; keep v1 reading.

# 19. Changes in v5.3 (as implemented, after implementation reviews 3 and 4) *(2026-10-03)*

None of these changes the MGT v1 wire format; the spec carries the normative text (§5, §6.8, §7.0, §7.1).

## 19.1 The resource contract as implemented (amends §14)

- **Per seed.** Budgets apply per seed: a request of n seeds can hold n × budget. A seed's stop does not depend on
  its position in the request.
- **Work is stated with its number, not a fixed W.** Work is charged when a fetch returns a row ("rows decoded": 8
  units per key, 1 per entry and per coordinate), plus label-state scans, pair evaluations, edge-reuse probes and
  enumerations. The budget is compared after every charge, so a stop exceeds it by at most what was charged since
  the previous comparison. Every work stop states the most its seed charged between two comparisons
  (`largest_charge`). `W` = 65 536 bounds the seed phase's interval and sizes fetch calls; it is no longer the
  overrun guarantee. The roots' rows of both arms are one charge, so that a stop still delivers a result complete
  to 0 bp. Rows the lookahead decodes without a fetch asking for them are not charged, and the spec says so.
- **Memory.** The account is a deterministic model, never measured RSS. Delivery is priced at worst-case escaped
  lengths, once per place a name is written. What a fetch, a validation or the depth-0 dictionary holds beyond the
  account is recorded before any check can stop the walk, and is stated as `memory_bound_soft.observed`.
  `memory_bound_soft` is stated on every result under a memory budget, failures included. A failed result echoes
  index-supplied names cut to 256 bytes under a memory budget; the request's own `seed_id` is stated, not cut.
- **Stops.** A `resource_stop` is created only when a head below the radius was actually censored. A memory bound
  is not a wall-clock bound, and delivery is not checked against the deadline; a hard overall deadline would need
  enforcement during delivery (staged).
- **The header index** (header → column, sequence) belongs to the `CoordToHeader` object and lives as long as it
  does; there is no address-keyed cache.

## 19.2 Library and tool decisions (amend §5, §5.1, §6; the reviewer's recommendations, adopted)

1. **Annotate claims.** A claim is a maximal end of a label's routes under the union rule of §5.1. Its displayed
   support starts at the last merge its witness route enters through a non-first parent. The witness is the first
   carrying parent in stored order. One witness does not enumerate all routes, and does not establish a
   contiguous occurrence in a source sequence.
2. **Route stamps.** `R.route_bp` is the earliest non-first-parent merge on the lineage's route (a run closed by a
   merge of three or more parents may carry `route_bp == to_bp`). The `T` route stamp is the latest entry. Displayed
   support is `[evidence_from, to_bp)`, derived from the partitions — not `[R.route_bp, to_bp)`.
3. **Names.** `next_request()` refuses names shared within the retrieval and names containing U+FFFD. A later
   capability could relax the U+FFFD rule, but only with its provenance saved in the retrieval; a server's
   capability today cannot vouch for an older body.
4. **Identity.** With a manifest on one side only, continue and replay refuse with `index_unverifiable` ("cannot
   be verified, not proven different"). A fallback needs the explicit opt-in `allow_unverified_index` (default
   off).
5. **Comparison at a depth.** The DAG is restricted to the comparison depth, and merges at or beyond it are left
   out (and named). Where that cannot be exact, the result is `qualified`; a side without bases gives `unknown`.
6. **Continuations.**
   - The loss budget is reduced by the largest terminal loss among the continued labels. This conserves the
     cumulative loss under unchanged costs: no continued route exceeds the original loss budget, and the result
     says so for the lower-loss labels. It is not equivalence with the uninterrupted walk: branch counts, edge
     history and other per-walk state restart in the new traversal (the notes state it), so a continuation can
     reach further than the uninterrupted walk would have (review of `d8d3ef86`, P7).
   - `labels.extra` is rebuilt around the new seed labels.
   - Continuations that cannot share one request raise an error, with `next_requests()` as the fallback.
   - Labels alive at a leaf but not seeded are named (`left_out`).
7. **Tool layer.**
   - Claims and walks filters both count `merge_entered`.
   - Every error honours an explicit `max_bytes`; the floor is 64 bytes.
   - The 16 KiB dry-run ceiling is stated and can be overridden.
   - A successful export, save or load always keeps its receipt (handle or file locator): optional fields are
     dropped first.
   - Local operations have no work budget yet (§14.1, stage L).

## 19.3 Decided after the review of `d8d3ef86` *(v5.4)*

- A continuation reduces `max_label_branches` by the largest terminal branch count of the seeded labels
  (preserving `"unlimited"`), states it, and offers an explicit reset.
- Extra labels are accepted when every label of the declared pool is reachable from the seed labels within the
  loss budget through the cost model (transitively, not in one switch); the walk enforces cumulative costs.
- The `seed_id` echo of a failed response is named in `memory_bound_soft`.

# 20. Stage 3 as implemented *(v5.4, 2026-10-03)*

Stage 3 of §14.1 is implemented for the row-diff annotation formats (`RowDiff<BRWT>`, `RowDiff<ColumnMajor>` and
the coordinate-aware `TupleRowDiff<...>`; the spec lists them).

- **Charging inside the decoder.** The traversal's annotation reads go through an opt-in, single-threaded decode
  path (`IRowDiff::decode_rows` / `decode_row_tuples` with a `DecodeBudget`). It charges the trace, each dependency
  row, each coordinate tuple and the reconstruction buffers *before* they are allocated, and refuses a read whole,
  without side effects, when it would not fit. The default decode path is unchanged.
- **Batch independence.** Every key a fetch returns is admitted against its standalone demand, so where a walk stops
  does not depend on `batch_kmers` or the lookahead. A stop states only facts that hold either way.
- **Stops.** A refused read in a level stops the seed in phase `annotation_decode` (actions `more_selective_seed`,
  `label_constrained_query`); in the seed phase, derivation or an annotate root it fails the seed. Labels that do not
  fit, and level lists that leave nothing for the read, are `traversal`-phase stops.
- **Work is deterministic logical work** (the review's answer). Each returned row is charged its row-diff dependency
  path (8 units per dependency row plus 1 per stored entry), whatever was shared in the actual decode; it is not
  measured decode effort, which depends on batching and caching. The physical decode counters in `timing` are kept
  separately. `capabilities.work_bound` says so.
- **Lower bounds.** A depth-0 failure from unread roots or unnamed labels states its need as "at least X"; raising the
  budget to X may still fail, because unread roots, labels or later state need more, and the statement says so.
- **What remains soft.** `memory_bound_soft` stays on every result under a memory budget and names what is still
  held uncharged: a level's lists and fetched rows until its heads are processed, the seed phase's intersection and
  hits, a failed seed's depth-0 dictionary, an index-wide header lookup, the dictionary's first table, and a failed
  response's `seed_id` echo. Stage 3 does not establish a hard request-wide memory bound.
- **Stage 3b** (not built): the same budget-aware path for formats without row-diff (BRWT, ColumnMajor and TupleCSC
  alone, the disk and flat row formats, `IntRowDiff`) and for `/resolve`; charging a level's rows and lists in the
  account.

# 21. Decisions for the next stages *(owner, 2026-10-03)*

The owner accepted the planners' recommendations (decision sheet C/R/B/L). Implementation order after stage 4:
per-graph identity → record coordinates → stage 3c → stage 3b → stage L.

**Record coordinates (§18).**
- C1: opt-in. `output.coordinates: true` is echoed only when given. With it, §18 applies exactly (an object, or
  `null` plus a reason); without it every output is unchanged.
- C2 *(provisional; the owner is not yet sure)*: the `coordinates` limitation (a cut occurrence list) belongs to no
  outcome class; the block states `complete: false`. To be revisited after the first review.
- C3: a dictionary mixing header and column labels gives kind `mixed`, and each label's kind says how its positions
  are read.
- C4: per-run `from_bp`/`to_bp`, plus `max_occurrences`, `complete` and `occurrences_total` on cut lists.
- C5: an occurrence is a chain that carries the whole run; a per-run `chains_ended` count appears when it is > 0.
- C6: support is advertised in `GET /traverse/capabilities` only. The owner notes that a general, server-wide
  capabilities document could be useful too (see the identity pass).
- C7: `graphlet_claims` rows always carry coordinates; `graphlet_walks` / `graphlet_labels` only with
  `coordinates=true`.
- C8: the library requests coordinates for trace strategies when the server supports them, and `next_request`
  carries the setting forward.
- C9: an inconsistent block is rejected eagerly.
- C10: column labels report global column positions in v1; mapping them to records is deferred.
- C11: the default cap is 16 occurrences per run.
- C12: new free tokens: K kind `coordinates` (extra `lists_cut`) and resource-stop action `drop_coordinates`.

**Stage 3c (selected-label decoding).**
- R2: `annotation.access_path` stays `tuples`; the selected mode is reported in timing and capabilities.
- R3: a selected read charges 8 per dependency row, 1 per BRWT probe and 1 per selected coordinate.
- R4: fall back to full-row decoding above about 64 distinct columns, or by a cost estimate.
- R5: the binary rd_direct of spec §8.2 is built in the same change.
- R6: no local refseq33m copy for now.
- R7: no pre-estimate warning for now; revisit after 3c.

**Stage 3b.**
- B1: every format except the four legacy ones (rbfish, rb_brwt, bin_rel_wt, row/EigenSpMat), which stay stated as
  uncovered.
- B2: ColumnMajor is read as column-wise charged pairs.
- B3: the slower BRWT single-row descent is accepted under a memory budget, to be measured on production indexes.
- B4: work-only budgets keep the default decode on formats without row-diff.
- B5: each head's rows are released once it is processed.
- B6: a head stop caused by a level's rows states the bytes, with no `more_selective_seed` action.
- B7: a stopped `/resolve` returns exactly the resolve of a query prefix.
- B8: new `/resolve` JSON fields and tokens as proposed.
- B9: discovery under a budget keeps every met label's runs.
- B10: an explicit `select.policy` that stops fails the request.
- B11: the default `--resolve-max-time-ms` is 0 (unlimited).
- B12: no server-side default memory budget for `/resolve`.
- B13: the memory bound stays soft for now.
- B14: the disk reader is charged as a measured constant.
- B15: the uncharged lookahead key list and the misleading `/resolve` error text are fixed inside 3b.
- B16: `annotation_reads` and `resolve_budgets` appear on `GET /traverse/capabilities` only.
- Also: a cap on `/resolve`'s explicit label list (`--resolve-max-labels`).

**Stage L.**
- L1: deterministic cold-price work units.
- L2: off by default.
- L3: a fetch whose parse stops stores the body as `parsed: false`.
- L4: no partial exports.
- L5: the agent may raise a budget up to the service's ceiling.
- L6: re-parse billing is service policy.
- L7: defaults are view 30 M / 1 GiB, heavy 300 M / 2 GiB, parse 200 M / 2 GiB, ceilings 10×.
- L8: the quadratic paths are separate items, except the one-line `graphlet_subtrie` fix.
- L9: memory is a modelled account (`memory_bound: "model"`).

**Capabilities.** `feature_level` is monotonic: each pass that adds capabilities fields or routes bumps it by one, and
SPEC §10.3 records the mapping (2 = attempts and the client-gone stop; 3 = pass 5, §23). `schema_version` stays the request schema
version (`strategy.schema_version`), which a bump would break.

# 22. Stage 4, backend half, as implemented *(v5.5, 2026-10-04)*

The ledger itself lives in the search service. The backend provides what a ledger needs; the SPEC (§5, §7.0, §7.3,
§10.3) is normative.

- **Attempts.** A request may carry `attempt_id` (`^[A-Za-z0-9._:-]{1,128}$`, unique per server process while
  running or retained), plus `budget_id` and `locus_id`, which are echoed and need an `attempt_id`. A reused id
  (running, retained or tombstoned) is refused with 409, whose body carries the other attempt's state; nothing
  runs twice. Requests without `attempt_id` keep their bytes.
- **Cancel and state.**
  - `POST /traverse/cancel {attempt_id, wait_ms}` answers:
    - 200 when it flagged a running attempt;
    - 404 when the attempt has finished;
    - 404 with `tombstone: true` when the id is unknown: the id is tombstoned, so a later request with it is
      refused (409);
    - 429 when the tombstone store is full (no promise).
  - `GET /traverse/attempt/{id}` returns `running`, `stopping` or `finished` with the reason (`completed`,
    `cancelled`, `client_gone`, `deadline`, `error`). It creates nothing, so its 404 is not completion. Only a
    tombstoned cancel answer, or a finished state, from the same `server_instance` says an attempt will not run
    there.
  - Finished attempts are kept for 3600 s / the last 10,000. Tombstones are kept by age, under their own cap.
- **Usage.** Every response to a request with an `attempt_id` carries a response-level `usage` block: the ids;
  `server_instance`; per seed and in total, work units, modelled memory (admitted, final, soft excess,
  `held_bound_bytes`, which is null without a memory budget because the account is then no bound on what is held),
  elapsed ms, and seeds started, finished and abandoned; and the bound actually used. This covers success,
  partial, failed and cancelled seeds, and 400/500 after registration. The 409, the 400s before registration and
  the loading 503 carry no usage; a client that is gone gets no response.
- **The attempt bound** = min(seeds × the effective per-seed `time_budget_ms` + allowance (10 s by default),
  ~899 s). The server enforces it itself:
  - seeds stop being walked at bound − allowance/2;
  - building, writing and compressing the response are checked against the bound;
  - an attempt at its bound answers 503 `{error, usage}` with reason `deadline`, delivering no seeds; its partial
    results come as a 200 instead (`resource_stop {scope: "attempt", resource: "attempt_deadline"}`, unstarted
    seeds failed with phase `not_started`).

  Enforcement is cooperative: it acts at the next checkpoint. Without `bounds.max_work_units`, a level's fetch
  and a seed's validation are single uninterruptible calls. Chunking them by time comes next, together with
  `not_after_ms` and the stated `max_uninterruptible_ms`.
- **Client gone.** Every `/traverse` (and `/resolve` between its phases) stops when its client has closed the
  connection, writing nothing. A half-closed connection counts as gone. The socket is polled at the walker's
  checkpoints at most every 100 ms, and the kernel's TCP state is read when bytes are pending.
- **Memory stops** state the smallest budget (whole MiB) whose account holds the need beside that budget's own
  cache allotments, so raising the budget to the stated value holds the state. A lower bound stays "at least".
- **New free-token values**, which MGT v1 admits:
  - resource `cancelled` and `attempt_deadline`;
  - scope `attempt`;
  - phase `not_started`;
  - action `retry_attempt`;
  - knob `attempt_id` on the stop's `walk_domain`. This is a new knob value, like `bounds.max_memory_mb` in stage 2;
    the kinds and value types are unchanged.
- **Capabilities.** `feature_level` 2 (§21). `GET /traverse/capabilities` adds an `attempts` block: fields, id
  pattern, routes, `server_instance`, retention, allowance, client check interval, and how the bound is computed.
- **Deployment.**
  - `-p` must exceed the number of concurrent traversals, so that the cancel and state routes find an io thread.
  - A reverse proxy must forward `/traverse/cancel` and `/traverse/attempt/*`, must close its upstream connection
    when the client goes, and must keep its read timeout at or above the largest bound.


# 23. Pass 5 as implemented *(v5.6, 2026-10-04)*

The identity pass of §21 with the owner's additions: a start deadline, the server-wide capabilities, deadlines
that reach into long annotation reads, and a delivery window that holds for large batches. MGT v1 is unchanged
(no record, field or token changed; the library's `J` envelope record reduces `usage`, below). The SPEC (§4.1, §5, §6.8, §7.3, §7.5.2, §10.3) is normative;
`feature_level` is 3.

- **`not_after_ms`.** An optional top-level request field (Unix epoch ms, an integer up to 2⁵³ − 1, with or without
  `attempt_id`): the instant after which the request must not be started. The server refuses, when the handler
  starts, a request whose instant has passed on its clock — 409 `{error, state: "expired", not_after_ms,
  server_time_ms, ids, server_instance}` — with nothing run and nothing registered (a later GET is a 404). The
  check is strict; the allowance a ledger adds is stated (`attempts.clock_skew_allowance_ms`, 2000 by default). A
  running, retained or tombstoned id is answered first, so `expired` always means that no attempt with the id
  exists there. Echoed in usage and in the state. The CLI applies the same check.
- **Per-graph identity.** The multi-graph list gains two optional columns, `manifest_path` and `index_ns`; each
  manifest is checked at start-up against the files its pair loads (sizes and the digest of the list, as
  `--index-manifest`), and a mismatch, two identities for one pair, more than five columns or a bad name refuses
  to start. Every response states its pair's identity. `GET /traverse/capabilities?graph=<name>&graph_path=<path>`
  describes one pair. A `graph_path` naming one graph listed with several annotations is refused (before, the
  first annotation answered and the others were never read: a silent narrowing). There is no union oracle; the
  SPEC no longer claims one, nor a parsed release column. `index_manifest.py --server-csv` writes one manifest
  per pair, every distinct file hashed once by parallel streams or taken from precomputed sha256 digests.
  *Review of pass 5*: the manifest check covers the sidecars the loader reads (row-diff anchors and fork
  successors, the dummy-edge mask, column coordinates, sequence headers), in both modes; a manifest that lists
  another graph or annotation is refused (one written for a directory lent one fingerprint to each of its
  annotations, swapped memberships included); two different pairs never state one `index_fp`, and lines are
  grouped by the pairs' real paths; a name listing one pair twice needs no `graph_path`, and one listing one
  graph with several annotations is refused as such. The tool matches digests by real path (relative paths
  from the digest file's directory), by base name only with `--match-base-names` and never one line for two
  files, and writes a pair's manifest to the path any of its lines names.
- **`GET /capabilities`.** A small, index-free document: routes and features, `feature_level`,
  `algorithm_version`, `release`, `server_instance`, `mode` and the graph names, the attempts block,
  `deadline_check`, `compression_level`. It answers while a single index loads (`ready: false`). The probe gains
  `algorithm_version`, `compression_level`, `deadline_check`, and in `attempts` `hard_cap_ms`,
  `clock_skew_allowance_ms`, the `not_after` rule and the `delivery_reserve`; `allowance_ms` is an integer.
- **Chunked deadlines.** Under any deadline, every `/traverse` annotation read (a level's fetch, the lookahead, a
  validation, a derivation's window) that the deadline may fall into is decoded in chunks, the deadline checked
  between them. *Review of pass 5*: a read is one piece, as before, when at the slowest per-row time the request
  has seen it would take less than 1/64 of the time left; otherwise its first chunk is at most 8 rows, each next
  at most 4 times the previous and sized at the rate the previous measured to take `chunk_target_ms` (50 by
  default), and the rest is one piece once predicted at that rate to take less than a quarter of the time left;
  chunks are taken in the walk's order. The first version sized every chunk from the slowest recent per-row
  time and cut sorted keys: on a row-diff annotation, whose rows share the decoding of their row-diff paths
  within a call, that decoded the paths again per chunk (0.75 ms a row against 0.016 ms whole), so a row_diff
  walk of 1.8 s took 30 s and came back partial, row_diff_brwt walks were 1.4-2.9 times slower, and a first
  chunk sized from another level's cheap rows took 650-880 ms of a 50 ms target. A read's counting, cache
  eviction and result are decided once for the whole read, so a request no deadline stops is byte-identical
  (checked on real requests, also with 1 ms chunks and with every read split into one-row chunks), and a
  stopped read censors the walk at the read. The rows a stopped read decoded are charged. A validation is
  stopped by the attempt only, never by its own time budget. One chunk (at least one row) stays
  uninterruptible, and so does a read decoded whole far from its deadline (a cancel is seen after it); no bound
  on one row exists before stage 3c, so `max_uninterruptible_ms` is null, and the longest single piece is
  stated as an observation (per process in `deadline_check`, per attempt in `usage`). Stated limits: a whole
  read overruns when its rows are more than 64 times slower than any seen before, the rest of a split read when
  more than 4 times slower than its chunk's. On the fan-out fixture budgets of 50-500 ms are kept to within 15
  ms with a column annotation and 60 ms with row_diff and row_diff_brwt (with and without a memory budget; the
  previous build overran by up to 0.9, 1.0 and 1.9 s), and walks that no deadline stops keep their bytes and
  time (60 real UHGG requests: 15,973 ms against 15,932 ms before and 16,021 ms unchunked).
- **The delivery window.** The traversal routes compress at zlib level 1 (630 MB/s against level 9's 153 MB/s on
  real responses, for 1.8 times the bytes). The server writes each seed's result as text once built (its tree
  freed at once) and assembles the response from the texts, byte for byte. Attempts stop walking at
  `bound − max(allowance/2, reserve)`: the reserve is the time to compress the text written so far and the
  walked seed's estimated text (its modelled account / an account per text byte) and to build the latter, times
  1.25, plus the time from the walk-until to the walk's end (`stop_ms`: the chunk target + 200 ms, or the longest
  the server measured recently) — *review of pass 5*: without the margin and the stop time a response the model
  fitted exactly was a 503 about half the time (a walk stopped 124 ms after its walk-until). The
  configured values (20 per text byte for JSON, 40 for a graphlet; compress 50, build 10 MB/s) apply until the
  server has measured its own on responses and seeds of 1 MiB or more: then the slowest rate and the smallest
  ratio (per detail) of its last 16 measurements, and the attempt's own where more conservative. A response is
  no longer lost whole for want of time; the reserve can stop a large attempt earlier than before (with the
  configured guesses alone, a 403 MB SRA response, 88 account bytes per text byte, stopped walking at 2.3 s on
  the M5 Max, which is why measurements replace them), and what it did not walk is stated. The configured
  values are uncalibrated on the staging server.
- **Library.** `standalone_text(body, result, response)` writes a seed's standalone `.mgt` text (H, J, body)
  without parsing the body, byte for byte what `save()` writes for a server body. The `J` envelope's usage is
  reduced to the totals plus this seed's `per_seed` entry (in `save()`, `dump()` and the store). `view` and
  `derived_from` are reserved envelope names. The client gains `not_after_ms`, `AttemptExpired`,
  `capabilities(graph=, graph_path=)` and `server_capabilities()`. A file saved by the stage-4 library with
  every seed's `per_seed` entry in `J` still loads as it is, but is not canonical under the reduced rule
  (`is_canonical()` false; `save()` rewrites it reduced): `J` is the library's record, MGT v1's records and
  tokens did not change.
- **Flags** (*review of pass 5*): `--traverse-clock-skew-ms`, `--traverse-chunk-target-ms` and
  `--traverse-attempt-allowance-ms` take integers in [0, 2⁵³ − 1]; a negative value had wrapped to 2⁶⁴ − 1 and
  was stated as such.
- **Deployment** (SPEC §10.3): `-p` above the number of concurrent traversals; a reverse proxy forwards the
  cancel, state and both capabilities routes with their query strings, closes its upstream connection when the
  client goes and keeps a read timeout of at least `hard_cap_ms` plus the transfer; the server's clock is
  synchronised within `clock_skew_allowance_ms`.

# 24. The review of pass 5, and the path cache's first reads *(v5.7, 2026-10-04)*

An external review of pass 5 (HEAD `5801aea1`) found the absolute deadline formula sound under its stated clock-skew
and uninterruptible-piece qualifications and MGT v1 intact, and four blockers for index identity and for
cancellation-based lease release (findings 1-4), two moderate findings (5, 6) and one in the library (7). A staging
benchmark of the efficiency pass (feature level 4) showed first reads slower than at feature level 3 (R10). Every
confirmed finding is fixed; MGT v1 is unchanged (no record, field or token). The SPEC (§4.2, §5, §6.8, §8.2, §8.4,
§10.3, T52) is normative; these changes are **feature level 5** (the constant is raised with the release).

- **Tombstones and early release (finding 1).** With `retention_s` 1, a request with `not_after_ms` = now + 60 s was
  half uploaded and cancelled (404, `tombstone: true`), a repeat cancel at 0.76 s gave the same promise, and the upload
  completed at 1.18 s ran: a tombstone held for `retention_s` alone does not justify releasing capacity. Now
  `POST /traverse/cancel` takes the cancelled request's `not_after_ms`, and a tombstone is held at least
  `retention_s` and until the server's clock reads `not_after_ms + clock_skew_allowance_ms`, capped by
  `--traverse-attempt-tombstone-max-s` (86,400 by default; `max(cap, retention_s)` as applied, stated). Every
  tombstone answer — the cancel's and a repeat's 404, the state's 404, and the 409 refusing a copy — states
  `suppressed_until_ms`, the `not_after_ms` it was judged against and `covers_admission` (true iff a `not_after_ms`
  was named and the expiry reaches it plus the skew; else `covers_admission_reason`). A repeat never shortens a
  tombstone and extends it to its own `not_after_ms`; a refused copy (the dispatch-time 409) extends it to the
  copy's own `not_after_ms`, so every later copy — a proxy's replay carries the same instant — is refused while it
  could still be admitted. The hold is kept while either clock says it is live (the steady clock for the computed
  duration, the wall clock until it reads the expiry), so a forward step of the wall clock cannot shorten it; a
  backward step by more than the skew after the expiry can revive a copy's instant (the ledger's skew assumption).
  Tombstones live in memory, so a restarted process holds none: a request may name `expect_server_instance`, and a
  process with another `server_instance` refuses it before anything runs (409 `instance_mismatch`). **The release
  rule** (normative, also in the capabilities): early release only on 404 + `tombstone: true` + `covers_admission:
  true` from the same `server_instance`, for an attempt sent with exactly that `not_after_ms` and with
  `expect_server_instance` naming that instance; an attempt sent without `not_after_ms` is never released early;
  otherwise only on a finished state or after `not_after_ms + clock_skew_allowance_ms + bound_ms`. Assumed: clocks
  within the skew, no rewriting of the request's ids, and no copy reaching another server that serves the same
  ledger. *(Completed at feature level 6, §27: the 409 refusing a copy as a finished-state source, the `expired`
  409 as a release ground, and the clock release's own assumption — the attempt past its bound apart from its
  run up to its next delivery check, whose length has no stated bound.)* The ledger's sentence is corrected: after `not_after_ms` + skew an unanswered attempt **cannot start
  subsequently** (it may already be running; `bound_ms` covers that), not "never started".
- **Retention settings (finding 4).** `--traverse-attempt-retention-s -1` started, stated 2⁶⁴ − 1 s, and expired every
  tombstone at once. The retention seconds, the count and the tombstone cap are bounded integers (a year; ten
  million), anything else refuses to start naming the option and its range. `retention_s` 0 means no suppression:
  a cancel of an unknown id is a 429 with `tombstone: false, reason: "no_suppression"` (a full table:
  `tombstones_full`), and the capabilities say so.
- **One loader dependency inventory (findings 2, 3).** Two pairs of symlinks to one graph and annotation with
  different `.seqs` beside them shared one manifest check (bundles were grouped by the main files' real paths), and
  stated one `index_fp` while answering with different accessions; a masked graph's `.bloom`, loaded by
  `DBGSuccinct::load` and able to turn `graph_runs` `[[0, 36]]` into `[]`, was in no manifest (the mask and the
  `.bloom` left the inventory again with the owner's decision #17, §29, this risk stated). `index_load_inventory`
  now lists every file the server's loaders open for a pair, derived from the listed spelling as the loaders derive
  it: graph; mask when it opens, then Bloom filter when it exists; annotation; anchors and fork successors beside the
  graph for `.row_diff`; `.seqs` for a coordinate annotation unless `--no-coord-mapping`. A unit test checks its
  annotation table against the types `initialize_annotation` builds. It decides the manifest check (a file of the
  inventory a manifest does not name refuses to start, naming it) and the bundle dedup (two lines are one index
  only when their whole inventories resolve to the same real paths in the same roles, so every distinct bundle is
  checked against its manifest). `metagraph traverse --index-inventory -i G -a A` prints it as JSON with its
  table; `scripts/traversal/index_manifest.py` (single mode, `--server-csv`, `--digests`) mirrors it, and an
  integration test compares the two on every sidecar kind. Not in it, because nothing opens them: a column
  annotation's `.coords`, the graph's `.weights`, anchors beside another annotation type, the in-memory header
  index. `index_fp` covers what it claims again — every file the pair loads, by size at start-up and by digest at
  build time (a same-size replacement is still found only by `--verify`); `index_meta_fp` is unchanged, a negative
  check.
- **Tombstone holds, refined.** A refused copy's suppression is judged against its own `not_after_ms`, a finished
  attempt's id is held like a tombstone, and the hold is live through `suppressed_until_ms` inclusive (the review of
  these fixes, below).
- **The budgeted lookahead (finding 5)** takes its runs and pieces in the walk's order, each sorted, as the default
  reads do (31,297 decoder charges against 161,679 on the reviewer's fixture). **Delivery checks (finding 6)** come
  every 64 KiB of a large token too (a piece is copied up to the next check); one token's preparation stays
  uninterruptible and its time is part of `observed_max_uninterruptible_ms`. **Finding 7** (the library's malformed
  continuation costs) is fixed in the library.
- **R10, the path cache's first reads.** On staging a first header-named walk ran 9,947 ms on a 1,000 ms budget with
  689 ms of reading; the time outside the read timers was the **sequence header index**, built by the first request
  naming a header (33 M headers: 4.6 s for a synthetic 33 M on the M5 Max, single-threaded, under a mutex), which the
  freshly deployed server had not built yet. It is built while the index loads now, and `timing` states the seed
  phase (`seed_phase_ms`, `seed_fetch_ms`, `label_resolve_ms`). The cache did slow first reads, though: keeping every
  row of every path copied each wide tuple row of a path (one allocation per column with more than two
  coordinates), shown on a synthetic index of 8,000 labels with 24,000-coordinate rows — 1.6-6.2 times slower first
  reads than keeping fewer rows on paths of up to 1,000 rows, and a walk reading one row per call 16.9 s against
  2.2 s with the cache off, its generations churned. A read now keeps the rows asked for, the 8 after each, the
  checkpoints (distance to the anchor a multiple of 16) and the narrow rows (under 4,096 bytes: UHGG keeps every
  row, as before), and keeps tuple rows flat. On the synthetic index no request is slower than the pass-5 binary or
  the cache off beyond noise, first reads are 1.1-2.9 times faster than with the cache off, one-row walks 6-15
  times; budgeted reads keep exact costs (tested across rules). Staging requests to rerun: `budgets.*`,
  `callsize.*`, `resolve.*`, `walks.*.first`, `lookahead.*` — on a server that has served one header-named request
  before (or the build with the header index at load), interleaved with the previous build.
- **R8, attributing an overrun.** A warm tuple walk ran 6,136 ms on 5,000 and nothing said which piece was late.
  `timing.deadline` states the seed's longest uninterruptible piece (kind, rows, coordinates mapped in it) and
  what stopped the walk how long after its deadline. A read maps its rows' coordinates inside its piece, so the rate
  a chunk is sized by includes the mapping.
- **The efficiency pass, audited (part B).** Path reuse: a row's logical charge is its standalone charge whether
  decoded, cached or shared (holds; tested across retention rules, capacities 0 / 64 KiB / 128 MiB, batches, chunks
  and a seed warm after a copy of itself); the physical work is in `timing.path_cache` only (fixed); the cache lives
  in the label cache's allotment under a memory budget, the lookahead admitted as without it (holds); reconstruction
  suffixes are kept as full rows at the requested rows, their successors, the checkpoints and the narrow rows, not
  as raw diffs (stated). `/resolve`'s single decode retained every present row until the profile was done; it now
  decodes in bounded batches, a discovery accumulating its profiles as it reads (fixed in its second version, below:
  the first kept a prefix of 256 MiB for a profile pass that decoded the rest again; the ranking and ties
  unchanged). Header ranges: the run mapper kept rank and select per run where every sequence
  matters; a query of named headers now filters by its own sequences' ranges instead of mapping every coordinate
  (fixed, as the review advised; pushing the ranges into decoding is 3c). `/resolve` with explicit labels, reported
  1-5% slower: not reproduced locally beyond noise (mini refseq and the synthetic index, three builds; the run
  mapper costs at most two extra selects per sequence change), to be rechecked on staging with repetitions.
  `WalkerDeadlineChunks` run on a virtual clock (exact bounds, no failure under load). The delivery reserve's
  calibration text is corrected: on a warm SRA server the first tree and full attempts were cut at 3-5 s of 40 s at
  levels 3 and 4 alike (the tree measuring 115-129 against the 30 assumed); per-detail starting ratios were not
  adopted on one locus.
- **`schema_version`** is the one request schema accepted (1, no window); a future request change adds
  `request_schema_versions` and keeps accepting 1 for a stated window; clients key features on `feature_level`.
- **Order of the later work, as the review recommends.** Stage 3b's `/resolve` budgets come before the remaining
  annotation formats. For 3c: the 64-column crossover is a tuning point, not a threshold (count distinct columns,
  BRWT ancestor sharing, dependency length, selected tuple volume); the charge model counts the coordinates examined,
  shifted or cancelled during reconstruction, not only the survivors; header-range filtering inside the decode
  shifts the ranges along the dependency path, never one absolute interval at every row; a tuple column's presence
  is not the XOR of Boolean nonemptiness; and interruption checkpoints inside dependency traversal, BRWT probing and
  tuple reconstruction come before a finite `max_uninterruptible_ms` is advertised. Measure first: first-touch
  batching and bounded prefetch of known dependency ranges, with major/minor fault counts; anchor density (path
  length percentiles, tuple volume, overlap; logical work changes with the representation and stays attributed to
  its identity); parallel decode only after reuse, with global concurrency bounds, commits in walk order and
  joined workers; cross-request caching with complete index identity, byte limits and unchanged logical charges;
  benchmarks that compare completed work, or certified reach and consumed work for deadline-limited requests.
- **The review of these fixes** (before release, of the working tree; seven findings, all fixed, each with a
  regression test from the reviewer's probe):
  - *A finished state did not hold against a replay (moderate).* The release rule lets a ledger release on a
    finished state, but a finished attempt's id was refused only while retained: with `retention_s` 1, or
    `retention_count` 1 and one later finish, the identical request (same `attempt_id`, `not_after_ms`,
    `expect_server_instance`) ran again. A finished attempt sent with `not_after_ms` is now held, its id refused,
    until its `not_after_ms` + skew (within `tombstone_max_s` of the finish, and extended by a refused copy's later
    `not_after_ms`), past its retention if need be, among the tombstones and on both clocks. Held attempts are never
    dropped early, so `retention_count` does not bound them; `tombstone_max_s` does (stated). *(Review of
    2026-10-06, C30: not quite — every refused copy, an identical replay included, extends the hold to
    `tombstone_max_s` from its own arrival, so the finishes and the refused copies within `tombstone_max_s` bound
    them; the texts say so from feature level 6, the behaviour, conservative for replay safety, is unchanged.)* Without `not_after_ms`
    nothing bounds a replay: the release rule states that a finished state then assumes no copy arrives after the
    attempt left retention (and with retention 0, after it finished).
  - *Off by one at the expiry (minor).* The tombstone ended when the wall clock read `suppressed_until_ms` exactly,
    while the strict `not_after_ms` check still admits in that millisecond: with skew 0, 27 of 40 covered copies
    ran. The hold is live through `suppressed_until_ms` inclusive (the steady hold one ms longer), which also makes
    the stated clock-step limit exact (a backward step of **more** than the skew).
  - *A copy without `not_after_ms` was told it was covered (low).* The 409 fell back to a cancel's `not_after_ms`;
    it is judged against the copy's own now, `covers_admission: false, no_not_after_ms` without one.
  - *Texts (nit).* The SPEC's "tombstones are kept `retention_s` each", the help of `--traverse-attempt-retention`
    (tombstones are never dropped oldest first: beyond the count a cancel is refused), and the start-up refusals,
    which now name each option's own range for every bad value.
  - *One manifest for a directory of bundles (major).* Loaded files are matched to manifest entries by base name,
    and nothing required base names to be distinct: a manifest listing `A/…` and `B/…` with shared names passed for
    both bundles, and two single-index servers stated one `index_fp` and `index_meta_fp` with different accessions
    (R3's dedup only covers pairs of one multi-graph server). Duplicate base names are refused, by the server, the
    CLI and `index_manifest.py --verify`.
  - *A manifest listing an optional file the pair does not load (minor).* The check ran one way (loaded ⊆ listed):
    a bundle without its `.seqs` beside the symlinks, or a server with `--no-coord-mapping`, stated the `index_fp`
    of one that loads it, and for a mask or a Bloom filter both fingerprints were equal. A manifest may not list an
    optional inventory file (`.edgemask`, `.bloom`, `.seqs`; since §29 the `.seqs` only, a mask or a Bloom filter
    never) the pair does not load; `--extra` files outside the inventory stay allowed.
  - *`/resolve` held every row anyway, and decoded most twice (moderate).* The first R6 fix kept a 256 MiB prefix of
    a discovery's rows and decoded the rest again in one call for the profile pass: the same peak (1.66 GB on the
    synthetic 3 kb locus of 8,000 wide columns) and 40–56% more CPU. Now every `/resolve` decodes in batches (64
    rows first, then about 64 MiB by the widest row of the batch before, at most twice its rows and 4,096) and holds
    one batch; a discovery accumulates each label's k-mers, runs and trace chain while it reads (no profile pass),
    keeping a repeated k-mer's row for its later occurrences within 256 MiB; explicit labels prime their query with
    each distinct row once. On the synthetic index the 0.3–3 kb queries take 3–17% less CPU, header presence −7% to
    +14%, at 144–167 MB instead of up to 1.66 GB; a 9.8 kb genome repeating the locus three times is faster for
    column presence, within 8% for explicit labels, and 24–48% slower for header and trace discovery (its repeats
    beyond 256 MiB decoded again) at about 370 MB instead of 4.8 GB (medians of 5 on a loaded machine). A discovery indexes its column labels by a 4-byte slot per column instead of an accumulator
    per column. Tested: the profiles equal across batch sizes, byte targets and kept bounds, on a query and on one
    repeating two thirds of it, and equal to the same labels given explicitly; no read exceeds its batch.

# 25. The review of levels 4–5 *(review 6, 2026-10-05; level 5 kept)*

An external review of levels 4 and 5 (pinned at `dcc0cebd`, the plan at `c9e05628`) confirmed the four blockers of
the review of pass 5 fixed and found no index-identity collision, traversal-evidence error or MGT v1
incompatibility. It found five implementation issues (findings 1–5, all P2), one qualification of the lease
contract (finding 6) and one design disagreement with the plan (X-R3). Findings 2, 4 and 5 are the library's
(`metagraph.traverse`) and are fixed by its half of this pass; this section notes the server's (1, 3, 6), the plan
revisions and a corrected benchmark statement. No request or response field is added and `feature_level` stays
5: what changes is physical work (finding 1), values in `timing` (finding 3) and two capabilities texts
(finding 6). *(Review of 2026-10-06, X1: true of these fixes, not of the builds they were made on: `7aaee760`
had put D3 and the first `output.coordinates` fields into the server, unconditionally, while it stated level 5,
so `7aaee760`..`67bef367` answer a derived seed whose time budget runs out mid-derivation with a walk where
`dcc0cebd` failed it, under the same level digit. None was deployed; none may be. Level 6 states both, §26.)* Each fix has a regression test made from the reviewer's probe (SPEC T53); untimed responses are
byte-identical to `7aaee760`.

- **The path cache allocated before admission (finding 1).** `RowDiffCache::insert` copied (or flattened) the
  whole row, then computed its size, read the shared bound and evicted. A cache bounded at 1 KiB given a row of
  1,048,576 columns allocated 4,210,688 bytes and refused the row, its bytes, peak and kept rows all 0; a small
  bound could thus make large, unreported copies over and over, and a full cache held its bound plus the copy
  (two 1.5 MiB rows under 4 MiB: 4.5 MiB during the insert of a third). Now the stored size is computed from the
  row (`StoredRow::bytes_of`: one buffer for a binary row, three exact buffers for a flat tuple row — the same
  value as the stored copy's `bytes()`, asserted), the bound is read and the older entries are evicted first, and
  the row is copied only when admitted; a refused row is counted (`rows_refused`, internal). The peak the cache
  states (`timing.path_cache.peak_bytes`) is what it held, insertion included: the copy is made after the
  evictions that make room for it, so no insert holds more than its end. What is kept, evicted or hit is
  unchanged (the same size decides). Measured with jemalloc's thread counters on the reviewer's probe: 0 bytes
  allocated by the refused insert (4,210,688 before); the kept third row adds at most 0 bytes to the cache's
  3,147,008 B (1,572,864 before, 4,719,872 B held at once under a 4,194,304 B bound).
- **The deadline record missed the seed phase's setup (finding 3).** A valid request naming 2,500 headers with a
  shared 1,024-character prefix stated, on a warm server, a seed phase of 126.37 ms, label resolution 3.00 ms, and a
  longest piece of 0.297 ms: the duplicate check of the names compared each with every earlier one (about 3 GB of
  `memcmp` on that prefix), neither it nor the resolution was a piece, and the first head's checkpoint did not count
  what came before it. Now:
  - the seed phase is cut into pieces wherever another piece (a read or chunk, a k-mer mapping, a derivation step)
    begins or ends, and every span between two of them is a **`setup`** piece (`DecodePacer::open_setup` /
    `close_setup`, the walker's `begin_setup` / `end_setup`); the phase ends at the walk's first checkpoint, or
    at its stop or end when it reaches none (a radius of 0, a D3 seed, no arm), and is flushed on every way out,
    a failed seed's included;
  - the walk after its last reading of the clock, to its stop or end, is a last `head` piece (it was no piece:
    the last head's remainder and the last level's merge and beam), so the pieces cover the seed's time from its
    start to its result;
  - `seed_phase_ms` runs from the seed's start to that end (before: the validation or derivation alone), and
    `seed_fetch_ms` counts the reads in it (an annotate root's read now among them);
  - the duplicate checks are hash sets in the request's order — the seed labels' by name (a resolved label
    carries the name it was given), the extra labels' by target (the same kind and column, and for a
    header its sequence) — refusing what the scans refused, the first duplicate in order and after resolving the
    names before it (tested, including an unknown name before a later duplicate).
  On the reviewer's request from the CLI: seed phase 13.9 ms against 139.1, elapsed 14.2 against 143.5, longest
  piece `setup` 9.8 ms against `read` 1.28 ms. The 9.8 ms span left is the seed id over the 2,500 names (their sort
  and FNV-1a over 2.6 MB, part of its definition), linear in the names' bytes and now stated rather than hidden.
  `usage.observed_max_uninterruptible_ms` stays the longest decoded piece (reads only, as stated in §6.8 of the
  SPEC; head pieces join them at feature level 6, §27); a setup piece is in the deadline record only. Timing only: no untimed byte changes, the walk's logic
  is untouched (the pieces are notes; `head_clock_ms_` is set at the seed phase's end instead of at the first
  checkpoint, which only adds a piece of about 0 ms there).
- **The finished-state release across a restart (finding 6).** The capabilities required instance pinning for
  early release on a tombstone but described a finished request within its hold as never run again. The hold
  lives in the process's memory: a request finished with `not_after_ms` 60 s ahead, without
  `expect_server_instance`, replayed after a restart of the same endpoint before then ran again (200, 170 work
  units), while the pinned copy was refused (409 `instance_mismatch`). The texts now say so, as the SPEC's §5 and
  §10.3: the hold is the process's; a finished state is replay-safe only for an attempt sent with
  `expect_server_instance` (equal to the instance that finished it) and with `not_after_ms` within
  `tombstone_max_s` of the finish, so a ledger releasing on a finished state before `not_after_ms + skew +
  bound_ms` pins the instance, as for a tombstone; an unpinned attempt's finished state assumes that no copy
  reaches a restarted process (or another server at the same address) while it could still be admitted there.
  The `instance` text names finished requests beside cancelled ones. The library's `release_verdict` already
  lists `sent_without_expect_server_instance` among a finished release's assumptions; its recorded answers keep
  the old texts (they are not compared).
- **The benchmark statement, corrected.** The review request's table (`REVIEW-REQUEST-traverse-level4-5.md`, not
  edited) says that walks cut off by their budget "reach about 2× more bases". The reviewer's reanalysis of the
  saved staging results (A after-pass5 `e7d99e0c`, level 3; B level 4 `f667d775`; 54 requests each, run about 8 h
  apart) limits that: completed walks with identical work, B/A latency 0.722; **deadline-limited walks, B/A
  certified reach 1.465** (13 further, 2 the same, 0 less; consumed work 1.441); **the short-budget suite, B/A
  certified reach 0.251** among the comparable nonzero pairs (1 further, 1 the same, 3 less; work 0.566); completed
  lookahead, B/A latency 0.093. About 2× holds for some walks (ec_rplJ 2.0–2.45), not in general, and runs hours
  apart do not separate the code from cache state and host load: the synthetic cache evidence shows the
  mechanism, the staging attribution is not established. Level 5's start-up (the header index built before
  readiness) is to be measured on its own.
- **Plan revisions for the owner** (`PLAN-traverse-next-stages.md`, marked there; not code):
  - *X-R3's reconstruction bound.* The lazy-shift argument does not bound reconstruction work: lazy shifting
    avoids re-adjusting inherited coordinates, not re-scanning them, and an ordinary merge at each nonempty diff
    copies the whole inherited window. A 256-row path with an anchor of 100,000 selected in-window coordinates and
    255 nonempty singleton diffs emits about 25.5 million coordinates in vector merges while the price counts
    100,255 — below the deep-term threshold D, and missed by the empty-diff test. Revised (§2.2 revision 2):
    diff-first reconstruction — per column the path's windowed diffs (down to the first clear marker) merged
    among themselves first, then once with the anchor's window — with its merge term Δ_w × ⌈log2 k⌉ charged (Δ_w
    the path's windowed diff coordinates per column, k their nonempty lists; not the slack D): in that case merge
    term 255 × 8 = 2,040, work bound and the price's coordinate terms 100,000 + 255 + 2,040 = 102,295; or as the
    fallback a deterministic merge-weighted charge; new tests with sparse nonempty diffs below and above D, with
    clear markers, and with cache and batch variations; a new question 3c-N10 (and whether the same applies to
    full reads, whose level-4 price counts stored entries).
  - *The 64-column crossover* stays provisional: the sweep covers selected-coordinate density, headers per
    column, path length and cold and full caches; the decision stays deterministic and `decode: "full"` stays.
  - *`max_uninterruptible_ms` stays null after level 8* (3c-N5): checkpoints bound index operations, page faults
    and scheduling leave the time open. The review request's "so `max_uninterruptible_ms` can become finite" and
    §24's "come before a finite `max_uninterruptible_ms` is advertised" are superseded; the SPEC's §6.8 and §10.3
    say it.
  - *Numbering*: these notes take §25, so the stages' sections in the plan move up by one (coordinates §26 …
    3b-2 §30), their versions unchanged.
- **Noted, not acted on here** (operations, the reviewer's recommendations): eager header indexing stays the
  default for header-capable servers (an explicit column-only opt-out could save 1–2 GB and start-up work if no
  header-dependent operation can trigger lazy construction; one-shot CLI invocations deserve their own policy);
  measure memory pressure (file and anonymous memory, refaults, reclaim, pressure stalls, faults, I/O per request)
  before raising the container's memory or changing read-ahead, then compare `MADV_RANDOM` against normal advice
  separately, in serial paced trials that abort on agreed sibling-latency or pressure thresholds.
- **After the verification of these fixes** (three low findings, all fixed; tests and documents only, no source):
  - *The GCC 13 gate covers `unit_tests`.* It compiled the 112 `src/` files only, so a gtest assertion under an
    unbraced `if` in the new `RowDiffPathCache.ARefusedRowIsNeverCopied` passed AppleClang and failed GCC's
    `-Wdangling-else` (an `EXPECT_*` expands to an if/else). The gate now compiles the 75 `tests/` sources too,
    with the target's own options as CMake gives them: gtest and gmock as system includes (googletest declares
    them so), `-Wno-uninitialized` after `-Wall`. It found 22 older sites in this branch's test files (19
    dangling-else, 3 range-for loops binding `const std::string &` to string literals) in `test_walker`,
    `test_graphlet_codec`, `test_label_oracle`, `test_label_oracle_budgeted`, `test_trie` and `test_trie_cases`;
    all are braced or iterate over `const char *`, and no test checks anything different. 112 + 75 files pass.
  - *SPEC §10.3 agrees with §6.8*: its `deadline_check` text now says `max_uninterruptible_ms` stays null after
    stage 3c-ii (3c-N5).
  - *The plan's merge term has its own symbol*: Δ_w × ⌈log2 k⌉ (Δ_w the path's windowed diff coordinates per
    selected column down to the first clear marker, k their nonempty lists), no longer D, which is the slack of
    revision 3. In the reviewer's case the merge term is 2,040; the work bound and the price's coordinate terms
    are 102,295 (corrected above, in the plan's revision 2, 3c-N10 and X-R3).

# 26. Record coordinates as built, and feature level 6 *(v5.8, 2026-10-05)*

Record coordinates (§18, the owner's decisions C1–C12 and the program plan's C-N1..C-N9, X-18.2, X-C8, X-C12;
`PLAN-traverse-next-stages.md` §2.1) are implemented in two chunks: W1 (increments C0–C2, `7aaee760`: the walker's
side tables and the wire contract, frozen at C2, with D3) and W2 (C3–C8: the delivery reserve's coordinate share,
the capabilities, level 6, the real-index tests and measurements, the library, the documents). SPEC §5, §6.8, §7.0,
§7.1, §7.5, §10.3 are normative; this section records what was built against §18, why, and what was measured.
`feature_level` is **6**.

## 26.1 The contract, as frozen at C2

- **Request** (`strategy.output`): `coordinates` (boolean, default false; false is byte-identical to omitting it)
  and `max_coordinate_occurrences` (an integer in [1, 2^64 − 2] or `"unlimited"`, default 16; 400 without
  `coordinates: true` — decision C-N6 —, inert under `support: kmer` and in annotate mode — C-N7). Echoed, the cap's
  default included, only when `coordinates` is true.
- **Response**, per seed in every detail (a graphlet's JSON summary, never its body): the block, or `null` with
  `coordinates_reason` (`index has no coordinates`, `support kmer`, `no traversal`, `partial derivation`, in that
  order). The block: `kind` (`record` | `column` | `mixed`, C3), `k`, `max_occurrences`, `complete`,
  `runs_lower_bound`?, `seed` (per seed label its occurrences of the seed) and `arms` (per requested arm one entry
  per run, in `runs` order: `run`, `label`, `from_bp`, `to_bp`, `occurrences`, `occurrences_total`? when cut,
  `chains_ended`? when > 0 (C5, C-N1), `lower_bound`? when true (C-N2, X-18.2)). Intervals 0-based, half-open,
  forward strand, `[c + k − L, c + k)` right and `[c, c + L)` left from the chain's coordinate c at the run's last
  node — the c₀ form of §18.2 rewritten, verified on 3,631 real runs before the build and on 10,696 occurrences of
  9,498 runs by `MiniRefSeq.CoordinatesAgainstTheSourceRecords` now (5,508 of them switch-entered, 48 lower
  bounds, 27 of them strictly below the string count; 326 lists cut; every recorded switch request's seed walked
  to its 3,000 bp — the review of W2 found the blaNDM and first repeat seeds' switch cells at the default 1,000 —,
  pinned: for each of the 7 seeds, a run of its switch cells ends past 1,000 bp).
- **Column positions** are the column's k-mer index space (record i's k-mer j is `offset_i + j`, `offset_{i+1} =
  offset_i + len_i − k + 1`), so **a record's last k − 1 bases share their numbers with the next record's first
  k − 1 positions** (review of W1, finding 6): an interval there is attributed to one record only with the record
  lengths, which are not stated (C10). Stated in SPEC §7.1, in the probe's `coordinates.rule` (C3) and exercised by
  the positional oracle (3 such intervals in its column switch cells, numbers in a record's last k − 1 bases whose
  bases are the next record's first ones, not the earlier record's).
- **Limitation** `coordinates` (seed level, only when a list was cut; in no outcome class, C2 kept), knob
  `output.max_coordinate_occurrences`, observed the largest true count, `lists_cut` the number of cut lists; MGT:
  exactly one `K * coordinates …` record more, Z one higher. **Action** `drop_coordinates` right after
  `use_graphlet` / `drop_sequences` wherever those are offered, for a request with coordinates (block or null form).
- **Memory**: each run's side-table entry charged once at its creation at min(chains, cap) occurrences, a seed
  label's with the depth-0 state, in the requested detail; `LabelRun`, `Entry` and `Item` unchanged (the layout
  trap: their sizes are part of the model). **Work**: none — every coordinate was charged one unit with its row.
- **D3** (the owner's decision R21 (3), all requests): a derived seed whose time budget runs out after 1 ≤ j < n of
  its k-mers is delivered as a partial walk at depth 0 with the set of the j k-mers read, a `derivation`
  limitation first (`cause: time_budget`, observed j) and `label_evidence: qualified`; j = 0, or a superset that
  would fail the seed any other way, fails as before (SPEC §7.0). The library accepts it from W2 (its `check_rules`
  flagged the qualified class against "complete" limitations before).

## 26.2 The delivery reserve's coordinate share (C3; revisions 3 and 8, decision C-N4)

- The walker's account has a **coordinate share** `ResourceAccount::coordinates`: the output's fixed part for
  coordinates (`DeliveryCosts::coordinate_fixed`, already in `fixed`: the block's skeleton with a cut list's
  limitation and K record, or the null form; counted again only in the share, nothing more is charged) plus every
  entry and occurrence as charged at its creation. It is published with the account at every level's end
  (`AttemptControl::progress(account, coordinates)`) and in the seed's meter.
- **The estimate**: E = ⌈(A − C) / ratio⌉ + ⌈C / 12⌉, `kCoordinateAccountPerTextByte` = 12 a **bound**, not an
  estimate: each part of the share is priced at least 12 times the most text it can write, whatever the digits
  (`GraphletCoordinates.CoordinateAccountBoundsItsText`, detail `full` / `graphlet`):

  | part | account (B) | widest text (B) | ratio |
  |---|---|---|---|
  | an occurrence (two 20-digit numbers, brackets, comma) | 872 | 44 | 19.8 |
  | a run's entry (every optional member at its widest) | 3,680 | 221 | 16.7 |
  | a seed label's entry | 1,728 | 79 | 21.9 |
  | the block's skeleton, a cut list's limitation (20-digit count, `"unlimited"`), `drop_coordinates` (JSON detail) | 11,628 | 597 | 19.5 |
  | the same with the K record, the Q token and the counts' digits (graphlet) | 15,151 | 960 | 15.8 |
  | the null form with its longest reason and the tokens (graphlet) | 1,584 | 106 | 14.9 |

  Without coordinates C = 0 and E is the formula before the split, to the bit.
- **The ratio sample** leaves both sides of the share out: (A − C) / (text − the coordinates' text), the latter
  counted **exactly** from digit counts and fixed punctuation (`coordinates_text_bytes`: the block or null form, a
  cut list's limitation, `drop_coordinates`, and in a graphlet the escaped K record, the Q token and the digits they
  add to Z, `graphlet_lines` and `graphlet_bytes`; `compact_json_size` counts compact JSON without writing it),
  and only on a rest of at least `measured_text_bytes`. So an attempt with coordinates measures exactly what the same
  walk without them does — `GraphletCoordinates.AttemptsMeasureTheSameRatioWithCoordinates` (the server's path on
  seeds of 1.3 MB: equal doubles, cut and uncut, full and graphlet), `.CoordinateTextIsExactAndBoundedByItsAccount`
  (142 results with 20-digit positions, memory stops and null forms: text − coordinates' text = the stripped
  result's text; A − C = the opt-out account), `MiniRefSeq.CoordinateShareIsExact` (60 real results) — and **W1's
  interim** (attempts with coordinates kept out of the server's measurements) is removed: they feed the measured
  ratios again, which no longer depend on whether such a request came first
  (`GraphletAttempt.CoordinatesLeaveTheServersRatioUnchanged`).

## 26.3 Capabilities and level 6 (revision 8)

The level-6 documents differ from level 5's exactly in: `feature_level` (both GET routes and every response's
capabilities); `attempts.delivery_reserve.rule` (the coordinate share and the exact sample) and the new numeric
`attempts.delivery_reserve.coordinate_account_per_text_byte: 12` (both GET routes); `deadline_check.rule`, one
sentence (both GET routes: "no time bound on one piece is stated (max_uninterruptible_ms: null, and it stays null:
checkpoints inside reads bound the index operations of a piece, not its wall time …)" in place of "no bound on one
row exists before stage 3c", which implied a finite value after 3c — decision 3c-N5); and the probe's new
`coordinates` block (`supported` = `supports_trace`, `knob`, `cap_knob`, `max_occurrences_default` 16, `kinds` this
index can report, `limitation`, `action`, `output_bound` — which names the probe's `max_memory_mb` and
`max_work_units` rather than a value, since a server maximum is the budget of every request without one —, `rule`).
The per-request capabilities change only in the digit. The compact probe grows from 15,908 to 19,490 bytes on
mini_refseq (the block 2,799, the reserve's rule 589 and the deadline rule 141 more; `GET /capabilities` 13,313 to
14,081), past the 16 KiB default ceiling the MCP tool `traverse_capabilities` had at level 5; the library raises
it to 32 KiB (W2's library part). Not changed, owner decision pending (review 6): whether `setup` pieces count toward
`observed_max_uninterruptible_ms` (they are in the deadline record only).

## 26.4 The owner's decisions, as built

| ID | Decision | As built |
|---|---|---|
| C1–C12 | §18 and §21 | as §18 with the amendments of §18.2 (v5.8) |
| C-N1 | `chains_ended` as drafted | chains continuing on no followed path of the lineage; partitions among clones not counted; clones inherit |
| C-N2 | (a) mark lower bounds | `lower_bound` per run, `complete: false`, `runs_lower_bound`; M1: 0.56% of runs, all switch-entered — (b) not needed now |
| C-N3 | `Coord*` subclasses | the library (W2) |
| C-N4 | fixed bound 12, exact sample | §26.2 |
| C-N5 | no server cap; true bound stated | the probe's `output_bound`; ops note O2 (`--traverse-max-memory-mb` on staging) |
| C-N6 | 400 for a cap without coordinates; the library strips it | server as frozen; library W2 |
| C-N7 | a cap under kmer support inert | as frozen |
| C-N8 | no coordinates in GFA | library |
| C-N9 | `compare()` notes, never clips | library |
| X-18.2 | §18.2 amended | §18.2 (v5.8) |
| X-C8 | auto-on under a memory budget only if the D4 gate passes | M1 by regime (§26.5): holds for header labels at refseq's chain counts and column labels below 16 chains a run (median −3 to −4%); does not hold for column labels with 16 or more chains a run (refseq33m's taxid columns: −13 to −16% projected, failures at depth 0 on the fixtures) nor for header-heavy seeds — returned to the owner (§26.9) |
| X-C12 | open token sets | §14.1 (v5.8), SPEC §7.5 |
| R21 (3) | D3 | §26.1 |
| R21 (4) | a displayed walk follows the merge parent carried by the most labels | §26.6: everywhere except `tree` / `full` under a memory budget, which keep the arrival order so that no budget's stop moves (owner's choice, §26.9) |
| R21 (5) | `compare()` above 64 labels: "qualified" can come from the per-node label limit | the library's docs (W2) |

## 26.5 Measurements

**M1 — the D4 gate, by regime, and C-N2.** Two measures per trace cell (`gate_cell` in `benchmarks/traversal/bench_traversal_measurements.cpp`): the
depth at the stop (`complete_to_bp` over both arms) under memory budgets of 50% and 75% of the cell's own opt-out
peak, exact bytes, with against without coordinates at cap 16; and — since that ratio is relative to the opt-out
walk and drops a cell whose opt-out walk fails at depth 0 too (the review of W2) — what **completing** takes: the
smallest whole-MiB budget that holds the walk with and without coordinates, and what the walk with coordinates does
at the budget that completes it without them. Cells, all to 3,000 bp, details `full` and `graphlet`:
mini_refseq (`BM_TraversalCoordinatesDepthAtTheStop`: 7 seeds — blaNDM both ways, three 200-bp windows of
its carriers, two repeat windows — × 8 strategies — branch limits 0 and 2, switch costs 0.5 and 1 at limits 0 and 2
within a loss budget of 2, column labels at limit 2 with and without a switch); its column cells again under the
**refseq33m projection** (every run's and seed label's list priced at the cap of 16 occurrences, as refseq33m's
taxid columns of many genomes give; mini_refseq's taxid columns reach 12 chains a run at most); and the fixtures
of revision 4 (`BM_TraversalCoordinatesDepthByRegime`): `scripts/traversal/make_column_coord_fixtures.sh`
— `coord_lockstep`, 200 columns of one shared 3,000-bp sequence, 16 records each (one path, 200 runs of 16 chains an
arm), and `coord_divcol`, 100 columns of their own sequences around a shared 600-bp core, 16 records each (seeds in
the core: a split into a branch a column at each end of it; and in a flank: one column) — and the wide fixture
(below: column labels one run of 5,000 chains an arm, header labels 5,000 seed labels, a run of one chain each).

| regime, cells (× full, graphlet) | chains a run | median depth with / without at 50% (full, graphlet) | at 75% | budget completing with / without (median; max) | with coordinates at the budget completing without |
|---|---|---|---|---|---|
| mini_refseq header labels, 42 | ≤ 6 | 0.95–0.97; 0.81, 0.52 single-lineage switch cells (limit 0) | 0.96–0.98; 0.94, 0.92 | all 56 cells: 1.00 full, 1.02 graphlet (1.33; 2.0, limit-0 switch cells) | all 56 cells: complete 29 / 24, stopped short 27 / 30, failed at depth 0 0 / 2 (full / graphlet, all strategies) |
| mini_refseq column labels, 14 | ≤ 12 | 0.96, 0.96 (switch 0.91, 0.94) | 0.98, 0.97 | 1.00–1.04 (1.07) | (in the row above) |
| the same, refseq33m projection | 16 | **0.87, 0.84** (col_b2 0.91, 0.90; with a switch 0.73, 0.68) | 0.91, 0.91 (switch 0.80; worst 0.28) | 1.06, 1.07 (1.13; 1.18) | complete 1 / 0, stopped 13 / 14 of 14 |
| `coord_divcol`, core seed | 16 | **0.85, 0.87** | 0.90, 0.90 | 1.06 (121 → 128 MiB full, 108 → 115 graphlet) | stopped at 2,757 / 2,743 of 6,000 bp |
| `coord_divcol`, flank seed (one column) | 16 | both at depth 0 | **0.66** full (398 → 262 bp); graphlet both 0 | 1.0 (3 MiB) | complete |
| `coord_lockstep` (200 labels, one path) | 16 | the opt-out walk itself fails at depth 0 at 50% and 75%: not measured by this ratio | — | **2.1 full (13 → 27 MiB), 3.8 graphlet (5 → 19 MiB)** | **fails at depth 0** |
| wide fixture, column labels | 5,000 (cap 16) | both fail at depth 0 | — | 1.0 (4 / 3 MiB) | complete |
| wide fixture, header labels (5,000 seed labels) | 1 | both fail at depth 0 | — | **1.28 full (235 → 301 MiB), 3.1 graphlet (38 → 119 MiB)** | **fails at depth 0** |

- **The gate (a median drop of at most 10%, X-C8) holds only where the coordinates' share of the account is small:
  header labels at refseq's chain counts and column labels below 16 chains a run** (mini_refseq: median −3 to −4%;
  every cell no deeper with coordinates, nearly every one a little shallower; the outliers are the single-lineage
  switch cells, limit 0, where two `full` cells fail at depth 0 at 50% with coordinates and not without).
- **It does not hold for column labels with 16 or more chains a run — refseq33m's taxid columns — nor for
  header-heavy seeds.** Projected to refseq33m's chain counts the median drop is 13–16% at 50% (27–32% with a
  switch); on the fixtures 13–15% (the core seed) and 34% (one column at 75%); and where many labels walk together
  the seed needs 2–4 times the memory to complete with coordinates and **fails at depth 0** at the budget that
  completes it without them (`coord_lockstep`; the wide fixture's 5,000 header labels). A library that turns
  coordinates on under a memory budget there turns a complete answer into a refused seed; `drop_coordinates` is
  offered, but only after the failure. X-C8 is returned to the owner with these numbers (§26.9).
- **C-N2**: unbudgeted, 21 of 3,780 runs (0.56%) are lower bounds, all switch-entered (21 of 1,971, 1.1%), in 10 of
  56 cells — 8 in `sw0.5b2`, 6 in `sw1b2`, 7 in `col_sw0.5b2`, none without a branch allowance. Marking (a) is
  enough; the exact chain set (b) is not needed now.

**M2 — the wide fixture** (`BM_TraversalCoordinateReserveOnTheWideFixture` on
`scripts/traversal/make_wide_coord_fixture.sh`: one 3,000-bp sequence in 5,000 records under one column, labelled
`wide.fa` — annotated on the relative name since the review of W2, which found the column's label, and with it the
column rows' bytes, to be the absolute path of wherever the fixture was built —, row-diff anchors every 1,000 rows;
a 100-bp seed in its middle walked to the records' ends; column labels give one run an arm with 5,000 chains,
header labels 5,000 runs an arm with one chain each; build times from the first, less loaded run):

| labels, cap, detail | text (coordinates) B | account (coordinate share) B | build ms with / without | coordinate estimate / their text | whole estimate / text | estimate before the split / text |
|---|---|---|---|---|---|---|
| column, 16, full | 14,608 (1,320) | 1,880,896 (62,572) | 0.15 / 0.24 | 3.95 | 4.51 | 0.94 |
| column, 16, graphlet | 8,990 (1,612) | 1,639,593 (66,095) | 0.10 / 0.16 | 3.42 | 4.11 | 0.86 |
| column, unlimited, full | 291,112 (277,824) | 14,919,040 (13,100,716) | 3.1 / 0.13 | 3.93 | 3.96 | 0.37 |
| column, unlimited, graphlet | 285,202 (277,824) | 14,677,737 (13,104,239) | 3.0 / 0.08 | 3.93 | 3.94 | 0.24 |
| header, 16 or unlimited, full | 5,076,957 (984,563)¹ | 222,416,410 (58,531,628) | 50 / 38 | 4.95 | 2.04 | 1.09 |
| header, 16 or unlimited, graphlet | 1,378,634 (984,563)¹ | 83,242,348 (58,535,151) | 14 / 2.7 | 4.95 | 3.90 | 0.96 |

¹ At cap 16. "Unlimited" writes 9 B more, all in the coordinates' text (the block's `max_occurrences`, `"unlimited"`
in place of `16`); the accounts are the same, since each run holds one chain and no list is cut at either cap.

- **D5 holds**: the coordinate share's estimate is 3.4–5.0 times the coordinates' real text, the whole estimate at
  the configured ratios (30 / 50) 2.0–4.5 times the real text. Estimated as before the split — the whole account at
  the ratio the rest of the output measures — the text of a coordinate-heavy seed was understated up to 4.1 times
  (column, "unlimited", graphlet: 0.24), which is what the split fixes.
- In every cell the account and the text less the coordinate share are exactly the opt-out walk's. The walk itself
  takes the same time (0.24 s column, 1.3–1.5 s header); building the block costs 3 ms for 10,000 occurrences.
- The reserve stays below the floor (allowance / 2) for these sizes, so the walk-until of a one-seed attempt at 30 s
  is the floor's (35,000 ms) with and without coordinates; the reserve itself grows with them (header, full: 2,551 ms
  against 1,819 ms; column, "unlimited": 1,173 against 1,009 ms).

**M3** (done before the build, plan §2.1): about 441 B a seed on the 12 recorded staging trace responses. **M4**
(staging, after the deploy, the owner): `bench_traverse.py --coordinates` on the refseq33m panel.

## 26.6 The displayed parent at a merge (the owner's decision R21 (4))

- **Where it was chosen**: `Walker::merge_level` made the first head to arrive at the node the merged segment's
  first parent (`parents[0]`); everything displayed follows first parents — a path's `segments` chain
  (`walk_path_leaf_first`), its spelled bases (`spell_path` of the tests) and continuation (`make_continuation`, the graphlet's
  `C`), the end labels' `route_bp` (a label taken from a later parent is routed from the merge) and the
  reconverge event's order. The graphlet writer, the JSON writer and the library read `parents[0]`; no other place
  chooses a displayed parent.
- **Now** the first parent is the one carried by the most labels: in constrain mode the labels whose lineages its
  head brings into the merge node (`Item::state`), in annotate mode the fewest labels present at a node of the
  parent's own segment (every head at the node holds the node's own labels, so the parents' own bases tell the
  routes apart); ties keep the arrival order. Under trace support nothing merges, so coordinates are untouched.
- **One exception keeps the arrival order: `tree` and `full` detail under a memory budget** (the request's or the
  server's). There every head reserves the delivery of its path's segment chain (`DeliveryCosts::chain_entry`,
  exactly what the JSON writes), whose length follows first parents, so the displayed parent decides the account
  and with it the memory stop. With the majority first everywhere (the first build of W2), a stop moved where the
  majority's chain was longer — on the cached real requests 52 of 6,552 budgeted responses, and in the review's
  sample at 8 MiB in `tree` detail 10 of 294 mini_refseq requests, 1–4 levels shallower on an arm (twice 30), never
  deeper — while the owner's decision changes what is displayed, never the depth a budget certifies. The other
  ways out were worse: charging below the output breaks the memory bound; charging the longest parent's chain at
  every merge makes the account independent of the display but moves more stops away from level 5. So where the
  account depends on the displayed chain, the display stays as it was at level 5; `graphlet` and `summary` detail
  charge no chain (the graphlet names a segment's first parent, not its chain) and follow the rule under every
  budget. A `tree` or `full` request under a memory budget therefore displays the arrival order at a merge where
  the same request without one displays the majority; which labels reach a leaf, with which loss and branches, is
  the same in both. A server maximum (`--traverse-max-memory-mb`, ops note O2) is the budget of every request
  without a smaller one, so on such a server every `tree` and `full` response keeps the arrival order and only
  `graphlet` and `summary` show the rule — the price of option (a), stated so that the owner weighs it (§26.9). Pinned: `WalkerTest.MergeKeepsTheArrivalOrderWhereTheBudgetChargesTheChain` (with the chain
  charged under a budget the same route is first whichever labels it carries; without either, the majority); on
  mini_refseq, win200_03 right at 8 and 2 MiB in `tree`, `full` and `graphlet` keeps level 5's `complete_to_bp`
  (2,782, 2,740, 3,779, 815, 800, 1,048). The owner may prefer the other
  trade (the rule everywhere, stops moving shallower in `tree` / `full` under a memory budget): it is the one
  condition in `Walker::merge_level` (`majority_first`), §26.9.
- **The annotate rule compares the parents' own segments, not the walks they display.** A parent that is itself a
  merged segment holds the union of its parents' labels, so after nested merges the rule can put first a short
  merged segment whose displayed walk upstream is carried by fewer labels than the other parent's. Pinned as the
  rule stands (`MiniRefSeq.AnnotateMergeRanksParentsByTheirOwnSegment`, the real cache's
  `mini_win200_02__annotate_merge`, right arm): at the merge at 288 bp the 5-bp merged segment from 283 (6 labels)
  is first before the 32-bp segment from 256 (5 labels), though its walk goes on through a 7-bp segment of 4 labels
  (276); the path's continuation is carried by label 3 alone where level 5's (through the 256 parent) was carried
  by 0, 1, 2 and 4. The review measured the net effect: over a fixed 10-bp window of the displayed walk, UHGG's
  annotate merges whose first parent is not a majority fall from 191 (arrival order) to 68 of 558, and of the
  continuations that changed 18 of 1,350 lose labels and 22 gain (mini_refseq 6 and 8 of 28). The library checks
  this rule as it stands (`derive.carried_labels`). A rule over equal stretches of the displayed walks (the last W
  bases through each parent and its first parents, labels present at every position) would rank these cases
  right; it changes the library's check with it and is left to the owner (§26.9).
- **What changes**, and only on requests with `on_reconverge: merge` where a merge's majority arrived later (and
  not in `tree` / `full` under a memory budget): the order of `parents`, `labels_via_parent` and the reconverge
  event's segments; the paths' chains, spelled bases and continuations through such a merge (a continuation's
  labels with them); the end labels' `route_bp` (the majority's labels now 0, the minority's the merge depth); and,
  on a tie of loss and branches, which parent's lineage continues (its run, its `labels_via_parent` entry, the
  closed run's segment) — in the graphlet the `G` parents and partition, `R` and `T` records. Without a memory
  budget, `tree` and `full` also state another account (`usage.memory`: the output they write is another one).
  **What does not**: which labels reach each leaf, with which loss and branches, the arms' `complete_to_bp` (under
  every budget), outcomes, limitations, `label_summary`'s direct and reach depths, work.
- **Fixtures**: `merge` (CLI) changed at its 62 merge — the G allele (b.fa, c.fa, both.fa) arrived second and is now
  first: a.fa is routed in at 62, b.fa only at 36 — and no longer shows `R:route_split`; `merge_ties` (new: the same
  locus with a.fa, b.fa, both.fa named, two labels on each allele: ties keep the arrival order) shows it; `annotate`
  changed at the same merge. The budget fixtures have no merge. The library's tests that encode the old `merge`
  document follow (W2's library part).
- **Checked** (T55): `WalkerTest.MergeIsSpelledThroughTheParentWithTheMostLabels` (the majority through either
  branch, constrain and annotate), the two pins above; every byte-identity difference against level 5 is a merge
  request and equal once the order-dependent fields above are canonicalised, and no budgeted stop moves (§26.8).

## 26.7 Shutdown on SIGTERM

`server_query` installed no signal handler; in a container it is PID 1, for which the kernel drops a signal with the
default action, so every `docker stop` of a staging deploy waited its timeout (120 s) and killed it. Now SIGTERM
and SIGINT write a byte to a pipe (the handler's only work, async-signal-safe); a thread reading it stops what is
in flight and the server (SPEC §10.3): every traversal — each `/traverse` has an `Attempt`, managed or not — is
stopped at its next poll as for a client gone (nothing is written for it, its connection closed: a partial result
would state a cause that is not true, "cancelled", and a ledger treats the attempt as unanswered), a `/traverse` or
`/resolve` arriving meanwhile finds its client gone, the HTTP server stops once they returned or after 3 s, and the
process exits as soon as it has stopped (`_Exit` from `run_server`, the log flushed), 3 s later at the latest even if
a handler (a `/search`, an uninterruptible read) has not returned. It does not return through `run_server`: its
destructors would join what still runs, and in single-index mode the server listens at once and answers 503 while
its index loads, so a signal can stop it with the load in progress — the first build of W2 then waited for the load
to end (the review: 19 s on SRA in RAM, on a cold staging disk the whole docker-stop timeout). Before the server
accepts connections it exits at once; a second signal exits immediately. Tested:
`TestTraverseAttempts.test_sigterm_stops_the_server_promptly` (idle: exit 0 within 5 s; walking a three-seed
attempt: the walk stopped, nothing written, exit 0 within 3.5 s) and `.test_sigterm_while_the_index_loads` (the
graph a FIFO nobody writes, so the load never ends: GET /capabilities `ready: false`, /traverse/capabilities 503,
then exit 0 within 5 s of the signal).

## 26.8 Verification

- **Opt-out byte identity** against level 5 (`67bef367`), decompressed, the `feature_level` digit and timing
  masked (the cached real requests, on the final build; each difference posted again to fresh servers and
  classified):

  | requests | identical | timing counters only | R21 (4) display only | memory stop moved | other |
  |---|---|---|---|---|---|
  | 2,184 unbudgeted, `full` and `graphlet` (mini_refseq 588, UHGG 792, SRA 804) | 1,807 | 26 | 346 | — | 5: 4 identical on the recheck, 1 the SRA attempt `sra_16s_PZ326290__beam20` (`full`) |
  | 6,552 budgeted, `full` and `graphlet` (300k work units, 4 MiB, 16 MiB × the three) | 6,054 | 6 | 492 (`graphlet` under the memory budgets, both details under the work budget) | 0 | 0 |
  | 2,184 unbudgeted, `tree` and `summary` | 1,929 | 36 | 215 | — | 4: 3 identical on the recheck, 1 the same SRA attempt (`tree`) |
  | 6,552 budgeted, `tree` and `summary` (2, 8, 16 MiB × the three) | 6,353 | 25 | 174 (all `summary`) | 0 | 0 |
  | 688 `/resolve` | 688 | 0 | 0 | — | 0 |

  Every differing request is `on_reconverge: merge` (most cached requests are: 502 of the 546 deterministic cells),
  and "display only" means equal once `parents` / `labels_via_parent` are paired and sorted (a label moved only
  between parents whose end sets hold it), the reconverge events' segments sorted, the paths' chains and
  continuations, the run tables, `label_summary`'s run lists, the end labels' run and `route_bp`, the graphlet's
  `R`, `C` and `T` extras and its partitions removed. **Under a memory budget the 5,460 `tree` and `full`
  responses are level 5's** apart from timing (R21 (4)'s exception, §26.6), and no stop moved in any detail; the
  first build of W2, with the rule everywhere, had moved 52 of the 6,552 `full` / `graphlet` budgeted stops, and 10
  of 294 in the review's `tree` / `summary` sample at 8 MiB on mini_refseq. The SRA attempt varies base against base (4 posts to fresh
  level-5 servers, 4 different responses: its walk-until depends on the load).
- **Opt-in** (the same requests with coordinates, stripped of what they add, against level 6's own opt-out):
  mini_refseq unbudgeted 586 of 588 identical (2 timing counters only), at 300k units, cap 1 and "unlimited" 588 of
  588; at 4 MiB 334 identical and 254 stopping earlier by memory, never deeper (`drop_coordinates` offered in 336);
  UHGG 779 identical and 13 timing only; SRA 776 identical, 27 timing only and the beam attempt above.
- **The frozen contract of W1** (the same requests with coordinates on both, level 5 against level 6): mini_refseq
  unbudgeted, at cap 1, at "unlimited" and at 4 MiB, and UHGG: every difference is a request of the opt-out run's
  R21 (4) set (112, 112, 112, 48 and 87 of its 87) but one, UHGG's `uhgg_rand100_05__beam20` attempt (`full`),
  which varies base against base too (4 of 4); the rest identical or timing only — no block, null form, echo,
  limitation or action changed.
- **Capabilities**: the diff against level 5 is exactly §26.3 (GET /capabilities: `feature_level`, the reserve's
  `rule` and `coordinate_account_per_text_byte`, `deadline_check.rule`; GET /traverse/capabilities: the same and
  the `coordinates` block).
- **Tests**: unit, the traversal filter 430 of 430 in Release and in Debug (416 before W2: + `ReserveCountsCoordinateText`,
  `CoordinatesLeaveTheServersRatioUnchanged` in place of W1's `CoordinatesRatioStaysTheAttemptsOwn`,
  `CoordinateTextIsExactAndBoundedByItsAccount`, `CoordinateAccountBoundsItsText`,
  `AttemptsMeasureTheSameRatioWithCoordinates`, `CapabilitiesBlockFollowsTheIndex`,
  `MiniRefSeq.CoordinateShareIsExact`, `.MemoryStopsDoNotDependOnTheDisplayedParent`,
  `.AnnotateMergeRanksParentsByTheirOwnSegment`, the typed `MergeIsSpelledThroughTheParentWithTheMostLabels` and
  `MergeKeepsTheArrivalOrderWhereTheBudgetChargesTheChain` × 3; M1, its regime test and M2 disabled, run on the
  final walker: what changed after them are comments and the positional oracle's depth pin); the annotation filter 574 of 574 (one more: the annotate pin matches it); integration *Traverse* 93
  of 93 (with `test_sigterm_while_the_index_loads`, which times out on the build before the fix) and *api* 171 (16
  expected failures, as before) and 21; the CLI fixtures and the hand-made ones up to date (R21 (4)'s exception
  changes none: no budget fixture has a merge); GCC 13.5 `-Werror` (`-Wno-error=stringop-overflow`, as staging
  builds) on the changed `src/` and `tests/` files, the one warning `walker.cpp`'s `-Wstringop-overflow` that was
  there before, and on the 112 `src/` and 75 `tests/` files with the headers as they are now.

## 26.9 Open (for the owner)

- **X-C8, returned** (§26.5): the D4 gate holds for header labels at refseq's chain counts and for column labels
  below 16 chains a run, and does not hold for column labels with 16 or more chains a run (refseq33m's taxid
  columns) nor for header-heavy seeds: there the depth at a stop drops 13–34% and, where many labels walk together,
  the seed needs 2–4 times the memory to complete and fails at depth 0 at the budget that completes it without
  coordinates. The library turns coordinates on under a memory budget (W2's library part:
  `AUTO_COORDINATES_UNDER_MEMORY_BUDGET`, which cites the gate). One option: auto-on under a memory budget only for
  header labels with few seed labels, off for column labels (the library's constant and its rule; the request can
  always ask for them).
- **R21 (4) under a memory budget** (§26.6): built as the review's option (a) — `tree` and `full` detail keep the
  arrival order under a memory budget, so every budget's stop is where level 5 put it and the rule changes the
  display only. Its price: on a server with a memory maximum (ops note O2, recommended for staging) every `tree`
  and `full` response is budgeted and keeps the arrival order; only `graphlet` and `summary` follow the rule. The
  alternative (b), the rule everywhere, moves stops shallower in those details (52 of 6,552 cached budgeted
  responses; 10 of 294 in the review's `tree` sample at 8 MiB; 1–4 levels on an arm, twice 30) and is one
  condition (`majority_first` in `Walker::merge_level`); (c), charging the longest parent's chain at every merge,
  moves more.
- **The annotate rule** (§26.6) compares the parents' own segments, so after nested merges it can display a walk
  carried by fewer labels; a rule over equal stretches of the displayed walks ranks those right but changes the
  library's check (`derive.carried_labels`) with it.
- The owner's pending decision of review 6: whether `setup` pieces count toward `observed_max_uninterruptible_ms`.
- M4 on staging; then C2 (an outcome class for `coordinates`) at the first external review of level 6.

# 27. The review of 2026-10-06: its P2 items, folded into level 6 *(v5.9, 2026-10-06)*

A three-day review of `278a53dd..67bef367` (report `CODE_REVIEW_metagraph_2026-10-06.md`, 189 items, no P1) found
ten P2 items. Level 6 was deployed nowhere, so the owner folded their fixes into it: `feature_level` stays 6, and
SPEC §10.3 states every text and output that changed (T58). The server's half:

- **X1, the D3 contract.** D3 went into the server at `7aaee760` with no contract while the level said 5 (§25's
  note). Level 6 had stated most of it; the sentences still false are corrected (SPEC §6.1 step 4, §7.0's
  `label_evidence` and `derivation` rows, T30c, a T-row for the D3 tests, §14 above, §25, the level 4–5 review
  request). **The owner's decision:** the `derivation` limitation's `observed` for `time_budget` keeps two units —
  the elapsed ms on a failed seed, j, the k-mers read, on a walked one — told apart by the result's shape and
  stated with both units; no wire change within level 6 (a field such as `kmers_read`, or a cause of its own,
  would be one).
- **X2, the release rule.** The bound is enforced cooperatively, and only the delivery checks compare it (every
  4096 objects of a result or of the MGT text, every 64 KiB of text, between compression blocks, before the
  transport); the walk's polls compare the walk-until, which lies below it, and only one poll in `poll_stride`
  (8) reads the clock, besides the forced ones (between seeds, before a paced read's chunk, in the lookahead). So
  an attempt runs past its bound until its next delivery check: the rest of the piece the bound fell into (a
  read's chunk, a head, the mapping of a seed's k-mers, the seed phase, a finalisation, a gap between delivery
  checks) and, when its walk had not stopped, the walk up to its next clock-reading poll (up to 7 heads), the
  stopped seed's finalisation and the building of its result up to the first check. No piece has a stated bound
  (`max_uninterruptible_ms` null, decision 3c-N5) — on a 6.94 Mbp seed GET answered `running` 2.7 s past the
  clock release's instant. The texts said "covers it" and "as stopped". Now `attempts.bound`, `not_after` and
  `release_rule`, `deadline_check.rule`, SPEC §5, §6.8 and §10.3 say that the attempt is past its bound apart
  from that run, whose length has no stated bound; that an answer of `running` or `stopping` past the instant
  shows it; and the Assumed list names that run and the other server that serves the same ledger. Text only.
  The library holds the clock release while the last answer said `running` or `stopping`, by default (LRG-R4).
  *(The first version of these texts said "one uninterruptible step"; the review of the fixes showed several
  pieces run past the bound, and that the walk-until is seen only by a poll that reads the clock.)*
- **The search service's release parity (LRG-R1, R2), text only.** The 409 refusing a copy carries the id's state
  as GET answers it, so its `finished` is a finished state (with the same replay conditions). An `expired` 409 is a
  release ground: the id's registration is checked before the expiry, under one lock, so no copy was running or
  registered on that `server_instance` and none can start there later; it settles nothing (an earlier copy may
  have run and left retention), and it assumes the server's clock does not step back below `not_after_ms`. With
  them, on the same strings: the refusal order (C24), the hold that refused copies extend (C30), the
  `expect_server_instance` pattern and its 400 (X4), the content-timeout cap (C16, also in the `attempt_deadline`
  statements and the 503 message) and "the lowest walk-until seen" (C20).
- **W1.** A derivation's windows are 64 k-mers whatever `annotation.batch_kmers` is. Each window is one work charge
  before its k-mers are consumed, so the later windows' width, `clamp(batch_kmers, 1, 64)` on a format whose reads
  are not budget-aware, placed the work comparisons: a work-budgeted derived seed walked partial at `batch_kmers`
  1 and failed at 64. An output change for `batch_kmers` below 64 with a derived seed, stated with level 6 —
  also, on a format whose reads are not budget-aware, the soft memory observation under a memory budget, since a
  later window holds 64 rows where it held `batch_kmers` (`soft_excess_bytes`, `held_bound_bytes`, the per-seed
  soft excess, `memory_bound_soft`'s `observed`; the review of the fixes found it unstated).
- **W2.** The derivation's coordinate bytes (`coords_at`, under trace) are a running total instead of a sum over
  every consumed k-mer at every observation, which made a trace derivation under any memory budget quadratic
  (51 s against 0.16 s for 100 kbp). The total equals the sum exactly — it decides `beside` and
  `memory_bound_soft` — which a debug build asserts and a test hook checks in Release, on a row-diff annotation
  through `beside` too (the one place where the total decides a stop). The timing test is skipped in a debug
  build, whose assert recounts the state at every observation.
- **W3.** The lookahead's chains ran between two checkpoints and read only the seed's own deadline; at
  `batch_kmers` 60,000 a cancel was seen 1.5–11.8 s late and a walk-until passed by 0.6–6 s. They now read the
  attempt's stop every 16 graph steps and before each chain's key mapping (`kLookaheadPollSteps`), a stop ends
  the lookahead, and each poll ends a head piece. The attempt's stops are read through the poll that reads the
  clock (`poll_now`), paced or not: it records a passed walk-until with the attempt, so the next checkpoint
  stops the walk (the first version, unpaced, stopped the lookahead on `ms_left()` without that poll, and the walk
  ran up to 7 more heads to the next clock-reading poll). `observed_max_uninterruptible_ms` counts head pieces as well as
  reads (and, in `deadline_check`, the delivery gaps), and `deadline_check.rule` names what stays unpolled: a
  chain's key mapping, a read's preparation before its first chunk (its cache lookups and the ordering of its
  keys: 90–260 ms for the 1 to 3.8 million keys of a lookahead in the reviewer's fixture, now the longest
  head piece there), the lookahead's clearing, the seed phase's own processing (D12). Whether setup pieces
  should count toward the observation stays the owner's open decision of review 6 (§26.9); they are not counted.

- **The review of the fixes (server half), folded into level 6 as well.** The overrun texts above (X2), the
  walk-until seen only by a poll that reads the clock (C20: SPEC §6.8 and §7.3, `attempts.bound`,
  `delivery_reserve.rule`, pinned with the server's stride 8 by a unit test of the reviewer's U11-02 driver), the
  unpaced lookahead's stop (W3 above), the W1 soft observation (above), the W2 timing test skipped in a debug
  build, and the duplicate's 409, which said "within the retention period: an attempt runs once" — false for an
  id held past its retention, and the D3 overclaim — and now says the id "is running, or retained or held after
  it finished, on this server (…): it is not run while so"; the cancel's 429 `tombstones_full`, which embeds the
  retention text, is named among the changed texts. An unbraced `EXPECT_EQ` failed GCC's
  `-Werror=dangling-else`: the changed test files are compiled with GCC 13 now, not only `src/`.

The library's half (L1, L2, O1–O3, the reader gaps LRG-G1..G7 and the release verdict, WORK_MODEL 2) is in the
library's own changes; the two halves state the release rule alike.

# 28. The review of 2026-10-06: its P3 server items, folded into level 6 *(v5.10, 2026-10-06)*

The owner approved the P3 server items of the same review ("Server (affects staging directly): fix them."). Level 6
is the level of the build deployed to staging (`ea285c2e`), so these are **level-6 corrections**: `feature_level`
stays 6, SPEC §10.3 states the three outputs of requests that ask for nothing new that change (C10, W9, C26) and
what else changes, and a client keyed on `feature_level` cannot tell a corrected build from `ea285c2e` (T59).

- **C10, `labels_from_seed` (X-DUP-01).** Four writers of a failed seed's `seed` block set it four ways; the
  not-started one said `seeds[].labels` empty, true for every annotate seed beside a walked one's false.
  `labels_derived_from_seed(seed, mode)` (walker.hpp) is the one rule — constrain mode and no labels named — for
  every result without a finished walk to read it from: a budget failure (`SeedBudgetError`), a failed
  derivation, a seed never started. An output change in annotate mode only.
- **W9, the format-(ii) dictionary stop (U03-03).** On a format whose reads are not budget-aware an annotate
  level's new labels are charged after its read (`charge_dictionary`), and the check after it threw a head's
  trip: cause HEAD, the dictionary's need as exact, without the levers that name fewer labels. When that charge
  is what crossed the budget (the account within it before) the trip is LABEL_NAMES at the level, as a lower
  bound (`ResourceStop::names_after_read`; the message names the labels and their bytes, not a row's demand or
  bytes left, which do not exist there). `observed`, `used` and the stop's depth are unchanged; the message, the
  actions and the `walk_domain` effect change (in a graphlet, its `Q` and `K` records and so `graphlet_bytes`).
  Charging each call's names inside the format-(ii) read, as the budget-aware read does, would stop such a level
  earlier and change where walks stop: not done.
- **C8, the measured stop latency (U11-01).** `seed_walked` is called again for a seed whose result was being
  built when the client left or a writer refused; the second call measured the time since the walk-until, the
  building included, as stop latency, which became the server's `measured_stop_ms` and shut out later attempts
  whose bound was no longer. Measured on the first call per seed only. No cap was added (the verdict's option):
  a long stop measured at a walk's end is a real stop.
- **C9, delivery memory (U12-02, U13-03, X-EFFICIENCY-02).** The per-seed texts are moved into the assembly,
  which reserves the response's exact size (the envelope's members after `results` written first) and frees each
  text once copied; the transport's buffer (the Response's `asio::streambuf`, which grows by doubling) is sized
  once through `rdbuf()` after the last check, so nothing reaches the transport before every check passed. On
  mini_refseq the delivery peak per byte of text fell from 2.60 to 1.52 (gzip) and 3.51 to 2.57 (identity).
  Feeding deflate from the texts, and sending the body in pieces, were left: the first saves little beside the
  freed texts, the second would commit a 200's header before the last check.
- **C26, zlib's 32-bit input (U13-04).** `compress_string` hands the text over in pieces of at most 2^32 − 1
  bytes (`Z_NO_FLUSH`, the last with `Z_FINISH`); a text of one piece makes the same calls as before. An output
  change: a compressed response of 4 GiB of text or more is the whole text where it was a stream of its first
  (size modulo 2^32) bytes answered 200, on every route (`process_request` is the server's one compressor).
- **C27, the server's own shutdown (U13-02).** With bytes waiting, every TCP state past ESTABLISHED reads gone (the
  content timeout's shutdown left Linux sockets in FIN_WAIT2, read connected); SYN_RECV (TCP Fast Open, not
  enabled here) stays connected; EBADF and ENOTSOCK read gone, as the header always said.
- **W11, the budget-aware lookahead (U05-01).** It decoded every run of a warm larger than the cache and evicted
  its own earlier runs (`cache_budgeted` evicted wholesale); it now stops at the first run the cache cannot keep
  beside its own runs (by count before decoding, by bytes after), evicting what was cached before it for its
  first run only. The level's fetch keeps its wholesale eviction. Physical work only: on the reviewer's walks
  36,275 → 33,015 tuple rows (batch 2,048, 16 MiB), 38,441 → 30,667 (8,192, 16 MiB), 1,149,725 → 258,101 (4,096,
  8 MiB), every untimed byte the same. Not uniformly fewer decodes, and not the whole excess (as the verification
  of this round found): of 180 walks on mini_refseq under 16-64 MiB the fix decoded fewer tuple rows in 24,
  more in 10 (up to 1.4%; 4.4% with the path cache's stored rows), the same in 146. What remains is a later
  warm's first run, which still evicts wholesale and with it rows earlier warms read ahead that the walk had not
  reached: up to +20% tuple rows over the walk under a work budget alone at 16 MiB, +52% at 32 MiB, +32% at
  64 MiB (`batch_kmers` 2,048-8,192; SPEC §6.8 states it). Evicting the oldest entries first (a list through the
  costs map, a warm keeping a quarter of the cache) was built and measured on a grid of 160 walks: 12% fewer
  tuple rows (the review's walk +0.4% over the work-budget walk), but the label cache stayed full and the
  row-diff path cache, which `make_room` gives only what the label cache leaves of the allotment, read 2.5 times
  the stored rows (`annotation_fetch_ms` +44%); evicting the oldest half read 2-7% more rows (tuple and path
  rows) in walks at the default `batch_kmers` 64, 2.8 times in one at 2,048. Not taken; the comment at
  `cache_budgeted` says why.
  `MiniRefSeq.LookaheadKeepsItsRunsUnderAMemoryBudget` holds the review's walk within 25% of the work-budget
  walk (+19%, +20%, +11% at 2,048, 4,096, 8,192; `ea285c2e` +31%, +33%, +39%). The structural lookahead's own
  re-warming of the chains it clears (the reviewer's "related mechanism" at `batch_kmers` near
  `max_lookahead_`) is not changed: still 5.4 million rows at 8,192 under 8 MiB.
- **W16, `/resolve` explicit labels (X-EFFICIENCY-04).** One scatter pass over the k-mers replaces the per-label
  scan (trace: 2,000 labels 3,457 → 262 ms; the same profiles, checked against the discovery's); presence uses
  the same accumulator instead of a labels × k-mers bitmap; the client is checked every 4,096 k-mers of the pass.
  The per-k-mer hit copies of `fetch`, the larger memory term with many labels, are the planned label cap's.

# 29. The owner's decision #17 (2026-10-08): derived data outside the identity *(v5.11)*

The graph's dummy-edge mask (`<graph without .dbg>.edgemask`) and its Bloom filter (`<graph without .dbg>.bloom`)
are **derived data** of the graph: computed from it alone, they decide which counts are exact (`/pattern` counts
exactly with the mask and states bounds with an estimate without it, `SPEC-pattern-search.md` §18) and how fast
k-mers are looked up, never what an exact answer is. They leave the identity inventory (`index_load_inventory`,
`index_bundle_files`); `index_derived_files` lists them apart, and `traverse --index-inventory` prints them under
`derived` with the rule (`derived_rule`). A manifest listing one (`*.edgemask`, `*.bloom`) is refused, naming the
entry and the rule; `index_manifest.py` never writes them, refuses them with `--extra` and flags them with
`--verify`. Consequence: adding a mask to a deployed index (`transform --mask-dummy`) leaves its `index_fp` and
`index_meta_fp` unchanged (shown on a copy of `build/mini_refseq`: the stored graphlets still compare as one
index, and every `/traverse` and `/resolve` answer replayed byte for byte apart from timing). Two spellings of one
graph, one with a mask beside it, are one index with one identity.

This reverses the pass-5 finding-3 fix (§24) for these two files, and its risk comes back, stated
(`SPEC-labeled-traversal-core.md`, "The index identity"): no fingerprint covers derived data, so a stale or
foreign derived file beside the graph is not detected by `index_fp`. The loader checks a mask's size and its
W = `$` edges (`mask_invalid` on `/pattern`) and a Bloom filter's k and mode only; a Bloom filter of another graph
with the same k and mode can hide k-mers under an unchanged `index_fp` (the pass-5 probe). Derived data is written
only by `build --mask-dummy`, `transform --mask-dummy` or `transform --initialize-bloom` on the graph it sits beside.
A spot check at load (a fixed-seed sample of the graph's k-mers looked up through the Bloom filter), or the graph's
size and digest stored in the derived files' headers, would close it; the owner's call. Manifests written between
the pass-5 review and this decision that list a mask or a Bloom filter refuse to start and are written again
without them (their `index_fp` changes once).
