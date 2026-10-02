# Design: the traversal graphlet — retrieve once, process locally

**Status:** v5.2 (2026-10-02; v5.2 = the owner's conservative outcome rule in §14) — **approved for implementation** by the fifth external review (no further architecture
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
O <walks c|p|f> <branch_diagnostics c|x> <label_evidence c|l> <delivery i|s|p>
                                       # (v5) the per-seed guarantee dimensions (§14), one per document, before the arms
Q <scope> <resource> <phase> <requested> <effective> <used> <remaining> <actions csv> <message>
                                       # (v5) resource_stop (§14), at most one per document, after O; message is the
                                       # free-text last field (percent-encoded like names)
K <arm l|r|*> <kind> <knob> <limit VALUE> <observed VALUE> <complete_to_bp|*> <extra name=VALUE,...|.> <effect>
                                       # (v5) one stated limitation; arm * = seed-level (seed_labels, server_clamp,
                                       # derivation). (v5.1) VALUE is typed by a one-letter prefix: i:<integer>,
                                       # f:<float, codec rule>, s:<string, percent-encoded incl. space and comma>, or
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

`walk` in `H`: every `G` stores its bases in walking order on **both** arms, so outward index `i ∈ [from_bp, from_bp+length_bp)` is `bases[i − from_bp]`, and every position field (`from_bp`, `to_bp`, `at_bp`, `route_bp`, `complete_to_bp`) indexes bases directly. Natural orientation: right flank = concat root→leaf; left flank = `reverse(concat root→leaf)` = concat leaf→root of per-segment reversals = `spell_path` (walker.cpp:3220-3232). Whole molecule, natural: `natural(left) + seed + natural(right)`. Seed coordinate of outward `i`: `|seed| + i` (right), `−(i+1)` (left). The serializer reverses `Segment::sequence` of the left arm back to walking order (finalize reversed it at walker.cpp:2877-2878).

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
  graph, annotation, sidecars (`.seqs`, `.anchors`, `.rd_succ`, `.edgemask`, `.coords`) — with size and sha256, plus
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

# 7. Size and memory (SRA case; estimates from per-record byte costs, to be re-measured by `scripts/traversal/graphlet_measure.py`)

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
| `label_evidence` | `complete` · `lower_bound` · `qualified` | `complete` only when no label-class limitation applies; `lower_bound`: evidence may be *missing* or understated — cut lists (`label_lists`, `inexact_counts`), dropped carriers (`seed_labels`), a cut switch-source list (`switch_sources`), greedy losses (`greedy_losses`); `qualified`: something reported may be *overstated* — a column-label trace across a record boundary (`trace_record_boundaries`); `qualified` wins when both apply |
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

*(v5.1)* **Charging work and releasing capacity are separate operations.** Consumed work is charged when the response (or the lease expiry) settles the attempt; occupied capacity — the worker slot and the memory reservation — is released only once the backend has actually stopped: on its response, or after the lease, which is only a valid release point because the backend enforces the same deadline itself (the request's `time_budget_ms` plus the server's hard request timeout, so an attempt cannot outlive its lease). A cancellation the backend has not acknowledged releases nothing.

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
`resource_limit` end reason; the four outcome dimensions; the `limitations` kinds and their effect templates; the
`resource_stop` structure; `attempt_id`/`budget_id`/`locus_id` in requests and usage in every response. Enforced in
stages, each stated in every response until it lands (`limitations` entries): (1) **now** — the existing caps, the
evidence boundary, stated limitations, the outcome dimensions (being implemented in the JSON); (2) the two-phase
head commit and lazy paths in the walker; (3) the budget-aware row-diff decode path in the annotation library
(`memory_bound_soft` until then); (4) the ledger in the search service. The promise at every stage is the same:
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
