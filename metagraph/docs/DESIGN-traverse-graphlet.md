# Design: the traversal graphlet — retrieve once, process locally

**Status.** MGT v1 is frozen: every change stays within its records, fields and tokens, and a change to any of them
is MGT v2 (§18.4 collects the reasons that would justify one). The server is at feature level 6
(`capabilities.feature_level`; SPEC §10.3 states what each level adds). `SPEC-labeled-traversal-core.md` is the
normative contract: its §7.5 carries §2 of this document, and where the two differ the SPEC is right. This document
states what the system is and why it is built that way. Companion document: `DESIGN-labeled-traversal-endpoint.md`
(the endpoint's design note).

All paths are relative to `metagraph/`.

**Contents, by topic.**

| topic | sections |
|---|---|
| purpose, architecture, the facts the format rests on | Purpose, §0 |
| the graphlet: request, MGT v1, JSON summary, index identity, the walker's run fields | §1–§4 |
| the Python library: queries, derivations, claims, continuations, local limits, files, MCP tools, sizes | §5–§7 |
| where the code and the conformance tests are; compatibility with the other details | §8–§10 |
| the walker: merges and the displayed parent; seeds, derivation and seed labels | §11, §12 |
| `/resolve` under a deadline | §13 |
| budgets, the guarantee contract, annotation reads, deadlines, stops | §14 (§14.1–§14.5) |
| attempts, cancellation and the release rule | §15 |
| the index identity as served | §16 |
| the server: capabilities, delivery, shutdown, deployment | §17 |
| record coordinates | §18 |
| open questions | §19 |

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

The backend owns graph work, the label oracle, caps and budgets, determinism, the completeness certificate
(`complete_to_bp`, `completeness_scope`, `walk_rule`, `cap_trigger`), the lossless dump and the small diagnostics.
The local library owns parsing, every derivation, ranking and pagination, spelling and orientation, FASTA/GFA,
comparisons, sub-tries, handle lifecycle, size caps on what the agent sees, and building the *next* backend request
(a continuation is a new traversal — the library prepares it, the backend runs it). Nothing in the local layer reads
the graph.

Why a custom text and not more JSON views: the measured 3.72 MB compact JSON (0.18 MB gzipped) of a 2,500-leaf
exhaustive trie is ~95 % key names, repeated id lists and structure derivable from a far smaller core; a line-based
text that stores only the primary data is ~0.3–0.5 MB, parses in tens of milliseconds, is self-contained (seed
sequence, provenance, orientation rule in the file), can be `grep`ped, and is written by a streaming writer rather
than a `Json::Value` tree (whose deep copies cost 35–40 s on 50 MB responses). GFA cannot carry per-node label runs,
label ends with reasons or `route_bp`; it is a local export instead.

Design choices: the text is embedded as a JSON string under `detail: graphlet` (transport, routing, CLI and error
path stay the same; it costs one `\n` escape per line); runs are primary and label ends derived (which halves the
biggest block) at the price of three walker fields (§4); orientation is stored in walking order on both arms; a
derivable field is written as `*` only when the writer has verified that the rule reproduces the walker's value, so
every production dump tests the derivations.

# 0. The facts the format rests on

`SeedResult` is fully materialised, so the retrieval format is a serializer in `src/cli/traverse.cpp`, not a walker
change, apart from the run fields of §4.

1. **Transport.** `process_request(..., compact=true)` writes compact JSON and compresses it (gzip first). The text
   is embedded in the JSON (`output.detail: "graphlet"`) rather than sent as a new content type: `server_utils.cpp`,
   the routing, the CLI (`json.loads(stdout)` in the integration tests) and the error path need nothing new, and
   multi-seed responses come for free. The price: the field separator must be a space (a tab would JSON-escape to
   two bytes), and one `\n` escape per line (~15 KB on the SRA case).
2. **Runs are primary, label-end events are derived, and the anchor needs a walker field.** Run *ids* cannot be
   reconstructed by a sweep over the output: clones created at splits (`commit_entries`) are identical rows whose
   ids depend on frontier order, and `end_labels[].run`, `label_summary[].runs` and `prev_run` reference those ids.
   Reconstructing them means replaying the walker's frontier order under `lowest_loss_first`/`most_supported_first`.
   The terminal `loss` and `branches` of a lineage that ends inside a walk likewise exist only in the walker's
   discarded state, and `Split::ambiguous` cannot be derived (§2.5). Three fields in the walker (§4) are cheaper.
3. **A run can end with no `label_end` event.** `end_run(arm, src, at, LABEL_LOST)` ends a switch *source* silently
   ("the lineage continues only under other names"); only the `SWITCH` event on the target documents it, and
   `growth[].label_ends` counts it. So "every ended run ↔ one label_end event" is false in one direction; the `R`
   record carries a distinct end token (`Lw`) for it, and the anchor is set there too (§4). `#label_end <= #runs`
   holds.
4. **Two `route_bp` values.** `LabelEnd.route_bp` is the *entry's* route_bp (set at the latest merge the entry came
   through), while `LabelRun.route_bp` is the *earliest* stamp (`stamp_route`); they differ on a run that passed two
   merges via non-first parents, so `end_labels` cannot take it from the run and the leaf record carries it
   explicitly (`T`).

Further facts the format relies on: splits always have ≥ 2 children (one followed successor continues the
segment), so grouping children by parent reproduces `splits[]` exactly; `revisit` carries `length_bp` and
`same_distance`, `hairpin` carries "followed" — the graphlet keeps them (and the JSON details carry them as additive
fields).

# 1. Requesting it

- `strategy.output.detail`: `summary | tree | full | graphlet`. `output.sequences` keeps its meaning (`false` → `G`
  lines without bases, first base kept, continuation sequence explicit); `continuation_bp`, `max_branch_events`,
  `profile_bin_bp` apply as in the other details.
- The normalized echo includes `output.detail` and `output.timing`; both keys parse, so the echo stays resubmittable
  verbatim.
- `capabilities` carries `"graphlet_format": 1` and `"detail_levels": ["summary","tree","full","graphlet"]` (in
  every response and `/traverse/capabilities`).
- The CLI is unchanged: `metagraph traverse` prints the JSON with the embedded strings. `Graphlet.save()` writes the
  standalone `.mgt` file (§2, the `J` line).
- `output.label_lists` / `output.label_runs` are not request fields: the graphlet's delta coding carries what they
  would have.

# 2. The graphlet text format, MGT v1

## 2.1 Lexical rules

- One record per line, `\n`-terminated, fields separated by one space, first field a capital letter. The only
  free-text field is the **last** field of `L`, `X`, `Q` and `K` (names, messages, effects) and runs to end of line;
  in it `%`, LF and CR are percent-encoded (`%25 %0A %0D`), nothing else. Every other field is ASCII without spaces.
  Seed ids live in the JSON summary, not in the body.
- Integers decimal. **Floats, one normative algorithm:** take the *shortest round-trip significant digits*
  `d₁…dₙ` and decimal exponent of the double — the digits of C++ `std::to_chars(first, last, x,
  std::chars_format::scientific)` without a precision, equal to those of Python `repr(x)` (both produce the unique
  shortest digit string, nearest on ties) — and **expand them positionally**: zeros padded up to the decimal point,
  a `.` only when there is a fractional part, no trailing zeros, `-0` written `0`, `inf` for +∞ (NaN never occurs;
  all floats in the format are costs, budgets and losses ≥ 0). Python: `format(decimal.Decimal(repr(x)), 'f')` with
  `-0`/trailing-`.` normalisation. So `1e23` → `100000000000000000000000` and `1.2345678901234568e20` →
  `123456789012345680000`. Why this rule: `to_chars` and `repr` in their default forms disagree (on `0.0001`), and
  `to_chars(..., fixed)` writes the exact binary value (`99999999999999991611392` for `1e23`), which round-trips but
  is not the bytes Python writes. A differential probe of the rule — 20,000 non-negative finite doubles drawn over
  raw bits plus the edge cases above, `DBL_MAX` and subnormals — gave 0 mismatches between the two languages.
  Booleans `0|1`. `*` = absent / "derive by the rule of this field". `.` = empty list. Lists comma-separated; a list
  of lists uses `|` (never contains spaces).
- **RANGES**: ascending label ids with consecutive runs collapsed: `0-5,7,9-12`; `.` = empty.
- **Strings** are UTF-8. `L.prefix_len` counts **bytes** of the previous name's UTF-8 encoding, and the encoder
  shortens the shared prefix to the nearest UTF-8 character boundary, so the suffix is always valid UTF-8 on its
  own; a reader reconstructs `prev_bytes[:prefix_len] + suffix_bytes` and decodes. A name may contain spaces and may
  equal `*` or `.`: it is the remainder of the line after the fixed fields, preserved exactly (after
  percent-decoding `%25 %0A %0D`).
- **SETEXPR** (against a *base set* defined per field; the bases are normative: `G.entry` → parent[0]'s end set
  (merges: union of the parents' end sets; root: none — explicit RANGES only); `G.end` → this segment's entry set;
  `P` → the previous `P` of the segment, the first `P` → the entry set; `C.labels` → the leaf segment's end set):
  `RANGES` (explicit) | `!` (equal to the base) | `!REMOVED` | `!+ADDED` | `!REMOVED+ADDED` (base minus REMOVED plus
  ADDED, both RANGES). The encoder emits whichever of explicit/delta is shorter, tie → explicit. Canonical in both
  encoders, so `dump(parse(x)) == x` byte-exact.
- **Codes.** End reasons (`traversal_types.cpp`): `D` dead_end, `L` label_lost, `B` loss_budget, `R` branch, `U`
  edge_reuse, `V` edge_reuse_rc, `J` rejoined_seed, `T` trace_break, `X` max_extension_bp, `S` max_steps, `P`
  max_live_paths, `N` max_paths, `O` max_output_bp, `M` time_budget, `W` beam_pruned, `Y` resource_limit (memory or
  work budget, §14); a second letter keeps the walker's text qualifier (the quorum texts at the BRANCH ends): `Rm`
  minority, `Rb` below_min_labels, `Rs` split_limit, `Dh` hairpin, `Ls` superseded, `Lx` switch_sources. Enum and
  qualifier both survive (the JSON details replace the enum by the text). Arm `l|r`; status `c|t|p`; scope `p|u`
  (per_path | united_history); mode `c|a`; support `k|t`; reconverge `m|k`; label kind `c|h`; hairpin `f|s`.

## 2.2 Records

Document order: `H [J] S X* L* O Q? K*`, then per requested arm, left before right: `A B* V*`, then per segment in id
order `G P* E* T? C?`, then `R*`; finally `Z`.

```
H mgt 1 <k> <regime> <alphabet> <mode c|a> <support k|t> <reconverge m|k> <cap> <continuation_bp> <seed_index> walk <index_ns|*> <index_fp|*> <index_meta_fp>
                                       # index identity (§3.1): <index_ns> = the deployment's name for the index
                                       # ([A-Za-z0-9._-]+, server flag --index-name; * = unset), <index_fp> = the
                                       # digest of the index bundle's manifest (* = no manifest: joins are
                                       # unverifiable), <index_meta_fp> = the metadata hash, a negative check only.
J <compact JSON>                       # standalone FILES only: the response envelope with results[] reduced to this
                                       # seed's summary (graphlet string removed; §5.5). Written by the library, never
                                       # by the server. At most once, directly after H.
S <validated_seed_id> <length_bp> <num_kmers> <num_seed_labels> <SEQUENCE>
                                       # the seed as validated (upper case, request orientation)
X <reason> <a-b,c-d> <name>            # dropped seed label with its k-mer runs on the seed (SeedResult::dropped_labels)
L <c|h> <column> <seq_id|*> <prefix_len> <suffix>
                                       # label id = ordinal (constrain: seed labels then extra; annotate: first seen).
                                       # name = previous L name[:prefix_len] + suffix (front coding). column/seq_id =
                                       # LabelRef (the stable cross-retrieval join key). seq_id * for column labels.
A <l|r> <status c|t|p> <complete_to_bp> <scope p|u> <live_paths> <live_labels> <exact> <max_seen> <nodes_truncated>
  <counters name=value,... >   # named and extensible, in the walker's order: steps, successor_enumerations,
                               # output_bp, pair_evaluations, edge_reuse_probes, reminimisation_rounds,
                               # max_reminimisation_rounds, refusal_scans, switch_sources_cut, ... — an unknown name
                               # is preserved by the reader, so a newer writer never breaks an older reader
  <cap_trigger reason,at_bp,segment,live_paths,live_labels,exact,demand | *> <branch_events_total>
  <evidence_complete_to_bp|*>   # branch-event boundary (* = no event dropped)
  <n_segments> <n_runs> <n_leaves> <n_splits> <n_merges> <bases>          # the last six validate the body
B <from_bp> <max_live_paths> <distinct_live_labels> <live_pairs> <exact> <steps> <divergences> <ambiguous_branches>
  <splits> <reconvergences> <bubbles> <tips> <blocked_repeat> <label_ends code:n,... | .>     # one per growth bin
V <at_bp> <segment> <chars> <labels_per_successor csv> <ambiguous RANGES> <dropped RANGES> <refused char:cause:RANGES;... | .>
                                       # the capped branch events, as stored (the only record of a not-followed successor)
G <parents csv|*> <from_bp> <length_bp> <entry SETEXPR|*> <entry_total|*> <end SETEXPR|*> <partition RANGES|RANGES|...|*> <split 0|1|*> <first_base|*> <bases|*>
                                       # segment id = ordinal within the arm (the walker creates parents first);
                                       # bases in WALKING order; first_base only when bases are absent (sequences:false)
                                       # <split>: Split::ambiguous of the split at this segment's end (1|0), * = the
                                       # segment does not end in a split. PRIMARY: the walker counts predecessor LINEAGES
                                       # before quorum filtering and excludes followed hairpins, so it cannot be derived
                                       # from the children's label sets (disjoint child sets can be ambiguous after a
                                       # switch; overlapping ones can be a divergence).
                                       # annotate merges: <partition> is one '.' per parent (labels_via_parent empty)
P <from_bp> <to_bp> <total> <SETEXPR>  # annotate: LabelSetRun; base = previous P of this segment, for the first P the entry set
E <at_bp> s <from> <to> <cost>         # switch
E <at_bp> b <char> <reason code> <total> <RANGES>      # blocked successor
E <at_bp> h <char> <total> <RANGES> <f|s>              # hairpin followed | skipped
E <at_bp> v <segment> <delta | =>      # revisit: distance delta, "=" = same_distance (keep mode)
E <at_bp> t <char> <length_bp>   /   E <at_bp> u <length_bp> <alleles>     # tip / bubble: reserved, never emitted
T <path_reason code|*> <extras label:loss:branches:route_bp,... | .>
                                       # this segment is a LEAF; extras only for labels where any of the three is non-zero
C <n> <loss_used> <branches_used> <labels SETEXPR against the leaf's end set> [<sequence>]
                                       # continuation (leaves ended by X or a resource reason): the walker's label logic is
                                       # PRIMARY (it is subtle), the spelling is derived, both arms stated:
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
                                       # <branches> <loss>: the lineage's TERMINAL values when the run ended or was
                                       # closed (LabelRun::branches/loss, §4). Valid at to_bp only: a claim cut at an
                                       # earlier depth reports them as unknown, never as the value at the cut.
O <walks c|p|f> <branch_diagnostics c|x> <label_evidence c|l|q> <delivery i|s|p>   # q = qualified
                                       # the per-seed guarantee dimensions (§14), one per document, before the arms
Q <scope> <resource> <phase> <requested> <effective> <used> <remaining> <actions csv> <message>
                                       # resource_stop (§14), at most one per document, after O; message is the
                                       # free-text last field (percent-encoded like names)
K <arm l|r|*> <kind> <knob> <limit VALUE> <observed VALUE> <complete_to_bp|*> <extra name=VALUE,...|.> <effect>
                                       # one stated limitation; arm * = seed-level (seed_labels, server_clamp,
                                       # derivation, coordinates). VALUE is typed by a one-letter prefix: i:<integer>,
                                       # f:<float, codec rule>, s:<string: every UTF-8 byte outside 0x21..0x7E, '%' and ','
                                       # as uppercase %XX; '=' stays raw — as in the golden vectors>, or u (unlimited) —
                                       # e.g. the scope limitation's limit is s:merge; <extra> holds further fields
                                       # (server_limit, demand, lists_cut). <effect> is STORED as the free-text last field
                                       # (percent-encoded like names): the C++ builds it from values, and templates would
                                       # have to substitute identically in two languages
Z <line count of the document including this line>
```

## 2.3 The `*` convention (lossless by construction, rules exercised on every dump)

A field marked `|*` has a *rule*. The C++ writer computes the rule's value from the walker's own data and emits `*`
**only if it equals the walker's value**, else the explicit value; the reader applies the rule on `*`. So the format
never depends on a derivation being right, and every production dump is a test of it. Rules:

- `G.parents *` = root. `G.entry *` (constrain only): root → ids `0..num_seed_labels-1` (the root's state is the seed
  labels); merge → union of the partition. In annotate mode `G.entry` is never `*`; its SETEXPR base is parent[0]'s
  end set (root: explicit RANGES). Split children in constrain mode: SETEXPR against parent[0]'s end set (which labels
  followed which branch leaves no event, so this is primary).
- `G.entry_total *` = |entry| (the only non-`*` cases are the annotate root's cut boundary list and cut children).
- `G.end *` (chronological): constrain → start from the entry set and apply, in increasing position (ties: ends
  before switch-ins), every run end anchored on this segment with `to_bp < from_bp+length_bp` (remove its label) and
  every switch event on it (add its `to`); the result is the end set. A set formula (entry ∪ switch targets − ended
  labels) would remove a label that re-entered after ending (A→B→A inside one segment). Runs ending exactly at the
  segment end are still in `labels_end` (SPEC §7.1: label sets change only through events). Annotate → last `P` set,
  or entry when the segment has no `P`. Any field may legitimately be written explicitly: the conformance test (T37)
  checks that every `*` the writer emitted reproduces the walker's value — not that `*` is always emitted.
- `G.partition *` = no merge (|parents| ≤ 1), or annotate mode, where `*` decodes to **one empty list per parent**
  (`labels_via_parent` is empty in annotate mode; stated so the reader reconstructs the right arity).
- **Canonical form.** The writer MUST emit `*` whenever the field's rule reproduces the walker's value, and the
  explicit value otherwise; every SETEXPR in its shortest form (tie → explicit); RANGES collapsed. `dump()` always
  writes this canonical form, so `dump(parse(x)) == x` holds for every canonical document — all writer output — while
  a valid non-canonical document (e.g. an explicit set where `*` would do) parses and dumps to its canonical form.
  `is_canonical(x)` is part of the library and of T37.
- `T.path_reason *` = none (semantic end). `R.from_label:cost *` = entered by seed. `R.prev_run *` = none.
  `R.structural_successors *` for `Lw`/`m`.

## 2.4 Orientation (the one rule)

`walk` in `H`: every `G` stores its bases in walking order on **both** arms, so outward index
`i ∈ [from_bp, from_bp+length_bp)` is `bases[i − from_bp]`, and every position field (`from_bp`, `to_bp`, `at_bp`,
`route_bp`, `complete_to_bp`) indexes bases directly. Natural orientation: right flank = concat root→leaf; left
flank = `reverse(concat root→leaf)` = concat leaf→root of per-segment reversals = `spell_path`
(`tests/graph/traversal/walker_paths_for_tests.hpp`). Whole molecule, natural: `natural(left) + seed +
natural(right)`. Seed coordinate of outward `i`: `|seed| + i` (right), `−(i+1)` (left). The serializer reverses
`Segment::sequence` of the left arm back to walking order (the walker's finalisation reverses it to natural order).

## 2.5 Stored vs derived

Stored (primary): bases; DAG (`G`); entry sets where not implied; merge partitions (`labels_via_parent`); annotate
presence runs (`P`); runs with anchor, end token, route_bp, from/cost, prev_run, structural_successors,
needed_budget (`R`); leaf path_reason and per-label `(loss, branches, route_bp)` (`T`); continuation
labels/loss/branches/n (`C`); switch/blocked/hairpin/revisit events (`E`); growth bins (`B`); capped branch events
incl. refusals (`V`); seed sequence (`S`); dropped labels (`X`); dictionary with column/seq_id (`L`).

Derived by the library (and by `to_json()` to reproduce `results[i]` of `detail: full`): `children`; leaves
(= segments with `T`); `paths[]` (leaf ordinal in segment-id order, chain via `parents[0]`, `length_bp = from_bp +
length_bp` — exactly the walker's finalisation, so path ids match the server's); `end_labels` (emitted and derived
in ascending label id) = `R` anchored at the leaf with `to_bp == length_bp` (every alive label at the leaf has an
ended run; every run ending at the leaf's last node is a label alive there) + `T` extras; `end_reasons`,
`n_labels`; `label_end` events from `R` (code ≠ `Lw`, `m`; at = to_bp); `reconverge` events from `G` with > 1
parent (at = from_bp, segments = parents); `splits[]` (children with one parent grouped by parent, ordered by
(at_bp, first child id); `kind` ambiguous ⇔ the parent's `G.split` is 1 — stored, see §2.2: "a label in ≥ 2
children's entry sets" is wrong in both directions; `labels_before` = |parent end| (constrain) / parent's last `P`
total or entry_total (annotate); branches: char = first base, `labels_distinct` = child entry_total, `labels` =
entry cut to `cap`); `labels_at_end` = `G.end`; `label_sets` = `P`; `runs[]` = `R`; `label_summary` (constrain:
`Walker::summarize()`, incl. "a merge does not clamp direct_bp"; annotate: the parents-first union/intersection
pass of `Walker::summarize_annotate()`; both normative in §5.1); `needed_budgets` (B-runs; a histogram, order not
information); `branch_events_truncated`; continuation sequence; `seed.labels`; every `*_truncated` flag
(`total > |list|`).

Not information, normalised by the conformance test: event order among equal `at_bp` (the walker sorts stably by
`at_bp` only), `needed_budgets` order, `growth[].label_ends` `null` vs `{}`.

## 2.6 Worked example (the `quorum` fixture)

`api/python/tests/data/traverse/documents/quorum.mgt`, written by the CLI: a right arm of 100 bp under
`min_successor_labels: 2` and `max_label_branches: 1`, four column labels, a three-way fork where the minority
successor is refused (counters shortened here):

```
H mgt 1 15 basic $ACGT c k m 64 1000 0 walk fixtures * 178cc93c330374f0
S 2e433ba3262fdd15 30 16 4 CGCAGGGGCGTGGGTCAGGCCAAAATCGGT
L c 0 * 0 x.fa
L c 1 * 0 y.fa
L c 2 * 0 z.fa
L c 3 * 0 w.fa
O c c c i
K r scope branching.on_reconverge s:merge i:0 * . completeness holds for the united-history rule, not per path; use "keep" for the per-path guarantee
A r c 100 u 0 0 1 3 0 steps=68,successor_enumerations=69,output_bp=68,…,switch_sources_cut=0 * 3 * 3 5 2 1 0 68
B 0 2 4 5 1 68 0 1 1 0 0 0 0 D:2,L:2,R:1
V 10 0 AG 1,4 0 . A:minority:0
V 12 0 CG 2,3 1 . .
V 40 2 G 1 . 1 G:minority:1
G * 0 12 * * * * 1 * AACTTCTGAAGT
G 0 12 28 1,3 * * * * * CCAAGACAATGGGCCGAGCAAATCCTCT
T * 1:0:1:0
G 0 12 28 !3 * * * * * GCTCGTTCCAGAGAACGAAACCCTACCT
T * 0:0:1:0,1:0:1:0
R 2 0 0 40 L 0 * * 1 1 0
R 1 1 0 40 D 0 * * 0 1 0
R 2 2 0 40 L 0 * * 1 0 0
R 1 3 0 40 D 0 * * 0 0 0
R 2 1 0 40 Rm 0 * * 1 1 0
Z 24
```

Reading: the index has a name (`fixtures`) and no manifest (`index_fp *`). The root's entry `*` is the four seed
labels {0,1,2,3}; at 12 it ends in an ambiguous split (`G.split` 1, stored): child 1's entry is `1,3`, child 2's
`!3` = the root's end set without 3 = {0,1,2}, so label 1 enters both children and has two runs (1 and 4, the
second ended `Rm`, minority, at 40 with one structural successor). The `V` at 10 keeps the refused successor A of
label 0 (`minority`), which no other record shows; it is why run 0 (label 0) ends with `branches` 1 while run 2
(label 2), ending at the same position for the same reason, has 0 — the terminal values of §4, which the JSON
details alone could not tell apart. `T` stores leaf extras only where a value is non-zero. `K r scope` states the
united-history rule of `on_reconverge: merge`; no merge happened (observed `i:0`), so `walks` is complete
(`O c c c i`).

## 2.7 Shared golden vectors — the freeze gate

Both implementations pass the same files byte-exactly, and these files are what holds MGT v1 fixed: a change to
either implementation that alters one of them is a format change.

- `api/python/tests/data/traverse/codec_vectors.tsv` (a gtest and a Python test read it): floats given as **raw
  IEEE bits** (hex) with the expected text — 0, -0, 1, 0.5, 0.1, 0.0001, 1e-05, 2.5e-07, 1/3, inf, 1e16, 1e22,
  1e23, 1.2345678901234568e20, 2⁵³, 2⁵³+2, every 10ᵏ for k = −20…25 and its two neighbouring doubles, `DBL_MAX`, the
  smallest normal and the smallest subnormal — plus a **differential section** of 10,000 non-negative finite doubles
  from a seeded generator over raw bits, expected strings produced once by the Python reference and committed;
  RANGES and SETEXPR against each normative base (shorter-wins, tie → explicit); percent escapes; front-coded names
  including multi-byte UTF-8 (`éfoo`→`ébar`), names with spaces, and the names `*`, `.`, `c:0`.
- **Whole-document fixtures** `api/python/tests/data/traverse/documents/*.{request.json,full.json,mgt}`: the C++
  writer produces each `.mgt` byte-exactly from its request (gtest + CLI), the library parses it, `dump()`s it back
  byte-exactly and reproduces `full.json` with `to_json()`. They include the counterexamples of §9 and the
  **clipped-merge evidence case** of §5.1 with its expected claims.

# 3. The JSON summary (per seed, `detail: graphlet`)

The envelope is that of the other details (`release`, `capabilities`, `strategy` with `output.detail/timing` and
`clamped`, `walk_rule`, `algorithm_version` `traverse-0.2`, `results[]`, `timing?`, and `usage` for an attempt,
§15). Per seed:

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
             "label_evidence": "complete|lower_bound|qualified", "delivery": "inline|spooled|paged"},   # §14
 "resource_stop"?: {…},  "limitations": [ {…} ],                     # §14; per arm too: "evidence", "limitations"
 "coordinates"?: {…} | null, "coordinates_reason"?: "…",             # §18, only with output.coordinates
 "graphlet": "<MGT text>", "graphlet_bytes": N, "graphlet_lines": N}
```

Derivation failures keep the `{seed, error}` shape of the other details, with no `graphlet`. No names, no per-leaf
rows, no `label_summary`, no `growth` (in `B`), no `label_dict` (in `L`). Size ≈ 2.5 KB envelope + ~0.7 KB per arm;
a 64-seed batch stays ~100 KB of summary.

## 3.1 Index identity

`release`, `k`, regime and alphabet do not identify an index: two different annotations over one graph report the
same four values, and their `c:0` would be joined. A hash of *metadata* (counts, ordered column names) is not
identity either: two indexes with identical counts and names that swap which sequences columns A and B annotate get
the same hash, and a per-ref name check agrees too. So:

- **`index_fp` is the digest of the immutable index bundle**, not of its metadata. The index build writes a
  **manifest** next to the files (`<prefix>.manifest.json`: every file of the identity inventory — the graph, the
  annotation and the sidecars the pair's loaders open (§16.2); never the graph's derived `.edgemask` and `.bloom`
  (§16.3) — with size and sha256, plus the builder's version and inputs), and `index_fp` = sha256 over the canonical
  manifest. Hashing hundreds of GB at server start is not an option; hashing at build time is free. For public
  indexes whose files were not built with a manifest, the deployment supplies one (computed once by
  `scripts/traversal/index_manifest.py`, or from immutable object-store digests such as S3 checksums) and the server
  loads it with `--index-manifest`.
- `index_ns` (`--index-name`) names the index for humans and routing; it is not identity.
- Without a manifest the server reports `index_fp: *`, and everything that joins labels across retrievals reports
  *unverifiable* — never equal. The metadata hash is kept as `index_meta_fp`: a cheap **negative** check (a mismatch
  proves different indexes), never a positive one.
- Tested: two indexes over the same records with identical counts and column names but swapped memberships get
  different `index_fp` and compare `unverifiable`/`different`, never equal; `scripts/traversal/build_mini_refseq.sh`
  writes the manifest so the fixture exercises the positive path.

How a server checks a manifest against what it loads, per graph, is §16.

# 4. The walker's run fields

`src/graph/traversal/walker.hpp` `LabelRun` carries three fields for the format, set at the sites where a run ends
or is closed, with no behavioural effect:

1. `uint32_t segment` — the segment on which the run ended or was closed (the anchor). `end_run()` takes it;
   `end_label()` passes the item's segment, and so does the silent switch-source end; `merge_level()`, in the
   `labels_via_parent` loop, sets it for an entry whose run is not the kept one (that run is closed by this merge).
2. `uint32_t branches` and `double loss` — set at the same sites from the `Entry` that ends or is closed. Branch
   counts increment *before* quorum filtering, so a label ending inside a walk can carry `branches = 1` or `0`
   depending on a successor that left no trace in the output (two annotations whose full JSON is identical can
   differ here; the `quorum` fixture of §2.6 shows it); `T` records these values only for leaf labels. Without the
   fields, `Claim.branches`/`loss` for interior ends would have to be reported as unavailable.

`Split::ambiguous` is already in `SeedResult`; it is serialized as `G.split` (§2.2).

Why the anchor is stored: the run↔segment association exists only in `Item`/`Entry` state and is dropped from
`SeedResult`. From the output alone it is ambiguous whenever two segments end the same label at the same depth
(clones at splits are identical rows), it is absent for the silent ends of §0.3 (no event at all), and for ≥
3-parent merges. Without it, `R` would need the events stored beside it (~+100 KB, the biggest block twice), or a
sweep that reconstructs run ids — not reproducible under non-breadth-first frontier orders. Also exposed as
`runs[].segment` in the JSON (additive) and asserted in `tests/graph/traversal/test_trie.cpp`: every ended run's
anchor holds its `label_end` event at `to_bp`; for a silent `Lw` end, at least one run has `prev_run` = this run,
starts at `to_bp` and was entered by switch (`from_label` = this run's label) — the switch events themselves may sit
on the *children* of the anchor when the source ends at a split; and every `m` run's anchor is a parent of the
segment created at `to_bp`.

# 5. The local library: `api/python/metagraph/traverse/` (stdlib only; pandas lazy)

`setup.py` uses `find_packages(include=['metagraph', 'metagraph.*'])` (a plain `include=['metagraph']` would exclude
the subpackage from a non-editable install). `metagraph/__init__.py` does not import `client`, so
`import metagraph.traverse` is pandas-free.

| module | contents |
|---|---|
| `_codec.py` | `encode_ranges/decode_ranges`, `encode_setexpr(ids, base)/decode_setexpr(tok, base)` (shorter-wins, tie explicit), `fmt_num/parse_num`, `pct_escape/unescape`, `REASON/QUAL` code tables, `LabelSetInterner` (content-hashed `array('I')`), `GraphletFormatError(line_no, msg)` |
| `model.py` | slotted dataclasses: `Label(id, kind, column, seq_id, name)`, `SeedInfo(…)`, `Segment(id, parents, from_bp, length_bp, walk, entry, entry_total, end, partition, split, first_base, presence: list[PresenceRun], events: list[Event], leaf: Leaf|None; derived: children, depth)`, `Event`, `PresenceRun(from_bp, to_bp, labels, total)`, `Leaf(path_reason, extras: dict[label,(loss,branches,route_bp)], continuation: Continuation|None)`, `Run(id, segment, label, from_bp, to_bp, end_code, reason, qualifier, silent, merged, route_bp, from_label, cost, prev_run, structural_successors, branches, loss, needed_budget)`, `Arm(side, status, complete_to_bp, scope, frontier, labels_per_node, cap_trigger, counters, growth, branch_events, segments, runs)`, `Graphlet(format, k, regime, alphabet, index_ns, index_fp, mode, support, reconverge, cap, continuation_bp, seed_index, seed: SeedInfo, labels, dropped, arms, summary: dict|None, derived_from)` |
| `parser.py` | `parse(text) -> Graphlet` (one pass; `A` counts and `Z` validated → raises on truncation), `dump(g) -> str` (canonical; `== text`), `Graphlet.from_response(result, response)`, `Graphlet.load(path)/save(path)` (the `J` line), `standalone_text` (§5.5) |
| `derive.py` | cached per arm: `leaves, paths, splits, label_end_events, end_labels(leaf), labels_at_end, label_summary, needed_budgets, continuation_sequence(leaf), reconverge_events`; `check_rules(g)` = the §2.3 rules recomputed and compared with explicit fields (the oracle that reports any `*`-vs-explicit divergence) |
| `ops.py` | the queries below |
| `export.py` | `to_json(g) -> dict` (`results[i]` of `detail: full`), `to_fasta`, `to_gfa` |
| `client.py` | `TraverseClient` (§5.5) |
| `attempts.py` | the answers about an attempt and the release rule (`release_verdict`, §15.5) |
| `coords.py` | record coordinates (§18.3) |
| `budget.py` | local limits (§5.4) |
| `store.py` | `GraphletStore` (handles, LRU, disk spool) for the MCP layer |
| `mcp_tools.py` | framework-agnostic tool functions returning dicts ≤ `max_bytes` (the MCP server, wherever it lives, registers them; §6) |
| `frames.py` | `frames(g)` → dict of DataFrames, lazy `import pandas` |

Public API (one line each):

```python
parse(text) -> Graphlet                       # BODY ONLY: no seed_id, derivation metadata, annotation counters or
                                              # replay strategy; to_json(), next_request(), deepen() and summary() raise
                                              # MissingEnvelope on such a graphlet instead of inventing defaults
                                              # raises GraphletFormatError(line_no, msg); validates A counts and Z
Graphlet.dump() -> str                        # canonical text; equals the input for an unmodified graphlet
Graphlet.from_response(result: dict, response: dict) -> Graphlet   # results[i] + envelope → summary attached
Graphlet.load(path) / .save(path)             # the .mgt file: H, J (envelope with this seed only), body
Graphlet.to_json() -> dict                    # results[i] of detail: full (natural orientation); the conformance oracle
Graphlet.summary() -> dict                    # agent-facing ≤ 2 KB: per-arm status/complete_to/counts, top labels, caveats
Graphlet.label(selector) -> Label             # TAGGED selectors: {'ref': 'h:<column>:<seq_id>' | 'c:<column>'} |
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
    # the SPEC §6.9 unit: one per RUN END (claims, not leaves). Claim(label, arm, path_id, segment, from_bp, to_bp,
    # route_bp, evidence_from (§5.1), kind∈{end, alive, diverged, merged, stretch, route_only, boundary}, reason,
    # qualifier, end_class∈{lost, dead_end, radius, capped, blocked, pruned, open}, entered_by, from_label, cost, loss,
    # branches, support, exact); strict=True raises IncompleteRecording when annotate nodes_truncated > 0 (a cut list
    # makes the oracle a lower bound)
Graphlet.label_walks(label, arm=None) -> list[LabelWalk]        # per run: (from, to, evidence_from, sequence, end, leaves_below, merged_into)
Graphlet.routes(label, arm) -> list[list[int]]                   # per-label routes through merges via G partitions
Graphlet.support_profile(arm, leaf, kind='displayed'|'route') -> list[SupportRun(from_bp, to_bp, labels, exact)]
Graphlet.support_changes(arm, leaf) -> list[Change(at_bp, added, removed, reasons)]   # where support changes and why
Graphlet.label_summary() -> dict                                 # direct_bp / reach_bp / reentries / runs, both modes
Graphlet.continuation(arm, leaf) -> Continuation(sequence, labels, loss_used, branches_used, seed_coord, as_seed())
Graphlet.next_request(arm, leaves, bp=None, reduce_budget=True, **overrides) -> dict   # resubmittable /traverse request (§5.3)
Graphlet.subgraph(selectors, arm=None, mode='any'|'all') -> GraphletView   # a VIEW over the backing graphlet:
                                              # every id is the backing graphlet's original id (no remapping, no new
                                              # record type); save() writes the UNCHANGED backing body plus, in J,
                                              # {view: {selectors, mode, arm, of: <entry handle or body digest>}}; load()
                                              # restores the view; its completeness is the backing complete_to_bp
                                              # qualified 'for the selected labels' (stated in the view's summary)
Graphlet.to_fasta(arm=None, leaves=None, with_seed=True, orientation='natural', width=None) -> str
Graphlet.to_gfa(with_seed=True) -> str          # S per segment on the seed strand, L with k-1 overlap, P per walk, LB/ER tags
Graphlet.compare(other, *, arm=None, labels=None, mode='claims'|'walks'|'labels'|'prefix_subset') -> Comparison
    # keyed by LabelRef (kind, column, seq_id) — names only for display; restricted to min(complete_to_bp); claims
    # cut at D get end_class 'open'; comparable only for the same index identity (§3.1: equal index_fp, and equal
    # index_ns and release when set; names agreeing per ref) and the same ORIENTED SEED SEQUENCE — not
    # validated_seed_id, which includes the permitted label names and so differs between a constrain and an
    # annotate retrieval of the same seed; results carry the support kind, both strategies and both completeness
    # scopes, and a per_path vs united_history pair is 'qualified', never 'equal' — cutting both to
    # min(complete_to_bp) does not make their termination comparable;
    # Comparison(comparable, reason, depth_used, equal, only_in_a, only_in_b, notes) — notes carry the merge qualification
Graphlet.memory_bytes() -> int
TraverseClient(host, port, api_path=None, *, session=None, timeout=900, release=None)
    .capabilities(graph=None, graph_path=None) -> dict;  .server_capabilities() -> dict
    .resolve(sequence, *, labels=None, discover=None, select=None, **opts) -> dict
    .traverse(seeds, strategy=None, *, detail='graphlet', timing=True, graph=None) -> TraverseResponse(envelope, graphlets, errors)
    .traverse_raw(request) -> dict;  .deepen(graphlet, arm, leaves, bp=None, **overrides) -> TraverseResponse
    # requests only; explicit Accept-Encoding: gzip, deflate; non-2xx → TraverseError(status, message); 503 →
    # ServerInitializing; an expired not_after_ms → AttemptExpired
GraphletStore(spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800, ttl_disk_s=7*86400)
    # ENTRY identity and BODY identity are separate: an entry handle is opaque ('g_' + 12 random hex) and owns
    # (normalized request, index identity, envelope, body digest); bodies are deduplicated by sha256 in the spool.
    # Two requests with different loss budgets can produce byte-identical bodies (no switch happened) but continue
    # differently, so a body hash must never be the handle. Spool bodies are verified against their digest before
    # reuse.
    .put(response, request) -> handle; .get(handle) -> Entry (lazy parse, LRU); .free/.list/.save/.load/.sweep
```

Evidence semantics in the round trip: `route_bp` lives on the run (`R`, the *earliest* stamp) and per leaf label
(`T`, the entry's *latest* stamp); `claims()` reports route support `[from_bp, to_bp)` and displayed-path support
`[evidence_from, to_bp)` separately, with `evidence_from` derived **per claim and displayed path, never from
`R.route_bp`** (§5.1); "contiguous occurrence" is reported only under `support == 't'` (from `H`); per-label end
reason = enum + qualifier; left-arm positions index `G` bases directly; completeness (`status`, `complete_to_bp`,
`scope`, `cap_trigger`) gates every comparison; label-list cuts carry totals everywhere and `A.nodes_truncated`.
Iterative deepening stays a backend call: the library derives the continuation and builds the request; it never
answers deepening locally.

## 5.1 Normative derivations

The library does not need the C++ to reproduce these; they are part of the format contract.

**Displayed-path evidence** (`Claim.evidence_from`, `Walk.labels_full`, `walks(by='support')` ranking,
`support_profile(kind='displayed')`). `max(from_bp, R.route_bp)` would be wrong after two non-first-parent merges at
`d1 < d2`: the run keeps `d1`, the entry and the leaf record `d2`, and `[d1, d2)` would be attributed to the label on
the displayed path. Provenance is derived **at the run's anchored endpoint** and then clipped — never by evaluating
the walk at a cut depth, which skips a later merge and attributes the displayed prefix to a label whose own route
spells something else (the clipped-merge fixture: B joins through a non-first parent at depth 14, the displayed path
begins `AACC`, B's own route `TCTA`; evaluated at `D = 4` it would give `evidence_from = 0`). For the run of a claim,
with `p` = the first-parent chain from the root to the run's **anchor** segment (also for merge-closed `m` runs) and
`t = to_bp`:

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
For leaf labels the merge-derived part (before the `max` with the run start) equals `T.route_bp` (asserted by T37);
interior claims use the same walk.

**Cuts.** *Run-start guard, applied first:* intersect the run's own interval `[from_bp, to_bp)` with the cut
`[0, D)`. If it is empty — the run starts at or after `D`, e.g. B entered by a switch at depth 10 and the cut is at
5 — the run makes **no claim** at `D` (it does not exist yet). Otherwise a claim cut at `D < to_bp` reports the
displayed-support interval `[evidence_from, to_bp) ∩ [0, D)`; when that is empty the claim is `route_only` at `D`:
the label carries a route of that length, but not the displayed bases; `routes()` reconstructs and spells the
label's own ancestral route through the partitions when the caller asks for it. A run of zero length (`from_bp ==
to_bp`: a label present at the boundary node only, e.g. ended `minority` at the seed boundary) is a **boundary
claim**: reported at its position with `kind: boundary`, no bases, and only when `from_bp ≤ D`.

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

**`label_summary`** (annotate), per arm, exactly `Walker::summarize_annotate()`:

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

## 5.2 Claims, routes, names and comparisons

1. **Annotate claims.** A claim is a maximal end of a label's routes under the union rule of §5.1. Its displayed
   support starts at the last merge its witness route enters through a non-first parent. The witness is the first
   carrying parent in stored order. One witness does not enumerate all routes, and does not establish a contiguous
   occurrence in a source sequence.
2. **Route stamps.** `R.route_bp` is the earliest non-first-parent merge on the lineage's route (a run closed by a
   merge of three or more parents may carry `route_bp == to_bp`). The `T` route stamp is the latest entry.
   Displayed support is `[evidence_from, to_bp)`, derived from the partitions — not `[R.route_bp, to_bp)`.
3. **Names.** `next_request()` refuses names shared within the retrieval and names containing U+FFFD (which servers
   write in place of bytes that are not UTF-8): a request built from them could name a different label. Relaxing
   the U+FFFD rule would need its provenance saved in the retrieval; a server's capability cannot vouch for an older
   body.
4. **Identity.** With a manifest on one side only, continue and replay refuse with `index_unverifiable` ("cannot be
   verified, not proven different"). A fallback needs the explicit opt-in `allow_unverified_index` (default off).
5. **Comparison at a depth.** The DAG is restricted to the comparison depth, and merges at or beyond it are left out
   (and named). Where that cannot be exact, the result is `qualified`; a side without bases gives `unknown`. Above
   64 labels, "qualified" can come from the per-node label limit (`max_labels_per_node`), which the docs of
   `compare()` state. Record coordinates are noted, never compared or clipped by `compare()`.
6. **The displayed parent at a merge** is the server's (§11); the library checks the rule as it stands
   (`derive.carried_labels`) and never re-ranks parents.

## 5.3 Continuations

A continuation is a new traversal (§14): its result is certified on its own, with the overlap stated.

- The loss budget is reduced by the largest terminal loss among the continued labels. This conserves the cumulative
  loss under unchanged costs: no continued route exceeds the original loss budget, and the result says so for the
  lower-loss labels. It is not equivalence with the uninterrupted walk: branch counts, edge history and other
  per-walk state restart in the new traversal (the notes state it), so a continuation can reach further than the
  uninterrupted walk would have.
- `max_label_branches` is reduced by the largest terminal branch count of the seeded labels (preserving
  `"unlimited"`); the request states it and offers an explicit reset (`reset_branches`).
- `labels.extra` is rebuilt around the new seed labels.
- Continuations that cannot share one request raise an error, with `next_requests()` as the fallback.
- Labels alive at a leaf but not seeded are named (`left_out`, with the reason: unreachable, or alive but not
  seeded).
- The request carries the record-coordinates setting forward (§18.3).

## 5.4 Local limits

Parsing (a compact MGT expands into large label sets), route enumeration, comparisons and exports can run under
work and allocation budgets (`budget.py`: `LocalLimits`, `LocalBudget`; `mcp_tools.ToolLimits` for the tools). The
rules:

- **Off by default.** An operation called without a budget charges nothing, and its output does not depend on the
  budget module. `max_bytes` bounds the bytes a tool returns, not its computation or peak allocation; a service that
  needs bounded tools sets `local_limits=ToolLimits(...)`.
- **Work is deterministic.** Local work units (lwu) are a weighted count of the model elements an algorithm visits
  and of the rows and text it produces (`WORK_MODEL` 2; 1 lwu is about 0.1 µs of CPython 3.11 on the reference
  machine). They are charged at the **cold price**: a derivation is charged its structural price on every call,
  whether a cache holds it or not, so the same call on the same graphlet with the same budget stops at the same row
  in any process and whatever earlier calls cached.
- **Memory is a modelled account** (`memory_bound: "model"`), not a measurement: soft in-process; hard bounds are
  the service's process limits. A budget's memory limit applies to each call on its own; its work accumulates across
  the calls that share it (an analysis allowance).
- **A stop says so.** Every charge is admitted before the work it pays for (with stated exceptions for what has a
  size known only once built). A stop raises `LocalBudgetExceeded` with the `LocalStop` that states it; lists whose
  order allows it return a page of whole rows with `complete: false`, `total_at_least` and a resume cursor;
  `compare()` answers `comparable: unknown`, never equality. No partial exports: an interrupted export or save writes
  nothing. Caches are stored only when complete, so a later unbudgeted call answers byte for byte as if the stop had
  not happened.
- **An interrupted parse is a local failure:** the original body (or the handle to it) is kept unchanged and the
  failure is reported; a fetch whose parse stops stores the body as `parsed: false`. A shallower graphlet cannot be
  manufactured from a cut-off parse.
- **Tool classes.** Every budgeted tool has a class (`TOOL_CLASS`: view, heavy, parse) with a default (view 30 M lwu
  / 1 GiB, heavy 300 M / 2 GiB, parse 200 M / 2 GiB) and a ceiling of 10× up to which the agent may raise its call's
  budget; every result carries `local {complete, usage, limits, work_model, memory_bound, clamped?, parsed?, stop?}`.
  Re-parse billing is the service's policy (`ToolLimits.budget_for`, `on_usage`).
- Quadratic paths in individual operations are fixed one by one (the cost model of `budget.py` states each
  non-obvious charge).

## 5.5 Files, the envelope and the client

- **The `J` line** carries the response envelope with `results[]` reduced to the seed's summary. Its `usage` is
  reduced to the totals plus this seed's `per_seed` entry (in `save()`, `dump()` and the store). `view` and
  `derived_from` are the library's own envelope names (`RESERVED_ENVELOPE_NAMES`): a response carrying one is
  refused. A `J` line that carries every seed's `per_seed` entry still loads, but is not canonical
  (`is_canonical()` false; `save()` rewrites it reduced): `J` is the library's record, not one of MGT v1's.
- `standalone_text(body, result, response)` writes a seed's standalone `.mgt` text (H, J, body) without parsing the
  body, byte for byte what `save()` writes for a server body.
- The client sends `not_after_ms` when given (§15.3), raises `AttemptExpired` on its 409, reads per-graph
  capabilities (`capabilities(graph=, graph_path=)`) and the server-wide document (`server_capabilities()`, §17.1),
  and keys features on `feature_level`.

# 6. MCP tool surface

Functions in `mcp_tools.py`; every return ≤ `max_bytes`; ids never leave the server, agents see `{name, ref}`; every
local answer carries an evidence block `{complete_to_bp, exact, support, reconverge, scope, outcome, limitations}`
(plus `view`/`qualified` on a derived handle; `exact` requires complete label evidence and no cut lists).

**Paging and sizes.** List tools (`graphlet_walks`, `_claims`, `_splits`, `_labels`, `_support`, `graphlet_list`)
take `cursor` and return `next_cursor` (opaque, bound to the entry handle, the tool and its normalized arguments; a
cursor presented with other arguments is rejected), plus `total`. `max_bytes` defaults to 2048 for list tools;
`graphlet_sequence` has its own ceiling (16 KiB, stated in its result and overridable) and is the way to get a long
sequence; `traverse_capabilities` has 32 KiB (the probe of a level-6 server with its coordinates block is 19,490
bytes on mini_refseq). Every error honours an explicit `max_bytes`; the floor is 64 bytes. A single row larger than
`max_bytes` is returned alone with `row_truncated: true` and the fields that were cut named (e.g. a walk's spelled
tail), never silently shortened. A filter never hides rows silently: walks and claims report `filtered` counts
(`route_only`, `merge_entered`). A successful export, save or load always keeps its receipt (handle or file
locator): optional fields are dropped first. `traverse_fetch(replay=<handle>)` re-runs the entry's stored normalized
request against the same index (same release required; a different release is an error, not a silent re-run).
*Oversize:* `max_graphlet_mb` is a RAM threshold — a larger body is spooled to disk complete and only the summary and
handle are returned; a hard transport/storage limit rejects the fetch explicitly. A truncated MGT document is never
stored or returned as a graphlet.

Backend-calling: `traverse_capabilities(index)`; `traverse_resolve(index, sequence, labels?|discover?, select?,
limit=10)`; `traverse_fetch(index, seed:{sequence, labels?, seed_id?}, strategy, keep=True, max_graphlet_mb=8) ->
summary + handle` (the probe is this with tight bounds); `traverse_continue(handle, arm, walk, overrides?,
execute=True) -> new summary + parent:{handle, arm, walk, overlap_bp}` (or the request when `execute=False`).

Local over a handle: `graphlet_summary(handle, arm?, detail?)` (`detail='limitations'`); `graphlet_walks(handle, arm,
rank=support|length|loss, n=10, min_bp=0, label?, spell=none|tail|full, tail_bp=60, route_consistent)`;
`graphlet_walk(handle, arm, walk, cursor?)` (segment chain, labels in/out per segment, continuation; the chain is
paged by `cursor` like a list); `graphlet_support(handle, arm, walk, step?)` (support runs with who left/joined and
the reason: a `K`-coded end or "split: N labels took C"); `graphlet_labels(handle, arm?, rank=direct_bp|reach_bp,
n=20, min_direct_bp=0, at_bp?, name?)` (table; `name` → that label's runs/routes incl. merge routes);
`graphlet_splits(handle, arm, n=10, min_labels_before=2)`; `graphlet_claims(handle, arm, min_bp=0,
route_consistent=True, labels?)`; `graphlet_sequence(handle, arm, walk|segment, from, to, orientation, with_seed)`
(≤ 16 KiB slice); `graphlet_export(handle, format=fasta|gfa|json|mgt, arm?, walks?, path) -> {path, bytes,
records}`; `graphlet_compare(a, b, arm, mode, cursor?) -> Comparison dict` (preconditions: same index identity
(§3.1) and the same oriented seed sequence, else `comparable: false`; differing completeness scopes →
`comparable: qualified`; returns the counts and the first page of `only_in_a` / `only_in_b` / `differ`, with
`next_cursor`; `graphlet_export(format=json, what=compare|walk, …)` writes the complete lists to a file — a truncated
field is always recoverable by cursor or export); `graphlet_subtrie(handle, labels, arm?) -> derived handle`;
`graphlet_list/free/save/load`. Unknown handle → "replay with traverse_fetch(replay=<handle>)" (the store keeps the
request). Rows carry `labels_in_total`/`labels_out_total`/`n_labels_total`. `graphlet_claims` rows always carry
record coordinates when the graphlet has them; `graphlet_walks` / `graphlet_labels` only with `coordinates=true`
(§18.3).

Errors are structured results with codes such as `unknown_label`, `backend_error`, `backend_unreachable`,
`io_error`, `path_not_allowed`, `not_replayable`, `view_unsupported`, `not_in_view`, `result_too_large`,
`receipt_too_large`, `local_budget_exceeded`. Exports and saves take file names confined to the export directory
(default `<spool>/exports`); the cursor secret is created in the spool when none is passed. There is no streaming
validator: a parse builds the model once and releases it, so the peak is the full model. Local limits (§5.4) apply
per tool class when the service sets them.

# 7. Size and memory

**Measured on the SRA primary index** (`scripts/traversal/graphlet_measure.py --server`):

| retrieval | body | gzip | full JSON (compact) | body / JSON | gzip / gzip | summary | parse | model |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 16S (1,396 bp), exhaustive constrain, 300 bp | 460 KB | 53 KB | 3.99 MB | 11.5 % | 27 % | 3.6 KB | 46 ms | 4.4 MB |
| 16S, label-free trie (annotate), 300 bp | 434 KB | 93 KB | 4.59 MB | 9.4 % | 28 % | 2.9 KB | 46 ms | 3.4 MB |
| 16S, default walk, 3 kb | 3 KB | 1 KB | 10 KB | 33 % (fixed overhead: seed, header) | 73 % | 2.4 KB | < 1 ms | < 0.1 MB |
| random 50-mer, label-free beam 20 (most supported), 10 kb | 3.5 MB | 196 KB | 75 MB | 4.7 % | 3 % | 3.6 KB | 0.9 s | 30 MB |

The body is below 15 % of the compact full JSON on the exhaustive cases and the summary below 10 KB. The interner
hit rate was 0.997 (constrain) and 0.83 (annotate, 2,264 distinct sets), which is why there is no per-arm alias
table. The small mini_refseq retrievals (18–29 %) are dominated by fixed per-document overhead.

**The model behind these numbers, an estimate and not a bound.** The model's memory depends on how many label sets
are distinct: 5,000 distinct 827-label arrays are ~17 MB for the arrays alone before model objects, interner keys,
names and caches, and a synthetic fragmented-set case produced ~17 MB of `SETEXPR` text. Per-record costs on the
constrained left arm of the 16S case (2839 segments / 1437 leaves / 1402 splits / 4160 ends / 69 kb bases): `G`
2839 × ~30 B + 69 KB bases ≈ 155 KB; `R` ~4200 × ~24 B ≈ 100 KB; `T` 1437 × ~8 B ≈ 12 KB; `E`/`V`/`B`/`A`/`S`/`L`
≈ 8 KB → ≈ 275 KB, ~60–90 KB gzipped. Rebuilt locally from it: paths 0.71 MB, splits 0.41 MB, labels_at_end,
label_summary, 1437 continuation strings. Label-free (annotate): `G` ≈ 155 KB; `P` ~5000 runs delta-coded ≈ 150 KB;
`L` 2836 file-path names front-coded ≈ 100–150 KB (0.42 MB raw: the intrinsic part); `T` 12 KB → ≈ 430–480 KB.
`labels.max_labels_per_node` bounds the label lists by construction, every cut keeping its true count. Server: the
body is one `std::string` appended straight from `SeedResult` (no `Json::Value` tree for it). Client: parse ~15k
lines in 50–100 ms; model ~5 MB per constrained arm, ~20 MB label-free (interned `array('I')` sets).

# 8. Where the code lives

C++:
- `src/cli/traverse.hpp/.cpp`: the request (`detail`, budgets, attempts, coordinates), the `GraphletWriter`
  (anonymous namespace: codec helpers `ranges/encode_setexpr/fmt_num/pct_escape/reason_code`, the records in
  document order, the `*` rule computations with their equality check), `graphlet_text(...)`, the per-seed JSON
  summary, `capabilities_to_json`, `kGraphletFormatVersion = 1`, `kTraverseFeatureLevel`; the additive JSON fields
  (`runs[].segment`, `segments[].labels_via_parent`, hairpin `followed`, revisit `length_bp`/`same_distance`,
  `label_dict[].column/seq_id`).
- `src/cli/traverse_attempts.hpp/.cpp`: attempts, the registry, tombstones and holds, the delivery reserve (§15,
  §17.3).
- `src/cli/server.cpp`, `server_utils.cpp`, `server_checks.hpp`: routes, per-graph identity and the manifest checks
  (§16), compression, the client-gone check, shutdown (§17).
- `src/graph/traversal/walker.hpp/.cpp`: the run fields (§4), head admission, budgets, deadlines, coordinates;
  `label_oracle.hpp/.cpp`: the budget-aware reads; `resolve.hpp/.cpp`: `/resolve`.
- `src/annotation/binary_matrix/base/decode_budget.hpp`, `row_diff/row_diff_budgeted.cpp`, `row_diff_cache.hpp`:
  the budget-aware decode path and the path cache (§14.3).

Python: `api/python/metagraph/traverse/` (§5), its tests `api/python/tests/test_traverse_*.py` and fixtures
`api/python/tests/data/traverse/` (codec vectors, whole documents, budget fixtures); the real-server suite
`api/python/tests/real/`.

Scripts (`scripts/traversal/`): `graphlet_fixtures.py` (regenerate + `--check` the CLI fixtures),
`make_codec_vectors.py`, `graphlet_measure.py`, `index_manifest.py`, `build_mini_refseq.sh`,
`make_column_coord_fixtures.sh`, `make_wide_coord_fixture.sh`. Measurements: `benchmarks/traversal/`.

Tests: `integration_tests/test_traverse.py` (§9), `tests/graph/traversal/test_trie.cpp` (anchor assert),
`tests/graph/traversal/test_trie_cases.cpp` (a bubble+switch+hairpin case whose dump feeds the Python fixtures),
`tests/cli/test_graphlet_codec.cpp` (codec, documents, budgets). User documentation: `docs/source/graphlets.rst`.

# 9. Conformance tests

C++ / integration (`integration_tests/test_traverse.py`, CLI and HTTP; the test venv installs `api/python` editable,
and the CI job installs the Python package so T37 cannot silently skip):
- T37 round trip: for constrain tuned (merge), constrain exhaustive (keep), annotate exhaustive, annotate beam,
  quorum, left-only, both arms, caps-tripped `max_steps`, `max_labels_per_node: 1`, 3-seed batch with a derivation
  failure: run `detail: full` and `detail: graphlet` on the same request; `Graphlet.from_response(r, out).to_json()`
  equals the full `results[i]` after normalisation (events sorted within equal `at_bp`, `needed_budgets` sorted,
  `timing` removed); `parse(text).dump() == text`; `A` counts match; `graphlet_lines == Z`; `check_rules()` reports
  no divergence; `#label_end events ≤ #runs`; every `Lw` run is the `prev_run` of at least one run entered by switch
  at its `to_bp` (§4; the switch events may sit on the anchor's children); `is_canonical(text)` holds for every
  writer output.
- T37 also covers these counterexamples: an ambiguous split with disjoint child sets after a switch, a divergence
  with overlapping child sets, a followed-hairpin split, a quorum-filtered split, an A→B→A re-entry inside one
  segment, a silent switch-source end at a split (switches on the children), two interior ends whose branch counts
  differ only through a quorum-rejected successor, two runs through two non-first-parent merges (`R.route_bp` ≠
  `T.route_bp`), and two same-name headers from different columns in annotate mode.
- T38 orientation: `spell(left) + seed + spell(right)` reproduces `acc1/acc2/acc3`; left
  `continuation(leaf).sequence == paths[].continuation.sequence` for continuations contained in the flank and for
  ones crossing into the seed, on both arms; resubmitting a locally derived continuation returns 200 and extends.
- T39 oracle through the library: `constrain.compare(annotate, mode='claims').equal` reproduces the constrain/annotate
  equivalence of the integration suite; `prefix_subset` on tuned-vs-exhaustive with `minority` as the omission
  reason.
- T40 protocol: summary key set; `strategy.output.detail/timing` echoed and resubmittable;
  `capabilities.graphlet_format == 1`; `sequences: false` → `G` without bases, `first_base` set, `C` with sequence,
  `splits[].char` reproduced; determinism (two runs byte-equal); gzip transport
  (`test_api_traversal_routes_are_compact_and_compressible`); body < ½ of compact full JSON on the fixture.
- `test_trie.cpp`: the §4 anchor assert on every trie case (covers merges/switches the integration fixture lacks);
  the `test_trie_cases.cpp` fixture (bubble under `merge`, constant-cost switch chain, hairpin) dumped via the CLI
  into the Python fixture set.

Python unit (`api/python/tests`): codec identities (random sets, shorter-wins, tie→explicit, floats, percent
escapes, front coding); fixture parse → `to_json()` ≡ `full.json` (normalised) and `dump()` ≡ text, byte-exact;
hand-built DAGs for `derive` (split, merge with partitions and `route_bp`, switch chain incl. a silent `Lw` end,
censored leaf, zero-length merged leaf); `claims()` kinds and `evidence_from`; `compare()` equal / prefix_subset /
incomparable; `subgraph(['acc2'])` keeps exactly RIGHT2's chain; FASTA/GFA parse back (S/L/P counts, k-1 overlaps,
paths spell the walks); `GraphletStore` LRU/TTL with a fake clock; truncated body raises;
`python -c "import sys; sys.modules['pandas']=None; import metagraph.traverse"` succeeds.

# 10. Compatibility with the other details

`detail: graphlet` is opt-in; `full/tree/summary` keep their shape (additive fields only; no test asserts exact key
sets of `segments`, `events`, `runs`, `label_dict`); the echo's two extra keys are accepted on resubmission;
`capabilities` checks are `assertIn`; `algorithm_version` stays `traverse-0.2` because the walk did not change for
the format — the format has its own version (`graphlet_format`). The right-arm JSON helpers of the integration suite
(`_walks`/`_structural_walks`) remain; the graphlet tests spell through the library.

# 11. The walker at a merge: the displayed parent

`Walker::merge_level` decides the merged segment's first parent (`parents[0]`), and everything displayed follows
first parents: a path's `segments` chain (`walk_path_leaf_first`), its spelled bases and continuation
(`make_continuation`, the graphlet's `C`), the end labels' `route_bp` (a label taken from a later parent is routed
from the merge) and the reconverge event's order. The graphlet writer, the JSON writer and the library read
`parents[0]`; no other place chooses a displayed parent.

**The rule: the first parent is the one carried by the most labels** (SPEC §7.1) — in constrain mode the labels
whose lineages its head brings into the merge node (`Item::state`); in annotate mode the fewest labels present at a
node of the parent's own segment (every head at the merge node holds the node's own labels, so only the parents' own
bases tell the routes apart). Ties keep the arrival order. Under trace support nothing merges, so record coordinates
are untouched.

**One exception keeps the arrival order: `tree` and `full` detail under a memory budget** (the request's or the
server's maximum). There every head reserves the delivery of its path's segment chain (`DeliveryCosts::chain_entry`,
exactly what the JSON writes), whose length follows first parents, so the displayed parent decides the account and
with it where a memory stop falls. The rule changes what is displayed, never the depth a budget certifies; applied
everywhere it would move stops shallower where the majority's chain is longer (on the cached real requests 52 of
6,552 budgeted responses; in a sample at 8 MiB in `tree` detail 10 of 294 mini_refseq requests, 1–4 levels shallower
on an arm, twice 30, never deeper). The other ways out are worse: charging below the output breaks the memory bound;
charging the longest parent's chain at every merge makes the account independent of the display but moves more
stops. `graphlet` and `summary` detail charge no chain (the graphlet names a segment's first parent, not its chain)
and follow the rule under every budget. So a `tree` or `full` request under a memory budget displays the arrival
order at a merge where the same request without one displays the majority; which labels reach a leaf, with which
loss and branches, is the same in both. On a server with a memory maximum (`--traverse-max-memory-mb`) every `tree`
and `full` response is budgeted and keeps the arrival order; only `graphlet` and `summary` show the rule. The
exception is the one condition `majority_first` in `Walker::merge_level` (§19 keeps the alternative open).

**The annotate rule compares the parents' own segments, not the walks they display.** A parent that is itself a
merged segment holds the union of its parents' labels, so after nested merges the rule can put first a short merged
segment whose displayed walk upstream is carried by fewer labels than the other parent's (pinned as the rule stands:
`MiniRefSeq.AnnotateMergeRanksParentsByTheirOwnSegment`, from the real cache's `mini_win200_02__annotate_merge`:
at the merge at 288 bp the 5-bp merged segment from 283 (6 labels) is first before the 32-bp segment from 256
(5 labels), though its walk goes on through a 7-bp segment of 4 labels). Net effect, over a fixed 10-bp window of the
displayed walk: UHGG's annotate merges whose first parent is not a majority are 68 of 558 under the rule against 191
under arrival order, and of the continuations it changes 18 of 1,350 lose labels and 22 gain (mini_refseq 6 and 8 of
28). A rule over equal stretches of the displayed walks would rank these cases right and would change the library's
check (`derive.carried_labels`) with it (§19).

**What the rule decides**, only on requests with `on_reconverge: merge` where a merge's majority arrived later (and
not in `tree` / `full` under a memory budget): the order of `parents`, `labels_via_parent` and the reconverge event's
segments; the paths' chains, spelled bases and continuations through such a merge (a continuation's labels with
them); the end labels' `route_bp` (the majority's labels 0, the minority's the merge depth); and, on a tie of loss
and branches, which parent's lineage continues (its run, its `labels_via_parent` entry, the closed run's segment) —
in the graphlet the `G` parents and partition, `R` and `T` records. **What it never decides**: which labels reach
each leaf, with which loss and branches, the arms' `complete_to_bp` (under every budget), outcomes, limitations,
`label_summary`'s direct and reach depths, work.

Tested (T55): `WalkerTest.MergeIsSpelledThroughTheParentWithTheMostLabels` (the majority through either branch,
constrain and annotate), `WalkerTest.MergeKeepsTheArrivalOrderWhereTheBudgetChargesTheChain` and the annotate pin
above; the CLI fixtures `merge` (the G allele, carried by b.fa, c.fa and both.fa, arrives second and is first) and
`merge_ties` (two labels on each allele: ties keep the arrival order). Against the cached real requests (2,184
unbudgeted and 6,552 budgeted, every detail) every difference the rule makes is a merge request and disappears once
the order-dependent fields above are canonicalised, and no budgeted stop moves.

# 12. Seeds: derivation and seed labels

**Partial delivery of a derived seed.** A derived seed (constrain mode with no labels named: its labels are the
intersection of its k-mers' rows, SPEC §6.1) whose time budget runs out after 1 ≤ j < n of its k-mers is delivered
as a result at depth 0 (its arms truncated at 0) with the set of the j k-mers read: a `derivation` limitation first
(`cause: time_budget`, observed j), `label_evidence: qualified` (the set can only shrink with more k-mers, so a
label in it may not carry the whole seed), and `coordinates: null` with the reason `partial derivation` (§18.5).
j = 0, or a superset that would fail the seed any other way, fails the seed (SPEC §7.0). The `derivation`
limitation's `observed` for `time_budget` has two units — the elapsed ms on a failed seed, j on a walked one — told
apart by the result's shape and stated with both units; a field such as `kmers_read`, or a cause of its own, would
be a wire change.

**Derivation windows.** A derivation reads its k-mers in windows of 64, whatever `annotation.batch_kmers` is. Each
window is one work charge before its k-mers are consumed, so a window's width places the work comparisons; were it
`batch_kmers`, a work-budgeted derived seed could walk partial at `batch_kmers` 1 and fail at 64. On a format whose
reads are not budget-aware a window holds 64 rows, and the soft memory observation states them (§14.2).

**Coordinate bytes of a derivation** (`coords_at`, under trace) are a running total. A sum over every consumed k-mer
at every observation would make a trace derivation under any memory budget quadratic (51 s against 0.16 s for 100
kbp). The total equals that sum exactly — it decides `beside` and `memory_bound_soft` — which a debug build asserts
and a test hook checks in Release, on a row-diff annotation through `beside` too.

**`labels_from_seed`** has one rule, `labels_derived_from_seed(seed, mode)` (walker.hpp): constrain mode and no
labels named. It is used for every result without a finished walk to read it from: a budget failure
(`SeedBudgetError`), a failed derivation, a seed never started.

**Seed labels and extra labels.**
- The duplicate checks are hash sets in the request's order — the seed labels' by name (a resolved label carries the
  name it was given), the extra labels' by target (the same kind and column, and for a header its sequence) —
  refusing the first duplicate in order, after resolving the names before it. A pairwise check would be quadratic in
  the names' bytes (2,500 header names with a shared 1,024-character prefix compare about 3 GB).
- Extra labels are accepted when every label of the declared pool is reachable from the seed labels within the loss
  budget through the cost model (transitively, not in one switch); the walk enforces cumulative costs.
- A seed that names headers maps only the coordinates in its own sequences' ranges; mapping every coordinate (rank
  and select per run) is kept for the queries where every sequence matters.
- The seed id over many names (their sort and FNV-1a, part of its definition) is linear in the names' bytes and runs
  as a `setup` piece (§14.4).

# 13. `/resolve` under a deadline

A `/resolve` with `bounds.time_budget_ms` reads its deadline between row batches and every `kResolveCheckKmers`
k-mers of the explicit labels' support pass. **A stop answers the resolve of a query prefix** (status 200, SPEC
§4.5): exactly the profile of the sequence's first x k-mers — a prefix, never a sample of the whole query, so a
stopped answer is an exact answer to a smaller question. Under a deadline the explicit labels' hits are fetched
`kResolveCheckKmers` k-mers at a time (one fetch of the whole query is a piece no clock read can end). The work after
the reads (a discovery's ranking and naming, the profiles, the grouping) checks every `kResolveCheckLabels` labels
whether the answer can still be built and written in time (else 503, deadline). The client is checked between the
phases and every 4,096 k-mers of the support pass.

**Memory.** Every `/resolve` decodes in batches (64 rows first, then about 64 MiB by the widest row of the batch
before, at most twice its rows and 4,096) and holds one batch; a discovery accumulates each label's k-mers, runs and
trace chain while it reads (no second profile pass), keeping a repeated k-mer's row for its later occurrences within
256 MiB; explicit labels prime their query with each distinct row once; a discovery indexes its column labels by a
4-byte slot per column. On a synthetic index of 8,000 wide columns the 0.3–3 kb queries run at 144–167 MB where
holding every row takes up to 1.66 GB, with 3–17% less CPU (header presence −7% to +14%); a 9.8 kb genome repeating
the locus three times runs at about 370 MB instead of 4.8 GB, faster for column presence, within 8% for explicit
labels and 24–48% slower for header and trace discovery (its repeats beyond 256 MiB decoded again). The profiles are
equal across batch sizes, byte targets and kept bounds, and equal to the same labels given explicitly.

**Explicit labels** are supported by one scatter pass over the k-mers instead of a scan per label (trace, 2,000
labels: 262 ms against 3,457 ms); presence uses the same accumulator instead of a labels × k-mers bitmap. The
per-k-mer hit copies of `fetch` remain the larger memory term with many labels.

Memory budgets for `/resolve` are not built (§14.3).

# 14. Resource and guarantee contract

**The principle.** Every response certifies what it covers and states every limitation: either an answer is
complete, or a valid shallower answer is returned together with exactly where it stops, what was cut, why, and which
knob or action would go further. Nothing is cut silently, and no resource measure may change the question being
answered without saying so.

**Outcome: four independent dimensions.** One value cannot say "partial traversal, spooled", and "all evidence
present" is stronger than any one boolean (an annotate run can reach its radius with every branch event kept while
capped label lists make its summaries lower bounds). Per seed result (JSON `outcome`, MGT `O`):

| dimension | values | guarantee |
|---|---|---|
| `walks` | `complete` · `partial` · `failed` | `complete` only when no walk-class limitation applies: every requested arm is complete to the radius **per path** (keep, or merge where no history was united) and no carrier of the seed was dropped; `partial`: any of `walk_domain`, `scope` with at least one merge, `seed_labels` applies — a valid certified prefix whose limits the `limitations` state; `failed`: no valid traversal exists (structured reason in `limitations`/`resource_stop`) |
| `branch_diagnostics` | `complete` · `cut` | every branch decision and refusal is reported; `cut`: only before each arm's `evidence.complete_to_bp` (`branch_events`) |
| `label_evidence` | `complete` · `lower_bound` · `qualified` | `complete` only when no label-class limitation applies; `lower_bound`: evidence may be *missing* or understated — cut lists (`label_lists`, `inexact_counts`), dropped carriers (`seed_labels`), a cut switch-source list (`switch_sources`), greedy losses (`greedy_losses`); `qualified`: something reported may be *overstated* — a column-label trace across a record boundary (`trace_record_boundaries`), or a permitted set derived from part of the seed (a `derivation` on a walked result, §12); `qualified` wins when both apply |
| `delivery` | `inline` · `spooled` · `paged` | the body is in the response; `spooled`: complete in the spool behind the handle, the response holds the summary; `paged`: delivered in pages. Independent of `walks`, so a partial traversal can be spooled |

**Conservative by construction**: each dimension is `complete` only when no limitation of its class applies, so a
reader can never find a limitation in `limitations` whose dimension still reads `complete`.

**Stated limitations** (SPEC §7.0): per arm `evidence: {complete, complete_to_bp}` and `limitations: [{kind, knob,
limit, observed, effect, complete_to_bp?}]` with kinds `walk_domain`, `branch_events` (the stored branch events stop
at `max_branch_events`; every branch decision and refusal before the first dropped event's depth is present — a run
cut there says so and names `output.max_branch_events`, which also accepts `"unlimited"`), `label_lists`,
`switch_sources`, `inexact_counts`, `scope`; per seed `seed_labels`, `server_clamp`, `derivation`, `coordinates`,
`memory_bound_soft`. MGT carries them: `A`'s `<evidence_complete_to_bp|*>` and one `K` record per limitation (§2.2),
so a saved graphlet keeps its own caveats.

**Budgets enforce the bound; depth is the graceful fallback.** A depth limit alone does not protect against a very
wide first level or one enormous annotation row (annotate mode materialises a full row, or every coordinate tuple,
before cutting the recorded list — label_oracle.cpp). So:

- *What counts.* **Memory**: every allocation the request retains at once — the walker's state (heads' label states,
  segments, runs, events, recorded sets), the request's `LabelQuery`/`LabelRecorder` caches, the decode scratch of
  the annotation reads it issues (row-diff dependency rows, coordinate tuples), and the output buffers. The index
  itself (mmapped graph and annotation) and process-wide shared caches are not charged to a request; a shared cache
  that grows because of a request is charged to it. **Work**: *charged work units*, a weighted sum of counters the
  walker already keeps — successor enumerations, rows decoded (weighted by their row-diff dependency path length),
  coordinates mapped, pair evaluations, refusal scans, edge-reuse probes, re-minimisation rounds — not CPU time and
  not elapsed time; `max_steps` alone undercounts (one step can perform many source–target comparisons). Elapsed
  time stays a separate deadline. The deadline and the budgets are checked at least every *W* work units (a bounded
  interval, stated in `capabilities`; §14.2 states what the comparison after every charge guarantees).
- *Row decoding is budget-aware inside the decoder.* Estimating the *final* row does not bound peak memory: row-diff
  reconstruction loads its dependency rows, whose differences can cancel into a tiny result (`row_diff.hpp`). The
  annotation reads used by traversal go through a decode path that charges each dependency row and tuple as it is
  materialised and aborts the batch, without side effects, when the request's remaining budget would be exceeded;
  the step then stops with `resource_stop.phase = annotation_decode` and the suggestion to use a more selective seed
  or a label-constrained query (narrowing the radius would not help). The path exists for the row-diff formats
  (§14.3); the memory bound stays *soft*, and every result under a memory budget says so (`memory_bound_soft`,
  observed = what was held beyond the account, §14.2).
- *Atomic commit per head.* "Check before each allocation" is not enough when earlier mutations already happened: the
  walker ends a switch source's run before it allocates the split's children, so a denial in between would leave
  inconsistent runs, events or split records, and lowering `complete_to_bp` does not repair them. Processing a head
  is therefore two-phase: **plan** (successors, derived states, label ends, splits, merges — no mutation of the
  result) → **admit** (reserve the plan's allocations *and* the cost of later stopping and delivering what it
  creates, see the next bullet) → **commit** (apply; cannot fail). A head that is not admitted is censored exactly
  like a head beyond a cap (its walk ends at the last committed level with `resource_limit`), and the level it
  belongs to does not count toward `complete_to_bp`.
- *Delivery is accounted per expansion, not reserved as a fraction.* A fixed reserve fails on a comb-shaped trie (one
  terminating branch per split): segments grow linearly, but materialising every leaf's ancestor chain is
  quadratic. So (a) paths are **lazy** — `PathResult` keeps its leaf; ancestor chains are produced while serialising,
  from parent pointers, through bounded buffers (MGT stores no chains at all; JSON `detail: full` streams them) — and
  (b) admission charges each committed head with the incremental delivery cost of what it adds in the requested
  output format. The guarantee: **every committed prefix remains deliverable within the resources already reserved
  for it.**

**Scopes.** *Locus*: peak working memory and cumulative work across both arms, retries and continuations of one
locus. *Analysis*: total work across loci, an overall deadline, retained graphlet memory and spool storage.
*Service*: concurrent workers and aggregate memory reservations. Requests carry a stable `budget_id` (analysis), a
`locus_id`, and a fresh `attempt_id` per dispatch; splitting a request or retrying never resets an allowance —
memory limits what is held at once, work accumulates across attempts.

**An authoritative ledger with reservations**, not reported usage: two requests that each receive the same remaining
allowance before either reports would both spend it. The MCP layer / search service keeps a ledger per `budget_id`
and `locus_id` (atomic operations, e.g. Redis) and **reserves** an allowance for an attempt before dispatching it;
the reservation is the request's bounds; the backend returns its usage on every response, successful or not; the
ledger **reconciles** reserved against used on the response, and **releases** a reservation on cancellation, or on
lease expiry when a response is lost (the attempt is then charged its full reservation, since its spending is
unknown). The backend enforces the locus scope per request; the ledger enforces the rest. The backend's half is §15.

**Charging work and releasing capacity are separate operations.** Consumed work is charged when the response (or the
lease expiry) settles the attempt; occupied capacity — the worker slot and the memory reservation — is released only
once the backend has actually stopped: on its response, or after the lease, which is only a valid release point
because the backend enforces the same deadline itself (the request's `time_budget_ms` plus the server's hard request
timeout, so an attempt cannot outlive its lease). A cancellation the backend has not acknowledged releases nothing.
The bound is compared only at the delivery checks, so an attempt outlives it by its run up to the next one, whose
length has no stated bound (§15.2); the lease is a valid release point under that stated assumption, and an answer
of `running` or `stopping` past it shows that run still going.

**No silent semantic changes.** Enabling a beam, dropping labels, raising a cut, switching `on_reconverge` or the
support kind changes the question or the evidence; none of them is ever applied as a resource measure. When the agent
chooses one, the response says what it changed (`strategy` echo + `limitations`/`scope`).

**The local library has budgets too** (§5.4): an interrupted comparison returns `comparable: unknown`, never
equality, and an interrupted parse is a **local failure** that keeps the original body.

**Continuation is a new traversal.** It does not resume the original exhaustive search (the backend keeps no frontier
or edge history between requests); its result is certified on its own, with the overlap stated (§5.3).

**Fixtures for this section.** Exhaustion during annotation decoding (including a row-diff row whose dependencies are
dense but whose result is tiny), mid-level traversal, finalisation and serialisation: a valid shallower result is
delivered whenever one exists (`walks: partial` with the right `complete_to_bp` and `resource_stop`), `failed` only
when none does. Allocation denial injected inside a head's commit around a switch, a split and a merge (the result
stays consistent: runs, events, splits and paths agree); a comb-shaped trie (finalisation stays within the admitted
reservation); concurrent attempts against one ledger entry (no double spending; a lost response charged its
reservation); a partial traversal delivered spooled (`walks: partial`, `delivery: spooled`); a non-zero
`refusal_scans` and an unknown future counter round-tripping through `A`; an interrupted Python parse keeping the
original body and reporting a local failure.

## 14.1 What is frozen with MGT v1, and how enforcement is staged

**Frozen with the format** (the wire contract): the `O`, `Q`, `K` records and the `A` counters/evidence fields; the
`resource_limit` end reason; the four outcome dimensions; the `limitations` kinds; the `resource_stop` structure;
`attempt_id`/`budget_id`/`locus_id` in requests and usage in every response. Frozen means that no kind a reader knows
changes its meaning, knob or value types. The kinds themselves are an **open `[a-z_]` token set**, as are a `K`
record's extra field names and a `Q` record's scope, resource, phase and action tokens: a later feature level may add
one, and a reader keeps a value it does not know, uninterpreted (examples: the kind `memory_bound_soft`; the kind
`coordinates` with the extra field `lists_cut` and the action `drop_coordinates`; the attempt tokens of §15.1). A
limitation's `effect` sentence is free text and may be reworded within MGT v1. No record, field or grammar rule
changes; that would be MGT v2.

**Enforcement is staged**, each part stated in every response until it is enforced (`limitations` entries). The
stages, numbered as the SPEC and the tests name them:

1. the caps, the evidence boundary, the stated limitations and the outcome dimensions — enforced;
2. the two-phase head commit and lazy paths in the walker, with the request budgets `bounds.max_memory_mb` and
   `bounds.max_work_units` — enforced (§14.2);
3. the budget-aware decode path in the annotation library — enforced for the row-diff formats (§14.3); the memory
   bound stays soft (`memory_bound_soft`);
4. the ledger — the search service's; the backend's half (attempt ids, usage, cancel, state, the release rule) is
   built (§15);
- (L) work and allocation budgets for the local library's operations — built, off by default (§5.4).

The promise at every stage is the same: deliver the largest certified prefix that can be safely committed and
delivered, and state the limiting resource and the useful next action — the stages only move hard bounds from
*stated* to *enforced*.

## 14.2 Work and memory as built

- **Per seed.** Budgets apply per seed: a request of n seeds can hold n × budget. A seed's stop does not depend on
  its position in the request.
- **Work is stated with its number.** Work is charged when a fetch returns a row ("rows decoded": 8 units per key, 1
  per entry and per coordinate), plus label-state scans, pair evaluations, edge-reuse probes and enumerations. The
  budget is compared after every charge, so a stop exceeds it by at most what was charged since the previous
  comparison, and every work stop states the most its seed charged between two comparisons (`largest_charge`).
  `W` = 65,536 (`kWorkCheckInterval`) bounds the seed phase's interval and sizes fetch calls. The roots' rows of both
  arms are one charge, so that a stop still delivers a result complete to 0 bp. Rows the lookahead decodes without a
  fetch asking for them are not charged, and the SPEC says so.
- **Work is deterministic logical work.** Each returned row is charged its row-diff dependency path (8 units per
  dependency row plus 1 per stored entry), whatever was shared in the actual decode; a row's logical charge is its
  standalone charge whether decoded, cached or shared. It is not measured decode effort, which depends on batching
  and caching; the physical decode counters are in `timing` (and `timing.path_cache`). `capabilities.work_bound` says
  so.
- **Memory is a deterministic model, never measured RSS.** Delivery is priced at worst-case escaped lengths, once per
  place a name is written. The model's sizes of `LabelRun`, `Entry` and `Item` are part of it (a larger struct
  changes every demand). What a fetch, a validation or the depth-0 dictionary holds beyond the account is recorded
  before any check can stop the walk, and is stated as `memory_bound_soft.observed`. `memory_bound_soft` is stated on
  every result under a memory budget, failures included, and names what is still held uncharged: a level's lists and
  fetched rows until its heads are processed, the seed phase's intersection and hits, a failed seed's depth-0
  dictionary, an index-wide header lookup, the dictionary's first table, and a failed response's `seed_id` echo. A
  failed result echoes index-supplied names cut to 256 bytes under a memory budget; the request's own `seed_id` is
  stated, not cut. There is no hard request-wide memory bound.
- **The header index** (header → column, sequence) belongs to the `CoordToHeader` object and lives as long as it
  does; there is no address-keyed cache. The server builds it while the index loads (§17.4).

## 14.3 Annotation reads under a budget

- **Charging inside the decoder.** On the row-diff formats (`RowDiff<BRWT>`, `RowDiff<ColumnMajor>` and the
  coordinate-aware `TupleRowDiff<...>`; SPEC §6.8 lists them) the traversal's annotation reads go through an opt-in,
  single-threaded decode path (`IRowDiff::decode_rows` / `decode_row_tuples` with a `DecodeBudget`). It charges the
  trace, each dependency row, each coordinate tuple and the reconstruction buffers *before* they are allocated, and
  refuses a read whole, without side effects, when it would not fit. The default decode path is unchanged.
- **Batch independence.** Every key a fetch returns is admitted against its standalone demand, so where a walk stops
  does not depend on `batch_kmers` or the lookahead. A stop states only facts that hold either way.
- **The row-diff path cache.** A read keeps decoded rows so that later reads stop their row-diff paths at them: the
  rows asked for, the 8 after each, the checkpoints (distance to the anchor a multiple of 16) and the narrow rows
  (under 4,096 bytes), tuple rows kept flat, as full rows rather than raw diffs. Keeping every row of every path would
  copy each wide tuple row of a path (one allocation per column with more than two coordinates): on a synthetic index
  of 8,000 labels with 24,000-coordinate rows that made first reads 1.6–6.2 times slower than keeping fewer rows, and
  a walk reading one row per call took 16.9 s against 2.2 s with the cache off. With the retention rule first reads
  are 1.1–2.9 times faster than with the cache off and one-row walks 6–15 times; budgeted reads keep exact costs. The
  stored size of a row is computed before it is copied (`StoredRow::bytes_of`, equal to the stored copy's `bytes()`,
  asserted); the older entries are evicted first, and the row is copied only when admitted, so a bound never lets the
  cache allocate a row it then refuses, and no insert holds more than its end (`timing.path_cache.peak_bytes` is what
  it held). Under a memory budget the cache lives in the label cache's allotment (`make_room` gives it what the label
  cache leaves).
- **The budget-aware lookahead** takes its runs and pieces in the walk's order, each sorted, as the default reads do
  (31,297 decoder charges against 161,679 in unsorted order on a stress fixture). It stops at the first run the cache
  cannot keep beside its own runs (by count before decoding, by bytes after), evicting what was cached before it for
  its first run only; the level's fetch keeps its wholesale eviction. A later warm's first run still evicts wholesale,
  and with it rows earlier warms read ahead that the walk had not reached: up to +20% tuple rows over the walk under a
  work budget alone at 16 MiB, +52% at 32 MiB, +32% at 64 MiB (`batch_kmers` 2,048–8,192; SPEC §6.8 states it).
  Evicting the oldest entries first instead (a list through the costs map, a warm keeping a quarter of the cache)
  reads 12% fewer tuple rows but keeps the label cache full, and the row-diff path cache then read 2.5 times the
  stored rows (`annotation_fetch_ms` +44%); evicting the oldest half read 2–7% more rows at the default `batch_kmers`
  64 and 2.8 times in one walk at 2,048. The comment at `cache_budgeted` says why wholesale eviction stays.
  `MiniRefSeq.LookaheadKeepsItsRunsUnderAMemoryBudget` holds a walk within 25% of the work-budget walk. The structural
  lookahead's own re-warming of the chains it clears is not bounded this way (5.4 million rows at `batch_kmers` 8,192
  under 8 MiB).
- **Stops.** A refused read in a level stops the seed in phase `annotation_decode` (actions `more_selective_seed`,
  `label_constrained_query`); in the seed phase, derivation or an annotate root it fails the seed. Labels that do not
  fit, and level lists that leave nothing for the read, are `traversal`-phase stops (§14.5).
- **Not covered:** the formats without row-diff (BRWT, ColumnMajor and TupleCSC alone, the disk and flat row formats,
  `IntRowDiff`) read through the default path, stated; `/resolve` has no memory budget; a level's rows and lists are
  not charged in the account. Agreed for when this is built: every format except the four legacy ones (rbfish,
  rb_brwt, bin_rel_wt, row/EigenSpMat), which stay stated as uncovered; ColumnMajor read as column-wise charged pairs;
  work-only budgets keep the default decode on formats without row-diff; each head's rows released once it is
  processed; a head stop caused by a level's rows states the bytes, without the `more_selective_seed` action; the disk
  reader charged as a measured constant; a stopped `/resolve` returns exactly the resolve of a query prefix, an
  explicit `select.policy` that stops fails the request, and discovery under a budget keeps every met label's runs; no
  server-side default memory budget for `/resolve`, and a cap on its explicit label list (`--resolve-max-labels`);
  `annotation_reads` and `resolve_budgets` stated on `GET /traverse/capabilities` only. Selected-label decoding
  (decoding only the selected columns) is agreed with these rules: the access path stays `tuples`, the selected mode
  reported in timing and capabilities; the binary `rd_direct` accessor (SPEC §8.2) built with it; a selected read
  charges 8 per dependency row, 1 per BRWT probe and 1 per selected coordinate (and the coordinates examined, shifted
  or cancelled during reconstruction, not only the survivors); a fall-back to full-row decoding decided
  deterministically (the 64-column crossover is a tuning point, not a threshold); header-range filtering inside the
  decode shifts the ranges along the dependency path; a tuple column's presence is not the XOR of Boolean
  nonemptiness. The program and its measurements are `PLAN-traverse-next-stages.md`.

## 14.4 Deadlines and uninterruptible pieces

- **Chunked reads.** Under any deadline, every `/traverse` annotation read (a level's fetch, the lookahead, a
  validation, a derivation's window) that the deadline may fall into is decoded in chunks, the deadline checked
  between them. A read is one piece when, at the slowest per-row time the request has seen, it would take less than
  1/64 of the time left; otherwise its first chunk is at most 8 rows, each next at most 4 times the previous and sized
  at the rate the previous measured to take `chunk_target_ms` (`--traverse-chunk-target-ms`, 50 by default), and the
  rest is one piece once predicted at that rate to take less than a quarter of the time left; chunks are taken in the
  walk's order. Why per read and not per row: on a row-diff annotation rows share the decoding of their row-diff paths
  within a call, so cutting sorted keys into chunks sized from the slowest recent row decodes the paths again per
  chunk (0.75 ms a row against 0.016 ms whole: a row_diff walk of 1.8 s took 30 s), and a first chunk sized from
  another level's cheap rows overshoots a 50 ms target many times over. A read's counting, cache eviction and result
  are decided once for the whole read, so a request no deadline stops is byte-identical (checked on real requests,
  also with 1 ms chunks and with every read split into one-row chunks), and a stopped read censors the walk at the
  read. The rows a stopped read decoded are charged. A validation is stopped by the attempt only, never by its own
  time budget. A read maps its rows' coordinates inside its piece, so the rate a chunk is sized by includes the
  mapping. Measured: on the fan-out fixture budgets of 50–500 ms are kept to within 15 ms with a column annotation
  and 60 ms with row_diff and row_diff_brwt (with and without a memory budget), and walks that no deadline stops keep
  their bytes and time (60 real UHGG requests: 15,973 ms chunked against 16,021 ms unchunked). Stated limits: a whole
  read overruns when its rows are more than 64 times slower than any seen before, the rest of a split read when more
  than 4 times slower than its chunk's.
- **Pieces.** One chunk (at least one row) is uninterruptible, and so is a read decoded whole far from its deadline
  (a cancel is seen after it). The seed phase is cut into pieces wherever another piece (a read or chunk, a k-mer
  mapping, a derivation step) begins or ends, and every span between two of them is a **`setup`** piece
  (`DecodePacer::open_setup` / `close_setup`, the walker's `begin_setup` / `end_setup`); the phase ends at the walk's
  first checkpoint, or at its stop or end when it reaches none (a radius of 0, a partially derived seed, no arm), and
  is flushed on every way out, a failed seed's included. The walk after its last reading of the clock, to its stop
  or end, is a last `head` piece, so the pieces cover the seed's time from its start to its result.
- **The deadline record.** `timing.deadline` states the seed's longest uninterruptible piece (kind, rows, coordinates
  mapped in it) and what stopped the walk how long after its deadline; `timing` states the seed phase
  (`seed_phase_ms` from the seed's start to the phase's end, `seed_fetch_ms` its reads, `label_resolve_ms`).
- **Polls.** The walk's polls compare the walk-until (§17.3); one poll in `poll_stride` (8) reads the clock, besides
  the forced ones (between seeds, before a paced read's chunk, in the lookahead). The lookahead's chains read the
  attempt's stop every 16 graph steps (`kLookaheadPollSteps`) and before each chain's key mapping, through the poll
  that reads the clock (`poll_now`), so a passed walk-until is recorded with the attempt and the next checkpoint stops
  the walk; a stop ends the lookahead, and each poll ends a head piece. Without these polls a lookahead at
  `batch_kmers` 60,000 saw a cancel 1.5–11.8 s late. The delivery checks come every 4096 objects of a result or of the
  MGT text, every 64 KiB of text (a large token too: one token's preparation stays a piece), between compression
  blocks and before the transport.
- **No time bound on one piece is stated:** `max_uninterruptible_ms` is null, and it stays null — checkpoints inside
  reads bound the index operations of a piece, not its wall time, which page faults and scheduling leave open (SPEC
  §6.8, §10.3). The longest piece is stated as an observation: `observed_max_uninterruptible_ms` per process in
  `deadline_check` (reads, head pieces and the delivery gaps) and per attempt in `usage` (reads and head pieces).
  `setup` pieces are in the deadline record only (§19). `deadline_check.rule` names what stays unpolled: a chain's
  key mapping, a read's preparation before its first chunk (its cache lookups and the ordering of its keys: 90–260 ms
  for the 1 to 3.8 million keys of a lookahead in a stress fixture), the lookahead's clearing, the seed phase's own
  processing.

## 14.5 Stops and what they state

- A `resource_stop` is created only when a head below the radius was actually censored. A memory bound is not a
  wall-clock bound, and delivery is not checked against a request's own deadline (the attempt bound is, §15.2).
- **Memory stops** state the smallest budget (whole MiB) whose account holds the need beside that budget's own cache
  allotments, so raising the budget to the stated value holds the state.
- **Lower bounds.** A depth-0 failure from unread roots or unnamed labels states its need as "at least X"; raising the
  budget to X may still fail, because unread roots, labels or later state need more, and the statement says so.
- **The dictionary stop on formats whose reads are not budget-aware.** An annotate level's new labels are charged
  after its read (`charge_dictionary`). When that charge is what crossed the budget (the account within it before),
  the trip is `LABEL_NAMES` at the level, as a lower bound (`ResourceStop::names_after_read`): its message names the
  labels and their bytes, its actions name the levers that name fewer labels. Charging each call's names inside the
  read, as the budget-aware read does, would stop such a level earlier and move where walks stop.
- **Failed seeds.** A seed a budget does not hold even to 0 bp (`SeedBudgetError`), and a seed the attempt never
  started (phase `not_started`, its action a new attempt), are failed in the shape of a failed derivation (no arms,
  an error, `outcome.walks: failed`) with a seed-level `walk_domain` naming the budget's knob or the attempt, and the
  `resource_stop`; the other seeds of the request are traversed. The `walk_domain` is stated although it is
  otherwise an arm's, because no dimension may be non-complete without the limitation that explains it.

# 15. Attempts, cancellation and the release rule

The ledger lives in the search service (§14). The backend provides what a ledger needs; SPEC §5, §7.0, §7.3 and
§10.3 are normative, and `GET /traverse/capabilities` names them (`attempts`: each rule a reference to the SPEC
section that states it).

## 15.1 Attempt ids and usage

- A request may carry `attempt_id` (`^[A-Za-z0-9._:-]{1,128}$`), plus `budget_id` and `locus_id`, which are echoed
  and need an `attempt_id`. An id is unique per server process while its attempt runs, is retained or is held
  (§15.4): a reused id is refused with 409, whose body carries the id's state as GET answers it, so nothing runs
  twice. Requests without `attempt_id` keep their bytes.
- **Usage.** Every response to a request with an `attempt_id` carries a response-level `usage` block: the ids;
  `server_instance`; per seed and in total, work units, modelled memory (admitted, final, soft excess,
  `held_bound_bytes`, which is null without a memory budget because the account is then no bound on what is held),
  elapsed ms, and seeds started, finished and abandoned; the bound actually used; and
  `observed_max_uninterruptible_ms` (§14.4). This covers success, partial, failed and cancelled seeds, and 400/500
  after registration. The 409, the 400s before registration and the loading 503 carry no usage; a client that is
  gone gets no response.
- **Free-token values** of the attempts, which MGT v1 admits (§14.1): resource `cancelled` and `attempt_deadline`;
  scope `attempt`; phase `not_started`; action `retry_attempt`; knob `attempt_id` on the stop's `walk_domain`.

## 15.2 The attempt bound

The attempt bound = min(seeds × the effective per-seed `time_budget_ms` + allowance
(`--traverse-attempt-allowance-ms`, 10 s by default), `hard_cap_ms`, derived from the HTTP server's content timeout,
~899 s). The server enforces it itself:

- seeds stop being walked at the walk-until, `bound − max(allowance/2, delivery reserve)` (§17.3);
- building, writing and compressing the response are checked against the bound (the delivery checks of §14.4);
- an attempt at its bound answers 503 `{error, usage}` with reason `deadline`, delivering no seeds; its partial
  results come as a 200 instead (`resource_stop {scope: "attempt", resource: "attempt_deadline"}`, unstarted seeds
  failed with phase `not_started`).

**Enforcement is cooperative**, and only the delivery checks compare the bound; the walk's polls compare the
walk-until, which lies below it, and only a poll that reads the clock sees it (§14.4). So an attempt runs past its
bound until its next delivery check: the rest of the piece the bound fell into (a read's chunk, a head, the mapping
of a seed's k-mers, the seed phase, a finalisation, a gap between delivery checks) and, when its walk had not
stopped, the walk up to its next clock-reading poll (up to 7 heads), the stopped seed's finalisation and the
building of its result up to the first check. No piece has a stated bound (`max_uninterruptible_ms` null) — on a
6.94 Mbp seed, GET answered `running` 2.7 s past the clock release's instant. `attempts.bound`, `not_after`,
`release_rule` and `deadline_check.rule` say that the attempt is past its bound apart from that run, whose length has
no stated bound, and that an answer of `running` or `stopping` past the instant shows it.

## 15.3 `not_after_ms` and `expect_server_instance`

- **`not_after_ms`**, an optional top-level request field (Unix epoch ms, an integer up to 2⁵³ − 1, with or without
  `attempt_id`): the instant after which the request must not be started. The server refuses, when the handler
  starts, a request whose instant has passed on its clock — 409 `{error, state: "expired", not_after_ms,
  server_time_ms, ids, server_instance}` — with nothing run and nothing registered (a later GET is a 404). The check
  is strict; the allowance a ledger adds is stated (`attempts.clock_skew_allowance_ms`, `--traverse-clock-skew-ms`,
  2000 by default). A running, retained, held or tombstoned id is answered first, so `expired` always means that no
  attempt with the id exists there. Echoed in usage and in the state. The CLI applies the same check.
- **`expect_server_instance`**: holds and tombstones live in the process's memory, so a restarted process holds none.
  A request may name the instance it expects; a process with another `server_instance` refuses it before anything
  runs (409 `instance_mismatch`; a malformed value is a 400).

## 15.4 Cancel, state, tombstones and holds

- `POST /traverse/cancel {attempt_id, wait_ms, not_after_ms?}` answers 200 when it flagged a running attempt; 404 when
  the attempt has finished; 404 with `tombstone: true` when the id is unknown — the id is tombstoned, so a later
  request with it is refused (409); 429 when no tombstone can be held (`tombstones_full`; with `retention_s` 0,
  `tombstone: false, reason: "no_suppression"`), which promises nothing.
- `GET /traverse/attempt/{id}` returns `running`, `stopping` or `finished` with the reason (`completed`, `cancelled`,
  `client_gone`, `deadline`, `error`). It creates nothing, so its 404 is not completion.
- **Retention.** Finished attempts stay queryable for `retention_s` (`--traverse-attempt-retention-s`, 3600 by
  default) and at most `retention_count` of them are kept (`--traverse-attempt-retention`, 10,000; oldest dropped
  first). A cancel of an unknown id is tombstoned only while fewer than `retention_count` tombstones and held
  attempts are held, else refused (429) — never by dropping the oldest. The retention seconds, the count and the
  tombstone cap are bounded integers (a year; ten million); any other value refuses to start, naming the option and
  its range.
- **Holds.** A tombstone is held at least `retention_s` and until the server's clock reads `not_after_ms +
  clock_skew_allowance_ms`, capped by `--traverse-attempt-tombstone-max-s` (86,400 by default; never less than
  `retention_s`). A tombstone held for `retention_s` alone would not justify releasing capacity: a request half
  uploaded when it is cancelled can complete after a short retention and run. Every tombstone answer — the cancel's
  and a repeat's 404, the state's 404, and the 409 refusing a copy — states `suppressed_until_ms`, the `not_after_ms`
  it was judged against and `covers_admission` (true iff a `not_after_ms` was named and the expiry reaches it plus
  the skew; else `covers_admission_reason`). A repeat cancel never shortens a tombstone and extends it to its own
  `not_after_ms`; a refused copy extends it to the copy's own `not_after_ms`, so every later copy — a proxy's replay
  carries the same instant — is refused while it could still be admitted. A refused copy is judged against its own
  `not_after_ms`, never a cancel's (without one: `covers_admission: false`, `no_not_after_ms`).
- **A finished attempt sent with `not_after_ms` is held like a tombstone**, its id refused, until its `not_after_ms`
  + skew (within `tombstone_max_s` of the finish, extended by a refused copy's later `not_after_ms`), past its
  retention if need be, among the tombstones and on both clocks. Held attempts are never dropped early, so
  `retention_count` does not bound them; every refused copy, an identical replay included, extends the hold to
  `tombstone_max_s` from its own arrival, so the finishes and the refused copies within `tombstone_max_s` bound them
  (conservative for replay safety).
- **Clocks.** A hold is live through `suppressed_until_ms` inclusive (the strict `not_after_ms` check still admits in
  that millisecond; the steady hold is one ms longer) and while either clock says it is live (the steady clock for
  the computed duration, the wall clock until it reads the expiry), so a forward step of the wall clock cannot
  shorten it; a backward step by **more** than the skew after the expiry can revive a copy's instant (the ledger's
  skew assumption).

## 15.5 The release rule

Normative, stated in SPEC §5 and §10.3 (the capabilities' `attempts.release_rule` refers there); the library applies it in
`attempts.release_verdict`, which lists the assumptions of each release (e.g. `sent_without_expect_server_instance`).
A ledger may release an attempt's capacity:

- **early, on a tombstone**: only on 404 + `tombstone: true` + `covers_admission: true` from the same
  `server_instance`, for an attempt sent with exactly that `not_after_ms` and with `expect_server_instance` naming
  that instance. An attempt sent without `not_after_ms` is never released early;
- **on a finished state** (GET, or the 409 refusing a copy, which carries the id's state): replay-safe for an attempt
  sent with `expect_server_instance` equal to the instance that finished it and with `not_after_ms` within
  `tombstone_max_s` of the finish. An unpinned attempt's finished state assumes that no copy reaches a restarted
  process (or another server at the same address) while it could still be admitted there; without `not_after_ms`
  nothing bounds a replay, and a finished state assumes that no copy arrives after the attempt left retention (with
  retention 0, after it finished);
- **on an `expired` 409**: the id's registration is checked before the expiry, under one lock, so no copy was running
  or registered on that `server_instance` and none can start there later. It settles nothing (an earlier copy may
  have run and left retention), and it assumes the server's clock does not step back below `not_after_ms`;
- **otherwise, by the clock**: after `not_after_ms + clock_skew_allowance_ms + bound_ms`. After `not_after_ms` +
  skew an unanswered attempt cannot start subsequently (it may already be running; `bound_ms` covers that, apart
  from the run past the bound of §15.2). The library holds the clock release while the last answer said `running`
  or `stopping`, by default.

Assumed throughout: clocks within the skew, no rewriting of the request's ids, and no copy reaching another server
that serves the same ledger.

## 15.6 The client gone

Every `/traverse` (and `/resolve` between its phases) stops when its client has closed the connection, writing
nothing. A half-closed connection counts as gone. The socket is polled at the walker's checkpoints at most every
100 ms, and the kernel's TCP state is read when bytes are pending: with bytes waiting, every TCP state past
ESTABLISHED reads gone (a content timeout's shutdown leaves Linux sockets in FIN_WAIT2, which reads connected
otherwise); SYN_RECV (TCP Fast Open, not enabled) stays connected; EBADF and ENOTSOCK read gone.

# 16. The index identity as served

## 16.1 Per-graph identity

The multi-graph list has two optional columns, `manifest_path` and `index_ns`; each manifest is checked at start-up
against the files its pair loads (sizes and the digest of the list, as `--index-manifest`), and a mismatch, two
identities for one pair, more than five columns or a bad name refuses to start. Every response states its pair's
identity; two different pairs never state one `index_fp`. `GET /traverse/capabilities?graph=<name>&graph_path=<path>`
describes one pair. A `graph_path` naming one graph listed with several annotations is refused: answering with the
first annotation would silently narrow the question. There is no union oracle over several annotations. A name
listing one pair twice needs no `graph_path`. `scripts/traversal/index_manifest.py --server-csv` writes one manifest
per pair, every distinct file hashed once by parallel streams or taken from precomputed sha256 digests (`--digests`).

## 16.2 The identity inventory

`index_load_inventory` lists every file the server's loaders open for a pair, derived from the listed spelling as the
loaders derive it: the graph; the annotation; the row-diff anchors and fork successors beside the graph for a
`.row_diff` annotation; `.seqs` for a coordinate annotation unless `--no-coord-mapping`. A unit test checks its
annotation table against the types `initialize_annotation` builds. Not in it, because nothing opens them: a column
annotation's `.coords`, the graph's `.weights`, anchors beside another annotation type, the in-memory header index.
The graph's derived files are listed apart (§16.3).

The inventory decides:

- **the manifest check**: a file of the inventory a manifest does not name refuses to start, naming it; a manifest
  may not list an optional inventory file the pair does not load (`.seqs`), since a bundle without its `.seqs`
  beside the symlinks, or a server with `--no-coord-mapping`, would otherwise state the `index_fp` of one that loads
  it; a manifest that lists another graph or annotation is refused (one written for a directory would lend one
  fingerprint to each of its annotations, swapped memberships included). Loaded files are matched to manifest
  entries by base name, so duplicate base names in a manifest are refused, by the server, the CLI and
  `index_manifest.py --verify` (otherwise one manifest listing `A/…` and `B/…` would pass for both bundles).
  `--extra` files outside the inventory stay allowed;
- **the bundle dedup**: two lines are one index only when their whole inventories resolve to the same real paths in
  the same roles, so every distinct bundle is checked against its manifest (two symlinked pairs of one graph and
  annotation with different `.seqs` beside them are two indexes).

`metagraph traverse --index-inventory -i G -a A` prints the inventory as JSON with its table;
`scripts/traversal/index_manifest.py` (single mode, `--server-csv`, `--digests`) mirrors it, and an integration test
compares the two on every sidecar kind. The tool matches digests by real path (relative paths from the digest file's
directory), by base name only with `--match-base-names` and never one line for two files, and writes a pair's
manifest to the path any of its lines names. `index_fp` covers every file the pair loads, by size at start-up and by
digest at build time (a same-size replacement is found only by `--verify`); `index_meta_fp` stays a negative check.

## 16.3 Derived data outside the identity

The graph's dummy-edge mask (`<graph without .dbg>.edgemask`) and its Bloom filter (`<graph without .dbg>.bloom`)
are **derived data** of the graph: computed from it alone, they decide which counts are exact (`/pattern` counts
exactly with the mask and states bounds with an estimate without it, `SPEC-pattern-search.md` §18) and how fast
k-mers are looked up, never what an exact answer is. They are not in the identity inventory (`index_load_inventory`,
`index_bundle_files`); `index_derived_files` lists them apart, and `traverse --index-inventory` prints them under
`derived` with the rule (`derived_rule`). A manifest listing one (`*.edgemask`, `*.bloom`) is refused, naming the
entry and the rule; `index_manifest.py` never writes them, refuses them with `--extra` and flags them with
`--verify`. Consequence: adding a mask to a deployed index (`transform --mask-dummy`) leaves its `index_fp` and
`index_meta_fp` unchanged (shown on a copy of `build/mini_refseq`: the stored graphlets still compare as one index,
and every `/traverse` and `/resolve` answer replayed byte for byte apart from timing). Two spellings of one graph,
one with a mask beside it, are one index with one identity.

The risk is stated (`SPEC-labeled-traversal-core.md`, "The index identity"): no fingerprint covers derived data, so
a stale or foreign derived file beside the graph is not detected by `index_fp`. The loader checks a mask's size and
its W = `$` edges (`mask_invalid` on `/pattern`) and a Bloom filter's k and mode only; a Bloom filter of another
graph with the same k and mode can hide k-mers under an unchanged `index_fp` (it can turn `graph_runs` `[[0, 36]]`
into `[]`). Derived data is written only by `build --mask-dummy`, `transform --mask-dummy` or `transform
--initialize-bloom` on the graph it sits beside. A manifest that lists a mask or a Bloom filter refuses to start and
is written again without them (its `index_fp` changes once). How to close the risk is open (§19).

# 17. The server

## 17.1 Capabilities and feature levels

- **`feature_level`** is what the server offers beyond the base contract, monotonic: a change that adds capabilities
  fields or routes raises it by one, and SPEC §10.3 states what each level adds. Clients key features on
  `feature_level >= n`. A correction that adds no field can stay within a level; SPEC §10.3 then lists every output
  it changes, since a client keyed on `feature_level` cannot tell the builds apart.
- **`schema_version`** is the request schema (`strategy.schema_version`), the one accepted (1, no window); a feature
  level never changes it. A future request change adds `request_schema_versions` and keeps accepting 1 for a stated
  window.
- **`GET /capabilities`** is a small, index-free document: routes and features, `feature_level`,
  `algorithm_version`, `release`, `server_instance`, `mode` and the graph names, the attempts block,
  `deadline_check`, `compression_level`. It answers while a single index loads (`ready: false`).
- **`GET /traverse/capabilities`** (with `?graph=&graph_path=` in multi-graph mode) is the probe of one pair: its
  identity, `algorithm_version`, `compression_level`, `deadline_check`, the `attempts` block (fields, id pattern,
  routes, `server_instance`, retention, `allowance_ms`, `hard_cap_ms`, `clock_skew_allowance_ms`, the client check
  interval, how the bound is computed, the `not_after` rule, the `release_rule` and the `delivery_reserve`),
  `work_bound`, and the `coordinates` block (§18.7). What depends on the index — record coordinates among it — is
  advertised there only.

## 17.2 Delivery

- The traversal routes compress at zlib level 1 (`--traverse-compression-level`; 630 MB/s against level 9's
  153 MB/s on real responses, for 1.8 times the bytes); the other routes use 9.
- The server writes each seed's result as text once built (its tree freed at once) and assembles the response from
  the texts, byte for byte. The assembly reserves the response's exact size (the envelope's members after `results`
  written first) and frees each text once copied; the transport's buffer (the Response's `asio::streambuf`, which
  grows by doubling) is sized once through `rdbuf()` after the last check, so nothing reaches the transport before
  every check passed. On mini_refseq the delivery peak per byte of text is 1.52 (gzip) and 2.57 (identity), against
  2.60 and 3.51 with the texts copied into a growing buffer and kept until the end. Feeding deflate from the texts
  would save little beside the freed texts; sending the body in pieces would commit a 200's header before the last
  check.
- `compress_string` hands the text to zlib in pieces of at most 2³² − 1 bytes (`Z_NO_FLUSH`, the last with
  `Z_FINISH`), because zlib's input length is 32-bit: a compressed response of 4 GiB of text or more is the whole
  text, on every route (`process_request` is the server's one compressor).

## 17.3 The delivery window

Attempts stop walking at the walk-until, `bound − max(allowance/2, reserve)`, so that a response is not lost whole for
want of time. The reserve is

```
reserve = 1.25 × ( (text written so far + walked seed's estimated text) / compress rate
                   + walked seed's estimated text / build rate )
          + stop time
```

- The walked seed's text is estimated from its modelled account: its record coordinates' share at its own bound
  (12 account bytes per text byte, §18.6), the rest at the ratio in use.
- The configured rates (`--traverse-delivery-compress-mbps` 50, `--traverse-delivery-build-mbps` 10) and ratios
  (30 account bytes per text byte for JSON, 50 for a graphlet) apply until the server has measured its own on
  responses and seeds of 1 MiB or more: then the slowest rate and the smallest ratio (per detail) of its last 16
  measurements, and the attempt's own where more conservative. The ratios start just below the smallest measured on
  real responses (33.5 for UHGG `full` to 1,344 for SRA `summary`; 58.4 and more for a graphlet); a server's first
  large attempt of a detail can still be cut early until it has measured that detail (a warm SRA server cut the
  first `tree` and `full` attempts of a 16S beam at 3–5 s of 40 s, the tree measuring 115–129).
- The stop time: the walk does not stop at its walk-until but at the first poll that reads the clock after it, and
  the stopped seed is finalised before its text is built. `chunk_target_ms` + 950 ms is assumed (SRA attempts
  stopped 352 ms after their walk-until on a quiet fresh server, 1,001 ms on a loaded one) until the server measures
  a longer one. Stop latency is measured on the first call per seed only: a seed whose result was being built when
  the client left or a writer refused is not measured again, since that would count its building as stop latency and
  shut out later attempts.
- The margin and the stop time are why a response the model fits exactly is delivered: without them it would be a
  503 about half the time. The reserve can stop a large attempt earlier, and what it did not walk is stated.

## 17.4 Start-up and shutdown

- **The sequence header index** is built while the index loads, for every header-capable server. Built by the first
  request that names a header, it would fall into that request's budget (33 M headers: 4.6 s single-threaded under a
  mutex). An explicit column-only opt-out could save 1–2 GB and start-up work if no header-dependent operation can
  trigger lazy construction; it does not exist.
- **Options** are bounded integers: `--traverse-clock-skew-ms`, `--traverse-chunk-target-ms` and
  `--traverse-attempt-allowance-ms` take integers in [0, 2⁵³ − 1], the retention options their ranges (§15.4); a bad
  value refuses to start, naming the option and its range.
- **Shutdown on SIGTERM** (SPEC §10.3). In a container `server_query` is PID 1, for which the kernel drops a signal
  with the default action, so without a handler every `docker stop` waits its timeout and kills it. SIGTERM and
  SIGINT write a byte to a pipe (the handler's only work, async-signal-safe); a thread reading it stops what is in
  flight and the server: every traversal — each `/traverse` has an `Attempt`, managed or not — is stopped at its next
  poll as for a client gone (nothing is written for it, its connection closed: a partial result would state a cause
  that is not true, "cancelled", and a ledger treats the attempt as unanswered), a `/traverse` or `/resolve` arriving
  meanwhile finds its client gone, the HTTP server stops once they returned or after 3 s, and the process exits as
  soon as it has stopped (`_Exit` from `run_server`, the log flushed), 3 s later at the latest even if a handler (a
  `/search`, an uninterruptible read) has not returned. It does not return through `run_server`: its destructors
  would join what still runs, and in single-index mode the server listens at once and answers 503 while its index
  loads, so a signal can stop it with the load in progress (19 s on SRA in RAM, on a cold disk the whole docker-stop
  timeout). Before the server accepts connections it exits at once; a second signal exits immediately. Tested:
  `TestTraverseAttempts.test_sigterm_stops_the_server_promptly` and `.test_sigterm_while_the_index_loads` (T56).

## 17.5 Deployment

SPEC §10.3: `-p` above the number of concurrent traversals, so that the cancel, state and capabilities routes find
an io thread; a reverse proxy forwards the cancel, state and both capabilities routes with their query strings,
closes its upstream connection when the client goes and keeps a read timeout of at least `hard_cap_ms` plus the
transfer; the server's clock is synchronised within `clock_skew_allowance_ms`. A server memory maximum
(`--traverse-max-memory-mb`) is the budget of every request without a smaller one (§11, §18.8). Before raising a
container's memory or changing read-ahead, measure memory pressure (file and anonymous memory, refaults, reclaim,
pressure stalls, faults, I/O per request); compare `MADV_RANDOM` against normal advice separately, in serial paced
trials that abort on agreed sibling-latency or pressure thresholds.

# 18. Record coordinates

MGT v1 is frozen, so record coordinates are not in the body. They are an **additive JSON field** of the response,
opt-in (`strategy.output.coordinates: true`), carried with the envelope the library already keeps (the `J` line of a
saved `.mgt`) and ignored by readers that do not know it. The size cost is small (about 2–3× the bytes of a body
record for this part before compression, ~1.2× after gzip, and only for trace retrievals; about 441 B a seed on 12
recorded staging trace responses), and the Python side parses JSON faster than MGT text. Reasons that would justify
coordinates in an MGT v2 are collected in §18.4. SPEC §5, §6.8, §7.0, §7.1, §7.5 and §10.3 are normative.

## 18.1 What is reported, and when

Coordinates are reported **only where they are well defined**: under `support: trace` on an index with k-mer
coordinates, where every claim is one coordinate-consecutive occurrence in one record. Under `support: kmer` a
stretch can be stitched from several occurrences, so no single position is honest; nothing is reported, and the
response says so (`coordinates: null` with `coordinates_reason`, §18.5).

For **header labels** (an accession via `CoordToHeader`), coordinates are positions **within the record**, 0-based,
in the record's own orientation; for **column labels**, positions in the column's global coordinate space (and the
`trace_record_boundaries` limitation applies). Trace runs only on basic (single-strand) graphs, so the strand is
always the record's forward strand on the right arm; on the left arm the walk runs backwards along the record, which
the interval below accounts for.

## 18.2 The field

Per seed result (JSON, every detail — `full`, `tree`, `summary` and a graphlet's JSON summary, never its body):

```
"coordinates": {
  "kind": "record" | "column" | "mixed",  // header labels: within the record; column labels: global; mixed: per label
  "k": 31,
  "max_occurrences": 16 | "unlimited",
  "complete": true | false,               // false when a list was cut or a run is a lower bound
  "runs_lower_bound"?: N,
  "seed": [ {"label": <label id>, "occurrences": [[start, end], ...]} , ... ],
  "arms": {
    "left"|"right": [                     // one entry per run, in runs order
      {"run": <run id (the R ordinal / runs[] index)>, "label": <label id>, "from_bp": ..., "to_bp": ...,
       "occurrences": [[start, end], ...],   // one per live coordinate chain of the run, ascending start
       "occurrences_total"?: N,              // when the list was cut
       "chains_ended"?: N,                   // when > 0
       "lower_bound"?: true}
    ]
  }
}
```

- An occurrence `[start, end)` is a **half-open base interval in record (or column) coordinates**. Normative rule:
  the run's first node (the k-mer entered by step `from_bp + 1`) has k-mer coordinate `c₀`, its last node `c₁`. The
  run's **own bases** — the base each of its steps adds, i.e. the walk's bases `[from_bp, to_bp)` — are record bases
  `[c₀ + k − 1, c₁ + k)` on the right arm, where coordinates increase outward, and `[c₁, c₀ + 1)` on the left arm,
  where they decrease outward (a left-arm step adds the first base of the k-mer it enters). The server writes the
  same intervals from the chain's coordinate c at the run's last node: `[c + k − L, c + k)` right, `[c, c + L)`
  left, L = `to_bp − from_bp`. A seed occurrence is `[c_first, c_last + k)` over the seed's first and last k-mer.
- An occurrence is a chain that carries the whole run. Several occurrences per run occur when the label's record
  repeats the walked sequence (two live coordinate chains). `chains_ended` counts the chains continuing on no
  followed path of the lineage (partitions among clones are not counted; clones inherit the count).
- Runs closed by a merge do not occur (trace and merging are mutually exclusive). Runs entered by a switch carry the
  new label's own coordinates from their start — except a run entered by a switch into a label whose own lineage is
  still live there: the walker's entry keeps only that label's chains that continue from before the switch, so
  chains starting at the switch node are missing; such a run carries the inherited chains only and is marked
  `lower_bound: true` (the block's `complete` false, counted in `runs_lower_bound`), the walk unchanged. Marking is
  enough: it arises only with a switch cost (the default is `forbid`), and unbudgeted on mini_refseq's 56 trace
  cells it affects 21 of 3,780 runs (0.56%), all switch-entered (21 of 1,971), in 10 cells, none without a branch
  allowance; the exact chain set is not needed.
- **Column positions** are the column's k-mer index space (record i's k-mer j is `offset_i + j`, `offset_{i+1} =
  offset_i + len_i − k + 1`), so **a record's last k − 1 bases share their numbers with the next record's first
  k − 1 positions**: an interval there is attributed to one record only with the record lengths, which are not
  stated. Mapping column positions to records is not done (§19). Stated in SPEC §7.1, in the probe's
  `coordinates.rule` and exercised by the positional oracle.
- Bounded like everything else: at most `output.max_coordinate_occurrences` (default 16, accepts `"unlimited"`) per
  run; a cut is stated by the list (`occurrences_total`), the block (`complete: false`) and a limitation
  `coordinates` (§18.5) — never silent.

## 18.3 The library

- `Graphlet.coordinates` (from the envelope; `None` with the reason when absent). Body-only graphlets have none and
  say so (`MissingEnvelope` for coordinate queries). An inconsistent block is rejected eagerly, when the graphlet is
  built.
- `Claim.coordinates`: the occurrences of the claimed run(s), **clipped to the claim's displayed interval** (and to a
  cut depth `D`) with the same arithmetic as §18.2; `route_only` claims carry the coordinates of the label's own
  route.
- `walks(...)`: per walk and label, the record intervals covering the walk's bases; `label_walks()` rows likewise.
- `to_fasta()`: headers gain `acc:start-end` (natural orientation, 1-based closed in the header for readability,
  stated in the docs) when coordinates are present. GFA carries no coordinates. `compare()` notes them, never clips
  or compares them.
- The client asks for coordinates for trace strategies when the server supports them (the probe's `coordinates`
  block), and `next_request` carries the setting forward; a cap given without coordinates is stripped (the server
  refuses it with 400). Under a request memory budget it does **not** ask automatically
  (`AUTO_COORDINATES_UNDER_MEMORY_BUDGET` false, the depth gate of §18.8): the answer says so and how to ask.
- MCP tools: `graphlet_claims` rows always carry `coordinates`; `graphlet_walks` / `graphlet_labels` only with
  `coordinates=true` (counted against the byte ceilings).
- **The positional oracle**: on mini_refseq the real-index suite fetches each record's sequence from the fixture
  FASTA and checks that the bases at every reported interval equal the spelled claim (natural orientation, both
  arms) — positional verification on top of the search-based one. `MiniRefSeq.CoordinatesAgainstTheSourceRecords`
  checks 10,696 occurrences of 9,498 runs this way (5,508 of them switch-entered, 48 lower bounds, 326 lists cut,
  and column intervals in a record's last k − 1 bases whose bases are the next record's), every recorded switch
  request's seed walked to 3,000 bp.

## 18.4 Reasons collected for an eventual MGT v2

1. Coordinates in the body (self-contained raw bodies, one conformance regime) — they are in the envelope (§18).
2. An extension mechanism (readers ignore unknown records of a reserved form), so later additions need no version
   bump — v1's strict reader rejects unknown records.

Cut a v2 only when more than one strong reason has accumulated; keep v1 reading.

## 18.5 The request and response contract

- **Request** (`strategy.output`): `coordinates` (boolean, default false; false is byte-identical to omitting it)
  and `max_coordinate_occurrences` (an integer in [1, 2^64 − 2] or `"unlimited"`, default 16). The cap without
  `coordinates: true` is refused (400), because accepting it would suggest a bound that bounds nothing; it is inert
  under `support: kmer` and in annotate mode. Echoed, the cap's default included, only when `coordinates` is true.
- **Response**: the block of §18.2, or `null` with `coordinates_reason` (`index has no coordinates`, `support kmer`,
  `no traversal`, `partial derivation`, checked in that order).
- **Limitation** `coordinates` (seed level, only when a list was cut), knob `output.max_coordinate_occurrences`,
  observed the largest true count among the cut lists, extra field `lists_cut` their number; in MGT exactly one
  `K * coordinates …` record more. It belongs to no outcome class: the block states `complete: false` (§19).
- **Action** `drop_coordinates`, right after `use_graphlet` / `drop_sequences` wherever those are offered, for a
  request with coordinates (block or null form).
- **Memory**: each run's side-table entry is charged once at its creation at min(chains, cap) occurrences, a seed
  label's with the depth-0 state, in the requested detail. `LabelRun`, `Entry` and `Item` do not grow (their sizes
  are part of the memory model, §14.2). **Work**: none beyond the row — every coordinate is charged one unit with
  its row.
- No server cap on the occurrences: the true bound is stated (the probe's `output_bound`, §18.7), and a server
  memory maximum is the operational bound.

## 18.6 The delivery reserve's coordinate share

- The walker's account has a **coordinate share** `ResourceAccount::coordinates`: the output's fixed part for
  coordinates (`DeliveryCosts::coordinate_fixed`, already in `fixed`: the block's skeleton with a cut list's
  limitation and K record, or the null form; counted again only in the share, nothing more is charged) plus every
  entry and occurrence as charged at its creation. It is published with the account at every level's end
  (`AttemptControl::progress(account, coordinates)`) and in the seed's meter.
- **The estimate** of a seed's text (§17.3): E = ⌈(A − C) / ratio⌉ + ⌈C / 12⌉. `kCoordinateAccountPerTextByte` = 12
  is a **bound**, not an estimate: each part of the share is priced at least 12 times the most text it can write,
  whatever the digits (`GraphletCoordinates.CoordinateAccountBoundsItsText`, detail `full` / `graphlet`):

  | part | account (B) | widest text (B) | ratio |
  |---|---|---|---|
  | an occurrence (two 20-digit numbers, brackets, comma) | 872 | 44 | 19.8 |
  | a run's entry (every optional member at its widest) | 3,680 | 221 | 16.7 |
  | a seed label's entry | 1,728 | 79 | 21.9 |
  | the block's skeleton, a cut list's limitation (20-digit count, `"unlimited"`), `drop_coordinates` (JSON detail) | 11,628 | 597 | 19.5 |
  | the same with the K record, the Q token and the counts' digits (graphlet) | 15,151 | 960 | 15.8 |
  | the null form with its longest reason and the tokens (graphlet) | 1,584 | 106 | 14.9 |

  Why a separate share: the measured ratio of the rest of the output (115–129 for a tree) would understate a
  coordinate-heavy seed's text up to 4.1 times (§18.8, M2). Without coordinates C = 0 and E is the estimate of the
  whole account at the ratio, to the bit.
- **The ratio sample** leaves both sides of the share out: (A − C) / (text − the coordinates' text), the latter
  counted **exactly** from digit counts and fixed punctuation (`coordinates_text_bytes`: the block or null form, a
  cut list's limitation, `drop_coordinates`, and in a graphlet the escaped K record, the Q token and the digits they
  add to Z, `graphlet_lines` and `graphlet_bytes`; `compact_json_size` counts compact JSON without writing it), and
  only on a rest of at least `measured_text_bytes`. So an attempt with coordinates measures exactly what the same
  walk without them does, and attempts with coordinates feed the server's measured ratios like any other
  (`GraphletCoordinates.AttemptsMeasureTheSameRatioWithCoordinates`,
  `.CoordinateTextIsExactAndBoundedByItsAccount`, `MiniRefSeq.CoordinateShareIsExact`,
  `GraphletAttempt.CoordinatesLeaveTheServersRatioUnchanged`).

## 18.7 Capabilities

The probe's `coordinates` block: `supported` (= `supports_trace`), `knob`, `cap_knob`, `max_occurrences_default` 16,
the `kinds` this index can report, `limitation`, `action`, `output_bound` — which names the probe's `max_memory_mb`
and `max_work_units` rather than a value, since a server maximum is the budget of every request without one — and
`rule`. `attempts.delivery_reserve` states the coordinate share and the exact sample in its `rule` and
`coordinate_account_per_text_byte: 12` (both GET routes). On mini_refseq the compact probe is 19,490 bytes (the
coordinates block 2,799) and `GET /capabilities` 14,081 — past a 16 KiB tool ceiling, hence the 32 KiB of
`traverse_capabilities` (§6).

## 18.8 When coordinates cost depth: measurements

Coordinates are charged when a run is created, so a walk with them stops shallower under a memory budget. **The depth
gate**: a library may ask for coordinates automatically under a memory budget only if the median depth at the
budget's stop with them stays within 10% of the depth without.

**M1, by regime** (`gate_cell` in `benchmarks/traversal/bench_traversal_measurements.cpp`). Two measures per trace
cell: the depth at the stop (`complete_to_bp` over both arms) under memory budgets of 50% and 75% of the cell's own
opt-out peak, exact bytes, with against without coordinates at cap 16; and, since that ratio drops a cell whose
opt-out walk fails at depth 0 too, what **completing** takes: the smallest whole-MiB budget that holds the walk with
and without coordinates, and what the walk with coordinates does at the budget that completes it without them.
Cells, all to 3,000 bp, details `full` and `graphlet`: mini_refseq (`BM_TraversalCoordinatesDepthAtTheStop`: 7 seeds
— blaNDM both ways, three 200-bp windows of its carriers, two repeat windows — × 8 strategies — branch limits 0 and
2, switch costs 0.5 and 1 at limits 0 and 2 within a loss budget of 2, column labels at limit 2 with and without a
switch); its column cells again under the **refseq33m projection** (every run's and seed label's list priced at the
cap of 16 occurrences, as refseq33m's taxid columns of many genomes give; mini_refseq's taxid columns reach 12 chains
a run at most); the column fixtures (`BM_TraversalCoordinatesDepthByRegime`,
`scripts/traversal/make_column_coord_fixtures.sh`) — `coord_lockstep`, 200 columns of one shared 3,000-bp sequence,
16 records each (one path, 200 runs of 16 chains an arm), and `coord_divcol`, 100 columns of their own sequences
around a shared 600-bp core, 16 records each (seeds in the core: a split into a branch a column at each end of it;
and in a flank: one column) — and the wide fixture (M2 below: column labels one run of 5,000 chains an arm, header
labels 5,000 seed labels, a run of one chain each).

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

- **The gate holds only where the coordinates' share of the account is small: header labels at refseq's chain counts
  and column labels below 16 chains a run** (mini_refseq: median −3 to −4%; every cell no deeper with coordinates,
  nearly every one a little shallower; the outliers are the single-lineage switch cells, limit 0, where two `full`
  cells fail at depth 0 at 50% with coordinates and not without).
- **It does not hold for column labels with 16 or more chains a run — refseq33m's taxid columns — nor for
  header-heavy seeds.** Projected to refseq33m's chain counts the median drop is 13–16% at 50% (27–32% with a
  switch); on the fixtures 13–15% (the core seed) and 34% (one column at 75%); and where many labels walk together
  the seed needs 2–4 times the memory to complete with coordinates and **fails at depth 0** at the budget that
  completes it without them (`coord_lockstep`; the wide fixture's 5,000 header labels). A library that turned
  coordinates on under a memory budget there would turn a complete answer into a refused seed, with
  `drop_coordinates` offered only after the failure. That is why the library does not ask automatically under a
  memory budget (§18.3).

**M2, the wide fixture** (`BM_TraversalCoordinateReserveOnTheWideFixture` on
`scripts/traversal/make_wide_coord_fixture.sh`: one 3,000-bp sequence in 5,000 records under one column, labelled
`wide.fa` — annotated on the relative name, so that the column's label, and with it the column rows' bytes, do not
depend on where the fixture was built —, row-diff anchors every 1,000 rows; a 100-bp seed in its middle walked to
the records' ends; column labels give one run an arm with 5,000 chains, header labels 5,000 runs an arm with one
chain each):

| labels, cap, detail | text (coordinates) B | account (coordinate share) B | build ms with / without | coordinate estimate / their text | whole estimate / text | estimate without the share / text |
|---|---|---|---|---|---|---|
| column, 16, full | 14,608 (1,320) | 1,880,896 (62,572) | 0.15 / 0.24 | 3.95 | 4.51 | 0.94 |
| column, 16, graphlet | 8,990 (1,612) | 1,639,593 (66,095) | 0.10 / 0.16 | 3.42 | 4.11 | 0.86 |
| column, unlimited, full | 291,112 (277,824) | 14,919,040 (13,100,716) | 3.1 / 0.13 | 3.93 | 3.96 | 0.37 |
| column, unlimited, graphlet | 285,202 (277,824) | 14,677,737 (13,104,239) | 3.0 / 0.08 | 3.93 | 3.94 | 0.24 |
| header, 16 or unlimited, full | 5,076,957 (984,563)¹ | 222,416,410 (58,531,628) | 50 / 38 | 4.95 | 2.04 | 1.09 |
| header, 16 or unlimited, graphlet | 1,378,634 (984,563)¹ | 83,242,348 (58,535,151) | 14 / 2.7 | 4.95 | 3.90 | 0.96 |

¹ At cap 16. "Unlimited" writes 9 B more, all in the coordinates' text (the block's `max_occurrences`, `"unlimited"`
in place of `16`); the accounts are the same, since each run holds one chain and no list is cut at either cap.

- The coordinate share's estimate is 3.4–5.0 times the coordinates' real text, the whole estimate at the configured
  ratios (30 / 50) 2.0–4.5 times the real text. Estimated without the share — the whole account at the ratio the rest
  of the output measures — the text of a coordinate-heavy seed would be understated up to 4.1 times (column,
  "unlimited", graphlet: 0.24).
- In every cell the account and the text less the coordinate share are exactly the opt-out walk's. The walk itself
  takes the same time (0.24 s column, 1.3–1.5 s header); building the block costs 3 ms for 10,000 occurrences.
- The reserve stays below the floor (allowance / 2) for these sizes, so the walk-until of a one-seed attempt at 30 s
  is the floor's (35,000 ms) with and without coordinates; the reserve itself grows with them (header, full: 2,551 ms
  against 1,819 ms; column, "unlimited": 1,173 against 1,009 ms).

**Verification against the walk without coordinates**: on the cached real requests with coordinates, stripped of
what they add, mini_refseq unbudgeted 586 of 588 are identical to the opt-out walk (2 differ in timing counters
only), at 300k work units, cap 1 and "unlimited" 588 of 588; at 4 MiB 334 identical and 254 stopping earlier by
memory, never deeper (`drop_coordinates` offered in 336); UHGG 779 identical and 13 timing only; SRA 776 identical
and 27 timing only, beside one beam attempt whose walk-until depends on the load.

# 19. Open questions

- **Coordinates under a memory budget.** The library does not ask for them automatically under a memory budget,
  because the depth gate fails for column labels with 16 or more chains a run and for header-heavy seeds (§18.8). An
  option: ask automatically only for header labels with few seed labels, never for column labels (the library's
  constant and its rule; a request can always ask).
- **The majority rule under a memory budget** (§11). `tree` and `full` detail keep the arrival order under a memory
  budget so that no stop moves; on a server with a memory maximum every `tree` and `full` response is therefore
  budgeted and only `graphlet` and `summary` follow the rule. The alternative, the rule everywhere, moves stops
  shallower in those details (52 of 6,552 cached budgeted responses; 1–4 levels on an arm, twice 30) and is the one
  condition `majority_first` in `Walker::merge_level`; charging the longest parent's chain at every merge moves more.
- **The annotate rule** (§11) compares the parents' own segments, so after nested merges it can display a walk
  carried by fewer labels; a rule over equal stretches of the displayed walks (the last W bases through each parent
  and its first parents, labels present at every position) ranks those right but changes the library's check
  (`derive.carried_labels`) with it.
- **`setup` pieces** (§14.4) are not counted in `observed_max_uninterruptible_ms`; whether they should be.
- **An outcome class for the `coordinates` limitation** (§18.5): it is in none; the block states `complete: false`.
- **Column coordinates mapped to records** (§18.2): positions stay in the column's space; mapping them needs the
  record lengths.
- **Derived data** (§16.3): a spot check at load (a fixed-seed sample of the graph's k-mers looked up through the
  Bloom filter), or the graph's size and digest stored in the derived files' headers, would close the risk.
- **Not measured yet**: coordinates on the refseq33m panel on staging (`bench_traverse.py --coordinates`); the
  start-up cost of building the header index before readiness on its own.
