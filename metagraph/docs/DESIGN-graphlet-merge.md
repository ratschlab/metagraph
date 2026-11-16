# Design: combining graphlets — the cross-pair union and stitching a continuation

**Owner-review edition, 2026-11-15, revised after the full Fable and GPT reviews.** Implemented since in train 7:
phase 1a-U (`union.py`, the source-supported mode), phase 1a-S (`_stitch.py`, experimental, every verdict
`qualified`) and a `StitchedGraphlet` as a union part (§2.3); the rest is design. Nothing is final. The owner asked to "write this design down first and I will review it end to end", and said of his
answers so far: "do not take my answers as the final answers yet". Every decision that rests on one of them is marked
provisional, and §8 lists them all.

The design went through eight revisions, four Fable reviews, an external review and a GPT review. Every finding was
checked against the code, reproduced where that was cheap, and settled; Appendix G maps each one to its verdict (the
external review in G.4, the full Fable review of `8e0083759` in G.5, the GPT review of `8486498fe` in G.6). The owner's
scope since the external review: two path modes for the union (§2.7) and the combination of two retrievals of one pair
as a later operation (§4.2). The body (§0–§8) reads end to end for the decision, and the appendices hold the exact
refusal texts, the data classes, the long tables, the later phases in detail, the full test lists and the reviews. Items
marked **[unsure]** are my reading and need a check.

Base: train 7, lane G1, `bc9684559` (train 6's head with W5's `93044c6ee`). Paths are relative to `metagraph/`;
library paths are under `api/python/metagraph/traverse/`. Line numbers are `bc9684559`'s.

## 0. The decision in brief

1. **Two combinations are proposed now, and a third is designed for later.**
   - **The cross-pair union**, `union_pairs(results)` (§2, phase 1a-U). It takes everything one seed's request
     returned across several pairs, such as r237's 11 graphs or two databases, and combines the graphlets by
     sequence: one trie per arm of the parts' spellings, with the pairs and labels that carry each segment. It is not
     exact and does not claim to be: it is "the best we can do today" (the owner's words, OQ3, provisional). Its
     paths come in two modes (§2.7): **source-supported**, the default, where every path is held end to end by one
     graph, and **inferred**, on request, where a path may switch graphs at a connection the input graphlets justify and
     names every switch. The first release builds the source-supported mode; the inferred mode is designed here and
     built in a later phase (OQ18, decided). Five statements come with every summary, printout and export, and
     `statements()` gives them for any arm and depth (§2.1): coverage per pair, provenance, what a turn-off means,
     labels per pair, and the scope of the paths.
   - **Stitching a continuation**, `stitch(parent, child)` (§3, phase 1a-S). It joins a continuation to the walk it
     continues, on one pair, and states at each junction and walk why the joined walks may differ from one
     uninterrupted walk (`qualified`, with named causes). The first release states only `qualified`, and its
     completeness and evidence carry every cause it finds (OQ19, decided). `exact`, on a narrow whitelist that
     compares every decision the child made (§3.4), is designed now and switched on later, once a seeded randomized
     comparison against uninterrupted traversals passes.
   - **Two retrievals of one pair** (§4.2): designed now, built later.
2. **The union ships first** (§7.1). It has a named user and a need today, it states nothing beyond its parts, so
   no part of it waits for an oracle, and it builds the shared machinery that stitching then reuses.
3. **Neither result is a graphlet.** Both hold their parts' bodies unchanged and compute every answer from them (§4.4).
   No answer reads more complete than the weakest part, junction or walk-scope cause it depends on (§2.6, §3.5).
4. **Not proposed:** combinations of different seeds, and an MGT body for any combination. The union of the two
   arms of one seed becomes a later validation helper (Appendix D.5).
5. **No engine change and no MCP tool.** The service keeps one record per graph and seed; a union is built from
   those records, not stored in their place (§5).
6. **What the owner decides:** whether to confirm the provisional answers (§8.1; OQ18, OQ19 and OQ20 are decided),
   whether two extras cut from 1a-U come back (OQ12, §8.2), and three questions the reviews raised (OQ15–OQ17, §8.3).
   The service session has answered the service's precondition for r237 (OQ21).

## 1. Use cases

### 1.1 The owner's case: one seed across r237's 11 graphs

A seed fetched on r237 comes back as 11 records, one per graph, and each graph is a taxonomic group. The question
is the release's: which graphs carry this flank, how far each was walked, which records in them carry it, and which
graphs did not hold the seed, failed or were not asked. Today the user reads the 11 records one by one (§1.3). The
cross-pair union answers it in one object; §2.2 shows what that looks like on r237's shape.

### 1.2 All use cases, and whether each needs a new operation

Sources: graphlets.rst, `DESIGN-traverse-graphlet.md` §5–§6, the library's tool functions (`mcp_tools.py`), the
service's tool descriptions (`traverse_fetch`, `traverse_results` and the seven reading tools `graphlet_summary`,
`graphlet_labels`, `graphlet_walk`, `graphlet_walks`, `graphlet_sequence`, `graphlet_compare`, `graphlet_export`), and
the owner's provisional answers of 2026-11-15.

| # | the user has | the user wants | where it arises | needs a new operation? |
|---|---|---|---|---|
| U1 | a retrieval to radius R, and the continuation of walk w (child) | one view of the long locus: one walk from the seed to R + E, support changes across the junction, one FASTA record, claims that do not stop at R | the intended loop "re-seed from a walk's continuation" (graphlets.rst); `next_request`, `deepen`, `traverse_continue` | **yes: stitching** (offered, §3) |
| U2 | one seed fetched with `direction: left`, then with `direction: right` | one graphlet with both arms; or a check that the two agree with a `direction: both` fetch | a budget spent on one arm first; validating the engine and the library | **no** for the answer: a fetch with `direction: both` serves it. **Later: a validation helper** (Appendix D.5) |
| U3 | one seed and one index retrieved twice (tuned and exhaustive, constrain and annotate, two radii, two `max_labels_per_node`) | which answer holds, and what one has that the other lacks | graphlets.rst "Comparing retrievals" | **no** today: `compare()` answers it. **Later: the same-pair combination** (§4.2), which holds both with their provenance but answers neither request |
| U4 | one seed constrained on label set A, then on set B | the walks of A ∪ B | label batches; a service request per labelled seed over several graphs | **no** today: one request naming both sets answers it, with its own certificate. **Later**: the same-pair combination (§4.2) |
| U5 | one seed on several pairs of one release (r237: 11 graphs, so 11 records), or in two databases | which graphs, and which records in them, carry this flank, and how far each graph was walked; one trie of the flanks; one FASTA or GFA | `traverse_fetch` over a database with named graphs ("a record per graph and seed, none merged"); the owner, today | **yes: the cross-pair union** (offered, §2, phase 1a-U): the best combination available today, with every limit stated; source-supported paths by default, inferred paths on request (§2.7) |
| U6 | two seeds of one locus (a gene tiled with seeds) | a map of the locus | the service's `graphlet_compare` description ("two seeds of one locus") | **not a combination of graphlets**: positions and certificates belong to each seed (§4.3); left out |

U5 and U1 need new operations. U3 and U4 are served today and get a later same-pair combination, U2 is served and
gets a later validation helper, and U6 does not combine graphlets.

### 1.3 What the library does today

A graphlet is one retrieval: one seed, one pair, one normalized strategy, up to two arms, a certificate per arm
(`complete_to_bp`, `completeness_scope`) and per seed (`outcome`, `limitations`, `resource_stop`). Every operation
reads one graphlet or compares two: `subgraph()` is a filtered view, and `compare()` compares two retrievals and
refuses nothing. A continuation (`next_request()`, `deepen()`, `traverse_continue`) is a new traversal whose result
is not joined to its parent; the tool functions store only a link. A `/traverse` request runs on exactly one pair,
so a seed sent to r237's 11 graphs gets 11 results; a graph that lacks the seed answers without a graphlet, and the
service keeps one record per graph and seed, none merged. graphlets.rst says "Merging graphlets is not supported."
Appendix F.4 gives the detail, with code references.

### 1.4 Is each combination worth carrying?

The questions follow the checklist in `tools/writing-rules/guides/STEWARDSHIP.md`.

- **The union.** The owner needs it today. Without it, a user joins the 11 records by hand and loses exactly what the
  five statements keep: which graphs were walked less far, which did not hold the seed or were not asked, that a graph
  not following a walk may still hold it, that equal label names in two graphs are not shown to be one record, and that
  a dead end in one graph is not a dead end in the release. A documented recipe is the cheapest alternative, and it is
  rejected for that reason. The smallest version is phase 1a-U after the cuts of §7.1, with no new file format: parts
  are saved as `.mgt` files today, and the union is pure, so building it again after `load()` is cheap. The library's
  maintainers keep it; the one-part property test (P2) holds it in step with `ops.py`. Removal deletes one module and
  one docs section.
- **Stitching.** Python users of the documented deepen loop (`TraverseClient.deepen()`, `next_request()`) and agents
  using `traverse_continue` would use it; no user request is on record. Without it, a user joins parent and child by
  hand (place the child at outward position L, §3.2, and drop the n-base overlap) and loses which labels the
  continuation did not carry, which knobs differ, which merges with other walks it could not see, and where the child
  walked edges the parent's walk had used. A recipe is the cheapest alternative; the owner's provisional answer is to
  build phase 1a (OQ1). The smallest version is phase 1a-S; claims across junctions (1b) and a bundle file (1c) come
  only if users need them. It ships with `qualified` verdicts only, and its API is marked experimental until `exact` is
  switched on (§7.1).

## 2. The cross-pair union (phase 1a-U)

A `/traverse` request runs on exactly one pair (`SPEC-multigraph-envelope.md` §1.1), labels are index-local (SPEC
§7.5.7: joined across retrievals "only when the `index_fp` agree"), and two different pairs never state one
`index_fp` (SPEC §10.3). None of this is a reason to refuse a combination across pairs. It is what the union has to
state, and it states it with every answer (§2.1):

- **Walks compare by their bases**, so different k, regime, release or database do not matter to the trie.
  Different alphabets do (§2.3).
- **Every certificate stays its pair's.** A pair is complete to a depth only as its own graphlet says, and the union
  is never more complete than its parts (§2.6).
- **A pair not following a walk did not follow it under its own request**: its mode, its permitted labels and its
  caps decide that, not only its graph (§2.4).
- **Labels are each pair's own**, shown with their pair (§2.5).
- **By default every path is held end to end by one graph** (§2.7). On r237, where each graph is a taxonomic group,
  a path across graphs would join two organisms. Paths that switch graphs are given only on request, as inferred
  paths that name every switch.

An earlier version of this note refused the combination and pointed to an engine feature ("the graph-less traverse
of a later release", envelope spec §1.1). The owner's answer (OQ3, provisional) is to offer it: "It's the only way I
can do it today", and "it won't be exact, I know, but the best we can do today, I think".

### 2.1 What it is, and the five statements

`union_pairs(results)` takes everything one seed's request returned across several pairs: the graphlets, and the
pairs that answered without one. It returns a `PairUnion`: per arm, one trie of every part's spellings (§2.4),
with the pairs and labels that carry each segment. The parts stay unchanged, and every answer is computed from them.

It is not exact, and it does not pretend to be. It is the best combination available today: everything one seed's
retrievals across a release say, in one place, with the limits next to the answers. Five statements make it
honest. They come with every summary and export: as sentences in `summary()`, `print(union)`, `coverage()`'s sentence,
`evidence()` and the GFA comment lines, and as tokens in the FASTA headers (§2.9). List queries (`walks()`,
`carriers()`, `segments()`, `labels()`) and `spell()` return plain lists and strings; `statements()` gives the five
for any arm and depth, structured (a `Statements` object, Appendix C.1):

1. **Coverage per pair** (§2.6): at the depth asked, which pairs are complete, which are complete with a named
   qualification, which are certified only to a smaller depth and why, which do not hold the seed, which were not
   walked or failed, and which were expected but not given.
2. **Provenance per segment and walk** (§2.4): which pairs spell each segment, how far each pair holds each union
   walk (as a displayed walk, or through a merge into its other walks), where it ends or turns off, and which of its
   labels carry the walk from where.
3. **What a turn-off means** (§2.4): a pair not following a walk did not follow it under its own request. Only a
   span marked `absence: graph` says the pair's graph lacks the walk.
4. **Labels per pair** (§2.5): every label is shown with its pair (`g1/NC_000913.3`), equal names in two pairs are not
   claimed to be one record, and header names on a pair whose record mapping was not checked are marked so.
5. **The scope of the paths** (§2.7). Source-supported, the default: "A walk that no single graph holds entirely is in
   no part, so not in this union." Inferred: "Joined paths connect graphs at the transitions they name, and no graph is
   shown to hold them whole; routes are held whole by their graph through its own merges; no inferred path is shown to
   be contained in a sample."

The statements lead every summary and are never cut (§2.8, `summary()`).

### 2.2 A worked example on r237's shape

**Illustrative.** No real r237 data was used: the graph names, lengths, counts and outcomes below are invented to
show one of each case on r237's shape (11 graphs, k 31, one seed of 40 bp, so 10 k-mers; constrain mode with
derived seed labels, `kmer` support, the default `on_reconverge: merge`, both arms, radius 1000 bp). The outputs
follow §2.6 and §2.9; the one-line-per-span listing is how `carriers()`'s result prints.

The request went to all 11 graphs. Six answered with a graphlet (g01, g03, g04, g06, g07, g08), three answered
without one (g02 and g09 do not hold the seed; g05 holds it, but no label carries all of it), g10's answer failed
with a deadline, and g11's result was not passed in.

```python
>>> u = union_pairs(results, expected_pairs=r237_graphs, seqs_order=r237_seqs_order)   # the 10 answers, any order
>>> print(u)
union of one seed (40 bp) over 11 expected pairs: 6 parts, 4 absent records (1 failed), 1 not included
right arm, to 1000 bp (constrain, kmer, in every pair): complete in g01, g04, g08; complete, qualified, in g03
  (complete for 1000 of 4312 seed labels), g06 (a merge united histories: complete for the united-history rule,
  not per path); incomplete in g07 (certified to 420 bp: max_steps); seed not in graph: g02 (0 of 10 k-mers),
  g09 (4 of 10 k-mers); not walked: g05 (seed in graph; derivation: no label carries the whole seed); failed:
  g10 (deadline exceeded); not included: g11. Seed of the absent records not checked: g02, g05, g09
  (unchecked).
left arm, to 1000 bp (constrain, kmer, in every pair): ...
walks: right 14 (3 held by two or more pairs; 2 stub walks not listed), left 6
labels: each pair's own, shown as <graph>/<name>; equal names in two pairs are not claimed to be one record
  (1 name occurs in two pairs)
A pair that does not follow a walk did not follow it under its own request.
A walk that no single graph holds entirely is in no part, so not in this union.

>>> c = u.coverage('right', at_bp=500)
>>> c.by_category
{'complete': ['g01', 'g04', 'g08'], 'complete_qualified': ['g03', 'g06'], 'incomplete': ['g07'],
 'seed_not_in_graph': ['g02', 'g09'], 'not_walked': ['g05'], 'failed': ['g10'], 'not_included': ['g11']}
>>> c.every_pair_certified_to_bp
None                       # g03 and g06 are qualified, and g05, g07, g10, g11 are not complete

>>> print(u.carriers('right', 3))               # walk 3: an 812 bp flank; one line per pair
g01  whole      to 812  labels 3, 2 of them from 0 bp           complete
g04  whole      to 812  labels 1, from 0 bp                     complete
g06  merges     at 85 into s31; route held to 812, goes on       labels 2: 1 from 0 bp (continues), 1 from 85 bp
g07  ends       at 420  max_steps                               absence: uncertified
g03  diverges   at 140                                          absence: under_strategy
g08  diverges   at 40                                           absence: under_strategy
header names of every pair rest on an unchecked record mapping (count_only)
```

And the FASTA record of the same walk starts with this header (one line in the file):

    >right.u3 kind=walk len=812 held=g01,g04,g06 whole=g01,g04 merges=g06@85>s31 ends=g07@420:max_steps/uncertified
     diverges=g03@140/under_strategy,g08@40/under_strategy certified=g01:812,g03:140,g04:812,g06:812,g07:420,g08:40
     labels=g01:3/2,g04:1/1,g06:2/1 modes=constrain absent=g02:seed_not_in_graph,g05:not_walked,g09:seed_not_in_graph,
     g10:failed not_included=g11 mapping=count_only scope=within_one_graph

In the later inferred mode (§2.7, OQ18), a request adds the paths that switch graphs where the parts allow it. One of
them:

```python
>>> p = u.walks('right', paths='inferred')[0]      # at most one switch (max_transitions=1)
>>> p.id, p.kind, p.length_bp
('i:4.31.52', 'joined', 860)                       # its union segments; one transition
>>> p.pieces                                       # (from_bp, to_bp, the pairs holding that stretch)
((0, 55, ('g06',)), (55, 700, ('g01', 'g06')), (700, 860, ('g01',)))
>>> p.transitions
(Transition(kind='merge', by=('g06',), at='s31', at_bp=85, from_bp=55, to_bp=700, before=('g06',), after=('g01',)),)
```

What the reader learns about this flank: g01, g04 and g06 hold all 812 bp. g06 holds it through an allele its
display does not show, so it is held but not displayed whole. g07 was cut at 420 bp by `max_steps`, so nothing is
known beyond that. g03 and g08 turned off under their own requests, and their graphs may still hold the flank. g02
and g09 lack the seed, g05 holds the seed but walked nothing, g10 failed, and g11 was not passed in. No flank that
runs from one graph into another is shown, and none is claimed absent. The inferred path runs along g06's other
allele, joins g06's displayed path at g06's merge at 85 bp, and leaves it at 700 bp along a branch only g01 carries.
g06 holds it to 700 bp and g01 from 55 bp (where its allele and g06's other one have the same bases), so the two
overlap on [55, 700) and neither holds it whole: it is a candidate, and each piece keeps its own graph's certificate.

### 2.3 Inputs, preconditions and refusals

**One list of results.** `union_pairs(results, *, seed=None, expected_pairs=None, seqs_order=None, budget=None)` takes
the results as they come back and sorts them into parts and absent records. The sorting is total:

- a `Graphlet` with its envelope, or a `TraverseResponse` holding one seed result (a part, or an absent record when the
  result has no graphlet);
- `{'response', 'index'}`: seed result `index` of a response with several seeds, its pair still stated by the
  envelope;
- `{'response', 'sequence'}`: a response with the caller's statement of the seed;
- `{'pair', 'graphlet'}`: a part on a pair the caller states;
- `{'pair', 'result'}` (optionally with `'error'` and `'sequence'`): an absent record from a seed result without a
  graphlet;
- `{'pair', 'error'}`: a pair whose answer was refused or never came;
- anything else is refused. Appendix A.1 gives the full decision table.

**Pairs and tags.** A pair is `{graph, graph_path, annotation_path}` (with `release` and `index_ns` where the caller
knows them) or a graph name. One rule decides whether two descriptions name one pair: fingerprints decide where both
sides state them, otherwise names, paths, release and namespace must agree wherever both state them. A description
that matches two pairs is refused `ambiguous_pair`. Every pair has a display tag, unique in the union, by which every
answer, selector, export, evidence block and provenance entry names it; caller-stated pairs are named as such in the
coverage sentence. Appendix A.4 gives the rule, the tags and the pair key.

**One retrieval per pair.** The same retrieval given twice (the same body and the same retrieval identity: the index
identity, the normalized strategy and the seed as requested) is a no-op, so a union can be rebuilt from overlapping
sets of records. A failed record followed by an answer for that pair is replaced by the answer, and provenance keeps
the failure. Two different retrievals of one pair, one body from two requests included, are refused `duplicate_pair`:
equal bodies are not equal retrievals, since two requests that differ in a knob nothing reached give byte-identical
bodies but continue differently. Combining two retrievals of one pair is the later operation of §4.2. Appendix A.4
defines retrieval identity.

**The seed of an absent record.** An absent result's `seed` object carries the request's `seed_id` (the caller's
text, echoed), `length_bp` and `labels_from_seed`, and nothing else: no sequence, no `validated_seed_id`, no labels
(`failed_seed_head`, `src/cli/traverse.cpp:4859`; the answer comes at SPEC §6.1 step 2, before labels are
validated). So each absent record states `seed_check`:

- `sent`: the record came in a `TraverseResponse` from the library's own client, which now keeps the seeds it sent
  (a client change in 1a-U, §7.1), and the sent sequence equals the union's seed;
- `stated`: the caller gave the seed's `sequence` with the record, and it equals the union's seed. This checks the
  caller's statement, not the server's record;
- `unchecked`: otherwise. A `length_bp` that differs from the union's seed is refused (`different_seed`).

`verified` is reserved for an id the server validates over the sequence in an absent result, which no server sends
today. Sequences are compared after the server's case mapping (upper case on DNA builds). The coverage sentence names
each absent record's `seed_check`.

**Required:** whole retrievals with their envelopes and their bases; one alphabet; one oriented seed sequence, base
for base, taken from the parts; `seed=`, when given, equal to it (`different_seed`), and required when no result
carries a graphlet ("the seed is in none of the 11 graphs" is a coverage answer); one retrieval per pair.

**Allowed, and stated per pair:** different k, regime, release, server or database; different strategies, including
mode (constrain, annotate), support (`kmer`, `trace`), radius, caps and `on_reconverge`; a part that walked one arm
only; a part whose identity is unverified (no `index_fp`) or whose pair the caller stated. Each pair's mode is
stated in the coverage sentence and on each span wherever a part is not in annotate mode (§2.4).

**Refusals.** The arguments are checked first (`expected_pairs` must be a list or tuple, not a string; `seqs_order` a
dict of pairs to `verified`, `count_only` or `None`; `seed` a string). Then each check of the table runs over all
elements in the order given before the next check runs (check-major order), and the first failing check decides; a
refusal names the first failing element. Commutativity (§2.10) covers successful unions only. `%s` names an element by
its tag or, before a tag exists, by its position. The last three codes refuse queries, not construction.

The codes, in that order: `bad_expected`, `bad_seqs_order`, `not_combinable`, `view_unsupported`, `not_one_seed`,
`result_has_graphlet`, `no_graphlet`, `no_seed`, MissingEnvelope, `pair_mismatch`, `no_pair`, `ambiguous_pair`,
`no_bases`, `different_alphabet`, `different_seed`, `duplicate_pair`, `not_expected`, `ambiguous_name`,
`summary_too_small`, `not_offered`. Their exact texts are in Appendix A.2.

- **A reverse-complement part is refused**, not flipped. On primary graphs a seed and its reverse complement are one
  seed; flipping a part's walks (each arm becomes the other, reverse complemented) is phase 2.
- **A `StitchedGraphlet`** (built in train 7, experimental) is one pair's answer: its root and continuations on that
  pair. Its walks are its stitched walks, and a span's `certified_to_bp` is the stitched one, with the junction
  verdicts named. The relation goes this way only: a union of stitched parts, never a stitch of a union, since a
  continuation runs on one pair. As built (the coordinator's decisions of 2026-11-16):
  - *Insertion.* The root goes into the trie as a plain part (its stubs, merges and routes unchanged), and each
    child's first-parent tree goes in at the trie node where the walk it continues ends (depth `at_bp`, positions
    shifted by the child's `offset_bp`); the continued leaf ends no union walk. So a union with `stitch(g)` answers
    every query byte for byte as with `g`, except the part's digest and retrieval identity (its `digest`), and the
    union's digest. Inserting walk by walk would drop the root's stubs and merges.
  - *Past a junction.* Every span read in a child (a span whose deepest spelling of the walk is a child's, or a route
    that crosses into one) names the junctions it reads (`PairSpan.junctions`, C.1), has `complete` false, and is
    certified to `min(certified_to_bp, the first junction's certified depth)`, which is the stitched walk's
    `certified_to_bp` for a walk the part holds whole. Its absence is `graph` at a dead end and `uncertified`
    otherwise, never `walk_rule` or `under_strategy`: the child's own seed re-entry or edge history may stop it where
    an uninterrupted walk goes on. A dead end where the child recorded a blocked base or an unfollowed hairpin for
    the next base (annotate mode records a skipped hairpin so, with path reason `D`) is not a dead end of the graph:
    `uncertified` there too. Each `SpanJunction` carries the junction's own cause codes, then the walk-scope causes
    (B.2, W8) of the stitched walk's pieces the span reads (to its `to_bp` for an `extends` span), in rank order;
    a `merges` or `diverges` span in a child names no walk there and adds none.
  - *A route that reaches a continued walk.* A route through a merge of the root (or of any part) that reaches the
    end of a walk a child continues goes on in that child's tree, from its root segment, along the union walk's
    remaining bases; it reads the junction, so the span is neither complete nor certified past it (the
    coordinator's decision of 2026-11-16). Stopping the route there would state that the pair's walk ends where the
    stitched graphlet goes on. The continued walk may end at a leaf with bases or at a zero-length merge leaf that
    the route reaches by its last step or by running out of bases at its parent; both cross. The tests build the
    case from `splits45` by making the reference's first parent at its merge the stub (no committed fixture has it),
    on the right arm and mirrored onto the left, and from a synthetic zero-length leaf on `limits`' model.
  - *Labels.* Not joined across a junction (D.1, phase 1b): a child's span labels start at the junction, and at a
    depth no child base reaches they are the parent's along the continued walk. The union therefore states less than
    `StitchedWalk.labels_full`, and `walks(labels=)`, which keeps the walks a label carries whole, never selects a
    walk with bases past a junction.
  - *Statements.* The pair's coverage category is its root's, with a note naming its junctions on the arm and the
    walk causes below them (`StitchedGraphlet._causes_below`); its evidence block is
    `StitchedGraphlet.evidence(arm)` on an arm with a junction (exact false), its root's otherwise, `exact` then
    judged on the root's identity statement alone; `summary()` lists its junctions with their `walk_causes`, as
    `StitchedGraphlet.summary()` does, `provenance()` each child's entry (`stitched`), FASTA the `junctions=` token
    with the span's cause codes (A.5) and GFA a comment line per stitched pair with the junction causes and the walk
    causes below.
  - *Order.* Two forms of one stitched graphlet with one digest (children that differ in a seed summary, say
    `seed_id`) are one retrieval; the form kept ranks on the canonical JSON of every sub-part's envelope and seed
    summary (A.4), so the answers do not depend on the order given (§2.10).
  - *Identity.* A stitched part beside a plain graphlet of its pair is `duplicate_pair` (`stitch(g)` and `g`
    included); the same stitched graphlet twice is a no-op; every part of it must carry bases (`no_bases`, the
    text naming the pair and its first part without them: `part A (its stitched part right:0) carries no bases`).

### 2.4 The union trie: spellings, segments, walks

**Spellings.** In each part, every segment's first parent makes a tree over the arm's segments. Its leaves are the
walks' leaves (a walk is the first-parent chain to its leaf, `derive.py:258`) and the ends of **stub** chains. A
merge consumes the heads of all its parents (`walker.cpp:4765-4772`), and the parser refuses a segment that ends in
both a split and a merge (`parser.py:694-697`), so a non-first parent of a merge has no other child: it is a leaf of
the first-parent tree, lies on no walk (`_first_paths`, `ops.py:641`), and ends a stub. Every segment of a part
lies on a walk or a stub, so the union holds every base of every part. The parser does not refuse a
segment that is a non-first parent of two merges, so `union.py` asserts one merge per stub end (a
`GraphletFormatError` naming the part and the segment otherwise). Positions are outward from the seed's boundary in
the seed's orientation, so spellings compare across k, regime and radius. In one graph a prefix from one seed names
at most one path (a de Bruijn path is fixed by its start and its bases), so a node of the union trie, a spelled
prefix, is at most one position on one path in each pair.

**Construction.** One depth-first pass over each part's first-parent tree inserts it into a radix trie, so the work is
linear in the part's bases. Inserting walk by walk would be quadratic on a comb, as `walk_rows()` warns (`ops.py:1473`).
Trie edges are compressed into union segments, written `s<id>`. A segment ends where the trie branches, where the set of
pairs following it changes, and at every part's merge, so that each merge is at a trie node: the start of a union
segment when the merge segment has bases, or the end of one when it has none. A merge segment can have none: a merge at
the radius is a zero-length leaf (`limits.mgt` line 21, `G 1,2 7 0`), and a merge node that splits at once is
zero-length with two children. A trie node is written `s31` (the start of union segment 31) or `s30$` (the end of union
segment 30). Children are ordered by their first base in the alphabet's order. A segment's bases are a reference into
the first part (in canonical order) that spells them, never copied.

Per segment and pair, a `SegmentSupport`: the part segments it covers with their offsets, the number of walks and of
stubs that pass it (the walks themselves are derived on demand, so the support stays linear on a comb), and the labels
alive there, as intervals: `LabelInterval(label, from_bp, to_bp)` within the segment. They come from the part segment's
own alive sets (`G.entry`, `G.end`) and its runs (`R`), not from the displayed profile, which is anchored at a walk's
leaf and misses labels alive at a cut depth (`derive.py:804-805`). A support that changes inside a segment keeps its
order and extent: in the `reentry` fixture one segment carries A, then B, then A again. A segment belongs to one path,
so its alive sets do not depend on which walk passes it; the intervals are built once per part segment from its entry
set and its label changes, in the order `alive_at()` applies them (`derive.label_changes()`, `derive.py:72`;
`ops.py:437`): a label's interval opens at the segment's start or at its switch-in and closes at its run's end, so the
work is the entry set plus the changes, with `labels_cut` where the part's lists are cut (`label_lists`,
`inexact_counts`). Annotate mode records one set per node.

**Walks and spans.** A union walk ends wherever some part's displayed walk or stub ends: the trie keeps a terminal
marker there, so a union walk need not be a leaf. Otherwise a pair whose walk is shorter
than another pair's would lose it: with A answering `AC` and B `ACGT`, a leaf-only union shows only `ACGT`, and
`walks(pairs=['A'])` returns nothing. A union walk's `kind` is `walk` when some pair's displayed walk ends exactly
there, and `stub` when only stubs end there. `walks()` and the exports leave stub walks out unless `stubs=True`,
and `summary()` counts them apart. For each pair a union walk has a `PairSpan`: how far that pair's spellings follow
it (`to_bp`) and how they stop:

| `how` | meaning |
|---|---|
| `whole` | the pair has a displayed walk spelled exactly so; its end reasons and `complete` are that walk's |
| `extends` | the union walk is a proper prefix of a displayed walk or a stub of the pair: the pair holds it and goes on, so its end is no end in that pair (no end reason; `certified_to_bp` is the smaller of the length and the part's `complete_to_bp`) |
| `ends` | a displayed walk of the pair ends at `to_bp`, before the union walk does, with its end reasons. A walk of the pair that ends at the seed (`to_bp` 0, a root leaf) is `ends` too, with its reason |
| `merges` | the pair's stub ends at `to_bp`, where its path enters its merge at trie node `into`; the span's `route` says how far the pair holds the rest (below) |
| `diverges` | the pair's walks and stubs go on at `to_bp` with other bases (another union walk) |
| `none` | the pair did not walk this arm (`to_bp` 0) |

**Routes through merges.** The part holds more than its displayed walks: past a merge, the bases below `into`
continue every parent's path. Under `on_reconverge: merge`, the constrain default, a graph that holds a walk through a
non-displayed allele would otherwise read as stopping at the merge, and "how far does each graph carry this flank" would
understate it. So a `merges` span carries a `route`: from `into`, the union walk's remaining bases are matched through
the part's segments, following children through any parent and crossing further merges (at most one child matches each
next base, since a prefix names one path per graph). The route states `to_bp`, how it stops (`held`: it holds the union
walk to its end, with `goes_on` when the part continues past that end; `ends`, at a part leaf, with its end reasons;
`diverges`), its `absence` where it stops short, and the merges it crossed. The work is O(the walk's remaining length)
per span; `carriers()` computes it on demand, and `walks(by='pairs')` and `pairs=` compute it for the walks they read,
charged to `budget=`. `whole` stays reserved for displayed walks spelled exactly so; `held_in` lists every pair that
holds the whole union walk, displayed (`whole`, `extends`) or through a route, and it is what ranking, filters and the
exports use. A route says the pair's graph holds these bases as one path its traversal explored: under `merge` the
history below a merge is the united one, which the coverage already states (`complete_qualified`, `scope`).

**What a turn-off means.** In constrain mode a walk is followed only while some permitted label supports
every node of it (SPEC §6.10), and each pair derives its own permitted set (SPEC §6.1 step 4). Beyond a pair's
certificate nothing is certified: "a capped or pruned result is never evidence of absence" (SPEC §6.7). So every
`ends`, `diverges` and `none` span, and every route that stops short, states `absence`:

| `absence` | when | it says |
|---|---|---|
| `graph` | where the pair stops, it recorded no event for c, and its record there is complete: the node is a `dead_end` (no successor, or only skipped hairpins, none of them c); or, in annotate mode with `per_path` scope, the pair turns off below its `complete_to_bp` and `evidence.complete_to_bp`, where every successor is either followed or recorded | the pair's graph, at its k and regime (stated), has no next k-mer that continues this walk |
| `walk_rule` | where the pair stops, it recorded a `blocked` event (`edge_reuse`, `edge_reuse_rc`, `rejoined_seed`) or a skipped `hairpin` event (`followed` false; SPEC §6.9, `walker.cpp:6161`, `Dh` in constrain mode) whose character is c (`Event.char`, `model.py:122`). A followed hairpin has a child, so it is no turn-off | the graph has the next k-mer, and the traversal rule (no edge twice, no seed entry, `hairpins`; SPEC §6.6) stopped the walk: another k, regime or `hairpins` setting can go on |
| `under_strategy` | inside the certificate otherwise: constrain mode, a `united_history` scope, or another end reason | the pair's request did not follow it, for example because no permitted label of that pair supports it; the graph may hold it |
| `uncertified` | otherwise, `to_bp` ≥ the part's `complete_to_bp` | nothing is known beyond `to_bp` in this pair |

Here c is the union walk's next base where the pair stops: after the last base of an `ends` or `diverges` span, or where
a route stops short. An event decides only for its own character: a pair that blocks A at a node by seed re-entry, where
another pair goes on with G, has not blocked G, and if it records nothing for G there its graph lacks that k-mer.
`graph` and `walk_rule` rest on what the walker recorded at that node for c, so they hold at any depth and come first;
`uncertified` comes next, and `under_strategy` last. A `dead_end` at the certificate's edge is therefore `graph`, as
g7's end at 120 bp in §2.9's example is.

A `none` span is `under_strategy` (the request named the other arm). Each span also carries the part's mode, and
`complete` and `certified_to_bp`: for `whole` and `ends` spans, and for a route that ends at a part leaf, the part
walk's own (`complete`, and `length_bp − beyond_certified_bp`, as `_walk_end` computes them, `ops.py:1463`);
otherwise `certified_to_bp` is the smaller of where the span or route stops and the part's `complete_to_bp`, and
`complete` says that it reaches that stop (`extends`, `merges` through the route) or is false (`diverges`, `none`).

**Which labels carry the span.** A label can be alive at a span's end without supporting all of its bases: a
lineage whose kept entry came through a non-first parent supports the spelled bases only from that merge on (SPEC
§7.1). So every span label is a `SpanLabel(label, from_bp)`, where `from_bp` is where its support of this span's
bases begins in this pair; a label with `from_bp` 0 carries the whole span.

- On a displayed walk (`whole`, `extends`, `ends`): the labels alive at the span's end (the part's alive set at that
  depth, not the leaf-anchored displayed profile, which misses labels alive at a cut depth: on the committed
  `merge.mgt`, `a.fa` is alive at 40 bp of walk 0 and carries [0, 40), but is not in the profile there). `from_bp`
  is the latest merge on the spelled chain, at or before the span's end, that the label did not enter through that
  chain's parent (`Segment.partition`), or the start of its run for a label switched in later; 0 when neither.
- On a `merges` span at the stub's end: the labels alive there, and `continues`, those whose kept entry entered
  `into` through this parent (`Segment.partition`, `model.py:176`). A lineage not kept at the merge ends its
  history there; its label may go on below the merge with another parent's history.
- Along a route: a label keeps its stub `from_bp` only if its kept entry entered every crossed merge through the
  route's parent; otherwise its `from_bp` is the last crossed merge where it did not.
- Annotate mode has no lineages: the labels recorded per node, `from_bp` where the label's run of presence along
  the span begins.

`carriers()` lists every span label with its `from_bp`; `walks(labels=)` keeps the walks those labels carry whole
(`from_bp` 0) in some pair.

Union walk ids number the terminal markers in trie preorder (a node's marker before its children, children by
base), walks and stub walks in one sequence; they are stable for one set of parts. Adding a part can renumber them;
the exports also name each walk by its length and pairs.

**The union of one part is that part.** The union walks of one graphlet with `kind: walk` are exactly the part's walks,
each `whole` with its spelling, length, end reasons, labels and completeness; the others are stub walks, each a `merges`
span whose route is `held` to the stub's end, with `goes_on` unless its merge is a zero-length leaf; and the coverage is
the part's arm certificate. This is the property test P2 (Appendix E.1).

### 2.5 Labels across pairs

Every label is a `PairLabel(pair, ref, name, kind)`, shown as `<tag>/<name>` (`g1/NC_000913.3`, `g4/562`). Refs are
index-local, so two pairs' refs never join, and names are never joined either. The server's certificates
cannot license a join: `headers_unique` true means no header repeats "in a pair or across pairs" (envelope spec
§9), and `columns_disjoint` true means no column name is a column of two pairs (SPEC §10.3). Wherever they hold, no
name occurs in two pairs, so there is nothing to join; where a name does occur in two pairs (a record in both
graphs), they do not hold. A join rule would only decide whether the tag is shown, so the tag is always shown and
the union reads no certificate.

- `label({'name': n})` resolves when exactly one `PairLabel` of the union has that name; otherwise it is refused
  `ambiguous_name`, listing the pairs. `{pair, ref}` always resolves within the pair. `{pair, name}` resolves as the
  part's own `label({'name': n})` does (`ops.py:238`): one label of that pair, otherwise `AmbiguousLabel` listing
  the refs, since one pair can hold two labels of one name (two headers ACC1 in two columns, the `same_name` fixture).
- `shared_names()` lists the names seen in two or more parts, with their pairs. The labels statement gives their
  count: "equal names in two pairs are not claimed to be one record (3 names occur in more than one pair)".

A header name rests on its pair's record mapping. On a pair whose `seqs_order` is `"count_only"`, as on every r237 pair
today, that mapping's column order was not checked at load, so a header name can name another record even inside one
pair (envelope spec §9.1). The tag does not repair that, and the library states it for no graphlet today. The owner
decided (OQ20) that the union warns per pair. `union_pairs(..., seqs_order=)` takes each pair's value as `GET
/capabilities` states it (`graph_summary.pairs[]` or `pair_summary`, envelope spec §9.1): `verified`, `count_only`, or
`None` where no `.seqs` is loaded. The fan-out helper fills it in (`PairResults.seqs_order`); a pair not given reads
`unknown`, as do the service's records unless the caller passes their pairs' values. Every header label (`h:` ref) of a
pair whose value is `count_only` or `unknown` is marked as resting on an unchecked record mapping: in the labels
statement, in `carriers()` and its printout, and in the exports (FASTA `mapping=`, a GFA comment line). Column labels,
`verified` pairs and pairs without `.seqs` carry no warning.

### 2.6 Coverage per pair: never more complete than the parts

`coverage(arm, at_bp=None)` puts every pair in one category at the depth `at_bp`. Without `at_bp`, the depth is the
largest radius (`bounds.max_extension_bp`) among the parts, for `by_category`, `summary()`, `print()` and the exports
alike, and `per_pair` also gives each pair's own certified depth.

| category | the statement | from |
|---|---|---|
| `complete` | the part holds every walk from this seed in this pair's graph, under the part's strategy, per path, to `at_bp` or to its end before it | the arm's `complete_to_bp` ≥ `at_bp`, no `seed_labels` limitation, and no `scope` limitation with `observed` > 0. A cap that cut the arm beyond `at_bp` does not matter: the certificate reaches `at_bp` |
| `complete_qualified` | as `complete`, with the qualification named: `seed_labels` (carriers of the seed were cut: "complete for 1000 of 4312 seed labels"; `max_seed_labels` is 1000 by default, SPEC §5) or `scope` (a merge united histories: complete for the united-history rule, not per path) | the arm's `complete_to_bp` ≥ `at_bp`, and a `seed_labels` limitation or a `scope` limitation with `observed` > 0 (SPEC §7.0) |
| `incomplete` | certified only to `x` < `at_bp`, with the reason | the arm's `complete_to_bp` < `at_bp`, and the reason from its `status` and `cap_trigger.reason` (truncated, pruned), else the part's radius |
| `seed_not_in_graph` | the graph does not hold the seed as one sequence, so no walk from it exists there | an absent record with `outcome.walks: "not_in_graph"`, and `not_in_graph.kmers_present` of `kmers` |
| `labels_not_in_graph` | none of the seed labels is in this pair, so nothing was walked | an absent record with `outcome.walks: "labels_not_in_graph"` |
| `not_walked` | a seed result without arms, with its cause; "the seed is in this graph" where the cause is a derivation, whose step follows the presence check (SPEC §6.1 steps 2 and 4: "nothing carries the seed in full") | an absent record with `outcome.walks: "failed"` |
| `arm_not_walked` | the part's request did not walk this arm | the part's `direction` |
| `failed` | the pair's answer is an error, so nothing is known | a `{pair, error}` record, its text cut at 60 bytes. On a single-graph server a seed missing from the graph fails the whole request with a 400 (SPEC §6.1 step 2), so in a union of two databases it arrives here, and the quoted text says so |
| `not_included` | the request went to this pair, and the union has no result for it | `expected_pairs=` |

`every_pair_certified_to_bp` is the depth to which every pair is `complete`, `seed_not_in_graph` or
`labels_not_in_graph`: the smallest such part's `complete_to_bp`. It is `None` when any pair is in another category,
`complete_qualified` included, and `None` without `expected_pairs=`, since a bare number would read release-wide.
Without `expected_pairs=`, every statement is about the pairs given and says so ("over the pairs given;
expected_pairs= was not passed"). A pair is never stated complete beyond its part's
certificate, and no union-level statement is stronger than the weakest pair it covers.

**A union without parts.** When no result carries a graphlet, the union has no arm, alphabet or radius. Its
coverage is arm-less: `coverage(arm=None)` puts each pair in its absent category, which holds for both arms, and the
sentence starts "both arms:"; `arms()` is empty and `walks()` returns nothing.

The sentence (one per arm and depth) lists the non-empty categories in the table's order, pairs in canonical order,
each pair's mode where any part is not in annotate mode, the pairs the caller stated, each absent record's seed
check, and ends with the scope of the paths. For example (nine illustrative pairs, g1–g9, to 300 bp):

    right arm, to 300 bp (constrain, kmer, in every pair): complete in g1, g4, g8; complete, qualified, in g6 (complete
    for 1000 of 4312 seed labels); incomplete in g7 (certified to 120 bp: max_steps); seed not in graph: g2 (0 of 10
    k-mers), g3 (4 of 10 k-mers); not walked: g5 (seed in graph; derivation: no label carries the whole seed);
    not included: g9. Pairs stated by the caller: g2, g3. Seed of the absent records: stated by the caller for g5,
    not checked for g2, g3. A pair that does not follow a walk did not follow it under its own request. A walk
    that no single graph holds entirely is in no part, so not in this union.

### 2.7 Path modes: source-supported and inferred

Should the union give source-supported paths only, or inferred paths across graphs as well (OQ13)? The owner's
provisional answer: "Expose both as separate modes"; and, decided: "Ship source-supported union first. Keep inferred
paths in the design and implement them in a later stage" (OQ18). Every query that returns or exports paths takes
`paths='supported'` (the default) or `paths='inferred'`; coverage, labels and the statements other than the fifth do not
depend on it.

**Source-supported (the default).** Every path returned or exported is held end to end by at least one source graph: a
union walk is a spelling from the seed that some pair holds whole, displayed or through a route of its own merges
(§2.4). Supported means supported by a source graph's recorded traversal; it does not show that one biological record
contains the whole sequence (under `merge`, a route below a merge rests on united histories, §2.4). The fifth statement
is the non-claim: **"A walk that no single graph holds entirely is in no part, so not in this union."** Every walk lies
inside one graph, and every end reason is that graph's: a `dead_end` in g7 without a skipped hairpin means g7 has no
next k-mer at its k and regime, not that the release has none (§2.4, `absence`). A seed only partly present in a graph
contributes nothing from it.

What it costs: on r237 each graph is a taxonomic group, so a real sequence (a chromosome, a plasmid) lies within one
graph, and a walk across graphs would join two organisms. A sequence shared by two groups (a plasmid, a transposon, a
conserved gene) is held entirely by both graphs, so its bases are in the union with both pairs and each pair's labels:
as a union walk where a pair displays it, otherwise as a stub walk whose merge leads into the walk below (`held_in`
through the route). The supported mode lists displayed walks and stub walks, not every route through the merges: where
both parts display the same allele of a shared bubble, the other allele is listed only as a stub walk, and the inferred
mode lists its full route. Where a release splits its graphs by something other than organism, this mode can miss real
walks, and it cannot measure that.

**Inferred (on request; designed now, built in a later phase).** A path may switch from one graph's spelling to
another's, but only at a connection the input graphlets justify. The connection rule of the first version admits two
kinds, and nothing else:

- **merge.** A pair P's stub ends where P's own path enters its merge, at trie node m (§2.4). P's graph joins both
  branches in one node there, so the path may go from the stub's end through m and on along any union segment below m,
  whichever pairs spell it (to each child of m where the merge segment has no bases). Where one pair holds the path on
  both sides, it is that pair's route and no transition is counted. Such routes can be many: the committed `merge.mgt`
  has two bubbles in a row, so its one part has one displayed walk but four routes from the seed to the leaf; the
  default mode lists the walk and two stub walks; the inferred mode lists four paths without a transition, the walk
  (`kind: walk`) and three routes. B bubbles in a row give up to 2^B routes, whatever `max_transitions` is.
- **seed** (GFA only, where both arms meet). The left flank of one pair and the right flank of another meet at the
  seed, which every part holds as the same oriented sequence.

No other junction is made: no join of sequences that lie at different outward positions or in different branches of
the trie, no join by k-mer equality across graphs (their k and regimes can differ), no join through a partial seed.
That would be sequence-overlap assembly, which the first version leaves out (provisional, OQ15).

**Pieces and transitions.** A pair holds a stretch of an inferred path when it spells those bases at those positions
as one path of its own: a displayed walk, a stub, or a route through its merges. A pair that spells another branch
can still hold a stretch: on the other side of a bubble, both alleles have the same bases from the first base after
the variant on, so the pair spelling the first allele holds the path from there. The **pieces** cut the path wherever
the set of holding pairs changes. A **transition** is where the path passes from the pairs that held it before to
pairs that hold it after, with no pair holding both sides; it has two positions: `from_bp`, where the after side
starts to hold the path, and `to_bp`, where the before side stops (for a merge, where P's route leaves the path,
which can lie well below the merge itself). Both sides hold `[from_bp, to_bp)`. An inferred path's id is `i:` and the
union segments it runs through (`i:4.31.52`): it does not depend on the listing, the ranking or `max_transitions`, and
`spell()` and the exports take it.

An inferred path follows three rules:

1. **Evidence stays per piece.** The path is a list of pieces, each held by named pairs. Support for each piece does not
   show that any sample holds the whole path. Each inferred path has a `kind`: `walk`, a union walk (listed in both
   modes); `route`, held whole by one or more pairs through their own merges but not a union walk; `joined`, with at
   least one transition. A walk or route names the pairs that hold it whole in `held_in`, as a `merges` span's route
   does; a joined path has no `held_in`. Labels are per piece and per pair, and no label is joined across a transition.
2. **Every transition is named**: its kind (`merge` or `seed`), the trie node and position of the merge (`at`, `at_bp`),
   every pair whose merge there permits it (`by`, in canonical order, so the witness does not depend on the
   implementation), `from_bp` and `to_bp`, and the pairs before and after.
3. **Completeness is stated narrowly.** Coverage stays per retrieval. Inferred mode lists the paths the connection rule
   admits over the parts given, within `max_transitions` (default 1) and `max_paths`, and states nothing about which
   biological sequences exist: each piece carries its pair's `certified_to_bp`, and the path has no `complete` of its
   own.

The fifth statement then reads: "Joined paths connect graphs at the transitions they name, and no graph is shown to hold
them whole; routes are held whole by their graph through its own merges; no inferred path is shown to be contained in a
sample." A union walk is listed in both modes, without a transition; a route only in inferred mode. Inferred paths can
be many even without a transition: routes multiply at every merge, and each transition multiplies the choices below a
merge. So the enumeration is bounded by `max_paths` (default 1000), independently of `max_transitions`, and by
`budget=`; it runs in path-id order, so a cut is the same in every run. The listing is an `InferredListing`: the paths,
`count` (how many were enumerated), `cut` (true when `max_paths` or the budget stopped the enumeration) and `resume` (a
token that continues it in path-id order, as `Partial` does, `budget.py`). A caller who receives exactly `max_paths`
paths reads `cut` to know whether more exist. `top=` and `by=` rank only the paths enumerated: with `cut`, the ranking
is within the first `count` paths in path-id order, never a global top, and the listing says so. Each FASTA record of an
inferred listing carries `of=<count>`, with `+` when it was cut. This contract belongs to the later inferred phase. What
neither mode gives is a walk of one traversal over several graphs, which no engine offers today.

### 2.8 Python API

New module `union.py`, exported from `metagraph.traverse`; the parts container, refusals, evidence and budgets are in
`composite.py`, shared with `_stitch.py` (Appendix F.2). Everything is phase 1a-U except what is marked `[1c]` and
`[inferred]` (the later inferred phase, §7.1).

```python
union_pairs(results, *, seed=None, expected_pairs=None, seqs_order=None, budget=None) -> PairUnion
    # results: Graphlet | TraverseResponse (one seed) | {response, index} | {response, sequence}
    #          | {pair, graphlet} | {pair, result[, error][, sequence]} | {pair, error}
    # seed: the oriented seed sequence; required only when no result carries a graphlet
    # expected_pairs: the pairs the request was sent to (graph names or pair objects); None: not known
    # seqs_order: pair -> 'verified' | 'count_only' | None, from GET /capabilities; a pair not given: 'unknown'
    # A new object every time; parts held by reference, never copied or modified.

class PairUnion:                          # every method also takes budget= (keyword-only, see below)
    parts -> tuple[PairPart]              # canonical order (§2.10)
    absent -> tuple[PairAbsent]
    pairs -> tuple[str]                   # every tag (unique, Appendix A.4): parts, absent records, missing pairs
    expected -> tuple[str] | None
    seed -> SeedInfo                      # the common oriented seed; each part's own seed record in parts
    alphabet; arms() -> list[str]         # the arms any part walked
    statements(arm=None, at_bp=None, *, paths='supported') -> Statements   # the five statements (§2.1)
    coverage(arm, at_bp=None) -> Coverage                 # §2.6
    labels(pair=None) -> list[PairLabel]
    label(selector) -> PairLabel          # {'pair', 'ref'} | {'pair', 'name'} | {'name'} when only one pair has it
    shared_names() -> dict                # name -> tags, for names seen in two or more parts
    segments(arm, *, budget=None, resume=None) -> list[UnionSegment]
    walks(arm, *, paths='supported'|'inferred', max_transitions=1, max_paths=1000, top=None, by='pairs'|'length'|'id',
          pairs=None, labels=None, min_bp=0, spell=False, stubs=False, budget=None, resume=None)
          -> list[UnionWalk] | InferredListing
    # [inferred] paths='inferred', max_transitions, max_paths and InferredListing, here and below
    spell(arm, walk, orientation='natural'|'walk', with_seed=False, *, budget=None) -> str  # a walk or inferred id
    carriers(arm, walk, *, at_bp=None, budget=None, resume=None) -> Carriers  # a list of PairSpan; str(): §2.2
    to_fasta(arm=None, walks=None, with_seed=True, orientation='natural', width=None, *, paths='supported',
             max_transitions=1, max_paths=1000, at_bp=None, stubs=False, budget=None) -> str
    to_gfa(*, paths='supported', with_seed=True, labels='counts'|'full', at_bp=None, budget=None) -> str
    evidence(arm=None) -> dict
    summary(max_bytes=4096, *, paths='supported', budget=None) -> dict   # the statements first and whole
    __str__() -> str                      # the statements as sentences: print(union)
    provenance() -> list[dict]
    save(path, *, budget=None)            # [1c] the bundle (Appendix D.3)
    digest -> str
    memory_bytes() -> int
    # Refused (CombineRefused 'not_offered', §2.3): dump(), to_json(). Not in 1a-U (§7.1): compare() of two
    # unions, support_by_depth(), carriers(sequence=), label routes (routes=), a per-pair FASTA,
    # next_request() per pair. Until the inferred phase (§7.1), paths='inferred' (with max_transitions and
    # max_paths) is refused 'not_offered' too.
    # the data classes (PairPart, PairAbsent, PairLabel, SpanLabel, LabelInterval, SegmentSupport, UnionSegment,
    # Route, PairSpan, Carriers, UnionWalk, InferredPath, Transition, Coverage, Statements) are in Appendix C.1
```

**Tags.** Every pair is named by its display tag, unique in the union (Appendix A.4); selectors take a tag, and
`parts[i].key` gives the full pair key.

**Budgets.** Every operation takes `budget=`, keyword-only, with the library's
semantics (`budget.py`): a stop raises `LocalBudgetExceeded`; the lists whose order allows it (`walks()`,
`segments()`, `carriers()`) carry a `Partial` of whole rows with a resume token for `resume=`, and the others stop
without one. `union_pairs()` charges construction (the trie pass) to its `budget=`.

**The queries.**

- `walks(arm, by=...)`: `'pairs'` ranks by the number of pairs in `held_in`, then in `whole_in`, then length, then id;
  `'length'` by length; `'id'` keeps trie order. With `paths='inferred'` it returns `InferredPath`s ranked by the number
  of transitions, then length, then id; `'pairs'` is then not offered. `pairs=` keeps the walks those pairs hold whole
  (`held_in`); `labels=` the walks those labels carry whole (`from_bp` 0) in some pair. Stub walks only with
  `stubs=True`.
- `carriers(arm, walk, at_bp=None)` answers "which pairs and labels carry this extension": per pair, how far it
  holds the walk (displayed or through its route), how it stops, with `absence` and mode, and its span labels with
  their `from_bp` (at `at_bp`, default each span's end).
- `spell()` spells a union walk from the parts' bases.
- No cross-pair label total is offered. Label counts are per pair, as label entries of that pair, and
  marked as lower bounds where the pair's lists are cut.
- `summary(max_bytes=4096)` puts the statements first and whole; if they alone exceed `max_bytes`, it refuses
  `summary_too_small` rather than cut them. For r237's 11 pairs and two arms they take about
  1–2 KB [estimate].

### 2.9 Exports

**FASTA.** One record per union walk, or per inferred path with `paths='inferred'` (the parts' own walks are exported
by each part's `to_fasta()`). Each header states the path's provenance, its coverage and the scope as `key=value`
tokens, empty ones left out, on one line in the file:

    >right.u12 kind=walk len=812 held=g1,g4,g6 whole=g1,g4 merges=g6@85>s31 ends=g7@120:dead_end/graph
     diverges=g8@40/under_strategy certified=g1:812,g4:812,g6:812,g7:120,g8:40 labels=g1:3/2,g4:1+/1,g6:2/0
     modes=constrain absent=g2:seed_not_in_graph,g3:seed_not_in_graph,g5:not_walked not_included=g9
     mapping=count_only scope=within_one_graph

and, in the later inferred mode, for the inferred path of §2.2:

    >right.i:4.31.52 kind=joined len=860 pieces=g06:0-55,g01+g06:55-700,g01:700-860
     transitions=merge:g06@s31:85[55-700] certified=g06:700,g01:860 scope=inferred

Appendix A.5 gives the token grammar and the escaping that makes every token decode back to the names it holds.

**GFA.** GFA 1.0 on the seed strand, with the statements as comment lines first. Its topology respects the path mode:
`0M` links concatenate segments for any GFA reader, whatever the tags say, so the file admits exactly the paths of its
mode. Segments are union segments, and links overlap `0M`: the union is a trie of spellings, not a de Bruijn graph,
and its parts' k can differ. Left-arm segments are written in natural orientation, as `to_gfa()` writes them.

- **Source-supported** (the default): one tree per arm. A union segment links to its children only, so every path the
  file admits runs down one branch of the trie and is spelled whole by each pair that spells its last segment. The seed
  is written once per arm (`seed_l`, `seed_r`, one sequence; the second comment line says so, and that no path of the
  file runs from a left flank through the seed into a right flank), because one seed segment joining both arms would
  admit a left flank of one graph followed by a right flank of another. A stub's merge is not a link: the stub's last
  segment carries `mg:Z:<tag>` and `mt:Z:<trie node>`, so a stub does not look like a dead-end tip and admits nothing
  more. A pair's route through its own merge is therefore not a path of the file; `carriers()`, the FASTA `held=` and
  the P lines' `ph` tag state it, and each part's own `to_gfa()` admits it.
- **Inferred** (the later phase): the same trees, one seed segment joining both arms, and an `L` link from each stub's
  last segment to the union segment that starts at its merge's trie node (to each child of that node where the merge
  segment has no bases; none where it is a leaf), tagged `mg:Z:<tag>`. The paths the file admits are exactly the paths
  of the connection rule (§2.7), with any number of transitions.

The default mode's file, abridged:

    H	VN:Z:1.0
    # union_pairs: seed 40 bp; parts g1 g4 g6 g7 g8; seed not in graph g2 g3; not walked g5; not included g9
    # paths: source-supported; seed_l and seed_r are one seed, written once per arm, so no path of this file runs
#   left flank, seed, right flank: each part's own GFA has those paths
    # right arm, to 300 bp (constrain, kmer, in every pair): complete in g1, g4, g8; ...
    # turn-offs: a pair not following a walk did not follow it under its own request (absence on each span)
    # labels: per pair, as <tag>/<name>; equal names in two pairs are not claimed to be one record (3 shared)
    # labels: header names of g1, g4, g6, g7, g8 rest on an unchecked record mapping (count_only)
    # scope: a walk that no single graph holds entirely is in no part, so not in this union
    S	seed_r	CGTACC...	LN:i:40
    S	s12	ACGTT...	LN:i:300	pp:Z:g1,g4	pl:Z:g1=3;g4=1+
    S	s30	GGATC...	LN:i:85	pp:Z:g6	mg:Z:g6	mt:Z:s31
    L	seed_r	+	s12	+	0M
    P	right.u12	seed_r+,s3+,s12+	*	pp:Z:g1,g4	ph:Z:g1,g4,g6	pe:Z:g7=120:dead_end/graph;g8=40/under_strategy

The lower-case tags are user tags: `pp` the pairs spelling a segment or following a path whole, `ph` the pairs
holding a path whole (displayed or through a route), `pl` the label entries per pair (counts, or with
`labels='full'` the names, `<tag>=<name>,<name>;…`, each tag once and each name escaped on its own), `pe` the other
spans of a path, as in FASTA, `mg` and `mt` a stub's merge. Tests check the
paths the file admits, not its lines (Appendix E.1, P9).

**Comparing.** A union-level `compare()` is not in 1a-U (§7.1). Each common pair's two retrievals compare with the
parts' own `compare()` (`ops.py:3778`).

### 2.10 Determinism, evidence and provenance

A union is a pure function of its inputs: the results, the stated pairs, `seed` and `expected_pairs`. Any order of
the results gives the same digest and the same answers, byte for byte (a refusal names the first failing element in
the order given). The union of one part answers as that part, and the same retrieval twice on one pair is a no-op,
the envelope-stated form kept. `evidence()` lists each part's own evidence block under its tag, and an answer that reads
several parts states the weakest outcome per dimension over them; `exact` holds only when every part read is exact
and verified. `provenance()` lists each part's replay request; replaying a union replays each part and builds the
union again. Appendix C.3 states these rules in full.

## 3. Stitching a continuation (phase 1a-S)

Stitching joins a continuation (a new traversal seeded with a walk's continuation, `next_request()`, `deepen()`,
`traverse_continue`) to the walk it continues, on one pair. The owner's provisional answer is to build it as phase
1a-S, after the union (OQ1). Its claims and support across junctions are phase 1b (Appendix D.1).

### 3.1 The contract

New private module `_stitch.py`, its names exported from `metagraph.traverse` (so `metagraph.traverse.stitch` is
only the function), on the composite machinery of `composite.py` that the
union builds first (Appendix F.2). Everything is phase 1a-S except what Appendix C.2 marks `[1b]` and `[1c]` (§7.1).

```python
stitch(parent, child=None, *, arm=None, walk=None, allow_unverified_index=False, budget=None) -> StitchedGraphlet
    # parent: Graphlet (with envelope) or StitchedGraphlet; child: Graphlet (with envelope) or StitchedGraphlet
    # rooted at the continuation; None -> the parent alone as a StitchedGraphlet with no junction.
    # arm: 'left' | 'right' | 'l' | 'r' (default: the child's one arm). walk: the continued walk, a path id on a
    # Graphlet parent, a walk key (tuple) on a StitchedGraphlet parent; default: inferred (§3.1). One junction
    # per call; a new object every time (parts held by reference, never copied or modified).

# class StitchedGraphlet (its methods) and the data classes Part, Junction, StitchedWalk, StitchedClaim
# and CombineRefused are in Appendix C.2
```

A walk key is the tuple of walk ids along the chain: `(3,)` is root walk 3, and `(3, 0)` is walk 0 of the child that
continues root walk 3. Positions in every answer are outward from the root seed's boundary. Each row also carries its
pieces in their own parts' coordinates. A label selector `{'id': n}` is refused on a stitched graphlet (BadSelector,
"label ids are local to each part of a stitched graphlet: select by ref or name"). Refs and names resolve over the
union.

The methods mirror a graphlet's: `walks()` and `spell()` over stitched walks, `continuation()` and `next_request()`
delegated to the part that owns the walk, `to_fasta()`, `completeness()`, `evidence()`, `summary()`, `provenance()`,
`digest` and `memory_bytes()`; `claims()`, `support_profile()` and `support_changes()` are phase 1b, and `save()` the
bundle of phase 1c.

The refusals, in the order they are checked: `not_combinable`, `view_unsupported`, MissingEnvelope, `different_index`,
`index_unverifiable`, `different_pair`, `different_graph`, `different_mode`, `not_one_arm`, `wrong_arm`,
`no_continuation`, `not_a_continuation`, `ambiguous_walk`, `already_continued`, `foreign_seed_labels`, and `not_offered`
for queries. Their exact texts are in Appendix B.1.

Checks run in that order, and the first failing check decides. The identity checks reuse the index rules of
`_identity` (`ops.py:3588`) **without** its seed rule, because parent and child seeds differ by design, and with the
one-sided-manifest rule of `_same_index`: refused unless opted in. When neither side has `index_fp`, the stitch is
not refused and the identity is stated as unverified, as `traverse_continue` does (OQ6, provisional).

- **Reachability.** Fingerprints are per pair (SPEC §10.3), so with manifests a different pair or graph is always
  refused `different_index` first. `different_pair` and `different_graph` are reached only when no fingerprint
  decides (no manifest on either side, and `index_meta_fp` equal or absent). The pair is not checked first on
  purpose: a fingerprint mismatch is proof, while pair tags can name one index under two names (a pair listed under
  two names counts once, envelope spec §9).
- **Walk inference (one input, one code).** Without `walk=`, the candidates are **all** walks of the arm whose
  continuation sequence equals the child's seed, continued or not. Several candidates: `ambiguous_walk`. One: that
  walk. The answer therefore does not depend on which children were stitched before.
- **An identical retrieval is a no-op.** A child with the body digest and the retrieval
  identity (§2.3) of the one already stitched at that walk returns the stitched graphlet unchanged, so rebuilding from
  a set of parts tolerates duplicates. A different body, or the same body from another request, at a continued walk
  is refused `already_continued`: the two would continue differently.

### 3.2 Coordinates, orientation and the overlap

Let the parent's walk w on arm a end at outward position L (its leaf's `end_bp`). Its continuation has n bases. The
walker chooses n (`make_continuation`, `src/graph/traversal/walker.cpp:3117-3128`; SPEC §7.1, `output.continuation_bp`
in SPEC §5):

    n = min(continuation_bp, max(k, L − b), |seed| + L)

Here b is the latest boundary among the labels alive at the leaf: where a run was entered by a switch, or a route
joined at a merge (`route_bp`); 0 when neither happened. `continuation_bp` defaults to 1000. So n ranges from k to
1000, and the child's seed is often a large part of the parent's walk. Without a switch or a route join and with
k ≤ L ≤ 1000, it is the **whole walk**. When L < k it reaches into the root seed. The committed fixtures show
both: `same_name` (L, n) = (40, 40), the whole walk; `fork` (10, 15), into the seed; `merge` (100, 38), cut at a
route join; `caps` (126, 23) and `switch_chain` (80, 20), cut at switches. The child's seed nodes are therefore all
nodes of the parent's walk from L − n on (with root-seed nodes when n > L). §3.4 uses this.

`continuation()` (`ops.py:2238`) places the n bases at seed coordinates `[|seed| + L − n, |seed| + L)` on the right
arm and `[−L, −L + n)` on the left. The child's seed is that sequence, in the request orientation (SPEC §7.5.2
`S`). Its one arm (direction a) starts at the child seed's outward boundary, which is the parent's outward position L
**on both arms**. So:

- child outward base i = composite outward base `L + i` (offset `L`; chains add offsets);
- walking order: stitched walk = parent's `walk(w)` + child's `walk(v)`; natural orientation: right flank =
  `natural_P(w) + natural_C(v)`, left flank = `natural_C(v) + natural_P(w)`. The whole molecule is
  `natural(left) + root seed + natural(right)`, as for one graphlet (SPEC §7.5.4);
- the child's seed is never spelled twice: its n bases are the parent's walk from L − n (and root-seed bases when
  n > L);
- a child seed coordinate x (child seed = `[0, n)`) is root seed coordinate `x + |root seed| + L − n` on the right
  arm and `x − L` on the left. `StitchedGraphlet.continuation()` maps the owning part's `seed_coord` this way,
  summing over a chain.

No reverse complement is involved. On primary graphs a seed and its reverse complement are one seed, but the child's
seed is written in the parent's request orientation, and its walks are spelled on that strand.

**Worked example (the committed `fork` fixture).** The seed is 40 bp, `CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATT`,
and right walk 0 stops at the radius, L = 10: `fork.spell('right', 0)` is `ATTGCTAAGA`. Its continuation is
`CTATTATTGCTAAGA` (n = 15 = k, since L < k), seed coordinates `(35, 50)`, labels `['acc2']`. Its first 5 bases,
`CTATT`, are the seed's last 5, and its last 10 are the walk. A child seeded with it and walked to the right has
outward base 0 at composite position 10. Its walk v becomes stitched walk `(0, v)` of length `10 + child length`,
spelled `ATTGCTAAGA + child.spell('right', v)`. `fork` is also a case for §3.4's `sibling_merge_unknown`: it runs
`merge`, and its right arm has a second walk with a continuation.

### 3.3 Walk identity and labels

A stitched graphlet's walks on arm a are the root's walks, except that every continued walk w is replaced by `(w, v)`
for each walk v of its child (recursively for chains). Nothing is renumbered: `(w, v)` names root walk w and child
walk v. Continued walks of one parent with identical tails are not deduplicated, because each child is its own trie;
under `merge` such children are a `sibling_merge_unknown` case (§3.4). The overlap is exactly the child's seed
(`overlap_bp = n`, the value `traverse_continue` stores), so nothing is double counted. The child's positions start
at `L`, where the parent's walk ends.

Labels join by `ref` across parts. This is valid because the identity check proved one index, or the caller accepted
an unverified one, which is then stated. Names must agree per ref (`different_index` otherwise, as `compare()`'s
names rule). In constrain mode the child's seed labels (its first `num_seed_labels` dictionary entries) must be alive
at the end of walk w (the leaf's end set, `G.end`), otherwise `foreign_seed_labels`. Labels alive at w that the
child does not seed are listed in the junction's `not_continued` (`alive_not_seeded` as `next_request` reports it,
or simply not named). In annotate mode the child has no seed labels (`ops.py:2769`), and none of this applies.

### 3.4 The junction statement: two verdicts

**The reference.** One uninterrupted walk of the parent's request with a larger radius and larger caps. A
continuation is a new traversal (DESIGN §5.3): branch state, loss budget, edge history, seed re-entry, merges with
other walks and (under `trace`) the followed occurrence all restart. Each junction derives, from the parent (its DAG,
the leaf's `T` and `C` records, the runs ending there, its normalized strategy) and the child (its normalized
strategy, `X` records and body), whether the child's walks below w are the reference's walks below w.

**Why only two verdicts.** The first version also had `lower_bound` ("the reference holds everything the stitched walk
shows"). It is not monotone (Appendix G.1): fewer labels or less budget can also *lengthen* a walk. Fewer paths mean
fewer merges and so fewer united edge histories (SPEC §6.4, §6.5). Fewer non-empty successors mean fewer splits, so the
reference can end with `split_limit` where the child goes on (SPEC §6.4; default 64). Under switching with a finite
branch allowance, a source that is ambiguous in the reference can be unambiguous in the child (SPEC §6.4; the
re-minimisation is greedy, SPEC §6.3). A direction can be stated only under conditions a comparison against
uninterrupted traversals has confirmed. So `exact` is limited to a narrow whitelist, and phase 1a-S ships `qualified`
only (OQ19, decided). `lower_bound` can come later, under the conditions of Appendix G.1 (`keep`, an unlimited split
limit, no quorum, minimum or beam, `forbid` or unlimited branches) pinned by T4.

**The whitelist.** Once `exact` is switched on (§7.1), a junction is `exact` when all of W1–W7 hold, and a stitched walk
is `exact` when all its junctions are and W8 holds for it. In phase 1a-S every rule is computed and every failing
condition is named, but nothing is `exact`: where W1–W7 hold at a junction, its verdict is `qualified` with the one
cause `exact_gated` (W8 is judged per walk below it; the junction rules found no difference, and they have not yet
passed the randomized comparison), and a walk reads `exact_gated` where its junctions do and W8 holds for it.

**What the rules read.** W2–W5 and the split check of W8 read the two normalized strategies. A field they read that
either envelope does not state (W5's knobs; in constrain mode also `labels.change_cost`, `labels.extra`,
`branching.max_label_branches` and `branching.max_splits_per_path`), or that one part states and the other does not, is
the junction cause `strategy_unknown`, not a pass: a missing `max_splits_per_path` is unknown, not unlimited. So is a
field W2 compares (every field the server's normalized strategy states, but the allowed ones) that neither envelope
states; a null field is not stated. It is no refusal, since every verdict is `qualified` anyway; nothing below it is
`exact`, the prospective result included. W2 compares the fields both state. W6's `keep` or `merge` is read from each
part's `H` record, so a difference there is `strategy_differs` even where the envelopes do not state
`branching.on_reconverge`.

- **W1.** `support: kmer` on both parts.
- **W2.** The two normalized strategies are equal except in the allowed fields: `direction`; the radius and the
  resource caps (`bounds.max_extension_bp`, `max_steps`, `max_paths`, `max_output_bp`, `time_budget_ms`,
  `max_memory_mb`, `max_work_units`, and `max_live_paths` where `frontier.on_overflow` is `stop` on both); `output.*`;
  `annotation.*`; `labels.max_labels_per_node` (it caps output lists only, SPEC §5); `labels.loss_budget`,
  `labels.extra` and `branching.max_label_branches` (judged by W3–W4); and `clamped`, the server's report of the
  bounds it clamped, not a knob. The child's `labels.extra` names no label outside the parent's permitted pool (its
  seed labels and `labels.extra`).
- **W3.** `labels.change_cost.model` is `forbid` on both parts. Under `forbid` no switch happens, every loss is 0, a
  lineage is one label, and `labels.extra` cannot be entered. Switching (any other model) is `qualified` in phase 1a-S.
  The uniform-shift argument (every continued label at one terminal loss and branch count, with the child's budget
  and allowance reduced by exactly those) is a candidate for promotion once its T4 fixture pins it.
- **W4.** `branching.max_label_branches` is `"unlimited"` on both, or every continued label ended walk w with one
  branch count b and the child's allowance is the parent's minus b. Under `forbid` a label's entries all descend from
  that label, so the tie-break and the merge minimum, which read branch counts (SPEC §6.3, §6.5), shift uniformly
  within each label.
- **W5.** No knob that reads other labels' counts: `min_successor_labels` ≤ 1, `min_successor_fraction` = 0,
  `bounds.min_live_labels` ≤ 1, `frontier.on_overflow: stop` (no beam).
- **W6.** Reconvergence. Either `keep` on both, or `merge` on both, and every other walk of the arm ended without a
  continuation at a length ≤ L. Under `merge` the reference joins paths that reach one node at one outward distance
  and unites their labels and edge histories (SPEC §6.5). A sibling walk with a continuation, or one longer than L,
  could meet w's descendants beyond L, and the child never sees it. Merges on w's own route before L are no
  obstacle: the parent's DAG holds them, and the rule check (W8) bars the edges of every DAG ancestor of w's leaf
  through any parent, as the walker does (SPEC §6.6, "for merged DAG nodes, any parent").
- **W7.** Constrain mode: every label alive at the end of w is a seed label of the child (`not_continued` empty), and
  the child's `X` drops none of them. Annotate mode: no condition.
- **W8 (per stitched walk): the child decided as the reference would.** Every part from the root to the child
  carries bases (the edge sets of B come from all of them; otherwise `no_bases` at the junction, and W8 is not run
  below it). Five checks, on the
  walk's route through the child (the DAG ancestors of its leaf in the child, through any parent):
  - *What the reference blocks.* B is the set of steps the reference blocks and the child cannot know of: the
    (k+1)-mers of the DAG ancestors of w's leaf, through any parent (canonicalized in primary and canonical regimes),
    and the steps into root-seed nodes that are not nodes of the child's seed (SPEC §6.6).
  - *Admitted successors.* Every successor the child admitted is checked against B: the step it took and the other
    children of a split (`G`), and every successor of a branch event (`V`: its characters, refused and ambiguous ones
    included). An admitted successor whose step is in B changed the child's decision, **whether or not the child took
    it**: the reference does not admit it, so a label the child found ambiguous there, or a split the child made, may
    not happen in the reference (SPEC §6.4).
  - *Blocked successors.* The successors the child blocked (`E`) are checked the other way: a child `rejoined_seed`
    into a node of its seed that is not a root-seed node changed the decision unless that step is in B too (§3.2: a
    node of the parent's walk from L − n; `same_stop` below), since the reference admits it.
  - *Split counts.* Splits are counted as the engine counts them, in parent and child alike: 0 at the root, one more
    on each child of a split, and the maximum over the parents at a merge (`walker.cpp:4774`, `4875`), so P, the
    parent's count at w's leaf, can exceed the splits on w's displayed chain. The reference refuses a split when its
    count plus one would exceed `max_splits_per_path` (`walker.cpp:5394`; finite by default, 64), and below the junction
    its count is P plus the child's, which starts from 0. So at every split the child made on the walk's route, P + c +
    1 ≤ `max_splits_per_path` must hold, c being the child's count at that node; otherwise `split_limit_restarted`.
  - *The evidence boundary.* A successor not taken at depth d has its recorded reason only when d <
    `evidence.complete_to_bp` (SPEC §7.0: "at or beyond it the reason may be among the events not kept"). The checks
    above are complete strictly below that depth, and certification stops there (§3.5). A walk whose end lies at or
    beyond it is `decisions_unrecorded`: certifying the bases before an end does not certify why the walk ended. The
    extreme case is a 0 bp child walk with a boundary of 0, the `loop_branch` child run with `max_branch_events` 0:
    its one branch event is not kept, and the reference goes on to 145 bp.

**Why the decisions, not the steps.** The first version checked only the steps the child took. Forgetting the parent's
used edges can change what the child sees before it takes any step. A counterexample, reproduced on the base binary
(Appendix G.4): the parent's walk ends at the radius, L = 115, back at a node it passed before, having used its one
branch allowance there; the reference blocks the edge the walk used then and has one exit; the child, seeded with an
allowance of 0, sees two exits, finds its label ambiguous and stops at 0 bp; the reference goes on to 145 bp. W1–W7
held, and the child took no forbidden step. Under the rule above the stitched walk is `qualified` (`candidate_changed`
at 115 bp), certified only to 115 bp, and not complete.

**The causes.** A condition that fails is stated as a named cause, with its scope (the junction or one walk) and the
modes where it can occur, and the verdict is `qualified`: `trace_restart`, `strategy_differs`, `extra_foreign`,
`switching`, `branches_unequal`, `label_count_knob`, `sibling_merge_unknown`, `not_continued`, `seed_dropped`,
`no_bases`, `strategy_unknown`, and per walk `rule_violation` (a step in B), `candidate_changed` (a recorded successor in B), `seed_reentry`
(`same_stop`, `earlier_stop`), `decisions_unrecorded` (the walk's end lies at or beyond the child's evidence boundary)
and `split_limit_restarted`; and, in phase 1a-S, `exact_gated`. Appendix B.2 states what each one means. Until `exact`
is switched on, every verdict is `qualified`, with `exact_gated` where no other cause applies; after that, a whitelist
entry or a promotion reads `exact` only with its T4 fixture and the randomized comparison passed. The verdict never
depends on notes `next_request()` returned: a child made another way is judged by what it ran.

### 3.5 Completeness and evidence: nothing is laundered

- **The arm.** `complete_to_bp` of a stitched arm is the root's, unchanged: stitching adds walks below continued
  walks only and never makes "every walk to depth D" true for a larger D. `completeness()` adds one extension per
  junction: `{walk, at_bp, extended_to_bp: offset + child.complete_to_bp, certified_to_bp, verdict, causes}`.
- **A stitched walk.** Pieces are checked in order. Piece i certifies its part's walk up to the earliest of its part's
  `complete_to_bp` and, inside a child, the first walk-scope cause on the walk's route (a changed decision, a step in B,
  the child's `evidence.complete_to_bp`): certification stops at the first relevant violation inside the child, not only
  at junctions. It passes to piece i + 1 only if piece i's walk ends within that depth and the junction between them is
  `exact`. In phase 1a-S no junction is `exact`, so a stitched walk is certified at most to its first junction, and
  every walk that reaches into a child reads `complete: false`; each piece still lists its own part's certificate. The
  rules still compute the **prospective result**: the same answer with `exact_gated` read as `exact`, that is, certified
  through such junctions, with `complete` and the evidence as they will read once `exact` is switched on. In 1a-S it is
  a private test hook (T4, T11) and never part of a public answer. `certified_to_bp` is where it stops,
  `beyond_certified_bp = length_bp − certified_to_bp`, and `complete` holds when that is 0, the last piece's walk is
  complete in its part (`Walk.complete`), and no walk-scope cause lies on the walk, at its end included. The
  `loop_branch` walk has `candidate_changed` at its end, 115 bp, after a child walk of 0 bp that its part calls
  complete: it is certified to 115 bp and not complete. `certified_to_bp` counts bases; whether the walk ends where it
  ends is `complete`. A parent walk cut by a cap beyond its arm's `complete_to_bp` therefore stays uncertified past that
  depth, however far its child reaches.
- **The evidence block** (`evidence()`, and with every answer): each part's own `evidence_block` (`ops.py:6006`) keyed
  by part, the junctions it depends on with verdicts and causes, the identity statement, and an `outcome` that takes the
  **weakest** value per dimension over every part and junction the answer reads. A `qualified` junction, and a
  walk-scope cause on any walk the answer reads, map to `walks: partial` and `label_evidence: qualified`. `exact` holds
  only when every part read is exact, every junction is `exact`, no walk the answer reads has a walk-scope cause, and,
  where the answer reads a junction, the identity is verified, so the CLI fixtures, which have no `index_fp`, can never
  read `exact` there (§6 T4). In phase 1a-S no answer that depends on a junction reads `exact`, and a `qualified`
  verdict never stands beside `complete: true` or exact evidence that the causes do not justify. An answer that depends
  on no junction, `stitch(g)` without a child or a root walk no child continues, states the root's own evidence, `exact`
  where the root's is, so `stitch(g)` still answers as `g` (T2). The verdict is `None` for a walk no junction joins,
  and its evidence is the root's. Each part keeps its `K` records, its `Q`
  (`resource_stop`), its dropped seed labels (`X`), its cut lists (`label_lists`, `inexact_counts`), its dropped branch
  events (`evidence.complete_to_bp`), its `coordinates` limitation and its delivery (`spooled`), all listed under its
  key. None is summed or dropped. A local-budget stop follows the library's rule (`budget.py`): `LocalBudgetExceeded`,
  with a `Partial` of whole rows and a resume token bound to the stitched graphlet's digest for the lists whose order
  allows it, and no partial otherwise. `completeness()` and `evidence()` run the walk rule (W8) over every walk they
  read, so both take `budget=` and charge it as `summary()` does, each walk's causes at their cold price whether the
  stitched graphlet holds them or not (F.1). The prospective hook is uncharged.

### 3.6 Exports, provenance and replay

- **Stitched FASTA** (phase 1a-S). One record per stitched walk, spelled from the parts without the overlap
  (§3.2), using the exports' streaming spelling. Headers are `to_fasta()`'s plus `walk=3.0`,
  `junctions=<at_bp,...>`, `verdict=<exact|qualified>`, `certified_to=<bp>`. A part without bases is refused as that
  part's own export would be.

`provenance()` returns one entry per part: `{key, role: root | child, junction: {arm, walk, at_bp, overlap_bp} | null,
pair: answered_pair(envelope), identity: index_identity(g), seed_id, digest, request}`. Here `request` is the part's
replay request as `GraphletStore.request_of()` (`store.py:1277`) builds it: the child's request names the
continuation seed and the parent's pair. No single request reproduces a stitched graphlet. Replaying one means
replaying each part (the identity rules apply as for `traverse_fetch(replay=)`) and stitching the results again.
Each child is then checked again against its replayed parent, so a changed index or walk numbering is refused rather
than joined.

### 3.7 Determinism and the algebra

Stitching is pure: junction verdicts are recomputed from the parts every time, never stored. A stitched graphlet is a
canonical tree with a digest. Stitching is commutative across junctions and associative along chains, but not between
parent and child: `stitch(C, P)` is refused in the usual cases, and where a circular sequence makes each one's
continuation the other's seed it passes the preconditions and is judged by W1–W8 like any stitch; an identical child is
a no-op and a different one at a continued walk is refused; `stitch(g)` answers every offered query as `g` does (T2).
Appendix C.3 states these rules in full.

## 4. What is offered later or not at all, and why

| operation | status | where |
|---|---|---|
| cross-pair union `union_pairs(results)` | phase 1a-U, first (built in train 7, source-supported mode) | §2 |
| continuation stitching `stitch(parent, child)` | phase 1a-S, after the union (built in train 7, experimental, `qualified` only) | §3 |
| two retrievals of one pair | designed, built later | §4.2 |
| the two arms of one seed against a `direction: both` fetch | later, a validation helper: as an answer a `direction: both` fetch serves it; as a check the three fetches should agree arm by arm (the owner's provisional answer to OQ2) | Appendix D.5 |
| different seeds, overlapping or not | left out: a different product | §4.3 |

Appendix A.3 gives every operation's preconditions in the order they are checked, and its refusal codes.

### 4.1 "Merge" is taken

In this codebase, merge already means reconvergence: `on_reconverge: merge`, a segment with two or more parents
(`Segment.is_merge`), a run closed by a merge (`m`, `Run.merged`), `completeness_scope: united_history`. A
`Graphlet.merge()` that combines retrievals would sit next to "a merge before the depth" in every compare note. The
operations are therefore named `union_pairs` and `stitch`. "Merging graphlets is not supported" stays true in the
docs.

### 4.2 Two retrievals of one pair: a later operation

The owner's provisional answer: same-pair combinations belong in the design, and their implementation is deferred.

**Why it is not the cross-pair union.** Which walks a traversal holds depends on the whole strategy: the permitted
set, the change cost and loss budget, the branch allowance, the quorum knobs, the beam and caps, and (under `merge`)
which edge histories were united and which parent is displayed (DESIGN §11). Two retrievals of one arm on one pair
are two answers to two requests. Their union is the answer to neither, and no `complete_to_bp` describes it. SPEC §6.4
states the one case where it is exact: "Equality with the union of independent single-label walks holds only under
`on_reconverge: keep`". In my reading it also needs no switching, quorum or beam **[unsure]**.

**The operation** (Appendix D.4 has it in full): `union_retrievals(results)`, whose sources are retrievals, told apart
by retrieval identity (§2.3), instead of pairs. It reuses the trie, spans, routes, absence, both path modes and the
exports of §2, states coverage per retrieval and never for the union ("This union answers none of its requests"), and
joins labels by ref only where the retrievals' `index_fp` agree. Its one exact case (the SPEC §6.4 conditions) comes
only with an oracle fixture. Until it is built, `union_pairs()` refuses two retrievals of one pair (§2.3), and
`compare(mode='prefix_subset')` shows whether one contains the other.

The two arms of one seed (disjoint arms) are well defined but not needed as an answer: a fetch with `direction:
both` serves it, and Appendix D.5 makes it a check. Different seeds are §4.3.

### 4.3 Two seeds of one locus are not one graphlet

Every position is outward from one seed boundary, every certificate is about walks from that seed, and the seed's
permitted set decides the walks in constrain mode. Placing seed B inside A's walks is an alignment, and seeds can
overlap A's walks more than once (repeats) or not at all. What U6 asks for, a k-mer graph of the union with labels
per node, is a different product (a local assembly graph). `to_gfa()` per graphlet plus a graph tool does it today.
Left out (OQ4; the owner's provisional answer agrees).

### 4.4 Neither result is an MGT body

A merged MGT v1 body would need one `H` (one `seed_index`, one identity), one `S`, an `L` dictionary whose first
`num_seed_labels` ids are the seed labels (the root's `G.entry *` rule), `A` counters and counts, `B` growth bins and
`V` branch events of one traversal, one `O` and at most one `Q`. Positions would also pass the arm's own certificate.
Summing counters or bins would invent values; dropping them would break "lossless" and the `*` rules. Across pairs
it would also need one index identity, which a union of pairs does not have. MGT v1 is frozen (SPEC §7.5). A
combination is therefore a **container of unchanged bodies** plus a stated relation (one seed across pairs; a walk
and its continuation), which is how views already work (`J.view`, the backing body unchanged).

## 5. Engine and service

### 5.1 The engine

**No engine change is needed or proposed.** Both combinations read only bodies and envelopes (the fan-out helper of §5.2
reads `GET /capabilities` for the pair list; the union itself reads no certificate). The union needs no engine feature:
inferred paths join graphs only where the parts already state a connection (§2.7), and a walk of one traversal over
several graphs is out of its reach. Two engine features would make continuations closer to exact. Both change the
`/traverse` request schema; the owner's provisional answer is "not now" (OQ7), and they are not part of this plan:

1. per-seed-label starting state in a seed (`seeds[].labels` as `{name, loss, branches}`), which removes the `switching`
   and `branches_unequal` causes and the `not_continued` gap for labels that switched in late. It is the more useful of
   the two;
2. an excluded-edge list (the ancestors' (k+1)-mers), which removes `rule_violation` at a request cost of one key per
   ancestor base.

Neither removes `sibling_merge_unknown` or `trace_restart`: those need the whole frontier or the followed
occurrence. A server-side frontier resume contradicts "the backend keeps no frontier between requests"
(graphlets.rst) and is not proposed. Under the project rule, an engine search-path change would go to master first,
as a PR Kulkov reviews.

### 5.2 MCP and the service

**No MCP exposure in phase 1.** The library's own tool functions (`mcp_tools.py`) gain nothing in 1a. A local tool
over a union or a stitched graphlet, returning a derived handle like `graphlet_subtrie`, would be phase 2 and only
if asked. The service's tools are the service session's design, and this note designs none. What the service
session would need to know:

1. **The library will offer `union_pairs()`** over one seed's records from several pairs of a release, or from
   several databases. The service's shape, one record per graph and seed with none merged, stays right: a union is
   built from those records, not stored in their place.
2. **What an honest union needs from the records**, per seed, all passed in one list (§2.3):
   - every record's body (`graphlet_export(format="mgt")`, then `load()`). The body must keep its `J` line with
     the envelope, or the union refuses it with MissingEnvelope; `{pair, graphlet}` can state a pair the envelope
     lacks, but it does not replace the envelope (`has_envelope` also needs the seed summary, `model.py:461`).
     The service session confirmed it for r237 (OQ21, §8.3): the service adds the `J` line to each stored record, and
     the REST `format=mgt` export and the MCP `graphlet_export(format="mgt")` both keep it. A multi-graph record carries
     its own pair's entry under a reduced envelope (that pair only; `envelope_version`, `route` and the other pairs
     dropped). The union reads the index fingerprints from the body's `H` line (`parser.py:228`), which the export
     copies unchanged (`parser.py:1365`), not from the envelope. A record from the master backend carries no pair
     fields: the union could name it only by its `index_ns`, which matches no graph name in `expected_pairs=`, so such
     records are passed as `{pair, graphlet}`;
   - every record without a body, with its pair and its seed result (`not_in_graph` with `kmers_present` and
     `kmers`, `labels_not_in_graph`, a failed derivation), as `{pair, result}`, with the seed's `sequence` where
     the service has it (otherwise the coverage says the record's seed was not checked; a seed id cannot be checked,
     §2.3); and a pair whose answer never came as `{pair, error}`. Without these, a graph that lacks the seed looks
     like a graph that was not asked;
   - the pairs the fetch was sent to, passed as `expected_pairs=`. Without them the union states its coverage "over
     the pairs given" only;
   - for a seed that no graph holds, the seed sequence as `seed=`, so the union can still say so.
3. **One orientation per union.** All records of one union share one oriented seed sequence. Records of one seed
   fetched in the other orientation are refused (`different_seed`, reverse-complement variant) until phase 2 flips
   them.
4. **Stitching in the service.** The service publishes no continuation tool today; what a stitched read in the
   service would need, how an agent can stitch locally, and a note on the wording of `graphlet_compare` are in
   Appendix F.5.

**For library users, one fan-out helper** (1a-U, small, no engine change):

```python
TraverseClient.traverse_pairs(seed, strategy=None, *, pairs=None, detail='graphlet', timing=True,
                              coordinates='auto', not_after_ms=None) -> PairResults
```

- *Dispatch.* One `/traverse` request per pair: `pairs` as graph names or `{graph, graph_path, annotation_path}`
  objects, by default every graph `GET /capabilities` lists.
- *Failures.* Every answer is kept as returned; a refused or failed request becomes `{pair, error}`.
- *Expected pairs.* `PairResults.expected_pairs` holds the full pair objects the requests were sent to, so two
  annotations of one graph name stay two pairs.
- *Seed and mappings.* `PairResults.seed` is the oriented seed sent, so a union whose every pair lacks the seed is still
  built (`seed=`, §2.3); `PairResults.seqs_order` holds each pair's record-mapping certificate from `GET /capabilities`
  (§2.5).
- *Kept seeds.* Each `TraverseResponse` keeps the seeds it sent (`seeds_sent`, beside the attempt fields it keeps as
  `.sent` today, `attempts.py:363`), so absent records read from it are `seed_check: 'sent'` (§2.3).
- *Not accepted.* `attempt_id` and `expect_server_instance`. Attempt ids share one server-wide namespace with request
  ids (envelope spec §1.1, §5), so one id sent to several pairs is refused with a 409 from the second request on. A
  ledger-managed fan-out, with a distinct id per pair and each receipt kept, is a later extension.

```python
>>> r = client.traverse_pairs(seed, strategy)                   # every graph of the release
>>> u = union_pairs(r.results, seed=r.seed, expected_pairs=r.expected_pairs, seqs_order=r.seqs_order)
```

## 6. Tests

Stdlib `unittest` with seeded random generation, as in the existing suites; no new dependency. Appendix E lists them
all: P1–P13 for the union and T1–T11 for stitching. Two of them are gates: the paths each exported GFA admits, checked
by enumerating them (P9), and, for switching stitching's `exact` on, a seeded randomized comparison against
uninterrupted traversals (T11; OQ19, decided). The oracle triples of T4 pin every cause, and they assert each walk's
`complete`, `certified_to_bp` and evidence beside the verdict.

## 7. Plan

### 7.1 Phases

- **Phase 0 (this note).** Four Fable reviews, an external review and a GPT review (Appendix G). The owner answered
  OQ1–OQ10, OQ13 and OQ14 provisionally and decided OQ18, OQ19 and OQ20, all on 2026-11-15 (§8.1); OQ12 and OQ15–OQ17
  are open (§8.2, §8.3); the service session answered OQ21.
- **Phase 1a-U, the cross-pair union, first.** `composite.py` and `union.py`, cut to what looking at one seed across
  r237's 11 graphs needs. Its first step is a round trip on one real r237 record: the service's `mgt` export, `load()`,
  `union_pairs()` (OQ21). Then: `union_pairs()` over one list of results, unique tags and retrieval identity, every
  refusal, the trie with terminal markers and stubs, spans, routes through merges and span labels, interval support,
  `absence`, coverage with the qualified and not-walked categories, labels per pair, walks, spell, carriers of a walk,
  FASTA and GFA with the statements, summary and `print()`, evidence and provenance, all in the source-supported mode
  only (OQ18, decided); per-pair warnings on header names whose record mapping was not checked, `unknown` where the
  certificate is not known (OQ20, decided, §2.5); the client's fan-out helper and its kept seeds (§5.2); P1–P12 and the
  integration check; the split-graph fixture; the docs. No new file format. One lane, its own review, then a merge train
  like every library change. **Why the union ships first.** It has a named user and a need today (the owner, on r237's
  11 graphs; OQ3). It states nothing beyond what its parts state, so no part of it waits for an oracle, while
  stitching's `exact` waits for a randomized comparison. And it builds the shared machinery (`composite.py`), which
  makes stitching smaller. **Cut from 1a-U** (Appendix G.2): `certificates=` and every name rule (only `seqs_order=`
  came back, for OQ20); `routes=` (replaced by stubs); a union-level `compare()` (each common pair's parts compare with
  their own `compare()`); `support_by_depth()`; `carriers(sequence=)`; a per-pair FASTA (it equals the parts' own
  exports); a `PairUnion` as an input; `walks(by='labels')`.
- **Phase 1a-S, stitching, after the union (the owner's provisional answer to OQ1).** `stitch()` and `StitchedGraphlet`
  with walks, spell, continuation and next_request delegation, FASTA, completeness, evidence, summary and provenance;
  every refusal; the junction and walk statement with every cause of §3.4, every verdict `qualified` (OQ19, decided);
  completeness, evidence and outcomes that carry every walk-scope cause (§3.5), so that no `qualified` verdict stands
  beside an unjustified `complete: true` or exact evidence; T1–T9 (T6 in its 1a form), with T4 and T5 asserting those
  fields and not only the verdicts; the CLI triples; the docs, marked experimental. A `StitchedGraphlet` is then
  accepted as a union part (§2.3).
- **Phase 1b (when a user needs label-level answers across junctions).** Stitched claims (constrain joined, annotate
  shifted), support profile and changes.
- **Phase 1c (when a user needs one file).** The bundle for both kinds (Appendix D.3) and T6's 1c form.
- **The inferred mode, a later phase of its own (OQ18, decided).** The inferred paths of §2.7 with their kinds, pieces
  and transitions, FASTA records and GFA links; P13 and the inferred cases of P9; after the source-supported union has
  been used on r237.
- **`exact` for stitching, a later step (OQ19, decided).** Gate: the seeded randomized comparison of T11 passes for
  W1–W8 as corrected. Then promotions to `exact` (`same_stop`; uniform-shift switching), each with its T4 fixture and
  T11 runs, and a `lower_bound` verdict under the conditions of Appendix G.1.
- **Later, a validation helper.** `validate_arms()` (Appendix D.5), when the owner wants it.
- **Phase 2 (each item only if wanted).** The union: the items cut from 1a-U; flipping a reverse-complement part on
  primary graphs; `next_request()` per pair (each part's own request for the walks it follows); a JSON view.
  Stitching: GFA; label walks, routes and label summary over stitched walks; annotate claims joined; trace
  contiguity from coordinates. Both: store entries and local tool functions if the service asks.
- **Later, the same-pair combination** (§4.2), when a user needs it; its exact case only with an oracle fixture.
- **Not planned.** Combinations of different seeds, any engine change.

### 7.2 Cost and size

- **Memory.** A union holds its parts' models: r237's 11 parts with two arms each take about 110 MB constrained, or
  440–660 MB label-free (about 5 MB per constrained arm, 20–30 MB label-free, DESIGN §7), plus a trie of 5–10 MB
  [estimate]. A stitched graphlet holds its parts' models plus about 1 MB per 10 kb of continued chain for the rule
  check [estimate]. `memory_bytes()` states both.
- **Work.** The union is built in one pass over each part's first-parent tree, linear in the parts' bases; routes,
  span labels and inferred paths are computed on demand. Stitching's decision check reads each child's recorded
  successors once per walk route, with the ancestors' (k+1)-mers in one hash set. Every operation of both
  combinations takes `budget=`.
- **Size** [estimate]: `composite.py` 300–400 lines; the union 900–1,300 lines with both path modes and about 850 lines
  of tests; stitching 1,200–1,800 lines and about 1,000 lines of tests; about 300 lines of docs. Appendix F has the
  details, the code plan and the docs plan.

### 7.3 Risks

Appendix F.6 lists sixteen risks with their mitigations. Two shaped the design most: a union read as one traversal of
the release (hence the five statements, a GFA topology per path mode and no `held_in` on inferred paths), and a stitched
walk read as one traversal (hence a whitelist that compares every decision of the child, certification that stops at the
first changed decision, and `qualified` only until a seeded randomized comparison passes).

## 8. Open decisions

### 8.1 The owner's answers so far

The owner: "do not take my answers as the final answers yet". Every answer below is **provisional, 2026-11-15**, except
OQ18, OQ19 and OQ20, which the owner **decided** on 2026-11-15; the body follows them all, and the provisional ones can
still change. The owner's own words are quoted where they were relayed.

| OQ | question | the owner's answer | status | where |
|---|---|---|---|---|
| OQ1 | stitching: a recipe, phase 1a, or phases 1a and 1b? | "Phase 1a": built as 1a-S, after the union | provisional, 2026-11-15 | §3, §7.1 |
| OQ2 | the union of the two arms of one seed | "lower priority, mostly useful for sanity checking/validation": a later validation helper | provisional, 2026-11-15 | Appendix D.5 |
| OQ3 | combining one seed across pairs | "why would I refuse it? It's the only way I can do it today. Referring to a non-existing and complex engine feature does not help", and "it won't be exact, I know, but the best we can do today, I think": offered, ships first, its limits stated with every answer | provisional, 2026-11-15 | §2 |
| OQ4 | a map of the locus from several seeds | left out | provisional, 2026-11-15 | §4.3 |
| OQ5 | stitching's names | `stitch`, `StitchedGraphlet` and `CombineRefused`; "merge" kept for reconvergence inside one traversal | provisional, 2026-11-15 | §4.1 |
| OQ6 | stitching without an index manifest on either side | allowed and marked unverified; the evidence then never reads `exact` | provisional, 2026-11-15 | §3.1 |
| OQ7 | an engine change for exact continuations | not now | provisional, 2026-11-15 | §5.1 |
| OQ8 | the union's names | `union_pairs`, `PairUnion`, `PairLabel`, `SpanLabel`, `PairSpan`, `Route`; the span values `whole`, `ends`, `merges`, `diverges`, `none`; the route values `held`, `ends`, `diverges`; the absence values `graph`, `under_strategy`, `uncertified`; `whole_in`, `held_in`; walk kinds `walk` and `stub`; the categories `complete_qualified` and `not_walked`; `every_pair_certified_to_bp`; labels shown as `<tag>/<name>`; the FASTA token `scope=within_one_graph` (names added since are OQ17) | provisional, 2026-11-15 | §2 |
| OQ9 | mixed settings across the graphs of one union (constrain with annotate, `kmer` with `trace`, different k) | allowed and stated: each pair's mode in the sentence, on each span and in `modes=`; no label total across pairs | provisional, 2026-11-15 | §2.3 |
| OQ10 | two answers on one pair | refused by the union; the same retrieval twice is a no-op; combining two retrievals is the later operation of OQ14 | provisional, 2026-11-15 | §2.3 |
| OQ13 | source-supported or inferred paths | "Expose both as separate modes": source-supported by default, inferred on request, each with a GFA topology that admits only its own paths; connections only where the input graphlets justify them, no sequence-overlap assembly | provisional, 2026-11-15 | §2.7, §2.9 |
| OQ14 | combining two retrievals of one pair | belongs in the design, implementation deferred; restated with OQ18 and OQ19: "combining separate retrievals from the same graph/annotation pair remains designed now and implemented later" | provisional, 2026-11-15 | §4.2 |
| OQ18 | when the inferred mode ships | "Ship source-supported union first. Keep inferred paths in the design and implement them in a later stage." | **decided, 2026-11-15** | §2.7, §7.1 |
| OQ19 | `exact` for stitching before a randomized comparison | "Ship stitching as qualified first. Defer exactness until the corrected rules pass a seeded randomized comparison against uninterrupted traversals. The qualified stitching release still needs the completeness fixes: walk-specific causes must propagate into `complete`, evidence and outcomes. A qualified verdict must not coexist with an unjustified `complete=True` or `exact=True`." | **decided, 2026-11-15** | §3.4, §3.5, §7.1, E.2 T11 |
| OQ20 | header names on `count_only` pairs | "add per-pair warnings for unchecked header-to-record mappings, using `unknown` when metadata is unavailable" (option (b); recorded from the owner's paste of the GPT follow-up review, for the owner to correct if it is not his decision) | **decided, 2026-11-15** | §2.5, §7.1 |

OQ11 (which rule joins label names across pairs) was withdrawn after the union review: under the server's
certificates no name occurs in two pairs, so no rule would ever join anything (§2.5).

### 8.2 Open: OQ12, two extras cut from 1a-U

The union review cut eight items from 1a-U (§7.1). Two of them were in the design brief, so they are the owner's to
confirm. Both can be added later without changing anything in 1a-U.

**(a) Support by depth across graphs**, `support_by_depth(arm, *, bin_bp=100)`. One table per arm: a row per depth
bin, a column per graph, and in each cell the walks alive at that depth (displayed or held through a route) and the
label entries carrying them, with a mark where the graph's certificate ends. For the example of §2.2, right arm,
every 100 bp (illustrative):

    depth (bp)   g01      g03       g04     g06       g07      g08
    0–100        5 / 9    4 / 1000  2 / 3   3 / 7     4 / 6    3 / 4
    100–200      4 / 7    3 / 980   2 / 3   3 / 6     3 / 5    2 / 2
    ...
    400–500      2 / 3    2 / 410   1 / 1   2 / 3     · cut    1 / 1
    ...
    900–1000     1 / 3    1 / 120   1 / 1   1 / 2     ·        —

    cell: walks alive / label entries; "· cut": the certificate ends in this bin; "·": beyond it; "—": no walk

- *What it adds:* where each graph's support thins out, on one screen, without reading walks; the depth at which a
  graph runs out of walks or of labels. Coverage gives one category per graph at one depth, and the spans give how
  far each graph follows one walk; the profile shows every depth at once.
- *What it costs:* one pass over the union segments per arm, O(segments × pairs); about 150 lines of code and
  tests [estimate]. Counts stay per graph, as label entries, with no total across graphs, because a record held by
  two graphs would count twice.
- *Without it:* `walks()` and the spans give the same numbers, binned by a few lines of user code.

**(b) Comparing two unions**, `union.compare(other)`. For example, the same seed on r237 and on a later release:
the graphs in both and in one only; per graph, both coverage statements side by side; the union walks spelled in
one only (with `mode='prefix_subset'`, a walk counts as contained when it is a prefix of a walk of the other); and,
for graphs of one index on both sides, each part's own `compare()`.

- *What it adds:* one call answers "which flanks appeared or vanished, and in which graph". Graphs of two releases
  have different indexes, so they are matched by name and the match is stated as unverified; their parts' own
  `compare()` would answer `comparable: false`, and the union-level comparison by spelling is what remains.
- *What it costs:* about 300 lines of code and tests [estimate]; the matching of graphs across releases needs its
  own rules and tests.
- *Without it:* two unions are built, and their walks are compared by spelling in a few lines of user code; two
  retrievals of one graph on one index compare with `compare()` today.

**Recommendation:** leave both out of 1a-U and add either when the owner wants it. Each is additive and changes nothing
that 1a-U states. Support by depth would be the natural first follow-up; comparing unions needs rules of its own for
matching graphs across releases. *Reviewers:* Fable and GPT recommend the same.

### 8.3 Open: three questions from the reviews, and OQ21 answered

Each question stands on its own; my recommendation follows it.

- **OQ15, where an inferred path may switch graphs.** In the first version, only where one graph's own answer rejoins
  two of its walks (a merge in that graph), and at the seed, where one graph's left flank can meet another's right flank
  in the GFA. No joining by overlapping sequence. *Recommendation: accept; joins by overlap only with a design of their
  own. Fable and GPT agree.*
- **OQ16, the default GFA writes the seed twice.** So that the file admits no path across graphs, the default-mode GFA
  writes the seed once per arm, and a graph's walk through its own merge is a tag, not a link. One graph's left flank,
  seed and right flank are then not one path of the file; each part's own GFA has that path. *Recommendation: accept;
  the inferred-mode GFA joins both, and the file says so in a comment line (§2.9). Fable and GPT agree.*
- **OQ17, the names added since OQ8.** `paths='supported' | 'inferred'`, `max_transitions`, `InferredPath` (ids
  `i:<segments>`), `Transition` (kinds `merge`, `seed`), the span value `extends`, the absence value `walk_rule`,
  `LabelInterval`, `Statements`, `Carriers`, `seed_check: 'sent'`, `TraverseClient.traverse_pairs()`, the stitch causes
  `candidate_changed` and `decisions_unrecorded`, and the later `union_retrievals()` and `RetrievalUnion`. The earlier
  `pair_id` is gone: unique display tags are the one identifier. Since the GPT review: the path kinds `walk`, `route`,
  `joined`, `max_paths`, `Transition.by` as a tuple, and `PairResults`. *Recommendation: accept. Fable and GPT agree,
  GPT after the contract fixes of Appendix G.6, which are made.*

**OQ21, the `J` line from the service: answered by the service session, 2026-11-15, for r237.** The service adds the `J`
line to every stored record, and both the REST `format=mgt` export and the MCP `graphlet_export(format="mgt")` keep it.
A multi-graph record carries its own pair's entry under a reduced envelope; the union reads the index fingerprints from
the record's `H` line, which the export copies unchanged. Records from the master backend carry no pair fields and are
passed as `{pair, graphlet}`. Nothing here was checked against a live exported record, so the first step of 1a-U is a
real round trip on one r237 record: export, `load()`, `union_pairs()` (§5.2, §7.1).

## Appendices

The body refers to these by letter. They hold the exact texts and the detail behind the body's statements; nothing
in them changes a decision of the body.

## Appendix A. The union: inputs, refusals and all operations

### A.1 How results are sorted

The sorting is total: every element becomes a part, an absent record, or a refusal.

| element | becomes |
|---|---|
| a `Graphlet` | a part, its pair from the envelope's tags (`pair_stated_by: envelope`); a view is refused `view_unsupported`, a body without its envelope MissingEnvelope |
| a `StitchedGraphlet` (alone or as `{'pair', 'graphlet'}`) | a part, as its root would be: its pair, envelope, seed and certificate are the root's; every part of it carries bases (else `no_bases`); its digest is its body and retrieval identity (§2.3) |
| a `TraverseResponse` | with exactly one seed result (else `not_one_seed`): its graphlet is a part; a result without a graphlet is read as `{result}` below, its pair from the response's envelope (`pair_stated_by: envelope`) |
| `{'response': TraverseResponse, 'index': int}` | seed result `index` of a response with any number of seeds, read as the response alone would read a one-seed response: its pair still from the envelope (`pair_stated_by: envelope`), and `seed_check: 'sent'` where the client kept the seeds it sent; an index out of range is `not_one_seed` |
| `{'response': TraverseResponse, 'sequence': str}` | as the response alone, with the caller's statement of the seed (R3) |
| `{'pair', 'graphlet'}` | a part on the pair the caller states (`pair_stated_by: caller`), refused `pair_mismatch` if the envelope states a pair that does not match it. The body still needs its `J` line: `has_envelope` needs the envelope and the seed summary (`model.py:461`) |
| `{'pair', 'result'}`, optionally with `'error'` and `'sequence'` | the result decides; an `'error'` beside it is the result's own text, as in `TraverseResponse.errors` entries (`client.py:909`). An absent record only when the result has no `graphlet`, no `arms`, and `outcome.walks` is `not_in_graph`, `labels_not_in_graph` or `failed`; a result with a graphlet is refused `result_has_graphlet`, and a walked result without one (`outcome.walks` `complete` or `partial`: fetched with another `detail`, `client.py:836`) `no_graphlet` |
| `{'pair', 'error'}` | a `failed` record: the pair's answer was refused or never came (`single_answer()` raises on a refused entry, `client.py:276`, and an HTTP failure raises). `error` is a string or an exception (its `str()`) |
| anything else | `not_combinable`: another type, a dict without `pair` (or `response`), with two of `graphlet`, `result` and `error` other than the pair above, or with unknown keys |

### A.2 Refusal texts

| code | text |
|---|---|
| `bad_expected` | `expected_pairs= takes a list of the pairs the request was sent to, as graph names or {graph, graph_path, annotation_path} objects; %s is not one` (a string: `expected_pairs= is a string: pass a list of graph names`) |
| `bad_seqs_order` | `seqs_order= maps pairs (graph names or pair objects) to "verified", "count_only" or None, as GET /capabilities states them; %s is not one` |
| `not_combinable` (TypeError) | `union_pairs() takes results: a Graphlet, a StitchedGraphlet, a TraverseResponse, {response, sequence}, or {pair, graphlet}, {pair, result} or {pair, error}; not %s` |
| `view_unsupported` | `union_pairs() reads whole retrievals, not views: pass view.backing, and select labels on the union (labels=)` |
| `not_one_seed` | `result %s holds %d seed results: pass {response, index} for this union's seed, or it as {pair, graphlet} or {pair, result}`; zero: `result %s holds no seed result`; out of range: `result %s has no seed result %d (it holds %d)` |
| `result_has_graphlet` | `result %s carries a graphlet: pass it as a Graphlet (Graphlet.from_response(result, response)) or as {pair, graphlet}` |
| `no_graphlet` | `result %s was walked (outcome.walks %s) but carries no graphlet: fetch with detail='graphlet'` |
| `no_seed` | `no result carries a graphlet, so the union has no seed sequence: pass seed= (the oriented seed the request sent)` |
| MissingEnvelope | today's `require_envelope('union_pairs()')` text, as for `stitch()` (§3.1) |
| `pair_mismatch` | `result %s was answered on pair %s, and the caller states pair %s: state the pair it was answered on, or none` |
| `no_pair` | `result %s names no pair (no graph tag, no index_ns): pass it as {pair, graphlet} or {pair, result}` |
| `ambiguous_pair` | `pair %s matches both %s and %s: state its paths, release or index_ns` |
| `no_bases` | `part %s carries no bases (output.sequences: false): a union by sequence needs every part's bases` |
| `different_alphabet` | `part %s walks the alphabet %s and part %s the alphabet %s: their bases are not comparable` (compared without `$`: succinct `$ACGT` and hash `ACGT` are one DNA alphabet) |
| `different_seed` | `part %s's seed is not part %s's (they first differ at base %d): a union by sequence needs one oriented seed, spelled base for base the same`; unequal lengths: `part %s's seed is %d bp and part %s's %d bp (equal for the first %d bases): a union by sequence needs one oriented seed`; reverse complement: `part %s's seed is the reverse complement of part %s's: a union by sequence needs one orientation; fetch part %s with the seed as part %s spells it`; `seed=`, an absent record's `length_bp` or a stated `sequence` against the parts' seed: the same texts with `seed=` or `the absent record of pair %s` in place of a part |
| `duplicate_pair` | `results %s and %s were answered on one pair (%s): a union holds one answer per pair; compare() reads two retrievals of one pair, and union_retrievals() will combine them (not yet offered)`; one body from two requests: `results %s and %s are byte-identical bodies on pair %s from different requests (%s differs): equal bodies are not equal retrievals, and a union holds one answer per pair`; when they cannot be told apart: `results %s and %s cannot be told apart (no index_fp on %s, and their names, paths, release, index_ns and index_meta_fp agree): they may be one pair; state each pair with {pair, graphlet} or {pair, result}` |
| `not_expected` | `pair %s (a part or an absent record) is not among expected_pairs=` |
| `ambiguous_name` | `label name %r is in pairs %s: equal names in two pairs are not claimed to be one record; select by {pair, name} or {pair, ref}`; a selector `{'id': n}` is BadSelector, `label ids are local to each part of a pair union: select by {pair, ref} or {pair, name}`; `{pair, name}` naming two labels of that pair is the part's own AmbiguousLabel, `%d labels named %r (%s): select by ref` (`ops.py:297`) |
| `summary_too_small` (ValueError) | `the statements alone take %d bytes, more than max_bytes=%d: they are never cut; pass max_bytes >= %d` |
| `not_offered` | `dump()`: `a pair union is not an MGT document: each part dumps on its own (parts[i].graphlet.dump())`; `to_json()`: `to_json() is detail "full" of one server result, and a pair union has none: export its parts`; any other: `%s() is not offered on a pair union` |

### A.3 All operations: preconditions and refusal codes

| operation | status | preconditions (all checked, in the order of the refusal tables) | refusal codes |
|---|---|---|---|
| cross-pair union `union_pairs(results)` | **offered**, phase 1a-U, **first** | results sorted into parts and absent records; every part a whole retrieval with an envelope and bases; one alphabet; one oriented seed, base for base; one answer per pair; expected pairs well formed (§2.3). Allowed and stated: different k, regime, release, database, strategies, radii, arm sets, unverified identity, caller-stated pairs | `bad_expected`, `not_combinable`, `view_unsupported`, `not_one_seed`, `result_has_graphlet`, `no_graphlet`, `no_seed`, MissingEnvelope, `pair_mismatch`, `no_pair`, `ambiguous_pair`, `no_bases`, `different_alphabet`, `different_seed`, `duplicate_pair`, `not_expected`; for queries `ambiguous_name`, `summary_too_small`, `not_offered` |
| continuation stitching `stitch(parent, child)` | **offered**, phase 1a-S, after the union (OQ1) | both carry an envelope; same index (§3.1: `index_meta_fp`, `index_fp`, `index_ns`, release, names per ref; one-sided `index_fp` only with `allow_unverified_index=True`); same pair where both state one; same k, regime and alphabet; same mode and support; the child walks one arm; its seed equals the continuation of one walk of that arm of the parent; that walk is not yet continued by another body; in constrain mode its seed labels are alive at the walk's end | `not_combinable`, `view_unsupported`, MissingEnvelope, `different_index`, `index_unverifiable`, `different_pair`, `different_graph`, `different_mode`, `not_one_arm`, `wrong_arm`, `no_continuation`, `not_a_continuation`, `ambiguous_walk`, `already_continued`, `foreign_seed_labels` |
| the two arms of one seed against a `direction: both` fetch | **later**, a validation helper (Appendix D.5) | as stitch's identity checks; one oriented seed | none designed yet |
| two retrievals of one pair `union_retrievals(results)` | **later**: designed (§4.2), built when a user needs it | one oriented seed; retrievals told apart by retrieval identity; labels joined by ref only where `index_fp` agree | until then `union_pairs()` refuses them `duplicate_pair`, and `compare()` reads them |
| different seeds (overlapping or not) | **left out** (a different product, §4.3) | n/a | `different_seed`, `not_a_continuation` |

### A.4 Pairs, tags and retrieval identity

**The pair key** is (release, `index_ns`, graph, graph_path, annotation_path, `index_fp`, `index_meta_fp`), the stated
pair taking the place of the envelope's tags where the caller states one, `None` first. It orders parts and absent
records (byte order). A pair described both by a bare name and in full keeps the fuller description.

**When two descriptions name one pair.** One rule decides `duplicate_pair`, `pair_mismatch`, `not_expected`,
`not_included` and whether an absent record and a part are one pair:

1. where both state an `index_meta_fp` and they differ, they are two pairs: a different graph or annotation build is
   proof, as `_identity` treats it (`ops.py:3593`);
2. where both state an `index_fp`, the `index_fp` decides;
3. otherwise the graph names must be equal, and `graph_path`, `annotation_path`, `release`, `index_ns` and
   `index_meta_fp` must be equal wherever both state them (a bare name states only the name).

Because a bare name can match two different pairs, a description that matches two distinct pairs of the union is
refused `ambiguous_pair`, so the outcome never depends on which comparison runs first.

**Tags.** A pair's display tag is the pair the caller stated, else its graph name on a multi-graph server, else its
`index_ns` on a single-graph one; qualified as `<release>/<graph>` when the parts span releases. Where two tags still
collide (one graph name and release, different paths, no `index_fp`), each gets the suffix `~` and the shortest
prefix, at least 6 hex digits, of the sha256 of its pair key that tells them apart. Tags are therefore unique in the
union, and they are the one identifier: every dict, selector, span, label (`<tag>/<name>`), evidence block and
provenance entry uses them, and each part keeps its pair key beside its tag (`PairPart.key`). A tag is stable for a
given set of pairs; it gains a suffix only when a colliding pair joins. A result with no name at all is refused
`no_pair`.

**Retrieval identity** is (the pair key; the index identity: `index_fp`, `index_meta_fp`, `index_ns`, release; the
normalized strategy from the envelope, all of it; the seed as requested: its case-mapped sequence, orientation and
labels). It leaves out the caller's `seed_id`, and it does not depend on `GraphletStore.request_of()`, which returns
`None` without an envelope and where the seed labels' names cannot be verified to resolve back, and adds the caller's
`seed_id` (`store.py:1289-1301`); so it is defined for every part with an envelope. An absent record's is (the pair key,
its result as canonical JSON without `seed_id`); a failed record's, (the pair key, the error text).

**No-ops and replacements.**

- The same body with the same retrieval identity on a matching pair is a no-op. When one retrieval comes both as a
  `Graphlet` and as `{pair, graphlet}`, the envelope-stated form is kept, whatever the order; of two envelope-stated
  copies, the one whose envelope's canonical JSON sorts first. Each part keeps the model it was given (results arrive
  parsed, from the client or `load()`): the union copies no model and merges none with an equal body. Sharing one parsed
  body between retrievals is a later optimization, because a model carries its envelope and caches beside the body; the
  store keeps body and entry identity apart for the same reason (`store.py:3-7`).
- Two identical absent records are one.
- A failed record (`{pair, error}`) and an answer for the same pair: the answer (a part or an absent record) is kept,
  and the failure is listed in `provenance()` as an earlier attempt. Of two failed records for one pair, the first in
  canonical order is kept and the other listed.
- Anything else on one pair (two bodies, or one body from two requests) is refused `duplicate_pair`.
- Two pairs with byte-identical bodies (possible without manifests, once the caller states the pairs) stay two parts.

### A.5 Export tokens and escaping

**FASTA header tokens** (`key=value`, empty ones left out):

- `kind=` `walk` or `stub`, and in inferred mode also `route` or `joined` (§2.7); `len=` the length in bp.
- `held=` the pairs that hold the whole walk (displayed or through a route); `whole=` those that display it exactly;
  `extends=` those whose displayed walk or stub goes on past it.
- A span token, under `ends=` or `diverges=`, is `<tag>@<to_bp>[:<end reason>]/<absence>`.
- A `merges` token is `<tag>@<to_bp>><trie node>` (`s31`, or `s30$` for a merge segment without bases), followed by
  `~<route to_bp>[:<end reason>]/<absence>` when the route stops short.
- `certified=` each pair's `certified_to_bp`, as `<tag>:<bp>`.
- `junctions=` the junctions of a stitched part that the walk's span reads, as
  `<tag>@<at_bp>:<verdict>[:<code>|<code>...]`, the codes those of `SpanJunction.causes` (`[a-z_]+`, so `|` needs
  no escape): `junctions=A@115:qualified:exact_gated|candidate_changed`; every verdict `qualified` in this phase,
  §2.3.
- `labels=` label entries per pair as `<alive>/<carrying the whole walk>`, `+` marking a lower bound (cut lists). The
  full lists with their `from_bp` are `carriers()`'s; `to_fasta()` writes the counts only.
- `modes=` whenever a part is not in annotate mode (`modes=g1:constrain,g7:annotate` when they differ).
- `absent=` the absent records as `<tag>:<category>`; `not_included=` the expected pairs not given.
- `mapping=` the pairs whose header labels rest on an unchecked record mapping, as `<tag>:count_only` or
  `<tag>:unknown`, or one value when every pair has it (`mapping=count_only` on r237).
- `scope=within_one_graph` (the non-claim of the default mode) or `scope=inferred`.
- Stub walks are exported only with `stubs=True`, as `kind=stub`.
- An inferred path states `pieces=` as `<tag>[+<tag>…]:<from_bp>-<to_bp>` and `transitions=` as
  `merge:<tag>[+<tag>…]@<trie node>:<at_bp>[<from_bp>-<to_bp>]` (every pair whose merge permits it) or
  `seed[<from_bp>-<to_bp>]`, instead of span tokens; a route also states `held=`, as a walk does; every record of an
  inferred listing states `of=<count>`, `of=<count>+` when the listing was cut. Its `certified=` gives each pair's
  certified depth on its piece.

**Escaping.** Names and tags in tokens are percent-escaped, so that every token decodes back to what it names: every
byte outside printable ASCII (`!` to `~`), which takes in the space, control characters and the UTF-8 bytes of every
non-ASCII character, and `%` and the token characters `,;=@/:>~+|[]$`, are written as `%XX`. GFA1 allows only printable
ASCII in names and `Z` tag values, so a label `café` is written `caf%C3%A9`; FASTA headers use the same rule, so one
decoder reads both. The MGT codec escapes only `%` and line breaks (`_codec.py:263`); a header token needs more. GFA tag
values use the same escaping.

## Appendix B. Stitching: refusal texts and junction causes

### B.1 Refusal texts

`%s` and `%d` are filled as stated; a label list prints at most 8 `{name, ref}` then "and N more", as `next_request`'s
notes do.

| code | text |
|---|---|
| `not_combinable` (TypeError) | `stitch() combines graphlets (Graphlet or StitchedGraphlet), not %s` |
| `view_unsupported` | `stitch() reads whole retrievals, not views: stitch view.backing instead` |
| MissingEnvelope | today's `require_envelope('stitch()')` text: `stitch() needs the response envelope (seed_id, the normalized strategy, annotation counters): build the graphlet with Graphlet.from_response() or load a saved .mgt file with its J line; a parsed body alone does not carry it` |
| `different_index` | `the child was retrieved from another index than its parent (%s): a continuation runs on its parent's index, so a child of another index continues none of its walks; union_pairs() combines one seed's graphlets across pairs`, `%s` one of `index_meta_fp differs`, `index_fp differs`, `index namespace %r vs %r`, `release %r vs %r`, `label %s is %r in one and %r in the other` |
| `index_unverifiable` | `the child states %s index manifest digest (index_fp) and its parent %s: whether they are the same index cannot be verified (not proven different), so the child is not stitched automatically; pass allow_unverified_index=True to stitch it anyway (the result then states the identity unverified)` (`an`/`no`, `one`/`none`, as `_same_index`, `mcp_tools.py:2578`); where a side's parts state it unevenly: `the parts state an index manifest digest (index_fp) unevenly (the parent's on %d of %d parts, the child's on %d of %d): whether they are the same index cannot be verified (not proven different), so the child is not stitched automatically; pass allow_unverified_index=True to stitch it anyway (the result then states the identity unverified)`. Both, and `different_index` and `different_pair`, are checked against every part on both sides |
| `different_pair` | `the child was answered on pair %s, its parent on pair %s: a continuation runs on its parent's pair` (each `{graph, graph_path, annotation_path}` as compact JSON) |
| `different_graph` | `the child's graph is not its parent's (%s %r vs %r): its walks cannot continue the parent's` (`k`, `regime` or `alphabet`) |
| `different_mode` | `the child was walked in %s mode with %s support, its parent in %s mode with %s support: its walks answer another question` |
| `not_one_arm` | `the child walked both arms: a continuation walks only the arm it continues (next_request() sets direction to it)` |
| `wrong_arm` | `the child walked the %s arm, not the %s arm named by arm=` |
| `no_continuation` | today's `_continuation` text: `walk %d of the %s arm has no continuation (it ended for a semantic reason, not at the radius or a cap)` |
| `not_a_continuation` | with `walk=`: `the child's seed is not the continuation of walk %s of the %s arm: a stitched child is seeded with that walk's continuation (continuation(), next_request())`; a `walk=` the arm does not have (out of range, or below a walk not continued): `walk= names no walk of the %s arm (%s): a stitched child is seeded with a walk's continuation (continuation(), next_request())`; without: `the child's seed is the continuation of no walk of the %s arm: a stitched child is seeded with a walk's continuation (continuation(), next_request())` |
| `ambiguous_walk` | `the child's seed is the continuation of walks %s of the %s arm: name the one it continues with walk=` |
| `already_continued` | `walk %s of the %s arm is already continued by another child (part %s, digest %s): one child per walk` (the digest of the subtree there, as `stitch()` of it alone states it; a plain graphlet is compared with the part at the walk only, a stitched child with the whole subtree there); the same bodies, a part's from another request: `walk %s of the %s arm is continued by this body from another request (%s differs): equal bodies are not equal retrievals, and the two would continue differently` (`%s` a field path, `… of part %s` for a part below the child's root). An identical subtree (every part's key, body and retrieval identity) is a no-op |
| `foreign_seed_labels` | `the child seeds %s, not alive at the end of walk %s of the %s arm: their runs would start new lineages at %d bp, not continue the walk's` |
| `summary_too_small` | `the summary's statements alone take %d bytes, more than max_bytes=%d: they are never cut; pass max_bytes >= %d` (`summary()` cuts junctions from the end first) |
| `not_offered` | `dump()`: `a stitched graphlet is not an MGT document: save() writes its parts unchanged in a bundle, and each part dumps on its own (parts[i].graphlet.dump())`; `to_json()`: `to_json() is detail "full" of one server result, and a stitched graphlet has none: export its parts`; `compare()`: `compare() reads one retrieval per side: compare the parts`; any other: `%s() is not offered on a stitched graphlet` |

### B.2 The junction causes

**The causes.** A condition that fails is stated as a cause, and the verdict is `qualified`. Each cause has a scope
(the junction, or one walk) and the modes where it can occur.

| cause | fails | scope | modes |
|---|---|---|---|
| `trace_restart` | W1: under `trace` a lineage follows one occurrence chain from the root seed (SPEC §6.3). The child validates its seed label against any coordinate-consecutive occurrence of the tail (SPEC §6.1), so it can follow an occurrence the parent's lineage never had | junction | both, where `trace` is accepted |
| `strategy_differs` | W2, naming the first differing field path in byte order | junction | both |
| `extra_foreign` | W2: the child's `labels.extra` names a label outside the parent's pool | junction | constrain |
| `switching` | W3 | junction | constrain |
| `branches_unequal` | W4 | junction | constrain |
| `label_count_knob` | W5, naming the knob | junction | both |
| `sibling_merge_unknown` | W6, naming the sibling walks | junction | both |
| `not_continued` | W7, naming the labels | junction | constrain |
| `seed_dropped` | W7, naming the labels and the `X` reasons | junction | constrain |
| `no_bases` | W8 cannot be checked: a part from the root to the child has `output.sequences: false`, named (`on the route from the root to this child, part right:0 carries no bases …`); W8 is not run below the junction | junction | both |
| `strategy_unknown` | W2–W5 or the split check of W8 cannot be judged: an envelope does not state a field they read (W5's knobs; in constrain mode also `labels.change_cost`, `labels.extra`, `branching.max_label_branches`, `branching.max_splits_per_path`), one part states it and the other does not, or a field W2 compares neither states (a null field is not stated); naming them (`(both)`, `(the parent)`, `(the child)`) | junction | both |
| `rule_violation` | W8: the walk takes a step in B (it reuses an ancestor's edge or enters a root-seed node); stated with `at_bp`. The walk is not cut and no end is invented | walk | both |
| `candidate_changed` | W8: a successor the child admitted at a node of the walk's route, without taking it, is a step in B: a branch event's successor or another child of a split. The child's ambiguity, split or refusal there may not happen in the reference (the external review's counterexample: the child stops at 0 bp, the reference goes on to 145 bp); stated with `at_bp` and the character | walk | both |
| `exact_gated` | at a junction: W1–W7 hold there, and W8 is judged per walk below it; on a walk: its junctions are `exact_gated` and W8 holds for it. `exact` is not switched on (phase 1a-S, OQ19): the rules found no difference and have not yet passed the randomized comparison (Appendix E.2, T11). `summary()` and `completeness()` state, per junction, the walk-scope causes on the walks below it (`walk_causes`) | junction, walk | both |
| `decisions_unrecorded` | W8 cannot be checked at or beyond the child's `evidence.complete_to_bp` (branch events cut by `output.max_branch_events`; SPEC §7.0 guarantees a recorded reason only strictly below it), and the walk's end lies there, a 0 bp walk with a boundary of 0 included; stated with that depth | walk | both |
| `seed_reentry` | W8: a child `rejoined_seed` on a node of the child's seed outside the root seed. Two sub-cases, told apart by the blocked event's character (`E … b`): `same_stop`, where the blocked (k+1)-mer is an ancestor's edge, so the reference stops there too (with `edge_reuse`, by the precedence of SPEC §6.6) and only the reason differs; `earlier_stop`, a new edge into a node of the walk, where the reference goes on. Common: with the whole walk as the child's seed, any repeat of the walk ends the child here (§3.2) | walk | both |
| `split_limit_restarted` | W8: the parent's split count at w's leaf (the maximum over a merge's parents, as the engine counts) plus the child's count after a split it made on the walk's route exceeds `max_splits_per_path` (the reference would have refused that split); stated with the position of that split | walk | both |

Every walk-scope cause also ends certification at its position (§3.5): a stitched walk is certified at most to the first
one on its route. `same_stop` is the first candidate for promotion: its walks are the reference's and only the end
reason differs, so once T4 pins it, it can be `exact` with the reference's reason restated, after the randomized
comparison passes. Until `exact` is switched on (§7.1), every verdict is `qualified`. The verdict never depends on notes
`next_request()` returned: a child made another way (a service job with preset and overrides, a request written by hand)
is judged by what it ran. Hairpins (`hairpins: skip | follow`) are decided on one (k+1)-mer and behave the same given
the same knob, which W2 requires. `tip_window_bp` and `bubble_window_bp` are rejected unless 0 (SPEC §6.5), so W2 keeps
them equal; if they are ever implemented, a non-zero window becomes a `qualified` cause.

## Appendix C. Data classes and the algebra

### C.1 The union

```python
@dataclass(frozen=True)
class PairPart:    tag: str; key: tuple; pair_stated_by: str         # tag: unique (A.4); 'envelope' | 'caller'
                   seqs_order: str                                  # 'verified' | 'count_only' | 'none' | 'unknown'
                   graphlet: Graphlet | StitchedGraphlet            # a stitched part: §2.3
                   identity: dict; retrieval: tuple; digest: str    # retrieval identity: A.4
@dataclass(frozen=True)
class PairAbsent:  tag: str; key: tuple; pair_stated_by: str; category: str; detail: dict; seqs_order: str
                   seed_check: str                                  # 'sent' | 'stated' | 'unchecked' (§2.3)
                   # category: 'seed_not_in_graph' | 'labels_not_in_graph' | 'not_walked' | 'failed'
@dataclass(frozen=True)
class PairLabel:   pair: str; ref: str; name: str; kind: str        # pair: a tag; str(): '<tag>/<name>'
@dataclass(frozen=True)
class SpanLabel:   label: PairLabel; from_bp: int                   # from_bp 0: carries the whole span
@dataclass(frozen=True)
class LabelInterval: label: PairLabel; from_bp: int; to_bp: int     # where a label is alive within a segment
@dataclass(slots=True)
class SegmentSupport:  pair: str; n_walks: int; n_stubs: int        # the walks and stubs passing, by count
                       pieces: tuple[tuple[int, int, int], ...]     # (part segment id, from_bp, to_bp)
                       labels: tuple[LabelInterval, ...]; labels_cut: bool   # from alive sets and runs, in order
                       part: tuple = ()                             # a stitched part's sub-part (Part.key)
@dataclass(slots=True)
class UnionSegment:    id: int; arm: str; from_bp: int; to_bp: int; parent: int | None; children: tuple[int, ...]
                       support: tuple[SegmentSupport, ...]          # one per pair spelling it, canonical order
@dataclass(slots=True)
class Route:       to_bp: int; how: str                             # 'held' | 'ends' | 'diverges'
                   goes_on: bool; end_reasons: dict; absence: str | None
                   merges: tuple[tuple[int, str], ...]              # (at_bp, trie node) of each merge crossed
                   labels: tuple[SpanLabel, ...]; complete: bool; certified_to_bp: int
@dataclass(slots=True)
class PairSpan:    pair: str; mode: str; how: str; to_bp: int       # whole | extends | ends | merges | diverges | none
                   walk: int | tuple | None                         # the part's path id (whole, extends, ends);
                                                                    # past a junction the stitched walk's key
                   stub_end: int | None; into: str | None           # merges: the part segment; trie node 's31' | 's30$'
                   route: Route | None                              # merges only, computed on demand
                   continues: tuple[PairLabel, ...]                 # merges: kept entries entering through it
                   end_reasons: dict; labels: tuple[SpanLabel, ...]; labels_cut: bool
                   complete: bool; certified_to_bp: int             # §2.4: defined for every how
                   absence: str | None                              # graph | walk_rule | under_strategy | uncertified;
                                                                    # None for whole, extends, merges
                   k: int; regime: str
                   junctions: tuple[SpanJunction, ...] = ()         # a stitched part's junctions the span reads
@dataclass(frozen=True)
class SpanJunction: at_bp: int; verdict: str; part: tuple; causes: tuple[str, ...]   # Part.key; the junction's
                                                                    # cause codes, then the walk-scope causes of
                                                                    # the child walk read (§2.3, B.2)
class Carriers(list):                                               # carriers(): a list of PairSpan
    def __str__(self) -> str: ...                                   # one line per pair, as in §2.2
@dataclass(slots=True)
class UnionWalk:   id: int; arm: str; kind: str                     # 'walk' | 'stub'; ends at a terminal marker
                   length_bp: int; sequence: str | None; segments: tuple[int, ...]
                   spans: tuple[PairSpan, ...]                      # one per pair, canonical order
                   whole_in: tuple[str, ...]                        # the pairs with a 'whole' span (displayed)
                   held_in: tuple[str, ...]                         # whole, extends, and routes held to the end
@dataclass(frozen=True)
class Transition:  kind: str; by: tuple[str, ...]; at: str | None   # 'merge' | 'seed'; every permitting pair; trie node
                   at_bp: int                                       # the merge's position (0 for the seed)
                   from_bp: int; to_bp: int                         # the after side holds from, the before side to
                   before: tuple[str, ...]; after: tuple[str, ...]  # tags holding each side
@dataclass(slots=True)
class InferredPath: id: str; arm: str; length_bp: int; sequence: str | None   # id 'i:' + its union segments
                   kind: str                                        # 'walk' | 'route' | 'joined' (§2.7)
                   held_in: tuple[str, ...]                         # walk, route: the pairs holding it whole
                   segments: tuple[int, ...]                        # union segments in path order
                   pieces: tuple[tuple[int, int, tuple[str, ...]], ...]   # cut wherever the holding set changes
                   transitions: tuple[Transition, ...]              # empty for a walk or a route
                   certified: dict                                  # tag -> certified_to_bp on its stretch
class InferredListing(list):                                        # [inferred] walks(paths='inferred'), §2.7
                   count: int; cut: bool; resume: str | None        # enumerated; stopped by max_paths or budget
@dataclass(frozen=True)
class Coverage:    arm: str; at_bp: int; per_pair: dict; by_category: dict   # at_bp: the depth used (§2.6)
                   every_pair_certified_to_bp: int | None; sentence: str     # None also without expected_pairs=
@dataclass(frozen=True)
class Statements:  arm: str | None; at_bp: int; paths: str           # 'supported' | 'inferred'
                   coverage: Coverage; provenance: str; turn_offs: str; labels: str; scope: str
                   def sentences(self) -> list[str]: ...            # the five, in order, as print() writes them
@dataclass(frozen=True)
class PairResults: results: list; seed: str                         # traverse_pairs() (§5.2): the seed sent
                   expected_pairs: tuple[dict, ...]                 # the pair objects sent
                   seqs_order: dict                                 # pair -> its record-mapping certificate (§2.5)
```

### C.2 Stitching

```python
class StitchedGraphlet:
    parts -> tuple[Part]                  # root first, then by part key (§3.7)
    junctions -> tuple[Junction]          # by (arm, walk key)
    k, regime, alphabet, mode, support    # the root's (equal in every part)
    seed -> SeedInfo                      # the root's
    labels -> tuple[Label]                # union by ref: the root's order, then labels new in later parts
    identity -> dict                      # {verified, index_fp | note, accepted?} as traverse_continue states it
    arms() -> list[str]
    stitch(child, *, arm=None, walk=None, allow_unverified_index=False, budget=None) -> StitchedGraphlet
    walks(arm, *, top=None, by='support'|'length'|'loss'|'id', labels=None, route_consistent=True, min_bp=0,
          budget=None, resume=None) -> list[StitchedWalk]
    spell(arm, walk, orientation='natural'|'walk', with_seed=False, *, budget=None) -> str
    claims(arm=None, labels=None, at_most_bp=None, strict=True, *, budget=None, resume=None)
          -> list[StitchedClaim]                                  # [1b]
    support_profile(arm, walk, kind='displayed'|'route', *, budget=None) -> list[SupportRun]   # [1b]
    support_changes(arm, walk, *, budget=None) -> list[Change]   # [1b] a junction change has reason 'junction'
    continuation(arm, walk, *, budget=None) -> Continuation       # the owning part's, seed_coord mapped (§3.2)
    next_request(arm, walks, bp=None, reduce_budget=True, reset_branches=False, **overrides) -> NextRequest
    next_requests(...)                                            # delegate to the part owning the walk(s)
    to_fasta(arm=None, walks=None, with_seed=True, orientation='natural', width=None, *, coordinates=None,
             budget=None) -> str        # coordinates as Graphlet.to_fasta (model.py:571), per piece (Appendix D.2)
    completeness(*, budget=None) -> dict  # per arm: the root's complete_to_bp and one extension per junction
    evidence(arm=None, walks=None, *, budget=None) -> dict   # §3.5; walks=: only what those walks (keys, on arm) read
    summary(max_bytes=2048, *, budget=None) -> dict
    provenance() -> list[dict]            # §3.6
    save(path, *, budget=None)            # [1c] the bundle, Appendix D.3; metagraph.traverse.load(path) returns it
    digest -> str                         # sha256 of the canonical tree: the part keys and the parts' digests,
                                          # each with its retrieval identity (§3.7, C.3)
    memory_bytes() -> int                 # sum of the parts' models + the junction indexes
    # Deliberately refused (CombineRefused 'not_offered', text in §3.1): dump(), to_json(), compare(),
    # subgraph(), to_gfa() (phase 2), label_walks(), routes(), label_summary(), splits() (phase 2)

@dataclass(frozen=True)
class Part:      key: tuple; graphlet: Graphlet; offset_bp: int; digest: str
                 # key () = the root; (('right', 3),) = the child of root walk 3; (('right', 3), ('right', 0)) ...
@dataclass(frozen=True)
class Junction:  arm: str; walk: tuple; at_bp: int; overlap_bp: int; part: tuple;
                 verdict: str                                   # 'exact' | 'qualified'
                 causes: tuple[tuple[str, str], ...]            # (code, sentence), in the order of §3.4's table
                 continued: tuple[Label, ...]; not_continued: tuple[Label, ...]
@dataclass(slots=True)
class StitchedWalk:  key: tuple; arm: str; length_bp: int; sequence: Optional[str]; path_reason: Optional[str];
                     end_reasons: dict; labels_full: list[Label]; n_alive: int; complete: bool;
                     verdict: Optional[str]                     # 'exact' | 'qualified': its junctions and own causes;
                                                                # None for a walk no junction joins (its evidence is the root's)
                     certified_to_bp: int; beyond_certified_bp: int; junctions: tuple[int, ...];
                     pieces: tuple[tuple[tuple, int, int, int], ...]   # (part key, path_id, from_bp, to_bp)
                     causes: tuple                                      # this walk's own causes (§3.4)
@dataclass(slots=True)
class StitchedClaim(Claim):  walk: tuple; pieces: tuple; junctions: tuple[int, ...]       # [1b]

class CombineRefused(ValueError):  code: str     # every refusal of §2.3 and §3.1 (composite.py); MissingEnvelope
                                                 # and BadSelector as today
```

### C.3 Determinism and the algebra

**The union.**

- **Pure.** `union_pairs()` is a function of the results' texts and envelopes, the stated pairs, `seed` and
  `expected_pairs`.
- **Canonical order.** Parts and absent records are ordered by the pair key of §2.3, byte order, `None` first. Every
  list answer has a canonical order: arms left then right, walks by id under `by='id'`, ties of other rankings by
  id, spans and supports in pair order.
- **Commutative** for successful unions: any order of the results gives the same digest and the same answers, byte
  for byte. A refusal names the first failing element in the order given.
- **Identity element.** The union of one part answers as that part (§2.4, P2).
- **Idempotent for identical retrievals on one pair.** Given twice (same body, same retrieval identity), a no-op,
  the envelope-stated form kept whatever the order; another body, or the same body from another request, on that
  pair is refused (§2.3).
- **Digest.** The sha256 of the canonical list of (pair key, tag, `pair_stated_by`, retrieval identity,
  part digest), the absent
  records as canonical JSON (an error as its text), `seed` and `expected_pairs`. With the no-op rule above it does
  not depend on the order of the results.
- **Evidence.** `evidence(arm)` lists each part's own `evidence_block` (`ops.py:6006`) under its tag, the absent
  records, the expected pairs not included, each part's identity and pair source, and the five statements. An
  answer that reads several parts states, as `outcome`, the weakest value per dimension over the parts it reads,
  and each part's outcome under its tag. `exact` holds only when every part read is exact and verified.
- **Provenance and replay.** `provenance()` lists each part (`{tag, pair, pair_stated_by, identity, seed_id,
  digest, request}`, the request built as `GraphletStore.request_of()` does, `store.py:1277`), each absent record
  and the expected pairs. Replaying a union replays each part (the identity rules apply as for
  `traverse_fetch(replay=)`) and builds the union again.

**Stitching.**

- **Pure.** `stitch()` is a function of the parts' texts and envelopes, the walk named or inferred, and
  `allow_unverified_index`. Junction verdicts are recomputed from the parts every time, never stored.
- **Canonical tree.** A stitched graphlet is a tree: root, then a child per continued `(arm, walk)`. Parts are ordered
  by key, comparing steps as `(arm: left < right, walk id)`. Every list answer has a canonical order: arms left then
  right, walks by key under `by='id'`, ties of the other rankings by key. `digest` is the sha256 of the canonical
  tree (the part keys in order with each part's sha256 and retrieval identity, as the union's digest states retrieval
  identity), which is also what the bundle's manifest lists (Appendix D.3): equal bodies from other requests, whose
  junctions can read otherwise, give another digest.
- **Commutative across junctions.** `stitch(stitch(P, C1), C2)` and `stitch(stitch(P, C2), C1)` have equal digests
  when C1 and C2 continue different walks. Inference does not depend on the order (§3.1).
- **Associative along chains.** If C2 continues walk v of C1, then `stitch(stitch(P, C1, walk=w), C2, walk=(w, v))`
  and `stitch(P, stitch(C1, C2, walk=v), walk=w)` have equal digests (a StitchedGraphlet child is accepted for this).
  A refusal is the same in both folds: the identity checks of §3.1 run against every part on both sides, so a proven
  index difference anywhere in the chain is `different_index` whichever pair is stitched first.
- **Not commutative between parent and child.** `stitch(C, P)` is usually refused, deterministically, by the first
  failing check of §3.1: a parent with both arms gives `not_one_arm`, a one-arm parent whose walks do not continue
  into P's seed `not_a_continuation`. It is not always refused: on a circular sequence, P's continuation can be C's
  seed and C's continuation P's seed (k 15, radius 25, `continuation_bp` 15, Appendix G.6). Then both orders pass
  the preconditions, give different stitched graphlets, and are judged by W1–W8; W8 qualifies a walk that re-enters
  the root seed.
- **Idempotent for identical retrievals.** Stitching the same retrieval again at the same walk is a no-op, whatever
  continues it there (for a StitchedGraphlet child: the same subtree, every part's key, body and retrieval identity);
  a different body or subtree, or the same bodies from another request, there is refused `already_continued` (§3.1).
- **Identity element.** `stitch(g)` (no child) answers every offered query as `g` does, with keys `(w,)` for path
  ids `w` (tested on every fixture, §6 T2).

## Appendix D. Later phases in detail

### D.1 Claims and support across a junction (phase 1b)

- **Constrain** (external review, Must 4: the first version joined a parent's claim to every descendant, switched
  ones included, and so made claims no label holds). A claim is one label's support over one interval, so only a
  claim of the **same label** that is **contiguous** across the junction is joined: for each continued label X, the
  parent's claim of X that is open at the end of walk w (end class `radius`, `capped` or `pruned`: a leaf ended by
  `beam_pruned` also carries a continuation, SPEC §7.1) is joined with the child's claim of X on X's seed run, the
  run that starts at the child's seed. The joined claim runs from the parent claim's `from_bp` to `L + child.to_bp`.
  Its `evidence_from` is the parent's, or `L + child.evidence_from` when the child's displayed support starts later
  (a merge in the child). Its end fields come from the child; `entered_by` is the parent's; `pieces` lists both
  runs. A child claim of another label, entered by a switch, is **not** joined: it stays that label's claim, shifted
  by `L`, with its switch (`from` label, position, cost) kept. Lineage is stated apart from the claims: each claim
  names the claim it descends from (`descends_from`, through the child's `prev_run` chain and across the junction),
  so A then B reads as A over `[0, 120)` and B over `[120, 140)` entered from A at 120, never one claim over
  `[0, 140)`. `loss` and `branches` of a claim are its lineage's at its end: the parent's terminal values of the
  descent's seed label plus the child's. Under an `exact` junction (W3: `forbid`) every loss is 0 and no switch
  occurs. Under a `qualified` one with switching, the loss is that of the route the child chose, which can exceed
  the reference's minimum (an upper bound when the reference has that route). A `qualified` verdict qualifies the
  claims; it does not repair a wrong label or interval, which is why the join rule itself must be right. The
  parent's claims of labels in `not_continued` stay open at `L`, with the junction named; they are never turned into
  ends. Claims of other walks are unchanged. Order: the root's claims in the root's order, each joined claim in place
  of the open one it extends, then each child's remaining claims in part-key order.
- **Annotate.** No lineages. The children's claims are shifted by `L` without joining (`joined: false` in the
  evidence). That understates continuity and is honest. Joining by presence at the junction is phase 2.
- **Support.** A stitched walk's profile is the parent's on w over `[0, L)` followed by the child's on v, shifted. At
  `L`, `support_changes()` reports the labels in `not_continued` as removed with reason `junction` (the continuation
  did not seed them, so the removal is not biology), and those the child's seed dropped with reason `junction` and
  the `X` reason. The junction never adds a label: none is reported as added at `L` because of it.

### D.2 Record coordinates across a junction

Coordinates are kept per piece and never joined. Under `support: kmer` a retrieval has none (DESIGN §18.1). Under
`trace` the junction is `qualified` (`trace_restart`), and the two pieces of a joined claim are two occurrences unless
their coordinates touch in one record (the parent's run ends at record position e, and the child's run's own bases
start at e in the same record and orientation). Phase 1b states `contiguous: null` for every joined claim; proving
it from the coordinates, and so lifting `trace_restart` for such a junction, is phase 2. `to_fasta(coordinates=)`
states each piece's coordinates in the header, as `Graphlet.to_fasta` does for one walk.

### D.3 Serialization: without a bundle, and the bundle (phase 1c)

The union's FASTA and GFA are in §2.9. This section holds what both combinations share, and stitching's exports.

- **Without a bundle (phase 1a).** Each part is saved as its own `.mgt` file, as today. Building the combination
  again from the loaded parts gives the same digest, because `union_pairs()` and `stitch()` are pure. A union also
  needs its absent records, stated pairs, `seed` and expected pairs, which the caller keeps as JSON in 1a.
- **The bundle** (phase 1c, one format for both kinds; `save()`, `load()`; `.mgtb`). One UTF-8 text file, written
  atomically like `save()`, read by `load()`, which dispatches on its first line as it already returns a view for a
  saved view. A stitched graphlet:

  ```
  #mgt-bundle 1 stitch
  {"allow_unverified_index":false,"parts":[{"key":[],"bytes":12345,"lines":210,"sha256":"..."},
   {"key":[["right",3]],"bytes":2345,"lines":40,"sha256":"..."}]}          <- one line, compact JSON
  <part 0: its standalone .mgt text (H, J, body, Z), byte for byte what save() writes>
  <part 1: ...>
  ```

  The manifest holds only what cannot be recomputed: part keys, sizes, digests and the opt-in. There is no field
  for a verdict, cause or certificate, and an unknown manifest member is refused, so a bundle cannot carry a
  laundered "exact". `load()` parses every part with the usual validation (`A` counts, `Z`, the `J` rules; a fault
  raises `GraphletFormatError` naming the part and its line), checks each digest, and rebuilds by calling `stitch()`
  in key order with the stored opt-in. A bundle whose tree does not hold gets the same codes as §3.1. `save(load(f))`
  is byte-identical to `f`.

  A pair union is written the same way, with first line `#mgt-bundle 1 union` and a manifest of `parts` (`{tag,
  pair, pair_stated_by, bytes, lines, sha256}` in canonical order), `absent` (the absent records as they were
  given), `seed` and `expected_pairs`. No coverage category, absence value or statement is stored: `load()` rebuilds
  by calling `union_pairs()`, so every statement is derived again, and a bundle cannot carry a stronger one.
- **Stitched GFA** (phase 2). Every part's segments under a part prefix, plus a link from the parent's leaf segment
  to each child root segment with the de Bruijn overlap (k − 1)M, which `to_gfa()` already uses inside a part. The
  child's seed is not a segment of its own. Since n ≥ k, the child's first step reads only its seed. **[unsure]**
  how the link is written when L < k − 1 and the overlap reaches into the root seed.
- **Not offered.** `dump()` / MGT (§4.4), `to_json()` (detail `full` of one server result), JSON views of a composite
  (phase 2: `to_dict()` with the parts' summaries, and the junctions or the statements).

### D.4 The same-pair combination (later, §4.2)

`union_retrievals(results, *, seed=None, budget=None) -> RetrievalUnion` (names provisional), designed now, built
when a user needs it:

- **Sources are retrievals.** Each input is one retrieval of one seed; retrievals are told apart by retrieval identity
  (§2.3) and each is named by a display tag made unique as pair tags are (Appendix A.4). A pair may repeat, retrievals
  of several pairs may be mixed, and each pair groups its retrievals. One oriented seed, as in §2.3; identical
  retrievals are no-ops.
- **The same machinery.** The trie, spans, routes, absence, both path modes and the exports of §2 apply with
  retrieval in place of pair. A path that switches between retrievals of one pair is an inferred path, with the
  connection rule of §2.7.
- **What it states.** Coverage per retrieval, never for the union: no depth is certified for a pair's retrievals
  together. Turn-offs per retrieval under its own strategy. The scope statement: "This union answers none of its
  requests: each retrieval's certificate is its own."
- **Labels.** Within one index labels can join by ref where the retrievals' `index_fp` agree (SPEC §7.5.7), which a
  cross-pair union can never do; without `index_fp` on both they stay per retrieval.
- **One exact case, later.** Where SPEC §6.4's conditions hold (`keep`, `forbid`, no quorum or beam, every
  retrieval complete to the depth asked), the union of the retrievals of label sets A and B equals one retrieval
  naming both, and the operation could say so with that request. Only after an oracle fixture pins it, as T4 does
  for stitching.
- **Until it is built**, `union_pairs()` refuses two retrievals of one pair (`duplicate_pair`, naming this operation
  as planned), and `compare(mode='prefix_subset')` shows whether one contains the other.

### D.5 The two arms of one seed: a later validation helper

The owner's provisional answer (OQ2, 2026-11-15): "lower priority, mostly useful for sanity checking/validation". As an
answer, the union of a `direction: left` fetch and a `direction: right` fetch adds nothing that a `direction: both`
fetch does not give. As a check it is useful: the two one-arm fetches and the two-arm fetch of one seed should agree arm
by arm. Where they do not, either a shared cap explains it, or something is wrong in the engine or the library.

A later item, `validate_arms(both, left, right, *, budget=None) -> ArmCheck` (name provisional):

- **Preconditions.** Stitch's identity checks (§3.1: one index, one pair; a one-sided `index_fp` refused), one
  oriented seed, `both` walking both arms, and `left` and `right` one arm each. Strategy differences other than
  `direction` are reported, not refused.
- **The check.** `compare(both, left, arm='left', mode='walks')`, the same for the right arm, and the seed-level
  records (`S` including `validated_seed_id`, `X`, seed-level `K` including `derivation`), which must be equal.
- **Expected differences.** In one request the arms share caps: `max_steps`, `max_paths`, `max_output_bp`, and the
  time, memory and work budgets. A cap reached in `both` can cut an arm that its one-arm fetch walked further. The
  check attributes a difference to a shared cap when `both` states a resource stop (`Q`) and the difference lies
  beyond that arm's `complete_to_bp` in `both`. Any other difference is a finding.
- **Not a combination.** It returns a report, not a composite, and needs nothing of §2 or §3 beyond the identity
  checks and `compare()`.

## Appendix E. Tests in full

Stdlib `unittest` with seeded random generation, as in the existing suites; no new dependency.

### E.1 The union (P1–P13)

- **P1 refusals.** One test per code of §2.3, with the exact text, including every variant of `different_seed` (first
  difference, unequal lengths, reverse complement; `seed=`, an absent record's `length_bp` and a stated `sequence`
  against the parts) and of `duplicate_pair` (with the remedy for absent records). Check-major order: `[no_bases part,
  not_combinable element]` gives `not_combinable`; arguments first (a string `expected_pairs`). An identical retrieval
  on the same pair is a no-op, and the envelope-stated form wins in either order; one body from two requests on one pair
  (`max_splits_per_path` 64 and 128, nothing reached) is refused `duplicate_pair` with the retrieval text; the same body
  on two stated pairs makes two parts. No parsed model is shared or merged between parts; `memory_bytes()` counts a
  model object given twice only once.
- **P2 the union of one part is that part.** On every committed fixture with bases (one without them is refused
  `no_bases`, P1), both arms: the union walks of `kind: walk` are exactly the part's walks (matched by spelling), each
  `whole`, with their lengths, end reasons, span labels (`from_bp` 0 exactly for `labels_full`) and `complete`; every
  other union walk is a stub walk whose `route` is `held`, with `goes_on` unless its merge is a zero-length leaf
  (`limits.mgt` line 21); the coverage is the part's arm certificate.
- **P3 order.** Seeded random permutations of the results give equal digests and byte-equal answers for every query
  and both exports.
- **P4 inputs.** Every row of the decision table of §2.3: a `Graphlet`; a one-seed `TraverseResponse` with a graphlet,
  with an error result, and fetched with another `detail` (`no_graphlet`); `{response, index}` on a two-seed response
  (the pair still stated by the envelope); `{response, sequence}`; `{pair, graphlet}` (CLI bodies under two stated
  pairs); `{pair, result}` with and without `sequence`, and with the `error` beside it as in `TraverseResponse.errors`;
  a result with a graphlet (`result_has_graphlet`); `{pair, error}` with a string and with an exception; zero and two
  seed results (`not_one_seed`); a union of absent records only, with `seed=` (arm-less coverage), and a fan-out whose
  every pair lacks the seed: `union_pairs(r.results, seed=r.seed, ...)` builds the arm-less union, and without `seed=`
  it is refused `no_seed`. The pair-match rule: bare names against envelope pairs (an absent `g2` and the part on `g2`
  are one pair; `expected_pairs=['g2']` is not `not_included`), two databases' same-named graphs (two pairs by
  `index_fp`; `ambiguous_pair` for a bare name that matches both). Two pairs with one graph name and release, different
  paths and no `index_fp` get two display tags with `~` suffixes (Appendix A.4), and no dict, label or span of one is
  read as the other's. Retrieval identity: the same retrieval twice is a no-op that keeps the envelope-stated form in
  either order; two identical absent records are one; a failed record and then an answer keep the answer, with the
  failure in `provenance()`; one body from two requests, and two requests that differ only in `seed_id` against a
  changed strategy, are told apart as §2.3 says. `seed_check` gives `sent` for a response from the client that kept its
  seeds, `stated` with a sequence and `unchecked` otherwise, sequences compared after case mapping, also when the
  record's `seed_id` happens to equal the union's seed id (a caller's text proves nothing).
- **P5 the trie.** Hand-built cases: two pairs that share a prefix and branch; a pair whose walk ends inside another's
  (`ends`, with its reason), and the same walk kept as a union walk by its terminal marker (A answers `AC`, B `ACGT`:
  `AC` is a walk, `whole` in A and `extends` in B, `walks(pairs=['A'])` returns it, and adding B to a union of A removes
  nothing); support that changes inside one segment (the `reentry` fixture: A, B, A as three `LabelInterval`s in order);
  a root leaf (`ends` at 0) against an arm not walked (`none`); a pair that turns off (`diverges`); a merge whose
  non-first parent is another pair's walk (`merges`, not `diverges`), whose route is `held` to the walk's end, so the
  pair is in `held_in` but not `whole_in`; a route that crosses a second merge, one that ends at a part leaf and one
  that diverges below the merge, each with its absence; a merge segment that lies only on a stub (nested merges); span
  labels with `from_bp` and `continues` on the bubble; union segments cut where a merge segment starts, and the GFA link
  from the stub; stub walks left out of `walks()` and FASTA unless `stubs=True` and counted apart; zero-length merge
  segments, a merge at the radius (`limits.mgt` line 21, `G 1,2 7 0`) and a merge node that splits at once (two
  children), where `into`, the `merges` token and `mt:Z:` name the trie node `s<id>$` and the inferred GFA links the
  stub to each child; span and segment labels from the alive sets, on the committed `merge.mgt`: at 40 bp of walk 0,
  `a.fa` is a span label with `from_bp` 0 and `b.fa` one with `from_bp` 36, and on the stub [0, 2] `b.fa` is in
  `continues`; a stub end that is a non-first parent of two merges refused (`GraphletFormatError`); segments that break
  where the set of pairs changes; children in base order. Every union walk is `whole` or `merges` in some pair.
  Construction stays linear on a comb (the pattern of the speed-up tests), supports included.
- **P6 the split-graph oracle and the non-claim.** One tiny index split by record into two graphs (two "taxonomic
  groups") with one record in both (a "plasmid"), and the unsplit index, under `forbid` and generous caps, once with
  `keep` and once with `merge` (r237's default), the shared record forming a bubble with each group's records. In both
  modes, every union spelling is a prefix of a spelling of the unsplit retrieval at the same radius, the shared record's
  walks are carried by both pairs, and a walk of the unsplit retrieval that no single graph holds entirely is not in the
  union, while the coverage sentence and both exports state the non-claim (the fixture is built so that such a walk
  exists). Annotate mode: every walk of the unsplit retrieval whose k-mers all lie in one graph is in the union, whole
  or as a prefix ended by that graph's `dead_end`. Constrain mode: every walk of the unsplit retrieval carried by a
  label whose record lies in graph g is in g's part; a walk whose k-mers are all in g but whose only carrier is a record
  of the other graph need not be.
- **P7 labels.** Every label is shown with its pair; `label({'name'})` resolves for a name in one pair and is refused
  `ambiguous_name` for the shared record's name; `shared_names()` lists it, and the statement counts it. `{pair, name}`
  on the `same_name` fixture (two labels ACC1 in one pair) raises the part's `AmbiguousLabel`; `{pair, ref}` resolves.
  Header-name warnings: `count_only`, `unknown` (a pair not in `seqs_order=`, and every pair without it) and `None`
  pairs; header labels marked and column labels not, in the statement, `carriers()`, FASTA `mapping=` and the GFA
  comment; a `verified` pair unmarked; `bad_seqs_order` for a value outside the three.
- **P8 coverage.** Every category of §2.6: absent records `not_in_graph` with none and with some k-mers present,
  `labels_not_in_graph`, a failed derivation (`not_walked`, seed in graph) and `{pair, error}`; a part with a
  `seed_labels` limitation (`max_seed_labels` below the carriers) and one with a `scope` limitation
  (`complete_qualified`); parts walked to different radii and with a resource stop; a one-arm part; `expected_pairs=`
  with a missing pair, and no `expected_pairs=`. The sentences are exact; the default depth is the largest radius, for
  `by_category` too; a part truncated beyond the depth asked (`complete_to_bp` 420, asked at 300) is `complete` at 300
  and `incomplete` at 500; and `every_pair_certified_to_bp` is never above any part's `complete_to_bp`, and is `None`
  with a qualified pair and without `expected_pairs=`.
- **P9 exports.** FASTA spellings equal the parts' `spell()`; headers carry every token (spans with their absence,
  `extends=`, `modes=`, lower-bound `+`) and leave empty ones out; names with every reserved character, and with
  non-ASCII characters (`café`), escape and decode back to themselves, and every GFA line is printable ASCII. GFA, **by
  the paths it admits**: on small fixtures every maximal path the exported file admits (oriented walks along its links,
  from any segment) is enumerated. Default mode: each admitted path is spelled whole by one pair (a substring of one
  part's walk or route, the seed included); the case of Appendix G.4, Must 2 (A: `AACGT`, and `TACGT` through a merge;
  B: `AACGA`) admits no `TACGA`; no path joins one pair's left flank to another's right flank. Inferred mode (the later
  phase): the admitted paths are exactly those the connection rule builds (§2.7), each transition at a stub link or the
  seed. Then the comment lines, tags and left-arm orientation, as the existing `to_gfa()` tests check them.
- **P10 absence.** Each value of §2.4 on fixtures: annotate mode with `keep` inside the certificate (`graph`); constrain
  mode, and a `united_history` scope (`under_strategy`); a turn-off at or beyond `complete_to_bp` (`uncertified`); a
  `dead_end` (`graph`), also at the certificate's edge, where it comes before `uncertified`. `walk_rule` on three cases
  of Appendix G.5: a canonical pair that ends at 0 bp with `dead_end` on a skipped hairpin inside its certificate,
  against a basic pair that spells 30 bp more, whose k-mers the canonical graph holds; `hairpins: skip` against `follow`
  on one graph, where `follow` has two walks more inside a `per_path` certificate; a k 11 pair that ends by `edge_reuse`
  at 56 bp against a k 21 pair that spells the record's 90 bp. The span states the rule's k and regime. An event decides
  only for its own character (the k 2 case of Appendix G.6): A blocks A by seed re-entry at a node where B goes on with
  G, and A's span there is `graph`, not `walk_rule`.
- **P11 mixed strategies.** Parts in constrain and annotate mode, and parts at two k: the union builds and states
  each pair's mode and k in the sentence, on each span and in `modes=`.
- **P12 budgets and memory.** Every operation under `budget=` completes or stops with the same row in every run.
  `memory_bytes()` bounds the traced peak on 11 parts (the pattern of `test_traverse_memory`). `summary()` refuses
  `summary_too_small` below the statements' size and never cuts them.

- **P13 inferred paths (the later inferred phase).** On hand-built parts with a stub of one pair merging into a segment
  another pair continues below: the inferred path is listed with its pieces, cut wherever the set of pairs holding the
  path changes (P's own route included, so a piece held by two pairs is listed as such), one `merge` transition (`by`,
  `at`, `at_bp`, `from_bp` and `to_bp`, the pairs before and after), no `held_in` and no `complete`; its id is `i:` and
  its union segments, the same whatever `max_transitions` lists it; a union walk is listed in both modes without a
  transition (`kind: walk`); on the committed `merge.mgt` the default mode lists one walk and two stub walks and the
  inferred mode one `walk` and three `route`s (`held_in` the part's pair, no transition); a chain of bubbles is cut at
  `max_paths` with `max_transitions` 0, the same cut in every run, with `cut` true, `count`, a `resume` that continues
  it, `of=<count>+` in FASTA, and `top=` ranking only the paths enumerated; exactly `max_paths` paths without a cut read
  `cut` false; two pairs whose merges at one trie node permit a connection are both in `by`; `max_transitions` bounds
  the transitions, which are stated; no junction is made at equal sequence in another branch or position (no overlap
  assembly); the fifth statement reads as §2.7 says in that mode.

Beside them, one integration check on a two-graph server (`integration_tests/test_multigraph.py`): a seed held by
one graph only, so the other answers `not_in_graph`, read as §2.3 and §2.6 expect from a real response. The new
graphlets.rst examples run in `run_examples.py` (CI's fixture-checks job).

### E.2 Stitching (T1–T11)

- **T1 refusals.** One test per code of §3.1, with the exact text, on committed fixtures. Identity cases edit `H`/`J`
  as W5's envelope tests do (a manifest digest on one side, on neither, different ones; a different pair; a different
  release). `different_pair` and `different_graph` are tested without manifests, where they are reachable.
- **T2 identity element.** On every committed fixture, for every offered query, `stitch(g)` answers as `g`, with keys
  mapped.
- **T3 algebra.** Commutativity across junctions, associativity along chains and canonical bytes, as digest equality
  over seeded random construction orders. `stitch(C, P)` gives `not_one_arm` for a two-arm parent and
  `not_a_continuation` for a one-arm parent whose walks do not continue into P's seed; on the circular fixture (Appendix
  G.6) both orders are accepted, differ, and carry the verdicts W1–W8 give them. An identical body at a continued walk
  is a no-op; a different body is `already_continued`, and so is the same body from another request. Inference over
  continued and uncontinued candidates gives one code in every order.
- **T4 oracle triples.** CLI triples on tiny indexes: a parent at radius R, its child from `next_request()`, and the
  uninterrupted retrieval at R + E. The CLI fixtures have no `index_fp`, so T4 asserts junction and walk verdicts, not
  the evidence's `exact`. (a) Each whitelist case is pinned (in phase 1a-S as `exact_gated`): `keep` and `forbid`;
  `merge` with w the arm's only open walk; an unlimited branch allowance; an equal finite one. In each, the stitched
  walks below w equal the uninterrupted retrieval's subtree of w's prefix, per path, and every junction says `exact`
  once it is switched on (`exact_gated` before). A whitelist entry without its passing fixture is not shipped as
  `exact`. (b) One fixture per cause, each asserting the stated cause and `qualified`: `sibling_merge_unknown` (two
  radius-ended walks under `merge` that reconverge beyond R); `trace_restart`; `strategy_differs` with
  `label_count_knob` (a changed quorum knob); `branches_unequal` (a reset branch allowance: W2 leaves
  `max_label_branches` to W3–W4); `extra_foreign`; `switching`; `not_continued`; both `seed_reentry` sub-cases (a
  circular record); `split_limit_restarted`, also with unequal split counts on a merge's parents (`splits45`, Appendix
  G.5: k 11, constrain, `forbid`, `merge`, `max_label_branches` 0, `max_splits_per_path` 2, radius 45; one split on w's
  displayed chain, two on the other parent's; the stitched walks are 85 bp, the reference stops at 50 bp with
  `split_limit`); `rule_violation`; `candidate_changed` with the counterexample of Appendix G.4 (`loop_branch`: basic, k
  15, one record `S P U C U Y` whose 15-mer U repeats; `forbid`, `max_label_branches` 1, `keep`, radius 115,
  `continuation_bp` 40; the child stops at 0 bp with `branch`, the reference goes on to 145 bp: the stitched walk is
  `qualified`, certified to 115 bp and not complete); `decisions_unrecorded` (the child's branch events cut at 1; the
  `loop_branch` child with `max_branch_events` 0, whose boundary is 0 and whose walk is 0 bp; a boundary equal to a
  walk's end). Every fixture asserts, beside the verdict and the causes, each stitched walk's `certified_to_bp` and
  `complete` and the evidence's outcome and `exact`, in the public result and in the prospective one (§3.5): the
  `loop_branch` walk is complete in neither. One fixture with the changed decision inside the child, past its first
  base, where the prospective certification stops at that depth and not at the junction, while the public result stops
  at the junction. Where the uninterrupted retrieval differs, the difference lies below a stated cause. (c) Promotion
  candidates (`same_stop`, uniform-shift switching) get fixtures that show their walks equal the reference's before they
  are moved to the whitelist.
- **T5 no laundering.** Over all fixture combinations: every answer's evidence is, per dimension, no stronger than any
  part or junction it reads. Every part's `K` kinds, `Q` and `X` appear under its key. A stitched walk is complete only
  if all its pieces and junctions are and no walk-scope cause lies on it; in phase 1a-S no stitched walk that reaches
  into a child is complete, and no answer that depends on a junction reads `exact`, while `stitch(g)` and the root walks
  no child continues keep the root's evidence (T2). With fixtures whose `H` carries an `index_fp`, the evidence's
  `exact` is checked too.
- **T6 round trip.** In 1a: saving the parts and stitching the loaded parts gives an equal digest. In 1c:
  `save`/`load` byte-identical with equal digests, and `load` refuses a truncated part (`GraphletFormatError` naming
  part and line), a digest mismatch, an unknown manifest member, a tree whose child is not a continuation (the §3.1
  code), and a one-sided identity without the stored opt-in.
- **T7 exports.** Stitched FASTA equals the parts' spellings joined at `L` without the overlap, on both arms in both
  orientations, with `L < n` included (`fork`). Headers carry junctions, verdicts and `certified_to`.
- **T8 budgets and memory.** Every operation under `budget=` completes or stops with the same row in every run. The
  account bounds the traced peak on a chain of comb continuations (the pattern of `test_traverse_memory`).
- **T9 docs.** The new graphlets.rst examples run in `run_examples.py` (CI's fixture-checks job).
- **T10 claims across a junction (phase 1b).** A child that carries the continued label A for 20 bases and then
  switches to B for 20: the stitched claims are A over `[0, L + 20)` and B over `[L + 20, L + 40)` entered from A,
  with `descends_from` linking them, and no claim spans both labels.
- **T11 randomized comparison (the gate for `exact`, OQ19).** Seeded random small graphs (generated records with
  repeats, loops and bubbles; basic, canonical and primary regimes; both label modes), each with parent, child and
  uninterrupted traversals, and chains of three or more parts. For every stitched walk the rules call `exact_gated`, the
  comparison checks the prospective result (§3.5) against the uninterrupted traversal: its endpoint and end reasons,
  `complete` and `certified_to_bp`, the evidence and the outcome per dimension, not only its spelling. A second check on
  the same runs asserts that the public result only reduces the prospective one: `certified_to_bp` no larger, `complete`
  only where the prospective result is complete, evidence and outcome no stronger in any dimension; runs with branch
  events cut (`max_branch_events` down to 0) are included. Every disagreement fails the gate and is kept, minimized, as
  a fixture. `exact` is switched on only when the comparison passes; it is not part of phase 1a-S.


## Appendix F. Cost, code plan, docs, today, service and risks

### F.1 Cost and memory bounds

**The union.**

- **Parts.** A union holds its parts' models: about 5 MB per constrained arm and about 20–30 MB label-free, measured
  (DESIGN §7). r237's 11 parts with two arms each hold about 110 MB constrained, or 440–660 MB label-free.
  `memory_bytes()` is their sum plus the trie.
- **The trie.** Built once, arm by arm, by one depth-first pass over each part's first-parent tree (walks and stubs,
  §2.4). Bases are read from the part's segments, not copied, so the work is O(sum of the parts' bases), linear on a
  comb too. A union segment costs about 200 bytes plus about 120 bytes per pair spelling it, and its bases are a
  reference into a part [estimate]; a support holds its pieces and two counts, never a list of walks. The trie
  has at most as many segments as the parts together, plus one per point where a pair's spelling ends, merges or
  turns off and one per merge segment's start: for 11 parts of about 2,000 segments each, about 22,000 union
  segments, 5–10 MB [estimate].
- **Labels.** Per (union segment, pair), the labels alive on the part segments it covers, as intervals. A part
  segment lies on one path, so its intervals do not depend on which walk passes it. In constrain mode they are built
  from the segment's entry set and its label changes (`derive.label_changes()`, `derive.py:72`), each change opening
  or closing one interval: O(entry labels + changes) per segment, O(the part's runs + segments) for a part. In
  annotate mode consecutive `presence` runs that hold a label extend its interval: O(the runs' sizes). The library's
  `_alive_pieces` (`ops.py:564`) is not used for this: it copies the whole alive set at every change, so a segment
  whose N entry labels end at N different bases costs N(N+1)/2 label ids (500,500 for N = 1,000, run against the
  library), where N intervals suffice. The leaf-anchored `_displayed_profile` (`ops.py:492`) is not
  used: it misses labels alive at a cut depth (`derive.py:804-805`). Counts are kept; full lists are built on demand
  (`carriers()`, `labels='full'`).
- **Routes and span labels.** A route costs O(the walk's remaining length) per `merges` span, through the
  part's segments; a span's labels with their `from_bp` come from the per-segment labels and the partitions of the
  merges crossed. Both are computed on demand: `carriers()` for one walk, `walks(by='pairs')`, `pairs=` and the
  exports for the walks they read, at worst O(union walks × arm length) for an arm, charged to `budget=`.
- **Inferred paths.** Enumerated on demand from the trie and the stubs' merge links, depth first in path-id order,
  within `max_transitions` and `max_paths`: O(paths listed × their length), charged to `budget=`. Their number grows
  with the merges on a pair's own path (up to 2^B routes over B bubbles in a row, with no transition) and with the
  transitions allowed, so `max_paths` bounds the enumeration whatever `max_transitions` is; the listing states the
  count it gave and whether it was cut.
- **Queries.** `walks()`: O(union walks × pairs), plus the routes it reads. `carriers()`: O(depth + pairs ×
  labels). FASTA and GFA stream as `to_fasta()` and `to_gfa()` do. Every operation takes `budget=` and charges the
  parts' cold prices plus the trie work, with prices in `budget.py` (the cost model states each one).
- **Bundle.** The parts' texts, the absent records and about 100 bytes per part.

**Stitching.**

- **Parts.** A stitched graphlet holds its parts' models: about 5 MB per constrained arm and about 20–30 MB
  label-free, measured (DESIGN §7). `memory_bytes()` is their sum plus the junction indexes. A parent continued at
  many walks holds as many child models. That is the caller's choice, bounded by `budget=`.
- **Junctions.** Per junction: the continuation check is O(n), with n up to 1000 (§3.2). The rule check builds the
  set of (k+1)-mers of the DAG ancestors of w's leaf (any parent; at a merge all parents end in the same oriented
  node, so any parent's last k bases give the window across the boundary), O(bases of those segments + |seed|),
  about 100 bytes per k-mer in a Python set [estimate], so about 1 MB for a 10 kb chain. It then scans the child's
  bases once per segment, O(bases of the child). The decision check (W8) also looks up, along each stitched walk's
  route through the child, every successor the child recorded (`G` children, `V` characters, `E` blocked and
  skipped): O(route bases + events). For chains, each ancestor's set is looked up in turn and not copied.
  The local check compares (k+1)-mers exactly. The engine compares 128-bit hashes where a (k+1)-mer does not fit 64
  bits (DNA k ≥ 32, protein k > 11; SPEC §6.6), so a hash collision can block a step the local check allows. That is
  ignored, and the docs say so.
- **Queries.** Walks: O(sum of the parts' walks). Spelled walks stream as the exports do (peak about twice the text).
  Claims: O(sum of the parts' claims + the continued labels × the child's runs). A profile is two profiles
  concatenated. Every operation takes `budget=` and charges the parts' cold prices plus the junction work, with prices
  in `budget.py` (the cost model states each one).
- **Bundle.** The sum of the parts' standalone texts plus about 100 bytes per part.

**Budgets, both classes.** A public method of `PairUnion` or `StitchedGraphlet` whose work grows with the walks or
segments takes `budget=` (keyword-only; `None` is the ambient `local_budget()`) and charges it. Every derivation a call
reads is charged at its cold price whether a cache holds it or not, the object's own derivations (a union's walk
figures and segment labels, a stitched graphlet's rows, rule-check contexts and walk causes) as the parts' are, so
a call charges the same, and stops at the same row, whatever earlier calls cached; warm calls pay more than they
spend.

- **`PairUnion`.** `statements()`, `segments()`, `walks()`, `spell()`, `carriers()`, `to_fasta()`, `to_gfa()`,
  `summary()`, `labels()`, `label()`, `shared_names()`, `coverage()` and `evidence()`; `print(u)` charges the ambient
  budget. A stitched part inside a union runs under the union call's budget as its ambient budget.
- **`StitchedGraphlet`.** `stitch()`, `walks()`, `walk()`, `spell()`, `continuation()`, `next_request()`,
  `next_requests()`, `to_fasta()`, `summary()`, `evidence()` and `completeness()`. `walk()` charges the walk it builds
  as `walks()` charges each; `summary()` reads the evidence under its own budget.
- **Exempt.** The accessors that return stored fields or measure the object: on `PairUnion` `arms()`,
  `provenance()`, `memory_bytes()` and the stored attributes (`parts`, `absent`, `pairs`, `expected`, `digest`,
  `seed`, `alphabet`); on `StitchedGraphlet` `arms()`, `arm_of()`, `label()`, `provenance()`, `memory_bytes()` and
  the stored attributes (`parts`, `root`, `junctions`, `labels`, `identity`, `digest`, `seed`, `k`, `regime`,
  `alphabet`, `mode`, `support`). The prospective test hook (§3.5) is uncharged.

### F.2 Code plan

| file | change |
|---|---|
| `composite.py` (new; 1a-U) | the parts container, `CombineRefused`, the shared preconditions (envelope, bases, alphabet, the index identity rules), the weakest-of evidence, digests, the budget plumbing, the streaming speller over pieces |
| `union.py` (new; 1a-U) | `union_pairs`, `PairUnion`, `PairPart`, `PairAbsent`, `PairLabel`, `SpanLabel`, `LabelInterval`, `SegmentSupport`, `UnionSegment`, `PairSpan`, `Carriers`, `UnionWalk`, `Coverage`, `Statements`; the input classification, the tags and retrieval identity, each pair's `seqs_order` and the header-name warnings (OQ20), the trie with terminal markers by one pass over each first-parent tree, the spans and their absence, coverage, the statements. The later inferred phase (OQ18) adds `InferredPath`, `InferredListing`, `Transition` and the inferred paths. The queries delegate to the parts' `ops` (no derivation reimplemented) |
| `_stitch.py` (new; 1a-S) | private, so that `metagraph.traverse.stitch` is only the function: `stitch`, `StitchedGraphlet`, `Part`, `Junction`, `StitchedWalk`, `StitchedClaim`; the preconditions; the whitelist and causes; the decision check (W8); the queries, which delegate to the parts' `ops` and shift |
| `ops.py` | split `_identity` into the index rule and the seed rule (the index rule shared, no behaviour change; 1a-U); a public per-segment label-interval accessor, built from the entry set and `derive.label_changes()` (constrain) or the `presence` runs (annotate), linear, priced like `_alive_prices` (Appendix F.1; 1a-U); expose the (k+1)-mer edge key of SPEC §6.6 (1a-S) |
| `client.py`, `attempts.py` | `TraverseClient.traverse_pairs()` and `PairResults` with the seed sent and each pair's `seqs_order` (§5.2), and `TraverseResponse.seeds_sent` (`attempts.py:363`), 1a-U. `TraverseClient.response()` already keeps results without a graphlet with their errors (`client.py:895-911`), which `union_pairs()` reads from a one-seed response |
| `parser.py` | (1c) `load()` dispatches on `#mgt-bundle`; the atomic writer reused |
| `export.py` | the union's FASTA and GFA in the source-supported mode, with the token escaping (1a-U); the inferred records and links (the later inferred phase); stitched FASTA (1a-S) |
| `budget.py` | prices for the trie, the union queries, stitching and the stitched queries |
| `__init__.py` | exports |
| `scripts/traversal/graphlet_fixtures.py` | `UNION_FIXTURES` (1a-U): one tiny index split by record into two graphs with one shared record, the unsplit index, and the split index at a second k (P6, P11); `STITCH_FIXTURES` (1a-S): CLI parent, child, uninterrupted triples (Appendix E.2 T4, `loop_branch` among them), the child's request committed as JSON and pinned to `parent.next_request()` by a test. The script still does not import the library |
| `api/python/tests/test_traverse_union.py`, `test_traverse_stitch.py` (new); `tests/data/traverse/documents/union/`, `stitch/` (new) | §6 |
| `integration_tests/test_multigraph.py` | one union check on a two-graph server (Appendix E.1) |
| `docs/source/graphlets.rst`, `docs/DESIGN-traverse-graphlet.md` | Appendix F.3 |

The CLI answers single-graph and has no `index_fp`, so the union fixtures' envelopes get their pair tags and an
`index_fp` added, as W5's envelope tests edit `H` and `J`, or are passed with stated pairs (P4); the absent records
follow SPEC §6.1's shapes (steps 2–4); the integration check reads a real one once.

Not touched: the engine, `mcp_tools.py`, `store.py` (phase 2), MGT v1, the service. Size [estimate]:
`composite.py` 300–400 lines; the union, after the cuts and with both path modes, 900–1,300 lines and about 850
lines of tests; stitching 1,200–1,800 lines and about 1,000 lines of tests; about 300 lines of docs.

### F.3 Docs

- graphlets.rst, the intro paragraph (lines 23–30), becomes: "**Graphlets are not merged.** Two combinations are
  offered. The graphlets of one seed answered on several pairs (the graphs of a release, or two databases) can be
  combined by sequence (`Combining one seed across pairs`_). The result keeps every retrieval unchanged and states, with
  every answer, how complete each pair is, which pairs carry each walk and how far, that a pair not following a walk did
  not follow it under its own request, and that labels are each pair's own. Every path it gives is held whole by one
  graph; inferred paths that switch graphs where the parts allow it, with every switch named, come in a later release. A
  continuation can be stitched to the walk it continues (`Stitching a continuation`_); the result states, at each
  junction and walk, why its walks may differ from one uninterrupted walk (`qualified`, with named causes); `exact`
  comes in a later release. Combining two retrievals of one pair is planned; every other combination is not offered,
  with the reason (`What is not combined`_)", followed by the existing sentences on `subgraph()`, `compare()` and
  `on_reconverge`.
- graphlets.rst, new sections. "Combining one seed across pairs", after "Comparing retrievals": the five statements
  first, then the inputs, the trie with its spans and absence values, the coverage categories, labels per pair, the
  source-supported mode with the r237 example (the inferred mode and its connection rule when that phase ships), the
  header-name warnings, one example on the split-graph fixture, FASTA and GFA (what the file admits). "Stitching a
  continuation", after "Continuations and the next request": offsets, orientation and the continuation's length, walk
  keys, the whitelist and the causes, the certified depth, one example on a committed triple, marked experimental. "What
  is not combined": the table of §4.
- `DESIGN-traverse-graphlet.md`: `composite.py`, `union.py` and `_stitch.py` rows in the §5 module table, a sentence
  in §5.3 pointing to this note, and the public-API block.
- No spec change: no wire, record or field changes.

### F.4 The library today, in detail

A graphlet is one retrieval: one seed, one index (pair), one normalized strategy, up to two arms, a certificate per
arm (`complete_to_bp`, `completeness_scope`) and per seed (`outcome`, `limitations`, `resource_stop`). Every operation
reads one graphlet or compares two:

- `subgraph()` (`ops.py:3502`, `GraphletView`, `ops.py:3290`) is a filtered view of one graphlet. Ids are not
  renumbered, `save()` writes the unchanged backing body plus `J.view`, and completeness is qualified "for the
  selected labels".
- `compare()` (`ops.py:3778`; identity rule `_identity`, `ops.py:3588`) compares two retrievals and refuses nothing.
  `comparable` is `False` for different indexes or different oriented seeds, `'unverifiable'` without `index_fp`,
  and `'qualified'` when scopes or label evidence differ.
- A continuation (`next_request()`, `model.py:547`; `TraverseClient.deepen()`, `client.py:912`;
  `GraphletTools.traverse_continue()`, `mcp_tools.py:1488`) is a **new traversal** whose result is not joined to its
  parent. The tool stores the link `parent: {handle, arm, walk, overlap_bp}` with the child's store entry
  (`store.py:129`, `Entry.parent`). A saved child keeps only `J.derived_from` (the parent's handle), not the walk
  (`parser.py:1205`, `j_object`).
- Across pairs, a `/traverse` request runs on exactly one pair (`SPEC-multigraph-envelope.md` §1.1). A seed sent to
  every graph of a release gets one result per graph. A graph that lacks the seed answers with `outcome.walks:
  "not_in_graph"` and no graphlet, and a pair holding none of the seed labels with `"labels_not_in_graph"` (SPEC
  §6.1 step 2, step 3; envelope spec §2). The library's `TraverseClient.response()` keeps such results as `None`
  graphlets with an error entry (`client.py:895-911`). The service keeps one record per graph and seed, none merged,
  so a user reads r237's 11 records one by one.
- `branching.on_reconverge: "merge"` merges paths inside one traversal.
- After W5 (`E23/T6/W5/DELTA.md`), continuations, `deepen()` and replays keep the pair `(graph, graph_path,
  annotation_path)` a graphlet was answered on, the MCP wrapper reads the multi-graph envelope, and
  `GraphletView.claims()` keeps the view's arms. graphlets.rst lines 23–30 say "**Merging graphlets is not
  supported.**"

Combining graphlets needs rules for overlaps, orientation, index-local label ids, provenance and incomplete
evidence. Concatenating records does not give a graphlet.

### F.5 The service: stitching, and a note on wording

Item 4 of the list in §5.2 in full, and one more item for the service session's information.

4. **The service publishes no continuation tool today.** Stitched reads in the service would need one, plus the
   parent link (parent record, arm, walk) stored with the child record, plus the bundle for downloads. An agent can
   stitch locally: download both bodies, provided the child was fetched on the parent's pair with the walk's
   continuation as its seed. A child made from a preset plus overrides differs from `next_request()`'s request: on
   the default path the overrides take only `bounds`, `direction`, `exhaustive`, `labels.max_labels_per_node` and
   `labels.seed_label_kind`, so there is no reduced loss budget or branch allowance and no rebuilt `labels.extra`.
   Its junction then states `strategy_differs` or `branches_unequal` where that matters (§3.4); the stitch itself
   works.
5. **For information.** The library's `compare()` answers `comparable: false` for retrievals of different indexes
   and for different seeds. The service's `graphlet_compare` description gives "the same seed in two databases or
   two seeds of one locus" as examples, and both are such pairs. For the first, a union of the two databases'
   records answers what each database carries (§2). The wording is the service's.

### F.6 Risks in full

**The union.**

1. **A union read as one traversal of the release.** Mitigated by the five statements leading every answer, export
   and `print()`, `scope=within_one_graph` on every FASTA record, and the coverage sentence ending with the
   non-claim.
2. **A pair read as complete when it was not asked or failed.** Without every result and `expected_pairs=`, a graph
   that was not asked looks like a graph without the seed. Mitigated by the categories `failed`, `not_walked` and
   `not_included`, by naming the absent records not checked against the seed, and by "over the pairs given" on
   every statement when `expected_pairs=` is missing. The caller must pass them (§5.2).
3. **A turn-off read as absence.** Mitigated by `absence` on every span that is not `whole`, `extends` or `merges`,
   the third statement, and each pair's mode in the sentence, on the spans and in `modes=` (review M1).
4. **Equal names read as one record.** Mitigated by showing every label with its pair, refusing a bare name held by
   two pairs (`ambiguous_name`), and counting shared names in the labels statement.
5. **A header name on a `count_only` pair names another record.** The tag does not repair it. Mitigated by the
   per-pair warning the owner decided (OQ20, §2.5): header labels of `count_only` and `unknown` pairs are marked in
   the statement, carriers and exports. A library-wide warning for every graphlet stays a later option.
6. **Reviewed four times.** Fable reviewed the union at `42023784f` (Appendix G.2), the stubs and the input list at
   `c607ae6df` (Appendix G.3) and the whole note at `8e0083759` (Appendix G.5), and an external reviewer the
   owner-review edition (Appendix G.4). The fixes follow
   the reviewers' proposals where they made one; the implementation lane's review checks them against the code.
7. **A GFA reader following a link across graphs.** `0M` links concatenate segments whatever the tags say (review
   G.4, Must 2). Mitigated by a topology per mode: the default file is one tree per arm with the seed written twice
   and merges as tags, so it admits only paths one graph holds whole; the inferred file links only what the
   connection rule admits. P9 tests the paths each file admits.
8. **An inferred path read as a sequence some sample holds.** Mitigated by pieces and named transitions on every row and
   record, no `held_in` on a joined path and no `complete` on any inferred path, `scope=inferred`, and the fifth
   statement in that mode. The connection rule makes no overlap joins, so every transition is one a part states.

**Stitching.**

1. **A stitched walk read as one traversal.** Mitigated by the verdict on every row and junction, the per-walk
   `certified_to_bp`, the summary leading with the junctions, and `qualified` for everything in phase 1a-S, with
   completeness and evidence that follow the causes (§3.5).
2. **An `exact` that is not exact.** Three were found in review: siblings under `merge` (Appendix G.1, B1), a changed
   decision without a forbidden step (G.4, Must 1), and splits counted on the displayed chain instead of as the maximum
   over a merge's parents (G.5, B1). Mitigated by the narrow whitelist, by a W8 that compares every recorded decision of
   the child, by certification that stops at the first changed decision inside a child, by T4 pinning every entry
   against the uninterrupted retrieval (the counterexample included), and by the rule that a whitelist entry or
   promotion ships only with its passing fixture. A fourth gap was found in the GPT review: a decision cut from the
   events at a walk's very end (Appendix G.6, M2), now fixed. So the first release states no `exact` at all, and `exact`
   is switched on only after the seeded randomized comparison of T11 passes (OQ19, decided).
3. **Mostly `qualified` in practice.** In phase 1a-S everything is `qualified`, by the owner's decision (OQ19). Later,
   under the constrain defaults (`merge`, finite splits, `forbid`, no branches) `exact` needs w to be the arm's only
   open walk, and `seed_reentry` qualifies every walk that meets a repeat of the parent's walk. That is the honest
   answer. Promotions (`same_stop`, uniform-shift switching) narrow it once pinned.
4. **Two query families drift apart.** The composite delegates every derivation to the parts' `ops`, and T2 compares
   `stitch(g)` with `g` on every fixture for every query, so a change in `ops.py` that the composite misses fails T2.
5. **Scope creep.** compare, subgraph, routes and GFA on composites are requested next. They are listed as phase 2 or
   refused with a reason, never half-done.
6. **Memory with many children.** Bounded by `budget=`, stated by `memory_bytes()`.
7. **A wrong claim across a junction (1b only).** Joining a parent's claim to a switched descendant made a claim no
   label holds (review G.4, Must 4). Mitigated by joining only the same label's contiguous claims, keeping switches
   as they are, and stating lineage apart (`descends_from`); T10 pins it. A `qualified` verdict does not repair a
   wrong interval, so the join rule must be right on its own.
8. **A second file format (1c only).** The bundle is the library's (like `J`), versioned in its first line, and holds
   only unchanged MGT texts plus keys and digests. It is deferred until a user needs one file.

## Appendix G. Review history

Every finding was checked against the code and specs at `bc9684559` before it was settled. The section numbers in
the tables are this edition's.

### G.1 Fable, of `dfd97d9cf` (stitching)

| finding | checked against | settled |
|---|---|---|
| B1 `exact` unsafe under `merge` with another open walk | SPEC §6.5 (merge at one node and one distance; histories united), §6.4 consequences | **accepted, made stricter**: W6 and the cause `sibling_merge_unknown`. Besides "no other walk with a continuation", W6 also requires no other walk longer than L (one that ended deeper than a capped w could also meet its descendants). The old `scope` cause is replaced: merges on w's own route are covered by W8's check over all DAG ancestors (SPEC §6.6, "any parent"). T4(b) fixture added |
| M1 `lower_bound` not monotone | SPEC §6.4 (splits, ambiguity, greedy re-minimisation), §6.5 (merges) | **accepted**: `lower_bound` deferred; two verdicts (§3.4); M1's conditions recorded for a later phase |
| M2 no cause for `trace` | SPEC §6.3 (`trace_break`, one occurrence chain), §6.1 (seed validated on any consecutive occurrence) | **accepted**: W1 and `trace_restart` |
| M3 strategy differences not judged | the strategy fields of SPEC §5 | **accepted, list adjusted**: W2 with the reviewer's allowed list, the resource caps spelled out, `max_live_paths` only without a beam, and `labels.max_labels_per_node` added (it caps output lists only, SPEC §5); `strategy_differs` and `extra_foreign` added |
| M4 the continuation's length | `walker.cpp:3117-3128`, SPEC §7.1; the fixtures' (L, n) | **accepted, wording narrowed**: the formula is in §3.2, and §3.4 W8 uses which nodes the child's seed covers. "Usually the whole flank" holds only without a switch or a route join; the fixtures show both (`same_name` (40, 40); `merge` (100, 38), `caps` (126, 23), `switch_chain` (80, 20)) |
| minor 1 two sub-cases of seed re-entry | SPEC §6.6 precedence | **accepted**: `same_stop` and `earlier_stop` (§3.4); both `qualified` in 1a; `same_stop` is the first promotion candidate |
| minor 2 `max_switch_sources`, branch-count tie-breaks | SPEC §6.3 (tie-break, source cut), §6.5 (merge minimum) | **accepted, settled by the whitelist**: W3 requires `forbid`, under which no source list or switch predecessor exists, and W4 makes branch counts shift uniformly within each label |
| minor 3 causes need a mode column | `ops.py:2769` (annotate seeds carry no labels) | **accepted**: the modes column in §3.4 |
| minor 4 `pruned` omitted | SPEC §7.1; `is_resource_stop` includes `BEAM_PRUNED` (`traversal_types.cpp:59`) | **accepted**: Appendix D.1 |
| minor 5 joined loss is an upper bound | SPEC §6.3 | **accepted**: Appendix D.1 (0 under an `exact` junction; the chosen route's loss otherwise) |
| minor 6 `stitch(C, P)` code | §3.1's order | **accepted**: §3.7 and T3 |
| minor 7 inference vs continued candidates | — | **accepted**: §3.1, all candidates count |
| minor 8 reachability of `different_pair` / `different_graph` | SPEC §10.3 | **accepted as stated reachability, order kept**: a fingerprint mismatch is proof; pair tags can name one index twice (§3.1) |
| minor 9 CLI fixtures have no `index_fp` | `fork.mgt` `H` (`*`) | **accepted**: T4 asserts verdicts only; T5 checks `exact` on fixtures with an `index_fp` |
| nit: hashed edges for DNA k ≥ 32 | SPEC §6.6 | **accepted**: Appendix F.1 |
| nit: `continuation()` shift formula | `ops.py:2275-2279` | **accepted**: §3.2 |
| nit: tip and bubble windows | SPEC §6.5 | **accepted**: §3.4 |
| nit: `to_fasta(coordinates=)` | `model.py:571` | **accepted**: §3.1, Appendix D.2 |
| TQ1 causes incomplete | see B1, M2, M3, minor 2 | **accepted**; hairpins need no cause (§3.4) |
| TQ2 seed re-entry | see minor 1, M4 | **accepted**; trace restart and lost sibling merges are the other artifacts (both causes now) |
| TQ3 the rule check's window | — | **accepted with one change**: any parent's window; collisions ignored and documented (Appendix F.1). The reviewer read W8's ancestor set as qualified by `scope`; here it is the united DAG-ancestor set, and only siblings beyond L qualify (B1) |
| TQ4 the bundle | — | **accepted**: one file, `load()` dispatching, only in 1c |
| TQ5 claims order | — | **accepted** as written |
| TQ6 annotate claims shifted | — | **accepted**; `support_changes()` adds no label at the junction (Appendix D.1) |
| TQ7 identical child | — | **accepted**: a no-op (§3.1, §3.7) |
| scope: union of two arms has no use case | U2 | **accepted** in the second version (not offered); the owner later made it a validation helper (Appendix D.5, OQ2) |
| scope: two verdicts in 1a | M1 | **adopted** (§3.4, §7.1) |

The reviewer agreed with refusing cross-index composition (OQ3). The owner reversed that on 2026-11-15, and the
union was reviewed separately (Appendix G.2).

### G.2 Fable, of `42023784f` (the union)

Verdict "implement after fixes": the container of unchanged parts is right; one blocker and three majors let a union
read stronger than its parts. The reviewer's checks that came out sound were read again too: positions compare
across k and regime (SPEC §6.1), the seed check matches `_identity` (`ops.py:3600`), per-walk `complete` matches
`_walk_end` (`ops.py:1463`), and a segment's displayed labels do not depend on the walk (`ops.py:492-545`), which
removes the second version's **[unsure]** in Appendix F.1.

| finding | checked against | settled |
|---|---|---|
| B1 `complete` while the part's `outcome.walks` is `partial` | SPEC §7.0 (a `seed_labels` limitation and a `scope` with `observed` > 0 make the walks partial and leave `complete_to_bp`); SPEC §5 (`max_seed_labels` 1000) | **accepted**: `complete` needs the certificate and no walk-class limitation on the arm or the seed; `complete_qualified` names `seed_labels` or `scope`; `every_pair_certified_to_bp` is `None` with a qualified pair; P8 (§2.6) |
| M1 a turn-off read as absence | SPEC §6.10 (constrain follows only label-supported walks), §6.1 step 4 (each pair's own permitted set), §6.7 ("never evidence of absence") | **accepted, one case added**: `absence` on every `ends`, `diverges` and `none` span; `ends` with `dead_end` is `graph` in either mode, since the graph has no next k-mer at all (`END_CLASS`, `ops.py:634`); the third statement; each pair's mode in the sentence, on the spans and in `modes=`; P10 (§2.4) |
| M2 merge parents dropped, a false `diverges` | `derive.chain` (`derive.py:258`), `_first_paths` (`ops.py:641`) | **accepted**: stubs by default from the first-parent tree, `how: merges` with `into`. Every segment lies on a walk or a stub: a segment on no walk reaches, through children whose first parent it is, a segment that is only a non-first parent, and lies on that segment's chain. `routes=` and `segments_off_spellings` are dropped; P2, P5 (§2.4) |
| M3 the certificate rule joins nothing | envelope spec §9 (`headers_unique`: no repeat in a pair or across pairs); SPEC §10.3 (`columns_disjoint`) | **accepted**: `certificates=` and the name rules dropped; labels always shown with their pair; `label({'name'})` resolves for one pair or is refused `ambiguous_name`; `shared_names()`; the labels statement reworded; P7 replaced; OQ11 withdrawn (§2.5) |
| M4 cross-pair totals count one record twice | the note's own text | **accepted, and mostly moot** after the cuts (`support_by_depth()`, `walks(by='labels')`): no cross-pair total is offered; counts are per pair, as label entries, `+` marking a lower bound (§2.8, §2.9) |
| M5 absent records cannot be checked; the parts/absent split | SPEC §7.1 (the `seed` object: `seed_id`, `length_bp`, labels, no sequence); SPEC §4.3 (a seed id hashes the labels); `single_answer()` raises on a refused entry (`client.py:274`); `TraverseResponse.graphlets` holds `None` (`attempts.py:363`, `client.py:908`) | **accepted**: one list of results that the library sorts (§2.3); `{pair, graphlet}` states a pair (CLI bodies, databases without a manifest); `{pair, error}` for an answer that never came; `length_bp` and an optional `sequence` checked, other absent records named as not checked; `not_expected` covers absent records; `bad_absent`, `absent_has_part` (now `duplicate_pair`) and `absent_other_seed` removed; §5.2 item 2 names the `J` line |
| minor 1 a failed derivation is not "nothing known" | SPEC §6.1 steps 2 and 4, "Per-seed failure" | **accepted as one category**: `not_walked`, stating "the seed is in this graph" for a derivation cause; other failures without arms (a budget that does not hold the seed) state their cause only |
| minor 2 the incomplete reason | SPEC §6.10 ("two boundaries, one per question"); `Arm.cap_trigger` (`model.py:276`) | **accepted** (§2.6) |
| minor 3 the depth asked | — | **accepted**: the default is the largest radius among the parts (§2.6) |
| minor 4 the non-claim's wording | — | **accepted** everywhere: "a walk that no single graph holds entirely" |
| minor 5 the no-op | — | **accepted**: same pair and same body (§2.3, §2.10) |
| minor 6 the pair key | `_identity` (`ops.py:3588`: `index_ns`, release, `index_meta_fp`, `index_fp`) | **accepted** (§2.3) |
| minor 7 refusal order; unequal seed lengths | — | **accepted**: commutativity covers successful unions; a variant for unequal lengths (§2.3, §2.10) |
| minor 8 the certificates digest | — | **moot** under M3 |
| minor 9 `max_bytes` against "never cut" | — | **accepted**: `summary_too_small`, default `max_bytes=4096` (§2.8) |
| minor 10 P6 holds only in annotate mode | SPEC §6.1 step 4 | **accepted**: P6 split by mode (Appendix E.1) |
| minor 11 per-leaf insertion is quadratic on a comb | `walk_rows()` (`ops.py:1473`) | **accepted**: one depth-first pass over each first-parent tree (§2.4, Appendix F.1) |
| minor 12 memory counted one arm | DESIGN §7 | **accepted** (Appendix F.1) |
| minor 13 a union of absent records only | — | **accepted**: allowed with `seed=`; `no_seed` otherwise (§2.3) |
| minor 14 `count_only` is not a tagging matter | envelope spec §9.1; the library reads no certificate (no `headers_unique` or `seqs_order` in `api/python/metagraph/traverse/`) | **accepted**: stated in §2.5 as each pair's own matter that the tag does not repair; a library-wide gap for the owner (Appendix F.6; now OQ20, §8.3) |
| nits | — | **accepted**: every union walk is `whole` or `merges` in some pair (stubs extend the reviewer's "whole"); a root leaf is `ends` at 0, `none` an arm not walked; the column extension is gone; `every_pair_certified_to_bp` |
| the cut list | the owner's use case: one seed across r237's 11 graphs | **accepted in full** (§7.1). Two cut items were in the design brief (`support_by_depth()`, a union-level `compare()`), so OQ12 asks the owner to confirm |
| the views on OQ8–OQ11 | — | **adopted**: OQ8 with the renames, OQ9 with the mode statements, OQ10 with the narrower no-op, OQ11 withdrawn (§8) |

### G.3 Fable, of `c607ae6df` (the merge stubs and the input list)

A targeted review, "implement after fixes", no blocker: neither a stub nor the input list can make the union state
something false. Its checks that came out sound were read again: a merge consumes the heads of all its parents
(`walker.cpp:4765-4772`) and the parser refuses a segment that ends in both a split and a merge
(`parser.py:694-697`), so a non-first parent is a leaf of the first-parent tree and every segment lies on a walk or
a stub; all parents of a merge end at its `from_bp`; a stub never spells a walk's prefix in its own part.

| finding | checked against | settled |
|---|---|---|
| S1 a `merges` span stops at the merge though the part holds the rest | `Segment.parents` and `children` (`model.py:165-185`): the bases below a merge segment continue every parent's path | **accepted, pinned** for the owner's question (§2.4): a `merges` span carries a `route` matched through the part's segments and further merges (`held` with `goes_on`, `ends`, `diverges`, with `absence`, the merges crossed); `whole` stays displayed-only; `held_in` drives ranking, `pairs=` and the exports (`held=`); routes are computed on demand; P2, P5 |
| S2 span labels do not say which labels carry the displayed bases | `_labels_full` (`ops.py:1133`), `_displayed_profile` (`ops.py:492-500`), SPEC §7.1, `Segment.partition` (`model.py:176`) | **accepted**: `SpanLabel(label, from_bp)` on every span and route; `continues` on a `merges` span; a lineage not kept at a merge ends its history there; `complete` and `certified_to_bp` defined for `merges` spans through the route; `walks(labels=)` means carried whole (§2.4, §2.8). Since G.5, M2: the labels come from the alive sets, not the displayed profile |
| S3 a stub exported as a dead-end tip | the note's §2.4 and §2.9 | **accepted**: union segments cut where a merge segment starts; `into` is a union segment; GFA links the stub to it (`mg:Z:`); FASTA `merges=g6@85>s31` (§2.9) |
| S4 stub walks not marked | the note's `UnionWalk` | **accepted**: `UnionWalk.kind`; `stubs=False` by default in `walks()` and FASTA; counted apart in `summary()` |
| S5 `SegmentSupport.walks` quadratic on a comb | `walk_rows()` (`ops.py:1473`) | **accepted**: pieces and two counts; walks derived on demand; a stub named by its end segment (§2.4, Appendix F.1) |
| nits (stubs) | SPEC §7.1 (nested merges), `parser.py:694-697` | **accepted**: "into segment `into`", which can itself lie only on a stub; the parser cited; one merge per stub end asserted |
| R1 pair equality undefined | `_identity` (`ops.py:3592-3599`: a missing value matches) | **accepted**: one match rule for every pair comparison, and `ambiguous_pair` where a description matches two pairs, so the result never depends on the comparison order; a remedy for absent records in `duplicate_pair` (§2.3) |
| R2 the classification is not total | `TraverseResponse.errors` (`client.py:909`), `detail=` (`client.py:836`), `response()` (`client.py:897-911`) | **accepted**: a total decision table; `{pair, result, error}` read as the errors entry; `no_graphlet`, `result_has_graphlet`, `not_one_seed` (zero or several) (§2.3); P4 |
| R3 `seed_checked` only checks the caller's string | SPEC §4.3 (the seed id formula), §7.1 | **accepted**: `seed_check` `verified` / `stated` / `unchecked`, with a recomputed seed id where the record allows it (**[unsure]** which members a `not_in_graph` seed object carries); a sequence beside a response (§2.3) |
| R4 `seed=` and the part-less union | — | **accepted**: `seed=` checked against the parts; an arm-less coverage for a union without parts (§2.3, §2.6) |
| R5 the caller's statement not visible | — | **accepted**: `PairAbsent.pair_stated_by`; caller-stated pairs named in the sentence; the envelope-stated form wins a no-op, so the digest does not depend on order; `error` is a string or an exception's text (§2.3, §2.10) |
| R6 refusal order | — | **accepted**: arguments first (a string `expected_pairs` refused), then check-major order (§2.3) |
| R7 a single-graph server's missing seed is a 400 | SPEC §6.1 step 2 | **accepted as stated**: it arrives as `{pair, error}`, and `failed` quotes its text (§2.6). Classifying a 400 by its text would be fragile, so it is not done |
| nits (inputs) | `client.py:276`, `model.py:461` | **accepted**: the line corrected; §5.2 says `{pair, graphlet}` states a pair but does not replace the envelope |

### G.4 External review, 2026-11-15 (of `59a07a185`, the owner-review edition)

The reviewer's verdict: keep the architecture (unchanged source graphlets, explicit provenance, coverage per source)
but revise before implementation, because several guarantees did not hold as written, the stitching `exact` among
them. Six findings were marked Must, three Should, and one open choice went to the owner, who answered it ("Expose
both as separate modes") and restated an earlier decision (same-pair combinations designed, implementation
deferred). Every finding was checked against the code at `bc9684559`; all are accepted, and none is argued.

**The reproduction of Must 1** (lane dir `repro-m1/`, a throwaway script, not repository code; base binary
`metagraph_DNA-bc9684559`, one thread). Index: basic, k 15, one record `S P U C U Y` (30, 40, 15, 45, 15 and 30 bp)
in one column, so the 15-mer U is one node u with two successors, C[0] = C and Y[0] = G. Strategy: right arm,
constrain, `forbid`, `max_label_branches` 1, `on_reconverge: keep`, `output.continuation_bp` 40.

- **Parent** (radius 115): walk 0 is `P U C U`, 115 bp, ended by the radius at u, its label's one branch spent at
  the first u; walk 1 is `P U Y`, 85 bp, `dead_end`.
- **Child**, from `parent.next_request('right', [0])`: seed = the walk's last 40 bp, `max_label_branches` 0
  (`branch_budget`: original 1, effective 0). It walks nothing: walk 0 is 0 bp, label end `branch`, and its one
  branch event is `V 0 0 CG 1,1 0 0 C:branch:0;G:branch:0`. Its arm reads complete (`A r c 885`).
- **Reference** (radius 1000): walk 0 is `P U C U Y`, 145 bp, `dead_end`: at the second u the step `U+C` was used
  by this path, so only G is admissible and nothing is ambiguous.

Under the first version's W1–W8 the stitched walk is `exact` and complete at 115 bp; the reference holds 145 bp.
Under the revised W8, the branch event's successor C is a step in B (`U+C` is a 16-mer of the parent's walk), so the
walk is `qualified` with `candidate_changed` at 115 bp and certified to 115 bp. With a `continuation_bp` that covers
the loop (the default 1000), the child's seed holds the node after `U+C`, so the child blocks C as `rejoined_seed`;
the step is in B too, both sides block it and decide alike (by reading, not run). The counterexample needs a loop
longer than the continuation, which a walk past 1000 bp through a repeat gives at the default.

| finding | checked against | settled |
|---|---|---|
| Must 1 the whitelist can call unequal results `exact` | SPEC §6.3 (a successor is admissible only if its edge is unused on the path), §6.4 (ambiguity is counted over the admissible successors; a source over its allowance ends with `branch`), §6.6; `check_structure` (`walker.cpp:5908-5964`: a seed node or a used edge blocks a successor before ambiguity is counted); `next_request` lowers `max_label_branches` by the largest terminal branch count (`ops.py:3056-3071`); reproduced at `bc9684559` (above) | **accepted, changed**: W8 compares every successor the child recorded on the walk's route (`G`, `V`, `E`) with the steps the reference blocks, not only the steps taken; new walk causes `candidate_changed` and `decisions_unrecorded`; certification stops at the first walk cause inside a child (§3.4, §3.5, B.2); T4 `loop_branch` and a violation past a child's first base (E.2) |
| Must 2 GFA links admit cross-source paths | §2.9 of `59a07a185` (a stub linked to the union segment its merge segment starts); GFA1 `L` lines join segments for every reader. Checking it found a second class: one seed segment joins any pair's left flank to any pair's right flank | **accepted, changed**: a topology per path mode (§2.9). Default: one tree per arm, the seed once per arm, a merge as tags (`mg`, `mt`). Inferred: stub links and one seed, exactly the connection rule. P9 tests the paths each file admits (OQ16) |
| Must 3 identical bodies are not identical retrievals | `store.py:3-7` (entry and body identity apart: "Two requests with different loss budgets can produce byte-identical bodies"), `request_of()` (`store.py:1277`) | **accepted, changed**: a no-op needs an equal body and an equal retrieval identity (index identity, normalized strategy, seed as requested); otherwise `duplicate_pair` or `already_continued`, each with a new text; identical bodies are shared for storage only (§2.3, §3.1, A.2, B.1, C.3) |
| Must 4 D.1 joins switched descendants | D.1 of `59a07a185` ("joined with each claim of the child's lineages descending from X's seed run") | **accepted, changed** (phase 1b): only the same label's contiguous claims join; a switched-in claim stays its label's with its switch; lineage is `descends_from`, apart from the claims (D.1); T10 |
| Must 5 a caller's `seed_id` cannot verify an absent record | `failed_seed_head` (`src/cli/traverse.cpp:4859-4867`: the request's `seed_id` echoed, `length_bp`, `labels_from_seed`, nothing else), `not_in_graph_seed_to_json` (`:5345`); `validated_seed_id` is computed after the labels (`walker.cpp:1168`; SPEC §6.1 step 3); `TraverseResponse` keeps no request (`attempts.py:363`) | **accepted, changed**: `seed_check` is `stated` or `unchecked`, `verified` reserved for a server field or a response that keeps its request (§2.3). This resolves the former **[unsure]** of §2.3: the seed object carries none of the members a recomputation needs |
| Must 6 pair tags can collide | §2.8 of `59a07a185` (tags told apart by an `index_fp` prefix, which two pairs without manifests lack) | **accepted, changed**: `pair_id` from the canonical pair key, unique in the union, keys every answer and selector; display tags are made unique with a `pair_id` suffix (§2.3, §2.8); P4. Since G.5, minor 1: the unique tag is the one identifier and `pair_id` is gone (Appendix A.4) |
| Should: shorter source results | §2.4 of `59a07a185` (a union walk is a leaf of the trie) | **accepted, changed**: terminal markers, so every part's walk end is a union walk; the span value `extends` for a pair that holds the walk and goes on; `held_in` includes it (§2.4); P5 |
| Should: support changes within a segment | the `reentry` fixture (labels A, B, A in one segment; `traverse_testlib.py`) | **accepted, changed**: interval support (`LabelSpan`, since G.5 `LabelInterval`) rather than splitting segments, so union segments stay defined by spelling and pairs (§2.4, C.1); P5 |
| Should: finish the API contracts | `ops.label` (`ops.py:238-298`: `{'name'}` raises `AmbiguousLabel` for two labels of one name in one graphlet); `budget.py` (a stop raises `LocalBudgetExceeded`; a `Partial` only for lists whose order allows it, with no `complete` member); `pct_escape` (`_codec.py:263`, `%` and line breaks only) | **accepted, changed**: `{pair, name}` resolves as the part's selector does; tokens percent-escaped and reversible; `budget=` on every method with the library's semantics, and the "`Partial` with `complete: false`" of §3.5 withdrawn (§2.5, §2.8, §2.9, §3.5) |
| Open choice: source-supported or inferred paths | the owner: "Expose both as separate modes"; the reviewer's three rules and "connections justified by the input graphlets only" | **designed** (§2.7, §2.9): two modes, source-supported by default; the first version's connection rule (a pair's own merge, the seed; no overlap assembly); pieces, named transitions, narrow completeness; OQ13 provisional, OQ15–OQ17 open |
| Scope: same-pair combinations | the owner's decision, restated | **designed** (§4.2): `union_retrievals()`, built later; the "not offered" statements changed (§0, §1.2, §4, A.2, A.3, F.3); OQ14 provisional |

### G.5 Fable, of `8e0083759` (the full review)

The reviewer's verdict: ready after fixes. The architecture holds; one stitching guarantee was false as written (B1)
and six union items were wrong on real bodies or left open (M1–M6), with a scope point (M7), eight minors, six nits,
five readability items and views on the open questions. The reviewer ran the reproductions on another worktree's
build whose `walker.cpp` is byte-identical to `bc9684559`. I re-ran every script this table relies on with the clean
`metagraph_DNA-bc9684559`, one thread, on the reviewer's fixtures (lane dir `fable-rerun/`, outputs in `out/`;
throwaway scripts, not repository code); every result matched the review. Every other finding was checked in the code
or in the note itself. The owner's decisions stand: where the reviewer proposed otherwise, the proposal is an owner
question in §8.3, not an edit.

| finding | checked against | settled |
|---|---|---|
| B1 W8 counts splits on the displayed chain, the engine the maximum over a merge's parents: a false `exact` | `walker.cpp:4774` and `4875` (a merge keeps `std::max` of its parents' counts), `5394` (a split is refused at `item.splits + 1 > max_splits_per_path`). Re-run of `splits45.py` and `splits_check.py` (k 11, constrain, `forbid`, `merge`, `max_label_branches` 0, `max_splits_per_path` 2): w, 45 bp, has one split on its displayed chain [0, 2, 5] and count 2 at its leaf (segment 5 merges 2 and 3); W1–W7 hold, the other walk is an 18 bp `dead_end`; W8 as written passes (1 + 1 ≤ 2, no step in B, no events); stitched walks 85 and 85 bp; reference walks 18 bp and 50 bp, the latter ended `branch` with the qualifier `split_limit` at 50 | **verified, fixed**: W8 counts splits as the engine does, in parent and child (0 at the root, one more on each child of a split, the maximum over a merge's parents), and requires P + c + 1 ≤ `max_splits_per_path` at every split the child made on the walk's route (§3.4); B.2 restated; T4(b) fixture `splits45`. The wider point, gating `exact` on a seeded random oracle, is put to the owner as OQ19 (§8.3), since stitching in 1a is the owner's provisional answer (OQ1) |
| M1 `absence: graph` and "`dead_end` means no next k-mer" are false on canonical and primary graphs and under mixed k | SPEC §7.0 (a successor that was only a skipped hairpin is a `dead_end` with its `hairpin` event), `walker.cpp:6161` (annotate), `5644` (constrain, `Dh`). Re-runs: `hairpin.py` (k 11, canonical, annotate, `hairpins: skip`: the walk ends at 0 bp with `dead_end` and a `hairpin` event, inside `complete_to_bp` 30, `per_path`; with `follow` it ends `rejoined_seed`); `hairpin_b.py` (the basic graph of the same record spells 30 bp past the seed, all of whose k-mers the canonical graph holds); `hairpin_div.py` (`follow` has two walks that `skip` lacks, both `per_path`, `complete_to_bp` 20); `more.py` part 2 (k 11 ends by `edge_reuse` at 56 bp; k 21 spells the record's 90 bp) | **verified, fixed**: the absence value `walk_rule` (edge reuse, seed re-entry, a `dead_end` on a skipped hairpin, or a turn-off at a base with a `blocked` or `hairpin` event); `graph` narrowed to a true tip, or an event-free turn-off below `evidence.complete_to_bp` in annotate mode; the span states the pair's k and regime; precedence stated: `graph` and `walk_rule`, then `uncertified`, then `under_strategy`, so a `dead_end` at the certificate's edge is `graph` (§2.4, §2.7); P10. Not reachable on r237 (every graph basic, k 31); OQ9 allows mixed settings, so the rule is needed |
| M2 span and segment labels from `_displayed_profile` are wrong at cut depths | `derive.py:804-805` (the profile is anchored at the leaf); `ops.py:437` `alive_at`, `ops.py:564` `_alive_pieces`. Re-run of `labels_cut.py` on the committed `merge.mgt`: at 40 bp of walk 0 the alive set is {`a.fa`, `b.fa`, `both.fa`, `c.fa`} and the displayed profile {`b.fa`, `both.fa`, `c.fa`}; from the partitions, `a.fa` carries from 0 and `b.fa` from 36; on the stub [0, 2], which ends at 36 and merges into segment 3, `b.fa` continues, while the profile at the stub's end holds only `both.fa` | **verified, fixed**: segment labels are `LabelInterval`s from the alive sets and runs, read through the library's existing `_alive_pieces` (annotate: `presence`); `SpanLabel.from_bp` is the latest merge on the spelled chain, at or before the span's end, that the label did not enter through that chain's parent, or the start of its run for a label switched in later (§2.4, F.1, F.2); P5 |
| M3 the coverage table leaves a truncated part uncategorised | §2.6 of `8e0083759`: `complete` required "no walk-class limitation on the arm", which a truncated arm has (a `walk_domain` limitation), and `complete_qualified` named only `seed_labels` and `scope`, so a part with `complete_to_bp` 420 asked at 300 fitted no row | **verified, fixed**: `complete` at `at_bp` means the certificate reaches it, with no `seed_labels` limitation and no `scope` limitation with `observed` > 0, whatever cut the arm beyond; without `at_bp` every view, `by_category` included, uses the largest radius; `every_pair_certified_to_bp` is `None` without `expected_pairs=` (§2.6, C.1); P8 |
| M4 zero-length merge segments | `limits.mgt` line 21 (`G 1,2 7 0`). Re-runs: `zero_merge.py` (a merge at the radius is a zero-length leaf: `G 1,2 36 0` at radius 36, `G 5,4 62 0` at radius 62); `more.py` part 1 (`G 1,2 23 0` with two children of 17 bp) | **verified, fixed**: every merge is at a trie node, the start of a union segment, or the end of one (`s30$`) for a zero-length merge segment; `into` names that trie node; the `merges` token and `mt:Z:` carry it; the inferred GFA links the stub to each child of a zero-length merge (§2.4, §2.9, A.5); P2, P5 |
| M5 inferred mode leaves pieces, `overlap_bp` and path ids open | §2.2 and §2.7 of `8e0083759`: the example gave g01 alone on [85, 860) although g06's route holds [85, 812); `overlap_bp` was 30 in the example and "at least P's merge node" in the text; `spell(arm, 'i<n>')` depended on the listing | **verified, fixed** (the owner's instruction: fix the open data in the design): pieces are cut wherever the set of pairs holding the path changes, P's own route included (the example now reads g06 on [0, 55), g01 and g06 on [55, 700), g01 on [700, 860)); a `Transition` carries `at_bp`, `from_bp` and `to_bp` in place of `overlap_bp`; ids are `i:` and the union segments (`i:4.31.52`), independent of the listing and of `max_transitions` (§2.2, §2.7, §2.9, C.1); P13 |
| M6 "every answer carries the five statements" is not implementable | `mcp_tools.py:1352` (the tool layer attaches the evidence block today); §2.8 of `8e0083759` (list queries return plain lists and strings); C.1 had no `Statements` | **verified, fixed**: the claim is narrowed to summaries, `print()`, exports, `coverage()`, `evidence()` and `statements()`; list queries stay plain lists; `Statements` and `Carriers` are defined (§0, §2.1, §2.8, C.1) |
| M7 scope for the r237 case: the `J` line, a fan-out helper, inferred mode | §5.2 of `8e0083759` (the `J` line unconfirmed); `client.py:836` (`traverse()` asks one pair); `attempts.py:363` | **partly**: the `J` line is a precondition the service session confirms, stated in §5.2 and raised as OQ21. **Fixed**: `TraverseClient.traverse_pairs()` asks every pair, wraps a failure as `{pair, error}` and returns `PairResults` with the pairs asked as `expected_pairs` (§5.2, C.1, F.2). **Not applied**: "ship source-supported first and gate inferred separately", since both modes are the owner's decision (OQ13); the phasing is put to the owner as OQ18, recommending (b), with the design of the inferred mode kept |
| minor 1 `pair_id` uneven and unstable | §2.2, §2.8, §2.10 and C.3 of `8e0083759` (examples, `evidence()` and `provenance()` keyed by tag) | **verified, fixed**: `pair_id` dropped; the display tag, made unique with a `~` key-digest suffix only on a collision, is the one identifier, and each part keeps its pair key (§2.3, A.4, C.1); P4. G.4's Must 6 row notes it |
| minor 2 the match rule ignores `index_meta_fp` | A.4; `ops.py:3593` (a `meta_fp` difference is proof) | **verified, fixed**: a differing `index_meta_fp` makes two pairs (A.4 rule 1), and rule 3 compares it where both state it |
| minor 3 retrieval identity "as `request_of()` builds it" | `store.py:1289-1301` (`None` without an envelope or with unverifiable label names; adds the caller's `seed_id`) | **verified, fixed**: retrieval identity is defined in A.4 without `request_of()` and without `seed_id`; two identical absent records are one; a failed record followed by an answer keeps the answer, the failure in `provenance()`; the envelope-stated copy survives a no-op, else the first in canonical order; P4 |
| minor 4 walk ids "leaves in trie order" | §2.4 (a union walk need not be a leaf) | **verified, fixed**: ids number the terminal markers in trie preorder, walks and stub walks in one sequence (§2.4) |
| minor 5 `PairSpan.complete` undefined for four `how`s | §2.4 | **verified, fixed**: `complete` and `certified_to_bp` defined for every `how` (§2.4, C.1) |
| minor 6 §5.1 says the union reads `GET /capabilities` | §5.1 against §2.5 | **verified, fixed**: the union reads bodies and envelopes only; the fan-out helper reads `GET /capabilities` (§5.1, §5.2) |
| minor 7 seed checks | `attempts.py:363` (`TraverseResponse.sent` keeps the attempt fields only), `cli/traverse.cpp:4859` (`failed_seed_head` sends no sequence) | **verified, fixed**: the client keeps the seeds it sent (`TraverseResponse.seeds_sent`, a client change in 1a-U), giving `seed_check: 'sent'`; the case mapping is stated; the input form `{response, index}` keeps `pair_stated_by: envelope` for multi-seed responses (§2.3, A.1, A.2); P4 |
| minor 8 tests | E.1 P2 and P6 of `8e0083759` | **verified, fixed**: P2 runs on fixtures with bases (others are refused `no_bases`), and a stub whose merge is a zero-length leaf is `held` without `goes_on`; P6 runs under `keep` and under `merge` |
| nits | §2.6, §2.9, §2.2, §2.8, F.3 of `8e0083759` | **verified, fixed**: g8 in the coverage sentence and g5 in the FASTA and GFA comments; union segments written `s<id>` and trie nodes `s31` or `s30$` throughout; docs 300 lines in §7.2 and F.3; the §2.8 comment lists every input form; `LabelSpan` renamed `LabelInterval`; the carriers listing is `print(u.carriers(...))` (`Carriers.__str__`) |
| readability | the body of `8e0083759` | **done**: the inline review tags are gone from §0–§8 (the verdicts are here); §6 and §7.3 are pointers to E and F.6; the pair-match and no-op rules moved to A.4, the token grammar to A.5, the old §4.5 to D.5; §8.1 is a table; `duplicate_pair` is defined once (§2.3) and `union_retrievals` once (§4.2), other mentions being the refusal list, the names list of OQ17 and Appendix D.4; the old §8.4 became OQ20 |
| OQ12, OQ15–OQ17 (the reviewer's views) | §8.2, §8.3 | **agreed**: OQ12 leave both out, support by depth the first follow-up (now possible with M2's alive sets); OQ15 accept, M5 being fixed; OQ16 accept; OQ17 accept, with `LabelInterval`, `pair_id` gone and `walk_rule` in the list. The recommendations say so |
| OQ1, OQ9, OQ13, the old §8.4 (risky provisional answers) | §8.1 | **raised, not edited**: OQ1's "`exact` off until a seeded random oracle exists" is OQ19; OQ9 stands, and M1's `walk_rule` states what mixed settings change; OQ13's "not inside the first lane" is OQ18; the old §8.4 note on `count_only` is now a decision, OQ20 |
| "checked and sound" (cited lines, refusal orders, the continuation formula, the `fork` example, `failed_seed_head`, trie invariants on 27 graphlets, the default GFA topology, the label non-join) | — | no change needed |

### G.6 GPT, of `8486498fe`, with the service's answer to OQ21 and the owner's decisions on OQ18 and OQ19

The reviewer's verdict: ready after fixes; no blocker, five majors, five minors, two nits, views on the open questions
and a cut list. It ran a MetaGraph 0.5.2 binary with the bundled parser; its files are in the lane dir
(`gpt-repro/`). I checked every finding again: M1, M2, M3 on the clean `metagraph_DNA-bc9684559`, one thread, on tiny
fixtures (lane dir `gpt-rerun/`; throwaway scripts, not repository code), M4 and M5 against the library, the rest in
the specs, the code and the reviewer's files. The owner decided OQ18 and OQ19 the same evening (§8.1); the service
session answered OQ21.

| finding | checked against | settled |
|---|---|---|
| M1 absence evidence must match the successor asked about | `model.py:122` (`Event.char`, `followed`). `m1.py`: basic, k 2, seed `AA`; graph A holds `AACAAT`, graph B `AACAG`. At outward depth 2 (prefix `CA`) A records `E 2 b A J` (the successor A blocked, `rejoined_seed`) and goes on with T; B goes on with G, and A has no k-mer `AG`. The rule of `8486498fe` read any `blocked` event at the turn-off as `walk_rule` | **verified, fixed**: absence is decided for c, the union walk's next base where the pair stops. `walk_rule` needs a `blocked` event or a skipped `hairpin` event (`followed` false) whose character is c; `graph` needs no event for c and a complete record there (§2.4); P10 gets this case |
| M2 W8 misses an unrecorded decision at the walk's end | SPEC §7.0 ("guaranteed its recorded reason when d < `evidence.complete_to_bp`; at or beyond it the reason may be among the events not kept"). `m2_m3.py`: the `loop_branch` child of Appendix G.4 run with `max_branch_events` 0 gives `evidence {complete: false, complete_to_bp: 0}`, a `K r branch_events` limitation, no `V` record, a 0 bp walk ended `branch`. Its end never passes the boundary, so `decisions_unrecorded` as written did not fire and W8 saw no event | **verified, fixed**: W8's checks are complete strictly below the boundary, and a walk whose end lies at or beyond it is `decisions_unrecorded`, a 0 bp walk with a boundary of 0 included (§3.4, B.2); T4 fixtures for that child and for a boundary equal to a walk's end |
| M3 a detected cause can still give `complete: true` and exact evidence | §3.5 of `8486498fe`. The same child's walk is 0 bp and `complete` in its part (`m2_m3.py`), so the stitched walk had `certified_to_bp` 115 = its length, `beyond_certified_bp` 0 and `complete: true`, while the reference reaches 145 bp; the evidence's `exact` read only parts and junctions | **verified, fixed**: `complete` also needs no walk-scope cause on the walk, its end included; a walk-scope cause maps to `walks: partial` and `label_evidence: qualified` and bars `exact` (§3.5). By the owner's decision (OQ19) this is a requirement of the qualified release: phase 1a-S states no `exact`, certifies a stitched walk at most to its first junction, and never shows `complete: true` or exact evidence beside an unjustified `qualified`. T4 and T5 assert `complete`, `certified_to_bp`, the evidence and the outcome beside each verdict |
| M4 zero-transition inferred paths differ from supported walks and can grow exponentially | `m4_m5.py` on the committed `merge.mgt`: 7 segments, two bubbles in a row (1 and 2 into 3; 4 and 5 into 6), one displayed walk, four routes from the seed to the leaf | **verified, fixed**: inferred paths have a kind (`walk`, a union walk, in both modes; `route`, held whole by its pair through its own merges, with `held_in`; `joined`, with a transition); the fifth statement is narrowed to joined paths; `max_paths` (default 1000) bounds the enumeration whatever `max_transitions` is, in path-id order, and the listing states its count and any cut (§2.7, §2.8, A.5, C.1, F.1); P13. The inferred mode is a later phase (OQ18, decided) |
| M5 label intervals from `_alive_pieces` are not linear | `ops.py:564-590` (a `frozenset` of the whole alive set per piece). `m4_m5.py` runs the library's `_iter_alive_pieces` on one segment whose N entry labels end at N bases: 55, 5,050 and 500,500 label ids copied for N 10, 100 and 1,000, N(N+1)/2. The proposal to read intervals from `_alive_pieces` was mine, in the G.5 fix of M2 | **verified, fixed**: intervals are built from the entry set and `derive.label_changes()` (`derive.py:72`), each change opening or closing one interval, O(entry labels + changes) per segment; annotate mode from the `presence` runs (§2.4, F.1, F.2) |
| minor 1 the helper's `**kw` forwards `attempt_id` | envelope spec §1.1 (a repeated `attempt_id` is a 409) and §5 (request and attempt ids are one namespace) | **verified, fixed**: an explicit keyword contract; `attempt_id` and `expect_server_instance` are not accepted; a ledger-managed fan-out with a distinct id per pair and each receipt kept is a later extension (§5.2) |
| minor 2 `expected_pairs` typed as names | C.1 (`tuple[str, ...]`) against §2.3 and A.2 (`expected_pairs=` takes pair objects); two annotations under one graph name | **verified, fixed**: `PairResults.expected_pairs` holds the full pair objects sent (§5.2, C.1) |
| minor 3 non-ASCII labels make invalid GFA | A.5 of `8486498fe` escaped only `%`, whitespace, control and token characters; GFA1 allows printable ASCII only in names and `Z` values | **verified, fixed**: every byte outside printable ASCII is percent-escaped, in GFA and FASTA alike (`café` → `caf%C3%A9`); P9 |
| minor 4 reverse stitching is not always refused | the reviewer's `stitch-reverse-review` files: P (seed `GGCCCCGCCGCAACT`, k 15, radius 25, `continuation_bp` 15) continues into `TGATTGGTCGTTAAT`, which is C's seed, and C continues into P's seed, on a circular record | **verified, fixed**: stitching is described as not commutative between parent and child; `stitch(C, P)` is refused in the usual cases and judged by W1–W8 in this one (§3.7, C.3); T3 gets the circular fixture. That both orders pass every check is the reviewer's reasoning and mine, not run |
| minor 5 `Transition.by` names one of several witnesses | C.1 (`by: str`) | **verified, fixed**: `by` lists every pair whose merge at that node permits the connection, in canonical order; the FASTA token lists them (§2.7, A.5, C.1); P13 |
| nit 1 W8 is one paragraph of five checks | §3.4 | **done**: W8 is five sub-items: what the reference blocks, admitted successors, blocked successors, split counts, the evidence boundary |
| nit 2 the helper paragraph mixes four things | §5.2 | **done**: a signature, five short items and an example |
| views on OQ12, OQ15–OQ17, OQ20 | §8.2, §8.3 | **recorded as supporting views**, the recommendations unchanged: OQ16's file now says in a comment line that no left–seed–right path is in it; OQ17's "supported" is defined in §2.7 (a source graph's recorded traversal, not one biological record); OQ20's option (b) gains the `unknown` state and the service's exports |
| views on OQ18, OQ19 | §8.1 | **decided by the owner**, 2026-11-15: OQ18 (b), as GPT and Fable recommended; OQ19 (b), stitching qualified first and `exact` after a seeded randomized comparison (T11), which covers what GPT asked of the oracle (endpoints, completeness, evidence and outcomes, event cutoffs, chains of three or more parts, minimized failures kept) |
| OQ21: confirm the `J` line by a real round trip | the service session's answer (lane dir `OQ21-service-answer.md`); `parser.py:228` (the index fingerprints come from the `H` line); `parser.py:1358-1366` (`standalone_text` copies `H` and the body and inserts `J`) | **answered for r237, with a first step**: the service adds the `J` line and both exports keep it; multi-graph records carry their own pair under a reduced envelope; master-backend records carry no pair fields and are passed as `{pair, graphlet}`; the first step of 1a-U is an export → `load()` → `union_pairs()` round trip on one r237 record (§5.2, §7.1, §8.3) |
| cut: mandatory sharing of parsed models | A.4 and P1 of `8486498fe` ("identical bodies are parsed and held once"); results arrive already parsed, and a model carries its envelope and caches | **agreed, cut**: each part keeps the model it was given, and no model is shared or merged between parts; sharing one parsed body is a later optimization (A.4, P1). Holding the caller's parts by reference, never copied, stays; it is not this optimization |
| cut: inferred enumeration and same-pair combination later | §7.1 | **agreed**: the inferred mode is a later phase (OQ18, decided); the same-pair combination stays designed and deferred (OQ14) |
| "checked and sound" (the default GFA topology, terminal markers and `extends`, interval support, pair-qualified labels, the revised W8 on the recorded example, the split recurrence, the root's arm-wide completeness, the same-pair decision) | — | no change needed |

### G.7 GPT, follow-up of `cd6a4245a`, with the owner's decision on OQ20

The reviewer's verdict: ready after minor fixes; no blocker or major, five minors and one nit, all reasoned from the
note (nothing run). Each was checked against the note and, where it touches them, the code and the specs; none needed
a run. The owner's paste of this review also carried the decision on OQ20, recorded in §8.1.

| finding | checked against | settled |
|---|---|---|
| minor 1 the fan-out example fails when no graph holds the seed | §2.3 (`seed=` is required when no result carries a graphlet), A.2 (`no_seed`); §5.2 of `cd6a4245a` called `union_pairs(r.results, expected_pairs=...)` without `seed=` | **verified, fixed**: `PairResults` carries the seed sent, and the example passes `seed=r.seed` (§5.2, C.1); P4 builds the arm-less union of an all-absent fan-out and checks `no_seed` without `seed=` |
| minor 2 qualification should apply only to answers that read a junction | §3.5 of `cd6a4245a` ("no stitched answer's evidence reads `exact`") against the identity element (C.3: `stitch(g)` answers as `g`) and T2 | **verified, fixed**: in phase 1a-S no answer that depends on a junction reads `exact`; `stitch(g)` without a child and root walks no child continues keep the root's own evidence; the identity condition applies where a junction is read (§3.5); T5 |
| minor 3 T11 must compare the prospective exact result, not the public one | T11 and T4 of `cd6a4245a`: the public 1a-S result certifies a walk only to its first junction, so a correct chain of a 100 bp parent and a 100 bp child reads 100 bp and incomplete while the uninterrupted traversal reads 200 bp and complete; the comparison as written would fail it | **verified, fixed**: §3.5 defines the prospective result (`exact_gated` read as `exact`), a private test hook in 1a-S; T11 compares it with the uninterrupted traversal and checks separately that the public result only reduces it; T4 asserts both, its interior-certification fixture on the prospective result |
| minor 4 the inferred listing needs a truncation contract and a ranking rule | §2.7 of `cd6a4245a` (a count and a cut promised, `walks()` returning a plain list, no FASTA token) | **verified, fixed for the later inferred phase**: `InferredListing` with `count`, `cut` and `resume`; `top=` and `by=` rank only the paths enumerated, within the path-id prefix when cut, never a global top; FASTA `of=<count>`, `+` when cut (§2.7, §2.8, A.5, C.1); P13 |
| minor 5 path claims contradict the walk/route distinction | §2.1 of `cd6a4245a` (the old fifth statement, "no graph is shown to hold them whole"), §2.7 ("it is in the union"), P13 (four results as `route`); `merge.mgt` has one displayed walk among its four routes (`gpt-rerun/m4_m5.py`) | **verified, fixed**: §2.1 quotes the fifth statement of §2.7; a shared sequence's bases are in the union, as walks where a pair displays them and as stub walks otherwise, and the supported mode does not list every route (two parts displaying the same allele list the other only as a stub walk); P13 expects one `walk` and three `route`s |
| nit availability descriptions | §2.8 ("everything is phase 1a-U except `[1c]`"), F.2 (`union.py`, `export.py`), F.3 (the graphlets.rst intro and sections) | **verified, fixed**: `[inferred]` marks in the API; the inferred classes, records and links listed under the later phase in F.2; the docs plan says inferred paths and stitching's `exact` come in a later release |
| OQ20 decided | the owner's paste: "add per-pair warnings for unchecked header-to-record mappings, using `unknown` when metadata is unavailable"; envelope spec §9.1 (`seqs_order`: `verified`, `count_only`, null; in `GET /capabilities`, not in `/traverse` answers) | **recorded as decided, 2026-11-15**, from the paste, for the owner to correct: option (b) is a requirement of 1a-U. `union_pairs(..., seqs_order=)` (the helper fills it; a pair not given reads `unknown`), refusal `bad_seqs_order`, `PairPart.seqs_order`, warnings on header labels in the labels statement, `carriers()`, FASTA `mapping=` and a GFA comment (§2.5, §2.8, §7.1, §8.1, A.2, A.5, C.1); P7. `seqs_order=` is the one input that comes back from the cut `certificates=` |
| views on OQ12, OQ15–OQ17 | §8.2, §8.3 | unchanged (defer; accept, accept, accept), as recorded |
| "checked and sound" (per-character absence, terminal decisions at the cutoff, walk-scope causes in completeness, linear intervals, zero-transition routes, every merge witness, ASCII escaping, full pair selectors, the keyword contract, no model sharing) | — | no change needed |
