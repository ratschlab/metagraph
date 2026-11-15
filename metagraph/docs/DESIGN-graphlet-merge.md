# Design: combining graphlets — the cross-pair union, stitching a continuation, and what is never merged

**Status.** Design note only (train 7, lane G1; base `bc9684559`, train 6's head with W5's `93044c6ee`). Nothing here
is implemented. The owner chose "design note first" on 2026-11-14 after the Python graphlet review
(`review-python-graphlets-6d97d6e8.md`: "no public merge/union/stitch of graphlets"; its suggested order ends with
"define a merge API separating same-index union, continuation stitching and cross-index composition"). The first
version (`dfd97d9cf`) was reviewed by Fable ("implement after fixes": one blocker, four majors, nine minors). The
second (`77243cc05`) settled every finding (§20.1). The third (`42023784f`) took the owner's answers of
2026-11-15: the cross-pair union, which the second version refused, is offered and ships first (OQ3: "why would I
refuse it? It's the only way I can do it today", and "it won't be exact, I know, but the best we can do today, I
think"), and the union of the two arms of one seed becomes a later validation helper (OQ2). Fable reviewed the union
("implement after fixes": one blocker, five majors, fourteen minors, and a cut list for 1a-U). This fourth version
settles every finding (§20.2), takes the cuts, and records the owner's answer to OQ1 ("Phase 1a": stitching is
built, as 1a-S after the union). Items marked **[unsure]** are my reading and need a check. §19 lists the open
questions.

All paths are relative to `metagraph/`; library paths are under `api/python/metagraph/traverse/`. Line numbers are
`bc9684559`'s.

## 0. The decision in brief

1. **Two combinations are offered.**
   - **The cross-pair union**, `union_pairs(results)` (§5). It combines the graphlets of one seed from several
     pairs (r237's 11 graphs, or two databases) by sequence: one trie of the parts' spellings per arm, with the
     pairs and labels that carry each segment. It is the best combination available today, and every answer says
     what it holds (§5.1): how complete each pair is at the depth asked, which pairs carry each segment and walk
     and how far, that a pair not following a walk did not follow it under its own request, that labels are each
     pair's own, and that a walk no single graph holds entirely is in no part, so not in the union.
   - **Stitching a continuation**, `stitch(parent, child)` (§6–§8). It joins a continuation to the walk it
     continues, on one pair, and states at each junction whether the joined walks are those of one uninterrupted
     walk.
2. **The union ships first, then stitching** (§18; the owner approved phase 1a, OQ1). The union has a named user
   and a need today (the owner, on r237). It states nothing beyond what its parts state, so no verdict waits for an
   oracle. It also builds the composite machinery (parts, refusals, evidence, budgets, FASTA over pieces) that
   stitching then reuses. Phase 1a-U is cut to what looking at one seed across r237's 11 graphs needs (§18).
3. **Neither result is a graphlet.** A `PairUnion` and a `StitchedGraphlet` hold their parts' bodies unchanged, and
   every answer is computed from the parts. MGT v1 is frozen, and per-traversal records (A counters, B growth bins,
   V branch events, one Q, seed-first label ids) cannot be merged without inventing values. This is the view's rule
   (a view saves its backing body unchanged) applied to several bodies (§3.5).
4. **Every stitch junction has one of two verdicts.** `exact` means the child's walks below the junction are the
   walks one uninterrupted walk of the parent's request would have. It is given only on a narrow whitelist of
   conditions (§7.4). Otherwise the verdict is `qualified`, with the named causes. A `lower_bound` verdict is
   deferred, because fewer labels or less budget can also lengthen a walk (§7.4, review M1). No answer of either
   combination reads more complete than the weakest part or junction it depends on.
5. **Not offered.** A union of two retrievals of one pair: they answer different requests, and `compare()` reads
   them (§3.2). Combinations of different seeds (§3.4). The union of the two arms of one seed is a later
   validation helper, not a combination (§10).
6. **No engine change, no MCP tool.** The name `merge` is not used: it already means reconvergence
   (`on_reconverge: merge`, `m` runs, merge segments).

## 1. Problem and today

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

## 2. Use cases: what a user combines, and whether it needs a new operation

Sources: graphlets.rst, `DESIGN-traverse-graphlet.md` §5–§6, the library's tool functions (`mcp_tools.py`), the
service's tool descriptions (`traverse_fetch`, `traverse_results` and the seven reading tools `graphlet_summary`,
`graphlet_labels`, `graphlet_walk`, `graphlet_walks`, `graphlet_sequence`, `graphlet_compare`, `graphlet_export`), and
the owner's answers of 2026-11-15.

| # | the user has | the user wants | where it arises | needs a new operation? |
|---|---|---|---|---|
| U1 | a retrieval to radius R, and the continuation of walk w (child) | one view of the long locus: one walk from the seed to R + E, support changes across the junction, one FASTA record, claims that do not stop at R | the intended loop "re-seed from a walk's continuation" (graphlets.rst); `next_request`, `deepen`, `traverse_continue` | **yes: stitching** (offered, §6–§8) |
| U2 | one seed fetched with `direction: left`, then with `direction: right` | one graphlet with both arms; or a check that the two agree with a `direction: both` fetch | a budget spent on one arm first; validating the engine and the library | **no** for the answer: a fetch with `direction: both` serves it. **Later: a validation helper** (§10) |
| U3 | one seed and one index retrieved twice (tuned and exhaustive, constrain and annotate, two radii, two `max_labels_per_node`) | which answer holds, and what one has that the other lacks | graphlets.rst "Comparing retrievals" | **no**: `compare()` answers it. The two are answers to different requests, so their union answers none |
| U4 | one seed constrained on label set A, then on set B | the walks of A ∪ B | label batches; a service request per labelled seed over several graphs | **no**: one request naming both sets answers it, with its own certificate (§3.2) |
| U5 | one seed on several pairs of one release (r237: 11 graphs, so 11 records), or in two databases | which graphs, and which records in them, carry this flank, and how far each graph was walked; one trie of the flanks; one FASTA or GFA | `traverse_fetch` over a database with named graphs ("a record per graph and seed, none merged"); the owner, today | **yes: the cross-pair union** (offered, §5, phase 1a-U): the best combination available today, with every limit stated |
| U6 | two seeds of one locus (a gene tiled with seeds) | a map of the locus | the service's `graphlet_compare` description ("two seeds of one locus") | **not a combination of graphlets**: positions and certificates belong to each seed (§3.4); left out |

U5 and U1 need new operations. U3 and U4 are served today, U2 is served and gets a later validation helper, and U6
does not combine graphlets.

### 2.1 Is each combination worth carrying?

The questions follow the checklist in `tools/writing-rules/guides/STEWARDSHIP.md`.

**The cross-pair union.**

- **Who needs it.** The owner, today: a seed fetched on r237 comes back as 11 records, and the question is the
  release's ("which graphs, which records, how far"). Also any library user who fetches one seed on several pairs
  or in two databases.
- **What happens without it.** The user reads 11 records and joins them by hand: spell each record's walks, line up
  equal prefixes, note which graph carries which branch. What the hand-made join loses is the honesty: which graphs
  were walked less far, which did not hold the seed, which were not asked, that a graph not following a walk may
  still hold it, that equal label names in two graphs are not shown to be one record, and that a dead end in one
  graph is not a dead end in the release.
- **The cheapest alternative.** A documented recipe. It is rejected because the five statements of §5.1 are
  exactly what a recipe user skips.
- **The smallest version.** Phase 1a-U (§18), after the review's cut list: `union_pairs()`, its refusals, the union
  trie, coverage per pair, labels per pair, walks, spelling, carriers of a walk, FASTA and GFA with the statements,
  summary and `print()`, evidence and provenance. No new file format: parts are saved as `.mgt` files today, and
  the union is pure, so building it again after `load()` is cheap.
- **Who keeps it.** The labeled-traversal library's maintainers (this project's lanes). The ongoing cost is the
  tests of §15.1 and keeping the composite in step with `ops.py`, which the one-part property test (P2) enforces.
- **How we know it is used, and how to remove it.** The owner uses it; beyond that, issues, questions and code
  search. Removal deletes one module and one docs section; nothing persisted depends on it until a bundle ships.

**Stitching a continuation.**

- **Who needs it.** Python users of the documented deepen loop (`TraverseClient.deepen()`, `next_request()`), and
  agents using the library's own tool functions (`traverse_continue`). No user request is on record. The need comes
  from the review, which found the gap, and from the docs, which prescribe the loop.
- **What happens without it.** A user reads parent and child separately and joins them by hand: spell the parent's
  walk, spell the child's walks, place the child at outward position L (§7.1) and drop the n-base overlap.
  Everything needed is public (`continuation().seed_coord`, the `overlap_bp` that `traverse_continue` returns). What
  the hand-made join loses is which labels the continuation did not carry, which knobs differ, which merges with
  other walks it could not see, and where the child walked edges the parent's walk had used. Without those, a
  stitched walk reads as one uninterrupted traversal.
- **The cheapest alternative.** A documented recipe with no API. The owner chose to build phase 1a instead (OQ1).
- **The smallest version.** Phase 1a-S (§18): `stitch()`, its refusals, the two-verdict junction statement, walks,
  spelling, FASTA, completeness, evidence and provenance. Claims and support across junctions (1b) and a bundle
  file (1c) come only if users need them.
- **Who keeps it, and how to remove it.** As for the union. The API is marked experimental in the docs until the
  oracle (T4) pins every whitelist entry, and is reviewed at the next library review.

## 3. Where the frame is wrong

### 3.1 "Merge" is taken

In this codebase, merge already means reconvergence: `on_reconverge: merge`, a segment with two or more parents
(`Segment.is_merge`), a run closed by a merge (`m`, `Run.merged`), `completeness_scope: united_history`. A
`Graphlet.merge()` that combines retrievals would sit next to "a merge before the depth" in every compare note. The
operations are therefore named `union_pairs` and `stitch`. "Merging graphlets is not supported" stays true in the
docs.

### 3.2 Same-pair union, "overlapping or not", is mostly ill-posed

Which walks a traversal holds depends on the whole strategy: the permitted set, the change cost and loss budget, the
branch allowance, the quorum knobs, the beam and caps, and (under `merge`) which edge histories were united and
which parent is displayed (DESIGN §11). Two retrievals of one arm on one pair are two answers to two requests. Their
union is the answer to neither, and no `complete_to_bp` describes it. SPEC §6.4 states the one case where a union is
exact: "Equality with the union of independent single-label walks holds only under `on_reconverge: keep`". In my
reading it also needs no switching, quorum or beam **[unsure]**. A request that names both label sets answers that
case with its own certificate. So:

- **one seed, one arm, two retrievals on one pair**: not offered. Keep the retrieval that
  `compare(mode='prefix_subset')` shows to contain the other, or fetch once with what both asked. `union_pairs()`
  refuses two retrievals of one pair (`duplicate_pair`, §5.2).
- **one seed, disjoint arms**: well defined, because arms share no used-edge sets and each arm certifies itself (SPEC
  §6.6). As an answer it is not needed: a fetch with `direction: both` serves it. As a check it is useful later
  (§10).
- **different seeds**: not offered (§3.4).

Across pairs the argument does not apply: each pair is walked by its own traversal, and the union keeps each pair's
answer and certificate apart (§3.3, §5).

### 3.3 Across pairs, a union by sequence states what it can

A `/traverse` request runs on exactly one pair (`SPEC-multigraph-envelope.md` §1.1); the service says "a traversal
stays inside its graph". Labels are index-local: `c:<column>` and `h:<column>:<seq_id>` are positions in one pair's
annotation (SPEC §7.5.7: joined across retrievals "only when the `index_fp` agree"), and two different pairs never
state one `index_fp` (SPEC §10.3). Different k makes two graphs' walks differ in what they can branch on, though
their bases still compare.

None of this is a reason to refuse. It is what a cross-pair union has to state, and the union states it (§5.1):

- **Walks are compared by their bases**, so different k, regime, release or database do not matter to the trie.
  Different alphabets do (§5.2).
- **Every certificate stays its pair's.** A pair is complete to a depth only as its own graphlet says, and the union
  is never more complete than its parts (§5.5).
- **A pair not following a walk did not follow it under its own request**: its mode, its permitted labels and its
  caps decide that, not only its graph (§5.3).
- **Labels are each pair's own**, shown with their pair. The server's certificates would license a bare name only
  where no name occurs in two pairs, so they never join anything (§5.4).
- **The union holds no walk that crosses graphs**: a walk that no single graph holds entirely is in no part (§5.6).
  On r237, where each graph is a taxonomic group, such a walk would join two organisms.

The second version of this note refused cross-pair composition and pointed to an engine feature ("the graph-less
traverse of a later release", envelope spec §1.1). That was wrong: the owner needs the combination now, and the
union is the only way to get it today.

### 3.4 Two seeds of one locus are not one graphlet

Every position is outward from one seed boundary, every certificate is about walks from that seed, and the seed's
permitted set decides the walks in constrain mode. Placing seed B inside A's walks is an alignment, and seeds can
overlap A's walks more than once (repeats) or not at all. What U6 asks for, a k-mer graph of the union with labels
per node, is a different product (a local assembly graph). `to_gfa()` per graphlet plus a graph tool does it today.
Left out (owner question OQ4).

### 3.5 The result must not be an MGT body

A merged MGT v1 body would need one `H` (one `seed_index`, one identity), one `S`, an `L` dictionary whose first
`num_seed_labels` ids are the seed labels (the root's `G.entry *` rule), `A` counters and counts, `B` growth bins and
`V` branch events of one traversal, one `O` and at most one `Q`. Positions would also pass the arm's own certificate.
Summing counters or bins would invent values; dropping them would break "lossless" and the `*` rules. Across pairs
it would also need one index identity, which a union of pairs does not have. MGT v1 is frozen (SPEC §7.5). A
combination is therefore a **container of unchanged bodies** plus a stated relation (one seed across pairs; a walk
and its continuation), which is how views already work (`J.view`, the backing body unchanged).

## 4. The operations

| operation | status | preconditions (all checked, in the order of the refusal tables) | refusal codes |
|---|---|---|---|
| cross-pair union `union_pairs(results)` | **offered**, phase 1a-U, **first** | results sorted into parts and absent records; every part a whole retrieval with an envelope and bases; one alphabet; one oriented seed, base for base; one answer per pair; expected pairs well formed (§5.2). Allowed and stated: different k, regime, release, database, strategies, radii, arm sets, unverified identity, caller-stated pairs | `not_combinable`, `view_unsupported`, `several_seeds`, `no_seed`, MissingEnvelope, `pair_mismatch`, `no_pair`, `no_bases`, `different_alphabet`, `different_seed`, `duplicate_pair`, `bad_expected`, `not_expected`; for queries `ambiguous_name`, `summary_too_small`, `not_offered` |
| continuation stitching `stitch(parent, child)` | **offered**, phase 1a-S, after the union (OQ1) | both carry an envelope; same index (§6.2: `index_meta_fp`, `index_fp`, `index_ns`, release, names per ref; one-sided `index_fp` only with `allow_unverified_index=True`); same pair where both state one; same k, regime and alphabet; same mode and support; the child walks one arm; its seed equals the continuation of one walk of that arm of the parent; that walk is not yet continued by another body; in constrain mode its seed labels are alive at the walk's end | `not_combinable`, `view_unsupported`, MissingEnvelope, `different_index`, `index_unverifiable`, `different_pair`, `different_graph`, `different_mode`, `not_one_arm`, `wrong_arm`, `no_continuation`, `not_a_continuation`, `ambiguous_walk`, `already_continued`, `foreign_seed_labels` |
| the two arms of one seed against a `direction: both` fetch | **later**, a validation helper (§10) | as stitch's identity checks; one oriented seed | none designed yet |
| union of two retrievals of one pair | **not offered** (§3.2) | n/a | `union_pairs()` refuses them `duplicate_pair`; `compare()` reads them |
| different seeds (overlapping or not) | **left out** (a different product, §3.4) | n/a | `different_seed`, `not_a_continuation` |

## 5. The cross-pair union

### 5.1 What it is, and the five statements

`union_pairs(results)` takes everything one seed's request returned across several pairs: the graphlets, and the
pairs that answered without one. It returns a `PairUnion`: per arm, one trie of every part's spellings (§5.3),
with the pairs and labels that carry each segment. The parts stay unchanged, and every answer is computed from them.

It is not exact, and it does not pretend to be. It is the best combination available today: everything one seed's
retrievals across a release say, in one place, with the limits next to the answers. Five statements make it
honest. Every answer carries them, structured (`statements()`) and as sentences (`summary()`, `print(union)`,
`evidence()`, the GFA comment lines), and the FASTA headers carry them as tokens (§5.8):

1. **Coverage per pair** (§5.5): at the depth asked, which pairs are complete, which are complete with a named
   qualification, which are certified only to a smaller depth and why, which do not hold the seed, which were not
   walked or failed, and which were expected but not given.
2. **Provenance per segment and walk** (§5.3): which pairs spell each segment, and how far each pair follows each
   union walk: to its end, ending earlier, merging into another of its walks, or turning off.
3. **What a turn-off means** (§5.3): a pair not following a walk did not follow it under its own request. Only a
   span marked `absence: graph` says the pair's graph lacks the walk.
4. **Labels per pair** (§5.4): every label is shown with its pair (`g1/NC_000913.3`), and equal names in two pairs
   are not claimed to be one record.
5. **The non-claim** (§5.6): "A walk that no single graph holds entirely is in no part, so not in this union."

The statements lead every summary and are never cut (§5.7, `summary()`).

### 5.2 Inputs, preconditions and refusals

**One list of results.** `union_pairs(results, *, seed=None, expected_pairs=None)` takes the results as they come
back, and the library sorts them into parts and absent records (review M5):

| element | becomes |
|---|---|
| a `Graphlet` with its envelope | a part; its pair from the envelope's tags |
| a `TraverseResponse` with one seed result | a part (its graphlet), or an absent record (its error result); the pair from the response's envelope |
| `{pair, graphlet}` | a part on the pair the caller states; refused `pair_mismatch` if the envelope states another. This unites CLI bodies and single-graph databases without a manifest |
| `{pair, result}`, optionally with `sequence` | an absent record from a seed result without a graphlet (SPEC §6.1 steps 2–4) |
| `{pair, error}` | a `failed` record: the pair's answer was refused or never came (`single_answer()` raises on a refused entry, `client.py:274`, and an HTTP failure raises), so there is no seed result |

A pair is `{graph, graph_path, annotation_path}` or a graph name. An absent record is classified by its
`outcome.walks` (§5.5). Its seed cannot be checked against the union's seed: a result's `seed` object carries
`seed_id`, `length_bp` and the labels, but no sequence (SPEC §7.1), and a seed id hashes the labels too (SPEC §4.3,
§6.1 step 3). So `length_bp` is checked where present, a given `sequence` is checked base for base, and every other
absent record is marked `seed_checked: false` and named in the coverage sentence.

**Required:** whole retrievals with their envelopes and their bases; one alphabet; one oriented seed sequence, base
for base (from the parts, or from `seed=` when no result carries a graphlet: "the seed is in none of the 11 graphs"
is a coverage answer, review minor 13); one retrieval per pair.

**Allowed, and stated per pair:** different k, regime, release, server or database; different strategies, including
mode (constrain, annotate), support (`kmer`, `trace`), radius, caps and `on_reconverge`; a part that walked one arm
only; a part whose identity is unverified (no `index_fp`) or whose pair the caller stated. Each pair's mode is
stated in the coverage sentence and on each span wherever a part is not in annotate mode (§5.3, review M1).

**Refusals.** The checks run in the table's order, element by element in the order given, and the first failing
check decides; a refusal names the first failing element in that order (commutativity, §5.9, covers successful
unions only). `%s` names an element by its tag (§5.7) or, before a tag exists, by its position. The last three codes
refuse queries, not construction.

| code | text |
|---|---|
| `not_combinable` (TypeError) | `union_pairs() takes results: a Graphlet, a TraverseResponse, or {pair, graphlet}, {pair, result} or {pair, error}; not %s` (once stitching ships, a StitchedGraphlet too) |
| `view_unsupported` | `union_pairs() reads whole retrievals, not views: pass view.backing, and select labels on the union (labels=)` |
| `several_seeds` | `result %s holds %d seed results: pass this union's seed as {pair, graphlet} or {pair, result}` |
| `no_seed` | `no result carries a graphlet, so the union has no seed sequence: pass seed= (the oriented seed the request sent)` |
| MissingEnvelope | today's `require_envelope('union_pairs()')` text, as for `stitch()` (§6.2) |
| `pair_mismatch` | `result %s was answered on pair %s, and the caller states pair %s: state the pair it was answered on, or none` |
| `no_pair` | `result %s names no pair (no graph tag, no index_ns): pass it as {pair, graphlet} or {pair, result}` |
| `no_bases` | `part %s carries no bases (output.sequences: false): a union by sequence needs every part's bases` |
| `different_alphabet` | `part %s walks the alphabet %s and part %s the alphabet %s: their bases are not comparable` (compared without `$`: succinct `$ACGT` and hash `ACGT` are one DNA alphabet) |
| `different_seed` | `part %s's seed is not part %s's (they first differ at base %d): a union by sequence needs one oriented seed, spelled base for base the same`; unequal lengths: `part %s's seed is %d bp and part %s's %d bp (equal for the first %d bases): a union by sequence needs one oriented seed`; reverse complement: `part %s's seed is the reverse complement of part %s's: a union by sequence needs one orientation; fetch part %s with the seed as part %s spells it`; an absent record: `the absent record of pair %s states a %d bp seed, and this union's seed is %d bp` (or, with `sequence`, the first variant) |
| `duplicate_pair` | `results %s and %s were answered on one pair (%s): a union holds one answer per pair; compare() reads two retrievals of one pair`; when they cannot be told apart: `results %s and %s cannot be told apart (no index_fp on %s, and their pair tags, index_ns, release and index_meta_fp agree): they may be one pair; state each pair with {pair, graphlet}` |
| `bad_expected` | `expected_pairs= takes the pairs the request was sent to, as graph names or {graph, graph_path, annotation_path} objects; %s is neither` |
| `not_expected` | `pair %s (a part or an absent record) is not among expected_pairs=` |
| `ambiguous_name` | `label name %r is in pairs %s: equal names in two pairs are not claimed to be one record; select by {pair, name} or {pair, ref}`; a selector `{'id': n}` is BadSelector, `label ids are local to each part of a pair union: select by {pair, ref} or {pair, name}` |
| `summary_too_small` (ValueError) | `the statements alone take %d bytes, more than max_bytes=%d: they are never cut; pass max_bytes >= %d` |
| `not_offered` | `dump()`: `a pair union is not an MGT document: each part dumps on its own (parts[i].graphlet.dump())`; `to_json()`: `to_json() is detail "full" of one server result, and a pair union has none: export its parts`; any other: `%s() is not offered on a pair union` |

- **One pair, or two.** The pair key is (release, `index_ns`, graph, graph_path, annotation_path, `index_fp`,
  `index_meta_fp`), the stated pair taking the place of the tags where the caller states one. Two results are one
  pair when both state an `index_fp` and the two agree (a pair listed under two names counts once, envelope spec
  §9), or when nothing tells them apart: no `index_fp` on one or both, and equal tags, `index_ns`, release and
  `index_meta_fp`, as `_identity` compares them (`ops.py:3588`).
- **An identical part is a no-op** only on the same pair with the same body digest, so a union can be rebuilt from
  overlapping sets of records. Two pairs with byte-identical bodies (possible without manifests, once the caller
  states the pairs) stay two parts.
- **A reverse-complement part is refused**, not flipped. On primary graphs a seed and its reverse complement are one
  seed; flipping a part's walks (each arm becomes the other, reverse complemented) is phase 2.
- **A `StitchedGraphlet`** (once stitching ships) is one pair's answer: its root and continuations on that pair. Its
  walks are its stitched walks, and a span's `certified_to_bp` is the stitched one, with the junction verdicts named.
  The relation goes this way only: a union of stitched parts, never a stitch of a union, since a continuation runs
  on one pair.

### 5.3 The union trie: spellings, segments, walks

**Spellings.** In each part, every segment's first parent makes a tree over the arm's segments. Its leaves are the
walks' leaves (a walk is the first-parent chain to its leaf, `derive.py:258`) and the ends of **stub** chains: a
segment that is only a non-first parent of a merge lies on no walk (`_first_paths`, `ops.py:641`), but its
first-parent chain ends at that merge. Every segment of a part lies on a walk or a stub, so the union holds every
base of every part (review M2). Positions are outward from the seed's boundary in the seed's orientation, so
spellings compare across k, regime and radius. In one graph a prefix from one seed names at most one path (a de
Bruijn path is fixed by its start and its bases), so a node of the union trie, a spelled prefix, is at most one
position on one path in each pair.

**Construction.** One depth-first pass over each part's first-parent tree inserts it into a radix trie, so the work
is linear in the part's bases. Inserting walk by walk would be quadratic on a comb, as `walk_rows()` warns
(`ops.py:1473`). Trie edges are compressed into union segments; a segment ends where the trie branches and where
the set of pairs following it changes. Children are ordered by their first base in the alphabet's order. A
segment's bases are a reference into the first part (in canonical order) that spells them, never copied.

Per segment and pair, a `SegmentSupport`: the part's walks or stubs passing it, the part segments it covers with
their offsets, and the labels the part records there. A part segment's displayed labels do not depend on which walk
passes it (a run counts on every walk that shares its displayed chain, `_displayed_profile`, `ops.py:492`), so they
are computed once per part segment by a new accessor (§14), with `labels_cut` where the part's lists are cut
(`label_lists`, `inexact_counts`). Annotate mode records one set per node.

**Walks and spans.** A union walk is a leaf of the trie. For each pair it has a `PairSpan`: how far that pair's
spellings follow it (`to_bp`) and how they stop:

| `how` | meaning |
|---|---|
| `whole` | the pair has a walk spelled exactly so; its end reasons, labels and `complete` are that walk's |
| `ends` | a walk of the pair ends at `to_bp`, before the union walk does, with its end reasons and labels. A walk of the pair that ends at the seed (`to_bp` 0, a root leaf) is `ends` too, with its reason |
| `merges` | the pair's stub ends at `to_bp`: there the pair's path merges into one of its displayed walks, and it continues as that walk does (`into`: the part's merge segment). The pair holds this walk's bases to `to_bp`; past it, its continuations are those of its walks through the merge segment |
| `diverges` | the pair's walks and stubs go on at `to_bp` with other bases (another union walk) |
| `none` | the pair did not walk this arm (`to_bp` 0) |

Every union walk is `whole` or `merges` in at least one pair, since every leaf ends some part's walk or stub.

**What a turn-off means** (review M1). In constrain mode a walk is followed only while some permitted label supports
every node of it (SPEC §6.10), and each pair derives its own permitted set (SPEC §6.1 step 4). Beyond a pair's
certificate nothing is certified: "a capped or pruned result is never evidence of absence" (SPEC §6.7). So every
`ends`, `diverges` and `none` span states `absence`:

| `absence` | when | it says |
|---|---|---|
| `graph` | `ends` with `dead_end`; or annotate mode, `per_path` scope, and `to_bp` < the part's `complete_to_bp` | the pair's graph holds no walk that continues this one under the traversal rule (no edge twice, no seed entry: SPEC §6.6) |
| `under_strategy` | inside the certificate otherwise: constrain mode, a `united_history` scope, or an end reason other than `dead_end` | the pair's request did not follow it, for example because no permitted label of that pair supports it; the graph may hold it |
| `uncertified` | `to_bp` ≥ the part's `complete_to_bp` | nothing is known beyond `to_bp` in this pair |

A `none` span is `under_strategy` (the request named the other arm). Each span also carries the part's mode, and
`complete` and `certified_to_bp`: for `whole` and `ends` spans the part walk's own (`complete`, and `length_bp −
beyond_certified_bp`, as `_walk_end` computes them, `ops.py:1463`), for `merges` and `diverges` spans the smaller of
`to_bp` and the part's `complete_to_bp`.

Union walk ids are the leaves in trie order (children by base), stable for one set of parts. Adding a part can
renumber them; the exports also name each walk by its spelling's length and pairs.

**The union of one part is that part.** The union of one graphlet has one union walk with a `whole` span per part
walk, with its spelling, length, end reasons, labels and completeness; every other union walk ends a stub
(`merges`); and the coverage is the part's arm certificate. This is the property test P2 (§15.1).

### 5.4 Labels across pairs

Every label is a `PairLabel(pair, ref, name, kind)`, shown as `<tag>/<name>` (`g1/NC_000913.3`, `g4/562`). Refs are
index-local, so two pairs' refs never join, and names are never joined either (review M3). The server's certificates
cannot license a join: `headers_unique` true means no header repeats "in a pair or across pairs" (envelope spec
§9), and `columns_disjoint` true means no column name is a column of two pairs (SPEC §10.3). Wherever they hold, no
name occurs in two pairs, so there is nothing to join; where a name does occur in two pairs (a record in both
graphs), they do not hold. A join rule would only decide whether the tag is shown, so the tag is always shown and
the union reads no certificate.

- `label({'name': n})` resolves when exactly one `PairLabel` of the union has that name; otherwise it is refused
  `ambiguous_name`, listing the pairs. `{pair, name}` and `{pair, ref}` always resolve within the pair.
- `shared_names()` lists the names seen in two or more parts, with their pairs. The labels statement gives their
  count: "equal names in two pairs are not claimed to be one record (3 names occur in more than one pair)".

A header name rests on its pair's record mapping. On a pair whose `seqs_order` is `"count_only"`, as on every r237
pair today, that mapping's column order was not checked at load, so a header name can name another record even
inside one pair (envelope spec §9.1). The tag does not repair that. It is not a union matter: the library states it
for no graphlet today, since it reads no `GET /capabilities` certificate (§17.1).

### 5.5 Coverage per pair: never more complete than the parts

`coverage(arm, at_bp=None)` puts every pair in one category. Without `at_bp`, `summary()`, `print()` and the exports
use the largest radius (`bounds.max_extension_bp`) among the parts, and `coverage()` gives each pair's own certified
depth (review minor 3).

| category | the statement | from |
|---|---|---|
| `complete` | the part holds every walk from this seed in this pair's graph, under the part's strategy, per path, to `at_bp` or to its end before it | the arm's `complete_to_bp` ≥ `at_bp`, and no walk-class limitation on the arm or the seed |
| `complete_qualified` | as `complete`, with the qualification named: `seed_labels` (carriers of the seed were cut: "complete for 1000 of 4312 seed labels"; `max_seed_labels` is 1000 by default, SPEC §5) or `scope` (a merge united histories: complete for the united-history rule, not per path) | the arm's `complete_to_bp` ≥ `at_bp`, and the walk-class limitations that make the part's `outcome.walks` `partial` (SPEC §7.0) |
| `incomplete` | certified only to `x` < `at_bp`, with the reason | the arm's `complete_to_bp`, and the reason from its `status` and `cap_trigger.reason` (truncated, pruned), else the part's radius (review minor 2) |
| `seed_not_in_graph` | the graph does not hold the seed as one sequence, so no walk from it exists there | an absent record with `outcome.walks: "not_in_graph"`, and `not_in_graph.kmers_present` of `kmers` |
| `labels_not_in_graph` | none of the seed labels is in this pair, so nothing was walked | an absent record with `outcome.walks: "labels_not_in_graph"` |
| `not_walked` | a seed result without arms, with its cause; "the seed is in this graph" where the cause is a derivation, whose step follows the presence check (SPEC §6.1 steps 2 and 4: "nothing carries the seed in full", review minor 1) | an absent record with `outcome.walks: "failed"` |
| `arm_not_walked` | the part's request did not walk this arm | the part's `direction` |
| `failed` | the pair's answer is an error, so nothing is known | a `{pair, error}` record, its text cut at 60 bytes |
| `not_included` | the request went to this pair, and the union has no result for it | `expected_pairs=` |

`every_pair_certified_to_bp` is the depth to which every pair is `complete`, `seed_not_in_graph` or
`labels_not_in_graph`: the smallest such part's `complete_to_bp`. It is `None` when any pair is in another category,
`complete_qualified` included. Without `expected_pairs=`, every statement is about the pairs given and says so
("over the pairs given; expected_pairs= was not passed"). A pair is never stated complete beyond its part's
certificate, and no union-level statement is stronger than the weakest pair it covers.

The sentence (one per arm and depth) lists the non-empty categories in the table's order, pairs in canonical order,
each pair's mode where any part is not in annotate mode, the unverified absent records, and ends with the
non-claim:

    right arm, to 300 bp (constrain, kmer, in every pair): complete in g1, g4; complete, qualified, in g6 (complete
    for 1000 of 4312 seed labels); incomplete in g7 (certified to 120 bp: max_steps); seed not in graph: g2 (0 of 10
    k-mers), g3 (4 of 10 k-mers); not walked: g5 (seed in graph; derivation: no label carries the whole seed);
    not included: g9. Absent records not checked against the seed: g2, g3, g5. A pair that does not follow a walk did
    not follow it under its own request. A walk that no single graph holds entirely is in no part, so not in this union.

### 5.6 What the union does not hold

Stated once, plainly, in every summary, evidence block and export: **a walk that no single graph holds entirely is
in no part, so not in this union.** Every walk the union holds lies inside one graph, and every end reason is that
graph's: `dead_end` in g7 means g7 has no next k-mer, not that the release has none. A seed only partly present in a
graph (`seed_not_in_graph` with some k-mers present) contributes nothing from that graph.

What it costs: on r237 each graph is a taxonomic group, so a real sequence (a chromosome, a plasmid) lies within one
graph, and a walk across graphs would join two organisms. A sequence shared by organisms of two groups (a plasmid, a
transposon, a conserved gene) is held entirely by both graphs, so it is in the union with both pairs and each pair's
labels. Where a release splits its graphs by something other than organism, the non-claim can cost real walks, and
the union cannot measure that cost.

What the union cannot give: walks that cross graphs. Those need one traversal over several graphs, which no engine
offers today.

### 5.7 Python API

New module `union.py`, exported from `metagraph.traverse`; the parts container, refusals, evidence and budgets are
in `composite.py`, shared with `stitch.py` (§14). Everything is phase 1a-U except the line marked `[1c]`.

```python
union_pairs(results, *, seed=None, expected_pairs=None, budget=None) -> PairUnion
    # results: Graphlet | TraverseResponse (one seed) | {pair, graphlet} | {pair, result[, sequence]} | {pair, error}
    # seed: the oriented seed sequence; required only when no result carries a graphlet
    # expected_pairs: the pairs the request was sent to (graph names or pair objects); None: not known
    # A new object every time; parts shared, never copied or modified.

class PairUnion:
    parts -> tuple[PairPart]              # canonical order (§5.9)
    absent -> tuple[PairAbsent]
    pairs -> tuple[str]                   # every tag: parts, absent records, missing expected pairs
    expected -> tuple[str] | None
    seed -> SeedInfo                      # the common oriented seed; each part's own seed record in parts
    alphabet; arms() -> list[str]         # the arms any part walked
    statements(arm=None, at_bp=None) -> Statements        # the five statements (§5.1), structured
    coverage(arm, at_bp=None) -> Coverage                 # §5.5
    labels(pair=None) -> list[PairLabel]
    label(selector) -> PairLabel          # {'pair', 'ref'} | {'pair', 'name'} | {'name'} when only one pair has it
    shared_names() -> dict                # name -> tags, for names seen in two or more parts
    segments(arm, *, budget=None) -> list[UnionSegment]
    walks(arm, *, top=None, by='pairs'|'length'|'id', pairs=None, labels=None, min_bp=0, spell=False,
          budget=None, resume=None) -> list[UnionWalk]
    spell(arm, walk, orientation='natural'|'walk', with_seed=False, *, budget=None) -> str
    carriers(arm, walk, *, at_bp=None, budget=None) -> list[PairSpan]   # which pairs and labels carry it
    to_fasta(arm=None, walks=None, with_seed=True, orientation='natural', width=None, *, at_bp=None,
             budget=None) -> str
    to_gfa(*, with_seed=True, labels='counts'|'full', at_bp=None, budget=None) -> str
    evidence(arm=None) -> dict
    summary(max_bytes=4096, *, budget=None) -> dict       # the statements first and whole (summary_too_small)
    __str__() -> str                      # the statements as sentences: print(union)
    provenance() -> list[dict]
    save(path, *, budget=None)            # [1c] the bundle (§9)
    digest -> str
    memory_bytes() -> int
    # Refused (CombineRefused 'not_offered', §5.2): dump(), to_json(). Not in 1a-U (§18): compare() of two
    # unions, support_by_depth(), carriers(sequence=), routes, a per-pair FASTA, next_request() per pair.

@dataclass(frozen=True)
class PairPart:    tag: str; pair: tuple | None; pair_stated_by: str     # 'envelope' | 'caller'
                   graphlet: Graphlet; identity: dict; digest: str
@dataclass(frozen=True)
class PairAbsent:  tag: str; pair: tuple; category: str; detail: dict; seed_checked: bool
                   # category: 'seed_not_in_graph' | 'labels_not_in_graph' | 'not_walked' | 'failed'
@dataclass(frozen=True)
class PairLabel:   pair: str; ref: str; name: str; kind: str        # str(): '<tag>/<name>'
@dataclass(slots=True)
class SegmentSupport:  pair: str; walks: tuple[int, ...]; kind: str  # 'walk' | 'stub'
                       pieces: tuple[tuple[int, int, int], ...]     # (part segment id, from_bp, to_bp)
                       labels: tuple[PairLabel, ...]; labels_cut: bool
@dataclass(slots=True)
class UnionSegment:    id: int; arm: str; from_bp: int; to_bp: int; parent: int | None; children: tuple[int, ...]
                       support: tuple[SegmentSupport, ...]          # one per pair spelling it, canonical order
@dataclass(slots=True)
class PairSpan:    pair: str; mode: str; how: str; to_bp: int; walk: int | None; into: int | None
                   end_reasons: dict; labels: tuple[PairLabel, ...]; labels_cut: bool
                   complete: bool; certified_to_bp: int; absence: str | None   # None for 'whole' and 'merges'
@dataclass(slots=True)
class UnionWalk:   id: int; arm: str; length_bp: int; sequence: str | None; segments: tuple[int, ...]
                   spans: tuple[PairSpan, ...]                      # one per pair, canonical order
                   whole_in: tuple[str, ...]                        # the pairs with a 'whole' span
@dataclass(frozen=True)
class Coverage:    arm: str; at_bp: int | None; per_pair: dict; by_category: dict
                   every_pair_certified_to_bp: int | None; sentence: str
```

**Tags.** A part's tag is the pair the caller stated, else its graph name on a multi-graph server, else its
`index_ns` on a single-graph one; qualified as `<release>/<graph>` when the parts span releases, and by an
`index_fp` prefix if two tags still collide. A result with none of these is refused `no_pair`. Tags are display
names; pairs are compared by their keys (§5.2).

**The queries.**

- `walks(arm, by=...)`: `'pairs'` ranks by the number of pairs with a `whole` span, then length, then id;
  `'length'` by length; `'id'` keeps trie order. `pairs=` keeps the walks those pairs follow whole; `labels=` the
  walks those labels carry at a span end.
- `carriers(arm, walk, at_bp=None)` answers "which pairs and labels carry this extension": per pair, how far it
  follows, how it stops, with `absence` and mode, and the labels at `at_bp` (default: at each span's end).
- `spell()` spells a union walk from the parts' bases.
- No cross-pair label total is offered (review M4). Label counts are per pair, as label entries of that pair, and
  marked as lower bounds where the pair's lists are cut.
- `summary(max_bytes=4096)` puts the statements first and whole; if they alone exceed `max_bytes`, it refuses
  `summary_too_small` rather than cut them (review minor 9). For r237's 11 pairs and two arms they take about
  1–2 KB [estimate].

### 5.8 Exports

**FASTA.** One record per union walk (the parts' own walks are exported by each part's `to_fasta()`). Each header
states the walk's provenance, its coverage and the non-claim as `key=value` tokens, empty ones left out:

    >right.u12 len=812 whole=g1,g4 merges=g6@85 ends=g7@120:dead_end/graph diverges=g8@40/under_strategy
     certified=g1:812,g4:812,g6:85,g7:120,g8:40 labels=g1:3,g4:1+ modes=constrain
     absent=g2:seed_not_in_graph,g3:seed_not_in_graph not_included=g9 scope=within_one_graph

(one line in the file). A span token is `<tag>@<to_bp>[:<end reason>]/<absence>`; `labels=` counts label entries
per pair, `+` marking a lower bound (cut lists; the full lists with `labels='full'`, as `<tag>/<name>`);
`modes=` appears whenever a part is not in annotate mode (`modes=g1:constrain,g7:annotate` when they differ);
`scope=within_one_graph` is the non-claim. A part without bases is refused at construction (§5.2).

**GFA.** GFA 1.0 on the seed strand, like `to_gfa()`, with the statements as comment lines first:

    H	VN:Z:1.0
    # union_pairs: seed 40 bp; parts g1 g4 g6 g7 g8; seed not in graph g2 g3; not included g9
    # right arm, to 300 bp (constrain, kmer, in every pair): complete in g1, g4; ...
    # turn-offs: a pair not following a walk did not follow it under its own request (absence on each span)
    # labels: per pair, as <tag>/<name>; equal names in two pairs are not claimed to be one record (3 shared)
    # scope: a walk that no single graph holds entirely is in no part, so not in this union
    S	seed	CGTACC...	LN:i:40
    S	r12	ACGTT...	LN:i:300	pp:Z:g1,g4	pl:Z:g1=3;g4=1+
    L	seed	+	r12	+	0M
    P	right.u12	seed+,r3+,r12+	*	pp:Z:g1,g4	pe:Z:g7=120:dead_end/graph;g8=40/under_strategy

Segments are union segments, and links overlap `0M`: the union is a trie of spellings, not a de Bruijn graph, and its
parts' k can differ. The lower-case tags are user tags: `pp` the pairs spelling a segment or following a path
whole, `pl` the label entries per pair (counts, or `<tag>/<name>` lists with `labels='full'`), `pe` the other spans
of a path, as in FASTA. Left-arm segments are written in natural orientation, as `to_gfa()` writes them.

**Comparing.** A union-level `compare()` is not in 1a-U (§18). Each common pair's two retrievals compare with the
parts' own `compare()` (`ops.py:3778`).

### 5.9 Determinism and the algebra

- **Pure.** `union_pairs()` is a function of the results' texts and envelopes, the stated pairs, `seed` and
  `expected_pairs`.
- **Canonical order.** Parts and absent records are ordered by the pair key of §5.2, byte order, `None` first. Every
  list answer has a canonical order: arms left then right, walks by id under `by='id'`, ties of other rankings by
  id, spans and supports in pair order.
- **Commutative** for successful unions: any order of the results gives the same digest and the same answers, byte
  for byte. A refusal names the first failing element in the order given.
- **Identity element.** The union of one part answers as that part (§5.3, P2).
- **Idempotent for identical bodies on one pair.** Given twice, a no-op; another body on that pair is refused
  (§5.2).
- **Digest.** The sha256 of the canonical list of (pair key, tag, `pair_stated_by`, part digest), the absent
  records as canonical JSON, `seed` and `expected_pairs`.
- **Evidence.** `evidence(arm)` lists each part's own `evidence_block` (`ops.py:6006`) under its tag, the absent
  records, the expected pairs not included, each part's identity and pair source, and the five statements. An
  answer that reads several parts states, as `outcome`, the weakest value per dimension over the parts it reads,
  and each part's outcome under its tag. `exact` holds only when every part read is exact and verified.
- **Provenance and replay.** `provenance()` lists each part (`{tag, pair, pair_stated_by, identity, seed_id,
  digest, request}`, the request built as `GraphletStore.request_of()` does, `store.py:1277`), each absent record
  and the expected pairs. Replaying a union replays each part (the identity rules apply as for
  `traverse_fetch(replay=)`) and builds the union again.

## 6. Stitching a continuation: the contract

### 6.1 Python API

New module `stitch.py`, exported from `metagraph.traverse`, on the composite machinery of `composite.py` that the
union builds first (§14). Everything is phase 1a-S except the lines marked `[1b]` and `[1c]` (§18).

```python
stitch(parent, child=None, *, arm=None, walk=None, allow_unverified_index=False, budget=None) -> StitchedGraphlet
    # parent: Graphlet (with envelope) or StitchedGraphlet; child: Graphlet (with envelope) or StitchedGraphlet
    # rooted at the continuation; None -> the parent alone as a StitchedGraphlet with no junction.
    # arm: 'left' | 'right' | 'l' | 'r' (default: the child's one arm). walk: the continued walk, a path id on a
    # Graphlet parent, a walk key (tuple) on a StitchedGraphlet parent; default: inferred (§6.2). One junction
    # per call; a new object every time (parts shared, never copied or modified).

class StitchedGraphlet:
    parts -> tuple[Part]                  # root first, then by part key (§8)
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
    continuation(arm, walk, *, budget=None) -> Continuation       # the owning part's, seed_coord mapped (§7.1)
    next_request(arm, walks, bp=None, reduce_budget=True, reset_branches=False, **overrides) -> NextRequest
    next_requests(...)                                            # delegate to the part owning the walk(s)
    to_fasta(arm=None, walks=None, with_seed=True, orientation='natural', width=None, *, coordinates=None,
             budget=None) -> str        # coordinates as Graphlet.to_fasta (model.py:571), per piece (§7.7)
    completeness() -> dict                # per arm: the root's complete_to_bp and one extension per junction
    evidence(arm=None) -> dict            # §7.6
    summary(max_bytes=2048, *, budget=None) -> dict
    provenance() -> list[dict]            # §7.8
    save(path, *, budget=None)            # [1c] the bundle, §9; metagraph.traverse.load(path) returns it
    digest -> str                         # sha256 of the canonical tree: the part keys and the parts' digests (§8)
    memory_bytes() -> int                 # sum of the parts' models + the junction indexes
    # Deliberately refused (CombineRefused 'not_offered', text in §6.2): dump(), to_json(), compare(),
    # subgraph(), to_gfa() (phase 2), label_walks(), routes(), label_summary(), splits() (phase 2)

@dataclass(frozen=True)
class Part:      key: tuple; graphlet: Graphlet; offset_bp: int; digest: str
                 # key () = the root; (('right', 3),) = the child of root walk 3; (('right', 3), ('right', 0)) ...
@dataclass(frozen=True)
class Junction:  arm: str; walk: tuple; at_bp: int; overlap_bp: int; part: tuple;
                 verdict: str                                   # 'exact' | 'qualified'
                 causes: tuple[tuple[str, str], ...]            # (code, sentence), in the order of §7.4's table
                 continued: tuple[Label, ...]; not_continued: tuple[Label, ...]
@dataclass(slots=True)
class StitchedWalk:  key: tuple; arm: str; length_bp: int; sequence: Optional[str]; path_reason: Optional[str];
                     end_reasons: dict; labels_full: list[Label]; n_alive: int; complete: bool;
                     verdict: str                               # 'exact' | 'qualified': its junctions and own causes
                     certified_to_bp: int; beyond_certified_bp: int; junctions: tuple[int, ...];
                     pieces: tuple[tuple[tuple, int, int, int], ...]   # (part key, path_id, from_bp, to_bp)
                     causes: tuple                                      # this walk's own causes (§7.4)
@dataclass(slots=True)
class StitchedClaim(Claim):  walk: tuple; pieces: tuple; junctions: tuple[int, ...]       # [1b]

class CombineRefused(ValueError):  code: str     # every refusal of §5.2 and §6.2 (composite.py); MissingEnvelope
                                                 # and BadSelector as today
```

A walk key is the tuple of walk ids along the chain: `(3,)` is root walk 3, and `(3, 0)` is walk 0 of the child that
continues root walk 3. Positions in every answer are outward from the root seed's boundary. Each row also carries its
pieces in their own parts' coordinates. A label selector `{'id': n}` is refused on a stitched graphlet (BadSelector,
"label ids are local to each part of a stitched graphlet: select by ref or name"). Refs and names resolve over the
union.

### 6.2 Refusals: exact texts

`%s` and `%d` are filled as stated; a label list prints at most 8 `{name, ref}` then "and N more", as `next_request`'s
notes do.

| code | text |
|---|---|
| `not_combinable` (TypeError) | `stitch() combines graphlets (Graphlet or StitchedGraphlet), not %s` |
| `view_unsupported` | `stitch() reads whole retrievals, not views: stitch view.backing instead` |
| MissingEnvelope | today's `require_envelope('stitch()')` text: `stitch() needs the response envelope (seed_id, the normalized strategy, annotation counters): build the graphlet with Graphlet.from_response() or load a saved .mgt file with its J line; a parsed body alone does not carry it` |
| `different_index` | `the child was retrieved from another index than its parent (%s): a continuation runs on its parent's index, so a child of another index continues none of its walks; union_pairs() combines one seed's graphlets across pairs`, `%s` one of `index_meta_fp differs`, `index_fp differs`, `index namespace %r vs %r`, `release %r vs %r`, `label %s is %r in one and %r in the other` |
| `index_unverifiable` | `the child states %s index manifest digest (index_fp) and its parent %s: whether they are the same index cannot be verified (not proven different), so the child is not stitched automatically; pass allow_unverified_index=True to stitch it anyway (the result then states the identity unverified)` (`an`/`no`, `one`/`none`, as `_same_index`, `mcp_tools.py:2578`) |
| `different_pair` | `the child was answered on pair %s, its parent on pair %s: a continuation runs on its parent's pair` (each `{graph, graph_path, annotation_path}` as compact JSON) |
| `different_graph` | `the child's graph is not its parent's (%s %r vs %r): its walks cannot continue the parent's` (`k`, `regime` or `alphabet`) |
| `different_mode` | `the child was walked in %s mode with %s support, its parent in %s mode with %s support: its walks answer another question` |
| `not_one_arm` | `the child walked both arms: a continuation walks only the arm it continues (next_request() sets direction to it)` |
| `wrong_arm` | `the child walked the %s arm, not the %s arm named by arm=` |
| `no_continuation` | today's `_continuation` text: `walk %d of the %s arm has no continuation (it ended for a semantic reason, not at the radius or a cap)` |
| `not_a_continuation` | with `walk=`: `the child's seed is not the continuation of walk %s of the %s arm: a stitched child is seeded with that walk's continuation (continuation(), next_request())`; without: `the child's seed is the continuation of no walk of the %s arm: a stitched child is seeded with a walk's continuation (continuation(), next_request())` |
| `ambiguous_walk` | `the child's seed is the continuation of walks %s of the %s arm: name the one it continues with walk=` |
| `already_continued` | `walk %s of the %s arm is already continued by another child (part %s, digest %s): one child per walk` |
| `foreign_seed_labels` | `the child seeds %s, not alive at the end of walk %s of the %s arm: their runs would start new lineages at %d bp, not continue the walk's` |
| `not_offered` | `dump()`: `a stitched graphlet is not an MGT document: save() writes its parts unchanged in a bundle, and each part dumps on its own (parts[i].graphlet.dump())`; `to_json()`: `to_json() is detail "full" of one server result, and a stitched graphlet has none: export its parts`; `compare()`: `compare() reads one retrieval per side: compare the parts`; any other: `%s() is not offered on a stitched graphlet` |

Checks run in the table's order, and the first failing check decides. The identity checks reuse the index rules of
`_identity` (`ops.py:3588`) **without** its seed rule, because parent and child seeds differ by design, and with the
one-sided-manifest rule of `_same_index`: refused unless opted in. When neither side has `index_fp`, the stitch is
not refused and the identity is stated as unverified, as `traverse_continue` does (OQ6).

- **Reachability.** Fingerprints are per pair (SPEC §10.3), so with manifests a different pair or graph is always
  refused `different_index` first. `different_pair` and `different_graph` are reached only when no fingerprint
  decides (no manifest on either side, and `index_meta_fp` equal or absent). The pair is not checked first on
  purpose: a fingerprint mismatch is proof, while pair tags can name one index under two names (a pair listed under
  two names counts once, envelope spec §9).
- **Walk inference (one input, one code).** Without `walk=`, the candidates are **all** walks of the arm whose
  continuation sequence equals the child's seed, continued or not. Several candidates: `ambiguous_walk`. One: that
  walk. The answer therefore does not depend on which children were stitched before.
- **An identical child is a no-op.** A child whose body digest equals the one already stitched at that walk returns
  the stitched graphlet unchanged, so rebuilding from a set of parts tolerates duplicates. A different body at a
  continued walk is refused `already_continued`.

## 7. Stitching: semantics

### 7.1 Coordinates, orientation and the overlap

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
nodes of the parent's walk from L − n on (with root-seed nodes when n > L). §7.4 uses this.

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
spelled `ATTGCTAAGA + child.spell('right', v)`. `fork` is also a case for §7.4's `sibling_merge_unknown`: it runs
`merge`, and its right arm has a second walk with a continuation.

### 7.2 Walk identity and overlap

A stitched graphlet's walks on arm a are the root's walks, except that every continued walk w is replaced by `(w, v)`
for each walk v of its child (recursively for chains). Nothing is renumbered: `(w, v)` names root walk w and child
walk v. Continued walks of one parent with identical tails are not deduplicated, because each child is its own trie;
under `merge` such children are a `sibling_merge_unknown` case (§7.4). The overlap is exactly the child's seed
(`overlap_bp = n`, the value `traverse_continue` stores), so nothing is double counted. The child's positions start
at `L`, where the parent's walk ends.

### 7.3 Labels

Labels join by `ref` across parts. This is valid because the identity check proved one index, or the caller accepted
an unverified one, which is then stated. Names must agree per ref (`different_index` otherwise, as `compare()`'s
names rule). In constrain mode the child's seed labels (its first `num_seed_labels` dictionary entries) must be alive
at the end of walk w (the leaf's end set, `G.end`), otherwise `foreign_seed_labels`. Labels alive at w that the
child does not seed are listed in the junction's `not_continued` (`alive_not_seeded` as `next_request` reports it,
or simply not named). In annotate mode the child has no seed labels (`ops.py:2769`), and none of this applies.

### 7.4 The junction statement: two verdicts

**The reference.** One uninterrupted walk of the parent's request with a larger radius and larger caps. A
continuation is a new traversal (DESIGN §5.3): branch state, loss budget, edge history, seed re-entry, merges with
other walks and (under `trace`) the followed occurrence all restart. Each junction derives, from the parent (its DAG,
the leaf's `T` and `C` records, the runs ending there, its normalized strategy) and the child (its normalized
strategy, `X` records and body), whether the child's walks below w are the reference's walks below w.

**Why only two verdicts.** The first version also had `lower_bound` ("the reference holds everything the stitched
walk shows"). The review's M1 showed it is not monotone: fewer labels or less budget can also *lengthen* a walk.
Fewer paths mean fewer merges and so fewer united edge histories (SPEC §6.4, §6.5). Fewer non-empty successors mean
fewer splits, so the reference can end with `split_limit` where the child goes on (SPEC §6.4; default 64). Under
switching with a finite branch allowance, a source that is ambiguous in the reference can be unambiguous in the
child (SPEC §6.4; the re-minimisation is greedy, SPEC §6.3). A direction can be stated only under conditions the
oracle has pinned, so phase 1a-S ships `exact` on a narrow whitelist and `qualified` with named causes for everything
else. `lower_bound` can come later, with M1's conditions (`keep`, an unlimited split limit, no quorum, minimum or
beam, `forbid` or unlimited branches) pinned by T4.

**The whitelist.** A junction is `exact` when all of W1–W7 hold; a stitched walk is `exact` when all its junctions
are and W8 holds for it.

- **W1.** `support: kmer` on both parts.
- **W2.** The two normalized strategies are equal except in the allowed fields: `direction`; the radius and the
  resource caps (`bounds.max_extension_bp`, `max_steps`, `max_paths`, `max_output_bp`, `time_budget_ms`,
  `max_memory_mb`, `max_work_units`, and `max_live_paths` where `frontier.on_overflow` is `stop` on both); `output.*`;
  `annotation.*`; `labels.max_labels_per_node` (it caps output lists only, SPEC §5); `labels.loss_budget`,
  `labels.extra` and `branching.max_label_branches` (judged by W3–W4). The child's `labels.extra` names no label
  outside the parent's permitted pool (its seed labels and `labels.extra`).
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
- **W8 (per stitched walk).** Both parts carry bases. The walk takes no step the reference would block and the child
  allowed: no (k+1)-mer of the DAG ancestors of w's leaf (canonicalized in primary and canonical regimes) and no
  entry into a root-seed node (SPEC §6.6). It took no step the child blocked and the reference would allow: no child
  `rejoined_seed` on a node of the child's seed that is not a root-seed node (§7.1: a node of the parent's walk
  from L − n; on a root-seed node both stop with `rejoined_seed`). It did not exceed
  `max_splits_per_path` counting the splits on w and on the child's walk together (the child counts from 0; finite by
  default, 64).

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
| `no_bases` | W8 cannot be checked: a part has `output.sequences: false` | junction | both |
| `rule_violation` | W8: the walk reuses an ancestor's edge or enters a root-seed node; stated with `at_bp`. The walk is not cut and no end is invented | walk | both |
| `seed_reentry` | W8: a child `rejoined_seed` on a node of the child's seed outside the root seed. Two sub-cases, told apart by the blocked event's character (`E … b`): `same_stop`, where the blocked (k+1)-mer is an ancestor's edge, so the reference stops there too (with `edge_reuse`, by the precedence of SPEC §6.6) and only the reason differs; `earlier_stop`, a new edge into a node of the walk, where the reference goes on. Common: with the whole walk as the child's seed, any repeat of the walk ends the child here (§7.1) | walk | both |
| `split_limit_restarted` | W8: the splits on w plus the child's exceed `max_splits_per_path`, stated with the position of the first split past it | walk | both |

`same_stop` is the first candidate for promotion: its walks are the reference's and only the end reason differs, so
once T4 pins it, it can be `exact` with the reference's reason restated. Until a whitelist entry or a promotion is
pinned by its T4 fixture, the verdict is `qualified`. The verdict never depends on notes `next_request()` returned: a
child made another way (a service job with preset and overrides, a request written by hand) is judged by what it ran.
Hairpins (`hairpins: skip | follow`) are decided on one (k+1)-mer and behave the same given the same knob, which W2
requires. `tip_window_bp` and `bubble_window_bp` are rejected unless 0 (SPEC §6.5), so W2 keeps them equal; if they
are ever implemented, a non-zero window becomes a `qualified` cause.

### 7.5 Claims and support across a junction (phase 1b)

- **Constrain.** For each continued label X, the parent's claim of X that is open at the end of walk w (end class
  `radius`, `capped` or `pruned`: a leaf ended by `beam_pruned` also carries a continuation, SPEC §7.1) is joined
  with each claim of the child's lineages descending from X's seed run. A switch-entered child run joins through its
  `prev_run` chain to the seed run it descends from. The joined claim runs from the parent claim's `from_bp` to `L +
  child.to_bp`. Its `evidence_from` is the parent's, or `L + child.evidence_from` when the child's displayed support
  starts later (a merge in the child). Its end fields come from the child. `loss` is the parent's terminal loss of
  the descent's seed label plus the child's loss, and `branches` adds the same way. Under an `exact` junction (W3:
  `forbid`) every loss is 0. Under a `qualified` one with switching, the loss is that of the route the child chose,
  which can exceed the reference's minimum (an upper bound when the reference has that route). `entered_by` is the
  parent's. `pieces` lists both runs. The parent's claims of labels in `not_continued` stay open at `L`, with the
  junction named; they are never turned into ends. Claims of other walks are unchanged. Order: the root's claims in
  the root's order, each joined claim in place of the open one it extends, then each child's remaining claims in
  part-key order.
- **Annotate.** No lineages. The children's claims are shifted by `L` without joining (`joined: false` in the
  evidence). That understates continuity and is honest. Joining by presence at the junction is phase 2.
- **Support.** A stitched walk's profile is the parent's on w over `[0, L)` followed by the child's on v, shifted. At
  `L`, `support_changes()` reports the labels in `not_continued` as removed with reason `junction` (the continuation
  did not seed them, so the removal is not biology), and those the child's seed dropped with reason `junction` and
  the `X` reason. The junction never adds a label: none is reported as added at `L` because of it.

### 7.6 Completeness and evidence: nothing is laundered

- **The arm.** `complete_to_bp` of a stitched arm is the root's, unchanged: stitching adds walks below continued
  walks only and never makes "every walk to depth D" true for a larger D. `completeness()` adds one extension per
  junction: `{walk, at_bp, extended_to_bp: offset + child.complete_to_bp, certified_to_bp, verdict, causes}`.
- **A stitched walk.** Pieces are checked in order. Piece i certifies its part's walk up to its part's
  `complete_to_bp`. Certification passes to piece i + 1 only if piece i's walk ends within its part's
  `complete_to_bp`, the junction between them is `exact`, and the walk has no walk-scope cause at or before that
  junction. `certified_to_bp` is where it stops, `beyond_certified_bp = length_bp − certified_to_bp`, and `complete`
  holds when that is 0 and the last piece's walk is complete in its part (`Walk.complete`). A parent walk cut by a
  cap beyond its arm's `complete_to_bp` therefore stays uncertified past that depth, however far its child reaches.
- **The evidence block** (`evidence()`, and with every answer): each part's own `evidence_block` (`ops.py:6006`)
  keyed by part, the junctions it depends on with verdicts and causes, the identity statement, and an `outcome` that
  takes the **weakest** value per dimension over every part and junction the answer reads. A `qualified` junction
  maps to `walks: partial` and `label_evidence: qualified`. `exact` holds only when every part read is exact, every
  junction is `exact`, and the identity is verified, so the CLI fixtures, which have no `index_fp`, can never read
  `exact` there (§15 T4). Each part keeps its `K` records, its `Q` (`resource_stop`), its dropped seed labels (`X`),
  its cut lists (`label_lists`, `inexact_counts`), its dropped branch events (`evidence.complete_to_bp`), its
  `coordinates` limitation and its delivery (`spooled`), all listed under its key. None is summed or dropped. A
  local-budget stop on a composite answer is a `Partial` with `complete: false` and a resume token bound to the
  stitched graphlet's digest, as for one graphlet.

### 7.7 Record coordinates

Coordinates are kept per piece and never joined. Under `support: kmer` a retrieval has none (DESIGN §18.1). Under
`trace` the junction is `qualified` (`trace_restart`), and the two pieces of a joined claim are two occurrences unless
their coordinates touch in one record (the parent's run ends at record position e, and the child's run's own bases
start at e in the same record and orientation). Phase 1b states `contiguous: null` for every joined claim; proving
it from the coordinates, and so lifting `trace_restart` for such a junction, is phase 2. `to_fasta(coordinates=)`
states each piece's coordinates in the header, as `Graphlet.to_fasta` does for one walk.

### 7.8 Provenance and replay

`provenance()` returns one entry per part: `{key, role: root | child, junction: {arm, walk, at_bp, overlap_bp} | null,
pair: answered_pair(envelope), identity: index_identity(g), seed_id, digest, request}`. Here `request` is the part's
replay request as `GraphletStore.request_of()` (`store.py:1277`) builds it: the child's request names the
continuation seed and the parent's pair. No single request reproduces a stitched graphlet. Replaying one means
replaying each part (the identity rules apply as for `traverse_fetch(replay=)`) and stitching the results again.
Each child is then checked again against its replayed parent, so a changed index or walk numbering is refused rather
than joined.

## 8. Stitching: determinism and the algebra

- **Pure.** `stitch()` is a function of the parts' texts and envelopes, the walk named or inferred, and
  `allow_unverified_index`. Junction verdicts are recomputed from the parts every time, never stored.
- **Canonical tree.** A stitched graphlet is a tree: root, then a child per continued `(arm, walk)`. Parts are ordered
  by key, comparing steps as `(arm: left < right, walk id)`. Every list answer has a canonical order: arms left then
  right, walks by key under `by='id'`, ties of the other rankings by key. `digest` is the sha256 of the canonical
  tree (the part keys in order with each part's sha256), which is also what the bundle's manifest lists (§9).
- **Commutative across junctions.** `stitch(stitch(P, C1), C2)` and `stitch(stitch(P, C2), C1)` have equal digests
  when C1 and C2 continue different walks. Inference does not depend on the order (§6.2).
- **Associative along chains.** If C2 continues walk v of C1, then `stitch(stitch(P, C1, walk=w), C2, walk=(w, v))`
  and `stitch(P, stitch(C1, C2, walk=v), walk=w)` have equal digests (a StitchedGraphlet child is accepted for this).
- **Not symmetric, by role.** `stitch(C, P)` is refused, deterministically, by the first failing check of §6.2: a
  parent with both arms gives `not_one_arm`, a one-arm parent `not_a_continuation`.
- **Idempotent for identical bodies.** Stitching the same body again at the same walk is a no-op; a different body
  there is refused `already_continued` (§6.2).
- **Identity element.** `stitch(g)` (no child) answers every offered query as `g` does, with keys `(w,)` for path
  ids `w` (tested on every fixture, §15 T2).

## 9. Serialization and exports

The union's FASTA and GFA are in §5.8. This section holds what both combinations share, and stitching's exports.

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
  in key order with the stored opt-in. A bundle whose tree does not hold gets the same codes as §6.2. `save(load(f))`
  is byte-identical to `f`.

  A pair union is written the same way, with first line `#mgt-bundle 1 union` and a manifest of `parts` (`{tag,
  pair, pair_stated_by, bytes, lines, sha256}` in canonical order), `absent` (the absent records as they were
  given), `seed` and `expected_pairs`. No coverage category, absence value or statement is stored: `load()` rebuilds
  by calling `union_pairs()`, so every statement is derived again, and a bundle cannot carry a stronger one.
- **Stitched FASTA** (phase 1a-S). One record per stitched walk, spelled from the parts without the overlap
  (§7.1), using the exports' streaming spelling. Headers are `to_fasta()`'s plus `walk=3.0`,
  `junctions=<at_bp,...>`, `verdict=<exact|qualified>`, `certified_to=<bp>`. A part without bases is refused as that
  part's own export would be.
- **Stitched GFA** (phase 2). Every part's segments under a part prefix, plus a link from the parent's leaf segment
  to each child root segment with the de Bruijn overlap (k − 1)M, which `to_gfa()` already uses inside a part. The
  child's seed is not a segment of its own. Since n ≥ k, the child's first step reads only its seed. **[unsure]**
  how the link is written when L < k − 1 and the overlap reaches into the root seed.
- **Not offered.** `dump()` / MGT (§3.5), `to_json()` (detail `full` of one server result), JSON views of a composite
  (phase 2: `to_dict()` with the parts' summaries, and the junctions or the statements).

## 10. The two arms of one seed: a later validation helper

The owner (2026-11-15): "lower priority, mostly useful for sanity checking/validation". As an answer, the union of a
`direction: left` fetch and a `direction: right` fetch adds nothing that a `direction: both` fetch does not give. As
a check it is useful: the two one-arm fetches and the two-arm fetch of one seed should agree arm by arm. Where they
do not, either a shared cap explains it, or something is wrong in the engine or the library.

A later item, `validate_arms(both, left, right, *, budget=None) -> ArmCheck` (name provisional):

- **Preconditions.** Stitch's identity checks (§6.2: one index, one pair; a one-sided `index_fp` refused), one
  oriented seed, `both` walking both arms, and `left` and `right` one arm each. Strategy differences other than
  `direction` are reported, not refused.
- **The check.** `compare(both, left, arm='left', mode='walks')`, the same for the right arm, and the seed-level
  records (`S` including `validated_seed_id`, `X`, seed-level `K` including `derivation`), which must be equal.
- **Expected differences.** In one request the arms share caps: `max_steps`, `max_paths`, `max_output_bp`, and the
  time, memory and work budgets. A cap reached in `both` can cut an arm that its one-arm fetch walked further. The
  check attributes a difference to a shared cap when `both` states a resource stop (`Q`) and the difference lies
  beyond that arm's `complete_to_bp` in `both`. Any other difference is a finding.
- **Not a combination.** It returns a report, not a composite, and needs nothing of §5 or §6 beyond the identity
  checks and `compare()`.

## 11. The engine

**No engine change is needed or proposed.** Both combinations read only bodies, envelopes and, for the union, the
`GET /capabilities` document. The union needs no engine feature: walks across graphs are out of its reach, and §5.6
says so. Two engine features would make continuations closer to exact. Both change the `/traverse` request schema;
they are for the owner (OQ7) and are not part of this plan:

1. per-seed-label starting state in a seed (`seeds[].labels` as `{name, loss, branches}`), which removes the
   `switching` and `branches_unequal` causes and the `not_continued` gap for labels that switched in late. The
   reviewer rates it the more useful of the two;
2. an excluded-edge list (the ancestors' (k+1)-mers), which removes `rule_violation` at a request cost of one key per
   ancestor base.

Neither removes `sibling_merge_unknown` or `trace_restart`: those need the whole frontier or the followed
occurrence. A server-side frontier resume contradicts "the backend keeps no frontier between requests"
(graphlets.rst) and is not proposed. Under the project rule, an engine search-path change would go to master first,
as a PR Kulkov reviews.

## 12. MCP and the service

**No MCP exposure in phase 1.** The library's own tool functions (`mcp_tools.py`) gain nothing in 1a. A local tool
over a union or a stitched graphlet, returning a derived handle like `graphlet_subtrie`, would be phase 2 and only
if asked. The service's tools are the service session's design, and this note designs none. What the service
session would need to know:

1. **The library will offer `union_pairs()`** over one seed's records from several pairs of a release, or from
   several databases. The service's shape, one record per graph and seed with none merged, stays right: a union is
   built from those records, not stored in their place.
2. **What an honest union needs from the records**, per seed, all passed in one list (§5.2):
   - every record's body (`graphlet_export(format="mgt")`, then `load()`). The body must keep its `J` line with
     the envelope and the pair tags, or the union refuses it with MissingEnvelope (or the caller states the pair
     with `{pair, graphlet}`). Whether the service's export carries the `J` line is for the service session to
     confirm;
   - every record without a body, with its pair and its seed result (`not_in_graph` with `kmers_present` and
     `kmers`, `labels_not_in_graph`, a failed derivation), as `{pair, result}`, with the seed's `sequence` where
     the service has it (otherwise the coverage says the record was not checked against the seed); and a pair
     whose answer never came as `{pair, error}`. Without these, a graph that lacks the seed looks like a graph that
     was not asked;
   - the pairs the fetch was sent to, passed as `expected_pairs=`. Without them the union states its coverage "over
     the pairs given" only;
   - for a seed that no graph holds, the seed sequence as `seed=`, so the union can still say so.
3. **One orientation per union.** All records of one union share one oriented seed sequence. Records of one seed
   fetched in the other orientation are refused (`different_seed`, reverse-complement variant) until phase 2 flips
   them.
4. **The service publishes no continuation tool today.** Stitched reads in the service would need one, plus the
   parent link (parent record, arm, walk) stored with the child record, plus the bundle for downloads. An agent can
   stitch locally: download both bodies, provided the child was fetched on the parent's pair with the walk's
   continuation as its seed. A child made from a preset plus overrides differs from `next_request()`'s request: on
   the default path the overrides take only `bounds`, `direction`, `exhaustive`, `labels.max_labels_per_node` and
   `labels.seed_label_kind`, so there is no reduced loss budget or branch allowance and no rebuilt `labels.extra`.
   Its junction then states `strategy_differs` or `branches_unequal` where that matters (§7.4); the stitch itself
   works.
5. **For information.** The library's `compare()` answers `comparable: false` for retrievals of different indexes
   and for different seeds. The service's `graphlet_compare` description gives "the same seed in two databases or
   two seeds of one locus" as examples, and both are such pairs. For the first, a union of the two databases'
   records answers what each database carries (§5). The wording is the service's.

## 13. Cost and memory bounds

### 13.1 The union

- **Parts.** A union holds its parts' models: about 5 MB per constrained arm and about 20–30 MB label-free, measured
  (DESIGN §7). r237's 11 parts with two arms each hold about 110 MB constrained, or 440–660 MB label-free.
  `memory_bytes()` is their sum plus the trie.
- **The trie.** Built once, arm by arm, by one depth-first pass over each part's first-parent tree (walks and stubs,
  §5.3). Bases are read from the part's segments, not copied, so the work is O(sum of the parts' bases), linear on a
  comb too. A union segment costs about 200 bytes plus about 120 bytes per pair spelling it, and its bases are a
  reference into a part [estimate]. The trie has at most as many segments as the parts together, plus one per point
  where a pair's spelling ends, merges or turns off: for 11 parts of about 2,000 segments each, about 22,000 union
  segments, 5–10 MB [estimate].
- **Labels.** Per (union segment, pair), the labels of the part segments it covers. A part segment's displayed
  labels do not depend on which walk passes it (`_displayed_profile`, `ops.py:492`), so one sweep over the part's
  first-parent tree computes them for every segment, O(the part's runs + segments). `support_profile()` itself scans
  every run per leaf, so the sweep is a new accessor (§14). Counts are kept; full lists are built on demand
  (`carriers()`, `labels='full'`).
- **Queries.** `walks()`: O(union walks × pairs). `carriers()`: O(depth + pairs × labels). FASTA and GFA stream as
  `to_fasta()` and `to_gfa()` do. Every operation takes `budget=` and charges the parts' cold prices plus the trie
  work, with prices in `budget.py` (the cost model states each one).
- **Bundle.** The parts' texts, the absent records and about 100 bytes per part.

### 13.2 Stitching

- **Parts.** A stitched graphlet holds its parts' models: about 5 MB per constrained arm and about 20–30 MB
  label-free, measured (DESIGN §7). `memory_bytes()` is their sum plus the junction indexes. A parent continued at
  many walks holds as many child models. That is the caller's choice, bounded by `budget=`.
- **Junctions.** Per junction: the continuation check is O(n), with n up to 1000 (§7.1). The rule check builds the
  set of (k+1)-mers of the DAG ancestors of w's leaf (any parent; at a merge all parents end in the same oriented
  node, so any parent's last k bases give the window across the boundary), O(bases of those segments + |seed|),
  about 100 bytes per k-mer in a Python set [estimate], so about 1 MB for a 10 kb chain. It then scans the child's
  bases once per segment, O(bases of the child). For chains, each ancestor's set is looked up in turn and not copied.
  The local check compares (k+1)-mers exactly. The engine compares 128-bit hashes where a (k+1)-mer does not fit 64
  bits (DNA k ≥ 32, protein k > 11; SPEC §6.6), so a hash collision can block a step the local check allows. That is
  ignored, and the docs say so.
- **Queries.** Walks: O(sum of the parts' walks). Spelled walks stream as the exports do (peak about twice the text).
  Claims: O(sum of the parts' claims + the continued labels × the child's runs). A profile is two profiles
  concatenated. Every operation takes `budget=` and charges the parts' cold prices plus the junction work, with prices
  in `budget.py` (the cost model states each one).
- **Bundle.** The sum of the parts' standalone texts plus about 100 bytes per part.

## 14. Code plan

| file | change |
|---|---|
| `composite.py` (new; 1a-U) | the parts container, `CombineRefused`, the shared preconditions (envelope, bases, alphabet, the index identity rules), the weakest-of evidence, digests, the budget plumbing, the streaming speller over pieces |
| `union.py` (new; 1a-U) | `union_pairs`, `PairUnion`, `PairPart`, `PairAbsent`, `PairLabel`, `SegmentSupport`, `UnionSegment`, `PairSpan`, `UnionWalk`, `Coverage`; the input classification, the trie by one pass over each first-parent tree, the spans and their absence, coverage, the statements. The queries delegate to the parts' `ops` (no derivation reimplemented) |
| `stitch.py` (new; 1a-S) | `stitch`, `StitchedGraphlet`, `Part`, `Junction`, `StitchedWalk`, `StitchedClaim`; the preconditions; the whitelist and causes; the rule check; the queries, which delegate to the parts' `ops` and shift |
| `ops.py` | split `_identity` into the index rule and the seed rule (the index rule shared, no behaviour change; 1a-U); a per-segment displayed-label accessor, one sweep per arm (§13.1; `_displayed_profile` scans every run per leaf; 1a-U); expose the (k+1)-mer edge key of SPEC §6.6 (1a-S) |
| `client.py` | none: `TraverseClient.response()` already keeps results without a graphlet with their errors (`client.py:895-911`), which `union_pairs()` reads from a one-seed response |
| `parser.py` | (1c) `load()` dispatches on `#mgt-bundle`; the atomic writer reused |
| `export.py` | the union's FASTA and GFA (1a-U); stitched FASTA (1a-S) |
| `budget.py` | prices for the trie, the union queries, stitching and the stitched queries |
| `__init__.py` | exports |
| `scripts/traversal/graphlet_fixtures.py` | `UNION_FIXTURES` (1a-U): one tiny index split by record into two graphs with one shared record, the unsplit index, and the split index at a second k (P6, P11); `STITCH_FIXTURES` (1a-S): CLI parent, child, uninterrupted triples (§15.2 T4), the child's request committed as JSON and pinned to `parent.next_request()` by a test. The script still does not import the library |
| `api/python/tests/test_traverse_union.py`, `test_traverse_stitch.py` (new); `tests/data/traverse/documents/union/`, `stitch/` (new) | §15 |
| `integration_tests/test_multigraph.py` | one union check on a two-graph server (§15.1) |
| `docs/source/graphlets.rst`, `docs/DESIGN-traverse-graphlet.md` | §16 |

The CLI answers single-graph and has no `index_fp`, so the union fixtures' envelopes get their pair tags and an
`index_fp` added, as W5's envelope tests edit `H` and `J`, or are passed with stated pairs (P4); the absent records
follow SPEC §6.1's shapes (steps 2–4); the integration check reads a real one once.

Not touched: the engine, `mcp_tools.py`, `store.py` (phase 2), MGT v1, the service. Size [estimate]:
`composite.py` 300–400 lines; the union, after the cuts, 700–1,000 lines and about 600 lines of tests; stitching
1,200–1,800 lines and about 1,000 lines of tests; about 250 lines of docs.

## 15. Tests

Stdlib `unittest` with seeded random generation, as in the existing suites; no new dependency.

### 15.1 The union (P1–P12)

- **P1 refusals.** One test per code of §5.2, with the exact text, including every variant of `different_seed`
  (first difference, unequal lengths, reverse complement, an absent record's `length_bp`) and of `duplicate_pair`.
  An identical body on the same pair is a no-op; the same body on two stated pairs makes two parts.
- **P2 the union of one part is that part.** On every committed fixture, both arms: the union walks with a `whole`
  span are exactly the part's walks (matched by spelling), with their lengths, end reasons, labels and `complete`;
  every other union walk ends a stub (`merges`); the coverage is the part's arm certificate.
- **P3 order.** Seeded random permutations of the results give equal digests and byte-equal answers for every query
  and both exports.
- **P4 inputs.** Every element form of §5.2: a `Graphlet`, a one-seed `TraverseResponse` with a graphlet and with an
  error result, `{pair, graphlet}` (CLI bodies under two stated pairs), `{pair, result}` with and without `sequence`,
  `{pair, error}`; a multi-seed response refused `several_seeds`; a union of absent records only, with `seed=`.
- **P5 the trie.** Hand-built cases: two pairs that share a prefix and branch; a pair whose walk ends inside
  another's (`ends`, with its reason); a root leaf (`ends` at 0) against an arm not walked (`none`); a pair that
  turns off (`diverges`); a merge whose non-first parent is another pair's walk (review M2's example: `merges`, not
  `diverges`); segments that break where the set of pairs changes; children in base order. Every union walk is
  `whole` or `merges` in some pair. Construction stays linear on a comb (the pattern of the speed-up tests).
- **P6 the split-graph oracle and the non-claim.** One tiny index split by record into two graphs (two "taxonomic
  groups") with one record in both (a "plasmid"), and the unsplit index, under `keep`, `forbid` and generous caps.
  In both modes, every union spelling is a prefix of a spelling of the unsplit retrieval at the same radius, the
  shared record's walks are carried by both pairs, and a walk of the unsplit retrieval that no single graph holds
  entirely is not in the union, while the coverage sentence and both exports state the non-claim (the fixture is
  built so that such a walk exists). Annotate mode: every walk of the unsplit retrieval whose k-mers all lie in one
  graph is in the union, whole or as a prefix ended by that graph's `dead_end`. Constrain mode: every walk of the
  unsplit retrieval carried by a label whose record lies in graph g is in g's part; a walk whose k-mers are all in g
  but whose only carrier is a record of the other graph need not be (review minor 10).
- **P7 labels.** Every label is shown with its pair; `label({'name'})` resolves for a name in one pair and is
  refused `ambiguous_name` for the shared record's name; `shared_names()` lists it, and the statement counts it.
- **P8 coverage.** Every category of §5.5: absent records `not_in_graph` with none and with some k-mers present,
  `labels_not_in_graph`, a failed derivation (`not_walked`, seed in graph) and `{pair, error}`; a part with a
  `seed_labels` limitation (`max_seed_labels` below the carriers) and one with a `scope` limitation
  (`complete_qualified`); parts walked to different radii and with a resource stop; a one-arm part;
  `expected_pairs=` with a missing pair, and no `expected_pairs=`. The sentences are exact, the default depth is the
  largest radius, and `every_pair_certified_to_bp` is never above any part's `complete_to_bp` and is `None` with a
  qualified pair.
- **P9 exports.** FASTA spellings equal the parts' `spell()`; headers carry every token (spans with their absence,
  `modes=`, lower-bound `+`) and leave empty ones out. GFA: the comment lines, the `pp`, `pl` and `pe` tags, `0M`
  links and left-arm orientation, checked line by line as the existing `to_gfa()` tests do.
- **P10 absence.** Each value of §5.3 on fixtures: annotate mode with `keep` inside the certificate (`graph`);
  constrain mode, and a `united_history` scope (`under_strategy`); a turn-off at or beyond `complete_to_bp`
  (`uncertified`); a `dead_end` (`graph`).
- **P11 mixed strategies.** Parts in constrain and annotate mode, and parts at two k: the union builds and states
  each pair's mode and k in the sentence, on each span and in `modes=`.
- **P12 budgets and memory.** Every operation under `budget=` completes or stops with the same row in every run.
  `memory_bytes()` bounds the traced peak on 11 parts (the pattern of `test_traverse_memory`). `summary()` refuses
  `summary_too_small` below the statements' size and never cuts them.

Beside them, one integration check on a two-graph server (`integration_tests/test_multigraph.py`): a seed held by
one graph only, so the other answers `not_in_graph`, read as §5.2 and §5.5 expect from a real response. The new
graphlets.rst examples run in `run_examples.py` (CI's fixture-checks job).

### 15.2 Stitching (T1–T9)

- **T1 refusals.** One test per code of §6.2, with the exact text, on committed fixtures. Identity cases edit `H`/`J`
  as W5's envelope tests do (a manifest digest on one side, on neither, different ones; a different pair; a different
  release). `different_pair` and `different_graph` are tested without manifests, where they are reachable.
- **T2 identity element.** On every committed fixture, for every offered query, `stitch(g)` answers as `g`, with keys
  mapped.
- **T3 algebra.** Commutativity across junctions, associativity along chains and canonical bytes, as digest
  equality over seeded random construction orders. `stitch(C, P)` gives `not_one_arm` for a two-arm parent and
  `not_a_continuation` for a one-arm parent. An identical body at a continued walk is a no-op; a different body is
  `already_continued`. Inference over continued and uncontinued candidates gives one code in every order.
- **T4 oracle (the semantic gate).** CLI triples on tiny indexes: a parent at radius R, its child from
  `next_request()`, and the uninterrupted retrieval at R + E. The CLI fixtures have no `index_fp`, so T4 asserts
  junction and walk verdicts, not the evidence's `exact`.
  (a) Each whitelist case is pinned: `keep` and `forbid`; `merge` with w the arm's only open walk; an unlimited
  branch allowance; an equal finite one. In each, the stitched walks below w equal the uninterrupted retrieval's
  subtree of w's prefix, per path, and every junction says `exact`. A whitelist entry without its passing fixture is
  not shipped as `exact`.
  (b) One fixture per cause, each asserting the stated cause and `qualified`: `sibling_merge_unknown` (two
  radius-ended walks under `merge` that reconverge beyond R; review B1); `trace_restart`; `strategy_differs` (a
  changed quorum knob, a reset branch allowance); `extra_foreign`; `switching`; `not_continued`; both `seed_reentry`
  sub-cases (a circular record); `split_limit_restarted`; `rule_violation`. Where the uninterrupted retrieval differs,
  the difference lies below a stated cause.
  (c) Promotion candidates (`same_stop`, uniform-shift switching) get fixtures that show their walks equal the
  reference's before they are moved to the whitelist.
- **T5 no laundering.** Over all fixture combinations: every answer's evidence is, per dimension, no stronger than any
  part or junction it reads. Every part's `K` kinds, `Q` and `X` appear under its key. A stitched walk is complete
  only if all its pieces and junctions are. With fixtures whose `H` carries an `index_fp`, the evidence's `exact` is
  checked too.
- **T6 round trip.** In 1a: saving the parts and stitching the loaded parts gives an equal digest. In 1c:
  `save`/`load` byte-identical with equal digests, and `load` refuses a truncated part (`GraphletFormatError` naming
  part and line), a digest mismatch, an unknown manifest member, a tree whose child is not a continuation (the §6.2
  code), and a one-sided identity without the stored opt-in.
- **T7 exports.** Stitched FASTA equals the parts' spellings joined at `L` without the overlap, on both arms in both
  orientations, with `L < n` included (`fork`). Headers carry junctions, verdicts and `certified_to`.
- **T8 budgets and memory.** Every operation under `budget=` completes or stops with the same row in every run. The
  account bounds the traced peak on a chain of comb continuations (the pattern of `test_traverse_memory`).
- **T9 docs.** The new graphlets.rst examples run in `run_examples.py` (CI's fixture-checks job).

## 16. Docs

- graphlets.rst, the intro paragraph (lines 23–30), becomes: "**Graphlets are not merged.** Two combinations are
  offered. The graphlets of one seed answered on several pairs (the graphs of a release, or two databases) can be
  combined by sequence (`Combining one seed across pairs`_). The result keeps every retrieval unchanged and states,
  with every answer, how complete each pair is, which pairs carry each walk and how far, that a pair not following a
  walk did not follow it under its own request, that labels are each pair's own, and that a walk no single graph
  holds entirely is in no part. A continuation can be stitched
  to the walk it continues (`Stitching a continuation`_); the result states, at each junction, whether its walks
  are those of one uninterrupted walk (`exact`) or why that is not established (`qualified`). Every other
  combination is refused or not offered, with the reason (`What is not combined`_)", followed by the existing
  sentences on `subgraph()`, `compare()` and `on_reconverge`.
- graphlets.rst, new sections. "Combining one seed across pairs", after "Comparing retrievals": the five statements
  first, then the inputs, the trie with its spans and absence values, the coverage categories, labels per pair, the
  non-claim with the r237 example, one example on the split-graph fixture, FASTA and GFA. "Stitching a
  continuation", after "Continuations and the next request": offsets, orientation and the continuation's length,
  walk keys, the whitelist and the causes, the certified depth, one example on a committed triple, marked
  experimental. "What is not combined": the table of §4.
- `DESIGN-traverse-graphlet.md`: `composite.py`, `union.py` and `stitch.py` rows in the §5 module table, a sentence
  in §5.3 pointing to this note, and the public-API block.
- No spec change: no wire, record or field changes.

## 17. Risks

### 17.1 The union

1. **A union read as one traversal of the release.** Mitigated by the five statements leading every answer, export
   and `print()`, `scope=within_one_graph` on every FASTA record, and the coverage sentence ending with the
   non-claim.
2. **A pair read as complete when it was not asked or failed.** Without every result and `expected_pairs=`, a graph
   that was not asked looks like a graph without the seed. Mitigated by the categories `failed`, `not_walked` and
   `not_included`, by naming the absent records not checked against the seed, and by "over the pairs given" on
   every statement when `expected_pairs=` is missing. The caller must pass them (§12).
3. **A turn-off read as absence.** Mitigated by `absence` on every span that is not `whole` or `merges`, the third
   statement, and each pair's mode in the sentence, on the spans and in `modes=` (review M1).
4. **Equal names read as one record.** Mitigated by showing every label with its pair, refusing a bare name held by
   two pairs (`ambiguous_name`), and counting shared names in the labels statement.
5. **A header name on a `count_only` pair names another record.** Not a union matter: the tag does not repair it,
   and the library states it for no graphlet, since it reads no `GET /capabilities` certificate. Reading each pair's
   `seqs_order` and stating it with header labels would be a library-wide change, outside this note (§19).
6. **Reviewed once.** Fable reviewed the union at `42023784f` (§20.2). The fixes for M2 (stubs and `merges` spans)
   and M5 (one list of results) are new design that no reviewer has seen.

### 17.2 Stitching

1. **A stitched walk read as one traversal.** Mitigated by the verdict on every row and junction, the per-walk
   `certified_to_bp`, the summary leading with the junctions, and `qualified` for anything not on the pinned
   whitelist.
2. **An `exact` that is not exact.** The first version had one (review B1: siblings under `merge`). Mitigated by the
   narrow whitelist, by T4(a) pinning every entry against the uninterrupted retrieval, and by the rule that a
   whitelist entry or promotion ships only with its passing fixture.
3. **Mostly `qualified` in practice.** Under the constrain defaults (`merge`, finite splits, `forbid`, no branches)
   `exact` needs w to be the arm's only open walk, and `seed_reentry` qualifies every walk that meets a repeat of the
   parent's walk. That is the honest answer. Promotions (`same_stop`, uniform-shift switching) narrow it once pinned.
4. **Two query families drift apart.** The composite delegates every derivation to the parts' `ops`, and T2 compares
   `stitch(g)` with `g` on every fixture for every query, so a change in `ops.py` that the composite misses fails T2.
5. **Scope creep.** compare, subgraph, routes and GFA on composites are requested next. They are listed as phase 2 or
   refused with a reason, never half-done.
6. **Memory with many children.** Bounded by `budget=`, stated by `memory_bytes()`.
7. **A second file format (1c only).** The bundle is the library's (like `J`), versioned in its first line, and holds
   only unchanged MGT texts plus keys and digests. It is deferred until a user needs one file.

## 18. Phased plan

- **Phase 0 (this note).** Fable reviewed stitching (§20.1) and the union (§20.2). The owner answered OQ1 ("Phase
  1a"), OQ2 and OQ3 (2026-11-15); OQ4–OQ10 and OQ12 are open.
- **Phase 1a-U, the cross-pair union, first.** `composite.py` and `union.py`, cut to what looking at one seed
  across r237's 11 graphs needs: `union_pairs()` over one list of results, every refusal, the trie with stubs and
  spans, `absence`, coverage with the qualified and not-walked categories, labels per pair, walks, spell, carriers
  of a walk, FASTA and GFA with the statements, summary and `print()`, evidence and provenance; P1–P12 and the
  integration check; the split-graph fixture; the docs. No new file format. One lane, its own review, then a merge
  train like every library change.
  **Why the union ships first.** It has a named user and a need today (the owner, on r237's 11 graphs; OQ3). It
  states nothing beyond what its parts state, so no part of it waits for an oracle, while stitching's `exact` waits
  for T4's fixtures. And it builds the shared machinery (`composite.py`), which makes stitching smaller.
  **Cut from 1a-U** (review cut list): `certificates=` and every name rule (M3); `routes=` (replaced by stubs, M2);
  a union-level `compare()` (each common pair's parts compare with their own `compare()`); `support_by_depth()`;
  `carriers(sequence=)`; a per-pair FASTA (it equals the parts' own exports); a `PairUnion` as an input;
  `walks(by='labels')`.
- **Phase 1a-S, stitching, after the union (approved, OQ1).** `stitch()` and `StitchedGraphlet` with walks, spell,
  continuation and next_request delegation, FASTA, completeness, evidence, summary and provenance; every refusal;
  the two-verdict junction statement with every cause of §7.4; T1–T9 (T6 in its 1a form); the CLI triples; the
  docs, marked experimental. A `StitchedGraphlet` is then accepted as a union part (§5.2).
- **Phase 1b (when a user needs label-level answers across junctions).** Stitched claims (constrain joined, annotate
  shifted), support profile and changes.
- **Phase 1c (when a user needs one file).** The bundle for both kinds (§9) and T6's 1c form.
- **Later, each only with its T4 fixtures.** Promotions to `exact` (`same_stop`; uniform-shift switching); a
  `lower_bound` verdict under M1's conditions.
- **Later, a validation helper.** `validate_arms()` (§10), when the owner wants it.
- **Phase 2 (each item only if wanted).** The union: the items cut from 1a-U; flipping a reverse-complement part on
  primary graphs; `next_request()` per pair (each part's own request for the walks it follows); a JSON view.
  Stitching: GFA; label walks, routes and label summary over stitched walks; annotate claims joined; trace
  contiguity from coordinates. Both: store entries and local tool functions if the service asks.
- **Not planned.** Unions of retrievals of one pair, combinations of different seeds, any engine change.

## 19. Open questions for the owner

**Answered on 2026-11-15.**

- **OQ1** (stitching: a recipe, phase 1a, or phases 1a and 1b): "Phase 1a". Taken: stitching is built as phase
  1a-S, after the union (§18).
- **OQ2** (the union of two arms): "lower priority, mostly useful for sanity checking/validation". Taken: a later
  validation helper (§10).
- **OQ3** (cross-index composition): "why would I refuse it? It's the only way I can do it today. Referring to a
  non-existing and complex engine feature does not help", and "it won't be exact, I know, but the best we can do
  today, I think". Taken: the cross-pair union is offered and ships first (§5, §18), stating its limits with every
  answer.

**Withdrawn.** OQ11 (which rule joins label names across pairs on r237) is moot after review M3: under the server's
certificates no name occurs in two pairs, so no rule joins anything, and the union shows every label with its pair
(§5.4).

**Open**, each with my recommendation, and the reviewer's view where it differs. OQ4–OQ7 are unchanged; the reviewer
of the union saw OQ8–OQ11.

- **OQ4.** Combinations of different seeds of one locus (a locus map) left out. **Recommendation: confirm.** The
  reviewer: confirm.
- **OQ5.** Names: `stitch`, `StitchedGraphlet`, `CombineRefused`, with "merge" kept for reconvergence (`union_arms`
  is moot). **Recommendation: accept.** The reviewer: the names are fine.
- **OQ6.** Without `index_fp` on either side, allow stitching and state it unverified (as `traverse_continue` does
  today), or refuse? The CLI fixtures and some deployments have no manifest. **Recommendation: allow and state
  unverified**; the evidence then never reads `exact`. The reviewer: the same.
- **OQ7.** The engine features of §11. **Recommendation: not now**; only if exact continuations become a requirement
  (master first, Kulkov's review). The reviewer: the same, adding that if ever, per-label starting state removes more
  causes than an excluded-edge list.
- **OQ8.** Names for the union: `union_pairs`, `PairUnion`, `PairLabel`, `PairSpan`; the span values `whole`,
  `ends`, `merges`, `diverges`, `none` and the absence values `graph`, `under_strategy`, `uncertified`; the
  categories `complete_qualified` and `not_walked`; `every_pair_certified_to_bp`; labels shown as `<tag>/<name>`; the
  FASTA token `scope=within_one_graph`. **Recommendation: accept.** "Union" is free in this codebase, unlike
  "merge", and "pairs" says what is united. The reviewer: accept, with these renames and `merges`.
- **OQ9.** Mixed strategies across parts (constrain with annotate, `kmer` with `trace`, different k) are allowed and
  stated: each pair's mode in the coverage sentence, on each span and in `modes=`, and no cross-pair label total at
  all. The alternative refuses mixed modes and supports. **Recommendation: allow and state.** Everything is per
  pair, and the owner listed strategies among what may differ. The reviewer: the same, with the mode statements.
- **OQ10.** Two answers on one pair are refused (`duplicate_pair`), and `compare()` reads two retrievals. The
  alternative keeps both under one pair. **Recommendation: refuse**, with the no-op limited to the same body on the
  same pair. Their union answers no request (§3.2), and a pair's coverage must come from one certificate. The
  reviewer: the same.
- **OQ12 (new).** The review's cut list removes two items the design brief named: "support by depth across pairs"
  (`support_by_depth()`) and a union-level `compare()`. For one seed across r237's 11 graphs, the coverage at a
  chosen depth and the spans of each walk say how far each graph carries each flank, and each pair's two
  retrievals compare with the parts' own `compare()`. **Recommendation: accept the cut**; both return in phase 2 if
  the owner wants a depth profile or to compare two unions (a rebuild, two databases) in one call.

**For the owner's information** (no decision needed for the union): a header label on a pair whose `seqs_order` is
`"count_only"`, as on every r237 pair today, rests on a column order the load did not check, so its name can name
another record even within one pair (envelope spec §9.1). The library states this for no graphlet, since it reads no
`GET /capabilities` certificate. Stating it would be a library-wide change (§17.1, risk 5).

## 20. The reviews and how each finding is settled

Every finding was checked against the code and specs at `bc9684559` before it was settled.

### 20.1 Fable, of `dfd97d9cf` (stitching)

| finding | checked against | settled |
|---|---|---|
| B1 `exact` unsafe under `merge` with another open walk | SPEC §6.5 (merge at one node and one distance; histories united), §6.4 consequences | **accepted, made stricter**: W6 and the cause `sibling_merge_unknown`. Besides "no other walk with a continuation", W6 also requires no other walk longer than L (one that ended deeper than a capped w could also meet its descendants). The old `scope` cause is replaced: merges on w's own route are covered by W8's check over all DAG ancestors (SPEC §6.6, "any parent"). T4(b) fixture added |
| M1 `lower_bound` not monotone | SPEC §6.4 (splits, ambiguity, greedy re-minimisation), §6.5 (merges) | **accepted**: `lower_bound` deferred; two verdicts (§7.4); M1's conditions recorded for a later phase |
| M2 no cause for `trace` | SPEC §6.3 (`trace_break`, one occurrence chain), §6.1 (seed validated on any consecutive occurrence) | **accepted**: W1 and `trace_restart` |
| M3 strategy differences not judged | the strategy fields of SPEC §5 | **accepted, list adjusted**: W2 with the reviewer's allowed list, the resource caps spelled out, `max_live_paths` only without a beam, and `labels.max_labels_per_node` added (it caps output lists only, SPEC §5); `strategy_differs` and `extra_foreign` added |
| M4 the continuation's length | `walker.cpp:3117-3128`, SPEC §7.1; the fixtures' (L, n) | **accepted, wording narrowed**: the formula is in §7.1, and §7.4 W8 uses which nodes the child's seed covers. "Usually the whole flank" holds only without a switch or a route join; the fixtures show both (`same_name` (40, 40); `merge` (100, 38), `caps` (126, 23), `switch_chain` (80, 20)) |
| minor 1 two sub-cases of seed re-entry | SPEC §6.6 precedence | **accepted**: `same_stop` and `earlier_stop` (§7.4); both `qualified` in 1a; `same_stop` is the first promotion candidate |
| minor 2 `max_switch_sources`, branch-count tie-breaks | SPEC §6.3 (tie-break, source cut), §6.5 (merge minimum) | **accepted, settled by the whitelist**: W3 requires `forbid`, under which no source list or switch predecessor exists, and W4 makes branch counts shift uniformly within each label |
| minor 3 causes need a mode column | `ops.py:2769` (annotate seeds carry no labels) | **accepted**: the modes column in §7.4 |
| minor 4 `pruned` omitted | SPEC §7.1; `is_resource_stop` includes `BEAM_PRUNED` (`traversal_types.cpp:59`) | **accepted**: §7.5 |
| minor 5 joined loss is an upper bound | SPEC §6.3 | **accepted**: §7.5 (0 under an `exact` junction; the chosen route's loss otherwise) |
| minor 6 `stitch(C, P)` code | §6.2's order | **accepted**: §8 and T3 |
| minor 7 inference vs continued candidates | — | **accepted**: §6.2, all candidates count |
| minor 8 reachability of `different_pair` / `different_graph` | SPEC §10.3 | **accepted as stated reachability, order kept**: a fingerprint mismatch is proof; pair tags can name one index twice (§6.2) |
| minor 9 CLI fixtures have no `index_fp` | `fork.mgt` `H` (`*`) | **accepted**: T4 asserts verdicts only; T5 checks `exact` on fixtures with an `index_fp` |
| nit: hashed edges for DNA k ≥ 32 | SPEC §6.6 | **accepted**: §13 |
| nit: `continuation()` shift formula | `ops.py:2275-2279` | **accepted**: §7.1 |
| nit: tip and bubble windows | SPEC §6.5 | **accepted**: §7.4 |
| nit: `to_fasta(coordinates=)` | `model.py:571` | **accepted**: §6.1, §7.7 |
| TQ1 causes incomplete | see B1, M2, M3, minor 2 | **accepted**; hairpins need no cause (§7.4) |
| TQ2 seed re-entry | see minor 1, M4 | **accepted**; trace restart and lost sibling merges are the other artifacts (both causes now) |
| TQ3 the rule check's window | — | **accepted with one change**: any parent's window; collisions ignored and documented (§13). The reviewer read W8's ancestor set as qualified by `scope`; here it is the united DAG-ancestor set, and only siblings beyond L qualify (B1) |
| TQ4 the bundle | — | **accepted**: one file, `load()` dispatching, only in 1c |
| TQ5 claims order | — | **accepted** as written |
| TQ6 annotate claims shifted | — | **accepted**; `support_changes()` adds no label at the junction (§7.5) |
| TQ7 identical child | — | **accepted**: a no-op (§6.2, §8) |
| scope: union of two arms has no use case | U2 | **accepted** in the second version (not offered); the owner later made it a validation helper (§10, OQ2) |
| scope: two verdicts in 1a | M1 | **adopted** (§7.4, §18) |

The reviewer agreed with refusing cross-index composition (OQ3). The owner reversed that on 2026-11-15, and the
union was reviewed separately (§20.2).

### 20.2 Fable, of `42023784f` (the union)

Verdict "implement after fixes": the container of unchanged parts is right; one blocker and three majors let a union
read stronger than its parts. The reviewer's checks that came out sound were read again too: positions compare
across k and regime (SPEC §6.1), the seed check matches `_identity` (`ops.py:3600`), per-walk `complete` matches
`_walk_end` (`ops.py:1463`), and a segment's displayed labels do not depend on the walk (`ops.py:492-545`), which
removes the second version's **[unsure]** in §13.1.

| finding | checked against | settled |
|---|---|---|
| B1 `complete` while the part's `outcome.walks` is `partial` | SPEC §7.0 (a `seed_labels` limitation and a `scope` with `observed` > 0 make the walks partial and leave `complete_to_bp`); SPEC §5 (`max_seed_labels` 1000) | **accepted**: `complete` needs the certificate and no walk-class limitation on the arm or the seed; `complete_qualified` names `seed_labels` or `scope`; `every_pair_certified_to_bp` is `None` with a qualified pair; P8 (§5.5) |
| M1 a turn-off read as absence | SPEC §6.10 (constrain follows only label-supported walks), §6.1 step 4 (each pair's own permitted set), §6.7 ("never evidence of absence") | **accepted, one case added**: `absence` on every `ends`, `diverges` and `none` span; `ends` with `dead_end` is `graph` in either mode, since the graph has no next k-mer at all (`END_CLASS`, `ops.py:634`); the third statement; each pair's mode in the sentence, on the spans and in `modes=`; P10 (§5.3) |
| M2 merge parents dropped, a false `diverges` | `derive.chain` (`derive.py:258`), `_first_paths` (`ops.py:641`) | **accepted**: stubs by default from the first-parent tree, `how: merges` with `into`. Every segment lies on a walk or a stub: a segment on no walk reaches, through children whose first parent it is, a segment that is only a non-first parent, and lies on that segment's chain. `routes=` and `segments_off_spellings` are dropped; P2, P5 (§5.3) |
| M3 the certificate rule joins nothing | envelope spec §9 (`headers_unique`: no repeat in a pair or across pairs); SPEC §10.3 (`columns_disjoint`) | **accepted**: `certificates=` and the name rules dropped; labels always shown with their pair; `label({'name'})` resolves for one pair or is refused `ambiguous_name`; `shared_names()`; the labels statement reworded; P7 replaced; OQ11 withdrawn (§5.4) |
| M4 cross-pair totals count one record twice | the note's own text | **accepted, and mostly moot** after the cuts (`support_by_depth()`, `walks(by='labels')`): no cross-pair total is offered; counts are per pair, as label entries, `+` marking a lower bound (§5.7, §5.8) |
| M5 absent records cannot be checked; the parts/absent split | SPEC §7.1 (the `seed` object: `seed_id`, `length_bp`, labels, no sequence); SPEC §4.3 (a seed id hashes the labels); `single_answer()` raises on a refused entry (`client.py:274`); `TraverseResponse.graphlets` holds `None` (`attempts.py:363`, `client.py:908`) | **accepted**: one list of results that the library sorts (§5.2); `{pair, graphlet}` states a pair (CLI bodies, databases without a manifest); `{pair, error}` for an answer that never came; `length_bp` and an optional `sequence` checked, other absent records named as not checked; `not_expected` covers absent records; `bad_absent`, `absent_has_part` (now `duplicate_pair`) and `absent_other_seed` removed; §12 item 2 names the `J` line |
| minor 1 a failed derivation is not "nothing known" | SPEC §6.1 steps 2 and 4, "Per-seed failure" | **accepted as one category**: `not_walked`, stating "the seed is in this graph" for a derivation cause; other failures without arms (a budget that does not hold the seed) state their cause only |
| minor 2 the incomplete reason | SPEC §6.10 ("two boundaries, one per question"); `Arm.cap_trigger` (`model.py:276`) | **accepted** (§5.5) |
| minor 3 the depth asked | — | **accepted**: the default is the largest radius among the parts (§5.5) |
| minor 4 the non-claim's wording | — | **accepted** everywhere: "a walk that no single graph holds entirely" |
| minor 5 the no-op | — | **accepted**: same pair and same body (§5.2, §5.9) |
| minor 6 the pair key | `_identity` (`ops.py:3588`: `index_ns`, release, `index_meta_fp`, `index_fp`) | **accepted** (§5.2) |
| minor 7 refusal order; unequal seed lengths | — | **accepted**: commutativity covers successful unions; a variant for unequal lengths (§5.2, §5.9) |
| minor 8 the certificates digest | — | **moot** under M3 |
| minor 9 `max_bytes` against "never cut" | — | **accepted**: `summary_too_small`, default `max_bytes=4096` (§5.7) |
| minor 10 P6 holds only in annotate mode | SPEC §6.1 step 4 | **accepted**: P6 split by mode (§15.1) |
| minor 11 per-leaf insertion is quadratic on a comb | `walk_rows()` (`ops.py:1473`) | **accepted**: one depth-first pass over each first-parent tree (§5.3, §13.1) |
| minor 12 memory counted one arm | DESIGN §7 | **accepted** (§13.1) |
| minor 13 a union of absent records only | — | **accepted**: allowed with `seed=`; `no_seed` otherwise (§5.2) |
| minor 14 `count_only` is not a tagging matter | envelope spec §9.1; the library reads no certificate (no `headers_unique` or `seqs_order` in `api/python/metagraph/traverse/`) | **accepted**: stated in §5.4 as each pair's own matter that the tag does not repair; a library-wide gap for the owner (§17.1, §19) |
| nits | — | **accepted**: every union walk is `whole` or `merges` in some pair (stubs extend the reviewer's "whole"); a root leaf is `ends` at 0, `none` an arm not walked; the column extension is gone; `every_pair_certified_to_bp` |
| the cut list | the owner's use case: one seed across r237's 11 graphs | **accepted in full** (§18). Two cut items were in the design brief (`support_by_depth()`, a union-level `compare()`), so OQ12 asks the owner to confirm |
| the views on OQ8–OQ11 | — | **adopted**: OQ8 with the renames, OQ9 with the mode statements, OQ10 with the narrower no-op, OQ11 withdrawn (§19) |
