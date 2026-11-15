# Design: combining graphlets — stitching a continuation, and what is never merged

**Status.** Design note only (train 7, lane G1; base `bc9684559`, train 6's head with W5's `93044c6ee`). Nothing here
is implemented. The owner chose "design note first" on 2026-11-14 after the Python graphlet review
(`review-python-graphlets-6d97d6e8.md`: "no public merge/union/stitch of graphlets"; its suggested order ends with
"define a merge API separating same-index union, continuation stitching and cross-index composition"). The first
version (`dfd97d9cf`) was reviewed by Fable ("implement after fixes": one blocker, four majors, nine minors). This
version settles every finding; §19 lists each one and how it was settled. Items marked **[unsure]** are my reading
and need a check; §18 lists the open questions for the owner.

All paths are relative to `metagraph/`; library paths are under `api/python/metagraph/traverse/`. Line numbers are
`bc9684559`'s.

## 0. The decision in brief

1. **One combination is offered: stitching a continuation to the walk it continues** (`stitch(parent, child)`). It is
   the only relation with a use case the docs already prescribe (explore, then re-seed from a walk's continuation,
   graphlets.rst "Choosing a retrieval") and with preconditions the library can check from the two bodies alone.
   The smallest version ships first: walks, spelling, FASTA and the evidence, with no new file format (§2.1, §17).
2. **The result is not a graphlet.** It is a `StitchedGraphlet`: the parts' bodies unchanged, a junction per stitched
   walk, every answer computed from the parts and shifted. MGT v1 is frozen and per-traversal records (A counters,
   B growth bins, V branch events, one Q, seed-first label ids) cannot be merged without inventing values. This is
   the view's rule (a view saves its backing body unchanged) applied to several bodies.
3. **Every junction has one of two verdicts.** `exact` means the child's walks below the junction are the walks one
   uninterrupted walk of the parent's request would have. It is given only on a narrow whitelist of conditions
   (§6.4). Otherwise the verdict is `qualified`, with the named causes. A `lower_bound` verdict is deferred, because
   fewer labels or less budget can also lengthen a walk (§6.4, review M1). Verdicts are derived each time from the
   parent's walk end and the child's own envelope, never stored. No answer reads more complete than the weakest part
   or junction it depends on.
4. **Unions are not offered.** Two retrievals of one arm answer different requests. A union of the two arms of one
   seed would be well defined, but a fetch with `direction: both` already serves it.
5. **Cross-index composition is refused.** A traversal stays inside its pair's graph, and label ids are local to an
   index. Aggregating per-pair answers is reporting work (the service already keeps one record per graph and seed)
   and does not combine graphlets.
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
- `branching.on_reconverge: "merge"` merges paths inside one traversal.
- After W5 (`E23/T6/W5/DELTA.md`), continuations, `deepen()` and replays keep the pair `(graph, graph_path,
  annotation_path)` a graphlet was answered on, the MCP wrapper reads the multi-graph envelope, and
  `GraphletView.claims()` keeps the view's arms. graphlets.rst lines 23–30 say "**Merging graphlets is not
  supported.**"

Combining graphlets needs rules for overlaps, orientation, index-local label ids, provenance and incomplete
evidence. Concatenating records does not give a graphlet.

## 2. Use cases: what a user combines, and whether it needs a merge

Sources: graphlets.rst, `DESIGN-traverse-graphlet.md` §5–§6, the library's tool functions (`mcp_tools.py`), and the
service's tool descriptions (`traverse_fetch`, `traverse_results` and the seven reading tools `graphlet_summary`,
`graphlet_labels`, `graphlet_walk`, `graphlet_walks`, `graphlet_sequence`, `graphlet_compare`, `graphlet_export`).

| # | the user has | the user wants | where it arises | needs a combination? |
|---|---|---|---|---|
| U1 | a retrieval to radius R, and the continuation of walk w (child) | one view of the long locus: one walk from the seed to R + E, support changes across the junction, one FASTA record, claims that do not stop at R | the intended loop "re-seed from a walk's continuation" (graphlets.rst); `next_request`, `deepen`, `traverse_continue` | **yes: stitching** (offered, §5–§8) |
| U2 | one seed fetched with `direction: left`, then with `direction: right` | one graphlet with both arms | a budget spent on one arm first | **no**: a fetch with `direction: both` serves it (§9) |
| U3 | one seed and one index retrieved twice (tuned and exhaustive, constrain and annotate, two radii, two `max_labels_per_node`) | which answer holds, and what one has that the other lacks | graphlets.rst "Comparing retrievals" | **no**: `compare()` answers it. The two are answers to different requests, so their union answers none |
| U4 | one seed constrained on label set A, then on set B | the walks of A ∪ B | label batches; a service request per labelled seed over several graphs | **no**: one request naming both sets answers it, with its own certificate (§3.2) |
| U5 | one seed on several pairs of one release (r237: 11 graphs → 11 records), or in two databases | "who carries this flank anywhere in the release" | `traverse_fetch` over a database with named graphs ("a record per graph and seed, none merged") | **no merge possible** (§3.3): a per-pair report; a name-level join only under the header certificates |
| U6 | two seeds of one locus (a gene tiled with seeds) | a map of the locus | the service's `graphlet_compare` description ("two seeds of one locus") | **not a merge of graphlets**: positions and certificates belong to each seed (§3.4); left out |

Only U1 needs a new operation. U2, U3 and U4 are served today. U5 and U6 do not combine graphlets.

### 2.1 Is stitching worth carrying?

The questions follow the checklist in `tools/writing-rules/guides/STEWARDSHIP.md`.

- **Who needs it.** Python users of the documented deepen loop (`TraverseClient.deepen()`, `next_request()`), and
  agents using the library's own tool functions (`traverse_continue`). No user request is on record. The need comes
  from the review, which found the gap, and from the docs, which prescribe the loop.
- **What happens without it.** A user reads parent and child separately and joins them by hand: spell the parent's
  walk, spell the child's walks, place the child at outward position L (§6.1) and drop the n-base overlap.
  Everything needed is public (`continuation().seed_coord`, the `overlap_bp` that `traverse_continue` returns). What
  the hand-made join loses is the honesty: which labels the continuation did not carry, which knobs differ, which
  merges with other walks it could not see, and where the child walked edges the parent's walk had used. Without
  those, a stitched walk reads as one uninterrupted traversal.
- **The cheapest alternative.** A documented recipe with no API. It is rejected as the only answer because the
  junction causes of §6.4 are exactly what a recipe user skips. It remains the fallback if the owner declines (OQ1).
- **The smallest version.** Phase 1a (§17): `stitch()`, its refusals, the two-verdict junction statement, walks,
  spelling, FASTA, completeness, evidence and provenance. No new file format: parts are saved as `.mgt` files today,
  and stitching them again after `load()` is cheap and pure. Claims and support across junctions (1b) and a bundle
  file (1c) come only if users need them.
- **Who keeps it.** The labeled-traversal library's maintainers (this project's lanes). The ongoing cost is the
  tests of §14 and keeping the composite in step with `ops.py`, which T2 enforces.
- **How we know it is used, and how to remove it.** The library has no telemetry, so the evidence is issues,
  questions and code search. The service would know if it ever adds a continuation tool. The API is marked
  experimental in the docs until the oracle (T4) pins every whitelist entry, and is reviewed at the next library
  review. Removal deletes one module and one docs section. Nothing persisted depends on it until a bundle format
  ships, which is one more reason to defer 1c.

## 3. Where the frame is wrong

### 3.1 "Merge" is taken

In this codebase, merge already means reconvergence: `on_reconverge: merge`, a segment with two or more parents
(`Segment.is_merge`), a run closed by a merge (`m`, `Run.merged`), `completeness_scope: united_history`. A
`Graphlet.merge()` that combines retrievals would sit next to "a merge before the depth" in every compare note. The
operation is therefore named `stitch`. "Merging graphlets is not supported" stays true in the docs.

### 3.2 Same-index union, "overlapping or not", is mostly ill-posed

Which walks a traversal holds depends on the whole strategy: the permitted set, the change cost and loss budget, the
branch allowance, the quorum knobs, the beam and caps, and (under `merge`) which edge histories were united and
which parent is displayed (DESIGN §11). Two retrievals of one arm are two answers to two requests. Their union is
the answer to neither, and no `complete_to_bp` describes it. SPEC §6.4 states the one case where a union is exact:
"Equality with the union of independent single-label walks holds only under `on_reconverge: keep`". In my reading it
also needs no switching, quorum or beam **[unsure]**. A request that names both label sets answers that case with its
own certificate. So:

- **one seed, one arm, two retrievals**: not offered. Keep the retrieval that `compare(mode='prefix_subset')` shows
  to contain the other, or fetch once with what both asked.
- **one seed, disjoint arms**: well defined, because arms share no used-edge sets and each arm certifies itself (SPEC
  §6.6). Not offered either: a fetch with `direction: both` serves the only use case (§9).
- **different seeds**: not offered (§3.4).

### 3.3 Cross-index composition cannot produce a graphlet

A `/traverse` request runs on exactly one pair (`SPEC-multigraph-envelope.md` §1.1). The service says "a traversal
stays inside its graph". A walk whose k-mers lie in two graphs of one release is in neither graphlet, so the union of
per-pair walks is not the walk set of the union graph, and no arm certificate covers it. Labels are index-local:
`c:<column>` and `h:<column>:<seq_id>` are positions in one pair's annotation (SPEC §7.5.7: joined across
retrievals "only when the `index_fp` agree"). Two different pairs never state one `index_fp` (SPEC §10.3), and there
is "no union oracle over several annotations" (DESIGN §16.1). Different k or alphabet makes the bases themselves
incomparable.

What users want in U5 is a per-pair report. The service already keeps one record per graph and seed, none merged. A
join by header name across pairs is valid only where `headers_unique`, `labels_disjoint` and every pair's `seqs_order:
"verified"` hold (envelope spec §9, §9.1, the rule the service's `/resolve` merge applies). That is the service's
report, not a library combination. The route to cross-pair walks is an engine feature ("the graph-less traverse of a
later release", envelope spec §1.1), not a client-side merge.

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
Summing counters or bins would invent values; dropping them would break "lossless" and the `*` rules. MGT v1 is
frozen (SPEC §7.5). A combination is therefore a **container of unchanged bodies** plus a stated relation, which is
how views already work (`J.view`, the backing body unchanged).

## 4. The operations

| operation | status | preconditions (all checked, in this order) | refusal codes |
|---|---|---|---|
| continuation stitching `stitch(parent, child)` | **offered** (phase 1a) | both carry an envelope; same index (§5.2: `index_meta_fp`, `index_fp`, `index_ns`, release, names per ref; one-sided `index_fp` only with `allow_unverified_index=True`); same pair where both state one; same k, regime and alphabet; same mode and support; the child walks one arm; its seed equals the continuation of one walk of that arm of the parent; that walk is not yet continued by another body; in constrain mode its seed labels are alive at the walk's end | `not_combinable`, `view_unsupported`, MissingEnvelope, `different_index`, `index_unverifiable`, `different_pair`, `different_graph`, `different_mode`, `not_one_arm`, `wrong_arm`, `no_continuation`, `not_a_continuation`, `ambiguous_walk`, `already_continued`, `foreign_seed_labels` |
| union of two arms of one seed | **not offered** (§9) | n/a | no function; `stitch()` answers such a pair with `not_one_arm` or `not_a_continuation` |
| union of two retrievals of one arm | **not offered** | n/a | as above |
| cross-index composition (pairs, releases, k, alphabet) | **refused** | n/a | `different_index`, `different_pair`, `different_graph` |
| different seeds (overlapping or not) | **left out** (a different product, §3.4) | n/a | `not_a_continuation` |

## 5. Contract (phase 1)

### 5.1 Python API

New module `stitch.py`, exported from `metagraph.traverse`. Everything is phase 1a except the lines marked `[1b]`
and `[1c]` (§17).

```python
stitch(parent, child=None, *, arm=None, walk=None, allow_unverified_index=False, budget=None) -> StitchedGraphlet
    # parent: Graphlet (with envelope) or StitchedGraphlet; child: Graphlet (with envelope) or StitchedGraphlet
    # rooted at the continuation; None -> the parent alone as a StitchedGraphlet with no junction.
    # arm: 'left' | 'right' | 'l' | 'r' (default: the child's one arm). walk: the continued walk, a path id on a
    # Graphlet parent, a walk key (tuple) on a StitchedGraphlet parent; default: inferred (§5.2). One junction
    # per call; a new object every time (parts shared, never copied or modified).

class StitchedGraphlet:
    parts -> tuple[Part]                  # root first, then by part key (§7)
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
    continuation(arm, walk, *, budget=None) -> Continuation       # the owning part's, seed_coord mapped (§6.1)
    next_request(arm, walks, bp=None, reduce_budget=True, reset_branches=False, **overrides) -> NextRequest
    next_requests(...)                                            # delegate to the part owning the walk(s)
    to_fasta(arm=None, walks=None, with_seed=True, orientation='natural', width=None, *, coordinates=None,
             budget=None) -> str        # coordinates as Graphlet.to_fasta (model.py:571), per piece (§6.7)
    completeness() -> dict                # per arm: the root's complete_to_bp and one extension per junction
    evidence(arm=None) -> dict            # §6.6
    summary(max_bytes=2048, *, budget=None) -> dict
    provenance() -> list[dict]            # §6.8
    save(path, *, budget=None)            # [1c] the bundle, §8; metagraph.traverse.load(path) returns it
    digest -> str                         # sha256 of the canonical tree: the part keys and the parts' digests (§7)
    memory_bytes() -> int                 # sum of the parts' models + the junction indexes
    # Deliberately refused (CombineRefused 'not_offered', text in §5.2): dump(), to_json(), compare(),
    # subgraph(), to_gfa() (phase 2), label_walks(), routes(), label_summary(), splits() (phase 2)

@dataclass(frozen=True)
class Part:      key: tuple; graphlet: Graphlet; offset_bp: int; digest: str
                 # key () = the root; (('right', 3),) = the child of root walk 3; (('right', 3), ('right', 0)) ...
@dataclass(frozen=True)
class Junction:  arm: str; walk: tuple; at_bp: int; overlap_bp: int; part: tuple;
                 verdict: str                                   # 'exact' | 'qualified'
                 causes: tuple[tuple[str, str], ...]            # (code, sentence), in the order of §6.4's table
                 continued: tuple[Label, ...]; not_continued: tuple[Label, ...]
@dataclass(slots=True)
class StitchedWalk:  key: tuple; arm: str; length_bp: int; sequence: Optional[str]; path_reason: Optional[str];
                     end_reasons: dict; labels_full: list[Label]; n_alive: int; complete: bool;
                     verdict: str                               # 'exact' | 'qualified': its junctions and own causes
                     certified_to_bp: int; beyond_certified_bp: int; junctions: tuple[int, ...];
                     pieces: tuple[tuple[tuple, int, int, int], ...]   # (part key, path_id, from_bp, to_bp)
                     causes: tuple                                      # this walk's own causes (§6.4)
@dataclass(slots=True)
class StitchedClaim(Claim):  walk: tuple; pieces: tuple; junctions: tuple[int, ...]       # [1b]

class CombineRefused(ValueError):  code: str     # every refusal of §5.2; MissingEnvelope and BadSelector as today
```

A walk key is the tuple of walk ids along the chain: `(3,)` is root walk 3, and `(3, 0)` is walk 0 of the child that
continues root walk 3. Positions in every answer are outward from the root seed's boundary. Each row also carries its
pieces in their own parts' coordinates. A label selector `{'id': n}` is refused on a stitched graphlet (BadSelector,
"label ids are local to each part of a stitched graphlet: select by ref or name"). Refs and names resolve over the
union.

### 5.2 Refusals: exact texts

`%s` and `%d` are filled as stated; a label list prints at most 8 `{name, ref}` then "and N more", as `next_request`'s
notes do.

| code | text |
|---|---|
| `not_combinable` (TypeError) | `stitch() combines graphlets (Graphlet or StitchedGraphlet), not %s` |
| `view_unsupported` | `stitch() reads whole retrievals, not views: stitch view.backing instead` |
| MissingEnvelope | today's `require_envelope('stitch()')` text: `stitch() needs the response envelope (seed_id, the normalized strategy, annotation counters): build the graphlet with Graphlet.from_response() or load a saved .mgt file with its J line; a parsed body alone does not carry it` |
| `different_index` | `the child was retrieved from another index than its parent (%s): a continuation runs on its parent's index, and graphlets of different indexes are never combined -- their label ids are local to each index and their walks to each graph`, `%s` one of `index_meta_fp differs`, `index_fp differs`, `index namespace %r vs %r`, `release %r vs %r`, `label %s is %r in one and %r in the other` |
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

## 6. Semantics of stitching

### 6.1 Coordinates, orientation and the overlap

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
nodes of the parent's walk from L − n on (with root-seed nodes when n > L). §6.4 uses this.

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
spelled `ATTGCTAAGA + child.spell('right', v)`. `fork` is also a case for §6.4's `sibling_merge_unknown`: it runs
`merge`, and its right arm has a second walk with a continuation.

### 6.2 Walk identity and overlap

A stitched graphlet's walks on arm a are the root's walks, except that every continued walk w is replaced by `(w, v)`
for each walk v of its child (recursively for chains). Nothing is renumbered: `(w, v)` names root walk w and child
walk v. Continued walks of one parent with identical tails are not deduplicated, because each child is its own trie;
under `merge` such children are a `sibling_merge_unknown` case (§6.4). The overlap is exactly the child's seed
(`overlap_bp = n`, the value `traverse_continue` stores), so nothing is double counted. The child's positions start
at `L`, where the parent's walk ends.

### 6.3 Labels

Labels join by `ref` across parts. This is valid because the identity check proved one index, or the caller accepted
an unverified one, which is then stated. Names must agree per ref (`different_index` otherwise, as `compare()`'s
names rule). In constrain mode the child's seed labels (its first `num_seed_labels` dictionary entries) must be alive
at the end of walk w (the leaf's end set, `G.end`), otherwise `foreign_seed_labels`. Labels alive at w that the
child does not seed are listed in the junction's `not_continued` (`alive_not_seeded` as `next_request` reports it,
or simply not named). In annotate mode the child has no seed labels (`ops.py:2769`), and none of this applies.

### 6.4 The junction statement: two verdicts

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
oracle has pinned, so phase 1a ships `exact` on a narrow whitelist and `qualified` with named causes for everything
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
  lineage is one label, and `labels.extra` cannot be entered. Switching (any other model) is `qualified` in phase 1a.
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
  `rejoined_seed` on a node of the child's seed that is not a root-seed node (§6.1: a node of the parent's walk
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
| `seed_reentry` | W8: a child `rejoined_seed` on a node of the child's seed outside the root seed. Two sub-cases, told apart by the blocked event's character (`E … b`): `same_stop`, where the blocked (k+1)-mer is an ancestor's edge, so the reference stops there too (with `edge_reuse`, by the precedence of SPEC §6.6) and only the reason differs; `earlier_stop`, a new edge into a node of the walk, where the reference goes on. Common: with the whole walk as the child's seed, any repeat of the walk ends the child here (§6.1) | walk | both |
| `split_limit_restarted` | W8: the splits on w plus the child's exceed `max_splits_per_path`, stated with the position of the first split past it | walk | both |

`same_stop` is the first candidate for promotion: its walks are the reference's and only the end reason differs, so
once T4 pins it, it can be `exact` with the reference's reason restated. Until a whitelist entry or a promotion is
pinned by its T4 fixture, the verdict is `qualified`. The verdict never depends on notes `next_request()` returned: a
child made another way (a service job with preset and overrides, a request written by hand) is judged by what it ran.
Hairpins (`hairpins: skip | follow`) are decided on one (k+1)-mer and behave the same given the same knob, which W2
requires. `tip_window_bp` and `bubble_window_bp` are rejected unless 0 (SPEC §6.5), so W2 keeps them equal; if they
are ever implemented, a non-zero window becomes a `qualified` cause.

### 6.5 Claims and support across a junction (phase 1b)

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

### 6.6 Completeness and evidence: nothing is laundered

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
  `exact` there (§14 T4). Each part keeps its `K` records, its `Q` (`resource_stop`), its dropped seed labels (`X`),
  its cut lists (`label_lists`, `inexact_counts`), its dropped branch events (`evidence.complete_to_bp`), its
  `coordinates` limitation and its delivery (`spooled`), all listed under its key. None is summed or dropped. A
  local-budget stop on a composite answer is a `Partial` with `complete: false` and a resume token bound to the
  stitched graphlet's digest, as for one graphlet.

### 6.7 Record coordinates

Coordinates are kept per piece and never joined. Under `support: kmer` a retrieval has none (DESIGN §18.1). Under
`trace` the junction is `qualified` (`trace_restart`), and the two pieces of a joined claim are two occurrences unless
their coordinates touch in one record (the parent's run ends at record position e, and the child's run's own bases
start at e in the same record and orientation). Phase 1b states `contiguous: null` for every joined claim; proving
it from the coordinates, and so lifting `trace_restart` for such a junction, is phase 2. `to_fasta(coordinates=)`
states each piece's coordinates in the header, as `Graphlet.to_fasta` does for one walk.

### 6.8 Provenance and replay

`provenance()` returns one entry per part: `{key, role: root | child, junction: {arm, walk, at_bp, overlap_bp} | null,
pair: answered_pair(envelope), identity: index_identity(g), seed_id, digest, request}`. Here `request` is the part's
replay request as `GraphletStore.request_of()` (`store.py:1277`) builds it: the child's request names the
continuation seed and the parent's pair. No single request reproduces a stitched graphlet. Replaying one means
replaying each part (the identity rules apply as for `traverse_fetch(replay=)`) and stitching the results again.
Each child is then checked again against its replayed parent, so a changed index or walk numbering is refused rather
than joined.

## 7. Determinism and the algebra

- **Pure.** `stitch()` is a function of the parts' texts and envelopes, the walk named or inferred, and
  `allow_unverified_index`. Junction verdicts are recomputed from the parts every time, never stored.
- **Canonical tree.** A stitched graphlet is a tree: root, then a child per continued `(arm, walk)`. Parts are ordered
  by key, comparing steps as `(arm: left < right, walk id)`. Every list answer has a canonical order: arms left then
  right, walks by key under `by='id'`, ties of the other rankings by key. `digest` is the sha256 of the canonical
  tree (the part keys in order with each part's sha256), which is also what the bundle's manifest lists (§8).
- **Commutative across junctions.** `stitch(stitch(P, C1), C2)` and `stitch(stitch(P, C2), C1)` have equal digests
  when C1 and C2 continue different walks. Inference does not depend on the order (§5.2).
- **Associative along chains.** If C2 continues walk v of C1, then `stitch(stitch(P, C1, walk=w), C2, walk=(w, v))`
  and `stitch(P, stitch(C1, C2, walk=v), walk=w)` have equal digests (a StitchedGraphlet child is accepted for this).
- **Not symmetric, by role.** `stitch(C, P)` is refused, deterministically, by the first failing check of §5.2: a
  parent with both arms gives `not_one_arm`, a one-arm parent `not_a_continuation`.
- **Idempotent for identical bodies.** Stitching the same body again at the same walk is a no-op; a different body
  there is refused `already_continued` (§5.2).
- **Identity element.** `stitch(g)` (no child) answers every offered query as `g` does, with keys `(w,)` for path
  ids `w` (tested on every fixture, §14 T2).

## 8. Serialization and exports

- **Without a bundle (phase 1a).** Each part is saved as its own `.mgt` file, as today. Stitching the loaded parts
  again rebuilds the same stitched graphlet (equal digest), because `stitch()` is pure.
- **The bundle** (phase 1c; `save()`, `load()`; `.mgtb`). One UTF-8 text file, written atomically like `save()`,
  read by `load()`, which dispatches on its first line as it already returns a view for a saved view:

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
  in key order with the stored opt-in. A bundle whose tree does not hold gets the same codes as §5.2. `save(load(f))`
  is byte-identical to `f`.
- **FASTA** (phase 1a). One record per stitched walk, spelled from the parts without the overlap (§6.1), using the
  exports' streaming spelling. Headers are `to_fasta()`'s plus `walk=3.0`, `junctions=<at_bp,...>`,
  `verdict=<exact|qualified>`, `certified_to=<bp>`. A part without bases is refused as that part's own export would
  be.
- **GFA** (phase 2). Every part's segments under a part prefix, plus a link from the parent's leaf segment to each
  child root segment with the de Bruijn overlap (k − 1)M, which `to_gfa()` already uses inside a part. The child's
  seed is not a segment of its own. Since n ≥ k, the child's first step reads only its seed. **[unsure]** how the
  link is written when L < k − 1 and the overlap reaches into the root seed.
- **Not offered.** `dump()` / MGT (§3.5), `to_json()` (detail `full` of one server result), JSON views of a composite
  (phase 2: `to_dict()` with the parts' summaries and the junctions).

## 9. Unions of retrievals: not offered

Neither union is built. For two retrievals of one arm (§3.2) no union is the answer to any request. For the two arms
of one seed the union is well defined (arms share no used-edge sets, SPEC §6.6), but a fetch with `direction: both`
serves the only use case (U2), and the reviewer agrees. If it is ever wanted, its preconditions would be: one index
(as §5.2), one oriented seed sequence, disjoint arm sets, normalized strategies equal except `direction`, and equal
seed-level records (`S` including `validated_seed_id`, `X`, seed-level `K` including `derivation`); each arm answered
from its own part, the outcome the weakest per dimension, and two `Q` records kept apart under their parts.

## 10. The engine

**No engine change is needed or proposed.** Everything above reads the two bodies and envelopes. Two engine features
would make continuations closer to exact. Both change the `/traverse` request schema; they are for the owner (OQ7)
and are not part of this plan:

1. per-seed-label starting state in a seed (`seeds[].labels` as `{name, loss, branches}`), which removes the
   `switching` and `branches_unequal` causes and the `not_continued` gap for labels that switched in late. The
   reviewer rates it the more useful of the two;
2. an excluded-edge list (the ancestors' (k+1)-mers), which removes `rule_violation` at a request cost of one key per
   ancestor base.

Neither removes `sibling_merge_unknown` or `trace_restart`: those need the whole frontier or the followed
occurrence. A server-side frontier resume contradicts "the backend keeps no frontier between requests"
(graphlets.rst) and is not proposed. Under the project rule, an engine search-path change would go to master first,
as a PR Kulkov reviews.

## 11. MCP and the service

**No MCP exposure in phase 1.** The library's own tool functions (`mcp_tools.py`) gain nothing. A local
`graphlet_stitch(parent, child)` returning a derived handle, like `graphlet_subtrie`, would be phase 2 and only if the
service asks. The service session would need to know:

1. The library will offer `stitch()` for a parent and its continuation **on the same pair**, and nothing that joins
   graphlets of different pairs, indexes, releases or seeds. The service's "one record per graph and seed, none
   merged" stays the right shape. A cross-pair join by header name needs `headers_unique`, `labels_disjoint` and
   `seqs_order: "verified"`, as its `/resolve` merge already requires.
2. The service publishes no continuation tool today. Stitched reads in the service would need one, plus the parent
   link (parent record, arm, walk) stored with the child record, plus the bundle for downloads. None of that is
   designed here.
3. An agent can stitch locally: download both bodies (`graphlet_export(format="mgt")`, `load()`), provided the child
   was fetched on the parent's pair with the walk's continuation as its seed. A child made from a preset plus
   overrides differs from `next_request()`'s request: on the default path the overrides take only `bounds`,
   `direction`, `exhaustive`, `labels.max_labels_per_node` and `labels.seed_label_kind`, so there is no reduced loss
   budget or branch allowance and no rebuilt `labels.extra`. Its junction then states `strategy_differs` or
   `branches_unequal` where that matters (§6.4); the stitch itself works.
4. For information: the library's `compare()` answers `comparable: false` for retrievals of different indexes and for
   different seeds. The service's `graphlet_compare` description gives "the same seed in two databases or two seeds
   of one locus" as examples, and both are such pairs. The wording is the service's.

## 12. Cost and memory bounds

- **Parts.** A stitched graphlet holds its parts' models: about 5 MB per constrained arm and about 20–30 MB
  label-free, measured (DESIGN §7). `memory_bytes()` is their sum plus the junction indexes. A parent continued at
  many walks holds as many child models. That is the caller's choice, bounded by `budget=`.
- **Stitching.** Per junction: the continuation check is O(n), with n up to 1000 (§6.1). The rule check builds the
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

## 13. Code plan (phase 1; `parser.py`'s dispatch and the bundle only in 1c)

| file | change |
|---|---|
| `stitch.py` (new) | `stitch`, `StitchedGraphlet`, `Part`, `Junction`, `StitchedWalk`, `StitchedClaim`, `CombineRefused`; the preconditions; the whitelist and causes; the rule check; the queries, which delegate to the parts' `ops` and shift (no derivation reimplemented); bundle save and load (1c) |
| `ops.py` | split `_identity` into the index rule and the seed rule (the index rule shared, no behaviour change); expose the (k+1)-mer edge key of SPEC §6.6 |
| `parser.py` | (1c) `load()` dispatches on `#mgt-bundle`; the atomic writer reused |
| `export.py` | `to_fasta` over a stitched graphlet (the streaming speller over pieces) |
| `budget.py` | prices for stitching and the composite queries |
| `__init__.py` | exports |
| `scripts/traversal/graphlet_fixtures.py` | `STITCH_FIXTURES`: CLI parent, child, uninterrupted triples (§14 T4). The child's request is committed as JSON, and a test pins it to `parent.next_request()`. The script still does not import the library |
| `api/python/tests/test_traverse_stitch.py` (new), `tests/data/traverse/documents/stitch/` (new) | §14 |
| `docs/source/graphlets.rst`, `docs/DESIGN-traverse-graphlet.md` | §15 |

Not touched: the engine, `mcp_tools.py`, `store.py` (phase 2), MGT v1, the service. Size [estimate]: 1,200–1,800
lines of library code, about 1,000 lines of tests, 6–10 CLI fixture triples, about 150 lines of docs.

## 14. Tests

Stdlib `unittest` with seeded random generation, as in the existing suites; no new dependency.

- **T1 refusals.** One test per code of §5.2, with the exact text, on committed fixtures. Identity cases edit `H`/`J`
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
  part and line), a digest mismatch, an unknown manifest member, a tree whose child is not a continuation (the §5.2
  code), and a one-sided identity without the stored opt-in.
- **T7 exports.** Stitched FASTA equals the parts' spellings joined at `L` without the overlap, on both arms in both
  orientations, with `L < n` included (`fork`). Headers carry junctions, verdicts and `certified_to`.
- **T8 budgets and memory.** Every operation under `budget=` completes or stops with the same row in every run. The
  account bounds the traced peak on a chain of comb continuations (the pattern of `test_traverse_memory`).
- **T9 docs.** The new graphlets.rst examples run in `run_examples.py` (CI's fixture-checks job).

## 15. Docs

- graphlets.rst, the intro paragraph (lines 23–30), becomes: "**Graphlets are not merged.** One combination is
  offered: a continuation can be stitched to the walk it continues (`Stitching a continuation`_); the result keeps
  both retrievals unchanged and states, at each junction, whether its walks are those of one uninterrupted walk
  (`exact`) or why that is not established (`qualified`). Every other combination is refused or not offered, with
  the reason (`What is not combined`_)", followed by the existing sentences on `subgraph()`, `compare()` and
  `on_reconverge`.
- graphlets.rst, a new section after "Continuations and the next request": "Stitching a continuation" (offsets,
  orientation and the continuation's length, walk keys, the whitelist and the causes, the certified depth, one
  example on a committed triple; the bundle once 1c ships; marked experimental), and "What is not combined" (the
  table of §4).
- `DESIGN-traverse-graphlet.md`: a `stitch.py` row in the §5 module table, a sentence in §5.3 pointing to this note,
  and the public-API block.
- No spec change: no wire, record or field changes.

## 16. Risks

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

## 17. Phased plan

- **Phase 0 (this note).** Fable review done (§19); the owner's answers to OQ1–OQ7.
- **Phase 1a, the smallest version (if OQ1 says build).** `stitch()` and `StitchedGraphlet` with walks, spell,
  continuation and next_request delegation, FASTA, completeness, evidence, summary and provenance; every refusal;
  the two-verdict junction statement with every cause of §6.4; T1–T9 (T6 in its 1a form); the CLI triples; the docs,
  marked experimental. No new file format. One lane, its own review, then a merge train like every library change.
- **Phase 1b (when a user needs label-level answers across junctions).** Claims (constrain joined, annotate shifted),
  support profile and changes.
- **Phase 1c (when a user needs one file).** The bundle (§8) and T6's 1c form.
- **Later, each only with its T4 fixtures.** Promotions to `exact` (`same_stop`; uniform-shift switching); a
  `lower_bound` verdict under M1's conditions.
- **Phase 2 (each item only if wanted).** GFA; label walks, routes and label summary over stitched walks; annotate
  claims joined; trace contiguity from coordinates; store entries and a local `graphlet_stitch` tool if the service
  asks.
- **Not planned.** Unions of retrievals (one arm or two), cross-index composition, combinations of different seeds,
  any engine change.

## 18. Open questions for the owner

Each with my recommendation, and the reviewer's view where it differs.

- **OQ1.** Which of: (a) a documented recipe only, no API (§2.1); (b) phase 1a; (c) phases 1a and 1b? No user has
  asked for stitching, and the service publishes no continuation tool, so its users would get stitching only through
  downloads and Python, or through the library's own tools (`traverse_continue`). **Recommendation: (b), with the
  two-verdict reduction.** The reviewer: the same.
- **OQ2.** A union of the two arms of one seed. **Recommendation: leave it out** (changed from "phase 2 if wanted"):
  `direction: both` serves the only use case. The reviewer: the same.
- **OQ3.** Cross-index composition refused in the library; per-pair aggregation (with the header certificates) stays
  the service's report. **Recommendation: confirm.** The reviewer: confirm.
- **OQ4.** Combinations of different seeds of one locus (a locus map) left out. **Recommendation: confirm.** The
  reviewer: confirm.
- **OQ5.** Names: `stitch`, `StitchedGraphlet`, `CombineRefused`, with "merge" kept for reconvergence (`union_arms`
  is moot). **Recommendation: accept.** The reviewer: the names are fine.
- **OQ6.** Without `index_fp` on either side, allow stitching and state it unverified (as `traverse_continue` does
  today), or refuse? The CLI fixtures and some deployments have no manifest. **Recommendation: allow and state
  unverified**; the evidence then never reads `exact`. The reviewer: the same.
- **OQ7.** The engine features of §10. **Recommendation: not now**; only if exact continuations become a requirement
  (master first, Kulkov's review). The reviewer: the same, adding that if ever, per-label starting state removes more
  causes than an excluded-edge list.

## 19. The review (Fable, of `dfd97d9cf`) and how each finding is settled

Every finding was checked against the code and specs at `bc9684559` before it was settled.

| finding | checked against | settled |
|---|---|---|
| B1 `exact` unsafe under `merge` with another open walk | SPEC §6.5 (merge at one node and one distance; histories united), §6.4 consequences | **accepted, made stricter**: W6 and the cause `sibling_merge_unknown`. Besides "no other walk with a continuation", W6 also requires no other walk longer than L (one that ended deeper than a capped w could also meet its descendants). The old `scope` cause is replaced: merges on w's own route are covered by W8's check over all DAG ancestors (SPEC §6.6, "any parent"). T4(b) fixture added |
| M1 `lower_bound` not monotone | SPEC §6.4 (splits, ambiguity, greedy re-minimisation), §6.5 (merges) | **accepted**: `lower_bound` deferred; two verdicts (§6.4); M1's conditions recorded for a later phase |
| M2 no cause for `trace` | SPEC §6.3 (`trace_break`, one occurrence chain), §6.1 (seed validated on any consecutive occurrence) | **accepted**: W1 and `trace_restart` |
| M3 strategy differences not judged | the strategy fields of SPEC §5 | **accepted, list adjusted**: W2 with the reviewer's allowed list, the resource caps spelled out, `max_live_paths` only without a beam, and `labels.max_labels_per_node` added (it caps output lists only, SPEC §5); `strategy_differs` and `extra_foreign` added |
| M4 the continuation's length | `walker.cpp:3117-3128`, SPEC §7.1; the fixtures' (L, n) | **accepted, wording narrowed**: the formula is in §6.1, and §6.4 W8 uses which nodes the child's seed covers. "Usually the whole flank" holds only without a switch or a route join; the fixtures show both (`same_name` (40, 40); `merge` (100, 38), `caps` (126, 23), `switch_chain` (80, 20)) |
| minor 1 two sub-cases of seed re-entry | SPEC §6.6 precedence | **accepted**: `same_stop` and `earlier_stop` (§6.4); both `qualified` in 1a; `same_stop` is the first promotion candidate |
| minor 2 `max_switch_sources`, branch-count tie-breaks | SPEC §6.3 (tie-break, source cut), §6.5 (merge minimum) | **accepted, settled by the whitelist**: W3 requires `forbid`, under which no source list or switch predecessor exists, and W4 makes branch counts shift uniformly within each label |
| minor 3 causes need a mode column | `ops.py:2769` (annotate seeds carry no labels) | **accepted**: the modes column in §6.4 |
| minor 4 `pruned` omitted | SPEC §7.1; `is_resource_stop` includes `BEAM_PRUNED` (`traversal_types.cpp:59`) | **accepted**: §6.5 |
| minor 5 joined loss is an upper bound | SPEC §6.3 | **accepted**: §6.5 (0 under an `exact` junction; the chosen route's loss otherwise) |
| minor 6 `stitch(C, P)` code | §5.2's order | **accepted**: §7 and T3 |
| minor 7 inference vs continued candidates | — | **accepted**: §5.2, all candidates count |
| minor 8 reachability of `different_pair` / `different_graph` | SPEC §10.3 | **accepted as stated reachability, order kept**: a fingerprint mismatch is proof; pair tags can name one index twice (§5.2) |
| minor 9 CLI fixtures have no `index_fp` | `fork.mgt` `H` (`*`) | **accepted**: T4 asserts verdicts only; T5 checks `exact` on fixtures with an `index_fp` |
| nit: hashed edges for DNA k ≥ 32 | SPEC §6.6 | **accepted**: §12 |
| nit: `continuation()` shift formula | `ops.py:2275-2279` | **accepted**: §6.1 |
| nit: tip and bubble windows | SPEC §6.5 | **accepted**: §6.4 |
| nit: `to_fasta(coordinates=)` | `model.py:571` | **accepted**: §5.1, §6.7 |
| TQ1 causes incomplete | see B1, M2, M3, minor 2 | **accepted**; hairpins need no cause (§6.4) |
| TQ2 seed re-entry | see minor 1, M4 | **accepted**; trace restart and lost sibling merges are the other artifacts (both causes now) |
| TQ3 the rule check's window | — | **accepted with one change**: any parent's window; collisions ignored and documented (§12). The reviewer read W8's ancestor set as qualified by `scope`; here it is the united DAG-ancestor set, and only siblings beyond L qualify (B1) |
| TQ4 the bundle | — | **accepted**: one file, `load()` dispatching, only in 1c |
| TQ5 claims order | — | **accepted** as written |
| TQ6 annotate claims shifted | — | **accepted**; `support_changes()` adds no label at the junction (§6.5) |
| TQ7 identical child | — | **accepted**: a no-op (§5.2, §7) |
| scope: union of two arms has no use case | U2 | **accepted**: not offered (§9); OQ2 changed |
| scope: two verdicts in 1a | M1 | **adopted** (§6.4, §17) |
