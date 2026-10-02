# Review request: the graphlet stage-1 implementation (MGT v1 writer, Python library, conformance)

**What I want from you:** an adversarial review of the *implementation* of the design you approved
(`docs/DESIGN-traverse-graphlet.md` v5.2) — before MGT v1 is declared frozen. Branch `gr/labeled-traversal`, head
`0d85b371` or later. The stage-1 commits are `3ecbfc47`…`5fbd9057` (read their messages: each states what it claims);
`f9e8be97`, `96cbc410`, `0876d366` are the stated-limitations work that precedes them. The staged enforcement of §14.1
(two-phase commit, budget-aware decoder, ledger) is **not** part of this stage and is not under review, except where
the implementation claims a guarantee it does not deliver.

The freeze criteria you set are met on paper: the golden vectors and the round-trip fixtures pass, and the size was
measured on the SRA index (design §7: 11.5 % / 9.4 % / 4.7 % of compact JSON). Your review decides whether that is true
in substance. **Please rank findings by whether they would force a change to the MGT v1 wire format** — those block
the freeze; everything else can follow.

## Priorities

1. **Wire-format conformance (blocks the freeze).** The C++ writer (`src/cli/traverse.cpp`, the `GraphletWriter` and
   the `mgt::` codec) against design §2 as amended through §17 and spec §7.5. Are records, fields, document order,
   the `*` rules and the canonical form exactly as specified? Is there output the reader would accept but another
   correct reader would interpret differently? Are the decisions recorded in the golden-vector header (canonical-only
   token decoding, RANGES collapsing runs of two, SETEXPR delta parts non-empty and base-consistent, percent-escape
   strictness, front coding on code points, `K s:` escaping beyond §2.2's wording, `i:` range) the right ones to
   freeze?
2. **Round-trip honesty.** The claims are `parse(text).dump() == text`, `is_canonical(text)`, and
   `Graphlet.to_json() ==` the server's `detail: full` result after the §2.5/§9 normalisation, on 29 real requests
   (two indexes) plus 10 CLI-generated fixtures (`api/python/tests/data/traverse/documents/`). Find a walker result
   (any regime, mode, cap, merge, switch, hairpin, trace, both arms, `sequences: false`) where one of them fails, or
   where the normalisation hides a real difference. The known exception: `time_budget` stops (walk_domain `observed`
   is elapsed ms, so T37's oracle cannot hold there) — is excluding them sound?
3. **The walker change** (`7be7c80e`, design §4): `LabelRun::segment`, `branches`, `loss`. Are they set at every
   place a run ends or is closed, with the right values (terminal, valid at `to_bp` only)? Does anything else in the
   walker now behave differently? (`EveryRunIsAnchoredWhereItEnds` is the test.)
4. **Evidence semantics in the library** (`api/python/metagraph/traverse/derive.py`, `ops.py`): §5.1's anchored
   evidence, the run-start guard, `route_only` and boundary claims, routes through merge partitions and `prev_run`,
   `label_summary` (both modes) exactly as the walker computes it, both continuation formulas, `compare()` under index
   identity and completeness scope, `prefix_subset` on the left arm. An internal reviewer checked these against an
   oracle built from the source FASTA k-mers; look for what that oracle could not see.
5. **The conservative outcome rule** (design v5.2 §14): each dimension `complete` only when no limitation of its
   class applies; `label_evidence: qualified` for possible overstatement. Implemented on both sides — find a result
   whose `limitations` contain an entry while its dimension still reads `complete`, or the reverse.
6. **The agent tool layer** (`mcp_tools.py`, `store.py`, `client.py`) — not yet registered in the MCP server (that
   lives in another repository). Size ceilings and cursors, `{name, ref}` and tagged selectors, the evidence block,
   error results, file-path confinement, spooled delivery, handle identity vs body deduplication. Security matters
   here: paths, cursor forgery (the cursor carries a MAC with a secret created in the spool), resource use on hostile
   bodies.
7. **Index identity** (`--index-name`, `--index-manifest`, `index_meta_fp`; the mini-refseq builder writes a
   manifest): is the manifest digest computed over the right bytes, canonically, and is "no manifest → unverifiable"
   enforced everywhere a join happens?

## Specific questions

- The writer refuses a label name that is not valid UTF-8 after sanitising names once after the walk (ill-formed
  sequences → U+FFFD, Python's `decode('utf-8','replace')` rule). Is replacing better than refusing the seed? Two
  distinct raw names can collide after replacement — is that detected?
- `G.first_base` under `sequences: false`: the walker drops the first base of the root and of merged segments, so
  the writer emits `*` there and the reader accepts it. Does that lose information the full JSON had?
- Label sets are bounds-checked against the number of `L` records before ranges are expanded (both languages). Is
  there still a way for a small hostile body to make the parser allocate a lot (e.g. many large SETEXPRs, deep
  merges, huge `n` in `C`)?
- The fixtures are regenerated from the CLI by `scripts/traversal/graphlet_fixtures.py --check`; hand-made ones remain
  where the CLI cannot produce the situation. Are the hand-made ones consistent with what the walker *would* produce,
  or could they encode a case the writer can never emit (and hide a writer bug)?
- `graphlet_measure.py` reported the SRA numbers in design §7. Is the comparison fair (compact full JSON vs body;
  gzip vs gzip)? Is the 30 MB Python model for a 3.5 MB beam body acceptable, or a sign of a representation problem?

## Known limitations — do not re-report

The staged enforcement (two-phase head commit, lazy paths, the budget-aware row-diff decoder, the reservation ledger)
is not built; responses state memory/work bounds as soft. There is no streaming parse validator (a parse builds the
model once). `/resolve` still emits non-UTF-8 names raw. `metagraph annotate --index-header-coords` crashes on a
FASTA header with a leading space (core MetaGraph, not this branch). The MCP registration lives outside this
repository.

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: configure with -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*'   # 233 tests (cwd = build/)
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'                    # 38 tests
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'               # 157 tests
python3 scripts/traversal/make_codec_vectors.py --check
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
# the real-format fixture: scripts/traversal/build_mini_refseq.sh ./build/mini_refseq ./build/metagraph scripts/traversal/mini_refseq_accessions.tsv
```

## How to report

For each finding: severity (blocker / major / minor), **whether it forces a wire-format change** (yes / no),
`file:line`, a concrete failure (input → wrong dump, wrong local answer, crash or unbounded resource use), and the fix
you propose. Separate *confirmed* (reproduced), *suspected* (say which experiment settles it) and *design
disagreement*. End with an explicit verdict: freeze MGT v1 now, freeze after named format fixes, or do not freeze.
