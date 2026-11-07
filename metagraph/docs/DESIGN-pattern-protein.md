# Design: pattern search on protein graphs (E16)

Status (2026-11-07, branch gr/mg-e16-protein-pattern, lane C of the every-route engine work): the engine, the
route and the capabilities are implemented behind `_PROTEIN_GRAPH`, the contract is SPEC-pattern-search.md §29, the
tests of section 9 pass locally on both builds; the branch is not rebased onto E15 yet (section 10: the CI job comes
with that rebase), not reviewed, not deployed.

The owner's decision of 2026-10-31: "patterns must work on UniParc too". The language is the minimum of the suite's plan row E16
(docs/PLAN-multigraph-implementation.md of metagraph-suite) and of the engine's note of 2026-10-31, section 6; this
note settles only what those two leave to the implementation, and lists those choices in section 8.

Protein patterns against DNA graphs (a peptide read through its codon automaton, SPEC §12.2) are unchanged. This
note is about a graph built with the Protein alphabet, served by the Protein build (`metagraph_Protein`,
`-DCMAKE_DBG_ALPHABET=Protein`, `_PROTEIN_GRAPH`): UniParc at k = 15. Every change is behind that compile-time
switch; the DNA build compiles the same code it compiled before (section 7).

## 1. The graph

- The BOSS alphabet is `$ABCDEFGHIJKLMNOPQRSTUVWYZX`: every letter A to Z. Twenty-two are residues read as
  themselves (the 20 amino acids A C D E F G H I K L M N P Q R S T V W Y, U selenocysteine, O pyrrolysine); four are
  ambiguity letters a record may hold (B, J, Z, X; any character outside the alphabet is stored as X).
- One strand: a protein graph is BASIC (the Protein build builds no canonical or primary graph), and a sequence has
  no reverse complement. Nothing reverse-complements on it: no reverse orientation, no palindrome, no mirror rows.
- k = 15 on UniParc: a k-mer is 15 residues; every offset, length and `*_bp` field counts residues.

## 2. A protein pattern

The request names it as before: `patterns[i].protein`. One position per character, read one residue per position
(no codon automaton). Case-insensitive; the answer's `pattern` is the text in upper case, brackets kept.

| in the text | the position's set |
|---|---|
| one of the 22 residues (20 amino acids, U, O) | that residue |
| `X` | any of the 22 residues |
| `B` | D or N |
| `Z` | E or Q |
| `J` | I or L |
| `[...]` (PROSITE-like class, e.g. `[STA]`) | the union of the letters inside, each read as above (a letter of the 22, or X, B, Z, J) |
| `[^...]` (negated class, e.g. `[^P]`) | the 22 residues less that union |

- **A stored ambiguity letter is matched by nothing.** No position's set holds B, J, Z or X as stored letters: a
  k-mer that holds one where the pattern has a position does not instantiate it, whatever the position (X
  included), and no substitution produces one (section 4). The free flanks of `any_offset` (the k-mer's residues
  before and after the pattern) are not pattern positions and admit every symbol but `$`, the stored ambiguity
  letters included, as a DNA5 graph's flanks admit N (SPEC §3).
- L, the pattern's length, is its number of positions (`[STA]` is one). A pattern of L <= k is searched inside one
  k-mer (`suffix`, `any_offset`); a longer one is a pattern longer than k (SPEC §7.7, §12.1, §20): its anchors, its
  paths (`long_search: "paths"`) and its supported paths (`"supported_paths"`), each on the one strand.
- **Refused in the pattern's slot with `bad_alphabet`** (the other patterns answered): the `dna` and `iupac` kinds
  (a protein graph holds no nucleotides); the stop `*`; any other character (digits, `-`, `.`, whitespace, `(`,
  `)`, `{`, `}`, `,`, `<`, `>`); an empty text; a class that is empty (`[]`, `[^]`), unclosed, nested, or whose
  negation is empty (`[^X]`); a `]` without its `[`; a `^` anywhere but right after `[`. The message names the first
  offending character and its 0-based position in the text.
- **Not in the language** (the plan row's "not in the minimum"): variable gaps (`x(2,4)`), repeats, anchors (`<`,
  `>`), mismatches in patterns longer than k. Their characters are refused as above.

## 3. The request on a protein graph

- `genetic_code`: named at all (any value), 400 `invalid_request` — a protein graph is read in residues, no genetic
  code applies. Checked where the field is read today (after `require_support` and `supported_paths_level`).
- `strands`: `"both"` (the default) and `"forward"` search the one strand; `"reverse"` is 400 `invalid_request`.
  The entry's `strands` is `["+"]`, `by_strand` has the key `+` alone; `palindromic` is false.
- `predicate_strands`: accepted with both values and evaluated as `"context"` (a context's own row: there is no
  reverse complement whose row "either" could read). Every answer with a predicate states the effective value:
  `predicate.strands` and `limits.predicate_strands` are `"context"`, and each selected context's
  `selection_strands` lists `"context"` per label.
- Every other field as on a DNA graph: modes, scopes, caps, `long_search`, labels (`output.labels`), predicates
  (SPEC §19, §20.9, §25), `max_mismatches`.

## 4. Matching, mismatches, bits

- **Matching.** The engine's range narrowing (DESIGN-pattern-search.md §4.1) is alphabet-generic over BOSS codes: a
  position's BOSS codes are those of its residues. Nothing else changes in discovery, the scans, the release, the
  extension (each node's last residue tested against the position's set) or the unmasked graph's check.
- **Mismatches** (SPEC §21) for patterns of L <= k, 0 to 3 per instance, the cap `--pattern-max-mismatches-protein`
  (default 3; at most `--pattern-max-mismatches`): a substitution at a position is any of the 22 residues outside
  its set (21 for a residue, 20 for B, Z, J, none for X). A stored ambiguity letter never substitutes. A pattern
  longer than k is answered in its slot with `mismatches_above_cap` (cap 0), as on DNA. A released context's
  `substitutions` use the DNA schema: `position` (0-based), `expected` (the position as written: a letter, or the
  class such as `[STA]`), `found` (the residue read).
- **Information.** A position of set size s carries log2(22 / s) bits (a residue 4.46, B, Z, J 3.46, X 0); with
  m substitutions the ball over the 22 residues, as SPEC §21.5 over 4 bases. The floor (`min_information_bits`,
  24 by default: about six residues) gates as on DNA. The low-complexity note (`low_complexity_pattern`) is not
  computed: sdust is a nucleotide measure.
- **Counts.** As on DNA: exact with the dummy-edge mask (`metagraph transform --mask-dummy` builds one for a
  protein graph too), upper bounds with an estimate without it (DESIGN-pattern-search.md §4.4). UniParc has no mask
  today: its counts are upper bounds, its lists exact.

## 5. The answer

As SPEC §8 states it, with: `kind` `"protein"`; `length` and `residues` both L (residues); `genetic_code` null; the
contexts' `kmer` and offsets in residues; the substitutions of section 4.

## 6. Capabilities

`GET /pattern/capabilities` (and its envelope on a multi-graph server) and `GET /capabilities`' `pattern` block on
the Protein build:

| field | Protein build | DNA build (unchanged) |
|---|---|---|
| `alphabet` (per graph) | `$ABCDEFGHIJKLMNOPQRSTUVWYZX` | `$ACGT` |
| `kinds` | `["protein"]` | `["dna", "iupac", "protein"]` |
| `protein_residues` | the 22 residues, then X B Z J | the 20, X B Z J, `*` |
| `protein_reading` (new, Protein build only) | `"residues"` | absent (codons, SPEC §12.2) |
| `protein_classes` (new, Protein build only) | `["[...]", "[^...]"]` | absent |
| `genetic_codes`, `default_genetic_code` | `[]`, null | the tables, 1 |
| `protein_rule` | `SPEC-pattern-search.md section 29` | `... section 12.2` |
| `strands`, `default_strands` | `["both", "forward"]`, `"both"` | three, `"both"` |
| `scopes_by_graph_mode` | `basic` only | basic, canonical, primary |
| `mismatches.max`, `.unit` | `{"protein": cap}`, `{"protein": "residue"}` | three kinds each |
| `mismatches.protein_substitutes` | `"residues"` | `"sense_codons"` |

The gate block of `GET /traverse/capabilities` (SPEC §23.2) is the projection onto its listed keys, unchanged: it
carries `kinds`, `protein_residues`, `genetic_codes` and `default_genetic_code`, which tell a client of that block
alone the residue reading (`kinds` without `dna`, `genetic_codes` empty). `alphabet` (per graph) names the unit.
`graph_summary` (multi-graph `GET /capabilities`) gains no field here: E17 adds the alphabet per pair for every
build.

## 7. What the DNA build keeps byte for byte

Every change in `pattern_search.{hpp,cpp}` and `pattern.cpp` is either inside `#if _PROTEIN_GRAPH` or
an identity in the DNA build: named constants that equal the old literals there (the 4 bases, their mask, the
estimate's 4 strings per position), a helper that returns the old table entry. The proof is the DNA build's own:
the DNA unit and integration pattern suites; `scripts/traversal/pattern_fixtures.py --check` with the new DNA binary
(every frozen fixture body of `pattern/`, answered on the mini index, timings blanked);
`TestPatternRegression` against the unmodified base binary; and the contract samples (the pattern samples 11, 12 and
16 regenerate byte for byte from their recorded inputs; this branch changes no sample).

## 8. Choices this note makes (the engine note and the plan row leave them open)

1. X is any of the 22 residues (U and O included, being literals); B, Z, J the usual pairs. Bits are log2(22 / s)
   so that X carries 0.
2. A class holds letters of the language (the 22 and X, B, Z, J, each as its set); `^` negates within the 22.
3. `genetic_code` named on a protein graph: 400 `invalid_request` (the plan row: "a 400"); `strands: "reverse"`:
   400 `invalid_request`; `predicate_strands` evaluated as `"context"`.
4. The entry's `genetic_code` is null (the key kept, so that a protein entry has the same keys on both builds).
5. Substitutions in the DNA schema, `expected` the position's text.
6. The capabilities' additions are the Protein build's only: the DNA build's documents, the contract samples among
   them, are unchanged.
7. The protein fixtures sit in `api/python/tests/data/traverse/pattern_protein/`, beside `pattern/` (whose
   directory list is the DNA generator's, `scripts/traversal/pattern_fixtures.py`, and checked as a whole).

## 9. Tests

- Unit (Protein build, `metagraph/build-protein`): `tests/graph/alignment/test_pattern_protein.cpp` (the parser,
  sets, bits, the engine on small protein graphs against a brute-force count, mismatches, long patterns, the
  unmasked counts) and `tests/cli/test_pattern_protein_route.cpp` (the route: kinds, refusals, the entry, the
  capabilities, labels and a predicate). The DNA-only pattern suites are compiled out under `_PROTEIN_GRAPH`.
- Integration: `integration_tests/test_pattern_protein.py`, skipped in DNA mode: the plan row's scenarios on a
  small protein index built in the test at k = 15 (an exact peptide of at most 15 residues, one with X, one with B,
  Z or J, a bracket and a negated class, mismatches 1 to 3 on a 15-mer, an exact peptide longer than 15, a label
  predicate, the refusals), the supported paths with and without a predicate, through `metagraph pattern` and the
  server; and the fixtures `api/python/tests/data/traverse/pattern_protein/` (a capabilities document and one
  answer), pinned byte for byte apart from timings (`TestPatternProtein.test_fixtures`).
- CI: one Protein build in CI, E15's job `Linux-Protein` of `.github/workflows/dev-branch.yml`: its unit filter
  already holds `*Pattern*` (PatternProtein, PatternProteinRoute, the predicate and mask suites a Protein build
  compiles); its integration filters gain `*(test_pattern_protein.*` (section 10). No job of its own.

## 10. Landing after E15 (protein traversal lands first)

E15 asserts the old refusal in places that flip with E16 (its review, 2026-11-07). At the rebase onto E15:

1. E15's commit "Leave pattern search out of a Protein build's unit tests; test its refusal" is in this branch as
   a cherry-pick (same patch): the rebase drops it; this branch's later commit removes
   `PatternRoute.ProteinGraphIsRefused` and rewords the guards of the DNA suites (they stay DNA-only: they send
   `dna` and `iupac` patterns, which a protein graph refuses).
2. `integration_tests/test_traverse.py` (E15's lines marked "until E16 (pattern on protein graphs)": two in
   `TestTraverseAPI`, one in `test_capabilities_answer_in_the_envelope`): flip to /pattern served on a protein
   graph (`available` true, `kinds` `["protein"]`).
3. `tests/cli/test_task_resources.cpp` (`PatternStatesCpuPerRequestAndPerEntry`, compiled out by E15 in a Protein
   build): drop E15's guard; this branch sends the build's kind (`protein` on a Protein build).
4. `tests/cli/test_server_multigraph.cpp` (`LoadSharesTheResidentRecordMapping`, the /pattern block E15 guards):
   a conflict; keep this branch's version (the build's kind), drop the guard. `RequestLevelChecksNeedNoGraph`:
   this branch's `#if` (genetic_code is 400 invalid_request on a Protein build).
5. `docs/SPEC-labeled-traversal-core.md` §6.12 ("Not offered": /pattern on a protein graph) and its T66 row: offered
   now, SPEC-pattern-search.md §29.
6. `.github/workflows/dev-branch.yml`: E15's comment on the Protein job (the pattern search left out) reworded; the
   job's `integration_test_filters` gain `*(test_pattern_protein.*`.
7. Then the Protein and DNA builds of every target, the suites of section 9 and the GCC 13 check again.
