# Pattern search on a protein graph: fixtures (SPEC-pattern-search.md §29)

What the Protein build (`metagraph_Protein`) answers on the protein index of
`integration_tests/test_pattern_protein.py` (k 15, five records `p1` to `p5` named by their headers, a column
annotation, the graph with its dummy-edge mask), served as `server_query --index-release fixtures-protein`. The
documents are written without their timings (`timing`, `elapsed_ms`); `TestPatternProtein.test_fixtures` compares
them byte for byte with what the build answers now (`METAGRAPH_WRITE_PROTEIN_FIXTURES=1` rewrites them).

They sit beside `../pattern/`, whose fixtures the DNA build answers on the mini index
(`scripts/traversal/pattern_fixtures.py`, which owns that directory's list).

- `capabilities` GET /pattern/capabilities: `alphabet`, `kinds` `["protein"]`, `protein_reading` `"residues"`,
  `protein_classes`, the 26 `protein_residues`, no genetic codes, one strand, the mismatches of the protein kind.
- `entry` POST /pattern with `long_search` `"paths"`: exact peptides of 11 and 15 residues, one with X, one with B
  and J, a class, negated classes, two patterns longer than k (their anchors and paths, one with a class), and two
  slot refusals (`bad_alphabet`: the `dna` kind, the stop `*`).
