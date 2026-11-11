# Parity with the labeled-traversal prototype

The responses of the prototype `metagraph traverse` (branch `gr/labeled-traversal`, tree of
commit `d5d26baa`) to the requests of `scripts/traversal/traverse_parity.py` on the mini RefSeq
index (`scripts/traversal/fetch_mini_refseq.sh`, `mini_refseq-bea44d60.tar.gz`), projected onto
the fields of the v1 response schema (`docs/traverse.md`).

Recorded on 2026-11-01 with a Release build of `d5d26baa` (GCC 16.2.1, Linux x86-64; binary
sha256 `c8f71225fe3d34877687ad3fde25102d2ad41fb94a811f731ebcc50d4164e5a1`):

    traverse_parity.py --metagraph <d5d26baa>/build/metagraph --index build/mini_refseq --record

The prototype was sent each request with the knobs v1 fixes stated explicitly
(`to_prototype()`): `on_reconverge: keep`, `hairpins: skip`, `max_splits_per_path: 64`,
`continuation_bp: 0`, `detail: full`.

v1 picks seeds to cover the labels, where the prototype picked them by support, so the
`/resolve` part of the parity check covers discovery only; the `/traverse` cases compare as
before (2026-11-09).

The `labels` of `resolve_discover_headers.json` are sorted by support, then name (2026-11-09): equally
supported labels are ranked differently by the prototype and v1.

v1 follows a header label along its record's coordinates, where the prototype follows it by
membership (2026-11-11). The eight requests whose header labels reach a record end or a repeat
were removed with their golden files. `test_mini_refseq` checks header walks against the
source records.
