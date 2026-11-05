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
