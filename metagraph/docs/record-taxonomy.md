# Record taxonomy for pattern predicates

An annotation can group sequences into balanced file bins while its coordinate mapping identifies individual records. To select taxa with the existing pattern predicate language, an optional record-taxonomy view translates those coordinates into virtual taxid columns. The stored graph and bin annotation remain unchanged.

The loader must be explicitly enabled with `--record-taxonomy`. Without that flag, taxonomy sidecars are not loaded, no taxonomy-view memory is allocated, and the original column namespace remains the default.

## Prepare complete mappings

Export the exact indexed headers for every graph:

```sh
metagraph stats --print-sequence-headers annotation.seqs > headers.tsv
```

Use the NCBI catalogue for the indexed release. Its TSV fields are taxonomy ID, species name, accession.version, release directories, status and length. The preparation script joins on the complete accession.version, never an unversioned accession or a representative organism:

```sh
python scripts/prepare_record_taxids.py \
  --catalogue RefSeq-release237.catalog.gz \
  --chunk bact_other=bact_other.headers.tsv \
  --chunk bact_pseudomonadota=bact_pseudomonadota.headers.tsv \
  --output prepared-taxonomy
```

List every chunk being deployed. The output directory must not exist. The SQLite join uses bounded cache memory and disk-backed temporary data. Only a `coverage.json` with `complete: true`, all expected chunks and zero missing accessions is deployable. Preserve its catalogue, inventory and sidecar SHA-256 values. Conflicting taxids are errors. Missing records never become taxid zero or an assumed absence.

Each `annotation.taxids.tsv` belongs beside its matching `annotation.seqs`. It contains the version line `METAGRAPH_RECORD_TAXIDS_V1`, followed by column, sequence ID, exact header and positive taxid, tab separated in record order. Loading verifies every record and rejects missing, extra, reordered or mismatched rows. A changed sidecar changes index identity: generate and verify manifests with `scripts/traversal/index_manifest.py --record-taxonomy`, matching the server's flag.

## Query contract

An enabled graph advertises `predicate.namespaces: ["column", "taxid"]`. Clients must check every selected graph. The request opts in:

```json
{"patterns":[{"dna":"AGCTTTTCATTCTGACTGCAACGGGCAATAT"}],"mode":"count","predicate_namespace":"taxid","predicate":{"any":["511145"]}}
```

`column` is the unchanged default. `taxid` requires a predicate and a complete loaded mapping; otherwise the request is refused with `taxonomy_unavailable`. All existing Boolean predicate operators, projections and supported-path rules apply to the taxid view. Taxids match exactly; there is no implicit descendant or merged-taxid expansion. Headers remain accession labels in record placements.

Each source record gets its own interval in the virtual taxid column. A long path can therefore have taxon-wide `label_intersection` support while failing `record_verified` support because its k-mers occur in different records. This distinction is preserved by coordinate translation. Duplicate record names retain separate record IDs and intervals.

The coordinate projection admits its transient/output buffers before allocation and charges remapping/sorting work. A resource stop keeps the existing partial-result and uncertainty rules. The complete taxonomy view itself is resident index metadata, prepared during loading, outside request budgets; measure its startup memory and query cost on a canary before enabling a large index.

Multi-graph answers keep the existing per-graph contract. A predicate is evaluated within one graph, counts may sum the same k-mer from several graphs, and motif verdicts do not establish a global comparison across graphs.
