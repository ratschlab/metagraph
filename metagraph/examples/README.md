# How to Use MetaGraph as a C++ Library

Examples demonstrating how to integrate MetaGraph into your C++ project using CMake.

## Integration Methods

Two approaches for integrating MetaGraph:

### 1. [add_subdirectory](subdirectory/)

Simple and direct approach using `add_subdirectory`.

- Best for: Active development, simpler setup
- See: [subdirectory/README.md](subdirectory/README.md)

### 2. [ExternalProject_Add](ExternalProject/)

Automatic download and isolated build.

- Best for: Stricter isolation
- See: [ExternalProject/README.md](ExternalProject/README.md)

## Example Code

[`src/basic_query.cpp`](src/basic_query.cpp) demonstrates:

- Loading a de Bruijn graph and annotations
- Querying sequences against the annotated graph
- Retrieving matching labels

## Test Data

Pre-built graphs and query sequences in [`data/`](data/):

- DNA and Protein graphs with annotations
- Sample query sequences

Both examples include a `make run_example` target for quick testing.

## Coordinate context exploration

`src/context_cursor.cpp` is an experimental JSON-lines driver for
`graph::traversal::ContextCursor`. Build the `context_cursor` target in the
`subdirectory` example. It uses the normal query loader, including `--mmap`,
coordinate annotations, row-diff graph attachment, and the adjacent `.seqs` file:

```sh
./context_cursor query --mmap -i graph.dbg -a annotation.column_coord.annodbg requests.jsonl
```

Each input line is one request, for example:

```json
{"seed":"ACGA","radius":1000,"columns":["sample.fa"],"steps":[{"side":"left","base":"C","radius":1000}]}
```

The cursor finds exact, consecutive occurrences of the entire seed within source
records. Closure extends a flank only while every retained occurrence supports
the same next base. Choosing a base filters that same occurrence set for both
flanks. For records `CACGAG` and `TACGAA`, seed `ACGA` has joint choices `(C,G)`
and `(T,A)`: choosing left `C` then resolves the right side to `G`.

Results include the current sequence, the original seed's offset in it, source
column/record/header and zero-based source offsets, and joint flank counts.
The `closure` object reports `left_added`, `right_added`, `left_stop` and
`right_stop` for each automatic extension. The stop reasons are:

| Stop | Meaning |
| --- | --- |
| `no_occurrences` | The complete seed has no indexed occurrence within a source record. |
| `radius_limit` | This call added the allowed number of bases on that flank. It does not establish that the source ends there. |
| `record_end` | At least one retained occurrence reaches its record boundary. |
| `unindexed` | At least one retained occurrence has an unindexed continuation. |
| `branch` | The retained occurrences have different indexed next bases. |

An empty occurrence set takes precedence over the radius; the radius takes
precedence over inspecting another base. Otherwise record ends take precedence
over unindexed continuations and then differing bases. The joint choices retain
the full breakdown, including mixtures of boundaries and continuing bases.
An empty flank character denotes a record end; `?` denotes an unindexed source
continuation (for example an ambiguous base). Neither is silently discarded.
A failed explicit choice leaves the cursor unchanged and reports `extended:false`
and `closure:null`. A successful choice reports `extended:true`; its closure
counts exclude the one explicitly chosen base. C++ callers can copy a
cursor to explore alternative branches without decoding the seed again. `columns` selects
annotation columns, not individual FASTA headers; omitting it selects all columns.
An empty `columns` array also selects all columns. Results keep separate records
even when their FASTA headers are identical; use column and record ID to identify them.

This first implementation supports immutable BASIC DNA graphs with `k >= 2`.
It requires both coordinates and record boundaries. Canonical/PRIMARY graphs,
reverse-strand traversal, approximate matching, label switching, and quorum
policies are not implemented. Ordinary k-mer match counts retain their existing
semantics; the cursor's whole-record occurrence requirement is stricter.

The example limits seeds to 10,000 bases, closure to 10,000 bases per side per
step, and accepted occurrences to 100,000 (`max_occurrences` may lower this to
any positive integer). Exceeding the occurrence limit fails the request; it never
silently samples occurrences. `radius` defaults to 1,000 and may be zero. Each
step inherits the request radius unless it specifies its own. At most 100 explicit
steps are accepted. The radius applies separately to each closure, not cumulatively
to the entire request.
These are result limits, not a hard decoding-memory or execution-time budget.
After closure the example separately inspects the joint choices, including at
a radius-limited boundary; the radius does not limit that inspection work.
It is an offline example, not a server endpoint. The core uses a bounded local
`NodeFirstCache`, decodes distinct seed k-mers once, and reuses an unchanged
opposite flank's boundary results. It does not build a whole-graph cache.

Unknown fields, invalid steps, duplicate JSON keys and trailing JSON content are
rejected. The complete request is validated before seed decoding. Each nonempty
input line produces one result or an `error` string; later lines are still processed
after an error, and the process exits nonzero if any line failed.

To run the native example tests, build `run_context_example` in a DNA
subdirectory-example build. It requires Python 3 with only the standard library,
creates temporary coordinate indexes, and checks successful and invalid requests
with heap/mmap loading. These checks also run as part of `run_example`, including
in DNA example CI jobs.
