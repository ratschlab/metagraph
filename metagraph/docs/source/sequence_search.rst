.. _sequence_search:

Sequence search
===============

.. attention:: This page is in development

MetaGraph allows for query sequences to be searched against the graph alone, returning
the closest path in the graph, or against the annotated graph, returning a set of associated
labels. These are referred to as the *align* and *query* regimes, respectively.

Align sequences to the graph
-----------------------------

Sequence alignment features can be accessed via :code:`metagraph align`.
By default, a sequence of graph paths is output which form a disjoint cover of the
query sequence. Depending on the desired level of sensitivity, alignment options range
from simply finding exact k-mer matches to performing an alignment to finding a
best-scoring path in the graph.

Exact k-mer matching
^^^^^^^^^^^^^^^^^^^^
Also referred to as pseudo-alignment, this feature is accessed via the additional :code:`--map` flag.
This mode extracts the sequence of k-mers from a query sequence and reports the indices
of the corresponding nodes in the graph. An example command may be::

    metagraph align --map -i MYGRAPH.dbg MYREADS.fa

Input sequences may be in FASTA or FASTQ format (uncompressed or gzipped).
The output is in TSV format with the first column being an input k-mer and the second
column being the corresponding node in the graph (or :code:`0` if not present).

For less verbose output, the additional :code:`--query-presence` and :code:`--count-kmers`
flags are available.

- :code:`--query-presence` outputs one line per sequence indicating whether the sequence is present (:code:`1`) or absent (:code:`0`). A sequence is considered to be present if its fraction of present k-mers is at least :code:`d`, as set by the :code:`--min-kmers-fraction-label` flag.
- :code:`--count-kmers` outputs one line per sequence in TSV format. The first column is the sequence header, while the second column is of the form :code:`a/b/c`, where :code:`a` is the number of matching k-mers, :code:`b` is the total number of k-mers, and :code:`c` is the total number of unique matching k-mers (where reverse complements are considered to be matching).

Sequence-to-graph alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the :code:`--map` flag is not used, this enables sequence-to-graph alignment, approximately finding the best-matching path in the graph to the query sequence. Alongside this path, this mode also returns an alignment score and a CIGAR string describing the edits needed to transform the spelling of the graph path to the query sequence. An example command may be::

    metagraph align -i MYGRAPH.dbg MYREADS.fa

The output of the query is in TSV format, with one line per query sequence, where the columns are as follows:

1. Query name
2. Query sequence
3. Strand
4. Reference sequence (the spelling of the matched path
5. Alignment score
6. Number of exact matches
7. CIGAR string-like alignment summary
8. Number of nucleotides trimmed from the prefix of the reference sequence
9. Reference label matches (:code:`;` separated, present if the :code:`-a` flag is passed)

An important parameter is the seed length, which can be set with :code:`--align-min-seed-length` and can be shorter than the value of k used to construct the graph.

If an annotator is provided with the :code:`-a` flag, the returned alignments will be label-consistent, meaning that there is at least one label that is shared by all nodes on the path.

Alignment search options
------------------------

``metagraph align`` defaults to independent seed extension. Full-k searches use
the column aligner and its covered-seed filtering, including label-constrained
searches. On succinct and SSHash graphs, shorter seeds also use wavefront
extension and suffix lookup to retain partial-node matches. The wavefront
aligner implements the opt-in modes below.

Use ``--align-chain`` to construct and connect seed chains, and
``--align-two-sided`` to extend both ends of each seed or chain. Both options are
off by default and can be enabled independently. They can improve sensitivity
for long, error-prone reads at a higher runtime cost. ``--align-only-forwards``
controls query orientation independently of the extension direction.
Two-sided extension also consumes unaligned characters inside the first and
last nodes before traversing graph edges.

``--align-post-chain`` joins completed local alignments and can insert graph
jumps (``G`` in the CIGAR). It uses the column aligner and cannot be combined
with seed chaining, two-sided extension, or disabled extension. Queries needing
the wavefront fallback and coordinate-annotated post-chaining are unsupported.
The CLI's JSON format supports only paths without graph jumps. If post-chaining
or disabled anchor connection produces a jump, use TSV output.

Alignment queries must contain ASCII characters. Case-insensitive DNA and
protein builds normalize lowercase letters before alignment.

For the sensitive pipeline with the standard X-drop threshold::

    metagraph align --align-chain --align-two-sided --align-xdrop 27 -i graph.dbg reads.fa

The alignment endpoint accepts the corresponding Boolean fields
``chain_alignments`` and ``two_sided`` (both default to false).
