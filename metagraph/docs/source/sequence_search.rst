.. _sequence_search:

Sequence search
===============

MetaGraph allows for query sequences to be searched against the graph alone, returning
the closest path in the graph, or against the annotated graph, returning a set of associated
labels. These are referred to as the *align* and *query* regimes, respectively.

Align sequences to the graph
-----------------------------

Sequence alignment features can be accessed via :code:`metagraph align`.
By default, a sequence of graph paths is output which form a disjoint cover of the
query sequence. Depending on the desired level of sensitivity, alignment options range
from simply finding exact k-mer matches to performing an alignment to find a
best-scoring path in the graph.

Exact k-mer matching
^^^^^^^^^^^^^^^^^^^^
Also referred to as pseudo-alignment, this feature is accessed via the additional :code:`--map` flag.
This mode extracts the sequence of k-mers from a query sequence and reports the indices
of the corresponding nodes in the graph. An example command may be::

    metagraph align --map -i MYGRAPH.dbg MYREADS.fa

Input sequences may be in FASTA or FASTQ format (both uncompressed and gzipped).
The output is in TSV format with the first column being an input k-mer and the second
column being the corresponding node in the graph (or :code:`0` if not present).

For less verbose output, the additional :code:`--query-presence` and :code:`--count-kmers`
flags are available.

- :code:`--query-presence` outputs one line per sequence indicating whether the sequence is present (:code:`1`) or absent (:code:`0`). A sequence is considered to be present if its fraction of present k-mers is at least :code:`d`, as set by the :code:`--discovery-fraction` flag.
- :code:`--count-kmers` outputs one line per sequence in TSV format. The first column is the sequence header, while the second column is of the form :code:`a/b/c`, where :code:`a` is the number of matching k-mers, :code:`b` is the total number of k-mers, and :code:`c` is the total number of unique matching k-mers (where reverse complements are considered to be matching).

Sequence-to-graph alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Additional Parameters
^^^^^^^^^^^^^^^^^^^^^

Query sequences against the index
---------------------------------
(Experiment discovery)

Parameters for exact k-mer matching
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameters for alignment
^^^^^^^^^^^^^^^^^^^^^^^^

