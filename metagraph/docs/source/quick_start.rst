.. _quick_start:

Quick start
===========

The MetaGraph indexes consist of two parts: k-mer index and graph annotation.

The k-mer index stores all k-mers from the input sequences and represents a de Bruijn graph. This index essentially serves as a dictionary mapping k-mers to their unique positive identifiers.

The graph annotation is a matrix encoding the relation between the k-mers and their attributes (e.g., k-mer ``i`` appears in SRA sample ``j``).

Construct index
---------------

The indexing workflow in MetaGraph consists of two major steps: graph and annotation construction.

Build graph
^^^^^^^^^^^

Simple example::

    metagraph build -v -p 4 -k 31 -o graph transcripts_1000.fa

where ``transcripts_1000.fa`` is a fasta/fastq file with input sequences.

It is possible to pass multiple input files::

    metagraph build -v -p 4 -k 31 -o graph file_1.fa file_2.fa ...

or pipe a list of input files::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph

which will construct a single de Bruijn graph, a k-mer index for all k-mers extracted from the input sequences. 

To see the list of all available flags, type ``metagraph build``.

There are various graph representations available in MetaGraph, which can be chosen with flag ``--graph``.
However, the default ``succinct`` representation is usually the best choice because of its great scalability (requires only 2-4 bits per k-mer).

Construct graph with disk swap
""""""""""""""""""""""""""""""

For very large inputs, graphs can be constructed with disk swap to limit the RAM usage (available only for succinct graphs)::

    metagraph build -v -k 31 -o graph --parallel 36 --disk-swap /tmp --disk-cap-gb 100 file.fa ...

Transform graph to other representations
""""""""""""""""""""""""""""""""""""""""

To transform a ``succinct`` graph to a more compressed and smaller representation, run::

    metagraph transform -v --state small -o graph_small graph.dbg

To transform a graph back to sequences, it can be traversed to extract all its contigs/unitigs::

    metagraph transform -v --to-fasta -o contigs -p 4 graph.dbg

These sequences contain exactly all k-mers indexed in the graph and can be used as their non-redundant representation.

Build graph in canonical mode
"""""""""""""""""""""""""""""

If the input sequences are raw reads where the direction (strand) is unknown, it makes sense to index along with each read its reverse complement sequence.

MetaGraph has a special graph mode where each k-mer indexed in the graph automatically adds its reverse complement k-mer to the index. To build a canonical graph from a set of reads/sequences, run ::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph --mode canonical

Build graph in primary mode
"""""""""""""""""""""""""""

Canonical graphs contain each k-mer in both of its forms (given and reverse complement), however, the same structure can be modeled by storing explicitly only one of them and implicitly modeling the other.
Often, different tools do this by only storing the lexicographically smallest of the two k-mers. However, this would not be possible to implement with the succinct graph representation.
Hence, we introduce relax this constraint and pick *any* of the two forms of each k-mer.
In particular, the canonical graph is fully traversed and a k-mer is marked as primary if it was reached in the traversal before its reverse complement.
The graph containing only primary k-mers is called *primary*.

The algorithm for primarization of a canonical graph is as follows.

1. First, we extract a set of primary contigs from the canonical graph::

    metagraph transform -v --to-fasta --primary-kmers -o primary_contigs -p 4 graph.dbg

2. Then, we construct a new graph from the primary contigs and indicate that this graph is *primary*::

    metagraph build -v -p 4 \
                    -k 31 \
                    -o graph_primary \
                    --mode primary \
                    primary_contigs.fasta.gz

Now this new graph ``graph_primary.dbg`` emulates the same canonical graph, while taking only half of its space.

Annotate graph
^^^^^^^^^^^^^^

Once a graph is constructed, there are multiple ways to annotate it to encode metadata.

Annotate sequence headers
"""""""""""""""""""""""""

For annotating each sequence with its header in the fasta/fastq file, run ::

    metagraph annotate -i graph.dbg --anno-header -o annotation transcripts_1000.fa

This is a common annotation scenario when indexing reference sequences or assembled genomes.

Annotate source filename
""""""""""""""""""""""""

For labeling all k-mers from a file with a single id, the commonly used command is ::

    metagraph annotate -v -i graph.dbg \
                       --anno-filename \
                       -o annotation \
                       file_1.fa file_2.fa ...

which will construct annotation labeling k-mers from the first file by label ``file_1.fa``, k-mers from the second file by label ``file_2.fa``, etc.

In this mode, usually it is preferred to independently construct a single annotation column for each input file, where each can be constructed in parallel by adding ``--separately -p <num_threads>``::

    metagraph annotate -v -i graph.dbg \
                       --anno-filename \
                       --separately -p 36 \
                       -o annotation \
                       file_1.fa file_2.fa ...

**Note** that it is recommended to run annotation from a set of long (primary) contigs/unitigs, where all k-mers have already been deduplicated, especially when annotating a (primary) graph in the ``succinct`` representation. This can be achieved by simply constructing individual (canonical) de Bruijn graphs from all read sets and transforming them to contigs. These contigs serve as an equivalent non-redundant representation of sets of k-mers from each read set and using them for annotation usually speeds up the process by one or two orders of magnitude.

Annotate graph with custom labels
"""""""""""""""""""""""""""""""""

To add a custom annotation label for all k-mers from an input file, add ``--anno-label <LABEL_NAME>`` when annotating the graph.


Transform annotation
^^^^^^^^^^^^^^^^^^^^

To enhance the query performance and reduce the memory footprint, annotations can be converted to other representations.

There are different annotation representations available in MetaGraph.
For instance, ``Rainbowfish`` can be used with relatively small instances to achieve a very fast query speed. In contrast, ``RowDiff<Multi-BRWT>`` typically achieves the smallest memory footprint while still providing a good query performance.

Convert annotation to Rainbowfish
"""""""""""""""""""""""""""""""""

The conversion to Rainbowfish consists of two steps.

1. First, convert the column-compressed annotation to the row-major representation::

    find . -name "*.column.annodbg" | metagraph transform_anno -v \
                                                 --anno-type row \
                                                 -o annotation ...

2. Then, transform the row-major annotation to the compressed Rainbowfish representation::

    metagraph transform_anno -v --anno-type rbfish \
                                -o annotation \
                                annotation.row.annodbg


Convert annotation to RowDiff<Multi-BRWT>
"""""""""""""""""""""""""""""""""""""""""

The conversion to ``RowDiff<Multi-BRWT>`` is done in two steps.

1. Transform annotation columns ``*.column.annodbg`` to ``row_diff`` in three stages::

    metagraph transform_anno -v --anno-type row_diff --row-diff-stage 0 ...
    metagraph transform_anno -v --anno-type row_diff --row-diff-stage 1 ...
    metagraph transform_anno -v --anno-type row_diff --row-diff-stage 2 ...

2. Transform the RowDiff-sparsified columns ``*.row_diff.annodbg`` to ``Multi-BRWT``::

    metagraph transform_anno -v --anno-type row_diff_brwt --greedy --fast ...
    metagraph relax_brwt -v -p 18 \
                         --relax-arity 32 \
                         -o annotation_relaxed \
                         annotation.row_diff_brwt.annodbg

Check stats
^^^^^^^^^^^

The stats for a constructed graph/annotation can always be checked with ::

    metagraph stats graph.dbg
    metagraph stats -a annotation.column.annodbg

Query index
-----------

Using Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To query a MetaGraph index (graph + annotation) using the command line interface (CLI), run ``metagraph query``, e.g.::

    metagraph query -i graph.dbg \
                    -a annotation.column.annodbg \
                    --count-kmers \
                    --discovery-fraction 0.1 \
                    transcripts_1000.fa

For alignment, see ``metagraph align``.

To load up a MetaGraph index in the server mode for querying it with Python API or HTTP requests, run::

    metagraph server_query -i graph.dbg \
                           -a annotation.column.annodbg \
                           --port <PORT> \
                           --parallel <NUM_THREADS>

Using Python API
^^^^^^^^^^^^^^^^
See :ref:`api`
