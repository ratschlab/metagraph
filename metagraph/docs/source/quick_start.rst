.. _quick_start:

Quick start
===========

MetaGraph constructs indexes composed of two main elements: a k-mer index and an annotation matrix.

The k-mer index stores all k-mers from the input sequences and represents a de Bruijn graph.
This index serves as a dictionary mapping k-mers to their unique positive identifiers.

.. It can also be used to map sub-k-mers (or spaced k-mers) to ranges of their identifiers (see TODO).

The second element is a matrix encoding the relation between the k-mers and their attributes.
These relations may represent, for instance:

* k-mer ``i`` is present in sequence/genome ``j``
* k-mer ``i`` is present in SRA sample ``j``
* k-mer ``i`` is present in file ``j``
* k-mer ``i`` is present in a collection of sequences/files marked by ``j``
* k-mer ``i`` is highly expressed in sample ``j``

.. TODO: Describe counts/coordinate annotation

In the following, you will find simplified instructions and examples for constructing a MetaGraph
index or querying it.

Construct index
---------------

The indexing workflow in MetaGraph consists of two major steps: graph construction and annotation construction.

Construct graph
^^^^^^^^^^^^^^^

Simple example::

    metagraph build -v -p 4 -k 31 -o graph transcripts_1000.fa

where ``transcripts_1000.fa`` is a fasta/fastq file (may be gzipped) with input sequences.

It is possible to pass multiple input files::

    metagraph build -v -p 4 -k 31 -o graph file_1.fa file_2.fa ...

or pipe a list of them::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph

which will construct a single de Bruijn graph (a k-mer index) from all k-mers extracted from the input sequences. 

To see the list of all available flags, type ``metagraph build``.

There are various graph representations available in MetaGraph, which can be chosen with flag ``--graph``.
However, the default ``succinct`` representation is usually the best choice because of its great scalability (requires only 2-4 bits per k-mer) and the ability to search sub-k-mers/k-mer ranges.

To check stats for a constructed graph, type::

    metagraph stats graph.dbg


Construct graph with disk swap
""""""""""""""""""""""""""""""

For very large inputs, graphs can be constructed with disk swap to limit the RAM usage (currently, available only for ``succinct`` graph representations).
For example, to restrict the buffer size to 4 GiB, the following build command can be used::

    metagraph build -v -k 31 -o graph -p 36 --disk-swap /tmp --disk-cap-gb 4 file.fa ...

Using larger buffers usually speeds up the construction. However, using a 50-100 GB buffer is always sufficient, even when constructing graphs with trillions of k-mers.

Transform graph to other representations
""""""""""""""""""""""""""""""""""""""""

To transform a ``succinct`` graph to a more compressed and smaller representation, run::

    metagraph transform -v --state small -p 4 -o graph_small graph.dbg

To transform a graph back to sequences, it can be traversed to extract all its contigs/unitigs::

    metagraph transform -v --to-fasta -o contigs -p 4 graph.dbg

These sequences contain exactly all k-mers indexed in the graph and can be used as their non-redundant (deduplicated) representation.

The assembled contigs are written to a compressed FASTA file, which can be inspected with::

    zless contigs.fasta.gz

Build graph in canonical mode
"""""""""""""""""""""""""""""

When the input sequences are raw reads from an unknown direction (strand), it is natural to index along with each sequence its reverse complement.

MetaGraph has a special graph mode where each k-mer indexed in the graph automatically adds its reverse complement k-mer to the index. To build a canonical graph from a set of reads/sequences, run ::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph --mode canonical

Build graph in primary mode
"""""""""""""""""""""""""""

Canonical graphs contain each k-mer in both of its forms (given and reverse complement), but the same data structure can be modeled by storing only one of them and implicitly modeling the other.
Often, different tools achieve this by only storing the lexicographically smallest of the two k-mers. However, this is not be possible to efficiently implement with the ``succinct`` graph representation.
Hence, we relax this constraint and pick *any* of the two forms of each k-mer.
In a nutshell, this representation is constructed by fully traversing the canonical graph and marking a k-mer as *primary* if it was reached before its reverse complement in the traversal.
The graph containing only primary k-mers is called *primary*.

The algorithm for primarization of a canonical graph is as follows:

1. First, extract a set of primary contigs (stretches of primary k-mers) from the canonical graph. ::

    metagraph transform -v --to-fasta --primary-kmers -o primary_contigs -p 4 graph.dbg

2. Then, construct a new graph from the primary contigs and mark this graph as *primary*. ::

    metagraph build -v -p 4 \
                    -k 31 \
                    -o graph_primary \
                    --mode primary \
                    primary_contigs.fasta.gz

Now, this new graph ``graph_primary.dbg`` emulates the original canonical graph (e.g., when querying or annotating), while taking only half of its space.

.. TODO: note that canonical graphs must not be used with row-diff<*> annotations and always must be primarized


Annotate graph
^^^^^^^^^^^^^^

Once a graph is constructed, there are multiple ways to annotate it to encode metadata.

Annotate sequence headers
"""""""""""""""""""""""""

For annotating each sequence with its header in the fasta/fastq file, run ::

    metagraph annotate -v -i graph.dbg --anno-header -o annotation transcripts_1000.fa

This is a common annotation scenario when indexing reference sequences or assembled genomes.

To check stats for the constructed annotation, type::

    metagraph stats -a annotation.column.annodbg

All annotation labels (column names) for an annotation matrix can be printed with::

    metagraph stats --print-col-names -a annotation.column.annodbg

Annotate source filename
""""""""""""""""""""""""

To label all k-mers from each file with the same id (for instance for the experiment discovery problem), the command is::

    metagraph annotate -v -i graph.dbg --anno-filename -o annotation file_1.fa file_2.fa ...

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
