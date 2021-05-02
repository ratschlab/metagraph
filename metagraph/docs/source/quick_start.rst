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

Simple example
""""""""""""""

De Bruijn graphs are constructed with ``metagraph build``::

    metagraph build --verbose --parallel 4 -k 31 --outfile-base graph transcripts_1000.fa

where ``transcripts_1000.fa`` is a fasta/fastq file (may be gzipped) with input sequences. The same
command can be shortened to::

    metagraph build -v -p 4 -k 31 -o graph transcripts_1000.fa

It is possible to pass multiple input files::

    metagraph build -v -p 4 -k 31 -o graph file_1.fa file_2.fa ...

or pipe a list of them::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph

All of the above commands will construct a single de Bruijn graph (a k-mer index) from all k-mers extracted from the input sequences. 

To see the list of all available flags, type ``metagraph build``.

There are various graph representations available in MetaGraph, which can be chosen with flag ``--graph``.
However, the default ``succinct`` representation is usually the best choice because of its great scalability (requires only 2-4 bits per k-mer) and the ability to search sub-k-mers/k-mer ranges.

To check stats for a constructed graph, type::

    metagraph stats graph.dbg

Construct from KMC counters
"""""""""""""""""""""""""""

Apart from the standard FASTA/FASTQ/VCF input formats (uncompressed or gzipped), a graph can be
constructed from k-mer counts produced by the `KMC <https://github.com/refresh-bio/KMC>`_ tool.

KMC is extremely efficient in counting k-mers and can be used to quickly pre-process the
input and deduplicate/count/filter the input k-mers.
For example, the following command can be used to construct a graph only from k-mers
occurring at least 5 times in the input::

    K=31
    ./KMC/kmc -ci5 -t4 -k$K -m5 -fm SRR403017.fasta.gz SRR403017.cutoff_5 ./KMC
    metagraph build -v -p 4 -k $K --mem-cap-gb 10 -o graph SRR403017.cutoff_5.kmc_pre

.. note:: In the above example, we use ``./KMC`` as the directory where KMC will store its
          intermediate caches. Depending on your input data, this directory should be at a location
          with a sufficient amount of free intermediate storage space.

Construct with disk swap
""""""""""""""""""""""""

For very large inputs, graphs can be constructed with disk swap to limit the RAM usage (currently, available only for ``succinct`` graph representations).
For example, to restrict the buffer size to 4 GiB, the following build command can be used::

    metagraph build -v -k 31 -o graph -p 36 --disk-swap /tmp --disk-cap-gb 4 file.fa ...

Using larger buffers usually speeds up the construction. However, using a 50-100 GB buffer is always sufficient, even when constructing graphs with trillions of k-mers.

Transform to other representations
""""""""""""""""""""""""""""""""""

To transform a ``succinct`` graph to a more compressed and smaller representation, run::

    metagraph transform -v --state small -p 4 -o graph_small graph.dbg

To transform a graph back to sequences, it can be traversed to extract all its contigs/unitigs::

    metagraph transform -v --to-fasta -o contigs -p 4 graph.dbg

These sequences contain exactly all k-mers indexed in the graph and can be used as their non-redundant (deduplicated) representation.

The assembled contigs are written to a compressed FASTA file, which can be inspected with::

    zless contigs.fasta.gz

Build graph in canonical mode
"""""""""""""""""""""""""""""

When the input sequences are raw reads of unknown directionality (strandedness), it is natural to index along with each sequence its reverse complement.

MetaGraph has a special graph mode where each k-mer indexed in the graph automatically adds its reverse complement k-mer to the index. To build a canonical graph from a set of reads/sequences, add ``--mode canonical`` to the build command::

    find . -name "*.fa" | metagraph build -v -p 4 -k 31 -o graph --mode canonical

Build graph in primary mode
"""""""""""""""""""""""""""

Canonical graphs contain each k-mer in both of its forms (forward and reverse complement), but the same data structure can be modeled by storing only one of them and implicitly modeling the other.
Often, different tools achieve this by only storing the lexicographically smallest of the two
k-mers. However, it is not possible to efficiently implement this with the ``succinct`` graph representation.
Hence, we relax this constraint and pick *any* of the two forms of each k-mer.
In a nutshell, this representation is constructed by fully traversing the canonical graph and marking a k-mer as *primary* if it was reached before its reverse complement in the traversal.
The graph containing only primary k-mers is called a *primary* graph.

The algorithm for primarization of a canonical graph is as follows:

1. First, extract a set of primary contigs (stretches of primary k-mers) from the canonical graph::

    metagraph transform -v --to-fasta --primary-kmers -o primary_contigs -p 4 graph.dbg

2. Then, construct a new graph from the primary contigs and mark this graph as *primary* by adding ``--mode primary`` to the build command::

    metagraph build -v -p 4 \
                    -k 31 \
                    -o graph_primary \
                    --mode primary \
                    primary_contigs.fasta.gz

Now, this new graph ``graph_primary.dbg`` emulates the original canonical graph (e.g., when querying
or annotating) containing the same information as the original canonical graph, while taking only
half of the space.

.. TODO: note that canonical graphs must not be used with row-diff<*> annotations and always must be primarized

Graph cleaning
""""""""""""""

For removing sequencing noise, there are graph cleaning and k-mer
filtering procedures implemented in MetaGraph. These are based on the assumption that
k-mers with a relatively low abundance in the input data are likely due to sequencing errors, and
hence should be dropped to keep the k-mer index free of the non-existent k-mers.

::

    K=31
    metagraph build -v -p 4 -k $K --count-kmers -o graph SRR403017.fasta.gz

    metagraph clean -v -p 4 --to-fasta --prune-tips $((2*$K)) --prune-unitigs 0 --fallback 2 \
                    -o SRR403017_clean_contigs graph.dbg

    zless SRR403017_clean_contigs.fasta.gz



Annotate graph
^^^^^^^^^^^^^^

Once a graph is constructed, there are multiple ways to construct the corresponding annotation to
encode its metadata.

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

which will annotate k-mers from the first file by label ``file_1.fa``, k-mers from the second file by label ``file_2.fa``, etc.

Annotate with disk swap
***********************
When the input files and the output annotation are very large, disk swap space can be used
by setting flags ``--disk-swap`` and ``--mem-cap-gb``, to limit the size of internal buffers
and reduce RAM usage during annotation construction::

    metagraph annotate -v -i graph.dbg --anno-filename --disk-swap /tmp --mem-cap-gb 1 \
                          -o annotation file_1.fa file_2.fa ...

Annotate files independently
****************************
It is recommended to independently construct a single annotation column per each input file.
To do this in parallel and avoid loading the same graph multiple times, run one annotation
command with flags ``--separately -p <num_threads>`` added::

    metagraph annotate -v -i graph.dbg --anno-filename --separately -p 36 \
                          -o annotation file_1.fa file_2.fa ...

This will create a new directory ``annotation/`` with individual annotation columns::

    file_1.fa.column.annodbg    file_2.fa.column.annodbg    ...

.. important:: It is recommended to run annotation from a set of long (primary) contigs/unitigs,
    where all k-mers have already been deduplicated, especially when annotating a (primary) graph
    in the ``succinct`` representation. In contrast, annotating a ``succinct`` graph from
    separate k-mers (especially not deduplicated) will take orders of magnitude longer.
    The contigs serve as an equivalent non-redundant representation of the k-mers sets and, thus,
    result in the same graph annotation.
    **Thus, in practice,** for large inputs, it is recommended to construct
    individual (canonical) de Bruijn graphs from all read sets, called sample graphs, and
    transform them to contigs. These contig sets are then used instead of the original read
    sets to construct and annotate the join (primary) graph.

Annotate graph with custom labels
"""""""""""""""""""""""""""""""""

To add a custom annotation label for all k-mers from an input file, add ``--anno-label <LABEL_NAME>`` when annotating the graph.


Transform annotation
^^^^^^^^^^^^^^^^^^^^

To enhance the query performance and reduce the memory footprint, annotations can be converted to other representations.

There are several different annotation representations available in MetaGraph (see the possible values for flag ``--anno-type`` in ``metagraph transform_anno``).
For instance, ``Rainbowfish`` can be used to achieve a very fast query speed, but it can
be applied only to relatively small problem instances (about 100 GB) because of the limited
compression performance and the complexity of the construction algorithm.
In contrast, ``RowDiff<Multi-BRWT>`` typically achieves
the best compression while still providing a good query performance, and thus, it is
recommended for very large problem instances.

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

To load up a MetaGraph index in the server mode for querying it with the Python API or via HTTP requests, run::

    metagraph server_query -i graph.dbg \
                           -a annotation.column.annodbg \
                           --port <PORT> \
                           --parallel <NUM_THREADS>

Using Python API
^^^^^^^^^^^^^^^^
See :ref:`api`
