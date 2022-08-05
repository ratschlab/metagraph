.. _api:

Python API
==========

The MetaGraph API provides a simple way to query indexes (running on a remote server or locally) in Python and supports both exact k-mer
matching as well as inexact search (alignment).

.. _install api:

Installation
------------

Install MetaGraph API in Python::

    pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/api/python"

Sequence search
---------------
The Python client has a ``search`` method for querying an index running on a server.
The method accepts a single sequence or a list of sequences represented with strings.

.. py:function:: metagraph.client.GraphClient.search(self, ...)

        :param      sequence:             Query sequence
        :type       sequence:             Union[str, Iterable[str]]
        :param      top_labels:           The maximum number of matched labels to retrieve [default: 100]
        :type       top_labels:           int
        :param      discovery_threshold:  The minimum fraction (between 0.0 and 1.0) of k-mers from the query required to match a label (occur in a sample) in order for that label to show up in the result [default: 0.0]
        :type       discovery_threshold:  float
        :param      with_signature:       Return the signature of k-mer matches
        :type       with_signature:       bool
        :param      abundance_sum:        Compute the sum of abundances for all k-mers matched
        :type       abundance_sum:        bool
        :param      query_counts:         Query k-mer counts
        :type       query_counts:         bool
        :param      query_coords:         Query k-mer coordinates
        :type       query_coords:         bool
        :param      align:                Align the query sequence to the joint graph and query labels for that alignment instead of the original sequence
        :type       align:                bool
        :param      align_params:         The parameters for alignment (see method align())
        :type       align_params:         dictionary
        :return:                          A data frame with query results
        :rtype:                           pandas.DataFrame

``GraphClient`` will return an instance of ``pandas.DataFrame`` storing the result.
You can also use the lower level ``GraphClientJson`` to instead receive a JSON response.


Search with alignment
^^^^^^^^^^^^^^^^^^^^^
It is possible to first align a sequence to the joint graph and use the aligned sequence to query the index.
This can be done by setting ``align=True`` (default False).
If the ``align`` flag is set, all the alignment options (explained in the :ref:`Sequence alignment <alignment>` section below) are accepted::

    metasub.search(query, discovery_threshold=0.0, top_labels=200,
                   align=True, min_exact_match=0.8, max_num_nodes_per_seq_char=10.0)

Querying k-mer abundance and coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MetaGraph supports k-mer count- and coordinate-aware indexes (see :ref:`Indexing k-mer counts <indexing counts>`
and :ref:`Indexing k-mer coordinates <indexing coordinates>`).
If you are querying such an index,
you can also set ``abundance_sum=True``, ``query_counts=True``, or ``query_coords=True``. If you try to query an
index which does not support either of these query types, the server will return an error.


.. _alignment:

Sequence alignment
------------------
The ``align`` method allows alignment of sequences to the graph.
The method accepts a single sequence or a list of sequences represented with strings.
Additionally, the method accepts the following keyword arguments:

.. py:function:: metagraph.client.GraphClient.align(self, ...)

        Align sequence(s) to the joint graph

        :param      sequence:                    The query sequence
        :type       sequence:                    Union[str, Iterable[str]]
        :param      min_exact_match:             The minimum fraction (between 0.0 and 1.0) of nucleotides covered by seeds required to align the sequence [default: 0]
        :type       min_exact_match:             float
        :param      max_alternative_alignments:  The number of different alignments to return [default: 1]
        :type       max_alternative_alignments:  int
        :param      max_num_nodes_per_seq_char:  The maximum number of nodes to consider per sequence character during extension [default: 10.0]
        :type       max_num_nodes_per_seq_char:  float

        :returns:   A data frame with alignments
        :rtype:     pandas.DataFrame


Examples
--------

.. _install metasub example:

Example of search in MetaSUB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    from metagraph.client import GraphClient

    metasub = GraphClient('dnaloc.ethz.ch', 80, api_path='/api/metasub19')

    lbls = metasub.column_labels()

    # >ENA|A14565|A14565.1 16S rRNA
    query = 'TCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAAT\
            GTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATA\
            ACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGG\
            GATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAG\
            GATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGG\
            GAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT\
            CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTG\
            ACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGG\
            GTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAG\
            ATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCG\
            TAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCG\
            GTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAA\
            ACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC\
            TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAA\
            GGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATT\
            CGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGAT\
            TGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAA\
            ATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCC\
            GGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCA\
            TCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA\
            CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACT\
            CGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTT\
            CCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTA\
            GCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAAC\
            AAGGTAACCGTAGGGGAAC'

    metasub.search(query, discovery_threshold=0.0, top_labels=200)

    metasub.align(query, min_exact_match=0.8)


Search multiple graphs in parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The API provides ``MultiGraphClient``, which can query multiple graph servers in parallel.
Both ``search`` and ``align`` have the keyword argument ``parallel`` [default: True].
If ``parallel=True``, the result will be a dictionary mapping the specified index names to instances
of ``concurrent.futures.Future``.
If ``parallel=False``, all graphs will simply be queried in sequence and the results will
be instances of ``pandas.DataFrame``.

::

    from metagraph.client import MultiGraphClient

    multi = MultiGraphClient()

    multi.add_graph('dnaloc.ethz.ch', 80, api_path='/api/metasub19', name='metasub')
    multi.add_graph('dnaloc.ethz.ch', 80, api_path='/api/uhgg', name='uhgg')

    multi.list_graphs()
    # {'metasub': ('dnaloc.ethz.ch', 80), 'uhgg': ('dnaloc.ethz.ch', 80)}

    # >ENA|A14565|A14565.1 16S rRNA
    query= 'TCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAAT\
            GTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATA\
            ACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGG\
            GATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAG\
            GATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGG\
            GAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT\
            CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTG\
            ACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGG\
            GTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAG\
            ATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCG\
            TAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCG\
            GTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAA\
            ACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC\
            TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAA\
            GGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATT\
            CGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGAT\
            TGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAA\
            ATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCC\
            GGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCA\
            TCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA\
            CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACT\
            CGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTT\
            CCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTA\
            GCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAAC\
            AAGGTAACCGTAGGGGAAC'

    # Search in parallel
    futures = multi.search(query, discovery_threshold=0.0, top_labels=100)
    # {'metasub': <Future at 0x116dbed10 state=running>,
       'uhgg': <Future at 0x116dad8d0 state=running>}

    # You can either handle the Future instances yourself
    # or block and wait for all of the results
    result = MultiGraphClient.wait_for_result(futures)


Query a locally hosted index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When an index is hosted locally, say on address ``localhost`` and ``5555``, the API client can
connect to it as follows::

    from metagraph.client import GraphClient

    graph_client = GraphClient('127.0.0.1', 5555, api_path='')

Since in this case requests directly go to the MetaGraph engine without being forwarded via an intermediate HTTP server,
the `api_path` flag should be omitted. (Compare this to the :ref:`example above <install metasub example>`).

Before initializing a client and initiating a connection, a search engine (the main MetaGraph app)
must be started to load up an index for query. This can be done, for instance, as follows::

    metagraph server_query -v -i graph.dbg -a annotation.row_diff_brwt.annodbg --port 5555 -p 10


Other examples
^^^^^^^^^^^^^^
Find more examples `here <https://github.com/ratschlab/metagraph_paper_resources/blob/master/notebooks/>`_.
