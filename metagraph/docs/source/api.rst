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
The Python client has a ``search`` method which allows querying an index running on a server.
The method accepts a single sequence or list of sequences represented with strings.

Example of search in MetaSUB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For basic usage, the method accepts the following keyword arguments:

- ``top_labels``: the maximum number of matched labels to retrieve [default: 100]
- ``discovery_fraction``: a value between 0.0 and 1.0 specifying the minimum fraction of
  matching k-mers needed to match a label [default: 0.0]

::

    from metagraph.client import GraphClient

    metasub = GraphClient('dnaloc.ethz.ch', 80, api_path='/api/metasub19')

    lbls = metasub.column_labels()

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

    metasub.search(query, discovery_threshold=0.0, top_labels=200)

``GraphClient`` will return an instance of ``pandas.DataFrame`` storing the result.
You can also use the lower level ``GraphClientJson`` to instead receive a JSON response.

Search with alignment
^^^^^^^^^^^^^^^^^^^^^
It is possible to first align a sequence and use the aligned sequence to query the index.
This can be done by setting ``align=True`` (default False).
If the ``align`` flag is set, all the alignment options (explained in the alignment section below) are accepted::

    metasub.search(query, discovery_threshold=0.0, top_labels=200,
                   align=True, min_exact_match=0.8, max_num_nodes_per_seq_char=10.0)

Searching count and coordinate indexes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Metagraph supports k-mer count-aware and coordinate-aware indexes. If you are querying such an index,
you can also set ``weighted_with_counts=True`` or ``query_coords=True``, respectively. If you try to query an
index which does not support either of these query types, the server will return an error.

Alignment
---------
The ``align`` method allows alignment of sequences to the the graph.
The method accepts a single sequence or list of sequences represented with strings.
Additionally, the method accepts the following keyword arguments:

- ``minimum_exact_match``: a value between 0.0 and 1.0 specifying the minimum fraction of
  matching nucleotides needed to align sequence [default: 0.0]
- ``max_alternative_alignments``: max number of different alignments to return [default: 1]
- ``max_num_nodes_per_seq_char``: maximum number of nodes to consider during extension [default: 10.0]

::

    metasub.align(query, min_exact_match=0.8)

Search multiple graphs in parallel
----------------------------------
The API provides ``MultiGraphClient``, which can query multiple graph servers in parallel.
Both ``search`` and ``align`` have the keyword argument ``parallel`` [default: True].
If ``parallel=True``, the result will be a dictionary mapping specified index names to instances
of ``concurrent.futures.Future``.
If ``parallel=False``, all graphs will simply be queried in sequence and the results will
be ``instances of pandas.DataFrame``.

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

Other examples
--------------

Find more examples `here <https://github.com/ratschlab/metagraph_paper_resources/blob/master/notebooks/>`_.
