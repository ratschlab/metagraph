.. _faq:

==========================
Frequently Asked Questions
==========================

General
=======

What is MetaGraph?
------------------

MetaGraph is a framework for scalable construction and querying of very large annotated de Bruijn graphs.
It uses succinct representations to index sequence collections efficiently, requiring only 2-4 bits per k-mer.

**Key capabilities:**

- Index trillions of k-mers with minimal memory
- Query sequences with exact matching or alignment
- Annotate k-mers with labels, counts, or coordinates
- Support DNA, RNA, and protein sequences

When should I use MetaGraph?
-----------------------------

**Good for:**

- Large-scale sequence search (millions of samples)
- Metagenomic classification
- Pangenome analysis
- Multi-sample comparison

**Not ideal for:**

- Single small datasets (MetaGraph might be overkill -- use BLAST)
- Full alignment profile or exact alignment (use dedicated aligners)

Installation
============

Which installation method?
--------------------------

1. **Conda** (easiest): ``conda install -c bioconda metagraph``
2. **Docker** (isolated): ``docker run ghcr.io/ratschlab/metagraph:master``
3. **Source** (for custom alphabets or latest features)

See :ref:`installation`.

Which alphabet?
---------------

- **DNA** (default): DNA/RNA sequences (A,C,G,T)
- **DNA5**: With ambiguous bases (A,C,G,T,N)
- **DNA_CASE_SENSITIVE**: When case matters (A,C,G,T,N,a,c,g,t)
- **Protein**: Amino acid sequences (A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,Y,Z,X)


Graph Construction
==================

What k-mer size?
----------------

**General rules:**

- DNA/RNA: k=21-41 (default: 31)
- Protein: k=7-15 (default: 10)
- Longer k = more specific, less sensitive
- Must be ≤ sequence/read length

How many threads to use?
------------------------

Most of the routines in MetaGraph are well parallelized, so the more threads, the better.
More than the number of cores is not recommended.

What buffer size?
-----------------
When constructing very large graphs, it is recommended to use disk swap
to limit the RAM usage::

   metagraph build --disk-swap <TMP_DIR> --mem-cap-gb 10

- For small graphs: 1 GB is enough (Default)
- For large graphs: the more, the better (50-80 GB is always sufficient, even for trillions of k-mers)
- Larger than 80 GB is **not** recommended
- ``<TMP_DIR>`` may be on a spinning disk (SSD might help but not necessary)

Canonical vs Primary graph?
----------------------------

**Workflow:**

1. Build **canonical** from reads (contains k-mer + reverse complement)
2. Extract primary contigs: ``metagraph transform --to-fasta --primary-kmers``
3. Build **primary** graph (half the size, same information)

**Important:** Use primary graphs with RowDiff annotations for best compression.

How to handle large datasets?
------------------------------

1. **Use disk swap:**

.. code-block:: bash

   metagraph build -k 31 --disk-swap /tmp --mem-cap-gb 50 input.fa

2. **Build sample graphs separately:**

.. code-block:: bash

   # Per sample
   for sample in samples/*.fasta.gz; do
      metagraph build -k 31 --mode basic -o $sample $sample
      metagraph transform --to-fasta --primary-kmers -o ${sample}.contigs ${sample}.dbg
   done

   # Joint graph from all contigs
   ls samples/*.contigs.fasta.gz | metagraph build -k 31 --mode canonical -o joint

   # Extract primary contigs from the joint graph
   metagraph transform --to-fasta --primary-kmers -o joint_contigs_primary joint.dbg

   # Joint primary graph from all contigs
   metagraph build -k 31 --mode primary -o joint_primary joint_contigs_primary.fasta.gz

3. **Extremely large graphs:**

Extremely large succinct graphs can be constructed by building their parts separately
and writing them to disk on the fly with flag ``--inplace``.
In such cases, don't forget to index suffix ranges afterwards with ``metagraph transform --index-ranges ...``.

Annotation
==========

Which annotation type?
----------------------

**By scale:**

- In most scenarios, start with: ``column`` (default)
- Fast and easy to construct: ``row_flat``
- Fast queries and small: ``row_diff_flat``
- Very large scale: ``row_diff_brwt`` (best compression)

**By feature:**

- K-mer counts: ``int_brwt`` or ``row_diff_int_brwt``
- Coordinates: ``brwt_coord`` or ``row_diff_brwt_coord``

How to annotate with counts?
-----------------------------

**Recommended (for large data):**

.. code-block:: bash

   # Build weighted sample graph
   metagraph build -k 31 --count-kmers -o sample.graph sample.fasta.gz

   # Extract contigs with counts
   metagraph transform --to-fasta -o contigs sample.graph.dbg
   # Creates: contigs.fasta.gz + contigs.kmer_counts.gz
   
   # Annotate joint graph
   metagraph annotate -i joint.dbg --count-kmers \
       --anno-filename -o annotation contigs.fasta.gz

**Simple (not recommended for large data):**

.. code-block:: bash

   metagraph annotate -i graph.dbg --count-kmers \
       --anno-filename -o annotation input.fa

**What if my k-mer counts are very large?**

Pass the ``--count-width`` flag to specify the number of bits used to represent the counts.
Default is 8 bits (max 255). E.g., with 12 bits, the max count is 4095:

.. code-block:: bash

   metagraph build -k 31 --count-kmers --count-width 12 -o sample.graph sample.fasta.gz

   metagraph transform --to-fasta -o contigs sample.graph.dbg

   metagraph annotate -i joint.dbg --count-kmers --count-width 12 \
       --anno-filename -o annotation contigs.fasta.gz

Querying
========

How to query?
-------------

**Command line:**

.. code-block:: bash

   # Presence/absence
   metagraph query -i graph.dbg -a annotation.annodbg query.fa
   
   # With threshold (80% k-mers must match)
   metagraph query --discovery-fraction 0.8 ...
   
   # K-mer counts
   metagraph query --query-mode counts -a annotation.int_brwt.annodbg ...
   
   # Coordinates
   metagraph query --query-mode coords -a annotation.brwt_coord.annodbg ...

**Python API:**

1. Start metagraph in the server mode:

.. code-block:: bash

   metagraph server_query -i graph.dbg -a annotation.annodbg --port 5555 -p 10 --mmap

2. Query the index:

.. code-block:: python

   from metagraph.client import GraphClient
   
   client = GraphClient('localhost', 5555, api_path='')
   results = client.search('ACGTACGT', discovery_fraction=0.8)

What is discovery_fraction?
----------------------------

Minimum fraction of k-mers that must match:

- **0.0**: Any match (at last 1 k-mer must match)
- **0.8**: ≥80% k-mers match (recommended for classification)
- **1.0**: All k-mers match (strict)

How to align?
-------------

**Exact k-mer matching:**

.. code-block:: bash

   metagraph align --map -i graph.dbg reads.fa

**Sequence-to-graph alignment:**

.. code-block:: bash

   metagraph align -i graph.dbg reads.fa
   
   # With annotation (label-consistent paths)
   metagraph align -i graph.dbg -a annotation.annodbg reads.fa
   
   # Adjust sensitivity
   metagraph align --align-min-seed-length 15 --min-exact-match 0.7 ...

Performance
===========

Memory requirements?
--------------------

**Graph:**

- Default construction: for k=31, at least 16 bytes per k-mer plus overhead (1B k-mers ≈ 18 GB)
- Construction with disk swap: buffer size plus output graph size
- Succinct (stat): ~4 bits/k-mer (10B k-mers ≈ 5 GB)
- Succinct (small): ~2 bits/k-mer (10B k-mers ≈ 2.5 GB)

**Large-scale indexing:**

When indexing large-scale datasets, everything depends on the chosen buffer sizes.
With carefully selected parameters, one can build and annotate graphs with trillion of k-mers.
For real examples, see https://github.com/ratschlab/metagraph/blob/master/metagraph/experiments/large_index_scripts.md.

How to speed up construction?
------------------------------

1. **Use more threads:** ``-p 64``
2. **Pre-process the input sequences with KMC:** ``kmc -k31 -m40 -sm input.fasta.gz output /tmp/``
3. **Parallelize samples:** Build sample graphs in parallel, with 4-8 threads (``-p ...``) on each.

How to reduce size?
-------------------

1. **Use primary graph:** 50% smaller than canonical (only applies to indexing reads)
2. **Use RowDiff<Multi-BRWT>** for annotation (10-20% of **Column**)
3. **Relax BRWT:** ``metagraph relax_brwt --relax-arity 32 ...`` (5-20% smaller)
4. **Transform graph:** ``metagraph transform --state small ...``
5. **Filter out low-abundance k-mers:** ``kmc -ci2 ...`` before building to remove singleton k-mers
6. **Use other graph cleaning techniques:** ``metagraph clean ...`` to remove sequencing errors

How to reduce RAM usage while querying?
--------------------------------

- **Use memory mapping (--mmap):** see :ref:`memory_mapping`
- **Reduce the batch size:** see ``./metagraph query --batch-size ...``
- **Alternative: Use disk-backed formats:** ``row_disk`` or ``row_diff_disk``

When to use memory mapping?
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- When the graphs are extremely large
- When the available RAM is limited

When not to use memory mapping?
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Never use it with slow (e.g., spinning) disks, unless it's for a single stats check

Troubleshooting
===============

Top issues and solutions:

Installation fails
------------------

**Conda:** Try installing in a fresh environment

**Docker:** Pull latest: ``docker pull ghcr.io/ratschlab/metagraph:master``

**Source:** Check compiler version and the dependencies

See :ref:`installation` and :ref:`troubleshooting`.

Out of memory
-------------

**Solution:** Use disk swap

.. code-block:: bash

   metagraph build --disk-swap /tmp --mem-cap-gb 50 ...
   metagraph annotate --disk-swap /tmp --mem-cap-gb 10 ...

RowDiff transform generates no output
-------------------------------------

**Solution:** Remember that for annotations with coordinates and counts, the output
files have the same name ``*.column.annodbg`` as the input files, hence, should be
written to a new directory. (See :ref:`transform_count_annotations`.)
