.. _troubleshooting:

===============
Troubleshooting
===============

Find a symptom below, check its cause, and use the matching recovery step.
For command syntax and prerequisites, see :ref:`quick_start` and :ref:`installation`.

Installation
============

See :ref:`installation` for platform-specific troubleshooting.

- **Conda dependency conflict:** Try a fresh environment with the channels and priority in :ref:`installation`.
- **Missing Docker image:** Pull it with ``docker pull ghcr.io/ratschlab/metagraph:master``.
- **Unsupported C++ features:** Use GCC 13 or newer, or Clang, with a standard library supporting C++20, as required by the source build.

Memory Issues
=============

Out of memory
-------------

**Symptoms:** Allocation fails, or the operating system kills the process. A killed process alone
does not establish an out-of-memory cause; check the system or scheduler logs.

**Solutions:**

1. **Use disk swap with a buffer that fits your memory budget:**

   .. code-block:: bash
   
      # Graph construction
      metagraph build --disk-swap /tmp --mem-cap-gb 4 -o graph input.fa
      
      # Annotation
      metagraph annotate -i graph.dbg --anno-filename --disk-swap /tmp \
          --mem-cap-gb 4 -o annotation input.fa

2. **Reduce simultaneous annotation jobs:**

   .. code-block:: bash
   
      metagraph annotate -i graph.dbg --anno-filename --separately \
          -p 1 --threads-each 4 -o annotation input_1.fa input_2.fa

3. **Allow for memory beyond the buffer:** The graph, annotation, worker threads, and
   compression stages also consume RAM. Their sizes depend on the representation and input;
   see :ref:`faq` for graph-size estimates and :ref:`to_multi_brwt` for BRWT construction estimates.

Choose a swap directory with enough disk space. If ``/tmp`` is backed by RAM, use a disk-backed
scratch directory instead.

Memory cap ignored
------------------

**Issue:** MetaGraph uses more RAM than ``--mem-cap-gb`` specifies

``--mem-cap-gb`` controls internal buffers, not a strict limit on the process's total RAM.
Lower it and leave room for the graph, annotation, and other allocations.

Annotation Issues
=================

Missing k-mer counts file
--------------------------

**Warning:** ``No k-mer counts found ...``

The matching ``.kmer_counts.gz`` file was not found. MetaGraph then tries counts in sequence
headers; if the first header has no count, it warns and uses count 1 for each sequence.
To preserve pre-counted abundances, generate and keep both contig files together:

**Solution:**

.. code-block:: bash

   # Build weighted graph
   metagraph build -k 31 --count-kmers -o graph input.fa
   
   # Transform creates both files
   metagraph transform --to-fasta -o contigs graph.dbg
   # Creates: contigs.fasta.gz + contigs.kmer_counts.gz

Annotation fails
----------------

**Issue:** Runs out of memory

**Solutions:**

1. **Use disk swap:**

   .. code-block:: bash
   
      metagraph annotate -i joint.dbg --anno-filename --separately \
          --disk-swap /tmp --mem-cap-gb 4 -o annotation genome.fasta.gz

2. **Annotate multiple files separately:**

   .. code-block:: bash
   
      find . -name '*.fasta.gz' | metagraph annotate -i joint.dbg \
          --anno-filename --separately --disk-swap /tmp \
          --mem-cap-gb 4 -o ./out/annotation

   Standard input takes one filename per line, not FASTA contents.

3. **For binary sample-presence annotations, annotate contigs instead of raw sequences:**

   .. code-block:: bash

      # Build sample graph
      metagraph build -k 31 -o sample.graph genome.fa

      # Extract contigs
      metagraph transform --to-fasta -o sample.contigs sample.graph.dbg

      # Label k-mers by their source file
      metagraph annotate -i joint.dbg --anno-filename -o annotation sample.contigs.fasta.gz

Transform Errors
================

RowDiff transform fails
-----------------------

**Error:** Graph required for row-diff transform

**Solution:** Pass the matching graph with ``-i graph.dbg`` at each row-diff stage.
Follow the complete :ref:`row-diff procedure <to_row_diff_brwt>`, including its output-directory requirements.

BRWT out-of-memory
------------------

Column clustering and BRWT construction have different memory requirements (see :ref:`to_multi_brwt`).
Reduce ``--subsample`` for clustering, or ``--parallel-nodes`` for the construction stage::

   metagraph transform_anno --anno-type brwt --parallel-nodes 1 \
       --greedy -o annotation input.column.annodbg

Other options:

.. code-block:: bash

   # Option 1: Reduce subsampling
   metagraph transform_anno --anno-type brwt --subsample 500000 \
       --greedy -o annotation input.column.annodbg

   # File-count issue: disable temporary files for very many columns (around 1M)
   metagraph transform_anno --anno-type brwt --disk-swap "" \
       --greedy -o annotation input.column.annodbg

Missing transform output
------------------------

**Issue:** ``kmer_counts.gz`` not created

**Cause:** Graph not weighted

**Solution:** Build with ``--count-kmers``:

.. code-block:: bash

   metagraph build -k 31 --count-kmers -o graph input.fa
   metagraph transform --to-fasta -o contigs graph.dbg
   # Now creates both .fasta.gz and .kmer_counts.gz

Alignment Issues
================

Empty sequence error
--------------------

**Error:** Assertion fails on empty sequence

**Solution:** Filter empty sequences:

.. code-block:: bash

   seqkit seq -m 1 input.fa > filtered.fa
   metagraph align -i graph.dbg filtered.fa

Reverse complement fails
-------------------------

**Issue:** RC alignment doesn't work

**Solution:** Use canonical or primary graph:

.. code-block:: bash

   # Build canonical
   metagraph build -k 31 --mode canonical -o graph input.fa
   
   # Transform to primary
   metagraph transform --to-fasta --primary-kmers -o contigs graph.dbg
   metagraph build -k 31 --mode primary -o graph_primary contigs.fasta.gz

API Issues
==========

Python API Connection Fails
----------------------------

**Error:** RuntimeError or connection refused

**Solution:**

1. **Ensure server is running:**

   .. code-block:: bash
   
      metagraph server_query -i graph.dbg -a annotation.column.annodbg \
          --port 5555 -p 10

2. **Use correct connection:**

   .. code-block:: python
   
      from metagraph.client import GraphClient
      
      # For local server
      client = GraphClient('localhost', 5555, api_path='')
      
      # Test connection
      labels = client.column_labels()
      print(f"Connected! {len(labels)} labels")

API not installed
-----------------

**Error:** ``ModuleNotFoundError: No module named 'metagraph'``

**Solution:**

.. code-block:: bash

   python -m pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/api/python"

Docker
======

**libhts.so.3 error:** Pull latest image:

.. code-block:: bash

   docker pull ghcr.io/ratschlab/metagraph:master

**Permission denied:** Add user to docker group or use sudo

Performance Tips
================

Speed up construction
---------------------

1. **Use your allocated cores:** Set ``-p`` within the CPU allocation for your job
2. **Pre-process with KMC:** Faster k-mer counting
3. **Use SSD for disk swap:** Faster I/O
4. **Process different samples in parallel:** GNU parallel or workflow manager

Reduce Index Size
-----------------

1. **Use primary graph:** 50% smaller than canonical
2. **Use RowDiff<BRWT>:** 10-20% of column size
3. **Filter k-mers:** with KMC (e.g., ``kmc -ci5``) or graph cleaning
4. **Compress graph:** ``metagraph transform --state small``

Getting Help
============

If your issue isn't listed:

1. Check :ref:`faq` and search the `GitHub issues <https://github.com/ratschlab/metagraph/issues>`_.
2. If no existing answer applies, open an issue with the command, ``metagraph --version``,
   platform, relevant error output, and a small input that reproduces the problem.
