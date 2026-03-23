.. _troubleshooting:

===============
Troubleshooting
===============

This page provides solutions to the most common MetaGraph issues.

Installation
============

See :ref:`installation` for platform-specific troubleshooting.

**Quick fixes:**

- **Conda fails:** ``conda install -c bioconda metagraph`` or fresh environment
- **Docker error:** pull the latest image: ``docker pull ghcr.io/ratschlab/metagraph:master``
- **Compiler too old:** Need GCC 9+ or Clang 9+

Memory Issues
=============

Out of memory
-------------

**Symptoms:** Process killed, allocation failed

**Solutions:**

1. **Use disk swap (recommended):**

   .. code-block:: bash
   
      # Graph construction
      metagraph build --disk-swap /tmp --mem-cap-gb 50 input.fa
      
      # Annotation
      metagraph annotate --disk-swap /tmp --mem-cap-gb 10 ...

2. **Process files separately:**

   .. code-block:: bash
   
      metagraph annotate --separately -p 4 --threads-each 9 ...

3. **Estimate memory needs:**

   - Graph: ~3 bits/k-mer
   - Annotation: ~0.5 bits/k-mer/label
   - RowDiff transform: graph + 2× annotation

Memory cap ignored
------------------

**Issue:** MetaGraph uses more RAM than ``--mem-cap-gb`` specifies

**Solution:** Lower the cap value -- it is rather a buffer size than a strict limit

Annotation Issues
=================

Missing k-mer counts file
--------------------------

**Error:** ``No k-mer counts found ...``

**Cause:** Missing ``.kmer_counts.gz`` file

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
   
      metagraph annotate --separately --disk-swap /tmp \
          --mem-cap-gb 10 -o annotation genome.fasta.gz

2. **Annotate multiple files separately:**

   .. code-block:: bash
   
      cat *.fasta.gz | metagraph annotate --separately \
          --disk-swap /tmp \
          --mem-cap-gb 10 -o ./out/annotation

3. **Annotate contigs instead of raw sequences:**

   .. code-block:: bash

      # Build sample graph
      metagraph build -k 31 -o sample.graph genome.fa

      # Extract contigs
      metagraph transform --to-fasta -o sample.contigs sample.graph.dbg

      # Annotate with coordinates
      metagraph annotate -i joint.dbg -o annotation sample.contigs.fasta.gz

Transform Errors
================

RowDiff transform fails
-----------------------

**Error:** Graph required for row-diff transform

**Solution:** Always provide graph with ``-i`` to the RowDiff transform:

BRWT out-of-memory
------------------

**Solution:** Reduce subsampling or use disk swap:

.. code-block:: bash

   # Option 1: Reduce subsampling
   metagraph transform_anno --anno-type brwt --subsample 500000 \
       --greedy -o annotation input.column.annodbg

   # Option 2: Disable temp files if there are too many columns (>1M)
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
   
      metagraph server_query -i graph.dbg -a annotation.annodbg \
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

   pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/api/python"

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

1. **Use all cores:** ``-p $(nproc)``
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

1. **See FAQ:** :ref:`faq`
2. **Search in** `GitHub issues <https://github.com/ratschlab/metagraph/issues>`_
3. **Prompt AI** with this documentation and ask for help
4. **Open an issue** on GitHub https://github.com/ratschlab/metagraph
