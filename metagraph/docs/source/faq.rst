.. _faq:

==========================
Frequently Asked Questions
==========================

General
=======

What is MetaGraph?
------------------

MetaGraph is a framework for scalable construction and querying of very largeannotated de Bruijn graphs.
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
