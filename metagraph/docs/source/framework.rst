Framework
=========

[TODO -- This page is in development]

Graph construction
------------------
Graph construction in external memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
construction using sorted set disk ...

Distributed graph construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
construction from chunks ...

Graph cleaning
--------------

Graph merging
-------------

Graph annotation
----------------
Column construction
^^^^^^^^^^^^^^^^^^^
get columns ...

Column compression
^^^^^^^^^^^^^^^^^^
diff ...

Annotation conversion
^^^^^^^^^^^^^^^^^^^^^
multi BRWT ...

Row-compressed construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Canonical indexes
-----------------
... for indexing raw read data

Canonical graphs
^^^^^^^^^^^^^^^^
Primary graphs
^^^^^^^^^^^^^^
Typical indexing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^




Supported alphabets
-------------------
Graph representations
---------------------
DBGSuccinct
^^^^^^^^^^^
- Properties
- Operation complexity

DBGHash
^^^^^^^
DBGBitmap
^^^^^^^^^

Annotation representations
--------------------------
Column-Compressed
^^^^^^^^^^^^^^^^^
Row-Compressed
^^^^^^^^^^^^^^
Rainbowfish
^^^^^^^^^^^
Multi-BRWT
^^^^^^^^^^
Rainbow-BRWT
^^^^^^^^^^^^
RowDiff annotation
^^^^^^^^^^^^^^^^^^
Benchmarks
----------
Benchmarking results [results of performing simple operations, like get_outgoing_edges, get_row, etc. Should be very important for developers]

Bit vectors
^^^^^^^^^^^
Graph routines
^^^^^^^^^^^^^^
Annotation routines
^^^^^^^^^^^^^^^^^^^
- Get 1 row
- Get batch of 1M random rows
- Get column

