.. _pattern:

Pattern search (experimental)
=============================

Pattern search finds the k-mers of a graph that contain a short pattern. A pattern can hold
ambiguity codes, for example ``ACGTNNRYACGT``. Use it to find every place of a motif in the
graph, also where its flanks are not known.

It has one server route:

.. list-table::
   :header-rows: 1

   * - Route
     - Method
     - Answer
   * - ``/pattern``
     - POST
     - the number of contexts of each pattern, and the contexts when they are few

The route is experimental. Its request and response fields can change in a later release
without notice.

A server on a single index serves the route:

.. code-block:: text

    metagraph server_query -i graph.dbg -a annotation.annodbg

The graph must be succinct (``.dbg``). The server refuses the route on another graph with
HTTP 400.

A request is a JSON object. The server refuses a field it does not know with HTTP 400, so a
misspelled field is never ignored.

Example
-------

.. code-block:: text

    curl -s http://localhost:5555/pattern -d '{
      "patterns": [{"id": "p1", "sequence": "ACGTNNRYACGT"}],
      "max_contexts": 1000
    }'

.. code-block:: json

    {
      "k": 31, "mode": "basic", "mask": false,
      "results": [
        {"id": "p1", "length": 12,
         "count": {"min": 1, "max": 1},
         "complete": true,
         "contexts": [
           {"kmer": "TTGACGTAAGCACGTACCGGTTAACCGGTAC", "offset": 3, "strand": "+"}
         ]}
      ]
    }

To get the labels of the k-mers found, send them to ``/search``.

Terms
-----

Pattern
    A sequence of 1 to k symbols. On a DNA graph a symbol is one of the 15 IUPAC codes, in
    upper or lower case: ``A``, ``C``, ``G``, ``T``, ``R``, ``Y``, ``S``, ``W``, ``K``, ``M``,
    ``B``, ``D``, ``H``, ``V`` and ``N``. On a protein graph a symbol is a residue, and ``X``
    stands for any residue.

Context
    A k-mer of the graph and an offset. The pattern matches the symbols of the k-mer from
    ``offset`` to ``offset + length``. A k-mer that contains the pattern at two offsets is two
    contexts.

Strand
    ``+`` when the stored k-mer contains the pattern. ``-`` when it contains the reverse
    complement of the pattern. The ``offset`` then refers to the reverse complement of the
    pattern.

Mask
    A graph stores dummy k-mers beside its real ones. They start at the beginning of an indexed
    sequence and hold fewer than k real symbols. A mask marks them. ``build --mask-dummy --in-ram``
    writes it to the file ``graph.edgemask``, and the server loads it with the graph.

``/pattern``
------------

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
     - Default
   * - ``patterns``
     - 1 to 16 patterns.
     - required
   * - ``patterns[].sequence``
     - The pattern, of 1 to k symbols.
     - required
   * - ``patterns[].id``
     - A name for the pattern. The answer repeats it.
     - the index in ``patterns``
   * - ``max_contexts``
     - List the contexts of a pattern when there are at most this many. At most 10,000.
     - 1000
   * - ``time_budget_ms``
     - The time for the whole request. The server lowers it to ``--pattern-max-time-ms``.
     - ``--pattern-max-time-ms``

The server refuses a pattern longer than k with HTTP 400.

Which strands are searched
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Graph
     - Searched
   * - basic
     - The pattern and its reverse complement. A pattern equal to its reverse complement is
       searched once.
   * - canonical
     - The pattern only. The graph stores both orientations of every k-mer, so every context is
       reported with the strand ``+``.
   * - primary
     - The stored k-mers, with the pattern and its reverse complement. The graph stores one
       orientation of every k-mer.
   * - protein
     - The pattern only.

The response
------------

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``k``
     - The k-mer length of the graph.
   * - ``mode``
     - ``basic``, ``canonical`` or ``primary``.
   * - ``mask``
     - Whether the graph has a mask.
   * - ``results[]``
     - One entry per pattern, in the order of the request.
   * - ``results[].id``, ``length``
     - The pattern's name and its number of symbols.
   * - ``results[].count``
     - ``min`` and ``max`` bound the number of contexts. The count is exact when they are equal.
   * - ``results[].complete``
     - ``false`` when a bound stopped the search.
   * - ``results[].stop``
     - The bound that stopped the search: ``time_budget`` or ``max_steps``. Absent when
       ``complete`` is ``true``.
   * - ``results[].contexts``
     - Every context, when the search is complete and ``count.max`` is at most
       ``max_contexts``. Absent otherwise.

A context holds its ``kmer``, ``offset`` and ``strand``. The contexts are ordered by strand,
then by offset.

How exact a count is
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Case
     - ``count``
   * - The contexts are listed
     - Exact. ``min`` and ``max`` equal the number of contexts. The server reads every listed
       k-mer from the graph and drops the dummy ones.
   * - Not listed, with a mask
     - Exact, with one exception. A pattern that ends at the last symbol of the k-mer can match
       dummy k-mers. The server checks up to 10,000 of them per group of k-mers and states a
       bound beyond that.
   * - Not listed, without a mask
     - ``max`` counts every k-mer that matches, dummy ones included. ``min`` counts the
       contexts with offset 0, where no dummy k-mer can match.
   * - The search was stopped
     - ``min`` counts what was found. ``max`` is ``null``.

So a mask only changes the counts above ``max_contexts``. A list is the same with and without
it, when the search completes in both cases. Counting with a mask can take more steps, so a
search near its budget can stop with a mask and complete without one.

For even k, a k-mer can equal its own reverse complement. On a basic or primary graph such a
k-mer contains the pattern on both strands. A list holds it with the strand ``+`` only. A
count without a list counts it on both strands.

Bounds
^^^^^^

A request holds at most 16 patterns. The server lists at most 10,000 contexts per pattern. A
request takes at most 10,000,000 search steps and at most its time budget. All patterns of a
request share both bounds. A step narrows the search by one symbol or reads one k-mer.

A short pattern, or one with many ``N``, matches many k-mers and can reach a bound. Its result
then has ``"complete": false``, and so do the results of the patterns after it.
