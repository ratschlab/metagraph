.. _traverse:

Labeled traversal (experimental)
================================

Labeled traversal extends a sequence in the graph along the labels that contain it. Use it to
get the flanks of a gene in each record or sample that carries the gene.

It has one server route:

.. list-table::
   :header-rows: 1

   * - Route
     - Method
     - Answer
   * - ``/traverse``
     - POST
     - the extension of one or more seeds to the left and to the right

The route is experimental. Its request and response fields can change in a later release
without notice.

A server on a single index serves the route:

.. code-block:: text

    metagraph server_query -i graph.dbg -a annotation.annodbg

So does a server started with a list of graphs. See `A list of graphs`_.

A request is a JSON object. The server refuses a field it does not know with HTTP 400, so a
misspelled field is never ignored.

Example
-------

Find the records that contain a query with ``/search`` and ``with_signature``:

.. code-block:: text

    curl -s http://localhost:5555/search -d '{
      "FASTA": ">query\nATGGAATTGCCC...",
      "with_signature": true,
      "discovery_fraction": 0
    }'

Each entry of ``results`` names a record in ``sample`` and gives its ``kmer_count``, the
number of query k-mers in that record. Its ``signature.presence_mask`` lists runs of k-mers
along the query, ``x`` for present and ``o`` for absent. For example, ``x423o31x329`` means
423 k-mers present, then 31 absent, then 329 present.

Then extend the query by 2,000 bp on each side. Use the query as the seed's ``sequence``. A
seed is at most 10,000 bp long, so send a part of a longer query. For the seed's ``labels``,
take the records that hold every k-mer of the seed. In their mask, the seed's k-mers lie
within one ``x`` run. When the seed is the whole query, the mask is a single ``x`` run:

.. code-block:: text

    curl -s http://localhost:5555/traverse -d '{
      "seeds": [{"sequence": "ATGGAATTGCCC...",
                 "labels": ["NZ_CM008882.1"]}],
      "strategy": {"bounds": {"max_extension_bp": 2000}}
    }'

One difference remains. ``/search`` checks the seed's k-mers one by one. ``/traverse`` keeps
a header label only when its record holds them at consecutive coordinates. ``/traverse`` drops
a record that holds them out of order and reports it in ``dropped_labels``. A seed left
without a label reports ``no_carrier``.

A list of graphs
----------------

A server started with a list of graphs walks every selected graph. The ``graphs`` field selects
them by name, as in ``/search``:

.. code-block:: json

    {
      "graphs": ["bacteria", "plants"],
      "seeds": [{"sequence": "ATGGAATTGCCC...", "labels": ["NZ_CM008882.1"]}]
    }

Without ``graphs``, the server walks every graph of a list with at most 10 names. A longer list
needs the field. A name can hold several graphs. All of them are walked.

Each entry of ``results`` then holds ``walks`` instead of one walk:

.. code-block:: json

    {
      "strategy": {...},
      "results": [
        {"walks": [
          {"graph": "bacteria", "seed": {...}, "label_dict": [...], "arms": {...}, "elapsed_ms": 12.5}
        ]}
      ]
    }

A walk has the fields of `The response`_ and the name of its ``graph``. ``walks`` holds one walk
per graph that has the seed and at least one of its labels. Each graph follows the labels it
has. A graph adds no walk in these cases:

- A k-mer of the seed is not in the graph.
- The graph has none of the named labels.
- No label is on every k-mer of the seed.

``walks`` is empty when no graph adds a walk. Any other error of a seed is listed as a walk with
an ``error``.

The server refuses a named label that no selected graph has with HTTP 400. The time budget
applies to each walk. The graphs are walked in parallel.

Terms
-----

Label
    What a walk follows. There are two kinds.

    A **column** label is an annotation column, for example a taxid or a file.

    A **header** label is one indexed sequence of a column, named by its FASTA header, for
    example an accession. Header labels need an index with k-mer coordinates and a sequence
    header index (``annotate --index-header-coords``). The server builds its index of header
    names on first use.

    A name is ambiguous when it is both a column and a header, or the header of several
    sequences. The server refuses a request that names it with HTTP 400.

Seed
    The sequence a walk starts from. Every k-mer of a seed must be in the graph.

Arm
    One direction of a walk. The right arm starts at the seed's last k-mer. The left arm
    starts at its first.

Path, segment
    A path is one walk from the seed to its end. Paths of an arm share their beginning, so an
    arm is a tree. A segment is one unbranched piece of that tree.

Interval
    A half-open range of k-mer positions ``[begin, end)`` in the seed.

``/traverse``
-------------

.. code-block:: json

    {
      "seeds": [
        {"sequence": "ATGGAATTGCCC...", "labels": ["NZ_CM008882.1", "562"]}
      ],
      "strategy": {
        "direction": "both",
        "labels": {"max_seed_labels": 1000, "seed_label_kind": "auto"},
        "branching": {"max_label_branches": 0},
        "bounds": {"max_extension_bp": 5000, "time_budget_ms": 30000}
      }
    }

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
     - Default
   * - ``seeds``
     - At most ``--traverse-max-seeds`` seeds. Each is at most 10,000 bp long.
     - required
   * - ``seeds[].labels``
     - The labels to follow, at most 100,000. Without it, the server takes the labels on every k-mer of the seed.
     - taken from the seed
   * - ``direction``
     - ``both``, ``left`` or ``right``.
     - ``both``
   * - ``labels.max_seed_labels``
     - Keep at most this many labels taken from the seed, at most 100,000. The server keeps the first ones in (column, sequence) order.
     - 1000
   * - ``labels.seed_label_kind``
     - The kind of the labels taken from the seed: ``column``, ``header`` or ``auto``. ``auto`` is ``header`` on an index with a header index.
     - ``auto``
   * - ``branching.max_label_branches``
     - How many times a column label can go on along more than one successor. An integer up to 1,000,000, or ``"unlimited"``.
     - 0
   * - ``bounds.max_extension_bp``
     - How far each arm extends.
     - 5000
   * - ``bounds.time_budget_ms``
     - The time for one seed, at least 1. A request can take this time once per seed. The server lowers it to ``--traverse-max-time-ms``.
     - ``--traverse-max-time-ms``

How a walk proceeds
^^^^^^^^^^^^^^^^^^^

The walk follows the named labels, or the labels on every k-mer of the seed. The server drops a
named label that is not on every k-mer of the seed. It reports the label with the runs of
k-mers it is on.

Each arm extends one base at a time. The two kinds of label are followed differently:

- A **column** label stays on a path while it is on every k-mer walked so far.
- A **header** label follows its record's k-mer coordinates. It goes on to the successor that
  holds the record's next coordinate. It ends at the record's first or last k-mer. A repeat
  inside the record is therefore walked through. If the record holds the seed several times,
  the walk starts from the place with the lowest coordinate.

Where the labels of a path go on along different successors, the path splits. Paths that meet
again stay apart, so each arm is a tree of segments.

A label leaves a path with one of these reasons:

.. list-table::
   :header-rows: 1

   * - Reason
     - Meaning
   * - ``dead_end``
     - The node has no successor. ``hairpin`` instead when its only successor was a hairpin step.
   * - ``label_lost``
     - Successors exist, and none carries the label.
   * - ``record_end``
     - The record of a header label ends here, and the graph goes on.
   * - ``branch``
     - A column label is on several successors, and its ``max_label_branches`` are used up.
   * - ``edge_reuse``
     - The next step of a column label takes an edge that this path took before.
   * - ``edge_reuse_rc``
     - The same, with the edge in the opposite orientation (an inverted repeat).
   * - ``rejoined_seed``
     - The next step of a column label goes back into the seed, for example in a circular record.
   * - ``max_extension_bp``
     - The arm reached its extension length.
   * - ``time_budget``, ``max_steps``, ``max_live_paths``, ``max_paths``, ``max_splits_per_path``
     - A bound stopped the walk. See below.

**Column labels on graphs with cycles.** A column label could circle forever. Three rules end
its walk. A path never takes an edge, that is a (k+1)-mer, twice. In canonical graphs this
holds in either orientation. A path never goes back into the seed. A path never takes a hairpin
step, which is a step whose (k+1)-mer equals its own reverse complement. These rules bind
column labels only. A header label ends because its record does. On a step that the rules
forbid, the column labels of the path end and its header labels go on.

Column labels have been walked on test data only. Every traversal run on the staging server so
far named records.

**Bounds.** Fixed limits bound the work for one seed, besides the extension length and the time
budget:

- 1,000,000 steps. Reaching this limit, or the time budget, stops the whole seed.
- 1000 live paths and 10,000 paths per arm. Reaching one of these limits stops that arm.
- 64 splits per path. Reaching this limit ends that path only.

Any of these limits marks the arm ``truncated``. ``complete_to_bp`` is then the extension up to
which every walk is present.

The response
^^^^^^^^^^^^

.. code-block:: json

    {
      "strategy": {...},
      "results": [
        {
          "seed": {
            "length_bp": 813, "num_kmers": 783,
            "labels": ["NZ_CM008882.1"],
            "labels_from_seed": false, "labels_supporting_total": 1,
            "labels_dropped": 0,
            "dropped_labels": [{"label": "562", "runs": [[0, 231]]}]
          },
          "label_dict": [{"name": "NZ_CM008882.1", "kind": "header", "column": 6, "seq_id": 3}],
          "arms": {
            "left": {...},
            "right": {
              "status": "complete",
              "complete_to_bp": 5000,
              "frontier_remaining": {"live_paths": 0},
              "segments": [
                {"id": 0, "parents": [], "from_bp": 0, "length_bp": 1200, "sequence": "GCAT...",
                 "labels": [0], "labels_at_end": [0],
                 "events": [{"at_bp": 1200, "type": "label_end", "label": 0, "reason": "record_end"}]}
              ],
              "splits": [],
              "paths": [{"id": 0, "segments": [0], "length_bp": 1200,
                         "end_labels": [{"label": 0, "reason": "record_end", "branches": 0}]}],
              "steps": 1200
            }
          },
          "elapsed_ms": 12.5
        }
      ]
    }

``results`` holds one entry per seed, in the order of the request.

The seed:

- ``seed.labels`` holds the names of the labels the walk follows.
- ``labels_dropped`` counts the labels taken from the seed beyond ``max_seed_labels``.
- ``dropped_labels`` lists the named labels that are not on every k-mer of the seed.

An arm:

- In an arm, a ``label`` field is an index into ``label_dict``. So are the entries of a
  segment's ``labels`` and ``labels_at_end``.
- A segment's ``sequence`` reads left to right on both arms. On the right arm, a segment comes
  before its children. On the left arm, it comes after them, because the children lie further
  from the seed.
- ``from_bp`` is the distance from the seed at which a segment starts. That is its first base on
  the right arm and its last base on the left arm.
- A path gives its ``segments`` from the seed outwards and its ``length_bp``. ``end_labels`` are
  the labels that reach its end, each with its reason and with ``branches``, the number of times
  it went on along several successors. ``path_reason`` is set when a bound ended the path.
- A split has the kind ``divergence`` when its labels part. It has the kind ``ambiguous`` when
  one label goes on along several children.
- ``cap_trigger`` gives the ``reason``, ``at_bp`` and ``demand`` of the bound that stopped a
  ``truncated`` arm.

A seed whose labels cannot be determined gets an error in its place. The other seeds are still
walked.

.. code-block:: json

    {"seed": {"length_bp": 813, "labels_from_seed": true},
     "error": {"cause": "no_carrier", "message": "...", "limit": 783, "observed": 412}}

.. list-table::
   :header-rows: 1

   * - ``cause``
     - Meaning
   * - ``no_carrier``
     - No label, named or taken from the seed, is on every k-mer of the seed.
   * - ``too_wide``
     - Too many labels are on the seed. The server checks the k-mer with the fewest labels among the first 64.
   * - ``time_budget``
     - The time budget ran out before the labels were determined.
   * - ``ambiguous_header``
     - A header name taken from the seed is ambiguous.

Any other error fails the whole request with HTTP 400. Examples are an unknown label, a seed
with a k-mer that is not in the graph, and too many seeds.
