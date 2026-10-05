.. _graphlets:

Graphlets: retrieve once, query locally
=======================================

The ``/traverse`` endpoint of the MetaGraph query server extends a seed sequence through
the de Bruijn graph, to the left and to the right, and reports which labels (samples,
genomes, accessions) support each extension. An agent or an analyst usually asks many
questions about one such locus: list the walks, show one sample's route, where does the
support change, compare a constrained with a label-free retrieval, give me FASTA. None of
these needs the index. The Python package ``metagraph.traverse`` therefore works on a
single retrieval, the **graphlet**:

* one ``POST /traverse`` with ``strategy.output.detail: "graphlet"`` returns, per seed, a
  small JSON summary and the *whole* traversal as one compact, line-based text (the
  format **MGT v1**). The text is lossless: the ``detail: "full"`` JSON of the same
  request is reproduced from it exactly;
* every follow-up question is answered locally by the library, which never reads the
  graph;
* going further is a **new traversal**: the library builds the request from a walk's
  continuation, the server runs it. Nothing about the first retrieval is resumed.

The traversal itself (seed validation, the walk rule, label modes, termination reasons)
is specified in ``metagraph/docs/SPEC-labeled-traversal-core.md`` (§5 to §7); the
format and the library in ``metagraph/docs/DESIGN-traverse-graphlet.md``. MGT v1 is
frozen: later versions of the library and the server keep its records, fields and tokens.

The package is part of the Python API (see :ref:`install api`) and is imported as
``metagraph.traverse``. It needs Python 3.10 or newer and nothing outside the standard
library; ``pandas`` is imported only by ``metagraph.traverse.frames``.

.. contents:: On this page
   :local:
   :depth: 2


Choosing a retrieval
--------------------

The strategy fields below are request fields of ``/traverse`` (spec §5); the graphlet
records whichever were used, and the normalized strategy comes back in the response.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Strategy
     - Use it when
   * - ``labels.mode: constrain`` (default)
     - You want **claims**: which labels carry the flank, and how far. A walk is followed
       only while some permitted label supports every node of it. The permitted set is
       the seed's own carriers (every label supporting all seed k-mers) unless you name
       ``seeds[].labels``. This is the mode for answers.
   * - ``labels.mode: annotate``
     - You want to see what the graph holds regardless of labels: every structural
       successor is followed and the labels present at each node are merely
       *recorded*. It shares no label logic with ``constrain`` and is therefore an
       independent oracle for it, but every node costs a full annotation row: a
       verification tool at small radius, not a search primitive.
   * - ``exhaustive: true``
     - You need the complete trie: unlimited branches, ``on_reconverge: keep``, no
       quorum, no split limit, no beam. Its walks are complete **per path**, which is
       what comparisons and verification need. A conflicting knob is refused, never
       overridden.
   * - tuned (the default knobs)
     - You want a small, fast answer. Reconvergent paths are merged and quorum or branch
       limits may prune; its walks are prefixes of the exhaustive ones, and
       ``compare(mode='prefix_subset')`` lists each walk it omits with the reason it
       recorded (or says that none was recorded).
   * - ``annotate`` + ``frontier.on_overflow: beam``
     - You are exploring beyond the depth a complete trie reaches (in read data the
       label-free trie explodes within tens to hundreds of bases). Set
       ``bounds.max_live_paths`` to the beam width and ``frontier.order:
       most_supported_first``. A beam path is a set of candidate bases with per-node
       support, **not** a claim that any sample contains it; an arm the beam cut is
       ``pruned``.

The intended loop is: explore with a beam, read where the support changes, then either
constrain (naming the carriers you saw) or re-seed from a walk's continuation.


Getting a graphlet
------------------

Every example on this page runs as shown against the fixtures committed in
``metagraph/api/python/tests/data/traverse/documents/``, from that directory, as one
session. Most are written by the CLI from tiny indexes; ``clipped_merge`` and the
``compare/oracle_*`` documents are hand-made. The outputs are real; see `Running the
examples`_.

From a response
^^^^^^^^^^^^^^^

A ``detail: "graphlet"`` response is ordinary JSON: ``results[i]`` holds the seed's
summary with the MGT text under ``graphlet``, and the rest of the response is the
*envelope* (release, capabilities, the normalized strategy, the walk rule).
``Graphlet.from_response()`` parses the text and attaches both. It also checks
``graphlet_bytes`` and ``graphlet_lines``, so a body cut in transport raises
``GraphletFormatError`` instead of parsing as a smaller locus.

.. graphlet-example: setup

.. code-block:: python

   import json
   from metagraph.traverse import Graphlet

   def fixture(name):
       """A committed fixture as a client receives it: the detail "graphlet" response."""
       with open(name + '.graphlet.json') as f:
           response = json.load(f)
       return Graphlet.from_response(response['results'][0], response)

   merge = fixture('merge')
   print(merge.mode, merge.k, list(merge.arms), merge.seed.length_bp,
         [label.name for label in merge.labels])

.. code-block:: text

   constrain 15 ['right'] 30 ['a.fa', 'b.fa', 'c.fa', 'both.fa']

From a server
^^^^^^^^^^^^^

``TraverseClient`` issues the requests; it sets ``output.detail`` to ``"graphlet"``,
asks for gzip, and parses each result. A seed whose permitted set could not be derived
is an error *result*, not an exception: its slot in ``graphlets`` is ``None`` and the
reason is in ``errors``. Non-2xx answers raise ``TraverseError`` (``status``,
``message``, ``body``, and ``usage`` for a request with ``attempt_id``); a server still
loading its index raises ``ServerInitializing`` (a 503 with ``retry_after``: nothing was
run), and an attempt stopped at its duration bound while its response was built raises
``AttemptAtBound`` (a 503 with ``usage``: it ran and its id is used up, so it is not
retried under the same ``attempt_id``). ``capabilities()`` and ``resolve()`` call the other two routes
(on a multi-graph server ``capabilities(graph=..., graph_path=...)`` describes one of its
indexes), ``server_capabilities()`` reads the server-wide document (``GET /capabilities``:
the routes and features, ``feature_level``, the mode and graph names, the attempts),
``traverse_raw(request)`` sends a request you built yourself, and ``api_path`` is for a
server behind a proxy prefix.

.. graphlet-example: client server

.. code-block:: python

   from metagraph.traverse import TraverseClient

   client = TraverseClient('localhost', 5555)       # a server_query serving /traverse
   with open('merge.request.json') as f:             # the seed and strategy of a fixture
       request = json.load(f)
   response = client.traverse(request['seeds'], request['strategy'], timing=False)
   g = response.graphlets[0]                         # one Graphlet (or None) per seed
   print(response.envelope['strategy']['output']['detail'], response.errors)
   print(g.arm('right').complete_to_bp, len(g.arm('right').segments), g.has_envelope)

.. code-block:: text

   graphlet []
   100 7 True

(The runner answers this request with a stand-in that replays the committed response.)

A request can be an *attempt* of a ledger that reserved an allowance for it (the
search service's): ``traverse(..., attempt_id=..., budget_id=..., locus_id=...)``. The
server then registers it, enforces a duration bound on it (the seeds' time budgets
plus an allowance) and states what it consumed in a response-level ``usage`` block
(``response.usage``: work units, the modelled memory peak, seeds started and finished,
elapsed ms, the bound and why the walk stopped); an error after the request was read
carries it too, raised: in the error's ``usage`` (``TraverseError.usage``, also
``AttemptAtBound``'s). ``cancel(attempt_id, wait_ms=None)`` stops it from anywhere (the running seed is
returned partial, its ``resource_stop`` naming ``cancelled``; seeds not started are
error results with ``resource_stop.phase`` ``not_started``), and ``attempt(attempt_id)``
reads its state (``running``, ``stopping`` or ``finished``, with the reason and when the
walk stopped) for the retention period after it finished. Both return the server's
answer as a dict, a 404 included (``state`` ``finished`` or ``unknown``): a cancel of an
id the server does not know tombstones it (``tombstone: True``), refusing that id for the
retention period, so nothing runs under it there in that time; when the server already
holds its maximum of tombstones it answers 429 with ``tombstone: False`` and promises
nothing (the cancel can be retried). An
id runs once per server process (a second request with it is a 409). A client that
closes its connection stops its traversal: nothing is written. ``not_after_ms=`` (Unix
epoch ms, with or without ``attempt_id``) is the instant after which the request must not
be started: a server whose clock has passed it when the request arrives refuses it, and the
client raises ``AttemptExpired`` (a 409 with ``state: "expired"``: nothing ran, nothing was
registered, no usage).

.. code-block:: python

   response = client.traverse(seeds, strategy, attempt_id='a-17', budget_id='b-3')
   response.usage['work_units'], response.usage['reason']      # e.g. 1234, 'completed'
   client.cancel('a-18', wait_ms=2000)   # {'cancelled': True, 'state': 'finished', ...}
   client.attempt('a-18')['reason']      # 'cancelled'
   client.traverse(seeds, strategy, attempt_id='a-19', not_after_ms=now_ms + 60_000)

**Feature level 5: releasing capacity on a tombstone.** A ledger that gives an attempt's
capacity back before the attempt finished must know that no copy of the request can
still start: a cancel can overtake its own request (queued, on the wire, half
uploaded), and a delayed copy can reach a server that was restarted in the meantime and
holds no tombstones. Two request fields close both gaps, and the capabilities'
``release_rule`` says when a release is safe:

* ``cancel(attempt_id, wait_ms=None, not_after_ms=None)``: ``not_after_ms`` is the
  ``not_after_ms`` of the request being cancelled. A tombstone is then held until no
  copy of that request can start -- this server's clock past ``not_after_ms`` plus its
  ``clock_skew_allowance_ms`` --, within the server's ``tombstone_max_s``, and every
  tombstone answer states it.
* ``traverse(..., attempt_id=..., expect_server_instance=...)`` (with ``attempt_id``
  only): the ``server_instance`` the attempt is meant for, a string or ``'auto'`` (the
  client fills in the server's own, ``client.server_instance()``, read from
  ``GET /capabilities`` and kept for that server process). A process with another
  instance refuses the request before anything runs: ``InstanceMismatch`` (409, with
  ``expect_server_instance`` and ``server_instance``).

Both fields are sent only to a server that states ``feature_level`` 5 or more, read from
``GET /capabilities`` (``client.server_feature_level()``; or pass ``feature_level=`` to
the client) and never for a request without them: a server below level 5 refuses the
fields with a 400, and a 400 to a request with ``attempt_id`` has already used that id
up. Below level 5 a cancel goes out without ``not_after_ms`` (it still stops the attempt)
and ``'auto'`` is left out; an explicit instance raises ``UnsupportedFeature`` (a
``ValueError``) before anything is sent.

What the client read is kept per server process. It is read again after an answer names
another ``server_instance`` (the server restarted at the same address), after an
``InstanceMismatch``, and after a 400 that refuses ``not_after_ms`` or
``expect_server_instance`` as an unknown field (an older binary now answers there). A
cancel refused that way is sent again at once without ``not_after_ms``: the 400 stopped
nothing. A request with ``expect_server_instance`` that is the first to reach a replaced
process cannot be saved: an older binary uses its ``attempt_id`` up with that 400, so a
ledger that rolls servers back calls ``client.refresh_capabilities()`` first.

What went out is stated: ``response.sent`` and the ``.sent`` of any exception
``traverse()`` raises once the request is built (a ``TraverseError``, and a transport
error such as a timeout or a reset connection, the unanswered attempt a ledger must
cancel) are an ``AttemptSent(attempt_id, not_after_ms, expect_server_instance)``;
``AttemptSent.of(request)`` reads one from a request you built, and a cancel's answer
carries the body it sent (``answer.sent``).

The answers are typed. ``cancel()`` and ``attempt()`` return an ``AttemptAnswer``: the
server's JSON (a ``dict``, as before) with its HTTP ``status``, ``state``,
``tombstone``, the 429's ``reason`` (``tombstones_full``: the server holds its maximum,
``answer.retryable``; ``no_suppression``: it keeps no tombstones) and ``suppression``, a
``Suppression(suppressed_until_ms, not_after_ms, covers_admission,
covers_admission_reason)`` on every tombstone answer of a level-5 server
(``covers_admission`` true: no copy of that request can start on this
``server_instance``; else the reason, ``no_not_after_ms`` or ``beyond_tombstone_max``).
A request refused because its ``attempt_id`` is held raises ``AttemptConflict`` (409:
running, retained, held after it finished, or tombstoned -- ``tombstoned``, with the
``suppression`` judged against the refused request's own ``not_after_ms``; the refusal
extends the tombstone), and ``AttemptExpired`` states ``not_after_ms``,
``server_time_ms`` and ``server_instance``.

``release_verdict(answer, sent, *, now_ms=None, clock_skew_allowance_ms=None,
bound_ms=None, attempts=None)`` applies the ``release_rule`` exactly and returns a
``ReleaseVerdict`` (``release``, ``early``, ``code``, ``why``, ``assumptions``; true as a
``bool`` when ``release``). ``attempts`` is the capabilities' attempts block; it supplies
``clock_skew_allowance_ms`` and, as the bound, ``hard_cap_ms`` when they are not given:

* **early** (``code: 'tombstone'``) only on a 404 with ``tombstone: true`` and
  ``covers_admission: true`` (a cancel's or ``GET /traverse/attempt``'s) from the
  ``server_instance`` the attempt was sent to with ``expect_server_instance``, judged
  against exactly the ``not_after_ms`` it was sent with. An attempt sent without either
  field is never released early on a tombstone.
* otherwise on a **finished** state (``code: 'finished'``: the response, an error that
  carries the attempt's ``usage``, a cancel's or ``attempt()``'s ``state: finished``),
  or once the ledger's clock (``now_ms``) passes ``not_after_ms +
  clock_skew_allowance_ms + bound_ms`` (``code: 'clock'``).
* a finished release is the rule's either way, but the server holds a finished id
  against a replay only for an attempt sent with ``not_after_ms``, until
  ``not_after_ms + clock_skew_allowance_ms`` and at most ``tombstone_max_s`` after the
  finish, on that process. Where the release instead assumes that no copy of the request
  arrives later, ``assumptions`` names why: ``sent_without_not_after_ms``,
  ``no_retention`` (``retention_s`` 0), ``beyond_tombstone_max``,
  ``sent_without_expect_server_instance`` (a copy reaching a restarted process). The two
  about the server are judged only when ``attempts`` is given. When the clock has passed
  as well, the ``clock`` release is returned, which assumes none of them.
* nothing else releases, and ``code`` names the first condition not met:
  ``not_tombstoned`` (a 429), ``stopping`` / ``running``, ``no_tombstone``,
  ``sent_without_not_after_ms``, ``sent_without_expect_server_instance``,
  ``no_suppression_stated`` (a server below level 5), ``covers_admission_false``,
  ``other_server_instance``, ``other_not_after_ms``, ``answer_for_other_attempt``,
  ``not_a_release_answer`` (a 409 -- held, tombstoned, expired, ``instance_mismatch``
  -- or an error before the attempt was registered), ``no_answer``.

.. code-block:: python

   from metagraph.traverse import AttemptConflict, AttemptSent, release_verdict

   naf = now_ms + 60_000
   request = client.build_request(seeds, strategy, attempt_id='a-20', not_after_ms=naf,
                                  expect_server_instance='auto')
   sent = AttemptSent.of(request)        # what the ledger records before it sends
   # ... the request is on its way (traverse_raw(request) elsewhere) and does not answer
   answer = client.cancel('a-20', not_after_ms=naf)
   verdict = release_verdict(answer, sent)
   verdict.release, verdict.early, verdict.code   # True, True, 'tombstone' (it never arrived)
   try:
       client.traverse_raw(request)      # the delayed request arrives after the cancel
   except AttemptConflict as e:
       e.tombstoned, e.suppression.covers_admission            # True, True
   other = AttemptSent('a-21', naf, sent.expect_server_instance)
   release_verdict(client.cancel('a-21'), other).why   # cancelled without not_after_ms
   # 'covers_admission is false (no_not_after_ms)'
   done = client.traverse(seeds, strategy, attempt_id='a-22', not_after_ms=naf)
   verdict = release_verdict(done, done.sent,
                             attempts=client.server_capabilities()['attempts'])
   verdict.code, verdict.assumptions   # 'finished', ('sent_without_expect_server_instance',)

From the text alone
^^^^^^^^^^^^^^^^^^^

``parse(text)`` reads an MGT body. A body is self-contained for every query on this page
(seed, dictionary, orientation rule and guarantees are in it), but it does not carry the
seed id, the normalized strategy or the annotation counters. The operations that need
them raise ``MissingEnvelope`` instead of inventing defaults: ``summary()``,
``to_json()`` and ``next_request()`` (and so ``TraverseClient.deepen()``).

.. graphlet-example: body

.. code-block:: python

   from metagraph.traverse import parse, MissingEnvelope

   with open('merge.mgt') as f:
       body = parse(f.read())
   print(body.has_envelope, len(body.walks('right')))
   for call in (body.summary, body.to_json, lambda: body.next_request('right', [0])):
       try:
           call()
       except MissingEnvelope:
           print('MissingEnvelope')

.. code-block:: text

   False 1
   MissingEnvelope
   MissingEnvelope
   MissingEnvelope

Files
^^^^^

``save(path)`` writes the standalone ``.mgt`` file: the body plus one ``J`` line after
the header that holds the envelope reduced to this seed — its ``results`` to this seed's
summary, and an attempt's ``usage`` to its totals plus this seed's ``per_seed`` entry (the
other seeds' entries describe graphlets this one is not; the store keeps the same reduced
envelope). ``load(path)`` restores both (and restores a saved view as a view, see
`Views`_). Writing is atomic. ``view`` and ``derived_from`` are the library's own names in
``J``: a response that carries either at its top level is refused rather than misread. A
file saved by an earlier version of the library, whose ``J`` kept every seed's ``per_seed``
entry, loads as it is; it is not in today's canonical form (``is_canonical()`` is false for
it), and ``save()`` rewrites it reduced.

.. graphlet-example: files

.. code-block:: python

   import os, tempfile
   from metagraph.traverse import load

   path = os.path.join(tempfile.mkdtemp(), 'merge.mgt')
   merge.save(path)
   with open(path) as f:
       print([line[0] for line in f][:4])
   again = load(path)
   print(again.has_envelope, again.dump(envelope=False) == body.dump())

.. code-block:: text

   ['H', 'J', 'S', 'L']
   True True

``standalone_text(body, result, response)`` writes the same text without parsing the body:
it splices the ``H`` line, the ``J`` line and the rest of the body (the ``Z`` count one
higher), checking only the summary's byte and line counts and the ``H`` and ``Z`` records,
so that a large retrieval can be saved or spooled before (or without) building its model;
``load()`` validates it. For a server's body it is byte for byte what ``save()`` writes.

.. graphlet-example: standalone

.. code-block:: python

   from metagraph.traverse import dump, from_response, standalone_text

   with open('merge.graphlet.json') as f:
       response = json.load(f)
   result = response['results'][0]
   text = standalone_text(result['graphlet'], result, response)
   print(text == dump(from_response(result, response), envelope=True),
         text.split('\n')[1][0], len(text) - len(result['graphlet']))

.. code-block:: text

   True J 3441

Integrity
^^^^^^^^^

The parser validates as it reads: every token has exactly one valid spelling, the ``A``
record of each arm states how many segments, runs, leaves, splits, merges and bases
follow, and the final ``Z`` record the line count. A truncated or corrupt document
raises ``GraphletFormatError(line_no, message)``; it is never read as a smaller
graphlet. ``dump()`` writes the canonical text, so ``dump(parse(x)) == x`` for every
document a server writes (``is_canonical(x)`` tests it).

.. graphlet-example: integrity

.. code-block:: python

   from metagraph.traverse import GraphletFormatError, is_canonical

   with open('merge.mgt') as f:
       text = f.read()
   print(is_canonical(text), body.dump() == text)
   damaged = text.replace('R 5 3 0 62 m 0 * * * 2 0\n', '')     # one run lost
   try:
       parse(damaged)
   except GraphletFormatError as e:
       print(str(e).split(' (')[0])

.. code-block:: text

   True True
   line 9: the right arm does not match its A counts


The model in brief
------------------

A ``Graphlet`` holds the header fields (``k``, ``regime``, ``alphabet``, ``mode``,
``support``, ``reconverge``, ``cap`` = ``labels.max_labels_per_node``,
``continuation_bp``, the index identity ``index_ns`` / ``index_fp`` /
``index_meta_fp``), the validated ``seed`` (``sequence``, ``length_bp``,
``num_seed_labels``, ...), the label dictionary ``labels``, the ``dropped`` seed labels,
the guarantees (``outcome``, ``limitations``, ``resource_stop``) and one ``Arm`` per
requested direction in ``arms``. Ids are the format's ordinals and never change.

Arms and orientation
^^^^^^^^^^^^^^^^^^^^

``g.arms`` holds ``'left'`` and/or ``'right'``; ``g.arm('l')`` and ``g.arm('right')``
select one (``g.arm()`` without an argument works when only one arm was requested).

Positions are **outward base indices** from the seed boundary: base ``i`` of an arm is
the ``(i+1)``-th base away from the seed, on both arms. Every position field
(``from_bp``, ``to_bp``, ``route_bp``, ``complete_to_bp``, ...) uses them, and every
segment stores its bases in **walking order** (``Segment.walk``), so outward base ``i``
of a segment is ``walk[i - from_bp]`` on both arms. The **natural orientation** is the
molecule as read along the seed: the right flank is the walking order, the left flank is
its reverse, and a whole molecule is ``natural(left) + seed + natural(right)``. In seed
coordinates outward base ``i`` is ``len(seed) + i`` on the right and ``-(i + 1)`` on the
left. ``spell(arm, walk, orientation='natural'|'walk', with_seed=False)`` spells either.

.. graphlet-example: orientation

.. code-block:: python

   fork = fixture('fork')                 # both arms, a split on each
   print(fork.seed.sequence)
   print(fork.spell('left', 0), fork.spell('left', 0, orientation='walk'))
   print(fork.spell('left', 0) + fork.seed.sequence + fork.spell('right', 1))

.. code-block:: text

   CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATT
   ATCCCTCACG GCACTCCCTA
   ATCCCTCACGCGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATTTTACACCACG

The left walk's natural spelling ends with the base next to the seed (``...CACG`` +
``CGTA...``); in walking order that base comes first.

Segments, walks and runs
^^^^^^^^^^^^^^^^^^^^^^^^

An arm is a trie (``on_reconverge: keep``) or a DAG (``merge``) of **segments**
(``arm.segments``): unbranched stretches with ``parents`` (none for the root, two or
more for a merge), ``from_bp`` / ``length_bp`` / ``end_bp``, the bases, the label sets
at their first and last node (``entry`` and ``end``, arrays of label ids, with
``entry_total`` the true count where a list was cut), the merge ``partition`` (which
labels came through which parent), the stored events, and, for a leaf, a ``leaf`` record
(the path-level end reason and the continuation).

A **walk** is a leaf. Its id (``path_id``, the ``walk`` of the tool functions) is the
leaf's ordinal in segment order, as in the server's ``paths[]``, and its segments are the
first-parent chain from the root: the **displayed** path. Behind a merge, a label may
have come through another parent; its *own* route is then not the displayed one (see
`Claims`_).

In ``constrain`` mode the arm also holds **runs** (``arm.runs``): one per stretch of a
label lineage, with ``label``, ``from_bp``, ``to_bp``, an end token (a reason code such
as ``X`` max_extension_bp or ``D`` dead_end, with a qualifier such as ``Rm`` =
``branch``/``minority``; ``Lw`` = it went on under another name after a switch; ``m`` =
closed by a merge, not ended), the switch it entered by (``from_label``, ``cost``,
``prev_run``) and its terminal ``loss`` and ``branches``. ``annotate`` mode tracks no
lineages and has no runs; it records the label sets along each segment instead
(``Segment.presence``).

The fixture ``merge`` used below: a 30 bp seed, the right arm to 100 bp, four column
labels, and two SNP bubbles that open at 20 and 46 and reconverge at 36 and 62
(``on_reconverge: merge``). ``a.fa`` takes the first allele of both bubbles, ``b.fa`` the
second of both, ``c.fa`` the first then the second, and ``both.fa`` is on every allele.

.. graphlet-example: segments

.. code-block:: python

   right = merge.arm('right')
   for s in right.segments:
       print(s.id, s.parents, s.from_bp, s.end_bp, s.walk[:6],
             [merge.labels[l].name for l in s.entry], 'leaf' if s.is_leaf else '')

.. code-block:: text

   0 () 0 20 TTTCCT ['a.fa', 'b.fa', 'c.fa', 'both.fa']
   1 (0,) 20 36 ATTACG ['a.fa', 'c.fa', 'both.fa']
   2 (0,) 20 36 TTTACG ['b.fa', 'both.fa']
   3 (1, 2) 36 46 ACCCAG ['a.fa', 'b.fa', 'c.fa', 'both.fa']
   4 (3,) 46 62 CCGGCG ['a.fa', 'both.fa']
   5 (3,) 46 62 GCGGCG ['b.fa', 'c.fa', 'both.fa']
   6 (4, 5) 62 100 GTCCGG ['a.fa', 'b.fa', 'c.fa', 'both.fa'] leaf

Labels
^^^^^^

A ``Label`` has ``kind`` (``'column'`` or ``'header'``), ``column``, ``seq_id`` (header
labels) and ``name``. Its ``ref`` -- ``c:<column>`` or ``h:<column>:<seq_id>`` -- is
the stable key that joins a label across retrievals of the same index; names are for
display, and need not be unique (an annotate retrieval can record two headers named
``ACC1`` from different columns, and a column may be *named* ``c:0``). Every result
therefore carries ``Label`` objects, and the tool functions emit ``{name, ref}``, never a
bare name or an id.

Arguments that select labels take **tagged selectors**: ``{'ref': 'h:1:0'}``,
``{'name': 'ACC2'}``, ``{'id': 2}``, or a ``Label``. A bare string is a convenience that
resolves only when exactly one label matches it as a ref *or* as a name; otherwise it
raises ``AmbiguousLabel`` (``UnknownLabel``, a subclass, when nothing matches).

.. graphlet-example: labels

.. code-block:: python

   from metagraph.traverse import AmbiguousLabel

   twins = fixture('same_name')          # annotate mode: two headers named ACC1
   for label in twins.labels:
       print(label.as_dict())
   print(twins.label({'ref': 'h:1:0'}).column, twins.label('ACC2').ref)
   try:
       twins.label('ACC1')
   except AmbiguousLabel as e:
       print(e)

.. code-block:: text

   {'name': 'ACC1', 'ref': 'h:0:0'}
   {'name': 'ACC1', 'ref': 'h:1:0'}
   {'name': 'ACC2', 'ref': 'h:0:1'}
   1 h:0:1
   'ACC1' matches 2 labels (h:0:0 'ACC1', h:1:0 'ACC1'): use {"ref": ...}


Guarantees
----------

Every response certifies what it covers and states every limit it hit, with the request
field to turn. A graphlet carries all of it in its own records (``A``, ``O``, ``K``),
so a saved file keeps its caveats, and every local answer can carry them too. Read these
before reading any walk.

What each arm certifies
^^^^^^^^^^^^^^^^^^^^^^^

``complete_to_bp``
   Every walk of at most this many bases that obeys the response's ``walk_rule`` is
   present, and every walk that ended before this depth ended for its reported reason.
   Walks reaching further may be present, but nothing is claimed about them. The arm's
   ``status`` is ``complete`` exactly when this equals the requested radius;
   otherwise ``truncated`` (a size cap; ``cap_trigger`` names it) or ``pruned`` (a
   beam).

``completeness_scope`` (``Arm.scope``)
   What "every walk" quantifies over. ``per_path`` (``on_reconverge: keep``): each walk
   is barred only from reusing its own (k+1)-mer edges. ``united_history``
   (``merge``): a merge unites the edge histories of the routes it joins, so a merged
   walk is also barred from edges another route used. Fewer walks are then present, and
   a structural block on a merged route need not be where the per-path trie ends a
   label. Route support itself is per label and unaffected.

``evidence.complete_to_bp`` (``Arm.evidence_complete_to_bp``)
   A second, independent boundary, for the *reasons*: below it every branch decision
   and every refusal (a successor not taken, and why) is in the arm's ``branch_events``;
   at or beyond it a reason may be among the events ``output.max_branch_events``
   dropped. ``None`` means no event was dropped.

The outcome: four dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``g.outcome`` (the ``O`` record) states, per seed, four independent guarantees:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Dimension
     - Values
     - Not ``complete`` when
   * - ``walks``
     - ``complete`` | ``partial`` | ``failed``
     - ``partial``: a ``walk_domain`` limitation (an arm was truncated or pruned), a
       ``scope`` limitation with at least one merge (complete for the united-history
       rule only), or ``seed_labels`` (walks only a cut carrier carries are missing).
       ``failed``: no traversal (the permitted set could not be derived; such a result
       has no graphlet).
   * - ``branch_diagnostics``
     - ``complete`` | ``cut``
     - ``cut``: a ``branch_events`` limitation; reasons at or beyond an arm's
       ``evidence.complete_to_bp`` may be missing.
   * - ``label_evidence``
     - ``complete`` | ``lower_bound`` | ``qualified``
     - ``lower_bound``: evidence may be *missing or understated* (``label_lists``,
       ``inexact_counts``, ``seed_labels``, ``switch_sources``, ``greedy_losses``).
       ``qualified``: something may be *overstated* (``trace_record_boundaries``); it
       wins when both apply.
   * - ``delivery``
     - ``inline`` | ``spooled`` | ``paged``
     - How the body reached you; independent of the other three. Servers write
       ``inline``; the tool layer marks a body it spooled to disk instead of keeping
       it in memory ``spooled``.

The rule is **conservative by construction**: a dimension reads ``complete`` only when no
limitation of its class applies, so you can never find a limitation whose dimension
still says ``complete``, nor a non-complete dimension without the limitation that
explains it.

Limitations
^^^^^^^^^^^

``g.limitations`` lists them (the ``K`` records; ``arm`` is ``None`` for a seed-level
one). Each has ``kind``, ``knob`` (the request field, relative to ``strategy``),
``limit`` (its value in this run), ``observed`` (what the run met against it),
``complete_to_bp`` where the limitation has a depth boundary, further fields in
``extra``, and ``effect``: one sentence saying what is missing and how to get it.

.. list-table::
   :header-rows: 1
   :widths: 22 10 20 48

   * - Kind
     - Level
     - Dimension
     - Knob, and when it is stated
   * - ``walk_domain``
     - arm
     - walks
     - The cap that stopped the arm: ``bounds.max_steps``, ``bounds.max_paths``,
       ``bounds.max_output_bp``, ``bounds.time_budget_ms`` or ``bounds.max_live_paths``
       (which is also a beam's width); ``complete_to_bp`` = the arm's.
   * - ``scope``
     - arm
     - walks, if merges > 0
     - ``branching.on_reconverge`` (limit ``merge``); ``observed`` = merges done. With
       0 merges it is informational and limits nothing.
   * - ``seed_labels``
     - seed
     - walks, label_evidence
     - ``labels.max_seed_labels``: the derived permitted set was cut.
   * - ``branch_events``
     - arm
     - branch_diagnostics
     - ``output.max_branch_events`` (accepts ``"unlimited"``); ``complete_to_bp`` =
       ``evidence.complete_to_bp``.
   * - ``label_lists``
     - arm
     - label_evidence
     - ``labels.max_labels_per_node``: recorded lists were cut (each keeps its true
       count).
   * - ``inexact_counts``
     - arm
     - label_evidence
     - ``labels.max_labels_per_node``: live-label counts flagged ``exact: false``.
   * - ``switch_sources``
     - arm
     - label_evidence
     - ``labels.max_switch_sources``: a ``table`` cost's source list was cut where a cut
       source could have switched; losses may be overestimated.
   * - ``greedy_losses``
     - arm
     - label_evidence
     - ``branching.max_label_branches``: losses were re-minimised greedily after a
       branch-limit exclusion and need not be optimal.
   * - ``trace_record_boundaries``
     - seed
     - label_evidence (qualified)
     - ``labels.seed_label_kind``: ``support: trace`` with column labels cannot see a
       record boundary.
   * - ``server_clamp``
     - seed
     - none of its own
     - The clamped field; what the clamp caused is stated by the ``walk_domain`` or
       ``seed_labels`` entry it led to.

Semantic ends (``dead_end``, ``label_lost``, ``max_extension_bp``, ...) are not
limitations: the requested domain is complete there.

.. graphlet-example: outcome

.. code-block:: python

   for name in ('merge', 'caps', 'annotate'):
       g = fixture(name)
       print(name, ' '.join('%s=%s' % kv for kv in g.outcome.as_dict().items()))
       for lim in g.limitations:
           print('    %-5s %-14s %-26s limit=%s observed=%s to=%s' % (
               lim.arm or 'seed', lim.kind, lim.knob, lim.limit, lim.observed,
               lim.complete_to_bp))

.. code-block:: text

   merge walks=partial branch_diagnostics=complete label_evidence=complete delivery=inline
       right scope          branching.on_reconverge    limit=merge observed=2 to=None
   caps walks=partial branch_diagnostics=cut label_evidence=complete delivery=inline
       right walk_domain    bounds.max_steps           limit=240 observed=241 to=126
       right branch_events  output.max_branch_events   limit=1 observed=2 to=95
   annotate walks=partial branch_diagnostics=complete label_evidence=lower_bound delivery=inline
       left  label_lists    labels.max_labels_per_node limit=2 observed=4 to=None
       left  inexact_counts labels.max_labels_per_node limit=2 observed=1 to=None
       left  scope          branching.on_reconverge    limit=merge observed=0 to=None
       right label_lists    labels.max_labels_per_node limit=2 observed=4 to=None
       right inexact_counts labels.max_labels_per_node limit=2 observed=2 to=None
       right scope          branching.on_reconverge    limit=merge observed=2 to=None

``merge`` walked its whole radius, but two merges united histories, so its walks are
complete only for the united-history rule. ``caps`` hit ``max_steps``: everything is
certified to 126 bp, and branch reasons only below 95 bp. ``annotate`` cut its recorded
lists at 2 labels per node, so everything derived from them is a lower bound.

.. graphlet-example: arm-certificate

.. code-block:: python

   import textwrap

   caps = fixture('caps')
   a = caps.arm('right')
   print(a.status, a.complete_to_bp, a.scope, a.evidence_complete_to_bp)
   t = a.cap_trigger                      # the cap that stopped the arm (S = max_steps)
   print(t.reason, t.at_bp, t.demand)
   print(textwrap.fill(caps.limitations[0].effect, 88))
   for w in caps.walks('right', by='id'):
       print(w.path_id, w.length_bp, w.path_reason, w.complete, w.beyond_certified_bp)

.. code-block:: text

   truncated 126 per_path 95
   S 126 241.0
   every walk of at most complete_to_bp bp is present, longer walks may be missing: the
   exploration stopped at this cap (max_steps); raise the knob
   0 126 max_steps False 0
   1 109 None True 0
   2 95 None True 0

A walk is ``complete`` when it lies within ``complete_to_bp`` and was not ended by a
cap; walk 0 was censored by ``max_steps`` at the certified depth. A capped or pruned
result is never evidence of absence.

``g.resource_stop`` (the ``Q`` record) says where a request budget stopped the walk
(``bounds.max_memory_mb`` or ``bounds.max_work_units``): the resource, the phase, what
was requested, used and left, and the actions that would help. The stopped arm's walks
are partial, and its ``walk_domain`` limitation names the budget as the knob. The phase
is ``annotation_decode`` when, on a row-diff annotation, reading the next level's
annotation with its row-diff dependency rows did not fit the memory left: such a read is
refused whole, and the levers are those that read fewer rows (``more_selective_seed``,
in annotate mode ``label_constrained_query``) -- a smaller radius would not help. In
annotate mode a row that fits but names more new labels than the memory left holds stops
in phase ``traversal``, with ``lower_max_labels_per_node`` among its levers; the message
always says which of the two did not fit.

These budgets bound the backend's walk. The library's own operations have budgets of
their own, which are off by default: called without one, ``compare()``, ``routes()`` and
the exports (``to_fasta()``, ``to_gfa()``, ``to_json()``) run locally with no work or
allocation budget, their time and peak memory follow the size of the graphlet (a
comparison reads both DAGs up to its depth, an export spells every chosen walk), and
nothing interrupts them. The exports and the comparison spell their walks as a stream --
each walk built from the prefixes it shares with the walks still to come and dropped
once its record or key is made -- so an export's peak is about twice its text (its
records and the joined text), not every spelling at once. With ``budget=`` every local
operation either completes or stops and says so (`Local budgets`_). A tool's
``max_bytes`` bounds the bytes it returns, not the computation behind them.

.. note::

   ``greedy_losses`` is an arm-level limitation: it is stated on the arm whose heads
   re-minimised their label sets, since every work counter belongs to the arm whose head
   did the work.

Every local answer carries the evidence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An answer read on its own must not look more certain than the retrieval it came from.
``metagraph.traverse.ops.evidence_block(g, arm)`` is the block the tool functions attach
to every result: ``complete_to_bp`` and ``scope`` per arm, ``support`` and
``reconverge``, the four outcome dimensions, the kinds of the limitations that apply
(seed level and per arm; an informational ``scope`` is left out) and ``exact``, which
holds for the arm only when its label evidence is complete *and* no recorded list on it
was cut -- a limitation of the other arm does not make it inexact. In the
query results, ``Walk.complete``, ``Walk.beyond_certified_bp`` and ``Claim.exact`` carry
the same information per row.

.. graphlet-example: evidence

.. code-block:: python

   import pprint
   from metagraph.traverse import ops

   pprint.pprint(ops.evidence_block(caps, 'right'), width=88, sort_dicts=False)

.. code-block:: text

   {'complete_to_bp': {'right': 126},
    'exact': True,
    'support': 'kmer',
    'reconverge': 'keep',
    'scope': {'right': 'per_path'},
    'outcome': {'walks': 'partial',
                'branch_diagnostics': 'cut',
                'label_evidence': 'complete',
                'delivery': 'inline'},
    'limitations': {'right': ['walk_domain', 'branch_events']}}

So: an answer is *complete* only when every dimension it depends on reads ``complete``
(walks for which walks exist, label_evidence for who carries them, branch_diagnostics
for why something is absent). Otherwise read the limitations of its class: they say
how far it holds and which knob would take it further.


Queries
-------

Walks
^^^^^

``walks(arm, *, top=None, by='support', labels=None, route_consistent=True, min_bp=0)``
returns ``Walk`` objects: ``path_id``, ``leaf``, ``segments``, ``length_bp``,
``sequence`` (natural orientation, ``None`` without bases), ``path_reason`` (``None``
for a semantic end), ``end_reasons`` (per-label end reasons at the leaf), ``claims`` (the
claims ending at the leaf), ``labels_full`` (labels with displayed support over the
whole walk), ``n_alive`` (labels alive at its end), ``complete`` and
``beyond_certified_bp``.

Rankings: ``'support'`` by (most ``labels_full``, most ``n_alive``, longest, id),
``'length'``, ``'loss'`` (lowest end-label loss first) and ``'id'``. ``labels=`` keeps
the walks that carry one of the selected labels; with ``route_consistent=True`` (the
default) the label must support the *displayed* walk, i.e. not have joined through a
non-first merge parent; ``False`` accepts any label alive at the walk's end.

.. graphlet-example: walks

.. code-block:: python

   for w in merge.walks('right'):
       print(w.path_id, w.length_bp, w.path_reason, w.n_alive,
             [label.name for label in w.labels_full], w.complete)
   print(merge.walks('right', labels=['b.fa']))
   print([w.path_id for w in merge.walks('right', labels=['b.fa'], route_consistent=False)])
   for by in ('support', 'length'):
       print(by, [(w.path_id, w.length_bp) for w in caps.walks('right', by=by)])
   print([w.path_id for w in caps.walks('right', by='length', min_bp=100, top=1)])

.. code-block:: text

   0 100 max_extension_bp 4 ['a.fa', 'both.fa'] True
   []
   [0]
   support [(1, 109), (2, 95), (0, 126)]
   length [(0, 126), (1, 109), (2, 95)]
   [0]

All four labels are alive at the end of ``merge``'s single walk, but only ``a.fa`` and
``both.fa`` spell it: ``b.fa`` and ``c.fa`` took the other allele of the second bubble
and rejoined at 62.

Claims
^^^^^^

``claims(arm=None, labels=None, at_most_bp=None, strict=True)`` returns the unit the
verification contract of the spec (§6.9) is stated in: **one claim per run end**, not
one per leaf. A label that ends inside a walk that other labels continue is a claim at
its own end position; comparing leaves would never see it.

A claim distinguishes two kinds of support:

* **route support** ``[from_bp, to_bp)``: the label's own lineage carries every node of
  *some* route of that length (its route may differ from the displayed walk behind a
  merge);
* **displayed support** ``[evidence_from, to_bp)``: the label carries the bases of the
  displayed walk. ``evidence_from`` is derived at the run's anchored end through the
  merge partitions; ``route_bp`` is the merge-derived part (0 when the label's route is
  the displayed one).

``kind`` is one of ``alive`` (still alive where the walk stopped at the radius, a cap
or a beam), ``end`` (the label ended), ``diverged`` (it went on under another name after a
switch), ``merged`` (closed by a merge: not an end, the lineage continues in the kept
entry), ``boundary`` (present at the boundary node only, no bases), and, for claims cut
at ``at_most_bp``, ``stretch`` and ``route_only``. ``end_class`` groups the reasons:
``lost``, ``dead_end``, ``radius``, ``capped``, ``blocked``, ``pruned``, and ``open``
(not ended: merged, or cut). ``reason`` and ``qualifier`` are the walker's, ``loss`` and
``branches`` the lineage's terminal values.

.. graphlet-example: claims

.. code-block:: python

   for c in merge.claims('right'):
       print(c.label.name, c.kind, (c.from_bp, c.to_bp), 'route_bp', c.route_bp,
             'evidence_from', c.evidence_from, c.end_class, c.branches)

.. code-block:: text

   a.fa alive (0, 100) route_bp 0 evidence_from 0 radius 0
   b.fa alive (0, 100) route_bp 62 evidence_from 62 radius 0
   c.fa alive (0, 100) route_bp 62 evidence_from 62 radius 0
   both.fa alive (0, 100) route_bp 0 evidence_from 0 radius 2
   both.fa merged (0, 36) route_bp 0 evidence_from 0 open 1
   both.fa merged (0, 62) route_bp 0 evidence_from 0 open 2

``b.fa`` and ``c.fa`` have route support over the whole 100 bp but displayed support
from 62 only. ``both.fa``'s lineage branched at both bubbles (``branches`` 2); the
branches on the second parents were closed at the merges at 36 and 62.

**Claims at a depth.** ``at_most_bp=D`` restricts every claim to ``[0, D)``, as a
comparison with a shallower retrieval must. A run that starts at or after ``D`` makes no
claim; a run that reaches beyond ``D`` is cut there (``end_class`` ``open``, ``loss`` and
``branches`` unknown, since they hold at the run's end only) and is a ``stretch`` when
some displayed support survives the cut, ``route_only`` when none does. The hand-made
``clipped_merge`` fixture is the case that motivated this: ``lineB`` joins the displayed
path through the second parent of a merge at 14, and its own route begins ``TCTA``,
not the displayed ``AACC``.

.. graphlet-example: claims-cut

.. code-block:: python

   clip = fixture('clipped_merge')
   print(clip.spell('right', 0), clip.routes('lineB', 'right', spell=True)[0][1])
   for depth in (None, 4):
       for c in clip.claims('right', at_most_bp=depth):
           print('at_most_bp', depth, c.label.name, c.kind, (c.from_bp, c.to_bp),
                 'evidence_from', c.evidence_from, c.end_class, 'loss', c.loss)

.. code-block:: text

   AACCGTAGTCATGCTTGAC TCTAAGGCACATGCTTGAC
   at_most_bp None lineA alive (0, 19) evidence_from 0 radius loss 0.0
   at_most_bp None lineB alive (0, 19) evidence_from 14 radius loss 0.0
   at_most_bp 4 lineA stretch (0, 4) evidence_from 0 open loss None
   at_most_bp 4 lineB route_only (0, 4) evidence_from None open loss None

**Ends inside walks, boundary claims and refusals.** In ``linear`` (quorum
``min_successor_labels: 2``) the minority label at each fork was refused and ended
``branch``/``minority``: on the left at 1 bp, inside the walk the two other labels carry
on, and on the right at the seed boundary itself, a zero-length ``boundary`` claim. The
refused successor is recorded in the arm's branch events (the only record of a
successor not taken).

.. graphlet-example: claims-linear

.. code-block:: python

   lin = fixture('linear')
   for side in ('left', 'right'):
       for c in lin.claims(side):
           print(side, c.label.name, c.kind, (c.from_bp, c.to_bp), c.end_class,
                 c.reason, c.qualifier)
   for be in lin.arm('right').branch_events:
       print(be.at_bp, be.chars, be.labels_per_successor,
             [(r.char, r.cause, [lin.labels[l].name for l in r.labels]) for r in be.refused])

.. code-block:: text

   left acc1 end (0, 40) dead_end dead_end None
   left acc2 end (0, 40) dead_end dead_end None
   left acc3 end (0, 1) pruned branch minority
   right acc1 end (0, 40) dead_end dead_end None
   right acc2 boundary (0, 0) pruned branch minority
   right acc3 end (0, 40) dead_end dead_end None
   0 AT [1, 2] [('A', 'minority', ['acc2'])]

**Switches.** Under a finite change cost a lineage can continue under another label.
In ``switch_chain`` ``A`` switches to ``B`` at 30 and ``B`` to ``C`` or ``D`` at the
split at 60; the switched-away runs are ``diverged``, the switched-in ones carry
``entered_by``, ``from_label`` and ``cost``, and ``loss`` accumulates along the
lineage. A path with switches is a candidate supported by a *chain* of labels, never
evidence that one sample contains it.

.. graphlet-example: claims-switch

.. code-block:: python

   sw = fixture('switch_chain')
   for c in sw.claims('right'):
       print(c.label.name, c.kind, (c.from_bp, c.to_bp), c.entered_by,
             c.from_label.name if c.from_label else None, c.cost, 'loss', c.loss)

.. code-block:: text

   A diverged (0, 30) seed None None loss 0.0
   B diverged (30, 60) switch A 1.0 loss 1.0
   D alive (60, 80) switch B 1.0 loss 2.0
   C alive (60, 80) switch B 1.0 loss 2.0

In ``annotate`` mode there are no lineages: ``claims()`` returns the oracle of the spec,
per label recorded at the seed boundary its maximal label-consistent stretches, filtered
from the recorded sets alone. If a recorded list was cut, those are lower bounds and
``strict=True`` (the default) raises ``IncompleteRecording``; pass ``strict=False`` to
accept them (each claim then has ``exact`` false).

Label walks and routes
^^^^^^^^^^^^^^^^^^^^^^

``label_walks(label, arm=None)`` lists, per run of a label, ``(from_bp, to_bp)``,
``evidence_from``, the bases of the label's **own** route over the run, the end (the
reason as the label-end event states it, ``'switched'`` or ``'merged'``), the leaves
below its end and, for a merged run, the segment it merged into. ``routes(label, arm,
spell=False)`` gives the segment chain of each run's own route, chosen at every merge
through the partition that holds the lineage.

.. graphlet-example: label-walks

.. code-block:: python

   (lw,) = merge.label_walks('b.fa')
   print(lw.run, (lw.from_bp, lw.to_bp), lw.evidence_from, lw.end, lw.leaves_below)
   shown = merge.spell('right', 0)
   print([(i, lw.sequence[i], shown[i]) for i in range(100) if lw.sequence[i] != shown[i]])
   print(merge.routes('b.fa', 'right'), merge.routes('a.fa', 'right'))
   for x in merge.label_walks('both.fa'):
       print(x.run, (x.from_bp, x.to_bp), x.end, 'merged_into', x.merged_into)

.. code-block:: text

   1 (0, 100) 62 max_extension_bp [6]
   [(20, 'T', 'A'), (46, 'G', 'C')]
   [[0, 2, 3, 5, 6]] [[0, 1, 3, 4, 6]]
   3 (0, 100) max_extension_bp merged_into None
   4 (0, 36) merged merged_into 3
   5 (0, 62) merged merged_into 6

``b.fa``'s own route differs from the displayed walk at the two SNPs; its route runs
through segments 2 and 5, the second alleles.

In ``annotate`` mode a label's route may pass through any parent of a merge (that is how
``label_summary()`` measures ``direct_bp``), and ``claims()``, ``label_walks()`` and
``routes()`` follow it there too: one claim per **maximal end** of the label's routes
under the union rule, each shown by **one witness route**, chosen deterministically by
the stored parent order (at every merge, the first parent that records the label). The
witness does not enumerate the other routes that reach the same end, and it does not
establish a contiguous source occurrence: the recorded sets say that every node of the
route carries the label, not that one indexed sequence spells it. For ``b.fa`` in the
``annotate`` fixture (the same locus, label-free) the result is its stretch to the radius
through the second alleles:

.. graphlet-example: annotate-routes

.. code-block:: python

   ann = fixture('annotate')
   for lw in ann.label_walks('b.fa', 'right'):
       print(lw.arm, (lw.from_bp, lw.to_bp), lw.end, lw.leaves_below)
   print(ann.routes('b.fa', 'right'))

.. code-block:: text

   right (0, 70) max_extension_bp [6]
   [[0, 2, 3, 5, 6]]

Support along a walk
^^^^^^^^^^^^^^^^^^^^

``support_profile(arm, walk, kind='displayed')`` returns maximal ``SupportRun`` pieces
(``from_bp``, ``to_bp``, ``labels``, ``exact``, ``total``) of the labels supporting the
displayed bases; ``kind='route'`` adds every run on the walk over its own route
interval. ``support_changes(arm, walk)`` says where the displayed support changes and
why: an end (with the walker's reason), a switch, a ``split`` (the label took another
branch, with its first base), a ``merge`` (it joined through another parent). In
``annotate`` mode both kinds are the recorded sets, with ``exact`` false where a list
was cut and ``total`` the true count.

.. graphlet-example: support

.. code-block:: python

   for r in merge.support_profile('right', 0):
       print((r.from_bp, r.to_bp), [label.name for label in r.labels])
   print([(r.from_bp, r.to_bp, len(r.labels))
          for r in merge.support_profile('right', 0, kind='route')])
   for ch in merge.support_changes('right', 0):
       for r in ch.reasons:
           detail = r.get('took') or r.get('via_parent') or r.get('to') or ''
           print(ch.at_bp, r['change'], r['label']['name'], r['why'], detail)

.. code-block:: text

   (0, 20) ['a.fa', 'b.fa', 'c.fa', 'both.fa']
   (20, 36) ['a.fa', 'c.fa', 'both.fa']
   (36, 46) ['a.fa', 'b.fa', 'c.fa', 'both.fa']
   (46, 62) ['a.fa', 'both.fa']
   (62, 100) ['a.fa', 'b.fa', 'c.fa', 'both.fa']
   [(0, 100, 4)]
   20 removed b.fa split ['T']
   36 added b.fa merge 2
   46 removed b.fa split ['G']
   46 removed c.fa split ['G']
   62 added b.fa merge 5
   62 added c.fa merge 5

Label summary
^^^^^^^^^^^^^

``label_summary()`` returns, per ``ref``, ``{name, ref, <arm>: {direct_bp, reach_bp,
reentries, runs}}`` exactly as the walker computes it. ``direct_bp`` is the longest
stretch from the seed boundary that the label carries along its **own route** (a merge
does not clamp it); ``reach_bp`` the end of its lineage including switches (credited to
the lineage's first label); ``reentries`` the switches into it. In ``annotate`` mode
``direct_bp`` is the continuous presence from the seed boundary along *some* recorded
route and ``reach_bp`` the furthest position the label was recorded at; there are no
runs.

``direct_bp`` is label-consistent **route** support. It is not support for the displayed
walk (that is ``evidence_from``), and it is not a contiguous occurrence in the label's
sequence: a label may hold several contigs or repeat copies, and only ``support:
trace`` or validation against the source record establishes contiguity.

.. graphlet-example: label-summary

.. code-block:: python

   for ref, row in merge.label_summary().items():
       print(ref, row['name'], row['right'])
   ann = fixture('annotate')               # lists cut at 2 labels per node
   for ref, row in ann.label_summary().items():
       print(ref, row['name'], 'direct_bp', row['right']['direct_bp'],
             'reach_bp', row['right']['reach_bp'])

.. code-block:: text

   c:0 a.fa {'direct_bp': 100, 'reach_bp': 100, 'reentries': 0, 'runs': [0]}
   c:1 b.fa {'direct_bp': 100, 'reach_bp': 100, 'reentries': 0, 'runs': [1]}
   c:2 c.fa {'direct_bp': 100, 'reach_bp': 100, 'reentries': 0, 'runs': [2]}
   c:3 both.fa {'direct_bp': 100, 'reach_bp': 100, 'reentries': 0, 'runs': [3, 4, 5]}
   c:0 a.fa direct_bp 70 reach_bp 70
   c:1 b.fa direct_bp 70 reach_bp 70
   c:2 c.fa direct_bp 0 reach_bp 61
   c:3 both.fa direct_bp 0 reach_bp 61

In the ``annotate`` retrieval of the same locus ``c.fa`` and ``both.fa`` read
``direct_bp`` 0: all four labels are at the seed boundary, but its recorded list was cut
to two. This is what ``label_evidence: lower_bound`` warns about.

Splits and successors not taken
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Splits are derived from the segments; ``metagraph.traverse.derive`` exposes them per
arm. A split is ``ambiguous`` when some lineage continued on more than one branch
(it counts against ``max_label_branches``) and a ``divergence`` when the branches'
sources are disjoint. A branch's label count is not a share of ``labels_before``: one
label may follow several branches. Successors the walker did **not** take, and why, are
in ``arm.branch_events`` (see the ``linear`` example above) and in the ``blocked`` and
``hairpin`` events of the segments.

.. graphlet-example: splits

.. code-block:: python

   from metagraph.traverse import derive

   for sp in derive.splits(right, merge.mode):
       branches = derive.split_branches(right, sp, merge.cap)
       print(sp.at_bp, 'segment', sp.segment, 'ambiguous' if sp.ambiguous else 'divergence',
             'labels_before', sp.labels_before,
             [(b['char'], b['labels_distinct']) for b in branches])

.. code-block:: text

   20 segment 0 ambiguous labels_before 4 [('A', 3), ('T', 2)]
   46 segment 3 ambiguous labels_before 4 [('C', 2), ('G', 3)]

Views
^^^^^

``subgraph(selectors, arm=None, mode='any'|'all')`` returns a ``GraphletView``: the
segments and walks that carry the selected labels, with the backing graphlet's original
ids (nothing is renumbered). A view needs at least one label (``None`` or ``[]`` is a
``BadSelector``). It answers ``walks()``, ``claims()``, ``to_fasta()`` and
``summary()`` for the selected labels -- ``walks(top=k)`` ranks the view's walks and
then keeps ``k`` of them -- and its completeness is the backing ``complete_to_bp``
*qualified "for the selected labels"*. Its ``of`` names the backing body (the digest of
its text, or the ``of=`` the caller gives, such as a store handle: then the body is not
dumped for its digest). ``save()`` writes the unchanged backing body plus the view's
selectors in the ``J`` line; ``load()`` returns the view again, and refuses a saved view
that names no label (``GraphletFormatError``) rather than reading it as every label or
none.

.. graphlet-example: views

.. code-block:: python

   v = fork.subgraph(['acc2'])
   print(v.segments, v.path_ids)
   print([w.path_id for w in v.walks('right')], {c.label.name for c in v.claims()})
   print(v.completeness()['right'], v.to_fasta().count('>'))

.. code-block:: text

   {'left': [0, 1], 'right': [0, 1]} {'left': [0], 'right': [0]}
   [0] {'acc2'}
   {'complete_to_bp': 10, 'qualified': 'for the selected labels', 'scope': 'united_history'} 2

Continuations and the next request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A walk that stopped at the radius or at a cap has a **continuation**: the last bases
before its end with the labels that cover all of them, valid ``/traverse`` input.
``continuation(arm, walk)`` returns ``sequence``, ``labels``, ``loss_used``,
``branches_used``, ``seed_coord`` (the half-open interval in seed coordinates) and
``as_seed()``. On the right arm the continuation is the *last* ``n`` bases of
``seed + natural(right)``, on the left arm the *first* ``n`` bases of
``natural(left) + seed``; it may reach into the seed when the radius is short (``fork``
walked 10 bp with k = 15). A walk with a semantic end has none (``ValueError``).

.. graphlet-example: continuation

.. code-block:: python

   c = fork.continuation('right', 0)
   print(c.sequence, [label.name for label in c.labels], c.seed_coord, c.loss_used)
   print(c.as_seed())
   left = fork.continuation('left', 1)
   print(left.sequence, left.seed_coord)

.. code-block:: text

   CTATTATTGCTAAGA ['acc2'] (35, 50) 0.0
   {'sequence': 'CTATTATTGCTAAGA', 'labels': ['acc2']}
   GTCGGACGGGCGTAC (-10, 5)

``loss_used`` is the *smallest* terminal loss among the continuation's labels;
``losses`` holds each label's own, and ``note`` says when one loss budget cannot serve
them exactly.

``next_request(arm, walks, bp=None, reduce_budget=True, reset_branches=False, **overrides)``
builds the resubmittable request: one seed per walk (its continuation, with its labels
named explicitly in ``constrain`` mode), the retrieval's normalized strategy with
``direction`` set to the arm, ``bounds.max_extension_bp = bp`` when given, and the
keyword overrides deep-merged into the strategy (``release``, ``graph`` and
``graph_path`` go to the request level). In ``constrain`` mode three rules hold:

* **No route exceeds its original loss budget.** The loss budget is reduced by the
  *largest* terminal loss of the continued labels. A request carries one loss budget
  and no per-label starting loss, so this is exact for the labels at that loss and
  conservative for the others: their continuation may stop earlier than one
  uninterrupted walk would. The returned request (a ``NextRequest``, a ``dict``) says so
  in ``notes`` and, where a loss budget applies, lists each label's terminal loss and
  remaining budget in ``loss_budget``; neither is sent to the server.
  ``reduce_budget=False`` keeps the original budget, and a note states that a route may
  then exceed it.
* **No lineage branches beyond its original allowance.** The server resets branch state
  for a new seed, so ``branching.max_label_branches`` is reduced the same way, by the
  *largest* terminal branch count of the continued labels (``"unlimited"`` stays
  unlimited): exact for the lineages at that count, conservative for those that used
  fewer. ``branch_budget`` (``{original, effective, largest_terminal_branches, reset}``)
  says how it was derived, and a note states where the reduction is conservative, for the
  lineages that used fewer branches (none when it is exact for every continued lineage,
  as for the loss budget). ``reset_branches=True`` keeps
  the original allowance, and a note states that the continuation may then branch further
  than one uninterrupted walk would.
* **The request is valid.** ``labels.extra`` is rebuilt around the new seed labels: the
  retrieval's permitted labels minus the seeds, each kept that a chain of switches from a
  seed label enters within the budget (what the server accepts: the walk enforces the
  cumulative loss switch by switch). Every label left out is listed in ``left_out``, per
  continued walk (``walk``), and the notes name the first few, with each that one
  uninterrupted walk could still have entered from a continued label at a lower loss (by a
  chain through any label of the retrieval, the kept ones included). An override section the
  rebuild reads that is not an object (``labels``, ``labels.change_cost``,
  ``branching``), or a malformed field of it, is refused with a ``ValueError`` naming it.
  The merged ``labels.change_cost`` is checked whole before any of it is used, by the
  server's rule: ``model`` a string; ``entries``, when present, a list under every model;
  a constant's ``value`` and a table's ``default`` a finite number >= 0 (the default may
  also be ``"forbid"``); a table's ``entries`` required, every entry ``[from, to, cost]``
  with two label names (strings) and a finite cost >= 0 -- the error names the entry
  (``labels.change_cost.entries[3]``). A ``labels.loss_budget`` override is a finite
  number >= 0, a ``branching.max_label_branches`` override an integer >= 0 or
  ``"unlimited"`` (not infinity or NaN). The MCP tool ``traverse_continue`` answers these
  with ``bad_argument``, also with ``execute=False`` and under local limits.

A label alive at a walk's leaf that does not cover the continuation's whole tail (it
switched in within the last bases) is not among the continuation's labels, so the request
does not seed it: the continuation can enter it again only by a switch from a seed label,
and may lack the lineage that one uninterrupted walk keeps. Each such label is named in
``notes`` and listed in ``left_out`` with ``why: 'alive_not_seeded'``, its walk and its
loss.

Walks share one strategy, so one ``labels.extra``: several walks are built into one
request only when that list is the same for all of them; otherwise
``IncompatibleContinuations`` is raised, and ``next_requests()`` builds one request per
walk. ``TraverseClient.deepen()`` sends the request and returns its ``notes`` with the
response. The result is a new traversal, certified on its own; the backend keeps no
frontier between requests (which is why the branch allowance is reduced, above).

.. graphlet-example: next-request

.. code-block:: python

   req = sw.next_request('right', [1], bp=500, bounds={'max_steps': 50000})
   print(req['seeds'])
   st = req['strategy']
   print(st['direction'], st['labels']['loss_budget'], st['bounds']['max_extension_bp'],
         st['bounds']['max_steps'], st['output']['detail'])

.. code-block:: text

   [{'sequence': 'TTACTCGTAGCCGGGCGTGA', 'labels': ['C']}]
   right 1.0 500 50000 graphlet

The request names labels, and ``/traverse`` resolves a name to a label. Where the
graphlet cannot verify that a name resolves back to the same label -- for instance when
two labels of the retrieval share it -- ``next_request()`` refuses to build the request
rather than constrain the new traversal to the wrong label:

.. graphlet-example: next-request-ambiguous

.. code-block:: python

   twin = fixture('fork')
   twin.labels[2].name = 'acc1'          # acc3 renamed: two labels now share a name
   try:
       print('built', twin.next_request('right', [1])['seeds'])
   except (LookupError, ValueError):
       print('refused')

.. code-block:: text

   refused

Comparing retrievals
^^^^^^^^^^^^^^^^^^^^

``compare(other, *, arm=None, labels=None, mode='claims')`` compares two retrievals of
the same locus, keyed by label ``ref`` and restricted to the smaller
``complete_to_bp`` of the two (``depth_used``); a claim reaching that depth is
``open`` on both sides. Both DAGs are **restricted to the depth** before anything is
keyed: the walks are the restricted DAG's leaves, and a claim reaching past the depth is
anchored on its own route there -- what a retrieval walked to that depth shows. A merge
beyond the depth therefore decides nothing before it: the branch that enters it through
a non-first parent, which no displayed walk of the deeper retrieval passes, is a walk of
its own at the depth. Everything at or beyond the depth is outside the restricted DAG,
a merge exactly at the depth included: a retrieval walked to a merge position records
that merge as a zero-length segment, and the runs it closes are open claims at the depth,
not merged ones. Where a route cannot be reconstructed, the comparison is
``'qualified'`` rather than a difference. It returns a ``Comparison`` with ``comparable``, ``reason``,
``depth_used``, ``equal``, ``only_in_a``, ``only_in_b``, ``differ`` and ``notes``.
Modes: ``claims`` (run ends, keyed by label, displayed start and prefix), ``walks``
(walks cut at the depth, with the labels supporting the whole cut prefix), ``labels``
(``label_summary`` cut at the depth) and ``prefix_subset`` (every walk of ``a`` under a
label is a prefix of a walk of ``b`` under it; ``b``'s walks that ``a`` omits are listed
with the reason ``a`` recorded for them).

``comparable`` is graded, and a comparison only moves to a weaker verdict:

* ``True``: the same index (equal ``index_fp``, and equal ``index_ns`` and release where
  set; names agreeing per ref) and the same **oriented seed sequence** (not
  ``validated_seed_id``, which differs between modes). Only then can ``equal`` be
  ``True``.
* ``'qualified'``: the scopes differ (``per_path`` against ``united_history``), or
  either side's label evidence is a lower bound or qualified (a cut list cannot tell
  absent from cut).
* ``'unverifiable'``: either side has no ``index_fp`` (the server loaded no index
  manifest), so label joins cannot be verified; ``index_meta_fp`` can only prove two
  indexes different.
* ``'unknown'``: nothing can be compared, e.g. no common certified depth.
* ``False``: different indexes, seeds or arms.

``equal`` is ``None`` whenever ``comparable`` is not ``True``; a note then says whether a
difference was found. The constrained trie over ``{acc1, acc3}`` equals the
label-free oracle filtered to those labels, and without the filter the oracle's
``acc2`` claim is reported:

.. graphlet-example: compare

.. code-block:: python

   def body_of(name):
       with open('compare/%s.mgt' % name) as f:
           return parse(f.read())

   constrained, oracle = body_of('oracle_constrain'), body_of('oracle_annotate')
   c = constrained.compare(oracle, labels=['acc1', 'acc3'])
   print(c.comparable, c.equal, c.depth_used, c.reason)
   c = constrained.compare(oracle)
   print(c.equal, [(x['label']['name'], x['prefix'], x['to_bp'], x['end_class'])
                   for x in c.only_in_b])

.. code-block:: text

   True True 20 same index and seed
   False [('acc2', 'TCCGGAAT', 8, 'dead_end')]

A tuned retrieval (quorum 2) against the exhaustive one: the tuned walks are a prefix
subset, and the walk it omits carries the refusal the tuned body recorded.

.. graphlet-example: compare-prefix

.. code-block:: python

   tuned, exhaustive = body_of('oracle_tuned'), body_of('oracle_exhaustive')
   c = tuned.compare(exhaustive, mode='prefix_subset')
   print(c.comparable, c.equal, c.only_in_a)
   print([(x['name'], x['walk'], x['at_bp'], x['reason']) for x in c.only_in_b])

.. code-block:: text

   True True []
   [('acc2', 'TCCGGAAT', 0, 'minority')]

The CLI-generated fixtures were made without an index manifest, so joins across them are
unverifiable: the differences are still computed and listed, but equality is not
claimed. A comparison with a lower-bound side or across scopes collects the
qualifications in ``notes``; a pair without a common certified depth is ``unknown``.

.. graphlet-example: compare-verdicts

.. code-block:: python

   qt, qe = body_of('cli_quorum_tuned'), body_of('cli_quorum_exhaustive')
   print(qt.index_ns, qt.index_fp, qt.index_meta_fp)
   c = qt.compare(qe, mode='prefix_subset')
   print(c.comparable, c.equal, [(x['name'], x['reason'], x['at_bp']) for x in c.only_in_b])
   c = merge.compare(ann, arm='right', labels=['a.fa'])
   print(c.comparable, c.equal)
   for note in c.notes:
       print('   ', note)
   nb = fixture('no_bins')                 # complete_to_bp 0 on both arms
   c = nb.compare(nb)
   print(c.comparable, c.equal, c.depth_used)

.. code-block:: text

   fixtures None 178cc93c330374f0
   unverifiable None [('x.fa', 'minority', 10), ('y.fa', 'minority', 40)]
   unverifiable None
       b: label_lists (labels.max_labels_per_node): the oracle is a lower bound
       b: inexact_counts (labels.max_labels_per_node): the oracle is a lower bound
       right arm of a: histories were united at merges
       no difference found, but equality is not claimed (unverifiable)
   unknown None 0

A retrieval made with ``output.sequences: false`` carries no bases, so walks and claim
prefixes cannot be matched: such a comparison is ``unknown``:

.. graphlet-example: compare-no-bases

.. code-block:: python

   for mode in ('claims', 'walks', 'prefix_subset'):
       c = caps.compare(fixture('caps'), mode=mode)
       print(mode, c.comparable, c.equal)

.. code-block:: text

   claims unknown None
   walks unknown None
   prefix_subset unknown None

Tables and size
^^^^^^^^^^^^^^^

``metagraph.traverse.frames.frames(g)`` returns pandas DataFrames (labels, segments,
runs, walks, splits, claims, events, presence; label columns hold refs with names
beside them); ``records(g)`` returns the same rows as plain dicts without pandas.
``memory_bytes()`` estimates the parsed model's footprint.

.. graphlet-example: frames

.. code-block:: python

   from metagraph.traverse.frames import records

   print(' '.join('%s=%d' % (table, len(rows)) for table, rows in records(merge).items()))
   print(merge.memory_bytes() > 0)

.. code-block:: text

   labels=4 segments=7 runs=6 walks=1 splits=2 claims=6 events=0 presence=0
   True


Exports
-------

FASTA
^^^^^

``to_fasta(arm=None, leaves=None, with_seed=True, orientation='natural', width=None)``
writes one record per walk, ``>{arm}_{path_id}`` with its length, reason, number of
labels and whether the seed is included.

.. graphlet-example: fasta

.. code-block:: python

   print(fork.to_fasta('right'), end='')
   print(merge.to_fasta(with_seed=False, width=60), end='')

.. code-block:: text

   >right_0 length_bp=10 reason=max_extension_bp labels=1 seed=yes
   CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATTATTGCTAAGA
   >right_1 length_bp=10 reason=max_extension_bp labels=2 seed=yes
   CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATTTTACACCACG
   >right_0 length_bp=100 reason=max_extension_bp labels=4 seed=no
   TTTCCTATTTAGCCTCTGTCATTACGTTTGACAATGACCCAGCCCTCCGGCGGGTCGACT
   TGGTCCGGACGATAGCACTTAGTTCCTCACTTCACAATAG

GFA
^^^

``to_gfa(with_seed=True)`` writes GFA 1.0 on the seed strand: an ``S`` line per segment
(and one for the seed) carrying k-1 bases of context, so that every ``L`` line overlaps
by ``(k-1)M``; a ``P`` line per walk with ``LB`` (the refs of the labels alive at its
end) and ``ER`` (its end reason) tags. Merges appear as nodes with two incoming links.
GFA cannot carry the per-node label runs or the merge partitions: it is an export, the
graphlet stays the record. (Tabs are shown as two spaces below.)

.. graphlet-example: gfa

.. code-block:: python

   print(fork.to_gfa().replace('\t', '  '), end='')

.. code-block:: text

   H  VN:Z:1.0
   S  seed  CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATT  LN:i:40
   S  l0  GCGTACCGTCGTAGC  LN:i:15
   S  l1  ATCCCTCACGCGTACCGTCGTAG  LN:i:23
   S  l2  GTCGGACGGGCGTACCGTCGTAG  LN:i:23
   S  r0  TTGCAGGTTCTATT  LN:i:14
   S  r1  TTGCAGGTTCTATTATTGCTAAGA  LN:i:24
   S  r2  TTGCAGGTTCTATTTTACACCACG  LN:i:24
   L  l0  +  seed  +  14M
   L  l1  +  l0  +  14M
   L  l2  +  l0  +  14M
   L  seed  +  r0  +  14M
   L  r0  +  r1  +  14M
   L  r0  +  r2  +  14M
   P  left_0  l1+,l0+,seed+  14M,14M  LB:Z:h:0:0,h:0:1  ER:Z:max_extension_bp
   P  left_1  l2+,l0+,seed+  14M,14M  LB:Z:h:0:2  ER:Z:max_extension_bp
   P  right_0  seed+,r0+,r1+  14M,14M  LB:Z:h:0:1  ER:Z:max_extension_bp
   P  right_1  seed+,r0+,r2+  14M,14M  LB:Z:h:0:0,h:0:2  ER:Z:max_extension_bp

JSON
^^^^

``to_json()`` rebuilds the server's ``detail: "full"`` result for the seed (natural
orientation, every derived field), which makes it the conformance oracle: compared after
``metagraph.traverse.export.normalize_result()`` (event order within one position,
``needed_budgets`` order and ``timing`` are not information), it equals what the server
returns for the same request with ``detail: "full"``. It needs the envelope.

.. graphlet-example: json

.. code-block:: python

   from metagraph.traverse.export import normalize_result

   with open('merge.full.json') as f:
       full = json.load(f)
   mine = merge.to_json()
   print(normalize_result(mine) == normalize_result(full['results'][0]))
   print(mine['arms']['right']['paths'][0]['end_labels'][1])     # b.fa, merged in at 62

.. code-block:: text

   True
   {'label': 1, 'loss': 0.0, 'branches': 0, 'run': 1, 'route_bp': 62}

MGT
^^^

``dump()`` returns the canonical text: the body as the server sent it for a graphlet
read with ``from_response()``; ``dump(envelope=True)`` adds the ``J`` line, which is what
``save()`` writes.

.. graphlet-example: mgt

.. code-block:: python

   text = merge.dump()
   print(text.splitlines()[0])
   print(len(text.encode()), merge.seed_summary['graphlet_bytes'],
         merge.dump(envelope=True).splitlines()[1][:30])

.. code-block:: text

   H mgt 1 15 basic $ACGT c k m 64 40 0 walk fixtures * fbfd396e439f0d22
   1225 1225 J {"algorithm_version":"traver


Local budgets
-------------

Every local operation takes ``budget=`` (keyword-only): a ``LocalBudget`` of **work
units** and a **memory account**, under which the call either completes or stops and
says so -- never a shorter answer that could be read as complete. Without one (the
default) nothing is budgeted and every answer is byte for byte what it was before.

* **Work** is counted in local work units (lwu): a deterministic, weighted count of the
  model elements an operation visits and of the rows and text it makes (work model 1;
  1 lwu is about 0.1 us of CPython 3.11 on the reference machine, so 30 M lwu is about
  3 s nominal). It is not CPU time. A derivation an operation uses (paths, splits, merge
  maps, evidence, label summaries, ...) is charged at its *cold price*, whether a cache
  holds it or not: the same call on the same graphlet with the same budget charges the
  same units and stops at the same row, in any process and whatever earlier calls cached.
  Warm calls pay more than they spend.
* **Memory** is a modelled account (``memory_bound: "model"``): the bytes of the rows,
  spellings, text and derived caches a call builds, each charged before it is made at a
  size modelled on CPython. It is not a measurement: Python can allocate more than the
  model says (interpreter internals, fragmentation, another CPython version, text beyond
  ASCII), so in-process it is a soft bound; a hard bound needs process limits. The
  memory limit holds for each call; the work of the calls that share one budget adds
  up (an allowance). One exception to "charged before": a record of output text is
  charged as soon as it is built (its size is known only then).
* ``LocalLimits(work_units=None, memory_mb=None, deadline_s=None)`` holds the limits
  (``None``: unlimited). ``deadline_s`` (elapsed seconds since the budget was made) and
  ``budget.cancel()`` (from another thread) stop a call at its next charge point; they
  are not reproducible and say so (``deterministic: false``).
* ``local_budget(...)`` sets an ambient budget for the code inside it, which every call
  without ``budget=`` charges; an explicit budget wins. A new thread starts without one.
* ``next_request()`` and ``next_requests()`` take a ``LocalBudget`` as ``budget=`` like
  every other operation; any other value of ``budget=`` is, as before, a keyword
  override like the others, deep-merged into the strategy (so the two together need the
  ambient ``local_budget()``).

A stop raises ``LocalBudgetExceeded`` (an ``Exception``, neither a ``ValueError`` nor a
``RuntimeError``) whose ``.stop`` states it: the resource (``work``, ``memory``,
``deadline``, ``cancelled``), the operation and its phase, the limit, what was used and
at least what would have been needed (in lwu or bytes), how far the call got
(``done``), and the levers (``raise_local_budget``, ``narrow_labels``,
``narrow_arm``, ``select_walks``, ``rank_id``, ``export_mgt``, ``process_locally``,
``retry_later``). A stop leaves only complete caches: an unbudgeted call afterwards
answers as on a fresh model.

* **Lists whose order allows it** -- ``claims()``, ``label_walks()``, ``routes()``,
  ``splits()``, ``walks(by='id')``, ``walks_at()`` -- put the whole rows they made, a
  prefix in their documented order, in ``e.partial.rows``, and a token in
  ``e.partial.resume`` that the same call takes as ``resume=`` to go on after the last
  row. A token of another call, other arguments or another graphlet is refused: it
  binds the digest of the graphlet's body (``Graphlet.body_digest``, made by the parse),
  which the same body parsed again, or saved and loaded, shares.
* **Ranked lists and summaries** -- ranked ``walks()``, ``rank_walks()``,
  ``label_summary()``, ``summary()``, ``support_profile()``, ``support_changes()``,
  ``subgraph()``, ``continuation()``, ``next_request()`` -- raise with no partial: a
  ranking or a summary of part of the arm would be another answer.
* **compare()** never raises for a budget: a stopped comparison returns
  ``comparable: 'unknown'``, ``equal: None``, no difference lists and ``local_stop``,
  because a list from a half-keyed side would be a false difference.
  ``ops.compare_cost(a, b, mode=...)`` estimates the charge beforehand: ``at_least`` (the
  structural phases, certain once both sides are keyed: never above what the completed
  comparison charges, for every mode, ``arm=`` and ``labels=``) and ``estimate``, with
  the phases whose size depends on the answer named in ``unpriced``.
* A budgeted call is admitted before the derived caches it reads are built (the arm's
  paths, indexes, price lists): a call stopped at any charge point -- a zero budget
  included -- leaves behind only derivations it was charged for.
* **Exports and saves** (``to_json()``, ``to_fasta()``, ``to_gfa()``, ``dump()``,
  ``save()``) return no text and write no file when they stop: never a partial export.
* **Parsing** (``parse()``, ``from_response()``, ``load()``): a body whose line count
  alone exceeds the work budget, or whose line list alone exceeds the account, is
  refused before any record is read; a parse that stops is a local failure -- no model,
  the text untouched, never a shallower graphlet.

.. graphlet-example: local-budgets

.. code-block:: python

   from metagraph.traverse import LocalBudget, LocalBudgetExceeded, local_budget

   b = LocalBudget()
   all_claims = merge.claims(budget=b)
   print(len(all_claims), b.usage()['work_units'] > 0, b.usage()['work_model'])
   try:
       merge.claims(budget=LocalBudget(work_units=b.used_work - 1))    # one unit short
   except LocalBudgetExceeded as e:
       print(e.stop.resource, e.stop.op, e.stop.phase, e.stop.deterministic,
             0 < len(e.partial.rows) < len(all_claims))
       rest = merge.claims(resume=e.partial.resume)
       print([c.as_dict() for c in e.partial.rows + rest]
             == [c.as_dict() for c in all_claims])
   whole = LocalBudget()
   merge.compare(merge, budget=whole)
   c = merge.compare(merge, budget=LocalBudget(work_units=whole.used_work // 2))
   print(c.comparable, c.equal, c.local_stop['resource'], c.local_stop['phase'])
   with local_budget(work_units=10 ** 9) as allowance:    # one allowance, two calls
       merge.to_fasta()
       merge.label_summary()
   print(allowance.used_work > 0, allowance.stops)

.. code-block:: text

   6 True 1
   work claims rows True True
   True
   unknown None work keys:a
   True []

Not charged: the store's bookkeeping (``refresh``, ``memory_bytes()``), the HTTP
client's decoding of a response, cursor MACs, ``check_rules()``, the tables of
``frames()`` (the library calls they make are charged), and the MCP framework's
serialisation of a result. The store runs its own parses under ``parse_limits`` (see
below), and the MCP tools take budgets per call (``local_limits``).


Handles, the store and the tool functions
-----------------------------------------

These pieces are for a service that keeps graphlets for an agent across many calls.

The store
^^^^^^^^^

``GraphletStore(spool_dir, max_ram_mb=512, max_handles=64, ttl_ram_s=1800,
ttl_disk_s=7*86400, *, max_body_mb=None, ttl_tomb_s=30*86400, clock=time.time,
parse_limits=None)`` keeps retrievals behind **opaque handles** (``g_`` and 12 random hex
digits):

* a handle names an *entry*: the request, the index identity, the envelope and a
  digest of the body. Bodies are deduplicated by sha256 in the spool, but a body hash is
  never a handle: two requests with different loss budgets can return byte-identical
  bodies that continue differently;
* ``put(response, request, index=0)`` parses the whole body first; an invalid or
  truncated body raises and nothing is stored. A body over ``max_body_mb`` raises
  ``StoreLimitExceeded``: rejected, never truncated. ``put_all()`` stores every seed;
* ``graphlet(handle)`` returns the parsed model from RAM, or parses it again from the
  spool. Parsed models are kept least-recently-used first out under ``max_handles``
  and the RAM budget ``max_ram_mb``; a model larger than the whole budget is never kept
  resident;
* entries expire ``ttl_disk_s`` after their last use (``sweep()`` applies the TTLs);
  an expired entry leaves a tombstone with its request, so ``UnknownHandle.replayable``
  tells a caller that the retrieval can be run again. ``free(handle)`` deletes for good;
* ``save(handle, path)`` / ``load(path)`` move entries through ``.mgt`` files
  (``load_graphlet(path)`` also returns the model it parsed, kept or not), and
  ``list()`` describes the entries. ``request_of(g, graph=..., graph_path=...)`` keeps
  the graph of a multi-graph server in the replay request (else the envelope's, when
  the server echoed it); a request without one is valid on a single-graph server only;
* the disk TTL counts from the last use also across restarts (the entry file is
  rewritten once its recorded use lags by a tenth of the TTL). Processes may share a
  spool: a body is deleted only when no entry file names it any more, a handle another
  process stored is read from its entry file on first use (by ``get``, ``free``, ``in``
  and ``list()``, which lists the spool's entries), one another process freed or
  expired answers as unknown, and an expiry re-reads the entry file's last use first.
  An entry whose body is missing expires (``UnknownHandle``, replayable when it kept a
  request), also when its body is copied without a parse (``body_text``,
  ``standalone_text``, ``save_body``); an entry file of another library version is read
  for the fields this one knows;
* with ``parse_limits=LocalLimits(...)`` every parse the store runs (``put``, a parse on
  demand, ``load``) has a fresh budget of those limits. A body whose parse stops is not
  lost: ``put`` stores it after the checks that need no parse (``graphlet_bytes``,
  ``graphlet_lines``, the ``H`` record, the ``Z`` line count -- a body cut in transport
  is still refused) as an entry marked ``parsed: false``, and raises the stop with
  ``e.handle`` naming it. No answer is derived from such an entry until a parse
  completes (``graphlet(handle, parse_budget=...)``); ``standalone_text(handle)`` and
  ``save_body(handle, path)`` copy any entry out without a parse, byte for byte what
  ``save()`` writes.

.. note::

   The RAM budget is enforced on ``memory_bytes()``: every structure a model retains,
   including the caches its queries add, with shared sets counted once. The store
   charges a model again when its caches grow, so a query can evict other entries.

.. graphlet-example: store

.. code-block:: python

   from metagraph.traverse import GraphletStore, UnknownHandle

   store = GraphletStore(tempfile.mkdtemp(), max_ram_mb=64, max_handles=8)
   with open('merge.graphlet.json') as f:
       stored_response = json.load(f)
   handle = store.put(stored_response, request)
   print(handle[:2], len(handle), store.in_ram(handle))
   print(store.graphlet(handle).arm('right').complete_to_bp,
         store.get(handle).request == request)
   store.free(handle)
   try:
       store.graphlet(handle)
   except UnknownHandle as e:
       print('UnknownHandle, replayable:', e.replayable)

.. code-block:: text

   g_ 14 True
   100 True
   UnknownHandle, replayable: False

The tool functions
^^^^^^^^^^^^^^^^^^

``metagraph.traverse.mcp_tools.GraphletTools(store, clients=None, *,
default_index=None, max_bytes=2048, sequence_max_bytes=16384, secret=None,
export_dir=None, local_limits=None)`` bundles framework-agnostic functions, one per tool,
each returning a JSON-serialisable dict. They are written for an MCP server to register (with whatever
MCP framework it uses); this package does not run a server itself. ``clients`` maps
index names to ``TraverseClient`` objects for the backend tools.

* **Backend tools:** ``traverse_capabilities``, ``traverse_resolve``,
  ``traverse_fetch`` (one retrieval: a handle plus a summary; ``replay=<handle>``
  re-runs a stored request against the same index and refuses a different release or
  index digest) and ``traverse_continue`` (a continuation as a new traversal; the new
  entry remembers its parent walk, and the result carries the request's ``notes``,
  where a loss budget applies ``loss_budget``, and where a branch allowance applies
  ``branch_budget`` (``reset_branches=true`` keeps the original allowance, stated) -- all
  part of its receipt -- and each label's terminal loss and remaining budget in
  ``loss_budget_labels``, the first optional field cut when the receipt would not fit). An index digest on one side only proves neither identity nor
  difference: it is refused as ``index_unverifiable``, and run only when the caller
  passes ``allow_unverified_index=true`` (the result then states the identity
  unverified). A body over ``max_graphlet_mb`` is spooled complete and only its summary
  and handle are returned (``delivery: spooled``).
* **Local tools over a handle:** ``graphlet_summary``, ``graphlet_walks``,
  ``graphlet_walk``, ``graphlet_support``, ``graphlet_labels``, ``graphlet_splits``,
  ``graphlet_claims``, ``graphlet_sequence``, ``graphlet_export``,
  ``graphlet_compare``, ``graphlet_subtrie`` (a derived handle that is a view) and
  ``graphlet_list`` / ``_free`` / ``_save`` / ``_load``.

The contract:

* labels leave as ``{name, ref}``, never as ids; label arguments are tagged selectors or
  a string unique as a ref or a name;
* every local answer carries the ``evidence`` block above;
* every result fits ``max_bytes`` (2 KB by default, 64 bytes to 16 MiB when given;
  ``graphlet_sequence`` has its own
  16 KB ceiling, stated in its result, and ``traverse_capabilities`` a 16 KB default, so
  that the server's description of itself comes whole; an explicit ``max_bytes`` still
  holds). List tools page, a page filled in time linear in its rows: they return
  ``total`` and an
  opaque ``next_cursor`` bound to the handle, the tool and its arguments, which stays
  valid across a restart while the handle does. A single row too large for a page comes
  alone with ``row_truncated`` and the cut fields named; an answer that cannot fit at
  all is the error ``result_too_large``, never a silently shortened one;
* an operation that makes something (a stored handle, an exported or saved file, a
  loaded handle, a view) always returns its receipt, the handle or the path: its other
  fields are cut first and named in ``fields_cut``, and an operation whose receipt could
  not fit is refused before it stores or writes anything (``receipt_too_large``, whose
  message names the lever the tool offers: a larger ``max_bytes``, or a shorter file name
  for an export or a save);
* ``max_bytes`` bounds the bytes returned, not the work. Without ``local_limits`` (the
  default) the local tools (``graphlet_compare``, ``graphlet_export``, the route
  listings) run with no work or allocation budget, and their time and peak memory follow
  the graphlet's size: a service sets ``local_limits`` (`Local budgets in the tools`_);
* a filter never hides rows silently: what it removed is counted (e.g.
  ``filtered: {merge_entered: N}`` with a hint how to include them);
* files are written to and read from the export directory (default
  ``<spool>/exports``) only;
* failures are results ``{error, message, ...}`` -- a malformed argument too (an arm
  that is not a string, a label list that is not a list, a malformed override of a
  continuation) -- with codes such as ``bad_argument``, ``bad_arm``, ``bad_cursor``, ``unknown_handle`` (with ``replayable``),
  ``unknown_label``, ``ambiguous_label``, ``missing_envelope``,
  ``incomplete_recording``, ``not_in_view``, ``view_unsupported``, ``no_bases``,
  ``path_not_allowed``, ``format_error``, ``io_error``, ``too_large``,
  ``result_too_large``, ``receipt_too_large``, ``backend_error`` (with the HTTP status,
  and ``retry_after_s`` while the server loads), ``backend_unreachable``,
  ``not_replayable``, ``seed_failed``, ``index_mismatch``, ``index_unverifiable``,
  ``release_mismatch``, ``unverifiable_label_name``.

.. graphlet-example: tools

.. code-block:: python

   from metagraph.traverse.mcp_tools import GraphletTools, tool_names

   print(textwrap.fill(', '.join(tool_names()), 88))
   tools = GraphletTools(store)            # local tools only: no clients configured
   h = store.put(stored_response, request)
   page = tools.graphlet_walks(h, 'right', spell='tail', tail_bp=24)
   print(page['evidence']['outcome']['walks'], page['evidence']['limitations'])
   row = page['rows'][0]
   print(row['walk'], row['end'], row['n_alive'], row['n_full'], row['sequence'])
   page = tools.graphlet_claims(h, 'right')
   print(page['total'], page['filtered'], [row['label']['name'] for row in page['rows']])
   page = tools.graphlet_claims(h, 'right', route_consistent=False, max_bytes=1400)
   rest = tools.graphlet_claims(h, 'right', route_consistent=False, max_bytes=1400,
                                cursor=page['next_cursor'])
   print(len(page['rows']), len(rest['rows']), 'next_cursor' in rest)
   print(tools.graphlet_claims(h, 'right', cursor=page['next_cursor'])['error'])
   print(tools.graphlet_labels(h, name='ACC9')['error'],
         tools.graphlet_export(h, format='fasta', path='/tmp/x.fa')['error'])

.. code-block:: text

   graphlet_claims, graphlet_compare, graphlet_export, graphlet_free, graphlet_labels,
   graphlet_list, graphlet_load, graphlet_save, graphlet_sequence, graphlet_splits,
   graphlet_subtrie, graphlet_summary, graphlet_support, graphlet_walk, graphlet_walks,
   traverse_capabilities, traverse_continue, traverse_fetch, traverse_resolve
   partial {'right': ['scope']}
   0 max_extension_bp 4 2 ACTTAGTTCCTCACTTCACAATAG
   4 {'merge_entered': 2} ['a.fa', 'both.fa', 'both.fa', 'both.fa']
   3 3 False
   bad_cursor
   unknown_label path_not_allowed


Local budgets in the tools
^^^^^^^^^^^^^^^^^^^^^^^^^^

``GraphletTools(..., local_limits=ToolLimits())`` runs every local tool under a budget
of its class:

==========  ==========================================================  ==================  =====================
class       tools                                                       default             ceiling
==========  ==========================================================  ==================  =====================
view        summary, walks, walk, support, labels, splits, claims,       30 M lwu, 1 GiB     300 M lwu, 10 GiB
            sequence
heavy       compare, export, subtrie, save                              300 M lwu, 2 GiB    3 G lwu, 20 GiB
parse       fetch, continue, load (and a parse on demand)               200 M lwu, 2 GiB    2 G lwu, 20 GiB
==========  ==========================================================  ==================  =====================

``traverse_capabilities``, ``traverse_resolve``, ``graphlet_list`` and ``graphlet_free``
are not budgeted (they compute nothing locally beyond handling JSON).

* An agent may raise or lower its call's budget with ``budget={"work_units": ...,
  "memory_mb": ..., "deadline_s": ...}``, up to the class ceiling; a value above it is
  clamped and the clamp stated (``local.clamped``). ``ToolLimits(per_tool={...})``
  overrides a tool's default, ``budget_for(tool, request)`` lets the service supply the
  budget object itself (its tiers, its ledger's reservation), and ``on_usage(tool,
  usage)`` receives the usage of every budgeted call, completed or stopped. Without
  ``local_limits`` a ``budget`` argument is refused (``bad_argument``).
* Every budgeted result carries ``local``: ``{complete, usage: {work_units,
  memory_bytes}, limits, work_model, memory_bound: "model", clamped?, parsed?, stop?}``
  (``parsed``: the usage of a parse the store ran for the call and the limits it ran
  under, which is not charged to it). A parse on demand runs under the parse class,
  raised field by field by the call's ``budget`` argument up to the parse class's
  ceiling (never lowered by it), so ``raise_local_budget`` is a lever for a stopped
  parse too; under the store's own ``parse_limits``, which no call raises, that lever is
  not offered. An interrupted answer also names the scope local in its evidence
  (``evidence.limitations.local: ["local_work"]``), so it never reads as complete.
* A stop is the error ``local_budget_exceeded`` with its ``stop`` (held to
  ``max_bytes``: the message is cut first, the code and the resource kept), except:
  the lists whose order allows it (``graphlet_claims``, ``graphlet_labels(name=...)``,
  ``graphlet_splits``, ``graphlet_walk``, and ``graphlet_walks`` once its ranking is
  made) answer the whole rows made so far with ``complete: false``,
  ``total_at_least`` and a ``next_cursor`` that resumes after the last row
  (``filtered`` becomes ``filtered_so_far``; ``total`` appears once the list is
  complete); ``graphlet_compare`` answers ``comparable: "unknown"``, ``equal: null``,
  ``counts: null`` and no rows; an export, save or view stopped makes nothing.
* A fetch or continuation whose parse stops keeps the body as an unparsed entry and
  returns its handle with the server's own per-seed summary (``server_summary``) and the
  stop; ``graphlet_export(format="mgt")`` and ``graphlet_save`` copy any entry out
  without a parse, so an agent can always take a graphlet it cannot afford to process
  here and process it with its own resources. The receipt of such a copy of an unparsed
  entry says ``parsed: false`` (and, cut first when it does not fit, ``validation``):
  only the body's frame was checked, and a load of the file validates its records.
* Budgeted results keep room for the ``local`` block: a page, a summary and a
  ``graphlet_sequence`` slice are fitted to the ceiling less that room.

.. graphlet-example: tools-local

.. code-block:: python

   from metagraph.traverse.mcp_tools import ToolLimits

   usage = []
   budgeted = GraphletTools(store, local_limits=ToolLimits(
       on_usage=lambda tool, u: usage.append((tool, u['complete']))))
   page = budgeted.graphlet_claims(h, 'right', route_consistent=False)
   print(page['local']['complete'], page['local']['limits']['work_units'],
         page['total'])
   small = budgeted.graphlet_claims(h, 'right', route_consistent=False,
                                    budget={'work_units': 1})
   print(small['error'], small['stop']['resource'], small['stop']['op'])
   big = budgeted.graphlet_claims(h, 'right', budget={'work_units': 10 ** 12})
   print(big['local']['clamped'], big['local']['limits']['work_units'])
   print(usage)

.. code-block:: text

   True 30000000 6
   local_budget_exceeded work claims
   ['work_units'] 300000000
   [('graphlet_claims', True), ('graphlet_claims', False), ('graphlet_claims', True)]


Pitfalls
--------

**Claims, not leaves.** A label can end inside a walk that other labels continue
(``acc3`` at 1 bp on the left arm of ``linear``). Leaves, FASTA records and walk lists
never show such an end; ``claims()`` does, and comparisons are made on claims for this
reason.

**The left arm is stored outward.** Positions and ``Segment.walk`` are in walking
order on both arms, so on the left arm the natural spelling is the reverse, an outward
prefix is a suffix of the natural string, and the continuation is the *first* ``n``
bases of ``natural(left) + seed``. Use ``spell()`` rather than concatenating
``Segment.walk`` yourself.

**Merged-in labels.** Behind a merge a label alive at a walk's end may have come
through another parent: it has route support over the whole walk but displayed support
only from ``evidence_from`` (``route_bp`` is the merge-derived part). Counting "labels
that share this exact flank" means counting claims with ``evidence_from == 0`` (or
``labels_full``), not labels alive at the leaf. ``walks(labels=...)`` and the tool
functions are route-consistent by default and say what they filtered.

**No bases.** With ``output.sequences: false`` the segments carry no bases:
``spell()``, ``to_fasta()`` and ``to_gfa()`` raise ``ValueError``, ``Walk.sequence`` is
``None``, and a comparison cannot match prefixes (``unknown``, see above). Continuations
still carry their sequence, so the next request can be built.

.. graphlet-example: no-bases

.. code-block:: python

   try:
       caps.spell('right', 0)
   except ValueError as e:
       print(e)
   print(caps.walks('right')[0].sequence, caps.continuation('right', 0).sequence)

.. code-block:: text

   this retrieval carries no bases (output.sequences: false)
   None CCGAGTGATATATCACCTTGCTA

**Cut lists make label evidence a lower bound.** When ``labels.max_labels_per_node``
cut a recorded list (``label_lists``), a label missing from it may still be there: the
``annotate`` example above reads ``direct_bp`` 0 for two labels that are on the seed.
``claims()`` refuses such an arm unless ``strict=False``, and ``compare()`` never calls
such a side equal. Raise the knob (every list keeps its true count, so you can see by
how much) before drawing conclusions from absence.

**Route support is not occurrence.** ``direct_bp``, ``reach_bp`` and claims mean that
a label's k-mers cover a graph route, not that its sequence contains those bases
contiguously; and a beam path, or a path with switches, is a candidate supported by a
chain of labels, not a sequence any one sample is claimed to contain.


Running the examples
--------------------

``metagraph/docs/source/_graphlet_examples/run_examples.py`` executes every example of
this page, in order and in one namespace, from the fixture directory with
``metagraph/api/python`` on ``sys.path``, and compares what each prints with the output
shown::

    python3 metagraph/docs/source/_graphlet_examples/run_examples.py            # check
    python3 metagraph/docs/source/_graphlet_examples/run_examples.py --write    # regenerate

An example whose directive carries the flag ``pending`` shows the intended behaviour
of a library fix that is not merged yet; the runner reports its difference without
failing, and ``--write --pending`` replaces its output once the fix is in. No example is
pending today. The client example runs against a stand-in server that
replays the committed fixture response.
