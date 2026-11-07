===================
Snakemake workflows
===================

Use ``metagraph-workflows`` to build and annotate an index with a single command.
The wrapper runs the graph construction, annotation, and compression stages through Snakemake.


Installation
------------

Install the CLI with Conda and the current workflow wrapper from a shallow
checkout. The same checkout supplies the example data; no C++ build or submodule
initialization is needed:

.. code-block:: bash

   conda create -n metagraph --override-channels -c conda-forge -c bioconda \
       --strict-channel-priority python=3.12 pip "metagraph>=0.5.2"
   conda activate metagraph
   git clone --depth 1 --branch master https://github.com/ratschlab/metagraph.git
   cd metagraph
   python -m pip install ./metagraph/workflows

For custom alphabets or features unavailable in the packaged CLI, build from
source as described in :ref:`installation` and install the wrapper from that checkout.


Creating graphs and annotations
-------------------------------

A single command runs the full pipeline — graph construction,
annotation, and all row-diff / BRWT transforms — for a set of samples::

    metagraph-workflows build samples.txt -o /tmp/mygraph --primary

``samples.txt`` is a text file listing input paths (one per line); a
directory of sample files also works. Process substitution is supported
too, so you can pipe a glob inline::

    metagraph-workflows build <(ls /data/samples/*.fa) -o /tmp/mygraph --primary


The same pipeline can be invoked from a Python script:

.. code-block:: python

    from metagraph_workflows import cli

    cli.run_workflow('/tmp/mygraph', samples='samples.txt', k=31, build_primary_graph=True)


The pipelines are written in the `Snakemake <https://snakemake.readthedocs.io/>`__ workflow management system and can also be directly invoked using the ``snakemake`` command line tool (see below).


Usage
~~~~~

Typically, the following steps would be performed:

1. Prepare a list of input files (or a directory).
2. Construct a MetaGraph index: invoke ``metagraph-workflows build``.
   Tell the workflow how much hardware is available and what kind of
   index you want; the per-stage memory caps, thread packing, and
   BRWT clustering parameters are derived automatically.

   Choose these parameters for your data and hardware:

   * ``-p N`` and ``--mem-gb GB`` for the hardware budget
   * ``-k`` for k-mer length (defaults to the selected MetaGraph executable's build default)
   * ``--primary`` for primary graph mode (recommended for most workloads)
   * ``--disk-swap-dir DIR`` to enable on-disk spill buffers
   * ``--anno-source`` (``header`` or ``filename``; default ``filename``)
   * ``--anno-type FMT`` to choose / add output annotation formats
   * ``--with-counts`` or ``--with-coords`` for count- / coordinate-aware
     annotation (mutually exclusive)
   * ``--graph EXISTING.dbg`` to reuse an already-built graph and run
     only the annotation + transform stages

   An example invocation:

   .. code-block:: bash

     metagraph-workflows build samples.txt -o /tmp/mygraph \
         -k 31 --primary \
         -p 34 --mem-gb 70 --disk-swap-dir /scratch/swap

3. Query the completed index with :ref:`metagraph query <query_index>`, or start
   ``metagraph server_query`` and connect through the :ref:`Python client <api>`.

The `end-to-end notebook <https://github.com/ratschlab/metagraph/blob/master/metagraph/workflows/notebooks/workflow_end_to_end_example.ipynb>`_
walks through indexing and API queries.

Count-aware annotations
^^^^^^^^^^^^^^^^^^^^^^^

The workflow supports these count-aware annotation formats:

* ``int_brwt``
* ``row_diff_int_brwt``
* ``row_diff_int_disk``

To enable counts explicitly, pass ``--with-counts``. If no annotation format is specified,
the default switches from ``relax.row_diff_brwt`` to ``row_diff_int_brwt``::

    metagraph-workflows build transcript_paths.txt -o [OUTPUT_DIR] \
        -k 31 --primary --with-counts --count-width 12

You can also select a count-aware format directly via ``--anno-type``; this
automatically enables count-aware mode::

    metagraph-workflows build transcript_paths.txt -o [OUTPUT_DIR] \
        -k 31 --primary --anno-type row_diff_int_brwt --count-width 12

Use ``--count-width`` to control the stored numeric range for counts
(valid range: ``2..32``, default: ``8``).

When reusing an output directory, the workflow keeps count and non-count intermediates in
separate mode-specific directories to avoid stale artifact reuse.

Coordinate-aware annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The workflow supports these coordinate-aware annotation formats:

* ``brwt_coord``
* ``row_diff_coord``
* ``row_diff_brwt_coord``
* ``row_diff_disk_coord``

To enable coordinates explicitly, pass ``--with-coords``. If no annotation format is specified,
the default switches from ``relax.row_diff_brwt`` to ``row_diff_brwt_coord``::

    metagraph-workflows build transcript_paths.txt -o [OUTPUT_DIR] \
        -k 31 --with-coords

You can also select a coordinate-aware format directly via ``--anno-type``; this
automatically enables coordinate-aware mode::

    metagraph-workflows build transcript_paths.txt -o [OUTPUT_DIR] \
        -k 31 --anno-type row_diff_brwt_coord

Coordinates are typically indexed for reference sequences, where preserving the original sequence context is important.
For this use case, primary graph mode is usually not recommended.

When coordinate-aware annotation is built with ``--anno-source filename``
(the default; one column per input file), the workflow additionally builds
a ``CoordToHeader`` sidecar at ``<graph>.seqs``. This is a second
``metagraph annotate --index-header-coords`` pass over the input fastas;
the pass uses the column order recovered from the final annotation via
``metagraph stats --print-col-names`` so that file labels and coord
offsets stay aligned even after BRWT clustering reorders columns. The
sidecar lets ``metagraph query --query-mode coords`` and ``metagraph
align`` report hits as ``<header>/<N>:<positions>`` instead of file-based
coordinates. In ``--anno-source header`` mode each sequence already has
its own column, so no sidecar is needed and the rule is skipped.

The loader always looks at ``<graph>.seqs`` (it strips the full
annotation extension and appends ``.seqs``), so a single sidecar is
emitted regardless of how many coord formats are requested. When
multiple are requested, the sidecar is derived from the first one in
the list -- pick that format when querying with coords.

Count-aware and coordinate-aware modes are mutually exclusive in this workflow.

Row-diff transform outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^

The annotate step writes ``columns.<mode>/<basename>.column.annodbg`` for each input sequence
file. Row-diff stages 0–2 then write under ``rd_cols.<mode>/`` (for example ``rd_cols.binary/``,
``rd_cols.coords/``, or ``rd_cols.counts.w8/`` depending on configuration):

* **Binary annotation mode** (no ``--with-counts`` / ``--with-coords``): stage 2 emits
  ``<basename>.row_diff.annodbg`` per column—the ``RowDiffColumnAnnotator`` format.
* **Count or coordinate mode**: stage 2 emits ``<basename>.column.annodbg`` and, when applicable,
  ``<basename>.column.annodbg.counts`` and/or ``<basename>.column.annodbg.coords``.

See ``metagraph-workflows build -h`` for more details.


Workflow management
~~~~~~~~~~~~~~~~~~~

The following snakemake options are exposed in the ``build`` subcommand

* ``--dryrun``: see what workflow steps would be done
* ``--force`` (corresponds to ``--forceall`` in snakemake): force run all steps

When building a new graph without ``-k``, a dry run calls the selected
MetaGraph executable's ``build --help`` to read its default k-mer length.
Pass ``-k`` explicitly to preview the workflow without that executable.
Reusing a graph with ``--graph`` does not need this lookup. A dry run does
not execute the workflow steps or validate the input alphabet.


Directly invoking Snakemake workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``metagraph-workflows`` command is only a wrapper around a snakemake workflow. You can also
directly invoke the snakemake workflow (assuming you checked out the `metagraph git repository <https://github.com/ratschlab/metagraph>`_):

.. code-block:: bash

    cat > build.yaml <<'YAML'
    k: 5
    seqs_file_list_path: transcript_paths.txt
    output_directory: /tmp/mygraph
    annotation_labels_source: header
    YAML
    snakemake --snakefile metagraph/workflows/snakemake/Snakefile \
        --configfile metagraph/workflows/snakemake/default.yml build.yaml --cores 2

The direct Snakemake invocation requires an explicit ``k`` in its configuration.
Use ``metagraph-workflows build`` to inherit the selected executable's default.
