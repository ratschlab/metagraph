===================
Snakemake workflows
===================

Since each indexing workflow in `MetaGraph <https://metagraph.ethz.ch>`_ comprises
several steps, we provide automated pipelines to make the process easier and more straightforward
for the most common scenarios.


Installation
------------

The ``metagraph-workflows`` CLI ships as a separate Python package. The
``metagraph`` conda recipe only installs the C++ binary, so the workflow
wrapper needs an extra ``pip install`` step alongside it:

.. code-block:: bash

   conda create -n metagraph-workflows python=3.8
   conda activate metagraph-workflows
   conda install -c bioconda -c conda-forge metagraph       # the metagraph binary
   pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/workflows"

After this, the ``metagraph`` binary and the ``metagraph-workflows``
command are both available on ``PATH``.


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

   Important parameters you may want to set:

   * ``--threads N`` and ``--mem-gb GB`` for the hardware budget
   * ``-k`` for k-mer length (default 31)
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
         --threads 34 --mem-gb 70 --disk-swap-dir /scratch/swap

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

3. Once a MetaGraph index has been created, it can be queried either by using the command line
   ``metagraph`` tool or by starting the MetaGraph server directly on a laptop or on another suitable
   machine and querying it using the python :ref:`API` client.


There is also a `jupyter notebook <https://github.com/ratschlab/metagraph/blob/master/metagraph/workflows/notebooks/workflow_end_to_end_example.ipynb>`_ showing the whole process: from indexing to api querying  on a simple example.



Workflow management
~~~~~~~~~~~~~~~~~~~

The following snakemake options are exposed in the ``build`` subcommand

* ``--dryrun``: see what workflow steps would be done
* ``--force`` (corresponds to ``--forceall`` in snakemake): force run all steps


Directly invoking Snakemake workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``metagraph-workflows`` command is only a wrapper around a snakemake workflow. You can also
directly invoke the snakemake workflow (assuming you checked out the `metagraph git repository <https://github.com/ratschlab/metagraph>`_):

.. code-block:: bash

    cd metagraph/workflows
    snakemake --forceall --configfile default.yml \
        --config k=5 seqs_file_list_path='transcript_paths.txt' output_directory=/tmp/mygraph \
        annotation_labels_source=header --cores 2
