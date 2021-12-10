=========
Workflows
=========

Since each indexing workflow in `MetaGraph <https://metagraph.ethz.ch>`_ comprises
several steps, here we provide pipelines to make the process easier and more straightforward
for the most common scenarios.


Installation
------------


Set up a conda environment and install the necessary packages using:

.. code-block:: bash

   conda create -n metagraph-workflows python=3.8
   conda activate metagraph-workflows
   conda install -c bioconda -c conda-forge metagraph
   pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/workflows"


Creating graphs and annotations
-------------------------------

Given some raw sequencing data and a few options like the k-mer length, graphs and annotations
are automatically built::

    metagraph-workflows build -k 5 transcript_paths.txt /tmp/mygraph


The same pipeline can be invoked from within a python script:

.. code-block:: python

    from metagraph_workflows import workflows

    workflows.run_build_workflow('/tmp/mygraph', seqs_file_list_path='transcript_paths.txt', k=5)



The pipelines are written in the `Snakemake <https://snakemake.readthedocs.io/>`__ workflow management system and can also be directly invoked using the ``snakemake`` command line tool (see below).


Usage
~~~~~

Typically, the following steps would be performed:

1. Prepare a list of files for indexing.
2. Construct a MetaGraph index: invoke a workflow using ``metagraph-workflows build``. Important parameters you may consider tuning are:

   * k-mer length
   * basic vs. primary graph mode
   * source of annotation labels: ``sequence_headers`` or ``sequence_file_names``

   An example invocation:

   .. code-block:: bash

     metagraph-workflows build -k 31 \
                               --seqs-dir-path [PATH_TO_FILES] \
                               --annotation-labels-source sequence_headers \
                               --build-primary-graph \
                               [OUTPUT_DIR]

   See ``metagraph-workflows build -h`` for more help.
3. Once a MetaGraph index has been created, it can be queried either by using the command line
   ``metagraph`` tool or by starting the metagraph server directly on a laptop or on another suitable
   machine and querying it using the python :ref:`API` client.


There is also a `jupyter notebook <https://github.com/ratschlab/metagraph/blob/master/metagraph/workflows/notebooks/workflow_end_to_end_example.ipynb>`_ showing the whole process: from indexing to api querying  on a simple example.



Workflow management
~~~~~~~~~~~~~~~~~~~

The following snakemake options are exposed in the ``build`` subcommand

* ``--dryrun``: see what workflow steps would be done
* ``--force`` (corresponds to ``--forceall`` in snakemake): force run all steps


Directly invoking Snakemake workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above command is only a wrapper around a snakemake workflow. You can also
directly invoke the snakemake workflow (assuming you checked out the `metagraph git repository <https://github.com/ratschlab/metagraph>`_):

.. code-block:: bash

    cd metagraph/workflows
    snakemake --forceall --configfile default.yml \
        --config k=5 seqs_file_list_path='transcript_paths.txt' output_directory=/tmp/mygraph \
        annotation_labels_source=sequence_headers --cores 2
