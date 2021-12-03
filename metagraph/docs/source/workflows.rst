=========
Workflows
=========

This package provides workflows for the `metagraph framework
<https://metagraph.ethz.ch>`_


Workflows for Creating Graphs and Annotations
---------------------------------------------

Since the creation of graph and indices comprises several steps, this package provides
some support to simplify these tasks - in particular for standard cases.

Given some raw sequence data and a few options like the kmer size (`k`) graphs and annotations
are automatically built:

.. code-block:: bash

    metagraph-workflows build -k 5 transcript_paths.txt /tmp/mygraph


If you prefer invoking the workflow from within a python script, the following is equivalent:

.. code-block:: python

    from metagraph_workflows import workflows
    workflows.run_build_workflow('/tmp/mygraph', seqs_file_list_path='transcript_paths.txt', k=5)



The workflow logic itself is expressed as a `Snakemake workflow
<https://snakemake.readthedocs.io/>`_ . You can also directly invoke the workflows
using the `snakemake` command line tool (see below).


Installation and Set up
~~~~~~~~~~~~~~~~~~~~~~~


Set up a conda environment and install the necessary packages using:

.. code-block:: bash

   conda create -n metagraph-workflows python=3.8
   conda activate metagraph-workflows
   conda install -c bioconda -c conda-forge metagraph
   pip install -U "git+https://github.com/ratschlab/metagraph.git#subdirectory=metagraph/workflows"




Usage Example
~~~~~~~~~~~~~

Typically, the following steps would be performed:

1. sequence file preparation: add your sequence files of interest into a directory.
2. running workflow: you can invoke the workflow using ``metagraph-workflows build``. Important parameters you may consider tuning are:

   * k
   * primary vs non primary graph creation
   * annotation label source: ``sequence_headers`` or ``sequence_file_names``

   An example invocation:

   .. code-block:: bash

     metagraph-workflows build -k 31 \
                               --seqs-dir-path [PATH_TO_SEQUENCES] \
                               --annotation-labels-source sequence_headers \
                               --build-primary-graph
                               [OUTPUT_DIR]

   see ``metagraph-workflows build --help`` for more help
3. do queries: once you created the indices you can query either by using the command line
   query tool or starting the metagraph server on your laptop or another suitable machine and access
   do queries using e.g. the python :ref:`API` client.


There is also a `jupyter notebook <https://github.com/ratschlab/metagraph/blob/master/metagraph/workflows/notebooks/workflow_end_to_end_example.ipynb>`_ walking you through an example from indexing to api querying.



Workflow Management
~~~~~~~~~~~~~~~~~~~

The following snakemake options are exposed in the ``build`` subcommand

   * ``--dryrun``: see what workflow steps would be done
   * ``--force`` (corresponds to ``--forceall`` in snakemake): force run all steps


Directly Invoking Snakemake Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above command is only a wrapper around a snakemake workflow. You can also
directly invoke the snakemake workflow (assuming you checked out the `metagraph git repository <https://github.com/ratschlab/metagraph>`_):

.. code-block:: bash

    cd metagraph/workflows
    snakemake --forceall --configfile default.yml \
        --config k=5 seqs_file_list_path='transcript_paths.txt' output_directory=/tmp/mygraph \
        annotation_labels_source=sequence_headers --cores 2
