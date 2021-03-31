===================
metagraph_workflows
===================

..
   
   .. image:: https://img.shields.io/pypi/v/metagraph_workflows.svg
           :target: https://pypi.python.org/pypi/metagraph_workflows

   .. image:: https://img.shields.io/travis/ratschlab/metagraph_workflows.svg
           :target: https://travis-ci.org/ratschlab/metagraph_workflows

   .. image:: https://readthedocs.org/projects/metagraph-workflows/badge/?version=latest
           :target: https://metagraph-workflows.readthedocs.io/en/latest/?badge=latest
           :alt: Documentation Status




Metagraph workflows


* Free software: MIT license
..
   * Documentation: https://metagraph-workflows.readthedocs.io.


Workflows for Creating Graphs and Annotations
---------------------------------------------

Since the creation of graph and indices comprises several steps, this package provides
some support to simplify these tasks - in particular for standard cases.

Given some raw sequence data and a few options like the kmer size (`k`) graphs and annotations
are automatically built:

.. code-block:: bash

    metagraph-utils build -k 5 transcript_paths.txt /tmp/mygraph


If you prefer invoking the workflow from within a python script, the following is equivalent:

.. code-block:: python

    from metagraph.cli import workflows
    workflows.run_build_workflow('transcript_paths.txt', '/tmp/mygraph', k=5)


The execution of the workflow is managed by https://snakemake.readthedocs.io/

Usage Example
~~~~~~~~~~~~~

Typically, the following steps would be performed:

0. installation: conda environemt
1. preparation files
2. running workflow
   - primary, non primary
   - TODO: other recommendations
3. do queries
    - command line
    - server and api



Workflow Management
~~~~~~~~~~~~~~~~~~~


The following options are exposed
 * `--force` (corresponds to `--forceall` in snakemake)
 * `--dryrun`
 * `--verbose`





Directly Invoking Snakemake Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above commands are only wrapper around a snakemake workflow. You can also
directly invoke the snakemake workflow:

.. code-block:: bash

    cd workflows
    snakemake --forceall --configfile default.yml \
        --config k=5 input_files_list_path='transcript_paths.txt' output_directory=/tmp/mygraph \
        annotation_labels_source=sequence_headers --cores 2




