.. _sequence_assembly:

Sequence assembly
=================

.. attention:: This page is in development

Basic assembly
--------------

See how to assemble simple contigs and unitigs in :ref:`to_sequences`.

Differential assembly
---------------------

MetaGraph allows for sequences to be assembled from a defined subgraph of an input annotated de Bruijn graph. In particular, a subset of labels may be defined as a *foreground* or **in-group**, while another may be defined as the *background* or **out-group**. From these, sequences can be assembled from a subgraph predominantly containing foreground k-mers and excluding background k-mers.

These subgraphs can be defined using a JSON configuration file. Examples are provided `here <https://github.com/ratschlab/metagraph/blob/master/metagraph/tests/data/example.diff.json>`__ and `here <https://github.com/ratschlab/metagraph/blob/master/metagraph/tests/data/example_simple.diff.json>`__.

Several differential assembly experiments may be defined for a single input annotated graph. These are organized into *groups* and *experiments*. A *group* is defined by a list of experiments and sets of foreground and background labels which are shared by these experiments. We recommend that labels with a substantially larger number of k-mers, such as reference genomes, are defined as shared. An *experiment* is defined by its own set of foreground and background labels, and by the assembly experiment parameters.

The assembly proceeds in four stages. In the first stage, the k-mers labeled with at least one of the foreground and background labels are enumerated to determine the initial subgraph. In the second stage, all labels associated with the subgraph unitigs are checked against the group's shared label set, updating the foreground and background label counts accordingly. In the third stage, the foreground and background label counts, and assembly parameters, are used to determine which k-mers are filtered out from the subgraph. In the fourth stage, sequences are assembled from the subgraph constructed in stage three.

Differential assembly parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the third stage of differential assembly, several parameters defined in the config JSON files are used to determine which k-mers from the initial subgraph are filtered out. A k-mer is considered to be an *in-k-mer* if it is annotated with ``in_min_fraction`` of the foreground labels (of the experiment plus the group's shared labels). Analogously, a k-mer is considered to be an *out-k-mer* if it is annotated with more than ``out_max_fraction`` of the background labels. Next, the numbers of in- and out-k-mers are considered in each unitig of the initial subgraph. An entire unitig is discarded from the subgraph if less than ``unitig_in_min_fraction`` of its constituent k-mers are in-k-mers, if more than ``unitig_out_max_fraction`` of these k-mers are out-k-mers, or if more than ``unitig_other_max_fraction`` of the constituent k-mers contain other (non-in- or non-out-) labels.

Each sequence which is extracted from an experiment is labeled with the ``name`` defined for that experiment.
