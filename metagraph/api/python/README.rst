====================
MetaGraph HTTP client
====================

Query a running MetaGraph index from Python. The client supports exact k-mer search
and sequence alignment and returns results as pandas data frames.

Install
-------

From this directory in a MetaGraph checkout, install into your active Python environment:

.. code-block:: bash

    python -m pip install .

Usage
-----

Connect to the hosted MetaSUB index:

.. code-block:: python

    from metagraph.client import GraphClient

    metasub = GraphClient('dnaloc.ethz.ch', 80, api_path='/api/metasub19')
    results = metasub.search('CAAATTCTTGTAAGTATTAAACATTGTATATGTATTTTGAA')

For a local index, first start ``metagraph server_query`` with your graph and annotation,
then connect with ``GraphClient('localhost', 5555)``. See the
`API manual <https://metagraph.ethz.ch/static/docs/api.html>`_ for server setup,
query options, and parallel searches, or explore the `example notebooks <./notebooks>`_.

The separate `MetaGraph Online API <https://metagraph.ethz.ch/help#api-cli>`_ provides
asynchronous searches across the public index collection.
