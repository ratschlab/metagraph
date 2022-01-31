=========
Metagraph
=========


Metagraph Python API


Install
--------
.. code-block:: bash

    make install


Usage
--------

.. code-block:: python

    from metagraph.client import GraphClient

    metagraph = GraphClient('dnaloc.ethz.ch', 80, api_path='/api/metasub19')

    metagraph.search('CAAATTCTTGTAAGTATTAAACATTGTATATGTATTTTGAA')
    metagraph.search('GAATGAAAGATATGTGTTTTTCA')
    metagraph.search('CATCTACTTGACTGGATTAAGAGACACACA')


For more examples, see `notebooks
<./notebooks>`_.

See more API docs in `<../../docs/source/api.rst>`_ on online: `<https://metagraph.ethz.ch/static/docs/api.html>`_.
