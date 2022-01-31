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
