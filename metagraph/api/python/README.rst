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

    from metagraph.client import Client

    metagraph = Client()

    metagraph.connect('localhost', 8555)

    metagraph.search('CAAATTCTTGTAAGTATTAAACATTGTATATGTATTTTGAA')
    metagraph.search('GAATGAAAGATATGTGTTTTTCA')
    metagraph.search('CATCTACTTGACTGGATTAAGAGACACACA')


For more examples, see `notebooks
<./notebooks>`_.

