.. _resources:

Resources	
=========

Preconstructed indexes
----------------------

The large preconstructed indexes are available from https://drive.google.com/drive/folders/1sOrgix1FWnH7SJpmY0P6gwl4wT5AEQ1B?usp=sharing.
(For more stable transmission, consider using `rclone <https://rclone.org/>`_, which reliably downloads files up to and larger 1 TB, from the command line.)

MetaGraph Online
----------------

`MetaGraph Online <https://metagraph.ethz.ch/search>`_ currently hosts various indexes constructed
from SRA runs and supports search and align queries.

To query these indexes using Python API, connect to their respective `endpoints <https://metagraph.ethz.ch/graphs>`_.
For example, write ``GraphClient('dnaloc.ethz.ch', 80, api_path='/api/metasub19')`` to connect to the MetaSUB index.

