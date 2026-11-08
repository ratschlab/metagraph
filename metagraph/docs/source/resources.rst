.. _resources:

Resources
=========

.. _download_indexes:

Preconstructed Indexes
----------------------

Access MetaGraph indexes via AWS S3 for local analysis using the MetaGraph command line tool.
All indexes are available for download from `AWS Open Data <https://registry.opendata.aws/metagraph/>`_. Install the AWS CLI following the `installation guide <https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html>`_ (supports Windows, macOS, Linux).

Example commands:

.. code-block:: bash

    # List available objects in a bucket
    aws s3 --no-sign-request ls s3://metagraph/refseq/

    # Replace FILE.dbg with an object name returned by the listing
    aws s3 --no-sign-request cp s3://metagraph/refseq/FILE.dbg .

    # Download an entire index directory (check its size before starting)
    aws s3 --no-sign-request sync s3://metagraph/refseq/ ./local-refseq/

The ``--no-sign-request`` flag indicates public access without AWS credentials.

Once downloaded, the indexes can be queried locally with the command line interface and Python API (see :ref:`query_index`).

For more information about the indexes and deployment manuals, see https://github.com/ratschlab/metagraph-open-data.

.. _metagraph_online:

MetaGraph Online
----------------

`MetaGraph Online <https://metagraph.ethz.ch/search>`_ provides a web interface for searching DNA, RNA,
and protein sequences across various public archives indexed with annotated de Bruijn graphs.

**Web Interface Features:**
    - Exact k-mer search for DNA, RNA, and protein sequences
    - Standard jobs: up to 10 sequences and 50,000 characters total
    - Large jobs: up to 10 MB total, with no sequence-count limit
    - Access to various indexes (for current coverage see `Databases <https://metagraph.ethz.ch/indexes>`_)
        - Microbe, Fungi, Plants, Metazoa (including Human/Mouse)
        - Reference/assembled sets (RefSeq, UHGG, Tara Oceans)
        - Protein databases (UniParc)
    - Export results as CSV or JSON
    - Shareable results: 48-hour expiration by default; locking guarantees retention for at least 6 months

**API Access:**
The production API endpoint is available at: ``https://metagraph.ethz.ch:8081/search``

For detailed API documentation and examples, see the `MetaGraph Online Help <https://metagraph.ethz.ch/help#api-cli>`_ page.

.. note:: Python API described in the :ref:`API documentation <api>` is used in the internal implementation of the service
    and can also be used to query indexes hosted locally. To query indexes hosted publicly, refer to
    `MetaGraph Online Help <https://metagraph.ethz.ch/help#api-cli>`_.

**Performance Tips:**
    - Use exact matching for high identity matches (fast)
    - For alignment of divergent or noisy sequences, download an index and use ``metagraph align`` locally
    - For larger jobs, use the API (rate limits apply for queued jobs) or download :ref:`download_indexes` and query locally
    - Split large inputs to avoid timeouts

Contact & Support
-----------------

For issues or feature requests:
    * Email: metagraph@inf.ethz.ch
    * GitHub: `ratschlab/metagraph <https://github.com/ratschlab/metagraph>`_
    * Web feedback: `Contact form <https://metagraph.ethz.ch/feedback>`_

