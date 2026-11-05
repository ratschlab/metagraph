Release notes
=============

Unreleased
----------

Case-sensitive DNA reverse complements now preserve case consistently during
construction and querying (``A`` complements ``T``, and ``a`` complements ``t``).
Rebuild existing case-sensitive canonical and primary graphs and their dependent
annotations. Basic graphs retain their previous behavior.

SSHash now supports DNA5 and case-sensitive DNA graphs, treating ``N`` as a
literal symbol. Rebuild SSHash graphs made with older DNA5 or case-sensitive
executables and their dependent annotations.

``metagraph build`` now chooses its default k-mer length at compile time from
the alphabet's packed 64-bit encoding, rounding down to an odd length for
nucleotide alphabets:

* DNA/RNA: 31
* DNA5: 21
* Case-sensitive DNA: 15
* Protein: 12

To reproduce graphs built with the previous default, pass ``-k 3``. Explicit
``-k`` values continue to override the default.

The workflow CLI now uses the selected MetaGraph executable's build default
instead of a separate default of 31. It records the resolved value in
``config.yaml`` and uses it throughout the build, including KMC stages. Pass
``-k 31`` to retain the previous workflow default. When reusing a graph with
``--graph``, its stored k-mer length is unchanged.
