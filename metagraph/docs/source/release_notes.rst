Release notes
=============

Unreleased
----------

Case-sensitive DNA string reverse complements now switch case, matching the
packed alphabet's strand convention: ``ACG`` becomes ``cgt``. This fixes
inconsistent graph construction, annotation, queries and sequence extraction.
Rebuild case-sensitive canonical and primary graphs and their annotations.
Also rebuild basic graphs created with ``--fwd-and-reverse`` and regenerate
outputs that used the previous string reverse-complement operation.
Ordinary DNA, DNA5 and protein builds keep their previous string behavior.

SSHash protein graphs with ``k=12`` now use 64-bit k-mer keys instead of 128-bit
keys. Rebuild existing SSHash protein graphs with ``k=12`` and their dependent
annotations before using this version. Other k-mer lengths retain their existing
storage width.

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
