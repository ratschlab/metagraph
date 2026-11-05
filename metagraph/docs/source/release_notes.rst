Release notes
=============

Unreleased
----------

``metagraph build`` now chooses its default k-mer length at compile time from
the alphabet's packed 64-bit encoding, rounding down to an odd length for
nucleotide alphabets:

* DNA/RNA: 31
* DNA5: 21
* Case-sensitive DNA: 15
* Protein: 12

To reproduce graphs built with the previous default, pass ``-k 3``. Explicit
``-k`` values continue to override the default. The workflow CLI's separate
default remains 31.
