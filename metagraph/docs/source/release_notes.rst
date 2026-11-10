Release notes
=============

Unreleased
----------

Annotations now load correctly with any build. Earlier versions saved the hash
table of label names as the writing build made it, and that hash differs between
C++ standard libraries (for example, libstdc++ on Linux and libc++ on macOS).
Label lookups then failed when a build with another standard library loaded the
annotation. This version rebuilds the table on load, which takes about 0.1
seconds per million labels. The file format is unchanged.

The HTTP ``/align`` endpoint now rejects invalid ``max_alternative_alignments``
values with status 400. Use a positive integer. Zero previously selected one
alignment; negative values could become very large unsigned counts.

A ``.seqs`` coordinate-to-header mapping is now checked when its annotation
loads. A mapping with a different number of columns is refused at load, and
``server_query`` exits before it serves requests. Previously the first query
that used the mapping stopped the process. ``annotate --index-header-coords``
now stores the label of each column in the ``.seqs`` file, and a load refuses a
mapping whose labels differ from the annotation's labels. With ``-a``, the
command spot-checks each input file against its annotation column first: the
column must hold the file's first and last k-mers in the graph, and a file with
no k-mer of the graph is refused, an empty file included. The check cannot tell
a file of the same sequences under other headers, or a copy cut short after its
last probed k-mer, so the records still rest on passing the files the annotation
was built from. Earlier versions ignore the stored labels.

Case-sensitive DNA packed reverse complements now preserve case, matching
string reverse complements: ``ACG`` becomes ``CGT``, and ``acg`` becomes ``cgt``.
Case can retain soft-masking information; uppercase and lowercase bases remain
distinct symbols during indexing and matching.
Rebuild existing case-sensitive canonical and primary graphs and their dependent
annotations. Basic graphs retain their previous behavior.

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
