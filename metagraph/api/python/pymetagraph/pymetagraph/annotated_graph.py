import math
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union

from . import _pymetagraph_core

Hits = List[Tuple[str, int]]
GraphAlignment = _pymetagraph_core.GraphAlignment
LabelMapping = _pymetagraph_core.LabelMapping


class AnnotatedGraph:
    """An annotated de Bruijn graph loaded into the current process.

    :param graph: Path to the graph (``.dbg``).
    :param annotation: Path to the annotation (e.g. ``.column.annodbg``).
    :param threads: Number of threads used by :meth:`query`, :meth:`align` and :meth:`map`.
    """

    def __init__(self, graph: Union[str, Path], annotation: Union[str, Path], threads: int = 1):
        # metagraph's loaders exit() the process on missing files, so check here first
        for path in (graph, annotation):
            if not Path(path).is_file():
                raise FileNotFoundError(f"{path} does not exist")
        if threads < 1:
            raise ValueError("threads must be at least 1")
        self._graph = _pymetagraph_core.AnnotatedGraph(str(graph), str(annotation), threads)

    @property
    def k(self) -> int:
        """k-mer length of the graph."""
        return self._graph.k

    @property
    def labels(self) -> List[str]:
        """All labels in the annotation."""
        return self._graph.labels

    def query(self, seqs: Union[str, Iterable[str]],
              num_top_labels: Optional[int] = None,
              discovery_fraction: float = 0.7,
              presence_fraction: float = 0.0) -> Union[Hits, List[Hits]]:
        """Top labels of each sequence with their numbers of matching k-mers.

        Same semantics as ``metagraph query``: a label is reported if at least
        ``discovery_fraction`` of the sequence's k-mers carry it, and only if at
        least ``presence_fraction`` of its k-mers are in the graph at all.
        The defaults are those of ``metagraph query``: all labels
        (``num_top_labels=None``), ``discovery_fraction=0.7`` and
        ``presence_fraction=0.0``.

        :param seqs: A sequence, or an iterable of sequences.
        :return: For each sequence, a list of ``(label, kmer_count)`` sorted by
            decreasing count, empty if nothing passes the thresholds or the
            sequence is shorter than ``k``. A single sequence gives a single list.
        """
        for fraction in (discovery_fraction, presence_fraction):
            if not (math.isfinite(fraction) and 0 <= fraction <= 1):
                raise ValueError("discovery_fraction and presence_fraction must be in [0, 1]")
        if isinstance(seqs, str):
            return self._graph.query([seqs], num_top_labels, discovery_fraction, presence_fraction)[0]
        return self._graph.query(list(seqs), num_top_labels, discovery_fraction, presence_fraction)

    def align(self, seqs: Union[str, Iterable[str]],
              labeled: bool = True,
              min_exact_match: Optional[float] = None,
              max_alternative_alignments: Optional[int] = None,
              min_seed_length: Optional[int] = None,
              max_seed_length: Optional[int] = None,
              max_nodes_per_seq_char: Optional[float] = None
              ) -> Union[List[GraphAlignment], List[List[GraphAlignment]]]:
        """Alignments of each sequence to the graph, as in ``metagraph align``.

        With ``labeled=True`` the alignments follow the annotation and report
        its labels (``metagraph align -a``). With ``labeled=False`` they are
        unlabeled, as from the server's ``/align`` endpoint. Arguments left at
        ``None`` keep the defaults of ``metagraph align``: ``min_exact_match=0.7``
        (``--align-min-exact-match``), one alignment per sequence
        (``--align-alternative-alignments``), seeds of at least ``min(k, 19)``
        characters (``--align-min-seed-length``) and
        ``max_nodes_per_seq_char=5`` (``--align-max-nodes-per-seq-char``).

        :param seqs: A sequence, or an iterable of sequences.
        :return: For each sequence, a list of :class:`GraphAlignment` (labels,
            score, max_score, cigar, orientation, query_begin, query_end), best
            first, empty if the sequence does not align. Query coordinates refer
            to the reverse complement if ``orientation`` is True. A single
            sequence gives a single list.
        """
        args = (labeled, min_exact_match, max_alternative_alignments,
                min_seed_length, max_seed_length, max_nodes_per_seq_char)
        if isinstance(seqs, str):
            return self._graph.align([seqs], *args)[0]
        return self._graph.align(list(seqs), *args)

    def map(self, seqs: Union[str, Iterable[str]],
            min_exact_match: Optional[float] = None,
            min_seed_length: Optional[int] = None,
            max_seed_length: Optional[int] = None
            ) -> Union[Optional[LabelMapping], List[Optional[LabelMapping]]]:
        """Labels of each sequence from its seeds alone, without alignment.

        Same as ``metagraph align -a --align-map-only``: seeds are found and
        filtered by label as for labeled :meth:`align`, but not extended, which
        is much faster. Unset arguments keep the defaults of ``metagraph align``.

        :param seqs: A sequence, or an iterable of sequences.
        :return: For each sequence, a :class:`LabelMapping` (the labels whose
            seeds cover the most query characters on the better strand,
            num_matches, orientation, and the query_begin/query_end they span on
            the forward query), or None. A single sequence gives a single result.
        """
        args = (min_exact_match, min_seed_length, max_seed_length)
        if isinstance(seqs, str):
            return self._graph.map([seqs], *args)[0]
        return self._graph.map(list(seqs), *args)
