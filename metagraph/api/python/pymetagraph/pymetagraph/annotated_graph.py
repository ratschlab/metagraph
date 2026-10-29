import math
from pathlib import Path
from typing import Iterable, List, Tuple, Union

from . import _pymetagraph_core

Hits = List[Tuple[str, int]]


class AnnotatedGraph:
    """An annotated de Bruijn graph loaded into the current process.

    :param graph: Path to the graph (``.dbg``).
    :param annotation: Path to the annotation (e.g. ``.column.annodbg``).
    :param threads: Number of threads used by :meth:`query`.
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
              num_top_labels: int = 1,
              discovery_fraction: float = 0.0,
              presence_fraction: float = 0.0) -> Union[Hits, List[Hits]]:
        """Top labels of each sequence with their numbers of matching k-mers.

        Same semantics as ``metagraph query``: a label is reported if at least
        ``discovery_fraction`` of the sequence's k-mers carry it, and only if at
        least ``presence_fraction`` of its k-mers are in the graph at all.

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
