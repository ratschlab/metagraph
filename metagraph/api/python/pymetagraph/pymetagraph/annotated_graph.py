import json
from pathlib import Path
from typing import Iterable, List, Tuple, Union, Optional

from . import _pymetagraph_core

Hits = List[Tuple[str, int]]


class AnnotatedGraph:
    """An annotated de Bruijn graph loaded into the current process.

    :param graph: Path to the graph (``.dbg``).
    :param annotation: Optional annotation path; required by :meth:`query`.
    :param threads: Number of workers used by :meth:`query` and :meth:`align`.
    """

    def __init__(self, graph: Union[str, Path], annotation: Optional[Union[str, Path]] = None, threads: int = 1):
        # metagraph's loaders exit() the process on missing files, so check here first
        for path in (graph,) if annotation is None else (graph, annotation):
            if not Path(path).is_file():
                raise FileNotFoundError(f"{path} does not exist")
        if threads < 1:
            raise ValueError("threads must be at least 1")
        self._graph = _pymetagraph_core.AnnotatedGraph(str(graph), str(annotation) if annotation is not None else "", threads)

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
              discovery_fraction: float = 0.7,
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
        if isinstance(seqs, str):
            return self._graph.query([seqs], num_top_labels, discovery_fraction, presence_fraction)[0]
        return self._graph.query(list(seqs), num_top_labels, discovery_fraction, presence_fraction)

    def align(self, seqs: Union[str, Iterable[str]], *, seed_length: int = 0,
              max_alternative_alignments: int = 1, min_exact_match: float = 0.7,
              xdrop: int = 27,
              rel_score_cutoff: float = 0.95,
              only_forwards: bool = False, connect_anchors: bool = True,
              extend_chains: bool = True, max_num_seeds_per_locus: int = 1000,
              max_num_nodes_per_seq_char: float = 5.0, max_ram_per_alignment: float = 200.0,
              chain_alignments: bool = False, two_sided: bool = False):
        """Align sequences, optionally constrained by annotation labels.

        Return one list of result dictionaries per query, or just that list
        for a string input. Results use the CLI JSON fields, plus ``strand``
        and zero-based, half-open ``query_start``/``query_end`` positions in
        the original input. Label ``nt_coords`` are reference coordinates in
        the CLI's one-based inclusive format; they are not query positions.
        ``seed_length=0`` selects min(19, k). No match gives an empty list.
        Chaining and two-sided extension are opt-in. ``only_forwards`` selects
        query orientation independently of these search options.
        """
        single = isinstance(seqs, str)
        batches = self._graph.align(
            [seqs] if single else list(seqs), seed_length, max_alternative_alignments,
            min_exact_match, xdrop, rel_score_cutoff, only_forwards, connect_anchors, extend_chains,
            max_num_seeds_per_locus, max_num_nodes_per_seq_char, max_ram_per_alignment,
            chain_alignments, two_sided)
        results = [[json.loads(result) for result in batch] for batch in batches]
        return results[0] if single else results
