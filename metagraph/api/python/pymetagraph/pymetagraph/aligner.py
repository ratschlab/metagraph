from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from readfish._config import Barcode, Region
from readfish.plugins.abc import AlignerABC
from readfish.plugins.utils import Result

from . import _pymetagraph_core


class Aligner(AlignerABC):
    """Readfish Aligner plugin backed by a metagraph graph + annotation index.

    Two unrelated query methods are available, picked via ``method``:

    - ``"query"`` (default): label-presence hits via
      ``AnnotatedDBG::get_top_labels``. ``Alignment.r_st``/``r_en``/``strand``
      are placeholder values (``0``, ``len(seq)``, ``+1``) since only the
      matched label is meaningful for accept/reject decisions in this mode.
    - ``"align"``: real seed-and-extend alignment via metagraph's
      alignment_redone module, the same code path backing the ``metagraph
      align`` CLI command and the server's ``/align`` endpoint (see the old
      client's ``GraphClient.align``). ``r_st``/``r_en``/``strand`` reflect
      the actual matched region and orientation of the read.

    :param debug_log: Filename (or ``"stdout"``/``"stderr"``) for debug logging.
    :keyword input: Path to the metagraph graph file (``.dbg``).
    :keyword annotator: Path to the metagraph annotation file. Required for
        ``method="query"``; optional for ``method="align"`` (omitting it runs
        plain graph alignment with no label lookups).
    :keyword threads: Number of threads to use for querying (default: 1).
    :keyword method: Either ``"query"`` (default) or ``"align"``.
    :keyword num_top_labels: Max number of labels to consider per query.
        Only used by ``method="query"``.
    :keyword discovery_fraction: Min fraction of k-mers required to be present
        in a label for it to count as a hit (default: 0.7). Only used by
        ``method="query"``.
    :keyword presence_fraction: Min fraction of k-mers required to be present
        in the graph at all before querying labels (default: 0.0). Only used
        by ``method="query"``.
    :keyword seed_length: Minimum seed length. Only used by ``method="align"``.
    :keyword max_alternative_alignments: Number of alternative paths to
        consider per seed (default: 1). Only used by ``method="align"``.
    :keyword max_num_nodes_per_seq_char: Max nodes to consider per sequence
        character during extension. Only used by ``method="align"``.
    :keyword min_exact_match: Min fraction of nucleotides covered by seeds
        required to align (default: 0.7). Only used by ``method="align"``.
    :keyword connect_anchors: If ``True`` (default), traverse the graph to
        align query regions falling between anchors. Only used by
        ``method="align"``.
    :keyword extend_chains: If ``True`` (default), perform ends-free
        extension from the first/last anchors in a chain. Only used by
        ``method="align"``.
    """

    def __init__(self, debug_log: Optional[str] = None, **kwargs):
        if debug_log:
            if debug_log == "stdout":
                self.logfile = sys.stdout
            elif debug_log == "stderr":
                self.logfile = sys.stderr
            else:
                self.logfile = open(debug_log, "w")
        else:
            self.logfile = None

        if "input" not in kwargs:
            raise AttributeError('Required argument "input" not found.')

        self.kwargs = kwargs
        self.validate()
        self.index = _pymetagraph_core.Index(**kwargs)

    def validate(self) -> None:
        index_path = Path(self.kwargs["input"])
        if not index_path.is_file():
            raise FileNotFoundError(f"{index_path} does not exist")

        # "query" (the default) needs an annotation to have any labels to
        # report; "align" can run against the bare graph (see _pymetagraph_core.Index).
        method = self.kwargs.get("method", "query")
        annotator = self.kwargs.get("annotator")
        if not annotator:
            if method != "align":
                raise AttributeError('Required argument "annotator" not found.')
            return
        if not Path(annotator).is_file():
            raise FileNotFoundError(f"{annotator} does not exist")

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: List[Region], barcodes: Dict[str, Barcode]) -> str:
        return (
            f"Using the pymetagraph plugin. Graph: {self.kwargs['input']}, "
            f"annotator: {self.kwargs.get('annotator', 'none')}."
        )

    def map_reads(self, basecall_results: Iterable[Result]) -> Iterable[Result]:
        skipped = []
        kept = []
        for result in basecall_results:
            if result.seq:
                kept.append(result)
            else:
                skipped.append(result)

        if kept:
            alignments = self.index.query_batch([r.seq for r in kept])
            for result, alignment in zip(kept, alignments):
                result.alignment_data = [] if alignment.ctg == "*" else [alignment]
                yield result

        for result in skipped:
            result.alignment_data = []
            yield result

    def disconnect(self) -> None:
        self.index = None
        if self.logfile and self.logfile not in (sys.stdout, sys.stderr):
            self.logfile.close()
