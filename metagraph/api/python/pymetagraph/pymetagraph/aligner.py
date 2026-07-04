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

    This reports label-presence hits (via ``AnnotatedDBG::get_top_labels``)
    rather than real base-level alignments: ``Alignment.r_st``/``r_en``/``strand``
    are placeholder values (``0``, ``len(seq)``, ``+1``) since only the matched
    label is meaningful for accept/reject decisions here.

    :param debug_log: Filename (or ``"stdout"``/``"stderr"``) for debug logging.
    :keyword input: Path to the metagraph graph file (``.dbg``).
    :keyword annotator: Path to the metagraph annotation file.
    :keyword threads: Number of threads to use for querying (default: 1).
    :keyword num_top_labels: Max number of labels to consider per query.
    :keyword discovery_fraction: Min fraction of k-mers required to be present
        in a label for it to count as a hit (default: 0.7).
    :keyword presence_fraction: Min fraction of k-mers required to be present
        in the graph at all before querying labels (default: 0.0).
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

        annotator = self.kwargs.get("annotator")
        if not annotator:
            raise AttributeError('Required argument "annotator" not found.')
        if not Path(annotator).is_file():
            raise FileNotFoundError(f"{annotator} does not exist")

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: List[Region], barcodes: Dict[str, Barcode]) -> str:
        return (
            f"Using the pymetagraph plugin. Graph: {self.kwargs['input']}, "
            f"annotator: {self.kwargs['annotator']}."
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
