"""metagraph.traverse: retrieve a traversal once, process it locally.

The MetaGraph query server returns a traversal as a compact, lossless text per seed --
the graphlet, MGT v1 (docs/DESIGN-traverse-graphlet.md §2) -- embedded in a small JSON
summary (`strategy.output.detail: "graphlet"`). This package parses it and answers the
follow-up questions an agent asks about one retrieval locally: walks, claims, label
routes, support profiles, comparisons, FASTA/GFA, continuations and the next request.
Nothing here reads the graph; deepening is a new backend traversal that the library
prepares (next_request) and the client runs.

Stdlib only. pandas is imported lazily by frames.frames() alone, so
`import metagraph.traverse` works without it.

    from metagraph.traverse import TraverseClient, Graphlet
    resp = TraverseClient('localhost', 5555).traverse([{'sequence': seq}], strategy)
    g = resp.graphlets[0]
    for w in g.walks('right', top=5):
        print(w.path_id, w.length_bp, [l.name for l in w.labels_full])
"""

from ._codec import CodecError, GraphletFormatError, UNLIMITED
from .model import (
    AmbiguousLabel, Arm, BadSelector, Change, Claim, Comparison, Continuation, Graphlet,
    IncompatibleContinuations, IncompleteRecording, Label, LabelWalk, MissingEnvelope,
    NextRequest, Segment, Run, SupportRun, UnknownLabel, Walk,
)
from .parser import FORMAT_VERSION, dump, from_response, is_canonical, load, parse, save
from .ops import GraphletView
from .client import (AttemptAtBound, ServerInitializing, TraverseClient, TraverseError,
                     TraverseResponse)
from .store import Entry, GraphletStore, StoreLimitExceeded, UnknownHandle

__all__ = [
    'FORMAT_VERSION', 'parse', 'dump', 'is_canonical', 'load', 'save', 'from_response',
    'Graphlet', 'GraphletView', 'Label', 'Arm', 'Segment', 'Run', 'Walk', 'Claim',
    'LabelWalk', 'SupportRun', 'Change', 'Continuation', 'NextRequest', 'Comparison',
    'GraphletFormatError', 'CodecError', 'MissingEnvelope', 'AmbiguousLabel',
    'UnknownLabel', 'BadSelector', 'IncompleteRecording', 'IncompatibleContinuations',
    'UNLIMITED',
    'TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
    'AttemptAtBound', 'GraphletStore', 'Entry', 'UnknownHandle', 'StoreLimitExceeded',
]
