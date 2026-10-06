"""metagraph.traverse: retrieve a traversal once, process it locally.

The MetaGraph query server returns a traversal as a compact, lossless text per seed --
the graphlet, MGT v1 (docs/DESIGN-traverse-graphlet.md §2) -- embedded in a small JSON
summary (`strategy.output.detail: "graphlet"`). This package parses it and answers the
follow-up questions an agent asks about one retrieval locally: walks, claims, label
routes, support profiles, comparisons, FASTA/GFA, continuations and the next request.
A retrieval made with strategy.output.coordinates (feature level 6) also says where each
run's bases lie in the indexed records (Graphlet.coordinates, coords.py).
Nothing here reads the graph; deepening is a new backend traversal that the library
prepares (next_request) and the client runs.

Stdlib only. pandas is imported lazily by frames.frames() alone, so
`import metagraph.traverse` works without it.

Every local operation takes budget= (stage L, metagraph.traverse.budget): a LocalBudget of
work units and a modelled memory account under which a call completes or stops and says
so (LocalBudgetExceeded; compare() answers comparable 'unknown'). Without one -- the
default -- nothing is budgeted and every answer is what it always was.

    from metagraph.traverse import TraverseClient, Graphlet
    resp = TraverseClient('localhost', 5555).traverse([{'sequence': seq}], strategy)
    g = resp.graphlets[0]
    for w in g.walks('right', top=5):
        print(w.path_id, w.length_bp, [l.name for l in w.labels_full])
"""

from ._codec import CodecError, GraphletFormatError, UNLIMITED
from .model import (
    AmbiguousLabel, Arm, BadSelector, Change, Claim, Comparison, Continuation, CoordClaim,
    CoordLabelWalk, CoordWalk, Graphlet, IncompatibleContinuations, IncompleteRecording,
    Label, LabelWalk, MissingEnvelope, NextRequest, Segment, Run, SupportRun, UnknownLabel,
    Walk,
)
from .coords import Coordinates, RunCoordinates, SeedOccurrences
from .parser import (FORMAT_VERSION, dump, from_response, is_canonical, load, parse, save,
                     seed_envelope, standalone_text)
from .ops import GraphletView
from .budget import (LocalBudget, LocalBudgetExceeded, LocalLimits, LocalStop, Partial,
                     local_budget)
from .client import (AttemptAnswer, AttemptAtBound, AttemptConflict, AttemptExpired,
                     AttemptSent, InstanceMismatch, ReleaseVerdict, ServerInitializing,
                     Suppression, TraverseClient, TraverseError, TraverseResponse,
                     UnsupportedFeature, release_verdict)
from .store import Entry, GraphletStore, StoreLimitExceeded, UnknownHandle

__all__ = [
    'FORMAT_VERSION', 'parse', 'dump', 'is_canonical', 'load', 'save', 'from_response',
    'standalone_text', 'seed_envelope',
    'Graphlet', 'GraphletView', 'Label', 'Arm', 'Segment', 'Run', 'Walk', 'Claim',
    'LabelWalk', 'SupportRun', 'Change', 'Continuation', 'NextRequest', 'Comparison',
    'CoordClaim', 'CoordWalk', 'CoordLabelWalk', 'Coordinates', 'RunCoordinates',
    'SeedOccurrences',
    'GraphletFormatError', 'CodecError', 'MissingEnvelope', 'AmbiguousLabel',
    'UnknownLabel', 'BadSelector', 'IncompleteRecording', 'IncompatibleContinuations',
    'UNLIMITED',
    'TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
    'AttemptAtBound', 'AttemptExpired', 'AttemptConflict', 'InstanceMismatch',
    'UnsupportedFeature', 'AttemptAnswer', 'AttemptSent', 'Suppression', 'ReleaseVerdict',
    'release_verdict', 'GraphletStore', 'Entry', 'UnknownHandle',
    'StoreLimitExceeded',
    'LocalBudget', 'LocalLimits', 'LocalStop', 'LocalBudgetExceeded', 'Partial',
    'local_budget',
]
