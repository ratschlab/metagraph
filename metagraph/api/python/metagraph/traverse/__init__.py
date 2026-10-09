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
`import metagraph.traverse` works without it. The package's names are loaded on first use
(PEP 562): `import metagraph.traverse.client` or `from metagraph.traverse import
release_verdict` loads the HTTP client and the attempt readers (attempts.py) only, not the
parser, the model or the local operations.

Every local operation takes budget= (local limits, metagraph.traverse.budget): a LocalBudget
of work units and a modelled memory account under which a call completes or stops and says
so (LocalBudgetExceeded; compare() answers comparable 'unknown'). Without one -- the
default -- nothing is budgeted and nothing is charged.

    from metagraph.traverse import TraverseClient, Graphlet
    resp = TraverseClient('localhost', 5555).traverse([{'sequence': seq}], strategy)
    g = resp.graphlets[0]
    for w in g.walks('right', top=5):
        print(w.path_id, w.length_bp, [l.name for l in w.labels_full])
"""

import importlib

# name -> the module that defines it. Loaded on first use (PEP 562): `from metagraph.traverse
# import release_verdict` (or TraverseClient) loads the light modules attempts.py and client.py
# only -- no parser, model or operation (a caller that only dispatches and releases attempts
# must not pay for the library's local processing) -- while every name below works as an
# attribute, in `from ... import`, and in `import *`.
_NAMES = {
    '_codec': ('CodecError', 'GraphletFormatError', 'UNLIMITED'),
    'model': ('AmbiguousLabel', 'Arm', 'BadSelector', 'Change', 'Claim', 'Comparison',
              'Continuation', 'CoordClaim', 'CoordLabelWalk', 'CoordWalk', 'Graphlet',
              'IncompatibleContinuations', 'IncompleteRecording', 'Label', 'LabelWalk',
              'MissingEnvelope', 'NextRequest', 'Segment', 'Run', 'SupportRun',
              'UnknownLabel', 'Walk'),
    'coords': ('Coordinates', 'RunCoordinates', 'SeedOccurrences'),
    'parser': ('FORMAT_VERSION', 'dump', 'from_response', 'is_canonical', 'load', 'parse',
               'save', 'seed_envelope', 'standalone_text'),
    'ops': ('GraphletView',),
    'budget': ('LocalBudget', 'LocalBudgetExceeded', 'LocalLimits', 'LocalStop', 'Partial',
               'local_budget'),
    'attempts': ('AttemptAnswer', 'AttemptAtBound', 'AttemptConflict', 'AttemptExpired',
                 'AttemptSent', 'InstanceMismatch', 'ReleaseVerdict', 'ServerInitializing',
                 'Suppression', 'TraverseError', 'TraverseResponse', 'classify_409',
                 'release_verdict'),
    'client': ('ResolveDeadline', 'TraverseClient', 'UnsupportedFeature'),
    'store': ('Entry', 'GraphletStore', 'StoreLimitExceeded', 'UnknownHandle'),
}
_HOME = {name: mod for mod, names in _NAMES.items() for name in names}
# the package's modules, also reachable as attributes before anything imported them
# (metagraph.traverse.ops right after `import metagraph.traverse`)
_MODULES = frozenset(('_codec', 'attempts', 'budget', 'client', 'coords', 'derive', 'export',
                      'frames', 'mcp_tools', 'model', 'ops', 'parser', 'store'))

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
    'UnsupportedFeature', 'ResolveDeadline', 'AttemptAnswer', 'AttemptSent', 'Suppression',
    'ReleaseVerdict',
    'release_verdict', 'classify_409', 'GraphletStore', 'Entry', 'UnknownHandle',
    'StoreLimitExceeded',
    'LocalBudget', 'LocalLimits', 'LocalStop', 'LocalBudgetExceeded', 'Partial',
    'local_budget',
]


def __getattr__(name):
    mod = _HOME.get(name)
    if mod is not None:
        value = getattr(importlib.import_module('.' + mod, __name__), name)
    elif name in _MODULES:
        value = importlib.import_module('.' + name, __name__)
    else:
        raise AttributeError('module %r has no attribute %r' % (__name__, name))
    globals()[name] = value             # once: later lookups never reach __getattr__
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
