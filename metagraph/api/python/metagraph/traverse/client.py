"""TraverseClient: the HTTP side of "retrieve once, process locally" (§5).

The client only issues requests (/traverse, /resolve, /traverse/capabilities); every
follow-up is answered by the graphlet locally, and deepening is a new traversal whose
request the graphlet builds (Graphlet.next_request). Stdlib only (urllib): a
requests-like `session` (anything with .request(method, url, data=, headers=,
timeout=) returning .status_code/.content/.headers) can be injected instead, e.g. a
requests.Session or a test double.

The traversal routes write compact JSON and compress it on request; the client asks
explicitly for `Accept-Encoding: gzip, deflate` (gzip preferred by the server) and
decodes either itself, so the transport does not depend on the HTTP library.

A request may be an ATTEMPT of a ledger that reserved an allowance for it (DESIGN §14):
`attempt_id` (unique per server process while it runs or is retained), with `budget_id`
and `locus_id` echoed. Every response to it carries a `usage` block
(TraverseResponse.usage), the server enforces a duration bound on it (usage.bound), and it
can be cancelled (cancel()) and its state read (attempt()) by id, also from another client.

Feature level 5 (SPEC §5, §10.3): a cancel names the cancelled request's `not_after_ms`,
so the tombstone it leaves is held until no copy of that request can start
(`covers_admission`); a request names the `server_instance` it is meant for
(`expect_server_instance`), so a copy that reaches a restarted process is refused (409
`instance_mismatch`); and the capabilities' `release_rule` says when a ledger may give an
attempt's capacity back -- release_verdict() applies it to an answer. Both fields are sent
only to a server that states `feature_level` >= 5 (read from GET /capabilities, or given to
the client): a server below it refuses them with a 400 after it registered the
attempt_id, so a request with that id is refused (409) while that server retains it
(retention_s, retention_count; below level 5 nothing holds it longer) and runs again after
that.

Feature level 6 (DESIGN §18): record coordinates. traverse() asks for them
(strategy.output.coordinates) by itself for a support: trace strategy that does not set
them, where GET /traverse/capabilities states the coordinates block for the index with
supported true (supports_coordinates()) -- never on an index that cannot report them,
where the null form would only cost depth. Under a request memory budget
(bounds.max_memory_mb) only while the depth gate passes
(AUTO_COORDINATES_UNDER_MEMORY_BUDGET): a run's coordinate account is charged with the
walk, so a budgeted walk stops a little shallower with them; when it does not pass the
response says so (coordinates_auto, notes) and coordinates=True asks for them, and when
coordinates asked for that way share a memory stop the notes say how to drop them. It never
retries silently without them: an older server's 400 is raised as it is. A request whose
output.coordinates is not true never carries output.max_coordinate_occurrences, which the
server refuses without it.

What GET /capabilities states is kept per server PROCESS: it is forgotten when an answer
names another server_instance than the one read (a restart at the same address), on a 409
instance_mismatch, and on a 400 that refuses a level-5 field as unknown (an older binary
now answers there: the level read is no longer the process's). A cancel refused that way
is sent again without not_after_ms, since that 400 stopped nothing. A request with
expect_server_instance that is the first to reach a replaced process cannot be saved that
way: an older binary registers its attempt_id before it answers the 400 (the id is
refused while that server retains it, and runs again after that), and a level-5 one
refuses it, 409 instance_mismatch, without registering it.

A multi-graph server answers /traverse and the capabilities probe of one pair in the envelope
(SPEC-multigraph-envelope.md section 1): the client reads its one entry as the single-graph
answer (single_answer()), with the attempt's usage the envelope states at its top level. Its
other answers (409, 503, 400) are not wrapped and read as before.
"""

import copy
import gzip
import hashlib
import http.client
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import uuid
import zlib
from dataclasses import dataclass
from typing import NamedTuple, Optional

# the answers' readers and the release rule live in the light module attempts.py (importing
# it loads no model): re-exported here as the same objects
from .attempts import (AttemptAnswer, AttemptAtBound, AttemptConflict, AttemptExpired,
                       AttemptSent, InstanceMismatch, ReleaseVerdict, ServerInitializing,
                       Suppression, TraverseError, TraverseResponse, classify_409,
                       release_verdict)

__all__ = ['TraverseClient', 'TraverseResponse', 'TraverseError', 'ServerInitializing',
           'AttemptAtBound', 'AttemptExpired', 'AttemptConflict', 'InstanceMismatch',
           'UnsupportedFeature', 'ResolveDeadline', 'ResponseFormatError', 'AttemptAnswer',
           'AttemptSent',
           'Suppression', 'ReleaseVerdict', 'release_verdict', 'classify_409', 'ACCEPT_ENCODING',
           'SUPPRESSION_LEVEL', 'COORDINATES_LEVEL', 'AUTO_COORDINATES_UNDER_MEMORY_BUDGET',
           'auto_coordinates', 'drop_coordinates_note', 'strip_coordinate_cap',
           'SELECTORS', 'pair_members', 'answered_pair', 'requested_pair', 'PairResults',
           'PAIR_ID_NAMESPACE', 'configured_pair_id', 'check_pair_id', 'Pair', 'PairIdError',
           'capabilities_pair_resolver', 'request_without_paths', 'entry_without_paths']

ACCEPT_ENCODING = 'gzip, deflate'
# the feature level of the tombstone's suppression (a cancel with not_after_ms,
# covers_admission) and of expect_server_instance (SPEC §10.3)
SUPPRESSION_LEVEL = 5
# the feature level whose GET /traverse/capabilities states the coordinates block (SPEC
# §10.3): supports_coordinates() reads the block itself, not the level
COORDINATES_LEVEL = 6
# The depth gate: under a request memory budget the library asks for
# coordinates automatically only if the median depth at the budget's stop with them stays
# within 10% of the depth without -- each run's coordinate entry and occurrences are
# charged when the run is created, so a budgeted walk with them stops shallower. Measured
# at level 6 (DESIGN §18.8): on mini_refseq's 56 trace cells the median depth is 0.959-0.971
# of the depth without (a 3-4% drop), but the gate FAILS where it matters most: column labels
# with 16 or more chains per run (refseq33m's taxid columns) and seeds with many header
# labels need 2-4x the memory and can fail at depth 0 at a budget that completes without
# them (the regime test). The gate is therefore not passed, and the library does not ask
# under a memory budget: the answer says so and how to ask (_MEMORY_BUDGET_NOTE,
# coordinates=True)
AUTO_COORDINATES_UNDER_MEMORY_BUDGET = False
_MEMORY_BUDGET_NOTE = (
    'record coordinates were not requested automatically: the request has a memory budget '
    '(bounds.max_memory_mb), under which their account makes the walk stop shallower; '
    'pass coordinates=True to ask for them (output.coordinates), or drop the memory budget')
# where coordinates requested automatically under a memory budget share a memory stop
# (DESIGN §18.8: the single-lineage switch cells are the outliers; at 50% of their peak two
# of them fail at depth 0 with coordinates and not without)
_DROP_COORDINATES_NOTE = (
    'record coordinates were requested automatically under a memory budget '
    '(bounds.max_memory_mb), and %d seed(s) stopped on it with drop_coordinates offered: '
    'their coordinate account is part of what stopped them; send the request again with '
    'coordinates=False (output.coordinates false) to walk further, or raise the budget')
# A server refuses a field it does not know with a 400 naming it ("request: unknown field
# 'not_after_ms'", server.cpp and traverse.cpp). For the fields the client sends by the
# server's level, that 400 says the level the client kept is not the process's any more.
_UNKNOWN_FIELD = re.compile(r"unknown field '([^']*)'")
_LEVELLED_FIELDS = frozenset(('not_after_ms', 'expect_server_instance'))


# the largest max_coordinate_occurrences a server accepts (2^64 - 2: the wire contract)
_MAX_CAP = (1 << 64) - 2

# What names the pair a /traverse request runs on (SPEC-multigraph-envelope.md section 6):
# request-level members, never strategy fields. `annotation_path` and `pair_id` are no members
# of a request: /traverse takes them inside a `pairs` object only (pair_members() writes them
# there); as overrides (deepen(), next_request()) they name the pair as pair_members() does
SELECTORS = ('graph', 'graph_path', 'annotation_path', 'pair_id', 'graphs', 'pairs')

# The n1 id of a configured pair (docs/DESIGN-pair-identity.md section 3.2): master's
# src/common/configured_pair_id.hpp (356fe1ee3), which the engine and the service compute too
# and which must stay byte-identical in all three. Its namespace:
PAIR_ID_NAMESPACE = uuid.UUID('d3565c73-e7e9-5bca-b6a9-64cb5c4c9f93')


def configured_pair_id(graph, graph_path, annotation_path):
    """The n1 id of the pair the graph list configures as (|graph|, |graph_path|,
    |annotation_path|): "n1-" and the lowercase UUIDv5 in PAIR_ID_NAMESPACE of the three
    fields, each as its UTF-8 byte length, ":" and its bytes, an absent one (None) as "-1:"
    (so None and "" differ). Nothing is normalized: the strings as the graph list spells
    them. A configured identity, not a content one (index_fp is that): the name is hashed,
    so a pair listed under two names has two ids. Public metadata, never a credential."""
    name = bytearray()
    for v in (graph, graph_path, annotation_path):
        if v is None:
            name += b'-1:'
        elif isinstance(v, str):
            b = v.encode('utf-8')
            name += b'%d:' % len(b) + b
        else:
            raise TypeError('a pair field is a string or None, not %r' % (v,))
    # uuid.uuid5() itself, over bytes (3.11's uuid5 takes a str name only)
    digest = hashlib.sha1(PAIR_ID_NAMESPACE.bytes + bytes(name), usedforsecurity=False).digest()
    return 'n1-' + str(uuid.UUID(bytes=digest[:16], version=5))


def _text(v):
    return v if isinstance(v, str) and v else None


class Pair(NamedTuple):
    """A pair of a multi-graph server as the library reads it (docs/DESIGN-pair-identity.md
    section 3.9): its graph name and its paths (a document stored before the engine named
    pairs by id; a path selector), or its n1 id (a document stored after; an id selector), or
    both -- |pair_id| is then the n1 id of the paths (configured_pair_id()), so that an old
    document and its id-form twin name one pair. A path or an id not stated is None.
    answered_pair(), requested_pair() and TraverseClient.selection() return one;
    pair_members(*pair) writes its selector."""
    graph: Optional[str]
    graph_path: Optional[str] = None
    annotation_path: Optional[str] = None
    pair_id: Optional[str] = None

    @property
    def identity(self):
        """What names the pair whatever its form: (graph, pair_id) where it has an id (stated,
        or computed from both paths), else (graph, graph_path, annotation_path)."""
        if self.pair_id:
            return self.graph, self.pair_id
        return self.graph, self.graph_path, self.annotation_path

    def named(self):
        """The pair as an output names it (docs/DESIGN-pair-identity.md section 2: no path in
        any response): {graph, pair_id} where it has an id -- an old document's computed from
        its paths -- else {graph} with the paths it was given (a caller's selector without
        the annotation path, which has no id to compute)."""
        if self.pair_id:
            return {'graph': self.graph, 'pair_id': self.pair_id}
        return {k: v for k, v in zip(('graph', 'graph_path', 'annotation_path'), self) if v}


# an n1 id as the engine states it: "n1-" and a lowercase UUID of version 5 (RFC 4122 variant)
_N1 = re.compile(r'n1-[0-9a-f]{8}-[0-9a-f]{4}-5[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}\Z')


class PairIdError(ValueError):
    """A pair id the library does not take, or a pair named by its id that cannot be sent.
    |code|: 'p2_pair_id' (the service's public p2 id where the engine's n1 id belongs: a
    document and a selector carry n1, never p2), 'bad_pair_id' (not an n1 id; an id beside
    paths in one selector), 'pair_id_mismatch' (paths stated with an id that is not theirs,
    or a resolver's answer whose id is not the one it resolved), 'pair_unresolved' (a pair
    named by its id only, to a server that does not state multigraph.pair_id "n1", with no
    resolver: TraverseClient(pair_resolver=)). Raised before anything is sent. A
    ValueError, so that a caller mapping ValueErrors to bad arguments still refuses it."""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


def _shown(v):
    r = repr(v)
    return r if len(r) <= 80 else r[:77] + '...'


def check_pair_id(v, where='pair_id'):
    """|v| when it is an n1 id, else PairIdError naming |where|: 'p2_pair_id' for the
    service's p2 id, 'bad_pair_id' for anything else."""
    if isinstance(v, str) and _N1.match(v):
        return v
    if isinstance(v, str) and v.startswith('p2-'):
        raise PairIdError('p2_pair_id', '%s %s is the service\'s public pair id (p2), not the '
                          'engine\'s n1 id: a document and a selector name a pair by its n1 '
                          'id; map the p2 id to its n1 id where it was issued' % (where, _shown(v)))
    raise PairIdError('bad_pair_id', '%s is not an n1 pair id ("n1-" and a lowercase '
                      'UUIDv5): %s' % (where, _shown(v)))


def _pair(graph, graph_path, annotation_path, pair_id, where):
    """Pair of tags as stated: an id checked (check_pair_id()), computed from both paths
    where none is stated, and refused (pair_id_mismatch) where one is stated beside both
    paths and is not theirs."""
    gp, ap = _text(graph_path), _text(annotation_path)
    pid = None if pair_id in (None, '') else check_pair_id(pair_id, where + ' pair_id')
    if gp is not None and ap is not None:
        own = configured_pair_id(graph, gp, ap)
        if pid is not None and pid != own:
            raise PairIdError('pair_id_mismatch', '%s states pair_id %s beside paths whose n1 '
                              'id is %s' % (where, pid, own))
        pid = own
    return Pair(graph, gp, ap, pid)


def pair_members(graph=None, graph_path=None, annotation_path=None, pair_id=None):
    """The request members that name one pair of a multi-graph server: `graph` (and
    `graph_path`), or -- with |annotation_path|, which /traverse reads inside a pair object
    only -- `pairs: [{graph, graph_path?, annotation_path}]` (`graph` beside `pairs` is
    refused, SPEC-multigraph-envelope.md section 6), or -- with |pair_id| and no path --
    `pairs: [{graph, pair_id}]` (docs/DESIGN-pair-identity.md section 3.3; TraverseClient
    sends it as its pair_resolver or the server allows). Both paths with their own id
    (Pair of an old document) write the paths: path selectors stay valid input on every
    engine. {} without a graph: the request names no pair (a single-graph server, or a
    multi-graph server's only pair). A path or an id without a graph is a ValueError: the
    server needs the name; an id beside one path, or beside paths that are not its,
    PairIdError."""
    if not graph:
        if graph_path or annotation_path or pair_id:
            raise ValueError('graph_path, annotation_path and pair_id name a pair of a graph: '
                             'give the graph too')
        return {}
    if pair_id:
        pid = check_pair_id(pair_id)
        if not graph_path and not annotation_path:
            return {'pairs': [{'graph': graph, 'pair_id': pid}]}
        if not (graph_path and annotation_path):
            raise PairIdError('bad_pair_id', 'name a pair by its pair_id or by its paths, '
                              'not by both')
        if configured_pair_id(graph, graph_path, annotation_path) != pid:
            raise PairIdError('pair_id_mismatch', 'pair_id %s is not the n1 id of the paths '
                              'given with it' % pid)
    if annotation_path:
        pair = {'graph': graph}
        if graph_path:
            pair['graph_path'] = graph_path
        pair['annotation_path'] = annotation_path
        return {'pairs': [pair]}
    out = {'graph': graph}
    if graph_path:
        out['graph_path'] = graph_path
    return out


def answered_pair(envelope):
    """Pair(graph, graph_path, annotation_path, pair_id) of the pair a retrieval was answered
    on, from the envelope its graphlet keeps: a multi-graph server's entry carries its tags
    (single_answer() keeps them) -- the paths before the engine named pairs by id, `pair_id`
    (n1) after; an old document's pair_id is computed from its paths. None for a
    single-graph server's answer, which names no pair. A tag the entry does not state is
    None; a pair_id that is no n1 id (a p2 id, say) is PairIdError."""
    env = envelope if isinstance(envelope, dict) else {}
    graph = _text(env.get('graph'))
    if graph is None:
        return None
    return _pair(graph, env.get('graph_path'), env.get('annotation_path'), env.get('pair_id'),
                 'the envelope')


def requested_pair(req):
    """Pair(graph, graph_path, annotation_path, pair_id) of the pair a /traverse request
    names -- its `graph` and `graph_path`, or its one `pairs` object (by paths or by
    pair_id) -- or None (no pair named, or `graphs`, or several pairs: the server
    decides). A pairs object with both paths gets their id, as answered_pair()."""
    req = req if isinstance(req, dict) else {}
    pairs = req.get('pairs')
    if pairs is not None:
        if isinstance(pairs, list) and len(pairs) == 1 and isinstance(pairs[0], dict) \
                and _text(pairs[0].get('graph')):
            p = pairs[0]
            return _pair(p['graph'], p.get('graph_path'), p.get('annotation_path'),
                         p.get('pair_id'), 'request.pairs[0]')
        return None
    graph = _text(req.get('graph'))
    return None if graph is None else Pair(graph, _text(req.get('graph_path')))


def request_without_paths(request):
    """|request| (a /traverse request, a replay or continuation request) as an output states
    it: each `pairs` object that names its pair by both paths names it by {graph, pair_id}
    instead (the n1 id of those paths), in a copy -- a request the library sends keeps the
    paths of an old document, which every engine takes, and sent again through a client this
    one goes through wire_request(). A selector without the annotation path (a caller's, with
    no id to compute) and any other request are returned as they are."""
    pairs = request.get('pairs') if isinstance(request, dict) else None
    if not isinstance(pairs, list):
        return request
    out = []
    for i, p in enumerate(pairs):
        if isinstance(p, dict) and _text(p.get('graph')) and _text(p.get('graph_path')) \
                and _text(p.get('annotation_path')):
            t = _pair(p['graph'], p['graph_path'], p['annotation_path'], p.get('pair_id'),
                      'request.pairs[%d]' % i)
            rest = {k: v for k, v in p.items() if k not in Pair._fields}
            out.append(dict(t.named(), **rest))
        else:
            out.append(p)
    return dict(request, pairs=out)


def entry_without_paths(entry):
    """A multi-graph entry's tags (an answer's, a probe's) as an output states them: graph
    and pair_id, no path -- the pair_id the n1 id of the paths an engine before ids states.
    Anything else is returned as it is."""
    if not isinstance(entry, dict) or not (_text(entry.get('graph_path'))
                                           or _text(entry.get('annotation_path'))):
        return entry
    pair = answered_pair(entry)
    out = {k: v for k, v in entry.items() if k not in ('graph_path', 'annotation_path')}
    if pair is not None and pair.pair_id:
        out['pair_id'] = pair.pair_id
    return out


def capabilities_pair_resolver(capabilities):
    """A pair_resolver (TraverseClient(pair_resolver=)) from a GET /capabilities document
    that lists its pairs with their paths (graph_summary.pairs[], a server before ids): n1
    id -> (graph, graph_path, annotation_path) as that server spells them, None for an id it
    does not list. A pair listed by its id only (a server after ids) is not in it: such a
    server takes the id itself."""
    summary = (capabilities or {}).get('graph_summary') if isinstance(capabilities, dict) \
        else None
    table = {}
    for i, p in enumerate((summary or {}).get('pairs') or []
                          if isinstance(summary, dict) else []):
        if not isinstance(p, dict) or not _text(p.get('graph')):
            continue
        t = _pair(p['graph'], p.get('graph_path'), p.get('annotation_path'), p.get('pair_id'),
                  'GET /capabilities graph_summary.pairs[%d]' % i)
        if t.graph_path and t.annotation_path:
            table[t.pair_id] = (t.graph, t.graph_path, t.annotation_path)
    return table.get


def strip_coordinate_cap(strategy):
    """|strategy| (in place) without output.max_coordinate_occurrences unless its
    output.coordinates is true: the server refuses the cap without coordinates (400), so
    the library never sends it alone -- build_request(), next_request() (whose request
    carries the retrieval's echo, cap included) and the tools' requests (so that dropping
    coordinates through a continuation works). Defined here, not in coords.py (which
    re-exports it), so that importing the client loads no model."""
    out = strategy.get('output') if isinstance(strategy, dict) else None
    if isinstance(out, dict) and out.get('coordinates') is not True:
        out.pop('max_coordinate_occurrences', None)
    return strategy


def _valid_cap(v):
    if v == 'unlimited':
        return True
    return isinstance(v, int) and not isinstance(v, bool) and 1 <= v <= _MAX_CAP


def _memory_budget(strategy):
    v = ((strategy or {}).get('bounds') or {}).get('max_memory_mb')
    return isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0


def auto_coordinates(client, strategy, graph=None, graph_path=None, annotation_path=None,
                     pair_id=None):
    """The automatic rule for record coordinates -> (value, statement):
    value True (ask), or None (leave the strategy as given); statement None when the rule
    does not apply (the strategy sets output.coordinates, or is not support: trace in
    constrain mode), else {requested, reason[, memory_budget][, note]} -- reason
    'supported' (asked; memory_budget true when under a request memory budget, the depth gate
    passing), 'unsupported' (the index states no coordinates block, or supported false) or
    'memory_budget' (a request memory budget with the gate not passing: not asked, |note|
    says how to ask). |client|: anything with supports_coordinates() (a TraverseClient);
    one without asks for none. |graph|, |graph_path|, |annotation_path|, |pair_id|: the
    pair whose capabilities decide (the client's own when none is given)."""
    s = strategy or {}
    out = s.get('output') if isinstance(s.get('output'), dict) else {}
    if 'coordinates' in out:
        return None, None
    labels = s.get('labels') if isinstance(s.get('labels'), dict) else {}
    if s.get('support', 'kmer') != 'trace' or labels.get('mode', 'constrain') != 'constrain':
        return None, None
    probe = getattr(client, 'supports_coordinates', None)
    # the probe gets the keywords it needs only (a test double may take fewer)
    kw = {k: v for k, v in (('annotation_path', annotation_path), ('pair_id', pair_id))
          if v is not None}
    if probe is None or not probe(graph, graph_path, **kw):
        return None, {'requested': False, 'reason': 'unsupported'}
    if _memory_budget(s):
        if not AUTO_COORDINATES_UNDER_MEMORY_BUDGET:
            return None, {'requested': False, 'reason': 'memory_budget',
                          'note': _MEMORY_BUDGET_NOTE}
        return True, {'requested': True, 'reason': 'supported', 'memory_budget': True}
    return True, {'requested': True, 'reason': 'supported'}


def drop_coordinates_note(results, statement):
    """The note for coordinates the automatic rule asked for under a memory budget
    (|statement|: auto_coordinates()'s) where seed results of |results| (the response's
    raw `results`, graphlet or failed) stopped on memory with drop_coordinates offered,
    else None: it says how to walk further without them (DESIGN §18.8)."""
    if not statement or not statement.get('requested') or not statement.get('memory_budget'):
        return None
    n = 0
    for r in results or ():
        stop = r.get('resource_stop') if isinstance(r, dict) else None
        if isinstance(stop, dict) and 'drop_coordinates' in (stop.get('actions') or ()):
            n += 1
    return _DROP_COORDINATES_NOTE % n if n else None


def single_answer(out):
    """|out| as a single-graph server answers it: a multi-graph server's envelope of ONE pair
    (SPEC-multigraph-envelope.md section 1: `envelope_version` and one entry in `answers`) is
    read as its entry -- the pair's answer with its tags (graph, graph_path, annotation_path,
    index_fp, outcome) -- plus the envelope's `usage` (on /traverse the attempt's, which the
    envelope states at its top level); an entry refused raises TraverseError with its
    refusal's http_status and error, a refusal that is not an object ResponseFormatError.
    Any other answer (a single-graph server's, an envelope over several pairs) is returned
    as it is."""
    if not isinstance(out, dict) or 'envelope_version' not in out:
        return out
    answers = out.get('answers')
    if not isinstance(answers, list) or len(answers) != 1 or not isinstance(answers[0], dict):
        return out
    entry = answers[0]
    if entry.get('outcome') == 'refused':
        refusal = entry.get('refusal') or {}
        if not isinstance(refusal, dict):
            # not shaped as a refusal (Opus fix review N3)
            raise ResponseFormatError('the refused entry\'s refusal is a JSON %s, not an '
                                      'object' % _json_kind(refusal), out)
        raise TraverseError(refusal.get('http_status'), refusal.get('error'), refusal)
    answer = dict(entry)
    if 'usage' in out:
        answer['usage'] = out['usage']
    return answer


def _json_kind(v):
    return {list: 'array', dict: 'object', str: 'string', bool: 'boolean',
            type(None): 'null'}.get(
        type(v), 'number' if isinstance(v, (int, float)) else type(v).__name__)


def _unknown_field(status, message):
    """The field a 400's |message| refuses as unknown, else None."""
    if status != 400 or not isinstance(message, str):
        return None
    m = _UNKNOWN_FIELD.search(message)
    return m.group(1) if m else None


class ResolveDeadline(TraverseError):
    """503 with `code: "deadline"` and no `usage`: a /resolve sent with bounds.time_budget_ms
    (resolve(time_budget_ms=)) whose answer could not be built and written within that budget
    -- its finalisation overran the reserve (SPEC §4.5). Nothing partial was sent and nothing
    was registered; send it again with a larger budget or a shorter sequence. Not a loading
    server (whose 503 has Retry-After and no code): never retry it as one."""


class ResponseFormatError(TraverseError):
    """A 2xx answer whose JSON is not shaped as a /traverse response (not an object, its
    `results` not a list of objects): |status| None, |message| what is wrong, |body| the
    answer. A TraverseError, so traverse_pairs() reports it as that pair's {pair, error}."""

    def __init__(self, message, body=None):
        super().__init__(None, message, body)
        self.args = ('not a /traverse response: %s' % message,)


class UnsupportedFeature(ValueError):
    """Raised before sending: the request asks for a field the server does not state (its
    |feature_level| is below |needed|). Nothing was sent -- a server below the level refuses
    the field with a 400 after it registered the attempt_id of the request, whose id it
    then refuses while it retains it."""

    def __init__(self, feature, feature_level, needed):
        self.feature = feature
        self.feature_level = feature_level
        self.needed = needed
        super().__init__('%s needs a server at feature_level >= %d; this one states %s'
                         % (feature, needed, feature_level))


class _ResponseCut(ConnectionError):
    """No whole answer was received: the response was cut in transfer -- the connection
    closed before its Content-Length (http.client.IncompleteRead: a server drops a
    connection at its content timeout, a process dies with its buffers unsent, a proxy cuts
    it) -- or its status line or a header could not be read (another
    http.client.HTTPException). A ConnectionError, so an OSError like every other transport
    failure: callers that map those (the tools' backend_unreachable) map this one too, never
    a bare HTTPException that no handler expects. Not a TraverseError of the status, even
    when one was read: a cut 2xx is no server answer, and a cut 409 or 503 body no longer
    tells its cases apart -- for an attempt the conservative reading is "unanswered" (cancel
    it, then judge it with release_verdict()). |reason| is the message (what the tools
    report)."""

    def __init__(self, message, status=None):
        super().__init__(message)
        self.reason = message
        self.status = status


def _cut(method, path, e, status=None, attempt_id=None):
    """The _ResponseCut of http.client.HTTPException |e| raised while |method| |path| was
    answered (|status|: the HTTP status, when its line was read; |attempt_id|: the
    request's, when it carried one). The message says the server was reached and may have
    run the request: the request went out before any answer was read, so a reader of
    "unreachable" alone would take a retry for free."""
    what = 'the response to %s %s' % (method, path) if status is None else \
        'the HTTP %d answer to %s %s' % (status, method, path)
    if isinstance(e, http.client.IncompleteRead):
        got = len(e.partial or b'')
        more = ', %d more expected' % e.expected if e.expected is not None else ''
        detail = 'cut in transfer (%d bytes read%s)' % (got, more)
    else:
        detail = 'unreadable (%s: %s)' % (type(e).__name__, e)
    then = 'cancel attempt %r, then judge it with release_verdict()' % attempt_id \
        if attempt_id is not None else 'a retry sends it again'
    return _ResponseCut('%s was %s: no whole answer was received, but the server was '
                        'reached and may have run the request (%s)' % (what, detail, then),
                        status)


def _decoded(data, encoding, status):
    """_decode_body(), a body its Content-Encoding does not decode (a truncated gzip or
    deflate stream) a TraverseError of the answer's status: zlib.error and EOFError are no
    error a caller of the client expects, and gzip's BadGzipFile, an OSError, read as a
    transport failure (backend_unreachable in the tools)."""
    try:
        return _decode_body(data, encoding)
    except (zlib.error, EOFError, OSError) as e:
        raise TraverseError(status, 'the response body could not be decoded (Content-Encoding '
                            '%s: %s)' % (encoding, type(e).__name__)) from None


def _decode_body(data, encoding):
    encoding = (encoding or '').strip().lower()
    if encoding == 'gzip':
        return gzip.decompress(data)
    if encoding == 'deflate':
        try:
            return zlib.decompress(data)            # zlib stream (what the server writes)
        except zlib.error:
            return zlib.decompress(data, -zlib.MAX_WBITS)   # raw deflate (RFC-lenient)
    return data


class TraverseClient:
    def __init__(self, host, port, api_path=None, *, session=None, timeout=900, release=None,
                 scheme='http', graph=None, feature_level=None, graph_path=None,
                 annotation_path=None, pair_id=None, pair_resolver=None):
        """|graph|, |graph_path|, |annotation_path|: the client's own pair of a multi-graph
        server, named by every request and probe that names none (selection()); a graph
        with several graphs needs graph_path, a graph whose pairs share one graph
        annotation_path (sent in a `pairs` object: pair_members()). |pair_id|: that pair
        named by its n1 id instead of its paths (docs/DESIGN-pair-identity.md section 3.3).

        |pair_resolver|: n1 id -> (graph, graph_path, annotation_path) of that pair as the
        server's graph list spells it, or None for an id it does not know -- the caller's
        binding (the service builds it from the engine's GET /capabilities;
        capabilities_pair_resolver() does that for a server that lists paths). Every request
        and probe this client sends that names a pair by its id only (the continuation,
        replay or re-ask of a document stored with a pair_id and no paths) passes
        wire_request() / wire_pair(): with a resolver the id is resolved, the n1 id of the
        answer recomputed and compared (PairIdError pair_id_mismatch otherwise), and the pair
        sent by its paths, which every server takes; without one it is sent by its id where
        the server states multigraph.pair_id "n1" in GET /capabilities, and refused
        (PairIdError pair_unresolved) where it does not. Nothing is sent before that.

        |feature_level|: the server's, when the caller knows it; else it is read from GET
        /capabilities the first time a feature-level-5 field is to be sent (and never for a
        request without one, so a request without level-5 fields reads no capabilities).
        Once GET /capabilities has been read (expect_server_instance='auto' reads it for the
        instance) the level it states wins over the given one, and a given level of 5 or
        more is dropped when the server refuses a level-5 field as unknown: either is the
        process's own statement."""
        self.host = host
        self.port = port
        self.base = '%s://%s:%s' % (scheme, host, port)
        if api_path:
            self.base += '/' + api_path.strip('/')
        self.session = session
        self.timeout = timeout
        self.release = release
        self.graph = graph
        self.graph_path = graph_path
        self.annotation_path = annotation_path
        self.pair_id = None if pair_id is None else check_pair_id(pair_id)
        self.pair_resolver = pair_resolver
        self.feature_level = feature_level
        self._server_doc = None
        # Pair (selection()) -> whether GET /traverse/capabilities states coordinates
        self._coord_support = {}

    # ---------------------------------------------------------------- transport

    def _request(self, method, path, payload=None, answers=(), with_status=False):
        """|answers|: non-2xx statuses whose JSON body is an answer, returned rather than
        raised (cancel and attempt: 404 states that nothing runs under the id).
        |with_status|: -> (HTTP status, body)."""
        url = self.base + path
        body = None if payload is None else json.dumps(payload).encode('utf-8')
        headers = {'Accept-Encoding': ACCEPT_ENCODING}
        if body is not None:
            headers['Content-Type'] = 'application/json'
        if self.session is not None:
            resp = self.session.request(method, url, data=body, headers=headers,
                                        timeout=self.timeout)
            status = resp.status_code
            hdrs = {k.lower(): v for k, v in dict(resp.headers or {}).items()}
            data = resp.content
            # requests decodes Content-Encoding itself; a raw double does not
            if getattr(resp, 'raw_encoded', False):
                data = _decoded(data, hdrs.get('content-encoding'), status)
        else:
            req = urllib.request.Request(url, data=body, headers=headers, method=method)
            status = None
            try:
                try:
                    with urllib.request.urlopen(req, timeout=self.timeout) as r:
                        status = r.status
                        hdrs = {k.lower(): v for k, v in r.headers.items()}
                        raw = r.read()
                    data = _decoded(raw, hdrs.get('content-encoding'), status)
                except urllib.error.HTTPError as e:
                    status = e.code
                    hdrs = {k.lower(): v for k, v in (e.headers or {}).items()}
                    try:
                        raw = e.read() or b''
                    finally:
                        e.close()
                    data = _decoded(raw, hdrs.get('content-encoding'), status)
            except http.client.HTTPException as e:
                # a body cut short (IncompleteRead) or an unreadable status line or header:
                # neither a TraverseError nor an OSError, so it would escape every caller's
                # mapping. The session= path needs no such case: requests raises its
                # ChunkedEncodingError (an OSError) for the same cut
                raise _cut(method, path, e, status,
                           payload.get('attempt_id') if isinstance(payload, dict) else None) \
                    from e
        out = None
        if isinstance(data, (bytes, bytearray)):
            try:
                text = data.decode('utf-8')
            except UnicodeDecodeError:
                # not UTF-8, so no JSON answer: an error page is still classified by its
                # status below, a 2xx is 'not JSON' -- never a bare UnicodeDecodeError, which
                # no caller expects from a server's answer
                text = data.decode('utf-8', 'replace')
                data = None
        else:
            text = data
        if data is not None:
            try:
                out = json.loads(text) if text else {}
            except ValueError:
                out = None
        if isinstance(out, dict) and self._server_doc and path != '/capabilities':
            self._check_instance(out)
        if status in answers and isinstance(out, dict):
            return (status, out) if with_status else out
        if not 200 <= status < 300:
            message = out.get('error') if isinstance(out, dict) else (text or '')[:500]
            refused = _unknown_field(status, message)
            if status == 400 and isinstance(message, str) and 'coordinate' in message:
                # a server that refuses the coordinates fields is not the one whose
                # capabilities were read (a rollback at the same address): read them again
                self._coord_support.clear()
            if refused in _LEVELLED_FIELDS and isinstance(payload, dict) and refused in payload:
                # a field sent by the server's level is unknown to the process answering:
                # an older binary replaced the one whose level the client kept (a rollback
                # at the same address; it cannot answer instance_mismatch). Kept, the level
                # would send the field again: a cancel that never stops anything, an 'auto'
                # request whose attempt_id that binary registers and answers with a 400
                # every time
                self._forget_level()
            if status == 503:
                # Two 503s mean opposite things to a ledger: a loading index or a route at its
                # cap (Retry-After, no usage: nothing registered; the latter with code "busy",
                # E19) and an attempt stopped at its bound (usage, no
                # Retry-After: registered and consumed -- a retry under its id is refused
                # while the server holds the id, and runs again after). Told apart by the body,
                # not by the status alone
                retry = hdrs.get('retry-after')
                if isinstance(out, dict) and 'usage' in out and not retry:
                    raise AttemptAtBound(status, message, out)
                if isinstance(out, dict) and out.get('code') == 'deadline' and not retry:
                    # a /resolve past its own bounds.time_budget_ms (SPEC §4.5): answered, not
                    # loading, so a ServerInitializing would invite a retry of the same budget
                    raise ResolveDeadline(status, message, out)
                raise ServerInitializing(status, message, out,
                                         int(retry) if retry and retry.isdigit() else None)
            if status == 409 and isinstance(out, dict):
                # three refusals, none of which ran anything: told apart by the body
                # (classify_409(): the attempt object's state first, an instance_mismatch
                # only to a request that named an instance)
                if out.get('state') == 'instance_mismatch':
                    # the process this client knew is gone: read its successor's instance
                    # before an 'auto' request names it again (forgetting it costs one read)
                    self._server_doc = None
                raise classify_409(out, payload if isinstance(payload, dict) else None,
                                   message)
            raise TraverseError(status, message, out)
        if out is None:
            raise TraverseError(status, 'the response is not JSON', text[:500])
        return (status, out) if with_status else out

    # ---------------------------------------------------------------- feature level

    def _server(self, refresh=False):
        if self._server_doc is None or refresh:
            self._server_doc = self.server_capabilities()
        return self._server_doc

    def _check_instance(self, out):
        """Forget the kept GET /capabilities when an answer names another server_instance:
        the process was restarted (or replaced by another binary) at the same address, and
        its instance and maybe its level are not the ones read. Every answer about an
        attempt names its process (a cancel's, a state's, a usage block, a 409's attempt),
        so a restart is noticed before the next 'auto' request names the old instance."""
        known = (self._server_doc.get('attempts') or {}).get('server_instance')
        if known is None:
            return
        inst = out.get('server_instance')
        if inst is None:
            for key in ('usage', 'attempt'):
                sub = out.get(key)
                if isinstance(sub, dict) and sub.get('server_instance') is not None:
                    inst = sub['server_instance']
                    break
        if inst is not None and inst != known:
            self._server_doc = None

    def _forget_level(self):
        self._server_doc = None
        if self.feature_level is not None and self.feature_level >= SUPPRESSION_LEVEL:
            self.feature_level = None

    def server_feature_level(self, refresh=False):
        """The server's feature_level: GET /capabilities' once read, else the one given to
        the client, else GET /capabilities' read now (|refresh| reads it again). Kept per
        server process: forgotten when an answer names another server_instance, on an
        InstanceMismatch and when the server refuses a level-5 field as unknown (see the
        module's docstring). A server without that route (before feature level 1) is 0."""
        if not refresh:
            if self._server_doc is not None:
                return int(self._server_doc.get('feature_level') or 0)
            if self.feature_level is not None:
                return self.feature_level
        try:
            doc = self._server(refresh=True)
        except TraverseError as e:
            if e.status != 404:
                raise
            self._server_doc = doc = {}         # kept: no route is the answer
        return int(doc.get('feature_level') or 0)

    def server_instance(self, refresh=False):
        """The server process's `server_instance` (16 hex, new at every start) from GET
        /capabilities' attempts block, read once and kept for that process (forgotten as
        server_feature_level() is, or read again with |refresh|); None where the server
        states none."""
        return ((self._server(refresh) or {}).get('attempts') or {}).get('server_instance')

    def refresh_capabilities(self):
        """Read GET /capabilities again (a restarted server: a new server_instance, maybe
        another feature level)."""
        return self._server(refresh=True)

    def _level(self, needed):
        level = self.server_feature_level()
        return level, level >= needed

    # ---------------------------------------------------------------- routes

    def selection(self, graph=None, graph_path=None, annotation_path=None, pair_id=None):
        """Pair(graph, graph_path, annotation_path, pair_id): the pair a request or probe
        names -- the arguments (|graph| the client's own when only the paths or the id are
        given), or, with none given, the client's own pair."""
        if graph is None and graph_path is None and annotation_path is None and pair_id is None:
            return Pair(self.graph, self.graph_path, self.annotation_path, self.pair_id)
        return Pair(graph or self.graph, graph_path, annotation_path, pair_id)

    def capabilities(self, graph=None, graph_path=None, annotation_path=None, pair_id=None):
        """GET /traverse/capabilities: what one index supports. On a multi-graph server name
        it: |graph| (the client's own pair by default, selection()) and, when the name spans
        several graphs, |graph_path|; |annotation_path| where its pairs share one graph; or
        |pair_id| (the n1 id) in place of the paths."""
        graph, graph_path, annotation_path, pair_id = self.selection(
            graph, graph_path, annotation_path, pair_id)
        sel = pair_members(graph, graph_path, annotation_path, pair_id)
        sel = sel['pairs'][0] if 'pairs' in sel else sel
        if sel.get('pair_id'):
            # a pair named by its id only: as wire_pair() sends it (resolved to its paths, or
            # by its id to a server that takes it)
            sel = self.wire_pair(sel['graph'], sel['pair_id'])
        query = ['%s=%s' % (k, urllib.parse.quote(sel[k], safe=''))
                 for k in ('graph', 'graph_path', 'annotation_path', 'pair_id') if sel.get(k)]
        # a multi-graph server's probe of one pair is an envelope of one entry
        return single_answer(self._request('GET', '/traverse/capabilities'
                                           + ('?' + '&'.join(query) if query else '')))

    def supports_coordinates(self, graph=None, graph_path=None, refresh=False, *,
                             annotation_path=None, pair_id=None):
        """Whether the index (|graph|, |graph_path|, |annotation_path| or |pair_id| on a
        multi-graph server; the client's own pair by default) reports record coordinates: its GET
        /traverse/capabilities states the coordinates block (feature level 6: on every
        index, `supported` saying whether this one records them) with supported true, and
        supports_trace. Read once per pair and kept (|refresh| reads it again; a 400 that
        refuses a coordinates field forgets it). A server whose capabilities route fails
        with an HTTP error is taken as not supporting them."""
        key = self.selection(graph, graph_path, annotation_path, pair_id)
        got = None if refresh else self._coord_support.get(key)
        if got is None:
            try:
                if graph is None and graph_path is None and annotation_path is None \
                        and pair_id is None:
                    caps = self.capabilities()
                else:
                    # the keywords given only (a subclass's capabilities() may take fewer)
                    kw = {k: v for k, v in (('annotation_path', annotation_path),
                                            ('pair_id', pair_id)) if v is not None}
                    caps = self.capabilities(graph, graph_path, **kw)
            except TraverseError:
                return False
            block = caps.get('coordinates') if isinstance(caps, dict) else None
            got = isinstance(block, dict) and block.get('supported') is True \
                and caps.get('supports_trace') is True
            self._coord_support[key] = got
        return got

    def auto_coordinates(self, strategy, graph=None, graph_path=None, annotation_path=None,
                         pair_id=None):
        """auto_coordinates(self, ...): the automatic rule for record coordinates."""
        return auto_coordinates(self, strategy, graph, graph_path, annotation_path, pair_id)

    def server_capabilities(self):
        """GET /capabilities: the server-wide document — routes and features, feature_level,
        algorithm_version, mode (single | multi) and the graph names, the attempts block and
        how deadlines are checked (deadline_check). Answered while a single index loads
        (`ready: false`)."""
        return self._request('GET', '/capabilities')

    def resolve(self, sequence, *, labels=None, discover=None, select=None, time_budget_ms=None,
                **opts):
        """POST /resolve -> the answer, a dict as the server wrote it (SPEC §4).

        |time_budget_ms|: the request's deadline, sent as bounds.time_budget_ms (merged into a
        `bounds` given in |opts|; giving it there too with another value is a ValueError). It
        must exceed the server's finalisation reserve (capabilities()['resolve']['time_budget']
        ['finalize_reserve_ms'], 250 ms) and is lowered to its max_time_ms, stated in the
        answer's limits.clamped. The answer then carries `limits` and `stop`: stop None when the
        work completed, else {phase, reason: 'time', resolved_kmers, query_kmers, resolved_bp,
        query_bp, remainder_from_bp, message}, the answer being exactly the resolve of the
        sequence's first resolved_kmers k-mers -- resolve sequence[remainder_from_bp:] for the
        rest. An answer that could not be built within the budget raises ResolveDeadline (503).
        A server without the field (no `resolve` block in its capabilities) refuses it: a
        TraverseError 400 naming it. None (the default): no bounds.time_budget_ms is sent.

        A discovery is bounded by the server (SPEC §4.6; capabilities()['resolve']['discovery']):
        discover.max_labels above its max_answer_labels is lowered to it and stated in
        limits.clamped (`limits` then present also without time_budget_ms, holding the clamp
        alone); its label cap's cut is stated as `labels_capped` (absent where it did not bind;
        ranking_exact true exactly when the kept labels number the effective discover.max_labels,
        the clamped one, and the last supports more than untracked_max_kmers), and its memory
        bound's stop as `stop` with reason 'memory' (the same members as a time stop's, with or
        without time_budget_ms).

        The explicit seeds of |select| are capped by the server (SPEC §4.1;
        capabilities()['resolve']['max_explicit_seeds'], 1024 by default, 0 for none): more raise
        TraverseError 400 with the body {error, code: 'too_many_explicit_seeds',
        max_explicit_seeds}; a select.max_seeds above the cap is lowered to it and stated in
        limits.clamped (SPEC §4.3). An explicit seed naming one label twice raises TraverseError
        400 {error: "request.select.seeds[i].labels: Duplicate seed label 'x'"} (SPEC §4.3: the
        words /traverse refuses such a seed with).

        The discovery's memory bound (capabilities()['resolve']['discovery']['max_memory_mb'])
        also bounds explicit labels on the row paths (the same stop with reason 'memory') and
        what the answer holds (SPEC §4.6, the answer account): an answer priced past it raises
        TraverseError 400 with the body {error, code: 'resolve_answer_too_large', max_memory_mb}
        (nothing partial; run_format='rle', a larger min_block_kmers, fewer labels or seeds, or a
        shorter sequence may fit)."""
        req = {'sequence': sequence}
        if labels is not None:
            req['labels'] = list(labels)
        if discover is not None:
            req['discover'] = discover
        if select is not None:
            req['select'] = select
        req.update(opts)
        if time_budget_ms is not None:
            bounds = dict(req.get('bounds') or {})
            if 'time_budget_ms' in bounds and bounds['time_budget_ms'] != time_budget_ms:
                raise ValueError('time_budget_ms=%r and bounds.time_budget_ms=%r differ'
                                 % (time_budget_ms, bounds['time_budget_ms']))
            bounds['time_budget_ms'] = time_budget_ms
            req['bounds'] = bounds
        if self.graph and 'graph' not in req:
            req['graph'] = self.graph
        return self._request('POST', '/resolve', req)

    def traverse_raw(self, request):
        """POST /traverse |request| as wire_request() sends it."""
        return self._request('POST', '/traverse', self.wire_request(request))

    # ---------------------------------------------------------------- pairs named by id

    def states_pair_ids(self, refresh=False):
        """Whether the server names pairs by their n1 id and takes `pairs: [{graph,
        pair_id}]` and ?pair_id=: its GET /capabilities states multigraph.pair_id "n1"
        (docs/DESIGN-pair-identity.md section 3.1). Read once per server process, as
        server_feature_level(); a server without the route states nothing."""
        try:
            doc = self._server(refresh)
        except TraverseError as e:
            if e.status != 404:
                raise
            self._server_doc = doc = {}
        mg = doc.get('multigraph') if isinstance(doc, dict) else None
        return isinstance(mg, dict) and mg.get('pair_id') == 'n1'

    def wire_pair(self, graph, pair_id):
        """The selector members that name the pair (|graph|, |pair_id|) on the wire -- the one
        place a pair named by its id only is resolved (TraverseClient(pair_resolver=)):
        {graph, graph_path, annotation_path} by the resolver, after its answer's n1 id was
        recomputed and found equal; else {graph, pair_id} where the server states
        multigraph.pair_id "n1" (states_pair_ids()); else PairIdError pair_unresolved.
        Never a selector without the pair: the path form would have named it."""
        pid = check_pair_id(pair_id)
        if self.pair_resolver is not None:
            got = self.pair_resolver(pid)
            if got is None:
                raise PairIdError('pair_unresolved', 'pair %s of %s: the pair_resolver knows '
                                  'no pair with this id; nothing was sent' % (pid, _shown(graph)))
            try:
                g, gp, ap = got
            except (TypeError, ValueError):
                raise PairIdError('pair_unresolved', 'pair %s of %s: the pair_resolver answered '
                                  'a %s, not (graph, graph_path, annotation_path); nothing was '
                                  'sent' % (pid, _shown(graph), type(got).__name__)) from None
            if not all(isinstance(v, str) and v for v in (g, gp, ap)):
                raise PairIdError('pair_unresolved', 'pair %s of %s: the pair_resolver answered '
                                  'no graph name, graph path or annotation path to send; '
                                  'nothing was sent' % (pid, _shown(graph)))
            own = configured_pair_id(g, gp, ap)
            if own != pid or g != graph:
                raise PairIdError('pair_id_mismatch', 'pair %s of %s: the pair_resolver resolved '
                                  'it to a pair of %s whose n1 id is %s; nothing was sent'
                                  % (pid, _shown(graph), _shown(g), own))
            return {'graph': g, 'graph_path': gp, 'annotation_path': ap}
        if self.states_pair_ids():
            return {'graph': graph, 'pair_id': pid}
        raise PairIdError('pair_unresolved', 'pair %s of %s is named by its n1 id only, and the '
                          'server does not state multigraph.pair_id "n1" in GET /capabilities '
                          '(it takes path selectors only): give the client a pair_resolver (n1 '
                          'id -> graph, graph_path, annotation_path), or send it to a server '
                          'that names pairs by id; nothing was sent' % (pid, _shown(graph)))

    def wire_request(self, request):
        """|request| as it is sent: each `pairs` object that names its pair by its id only
        replaced by wire_pair()'s members, in a copy (the request itself -- as a store keeps
        it, or a NextRequest -- keeps the id); one with both paths and their own id by its
        paths. Any other request is returned as it is. PairIdError where an object states an
        id beside one path, or beside paths that are not its."""
        pairs = request.get('pairs') if isinstance(request, dict) else None
        if not isinstance(pairs, list) or not any(
                isinstance(p, dict) and p.get('pair_id') is not None for p in pairs):
            return request
        out = []
        for i, p in enumerate(pairs):
            if not isinstance(p, dict) or p.get('pair_id') is None:
                out.append(p)
                continue
            if not _text(p.get('graph')):
                raise PairIdError('bad_pair_id', 'request.pairs[%d]: a pair_id names a pair of a '
                                  'graph: give the graph too' % i)
            sel = pair_members(p['graph'], p.get('graph_path'), p.get('annotation_path'),
                               p['pair_id'])['pairs'][0]
            if sel.get('pair_id'):
                sel = self.wire_pair(sel['graph'], sel['pair_id'])
            rest = {k: v for k, v in p.items() if k not in Pair._fields}
            out.append(dict(sel, **rest))
        return dict(request, pairs=out)

    def cancel(self, attempt_id, wait_ms=None, not_after_ms=None):
        """POST /traverse/cancel -> AttemptAnswer (the server's answer, a dict, with its HTTP
        status): `cancelled: True` and `state` stopping (or finished, within |wait_ms|) when
        the attempt was asked to stop; `cancelled: False` with `state` finished or unknown
        (HTTP 404) when nothing runs under the id. An unknown id is tombstoned (`tombstone:
        True`): a request arriving later with it is refused (409, AttemptConflict) and never
        runs there while the tombstone is held. When the server holds its maximum of
        tombstones (429, reason `tombstones_full`) or keeps none (429, `no_suppression`) the
        id is NOT tombstoned: nothing is promised.

        |not_after_ms| (feature level 5): the not_after_ms of the request being cancelled.
        The tombstone is then held until no copy of that request can start, within the
        server's tombstone_hold_s (900 s by default; tombstone_max_s before train 3), and the
        answer states so (`suppression`: covers_admission). One further ahead of the server's
        clock than its attempts.not_after_horizon_ms (hard_cap_ms + the skew, 901 s by
        default; H3) is refused 400 `not_after_beyond_horizon`, as on /traverse;
        release_verdict() says whether the attempt's capacity may be released on it. Sent
        only to a server at feature level 5 or above: to an older one the cancel goes out
        without it (the cancel still stops the attempt; |sent| shows what went out, and the
        answer states no suppression, which release_verdict() never releases on). A process
        that refuses it as unknown although the level read was 5 (an older binary now
        answers at the address) gets the cancel again without it, once."""
        payload = {'attempt_id': attempt_id}
        if wait_ms is not None:
            payload['wait_ms'] = wait_ms
        if not_after_ms is not None and self._level(SUPPRESSION_LEVEL)[1]:
            payload['not_after_ms'] = not_after_ms
        try:
            status, out = self._request('POST', '/traverse/cancel', payload,
                                        answers=(404, 429), with_status=True)
        except TraverseError as e:
            if 'not_after_ms' not in payload or \
                    _unknown_field(e.status, e.message) != 'not_after_ms':
                raise
            # the 400 refused the cancel before it was read: nothing was stopped and nothing
            # tombstoned, so sending it again is safe, and stopping the attempt is the
            # cancel's first job (the request it cancels would otherwise run). _request has
            # forgotten the level that let not_after_ms out
            payload = {k: v for k, v in payload.items() if k != 'not_after_ms'}
            status, out = self._request('POST', '/traverse/cancel', payload,
                                        answers=(404, 429), with_status=True)
        return AttemptAnswer(out, status, 'cancel', payload)

    def attempt(self, attempt_id):
        """GET /traverse/attempt/{attempt_id} -> AttemptAnswer: the attempt's state (running |
        stopping | finished, the reason, when it stopped, the response written, its usage), or
        `state: unknown` (HTTP 404) once it is no longer retained -- for a tombstoned id with
        `tombstone: True` and, at feature level 5, its suppression (reading it does not extend
        it)."""
        status, out = self._request('GET', '/traverse/attempt/' + attempt_id, answers=(404,),
                                    with_status=True)
        return AttemptAnswer(out, status, 'attempt')

    def build_request(self, seeds, strategy=None, *, detail='graphlet', timing=True,
                      graph=None, attempt_id=None, budget_id=None, locus_id=None,
                      not_after_ms=None, expect_server_instance=None, coordinates=None,
                      max_coordinate_occurrences=None, graph_path=None, annotation_path=None,
                      pair_id=None):
        """The /traverse request (a dict). |graph|, |graph_path|, |annotation_path|,
        |pair_id|: the pair of a multi-graph server it runs on (the client's own when none is
        given, selection()), written by pair_members(): `graph` and `graph_path`, or a `pairs`
        object with the annotation_path or the pair_id. |expect_server_instance| (with
        |attempt_id| only, feature level 5): the server_instance the attempt is meant for -- a
        string, or 'auto' for the server's own (server_instance(), kept for that process). An
        explicit one to a server below feature level 5 raises UnsupportedFeature (nothing is
        sent: that server would register the attempt_id and refuse the field with a 400);
        'auto' is then left out. Either way AttemptSent.of(request) states what is sent.

        |coordinates| (feature level 6): True sets strategy.output.coordinates, False
        removes it (the same as false to the server, and valid to one that does not know
        the field), None leaves the strategy as given (pure: no capabilities are read here;
        traverse() applies the automatic rule). |max_coordinate_occurrences|: the cap per
        occurrence list (an int >= 1 or "unlimited"). The cap is removed whenever the
        resulting output.coordinates is not true: the server refuses it alone."""
        if expect_server_instance is not None and attempt_id is None:
            raise ValueError('expect_server_instance is sent with attempt_id only')
        if coordinates is not None and not isinstance(coordinates, bool):
            # by type, not by ==: 0 and 1 equal False and True, but the dispatch below is
            # by identity, so they would be taken for None (coordinates=0 would keep
            # coordinates, coordinates=1 would drop an explicit cap)
            raise ValueError('coordinates is True, False or None, not %r' % (coordinates,))
        if max_coordinate_occurrences is not None and not _valid_cap(max_coordinate_occurrences):
            # refused here, before anything is sent: the server's 400 would come after it
            # registered an attempt_id, refused then while the server holds it (the frozen
            # range of strategy.output.max_coordinate_occurrences)
            raise ValueError('max_coordinate_occurrences is an integer in [1, %d] or '
                             '"unlimited", not %r' % (_MAX_CAP, max_coordinate_occurrences))
        norm = []
        for s in seeds:
            norm.append({'sequence': s} if isinstance(s, str) else dict(s))
        strategy = copy.deepcopy(strategy) if strategy else {}
        output = strategy.setdefault('output', {})
        output['detail'] = detail
        output['timing'] = timing
        if coordinates is True:
            output['coordinates'] = True
        elif coordinates is False:
            output.pop('coordinates', None)
        if max_coordinate_occurrences is not None:
            output['max_coordinate_occurrences'] = max_coordinate_occurrences
        strip_coordinate_cap(strategy)
        req = {'seeds': norm, 'strategy': strategy}
        if self.release:
            req['release'] = self.release
        req.update(pair_members(*self.selection(graph, graph_path, annotation_path, pair_id)))
        for key, value in (('attempt_id', attempt_id), ('budget_id', budget_id),
                           ('locus_id', locus_id), ('not_after_ms', not_after_ms)):
            if value is not None:
                req[key] = value
        if expect_server_instance is not None:
            level, ok = self._level(SUPPRESSION_LEVEL)
            if expect_server_instance == 'auto':
                instance = self.server_instance() if ok else None
                # reading the instance read the level the process states, which wins over
                # a level given to the client
                if instance is not None and self._level(SUPPRESSION_LEVEL)[1]:
                    req['expect_server_instance'] = instance
            elif not ok:
                raise UnsupportedFeature('expect_server_instance', level, SUPPRESSION_LEVEL)
            else:
                req['expect_server_instance'] = expect_server_instance
        return req

    def traverse(self, seeds, strategy=None, *, detail='graphlet', timing=True, graph=None,
                 attempt_id=None, budget_id=None, locus_id=None, not_after_ms=None,
                 expect_server_instance=None, coordinates='auto',
                 max_coordinate_occurrences=None, graph_path=None, annotation_path=None,
                 pair_id=None):
        """-> TraverseResponse(envelope, graphlets, errors). A seed whose permitted set
        could not be derived is an error result ({seed, error}), not an exception: the
        other seeds' graphlets are still returned (so is a seed an attempt never started,
        cancelled or past its bound: resource_stop phase not_started). A request error is
        raised (TraverseError; with attempt_id its usage is in the error's .usage), an
        attempt whose response was stopped at its bound raises AttemptAtBound, a request
        whose |not_after_ms| (Unix epoch ms) had passed on the server's clock raises
        AttemptExpired (not started), an attempt_id the server holds (running, retained, held
        or tombstoned) AttemptConflict, and an |expect_server_instance| (build_request())
        that is not the server's InstanceMismatch. The response's .sent, and the .sent of
        any exception raised once the request was built (a TraverseError, and a transport
        error too: a timeout, a reset connection), state the attempt fields as sent (for
        release_verdict()); an error raised while building it (UnsupportedFeature) sent
        nothing. |graph|, |graph_path|, |annotation_path|, |pair_id|: the pair, as in
        build_request();
        each graphlet keeps the pair a multi-graph server answered on (its envelope), and
        its continuations run there (next_request(), deepen()).

        |coordinates| (feature level 6): 'auto' (the default) asks for record coordinates
        for a support: trace strategy that does not set output.coordinates, where the index
        reports them (supports_coordinates()) -- under a request memory budget only while
        the depth gate passes (AUTO_COORDINATES_UNDER_MEMORY_BUDGET); what it decided is
        .coordinates_auto, and .notes states a decision not to ask under a memory budget,
        or, where coordinates it asked for under one share a memory stop, how to drop them
        (drop_coordinates_note()). True / False / None as in build_request(). A server that
        refuses the field (an older one: 400) raises; nothing is retried without it."""
        auto = None
        if coordinates == 'auto':
            coordinates, auto = auto_coordinates(self, strategy, graph, graph_path,
                                                 annotation_path, pair_id)
        req = self.build_request(seeds, strategy, detail=detail, timing=timing, graph=graph,
                                 graph_path=graph_path, annotation_path=annotation_path,
                                 pair_id=pair_id,
                                 attempt_id=attempt_id, budget_id=budget_id, locus_id=locus_id,
                                 not_after_ms=not_after_ms,
                                 expect_server_instance=expect_server_instance,
                                 coordinates=coordinates,
                                 max_coordinate_occurrences=max_coordinate_occurrences)
        # a pair named by its id is resolved (or refused) before the attempt is stated as
        # sent: a PairIdError sent nothing
        req = self.wire_request(req)
        sent = AttemptSent.of(req)
        try:
            resp = self.response(self.traverse_raw(req))
        except Exception as e:
            # an unanswered attempt (the socket timed out, the connection was reset) is the
            # one a ledger must cancel and judge with release_verdict(), and with 'auto' only
            # the client knows whether an instance went out
            try:
                e.sent = sent
            except (AttributeError, TypeError):     # an exception type without attributes
                pass
            raise
        resp.sent = sent
        resp.seeds_sent = copy.deepcopy(req['seeds'])
        resp.coordinates_auto = auto
        if auto is not None and auto.get('note'):
            resp.notes.append(auto['note'])
        dropped = drop_coordinates_note((resp.raw or {}).get('results'), auto)
        if dropped is not None:
            resp.notes.append(dropped)
        return resp

    @staticmethod
    def response(out):
        from .parser import from_response
        # a multi-graph server's answer: its one entry, the attempt's usage beside it
        out = single_answer(out)
        # the answer's shape, checked before it is read (review m4): an object whose results
        # are objects; anything else is this answer's ResponseFormatError
        if not isinstance(out, dict):
            raise ResponseFormatError('the answer is a JSON %s, not an object'
                                      % _json_kind(out), out)
        if not isinstance(out.get('results', []), list) or not all(
                isinstance(r, dict) for r in out.get('results', [])):
            raise ResponseFormatError('the answer\'s results are not a list of objects', out)
        for i, r in enumerate(out.get('results', [])):
            # a graphlet is MGT text (Opus fix review N3)
            if 'graphlet' in r and not isinstance(r['graphlet'], str):
                raise ResponseFormatError('results[%d].graphlet is a JSON %s, not a string'
                                          % (i, _json_kind(r['graphlet'])), out)
        envelope = {k: v for k, v in out.items() if k != 'results'}
        graphlets, errors = [], []
        for i, r in enumerate(out.get('results', [])):
            if 'graphlet' in r:
                graphlets.append(from_response(r, out))
            else:
                graphlets.append(None)
                errors.append({'index': i, 'error': r.get('error'), 'result': r})
        return TraverseResponse(envelope, graphlets, errors, out)

    def deepen(self, graphlet, arm, leaves, bp=None, **overrides):
        """Iterative deepening stays a backend call: the graphlet derives the
        continuations and builds the request (Graphlet.next_request), the server runs it.
        The response's notes carry the request's: how the continuation may differ from
        one uninterrupted walk (a loss budget conservative for the labels that ended at a
        lower loss, ...). Walks whose continuations carry different labels cannot share a
        request (IncompatibleContinuations): deepen them one walk at a time.

        The continuation runs on the pair the graphlet was answered on (a multi-graph
        server's envelope names it, and next_request() keeps it), else on the client's own
        pair (selection()); a selection in |overrides| (graph, graph_path, annotation_path,
        pair_id, graphs, pairs) replaces either whole. A graphlet answered on a pair the
        envelope names by its id only (a document stored after the engine named pairs by id)
        is continued by that id: `pairs: [{graph, pair_id}]`."""
        if answered_pair(getattr(graphlet, 'envelope', None)) is None \
                and not any(k in overrides for k in SELECTORS):
            overrides.update((k, v) for k, v in zip(Pair._fields, self.selection()) if v)
        req = graphlet.next_request(arm, leaves, bp, **overrides)
        out = self.response(self.traverse_raw(req))
        out.notes = list(getattr(req, 'notes', ()))
        out.seeds_sent = copy.deepcopy(req.get('seeds')) if isinstance(req, dict) else None
        return out

    def traverse_pairs(self, seed, strategy=None, *, pairs=None, detail='graphlet',
                       timing=True, coordinates='auto', not_after_ms=None):
        """One /traverse request per pair for one seed (docs/DESIGN-graphlet-merge.md §5.2):
        the fan-out whose results union_pairs() combines. |pairs|: graph names, {graph,
        graph_path, annotation_path} objects or {graph, pair_id} objects; by default (None)
        every pair GET /capabilities lists (graph_summary.pairs, by their paths or by their
        n1 ids), or the client's own pair on a single-graph server; an empty list selects no
        pair and sends nothing. A pair named by its id that GET /capabilities lists with its
        paths (an engine before ids) is sent by those paths, the engine's own spelling (their
        id is the one given); one listed by its id only is sent by its id. Every answer is
        kept as returned (a TraverseResponse of one seed, which keeps the seed it sent); a
        refused request (TraverseError, ResponseFormatError for an answer not shaped as a
        response), a transport failure (OSError) or a body that does not parse
        (GraphletFormatError) becomes {pair, error}, and any other exception is raised.
        attempt_id and expect_server_instance are not accepted: attempt ids share one
        server-wide namespace with request ids, so one id sent to several pairs is refused
        from the second request on (a ledger-managed fan-out is a later extension).

        -> PairResults(results, seed, expected_pairs, seqs_order): |seed| the oriented
        seed sent, |expected_pairs| the pair objects the requests were sent to, |seqs_order|
        each pair's record-mapping certificate as GET /capabilities states it (keys:
        (graph, pair_id), or the index namespace on a single-graph server; values
        'verified', 'count_only' or None). A pair is stated by its graph and n1 id, no path
        (Pair.named()): an engine's listed paths, or a caller's, become their id; only a
        caller's pair without the annotation path, which has no id, keeps what it gave
        ((graph, graph_path, annotation_path) keys)."""
        seq = seed.get('sequence') if isinstance(seed, dict) else seed
        if not isinstance(seq, str) or not seq:
            raise ValueError('traverse_pairs() takes one seed: a sequence or {sequence, ...}, '
                             'not %r' % (seed,))
        if pairs is not None:
            if isinstance(pairs, (str, dict)) or not isinstance(pairs, (list, tuple)):
                raise ValueError('pairs= is a list of graph names, {graph, graph_path, '
                                 'annotation_path} or {graph, pair_id} objects, not a %s'
                                 % type(pairs).__name__)
            if not pairs:
                # an explicit empty selection: no pair, no request (only pairs=None
                # discovers the pairs)
                return PairResults([], seq, (), {})
        caps = self.server_capabilities()
        caps = caps if isinstance(caps, dict) else {}
        summary = caps.get('graph_summary') if isinstance(caps.get('graph_summary'), dict) \
            else {}
        # each listed pair with its Pair: its id stated (an engine after ids) or computed from
        # its paths (one before)
        listed = [(p, _pair(p['graph'], p.get('graph_path'), p.get('annotation_path'),
                            p.get('pair_id'), 'GET /capabilities graph_summary.pairs[%d]' % i))
                  for i, p in enumerate(summary.get('pairs') or [])
                  if isinstance(p, dict) and _text(p.get('graph'))]
        targets = []
        if pairs is None:
            seen = set()
            for _, t in listed:
                if t.identity not in seen:
                    seen.add(t.identity)
                    targets.append(t)
        else:
            for i, p in enumerate(pairs):
                if isinstance(p, str) and p:
                    targets.append(Pair(p))
                elif isinstance(p, dict) and _text(p.get('graph')):
                    targets.append(_pair(p['graph'], p.get('graph_path'),
                                         p.get('annotation_path'), p.get('pair_id'),
                                         'pairs[%d]' % i))
                else:
                    raise ValueError('pairs[%d] is a graph name, a {graph, graph_path, '
                                     'annotation_path} or a {graph, pair_id} object, not a %s'
                                     % (i, type(p).__name__))
        seqs_order = {}
        if targets:
            sent = []
            for t in targets:
                hit = [(p, lp) for p, lp in listed if lp.graph == t.graph and _names(t, lp)]
                vals = {p.get('seqs_order') for p, _ in hit if 'seqs_order' in p}
                if len(vals) == 1:
                    seqs_order[_seqs_key(t)] = vals.pop()
                if t.pair_id and not t.graph_path and not t.annotation_path and len(hit) == 1 \
                        and hit[0][1].graph_path and hit[0][1].annotation_path:
                    # named by its id, listed with its paths: sent by the listed paths (the
                    # engine's own spelling; their id is t's, _names())
                    sent.append(hit[0][1])
                else:
                    sent.append(t)
            expected = tuple(_pair_object(t) for t in targets)
        else:
            # a single-graph server (or no pair listed): the client's own pair, named by its
            # index namespace
            sel = Pair(*self.selection())
            targets = sent = [sel if sel.graph else Pair(None)]
            ns = None
            try:
                ns = _text((self.capabilities() or {}).get('index_ns'))
            except TraverseError:
                pass
            ps = caps.get('pair_summary')
            if ns is not None and isinstance(ps, dict) and 'seqs_order' in ps:
                seqs_order[ns] = ps.get('seqs_order')
            if sel.graph:
                expected = (_pair_object(sel),)
            else:
                expected = ({'index_ns': ns},) if ns is not None else ()
        # traverse() parses each pair's graphlets (the parser loads _codec): imported here, so
        # that importing the client stays light
        from ._codec import GraphletFormatError
        results = []
        for t, w in zip(targets, sent):
            # the keywords given only (a subclass's traverse() may take fewer)
            kw = {} if w.pair_id is None or (w.graph_path and w.annotation_path) \
                else {'pair_id': w.pair_id}
            try:
                results.append(self.traverse([seed], strategy, detail=detail, timing=timing,
                                             graph=w.graph, graph_path=w.graph_path,
                                             annotation_path=w.annotation_path,
                                             coordinates=coordinates,
                                             not_after_ms=not_after_ms, **kw))
            except (TraverseError, OSError, GraphletFormatError) as e:
                # a refused request (TraverseError), a transport failure (an OSError: a
                # connection refused, a timeout, a cut response) or a body that does not
                # parse (GraphletFormatError: this pair's answer only); anything else -- a
                # bad argument, a library bug -- is raised, not reported as every pair
                # failing
                pair = _pair_object(t) or (expected[0] if expected else {})
                results.append({'pair': pair, 'error': e})
        return PairResults(results, seq, expected, seqs_order)


def _names(t, lp):
    """Whether selector Pair |t| names listed Pair |lp| (of the same graph): by id where both
    have one (an id computed from both paths included), else by the paths |t| states."""
    if t.pair_id and lp.pair_id:
        return t.pair_id == lp.pair_id
    if t.pair_id:
        return False
    return (t.graph_path is None or lp.graph_path == t.graph_path) \
        and (t.annotation_path is None or lp.annotation_path == t.annotation_path)


def _seqs_key(t):
    """traverse_pairs()' seqs_order key of target |t|: (graph, pair_id), no path."""
    if t.pair_id:
        return t.graph, t.pair_id
    return t.graph, t.graph_path, t.annotation_path


def _pair_object(t):
    """The pair object traverse_pairs() states for target |t| (Pair.named())."""
    return t.named()


@dataclass(frozen=True)
class PairResults:
    """traverse_pairs()' answer: |results| one per pair (a TraverseResponse, or {pair,
    error}), |seed| the oriented seed sent, |expected_pairs| the pair objects the requests
    were sent to, |seqs_order| each pair's record-mapping certificate (GET /capabilities).
    union_pairs(r.results, seed=r.seed, expected_pairs=r.expected_pairs,
    seqs_order=r.seqs_order) combines them."""
    results: list
    seed: str
    expected_pairs: tuple
    seqs_order: dict
