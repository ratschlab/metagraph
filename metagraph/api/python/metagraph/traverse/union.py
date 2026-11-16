"""The cross-pair union of one seed's graphlets (docs/DESIGN-graphlet-merge.md §2, phase 1a-U).

union_pairs(results) takes everything one seed's request returned across several pairs (the
graphs of a release such as r237's 11, or two databases): the graphlets, and the pairs that
answered without one. It returns a PairUnion: per arm, one trie of every part's spellings,
with the pairs and labels that carry each segment. The parts stay unchanged and every answer
is computed from them (§4.4). It is not exact and does not claim to be; five statements come
with every summary, printout and export (§2.1):

  1. coverage per pair (§2.6);
  2. provenance per segment and walk (§2.4);
  3. what a turn-off means: a pair not following a walk did not follow it under its own
     request (only a span marked absence 'graph' says the pair's graph lacks it);
  4. labels per pair: shown as <tag>/<name>, equal names in two pairs not claimed to be one
     record, header names on pairs whose record mapping was not checked marked so (§2.5);
  5. the scope of the paths: source-supported only in this phase (§2.7) -- "A walk that no
     single graph holds entirely is in no part, so not in this union."

Inferred paths (paths='inferred') are a later phase and are refused 'not_offered' here.
Positions are outward from the seed's boundary in the seed's orientation; bases are compared
in walking order (outward), so spellings compare across k, regime and radius.

A StitchedGraphlet (_stitch.py, experimental) is one pair's answer (§2.3): its root goes
into the trie as a plain part, and each child's first-parent tree at the trie node where the
walk it continues ends (the continued walk ends no union walk there). Where no junction is
read, the part answers as its root alone. A span read past a junction names the junctions it
reads (PairSpan.junctions), is never complete, is certified at most to the first junction
(the stitched walk's certified_to_bp), states absence 'graph' (a dead end) or 'uncertified'
only, and carries the child's own labels from the junction on: labels are not joined across a
junction (DESIGN-graphlet-merge.md D.1, phase 1b), so the union states less than
StitchedWalk.labels_full.
"""

import collections
import hashlib
from dataclasses import dataclass
from typing import Any, Optional

from . import budget as _B
from . import composite as C
from . import derive
from . import ops as _ops
from ._codec import REASON, GraphletFormatError
from .attempts import TraverseResponse
from .budget import (BASE_SHIFT, DICT, DICT_KEY, INT, LIST_ITEM, SET, SET_ITEM, STR, TUPLE,
                     W_ELEM, W_ROW, W_STEP, LocalBudgetExceeded, Partial)
from ._stitch import _RANK as _CAUSE_RANK
from ._stitch import StitchedGraphlet, part_name
from .client import answered_pair
from .model import ARM_SIDES, BadSelector, Graphlet, UnknownLabel

__all__ = [
    'union_pairs', 'PairUnion', 'PairPart', 'PairAbsent', 'PairLabel', 'SpanLabel',
    'LabelInterval', 'SegmentSupport', 'UnionSegment', 'Route', 'PairSpan', 'SpanJunction',
    'Carriers', 'UnionWalk', 'Coverage', 'Statements', 'TEXTS', 'CATEGORIES',
]

OP = 'union_pairs()'

# the refusal texts (Appendix A.2), filled with % by the check that refuses
TEXTS = {
    'bad_expected': 'expected_pairs= takes a list of the pairs the request was sent to, as '
                    'graph names or {graph, graph_path, annotation_path} objects; %s is not '
                    'one',
    'bad_expected_str': 'expected_pairs= is a string: pass a list of graph names',
    'bad_seqs_order': 'seqs_order= maps pairs (graph names or pair objects) to "verified", '
                      '"count_only" or None, as GET /capabilities states them; %s is not one',
    'not_combinable': 'union_pairs() takes results: a Graphlet, a StitchedGraphlet, a '
                      'TraverseResponse, {response, sequence}, or {pair, graphlet}, '
                      '{pair, result} or {pair, error}; not %s',
    'view_unsupported': 'union_pairs() reads whole retrievals, not views: pass view.backing, '
                        'and select labels on the union (labels=)',
    'not_one_seed': 'result %s holds %d seed results: pass {response, index} for this '
                    'union\'s seed, or it as {pair, graphlet} or {pair, result}',
    'not_one_seed_zero': 'result %s holds no seed result',
    'not_one_seed_range': 'result %s has no seed result %d (it holds %d)',
    'result_has_graphlet': 'result %s carries a graphlet: pass it as a Graphlet '
                           '(Graphlet.from_response(result, response)) or as {pair, graphlet}',
    'no_graphlet': 'result %s was walked (outcome.walks %s) but carries no graphlet: fetch '
                   'with detail=\'graphlet\'',
    'no_seed': 'no result carries a graphlet, so the union has no seed sequence: pass seed= '
               '(the oriented seed the request sent)',
    'pair_mismatch': 'result %s was answered on pair %s, and the caller states pair %s: state '
                     'the pair it was answered on, or none',
    'no_pair': 'result %s names no pair (no graph tag, no index_ns): pass it as {pair, '
               'graphlet} or {pair, result}',
    'ambiguous_pair': 'pair %s matches both %s and %s: state its paths, release or index_ns',
    'no_bases': 'part %s carries no bases (output.sequences: false): a union by sequence '
                'needs every part\'s bases',
    'different_alphabet': 'part %s walks the alphabet %s and part %s the alphabet %s: their '
                          'bases are not comparable',
    'different_seed': '%s is not %s (they first differ at base %d): a union by sequence needs '
                      'one oriented seed, spelled base for base the same',
    'different_seed_length': '%s is %d bp and %s %d bp (equal for the first %d bases): a '
                             'union by sequence needs one oriented seed',
    # an absent record states the seed's length, not its bases
    'different_seed_length_absent': '%s is %d bp and %s %d bp: a union by sequence needs one '
                                    'oriented seed',
    'different_seed_revcomp': '%s is the reverse complement of %s: a union by sequence needs '
                              'one orientation; fetch %s with the seed as %s spells it',
    'duplicate_pair': 'results %s and %s were answered on one pair (%s): a union holds one '
                      'answer per pair; compare() reads two retrievals of one pair, and '
                      'union_retrievals() will combine them (not yet offered)',
    'duplicate_pair_request': 'results %s and %s are byte-identical bodies on pair %s from '
                              'different requests (%s differs): equal bodies are not equal '
                              'retrievals, and a union holds one answer per pair',
    'duplicate_pair_unknown': 'results %s and %s cannot be told apart (no index_fp on %s, and '
                              'their names, paths, release, index_ns and index_meta_fp agree): '
                              'they may be one pair; state each pair with {pair, graphlet} or '
                              '{pair, result}',
    'not_expected': 'pair %s (a part or an absent record) is not among expected_pairs=',
    'ambiguous_name': 'label name %r is in pairs %s: equal names in two pairs are not claimed '
                      'to be one record; select by {pair, name} or {pair, ref}',
    'label_id': 'label ids are local to each part of a pair union: select by {pair, ref} or '
                '{pair, name}',
    # %d: the least answer, the statements with the notes of what is left out
    'summary_too_small': 'the statements alone take %d bytes, more than max_bytes=%d: they '
                         'are never cut; pass max_bytes >= %d',
    'not_offered_dump': 'a pair union is not an MGT document: each part dumps on its own '
                        '(parts[i].graphlet.dump())',
    'not_offered_json': 'to_json() is detail "full" of one server result, and a pair union '
                        'has none: export its parts',
    'not_offered': '%s() is not offered on a pair union',
    'not_offered_inferred': 'paths=\'inferred\' (and max_transitions, max_paths) is not '
                            'offered on a pair union yet: inferred paths come in a later '
                            'release; paths=\'supported\' lists the walks one graph holds whole',
}

# the coverage categories, in the order of the statement (§2.6)
CATEGORIES = ('complete', 'complete_qualified', 'incomplete', 'seed_not_in_graph',
              'labels_not_in_graph', 'not_walked', 'arm_not_walked', 'failed', 'not_included')
_CATEGORY_WORDS = {
    'complete': 'complete in', 'complete_qualified': 'complete, qualified, in',
    'incomplete': 'incomplete in', 'seed_not_in_graph': 'seed not in graph:',
    'labels_not_in_graph': 'labels not in graph:', 'not_walked': 'not walked:',
    'arm_not_walked': 'arm not walked:', 'failed': 'failed:', 'not_included': 'not included:',
}
_ABSENT_OF_WALKS = {'not_in_graph': 'seed_not_in_graph',
                    'labels_not_in_graph': 'labels_not_in_graph', 'failed': 'not_walked'}
_PAIR_FIELDS = ('graph', 'graph_path', 'annotation_path', 'release', 'index_ns', 'index_fp',
                'index_meta_fp')
# the pair key's order (Appendix A.4)
_KEY_FIELDS = ('release', 'index_ns', 'graph', 'graph_path', 'annotation_path', 'index_fp',
               'index_meta_fp')
_SEQS_ORDER = ('verified', 'count_only', None)

# The memory prices of the union's own objects (bytes): fitted to tracemalloc's traced peak
# of each operation, cold, on CPython 3.11 and 3.14 (the split fixtures' 11 parts, r237's
# export, a comb), so that the budget's account bounds the call's peak (budget.py).
_M_NODE = 300           # a raw trie node: its children dict, its depth
_M_EDGE = 420           # a raw edge: its five list slots, its reference, its cover dict
_M_COVER = 110          # a cover entry of an edge (part -> (segment, offset))
_M_SEG = 200            # a part segment's entries: its nodes, the walks and stubs below it
_M_USEG = 1280          # a union segment: the object, its dicts and lists, node_in/out
_M_PIECE = 200          # a piece of a union segment (part, segment, from, to)
_M_PART = 4096          # a part's record: its identity, retrieval identity, digests
_M_RECORD = 1024        # an absent or failed record
_M_SPAN = 900           # a PairSpan with its Route, merges and continues tuples
_M_LABEL = 320          # a span label: its SpanLabel, PairLabel (with its ref) and slot
_M_INTERVAL = 360       # a LabelInterval with its PairLabel, its ints and slot
_M_NAME = 300           # shared_names() per label: the name sets, the dict and its lists
_M_PAIR_LABEL = 260     # a PairLabel: the record and its ref string
_M_COV_PAIR = 700       # a coverage entry per pair (its dict, note and category lists)
_M_COV_NOTE = 120       # a coverage note per pair beside its tag (certified to ..., k-mers)
_M_SENTENCE = 1500      # a sentence's fixed words (coverage, the turn-offs and the scope)
_M_BLOCK = 2500         # a part's evidence block and identity entry
_M_META = 6             # copies of the inputs' metadata JSON at once: an absent record's
#                         detail JSON, the digest's JSON (that detail escaped again, at most
#                         twice its size), its UTF-8 bytes and the encoder's chunks (fitted:
#                         a detail of 100,000 quotes peaks at 1.08x five copies on 3.14)
_M_LABEL_WORK = 260     # what _labels_for() holds per label while it runs
_M_FIG = 160            # a walk's figures: the record, whole_in and held_in, its slot
_DERIVATION_CAUSES = {'no_carrier': 'no label carries the whole seed',
                      'no_trace_carrier': 'no label carries the whole seed by trace',
                      'too_wide': 'the seed\'s rows are too wide to derive from',
                      'over_seed_label_cap': 'more carriers than max_seed_labels',
                      'ambiguous_header': 'a header names two records',
                      'time_budget': 'the time budget ran out'}

TURN_OFFS = 'A pair that does not follow a walk did not follow it under its own request.'
SCOPE = 'A walk that no single graph holds entirely is in no part, so not in this union.'
OVER_GIVEN = 'Over the pairs given; expected_pairs= was not passed.'


def _refuse(code, key=None, *args):
    text = TEXTS[key or code]
    return C.CombineRefused(code, text % args if args else text)


# ------------------------------------------------------------------ data classes

@dataclass(frozen=True)
class PairLabel:
    pair: str                   # the pair's tag
    ref: str                    # index-local: never joined across pairs
    name: str
    kind: str                   # 'column' | 'header'

    def __str__(self):
        return '%s/%s' % (self.pair, self.name)


@dataclass(frozen=True)
class SpanLabel:
    label: PairLabel
    from_bp: int                # where its support of this span's bases begins; 0: whole


@dataclass(frozen=True)
class LabelInterval:
    label: PairLabel
    from_bp: int                # where the label is alive within a segment (outward bp)
    to_bp: int


@dataclass(slots=True)
class SegmentSupport:
    pair: str
    n_walks: int                # the pair's displayed walks passing, by count
    n_stubs: int                # the pair's stubs passing
    pieces: tuple               # (part segment id, from_bp, to_bp), in order
    labels: tuple               # LabelInterval, from the alive sets and runs, in order
    labels_cut: bool            # the lists read were cut (a lower bound)
    part: tuple = ()            # a stitched part's sub-part key (Part.key); () for a plain part


@dataclass(slots=True)
class UnionSegment:
    id: int
    arm: str
    from_bp: int
    to_bp: int
    parent: Optional[int]
    children: tuple
    support: tuple              # SegmentSupport per pair spelling it, canonical order


@dataclass(slots=True)
class Route:
    to_bp: int
    how: str                    # 'held' | 'ends' | 'diverges'
    goes_on: bool
    end_reasons: dict           # {'path_reason': name | None, 'labels': {reason: n}}
    absence: Optional[str]
    merges: tuple               # (at_bp, trie node) of each merge crossed, the first included
    labels: tuple               # SpanLabel (carriers(), labels='full'); () otherwise
    complete: bool
    certified_to_bp: int
    labels_cut: bool = False    # the route's labels read a cut list (a lower bound)


@dataclass(slots=True)
class PairSpan:
    pair: str
    mode: str
    how: str                    # whole | extends | ends | merges | diverges | none
    to_bp: int
    walk: Any                   # the part's path id (whole, extends, ends); past a junction
                                # the stitched walk key (a tuple, StitchedWalk.key)
    stub_end: Optional[int]     # merges: the part segment ending the stub
    into: Optional[str]         # merges: the trie node of the merge, 's31' | 's30$'
    route: Optional[Route]      # merges only
    continues: tuple            # merges: PairLabels whose kept entry enters through it
    end_reasons: dict
    labels: Optional[tuple]     # SpanLabel; None where the query did not compute them
    labels_cut: bool
    complete: bool
    certified_to_bp: int
    absence: Optional[str]      # graph | walk_rule | under_strategy | uncertified | None
    k: int = 0
    regime: str = ''
    junctions: tuple = ()       # SpanJunction per junction of a stitched part the span reads


@dataclass(frozen=True)
class SpanJunction:
    """A junction of a stitched part that a span reads (§2.3): where the child continues a
    walk (|at_bp|, outward from the seed), the junction's |verdict| ('qualified': every verdict
    in this phase), the child's part key (Part.key), and the cause codes: the junction's own,
    then the walk-scope causes (W8) of the child's walk the span reads (none where the span
    names no walk there: a merges or diverges span in the child)."""
    at_bp: int
    verdict: str
    part: tuple
    causes: tuple


class Carriers(list):
    """carriers(): a list of PairSpan (canonical pair order), one line per pair in str()."""

    mapping_note = None

    def __str__(self):
        lines = []
        for s in self:
            lines.append(_carrier_line(s))
        if self.mapping_note:
            lines.append(self.mapping_note)
        return '\n'.join(lines)


@dataclass(slots=True)
class UnionWalk:
    id: int
    arm: str
    kind: str                   # 'walk' | 'stub'
    length_bp: int
    sequence: Optional[str]
    segments: tuple             # union segment ids, root to end
    spans: tuple                # PairSpan per part, canonical order
    whole_in: tuple             # the pairs with a 'whole' span
    held_in: tuple              # whole, extends, and routes held to the end


_WalkFigures = collections.namedtuple('_WalkFigures', 'kind length_bp whole_in held_in')


@dataclass(frozen=True)
class Coverage:
    arm: Optional[str]
    at_bp: int
    per_pair: dict              # tag -> {category, certified_to_bp, detail, ...}
    by_category: dict           # category -> [tags], non-empty categories in table order
    every_pair_certified_to_bp: Optional[int]
    sentence: str


@dataclass(frozen=True)
class Statements:
    arm: Optional[str]
    at_bp: int
    paths: str
    coverage: Coverage
    provenance: str
    turn_offs: str
    labels: str
    scope: str

    def sentences(self):
        """The five, in order: coverage, provenance, turn-offs, labels, scope."""
        return [_core_sentence(self.coverage), self.provenance, self.turn_offs, self.labels,
                self.scope]


@dataclass(frozen=True)
class PairPart:
    tag: str
    key: tuple
    pair_stated_by: str         # 'envelope' | 'caller'
    seqs_order: str             # 'verified' | 'count_only' | 'none' | 'unknown'
    graphlet: Any               # a Graphlet, or a StitchedGraphlet (experimental)
    identity: dict              # {verified, index_fp} | {verified: False, note}
    retrieval: tuple            # (pair key, retrieval identity as canonical JSON)
    digest: str                 # the body's digest (a stitched part's: StitchedGraphlet.digest)


@dataclass(frozen=True)
class PairAbsent:
    tag: str
    key: tuple
    pair_stated_by: str
    category: str               # seed_not_in_graph | labels_not_in_graph | not_walked | failed
    detail: dict
    seqs_order: str
    seed_check: Optional[str]   # 'sent' | 'stated' | 'unchecked'; None for a failed record


# ------------------------------------------------------------------ inputs

class _Elem:
    __slots__ = ('pos', 'form', 'graphlet', 'result', 'error', 'response', 'index',
                 'sequence', 'sent', 'stated', 'desc', 'env', 'cls', 'kind', 'by')

    def __init__(self, pos, form):
        self.pos = pos
        self.form = form
        self.graphlet = self.result = self.error = self.response = self.index = None
        self.sequence = self.sent = self.stated = self.desc = self.env = None
        self.cls = None
        self.kind = None           # 'part' | 'absent' | 'failed'
        self.by = None             # 'envelope' | 'caller'


def _root(g):
    """The graphlet a part's pair, envelope, seed, request and certificate are read from: a
    stitched part's root (its children run on the same pair, §3.1), else the part itself."""
    return g.root if isinstance(g, StitchedGraphlet) else g


def _sub_graphlets(g):
    """The retrievals a part holds: a stitched part's parts, root first, else the part."""
    return [p.graphlet for p in g.parts] if isinstance(g, StitchedGraphlet) else [g]


class _Src:
    """One retrieval the trie reads: part |pi|'s graphlet, or one part of a stitched part
    (|key| its Part.key, |offset| its Part.offset_bp, |parent| the source it continues,
    |leaf| the parent's segment ending the walk it continues, |junction| that Junction)."""
    __slots__ = ('pi', 'tag', 'graphlet', 'offset', 'key', 'parent', 'leaf', 'junction',
                 'wcauses', 'below')

    def __init__(self, pi, tag, graphlet, offset, key, parent, leaf, junction):
        self.pi, self.tag, self.graphlet, self.offset = pi, tag, graphlet, offset
        self.key, self.parent, self.leaf, self.junction = key, parent, leaf, junction
        # a child: its walks' walk-scope causes (path id -> [(code, sentence, local bp)],
        # StitchedGraphlet's W8) and the codes on the stitched walks through its junction
        self.wcauses, self.below = {}, ()


def _describe(x):
    if isinstance(x, dict):
        return 'a dict with keys %s' % ', '.join(sorted(str(k) for k in x)) if x else 'an empty dict'
    return 'a %s' % type(x).__name__


def _pair_desc(p):
    """A caller's pair (a graph name, or a dict of _PAIR_FIELDS naming a graph or, for a
    single-graph server, an index_ns) as a description, or None when it is malformed."""
    if isinstance(p, str):
        return {'graph': p} if p else None
    if isinstance(p, dict) and p and set(p) <= set(_PAIR_FIELDS) \
            and (p.get('graph') or p.get('index_ns')) \
            and all(v is None or isinstance(v, str) for v in p.values()):
        return {k: v for k, v in p.items() if v not in (None, '')}
    return None


def _env_desc(env, g=None):
    """The pair an envelope (and a body's H line) states."""
    env = env if isinstance(env, dict) else {}
    caps = env.get('capabilities') if isinstance(env.get('capabilities'), dict) else {}
    d = {}
    ap = answered_pair(env)
    if ap is not None:
        d['graph'], d['graph_path'], d['annotation_path'] = ap
    d['release'] = env.get('release') or caps.get('release') or None
    if g is not None:
        d['index_ns'], d['index_fp'], d['index_meta_fp'] = g.index_ns, g.index_fp, g.index_meta_fp
    else:
        d['index_ns'] = env.get('index_ns') or caps.get('index_ns')
        d['index_fp'] = env.get('index_fp') or caps.get('index_fp')
        d['index_meta_fp'] = env.get('index_meta_fp') or caps.get('index_meta_fp')
    return {k: v for k, v in d.items() if v not in (None, '')}


def _matches(a, b):
    """The one rule of Appendix A.4: do descriptions |a| and |b| name one pair? It decides
    every match: duplicate_pair, pair_mismatch, not_expected, not_included, an absent record
    against a part, and seqs_order=. A single-graph server's pair states no graph name: its
    tag is its index_ns, and a bare name equal to it names that pair."""
    if a.get('index_meta_fp') and b.get('index_meta_fp') \
            and a['index_meta_fp'] != b['index_meta_fp']:
        return False
    if a.get('index_fp') and b.get('index_fp'):
        return a['index_fp'] == b['index_fp']
    if a.get('graph') != b.get('graph'):
        return _names_ns(a, b) or _names_ns(b, a)
    for k in ('graph_path', 'annotation_path', 'release', 'index_ns', 'index_meta_fp'):
        if a.get(k) and b.get(k) and a[k] != b[k]:
            return False
    return True


def _names_ns(bare, d):
    """Bare name |bare| names the single-graph pair |d| (no graph name) by its index_ns."""
    return set(bare) == {'graph'} and not d.get('graph') and d.get('index_ns') == bare['graph']


def _key(d):
    return tuple(d.get(k) for k in _KEY_FIELDS)


def _sort_key(key):
    return tuple((0, b'') if v is None else (1, str(v).encode('utf-8')) for v in key)


def _show(d):
    return C.canonical_json({k: v for k, v in d.items() if v is not None}) \
        if set(d) != {'graph'} else d['graph']


def _seed_meta(result):
    s = result.get('seed') if isinstance(result, dict) else None
    return s if isinstance(s, dict) else {}


def _classify(results):
    elems = []
    for i, x in enumerate(results):
        if C.is_view(x) or isinstance(x, (Graphlet, StitchedGraphlet)):
            e = _Elem(i, 'graphlet')
            e.graphlet = x
        elif isinstance(x, TraverseResponse):
            e = _Elem(i, 'response')
            e.response = x
        elif isinstance(x, dict):
            keys = set(x)
            e = None
            if keys == {'response', 'index'} and isinstance(x['response'], TraverseResponse) \
                    and isinstance(x['index'], int) and not isinstance(x['index'], bool):
                e = _Elem(i, 'response')
                e.response, e.index = x['response'], x['index']
            elif keys == {'response', 'sequence'} \
                    and isinstance(x['response'], TraverseResponse) \
                    and isinstance(x['sequence'], str):
                e = _Elem(i, 'response')
                e.response, e.sequence = x['response'], x['sequence']
            elif 'pair' in keys and _pair_desc(x['pair']) is not None:
                rest = keys - {'pair'}
                if rest == {'graphlet'} and (isinstance(x['graphlet'],
                                                        (Graphlet, StitchedGraphlet))
                                             or C.is_view(x['graphlet'])):
                    e = _Elem(i, 'pair_graphlet')
                    e.graphlet = x['graphlet']
                elif 'result' in rest and rest <= {'result', 'error', 'sequence'} \
                        and isinstance(x['result'], dict) \
                        and (x.get('sequence') is None or isinstance(x['sequence'], str)):
                    e = _Elem(i, 'pair_result')
                    e.result = x['result']
                    e.error = _err_text(x.get('error')) if x.get('error') is not None else None
                    e.sequence = x.get('sequence')
                elif rest == {'error'} and isinstance(x['error'], (str, BaseException)):
                    e = _Elem(i, 'pair_error')
                    e.error = _err_text(x['error'])
                if e is not None:
                    e.stated = _pair_desc(x['pair'])
            if e is None:
                raise C.NotCombinable(TEXTS['not_combinable'] % ('result %d (%s)'
                                                                 % (i, _describe(x))))
        else:
            raise C.NotCombinable(TEXTS['not_combinable'] % ('result %d (%s)' % (i, _describe(x))))
        elems.append(e)
    return elems


def _joined_bytes(texts):
    """The bytes |texts| joined into one str hold beyond its header: one a character for
    ASCII; a single other character widens every character of the result, so then four
    (CPython's widest form)."""
    texts = list(texts)
    n = sum(map(len, texts))
    return n if all(t.isascii() for t in texts) else 4 * n


def _json_bound(x):
    """An upper bound of len(composite.canonical_json(x)), read without building it: a
    printable ASCII string its length, quotes and backslashes doubled; any other string 12
    characters per character (an escape per UTF-16 unit of an astral character)."""
    if isinstance(x, str):
        n = len(x)
        if x.isascii() and x.isprintable():
            return n + 2 + x.count('"') + x.count('\\')
        return 12 * n + 2
    if x is None or isinstance(x, bool):
        return 5
    if isinstance(x, (int, float)):
        return 24 if isinstance(x, float) or -10 ** 20 < x < 10 ** 20 else len(str(x))
    if isinstance(x, dict):
        return 2 + sum(_json_bound(k if isinstance(k, str) else str(k)) + _json_bound(v) + 2
                       for k, v in x.items())
    if isinstance(x, (list, tuple)):
        return 2 + sum(_json_bound(v) + 1 for v in x)
    return _json_bound(str(x))


def _elem_meta(e):
    """The canonical-JSON bound of an input's variable-length metadata."""
    n = _json_bound(e.error) + _json_bound(e.result) + _json_bound(e.stated)
    n += len(e.sequence or '') + len(e.sent or '')
    if e.graphlet is not None and not C.is_view(e.graphlet):
        n += _json_bound(_root(e.graphlet).envelope)
    if e.env is not None:
        n += _json_bound(e.env)
    return n


def _err_text(e):
    return e if isinstance(e, str) else str(e)


def _resolve_response(e):
    """not_one_seed; then the element's graphlet or result."""
    r = e.response
    n = len(r.graphlets)
    if e.index is None:
        if n == 0:
            raise _refuse('not_one_seed', 'not_one_seed_zero', e.pos)
        if n != 1:
            raise _refuse('not_one_seed', None, e.pos, n)
        idx = 0
    else:
        if not 0 <= e.index < n:
            raise _refuse('not_one_seed', 'not_one_seed_range', e.pos, e.index, n)
        idx = e.index
    e.env = r.envelope
    g = r.graphlets[idx]
    if g is not None:
        e.graphlet = g
    else:
        err = next((x for x in r.errors if x.get('index') == idx), None)
        res = (err or {}).get('result')
        if res is None:
            raw = (r.raw or {}).get('results') or []
            res = raw[idx] if idx < len(raw) else {}
        e.result = res if isinstance(res, dict) else {}
        e.error = _err_text((err or {}).get('error')) if (err or {}).get('error') else None
    sent = getattr(r, 'seeds_sent', None)
    if sent and idx < len(sent):
        s = sent[idx]
        e.sent = s.get('sequence') if isinstance(s, dict) else s


# ------------------------------------------------------------------ construction

def union_pairs(results, *, seed=None, expected_pairs=None, seqs_order=None, budget=None):
    """The cross-pair union of one seed's results (§2): a new PairUnion every time, its parts
    held by reference, never copied or modified.

    results: a list of Graphlet | StitchedGraphlet (experimental, §2.3) | TraverseResponse
    (one seed) | {response, index} | {response, sequence} | {pair, graphlet} | {pair,
    result[, error][, sequence]} | {pair, error}. seed: the oriented seed sequence;
    required only when no result carries a graphlet. expected_pairs: the pairs the request
    was sent to (graph names or pair objects); None: not known, and every statement is then
    over the pairs given.
    seqs_order: pair (a graph name, or a (graph, graph_path, annotation_path) tuple) ->
    'verified' | 'count_only' | None, as GET /capabilities states it; a pair not given
    reads 'unknown'. budget= (local limits): the trie pass is charged; a stop raises with no
    partial. Refusals (CombineRefused, .code) in the order of Appendix A.2."""
    return _B.run_scoped(budget, 'union_pairs', (),
                         lambda b: PairUnion._build(results, seed, expected_pairs, seqs_order, b))


def _check_args(results, seed, expected_pairs, seqs_order):
    if isinstance(expected_pairs, str):
        raise _refuse('bad_expected', 'bad_expected_str')
    exp = None
    if expected_pairs is not None:
        if not isinstance(expected_pairs, (list, tuple)):
            raise _refuse('bad_expected', None, repr(expected_pairs))
        exp = []
        for i, p in enumerate(expected_pairs):
            d = _pair_desc(p)
            if d is None:
                raise _refuse('bad_expected', None, 'expected_pairs[%d] (%r)' % (i, p))
            exp.append(d)
    so = []
    if seqs_order is not None:
        if not isinstance(seqs_order, dict):
            raise _refuse('bad_seqs_order', None, repr(seqs_order))
        for k, v in seqs_order.items():
            if isinstance(k, tuple) and 1 <= len(k) <= 3 and isinstance(k[0], str) and k[0] \
                    and all(x is None or isinstance(x, str) for x in k):
                d = {kk: vv for kk, vv in zip(('graph', 'graph_path', 'annotation_path'), k)
                     if vv}
            else:
                d = _pair_desc(k) if isinstance(k, str) else None
            if d is None or v not in _SEQS_ORDER:
                raise _refuse('bad_seqs_order', None, 'seqs_order[%r] = %r' % (k, v))
            so.append((d, 'none' if v is None else v, k))
        # two values for one pair, seen among the keys themselves (an argument error, so
        # before any result is read); a conflict only a result's pair shows is found once
        # the pairs are known (_seqs_of)
        for i in range(len(so)):
            for j in range(i + 1, len(so)):
                if so[i][1] != so[j][1] and _matches(so[i][0], so[j][0]):
                    raise _refuse('bad_seqs_order', None,
                                  'seqs_order[%r] and seqs_order[%r] (two values for one pair)'
                                  % (so[i][2], so[j][2]))
    if seed is not None and not isinstance(seed, str):
        raise TypeError('seed= is the oriented seed sequence (a string), not %r' % (seed,))
    if not isinstance(results, (list, tuple)):
        raise C.NotCombinable(TEXTS['not_combinable'] % _describe(results))
    return exp, so


def _absent_category(e):
    """The absent record's category, or a refusal: result_has_graphlet / no_graphlet."""
    r = e.result
    if 'graphlet' in r:
        raise _refuse('result_has_graphlet', None, e.pos)
    walks = (r.get('outcome') or {}).get('walks') if isinstance(r.get('outcome'), dict) \
        else None
    if r.get('arms') or walks not in _ABSENT_OF_WALKS:
        raise _refuse('no_graphlet', None, e.pos, walks)
    return _ABSENT_OF_WALKS[walks]


def _seed_subject(kind, tag):
    return {'part': "part %s's seed" % tag, 'seed': 'seed=',
            'absent': 'the seed of the absent record of pair %s' % tag}[kind]


def _seed_object(tag):
    return "part %s's" % tag if tag is not None else 'seed='


def _check_seed(ref_seq, ref_tag, seq, kind, tag, alphabet):
    """different_seed for |seq| (case mapped) against the reference: the first part's
    seed (|ref_tag|), or seed= when no part carries a graphlet (|ref_tag| None)."""
    diff = C.seed_difference(ref_seq, seq, alphabet)
    if diff is None:
        return
    subj = _seed_subject(kind, tag)
    obj = _seed_object(ref_tag)
    if diff[0] == 'revcomp':
        raise _refuse('different_seed', 'different_seed_revcomp', subj, obj,
                      {'part': 'part %s', 'absent': 'pair %s'}.get(kind, 'seed=%s') % (
                          tag if kind != 'seed' else ''),
                      'part %s' % ref_tag if ref_tag is not None else 'seed=')
    if diff[0] == 'length':
        raise _refuse('different_seed', 'different_seed_length', subj, diff[2], obj, diff[1],
                      diff[3])
    raise _refuse('different_seed', None, subj, obj, diff[1])


class PairUnion:
    """The union of one seed's results across pairs (§2.8). Every method takes budget=
    (keyword-only); the lists walks(), segments() and carriers() also resume=."""

    # ------------------------------------------------------------ build

    @classmethod
    def _build(cls, results, seed, expected_pairs, seqs_order, b):
        self = cls.__new__(cls)
        exp, so = _check_args(results, seed, expected_pairs, seqs_order)
        elems = _classify(results)
        # view_unsupported
        for e in elems:
            if e.graphlet is not None and C.is_view(e.graphlet):
                raise _refuse('view_unsupported')
        # not_one_seed (and each response's graphlet or result)
        for e in elems:
            if e.form == 'response':
                _resolve_response(e)
        # result_has_graphlet / no_graphlet; the kind of every element
        for e in elems:
            if e.form == 'pair_error':
                e.kind = 'failed'
            elif e.graphlet is not None:
                e.kind = 'part'
            else:
                e.kind = 'absent'
                if 'graphlet' in e.result:
                    raise _refuse('result_has_graphlet', None, e.pos)
        for e in elems:
            if e.kind == 'absent':
                e.cls = _absent_category(e)
        # no_seed
        if seed is None and not any(e.kind == 'part' for e in elems):
            raise _refuse('no_seed')
        # MissingEnvelope
        for e in elems:
            if e.kind == 'part':
                C.require_envelope(_root(e.graphlet), OP)
        # pair_mismatch: the stated pair against the envelope's (a stitched part's: its
        # root's; stitch() checked its parts against each other)
        for e in elems:
            if e.kind == 'part':
                g0 = _root(e.graphlet)
                env = _env_desc(g0.envelope, g0)
                if e.stated is not None:
                    if answered_pair(g0.envelope) is not None or any(
                            e.stated.get(k) and env.get(k) and e.stated[k] != env[k]
                            for k in ('index_fp', 'index_meta_fp', 'index_ns')):
                        if not _matches(e.stated, env):
                            raise _refuse('pair_mismatch', None, e.pos, _show(env),
                                          _show(e.stated))
                    d = dict(env)
                    for k in ('graph', 'graph_path', 'annotation_path'):
                        d.pop(k, None)
                    d.update(e.stated)
                    e.desc, e.by = d, 'caller'
                else:
                    e.desc, e.by = env, 'envelope'
            elif e.form == 'response':
                e.desc, e.by = _env_desc(e.env), 'envelope'
            else:
                e.desc, e.by = dict(e.stated), 'caller'
        # no_pair
        for e in elems:
            if not e.desc.get('graph') and not e.desc.get('index_ns'):
                raise _refuse('no_pair', None, e.pos)
        # ambiguous_pair: a description matching two descriptions that do not match
        descs = [(('result %d' % e.pos), e.desc) for e in elems] + \
            [('expected_pairs[%d]' % i, d) for i, d in enumerate(exp or [])]
        for name, d in descs:
            hits = [x for _, x in descs if x is not d and _matches(d, x)]
            for i in range(len(hits)):
                for j in range(i + 1, len(hits)):
                    if not _matches(hits[i], hits[j]):
                        raise _refuse('ambiguous_pair', None, _show(d), _show(hits[i]),
                                      _show(hits[j]))
        classes = _classes([d for _, d in descs])
        cls_of = {id(d): c for c in classes for d in c['members']}
        tags = _tags(classes)
        seqs = _seqs_of(classes, so, tags)
        self._b = b
        # no_bases, different_alphabet
        parts = [e for e in elems if e.kind == 'part']
        for e in parts:
            if isinstance(e.graphlet, StitchedGraphlet):
                bare = next((q for q in e.graphlet.parts if not C.has_bases(q.graphlet)), None)
                if bare is not None:
                    raise _refuse('no_bases', None, '%s (its stitched part %s)' % (
                        tags[id(cls_of[id(e.desc)])], part_name(bare.key)))
            elif not C.has_bases(e.graphlet):
                raise _refuse('no_bases', None, tags[id(cls_of[id(e.desc)])])
        if parts:
            ref = parts[0]
            ref_tag = tags[id(cls_of[id(ref.desc)])]
            for e in parts[1:]:
                if C.alphabet_key(_root(e.graphlet)) != C.alphabet_key(_root(ref.graphlet)):
                    raise _refuse('different_alphabet', None, tags[id(cls_of[id(e.desc)])],
                                  C.alphabet_key(_root(e.graphlet)), ref_tag,
                                  C.alphabet_key(_root(ref.graphlet)))
        if b is not None:
            # the inputs' variable-length metadata (error strings, absent results, stated
            # pairs, envelopes, the seed checks), before the seed checks case-map it and the
            # identities, sort keys and the digest serialize it (GPT review M4): _M_META
            # copies of its canonical JSON at its bound, and 1 lwu per 8 bytes of it; the
            # seed (seed=, else the first part's): its case-mapped copy, kept, and the
            # digest's JSON of it (the encoder's chunk and their join) and its UTF-8 bytes
            # (Opus fix review N2)
            meta = sum(map(_elem_meta, elems)) + _json_bound(exp) + _json_bound(
                [k for _, _, k in so])
            n = len(seed) if seed is not None else len(parts[0].graphlet.seed.sequence)
            b.charge((meta >> 3) + (n >> 6), _M_META * meta + 4 * (STR + n))
        # different_seed
        alphabet = parts[0].graphlet.alphabet if parts else None
        if parts:
            useed = parts[0].graphlet.seed.sequence
            ref_tag = tags[id(cls_of[id(parts[0].desc)])]
            if seed is not None:
                _check_seed(useed, ref_tag, C.case_mapped(seed, alphabet), 'seed', None,
                            alphabet)
        else:
            # no part: the alphabet is not known, so the case mapping of every DNA and
            # protein build (upper case) applies
            useed, ref_tag = seed.upper(), None
        for e in elems:
            tag = tags[id(cls_of[id(e.desc)])]
            if e.kind == 'part':
                if e is not parts[0]:
                    _check_seed(useed, ref_tag, e.graphlet.seed.sequence, 'part', tag, alphabet)
            elif e.kind == 'absent':
                n = _seed_meta(e.result).get('length_bp')
                if isinstance(n, int) and not isinstance(n, bool) and n != len(useed):
                    raise _refuse('different_seed', 'different_seed_length_absent',
                                  _seed_subject('absent', tag), n, _seed_object(ref_tag),
                                  len(useed))
                for sq in (e.sent, e.sequence):
                    if sq is not None:
                        sq = C.case_mapped(sq, alphabet) if alphabet else sq.upper()
                        _check_seed(useed, ref_tag, sq, 'absent', tag, alphabet)
        # duplicate_pair (and the no-ops)
        kept = _dedupe(elems, cls_of, tags)
        # not_expected
        if exp is not None:
            exp_ids = {id(d) for d in exp}
            for e in kept['order']:
                c = cls_of[id(e.desc)]
                if not any(id(m) in exp_ids for m in c['members']):
                    raise _refuse('not_expected', None, tags[id(c)])
        self._finish(kept, classes, cls_of, tags, exp, seqs, useed, alphabet, seed,
                     expected_pairs, b)
        return self

    def _finish(self, kept, classes, cls_of, tags, exp, seqs, useed, alphabet, seed,
                expected_pairs, b):
        if b is not None:
            n = sum(1 for e in kept['by_class'].values() if e.kind == 'part')
            # a part record per part, and per sub-part of a stitched part (its source)
            ns = sum(len(_sub_graphlets(e.graphlet)) for e in kept['by_class'].values()
                     if e.kind == 'part')
            b.charge(W_ROW * len(classes), _M_PART * ns + _M_RECORD * (
                len(classes) - n + sum(map(len, kept['failures'].values()))))
        parts, absent = [], []
        self._failures = []
        for c in sorted(classes, key=lambda c: _sort_key(_key(c['desc']))):
            tag = tags[id(c)]
            e = kept['by_class'].get(id(c))
            key = _key(c['desc'])
            if e is not None and e.kind == 'part':
                g = e.graphlet
                if isinstance(g, StitchedGraphlet):
                    # its digest is its retrieval identity (§2.3, C.2); stitch(g), no
                    # junction, states g's identity
                    ident = dict(g.identity) if g.junctions else _identity_statement(g.root)
                    parts.append(PairPart(tag, key, e.by, seqs[id(c)], g, ident,
                                          (key, C.canonical_json(g.digest)), g.digest))
                else:
                    ident = _identity_statement(g)
                    ri = C.retrieval_identity(g)
                    parts.append(PairPart(tag, key, e.by, seqs[id(c)], g, ident,
                                          (key, C.canonical_json(ri)), C.part_digest(g)))
            elif e is not None:
                if e.kind == 'absent':
                    detail = _absent_detail(e)
                    absent.append(PairAbsent(tag, key, e.by, e.cls, detail, seqs[id(c)],
                                             _seed_check(e)))
                else:
                    absent.append(PairAbsent(tag, key, e.by, 'failed',
                                             {'error': e.error}, seqs[id(c)], None))
            for f in sorted(kept['failures'].get(id(c), ()),
                            key=lambda f: f.error.encode('utf-8')):
                if f is not e:
                    self._failures.append({'tag': tag, 'error': f.error,
                                           'kept': 'answer' if e is not None and e.kind != 'failed'
                                           else 'failure'})
        self.parts = tuple(parts)
        self.absent = tuple(absent)
        given = {id(c) for c in classes if kept['by_class'].get(id(c)) is not None}
        exp_ids = {id(d) for d in (exp or [])}
        self._not_included = tuple(tags[id(c)] for c in sorted(
            classes, key=lambda c: _sort_key(_key(c['desc'])))
            if id(c) not in given and any(id(m) in exp_ids for m in c['members']))
        self.expected = None if exp is None else tuple(
            tags[id(c)] for c in sorted(classes, key=lambda c: _sort_key(_key(c['desc'])))
            if any(id(m) in exp_ids for m in c['members']))
        self.pairs = tuple(p.tag for p in self.parts) + tuple(a.tag for a in self.absent) \
            + self._not_included
        # the seed and alphabet as the canonical first part states them (the parts agree
        # base for base after case mapping and on the alphabet's key, not on its spelling)
        self._seed_seq = self.parts[0].graphlet.seed.sequence if self.parts else useed
        self.alphabet = self.parts[0].graphlet.alphabet if self.parts else alphabet
        self._order = {ch: i for i, ch in enumerate(C._letters(alphabet or ''))}
        self._part_of = {p.tag: p for p in self.parts}
        self._tag_index = {p.tag: i for i, p in enumerate(self.parts)}
        self._exp_given = expected_pairs is not None
        self._expected_raw = [C.canonical_json(d) for d in (exp or [])]
        self._sources(b)
        self._arms = {}
        for side in ARM_SIDES:
            if any(side in _root(p.graphlet).arms for p in self.parts):
                self._arms[side] = _Trie(self, side, b)
        self._number_segments()
        # each pair's record-mapping status (seqs_order) is in it: the warnings, carriers()
        # and evidence() read it, so a resume token does not cross a change of it
        self.digest = C.combined_digest([
            [[list(p.key), p.tag, p.pair_stated_by, p.retrieval[1], p.digest, p.seqs_order]
             for p in self.parts],
            [[list(a.key), a.tag, a.pair_stated_by, a.category, C.canonical_json(a.detail),
              a.seed_check, a.seqs_order] for a in self.absent],
            [[f['tag'], f['error'], f['kept']] for f in self._failures],
            self._seed_seq, sorted(self._expected_raw) if exp is not None else None])
        self._fig_cache = {}
        self._seg_cache = {}

    # ------------------------------------------------------------ basic properties

    @property
    def seed(self):
        """The common oriented seed: the first part's SeedInfo, or a stand-in built from
        seed= when no part carries a graphlet ({sequence, length_bp})."""
        if self.parts:
            return self.parts[0].graphlet.seed
        return _SeedStandIn(self._seed_seq)

    def arms(self):
        """The arms any part walked, left before right."""
        return [s for s in ARM_SIDES if s in self._arms]

    def __repr__(self):
        return '<PairUnion %d parts, %d absent, %d not included, digest %s>' % (
            len(self.parts), len(self.absent), len(self._not_included), self.digest[:12])

    # refused (not_offered)
    def dump(self, *a, **kw):
        raise _refuse('not_offered', 'not_offered_dump')

    def to_json(self, *a, **kw):
        raise _refuse('not_offered', 'not_offered_json')

    def compare(self, *a, **kw):
        raise _refuse('not_offered', None, 'compare')

    def save(self, *a, **kw):
        raise _refuse('not_offered', None, 'save')

    def support_by_depth(self, *a, **kw):
        raise _refuse('not_offered', None, 'support_by_depth')

    def next_request(self, *a, **kw):
        raise _refuse('not_offered', None, 'next_request')

    def memory_bytes(self):
        """What the union holds in the heap: the parts' models (a model given twice counted
        once), the absent records, its own index and what its queries keep (the walk
        figures per arm, the part segments' label intervals)."""
        return C.models_memory_bytes([self])

    # ------------------------------------------------------------ helpers

    def _side(self, arm):
        s = {'l': 'left', 'r': 'right'}.get(arm, arm)
        if s not in ARM_SIDES:
            raise ValueError("arm is 'left' or 'right' ('l', 'r'), not %r" % (arm,))
        return s

    def _trie(self, arm):
        return self._arms.get(self._side(arm))

    def _sources(self, b=None):
        """The sources the trie reads, one per retrieval: a plain part's graphlet, or each
        part of a stitched part, root first (Part order: a parent before its children). For
        a union without a stitched part source i is part i."""
        self._src, self._root_src = [], []
        for pi, p in enumerate(self.parts):
            g = p.graphlet
            self._root_src.append(len(self._src))
            if not isinstance(g, StitchedGraphlet):
                self._src.append(_Src(pi, p.tag, g, 0, (), None, None, None))
                continue
            index = {}
            js = {j.part: j for j in g.junctions}
            for sp in g.parts:
                parent = leaf = None
                if sp.key:
                    parent = index[sp.key[:-1]]
                    side, pid = sp.key[-1]
                    pg = self._src[parent].graphlet
                    pa = pg.arms[side]
                    derive.uses(b, pg, pa, 'paths')
                    leaf = derive.paths(pa)[pid].leaf
                index[sp.key] = len(self._src)
                src = _Src(pi, p.tag, sp.graphlet, sp.offset_bp, sp.key, parent, leaf,
                           js.get(sp.key))
                if sp.key:
                    # the walk-scope causes the stitched graphlet states (review minor 2),
                    # computed once here (W8 per walk, charged when computed)
                    ca = sp.graphlet.arms[sp.key[-1][0]]
                    src.wcauses = {q.id: tuple(g._walk_causes(sp.key, q.id, b))
                                   for q in derive.paths(ca)}
                    src.below = tuple(g._causes_below(src.junction, b))
                self._src.append(src)

    def _below_codes(self, tag):
        """A stitched part's junctions (by the child's part key) -> the walk-scope cause codes
        on the stitched walks through each (StitchedGraphlet.summary()'s walk_causes)."""
        return {s.key: s.below for s in self._src if s.key and s.tag == tag}

    def _resolve(self, t, pi, node):
        """(source, stop, below) of part |pi| on the union walk ending at trie node |node|:
        t.stop_of() of the part's root; where that stop is where a child continues the
        walk (a graft), the child's own stop below it, and so on down the chain."""
        si = self._root_src[pi]
        stop, below = t.stop_of(si, node)
        while True:
            c = t.graft.get((si, stop))
            if c is None:
                return si, stop, below
            stop, below = t.stop_of(c, node, stop)
            si = c

    def _number_segments(self):
        n = 0
        for side in ARM_SIDES:
            t = self._arms.get(side)
            if t is not None:
                n = t.number(n)

    def _seg_labels(self, part, a, b=None):
        """seg -> (label_intervals(), the end set as a frozenset, label_intervals_cut()) of
        a part segment, built once per union (counted in memory_bytes()). |b|: each segment
        read is charged once per call at its cold price (label_intervals()' own, and its
        end set), whether the union holds it or not."""
        # |part|: a source (a stitched part's sub-parts each have their own segments)
        key = (part.tag, a.side) if not getattr(part, 'key', ()) else (part.tag, a.side,
                                                                         part.key)
        cache = self._seg_cache.get(key)
        if cache is None:
            cache = self._seg_cache[key] = {}
        g = part.graphlet

        def info(sid):
            if b is not None:
                pw, pm = _ops._interval_prices(b, g, a)
                b.uses_of(a, 'union_seg_labels', pw[sid] + W_ELEM,
                          pm[sid] + SET + SET_ITEM * len(a.segments[sid].end) + TUPLE,
                          depth=sid)
            got = cache.get(sid)
            if got is None:
                # built without a budget: its price is charged above, to the union call's
                # budget, and an ambient one is not charged again (Opus fix review N1)
                with _B.unbudgeted():
                    got = cache[sid] = (_ops.label_intervals(g, a.side, sid),
                                        frozenset(a.segments[sid].end),
                                        _ops.label_intervals_cut(g, a.side, sid))
            return got
        return info

    def _pair_label(self, part, lab_id):
        lab = part.graphlet.labels[lab_id]
        return PairLabel(part.tag, lab.ref, lab.name, lab.kind)

    @staticmethod
    def _check_paths(paths, max_transitions=1, max_paths=1000):
        if paths != 'supported' or max_transitions != 1 or max_paths != 1000:
            if paths not in ('supported', 'inferred'):
                raise ValueError("paths is 'supported' or 'inferred', not %r" % (paths,))
            raise _refuse('not_offered', 'not_offered_inferred')

    # ------------------------------------------------------------ labels

    def labels(self, pair=None, *, budget=None):
        """Every label of the parts as PairLabel (canonical pair order, label ids in
        order), or one pair's (its tag). budget=: per label read, and each PairLabel."""
        if pair is not None and pair not in self._part_of:
            raise KeyError('no part has the tag %r (parts: %s)'
                           % (pair, ', '.join(self._part_of) or 'none'))
        return _B.run_scoped(budget, 'union_labels', (), lambda b: self._labels(pair, b))

    def _labels(self, pair, b):
        out = []
        for p in self.parts:
            if pair is not None and p.tag != pair:
                continue
            if b is not None:
                n = len(p.graphlet.labels)
                b.charge(W_ROW + W_ELEM * n, n * (_M_PAIR_LABEL + LIST_ITEM))
            out.extend(self._pair_label(p, l.id) for l in p.graphlet.labels)
        return out

    def label(self, selector, *, budget=None):
        """{'pair', 'ref'} | {'pair', 'name'} | {'name'} (when only one pair has it) ->
        PairLabel. {'id'} is refused: ids are local to each part. budget=: per label of
        the parts searched (every part's for {'name'})."""
        return _B.run_scoped(budget, 'union_label', (), lambda b: self._label(selector, b))

    def _label(self, selector, b):
        if isinstance(selector, PairLabel):
            selector = {'pair': selector.pair, 'ref': selector.ref}
        if not isinstance(selector, dict):
            raise BadSelector('a pair union selects labels by {pair, ref}, {pair, name} or '
                              '{name}, not %r' % (selector,))
        if 'id' in selector:
            raise BadSelector(TEXTS['label_id'])
        keys = set(selector)
        if keys == {'name'}:
            name = selector['name']
            if b is not None:
                n = sum(len(p.graphlet.labels) for p in self.parts)
                b.charge(W_ELEM * n, _M_PAIR_LABEL * len(self.parts))
                for p in self.parts:
                    _label_index_price(b, p.graphlet)
            hits = [(p, l) for p in self.parts for l in p.graphlet.labels if l.name == name]
            tags = list(dict.fromkeys(p.tag for p, _ in hits))
            if not hits:
                raise UnknownLabel('no label named %r in any part' % (name,))
            if len(tags) > 1:
                raise _refuse('ambiguous_name', None, name, ', '.join(tags))
            p = hits[0][0]
            return self._pair_label(p, _label_in(p.graphlet, {'name': name}).id)
        if keys in ({'pair', 'ref'}, {'pair', 'name'}):
            p = self._part_of.get(selector['pair'])
            if p is None:
                raise UnknownLabel('no part has the tag %r' % (selector['pair'],))
            k = 'ref' if 'ref' in keys else 'name'
            if b is not None:
                b.charge(W_ELEM * len(p.graphlet.labels), _M_PAIR_LABEL)
                _label_index_price(b, p.graphlet)
            return self._pair_label(p, _label_in(p.graphlet, {k: selector[k]}).id)
        raise BadSelector('a pair union selects labels by {pair, ref}, {pair, name} or {name}, '
                          'not %r' % (selector,))

    def shared_names(self, *, budget=None):
        """name -> the tags of the parts where it occurs, for names seen in two or more
        parts (names in byte order). budget=: per label, its name set and the dicts."""
        return _B.run_scoped(budget, 'union_shared_names', (),
                             lambda b: self._shared_names(b))

    def _shared_names(self, b):
        """shared_names(); with |b| its charge: per label an element, and at its bound a
        name set entry, a sort key and the dict's lists (_M_NAME)."""
        if b is not None:
            n = sum(len(p.graphlet.labels) for p in self.parts)
            b.charge(W_ELEM * n, n * _M_NAME + DICT)
        seen = collections.defaultdict(list)
        for p in self.parts:
            for n in sorted({l.name for l in p.graphlet.labels}):
                seen[n].append(p.tag)
        # UTF-8 byte order is code point order: the names sort as they are, with no copy of
        # each as its bytes (Opus fix review N4)
        return {n: tuple(t) for n, t in sorted(seen.items(), key=lambda x: x[0])
                if len(t) > 1}

    def _mapping_warnings(self):
        """[(tag, 'count_only' | 'unknown')] of the parts whose header labels rest on an
        unchecked record mapping (§2.5, OQ20)."""
        out = []
        for p in self.parts:
            if p.seqs_order in ('count_only', 'unknown') \
                    and any(l.kind == 'header' for l in p.graphlet.labels):
                out.append((p.tag, p.seqs_order))
        return out

    def _mapping_sentence(self):
        w = self._mapping_warnings()
        if not w:
            return None
        groups = collections.OrderedDict()
        for t, v in w:
            groups.setdefault(v, []).append(t)
        if len(w) == len(self.parts) and len(groups) == 1:
            v = next(iter(groups))
            return 'header names of every pair rest on an unchecked record mapping (%s)' % v
        return 'header names rest on an unchecked record mapping in ' + ' and '.join(
            '%s (%s)' % (', '.join(ts), v) for v, ts in groups.items())

    # ------------------------------------------------------------ coverage

    def _default_depth(self):
        r = 0
        for p in self.parts:
            # a stitched part: its root's request (stitching certifies nothing past the
            # root's complete_to_bp, §3.5)
            g = _root(p.graphlet)
            bounds = ((g.envelope or {}).get('strategy') or {}).get('bounds') or {}
            v = bounds.get('max_extension_bp')
            if isinstance(v, int) and not isinstance(v, bool):
                r = max(r, v)
            else:
                r = max([r] + [a.complete_to_bp for a in g.arms.values()])
        return r

    def coverage(self, arm=None, at_bp=None, *, budget=None):
        """Every pair in one category at depth |at_bp| (default: the largest radius among
        the parts), §2.6. arm None: an arm-less union only (no part) -- or both arms' absent
        categories. budget=: per pair and limitation read, the entries and the sentence."""
        side = None if arm is None else self._side(arm)
        if side is None and self.parts:
            raise ValueError('name an arm (left or right): the union has parts')
        self._depth_of(at_bp)
        return _B.run_scoped(budget, 'union_coverage', (),
                             lambda b: self._coverage(side, at_bp, b))

    def _depth_of(self, at_bp):
        depth = self._default_depth() if at_bp is None else at_bp
        if not isinstance(depth, int) or isinstance(depth, bool) or depth < 0:
            raise ValueError('at_bp is a depth in bp (an int >= 0), not %r' % (at_bp,))
        return depth

    def _coverage(self, side, at_bp, b):
        depth = self._depth_of(at_bp)
        if b is not None:
            # per pair its entry and the limitations read; the sentence and the strings it
            # is joined from, at their bound (every tag at most three times, each note)
            n = len(self.pairs)
            nl = sum(len(_root(p.graphlet).limitations) for p in self.parts)
            b.charge(W_ROW * n + W_ELEM * nl, n * _M_COV_PAIR + 4 * (
                STR + _M_SENTENCE + 3 * _joined_bytes(self.pairs) + n * _M_COV_NOTE))
        per = collections.OrderedDict()
        for p in self.parts:
            per[p.tag] = _part_category(p, side, depth, self._below_codes(p.tag))
        for a in self.absent:
            per[a.tag] = _absent_entry(a)
        for t in self._not_included:
            per[t] = {'category': 'not_included', 'certified_to_bp': None}
        by = collections.OrderedDict((c, []) for c in CATEGORIES)
        for t, v in per.items():
            by[v['category']].append(t)
        by = collections.OrderedDict((c, ts) for c, ts in by.items() if ts)
        every = None
        if self.expected is not None and all(
                v['category'] in ('complete', 'seed_not_in_graph', 'labels_not_in_graph')
                for v in per.values()):
            ds = [v['certified_to_bp'] for t, v in per.items() if v['category'] == 'complete']
            every = min(ds) if ds else None
        sentence = self._sentence(side, depth, per, by)
        return Coverage(side, depth, dict(per), dict(by), every, sentence)

    def _settings(self, p):
        g = p.graphlet
        return '%s, %s, k %d' % (g.mode, g.support, g.k)

    def _sentence(self, side, depth, per, by):
        head = ('%s arm' % side) if side else 'both arms'
        settings = {self._settings(p) for p in self.parts}
        uniform = len(settings) == 1
        head += ', to %d bp' % depth if side else ''
        if uniform:
            head += ' (%s, in every pair)' % next(iter(settings))
        parts_txt = []
        for cat, tags in by.items():
            items = []
            for t in tags:
                v = per[t]
                s = t
                if not uniform and t in self._part_of:
                    s += ' [%s]' % self._settings(self._part_of[t])
                q = v.get('note')
                if q:
                    s += ' (%s)' % q
                items.append(s)
            parts_txt.append('%s %s' % (_CATEGORY_WORDS[cat], ', '.join(items)))
        out = head + ': ' + ('; '.join(parts_txt) if parts_txt else 'no pair') + '.'
        stated = [t for t in self.pairs if self._stated_by(t) == 'caller']
        if stated:
            out += ' Pairs stated by the caller: %s.' % ', '.join(stated)
        checks = collections.OrderedDict((('sent', []), ('stated', []), ('unchecked', [])))
        for a in self.absent:
            if a.seed_check is not None:
                checks[a.seed_check].append(a.tag)
        words = {'sent': 'sent by this client for', 'stated': 'stated by the caller for',
                 'unchecked': 'not checked for'}
        got = ['%s %s' % (words[k], ', '.join(v)) for k, v in checks.items() if v]
        if got:
            out += ' Seed of the absent records: %s.' % ', '.join(got)
        if not self._exp_given:
            out += ' ' + OVER_GIVEN
        return out + ' ' + TURN_OFFS + ' ' + SCOPE

    def _stated_by(self, tag):
        p = self._part_of.get(tag)
        if p is not None:
            return p.pair_stated_by
        for a in self.absent:
            if a.tag == tag:
                return a.pair_stated_by
        return None

    # ------------------------------------------------------------ statements

    def statements(self, arm=None, at_bp=None, *, paths='supported', budget=None):
        """The five statements (§2.1) for |arm| (None: an arm-less union, or the first arm
        any part walked) at |at_bp|."""
        self._check_paths(paths)
        side = None if arm is None else self._side(arm)
        return _B.run_scoped(budget, 'union_statements', (),
                             lambda b: self._statements(side, at_bp, paths, b))

    def _statements(self, side, at_bp, paths, b):
        if side is None and self.parts:
            side = self.arms()[0]
        cov = self._coverage(side, at_bp, b)
        return Statements(side, cov.at_bp, paths, cov, self._provenance_sentence(side, b),
                          TURN_OFFS, self._labels_sentence(b), SCOPE)

    def _provenance_sentence(self, side, b=None):
        sides = [side] if side else []
        items = []
        for s in sides:
            t = self._arms.get(s)
            if t is None:
                items.append('%s 0' % s)
                continue
            ws = self._figures(s, b)
            walks = [w for w in ws if w.kind == 'walk']
            stubs = len(ws) - len(walks)
            multi = sum(1 for w in walks if len(w.held_in) > 1)
            extra = ['%d held by two or more pairs' % multi]
            if stubs:
                extra.append('%d stub walks not listed' % stubs)
            items.append('%s %d (%s)' % (s, len(walks), '; '.join(extra)))
        if not items:
            return 'walks: none (no part carries a graphlet)'
        return 'walks: ' + ', '.join(items)

    def _labels_sentence(self, b=None):
        n = len(self._shared_names(b))
        if b is not None:
            # shared_names()' sets and dicts, dropped once counted; the sentence, which names
            # the tags of the parts its mapping warning covers, twice while built
            b.release(sum(len(p.graphlet.labels) for p in self.parts) * _M_NAME + DICT)
            b.charge(0, 2 * (STR + 600 + _joined_bytes(p.tag for p in self.parts)))
        s = ('labels: each pair\'s own, shown as <tag>/<name>; equal names in two pairs are not '
             'claimed to be one record (%s)' % (
                 'no name occurs in more than one pair' if n == 0 else
                 '1 name occurs in more than one pair' if n == 1 else
                 '%d names occur in more than one pair' % n))
        m = self._mapping_sentence()
        if m:
            s += '; ' + m
        return s

    def __str__(self):
        """The statements as sentences: print(union). Under an ambient budget
        (local_budget()) it is charged as statements() is."""
        return _B.run_scoped(None, 'union_str', (), self._text)

    def _text(self, b):
        n_failed = sum(1 for a in self.absent if a.category == 'failed')
        if self.expected is not None:
            head = 'union of one seed (%d bp) over %d expected pairs: ' % (
                len(self._seed_seq), len(self.expected))
        else:
            head = 'union of one seed (%d bp) over %d pairs given: ' % (
                len(self._seed_seq), len(self.pairs))
        head += '%d parts, %d absent records' % (len(self.parts), len(self.absent))
        if n_failed:
            head += ' (%d failed)' % n_failed
        if self._not_included:
            head += ', %d not included' % len(self._not_included)
        lines = [head]
        sides = self.arms() or [None]
        for s in sides:
            lines.append(_core_sentence(self._coverage(s, None, b)))
        if self.parts:
            prov = []
            for s in self.arms():
                prov.append(self._provenance_sentence(s, b)[len('walks: '):])
            lines.append('walks: ' + ', '.join(prov))
        lines.append(self._labels_sentence(b))
        lines.append(TURN_OFFS)
        lines.append(SCOPE)
        return '\n'.join(lines)

    # ------------------------------------------------------------ segments

    def segments(self, arm, *, budget=None, resume=None):
        """The union segments of |arm| in id order (trie preorder, children by base), each
        with its support per pair: the part segments it covers, the walks and stubs passing
        and the labels alive there as intervals (from the parts' alive sets, linear)."""
        side = self._side(arm)
        t = self._arms.get(side)
        args = [side]
        start = C.resume_pos(self.digest, 'segments', args, resume, 1)[0]
        if t is None:
            return []
        b = _B.resolve(budget)
        out = []
        try:
            if b is not None:
                with b.scope('segments'):
                    for sid in t.order[start:]:
                        self._charge_segment(b, t, sid)
                        out.append(self._segment(t, sid, b))
            else:
                for sid in t.order[start:]:
                    out.append(self._segment(t, sid))
        except LocalBudgetExceeded as e:
            raise e.restate(partial=Partial(out, C.token(self.digest, 'segments', args,
                                                           start + len(out)), len(t.order)))
        return out

    def _charge_segment(self, b, t, sid):
        """A union segment's row: per covering part its support and a LabelInterval per
        interval of the part segments it covers (at most their entry labels and changes,
        or P-run labels); the intervals read are charged by _seg_labels()."""
        us = t.useg[sid]
        n = 0
        for pi in us.cover:
            part = self._src[pi]
            info = self._seg_labels(part, part.graphlet.arms[t.side], b)
            for x in us.pieces[pi]:
                n += len(info(x[0])[0])
        b.charge(W_ROW + W_ELEM * n, 400 + 300 * len(us.cover) + n * _M_INTERVAL)

    def _segment(self, t, sid, b=None):
        us = t.useg[sid]
        sup = []
        # sources in order: a part's sources are consecutive and at most one of them
        # spells a union segment, so this is the canonical pair order
        for pi in sorted(us.cover):
            part = self._src[pi]
            off = part.offset
            info = self._seg_labels(part, part.graphlet.arms[t.side], b)
            ivs = []
            cut = False
            for seg, f, to in us.pieces[pi]:
                for lab, a, z in info(seg)[0]:
                    lo, hi = max(a + off, f), min(z + off, to)
                    if lo < hi:
                        ivs.append(LabelInterval(self._pair_label(part, lab), lo, hi))
                cut = cut or info(seg)[2]
            last = us.pieces[pi][-1][0]
            sup.append(SegmentSupport(part.tag, t.below_walks[pi][last], t.below_stubs[pi][last],
                                      tuple(us.pieces[pi]), tuple(ivs), cut, part.key))
        return UnionSegment(us.id, t.side, us.from_bp, us.to_bp,
                            None if us.parent is None else t.useg[us.parent].id,
                            tuple(t.useg[c].id for c in us.children), tuple(sup))

    # ------------------------------------------------------------ walks

    def _figures(self, side, b=None):
        """Per union walk of |side|, in id order: (kind, length_bp, whole_in, held_in), from
        one pass over the trie. Each part's stop is carried down the union segments, so a
        walk's figures cost one step per part; only a merges span follows its route (the
        walk's remaining length) to tell whether the part holds the walk. Kept per side
        (memory_bytes() counts it); charged at its cold price on every call."""
        t = self._arms[side]
        n = len(self.parts)
        ns = len(self._src)
        if b is not None:
            # per union segment a stop per source, per walk a step per source and, for a
            # part that merges, its route along at most the walk's length (one span at a
            # time); the derivations a span reads; the stops are dropped at the end, the
            # figures kept
            w = W_ELEM * ns * len(t.order)
            for node in t.markers:
                w += W_ROW + ns * (W_ELEM + (t.depth[node] >> 4))
            for p in self._src:
                a = p.graphlet.arms.get(side)
                if a is not None:
                    derive.uses(b, p.graphlet, a, 'first_paths', *(
                        ('end_labels',) if p.graphlet.mode == 'constrain' else ()))
            tmp = len(t.order) * (TUPLE + 8 * ns + LIST_ITEM) + _M_SPAN
            b.charge(w, tmp + len(t.markers) * (_M_FIG + 16 * n))
        got = self._fig_cache.get(side)
        if got is None:
            got = self._build_figures(t, side, n)
        if b is not None:
            b.release(tmp)
        return got

    def _build_figures(self, t, side, n):
        has_arm = [side in _root(p.graphlet).arms for p in self.parts]
        ns = len(self._src)
        # a stop per source; a child source's is None until it spells a segment (its stop
        # is then where it continues the walk, PairUnion._resolve())
        root = tuple(0 if s.parent is None else None for s in self._src)
        stops = [None] * len(t.useg)
        for u in t.order:
            us = t.useg[u]
            up = root if us.parent is None else stops[us.parent]
            stops[u] = tuple(us.end if si in us.cover else up[si] for si in range(ns))
        out = []
        for wid, node in enumerate(t.markers):
            u = t.node_in.get(node)
            st = root if u is None else stops[u]
            length = t.depth[node]
            whole, held = [], []
            for pi in range(n):
                if not has_arm[pi]:
                    continue
                si = self._root_src[pi]
                stop = st[si]
                while True:
                    c = t.graft.get((si, stop))
                    if c is None:
                        break
                    si, stop = c, (stop if st[c] is None else st[c])
                how, _ = _case(t, si, stop, length)
                if how == 'whole':
                    whole.append(self.parts[pi].tag)
                    held.append(self.parts[pi].tag)
                elif how == 'extends':
                    held.append(self.parts[pi].tag)
                elif how == 'merges':
                    sp = _span(self, t, pi, self.parts[pi], node, length, False)
                    if sp.route.how == 'held':
                        held.append(self.parts[pi].tag)
            kind = 'walk' if any(t.walk_end.get(node, ())) else 'stub'
            out.append(_WalkFigures(kind, length, tuple(whole), tuple(held)))
        self._fig_cache[side] = out
        return out

    def _charge_walk(self, b, t, node, with_labels):
        """A union walk built whole: per part its span (the segments past its stop, at most
        the walk's) and route (at most the walk's length) and, with labels, the alive sets
        read along its chain."""
        d = t.nsegs(node)
        n = len(self.parts)
        ns = len(self._src)
        w = W_ROW + ns * (W_ELEM * (1 + d) + (t.depth[node] >> 4))
        if with_labels:
            w += ns * W_ELEM * 4 * (1 + d)
        # the walk, its path (a slot and an id per union segment) and the spans' tuple; each
        # span charges itself (_span)
        b.charge(w, _walk_bytes(d, n))

    def _walk(self, t, wid, node, with_labels=False, spell=False, at_bp=None, b=None):
        path = t.path_usegs(node)
        length = t.depth[node]
        kind = 'walk' if any(t.walk_end.get(node, ())) else 'stub'
        spans = []
        for pi, part in enumerate(self.parts):
            spans.append(_span(self, t, pi, part, node, length, with_labels, at_bp, b))
        whole = tuple(s.pair for s in spans if s.how == 'whole')
        held = tuple(s.pair for s in spans if s.how in ('whole', 'extends')
                     or (s.how == 'merges' and s.route is not None and s.route.how == 'held'))
        seq = self._spell_nodes(t, path, 'natural', True) if spell else None
        return UnionWalk(wid, t.side, kind, length, seq, tuple(t.useg[u].id for u in path),
                         tuple(spans), whole, held)

    def _walks_full(self, side, b, with_labels, ids=None, stubs=True, at_bp=None):
        """The union walks of |side| built whole on demand, one at a time in id order (only
        |ids| when given; stub walks only with |stubs|; span labels at |at_bp|, default each
        span's end), each charged before it is built; nothing is kept."""
        t = self._arms[side]
        figs = self._figures(side, b)
        want = None if ids is None else set(ids)
        prev = 0
        for wid, node in enumerate(t.markers):
            if want is not None and wid not in want:
                continue
            if not stubs and figs[wid].kind != 'walk':
                continue
            if b is not None:
                self._charge_walk(b, t, node, with_labels)
            w = self._walk(t, wid, node, with_labels, at_bp=at_bp, b=b)
            yield w
            # the consumer has moved on: the walk before this one is gone (this one is
            # still bound to its loop variable while the next is built)
            if b is not None:
                b.release(prev)
                prev = _walk_bytes(len(w.segments), len(w.spans)) + sum(map(_span_bytes,
                                                                            w.spans))

    def walks(self, arm, *, paths='supported', max_transitions=1, max_paths=1000, top=None,
              by='pairs', pairs=None, labels=None, min_bp=0, spell=False, stubs=False,
              budget=None, resume=None):
        """The union walks of |arm| (§2.8). by: 'pairs' (the pairs holding it whole, then
        displaying it, then length, then id), 'length' or 'id'. pairs=: the walks every
        named pair holds whole (held_in); labels=: the walks those labels carry whole
        (from_bp 0) in some pair. Stub walks only with stubs=True. Every span's labels are
        None here (carriers() lists them)."""
        self._check_paths(paths, max_transitions, max_paths)
        if by not in ('pairs', 'length', 'id'):
            raise ValueError("by is 'pairs', 'length' or 'id', not %r" % (by,))
        side = self._side(arm)
        b = _B.resolve(budget)
        if b is None:
            return self._walks(side, by, top, pairs, labels, min_bp, spell, stubs, None,
                               resume)
        with b.scope('walks'):
            return self._walks(side, by, top, pairs, labels, min_bp, spell, stubs, b, resume)

    def _walks(self, side, by, top, pairs, labels, min_bp, spell, stubs, b, resume):
        if b is not None and (pairs or labels):
            # the selectors, caller-sized, bound into the arguments a resume token is
            # checked against: their repr, its JSON and UTF-8 bytes in the token's digest, at
            # twice their JSON bound (Opus fix review N4)
            n = 2 * (_json_bound(list(pairs or ())) + _json_bound(list(labels or ())))
            b.charge(n >> 6, 4 * (STR + n))
        args = [side, by, top, sorted(pairs) if pairs else None,
                [repr(x) for x in labels] if labels else None, min_bp, spell, stubs]
        start = C.resume_pos(self.digest, 'walks', args, resume, 1)[0]
        if side not in self._arms:
            return []
        out = []
        try:
            return self._walks_body(self._arms[side], by, top, pairs, labels, min_bp, spell,
                                    stubs, b, start, out)
        except LocalBudgetExceeded as e:
            raise e.restate(partial=Partial(out, C.token(self.digest, 'walks', args,
                                                           start + len(out))))

    def _walks_body(self, t, by, top, pairs, labels, min_bp, spell, stubs, b, start, out):
        # the walks are selected and ranked on their figures (one trie pass); only the
        # walks returned are built, with their spans and path
        figs = self._figures(t.side, b)
        if b is not None:
            # the selection, and the sort's key per walk (a tuple of four, its ints)
            b.charge(W_ELEM * len(figs), len(figs) * (LIST_ITEM + TUPLE + 4 * 8 + 2 * INT))
        sel = [wid for wid, f in enumerate(figs)
               if (stubs or f.kind == 'walk') and f.length_bp >= min_bp]
        if pairs:
            for p in pairs:
                if p not in self._part_of:
                    raise KeyError('no part has the tag %r' % (p,))
            sel = [wid for wid in sel if all(p in figs[wid].held_in for p in pairs)]
        if labels is not None:
            want = {self._label(x, b) for x in labels}
            # only a span that holds the walk carries it: its own labels (whole, extends) or
            # its route's (merges, held to the walk's end); a stub's labels before its merge
            # and the labels of a span that ends or turns off carry nothing to the end
            keep = []
            for wid in sel:
                node = t.markers[wid]
                held = figs[wid].held_in
                if b is not None:
                    self._charge_walk(b, t, node, True)
                hit = False
                for pi, part in enumerate(self.parts):
                    if part.tag not in held:
                        continue
                    sp = _span(self, t, pi, part, node, figs[wid].length_bp, True, None, b)
                    hit = any(sl.from_bp == 0 and sl.label in want
                              for sl in _carrying_labels(sp))
                    if b is not None:
                        b.release(_span_bytes(sp))
                    if hit:
                        break
                if b is not None:
                    b.release(_walk_bytes(t.nsegs(node), len(self.parts)))
                if hit:
                    keep.append(wid)
            sel = keep
        if by == 'pairs':
            sel.sort(key=lambda i: (-len(figs[i].held_in), -len(figs[i].whole_in),
                                    -figs[i].length_bp, i))
        elif by == 'length':
            sel.sort(key=lambda i: (-figs[i].length_bp, i))
        if b is not None:
            b.release(len(figs) * (TUPLE + 4 * 8 + 2 * INT))
        if top is not None:
            sel = sel[:top]
        for wid in sel[start:]:
            node = t.markers[wid]
            length = figs[wid].length_bp
            if b is not None:
                self._charge_walk(b, t, node, False)
                b.charge(W_ROW + (length >> BASE_SHIFT if spell else 0),
                         400 + (STR + length if spell else 0))
            out.append(self._walk(t, wid, node, False, spell, b=b))
        return out

    def _walk_by_id(self, side, walk):
        if isinstance(walk, str) and walk.startswith('i:'):
            raise _refuse('not_offered', 'not_offered_inferred')
        if isinstance(walk, UnionWalk):
            walk = walk.id
        if isinstance(walk, bool) or not isinstance(walk, int):
            raise BadSelector('a union walk is selected by its id (an int), not %r' % (walk,))
        t = self._arms.get(side)
        if t is None or not 0 <= walk < len(t.markers):
            raise KeyError('the %s arm has no union walk %r' % (side, walk))
        return t, t.markers[walk]

    def _spell_nodes(self, t, path, orientation, with_seed):
        pieces = []
        for u in path:
            us = t.useg[u]
            pi = min(us.cover)
            g = self._src[pi].graphlet
            off = self._src[pi].offset
            for seg, f, to in us.pieces[pi]:
                pieces.append((g, t.side, seg, f - off, to - off))
        w = C.spell_pieces(pieces, t.side, 'walk')
        seed = self._seed_seq
        if orientation == 'natural':
            flank = w if t.side == 'right' else w[::-1]
            if not with_seed:
                return flank
            return seed + flank if t.side == 'right' else flank + seed
        if orientation != 'walk':
            raise ValueError("orientation is 'natural' or 'walk'")
        if not with_seed:
            return w
        return (seed if t.side == 'right' else seed[::-1]) + w

    def spell(self, arm, walk, orientation='natural', with_seed=False, *, budget=None):
        """The bases of union walk |walk| of |arm| from the parts' bases (natural
        orientation: the left arm reversed as Graphlet.spell() writes it)."""
        side = self._side(arm)
        t, node = self._walk_by_id(side, walk)
        b = _B.resolve(budget)
        if b is not None:
            with b.scope('spell'):
                n = t.depth[node] + len(self._seed_seq)
                b.charge(W_ROW + W_STEP * t.nsegs(node) + (n >> 8), 3 * (STR + n))
                return self._spell_nodes(t, t.path_usegs(node), orientation, with_seed)
        return self._spell_nodes(t, t.path_usegs(node), orientation, with_seed)

    def carriers(self, arm, walk, *, at_bp=None, budget=None, resume=None):
        """Which pairs and labels carry union walk |walk| (§2.8): per part, how far it
        holds the walk (displayed or through its route), how it stops, with its absence and
        mode, and its span labels with their from_bp (at |at_bp|, default each span's
        end). A Carriers list; str() prints one line per pair."""
        side = self._side(arm)
        t, node = self._walk_by_id(side, walk)
        args = [side, walk if isinstance(walk, int) else getattr(walk, 'id', None), at_bp]
        start = C.resume_pos(self.digest, 'carriers', args, resume, 1)[0]
        length = t.depth[node]
        b = _B.resolve(budget)
        out = Carriers()
        out.mapping_note = self._mapping_sentence()
        try:
            if b is not None:
                with b.scope('carriers'):
                    for pi in range(start, len(self.parts)):
                        b.charge(W_ROW + W_STEP * length // 16 + W_ELEM * t.nsegs(node),
                                 LIST_ITEM)
                        out.append(_span(self, t, pi, self.parts[pi], node, length, True,
                                         at_bp, b))
            else:
                for pi in range(start, len(self.parts)):
                    out.append(_span(self, t, pi, self.parts[pi], node, length, True, at_bp))
        except LocalBudgetExceeded as e:
            raise e.restate(partial=Partial(list(out), C.token(self.digest, 'carriers', args,
                                                                 start + len(out)),
                                            len(self.parts)))
        return out

    # ------------------------------------------------------------ exports

    def to_fasta(self, arm=None, walks=None, with_seed=True, orientation='natural', width=None,
                 *, paths='supported', max_transitions=1, max_paths=1000, at_bp=None,
                 stubs=False, budget=None):
        """One record per union walk (§2.9), its header stating its provenance, coverage and
        scope as key=value tokens (Appendix A.5). at_bp: the labels= counts are the span
        labels at that depth, as carriers(at_bp=) gives them (default each span's end)."""
        self._check_paths(paths, max_transitions, max_paths)
        from . import export
        return _B.run_scoped(budget, 'union_to_fasta', (), lambda b: export.union_fasta(
            self, arm, walks, with_seed, orientation, width, at_bp, stubs, b))

    def to_gfa(self, *, paths='supported', with_seed=True, labels='counts', at_bp=None,
               budget=None):
        """GFA 1.0 on the seed strand (§2.9), the statements first as comment lines; in the
        source-supported mode one tree per arm, the seed written once per arm, a stub's
        merge a tag (mg, mt), not a link: every path the file admits is spelled whole by
        one pair."""
        self._check_paths(paths)
        if labels not in ('counts', 'full'):
            raise ValueError("labels is 'counts' or 'full', not %r" % (labels,))
        from . import export
        return _B.run_scoped(budget, 'union_to_gfa', (), lambda b: export.union_gfa(
            self, with_seed, labels, at_bp, b))

    # ------------------------------------------------------------ summary, evidence

    def summary(self, max_bytes=4096, *, paths='supported', budget=None):
        """Agent-facing, <= max_bytes of JSON: the statements first and whole (refused
        summary_too_small rather than cut), then the parts and absent records, then per arm
        the walk counts and the walks held by the most pairs, as long as they fit."""
        self._check_paths(paths)
        return _B.run_scoped(budget, 'union_summary', (),
                             lambda b: self._summary(max_bytes, b))

    def _summary(self, max_bytes, b):
        if b is not None:
            # the answer (at most max_bytes of JSON as dicts and lists) and the JSON text
            # each fit test writes, two at a time; a record per pair
            b.charge(W_ROW * len(self.pairs), 8 * max_bytes + 600 * len(self.pairs))
        sts = []
        for s in (self.arms() or [None]):
            st = self._statements(s, None, 'supported', b)
            sts.append({'arm': s, 'coverage': _core_sentence(st.coverage),
                        'provenance': st.provenance})
        base = {'statements': {'per_arm': sts, 'turn_offs': TURN_OFFS,
                               'labels': self._labels_sentence(b), 'scope': SCOPE}}

        # every fit test is of the answer as it would be returned, its notes included: a
        # section is added only where the answer still fits with the notes of what may be
        # left out after it (GPT review M6)
        def size(o, notes):
            if notes:
                o = dict(o, notes=list(notes))
            return len(C.canonical_json(o))

        pairs = {p.tag: _summary_part(p, self._below_codes(p.tag)) for p in self.parts}
        pairs.update({a.tag: {'absent': a.category, 'seed_check': a.seed_check}
                      for a in self.absent})
        pairs.update({t: {'not_included': True} for t in self._not_included})
        arms = {}
        ranked = []
        for s in self.arms():
            figs = self._figures(s, b)
            walks = [w for w in figs if w.kind == 'walk']
            arms[s] = {'walks': len(walks), 'stub_walks': len(figs) - len(walks),
                       'segments': len(self._arms[s].order), 'top': []}
            order = sorted((i for i, w in enumerate(figs) if w.kind == 'walk'),
                           key=lambda i: (-len(figs[i].held_in), -len(figs[i].whole_in),
                                          -figs[i].length_bp, i))
            ranked += [(s, {'id': i, 'length_bp': figs[i].length_bp,
                            'held_in': list(figs[i].held_in)}) for i in order[:10]]
        no_pairs, no_arms = 'pairs left out (max_bytes)', 'arms left out (max_bytes)'
        cut = ['top walks cut (max_bytes)'] if ranked else []

        def least_after(o, notes):
            # the smallest answer that goes on from |o|: the arms with no top walk listed
            # (and the note of the cut, where a walk is ranked), or the note that they are
            # left out -- for a union without arms, its empty arms (Opus fix review N7)
            return min(size(dict(o, arms=arms), notes + cut), size(o, notes + [no_arms]))

        with_pairs = dict(base, pairs=pairs)
        if least_after(with_pairs, []) <= max_bytes:
            out, notes = with_pairs, []
        elif least_after(base, [no_pairs]) <= max_bytes:
            out, notes = dict(base), [no_pairs]
        else:
            # the smallest answer: the statements and the notes of what is left out
            least = min(least_after(with_pairs, []), least_after(base, [no_pairs]))
            raise _refuse('summary_too_small', None, least, max_bytes, least)
        cand = dict(out, arms=arms)
        if size(cand, notes + cut) <= max_bytes:
            out = cand
            for k, (s, row) in enumerate(ranked):
                arms[s]['top'].append(row)
                if size(out, notes + (cut if k + 1 < len(ranked) else [])) > max_bytes:
                    arms[s]['top'].pop()
                    notes += cut
                    break
        else:
            notes.append(no_arms)
        if notes:
            out['notes'] = notes
        return out

    def evidence(self, arm=None, *, budget=None):
        """Each part's own evidence block under its tag, the absent records, the expected
        pairs not included, each part's identity and pair source, the weakest outcome over
        the parts, exact (every part exact and its identity verified) and the five
        statements (without an arm: each arm's coverage sentence, then the provenance over
        every arm, the turn-offs, the labels and the scope). budget=: the blocks, and the
        coverage and statements as coverage() and statements() charge them; a stitched
        part's own evidence() runs under it as the ambient budget."""
        side = None if arm is None else self._side(arm)
        return _B.run_scoped(budget, 'union_evidence', (), lambda b: self._evidence(side, b))

    def _evidence(self, side, b):
        if b is not None:
            # per part its block (and its identity entry), per absent record its entry
            b.charge(W_ROW * (len(self.parts) + len(self.absent)) + W_ELEM * sum(
                len(_root(p.graphlet).limitations) for p in self.parts),
                _M_BLOCK * len(self.parts) + _M_COV_PAIR * len(self.absent))
        ev = _combined_evidence(self.parts, side, b)
        if side is not None or not self.parts:
            sentences = self._statements(side, None, 'supported', b).sentences()
        else:
            # every arm: the coverage sentence of each, then the provenance over them all
            sides = self.arms()
            sentences = [_core_sentence(self._coverage(s, None, b)) for s in sides] + [
                'walks: ' + ', '.join(self._provenance_sentence(s, b)[len('walks: '):]
                                      for s in sides),
                TURN_OFFS, self._labels_sentence(b), SCOPE]
        return {
            'parts': ev['parts'],
            'outcome': ev['outcome'],
            'exact': ev['exact'],
            'identity': {p.tag: dict(p.identity, pair_stated_by=p.pair_stated_by)
                         for p in self.parts},
            'absent': [{'pair': a.tag, 'category': a.category, 'seed_check': a.seed_check}
                       for a in self.absent],
            'not_included': list(self._not_included),
            'statements': sentences,
        }

    def provenance(self):
        """One entry per part ({tag, pair, pair_stated_by, identity, seed_id, digest,
        request}: the replay request as GraphletStore.request_of() builds it), per absent
        record, per earlier failed attempt, and the expected pairs."""
        from .store import GraphletStore
        out = []
        for p in self.parts:
            g = _root(p.graphlet)
            pair = dict(zip(_KEY_FIELDS, p.key))
            entry = {'tag': p.tag, 'kind': 'part', 'pair': pair,
                     'pair_stated_by': p.pair_stated_by, 'identity': p.identity,
                     'seed_id': ((g.seed_summary or {}).get('seed') or {}).get('seed_id'),
                     'digest': p.digest,
                     'request': GraphletStore.request_of(
                         g, graph=pair['graph'], graph_path=pair['graph_path'],
                         annotation_path=pair['annotation_path'])
                     if p.pair_stated_by == 'caller' and pair['graph'] else
                     GraphletStore.request_of(g)}
            if isinstance(p.graphlet, StitchedGraphlet) and p.graphlet.junctions:
                # a stitched part: its root's request above, and each child's entry
                # (StitchedGraphlet.provenance()): no single request reproduces it
                entry['stitched'] = p.graphlet.provenance()[1:]
            out.append(entry)
        for a in self.absent:
            out.append({'tag': a.tag, 'kind': 'absent', 'pair': dict(zip(_KEY_FIELDS, a.key)),
                        'pair_stated_by': a.pair_stated_by, 'category': a.category,
                        'detail': a.detail, 'seed_check': a.seed_check})
        for f in self._failures:
            out.append({'tag': f['tag'], 'kind': 'earlier_failure', 'error': f['error']})
        if self.expected is not None:
            out.append({'kind': 'expected_pairs', 'tags': list(self.expected)})
        return out


class _SeedStandIn:
    __slots__ = ('sequence', 'length_bp')

    def __init__(self, seq):
        self.sequence = seq
        self.length_bp = len(seq)


# ------------------------------------------------------------------ pairs, tags, dedupe

def _seqs_of(classes, so, tags):
    """id(class) -> its seqs_order value ('unknown' where seqs_order= names it not), by the
    one matching rule; bad_seqs_order where two keys give one pair two values."""
    seqs = {}
    for c in classes:
        vals = {v for d, v, _ in so if _matches(d, c['desc'])}
        if len(vals) > 1:
            k = next(k for d, v, k in so if _matches(d, c['desc']))
            raise _refuse('bad_seqs_order', None, 'seqs_order[%r] (two values for pair %s)'
                          % (k, tags[id(c)]))
        seqs[id(c)] = vals.pop() if vals else 'unknown'
    return seqs


def _classes(descs):
    """The pair classes of |descs| (match is transitive once ambiguous_pair passed):
    [{'members': [desc...], 'desc': the fullest description}]."""
    classes = []
    for d in descs:
        hit = None
        for c in classes:
            if _matches(d, c['members'][0]):
                hit = c
                break
        if hit is None:
            classes.append({'members': [d]})
        else:
            hit['members'].append(d)
    for c in classes:
        c['desc'] = min(c['members'], key=lambda d: (-len(d), _sort_key(_key(d))))
    return classes


def _tags(classes):
    """id(class) -> its display tag, unique (Appendix A.4)."""
    releases = {c['desc'].get('release') for c in classes if c['desc'].get('release')}
    base = {}
    for c in classes:
        d = c['desc']
        name = d.get('graph') or d.get('index_ns')
        if len(releases) > 1 and d.get('release'):
            name = '%s/%s' % (d['release'], name)
        base[id(c)] = name
    groups = collections.defaultdict(list)
    for c in classes:
        groups[base[id(c)]].append(c)
    out = {}
    for name, cs in groups.items():
        if len(cs) == 1:
            out[id(cs[0])] = name
            continue
        hexes = {id(c): hashlib.sha256(C.canonical_json(list(_key(c['desc']))).encode())
                 .hexdigest() for c in cs}
        n = 6
        while len({h[:n] for h in hexes.values()}) < len(cs):
            n += 1
        for c in cs:
            out[id(c)] = '%s~%s' % (name, hexes[id(c)][:n])
    return out


def _absent_identity(e):
    r = dict(e.result)
    if isinstance(r.get('seed'), dict):
        r['seed'] = {k: v for k, v in r['seed'].items() if k != 'seed_id'}
    return C.canonical_json(r)


def _dedupe(elems, cls_of, tags):
    """duplicate_pair, and the no-ops of Appendix A.4 -> {'by_class': id(class) -> the
    element kept (an answer, else the first failure in canonical order), 'failures':
    id(class) -> [failed elements], 'order': kept elements in the order given}."""
    seen = {}
    failures = collections.defaultdict(list)
    for e in elems:
        c = cls_of[id(e.desc)]
        if e.kind == 'failed':
            failures[id(c)].append(e)
            continue
        prev = seen.get(id(c))
        if prev is None:
            seen[id(c)] = e
            continue
        tag = tags[id(c)]
        if prev.kind == 'part' and e.kind == 'part':
            # a stitched part: its digest is its body and its retrieval identity, so two
            # equal ones are a no-op, and one beside a plain graphlet of its pair is
            # duplicate_pair (stitch(g) and g included)
            same_body = _body_digest(prev.graphlet) == _body_digest(e.graphlet)
            ra, rb = _retrieval_of(prev.graphlet), _retrieval_of(e.graphlet)
            if same_body and ra == rb:
                # a no-op: the envelope-stated form wins; of two, the smaller envelope
                if _prefer(e, prev):
                    seen[id(c)] = e
                continue
            if same_body:
                raise _refuse('duplicate_pair', 'duplicate_pair_request', prev.pos, e.pos, tag,
                              C.first_difference(ra, rb))
        elif prev.kind == 'absent' and e.kind == 'absent' \
                and _absent_identity(prev) == _absent_identity(e):
            # a no-op: the envelope-stated form wins, then the stronger seed check, so the
            # form kept does not depend on the order given
            if _absent_rank(e) < _absent_rank(prev):
                seen[id(c)] = e
            continue
        if not prev.desc.get('index_fp') and not e.desc.get('index_fp') \
                and all(prev.desc.get(k) == e.desc.get(k)
                        for k in ('graph', 'graph_path', 'annotation_path', 'release',
                                  'index_ns', 'index_meta_fp')):
            raise _refuse('duplicate_pair', 'duplicate_pair_unknown', prev.pos, e.pos, tag)
        raise _refuse('duplicate_pair', None, prev.pos, e.pos, tag)
    by_class = dict(seen)
    for cid, fs in failures.items():
        if cid not in by_class:
            by_class[cid] = min(fs, key=lambda f: f.error.encode('utf-8'))
    order = [e for e in elems if by_class.get(id(cls_of[id(e.desc)])) is e]
    return {'by_class': by_class, 'failures': failures, 'order': order}


def _seed_check(e):
    return 'sent' if e.sent is not None else ('stated' if e.sequence is not None
                                              else 'unchecked')


def _absent_rank(e):
    """Which of two forms of one absent record is kept: the envelope-stated one, then the
    stronger seed check, then the whole record as canonical JSON (its seed_id included)
    and the stated pair: a total order, so the form kept never depends on the order
    given."""
    return (e.by != 'envelope', ('sent', 'stated', 'unchecked').index(_seed_check(e)),
            C.canonical_json(e.result), C.canonical_json(e.stated or {}), e.error or '')


def _part_rank(e):
    """Which of two forms of one retrieval is kept: the envelope-stated one, then the
    envelope's canonical JSON, then the seed summary's (the caller's seed_id) and the
    stated pair: a total order over what the union reports of a part."""
    g = e.graphlet
    if isinstance(g, StitchedGraphlet):
        # every part's envelope and seed summary (review minor 3): two forms of one
        # stitched graphlet may differ in a child's seed_id
        gs = _sub_graphlets(g)
        return (e.by != 'envelope', C.canonical_json([x.envelope for x in gs]),
                C.canonical_json([x.seed_summary or {} for x in gs]),
                C.canonical_json(e.stated or {}))
    return (e.by != 'envelope', C.canonical_json(g.envelope),
            C.canonical_json(g.seed_summary or {}), C.canonical_json(e.stated or {}))


def _prefer(e, prev):
    return _part_rank(e) < _part_rank(prev)


def _body_digest(g):
    return g.digest if isinstance(g, StitchedGraphlet) else C.part_digest(g)


def _retrieval_of(g):
    return {'stitched': g.digest} if isinstance(g, StitchedGraphlet) else \
        C.retrieval_identity(g)


def _label_index_price(b, g):
    """A lookup by name or ref in part |g| at its cold price: the name/ref index of each of
    its graphlets (ops.label() scans or builds it; the index's price either way)."""
    for x in _sub_graphlets(g):
        _ops._uses_g(b, x, 'label_index')


def _label_in(g, selector):
    """A label of a part by {'ref'} or {'name'}: a stitched part's over its parts' labels
    (StitchedGraphlet.label(); ids in its own numbering)."""
    return g.label(selector) if isinstance(g, StitchedGraphlet) else _ops.label(g, selector)


def _arm_junctions(g, side):
    """The junctions of a stitched part on |side| (every arm: None); () for a plain part."""
    if not isinstance(g, StitchedGraphlet):
        return ()
    return tuple(j for j in g.junctions if side is None or j.arm == side)


def _junction_note(g, side, below):
    """The coverage note of a stitched part's junctions on |side|, or None: its category is
    its root's (stitching certifies nothing past a junction). |below|: part key -> the
    walk-scope cause codes on the stitched walks through that junction."""
    js = _arm_junctions(g, side)
    if not js:
        return None
    depths = sorted({j.at_bp for j in js})
    verdicts = sorted({j.verdict for j in js})
    codes = sorted({c for j in js for c in below.get(j.part, ())}, key=_CAUSE_RANK.get)
    return 'stitched (experimental): junction%s at %s bp, %s%s; nothing past a junction is ' \
        'certified' % ('s' if len(depths) > 1 else '', ', '.join(map(str, depths)),
                       ' or '.join(verdicts),
                       '; walk causes below: %s' % ', '.join(codes) if codes else '')


def _summary_part(p, below):
    out = {'part': True, 'mode': p.graphlet.mode, 'k': p.graphlet.k,
           'stated_by': p.pair_stated_by}
    js = _arm_junctions(p.graphlet, None)
    if js:
        # as StitchedGraphlet.summary() states them: the junction's causes, and the
        # walk-scope causes on the stitched walks through it
        out['junctions'] = [{'arm': j.arm, 'at_bp': j.at_bp, 'verdict': j.verdict,
                             'causes': [c[0] for c in j.causes],
                             'walk_causes': list(below.get(j.part, ()))} for j in js]
    return out


def _part_block(p, side, b=None):
    """A part's evidence block: its own (ops.evidence_block()); a stitched part's
    StitchedGraphlet.evidence() where it has a junction on |side| (every arm: None), run
    under the union call's budget |b| as the ambient one and stated as the union's call (a
    stop there names the union's operation and levers: Opus fix review N5), its root's
    otherwise (no junction is read there)."""
    g = p.graphlet
    if _arm_junctions(g, side):
        with _B.part_of(b):
            return g.evidence(side)
    g = _root(g)
    return _ops.evidence_block(g, side if side is None or side in g.arms else None)


def _combined_evidence(parts, side, b=None):
    """composite.combined_evidence() over the parts, a stitched part's block by
    _part_block() (|b|: the union call's budget)."""
    if not any(_arm_junctions(p.graphlet, side) for p in parts):
        return C.combined_evidence([(p.tag, _root(p.graphlet), _part_verified(p, side))
                                    for p in parts], side)
    blocks = {}
    exact = True
    for p in parts:
        blk = _part_block(p, side, b)
        blocks[p.tag] = blk
        exact = exact and bool(blk.get('exact')) and bool(_part_verified(p, side))
    return {'parts': blocks,
            'outcome': C.weakest_outcome([x.get('outcome') for x in blocks.values()]),
            'exact': exact and bool(blocks)}


def _part_verified(p, side):
    """Whether part |p|'s identity is verified for an answer on |side| (every arm: None):
    a stitched part's own identity where it has a junction there, its root's otherwise (no
    junction is read, so the root answers alone; review minor 4)."""
    if isinstance(p.graphlet, StitchedGraphlet) and not _arm_junctions(p.graphlet, side):
        return _identity_statement(p.graphlet.root)['verified']
    return p.identity.get('verified')


def _identity_statement(g):
    if g.index_fp:
        return {'verified': True, 'index_fp': g.index_fp}
    return {'verified': False, 'note': 'no index manifest digest (index_fp): the index is '
                                       'not verified'}


def _absent_detail(e):
    r = e.result
    out = {'outcome': r.get('outcome')}
    for k in ('not_in_graph', 'unknown_labels', 'limitations'):
        if r.get(k) is not None:
            out[k] = r[k]
    err = r.get('error') or e.error
    if err:
        out['error'] = err
    meta = _seed_meta(r)
    for k in ('seed_id', 'length_bp', 'labels_from_seed'):
        if k in meta:
            out.setdefault('seed', {})[k] = meta[k]
    return out


# ------------------------------------------------------------------ coverage helpers

def _cut60(text):
    # 61 characters encode to more than 60 bytes, so a longer text is never encoded whole
    text = text or ''
    b = text[:61].encode('utf-8')
    if len(b) <= 60:
        return text
    return b[:57].decode('utf-8', 'ignore') + '...'


def _part_category(p, side, depth, below=None):
    """A part's coverage entry; a stitched part's is its root's, with a note naming its
    junctions on |side| (|below|: _junction_note())."""
    got = _root_category(_root(p.graphlet), side, depth)
    note = _junction_note(p.graphlet, side, below or {}) if side is not None else None
    if note:
        got['note'] = note if not got.get('note') else '%s; %s' % (got['note'], note)
    return got


def _root_category(g, side, depth):
    if side is None or side not in g.arms:
        return {'category': 'arm_not_walked', 'certified_to_bp': 0,
                'note': None if side is None else 'direction %s' % (
                    ((g.envelope or {}).get('strategy') or {}).get('direction'))}
    a = g.arms[side]
    ctb = a.complete_to_bp
    if ctb >= depth:
        quals = []
        for lim in g.limitations:
            if lim.arm not in (None, side):
                continue
            if lim.kind == 'seed_labels':
                quals.append('complete for %s of %s seed labels' % (lim.limit, lim.observed))
            elif lim.kind == 'scope' and isinstance(lim.observed, (int, float)) \
                    and not isinstance(lim.observed, bool) and lim.observed > 0:
                quals.append('a merge united histories: complete for the united-history '
                             'rule, not per path')
        if quals:
            return {'category': 'complete_qualified', 'certified_to_bp': ctb,
                    'note': '; '.join(dict.fromkeys(quals))}
        return {'category': 'complete', 'certified_to_bp': ctb, 'note': None}
    reason = None
    if a.status != 'complete' and a.cap_trigger is not None:
        reason = REASON.get(a.cap_trigger.reason, a.cap_trigger.reason)
    elif a.status != 'complete':
        reason = a.status
    else:
        reason = 'radius'
    return {'category': 'incomplete', 'certified_to_bp': ctb,
            'note': 'certified to %d bp: %s' % (ctb, reason)}


def _absent_entry(a):
    d = a.detail
    note = None
    if a.category == 'seed_not_in_graph':
        nig = d.get('not_in_graph') or {}
        if 'kmers' in nig:
            note = '%s of %s k-mers' % (nig.get('kmers_present'), nig.get('kmers'))
    elif a.category == 'not_walked':
        der = next((l for l in (d.get('limitations') or [])
                    if isinstance(l, dict) and l.get('kind') == 'derivation'), None)
        if der is not None:
            cause = der.get('cause')
            note = 'seed in graph; derivation: %s' % _DERIVATION_CAUSES.get(cause, cause)
        elif d.get('error'):
            note = _cut60(d['error'])
    elif a.category == 'failed':
        note = _cut60(d.get('error'))
    return {'category': a.category, 'certified_to_bp': None, 'note': note,
            'seed_check': a.seed_check}


def _core_sentence(cov):
    """The coverage sentence without its trailing turn-off and scope sentences."""
    s = cov.sentence
    tail = ' ' + TURN_OFFS + ' ' + SCOPE
    return s[:-len(tail)] if s.endswith(tail) else s


# ------------------------------------------------------------------ the trie

class _USeg:
    __slots__ = ('id', 'from_bp', 'to_bp', 'parent', 'children', 'cover', 'pieces', 'end',
                 'first', 'n')

    def __init__(self):
        self.id = -1
        self.children = []
        self.parent = None


class _Trie:
    """One arm's radix trie of every part's spellings (§2.4), built by one depth-first pass
    over each part's first-parent tree. Raw nodes sit at every part segment boundary;
    union segments compress the raw edges between nodes where nothing happens."""

    def __init__(self, union, side, b):
        self.side = side
        # the sources (PairUnion._sources()): part i for a union without a stitched part
        self.parts = union._src
        order = union._order
        self.kids = [{}]                # node -> {char: edge}
        self.depth = [0]
        self.e_from, self.e_to, self.e_len, self.e_ref, self.e_cover = [], [], [], [], []
        self.bases_read = 0             # bases the insertion compared or indexed (a count)
        self.walk_end = {}              # node -> [(pi, leaf seg, path id)]
        self.stub_end = {}              # node -> [(pi, stub seg, merge seg)]
        self.merge_at = {}              # node -> [(pi, merge seg)]
        self.seg_nodes = []             # pi -> {seg: (start node, end node)}
        self.zero_at = {}               # node -> [(pi, zero-length seg)]
        self.below_walks, self.below_stubs = [], []
        # a stitched part (§2.3): a child's tree goes in where the walk it continues ends,
        # (parent source, that node) -> the child's source; the continued leaf ends no walk
        self.graft = {}
        self.graft_nodes = set()
        self.continued = {}             # (parent source, leaf seg) -> the child's source
        for si, src in enumerate(self.parts):
            if src.parent is not None and side in src.graphlet.arms:
                self.continued[(src.parent, src.leaf)] = si
        fp_kids_of = []
        for pi, part in enumerate(self.parts):
            g = part.graphlet
            a = g.arms.get(side)
            self.seg_nodes.append({})
            self.below_walks.append({})
            self.below_stubs.append({})
            fp_kids_of.append(None)
            if a is None:
                continue
            start = 0
            if part.parent is not None:
                start = self.seg_nodes[part.parent][part.leaf][1]
                self.graft[(part.parent, start)] = pi
                self.graft_nodes.add(start)
            fp_kids_of[pi] = self._insert_part(pi, part, g, a, b, start)
        # the walks and stubs below each segment (first-parent subtree), children before
        # their parents: a continued leaf counts its child's
        for pi in range(len(self.parts) - 1, -1, -1):
            if fp_kids_of[pi] is not None:
                self._below(pi, self.parts[pi].graphlet.arms[side], fp_kids_of[pi])
        if b is not None:
            # the union segments at their bound: one per raw edge, a piece per cover entry
            b.charge(W_ELEM * len(self.e_from),
                     _M_USEG * len(self.e_from) + _M_PIECE * sum(map(len, self.e_cover)))
        self._compress(order)

    # ---------------------------------------------------------- insertion

    def _bases(self, e):
        pi, seg, off = self.e_ref[e]
        s = self.parts[pi].graphlet.arms[self.side].segments[seg]
        return s.walk[off:off + self.e_len[e]]

    def _ref_walk(self, e):
        """(the walk edge |e|'s bases are read from, their offset in it): an edge's bases
        are compared and indexed there, never sliced whole (GPT review M5)."""
        pi, seg, off = self.e_ref[e]
        return self.parts[pi].graphlet.arms[self.side].segments[seg].walk, off

    def _new_node(self, depth):
        self.kids.append({})
        self.depth.append(depth)
        return len(self.kids) - 1

    def _new_edge(self, f, t, n, ref, cover):
        self.e_from.append(f)
        self.e_to.append(t)
        self.e_len.append(n)
        self.e_ref.append(ref)
        self.e_cover.append(cover)
        return len(self.e_from) - 1

    def _split(self, e, m):
        """Split edge |e| after |m| bases -> the new middle node."""
        f, t, n = self.e_from[e], self.e_to[e], self.e_len[e]
        mid = self._new_node(self.depth[f] + m)
        pi, seg, off = self.e_ref[e]
        lower_cover = {k: (s, o + m) for k, (s, o) in self.e_cover[e].items()}
        low = self._new_edge(mid, t, n - m, (pi, seg, off + m), lower_cover)
        src, o = self._ref_walk(low)
        first = src[o]
        self.bases_read += 1
        self.kids[mid][first] = low
        self.e_to[e] = mid
        self.e_len[e] = m
        return mid

    def _insert_seg(self, pi, s, start):
        """Insert part |pi|'s segment |s| at trie node |start| -> its end node."""
        w = s.walk
        n = len(w)
        node = start
        i = 0
        while i < n:
            ch = w[i]
            e = self.kids[node].get(ch)
            if e is None:
                t = self._new_node(self.depth[node] + n - i)
                e = self._new_edge(node, t, n - i, (pi, s.id, i), {pi: (s.id, i)})
                self.kids[node][ch] = e
                return t
            # the edge's bases against w from i, in place in the walk they are read from:
            # a full match costs the bases it consumes, a mismatch at most the segment's
            # remaining bases (and it ends the insertion), so a segment costs O(its bases)
            # however long the edges it meets (GPT review M5)
            src, o = self._ref_walk(e)
            el = self.e_len[e]
            lim = min(el, n - i)
            self.bases_read += lim
            if src.startswith(w[i:i + lim], o):
                m = lim
            else:
                m = 1
                while m < lim and src[o + m] == w[i + m]:
                    m += 1
            if m < el:
                self._split(e, m)
            self.e_cover[e][pi] = (s.id, i)
            node = self.e_to[e]
            i += m
        return node

    def _insert_part(self, pi, part, g, a, b, start=0):
        """Insert source |pi|'s first-parent tree at trie node |start| (the root, or where a
        stitched child continues its walk) -> the first-parent children per segment."""
        segs = a.segments
        if b is not None:
            # one pass over the part's first-parent tree: an element per segment, its bases
            # compared, and at most two trie nodes, an edge and a cover entry per segment.
            # Memory at its bound: the first part's segments each add one node and one edge
            # (siblings start with distinct bases); a later part's segment at most two of
            # each (a split) and a copy of the split edge's cover (one entry per part before
            # it)
            n = sum(s.length_bp for s in segs)
            per = (_M_NODE + _M_EDGE + _M_COVER) if pi == 0 else \
                2 * (_M_NODE + _M_EDGE) + _M_COVER * (pi + 1)
            b.charge(W_ELEM * 4 * len(segs) + (n >> 4),
                     3 * DICT + len(segs) * (per + _M_SEG))
            derive.uses(b, g, a, 'leaves', 'paths', 'path_of_leaf')
        fp_kids = [[] for _ in segs]
        roots = []
        for s in segs:
            if s.parents:
                fp_kids[s.parents[0]].append(s.id)
            else:
                roots.append(s.id)
        pol = derive.path_of_leaf(a)
        nodes = self.seg_nodes[pi]
        stack = [(r, start) for r in reversed(roots)]
        first = {}                      # (trie node, first base) -> the segment starting so
        while stack:
            sid, start = stack.pop()
            s = segs[sid]
            if s.length_bp:
                # a prefix names one path per graph: two segments of one part that start at
                # one node with one base would share the trie's edge, one overwriting the
                # other's cover
                other = first.setdefault((start, s.walk[0]), sid)
                if other != sid:
                    raise GraphletFormatError(
                        0, 'part %s: segments %d and %d of the %s arm both start with %r at %d '
                        'bp (sibling segments start with distinct bases)'
                        % (part.tag, other, sid, self.side, s.walk[0], s.from_bp))
            end = self._insert_seg(pi, s, start)
            nodes[sid] = (start, end)
            if s.length_bp == 0:
                self.zero_at.setdefault(end, []).append((pi, sid))
            if s.is_merge:
                self.merge_at.setdefault(start, []).append((pi, sid))
            kids = fp_kids[sid]
            if s.leaf is not None:
                if (pi, sid) not in self.continued:
                    self.walk_end.setdefault(end, []).append((pi, sid, pol[sid]))
            elif not kids:
                merges = [c for c in s.children if segs[c].parents[0] != sid]
                if len(merges) != 1:
                    raise GraphletFormatError(
                        0, 'part %s: segment %d of the %s arm ends a stub and is a non-first '
                        'parent of %d merges (one is expected)'
                        % (part.tag, sid, self.side, len(merges)))
                self.stub_end.setdefault(end, []).append((pi, sid, merges[0]))
            elif any(segs[c].parents[0] != sid for c in s.children):
                raise GraphletFormatError(
                    0, 'part %s: segment %d of the %s arm ends in both a split and a merge'
                    % (part.tag, sid, self.side))
            for c in reversed(kids):
                stack.append((c, end))
        return fp_kids

    def _below(self, pi, a, fp_kids):
        """The walks and stubs below each segment of source |pi| (first-parent subtree); a
        continued leaf ends no walk and counts those of the child continuing it."""
        bw, bs = self.below_walks[pi], self.below_stubs[pi]
        for s in reversed(a.segments):
            c = self.continued.get((pi, s.id)) if s.leaf is not None else None
            if c is not None:
                ca = self.parts[c].graphlet.arms[self.side]
                w = sum(self.below_walks[c][r.id] for r in ca.segments if not r.parents)
                st = sum(self.below_stubs[c][r.id] for r in ca.segments if not r.parents)
            else:
                w = 1 if s.leaf is not None else 0
                st = 1 if (s.leaf is None and not fp_kids[s.id]) else 0
            for c in fp_kids[s.id]:
                w += bw[c]
                st += bs[c]
            bw[s.id], bs[s.id] = w, st

    # ---------------------------------------------------------- compression

    def _compress(self, order):
        """Union segments: maximal chains of raw edges whose inner nodes have one child,
        no marker and no merge start, and the same pairs."""
        def rank(ch):
            return (order.get(ch, len(order)), ch)

        self.useg = []
        self.node_in = {}               # node -> union segment ending there
        self.node_out = {}              # node -> {char: union segment starting there}
        self.markers = []               # terminal nodes, preorder (walk ids)
        stack = [(0, None)]
        while stack:
            node, parent = stack.pop()
            if node in self.walk_end or node in self.stub_end:
                self.markers.append(node)
            kids = sorted(self.kids[node].items(), key=lambda kv: rank(kv[0]))
            outs = {}
            todo = []
            for ch, e in kids:
                u = _USeg()
                u.parent = parent
                u.from_bp = self.depth[node]
                u.first = ch
                u.n = 1 if parent is None else self.useg[parent].n + 1
                u.cover = dict.fromkeys(self.e_cover[e])
                edges = [e]
                cur = self.e_to[e]
                while True:
                    k = self.kids[cur]
                    if len(k) != 1 or cur in self.walk_end or cur in self.stub_end \
                            or cur in self.merge_at or cur in self.zero_at \
                            or cur in self.graft_nodes:
                        break
                    (e2,) = k.values()
                    if set(self.e_cover[e2]) != set(u.cover):
                        break
                    edges.append(e2)
                    cur = self.e_to[e2]
                u.end = cur
                u.to_bp = self.depth[cur]
                u.pieces = {}
                for pi in u.cover:
                    ps = []
                    # positions outward from the seed: a stitched child's at its offset
                    base = self.parts[pi].offset
                    for x in edges:
                        seg, off = self.e_cover[x][pi]
                        s = self.parts[pi].graphlet.arms[self.side].segments[seg]
                        f = base + s.from_bp + off
                        t = f + self.e_len[x]
                        if ps and ps[-1][0] == seg and ps[-1][2] == f:
                            ps[-1] = (seg, ps[-1][1], t)
                        else:
                            ps.append((seg, f, t))
                    u.pieces[pi] = ps
                idx = len(self.useg)
                self.useg.append(u)
                if parent is not None:
                    self.useg[parent].children.append(idx)
                self.node_in[cur] = idx
                outs[ch] = idx
                todo.append((cur, idx))
            self.node_out[node] = outs
            for item in reversed(todo):
                stack.append(item)

    def number(self, n):
        """Number the union segments in preorder from |n| -> the next free id."""
        self.order = []
        stack = [u for u in reversed(self.node_out.get(0, {}).values())]
        while stack:
            u = stack.pop()
            self.useg[u].id = n
            self.order.append(u)
            n += 1
            stack.extend(reversed(self.useg[u].children))
        return n

    def path_usegs(self, node):
        """The union segments from the root to |node|, root first: built from the parent
        pointers on each call (nothing kept per node)."""
        out = []
        u = self.node_in.get(node)
        while u is not None:
            out.append(u)
            u = self.useg[u].parent
        out.reverse()
        return out

    def nsegs(self, node):
        """How many union segments lie between the root and |node|."""
        u = self.node_in.get(node)
        return 0 if u is None else self.useg[u].n

    def stop_of(self, pi, node, floor=None):
        """(stop, below): the deepest node on the path to |node| whose spellings part |pi|
        follows (covers are prefix-closed), and the union segments past it to |node|, root
        first. Costs the segments past the stop only. |floor|: a stitched child's graft
        node, where it starts (its stop when it spells nothing of the path below it)."""
        below = []
        u = self.node_in.get(node)
        while u is not None and pi not in self.useg[u].cover \
                and (floor is None or self.useg[u].end != floor):
            below.append(u)
            u = self.useg[u].parent
        below.reverse()
        return (0 if u is None else self.useg[u].end), below

    def last_seg(self, pi, a, stop):
        """The segment of part |pi| (arm |a|) that holds the stop's depth: the last piece
        of the union segment ending at |stop|, or the arm's root at the trie root."""
        u = self.node_in.get(stop)
        if u is None or pi not in self.useg[u].cover:
            # the trie root, or a stitched child's graft node: the arm's root segment
            return next((s.id for s in a.segments if not s.parents), None)
        return self.useg[u].pieces[pi][-1][0]

    def node_name(self, node, merge_seg_len, first=None):
        """The trie node of a merge: 's<id>' (the start of the union segment holding the
        merge segment's first base) or 's<id>$' (the end of the union segment ending at the
        node, for a merge segment without bases); 'seed' at the root."""
        if merge_seg_len > 0 and first is not None:
            u = self.node_out.get(node, {}).get(first)
            if u is not None:
                return 's%d' % self.useg[u].id
        u = self.node_in.get(node)
        return 'seed' if u is None else 's%d$' % self.useg[u].id



# ------------------------------------------------------------------ spans

def _leaf_reasons(g, a, seg):
    """(end reasons, code) of a part leaf. In constrain mode a structural dead end has no
    path reason: the walker ends every label there with a bare dead end (D), so a leaf whose
    every label ended with D reads code 'D' (a hairpin dead end, Dh, does not count: a
    hairpin successor exists, and successors no label carries may too)."""
    s = a.segments[seg]
    code = s.leaf.path_reason if s.leaf is not None else None
    labels = {}
    if g.mode == 'constrain' and s.leaf is not None:
        ends = derive.end_labels(a, seg)
        for e in ends:
            name = a.runs[e.run].reason
            labels[name] = labels.get(name, 0) + 1
        if code is None and ends and all(a.runs[e.run].end == 'D' for e in ends):
            code = 'D'
    path_reason = s.leaf.path_reason if s.leaf is not None else None
    return {'path_reason': REASON[path_reason] if path_reason else None,
            'labels': labels}, code


def _leaf_cert(a, seg):
    """(complete, certified_to_bp) of the part walk ending at leaf |seg|: ops._walk_end()'s
    own (no resource stop at the leaf, within complete_to_bp)."""
    p = derive.paths(a)[derive.path_of_leaf(a)[seg]]
    _, complete, beyond = _ops._walk_end(a, p)
    return complete, p.length_bp - beyond


def _events_at(t, pi, a, stop, last):
    """The part segments whose node at the stop's depth is the stop: |last|, the one ending
    there or holding it, and the zero-length ones there."""
    out = [a.segments[last]] if last is not None else []
    for pi2, z in t.zero_at.get(stop, ()):
        if pi2 == pi:
            out.append(a.segments[z])
    return out


def _absence(g, a, segs_here, d, c, code=None):
    """§2.4's absence where a pair stops at depth |d| and the union walk goes on with |c|."""
    if c is not None:
        for s in segs_here:
            for ev in s.events:
                if ev.at_bp != d or ev.char != c:
                    continue
                if ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed):
                    return 'walk_rule'
    if code is not None and code[0] == 'D':
        return 'graph'
    ctb = a.complete_to_bp
    ev_ctb = a.evidence_complete_to_bp
    if g.mode == 'annotate' and a.scope == 'per_path' and d < ctb \
            and (ev_ctb is None or d < ev_ctb) and c is not None:
        return 'graph'
    if d >= ctb:
        return 'uncertified'
    return 'under_strategy'


def _case(t, pi, stop, length):
    """How part |pi| stops on a union walk of |length| bp whose spellings it follows to
    trie node |stop|: ('whole' or 'ends', its walk end there), ('extends', None), ('merges',
    its stub end there) or ('diverges', None). One decision for the walk figures and the
    spans."""
    d = t.depth[stop]
    wh = next((x for x in t.walk_end.get(stop, ()) if x[0] == pi), None)
    if wh is not None:
        return ('whole' if d == length else 'ends'), wh
    if d == length and any(pi in t.useg[u].cover for u in t.node_out.get(stop, {}).values()):
        return 'extends', None
    sh = next((x for x in t.stub_end.get(stop, ()) if x[0] == pi), None)
    if sh is not None:
        return 'merges', sh
    return 'diverges', None


def _labels_to(un, part, a, seg, end, with_labels, at_bp, b=None):
    """_labels_for() along the first-parent chain to |seg| (the chain built only when the
    query computes labels). |part|: the source (_labels_src())."""
    if not with_labels:
        return None, False
    if seg is None:
        return (), False
    return _labels_src(un, part, a, lambda: _fp_chain(a, seg), end, at_bp, b)


def _labels_src(un, src, a, chain, end, at_bp, b=None):
    """_labels_for() of source |src| at outward depth |end| (or |at_bp| when smaller) along
    chain() -> (SpanLabels, cut), positions outward from the seed. A stitched child's labels
    are its own from the junction on (from_bp at least its offset: not joined across a
    junction, D.1); at a depth no base of the child reaches, its parent's along the walk it
    continues."""
    e = end if at_bp is None else min(end, at_bp)
    if src.key and e <= src.offset:
        p = un._src[src.parent]
        pa = p.graphlet.arms[a.side]
        return _labels_src(un, p, pa, lambda: _fp_chain(pa, src.leaf), e, None, b)
    if not src.offset:
        return _labels_for(un, src, a, chain(), end, True, at_bp, b)
    got, cut = _labels_for(un, src, a, chain(), e - src.offset, True, None, b)
    return tuple(SpanLabel(x.label, x.from_bp + src.offset) for x in got), cut


def _absence_at(src, g, a, segs_here, d, c, code=None):
    """_absence() at local depth |d| of source |src|; past a junction (a stitched child)
    only a dead end is the graph's, every other turn-off 'uncertified' (§2.3): the child's
    own seed re-entry and edge history may block where the reference goes on."""
    if src.key:
        # a base the child blocked or a hairpin it did not follow is there in the graph
        # (review M1): uncertified, never graph
        if c is not None:
            for s in segs_here:
                for ev in s.events:
                    if ev.at_bp == d and ev.char == c and (
                            ev.type == 'blocked' or (ev.type == 'hairpin' and not ev.followed)):
                        return 'uncertified'
        return 'graph' if code is not None and code[0] == 'D' else 'uncertified'
    return _absence(g, a, segs_here, d, c, code)


def _past_junction(un, src, span):
    """A span read in a stitched child (§2.3): its walk named by its stitched key, never
    complete, certified at most to the first junction (the stitched walk's certified_to_bp:
    no junction is exact in this phase), and the junctions it reads, root first."""
    if not src.key:
        return span
    js = []
    s = src
    # the walk read in each child: the span's own in the last (whole, ends, extends; to its
    # end for extends), the continued walk in the others
    pid = span.walk if span.how in ('whole', 'ends', 'extends') else None
    limit = span.to_bp - src.offset if span.how == 'extends' else None
    while s.key:
        j = s.junction
        own = [c[0] for c in s.wcauses.get(pid, ()) if limit is None or c[2] <= limit] \
            if pid is not None else []
        codes = tuple(c[0] for c in j.causes)
        js.append(SpanJunction(j.at_bp, j.verdict, j.part, codes + tuple(
            c for c in sorted(set(own), key=_CAUSE_RANK.get) if c not in codes)))
        first = s
        pid, limit = s.key[-1][1], None
        s = un._src[s.parent]
    js.reverse()
    # the root's certificate on the walk the first child continues: its certified_to_bp
    side, pid = first.key[0]
    ra = s.graphlet.arms[side]
    cert = min(derive.paths(ra)[pid].length_bp, ra.complete_to_bp)
    if span.walk is not None:
        span.walk = tuple(w for _, w in src.key) + (span.walk,)
    span.complete = False
    span.junctions = tuple(js)
    r = span.route
    if r is not None:
        # a merges span: certified as its route (which may cross the junction)
        r.complete = False
        r.certified_to_bp = min(r.certified_to_bp, cert)
        span.certified_to_bp = r.certified_to_bp
    else:
        span.certified_to_bp = min(span.certified_to_bp, cert)
    return span


def _span(un, t, pi, part, node, length, with_labels, at_bp=None, b=None):
    """Part |pi|'s PairSpan on the union walk ending at trie node |node|. |b|: its memory
    (the span and route, _M_SPAN; the labels, by _labels_for()) and the part's derivations
    it reads, at their cold price. A stitched part answers from the source that holds the
    deepest of its spellings on the walk (PairUnion._resolve())."""
    g0 = _root(part.graphlet)
    if t.side not in g0.arms:
        if b is not None:
            b.charge(0, _M_SPAN)
        return PairSpan(part.tag, g0.mode, 'none', 0, None, None, None, None, (), {},
                        None if not with_labels else (), False, False, 0, 'under_strategy',
                        g0.k, g0.regime)
    # how far the pair's spellings follow the walk: to trie node |stop| of source |si|;
    # |below|: the walk's union segments past it (what a route reads), never the whole path
    si, stop, below = un._resolve(t, pi, node)
    src = un._src[si]
    g = src.graphlet
    a = g.arms[t.side]
    off = src.offset
    mode = g.mode
    if b is not None:
        b.charge(0, _M_SPAN)
        derive.uses(b, g, a, 'first_paths', *(('end_labels',) if mode == 'constrain'
                                              else ()))
    d = t.depth[stop]
    how, item = _case(t, si, stop, length)
    if how == 'whole':
        _, leaf, pid = item
        reasons, code = _leaf_reasons(g, a, leaf)
        complete, cert = _leaf_cert(a, leaf)
        labels, cut = _labels_to(un, src, a, leaf, length, with_labels, at_bp, b)
        return _past_junction(un, src, PairSpan(
            part.tag, mode, 'whole', length, pid, None, None, None, (), reasons, labels, cut,
            complete, off + cert, None, g.k, g.regime))
    last = t.last_seg(si, a, stop)
    if how == 'extends':
        pid = _ops._first_paths(a)[last] if last is not None else None
        cert = off + min(length - off, a.complete_to_bp)
        labels, cut = _labels_to(un, src, a, last, length, with_labels, at_bp, b)
        return _past_junction(un, src, PairSpan(
            part.tag, mode, 'extends', length, pid, None, None, None, (),
            {'path_reason': None, 'labels': {}}, labels, cut,
            length - off <= a.complete_to_bp, cert, None, g.k, g.regime))
    # the union walk's base after the stop (none where the walk ends there)
    c = t.useg[below[0]].first if below else None
    if how == 'ends':
        _, leaf, pid = item
        reasons, code = _leaf_reasons(g, a, leaf)
        complete, cert = _leaf_cert(a, leaf)
        absence = _absence_at(src, g, a, _events_at(t, si, a, stop, last), d - off, c, code)
        labels, cut = _labels_to(un, src, a, leaf, d, with_labels, at_bp, b)
        return _past_junction(un, src, PairSpan(
            part.tag, mode, 'ends', d, pid, None, None, None, (), reasons, labels, cut,
            complete, off + cert, absence, g.k, g.regime))
    if how == 'merges':
        _, stub, m = item
        ms = a.segments[m]
        into = t.node_name(t.seg_nodes[si][m][0], ms.length_bp,
                           ms.walk[0] if ms.length_bp and ms.walk else None)
        j = ms.parents.index(stub)
        cont = tuple(un._pair_label(src, x) for x in (ms.partition[j] if j < len(
            ms.partition) else ()))
        route, rsrc = _route(un, t, si, src, a, stub, m, j, d, below, into, with_labels,
                             at_bp, b)
        labels, cut = _labels_to(un, src, a, stub, d, with_labels, at_bp, b)
        # the source the route ends in: a stitched child when it crossed a junction
        return _past_junction(un, rsrc, PairSpan(
            part.tag, mode, 'merges', d, None, stub, into, route, cont,
            {'path_reason': None, 'labels': {}}, labels, cut, route.complete,
            route.certified_to_bp, None, g.k, g.regime))
    # diverges: the pair goes on with other bases (or a stub/walk elsewhere)
    absence = _absence_at(src, g, a, _events_at(t, si, a, stop, last), d - off, c)
    labels, cut = _labels_to(un, src, a, last, d, with_labels, at_bp, b)
    return _past_junction(un, src, PairSpan(
        part.tag, mode, 'diverges', d, None, None, None, None, (),
        {'path_reason': None, 'labels': {}}, labels, cut, False,
        off + min(d - off, a.complete_to_bp), absence, g.k, g.regime))


def _fp_chain(a, seg):
    """[(seg, parent index used)] root -> |seg| along first parents."""
    if seg is None:
        return []
    return [(x, 0) for x in derive.chain(a, seg)]


def _route(un, t, pi, part, a, stub, m, j, d, below, into, with_labels, at_bp, b=None):
    """The pair's route from its merge |m| (entered from the stub through parent |j|)
    along the union walk's remaining bases, those of the union segments |below| the stop
    (§2.4 'Routes through merges') -> (Route, the source it ends in). |pi|, |part|: the
    source; |d| and the route's positions outward from the seed. A route that reaches the
    end of a walk a stitched child continues goes on in the child (§2.3): the junction is
    read, and the route is certified at most to the first junction."""
    # the union walk's remaining bases, walking order
    rest = []
    for u in below:
        us = t.useg[u]
        pi0 = min(us.cover)
        g0 = un._src[pi0].graphlet
        off0 = un._src[pi0].offset
        for seg, f, to in us.pieces[pi0]:
            lo = max(f, d)
            if lo < to:
                s0 = g0.arms[t.side].segments[seg]
                rest.append(s0.walk[lo - off0 - s0.from_bp:to - off0 - s0.from_bp])
    rest = ''.join(rest)
    chain = [(x, 0) for x in derive.chain(a, stub)]
    route = []
    # per source the route reads: (source index, source, arm, its chain and route)
    legs = [(pi, part, a, chain, route)]
    g, segs, off = part.graphlet, a.segments, part.offset
    merges = [(d, into)]
    cur, cur_j = m, j
    n = 0                       # bases of |rest| matched
    none = {'path_reason': None, 'labels': {}}
    reasons, code = none, None
    held_leaf = None
    while True:
        s = segs[cur]
        route.append((cur, cur_j))
        w = s.walk or ''
        i = 0
        while i < len(w) and n < len(rest) and w[i] == rest[n]:
            i += 1
            n += 1
        kid = via = None
        if i == len(w) and t.continued:
            if s.leaf is not None:
                kid = t.continued.get((pi, cur))
            elif (not any(segs[c].length_bp or segs[c].children for c in s.children)
                  if n == len(rest) else _step(segs, cur, rest[n]) is None):
                # the walk a child continues ends at a zero-length leaf child of this
                # segment (a merge at the radius): nothing else here spells the next base
                via = next((c for c in s.children if not segs[c].length_bp
                            and segs[c].leaf is not None and (pi, c) in t.continued), None)
                if via is not None:
                    kid = t.continued[(pi, via)]
        if kid is not None:
            # the route reaches the end of a walk a stitched child continues: it goes on
            # in the child's tree, from its root segment (its offset is this leaf's end),
            # as a span there does (PairUnion._resolve())
            if via is not None:
                route.append((via, segs[via].parents.index(cur)))
                _note_merge(t, pi, segs[via], merges, off)
            pi, part = kid, un._src[kid]
            g, a, off = part.graphlet, part.graphlet.arms[t.side], part.offset
            segs = a.segments
            route = []
            legs.append((pi, part, a, [], route))
            cur, cur_j = next(x.id for x in segs if not x.parents), 0
            continue
        if n == len(rest):
            how = 'held'
            goes_on = i < len(w) or any(segs[c].length_bp or segs[c].children
                                        for c in s.children)
            if not goes_on:
                # the route ends at a part leaf (this segment's, or a zero-length leaf
                # child's: a merge node there): it takes that leaf's reasons and certificate
                leaf = cur if s.leaf is not None else next(
                    (c for c in s.children if segs[c].leaf is not None), None)
                if leaf is not None:
                    reasons, code = _leaf_reasons(g, a, leaf)
                    held_leaf = leaf
            break
        if i < len(w):
            how, goes_on = 'diverges', True
            break
        step = _step(segs, cur, rest[n])
        if step is None:
            leaf = cur if s.leaf is not None else next(
                (c for c in s.children if not segs[c].length_bp and segs[c].leaf is not None),
                None)
            if leaf is not None:
                reasons, code = _leaf_reasons(g, a, leaf)
                how, goes_on = 'ends', False
                if leaf != cur:
                    route.append((leaf, segs[leaf].parents.index(cur)))
                    cur = leaf
            else:
                how, goes_on = 'diverges', True
            break
        for x, xj in step[:-1]:
            route.append((x, xj))
            _note_merge(t, pi, segs[x], merges, off)
        cur, cur_j = step[-1]
        _note_merge(t, pi, segs[cur], merges, off)
    pos = d + n
    # certification
    if how == 'ends':
        complete, cert = _leaf_cert(a, cur)
        cert += off
    elif held_leaf is not None:
        complete, cert = _leaf_cert(a, held_leaf)
        cert += off
    else:
        cert = off + min(pos - off, a.complete_to_bp)
        complete = how == 'held' and pos - off <= a.complete_to_bp
    absence = None
    if how != 'held':
        c = rest[n] if n < len(rest) else None
        lpos = pos - off
        here = [segs[x] for x, _ in route if segs[x].from_bp <= lpos <= segs[x].end_bp]
        absence = _absence_at(part, g, a, here, lpos, c, code)
    labels, lcut = (), False
    if with_labels:
        # the labels at the route's end (or |at_bp|) are those of the leg holding that
        # depth: never joined across a junction (D.1); their cut flag is the route's (GPT
        # review M2), not the stub's
        e = pos if at_bp is None else min(pos, at_bp)
        k = len(legs) - 1
        while k > 0 and e <= legs[k][1].offset:
            k -= 1
        _, lsrc, la, lchain, lroute = legs[k]
        labels, lcut = _labels_src(un, lsrc, la, lambda: lchain + lroute, pos, at_bp, b)
    return Route(pos, how, goes_on, reasons, absence, tuple(merges), labels, complete,
                 cert, lcut), part


def _step(segs, cur, ch):
    """[(segment, parent index)] from |cur| to its child spelling |ch| first, through
    zero-length segments (a merge node that splits at once), or None: at most one child of
    a node spells each next base, since a prefix names one path per graph."""
    for c in segs[cur].children:
        cs = segs[c]
        if cs.length_bp:
            if cs.walk and cs.walk[0] == ch:
                return [(c, cs.parents.index(cur))]
        elif cs.children:
            sub = _step(segs, c, ch)
            if sub:
                return [(c, cs.parents.index(cur))] + sub
    return None


def _note_merge(t, pi, s, merges, off=0):
    if s.is_merge:
        first = s.walk[0] if s.length_bp and s.walk else None
        merges.append((off + s.from_bp, t.node_name(t.seg_nodes[pi][s.id][0], s.length_bp,
                                                    first)))


# ------------------------------------------------------------------ span labels

def _labels_for(un, part, a, chain, end, with_labels, at_bp=None, b=None):
    """(SpanLabels, cut) at outward depth |end| (or |at_bp| when smaller) along |chain|
    [(segment, the parent index it was entered through)], root first; (None, False) when
    the query does not compute labels. |b|: the segments read at their cold price, the
    labels returned and, while it runs, its working sets (released)."""
    if not with_labels:
        return None, False
    if at_bp is not None:
        end = min(end, at_bp)
    g = part.graphlet
    segs = a.segments
    # trim the chain to the segments at or before |end|: a segment with bases from |end| on,
    # and a zero-length one beyond it (a merge node or leaf deeper than the cut: its empty
    # intervals would hide the labels alive at |end|, GPT review m1)
    k = len(chain)
    while k > 1:
        s = segs[chain[k - 1][0]]
        if not (s.from_bp >= end if s.length_bp else s.from_bp > end):
            break
        k -= 1
    chain = chain[:k]
    if not chain:
        return (), False
    info = un._seg_labels(part, a, b)
    last = segs[chain[-1][0]]
    cut = any(info(x)[2] for x, _ in chain)
    if b is not None:
        # at most the labels alive anywhere in the last segment: those returned, and the
        # working sets (out, open, starts) held while the chain is read
        li, le = info(last.id)[:2]
        n = max(len(li), len(le))
        b.charge(W_ELEM * len(chain), (_M_LABEL + _M_LABEL_WORK) * n)
    if end >= last.end_bp:
        alive = set(info(last.id)[1])
    else:
        alive = {lab for lab, f, to in info(last.id)[0] if f <= end - 1 < to}
    out = {}
    open_ = alive
    pos = end
    constrain = g.mode == 'constrain'
    # linear in the labels of each segment of the chain and in the merges' partitions: the
    # intervals and the end set of a segment are built once per union (_seg_labels), each
    # boundary tests every open label once against a set
    for idx in range(len(chain) - 1, -1, -1):
        sid, j = chain[idx]
        s = segs[sid]
        if not open_:
            break
        if s.length_bp:
            starts = {}
            for lab, f, to in info(sid)[0]:
                if f < pos <= to and lab in open_:
                    starts[lab] = f
            for lab in list(open_):
                f = starts.get(lab)
                if f is None or f > s.from_bp:
                    out[lab] = f if f is not None else pos
                    open_.discard(lab)
        if not open_:
            break
        if idx == 0:
            for lab in open_:
                out[lab] = s.from_bp
            open_ = set()
            break
        # the boundary into |s| from its parent chain[idx - 1]
        par = segs[chain[idx - 1][0]]
        keep = info(par.id)[1]
        if constrain and s.is_merge:
            keep = keep & frozenset(s.partition[j] if j < len(s.partition) else ())
        for lab in [x for x in open_ if x not in keep]:
            out[lab] = s.from_bp
            open_.discard(lab)
        pos = par.end_bp
    labs = sorted(out.items(), key=lambda x: (x[1], x[0]))
    got = tuple(SpanLabel(un._pair_label(part, lab), f) for lab, f in labs)
    if b is not None:
        b.release(_M_LABEL_WORK * n)
    return got, cut


def _walk_bytes(d, n):
    """A UnionWalk's own bytes: the record, its path (the list of union segments built, then
    a tuple of their ids) and its tuples of n spans and pairs."""
    return 300 + (LIST_ITEM + 8) * d + 3 * 8 * n


def _span_bytes(s):
    """What a built PairSpan holds (at most what _span() and _labels_for() charged)."""
    n = len(s.labels or ()) + (len(s.route.labels or ()) if s.route is not None else 0)
    return _M_SPAN + _M_LABEL * n


def _carrying_labels(s):
    """The span labels of a span that holds its walk: its own (whole, extends) or its
    route's (merges, the route held to the walk's end); () for any other span."""
    if s.how in ('whole', 'extends'):
        return s.labels or ()
    if s.how == 'merges' and s.route is not None and s.route.how == 'held':
        return s.route.labels or ()
    return ()


def _shown_labels(s):
    """(labels, cut) a span shows: a merges span its route's (the labels where the route
    ends, with the route's cut flag), any other its own."""
    if s.how == 'merges':
        r = s.route
        return (r.labels, r.labels_cut) if r is not None else ((), False)
    return s.labels, s.labels_cut


def _carrier_line(s):
    if s.how in ('whole', 'extends'):
        mid = 'to %d' % s.to_bp
    elif s.how == 'merges':
        r = s.route
        mid = 'at %d into %s; route %s to %d%s' % (
            s.to_bp, s.into, r.how, r.to_bp,
            ', goes on' if r.how == 'held' and r.goes_on else
            '' if r.how == 'held' else
            (' %s' % (r.end_reasons.get('path_reason') or '') if r.how == 'ends' else ''))
    elif s.how == 'none':
        mid = 'arm not walked'
    else:
        mid = 'at %d' % s.to_bp
        if s.how == 'ends' and s.end_reasons.get('path_reason'):
            mid += '  %s' % s.end_reasons['path_reason']
    labs, cut = _shown_labels(s)
    lab = ''
    if labs is not None and s.how in ('whole', 'extends', 'merges', 'ends'):
        n = len(labs)
        whole = sum(1 for x in labs if x.from_bp == 0)
        lab = 'labels %d%s' % (n, '+' if cut else '')
        if n:
            lab += ', %d of them from 0 bp' % whole if whole != n else ', from 0 bp'
    tail = []
    if s.how in ('whole', 'extends'):
        tail.append('complete' if s.complete else 'certified to %d' % s.certified_to_bp)
    absence = s.absence if s.how != 'merges' else (s.route.absence if s.route else None)
    if absence:
        tail.append('absence: %s' % absence)
    if s.mode != 'annotate':
        tail.append(s.mode)
    return '%-6s %-9s %-34s %-34s %s' % (s.pair, s.how, mid, lab, '  '.join(tail))
