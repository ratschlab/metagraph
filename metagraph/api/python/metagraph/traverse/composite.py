"""Combining graphlets: what the cross-pair union (union.py) and stitching (_stitch.py) share
(docs/DESIGN-graphlet-merge.md, Appendix F.2).

A combination is a container of unchanged parts plus a stated relation. It never builds an
MGT body: every answer is computed from the parts' own bodies, which are held by reference
and never copied or modified (§4.4). This module holds the pieces both combinations use:

  * CombineRefused (a ValueError with a .code) and NotCombinable (also a TypeError), the
    refusals of both combinations (Appendix A.2, B.1);
  * the shared preconditions: whole retrievals (no views), the envelope, the bases, one
    alphabet, the oriented seed, and the index rule of compare() without its seed rule
    (ops.index_rule_check(), which stitching words as different_index);
  * retrieval identity (Appendix A.4): what tells two retrievals of one pair apart, and
    the first field path where two of them differ;
  * the weakest-of evidence over several parts (Appendix C.3);
  * digests (canonical JSON, sha256) and resume tokens bound to a combination's digest;
  * the budget plumbing (budget.py's rules: charged before it is built, a stop raises);
  * the streaming speller over pieces of the parts' segments.

Nothing here reads a certificate or derives a walk: those stay the parts' (ops, derive).
"""

import hashlib
import json

from . import budget as _B
from . import derive
from . import ops as _ops
from .budget import BASE_SHIFT, STR, TUPLE, W_ELEM, W_ROW, list_bytes, str_bytes

__all__ = [
    'CombineRefused', 'NotCombinable',
    'is_view', 'require_whole', 'require_envelope', 'has_bases', 'alphabet_key',
    'case_mapped', 'reverse_complement', 'seed_difference',
    'index_check', 'different_index_clause', 'INDEX_CLAUSES',
    'retrieval_identity', 'first_difference',
    'OUTCOME_RANKS', 'weakest_outcome', 'combined_evidence',
    'canonical_json', 'sha256_hex', 'part_digest', 'combined_digest',
    'token', 'resume_pos', 'run_scoped', 'models_memory_bytes',
    'iter_piece_bases', 'spell_pieces',
]


# ------------------------------------------------------------------ refusals

class CombineRefused(ValueError):
    """A combination refused: .code names the rule (Appendix A.2 for the union, B.1 for
    stitching), the message is its text. MissingEnvelope and BadSelector are raised as
    today, not as this class."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


class NotCombinable(CombineRefused, TypeError):
    """not_combinable: an input that is no result a combination reads (a TypeError, so a
    caller catching either class sees it)."""

    def __init__(self, message):
        super().__init__('not_combinable', message)


# ------------------------------------------------------------------ preconditions

def is_view(g):
    """Whether |g| is a GraphletView (a filtered view of one retrieval), not a whole one."""
    return isinstance(g, _ops.GraphletView)


def require_whole(g, text):
    """view_unsupported with |text| (each combination states its own remedy) for a view."""
    if is_view(g):
        raise CombineRefused('view_unsupported', text)


def require_envelope(g, what):
    """MissingEnvelope with the library's text (Graphlet.require_envelope) for a body
    without its J line: |what| names the operation ('union_pairs()', 'stitch()')."""
    g.require_envelope(what)


def has_bases(g):
    """Whether every arm of |g| carries its bases (output.sequences was not false)."""
    return all(_ops._has_bases(g, side) for side in g.arms)


def _letters(alphabet):
    return (alphabet or '').replace('$', '')


def alphabet_key(g):
    """The alphabet two parts compare by: without the sentinel '$', so that a succinct
    graph's '$ACGT' and a hash graph's 'ACGT' are one DNA alphabet."""
    return _letters(g.alphabet)


def case_mapped(seq, alphabet):
    """|seq| after the server's case mapping: upper case where the alphabet's letters are
    upper case (every DNA build), unchanged otherwise."""
    letters = _letters(alphabet)
    if letters and letters == letters.upper():
        return seq.upper()
    return seq


_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reverse_complement(seq):
    """The reverse complement of a DNA sequence (A, C, G, T, N in either case)."""
    return seq.translate(_COMPLEMENT)[::-1]


def _is_dna(alphabet):
    return set(_letters(alphabet).upper()) <= set('ACGTN')


def seed_difference(a, b, alphabet=None):
    """How oriented seed |b| differs from |a| (both already case mapped) -> None when they
    are equal, else one of ('revcomp',) -- |b| is |a|'s reverse complement on a DNA
    alphabet --, ('length', len a, len b, the common prefix) or ('first', the first
    position where they differ), in that order of precedence: the texts of
    different_seed (Appendix A.2) are worded from it."""
    if a == b:
        return None
    if alphabet is not None and _is_dna(alphabet) and len(a) == len(b) \
            and reverse_complement(a) == b:
        return ('revcomp',)
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    if len(a) != len(b):
        return ('length', len(a), len(b), i)
    return ('first', i)


# ------------------------------------------------------------------ index identity

# stitching's different_index clauses (Appendix B.1), keyed by ops.index_rule_check()'s keys
INDEX_CLAUSES = {
    'meta_fp': 'index_meta_fp differs',
    'fp': 'index_fp differs',
    'ns': 'index namespace %r vs %r',
    'release': 'release %r vs %r',
    'names': 'label %s is %r in one and %r in the other',
}


def index_check(a, b):
    """compare()'s index rule without its seed rule (ops.index_rule_check()) -> (key,
    values): proven different ('meta_fp', 'fp', 'ns', 'release', 'names'), unverifiable
    ('unverifiable_both', 'unverifiable_one': no index manifest digest on both sides or on
    one), or 'same' (verified by index_fp)."""
    return _ops.index_rule_check(a, b)


def different_index_clause(key, values):
    """The different_index clause of a proven difference (INDEX_CLAUSES), or None for an
    unverifiable or verified identity."""
    clause = INDEX_CLAUSES.get(key)
    return None if clause is None else (clause % values if values else clause)


# ------------------------------------------------------------------ retrieval identity

def retrieval_identity(g):
    """What makes |g| one retrieval (Appendix A.4), without the pair: {index: {fp,
    meta_fp, ns, release}, strategy: the normalized strategy from the envelope, all of it,
    seed: {sequence (case mapped, in the request orientation), labels_from_seed, labels:
    the seed labels' refs (constrain, stated seed labels)}}. It leaves out the caller's
    seed_id, and needs only the envelope (GraphletStore.request_of() may return None
    where seed label names do not resolve back; this does not). |g| carries its
    envelope (require_envelope() first)."""
    meta = (g.seed_summary or {}).get('seed') or {}
    from_seed = bool(meta.get('labels_from_seed'))
    seed = {'sequence': case_mapped(g.seed.sequence, g.alphabet),
            'labels_from_seed': from_seed}
    if g.mode == 'constrain' and not from_seed:
        seed['labels'] = [l.ref for l in g.labels[:g.seed.num_seed_labels]]
    return {'index': _ops.index_identity(g),
            'strategy': (g.envelope or {}).get('strategy') or {},
            'seed': seed}


def first_difference(a, b, path=''):
    """The first field path, in byte order of the keys, where JSON-like values |a| and |b|
    differ ('strategy.branching.max_splits_per_path'; a list index as '[3]'), or None when
    they are equal. A key on one side only differs at that key."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=lambda x: str(x).encode('utf-8')):
            sub = '%s.%s' % (path, k) if path else str(k)
            if k not in a or k not in b:
                return sub
            got = first_difference(a[k], b[k], sub)
            if got is not None:
                return got
        return None
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        for i in range(min(len(a), len(b))):
            got = first_difference(a[i], b[i], '%s[%d]' % (path, i))
            if got is not None:
                return got
        return None if len(a) == len(b) else '%s[%d]' % (path, min(len(a), len(b)))
    if type(a) is not type(b) or a != b:
        return path or '.'
    return None


# ------------------------------------------------------------------ evidence

# per outcome dimension, each value's rank: the weakest value has the largest rank
OUTCOME_RANKS = {
    'walks': {'complete': 0, 'partial': 1, 'failed': 2},
    'branch_diagnostics': {'complete': 0, 'cut': 1},
    'label_evidence': {'complete': 0, 'lower_bound': 1, 'qualified': 2},
    'delivery': {'inline': 0, 'spooled': 1, 'paged': 2},
}


def weakest_outcome(outcomes):
    """The weakest value per dimension over |outcomes| (Outcome.as_dict() dicts, or None
    for a part without one, which is skipped) -> a dict in OUTCOME_RANKS' order, or None
    when no outcome was given. A value outside the ranks is weaker than every known one."""
    got = {}
    for o in outcomes:
        if o is None:
            continue
        for dim, ranks in OUTCOME_RANKS.items():
            v = o.get(dim)
            r = ranks.get(v, len(ranks))
            if dim not in got or r > got[dim][0]:
                got[dim] = (r, v)
    if not got:
        return None
    return {dim: got[dim][1] for dim in OUTCOME_RANKS if dim in got}


def combined_evidence(parts, side=None):
    """The evidence of an answer that reads several parts (Appendix C.3): |parts| is a
    list of (key, graphlet, verified) in canonical order -> {parts: {key: the part's own
    evidence_block (ops.evidence_block(), for |side| where the part walked it)}, outcome:
    the weakest per dimension, exact: every part read is exact and its identity verified}.
    A part that did not walk |side| contributes its block over the arms it walked."""
    blocks = {}
    exact = True
    for key, g, verified in parts:
        s = side if side is None or side in g.arms else None
        blk = _ops.evidence_block(g, s)
        blocks[key] = blk
        exact = exact and bool(blk.get('exact')) and bool(verified)
    return {'parts': blocks,
            'outcome': weakest_outcome([b.get('outcome') for b in blocks.values()]),
            'exact': exact and bool(blocks)}


# ------------------------------------------------------------------ digests

def canonical_json(obj):
    """The canonical JSON text of |obj|: keys sorted, no spaces, ASCII escapes."""
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                      default=str)


def sha256_hex(text):
    return hashlib.sha256(text.encode('utf-8') if isinstance(text, str) else text).hexdigest()


def part_digest(g):
    """The digest of a part's body: the sha-256 of its records without the J line
    (Graphlet.body_digest, made by the parse; the canonical text's otherwise)."""
    return _ops._graphlet_key(g)


def combined_digest(items):
    """The sha-256 of the canonical JSON of |items| (a list in canonical order)."""
    return sha256_hex(canonical_json(items))


# ------------------------------------------------------------------ budgets

def run_scoped(budget, op, fn, levers=()):
    """budget.run_scoped(): fn(None) without a budget, else fn(b) inside b.scope(op)."""
    return _B.run_scoped(budget, op, levers, fn)


def _token_digest(digest, op, args):
    return hashlib.sha256(json.dumps([op, digest, args], sort_keys=True,
                                     default=str).encode()).hexdigest()[:16]


def token(digest, op, args, *pos):
    """A resume token of a stopped list of a combination: bound to the operation, its
    arguments and the combination's |digest|, so it is refused on any other."""
    return (op, _token_digest(digest, op, args)) + tuple(pos)


def resume_pos(digest, op, args, resume, n):
    """The |n| ints a resume token carries; zeros without one. A token of another
    operation, other arguments or another combination is a ValueError."""
    if resume is None:
        return (0,) * n
    if not (isinstance(resume, (tuple, list)) and len(resume) == 2 + n and resume[0] == op
            and all(isinstance(x, int) and not isinstance(x, bool) and x >= 0
                    for x in resume[2:])):
        raise ValueError('resume= takes the Partial.resume of a stopped %s() call, not %r'
                         % (op, resume))
    if resume[1] != _token_digest(digest, op, args):
        raise ValueError('this resume token belongs to a %s() call with other arguments or '
                         'on another combination: repeat that call\'s arguments with it' % op)
    return tuple(resume[2:])


def models_memory_bytes(models, extra=0):
    """The heap footprint of the parts' models (ops' deduplicating deep sizeof, so a model
    object given twice counts once) plus |extra| (a combination's own index)."""
    return _ops._deep_bytes(list(models)) + extra


# ------------------------------------------------------------------ spelling

def _piece_bases(g, side, seg, lo, hi):
    s = g.arms[side].segments[seg]
    if s.walk is None:
        raise ValueError(derive.NO_BASES)
    if not s.from_bp <= lo <= hi <= s.end_bp:
        raise ValueError('piece [%d, %d) lies outside segment %d [%d, %d)'
                         % (lo, hi, seg, s.from_bp, s.end_bp))
    return s.walk[lo - s.from_bp:hi - s.from_bp]


def iter_piece_bases(pieces, b=None):
    """The bases of |pieces|, (graphlet, side, segment id, from_bp, to_bp) in outward
    positions of that segment, one piece at a time in walking order (outward from the
    seed), read from the parts' segments, never copied beyond the piece. |b| (a resolved
    budget or None): each piece charged before it is read."""
    for g, side, seg, lo, hi in pieces:
        if b is not None:
            b.charge(W_ELEM + ((hi - lo) >> BASE_SHIFT), str_bytes(hi - lo))
        yield _piece_bases(g, side, seg, lo, hi)


def spell_pieces(pieces, side, orientation='natural', b=None):
    """The bases of |pieces| as one string: in walking order (orientation 'walk'), or in
    natural orientation, the left arm's reversed as ops.spell() writes it ('natural')."""
    if orientation not in ('natural', 'walk'):
        raise ValueError("orientation is 'natural' or 'walk'")
    if b is not None:
        b.charge(W_ROW, TUPLE + list_bytes(0) + STR)
    w = ''.join(iter_piece_bases(pieces, b))
    return w if orientation == 'walk' or side == 'right' else w[::-1]
