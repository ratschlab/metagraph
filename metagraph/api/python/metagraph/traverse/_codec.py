"""The token codecs of MGT v1 (DESIGN-traverse-graphlet.md §2.1, §2.2).

Every decoder here is canonical-only at the token level: a token is accepted iff it is
the encoding of the value it denotes (the rule of tests/data/traverse/codec_vectors.tsv).
A decoder built on a native parser would accept language-specific extras (Python's
float() takes '1_0', ' 1', 'nan'), and two readers that accept different sloppy
spellings disagree about which documents are valid. The format's only freedoms are
field-level -- '*' vs an explicit value, SETEXPR explicit vs delta, a front-coding
prefix shorter than the longest -- and those decode and dump canonical (§2.3).

Stdlib only.
"""

import decimal
import math
import re
from array import array

__all__ = [
    'GraphletFormatError', 'CodecError', 'UNLIMITED',
    'fmt_int', 'parse_int', 'fmt_bool', 'parse_bool', 'fmt_num', 'parse_num',
    'encode_ranges', 'decode_ranges', 'encode_setexpr', 'decode_setexpr',
    'pct_escape', 'pct_unescape', 'front_encode', 'front_decode',
    'encode_kvalue', 'decode_kvalue', 'encode_pairs', 'decode_pairs',
    'REASON', 'REASON_CODE', 'QUAL', 'QUAL_CODE', 'RESOURCE_CODES',
    'end_token', 'parse_end_token', 'LabelSetInterner',
]

MAX_U64 = (1 << 64) - 1


class GraphletFormatError(ValueError):
    """A document that is not valid MGT v1. |line_no| is 1-based (0: the document as a
    whole, e.g. a missing trailer)."""

    def __init__(self, line_no, msg):
        self.line_no = line_no
        self.msg = msg
        super().__init__('line %d: %s' % (line_no, msg) if line_no else msg)


class CodecError(ValueError):
    """A token that is not the canonical encoding of any value; the parser turns it into
    a GraphletFormatError naming the line."""


# how much of a token an error message repeats: a corrupt or hostile body must not make
# its own error message huge (a 5,000-digit token would turn the tool layer's
# format_error into result_too_large and hide the error class)
TOKEN_ECHO = 32


def tok(x, n=TOKEN_ECHO):
    """|x| as an error message shows it: the repr of its first |n| characters, and its
    length when it was cut. An int beyond 64 bits is described, not printed (its decimal
    text may exceed sys.get_int_max_str_digits() and raise while formatting)."""
    if isinstance(x, str):
        if len(x) <= n:
            return repr(x)
        return '%r... (%d characters)' % (x[:n], len(x))
    if isinstance(x, int) and not isinstance(x, bool) and x.bit_length() > 64:
        return 'an integer of %d bits' % x.bit_length()
    s = repr(x)
    return s if len(s) <= n else '%s... (%d characters)' % (s[:n], len(s))


class _Unlimited:
    """The K value 'u'. A sentinel rather than the string 'unlimited', which is a legal
    string value of its own ('s:unlimited')."""

    __slots__ = ()

    def __repr__(self):
        return 'UNLIMITED'

    def __reduce__(self):
        return 'UNLIMITED'


UNLIMITED = _Unlimited()


# ------------------------------------------------------------------- integers, booleans

_UINT_RE = re.compile(r'(?:0|[1-9][0-9]*)\Z')
# every u64 has at most 20 digits: a longer digit string is out of range before int()
# sees it (int() raises a bare ValueError beyond sys.get_int_max_str_digits(), which
# escaped parse() as no CodecError)
_U64_DIGITS = 20


def fmt_int(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 0 or n > MAX_U64:
        raise CodecError('not an unsigned 64-bit integer: %s' % tok(n))
    return str(n)


def parse_int(t):
    if not _UINT_RE.match(t):
        raise CodecError('integer %s' % tok(t))
    if len(t) > _U64_DIGITS:
        raise CodecError('integer %s does not fit 64 bits' % tok(t))
    n = int(t)
    if n > MAX_U64:
        raise CodecError('integer %s does not fit 64 bits' % tok(t))
    return n


def fmt_bool(b):
    return '1' if b else '0'


def parse_bool(t):
    if t == '1':
        return True
    if t == '0':
        return False
    raise CodecError('boolean %s' % tok(t))


# --------------------------------------------------------------------------- floats

def fmt_num(x):
    """§2.1 (v3): the shortest round-trip significant digits of the double, expanded
    positionally. repr() gives the unique shortest digits (nearest on ties), the same
    digits as C++ std::to_chars(..., scientific); Decimal 'f' expands them exactly, so
    1e23 is '100000000000000000000000', not to_chars(fixed)'s exact binary value."""
    x = float(x)
    if math.isnan(x) or x < 0:   # -0 < 0 is false: -0 is in the domain, written 0
        raise CodecError('outside the MGT float domain: %s' % tok(x))
    if x == math.inf:
        return 'inf'
    s = format(decimal.Decimal(repr(x)), 'f')
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return '0' if s == '-0' else s


_FLOAT_RE = re.compile(r'(?:0|[1-9][0-9]*)(?:\.[0-9]*[1-9])?\Z')


def parse_num(t):
    if t != 'inf' and not _FLOAT_RE.match(t):
        raise CodecError('float %s' % tok(t))
    x = float(t)
    # re-encoding is the canonical check: it rejects digits that are not the shortest
    # (the exact expansion of 1e23, '0.10000000000000001'), overflow to inf and
    # underflow to 0
    if fmt_num(x) != t:
        raise CodecError('float %s is not the canonical text of %s' % (tok(t), tok(x)))
    return x


# --------------------------------------------------------------------------- RANGES

def encode_ranges(ids):
    """Ascending ids with every maximal run of >= 2 consecutive ids written 'a-b' (the
    golden vectors' decision: {5,6} is '5-6'); '.' is the empty set."""
    out = []
    start = prev = None
    for i in ids:
        if prev is not None and i <= prev:
            raise CodecError('RANGES ids not strictly ascending: %s after %s' % (tok(i), tok(prev)))
        if start is None:
            start = i
        elif i != prev + 1:
            out.append(str(start) if start == prev else '%d-%d' % (start, prev))
            start = i
        prev = i
    if start is None:
        return '.'
    if start is not None and start < 0:
        raise CodecError('negative id')
    out.append(str(start) if start == prev else '%d-%d' % (start, prev))
    return ','.join(out)


_RANGES_ITEM_RE = re.compile(r'(0|[1-9][0-9]*)(?:-(0|[1-9][0-9]*))?\Z')


def decode_ranges(t, bound=None):
    """-> ascending ids. |bound|: every id must be below it, checked BEFORE a run is
    expanded, so that a few bytes of a corrupt or hostile token ('0-3000000000') cannot
    allocate in proportion to the range (§14: a parse that cannot complete fails
    locally); the reader passes the number of L records."""
    if t == '.':
        return []
    ids = []
    last_end = None
    for item in t.split(','):
        m = _RANGES_ITEM_RE.match(item)
        if not m:
            raise CodecError('RANGES item %s' % tok(item))
        if max(len(m.group(1)), len(m.group(2) or '')) > _U64_DIGITS:
            raise CodecError('RANGES id %s does not fit 64 bits' % tok(item))
        a = int(m.group(1))
        if m.group(2) is not None:
            b = int(m.group(2))
            if b <= a:
                raise CodecError('RANGES run %s' % tok(item))
        else:
            b = a
        if last_end is not None and a <= last_end + 1:
            # descending, overlapping, or adjacent items that collapse into one run
            raise CodecError('RANGES %s is not ascending and collapsed' % tok(t))
        if b > MAX_U64:
            raise CodecError('RANGES id does not fit 64 bits')
        if bound is not None and b >= bound:
            raise CodecError('label id %d beyond the %d L records' % (b, bound))
        ids.extend(range(a, b + 1))
        last_end = b
    # a lone pair 'a,a+1' is adjacent and was rejected above; one more spelling
    # check keeps the decoder honest against the encoder
    return ids


# --------------------------------------------------------------------------- SETEXPR

def encode_setexpr(ids, base):
    """The shorter of the explicit RANGES and the delta against |base| ('!', '!R', '!+A',
    '!R+A'), a tie explicit; no base (the root's G.entry) -> explicit only."""
    explicit = encode_ranges(ids)
    if base is None:
        return explicit
    s = set(ids)
    b = set(base)
    removed = sorted(b - s)
    added = sorted(s - b)
    delta = '!' + (encode_ranges(removed) if removed else '') \
        + ('+' + encode_ranges(added) if added else '')
    return delta if len(delta) < len(explicit) else explicit


def decode_setexpr(t, base, bound=None):
    """-> sorted list of ids. A delta whose parts do not fit |base| is an error, never a
    different set: it means the reader derived another base than the writer. |bound| as
    in decode_ranges."""
    if not t.startswith('!'):
        return decode_ranges(t, bound)
    if base is None:
        raise CodecError('a SETEXPR delta needs a base: %s' % tok(t))
    body = t[1:]
    rem_tok, plus, add_tok = body.partition('+')
    # an empty part is omitted, never written '.': one spelling per delta
    if rem_tok == '.' or (plus and add_tok in ('', '.')):
        raise CodecError('empty SETEXPR part in %s' % tok(t))
    removed = decode_ranges(rem_tok, bound) if rem_tok else []
    added = decode_ranges(add_tok, bound) if plus else []
    b = set(base)
    if not b.issuperset(removed) or not b.isdisjoint(added):
        raise CodecError('SETEXPR %s does not fit its base' % tok(t))
    if removed:
        b.difference_update(removed)
    if added:
        b.update(added)
    return sorted(b)


# ---------------------------------------------------------------- free-text fields

_PCT = {'%': '%25', '\n': '%0A', '\r': '%0D'}
_UNPCT = {'%25': '%', '%0A': '\n', '%0D': '\r'}


def pct_escape(s):
    """The free-text last field (L/X names, Q message, K effect): '%', LF and CR are
    escaped, nothing else -- TAB, NUL, other controls and non-ASCII stay raw."""
    if '%' not in s and '\n' not in s and '\r' not in s:
        return s
    return ''.join(_PCT.get(c, c) for c in s)


def pct_unescape(t):
    if '\n' in t or '\r' in t:
        raise CodecError('raw LF/CR in a free-text field')
    if '%' not in t:
        return t
    out = []
    i = 0
    n = len(t)
    while i < n:
        c = t[i]
        if c == '%':
            esc = t[i:i + 3]
            if esc not in _UNPCT:
                raise CodecError('percent escape %s' % tok(esc))
            out.append(_UNPCT[esc])
            i += 3
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def _front_prefix(pb, nb):
    p = 0
    m = min(len(pb), len(nb))
    while p < m and pb[p] == nb[p]:
        p += 1
    # back off to a code-point boundary: a continuation byte (10xxxxxx) at p means p
    # is inside a character, and the suffix must be valid UTF-8 on its own
    while 0 < p < len(nb) and (nb[p] & 0xC0) == 0x80:
        p -= 1
    return p


def front_encode(prev, name):
    """-> '<prefix_len> <suffix>': the L record's tail after <seq_id|*>. prefix_len counts
    BYTES of the previous name's UTF-8, cut back to a character boundary. The separating
    space is always written, so an empty suffix leaves the record ending in a space."""
    pb = prev.encode('utf-8')
    nb = name.encode('utf-8')
    p = _front_prefix(pb, nb)
    return '%d %s' % (p, pct_escape(nb[p:].decode('utf-8')))


def front_prefix_len(prev, name):
    return _front_prefix(prev.encode('utf-8'), name.encode('utf-8'))


def front_decode_parts(prev, p_tok, suffix_tok):
    p = parse_int(p_tok)
    pb = prev.encode('utf-8')
    if p > len(pb):
        raise CodecError('front-coding prefix %d beyond the previous name' % p)
    try:
        return (pb[:p] + pct_unescape(suffix_tok).encode('utf-8')).decode('utf-8')
    except UnicodeDecodeError:
        raise CodecError('front-coding prefix %d ends inside a character' % p) from None


_FRONT_RE = re.compile(r'(0|[1-9][0-9]*) (.*)\Z', re.S)


def front_decode(prev, t):
    m = _FRONT_RE.match(t)
    if not m:
        raise CodecError('front-coded name %s' % tok(t))
    return front_decode_parts(prev, m.group(1), m.group(2))


# --------------------------------------------------------------------------- K values

def _s_raw(byte):
    return 0x21 <= byte <= 0x7E and byte not in (0x25, 0x2C)


def encode_kvalue(v):
    """Typed VALUE (§2.2 v5.1): int -> 'i:', float -> 'f:', str -> 's:' (every byte
    outside 0x21..0x7E, '%' and ',' as uppercase %XX), UNLIMITED -> 'u'."""
    if v is UNLIMITED:
        return 'u'
    if isinstance(v, bool):
        raise CodecError('a K value is not a boolean: %s' % tok(v))
    if isinstance(v, int):
        return 'i:' + fmt_int(v)
    if isinstance(v, float):
        return 'f:' + fmt_num(v)
    if isinstance(v, str):
        return 's:' + ''.join(chr(b) if _s_raw(b) else '%%%02X' % b
                              for b in v.encode('utf-8'))
    raise CodecError('not a K value: %s' % tok(v))


_HEX2_RE = re.compile(r'[0-9A-F]{2}\Z')


def decode_kvalue(t):
    if t == 'u':
        return UNLIMITED
    kind, colon, body = t.partition(':')
    if not colon or kind not in ('i', 'f', 's'):
        raise CodecError('K value %s' % tok(t))
    if kind == 'i':
        return parse_int(body)
    if kind == 'f':
        return parse_num(body)
    raw = bytearray()
    i = 0
    n = len(body)
    while i < n:
        c = body[i]
        if c == '%':
            if not _HEX2_RE.match(body[i + 1:i + 3]):
                raise CodecError('K string escape in %s' % tok(t))
            b = int(body[i + 1:i + 3], 16)
            if _s_raw(b):
                raise CodecError('K string escapes a raw byte: %s' % tok(t))
            raw.append(b)
            i += 3
        else:
            if ord(c) > 0x7F or not _s_raw(ord(c)):
                raise CodecError('K string holds a byte that must be escaped: %s' % tok(t))
            raw.append(ord(c))
            i += 1
    try:
        return raw.decode('utf-8')
    except UnicodeDecodeError:
        raise CodecError('K string is not UTF-8: %s' % tok(t)) from None


# ------------------------------------------------------------------ interval pairs

def encode_pairs(pairs):
    """X <runs>: 'a-b,c-d' (k-mer runs on the seed, stored as given), '.' when empty."""
    if not pairs:
        return '.'
    return ','.join('%s-%s' % (fmt_int(a), fmt_int(b)) for a, b in pairs)


def decode_pairs(t):
    if t == '.':
        return []
    out = []
    for item in t.split(','):
        a, dash, b = item.partition('-')
        if not dash:
            raise CodecError('interval %s' % tok(item))
        out.append((parse_int(a), parse_int(b)))
    return out


# ------------------------------------------------------------------ end-reason codes

# traversal_types.cpp:37-55, plus (v4) resource_limit
REASON = {
    'D': 'dead_end', 'L': 'label_lost', 'B': 'loss_budget', 'R': 'branch',
    'U': 'edge_reuse', 'V': 'edge_reuse_rc', 'J': 'rejoined_seed', 'T': 'trace_break',
    'X': 'max_extension_bp', 'S': 'max_steps', 'P': 'max_live_paths', 'N': 'max_paths',
    'O': 'max_output_bp', 'M': 'time_budget', 'W': 'beam_pruned', 'Y': 'resource_limit',
}
REASON_CODE = {v: k for k, v in REASON.items()}
# the walker's EndReason order (traversal_types.hpp): the canonical order of the
# per-reason maps (B label_ends), as the C++ iterates them
REASON_ORDER = 'DLBRUVJTXSPNOMWY'
# the second letter keeps the walker's text qualifier (walker.cpp: quorum texts at the
# BRANCH ends, the hairpin dead end, the label_lost texts)
QUAL = {
    'Rm': 'minority', 'Rb': 'below_min_labels', 'Rs': 'split_limit', 'Dh': 'hairpin',
    'Ls': 'superseded', 'Lx': 'switch_sources',
}
QUAL_CODE = {(tok[0], text): tok for tok, text in QUAL.items()}
# exploration of the requested domain incomplete (is_resource_stop) plus resource_limit
RESOURCE_CODES = frozenset('SPNOMWY')


def end_token(code, qualifier=None):
    if qualifier is None:
        return code
    t = QUAL_CODE.get((code, qualifier))
    if t is None:
        raise CodecError('no code for qualifier %s of %s' % (tok(qualifier), tok(code)))
    return t


def parse_end_token(t):
    """R <end>: -> (code, qualifier, silent, merged); code None for 'm'."""
    if t == 'm':
        return None, None, False, True
    if t == 'Lw':
        return 'L', None, True, False
    if len(t) == 1 and t in REASON:
        return t, None, False, False
    if t in QUAL:
        return t[0], QUAL[t], False, False
    raise CodecError('end token %s' % tok(t))


# -------------------------------------------------------------------- label sets

class LabelSetInterner:
    """Content-hashed label sets: equal sets share one array('I') (4 bytes per id instead
    of a list of ints), the reason a label-free arm with thousands of segments over a
    few hundred distinct sets stays small (§7). The arrays must be treated as
    immutable. hits/misses are what graphlet_measure.py reports."""

    __slots__ = ('_table', 'hits', 'misses', 'ids_stored')

    _EMPTY = array('I')

    def __init__(self):
        self._table = {}
        self.hits = 0
        self.misses = 0
        self.ids_stored = 0

    def intern(self, ids):
        if not ids:
            self.hits += 1
            return self._EMPTY
        a = ids if isinstance(ids, array) else array('I', ids)
        key = a.tobytes()
        got = self._table.get(key)
        if got is not None:
            self.hits += 1
            return got
        self.misses += 1
        self.ids_stored += len(a)
        self._table[key] = a
        return a

    def __len__(self):
        return len(self._table)

    def nbytes(self):
        return sum(len(k) * 2 + 64 for k in self._table)
