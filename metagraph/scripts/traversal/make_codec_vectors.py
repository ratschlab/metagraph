#!/usr/bin/env python3
"""Write (or --check) the shared MGT v1 codec golden vectors.

DESIGN-traverse-graphlet.md §2.7: before MGT v1 is frozen, the C++ writer and the
Python library must both reproduce every vector of
api/python/tests/data/traverse/codec_vectors.tsv byte-exactly. This script is the
reference that produced the expected strings. It is deliberately self-contained
(stdlib only, no import of metagraph.traverse): the library is tested against it,
not against itself. The file header it writes (HEADER below) is the format contract
both readers follow; the codec rules it fixes are the design's §2.1/§2.2, with the
points the design leaves open decided there and marked "decision".

usage: make_codec_vectors.py [--out PATH] [--check]
"""

import argparse
import decimal
import json
import math
import os
import re
import struct
import sys
import unicodedata

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_OUT = os.path.join(REPO, 'api', 'python', 'tests', 'data', 'traverse',
                           'codec_vectors.tsv')

XORSHIFT_SEED = 0x9E3779B97F4A7C15
N_DIFFERENTIAL = 10000
MASK64 = (1 << 64) - 1


# --------------------------------------------------------------------------------------
# The reference codec. Every decoder is canonical-only at the token level: a token is
# accepted iff it is the encoding of the value it denotes. The only valid non-canonical
# spellings are field-level choices (SETEXPR explicit vs delta, a shorter front-coding
# prefix), which decode and re-encode to the canonical form (§2.3).
# --------------------------------------------------------------------------------------

class Reject(ValueError):
    pass


def bits_to_float(h):
    return struct.unpack('<d', struct.pack('<Q', int(h, 16)))[0]


def float_to_bits(x):
    return '%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]


def fmt_float(x):
    """§2.1 (v3): shortest round-trip digits, expanded positionally."""
    if math.isnan(x) or x < 0:   # -0 < 0 is false: -0 is in the domain
        raise ValueError('outside the MGT float domain: %r' % x)
    if x == math.inf:
        return 'inf'
    # repr() gives the unique shortest round-trip digits (nearest on ties), the same
    # digits as std::to_chars(..., scientific); Decimal 'f' expands them exactly.
    s = format(decimal.Decimal(repr(x)), 'f')
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return '0' if s == '-0' else s


_FLOAT_TOKEN = re.compile(r'(?:0|[1-9][0-9]*)(?:\.[0-9]*[1-9])?\Z')


def parse_float(tok):
    if tok != 'inf' and not _FLOAT_TOKEN.match(tok):
        raise Reject('float syntax')
    x = float(tok)
    # Re-encoding is the canonical check: it rejects digits that are not the shortest
    # (the exact binary expansion of 1e23, "0.10000000000000001") and overflow to inf.
    if fmt_float(x) != tok:
        raise Reject('not the canonical float text')
    return x


def enc_ranges(ids):
    ids = list(ids)
    assert all(a < b for a, b in zip(ids, ids[1:])) and all(i >= 0 for i in ids)
    if not ids:
        return '.'
    out = []
    i = 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if i == j else '%d-%d' % (ids[i], ids[j]))
        i = j + 1
    return ','.join(out)


_UINT = r'(?:0|[1-9][0-9]*)'
_RANGES_ITEM = re.compile(r'(%s)(?:-(%s))?\Z' % (_UINT, _UINT))


def dec_ranges(tok):
    if tok == '.':
        return []
    ids = []
    for item in tok.split(','):
        m = _RANGES_ITEM.match(item)
        if not m:
            raise Reject('RANGES item %r' % item)
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) is not None else a
        if m.group(2) is not None and b <= a:
            raise Reject('RANGES run %r' % item)
        if ids and a <= ids[-1]:
            raise Reject('RANGES not ascending')
        ids.extend(range(a, b + 1))
    if enc_ranges(ids) != tok:
        raise Reject('RANGES not collapsed')
    return ids


def enc_setexpr(ids, base):
    explicit = enc_ranges(ids)
    if base is None:
        return explicit
    s, b = set(ids), set(base)
    removed, added = sorted(b - s), sorted(s - b)
    delta = '!' + (enc_ranges(removed) if removed else '') + \
        ('+' + enc_ranges(added) if added else '')
    return delta if len(delta) < len(explicit) else explicit


def dec_setexpr(tok, base):
    if not tok.startswith('!'):
        return dec_ranges(tok)
    if base is None:
        raise Reject('delta without a base')
    body = tok[1:]
    rem_tok, plus, add_tok = body.partition('+')
    # an empty part is omitted, never written '.': one spelling per delta
    if rem_tok == '.' or (plus and add_tok in ('', '.')):
        raise Reject('empty delta part')
    removed = dec_ranges(rem_tok) if rem_tok else []
    added = dec_ranges(add_tok) if plus else []
    b = set(base)
    # a part that does not fit its base means the reader computed another base
    # than the writer: make that an error rather than a silently different set
    if not set(removed) <= b or set(added) & b:
        raise Reject('delta does not fit its base')
    return sorted((b - set(removed)) | set(added))


_PCT_TEXT = {'%': '%25', '\n': '%0A', '\r': '%0D'}
_PCT_TEXT_DEC = {'%25': '%', '%0A': '\n', '%0D': '\r'}


def pct_text(s):
    return ''.join(_PCT_TEXT.get(c, c) for c in s)


def unpct_text(tok):
    out = []
    i = 0
    while i < len(tok):
        c = tok[i]
        if c in '\n\r':
            raise Reject('raw LF/CR')
        if c == '%':
            esc = tok[i:i + 3]
            if esc not in _PCT_TEXT_DEC:
                raise Reject('percent escape %r' % esc)
            out.append(_PCT_TEXT_DEC[esc])
            i += 3
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def enc_front(prev, name):
    pb, nb = prev.encode('utf-8'), name.encode('utf-8')
    p = 0
    while p < min(len(pb), len(nb)) and pb[p] == nb[p]:
        p += 1
    # back off to a code-point boundary: a continuation byte (10xxxxxx) at p means p is
    # inside a character, and the suffix must be valid UTF-8 on its own
    while 0 < p < len(nb) and (nb[p] & 0xC0) == 0x80:
        p -= 1
    return '%d %s' % (p, pct_text(nb[p:].decode('utf-8')))


_FRONT = re.compile(r'(%s) (.*)\Z' % _UINT, re.S)


def dec_front(prev, tok):
    m = _FRONT.match(tok)
    if not m:
        raise Reject('front-coded name syntax')
    p = int(m.group(1))
    pb = prev.encode('utf-8')
    if p > len(pb):
        raise Reject('prefix beyond the previous name')
    try:
        return (pb[:p] + unpct_text(m.group(2)).encode('utf-8')).decode('utf-8')
    except UnicodeDecodeError:
        raise Reject('prefix ends inside a character')


def _s_raw(byte):
    return 0x21 <= byte <= 0x7E and byte not in (0x25, 0x2C)


def enc_kvalue(kind, value):
    if kind == 'i':
        n = int(value)
        assert 0 <= n <= MASK64
        return 'i:%d' % n
    if kind == 'f':
        return 'f:' + fmt_float(bits_to_float(value))
    if kind == 's':
        return 's:' + ''.join(chr(b) if _s_raw(b) else '%%%02X' % b
                              for b in value.encode('utf-8'))
    if kind == 'u':
        assert value is None
        return 'u'
    raise AssertionError(kind)


_HEX2 = re.compile(r'[0-9A-F]{2}\Z')


def dec_kvalue(tok):
    if tok == 'u':
        return ('u', None)
    kind, colon, body = tok.partition(':')
    if not colon or kind not in ('i', 'f', 's'):
        raise Reject('K value type')
    if kind == 'i':
        if not re.match(_UINT + r'\Z', body) or int(body) > MASK64:
            raise Reject('K integer')
        return ('i', body)
    if kind == 'f':
        return ('f', float_to_bits(parse_float(body)))
    raw = bytearray()
    i = 0
    while i < len(body):
        c = body[i]
        if c == '%':
            if not _HEX2.match(body[i + 1:i + 3]):
                raise Reject('K string escape')
            b = int(body[i + 1:i + 3], 16)
            if _s_raw(b):
                raise Reject('K string: escaped a raw byte')
            raw.append(b)
            i += 3
        else:
            if ord(c) > 0x7F or not _s_raw(ord(c)):
                raise Reject('K string: raw byte that must be escaped')
            raw.append(ord(c))
            i += 1
    try:
        return ('s', raw.decode('utf-8'))
    except UnicodeDecodeError:
        raise Reject('K string: not UTF-8')


# --------------------------------------------------------------------------------------
# The vectors
# --------------------------------------------------------------------------------------

def xorshift64_bits(seed, n):
    """Non-negative finite doubles over raw bits (Marsaglia's xorshift64, 13/7/17)."""
    x = seed
    out = []
    while len(out) < n:
        x ^= (x << 13) & MASK64
        x ^= x >> 7
        x ^= (x << 17) & MASK64
        bits = x & 0x7FFFFFFFFFFFFFFF
        if bits >> 52 == 0x7FF:   # inf / NaN: outside the domain, draw again
            continue
        out.append('%016x' % bits)
    return out


def named_floats():
    """The named edge cases of §2.7 (v3), plus a few that pin the Python reference's
    two repr() styles (positional below 1e16, scientific from 1e16 / below 1e-4)."""
    v = [
        (0.0, '0'), (-0.0, '-0 (written 0; decodes to +0)'), (1.0, '1'), (0.5, '0.5'),
        (0.1, '0.1'), (0.0001, '0.0001 (v1: to_chars and repr disagreed)'),
        (1e-05, '1e-05 (repr switches to scientific)'), (2.5e-07, '2.5e-07'),
        (1 / 3, '1/3'), (math.inf, 'inf'), (1e16, '1e16 (repr switches to scientific)'),
        (1e22, '1e22 (largest exact power of ten)'),
        (1e23, '1e23 (v2: to_chars fixed wrote 99999999999999991611392)'),
        (1.2345678901234568e20, '1.2345678901234568e20'), (2.0 ** 53, '2^53'),
        (2.0 ** 53 + 2, '2^53+2'),
        (2.0 ** 53 - 1, '2^53-1'), (0.1 + 0.2, '0.1+0.2'), (0.3, '0.3'), (1.5, '1.5'),
        (0.75, '0.75'), (2.0, '2'), (10.0, '10 (only fractional zeros are stripped)'),
        (64.0, '64'), (100.0, '100'), (30000.0, '30000 (time_budget_ms default)'),
        (123456.789, '123456.789'), (1e15, '1e15'),
        (9999999999999998.0, 'largest double below 1e16'), (1e21, '1e21'),
        (2.0 ** 63, '2^63: shortest digits, not the exact 9223372036854775808'),
        (2.0 ** 64, '2^64: shortest digits, not the exact 18446744073709551616'),
        (2.0 ** 1023, '2^1023'),
        (sys.float_info.max, 'DBL_MAX'),
        (math.nextafter(sys.float_info.max, 0), 'DBL_MAX predecessor'),
        (sys.float_info.min, 'smallest normal 2^-1022'),
        (bits_to_float('000fffffffffffff'), 'largest subnormal'),
        (bits_to_float('0000000000000002'), 'second smallest subnormal'),
        (bits_to_float('0000000000000001'), 'smallest subnormal 2^-1074'),
    ]
    for k in range(-20, 26):
        x = float('1e%d' % k)
        v.append((x, '10^%d' % k))
        v.append((math.nextafter(x, 0.0), '10^%d predecessor' % k))
        v.append((math.nextafter(x, math.inf), '10^%d successor' % k))
    seen, out = set(), []
    for x, note in v:
        h = float_to_bits(x)
        if h not in seen:
            seen.add(h)
            out.append((h, note))
    return out


def setexpr_vectors():
    def r(a, b):
        return list(range(a, b + 1))

    evens = list(range(0, 19, 2))
    return [
        # no base: the root's entry (annotate mode; constrain roots are '*' or explicit)
        ([0, 1, 2], None, 'G.entry of a root: no base, explicit RANGES only'),
        ([], None, 'G.entry of a root: no base, empty'),
        ([5], None, 'G.entry of a root: no base, single'),
        # equal to the base
        ([0, 1, 2], [0, 1, 2], "G.end, no label ended (base = entry set): '!'"),
        ([12], [12], "equal, explicit '12' longer than '!'"),
        ([0], [0], "equal, tie '0' vs '!': explicit"),
        ([], [], "equal and empty, tie '.' vs '!': explicit"),
        # removals only
        ([0, 1], [0, 1, 2], "split child (base = parent[0]'s end set): '!2' (§2.6)"),
        ([2], [0, 1, 2], "split child (base = parent[0]'s end set): explicit wins"),
        ([0, 2], [0, 1, 2], 'removal inside a run'),
        ([0, 1, 2, 3], [0, 1, 2, 3, 4], 'removal at the end of a run'),
        ([], [0, 1, 2], "G.end, every label ended: '.' beats '!0-2'"),
        ([i for i in range(1000) if i != 500], list(range(1000)),
         'P after P (base = previous P): one label leaves'),
        (evens, evens + [20], "removal of a lone id: '!20'"),
        ([0, 1, 5], [0, 1, 2, 3, 4, 5], "removal of a run: '!2-4'"),
        # additions only
        ([0, 1], [0], "addition, tie '0-1' vs '!+1': explicit"),
        ([3], [], "addition to an empty base: explicit"),
        (r(0, 15) + [17], r(0, 15), "G.end with a switch-in (base = entry set): '!+17'"),
        (r(0, 20), r(0, 9) + r(11, 20), "addition closing a gap, tie '0-20' vs '!+10': explicit"),
        (r(0, 9) + r(20, 29), r(0, 9), "addition of a run: '!+20-29'"),
        # removals and additions
        (r(0, 4) + r(6, 9) + [12], r(0, 9), "first P (base = entry set): '!5+12'"),
        (r(0, 9) + [11], r(0, 10), "tie '0-9,11' vs '!10+11': explicit"),
        ([1, 3], [1, 2], "explicit '1,3' shorter than '!2+3'"),
        ([5], [3], "disjoint singles: explicit"),
        (evens + [21], evens + [20], "'!20+21' beats the long explicit list"),
        ([i for i in range(41) if i != 17] + [41], r(0, 40),
         "P after P (base = previous P): one leaves, one joins"),
        ([0, 2] + r(10, 19), [0, 1, 2, 3] + r(10, 18), "'!1,3+19'"),
        # the normative bases, named
        ([0, 1, 3, 4], [0, 1, 3, 4],
         "G.entry of a merge (base = union of the parents' end sets {0,1}|{3,4}): '!'"),
        ([0, 1, 3], [0, 1, 3, 4], "G.entry of a merge (base = union {0,1}|{3,4}): '!4'"),
        ([0, 2], [0, 1, 2], "C.labels (base = the leaf's end set): '!1'"),
        ([0, 1, 2], [0, 1, 2], "C.labels equal to the leaf's end set: '!'"),
        ([7], [7, 8], "C.labels: explicit '7' shorter than '!8'"),
    ]


def setexpr_decode_vectors():
    return [
        ('0-2', [0, 1, 2], '!', "valid non-canonical: explicit where '!' is shorter"),
        ('!', [0], '0', 'valid non-canonical: a tie written as the delta'),
        ('!+1', [0], '0-1', 'valid non-canonical: a tie written as the delta'),
        ('!', [], '.', "valid non-canonical: '!' against an empty base"),
        ('!0-1+5', [0, 1, 2], '2,5', 'valid non-canonical: delta longer than explicit'),
        ('0-2', None, '0-2', 'explicit without a base'),
        ('!', None, None, 'a delta needs a base (the root has none)'),
        ('!0', None, None, 'a delta needs a base'),
        ('!.', [0, 1], None, "empty part written '.'"),
        ('!+.', [0, 1], None, "empty part written '.'"),
        ('!0+.', [0, 1], None, "empty part written '.'"),
        ('!.+5', [0, 1], None, "empty part written '.'"),
        ('!0+', [0, 1], None, 'empty ADDED after +'),
        ('!5', [0, 1], None, 'REMOVED not in the base'),
        ('!+1', [0, 1], None, 'ADDED already in the base'),
        ('!0+0', [0, 1], None, 'the same id removed and added'),
        ('!!', [0, 1], None, 'malformed'),
        ('+1', [0], None, 'malformed: ADDED without !'),
        ('!0+1+2', [0], None, 'malformed: two + parts'),
        ('!1,0', [0, 1], None, 'malformed RANGES in a part'),
        ('!0,1', [0, 1, 2], None, 'non-collapsed RANGES in a part'),
        ('', [0], None, 'empty token'),
    ]


FRONT = [
    ('', 'acc1', 'first L of a document: previous name is ""'),
    ('acc1', 'acc2', '§2.6'),
    ('éfoo', 'ébar', 'prefix in bytes: é is 2 bytes'),
    ('éfoo', 'èfoo', 'common byte C3 is inside é/è: shortened to 0'),
    ('aéb', 'aèb', 'common bytes a,C3: shortened to the boundary after a'),
    ('日本', '日本語', '3-byte characters'),
    ('日x', '旦x', 'E6 97 A5 vs E6 97 A6: two common bytes inside a character'),
    ('🧬a', '🧬b', '4-byte character kept'),
    ('🧬a', '🧭a', 'F0 9F A7 AC vs F0 9F A7 AD: three common bytes inside a character'),
    ('e\u0301x', 'e\u0302y',
     'code-point boundary, not grapheme: suffix starts with a combining mark'),
    ('ACC1', 'ACC1',
     'same name twice (annotate: two columns): empty suffix, line ends in a space'),
    ('日本語', '日本', 'name is a prefix of the previous one: empty suffix'),
    ('abc', 'ab', 'empty suffix'),
    ('abc', '', 'empty name'),
    ('', '', 'empty name first'),
    ('abc', 'abd', ''),
    ('acc 1', 'acc 2', 'spaces in names'),
    ('x', ' lead', 'leading space: two spaces after prefix_len'),
    ('trail ', 'trail  ', 'trailing spaces preserved'),
    ('', '*', "name '*' is not 'absent'"),
    ('*', '.', "name '.' is not 'empty list'"),
    ('.', 'c:0', "name 'c:0' is not a ref"),
    ('c:0', 'c:1', ''),
    ('c:1', '*', ''),
    ('a%b', 'a%c', 'prefix counts raw bytes, not escaped text'),
    ('a%', 'a%25', 'a literal %25 in a name'),
    ('50%', '50%x', ''),
    ('a\nb', 'a\nc', 'LF is one raw byte of the prefix'),
    ('x', 'a\nb', 'LF in the suffix is escaped'),
    ('a\r', 'a\rb', 'CR in the prefix'),
    ('tab\there', 'tab\tthere', 'TAB stays raw'),
    ('C:\\data\\a', 'C:\\data\\b', 'backslashes stay raw'),
    ('NZ_STEQ01000045.1', 'NZ_STEQ01000046.1', 'accession headers'),
]

FRONT_DECODE = [
    ('abc', '0 abd', '2 d', 'valid non-canonical: prefix shorter than possible'),
    ('abc', '1 bd', '2 d', 'valid non-canonical'),
    ('éfoo', '0 ébar', '2 bar', 'valid non-canonical'),
    ('abc', '3 ', '3 ', 'empty suffix'),
    ('éfoo', '1 bar', None, 'prefix ends inside é'),
    ('abc', '4 x', None, 'prefix beyond the previous name'),
    ('', '1 a', None, 'first L must have prefix 0'),
    ('abc', '02 d', None, 'leading zero'),
    ('abc', '-1 d', None, 'negative'),
    ('abc', '2', None, 'separator missing'),
    ('abc', ' 2 d', None, 'leading space'),
    ('abc', 'x d', None, 'not a number'),
    ('abc', '2 d%41', None, 'percent escape of a byte that is written raw'),
    ('abc', '2 d%0a', None, 'lowercase escape'),
]

PCT = [
    ('', 'empty'), ('acc1', ''), ('NZ_STEQ01000045.1', ''),
    ('50%', ''), ('%', ''), ('%25', 'a literal %25'), ('%0A', 'a literal %0A'),
    ('a\nb', 'LF'), ('a\rb', 'CR'), ('\r\n', 'CR LF'), ('a\r', 'trailing CR'),
    ('a\tb', 'TAB stays raw'), ('with space', ''), (' lead and trail ', 'spaces kept'),
    ('*', ''), ('.', ''), ('c:0', ''), ('h:0:1', ''), ('éfoo', 'UTF-8 stays raw'),
    ('日本語 🧬', ''), ('C:\\data\\x', 'backslash stays raw'), ('a\x00b', 'NUL stays raw'),
    ('a\x7fb', 'DEL stays raw'),
    ('v\x0bf\x0cs\x1cg\x1dr\x1en\x85l\u2028p\u2029',
     'VT FF FS GS RS NEL LS PS stay raw: split records on LF only (not str.splitlines())'),
    ('Branch decisions at or beyond complete_to_bp are not reported; '
     'raise "output.max_branch_events".', 'a K effect / Q message'),
]

PCT_DECODE = [
    ('%25%0A%0D', '%25%0A%0D', 'all three escapes'),
    ('%250', '%250', "'%' then '0'"),
    ('%41', None, 'escape of a byte that is written raw'),
    ('%0a', None, 'lowercase escape'),
    ('%2', None, 'truncated escape'),
    ('%', None, 'truncated escape'),
    ('a%', None, 'truncated escape'),
    ('%%', None, 'malformed escape'),
    ('%2G', None, 'malformed escape'),
    ('%09', None, 'TAB is written raw'),
    ('a\nb', None, 'raw LF'),
    ('a\rb', None, 'raw CR'),
]

KVALUE = [
    (['i', '0'], ''), (['i', '1'], ''), (['i', '64'], 'max_labels_per_node default'),
    (['i', '1000000'], ''), (['i', '18446744073709551615'], '2^64-1'),
    (['f', float_to_bits(0.0)], ''), (['f', float_to_bits(-0.0)], '-0'),
    (['f', float_to_bits(30000.0)], 'time_budget_ms default'),
    (['f', float_to_bits(0.5)], ''), (['f', float_to_bits(1e23)], ''),
    (['f', float_to_bits(2.5e-07)], ''), (['f', float_to_bits(math.inf)], ''),
    (['s', 'merge'], 'the scope limitation (§2.2 v5.1)'),
    (['s', 'unlimited'], "the STRING 'unlimited', not u"),
    (['s', ''], 'empty string'), (['s', 'a b'], 'space'),
    (['s', 'a,b'], 'comma: the extra-list separator'), (['s', '50%'], ''),
    (['s', 'x=y'], "'=' stays raw: extra items split at the first '='"),
    (['s', 'a\tb'], 'TAB'), (['s', 'a\nb\r'], 'LF CR'),
    (['s', 'é'], 'non-ASCII: fixed fields are ASCII'),
    (['s', '*'], ''), (['s', '.'], ''), (['s', '|'], ''), (['s', 'c:0'], ''),
    (['s', 'header 1, plasmid'], ''), (['s', '\x7f'], 'DEL'),
    (['s', '!~'], '0x21 and 0x7E stay raw'),
    (['s', 'trace'], 'derivation limit (spec §7.0)'), (['s', 'NZ_STEQ01000045.1'], 'a header'),
    (['u', None], 'unlimited'),
]

KVALUE_DECODE = [
    ('U', 'unknown'), ('u:', 'unknown'), ('unlimited', 'unknown'), ('', 'empty'),
    ('x:1', 'unknown type'), ('S:merge', 'type is lower case'), ('i', 'missing colon'),
    ('i:', 'empty integer'), ('i:-1', 'negative'), ('i:+1', 'sign'), ('i:01', 'leading zero'),
    ('i:1.5', 'not an integer'), ('i:18446744073709551616', '2^64 does not fit'),
    ('f:', 'empty float'), ('f:1.0', 'non-canonical float'), ('f:-0', '-0 is written 0'),
    ('f:-1', 'negative'), ('f:nan', 'NaN'), ('f:1e5', 'scientific'),
    ('s:a b', 'raw space'), ('s:a,b', 'raw comma'), ('s:%41', 'escape of a raw byte'),
    ('s:%c3%a9', 'lowercase hex'), ('s:é', 'raw non-ASCII'), ('s:%FF', 'not UTF-8'),
    ('s:%2', 'truncated escape'), ('s:a\tb', 'raw TAB'),
]

FLOAT_DECODE = [
    ('', 'empty'), ('1.0', 'trailing zero'), ('1.', 'trailing point'), ('.5', 'no integer part'),
    ('0.0', 'trailing zero'), ('0.50', 'trailing zero'), ('01', 'leading zero'),
    ('00.5', 'leading zero'), ('-0', '-0 is written 0'), ('-1', 'negative'), ('+1', 'sign'),
    ('1e5', 'scientific'), ('1E5', 'scientific'), ('1e+23', 'scientific'), ('nan', 'NaN'),
    ('NaN', 'NaN'), ('-inf', 'negative'), ('+inf', 'sign'), ('Infinity', 'spelling'),
    ('INF', 'spelling'), ('1_0', 'Python float() accepts it'), (' 1', 'whitespace'),
    ('0x1p0', 'hex float'), ('99999999999999991611392', '1e23 by to_chars(fixed): not shortest'),
    ('0.1000000000000000055511151231257827021181583404541015625', 'exact value of 0.1'),
    ('0.10000000000000001', 'parses to 0.1, not the shortest digits'),
    ('1' + '0' * 400, 'overflows to inf'),
    ('0.' + '0' * 400 + '1', 'underflows to 0'),
    ('9007199254740993', '2^53+1 rounds to 2^53'),
]

RANGES = [
    [], [0], [7], [0, 1], [5, 6], [9, 10], [0, 1, 2], [0, 2], [1, 3, 5, 7],
    [0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12], [0, 1, 3, 4], [2, 3, 5], [0, 1, 2, 5, 6, 9],
    [99, 100, 101], [8, 9, 10, 11], [1000000], [4294967293, 4294967294],
    list(range(0, 21, 2)), list(range(1000)),
]
RANGES_NOTES = {
    (0, 1): 'decision: every maximal run of >= 2 ids is collapsed',
    (5, 6): 'a run of two',
    (0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12): '§2.1 example',
}

RANGES_DECODE = [
    ('', 'empty token'), (',', ''), ('0,', ''), (',0', ''), ('0,,1', ''), ('1,0', 'descending'),
    ('0,0', 'duplicate'), ('0-0', 'run of one'), ('1-0', 'reversed run'),
    ('0-1,1', 'overlap'), ('0-1,2', "adjacent runs: '0-2'"), ('0,1', "run of two: '0-1'"),
    ('0,1,2', "run: '0-2'"), ('00', 'leading zero'), ('01', 'leading zero'),
    ('-1', 'negative'), ('+1', 'sign'), ('0--1', ''), ('0-1-2', ''), ('a', ''),
    ('.,0', ''), ('0,.', ''), ('..', ''),
]


def build():
    """[(tag, input, expected, note)]; the encoders produce every expected string."""
    rows = []

    def add(tag, inp, exp, note=''):
        rows.append((tag, inp, exp, note))

    for h, note in named_floats():
        add('float', h, fmt_float(bits_to_float(h)), note)
    for tok, note in FLOAT_DECODE:
        add('float_decode', tok, _canon(lambda: fmt_float(parse_float(tok))), note)
    for ids in RANGES:
        add('ranges', ids, enc_ranges(ids), RANGES_NOTES.get(tuple(ids), ''))
    for tok, note in RANGES_DECODE:
        add('ranges_decode', tok, _canon(lambda: enc_ranges(dec_ranges(tok))), note)
    for ids, base, note in setexpr_vectors():
        add('setexpr', [ids, base], enc_setexpr(ids, base), note)
    for tok, base, exp, note in setexpr_decode_vectors():
        add('setexpr_decode', [tok, base],
            _canon(lambda: enc_setexpr(dec_setexpr(tok, base), base)), note)
        assert rows[-1][2] == exp, (tok, base, rows[-1][2], exp)
    for s, note in PCT:
        add('pct', s, pct_text(s), note)
    for tok, exp, note in PCT_DECODE:
        add('pct_decode', tok, _canon(lambda: pct_text(unpct_text(tok))), note)
        assert rows[-1][2] == exp, (tok, rows[-1][2], exp)
    for prev, name, note in FRONT:
        add('front', [prev, name], enc_front(prev, name), note)
    for prev, tok, exp, note in FRONT_DECODE:
        add('front_decode', [prev, tok], _canon(lambda: enc_front(prev, dec_front(prev, tok))),
            note)
        assert rows[-1][2] == exp, (prev, tok, rows[-1][2], exp)
    for inp, note in KVALUE:
        add('kvalue', inp, enc_kvalue(*inp), note)
    for tok, note in KVALUE_DECODE:
        add('kvalue_decode', tok, _canon(lambda: enc_kvalue(*dec_kvalue(tok))), note)
    for h in xorshift64_bits(XORSHIFT_SEED, N_DIFFERENTIAL):
        add('float', h, fmt_float(bits_to_float(h)), None)

    _self_check(rows)
    return rows


def _canon(f):
    try:
        return f()
    except Reject:
        return None


def _self_check(rows):
    """Every encode vector must decode back to its input (the obligation stated in the
    header), and the decode vectors that are meant to fail must fail."""
    for tag, inp, exp, _ in rows:
        if tag == 'float':
            want = inp if inp != '8000000000000000' else '0000000000000000'
            assert float_to_bits(parse_float(exp)) == want, (inp, exp)
        elif tag == 'ranges':
            assert dec_ranges(exp) == inp, (inp, exp)
        elif tag == 'setexpr':
            assert dec_setexpr(exp, inp[1]) == inp[0], (inp, exp)
        elif tag == 'pct':
            assert unpct_text(exp) == inp, (inp, exp)
        elif tag == 'front':
            assert dec_front(inp[0], exp) == inp[1], (inp, exp)
        elif tag == 'kvalue':
            kind, value = dec_kvalue(exp)
            assert kind == inp[0] and (value == inp[1] or (
                kind == 'f' and inp[1] == '8000000000000000' and value == '0' * 16)), (inp, exp)
    expected_rejects = {'float_decode', 'ranges_decode', 'kvalue_decode'}
    for tag, inp, exp, note in rows:
        if tag in expected_rejects:
            assert exp is None, ('meant to be rejected', tag, inp, exp)


# --------------------------------------------------------------------------------------
# The file
# --------------------------------------------------------------------------------------

HEADER = """\
MGT v1 codec golden vectors: the freeze gate of docs/DESIGN-traverse-graphlet.md §2.7 (v5.1).
GENERATED by scripts/traversal/make_codec_vectors.py -- do not edit by hand.
  regenerate: python3 scripts/traversal/make_codec_vectors.py
  verify:     python3 scripts/traversal/make_codec_vectors.py --check
Read by a gtest (the C++ writer) and by api/python/tests (the Python library); MGT v1 is
frozen when both pass every line byte-exactly.

FILE FORMAT
  UTF-8 without BOM, LF line ends; LF is the only line separator in the file (no raw CR, VT,
  FF, FS, GS, RS, NEL, LS or PS anywhere: inside the JSON they are \\u-escaped). A line that
  starts with '#' is a comment; an empty line is ignored. Every other line is one vector,
  TAB-separated:
      <tag> TAB <input> TAB <expected> [TAB <note>]
  <tag>       one of the tags below, [a-z_]+.
  <input>     one JSON text (RFC 8259); its shape depends on the tag.
  <expected>  one JSON text: a string = the exact expected token (every byte, after JSON
              unescaping), or null = decoding the input MUST fail (decode tags only).
  <note>      optional JSON string for humans; readers ignore it.
  Why JSON cells: tokens can be empty or end in a space ('ACC1' after 'ACC1' front-codes
  to "4 "), and names hold TAB, LF, CR and backslashes; a bare TSV cell would lose them to
  editors or need an escape layer of its own. The JSON contains no raw TAB, so splitting a
  line on TAB is safe. Non-ASCII is written raw (UTF-8), except that control and format
  characters, DEL, separators other than ' ' and combining marks are \\u-escaped so that
  every character is visible. Hex digits are lowercase.
  Columns are always 3 or 4; nothing else follows.

OBLIGATIONS
  encode tags (float, ranges, setexpr, pct, front, kvalue):
      encode(input) == expected, byte for byte, AND decode(expected) == input.
      (float and kvalue f: decoding gives the same bits, except that "0" decodes to +0.)
  decode tags (float_decode, ranges_decode, setexpr_decode, pct_decode, front_decode,
  kvalue_decode):
      expected null    -> decoding the input fails (an error, never a value);
      expected string  -> decoding succeeds and encode(decode(input)) == expected: a
                          valid non-canonical token dumps to its canonical form (§2.3).

TAGS
  float           input  "<16 hex digits>": IEEE-754 binary64 bits, most significant first
                         (Python: struct.unpack('<d', struct.pack('<Q', int(h, 16)))[0];
                          C++: std::bit_cast<double>(std::stoull(h, nullptr, 16))).
                  expect the float text of §2.1 (v3).
  float_decode    input  a token string.
  ranges          input  [ids...]: strictly ascending non-negative integers.
                  expect the RANGES token.
  ranges_decode   input  a token string.
  setexpr         input  [[ids...], [base...] | null]: the set and its base, both strictly
                         ascending; null = no base (the root's G.entry).
                  expect the SETEXPR token.
  setexpr_decode  input  [token, [base...] | null].
  pct             input  a string: the raw text of a free-text last field (L/X name, the
                         L suffix, Q message, K effect).
                  expect the escaped text.
  pct_decode      input  a token string.
  front           input  [previous name, name]: the previous L record's full name ("" for
                         the first L of a document) and this one's.
                  expect "<prefix_len> <suffix>": the L record's text after <seq_id|*> and
                         its separating space (suffix pct-escaped).
  front_decode    input  [previous name, token].
  kvalue          input  ["i", "<decimal>"] | ["f", "<16 hex digits>"] | ["s", "<string>"]
                         | ["u", null]: a typed K value (§2.2 v5.1) as the writer holds it.
                  expect the K VALUE token.
  kvalue_decode   input  a token string.

CODEC RULES (the design's §2.1/§2.2; "decision" marks what these vectors settle)
  Tokens: each value has exactly one valid spelling and decoders reject every other
    (decision): a decoder built on a native parser accepts language-specific extras
    (Python float() takes '1_0', ' 1', 'nan', 'Infinity'), and two decoders that accept
    different sloppy spellings disagree about which documents are valid. The format's
    only freedoms are field-level: '*' vs an explicit value, SETEXPR explicit vs delta,
    and a front-coding prefix shorter than the longest; those decode, and dump canonical.
  Integers: decimal, no sign, no leading zeros ('0' itself excepted).
  Floats: the shortest round-trip significant digits of the double (C++ std::to_chars(...,
    std::chars_format::scientific) without a precision == Python repr()), expanded
    positionally: zeros up to the decimal point, '.' only before a non-empty fraction, no
    trailing zeros; -0 -> "0"; +inf -> "inf". Domain: x >= 0 (and -0), finite or +inf;
    NaN, -inf and negatives are never encoded. Python reference: format(Decimal(repr(x)),
    'f'), fractional zeros and a bare '.' stripped, '-0' -> '0'. Decoding: the token must
    match (0|[1-9][0-9]*)(\\.[0-9]*[1-9])? or be 'inf', AND equal the encoding of the double
    it parses to (so '99999999999999991611392', the exact value of 1e23, is rejected).
  RANGES: '.' = the empty set; else comma-separated items in ascending order, every
    maximal run of >= 2 consecutive ids written 'a-b' (decision: {5,6} is '5-6', never
    '5,6'; both spellings have the same length, §2.1 says only "consecutive runs
    collapsed"), a lone id 'a'.
  SETEXPR against a base (§2.1 v2: G.entry -> parent[0]'s end set, a merge the union of
    its parents' end sets, the root none; G.end -> this segment's entry set; P -> the
    previous P of the segment, the first P the entry set; C.labels -> the leaf's end set).
    Forms: RANGES | '!' | '!R' | '!+A' | '!R+A', R = base - ids, A = ids - base, each a
    non-empty RANGES; an empty part is omitted, never written '.' (decision). The encoder
    writes the shorter (in bytes) of the explicit RANGES and the delta, a tie -> explicit;
    with no base, explicit only. Decoders reject a delta without a base, and a delta whose
    R is not inside the base or whose A meets it (decision: a mismatch means the reader
    derived another base than the writer, which must be an error, not a different set).
  pct (free-text last fields): '%' -> %25, LF -> %0A, CR -> %0D, nothing else (TAB, NUL,
    other controls, backslash, non-ASCII stay raw). Decoders reject any other '%'
    sequence, lowercase hex included, and a raw LF/CR. A record ends at LF only: names may
    hold VT, FF, FS, GS, RS, NEL, LS, PS raw, so Python must not split with str.splitlines().
  front: prefix_len = the longest common prefix of the two names' UTF-8 encodings, in
    BYTES, shortened to the nearest code-point boundary at or below it (code points, not
    grapheme clusters), so the suffix is valid UTF-8 on its own; the suffix is pct-escaped.
    The separating space is always written: an empty suffix leaves the record ending in a
    space (decision). Decoders accept a shorter prefix (valid, non-canonical) and reject a
    prefix that ends inside a character or beyond the previous name.
  kvalue: 'i:' integer (0 .. 2^64-1) | 'f:' float | 's:' string | 'u' (unlimited). 's:'
    writes every UTF-8 byte outside 0x21..0x7E, and '%' and ',', as %XX with uppercase
    hex, every other byte raw (decision: §2.2 says "percent-encoded incl. space and
    comma"; §2.1 makes every fixed field ASCII without spaces, so control and non-ASCII
    bytes are escaped too; '=' stays raw because an <extra> item splits at its first '=').
    The string 'unlimited' is 's:unlimited', distinct from 'u'.

DIFFERENTIAL SECTION
  The last {n_diff} float vectors (note-less) come from xorshift64 (Marsaglia 2003, shifts
  13/7/17) over raw bits: x = seed = 0x{seed:016x}; per draw, all mod 2^64:
      x ^= x << 13;  x ^= x >> 7;  x ^= x << 17;  bits = x & 0x7fffffffffffffff
  a draw whose exponent field (bits >> 52) is 0x7ff (inf/NaN) is skipped; the first
  {n_diff} accepted draws, in order. Expected strings by the Python reference above. When
  this file was generated, an independent C++ program (std::to_chars scientific digits
  expanded positionally, AppleClang/libc++) reproduced every float vector: 0 mismatches.

COUNTS (a reader can assert it saw every vector)
  {counts}
"""

TAG_ORDER = ['float', 'float_decode', 'ranges', 'ranges_decode', 'setexpr', 'setexpr_decode',
             'pct', 'pct_decode', 'front', 'front_decode', 'kvalue', 'kvalue_decode']

SECTION_TITLES = {
    'float': 'float: named edge cases (§2.7 v3), then 10^k and neighbours',
    'float_decode': 'float_decode: spellings that are not the canonical float text',
    'ranges': 'ranges', 'ranges_decode': 'ranges_decode',
    'setexpr': 'setexpr: each normative base, shorter wins, tie -> explicit',
    'setexpr_decode': 'setexpr_decode: valid non-canonical forms, and rejects',
    'pct': 'pct: free-text last fields', 'pct_decode': 'pct_decode',
    'front': 'front: front-coded names', 'front_decode': 'front_decode',
    'kvalue': 'kvalue: typed K values', 'kvalue_decode': 'kvalue_decode',
    'float_differential':
        'float: differential section (seeded xorshift64 over raw bits, see the header)',
}


def _visible(c):
    # controls, format characters, separators other than ' ', and combining marks are
    # \u-escaped so that every character in the file can be seen (a raw DEL or U+0301
    # is invisible) and no raw NEL/LS/PS can split a line for str.splitlines()
    cat = unicodedata.category(c)
    if c == ' ' or not (cat[0] in 'CM' or cat in ('Zl', 'Zp', 'Zs')):
        return c
    o = ord(c)
    if o <= 0xFFFF:
        return '\\u%04x' % o
    o -= 0x10000
    return '\\u%04x\\u%04x' % (0xD800 + (o >> 10), 0xDC00 + (o & 0x3FF))


def _json(v):
    s = json.dumps(v, ensure_ascii=False, separators=(',', ':'))
    # json.dumps has already escaped C0 controls, '"' and '\\'; nothing that remains
    # needs escaping outside a string, so a per-character pass is safe
    s = ''.join(_visible(c) if ord(c) >= 0x7F else c for c in s)
    assert '\t' not in s and '\n' not in s and '\r' not in s
    assert json.loads(s) == v
    return s


def render(rows):
    counts = {}
    for tag, *_ in rows:
        counts[tag] = counts.get(tag, 0) + 1
    count_line = ' '.join('%s=%d' % (t, counts[t]) for t in TAG_ORDER) + \
        ' total=%d' % len(rows)
    header = HEADER.replace('{n_diff}', str(N_DIFFERENTIAL)) \
                   .replace('{seed:016x}', '%016x' % XORSHIFT_SEED) \
                   .replace('{counts}', count_line)
    out = ['#' + (' ' + line if line else '') for line in header.rstrip('\n').split('\n')]
    section = None
    for tag, inp, exp, note in rows:
        # the differential floats are the only note-less float vectors
        key = 'float_differential' if tag == 'float' and note is None else tag
        if key != section:
            out += ['', '# ---- ' + SECTION_TITLES[key]]
            section = key
        cells = [tag, _json(inp), _json(exp)]
        if note:
            cells.append(_json(note))
        out.append('\t'.join(cells))
    return '\n'.join(out) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--out', default=DEFAULT_OUT)
    ap.add_argument('--check', action='store_true',
                    help='compare with the existing file instead of writing it')
    args = ap.parse_args()
    text = render(build())
    if args.check:
        with open(args.out, encoding='utf-8', newline='') as f:
            if f.read() != text:
                sys.exit('%s differs from the reference output' % args.out)
        print('%s: up to date' % args.out)
        return
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w', encoding='utf-8', newline='') as f:
        f.write(text)
    print('wrote %s (%d bytes)' % (args.out, len(text.encode('utf-8'))))


if __name__ == '__main__':
    main()
