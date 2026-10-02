"""MGT v1 codecs against the shared golden vectors (DESIGN-traverse-graphlet.md §2.7).

codec_vectors.tsv is the freeze gate both the C++ writer and this library must pass
byte-exactly; its header (OBLIGATIONS) is the spec this test implements.
"""

import json
import os
import random
import re
import struct
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from metagraph.traverse import _codec as C  # noqa: E402

VECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'traverse',
                       'codec_vectors.tsv')


def bits_to_float(h):
    return struct.unpack('<d', struct.pack('<Q', int(h, 16)))[0]


def float_to_bits(x):
    return '%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]


def load_vectors():
    counts = None
    rows = []
    with open(VECTORS, 'rb') as f:
        data = f.read()
    text = data.decode('utf-8')
    # LF is the only line separator (names may hold VT, FF, NEL, LS, PS raw)
    for line in text.split('\n'):
        if line.startswith('#'):
            m = re.match(r'#\s+(float=\d+.*total=\d+)\s*$', line)
            if m:
                counts = dict(kv.split('=') for kv in m.group(1).split())
            continue
        if not line:
            continue
        cells = line.split('\t')
        assert len(cells) in (3, 4), line
        rows.append((cells[0], json.loads(cells[1]), json.loads(cells[2])))
    return rows, {k: int(v) for k, v in counts.items()}


def kvalue_of(spec):
    kind, value = spec
    if kind == 'i':
        return int(value)
    if kind == 'f':
        return bits_to_float(value)
    if kind == 's':
        return value
    assert kind == 'u' and value is None
    return C.UNLIMITED


def same_kvalue(a, b):
    if a is C.UNLIMITED or b is C.UNLIMITED:
        return a is b
    if type(a) is not type(b):
        return False
    if isinstance(a, float):
        # the -0 vectors decode "0" to +0
        return float_to_bits(a) == float_to_bits(b) or (a == 0 and b == 0)
    return a == b


class TestGoldenVectors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.counts = load_vectors()

    def test_every_vector_is_present(self):
        seen = {}
        for tag, _, _ in self.rows:
            seen[tag] = seen.get(tag, 0) + 1
        expected = dict(self.counts)
        total = expected.pop('total')
        self.assertEqual(expected, seen)
        self.assertEqual(total, len(self.rows))

    def _check(self, tag, encode, decode, same=lambda a, b: a == b):
        n = 0
        for t, inp, expected in self.rows:
            if t == tag:
                n += 1
                got = encode(inp)
                self.assertEqual(expected, got, (tag, inp))
                self.assertTrue(same(decode(expected), inp), (tag, inp, decode(expected)))
            elif t == tag + '_decode':
                n += 1
                if expected is None:
                    with self.assertRaises(C.CodecError, msg=(t, inp)):
                        decode(inp)
                else:
                    self.assertEqual(expected, encode(decode(inp)), (t, inp))
        self.assertGreater(n, 0)

    def test_float(self):
        def same(x, h):
            return float_to_bits(x) == h or (x == 0 and bits_to_float(h) == 0)
        self._check('float', lambda h: C.fmt_num(bits_to_float(h)), C.parse_num, same)
        # float_decode inputs are tokens, encode() takes a value
        for t, inp, expected in self.rows:
            if t == 'float_decode' and expected is not None:
                self.assertEqual(expected, C.fmt_num(C.parse_num(inp)))

    def test_ranges(self):
        self._check('ranges', C.encode_ranges, C.decode_ranges)

    def test_setexpr(self):
        n = 0
        for t, inp, expected in self.rows:
            if t == 'setexpr':
                ids, base = inp
                self.assertEqual(expected, C.encode_setexpr(ids, base), inp)
                self.assertEqual(ids, C.decode_setexpr(expected, base), inp)
                n += 1
            elif t == 'setexpr_decode':
                tok, base = inp
                if expected is None:
                    with self.assertRaises(C.CodecError, msg=inp):
                        C.decode_setexpr(tok, base)
                else:
                    self.assertEqual(expected,
                                     C.encode_setexpr(C.decode_setexpr(tok, base), base), inp)
                n += 1
        self.assertEqual(self.counts['setexpr'] + self.counts['setexpr_decode'], n)

    def test_pct(self):
        self._check('pct', C.pct_escape, C.pct_unescape)

    def test_front(self):
        n = 0
        for t, inp, expected in self.rows:
            if t == 'front':
                prev, name = inp
                self.assertEqual(expected, C.front_encode(prev, name), inp)
                self.assertEqual(name, C.front_decode(prev, expected), inp)
                n += 1
            elif t == 'front_decode':
                prev, tok = inp
                if expected is None:
                    with self.assertRaises(C.CodecError, msg=inp):
                        C.front_decode(prev, tok)
                else:
                    self.assertEqual(expected, C.front_encode(prev, C.front_decode(prev, tok)),
                                     inp)
                n += 1
        self.assertEqual(self.counts['front'] + self.counts['front_decode'], n)

    def test_kvalue(self):
        n = 0
        for t, inp, expected in self.rows:
            if t == 'kvalue':
                v = kvalue_of(inp)
                self.assertEqual(expected, C.encode_kvalue(v), inp)
                self.assertTrue(same_kvalue(C.decode_kvalue(expected), v), inp)
                n += 1
            elif t == 'kvalue_decode':
                if expected is None:
                    with self.assertRaises(C.CodecError, msg=inp):
                        C.decode_kvalue(inp)
                else:
                    self.assertEqual(expected, C.encode_kvalue(C.decode_kvalue(inp)), inp)
                n += 1
        self.assertEqual(self.counts['kvalue'] + self.counts['kvalue_decode'], n)


class TestCodecIdentities(unittest.TestCase):
    """Random identities beyond the fixed vectors (§9: random sets, shorter-wins,
    tie -> explicit)."""

    def test_random_sets_round_trip(self):
        rng = random.Random(7)
        for _ in range(3000):
            universe = rng.randint(0, 60)
            ids = sorted(rng.sample(range(universe), rng.randint(0, universe)))
            base = sorted(rng.sample(range(universe + 5), rng.randint(0, universe)))
            tok = C.encode_ranges(ids)
            self.assertEqual(ids, C.decode_ranges(tok))
            for b in (None, base):
                expr = C.encode_setexpr(ids, b)
                self.assertEqual(ids, C.decode_setexpr(expr, b))
                if b is not None:
                    explicit = C.encode_ranges(ids)
                    if expr != explicit:
                        # a delta is chosen only when strictly shorter
                        self.assertLess(len(expr), len(explicit))
                        self.assertTrue(expr.startswith('!'))

    def test_random_floats_round_trip(self):
        rng = random.Random(11)
        for _ in range(3000):
            x = abs(rng.uniform(0, 1e6) * 10 ** rng.randint(-30, 30))
            tok = C.fmt_num(x)
            self.assertEqual(x, C.parse_num(tok))

    def test_end_tokens(self):
        for tok in ['D', 'Rm', 'Rb', 'Rs', 'Dh', 'Ls', 'Lx', 'Lw', 'm', 'X', 'Y', 'W']:
            code, qual, silent, merged = C.parse_end_token(tok)
            if merged:
                self.assertEqual('m', tok)
            elif silent:
                self.assertEqual('Lw', tok)
            else:
                self.assertEqual(tok, C.end_token(code, qual))
        for bad in ['', 'Q', 'Rx', 'DD', 'w', 'M ']:
            with self.assertRaises(C.CodecError):
                C.parse_end_token(bad)

    def test_interner_shares_equal_sets(self):
        it = C.LabelSetInterner()
        a = it.intern([1, 2, 3])
        b = it.intern([1, 2, 3])
        self.assertIs(a, b)
        self.assertEqual(1, it.misses)
        self.assertEqual(1, it.hits)
        self.assertEqual([], list(it.intern([])))


if __name__ == '__main__':
    unittest.main()
