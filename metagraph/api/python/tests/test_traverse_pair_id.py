"""Pair identity in both forms (docs/DESIGN-pair-identity.md sections 3.2 and 3.9): the n1 id
of a configured pair, tested with master's five vectors (tests/common/test_configured_pair_id.cpp
on master 356fe1ee3), and the library reading a pair named by its paths (a document stored
before the engine change) or by its id (one stored after)."""

import os
import sys
import unittest
import uuid

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from metagraph.traverse.client import PAIR_ID_NAMESPACE, configured_pair_id  # noqa: E402

# master's tests/common/test_configured_pair_id.cpp, PythonUuid5Vectors, byte for byte
MASTER_VECTORS = [
    (('refseq', '/graphs/a.dbg', '/annotations/a.annodbg'),
     'n1-0c8ac32d-6ec3-5dc6-85c0-6993b9c02020'),
    (('refseq', '/graphs/a.dbg', '/annotations/b.annodbg'),
     'n1-f3db2c00-cc7c-5044-a87d-a679fd98328e'),
    (('réf', '/图/a', ''), 'n1-b226283c-5004-52d2-a5f4-ccc494acf745'),
    (('refseq', None, None), 'n1-d600c46d-30c3-5f3d-a520-7d231fcb1141'),
    (('g:12', '12:/a:b', '/b:34'), 'n1-5ac59ea8-6d4a-5b01-aa08-e763a91371d2'),
]


class TestConfiguredPairId(unittest.TestCase):
    def test_master_vectors(self):
        for fields, want in MASTER_VECTORS:
            self.assertEqual(configured_pair_id(*fields), want, fields)

    def test_python_uuid5(self):
        # the same id as uuid.uuid5 over the length-prefixed name (lengths in UTF-8 bytes)
        for fields, want in MASTER_VECTORS:
            name = ''.join('-1:' if v is None else '%d:%s' % (len(v.encode('utf-8')), v)
                           for v in fields)
            self.assertEqual('n1-' + str(uuid.uuid5(PAIR_ID_NAMESPACE, name)), want)

    def test_fields_and_absence_differ(self):
        # master's DistinguishesFieldsAndAbsence
        self.assertNotEqual(configured_pair_id('g', 'ab', 'c'), configured_pair_id('g', 'a', 'bc'))
        self.assertNotEqual(configured_pair_id('g', 'a', 'b'), configured_pair_id('h', 'a', 'b'))
        self.assertNotEqual(configured_pair_id('g', 'a', 'b'), configured_pair_id('g', 'c', 'b'))
        self.assertNotEqual(configured_pair_id('g', '', 'b'), configured_pair_id('g', None, 'b'))

    def test_nothing_is_normalized(self):
        a = configured_pair_id('refseq', '/graphs/a.dbg', '/annotations/a.annodbg')
        for spelled in ('/graphs//a.dbg', '/graphs/a.dbg/', 'graphs/a.dbg', '/graphs/./a.dbg'):
            self.assertNotEqual(configured_pair_id('refseq', spelled, '/annotations/a.annodbg'),
                                a, spelled)
        # NFC and NFD spell one name two ways: two ids
        self.assertNotEqual(configured_pair_id('réf', 'g', 'a'),
                            configured_pair_id('réf', 'g', 'a'))

    def test_a_field_is_a_string_or_none(self):
        for bad in (b'/graphs/a.dbg', 1, ['a']):
            with self.assertRaises(TypeError):
                configured_pair_id('refseq', bad, None)


if __name__ == '__main__':
    unittest.main()
