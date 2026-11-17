"""Pair identity in both forms (docs/DESIGN-pair-identity.md sections 3.2 and 3.9): the n1 id
of a configured pair, tested with master's five vectors (tests/common/test_configured_pair_id.cpp
on master 356fe1ee3), and the library reading a pair named by its paths (a document stored
before the engine change) or by its id (one stored after)."""

import copy
import gzip
import json
import os
import sys
import unittest
import uuid

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

import metagraph.traverse as MT  # noqa: E402
from metagraph.traverse import CombineRefused, GraphletStore, parse, stitch  # noqa: E402
from metagraph.traverse.client import (PAIR_ID_NAMESPACE, Pair, PairIdError,  # noqa: E402
                                       answered_pair, configured_pair_id, pair_members,
                                       requested_pair)

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


# ---------------------------------------------------------------- documents in both forms

GRAPH, GP, AP = 'refseq', '/data/refseq/bact_other.dbg', '/data/refseq/bact_other.annodbg'
GP2, AP2 = '/data/refseq/bact_pseudo.dbg', '/data/refseq/bact_pseudo.annodbg'
ID = configured_pair_id(GRAPH, GP, AP)
ID2 = configured_pair_id(GRAPH, GP2, AP2)
P2 = 'p2-0123456789abcdef0123456789abcdef'


def tagged(resp, form, graph=GRAPH, gp=GP, ap=AP):
    """|resp| (a detail graphlet response) as a multi-graph server's entry states its pair:
    'paths' (before the engine named pairs by id: a document stored before), 'id' (after:
    graph and pair_id only), 'both' (paths and their own id)."""
    resp = copy.deepcopy(resp)
    for k in ('graph_path', 'annotation_path', 'pair_id'):
        resp.pop(k, None)
    resp['graph'] = graph
    if form in ('paths', 'both'):
        resp['graph_path'], resp['annotation_path'] = gp, ap
    if form in ('id', 'both'):
        resp['pair_id'] = configured_pair_id(graph, gp, ap)
    return resp


def fixture(name, form, **kw):
    resp = tagged(T.doc_json(name, 'graphlet'), form, **kw)
    return MT.from_response(resp['results'][0], resp)


R237 = os.path.join(T.DOCS, 'union', 'r237')


def r237_text(name):
    with gzip.open(os.path.join(R237, name + '.mgt.gz'), 'rt', encoding='utf-8',
                   newline='') as f:
        return f.read()


def id_form(text):
    """A stored document (H, J, body) as the service stores it after the engine change: its
    J line names the pair by graph and pair_id (the id of the paths it stated), no path."""
    lines = text.split('\n')
    i = next(i for i, line in enumerate(lines) if line.startswith('J '))
    j = json.loads(lines[i][2:])
    gp, ap = j.pop('graph_path'), j.pop('annotation_path')
    j['pair_id'] = configured_pair_id(j['graph'], gp, ap)
    lines[i] = 'J ' + json.dumps(j, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    return '\n'.join(lines)


class TestReadingBothForms(unittest.TestCase):
    """docs/DESIGN-pair-identity.md section 3.9: a document's pair id is its pair_id, or the
    n1 id of its paths; the library reads both forms wherever it reads pair tags."""

    def test_an_old_document_and_its_id_form_twin_are_one_pair(self):
        text = r237_text('bact_other')
        old, new = parse(text), parse(id_form(text))
        gp, ap = old.envelope['graph_path'], old.envelope['annotation_path']
        want = configured_pair_id('bact_other', gp, ap)
        self.assertEqual(answered_pair(old.envelope), Pair('bact_other', gp, ap, want))
        self.assertEqual(answered_pair(new.envelope), Pair('bact_other', None, None, want))
        self.assertEqual(answered_pair(old.envelope).identity,
                         answered_pair(new.envelope).identity)
        self.assertNotIn('graph_path', new.envelope)
        # the body is the same: the body digest leaves out the J line
        self.assertEqual(old.dump(envelope=False), new.dump(envelope=False))

    def test_a_single_graph_answer_names_no_pair(self):
        self.assertIsNone(answered_pair(T.graphlet('fork').envelope))
        self.assertIsNone(answered_pair({'graph': None, 'pair_id': None}))

    def test_ids_the_library_does_not_take(self):
        env = {'graph': GRAPH}
        for bad, code in ((P2, 'p2_pair_id'), ('n1-NOT-AN-ID', 'bad_pair_id'),
                          (ID.upper().replace('N1-', 'n1-'), 'bad_pair_id'),
                          ('n1-' + str(uuid.uuid4()), 'bad_pair_id'), (17, 'bad_pair_id')):
            with self.assertRaises(PairIdError) as cm:
                answered_pair(dict(env, pair_id=bad))
            self.assertEqual(cm.exception.code, code, bad)
            self.assertIsInstance(cm.exception, ValueError)
        with self.assertRaises(PairIdError) as cm:
            answered_pair(dict(env, pair_id=P2))
        self.assertIn("is the service's public pair id (p2), not the engine's n1 id",
                      str(cm.exception))
        # an id beside both paths must be theirs; the error names no path
        with self.assertRaises(PairIdError) as cm:
            answered_pair(dict(env, graph_path=GP, annotation_path=AP, pair_id=ID2))
        self.assertEqual(cm.exception.code, 'pair_id_mismatch')
        self.assertNotIn(GP, str(cm.exception))
        self.assertEqual(answered_pair(dict(env, graph_path=GP, annotation_path=AP, pair_id=ID)),
                         Pair(GRAPH, GP, AP, ID))

    def test_selectors_in_both_forms(self):
        self.assertEqual(pair_members(GRAPH, pair_id=ID),
                         {'pairs': [{'graph': GRAPH, 'pair_id': ID}]})
        self.assertEqual(requested_pair(pair_members(GRAPH, pair_id=ID)),
                         Pair(GRAPH, None, None, ID))
        # an old document's Pair (paths and their id) is written by its paths
        self.assertEqual(pair_members(*Pair(GRAPH, GP, AP, ID)),
                         {'pairs': [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP}]})
        self.assertEqual(requested_pair(pair_members(GRAPH, GP, AP)), Pair(GRAPH, GP, AP, ID))
        self.assertEqual(requested_pair({'graph': GRAPH, 'graph_path': GP}), Pair(GRAPH, GP))
        with self.assertRaises(PairIdError) as cm:
            pair_members(GRAPH, GP, None, ID)
        self.assertEqual(cm.exception.code, 'bad_pair_id')
        with self.assertRaises(PairIdError) as cm:
            pair_members(GRAPH, GP, AP, ID2)
        self.assertEqual(cm.exception.code, 'pair_id_mismatch')
        with self.assertRaises(PairIdError) as cm:
            pair_members(GRAPH, pair_id=P2)
        self.assertEqual(cm.exception.code, 'p2_pair_id')
        with self.assertRaises(ValueError) as cm:
            pair_members(None, pair_id=ID)
        self.assertNotIsInstance(cm.exception, PairIdError)
        with self.assertRaises(PairIdError):
            requested_pair({'pairs': [{'graph': GRAPH, 'pair_id': P2}]})

    def test_replay_and_continuation_of_both_forms(self):
        old, new = fixture('fork', 'paths'), fixture('fork', 'id')
        by_paths = {'pairs': [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP}]}
        by_id = {'pairs': [{'graph': GRAPH, 'pair_id': ID}]}
        a, b = GraphletStore.request_of(old), GraphletStore.request_of(new)
        self.assertEqual(a['pairs'], by_paths['pairs'])
        self.assertEqual(b['pairs'], by_id['pairs'])
        self.assertEqual(dict(a, pairs=None), dict(b, pairs=None))
        leaf = [old.walks('right')[0].path_id]
        a, b = old.next_request('right', leaf), new.next_request('right', leaf)
        self.assertEqual((a['pairs'], b['pairs']), (by_paths['pairs'], by_id['pairs']))
        self.assertEqual(dict(a, pairs=None), dict(b, pairs=None))
        # a selector override replaces the envelope's pair whole, by id too
        c = old.next_request('right', leaf, graph=GRAPH, pair_id=ID2)
        self.assertEqual(c['pairs'], [{'graph': GRAPH, 'pair_id': ID2}])
        with self.assertRaises(ValueError):
            old.next_request('right', leaf, graphs=[GRAPH], pair_id=ID2)
        # a stored request of the id form reads back as its pair
        self.assertEqual(requested_pair(GraphletStore.request_of(new)),
                         Pair(GRAPH, None, None, ID))

    def test_stitching_an_old_parent_and_an_id_form_child(self):
        def load(role, form, **kw):
            with open(os.path.join(T.STITCH, 'keep_forbid.%s.json' % role),
                      encoding='utf-8') as f:
                resp = tagged(json.load(f), form, **kw)
            return MT.from_response(resp['results'][0], resp)
        # one pair in two forms: one stitch identity
        sg = stitch(load('parent', 'paths'), load('child', 'id'))
        self.assertEqual(len(sg.parts), 2)
        stitch(load('parent', 'id'), load('child', 'paths'))
        stitch(load('parent', 'both'), load('child', 'id'))
        # another pair of the same name, in either form: different_pair
        for pf, cf in (('paths', 'id'), ('id', 'paths'), ('id', 'id')):
            with self.assertRaises(CombineRefused) as cm:
                stitch(load('parent', pf), load('child', cf, gp=GP2, ap=AP2))
            self.assertEqual(cm.exception.code, 'different_pair', (pf, cf))


class _FanOut(MT.TraverseClient):
    """A client whose GET /capabilities is |caps| and whose traverse() records the pair it
    was asked to run on."""

    def __init__(self, caps, **kw):
        super().__init__('localhost', 1, **kw)
        self.caps, self.sent = caps, []

    def server_capabilities(self):
        return self.caps

    def capabilities(self, *a, **kw):
        return {}

    def traverse(self, seeds, strategy=None, **kw):
        self.sent.append({k: v for k, v in kw.items()
                          if k in Pair._fields and v is not None})
        return 'answer'


class TestDiscoveryInBothForms(unittest.TestCase):
    """traverse_pairs() on GET /capabilities of an engine that lists pairs by id (one name,
    several pairs: production RefSeq's form) and of one that lists their paths."""

    BY_ID = {'graph_summary': {'pairs': [
        {'graph': GRAPH, 'pair_id': ID, 'seqs_order': 'verified'},
        {'graph': GRAPH, 'pair_id': ID2, 'seqs_order': 'count_only'}]},
        'multigraph': {'pair_id': 'n1'}}
    BY_PATHS = {'graph_summary': {'pairs': [
        {'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP, 'seqs_order': 'verified'},
        {'graph': GRAPH, 'graph_path': GP2, 'annotation_path': AP2}]}}

    def test_every_pair_an_engine_lists_by_id(self):
        c = _FanOut(self.BY_ID)
        r = c.traverse_pairs('ACGT', coordinates=None)
        self.assertEqual(c.sent, [{'graph': GRAPH, 'pair_id': ID},
                                  {'graph': GRAPH, 'pair_id': ID2}])
        self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': ID},
                                            {'graph': GRAPH, 'pair_id': ID2}))
        self.assertEqual(r.seqs_order, {(GRAPH, ID): 'verified', (GRAPH, ID2): 'count_only'})
        # one of them, by its id
        c = _FanOut(self.BY_ID)
        r = c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': ID2}], coordinates=None)
        self.assertEqual((c.sent, r.seqs_order), ([{'graph': GRAPH, 'pair_id': ID2}],
                                                  {(GRAPH, ID2): 'count_only'}))

    def test_a_pair_named_by_id_on_an_engine_that_lists_paths(self):
        # sent by the paths the engine lists (its own spelling), whose id is the one given
        c = _FanOut(self.BY_PATHS)
        r = c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': ID}], coordinates=None)
        self.assertEqual(c.sent, [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP}])
        self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': ID},))
        self.assertEqual(r.seqs_order, {(GRAPH, ID): 'verified'})
        # an id the listing does not have is sent as it is (the client decides, below)
        other = configured_pair_id(GRAPH, '/elsewhere.dbg', AP)
        c = _FanOut(self.BY_PATHS)
        r = c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': other}],
                             coordinates=None)
        self.assertEqual((c.sent, r.seqs_order), ([{'graph': GRAPH, 'pair_id': other}], {}))

    def test_a_p2_id_is_refused_before_anything_is_sent(self):
        c = _FanOut(self.BY_ID)
        with self.assertRaises(PairIdError) as cm:
            c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': P2}], coordinates=None)
        self.assertEqual((cm.exception.code, c.sent), ('p2_pair_id', []))
        bad = copy.deepcopy(self.BY_ID)
        bad['graph_summary']['pairs'][1]['pair_id'] = P2
        with self.assertRaises(PairIdError):
            _FanOut(bad).traverse_pairs('ACGT', coordinates=None)


if __name__ == '__main__':
    unittest.main()
