"""Pair identity in both forms (docs/DESIGN-pair-identity.md sections 3.2 and 3.9): the n1 id
of a configured pair, tested with master's five vectors (tests/common/test_configured_pair_id.cpp
on master 356fe1ee3), and the library reading a pair named by its paths (a document stored
before the engine change) or by its id (one stored after)."""

import copy
import gzip
import json
import os
import sys
import tempfile
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
        # an id the listing does not have: sent by id where the client's resolver knows it
        # (traverse() resolves it), else refused before any pair's request is sent -- also
        # the listed pair before it
        other = configured_pair_id(GRAPH, '/elsewhere.dbg', AP)
        c = _FanOut(self.BY_PATHS,
                    pair_resolver={other: (GRAPH, '/elsewhere.dbg', AP)}.get)
        r = c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': other}],
                             coordinates=None)
        self.assertEqual((c.sent, r.seqs_order), ([{'graph': GRAPH, 'pair_id': other}], {}))
        c = _FanOut(self.BY_PATHS)
        with self.assertRaises(PairIdError) as cm:
            c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': ID},
                                            {'graph': GRAPH, 'pair_id': other}],
                             coordinates=None)
        self.assertEqual((cm.exception.code, c.sent), ('pair_unresolved', []))

    def test_a_p2_id_is_refused_before_anything_is_sent(self):
        c = _FanOut(self.BY_ID)
        with self.assertRaises(PairIdError) as cm:
            c.traverse_pairs('ACGT', pairs=[{'graph': GRAPH, 'pair_id': P2}], coordinates=None)
        self.assertEqual((cm.exception.code, c.sent), ('p2_pair_id', []))
        bad = copy.deepcopy(self.BY_ID)
        bad['graph_summary']['pairs'][1]['pair_id'] = P2
        with self.assertRaises(PairIdError):
            _FanOut(bad).traverse_pairs('ACGT', coordinates=None)


class _Engine(MT.TraverseClient):
    """A client whose server is |caps| (GET /capabilities), an empty probe answer and
    |answer| for every POST /traverse; |log| records what was sent."""

    def __init__(self, caps, answer=None, **kw):
        super().__init__('localhost', 1, **kw)
        self.caps, self.answer, self.log = caps, answer, []

    def _request(self, method, path, payload=None, answers=(), with_status=False):
        self.log.append((method, path, copy.deepcopy(payload)))
        if path == '/capabilities':
            out = copy.deepcopy(self.caps)
        elif path.startswith('/traverse/capabilities'):
            out = copy.deepcopy((self.answer or {}).get('capabilities') or {'index_ns': 'ns'})
        else:
            out = copy.deepcopy(self.answer)
        return (200, out) if with_status else out

    def posted(self):
        return [p for m, path, p in self.log if (m, path) == ('POST', '/traverse')]

    def probes(self):
        return [path for _, path, _ in self.log if path.startswith('/traverse/capabilities')]


BEFORE_IDS = TestDiscoveryInBothForms.BY_PATHS             # an engine before M6: paths only
AFTER_IDS = TestDiscoveryInBothForms.BY_ID                 # after: multigraph.pair_id "n1"
# an engine before M6 that does not serve the pair (graph, GP, AP): it lists GP2's only
UNLISTED = {'graph_summary': {'pairs': [BEFORE_IDS['graph_summary']['pairs'][1]]}}
BY_PATHS = [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP}]
BY_ID = [{'graph': GRAPH, 'pair_id': ID}]


class TestTheResolverHook(unittest.TestCase):
    """A continuation, replay or re-ask of a document that names its pair by id only reaches
    an engine through TraverseClient.wire_request(): by the caller's resolver (paths sent,
    n1 recomputed); where it knows no such id, by an engine that states
    multigraph.pair_id "n1" (the id sent); else by the paths that engine's own GET
    /capabilities lists (n1 recomputed); refused where none of these names it. An old
    document is sent by its paths whatever the engine."""

    def setUp(self):
        self.old, self.new = fixture('fork', 'paths'), fixture('fork', 'id')
        self.answer = tagged(T.doc_json('fork', 'graphlet'), 'id')
        self.leaf = [self.old.walks('right')[0].path_id]

    def engines(self):
        resolver = MT.client.capabilities_pair_resolver(BEFORE_IDS)
        # the service's case once engines run M6: its resolver built from the engine's own
        # listing, which names pairs by id (review M1)
        m6 = MT.client.capabilities_pair_resolver(AFTER_IDS)
        return {'resolver': _Engine(BEFORE_IDS, self.answer, pair_resolver=resolver),
                'resolver_after': _Engine(AFTER_IDS, self.answer, pair_resolver=resolver),
                'm6_resolver': _Engine(AFTER_IDS, self.answer, pair_resolver=m6),
                'after_ids': _Engine(AFTER_IDS, self.answer),
                'listing': _Engine(BEFORE_IDS, self.answer),
                'neither': _Engine(UNLISTED, self.answer),
                'neither_resolver': _Engine(UNLISTED, self.answer, pair_resolver={}.get)}

    def ask(self, how, c):
        """One re-ask of the id-form document |how|: deepen(), a replay of request_of()
        through traverse_raw(), traverse() of its pair, or the probe."""
        if how == 'deepen':
            c.deepen(self.new, 'right', self.leaf)
        elif how == 'replay':
            c.traverse_raw(GraphletStore.request_of(self.new))
        elif how == 'traverse':
            c.traverse(['ACGTACGTACGTACGTACGTACGTACGTACGT'], graph=GRAPH, pair_id=ID,
                       coordinates=None)
        else:
            c.capabilities(GRAPH, pair_id=ID)

    def test_each_reask_by_resolver_engine_or_refused(self):
        want = {'resolver': BY_PATHS, 'resolver_after': BY_PATHS, 'm6_resolver': BY_ID,
                'after_ids': BY_ID, 'listing': BY_PATHS}
        probe = {'resolver': '/traverse/capabilities?graph=refseq&graph_path=%2Fdata%2Frefseq'
                             '%2Fbact_other.dbg&annotation_path=%2Fdata%2Frefseq%2F'
                             'bact_other.annodbg',
                 'after_ids': '/traverse/capabilities?graph=refseq&pair_id=' + ID}
        probe['resolver_after'] = probe['listing'] = probe['resolver']
        probe['m6_resolver'] = probe['after_ids']
        for how in ('deepen', 'replay', 'traverse', 'probe'):
            for name, c in self.engines().items():
                if name.startswith('neither'):
                    with self.assertRaises(PairIdError) as cm:
                        self.ask(how, c)
                    self.assertEqual(cm.exception.code, 'pair_unresolved', how)
                    self.assertIsNone(getattr(cm.exception, 'sent', None), how)
                    self.assertIn('give the client a pair_resolver', str(cm.exception))
                    self.assertEqual((c.posted(), c.probes()), ([], []), how)
                    continue
                self.ask(how, c)
                if how == 'probe':
                    self.assertEqual(c.probes(), [probe[name]], (how, name))
                else:
                    (req,) = c.posted()
                    self.assertEqual(req['pairs'], want[name], (how, name))
                    self.assertNotIn('graph', req)

    def test_an_old_document_is_sent_by_its_paths_on_every_engine(self):
        for name, c in self.engines().items():
            c.deepen(self.old, 'right', self.leaf)
            c.traverse_raw(GraphletStore.request_of(self.old))
            self.assertEqual([r['pairs'] for r in c.posted()], [BY_PATHS] * 2, name)
            # nothing about ids is read: no GET /capabilities
            self.assertNotIn('/capabilities', [path for _, path, _ in c.log], name)

    def test_the_request_kept_names_the_pair_by_id(self):
        # the request a store keeps (and a NextRequest) is not changed by sending it
        c = self.engines()['resolver']
        req = GraphletStore.request_of(self.new)
        c.traverse_raw(req)
        self.assertEqual(req['pairs'], BY_ID)
        self.assertEqual(c.posted()[0]['pairs'], BY_PATHS)

    def test_a_resolver_that_answers_another_pair(self):
        bad = {ID: (GRAPH, GP2, AP2)}.get                 # spelled for the other pair
        renamed = {ID: ('other', GP, AP)}.get            # the paths under another name
        for resolver, code in ((bad, 'pair_id_mismatch'), (renamed, 'pair_id_mismatch'),
                               (lambda pid: 'x', 'pair_unresolved'),
                               ({ID: (GRAPH, None, AP)}.get, 'pair_unresolved')):
            for caps in (AFTER_IDS, BEFORE_IDS):
                c = _Engine(caps, self.answer, pair_resolver=resolver)
                with self.assertRaises(PairIdError) as cm:
                    c.deepen(self.new, 'right', self.leaf)
                self.assertEqual(cm.exception.code, code)
                self.assertEqual(c.posted(), [])
                for path in (GP, GP2, AP, AP2):
                    self.assertNotIn(path, str(cm.exception))
        # a resolver that does not know the id is not the last word (review M1): the engine
        # that states ids takes the id, one before ids is asked by the paths it lists
        for caps, sent in ((AFTER_IDS, BY_ID), (BEFORE_IDS, BY_PATHS)):
            c = _Engine(caps, self.answer, pair_resolver={}.get)
            c.deepen(self.new, 'right', self.leaf)
            self.assertEqual(c.posted()[0]['pairs'], sent)

    def test_the_listing_is_read_once_and_goes_back_to_its_engine(self):
        # the default fallback: the paths of an engine's own GET /capabilities, read once per
        # server process (with states_pair_ids()), sent only to that engine
        c = _Engine(BEFORE_IDS, self.answer)
        for _ in range(3):
            c.deepen(self.new, 'right', self.leaf)
        self.assertEqual([r['pairs'] for r in c.posted()], [BY_PATHS] * 3)
        self.assertEqual([path for _, path, _ in c.log].count('/capabilities'), 1)
        other = _Engine(UNLISTED, self.answer)
        with self.assertRaises(PairIdError) as cm:
            other.deepen(self.new, 'right', self.leaf)
        self.assertEqual((cm.exception.code, other.posted()), ('pair_unresolved', []))
        # a listing whose entry states an id that is not its paths' is refused, not sent
        lying = copy.deepcopy(BEFORE_IDS)
        lying['graph_summary']['pairs'][0]['pair_id'] = ID2
        c = _Engine(lying, self.answer)
        with self.assertRaises(PairIdError) as cm:
            c.deepen(self.new, 'right', self.leaf)
        # the entry the library does not take resolves nothing (fix review N3)
        self.assertEqual((cm.exception.code, c.posted()), ('pair_unresolved', []))

    def test_a_bad_listing_entry_leaves_the_others(self):
        # fix review N3: one malformed entry (a p2 id, an id that is not its paths') is
        # left out of the listing's resolver; the others still resolve, also as fallback
        caps = copy.deepcopy(BEFORE_IDS)
        caps['graph_summary']['pairs'] += [
            {'graph': 'other', 'graph_path': '/x.dbg', 'pair_id': P2},
            {'graph': GRAPH, 'graph_path': GP2, 'annotation_path': AP2, 'pair_id': ID}]
        r = MT.client.capabilities_pair_resolver(caps)
        self.assertEqual(r(ID), (GRAPH, GP, AP))
        c = _Engine(caps, self.answer)
        c.deepen(self.new, 'right', self.leaf)
        self.assertEqual(c.posted()[0]['pairs'], BY_PATHS)

    def test_selectors_the_wire_refuses(self):
        c = _Engine(AFTER_IDS, self.answer)
        for pairs, code in (([{'graph': GRAPH, 'pair_id': ID, 'graph_path': GP}], 'bad_pair_id'),
                            ([{'graph': GRAPH, 'pair_id': P2}], 'p2_pair_id'),
                            ([{'pair_id': ID}], 'bad_pair_id'),
                            ([{'graph': GRAPH, 'pair_id': ID2, 'graph_path': GP,
                               'annotation_path': AP}], 'pair_id_mismatch')):
            with self.assertRaises(PairIdError) as cm:
                c.traverse_raw({'seeds': [], 'pairs': pairs})
            self.assertEqual(cm.exception.code, code, pairs)
        self.assertEqual(c.posted(), [])
        # both paths with their own id: sent by the paths
        c.traverse_raw({'seeds': [], 'pairs': [dict(BY_PATHS[0], pair_id=ID)]})
        self.assertEqual(c.posted()[0]['pairs'], BY_PATHS)

    def test_the_capabilities_resolver(self):
        r = MT.client.capabilities_pair_resolver(BEFORE_IDS)
        self.assertEqual((r(ID), r(ID2)), ((GRAPH, GP, AP), (GRAPH, GP2, AP2)))
        self.assertIsNone(r(configured_pair_id(GRAPH, GP, AP2)))
        # an engine after ids lists no paths: nothing to resolve to
        self.assertIsNone(MT.client.capabilities_pair_resolver(AFTER_IDS)(ID))
        self.assertIsNone(MT.client.capabilities_pair_resolver(None)(ID))


class TestToolsInBothForms(unittest.TestCase):
    """GraphletTools: a loaded id-form document is replayed and continued through its
    client's hook; a refusal is answered with the PairIdError's code."""

    def test_replay_and_continue_an_id_form_entry(self):
        from metagraph.traverse.mcp_tools import GraphletTools
        base = TestTheResolverHook('setUp')
        base.setUp()
        want = {'resolver': BY_PATHS, 'resolver_after': BY_PATHS, 'm6_resolver': BY_ID,
                'after_ids': BY_ID, 'listing': BY_PATHS}
        for name, client in base.engines().items():
            with tempfile.TemporaryDirectory() as d:
                store = GraphletStore(d)
                tools = GraphletTools(store, {'mini': client})
                path = os.path.join(d, 'new.mgt')
                base.new.save(path)
                h = tools.graphlet_load(path)['handle']
                self.assertEqual(store.get(h).request['pairs'], BY_ID)
                fetched = tools.traverse_fetch(replay=h, index='mini')
                cont = tools.traverse_continue(h, 'right', base.leaf[0])
                if name.startswith('neither'):
                    for out in (fetched, cont):
                        self.assertEqual(out['error'], 'pair_unresolved', out)
                    self.assertEqual((client.posted(), client.probes()), ([], []))
                    continue
                self.assertNotIn('error', fetched, (name, fetched))
                self.assertNotIn('error', cont, (name, cont))
                self.assertEqual([r['pairs'] for r in client.posted()], [want[name]] * 2, name)
                # the entries keep the id form; the probe named the same pair
                self.assertEqual(store.get(fetched['handle']).request['pairs'], BY_ID)
                self.assertEqual(len(client.probes()), 2)

    def test_a_p2_id_is_its_own_error(self):
        from metagraph.traverse.mcp_tools import GraphletTools
        with tempfile.TemporaryDirectory() as d:
            client = _Engine(AFTER_IDS, tagged(T.doc_json('fork', 'graphlet'), 'id'))
            tools = GraphletTools(GraphletStore(d), {'mini': client})
            out = tools.traverse_continue('g_000000000000', 'right', 0)
            self.assertEqual(out['error'], 'unknown_handle')
            h = tools.traverse_fetch(seed={'sequence': 'ACGT'}, index='mini')['handle']
            out = tools.traverse_continue(h, 'right', 0, overrides={'graph': GRAPH,
                                                                    'pair_id': P2})
            self.assertEqual(out['error'], 'p2_pair_id', out)


# ---------------------------------------------------------------- the union of both forms

R237_PARTS = ('bact_other', 'bact_pseudomonadota')


def r237_inputs(parts='paths', absent='paths', expected='paths', graph=None):
    """The r237 union's inputs (test_traverse_union.r237_union()) with each group in the form
    asked: 'paths' as the service delivered them, 'id' as it stores them after the engine
    change (graph and pair_id, no path), or a list of forms per part. |graph|: every pair
    renamed to one name (production RefSeq's shape: one name, several pairs)."""
    from metagraph.traverse import union_pairs  # noqa: F401
    forms = [parts] * 2 if isinstance(parts, str) else parts
    texts = [r237_text(n) for n in R237_PARTS]
    if graph is not None:
        texts = [t.replace('"graph":"%s"' % n, '"graph":"%s"' % graph)
                 for t, n in zip(texts, R237_PARTS)]
    gs = [parse(id_form(t) if f == 'id' else t) for t, f in zip(texts, forms)]
    with open(os.path.join(R237, 'absent.json'), encoding='utf-8') as f:
        recs = json.load(f)
    with open(os.path.join(R237, 'seed.fa'), encoding='utf-8') as f:
        seed = ''.join(line.strip() for line in f if not line.startswith('>'))
    exp = [{k: g0.envelope[k] for k in ('graph', 'graph_path', 'annotation_path')}
           for g0 in (parse(t) for t in texts)] + [dict(a['pair']) for a in recs]
    if graph is not None:
        for r in recs:
            r['pair']['graph'] = graph
        for d in exp:
            d['graph'] = graph

    def to_id(d):
        d = dict(d)
        d['pair_id'] = configured_pair_id(d['graph'], d.pop('graph_path'),
                                          d.pop('annotation_path'))
        return d
    if absent == 'id':
        for r in recs:
            r['pair'] = to_id(r['pair'])
    if expected == 'id':
        exp = [to_id(d) for d in exp]
    return gs + recs, seed, exp


def r237_union(**kw):
    from metagraph.traverse import union_pairs
    so = kw.pop('seqs_order', None)
    results, seed, exp = r237_inputs(**kw)
    return union_pairs(results, seed=seed, expected_pairs=exp, seqs_order=so)


class TestUnionOfBothForms(unittest.TestCase):
    """One identity, form-independent (the rule of union.py's _KEY_FIELDS): an old document
    and its id-form twin give one part key, one tag, one union part and the same digest, in
    a union of one form or of both."""

    def same(self, u, v, what):
        self.assertEqual(u.digest, v.digest, what)
        self.assertEqual([(p.tag, p.key) for p in u.parts], [(p.tag, p.key) for p in v.parts],
                         what)
        self.assertEqual([(a.tag, a.key) for a in u.absent], [(a.tag, a.key) for a in v.absent],
                         what)
        self.assertEqual(u.summary(), v.summary(), what)
        self.assertEqual(u.to_fasta(), v.to_fasta(), what)
        self.assertEqual(u.to_gfa(), v.to_gfa(), what)
        self.assertEqual(u.provenance(), v.provenance(), what)

    def test_old_new_and_mixed_unions_are_one(self):
        old = r237_union()
        for kw in ({'parts': 'id'}, {'parts': 'id', 'absent': 'id', 'expected': 'id'},
                   {'parts': ['paths', 'id']}, {'parts': ['id', 'paths'], 'absent': 'id'},
                   {'expected': 'id'}):
            self.same(old, r237_union(**kw), kw)
        # the key holds the id, never a path
        for p in old.parts + old.absent:
            self.assertTrue(p.key[3].startswith('n1-'), p.key)
            self.assertEqual(p.key[4:6], (None, None))

    def test_one_name_with_several_pairs(self):
        # production RefSeq's shape: the ~ suffixes come from the keys, so they do not move
        # with the form either
        old = r237_union(graph='refseq')
        self.assertEqual(len(old.parts), 2)
        self.assertTrue(all(p.tag.startswith('refseq~') for p in old.parts + old.absent))
        self.assertEqual(len(set(old.pairs)), 11)
        for kw in ({'parts': 'id', 'absent': 'id', 'expected': 'id'},
                   {'parts': ['id', 'paths'], 'expected': 'id'}):
            self.same(old, r237_union(graph='refseq', **kw), kw)
        # a bare name is ambiguous there, in either form
        results, seed, exp = r237_inputs(parts='id', graph='refseq')
        from metagraph.traverse import union_pairs
        with self.assertRaises(CombineRefused) as cm:
            union_pairs(results, seed=seed, expected_pairs=exp + ['refseq'])
        self.assertEqual(cm.exception.code, 'ambiguous_pair')

    def test_an_old_document_and_its_twin_given_together_are_one_part(self):
        from metagraph.traverse import union_pairs
        text = r237_text('bact_other')
        old, new = parse(text), parse(id_form(text))
        u = union_pairs([old, new])
        self.assertEqual(len(u.parts), 1)
        self.assertEqual(u.digest, union_pairs([old]).digest)
        self.assertEqual(u.digest, union_pairs([new]).digest)
        self.assertEqual(u.digest, union_pairs([new, old]).digest)

    def test_ids_decide_without_index_fp(self):
        # the fork fixture has no index_fp: the ids (computed from paths, or stated) decide
        from metagraph.traverse import union_pairs
        a_old, a_new = fixture('fork', 'paths'), fixture('fork', 'id')
        b_new = fixture('fork', 'id', gp=GP2, ap=AP2)
        u = union_pairs([a_old, a_new])
        self.assertEqual((len(u.parts), u.digest), (1, union_pairs([a_new]).digest))
        u = union_pairs([a_old, b_new])
        self.assertEqual(len(u.parts), 2)
        self.assertEqual(sorted(p.key[3] for p in u.parts), sorted((ID, ID2)))
        self.assertEqual(u.digest, union_pairs([fixture('fork', 'paths', gp=GP2, ap=AP2),
                                                a_new]).digest)
        # an expected pair by id or by paths
        digests = {union_pairs([g], expected_pairs=exp).digest
                   for g in (a_old, a_new)
                   for exp in ([{'graph': GRAPH, 'pair_id': ID}],
                               [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP}])}
        self.assertEqual(len(digests), 1)
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([a_new], expected_pairs=[{'graph': GRAPH, 'pair_id': ID2}])
        self.assertEqual(cm.exception.code, 'not_expected')
        # a caller's partial paths (no annotation) name the pair by name and graph path
        u = union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP}, 'graphlet': a_old}])
        self.assertEqual(u.parts[0].key[3:6], (None, GP, None))

    def test_aliases_with_manifests_are_one_pair(self):
        # review m2: one pair listed under two names has two ids, and one index_fp, which
        # decides first
        from metagraph.traverse import union_pairs
        text = r237_text('bact_other')
        alias = parse(id_form(text.replace('"graph":"bact_other"', '"graph":"alias"')))
        old = parse(text)
        self.assertNotEqual(answered_pair(alias.envelope).pair_id,
                            answered_pair(old.envelope).pair_id)
        self.assertEqual(len(union_pairs([old, alias]).parts), 1)

    def test_seqs_order_by_id_or_by_paths(self):
        text = r237_text('bact_other')
        g = parse(text)
        pid = answered_pair(g.envelope).pair_id
        gp, ap = g.envelope['graph_path'], g.envelope['annotation_path']
        a = r237_union(seqs_order={('bact_other', gp, ap): 'verified'})
        for so in ({('bact_other', pid): 'verified'},
                   {Pair('bact_other', None, None, pid): 'verified'},
                   {Pair('bact_other', gp, ap, pid): 'verified'}):
            b = r237_union(parts='id', seqs_order=so)
            self.assertEqual(a.digest, b.digest, so)
            self.assertEqual([p.seqs_order for p in b.parts], ['verified', 'unknown'])
        from metagraph.traverse import union_pairs
        with self.assertRaises(CombineRefused) as cm:
            r237_union(seqs_order={('bact_other', pid): 'verified',
                                   ('bact_other', gp, ap): 'count_only'})
        self.assertEqual(cm.exception.code, 'bad_seqs_order')
        self.assertTrue(union_pairs)

    def test_p2_ids_are_refused(self):
        from metagraph.traverse import union_pairs
        g = fixture('fork', 'id')
        for call, code in (
                (lambda: union_pairs([g], expected_pairs=[{'graph': GRAPH, 'pair_id': P2}]),
                 'bad_expected'),
                (lambda: union_pairs([{'pair': {'graph': GRAPH, 'pair_id': P2}, 'graphlet': g}]),
                 'not_combinable'),
                (lambda: union_pairs([g], seqs_order={(GRAPH, P2): 'verified'}),
                 'bad_seqs_order')):
            with self.assertRaises(CombineRefused) as cm:
                call()
            self.assertEqual(cm.exception.code, code)
            self.assertIn("service's public pair id (p2)", str(cm.exception))
        # a document whose J line states a p2 id: TestMalformedPairStatements


def _doc(form='id', **tags):
    """fork's document in |form| with its J line's tags changed (None removes one)."""
    resp = tagged(T.doc_json('fork', 'graphlet'), form)
    for k, v in tags.items():
        if v is None:
            resp.pop(k, None)
        else:
            resp[k] = v
    return MT.from_response(resp['results'][0], resp)


class TestMalformedPairStatements(unittest.TestCase):
    """One rule for a pair statement the library does not take (review m2, m3), whoever
    reads it: a pair_id that is no n1 id ('', 0, False and a p2 id included), an id beside
    one path only, an id beside paths that are not its. The client raises PairIdError,
    union_pairs() and stitch() a CombineRefused with its code and words, GraphletTools the
    tool error with its code."""

    BAD = (({'pair_id': P2}, 'p2_pair_id'), ({'pair_id': ''}, 'bad_pair_id'),
           ({'pair_id': 0}, 'bad_pair_id'), ({'pair_id': False}, 'bad_pair_id'),
           ({'pair_id': ID, 'graph_path': GP}, 'bad_pair_id'),
           ({'pair_id': ID, 'annotation_path': AP}, 'bad_pair_id'),
           ({'pair_id': ID2, 'graph_path': GP, 'annotation_path': AP}, 'pair_id_mismatch'))

    def test_a_falsy_pair_id_is_refused_not_a_key_error(self):
        c = _Engine(AFTER_IDS, None)
        for pid in ('', 0, False):
            with self.assertRaises(PairIdError) as cm:
                c.wire_request({'seeds': [], 'pairs': [{'graph': GRAPH, 'pair_id': pid}]})
            self.assertEqual(cm.exception.code, 'bad_pair_id', pid)
            with self.assertRaises(PairIdError) as cm:
                pair_members(GRAPH, pair_id=pid)
            self.assertEqual(cm.exception.code, 'bad_pair_id', pid)
        with self.assertRaises(PairIdError):
            MT.TraverseClient('localhost', 1, graph=GRAPH, pair_id='')
        self.assertEqual(c.posted(), [])

    def test_a_falsy_pair_id_through_the_tools(self):
        from metagraph.traverse.mcp_tools import GraphletTools
        with tempfile.TemporaryDirectory() as d:
            client = _Engine(AFTER_IDS, tagged(T.doc_json('fork', 'graphlet'), 'id'))
            tools = GraphletTools(GraphletStore(d), {'mini': client})
            h = tools.traverse_fetch(seed={'sequence': 'ACGT'}, index='mini')['handle']
            for pid in ('', 0, False):
                out = tools.traverse_continue(h, 'right', 0, overrides={
                    'pairs': [{'graph': GRAPH, 'pair_id': pid}]})
                self.assertEqual(out['error'], 'bad_pair_id', (pid, out))

    def test_every_reader_refuses_with_one_code(self):
        from metagraph.traverse import union_pairs
        from metagraph.traverse.mcp_tools import GraphletTools
        good = fixture('fork', 'paths')
        for tags, code in self.BAD:
            doc = _doc('id', **tags)
            with self.assertRaises(PairIdError) as cm:
                answered_pair(doc.envelope)
            self.assertEqual(cm.exception.code, code, tags)
            words = str(cm.exception)
            with self.assertRaises(PairIdError) as cm:
                GraphletStore.request_of(doc)
            self.assertEqual(cm.exception.code, code, tags)
            for call in (lambda: union_pairs([doc]), lambda: union_pairs([good, doc]),
                         lambda: union_pairs([{'pair': GRAPH, 'graphlet': doc}])):
                with self.assertRaises(CombineRefused) as cm:
                    call()
                self.assertEqual(cm.exception.code, code, tags)
                self.assertIn(words, str(cm.exception))
                self.assertRegex(str(cm.exception), r'^result \d+: ')
            for call, who in ((lambda: stitch(doc), 'the parent'),
                              (lambda: stitch(good, doc), 'the child')):
                with self.assertRaises(CombineRefused) as cm:
                    call()
                self.assertEqual(cm.exception.code, code, (tags, who))
                self.assertEqual(str(cm.exception), '%s: %s' % (who, words))
            with tempfile.TemporaryDirectory() as d:
                path = os.path.join(d, 'bad.mgt')
                doc.save(path)
                out = GraphletTools(GraphletStore(d), {}).graphlet_load(path)
                self.assertEqual(out.get('error'), code, (tags, out))
            for path in (GP, AP):
                self.assertNotIn(path, words)

    def test_an_empty_path_is_a_path(self):
        # master's vector 3 hashes an empty field as "0:"; a document stating an empty
        # annotation path names the pair the engine names by that id
        want = configured_pair_id(GRAPH, GP, '')
        self.assertNotEqual(want, configured_pair_id(GRAPH, GP, None))
        old = _doc('paths', annotation_path='')
        self.assertEqual(answered_pair(old.envelope), Pair(GRAPH, GP, '', want))
        twin = _doc('id', pair_id=want)
        self.assertEqual(answered_pair(twin.envelope).identity, (GRAPH, want))
        from metagraph.traverse import union_pairs
        u1, u2 = union_pairs([old]), union_pairs([twin])
        self.assertEqual((u1.parts[0].key, u1.digest), (u2.parts[0].key, u2.digest))
        self.assertEqual(u1.parts[0].key[3], want)
        # the engine's listing with the empty path resolves the id (n1 recomputed)
        caps = {'graph_summary': {'pairs': [{'graph': GRAPH, 'graph_path': GP,
                                             'annotation_path': ''}]}}
        self.assertEqual(MT.client.capabilities_pair_resolver(caps)(want), (GRAPH, GP, ''))

    def test_an_empty_path_is_written_as_it_is_read(self):
        # fix review N1: one rule for '' -- a value, as configured_pair_id() hashes it -- in
        # the writers too. An old document whose annotation path is "" loads, is re-asked,
        # continued and deepened by its paths, and unions with provenance by its id
        from metagraph.traverse import union_pairs
        from metagraph.traverse.client import request_without_paths
        from metagraph.traverse.mcp_tools import GraphletTools
        want = configured_pair_id(GRAPH, GP, '')
        old = _doc('paths', annotation_path='')
        sel = [{'graph': GRAPH, 'graph_path': GP, 'annotation_path': ''}]
        req = GraphletStore.request_of(old)
        self.assertEqual(req['pairs'], sel)
        self.assertEqual(request_without_paths(req)['pairs'],
                         [{'graph': GRAPH, 'pair_id': want}])
        leaf = [old.walks('right')[0].path_id]
        self.assertEqual(old.next_request('right', leaf)['pairs'], sel)
        for caps in (BEFORE_IDS, AFTER_IDS):
            c = _Engine(caps, tagged(T.doc_json('fork', 'graphlet'), 'paths', ap=''))
            c.deepen(old, 'right', leaf)
            self.assertEqual(c.posted()[0]['pairs'], sel)
        (part,) = [e for e in union_pairs([old]).provenance() if e['kind'] == 'part']
        self.assertEqual(part['pair']['pair_id'], want)
        self.assertEqual(part['request']['pairs'], [{'graph': GRAPH, 'pair_id': want}])
        self.assertNotIn(GP, json.dumps(part))
        with tempfile.TemporaryDirectory() as d:
            tools = GraphletTools(GraphletStore(os.path.join(d, 's')), {}, export_dir=d)
            path = os.path.join(d, 'old.mgt')
            old.save(path)
            out = tools.graphlet_load(path)
            self.assertNotIn('error', out, out)
            self.assertEqual(tools.store.get(out['handle']).request['pairs'], sel)
        # a probe of that pair names its empty path
        c = _Engine(BEFORE_IDS, None)
        c.capabilities(GRAPH, GP, '')
        self.assertEqual(c.probes(), ['/traverse/capabilities?graph=refseq&graph_path='
                                      '%2Fdata%2Frefseq%2Fbact_other.dbg&annotation_path='])

    def test_a_fan_out_sends_and_expects_one_pair(self):
        # fix review N1: a pair object of empty paths is sent as given and expected by its
        # id, so the fan-out's union reads its own results, before and after ids
        from metagraph.traverse import union_pairs
        want = configured_pair_id(GRAPH, '', '')
        seed = 'ACGTACGTACGTACGTACGTACGTACGTACGTAC'
        for caps, form in ((BEFORE_IDS, 'paths'), (AFTER_IDS, 'id')):
            c = _Engine(caps, tagged(T.doc_json('fork', 'graphlet'), form, gp='', ap=''))
            r = c.traverse_pairs(seed, pairs=[{'graph': GRAPH, 'graph_path': '',
                                               'annotation_path': ''}], coordinates=None)
            self.assertEqual(c.posted()[0]['pairs'],
                             [{'graph': GRAPH, 'graph_path': '', 'annotation_path': ''}])
            self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': want},))
            u = union_pairs(r.results, expected_pairs=r.expected_pairs)
            self.assertEqual((u.expected, u.parts[0].key[3]), ((GRAPH,), want))

    def test_a_caller_empty_path_is_a_value(self):
        # fix review N1: a caller's '' in a union pair object is a path, not refused
        from metagraph.traverse import union_pairs
        old = _doc('paths', annotation_path='')
        full = {'graph': GRAPH, 'graph_path': GP, 'annotation_path': ''}
        u = union_pairs([{'pair': full, 'graphlet': old}], expected_pairs=[dict(full)])
        self.assertEqual((u.parts[0].pair_stated_by, u.expected), ('caller', (GRAPH,)))
        self.assertEqual(u.parts[0].key[3], configured_pair_id(GRAPH, GP, ''))
        # and '' is not the old document's annotation path
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': full, 'graphlet': fixture('fork', 'paths')}])
        self.assertEqual(cm.exception.code, 'pair_mismatch')

    def test_an_empty_pair_id_is_refused_in_the_union(self):
        # fix review N2: pair_id '' is an id to check, never read as a bare name
        from metagraph.traverse import union_pairs
        a = fixture('fork', 'paths')
        for call, code in (
                (lambda: union_pairs([{'pair': {'graph': GRAPH, 'pair_id': ''},
                                       'graphlet': a}]), 'not_combinable'),
                (lambda: union_pairs([a], expected_pairs=[{'graph': GRAPH, 'pair_id': ''}]),
                 'bad_expected'),
                (lambda: union_pairs([a], seqs_order={Pair(GRAPH, None, None, ''): 'verified'}),
                 'bad_seqs_order')):
            with self.assertRaises(CombineRefused) as cm:
                call()
            self.assertEqual(cm.exception.code, code)
            self.assertIn('is not an n1 pair id', str(cm.exception))
        # (graph, '') names the empty graph path: no pair of these has it
        b = fixture('fork', 'paths', gp=GP2, ap=AP2)
        u = union_pairs([a, b], seqs_order={(GRAPH, ''): 'verified'})
        self.assertEqual({p.seqs_order for p in u.parts}, {'unknown'})


class TestPartialDescriptions(unittest.TestCase):
    """A description that names its pair by a path and no id (a caller's {graph,
    graph_path}, a (graph, graph_path[, None]) seqs_order key) is compared with a pair named
    by its id by the paths the call knows for that id; with none it does not match (review
    M2). One set of inputs gives one union in every order."""

    def setUp(self):
        self.old, self.new = fixture('fork', 'paths'), fixture('fork', 'id')

    def test_a_caller_pair_by_path_is_compared_or_refused(self):
        from metagraph.traverse import union_pairs
        # another pair's path: refused in both forms
        for g in (self.old, self.new):
            with self.assertRaises(CombineRefused) as cm:
                union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP2}, 'graphlet': g}])
            self.assertEqual(cm.exception.code, 'pair_mismatch')
        # its own path: compared with the old document's paths; the id-form document alone
        # cannot be compared with it, and is refused, never accepted unchecked
        u = union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP}, 'graphlet': self.old}])
        self.assertEqual(u.parts[0].pair_stated_by, 'caller')
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP}, 'graphlet': self.new}])
        self.assertEqual(cm.exception.code, 'pair_mismatch')
        self.assertIn('the two cannot be compared', str(cm.exception))
        # an expected pair by path: matched against the id-form answer by the paths its old
        # twin states in the same call, else not matched (not_expected)
        exp = [{'graph': GRAPH, 'graph_path': GP}]
        self.assertEqual(union_pairs([self.old], expected_pairs=exp).expected, (GRAPH,))
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([self.new], expected_pairs=exp)
        self.assertEqual(cm.exception.code, 'not_expected')

    def test_one_union_in_every_order(self):
        import itertools
        from metagraph.traverse import union_pairs
        exp = [{'graph': GRAPH, 'graph_path': GP}, {'graph': GRAPH, 'pair_id': ID}]
        for so in ({(GRAPH, GP): 'verified'}, {(GRAPH, GP2): 'verified'},
                   {(GRAPH, GP, None): 'count_only'}, {(GRAPH, ID): 'verified'}):
            got = set()
            for order in itertools.permutations([self.old, self.new]):
                u = union_pairs(list(order), expected_pairs=exp, seqs_order=so)
                got.add((u.digest, u.parts[0].seqs_order, u.parts[0].key))
            self.assertEqual(len(got), 1, so)
            # GP is this pair's path, GP2 another's
            self.assertEqual(got.pop()[1], 'unknown' if GP2 in next(iter(so)) else so[
                next(iter(so))], so)

    def test_a_certificate_attaches_to_its_own_pair_only(self):
        from metagraph.traverse import union_pairs
        so = {(GRAPH, GP2, None): 'verified'}
        old = [fixture('fork', 'paths'), fixture('fork', 'paths', gp=GP2, ap=AP2)]
        new = [fixture('fork', 'id'), fixture('fork', 'id', gp=GP2, ap=AP2)]
        u = union_pairs(old, seqs_order=so)
        self.assertEqual({p.key[3]: p.seqs_order for p in u.parts},
                         {ID: 'unknown', ID2: 'verified'})
        # by id only, the path key cannot be compared: no pair gets it, never the other one
        u = union_pairs(new, seqs_order=so)
        self.assertEqual({p.key[3]: p.seqs_order for p in u.parts},
                         {ID: 'unknown', ID2: 'unknown'})
        # with the old documents given too, the paths are known: as for the old form
        u = union_pairs(new + old, seqs_order=so)
        self.assertEqual({p.key[3]: p.seqs_order for p in u.parts},
                         {ID: 'unknown', ID2: 'verified'})

    def test_a_malformed_id_in_a_seqs_order_key_is_refused(self):
        from metagraph.traverse import union_pairs
        for bad in (ID.replace('n1-', 'N1-'), 'n2-' + ID[3:], ID.upper(),
                    'P2-0123456789ABCDEF0123456789ABCDEF', P2):
            for key in ((GRAPH, bad), (GRAPH, bad, None)):
                for g in (self.old, self.new):
                    with self.assertRaises(CombineRefused) as cm:
                        union_pairs([g], seqs_order={key: 'verified'})
                    self.assertEqual(cm.exception.code, 'bad_seqs_order', key)
                    self.assertIn('pair id', str(cm.exception))
        # a relative graph path with a dash is a path
        u = union_pairs([self.old], seqs_order={(GRAPH, 'graph-1.dbg'): 'verified'})
        self.assertEqual(u.parts[0].seqs_order, 'unknown')

    def test_a_caller_path_is_echoed_where_there_is_no_id(self):
        # review nits n4, n5: a pair the caller stated by a path without its id is shown in
        # its own words -- in refusal texts (two such pairs read apart) and in provenance,
        # pair and replay request alike; a pair with an id never shows a path
        from metagraph.traverse import union_pairs
        with self.assertRaises(CombineRefused) as cm:
            union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP}, 'error': 'x'},
                         {'pair': {'graph': GRAPH, 'graph_path': GP2}, 'error': 'y'}],
                        expected_pairs=[GRAPH], seed='ACGT')
        self.assertEqual(cm.exception.code, 'ambiguous_pair')
        self.assertIn(GP, str(cm.exception))
        self.assertIn(GP2, str(cm.exception))
        u = union_pairs([{'pair': {'graph': GRAPH, 'graph_path': GP}, 'graphlet': self.old}])
        (part,) = [e for e in u.provenance() if e['kind'] == 'part']
        self.assertEqual((part['pair']['graph_path'], part['request']['graph_path']), (GP, GP))
        self.assertIsNone(part['pair']['annotation_path'])
        u = union_pairs([{'pair': {'graph': GRAPH}, 'graphlet': self.old}])
        self.assertNotIn(GP, json.dumps(u.provenance()))

    def test_a_fan_out_pair_by_path_takes_the_answering_pair_id(self):
        from metagraph.traverse import union_pairs
        seed = 'ACGTACGTACGTACGTACGTACGTACGTACGTAC'
        for caps, form in ((AFTER_IDS, 'id'), (BEFORE_IDS, 'paths')):
            c = _Engine(caps, tagged(T.doc_json('fork', 'graphlet'), form))
            r = c.traverse_pairs(seed, pairs=[{'graph': GRAPH, 'graph_path': GP}],
                                 coordinates=None)
            self.assertEqual(c.posted()[0]['graph_path'], GP)
            self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': ID},), form)
            self.assertEqual(r.seqs_order, {(GRAPH, ID): 'verified'}, form)
            u = union_pairs(r.results, expected_pairs=r.expected_pairs,
                            seqs_order=r.seqs_order)
            self.assertEqual((u.expected, u.parts[0].seqs_order), ((GRAPH,), 'verified'))
        # an answer on another pair than the path given: kept as given
        c = _Engine(BEFORE_IDS, tagged(T.doc_json('fork', 'graphlet'), 'paths', gp=GP2,
                                       ap=AP2))
        r = c.traverse_pairs(seed, pairs=[{'graph': GRAPH, 'graph_path': GP}],
                             coordinates=None)
        self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'graph_path': GP},))


def configured_paths(*docs):
    """Every configured path the documents |docs| (texts) state."""
    out = set()
    for text in docs:
        j = json.loads(next(line for line in text.split('\n') if line.startswith('J '))[2:])
        out.update(v for v in (j.get('graph_path'), j.get('annotation_path')) if v)
    with open(os.path.join(R237, 'absent.json'), encoding='utf-8') as f:
        for a in json.load(f):
            out.update((a['pair']['graph_path'], a['pair']['annotation_path']))
    return out


class TestOutputsWithoutPaths(unittest.TestCase):
    """The owner's decision (docs/DESIGN-pair-identity.md section 2): no path in any
    response. What the library answers names a pair by graph and pair_id -- computed for an
    old document -- in both forms; inputs keep taking paths."""

    def test_the_union_of_old_documents_states_no_path(self):
        u = r237_union()
        paths = configured_paths(*(r237_text(n) for n in R237_PARTS))
        self.assertEqual(len(paths), 22)
        outs = [u.summary(), u.provenance(), u.evidence(), u.to_fasta(), u.to_gfa(),
                str(u), [(p.tag, p.key) for p in u.parts + u.absent]]
        text = json.dumps(outs, default=str)
        self.assertEqual([p for p in paths if p in text], [])
        # the replay request names the pair by id
        (prov,) = [e for e in u.provenance() if e.get('tag') == 'bact_other']
        self.assertEqual(prov['request']['pairs'], [{'graph': 'bact_other',
                                                     'pair_id': prov['pair']['pair_id']}])
        self.assertEqual((prov['pair']['graph_path'], prov['pair']['annotation_path']),
                         (None, None))

    def test_refusal_texts_state_no_path(self):
        from metagraph.traverse import union_pairs
        old = parse(r237_text('bact_other'))
        paths = configured_paths(r237_text('bact_other'))
        texts = []
        for call in (lambda: union_pairs([{'pair': 'other', 'graphlet': old}]),
                     lambda: union_pairs([old], seqs_order={
                         ('bact_other', old.envelope['graph_path'],
                          old.envelope['annotation_path']): 'verified',
                         ('bact_other', answered_pair(old.envelope).pair_id): 'count_only'})):
            with self.assertRaises(CombineRefused) as cm:
                call()
            texts.append(str(cm.exception))
        self.assertIn(answered_pair(old.envelope).pair_id, texts[0])
        self.assertEqual([p for p in paths if p in ' '.join(texts)], [], texts)

    def test_stitch_provenance_and_refusal_state_no_path(self):
        def load(role, form, **kw):
            with open(os.path.join(T.STITCH, 'keep_forbid.%s.json' % role),
                      encoding='utf-8') as f:
                resp = tagged(json.load(f), form, **kw)
            return MT.from_response(resp['results'][0], resp)
        sg = stitch(load('parent', 'paths'), load('child', 'id'))
        prov = sg.provenance()
        self.assertEqual([e['pair'] for e in prov], [{'graph': GRAPH, 'pair_id': ID}] * 2)
        self.assertEqual([e['request']['pairs'] for e in prov], [BY_ID] * 2)
        text = json.dumps(prov, default=str)
        self.assertNotIn(GP, text)
        with self.assertRaises(CombineRefused) as cm:
            stitch(load('parent', 'paths'), load('child', 'paths', gp=GP2, ap=AP2))
        self.assertEqual(str(cm.exception),
                         'the child was answered on pair {"graph":"refseq","pair_id":"%s"}, its '
                         'parent on pair {"graph":"refseq","pair_id":"%s"}: a continuation runs '
                         'on its parent\'s pair' % (ID2, ID))

    def test_discovery_states_pairs_by_id(self):
        c = _FanOut(TestDiscoveryInBothForms.BY_PATHS)
        r = c.traverse_pairs('ACGT', coordinates=None)
        # sent by the listed paths, stated by their ids
        self.assertEqual(c.sent, BY_PATHS + [{'graph': GRAPH, 'graph_path': GP2,
                                              'annotation_path': AP2}])
        self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': ID},
                                            {'graph': GRAPH, 'pair_id': ID2}))
        self.assertEqual(r.seqs_order, {(GRAPH, ID): 'verified'})
        # a caller's pair by its paths too; one without the annotation path as given
        c = _FanOut(TestDiscoveryInBothForms.BY_PATHS)
        r = c.traverse_pairs('ACGT', pairs=[dict(BY_PATHS[0]), {'graph': GRAPH,
                                                                'graph_path': GP2}],
                             coordinates=None)
        self.assertEqual(r.expected_pairs, ({'graph': GRAPH, 'pair_id': ID},
                                            {'graph': GRAPH, 'graph_path': GP2}))
        self.assertEqual(r.seqs_order, {(GRAPH, ID): 'verified'})

    def test_request_and_entry_without_paths(self):
        from metagraph.traverse.client import entry_without_paths, request_without_paths
        req = {'seeds': [], 'pairs': [dict(BY_PATHS[0])]}
        self.assertEqual(request_without_paths(req), {'seeds': [], 'pairs': BY_ID})
        self.assertEqual(req['pairs'], BY_PATHS)                   # not changed in place
        for same in ({'seeds': [], 'graph': GRAPH, 'graph_path': GP}, {'seeds': []},
                     {'seeds': [], 'pairs': BY_ID}, None):
            self.assertEqual(request_without_paths(same), same)
        entry = {'graph': GRAPH, 'graph_path': GP, 'annotation_path': AP, 'index_fp': 'f',
                 'k': 31}
        self.assertEqual(entry_without_paths(entry), {'graph': GRAPH, 'pair_id': ID,
                                                      'index_fp': 'f', 'k': 31})
        for same in ({'k': 31}, {'graph': GRAPH, 'pair_id': ID},
                     {'graph': None, 'graph_path': None, 'annotation_path': None}):
            self.assertEqual(entry_without_paths(same), same)


class TestDocumentsHandedOut(unittest.TestCase):
    """A document the library hands out names its pair by {graph, pair_id} in its J line,
    an old one too (the owner's rule, review m5); the body and its digest do not change, and
    the store's own files keep what they were stored with."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__('shutil').rmtree(self.d, ignore_errors=True))
        self.store = GraphletStore(os.path.join(self.d, 'store'))
        path = os.path.join(self.d, 'old.mgt')
        fixture('fork', 'paths').save(path)
        self.h = self.store.load(path)

    def j(self, text):
        return json.loads(next(line for line in text.split('\n') if line.startswith('J '))[2:])

    def body(self, text):
        return [line for line in text.split('\n') if not line.startswith('J ')]

    def test_texts_and_files(self):
        from metagraph.traverse.client import document_without_paths
        kept = self.store.standalone_text(self.h)
        self.assertEqual(self.j(kept)['graph_path'], GP)          # as stored
        out = self.store.standalone_text(self.h, paths=False)
        self.assertEqual(out, document_without_paths(kept))
        j = self.j(out)
        self.assertEqual((j['graph'], j['pair_id']), (GRAPH, ID))
        self.assertNotIn('graph_path', j)
        self.assertNotIn('annotation_path', j)
        self.assertEqual(self.body(out), self.body(kept))
        self.assertTrue(MT.is_canonical(out))
        for write in (lambda p: self.store.save(self.h, p, paths=False),
                      lambda p: self.store.save_body(self.h, p, paths=False),
                      lambda p: self.store.save_body(self.h, p, paths=False,
                                                     budget=MT.LocalBudget(MT.LocalLimits()))):
            p = os.path.join(self.d, 'out.mgt')
            write(p)
            with open(p, encoding='utf-8', newline='') as f:
                text = f.read()
            for path in (GP, AP):
                self.assertNotIn(path, text)
            g = parse(text)
            self.assertEqual(answered_pair(g.envelope), Pair(GRAPH, None, None, ID))
            self.assertEqual(g.body_digest, parse(kept).body_digest)
        # an id-form text, or one without a J line, is returned as it is
        self.assertIs(document_without_paths(out), out)
        bare = '\n'.join(self.body(kept))
        self.assertIs(document_without_paths(bare), bare)


class TestToolOutputsWithoutPaths(unittest.TestCase):
    """GraphletTools' answers name a pair by graph and pair_id, no path."""

    def test_capabilities_and_a_continuation_request(self):
        from metagraph.traverse.mcp_tools import GraphletTools

        class Before(_Engine):
            # a probe of one pair on an engine before ids: its entry with the paths
            def capabilities(self, *a, **kw):
                return dict(super().capabilities(*a, **kw), graph=GRAPH, graph_path=GP,
                            annotation_path=AP)
        old = fixture('fork', 'paths')
        with tempfile.TemporaryDirectory() as d:
            client = Before(BEFORE_IDS, tagged(T.doc_json('fork', 'graphlet'), 'paths'))
            tools = GraphletTools(GraphletStore(d), {'mini': client})
            caps = tools.traverse_capabilities('mini')['capabilities']
            self.assertEqual((caps['graph'], caps['pair_id']), (GRAPH, ID))
            self.assertNotIn('graph_path', caps)
            path = os.path.join(d, 'old.mgt')
            old.save(path)
            h = tools.graphlet_load(path)['handle']
            out = tools.traverse_continue(h, 'right', old.walks('right')[0].path_id,
                                          execute=False)
            self.assertEqual(out['request']['pairs'], BY_ID)
            text = json.dumps(out)
            self.assertNotIn(GP, text)
            # executed, it is sent by the old document's paths
            done = tools.traverse_continue(h, 'right', old.walks('right')[0].path_id)
            self.assertNotIn('error', done, done)
            self.assertEqual(client.posted()[-1]['pairs'], BY_PATHS)


    def test_exports_and_saves_state_no_path(self):
        # graphlet_export(format=mgt) and graphlet_save write the J line by {graph, pair_id},
        # an old document's too, with and without local limits; the body is the stored one
        from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits
        old = fixture('fork', 'paths')
        for limits in (None, ToolLimits()):
            with tempfile.TemporaryDirectory() as d:
                ex = os.path.join(d, 'exports')
                tools = GraphletTools(GraphletStore(os.path.join(d, 's')), {},
                                      export_dir=ex, local_limits=limits)
                src = os.path.join(ex, 'old.mgt')
                os.makedirs(ex)
                old.save(src)
                h = tools.graphlet_load(src)['handle']
                outs = [tools.graphlet_export(h, format='mgt', path='e.mgt'),
                        tools.graphlet_save(h, 's.mgt')]
                for out in outs:
                    self.assertNotIn('error', out, out)
                    with open(out['path'], encoding='utf-8', newline='') as f:
                        text = f.read()
                    for path in (GP, AP):
                        self.assertNotIn(path, text, limits)
                    g = parse(text)
                    self.assertEqual(answered_pair(g.envelope), Pair(GRAPH, None, None, ID))
                    self.assertEqual(g.body_digest, old.body_digest)
                    # loaded again, it is the same pair, re-asked by its id
                    h2 = tools.graphlet_load(out['path'])['handle']
                    self.assertEqual(GraphletStore.request_of(tools.store.graphlet(h2))
                                     ['pairs'], BY_ID)

if __name__ == '__main__':
    unittest.main()
