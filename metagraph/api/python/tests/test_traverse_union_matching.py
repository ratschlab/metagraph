"""The union's pair matching by candidates (union.py's _Candidates, _check_keys(),
_pair_classes(), _seqs_of()) against the all-pairs code it replaced (3e9a2c131, kept here
verbatim as the reference): the same matches, classes and class descriptions in the same
order, the same seqs_order values, and the same refusal -- code, the descriptions it names,
its text -- for every input. Randomized with fixed seeds over descriptions that mix names,
several pairs of one name, ids, both paths, one path, '' paths, known paths, index_fp,
index_meta_fp, index_ns, release, single-graph pairs, and seqs_order keys of every form
(every value a string or None: a description holds no other, union.py's _facets()); then
whole union_pairs() calls (failed and absent records, caller statements with a value that is
not a string, parts in both document forms, expected pairs, keys) against union_pairs() run
with the reference code; then the charges: each pass is charged before it runs, for what it
runs."""

import copy
import random
import unittest
from unittest import mock

from tests import traverse_testlib as T  # noqa: F401  (puts the library on sys.path)

from metagraph.traverse import LocalBudget, LocalBudgetExceeded, union_pairs
from metagraph.traverse import union as U
from metagraph.traverse import from_response
from metagraph.traverse.client import Pair, PairIdError, configured_pair_id


# ---------------------------------------------------------------- the all-pairs reference
# union.py at 3e9a2c131: _check_args()' loop over the keys, _build()'s ambiguous_pair loop,
# _classes() and _seqs_of(), unchanged but for the U. prefixes (and no budget).

def old_check_keys(so, known, b=None):
    for i in range(len(so)):
        for j in range(i + 1, len(so)):
            if so[i][1] != so[j][1] and U._matches(so[i][0], so[j][0], known):
                raise U._refuse('bad_seqs_order', None,
                                'seqs_order[%s] and seqs_order[%s] (two values for one pair)'
                                % (U._key_text(so[i][2], so[i][0]),
                                   U._key_text(so[j][2], so[j][0])))


def old_classes(descs, known=None):
    classes = []
    for d in descs:
        hit = None
        for c in classes:
            if U._matches(d, c['members'][0], known):
                hit = c
                break
        if hit is None:
            classes.append({'members': [d]})
        else:
            hit['members'].append(d)
    for c in classes:
        c['desc'] = min(c['members'], key=lambda d: (-U._fullness(d), U._sort_key(U._key(d)),
                                                     not U._has_path(d),
                                                     U.C.canonical_json(d)))
    return classes


def old_pair_classes(descs, known=None, b=None, check=True):
    if check:
        for d in descs:
            hits = [x for x in descs if x is not d and U._matches(d, x, known)]
            for i in range(len(hits)):
                for j in range(i + 1, len(hits)):
                    if not U._matches(hits[i], hits[j], known):
                        raise U._refuse('ambiguous_pair', None, U._show(d), U._show(hits[i]),
                                        U._show(hits[j]))
    return old_classes(descs, known)


def old_seqs_of(classes, so, tags, known=None, b=None):
    seqs = {}
    for c in classes:
        vals = {v for d, v, _ in so if U._matches(d, c['desc'], known)}
        if len(vals) > 1:
            d, k = next((d, k) for d, v, k in so if U._matches(d, c['desc'], known))
            raise U._refuse('bad_seqs_order', None, 'seqs_order[%s] (two values for pair %s)'
                            % (U._key_text(k, d), tags[id(c)]))
        seqs[id(c)] = vals.pop() if vals else 'unknown'
    return seqs


def with_reference():
    """union_pairs() as it was: the reference passes patched in."""
    return mock.patch.multiple(U, _check_keys=old_check_keys, _pair_classes=old_pair_classes,
                               _seqs_of=old_seqs_of)


def outcome(fn, *args, **kw):
    """('ok', value) or ('raised', type, code, text); a stop and a failed assertion (a
    test's own, inside a mock) are raised."""
    try:
        return ('ok', fn(*args, **kw))
    except (LocalBudgetExceeded, AssertionError):
        raise
    except Exception as e:  # noqa: BLE001 -- every refusal and error is compared
        return ('raised', type(e).__name__, getattr(e, 'code', None), str(e))


def class_shape(descs, classes):
    """The classes by positions: members' positions (by identity) and the description's."""
    where = {}
    for i, d in enumerate(descs):
        where.setdefault(id(d), i)
    return [([where[id(m)] for m in c['members']], where[id(c['desc'])]) for c in classes]


# ---------------------------------------------------------------- random descriptions

GRAPHS = ('g1', 'g2', 'g3')
NS = ('fixtures', 'ns2', 'g1')        # 'g1': a bare name that is also an index_ns
GPS = ('/a.dbg', '/b.dbg', '')
APS = ('/a.anno', '/b.anno', '')
RELEASES = ('r1', 'r2')
FPS = ('fpA', 'fpB', 'fpC')
MFPS = ('m1', 'm2')
LONE_ID = configured_pair_id('g1', '/elsewhere.dbg', '/elsewhere.anno')   # paths never stated


def universe(rng):
    """A few pairs, each with its tags (graph or none, paths, release, index_ns, index_fp,
    index_meta_fp): what the descriptions of a case state, wholly or in part."""
    out = []
    for _ in range(rng.randint(1, 5)):
        g = rng.choice(GRAPHS) if rng.random() < 0.85 else None
        p = {'graph': g, 'graph_path': rng.choice(GPS), 'annotation_path': rng.choice(APS),
             'release': rng.choice(RELEASES), 'index_ns': rng.choice(NS),
             'index_fp': rng.choice(FPS), 'index_meta_fp': rng.choice(MFPS)}
        # a single-graph pair's id too: a caller may state one beside its index_ns
        p['pair_id'] = configured_pair_id(g or p['index_ns'], p['graph_path'],
                                          p['annotation_path'])
        out.append(p)
    return out


def view(rng, p, wild=False):
    """A description of pair |p| in one of the forms a union reads, or of a nearby pair (a
    field changed), through the library's own constructors (_pair_desc(), _with_id());
    with |wild|, now and then one with a value that is not a string: a caller's pair
    statement, which union_pairs() refuses before any pass reads it (no description holds
    such a value: _pair_desc(), and an envelope's bad_pair_tag)."""
    p = dict(p)
    if rng.random() < 0.25:     # a nearby pair: one field differs
        k = rng.choice(('graph', 'graph_path', 'annotation_path', 'release', 'index_ns',
                        'index_fp', 'index_meta_fp'))
        pool = {'graph': GRAPHS, 'graph_path': GPS, 'annotation_path': APS,
                'release': RELEASES, 'index_ns': NS, 'index_fp': FPS,
                'index_meta_fp': MFPS}[k]
        p[k] = rng.choice(pool)
        p['pair_id'] = configured_pair_id(p.get('graph') or p['index_ns'], p['graph_path'],
                                          p['annotation_path'])
    form = rng.random()
    d = {}
    if p.get('graph'):
        if form < 0.15:
            return {'graph': p['graph']}                       # a bare name
        d['graph'] = p['graph']
        r = rng.random()
        if r < 0.3:
            d['pair_id'] = p['pair_id']                        # by its id (a new document)
        elif r < 0.55:
            d['graph_path'], d['annotation_path'] = p['graph_path'], p['annotation_path']
        elif r < 0.65:
            d[rng.choice(U._PATHS)] = p[rng.choice(U._PATHS)]  # one path: no id
        elif r < 0.7:
            d['pair_id'] = LONE_ID                             # an id no path is known for
    elif rng.random() < 0.2 and p.get('index_ns') in GRAPHS:
        return {'graph': p['index_ns']}                        # a bare name for its index_ns
    else:                                                      # no graph: by index_ns, and
        r = rng.random()
        if r < 0.25:
            d['pair_id'] = p['pair_id']                        # an id
        elif r < 0.4:
            d[rng.choice(U._PATHS)] = p[rng.choice(U._PATHS)]  # a path
        elif r < 0.5:
            d['pair_id'], d['graph_path'] = p['pair_id'], p['graph_path']  # an id and a path
    for k, q in (('release', 0.35), ('index_ns', 0.35), ('index_fp', 0.3),
                 ('index_meta_fp', 0.3)):
        if rng.random() < q:
            d[k] = p[k]
    if not d.get('graph') and not d.get('index_ns'):
        d['index_ns'] = p['index_ns']
    try:
        out = U._pair_desc(d)
    except PairIdError:
        return None
    if wild and out is not None and rng.random() < 0.04:   # not a string: refused
        out = dict(out, **{rng.choice(('release', 'index_fp', 'index_meta_fp',
                                       'index_ns')): rng.choice((5, 0.5, ['x'], True,
                                                                 float('nan')))})
    return out


def key_form(rng, p):
    """A seqs_order key of pair |p| (or of a nearby one) in one of the forms
    union_pairs() takes."""
    g = p.get('graph') or rng.choice(GRAPHS)
    r = rng.random()
    gp, ap = p['graph_path'], p['annotation_path']
    if rng.random() < 0.2:
        gp, ap = rng.choice(GPS), rng.choice(APS)
    if r < 0.25:
        return g
    if r < 0.3:
        return (g,)
    if r < 0.5:
        return (g, configured_pair_id(g, gp, ap) if rng.random() < 0.9 else LONE_ID)
    if r < 0.6:
        return (g, gp)
    if r < 0.8:
        return (g, gp, ap)
    if r < 0.85:
        return (g, None, ap)
    if r < 0.95:
        return Pair(g, gp, ap, None)
    return Pair(g, None, None, configured_pair_id(g, gp, ap))


def key_list(rng, pairs):
    """seqs_order= of up to six keys, mostly one value (a conflict now and then)."""
    so, seen = [], set()
    own = rng.choice(('verified', 'count_only', None))
    for _ in range(rng.randint(0, 6)):
        k = key_form(rng, rng.choice(pairs))
        if k in seen:
            continue
        seen.add(k)
        so.append((k, own if rng.random() < 0.75 else rng.choice(('verified', 'count_only',
                                                                     None))))
    return dict(so)


def parsed_keys(seqs_order):
    """The (description, value, key) list _check_args() builds, unchecked among
    themselves; None where a key is refused on its own."""
    with mock.patch.object(U, '_check_keys', lambda *a: None):
        try:
            return U._check_args([], None, None, seqs_order)[1]
        except U.C.CombineRefused:
            return None


def case(rng):
    pairs = universe(rng)
    descs = []
    for _ in range(rng.randint(0, 14)):
        r = rng.random()
        if r < 0.12 and descs:
            descs.append(dict(rng.choice(descs)))       # the same content again
        elif r < 0.15 and descs:
            descs.append(rng.choice(descs))             # the same object again
        else:
            d = view(rng, rng.choice(pairs))
            if d is not None:
                descs.append(d)
    so = parsed_keys(key_list(rng, pairs)) or []
    return descs, so


class TestPasses(unittest.TestCase):
    """Each pass against its reference, on 6,000 random cases (seeds 0..5999)."""

    def _compare(self, seed):
        rng = random.Random(seed)
        descs, so = case(rng)
        # as _build() reads them: every description of the call
        known = U._known_paths(descs + [d for d, _, _ in so])
        # the keys among themselves
        self.assertEqual(outcome(U._check_keys, so, known),
                         outcome(old_check_keys, so, known), seed)
        # the classes alone (greedy, whatever ambiguous_pair says)
        self.assertEqual(class_shape(descs, U._classes(descs, known)),
                         class_shape(descs, old_classes(descs, known)), seed)
        # ambiguous_pair, then the classes
        new, old = (outcome(f, descs, known) for f in (U._pair_classes, old_pair_classes))
        if old[0] == 'ok':
            self.assertEqual(new[0], 'ok', (seed, new))
            self.assertEqual(class_shape(descs, new[1]), class_shape(descs, old[1]), seed)
            classes = new[1]
            tags = {id(c): 't%d' % i for i, c in enumerate(classes)}
            a, b = outcome(U._seqs_of, classes, so, tags, known), \
                outcome(old_seqs_of, classes, so, tags, known)
            if b[0] == 'ok':
                self.assertEqual(a[0], 'ok', (seed, a))
                self.assertEqual([a[1][id(c)] for c in classes],
                                 [b[1][id(c)] for c in classes], seed)
            else:
                self.assertEqual(a, b, seed)
            return old[0], b[0]
        self.assertEqual(new, old, seed)
        return old[0], None

    def test_every_pass_gives_the_reference_answer(self):
        seen = {}
        for seed in range(6000):
            r = self._compare(seed)
            seen[r] = seen.get(r, 0) + 1
        # the cases reach every outcome: refused, classes with seqs_order values, a conflict
        self.assertGreater(seen.get(('raised', None), 0), 300, seen)
        self.assertGreater(seen.get(('ok', 'ok'), 0), 300, seen)
        self.assertGreater(seen.get(('ok', 'raised'), 0), 30, seen)

    def test_match_is_symmetric(self):
        # the passes compare a pair from either end (ambiguous_pair's cluster check counts
        # a match at its later group, the keys' pairs are counted from both ends)
        for seed in range(1500):
            rng = random.Random(seed)
            descs, so = case(rng)
            ds = descs + [d for d, _, _ in so]
            known = U._known_paths(ds)
            for a in ds:
                for b in ds:
                    self.assertEqual(U._matches(a, b, known), U._matches(b, a, known), seed)

    def test_the_candidates_hold_every_match(self):
        # _Candidates' rule: every match is among the candidates, once
        for seed in range(3000):
            rng = random.Random(seed)
            descs, so = case(rng)
            ds = descs + [d for d, _, _ in so]
            known = U._known_paths(ds)
            idx = U._Candidates(ds, known)
            for a in ds:
                cand = [i for lst in idx.query(a) for i in lst]
                self.assertEqual(len(cand), len(set(cand)), seed)
                self.assertEqual(len(cand), idx.count(a), seed)
                got = set(cand)
                for i, b in enumerate(ds):
                    if U._matches(a, b, known):
                        self.assertIn(i, got, (seed, a, b))


# ---------------------------------------------------------------- whole union_pairs() calls

def _fork_variant(rng, p):
    """fork's graphlet as a multi-graph server answered it on pair |p| (or a single-graph
    server, |p| naming no graph): in the form of a document stored before the engine named
    pairs by id (paths), after (id), or both."""
    resp = copy.deepcopy(T.doc_json('fork', 'graphlet'))
    for k in ('graph', 'graph_path', 'annotation_path', 'pair_id'):
        resp.pop(k, None)
    if p.get('graph'):
        resp['graph'] = p['graph']
        form = rng.choice(('paths', 'id', 'both'))
        if form in ('paths', 'both'):
            resp['graph_path'], resp['annotation_path'] = p['graph_path'], p['annotation_path']
        if form in ('id', 'both'):
            resp['pair_id'] = configured_pair_id(p['graph'], p['graph_path'],
                                                 p['annotation_path'])
    caps = resp['capabilities']
    resp['release'] = p['release'] if rng.random() < 0.6 else ''
    caps['release'] = resp['release']
    caps['index_ns'] = p['index_ns']
    caps['index_fp'] = p['index_fp'] if rng.random() < 0.5 else None
    caps['index_meta_fp'] = p['index_meta_fp'] if rng.random() < 0.6 else None
    return from_response(resp['results'][0], resp)


class TestWholeUnions(unittest.TestCase):
    """union_pairs() against union_pairs() with the reference passes, on 2,500 random
    calls: every refusal and every union (its digest, pairs, keys, tags, seqs_order values,
    expected and not-included pairs) the same."""

    @classmethod
    def setUpClass(cls):
        rng = random.Random(7)
        cls.pairs = []
        for _ in range(6):
            cls.pairs += universe(rng)
        cls.parts = [(_fork_variant(rng, p), p) for p in cls.pairs[:24]]
        cls.seed = T.graphlet('fork').seed.sequence

    def _inputs(self, rng):
        parts = rng.sample(self.parts, rng.randint(0, 3))
        pairs = [q for _, q in parts] + rng.sample(self.pairs, rng.randint(1, 3))
        wild = rng.random() < 0.1
        results = []
        for _ in range(rng.randint(1, 8)):
            r = rng.random()
            if r < 0.45 and parts:
                g, q = rng.choice(parts)
                stmt = view(rng, q if rng.random() < 0.8 else rng.choice(pairs), wild)
                if rng.random() < 0.4 or stmt is None:
                    results.append(g)          # stated by its envelope
                else:
                    results.append({'pair': stmt, 'graphlet': g})
                continue
            stmt = view(rng, rng.choice(pairs), wild) or rng.choice(GRAPHS)
            if rng.random() < 0.02:
                stmt = {'x': 1}                # not a pair statement: refused
            if r < 0.8:
                results.append({'pair': stmt,
                                'error': rng.choice(('request failed', 'timeout', 'x'))})
            else:
                res = {'seed': {'seed_id': 'x', 'length_bp': len(self.seed)},
                       'outcome': {'walks': 'not_in_graph'}, 'limitations': []}
                results.append({'pair': stmt, 'result': res})
        kw = {}
        if rng.random() < 0.5:
            exp = [view(rng, rng.choice(pairs), wild) for _ in range(rng.randint(0, 5))]
            kw['expected_pairs'] = [d['graph'] if set(d) == {'graph'} and rng.random() < 0.5
                                    else d for d in exp if d is not None]
        if rng.random() < 0.6:
            kw['seqs_order'] = key_list(rng, pairs)
        if rng.random() < 0.3 or not any(not isinstance(x, dict) or 'graphlet' in x
                                         for x in results):
            kw['seed'] = self.seed
        return results, kw

    def _union(self, results, kw):
        def run():
            u = union_pairs(results, **kw)
            return (u.digest, u.pairs, u.expected,
                    [(p.key, p.tag, p.pair_stated_by, p.seqs_order) for p in u.parts],
                    [(a.key, a.tag, a.category, a.seqs_order) for a in u.absent],
                    u._failures, u._not_included)
        return outcome(run)

    def test_every_call_gives_the_reference_answer(self):
        seen = {}
        for seed in range(2500):
            rng = random.Random(10000 + seed)
            results, kw = self._inputs(rng)
            new = self._union(results, kw)
            with with_reference():
                old = self._union(results, kw)
            self.assertEqual(new, old, seed)
            key = new[0] if new[0] == 'ok' else new[2] or new[1]
            seen[key] = seen.get(key, 0) + 1
        for code in ('ok', 'ambiguous_pair', 'duplicate_pair', 'not_expected',
                     'bad_seqs_order'):
            self.assertGreater(seen.get(code, 0), 20, seen)

    def test_the_service_ceiling(self):
        # the service's _ceiling_inputs (tests/test_traversal_graphlet_union_op.py), scaled
        # down and perturbed: one entry per task, its own graph name, index_fp and n1 id
        # (no paths), one release; each a {pair, error} input and an expected pair; one
        # (graph, n1) key per pair. Then the shapes real seeds add: a few names with two or
        # more pairs that differ by id and index_fp, a pair whose row has no path column
        # (no pair_id: matched by name alone), parts in both document forms
        rng = random.Random(3)
        for trial in range(120):
            n = rng.randint(1, 40)
            ents = []
            for i in range(n):
                g = 'g%d' % (i if rng.random() < 0.85 else rng.randrange(max(1, i)))
                ents.append({'graph': g, 'pair_id': configured_pair_id(g, '/p%d' % i, '/a%d' % i),
                             'index_fp': 'fp%d' % i, 'release': 'r237'})
            for e in rng.sample(ents, min(len(ents), rng.randint(0, 3))):
                del e['pair_id']               # no path column: by name alone
                if rng.random() < 0.5:
                    del e['index_fp']
            results = [{'pair': dict(e), 'error': 'request failed'} for e in ents]
            exp = [dict(e) for e in ents]
            so = {(e['graph'], e['pair_id']): 'verified' for e in ents if 'pair_id' in e}
            for g, p in rng.sample(self.parts, rng.randint(0, 3)):
                results.append(g)
                if rng.random() < 0.5:
                    exp.append(U._env_desc(g.envelope, g))
            if rng.random() < 0.3:
                so[rng.choice(ents)['graph']] = rng.choice(('count_only', 'verified'))
            rng.shuffle(results)
            kw = {'expected_pairs': exp, 'seqs_order': so, 'seed': self.seed}
            new = self._union(results, kw)
            with with_reference():
                old = self._union(results, kw)
            self.assertEqual(new, old, trial)


# ---------------------------------------------------------------- the charges

class _Recorder:
    """A stand-in budget: the work charged so far, and each charge in |events|."""

    def __init__(self, events=None):
        self.work = 0
        self.events = events if events is not None else []

    def charge(self, work, nbytes=0):
        self.work += work
        self.events.append(('charge', work))

    def release(self, nbytes):
        pass


class TestCharges(unittest.TestCase):

    def test_each_pass_is_charged_before_it_runs(self):
        """Admitted before it runs (budget.py), read from what each pass did, in order:
        before it groups, its groups are charged; before it builds an index, that index and
        its queries are charged: one query per description the pass then counts, and no
        more scans than that (Opus check of dada7f430, N1: a pass charging no queries
        passed); in the refusal search (_ambiguous()) each group's count follows its own
        query charge (N3); at every comparison it has charged _W_MATCH for at least as
        many; a pass that refuses nothing charges exactly the comparisons it makes (the
        keys' pairs with equal values are charged and not compared)."""
        real, Index, groups, ambiguous = U._matches, U._Candidates, U._groups, U._ambiguous

        class Spy:
            """The index, its counts and scans recorded."""

            def __init__(self, ix, events):
                self.ix, self.events = ix, events

            def count(self, d):
                self.events.append(('count',))
                return self.ix.count(d)

            def query(self, d):
                self.events.append(('scan',))
                return self.ix.query(d)
        for seed in range(3000):
            rng = random.Random(seed)
            descs, so = case(rng)
            known = U._known_paths(descs + [d for d, _, _ in so])
            classes = U._classes(descs, known)
            tags = {id(c): 't%d' % i for i, c in enumerate(classes)}
            for name, run in (
                    ('keys', lambda b: U._check_keys(so, known, b)),
                    ('classes', lambda b: U._pair_classes(descs, known, b)),
                    ('seqs', lambda b: U._seqs_of(classes, so, tags, known, b))):
                events = []
                b = _Recorder(events)

                def index_price(b_, grouped, indexed, queried, held):
                    events.append(('price', grouped, indexed, queried))
                    return 0

                def compare_price(b_, compared, longest):
                    b_.charge(U._W_MATCH * compared)
                    return 0

                def grouping(ds):
                    events.append(('groups', len(ds)))
                    return groups(ds)

                def index(ds, k):
                    events.append(('index', len(ds)))
                    return Spy(Index(ds, k), events)

                def searched(*a):
                    events.append(('ambiguous',))
                    return ambiguous(*a)

                def counted(x, y, k=None):
                    events.append(('match', b.work))
                    return real(x, y, k)
                with mock.patch.multiple(U, _matches=counted, _compare_price=compare_price,
                                         _index_price=index_price, _Candidates=index,
                                         _groups=grouping, _ambiguous=searched):
                    r = outcome(run, b)
                calls, priced, queried = 0, [], None
                counts = scans = 0
                search, last = False, None
                for e in events:
                    if e[0] == 'price':
                        priced.append(e)
                    elif e[0] == 'groups':
                        # charged since the last structure was built
                        self.assertIn(e[1], [p[1] for p in priced], (seed, name, events))
                        priced = []
                    elif e[0] == 'index':
                        self.assertIn(e[1], [p[2] for p in priced], (seed, name, events))
                        # its queries, priced with it
                        queried = sum(p[3] for p in priced)
                        priced = []
                    elif e[0] == 'ambiguous':
                        search = True
                    elif e[0] in ('count', 'scan') and search:
                        if e[0] == 'count':
                            self.assertEqual(last, ('charge', 2 * U._W_QUERY), (seed, name))
                    elif e[0] == 'count':
                        counts += 1
                    elif e[0] == 'scan':
                        scans += 1
                    elif e[0] == 'match':
                        calls += 1
                        self.assertLessEqual(U._W_MATCH * calls, e[1], (seed, name))
                    last = e
                if queried is not None:
                    # every description the pass queries was charged a count and a scan
                    self.assertEqual(counts, queried, (seed, name))
                    self.assertLessEqual(scans, queried, (seed, name))
                if r[0] == 'ok' and name != 'keys':
                    self.assertEqual(U._W_MATCH * calls, b.work, (seed, name))

    def test_each_pass_holds_what_it_is_charged(self):
        """A pass's own account (a budget outside any call: its charges only) bounds the
        memory it holds while it runs beyond what it returns (review M1: a _Paths per
        description where each names its own graph), warm free lists as above; every
        seqs_order key form, by name or by path, distinct names and one name, and the
        description shapes that index most (all fields stated, each its own name)."""
        import gc
        import tracemalloc
        from metagraph.traverse.client import Pair

        def key_forms(i):
            g, gp, ap = 'g%d' % i, '/p%d' % i, '/a%d' % i
            return [g, (g,), (g, configured_pair_id(g, gp, ap)), (g, gp), (g, gp, ap),
                    (g, None, ap), Pair(g, gp, ap, None),
                    Pair(g, None, None, configured_pair_id(g, gp, ap)),
                    ('g', configured_pair_id('g', gp, ap)), ('g', gp), ('g', gp, ap)]

        def desc_forms(i):
            g, gp, ap = 'g%d' % i, '/p%d' % i, '/a%d' % i
            return [{'graph': g, 'pair_id': configured_pair_id(g, '/p', '/a'),
                     'index_fp': 'f%d' % i, 'release': 'r'},
                    {'graph': g, 'graph_path': gp, 'annotation_path': ap, 'index_fp': 'f%d' % i,
                     'index_meta_fp': 'm%d' % i, 'release': 'r', 'index_ns': 'n%d' % i},
                    {'graph': g, 'graph_path': gp, 'annotation_path': ap, 'release': 'r',
                     'index_meta_fp': 'm%d' % i, 'index_ns': 'n%d' % i},
                    {'graph': g, 'graph_path': gp}, {'graph': g},
                    {'index_ns': 'n%d' % i, 'pair_id': configured_pair_id(g, gp, ap)},
                    {'index_ns': 'n%d' % i, 'graph_path': gp}]

        def held(fn):
            fn()
            gc.collect()
            gc.disable()
            try:
                tracemalloc.start()
                out = fn()
                cur, peak = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
                gc.enable()
            del out
            return peak - cur

        def account(fn):
            b = LocalBudget()
            fn(b)
            return b.peak_bytes

        # n = 2,731: past a resize of the index's tables, where a pass holds 1,834-1,867 B a
        # description (1,741 at 2,000)
        n = 2731
        for f in range(len(key_forms(0))):
            so = U._check_args([], None, None,
                               {key_forms(i)[f]: 'verified' for i in range(n)})[1]
            ks = [d for d, _, _ in so]
            known = U._known_paths(ks)
            classes = U._classes(ks, known)
            tags = {id(c): 't' for c in classes}
            for name, fn in (('keys', lambda b=None: U._check_keys(so, known, b)),
                             ('seqs', lambda b=None: U._seqs_of(classes, so, tags, known, b))):
                with self.subTest(form=repr(key_forms(0)[f]), op=name):
                    self.assertLessEqual(held(fn), account(fn))
        for f in range(len(desc_forms(0))):
            ds = [U._pair_desc(desc_forms(i)[f]) for i in range(n)]
            known = U._known_paths(ds)
            with self.subTest(form=desc_forms(0)[f], op='classes'):
                fn = (lambda b=None: U._pair_classes(ds, known, b))
                self.assertLessEqual(held(fn), account(fn))

    def test_the_service_ceiling_is_admitted(self):
        """The service's 2,000-entry seed (each entry an input and an expected pair, one
        (graph, n1) key each: 4,000 descriptions) builds under its 200M lwu union budget,
        at a linear price; so do 2,000 inputs of one name with distinct ids, and 1,000
        inputs naming one pair."""
        def ceiling(n):
            ents = [{'graph': 'g%d' % i, 'pair_id': configured_pair_id('g%d' % i, '/p', '/a'),
                     'index_fp': '%016x' % i, 'release': 'r237'} for i in range(n)]
            return ([{'pair': dict(e), 'error': 'request failed'} for e in ents],
                    {'expected_pairs': [dict(e) for e in ents],
                     'seqs_order': {(e['graph'], e['pair_id']): 'verified' for e in ents}})

        def one_name(n):
            ids = [configured_pair_id('g', '/p%d' % i, '/a') for i in range(n)]
            return ([{'pair': {'graph': 'g', 'pair_id': x}, 'error': 'x'} for x in ids],
                    {'seqs_order': {('g', x): 'verified' for x in ids}})

        def one_pair(n):
            return ([{'pair': 'g', 'error': 'request failed %d' % i} for i in range(n)], {})

        for name, make, n, most in (('ceiling', ceiling, 2000, 2_000_000),
                                    ('one name', one_name, 2000, 1_000_000),
                                    ('one pair', one_pair, 1000, 500_000)):
            with self.subTest(name):
                r, kw = make(n)
                b = LocalBudget(work_units=200_000_000)
                u = union_pairs(r, seed='ACGT', budget=b, **kw)
                self.assertLess(b.used_work, most)
                # linear: twice the inputs, about twice the charge
                r2, kw2 = make(n // 2)
                b2 = LocalBudget(work_units=200_000_000)
                union_pairs(r2, seed='ACGT', budget=b2, **kw2)
                self.assertLess(b.used_work, 2.2 * b2.used_work)
                self.assertEqual(u.digest, union_pairs(r, seed='ACGT', **kw).digest)


if __name__ == '__main__':
    unittest.main()
