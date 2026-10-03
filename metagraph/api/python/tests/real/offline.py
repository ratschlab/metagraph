"""The real-data suite, offline: a snapshot of small real retrievals committed as fixtures
(api/python/tests/data/traverse/real/), checked by the suite's own conformance and
search-oracle code without any server.

    Harness(fixtures)   realdata.py, test_real_conformance.py and test_real_oracle.py loaded
                        as PRIVATE module copies whose harness points at the fixtures:
                        servers off (no socket is ever opened), retrievals from
                        <fixtures>/<index>/<cell>.{request,graphlet,full}.json.gz, the
                        recorded /resolve answers from <fixtures>/<index>/<cell>.oracle.json.gz
                        only, the fixtures' manifest, catalog and capabilities. A /resolve
                        payload without a recorded answer is a FAILURE (OfflineMiss), not a
                        skip: it means the oracle plan changed, i.e. the library's output did.
    snapshot()          (re)build the fixtures from the retrieval cache of the live suite
                        (realdata.py fill) plus the DERIVED cells, which are fetched live; the
                        oracle answers come from the oracle cache, missing ones from the live
                        servers. Every check is then run offline on the new fixtures and its
                        outcome (ran / skipped and why) is recorded in the manifest; a failing
                        check aborts the snapshot and leaves the old fixtures in place.

    python3 offline.py snapshot     # needs the cache; the servers for derived cells / misses
    python3 offline.py check        # run every recorded check offline, print the table
    python3 offline.py situations   # what the fixture bodies cover (the README's claims)

The unit test api/python/tests/test_traverse_real_offline.py is the CI side of this: it
generates one test per (fixture cell, recorded check) and turns a skip of a check that ran
when the fixtures were made into a failure.
"""

import argparse
import collections
import copy
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
API = os.path.dirname(TESTS)
FIXTURES = os.path.join(TESTS, 'data', 'traverse', 'real')
FORMAT = 1
MAX_FILE_BYTES = 60 * 1024          # every fixture file (gzip) stays below this
if API not in sys.path:                 # metagraph.traverse, for the private module copies
    sys.path.insert(0, API)

# ------------------------------------------------------------------ the selection

# (index, cell, why): cells copied from the retrieval cache
CELLS = [
    ('sra', 'sra_hub__strict',
     'primary regime, the default strategy on a 50-mer carried by 508 SRA runs: 508 derived '
     'seed labels, merges (united history), runs whose displayed support starts at a merge'),
    ('sra', 'sra_rand50_00__limit1',
     'primary regime, branch limit 1 on a random gut-genome 50-mer: merges, merge-closed (m) '
     "runs, and runs cloned at a split that carry their source's route_bp (8 runs that "
     'the server fix of 2026-10-03 changed from 0)'),
    ('sra', 'sra_rand50_07__switch1',
     'primary regime, one switch (E s on the right arm): a run entered by switch that ends '
     'at the loss budget (B), the silent Lw end of its source, a superseded (Ls) end'),
    ('sra', 'sra_rand50_00__annotate_merge',
     'primary regime, annotate mode with on_reconverge merge: P runs, merges with one empty '
     'partition per parent, the section 5.1 union rule (b_maximal)'),
    ('sra', 'sra_rand50_00__cut_lists',
     'cut lists: annotate with max_labels_per_node 2, label_evidence lower_bound; /resolve '
     'confirms each cut list is a subset carrying the true count'),
    ('sra', 'sra_rand50_01__cut_events',
     'cut evidence: max_branch_events 1, branch_diagnostics cut with a branch_events '
     'limitation'),
    ('sra', 'sra_16s_OK510341__time_tight',
     'a partial, time-budgeted result: time_budget_ms 200 stops both arms (M ends, '
     'cap_trigger time_budget), walks partial; the graphlet and the full fetch ran '
     'separately, so each is checked on its own'),
    ('sra', 'sra_batch3__batch_refuse',
     'a 3-seed batch with max_seed_labels 3: the middle seed is refused per seed '
     '({seed, error}, no body), the other two are walked'),
    ('uhgg', 'uhgg_conserved__limit2',
     'basic regime, a 100 bp window carried by 7 genomes, branch limit 2 with merge: merges, '
     'merge-closed runs, routes through non-first parents, 21 split clones carrying their '
     "source's route_bp"),
    ('uhgg', 'uhgg_contig_end__annotate_merge',
     'basic regime, annotate with merge on the last 100 bp of a contig: the genome labels '
     'recorded past the contig end'),
    ('mini_refseq', 'mini_ndm1__trace',
     'header (accession) labels from the .seqs sidecar under support trace '
     '(coordinate-consecutive), keep; verifiable index identity (index_fp from the '
     'manifest); edge-reuse (U) ends'),
    ('mini_refseq', 'mini_ndm1_rc__switch1',
     'the reverse complement of blaNDM-1 in a basic graph (the other strand): switches '
     'between header labels, routes through merges'),
    ('mini_refseq', 'mini_ndm1__annotate_exh',
     'annotate exhaustive on header labels: the recorded sets against /resolve discover '
     '(kind header)'),
]

# cells NOT in the cache matrix, fetched live by snapshot() (a situation the cache has only
# in a body far over MAX_FILE_BYTES, or not at all)
DERIVED = [
    {'index': 'sra', 'seed': 'sra_hairpin', 'strategy': 'left_only',
     'why': 'primary regime, a skipped hairpin (E h .. s) and its dead_end/hairpin (Dh) run '
            'ends, on a 50-mer cut 8 bp before the first hairpin of the 2.6 MB body '
            'sra_16s_PZ326290__beam20 (a poly-A/poly-T run)'},
    {'index': 'sra', 'seed': 'sra_hairpin', 'strategy': 'left_only_follow',
     'why': 'primary regime, the same seed with hairpins: follow: the followed hairpin '
            '(E h .. f), the reverse-complement retrace ending edge_reuse_rc (V)'},
    {'index': 'uhgg', 'seed': 'uhgg_contig_end', 'strategy': 'switch1_limit1_keep',
     'why': 'basic regime, the last 100 bp of contig MGYG000000041_1 with one switch allowed '
            '(the genome that continues past the contig end as labels.extra), branch limit '
            '1, keep: the lineages switch into that genome past the end'},
]

HAIRPIN_SOURCE = ('sra', 'sra_16s_PZ326290__beam20')
HAIRPIN_OFFSET = 8           # bases from the hairpin node into the molecule
HAIRPIN_LENGTH = 50


def _derived_strategies(R):
    S = R.Strategy
    return [
        S('left_only_follow', 'strict, direction left, branching.hairpins follow (a followed '
          'hairpin: the reverse-complement retrace ends edge_reuse_rc)',
          lambda i, s: {'direction': 'left', 'branching': {'hairpins': 'follow'}},
          tags=('constrain', 'direction', 'hairpin_follow')),
        S('switch1_limit1_keep', 'switch1 (constant cost 1, loss budget 1, labels.extra = the '
          "seed's switch targets) with branch limit 1 and on_reconverge keep",
          lambda i, s: dict(R._switch(i, s), branching={'max_label_branches': 1,
                                                          'on_reconverge': 'keep'}),
          tags=('constrain', 'switch', 'keep')),
    ]


# ------------------------------------------------------------------ loading private copies

class OfflineMiss(AssertionError):
    """The offline harness was asked for something only a server could answer."""


_SERIAL = [0]


def _load_module(stem, realdata=None):
    """tests/real/<stem>.py executed as a private module (a fresh name per load), with
    |realdata| bound as its 'realdata' import. sys.path and sys.modules['realdata'] are
    restored afterwards, so a live run of the suite in the same process is unaffected."""
    _SERIAL[0] += 1
    name = '_mgt_offline_%s_%d' % (stem, _SERIAL[0])
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, stem + '.py'))
    mod = importlib.util.module_from_spec(spec)
    saved_path = list(sys.path)
    had = 'realdata' in sys.modules
    saved_rd = sys.modules.get('realdata')
    sys.modules[name] = mod
    try:
        if realdata is not None:
            sys.modules['realdata'] = realdata
        spec.loader.exec_module(mod)
    finally:
        sys.path[:] = saved_path
        if had:
            sys.modules['realdata'] = saved_rd
        else:
            sys.modules.pop('realdata', None)
    return mod


def _gz_json(path):
    with gzip.open(path, 'rb') as f:
        return json.loads(f.read())


class Harness:
    """The suite's harness and check modules pointed at |fixtures| (see the module doc)."""

    def __init__(self, fixtures=FIXTURES):
        self.fixtures = fixtures
        if not os.path.exists(os.path.join(fixtures, 'manifest.json')):
            raise FileNotFoundError('no fixture manifest in %s (python3 %s snapshot)'
                                    % (fixtures, os.path.relpath(__file__)))
        self.R = R = _load_module('realdata')
        # a private temporary directory holds what the harness reads as plain files: the
        # catalog (stored gzipped) and the oracle answers (stored as one bundle per cell)
        self._tmp = tempfile.TemporaryDirectory(prefix='mgt_offline_')
        R.SCRATCH = R.SERVERS_DIR = R.CACHE_DIR = fixtures
        R.CAPABILITIES_PATH = os.path.join(fixtures, 'capabilities.json')
        R.MANIFEST_PATH = os.path.join(fixtures, 'manifest.json')
        R.CATALOG_PATH = os.path.join(self._tmp.name, 'catalog.json')
        R.ORACLE_DIR = os.path.join(self._tmp.name, 'oracle')
        with open(R.CATALOG_PATH, 'w') as f:
            json.dump(_gz_json(os.path.join(fixtures, 'catalog.json.gz')), f)
        R.SERVERS = collections.OrderedDict((i, None) for i in R.INDEXES)
        R._UP.clear()
        R._MANIFEST_MEMO.update(mtime=None, data=None)

        def no_server(index, *a, **kw):
            raise OfflineMiss(
                '%s: offline, and the fixtures hold no recorded answer for this request: the '
                "oracle plan changed (the library's output changed) or the fixtures are "
                'stale; rebuild them with tests/real/offline.py snapshot while the servers '
                'run' % index)

        def no_http(index, method, route, *a, **kw):
            raise OfflineMiss('%s %s %s: no network access offline' % (index, method, route))
        R.require_server = no_server
        R.http = no_http
        R.server_up = lambda index, *, refresh=False: False
        self.rows = R.cells(heavy=None)
        self.oracle_answers = 0
        for row in self.rows:
            self.oracle_answers += self._materialise_oracle(row)
        self.C = _load_module('test_real_conformance', R)
        self.O = _load_module('test_real_oracle', R)

    def _materialise_oracle(self, row):
        """The cell's recorded /resolve answers, written where realdata.oracle_resolve
        looks for its disk cache (a private temporary directory)."""
        path = bundle_path(self.fixtures, row['index'], row['cell'])
        if not os.path.exists(path):
            return 0
        answers = _gz_json(path)['answers']
        for key, entry in answers.items():
            self.R._write_json_gz(os.path.join(self.R.ORACLE_DIR, row['index'], key[:2],
                                               key + '.json.gz'), entry)
        return len(answers)

    def manifest(self):
        return self.R.load_manifest()

    def conformance_case(self):
        return conformance_base(self)('runTest')

    def oracle_case(self):
        return oracle_base(self)('runTest')

    def checks(self, row):
        """[(kind, check)] that apply to |row|: every conformance check (plus
        full_invariants on a non-deterministic cell) and the oracle areas of
        test_real_oracle._applies."""
        out = [('conformance', c) for c in self.C.CHECKS]
        if not row.get('deterministic'):
            out.append(('conformance', 'full_invariants'))
        out += [('oracle', a) for a in self.O.AREAS if self.O._applies(a, row)]
        return out


def conformance_base(h):
    """The conformance TestCase base of harness |h| (bound to its private modules, kept on
    the harness) with the offline-only check."""
    if '_conformance_base' not in h.__dict__:
        class OfflineConformance(h.C._ConformanceBase):
            def runTest(self):
                pass

            def check_full_invariants(self, row):
                """A non-deterministic cell (time budget): its full fetch stopped elsewhere
                than its graphlet fetch, so the two are not compared, but the full JSON's
                own invariants hold."""
                self.check_invariants(dict(row, deterministic=True))
        h._conformance_base = OfflineConformance
    return h._conformance_base


def oracle_base(h):
    """The oracle TestCase base of harness |h|."""
    if '_oracle_base' not in h.__dict__:
        class OfflineOracle(h.O._OracleBase):
            def runTest(self):
                pass
        h._oracle_base = OfflineOracle
    return h._oracle_base


def write_gz_json(path, obj):
    """Compact JSON, gzip with a fixed header (no file name, mtime 0): the same content
    gives the same bytes, so a rebuild without changes leaves git clean."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = json.dumps(obj, separators=(',', ':')).encode()
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw,
                                                 compresslevel=9, mtime=0) as f:
        f.write(data)


def bundle_path(fixtures, index, cell):
    return os.path.join(fixtures, index, cell + '.oracle.json.gz')


def run_check(case, check, row):
    """('ran', None) | ('skipped', reason); a failure raises."""
    try:
        getattr(case, 'check_' + check)(row)
    except unittest.SkipTest as e:
        return 'skipped', str(e)
    return 'ran', None


# ------------------------------------------------------------------ situations (coverage)

def situations(h):
    """{situation: sorted [index/cell]} found in the fixture BODIES (raw records and the
    JSON summaries; not the manifest's notes): what the README promises the fixtures cover."""
    R = h.R
    found = collections.defaultdict(set)
    seeds = {}
    for row in h.rows:
        cell = R.load_cell(row['index'], row)
        key = '%s/%s' % (row['index'], row['cell'])
        req = cell.request
        extra = set((req['strategy'].get('labels') or {}).get('extra') or ())
        for i, res in enumerate(cell.graphlet_response['results']):
            if 'graphlet' not in res:
                found['derivation_error'].add(key)
                continue
            seeds.setdefault(row['index'], {}).setdefault(res['seed']['length_bp'], set()).add(
                (cell.seed_sequence(i), key))
            for side, a in res['arms'].items():
                cap = a.get('cap_trigger') or {}
                if cap.get('reason') == 'time_budget':
                    found['time_budget'].add(key)
                if (a.get('labels_per_node') or {}).get('nodes_truncated'):
                    found['cut_list_truncated_nodes'].add(key)
            if res['outcome']['walks'] == 'partial':
                found['walks_partial'].add(key)
            if res['outcome']['branch_diagnostics'] == 'cut':
                found['cut_events'].add(key)
            if res['outcome']['label_evidence'] == 'lower_bound':
                found['label_evidence_lower_bound'].add(key)
            text = res['graphlet']
            for line in text.split('\n'):
                f = line.split(' ')
                t = f[0]
                if t == 'H':
                    found['regime_' + f[4]].add(key)
                    found['mode_' + {'c': 'constrain', 'a': 'annotate'}[f[6]]].add(key)
                    found['support_' + {'k': 'kmer', 't': 'trace'}[f[7]]].add(key)
                    found['reconverge_' + {'m': 'merge', 'k': 'keep'}[f[8]]].add(key)
                    found['identity_' + ('unverifiable' if f[14] == '*' else 'verifiable')].add(key)
                elif t == 'L':
                    found['labels_' + {'c': 'column', 'h': 'header'}[f[1]]].add(key)
                elif t == 'G' and ',' in f[1]:
                    found['merge'].add(key)
                elif t == 'E' and f[2] == 'h':
                    found['hairpin_' + {'f': 'followed', 's': 'skipped'}[f[-1]]].add(key)
                elif t == 'E' and f[2] == 's':
                    found['switch'].add(key)
                elif t == 'R':
                    end = f[5]
                    if end == 'Dh':
                        found['end_dead_end_hairpin'].add(key)
                    elif end == 'V':
                        found['end_edge_reuse_rc'].add(key)
                    elif end == 'U':
                        found['end_edge_reuse'].add(key)
                    elif end == 'm':
                        found['merge_closed_run'].add(key)
                    elif end == 'B':
                        found['end_loss_budget'].add(key)
                elif t == 'K' and f[2] == 'label_lists':
                    found['cut_list'].add(key)
            if extra:
                g = cell.graphlet(i)
                for a in g.arms.values():
                    for s in a.segments:
                        for e in s.events:
                            if e.type == 'switch' and g.labels[e.to_label].name in extra:
                                found['switch_into_extra'].add(key)
            if res['seed']['num_seed_labels'] >= 5 and row['index'] == 'uhgg':
                found['multi_genome_seed'].add(key)
    for index, by_len in seeds.items():
        for n, items in by_len.items():
            seqs = {s for s, _ in items}
            for s, key in items:
                if R.revcomp(s) in seqs and R.revcomp(s) != s:
                    found['reverse_complement_seed_pair'].add(key)
    return {k: sorted(v) for k, v in found.items()}


# ------------------------------------------------------------------ snapshot

def _live():
    """The live suite's harness (scratch cache, oracle cache, servers) and its oracle
    module, imported normally."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import realdata as R
    import test_real_oracle as O
    return R, O


def _redirect(R, **kw):
    old = {k: getattr(R, k) for k in kw}
    for k, v in kw.items():
        setattr(R, k, v)
    R._MANIFEST_MEMO.update(mtime=None, data=None)
    return old


def _hairpin_seed(R, fixtures, log):
    """The derived hairpin seed: from the fixtures' catalog when present (stable), else cut
    from HAIRPIN_SOURCE: the molecule of the left arm spelled to the first segment (id
    order) that records a hairpin event, bases [8, 58) from its far end, so the hairpin
    node sits at depth 8 of the new seed's left arm."""
    try:
        for s in _gz_json(os.path.join(fixtures, 'catalog.json.gz'))['indexes']['sra']:
            if s['name'] == 'sra_hairpin':
                return s
    except (OSError, KeyError, ValueError):
        pass
    index, cell = HAIRPIN_SOURCE
    log('deriving the hairpin seed from %s/%s' % HAIRPIN_SOURCE)
    g = R.load_cell(index, cell).graphlet(0)
    a = g.arms['left']
    seg = next(s for s in a.segments if any(e.type == 'hairpin' for e in s.events))
    ev = next(e for e in seg.events if e.type == 'hairpin')
    mol = R.spell_with_seed(g, 'left', seg.id)
    seq = mol[HAIRPIN_OFFSET:HAIRPIN_OFFSET + HAIRPIN_LENGTH]
    res = R.oracle_profile(index, seq)
    return R._seed_record(
        index, 'sra_hairpin', 'hairpin', 'core', seq,
        {'derived_from': '%s/%s' % HAIRPIN_SOURCE, 'arm': 'left', 'segment': seg.id,
         'hairpin_at_bp': ev.at_bp, 'molecule_bp': len(mol),
         'start': HAIRPIN_OFFSET, 'end': HAIRPIN_OFFSET + HAIRPIN_LENGTH}, res,
        'a 50-mer whose left arm meets the first hairpin of %s/%s at depth %d'
        % (index, cell, HAIRPIN_OFFSET))


# what a recorded /resolve answer leaves out: the oracle checks never read 'candidates' (the
# bulk of an answer); the request keeps its labels / discover / support, and the queried
# sequence only as its length and sha256 (the answer's key already binds the full payload)
SLIMMED = ['response.candidates', 'request.sequence (kept: sequence_bp, sequence_sha256)']


def _slim(entry):
    req = dict(entry['request'])
    seq = req.pop('sequence', '')
    req['sequence_bp'] = len(seq)
    req['sequence_sha256'] = hashlib.sha256(seq.encode()).hexdigest()
    resp = {k: v for k, v in entry['response'].items() if k != 'candidates'}
    return {'request': req, 'response': resp}


def _copy_cell(src_dir, dst_dir, index, cell):
    for what in ('request', 'graphlet', 'full'):
        src = os.path.join(src_dir, index, '%s.%s.json.gz' % (cell, what))
        dst = os.path.join(dst_dir, index, '%s.%s.json.gz' % (cell, what))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)


def snapshot(fixtures=FIXTURES, *, log=print, refetch=False):
    R, O = _live()
    R.require_cache()
    for st in _derived_strategies(R):
        R.STRATEGIES.setdefault(st.name, st)
    build = os.path.join(R.SCRATCH, 'offline_build')
    staging = os.path.join(R.SCRATCH, 'offline_derived')
    if os.path.exists(build):
        shutil.rmtree(build)
    os.makedirs(build)
    rows, seeds = collections.OrderedDict(), collections.OrderedDict()

    def add_seed(index, rec):
        seeds.setdefault(index, collections.OrderedDict())[rec['name']] = rec

    for index, cell, why in CELLS:
        meta = R.cell_meta(index, cell)
        if not meta or meta.get('status') != 'ok':
            raise RuntimeError('cache cell %s/%s is missing or not ok' % (index, cell))
        _copy_cell(R.CACHE_DIR, build, index, cell)
        row = copy.deepcopy(meta)
        row['why'] = why
        rows['%s/%s' % (index, cell)] = row
        add_seed(index, R.catalog_seed(index, meta['seed']))

    derived_seeds = {}
    for n, d in enumerate(DERIVED):
        index = d['index']
        if d['seed'] not in derived_seeds:
            derived_seeds[d['seed']] = (_hairpin_seed(R, fixtures, log)
                                        if d['seed'] == 'sra_hairpin'
                                        else R.catalog_seed(index, d['seed']))
        rec = derived_seeds[d['seed']]
        add_seed(index, rec)
        cid = R.cell_id(rec['name'], d['strategy'])
        old = _redirect(R, CACHE_DIR=staging,
                        MANIFEST_PATH=os.path.join(staging, 'manifest.json'))
        try:
            meta = R.cell_meta(index, cid)
            want = R.build_request(index, rec, d['strategy'])
            have = None
            if meta and meta.get('status') == 'ok' and not refetch:
                have = json.loads(R._read_gz_bytes(os.path.join(
                    staging, index, cid + '.request.json.gz')))
            if have != want:
                log('fetching the derived cell %s/%s live' % (index, cid))
                R.require_server(index)
                meta = R.fetch_cell(index, rec, d['strategy'], seed_order=100 + n)
                R._save_manifest([meta])
            if meta.get('status') != 'ok':
                raise RuntimeError('derived cell %s/%s: %s' % (index, cid, meta.get('error')))
        finally:
            _redirect(R, **old)
        _copy_cell(staging, build, index, cid)
        row = copy.deepcopy(meta)
        row['why'] = d['why']
        row['derived'] = {'seed': rec['name'], 'strategy': d['strategy'],
                          'strategy_doc': R.STRATEGIES[d['strategy']].doc,
                          'seed_source': rec.get('source')}
        rows['%s/%s' % (index, cid)] = row

    scratch_manifest = R.load_manifest()
    caps = {}                       # per index, without the local process facts
    for index in R.INDEXES:
        rec = R.recorded_capabilities(index)
        if rec:
            caps[index] = {k: v for k, v in rec.items()
                           if k not in ('pid', 'log', 'host', 'port', 'binary')}
    with open(os.path.join(build, 'capabilities.json'), 'w') as f:
        json.dump(caps, f, indent=1, sort_keys=True)
    write_gz_json(os.path.join(build, 'catalog.json.gz'), {
        'format': 1, 'note': 'the seeds of the offline fixture cells (a subset of the live '
        'suite catalog, plus derived seeds)',
        'indexes': {i: list(s.values()) for i, s in seeds.items()}})
    manifest = collections.OrderedDict([
        ('format', FORMAT),
        ('description', 'offline regression fixtures: small real retrievals of the live '
                        'real-data suite (tests/real) with the /resolve answers their oracle '
                        'checks need; built by tests/real/offline.py snapshot'),
        ('servers', {i: (scratch_manifest.get('servers') or {}).get(i)
                     for i in sorted({r['index'] for r in rows.values()})}),
        ('cells', rows),
    ])

    def write_manifest():
        with open(os.path.join(build, 'manifest.json'), 'w') as f:
            json.dump(manifest, f, indent=1)
    write_manifest()

    # the /resolve answers every oracle plan of every cell asks for (misses go live)
    for index in {r['index'] for r in rows.values()}:
        live, rec = R.identity(index), (caps.get(index) or {}).get('capabilities') or {}
        for k in ('index_ns', 'index_fp', 'index_meta_fp', 'k'):
            if live.get(k) != rec.get(k):
                raise RuntimeError('%s: live %s %r != recorded %r: the oracle keys would not '
                                   'match offline' % (index, k, live.get(k), rec.get(k)))
    old = _redirect(R, CACHE_DIR=build, MANIFEST_PATH=os.path.join(build, 'manifest.json'))
    orig_key = R._oracle_key
    seen = []

    def recording_key(index, payload):
        key = orig_key(index, payload)
        seen.append((index, key))
        return key
    R._oracle_key = recording_key
    try:
        for key, row in rows.items():
            del seen[:]
            O._PlanCache.key = None
            t0 = time.time()
            for p in O._PlanCache.get(row):
                if p is not None:
                    p.run()
            answers = collections.OrderedDict()
            for index, k in seen:
                path = os.path.join(R.ORACLE_DIR, index, k[:2], k + '.json.gz')
                answers[k] = _slim(R._read_json_gz(path))
            write_gz_json(bundle_path(build, row['index'], row['cell']),
                          {'format': 1, 'index': row['index'], 'cell': row['cell'],
                           'slimmed': SLIMMED, 'answers': answers})
            log('  %-50s %3d /resolve answers (%.1f s)' % (key, len(answers), time.time() - t0))
    finally:
        R._oracle_key = orig_key
        _redirect(R, **old)
        O._PlanCache.key = None

    # every check, offline, on the new fixtures: what ran is what CI will require
    h = Harness(build)
    cc, oc = h.conformance_case(), h.oracle_case()
    failures = []
    for key, row in rows.items():
        hrow = h.R.cell_meta(row['index'], row['cell'])
        ran, skipped = {'conformance': [], 'oracle': []}, {}
        for kind, check in h.checks(hrow):
            try:
                status, why = run_check(cc if kind == 'conformance' else oc, check, hrow)
            except Exception as e:                          # noqa: BLE001 -- reported below
                failures.append('%s %s: %s: %s' % (key, check, type(e).__name__, str(e)[:300]))
                continue
            if status == 'ran':
                ran[kind].append(check)
            else:
                skipped[check] = why
        row['expect'] = {'ran': ran, 'skipped': skipped}
    if failures:
        raise RuntimeError('checks fail on the new fixtures (the old ones are kept):\n  '
                           + '\n  '.join(failures))
    write_manifest()

    big = []
    for dirpath, _, files in os.walk(build):
        for fn in files:
            p = os.path.join(dirpath, fn)
            if fn.endswith('.gz') and os.path.getsize(p) >= MAX_FILE_BYTES:
                big.append((os.path.relpath(p, build), os.path.getsize(p)))
    if big:
        raise RuntimeError('fixture files of %d bytes or more: %s' % (MAX_FILE_BYTES, big))
    if os.path.exists(fixtures):
        shutil.rmtree(fixtures)
    shutil.copytree(build, fixtures)
    shutil.rmtree(build)
    total = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(fixtures)
                for f in fs)
    log('wrote %d cells to %s (%d bytes)' % (len(rows), fixtures, total))
    return manifest


# ------------------------------------------------------------------ check (CLI)

def check(fixtures=FIXTURES, log=print):
    h = Harness(fixtures)
    cc, oc = h.conformance_case(), h.oracle_case()
    n_ran = n_skip = 0
    bad = []
    for row in h.rows:
        exp = row.get('expect') or {}
        t0 = time.time()
        got = []
        for kind in ('conformance', 'oracle'):
            for c in (exp.get('ran') or {}).get(kind) or []:
                try:
                    status, why = run_check(cc if kind == 'conformance' else oc, c, row)
                except Exception as e:                      # noqa: BLE001 -- reported
                    bad.append('%s/%s %s: %s' % (row['index'], row['cell'], c, e))
                    got.append('!' + c)
                    continue
                if status != 'ran':
                    bad.append('%s/%s %s now skips: %s' % (row['index'], row['cell'], c, why))
                    got.append('-' + c)
                    continue
                n_ran += 1
                got.append(c)
        n_skip += len(exp.get('skipped') or {})
        log('%-12s %-42s %5.2f s  %s' % (row['index'], row['cell'], time.time() - t0,
                                          ' '.join(got)))
    log('%d checks ran, %d recorded skips, %d problems; %d /resolve answers'
        % (n_ran, n_skip, len(bad), h.oracle_answers))
    for b in bad:
        log('  ' + b)
    return not bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('cmd', choices=['snapshot', 'check', 'situations'])
    ap.add_argument('--fixtures', default=FIXTURES)
    ap.add_argument('--refetch', action='store_true',
                    help='fetch the derived cells again even when staged')
    a = ap.parse_args(argv)

    def log(*x):
        print(*x, flush=True)
    if a.cmd == 'snapshot':
        snapshot(a.fixtures, log=log, refetch=a.refetch)
        return 0 if check(a.fixtures, log=log) else 1
    if a.cmd == 'check':
        return 0 if check(a.fixtures, log=log) else 1
    for k, v in sorted(situations(Harness(a.fixtures)).items()):
        log('%-30s %s' % (k, ', '.join(v)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
