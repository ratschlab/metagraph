import copy
import json
import os
import random
import re
import shlex
import subprocess
import sys
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
try:
    import metagraph.traverse as graphlet_lib
except ImportError:
    sys.path.insert(0, os.path.join(REPO, 'api', 'python'))
    import metagraph.traverse as graphlet_lib
from metagraph.traverse.export import normalize_result  # noqa: E402

"""
End-to-end tests of the label-consistent traversal on a protein graph
(SPEC-labeled-traversal-core.md §6.12): `metagraph traverse [--resolve]` and the server
routes /resolve, /traverse and /traverse/capabilities on an index of UniParc's shape, a
Protein build's basic succinct graph at k = 15 with a coordinate-aware annotation
(row_diff_brwt_coord) and a CoordToHeader (.seqs) sidecar, so labels are protein headers.

Sequence layout (random residues, fixed seed), two families in two files (two columns):
    family_a.fa  prot1 = N1 + DOMAIN + C1
                 prot2 = N1 + DOMAIN + C2   (same N-terminus, another C-terminus)
    family_b.fa  prot3 = N2 + DOMAIN + C1   (another N-terminus, prot1's C-terminus)
The proteins share every 15-mer of the blocks they share, so the graph branches where
they diverge; a label-consistent walk reads one protein (or one family) through the
segments it shares with the others. A protein graph has one strand: no reverse
complement, no canonical or primary mode.

test_traverse.py's fixtures run in a Protein build too, as protein graphs over the residues
A, C, G and T; this module is the protein graph over amino acids, at UniParc's k and forms.
"""

K = 15
BLOCK = 60
RESIDUES = 'ACDEFGHIKLMNPQRSTVWY'


def _supports_traverse():
    res = subprocess.run([METAGRAPH, 'traverse'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = (res.stdout + res.stderr).decode()
    return 'Usage' in out and 'traverse' in out


@unittest.skipUnless(PROTEIN_MODE, "the fixture is a protein graph (a Protein build)")
@unittest.skipUnless(_supports_traverse(), "`metagraph traverse` is not available in this build")
class TestTraverseProtein(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        rng = random.Random(20261031)

        def seq(n):
            return ''.join(rng.choice(RESIDUES) for _ in range(n))

        cls.n1, cls.n2 = seq(BLOCK), seq(BLOCK)
        cls.domain = seq(BLOCK)
        cls.c1, cls.c2 = seq(BLOCK), seq(BLOCK)
        cls.records = {
            'prot1': cls.n1 + cls.domain + cls.c1,
            'prot2': cls.n1 + cls.domain + cls.c2,
            'prot3': cls.n2 + cls.domain + cls.c1,
        }
        cls.files = [cls.tempdir.name + '/family_a.fa', cls.tempdir.name + '/family_b.fa']
        for path, names in zip(cls.files, (('prot1', 'prot2'), ('prot3',))):
            with open(path, 'w') as f:
                for name in names:
                    f.write(f'>{name}\n{cls.records[name]}\n')

        cls.graph = cls.tempdir.name + '/graph_k15.dbg'
        cls.anno_base = cls.tempdir.name + '/annotation'
        cls.anno = cls.anno_base + '.row_diff_brwt_coord.annodbg'
        cls._build_graph(cls.files, cls.graph, K, 'succinct', mode='basic')
        cls._annotate_graph(cls.files, cls.graph, cls.anno_base, 'row_diff_brwt_coord',
                            anno_type='filename')
        # The CoordToHeader sidecar, next to the annotation's basename (the loader strips the
        # whole extension). It is written one column per input file in the order the files
        # are given, while a parallel annotation orders its columns as its threads finish
        # them: the files are therefore listed in the annotation's own column order.
        res = cls._run_command(f'{METAGRAPH} stats --print-col-names -a {cls.anno}', 'stats')
        lines = res.stdout.decode().split('\n')
        columns = lines[lines.index('Number of columns: 2') + 1:][:2]
        assert sorted(columns) == sorted(cls.files), columns
        cls.columns = columns
        cls._run_command(
            f'{METAGRAPH} annotate -i {cls.graph} --anno-filename --index-header-coords '
            f'-o {cls.anno_base} {" ".join(columns)}',
            'index-header-coords')
        assert os.path.exists(cls.anno_base + '.seqs'), 'CoordToHeader sidecar was not created'

        # the server, on a port of its own choosing (it states it in its log)
        cls.log = tempfile.NamedTemporaryFile(mode='w+', dir=cls.tempdir.name, delete=False)
        cls.server = subprocess.Popen([METAGRAPH, 'server_query', '-i', cls.graph, '-a', cls.anno,
                                       '--address', '127.0.0.1', '--port', '0', '-p', '2', '-v'],
                                      stdout=cls.log, stderr=cls.log)
        deadline = time.monotonic() + 60
        cls.url = None
        while time.monotonic() < deadline and cls.server.poll() is None:
            with open(cls.log.name) as f:
                match = re.search(r'\[Server\] Listening on port (\d+)', f.read())
            if match:
                cls.url = f'http://127.0.0.1:{match.group(1)}'
                break
            time.sleep(0.2)
        if cls.url is None:
            cls.server.kill()
            with open(cls.log.name) as f:
                raise AssertionError('the server did not start:\n' + f.read()[-2000:])
        cls.session = requests.Session()
        cls.session.trust_env = False

    @classmethod
    def tearDownClass(cls):
        cls.server.kill()
        cls.server.wait()
        cls.session.close()
        cls.log.close()
        super().tearDownClass()

    @classmethod
    def _traverse(cls, request, resolve=False, json_flag=False):
        """The CLI on one request dict; -> (the parsed JSON, the return code, the raw bytes)."""
        path = cls.tempdir.name + '/request.json'
        with open(path, 'w') as f:
            json.dump(request, f)
        cmd = [METAGRAPH, 'traverse'] + (['--resolve'] if resolve else []) \
            + (['--json'] if json_flag else []) + ['-i', cls.graph, '-a', cls.anno, path]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = res.stdout.decode()
        assert out.strip(), f'no output\nstderr: {res.stderr.decode()}'
        return json.loads(out), res.returncode, res.stdout

    def _post(self, route, payload):
        return self.session.post(f'{self.url}/{route}', data=json.dumps(payload), timeout=60)

    @staticmethod
    def _names(result):
        return {i: l['name'] for i, l in enumerate(result['label_dict'])}

    @staticmethod
    def _flank(arm, path, side):
        """The residues one path of an arm spells, in natural orientation: segments are
        reported in natural orientation (SPEC §6.2), so a left path's segments, listed from
        the seed outwards, are joined leaf first."""
        segments = {s['id']: s for s in arm['segments']}
        order = path['segments'] if side == 'right' else reversed(path['segments'])
        return ''.join(segments[sid]['sequence'] for sid in order)

    def _leaves(self, result, side):
        """-> {flank: the names of the labels that travelled it to its end} of one arm."""
        names = self._names(result)
        arm = result['arms'][side]
        return {self._flank(arm, path, side): {names[e['label']] for e in path['end_labels']
                                               if e['route_bp'] == 0}
                for path in arm['paths']}

    def _column_names(self):
        """The two column labels (the family files), in column order, as a derivation
        of the column kind states them."""
        out, rc, _ = self._traverse({'seeds': [{'sequence': self.domain}],
                                     'strategy': {'direction': 'right',
                                                  'bounds': {'max_extension_bp': 1},
                                                  'labels': {'seed_label_kind': 'column'}}})
        self.assertEqual(0, rc, out.get('error'))
        return sorted(l['name'] for l in out['results'][0]['label_dict'])

    def test_capabilities_state_the_alphabet(self):
        """The probe and every response state the graph's alphabet, which tells a client
        that the unit of every *_bp length is a residue; one strand: the basic regime,
        trace support, no reverse complement anywhere."""
        caps = self.session.get(f'{self.url}/traverse/capabilities', timeout=30).json()
        self.assertEqual((K, 'basic'), (caps['k'], caps['regime']))
        alphabet = caps['alphabet']
        for residue in RESIDUES:
            self.assertIn(residue, alphabet)
        self.assertTrue(caps['has_coordinates'])
        self.assertTrue(caps['has_coord_to_header'])
        self.assertTrue(caps['supports_trace'])
        out = self._post('traverse', {'seeds': [{'sequence': self.domain}],
                                      'strategy': {'bounds': {'max_extension_bp': 5}}}).json()
        self.assertEqual(alphabet, out['capabilities']['alphabet'])
        cli, rc, _ = self._traverse({'seeds': [{'sequence': self.domain}],
                                     'strategy': {'bounds': {'max_extension_bp': 5}}})
        self.assertEqual(0, rc, cli.get('error'))
        self.assertEqual(alphabet, cli['capabilities']['alphabet'])

    def test_resolve_header_labels(self):
        """resolve reports per-protein support in residues and k-mers and groups identical
        runs; a protein's own sequence is supported end to end only by itself."""
        out, rc, _ = self._traverse({'sequence': self.domain,
                                     'discover': {'max_labels': 10, 'kind': 'header'}},
                                    resolve=True)
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual((K, 'basic', BLOCK - K + 1), (out['k'], out['regime'], out['num_kmers']))
        self.assertEqual(['prot1', 'prot2', 'prot3'], sorted(l['label'] for l in out['labels']))
        for label in out['labels']:
            self.assertEqual('header', label['kind'])
            self.assertEqual(BLOCK - K + 1, label['kmers_supported'])
        self.assertEqual(1, len(out['candidates']))
        self.assertEqual(3, len(out['candidates'][0]['labels']))

        whole = self.records['prot1']
        out, rc, _ = self._traverse({'sequence': whole, 'labels': ['prot1', 'prot2', 'prot3'],
                                     'select': {'policy': 'max_support', 'max_seeds': 1,
                                                'min_block_bp': K}}, resolve=True)
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual(len(whole) - K + 1, out['num_kmers'])
        supported = {l['label']: l['kmers_supported'] for l in out['labels']}
        # prot1 carries every k-mer; prot2 shares N1 + DOMAIN, prot3 DOMAIN + C1
        self.assertEqual({'prot1': 3 * BLOCK - K + 1, 'prot2': 2 * BLOCK - K + 1,
                          'prot3': 2 * BLOCK - K + 1}, supported)
        # max_support merges the overlapping runs: the block all three carry is the domain
        best = out['selection']['seeds'][0]
        self.assertEqual((self.domain, ['prot1', 'prot2', 'prot3']),
                         (best['sequence'], sorted(best['labels'])))

    def test_walk_branches_where_proteins_diverge(self):
        """Seeded on the shared domain: the right arm splits into C1 (prot1, prot3) and C2
        (prot2), the left arm into N1 (prot1, prot2) and N2 (prot3); each leaf reads its
        terminus residue for residue, and the lengths named *_bp count residues."""
        out, rc, _ = self._traverse({
            'seeds': [{'sequence': self.domain, 'labels': ['prot1', 'prot2', 'prot3']}],
            'strategy': {'bounds': {'max_extension_bp': BLOCK + 10}}})
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        self.assertEqual((BLOCK, BLOCK - K + 1),
                         (result['seed']['length_bp'], result['seed']['num_kmers']))
        self.assertEqual({self.c1: {'prot1', 'prot3'}, self.c2: {'prot2'}},
                         self._leaves(result, 'right'))
        self.assertEqual({self.n1: {'prot1', 'prot2'}, self.n2: {'prot3'}},
                         self._leaves(result, 'left'))
        # every label reaches its protein's end: BLOCK residues on each arm
        for entry in result['label_summary']:
            self.assertEqual((BLOCK, BLOCK), (entry['left']['reach_bp'],
                                              entry['right']['reach_bp']), entry)
        # a divergence is not an ambiguous branch
        for side in ('left', 'right'):
            for seg in result['arms'][side]['segments']:
                for ev in seg['events']:
                    if ev['type'] == 'label_end':
                        self.assertNotEqual('branch', ev['reason'])

    def test_a_label_consistent_walk_reads_one_protein(self):
        """Under one protein's header label, or under the column of a family of one, the
        walk reads that protein end to end through the segments it shares."""
        columns = self._column_names()
        family_b = [c for c in columns if 'family_b' in c]
        self.assertEqual(1, len(family_b), columns)
        for labels, name in ((['prot2'], 'prot2'), (['prot3'], 'prot3'), (family_b, 'prot3')):
            out, rc, _ = self._traverse({
                'seeds': [{'sequence': self.domain, 'labels': labels}],
                'strategy': {'bounds': {'max_extension_bp': 2 * BLOCK}}})
            self.assertEqual(0, rc, out.get('error'))
            result = out['results'][0]
            flanks = {}
            for side in ('left', 'right'):
                arm = result['arms'][side]
                self.assertEqual(1, len(arm['paths']), (labels, side))
                flanks[side] = self._flank(arm, arm['paths'][0], side)
            self.assertEqual(self.records[name], flanks['left'] + self.domain + flanks['right'],
                             labels)

    def test_derived_labels(self):
        """Without labels the permitted set is derived from the seed: the proteins (the index
        has a CoordToHeader), resubmittable verbatim with the same result; the column kind
        derives the families."""
        request = {'seeds': [{'sequence': self.domain}],
                   'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': BLOCK},
                                'output': {'timing': False}}}
        out, rc, _ = self._traverse(request)
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        seed = result['seed']
        self.assertTrue(seed['labels_from_seed'])
        self.assertEqual(['prot1', 'prot2', 'prot3'], sorted(seed['labels']))
        self.assertEqual(['header'] * 3, [l['kind'] for l in result['label_dict']])

        named = copy.deepcopy(request)
        named['seeds'][0]['labels'] = seed['labels']
        again, rc, _ = self._traverse(named)
        self.assertEqual(0, rc, again.get('error'))
        mine, theirs = copy.deepcopy(result), copy.deepcopy(again['results'][0])
        for r in (mine, theirs):
            r['seed'].pop('labels_from_seed')
            r.pop('annotation')
            r.pop('work', None)
        self.assertEqual(json.dumps(theirs, sort_keys=True), json.dumps(mine, sort_keys=True))

        columns = self._column_names()
        self.assertEqual(2, len(columns))
        request['strategy']['labels'] = {'seed_label_kind': 'column'}
        out, rc, _ = self._traverse(request)
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        self.assertEqual(columns, sorted(result['seed']['labels']))
        self.assertEqual(['column'] * 2, [l['kind'] for l in result['label_dict']])

    def test_one_strand_and_the_alphabet(self):
        """There is one strand: the reversed domain is not in the graph (nor is any other
        orientation); a character outside the alphabet is refused naming its position; a
        lower-case seed is the upper-case one."""
        out, rc, _ = self._traverse({'seeds': [{'sequence': self.domain[::-1]}],
                                     'strategy': {}})
        self.assertNotEqual(0, rc)
        self.assertIn('Seed is not fully present in the graph', out['error'])
        bad = self.domain[:20] + '*' + self.domain[21:]
        out, rc, _ = self._traverse({'seeds': [{'sequence': bad}], 'strategy': {}})
        self.assertNotEqual(0, rc)
        self.assertIn("Invalid character '*' at position 20 of the seed", out['error'])
        upper, rc, _ = self._traverse({'seeds': [{'sequence': self.domain}],
                                       'strategy': {'bounds': {'max_extension_bp': 20},
                                                    'output': {'timing': False}}})
        self.assertEqual(0, rc, upper.get('error'))
        # edges are the forward (k+1)-mers, and at k = 15 (16 residues of 5 bits each) they no
        # longer fit 64 bits: hashed, as walk_rule states (SPEC §6.6, §6.12); no other strand
        rule = upper['walk_rule']
        self.assertIn('(edges are compared as (k+1)-mers, identified by a 128-bit FNV-1a hash',
                      rule)
        self.assertNotIn('canonical', rule)
        self.assertNotIn('either strand', rule)
        lower, rc, _ = self._traverse({'seeds': [{'sequence': self.domain.lower()}],
                                       'strategy': {'bounds': {'max_extension_bp': 20},
                                                    'output': {'timing': False}}})
        self.assertEqual(0, rc, lower.get('error'))
        self.assertEqual(upper['results'][0]['seed']['validated_seed_id'],
                         lower['results'][0]['seed']['validated_seed_id'])
        self.assertEqual(upper['results'][0]['arms'], lower['results'][0]['arms'])

    def test_record_coordinates(self):
        """Under support trace the seed's occurrences and every run's own residues are
        located in the proteins: positions count residues, and the records confirm them."""
        req = {'seeds': [{'sequence': self.domain}],
               'strategy': {'support': 'trace',
                            'branching': {'on_reconverge': 'keep', 'max_label_branches': 2},
                            'bounds': {'max_extension_bp': BLOCK},
                            'output': {'coordinates': True}}}
        out, rc, _ = self._traverse(req)
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        block = result['coordinates']
        self.assertEqual(('record', K, True), (block['kind'], block['k'], block['complete']))
        names = self._names(result)
        self.assertEqual(3, len(block['seed']))
        checked = 0
        for e in block['seed']:
            record = self.records[names[e['label']]]
            for s, t in e['occurrences']:
                self.assertEqual((BLOCK, BLOCK), (s, t - BLOCK))
                self.assertEqual(self.domain, record[s:t])
                checked += 1
        for side, entries in block['arms'].items():
            arm = result['arms'][side]
            segments = {g['id']: g for g in arm['segments']}
            self.assertEqual(len(arm['runs']), len(entries))
            for entry, run in zip(entries, arm['runs']):
                chain = []
                sid = run['segment']
                while True:
                    chain.append(sid)
                    if not segments[sid]['parents']:
                        break
                    sid = segments[sid]['parents'][0]
                if side == 'right':
                    flank = ''.join(segments[x]['sequence'] for x in reversed(chain))
                    own = flank[run['from_bp']:run['to_bp']]
                else:
                    flank = ''.join(segments[x]['sequence'] for x in chain)
                    own = flank[len(flank) - run['to_bp']:len(flank) - run['from_bp']]
                record = self.records[names[run['label']]]
                for s, t in entry['occurrences']:
                    self.assertEqual(own, record[s:t], (side, s, t))
                    checked += 1
        # the seed's three occurrences (one in each protein) and the six runs' own (the right
        # arm's C1 under prot1 and prot3 and C2 under prot2, the left arm's N1 under prot1 and
        # prot2 and N2 under prot3), one occurrence each
        self.assertEqual(9, checked)

    def test_graphlet_round_trip(self):
        """The graphlet (MGT v1) of a protein walk: its H record names the protein alphabet
        and k, the body parses, dumps back byte for byte, is canonical and reproduces the
        detail: full result."""
        caps = self.session.get(f'{self.url}/traverse/capabilities', timeout=30).json()
        radius = BLOCK + 10
        both = {'direction': 'both', 'bounds': {'max_extension_bp': radius}}
        cases = {
            'tuned (merge), both arms': {'strategy': both},
            'exhaustive (keep)': {'strategy': dict(both, exhaustive=True)},
            'annotate': {'strategy': dict(both, exhaustive=True, labels={'mode': 'annotate'})},
            'trace, column labels': {'strategy': dict(both, support='trace',
                                                      branching={'on_reconverge': 'keep'},
                                                      labels={'seed_label_kind': 'column'})},
            'switch': {'seeds': [{'sequence': self.n2[-30:] + self.domain, 'labels': ['prot3']}],
                       'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': radius},
                                    'branching': {'max_label_branches': 1},
                                    'labels': {'extra': ['prot2'], 'loss_budget': 1,
                                               'change_cost': {'model': 'constant',
                                                               'value': 0.5}}}},
        }
        for case, request in cases.items():
            request = copy.deepcopy(request)
            request.setdefault('seeds', [{'sequence': self.domain}])
            outs = {}
            for detail in ('full', 'graphlet'):
                r = copy.deepcopy(request)
                r['strategy'].setdefault('output', {})['detail'] = detail
                outs[detail], rc, _ = self._traverse(r)
                self.assertEqual(0, rc, (case, outs[detail].get('error')))
            full, out = outs['full']['results'][0], outs['graphlet']['results'][0]
            text = out['graphlet']
            header = text.split('\n')[0].split(' ')
            self.assertEqual(['H', 'mgt', '1', str(K), 'basic', caps['alphabet']], header[:6], case)
            g = graphlet_lib.parse(text)
            self.assertEqual(text, g.dump(), case)
            self.assertTrue(graphlet_lib.is_canonical(text), case)
            gg = graphlet_lib.from_response(out, outs['graphlet'])
            self.assertEqual(normalize_result(full), normalize_result(gg.to_json()), case)
            self.assertEqual([], gg.check_rules(full), case)

    @classmethod
    def _untimed(cls, value):
        """|value| (a parsed body) without its `timing` members, which the server and the CLI
        each measure for themselves."""
        if isinstance(value, dict):
            return {k: cls._untimed(v) for k, v in value.items() if k != 'timing'}
        if isinstance(value, list):
            return [cls._untimed(v) for v in value]
        return value

    def test_api_bodies_are_the_cli_output(self):
        """The server's /traverse answers a protein request with the bytes of `metagraph
        traverse --json` (timing off), and its /resolve with the body of `metagraph traverse
        --resolve --json` but for the timing, which each measures for itself."""
        req = {'seeds': [{'sequence': self.domain},
                         {'sequence': self.n1[-20:] + self.domain, 'labels': ['prot1']}],
               'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': BLOCK},
                            'output': {'timing': False}}}
        ret = self._post('traverse', req)
        self.assertEqual(200, ret.status_code, ret.text)
        _, rc, raw = self._traverse(req, json_flag=True)
        self.assertEqual(0, rc)
        self.assertEqual(ret.content, raw.rstrip(b'\n'))

        req = {'sequence': self.records['prot2'], 'discover': {'max_labels': 10, 'kind': 'header'}}
        ret = self._post('resolve', req)
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        cli, rc, _ = self._traverse(req, resolve=True, json_flag=True)
        self.assertEqual(0, rc, cli.get('error'))
        self.assertIn('timing', out)
        self.assertEqual(json.dumps(self._untimed(cli), sort_keys=True),
                         json.dumps(self._untimed(out), sort_keys=True))
        supported = {l['label']: l['kmers_supported'] for l in out['labels']}
        self.assertEqual({'prot1': 2 * BLOCK - K + 1, 'prot2': 3 * BLOCK - K + 1,
                          'prot3': BLOCK - K + 1}, supported)

    def test_the_record_mapping_is_checked_at_load(self):
        """The load's check of the .seqs (SPEC-multigraph-envelope.md section 10, H4) on a protein
        coordinate index, as UniParc's: the fixture's mapping states its labels (verified), one
        written for the files in the other order is refused at load, and annotate -a spot-checks
        the files against the annotation's columns (the coordinates of their first and last k-mers
        in the graph; a file with none is refused)."""
        caps = self.session.get(f'{self.url}/capabilities', timeout=30).json()
        self.assertEqual('verified', caps['pair_summary']['seqs_order'])
        self.assertEqual({'value': True, 'collisions': 0, 'examples': []}, caps['labels_disjoint'])
        d = os.path.join(self.tempdir.name, 'mapping_check')
        os.makedirs(d, exist_ok=True)
        anno = os.path.join(d, 'annotation.row_diff_brwt_coord.annodbg')
        with open(self.anno, 'rb') as src, open(anno, 'wb') as dst:
            dst.write(src.read())

        def annotate(*args):
            return subprocess.run([METAGRAPH, 'annotate', '-i', self.graph, '--anno-filename',
                                   '--index-header-coords', '-o', os.path.join(d, 'annotation'),
                                   *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        def query():
            return subprocess.run([METAGRAPH, 'query', '--query-mode', 'coords', '-i', self.graph,
                                   '-a', anno, self.files[0]],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        swapped = list(reversed(self.columns))
        res = annotate(*swapped)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        res = query()
        self.assertEqual(1, res.returncode)
        self.assertIn(f"The record mapping '{os.path.join(d, 'annotation.seqs')}' does not belong "
                      f"to the annotation '{anno}': its column 0 is '{swapped[0]}', the "
                      f"annotation's column 0 '{self.columns[0]}'", res.stderr.decode())
        # -a: the swapped files are refused and nothing is written; in column order they pass
        os.remove(os.path.join(d, 'annotation.seqs'))
        res = annotate('-a', anno, *swapped)
        self.assertEqual(1, res.returncode)
        self.assertIn(f"The file '{swapped[0]}' is not column 0 ('{self.columns[0]}') of the "
                      f"annotation '{anno}'", res.stderr.decode())
        self.assertFalse(os.path.exists(os.path.join(d, 'annotation.seqs')))
        res = annotate('-a', anno, *self.columns)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        res = query()
        self.assertEqual(0, res.returncode, res.stderr.decode())
        # a file with k-mers none of which the graph has is refused, nothing written (GPT review
        # 6, finding 12: it passed, nothing checked)
        foreign = os.path.join(d, 'foreign.fa')
        rng = random.Random(20261109)
        with open(foreign, 'w') as f:
            f.write('>f1\n' + ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(300)) + '\n')
        os.remove(os.path.join(d, 'annotation.seqs'))
        res = annotate('-a', anno, foreign, *self.columns[1:])
        self.assertEqual(1, res.returncode, res.stderr.decode())
        self.assertIn(f"The file '{foreign}' has no k-mer of the graph: it cannot be checked "
                      f"against column 0 ('{self.columns[0]}') of the annotation '{anno}'",
                      res.stderr.decode())
        self.assertFalse(os.path.exists(os.path.join(d, 'annotation.seqs')))


if __name__ == '__main__':
    unittest.main()
