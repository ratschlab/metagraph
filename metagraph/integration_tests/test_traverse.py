import json
import os
import random
import shlex
import subprocess
import time
import unittest

import requests
from parameterized import parameterized_class

from base import PROTEIN_MODE, TestingBase, METAGRAPH, start_server
from metagraph.client import GraphClient

GRAPH_MODES = ['basic'] if PROTEIN_MODE else ['basic', 'canonical', 'primary']
K = 31


def random_seq(n, rng, alphabet='ACGT'):
    return ''.join(rng.choice(alphabet) for _ in range(n))


class TraverseFixture(TestingBase):
    """Two records sharing a 100 bp prefix X: A = X + P, B = X + Q (P, Q start differently)."""

    @classmethod
    def _records(cls):
        return f'>A\n{cls.X + cls.P}\n>B\n{cls.X + cls.Q}\n'

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(42)
        alphabet = 'ACDEFGHIKLMNPQRSTVWY' if PROTEIN_MODE else 'ACGT'
        cls.X = random_seq(100, rng, alphabet)
        cls.P = random_seq(80, rng, alphabet)
        cls.Q = random_seq(60, rng, alphabet)
        while cls.Q[0] == cls.P[0]:
            cls.Q = random_seq(60, rng, alphabet)
        fasta = cls.tempdir.name + '/records.fa'
        with open(fasta, 'w') as f:
            f.write(cls._records())
        cls.graph = cls.tempdir.name + '/graph.dbg'
        anno_base = cls.tempdir.name + '/annotation'
        cls.annotation = anno_base + '.column.annodbg'
        cls._build_graph(fasta, cls.graph, K, 'succinct', mode=cls.mode)
        cls._annotate_graph(fasta, cls.graph, anno_base, 'column')


class TraverseCLI:
    def _run(self, request, *flags):
        path = self.tempdir.name + '/request.json'
        with open(path, 'w') as f:
            json.dump(request, f)
        command = f'{METAGRAPH} traverse -i {self.graph} -a {self.annotation} {" ".join(flags)} {path}'
        res = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return res

    def _traverse(self, request, *flags):
        res = self._run(request, *flags)
        self.assertEqual(res.returncode, 0, res.stderr.decode())
        return json.loads(res.stdout.decode())


@parameterized_class(('mode',), input_values=[(mode,) for mode in GRAPH_MODES])
class TestTraverseCLI(TraverseCLI, TraverseFixture):
    def test_fork(self):
        out = self._traverse({'seeds': [{'seed_id': 'x', 'sequence': self.X, 'labels': ['A', 'B']}]})
        self.assertEqual(1, out['schema_version'])
        result = out['results'][0]
        self.assertEqual(['A', 'B'], result['seed']['labels'])
        self.assertFalse(result['seed']['labels_from_seed'])
        self.assertEqual(16, len(result['seed']['validated_seed_id']))
        self.assertTrue(result['seed']['seed_id_mismatch'])

        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(1, len(right['splits']))
        self.assertEqual('divergence', right['splits'][0]['kind'])
        flanks = {}
        for path in right['paths']:
            flank = ''.join(right['segments'][s]['sequence'] for s in path['segments'])
            self.assertEqual(len(flank), path['length_bp'])
            self.assertEqual(1, len(path['end_labels']))
            self.assertEqual('dead_end', path['end_labels'][0]['reason'])
            flanks[result['label_dict'][path['end_labels'][0]['label']]['name']] = flank
        self.assertEqual({'A': self.P, 'B': self.Q}, flanks)

        # the left arm: both records start with the seed
        left = result['arms']['left']
        self.assertEqual(1, len(left['paths']))
        self.assertEqual(0, left['paths'][0]['length_bp'])

    def test_derived_labels_direction_and_bounds(self):
        out = self._traverse({
            'seeds': [{'sequence': self.X}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10}},
        })
        result = out['results'][0]
        self.assertTrue(result['seed']['labels_from_seed'])
        self.assertEqual(['A', 'B'], sorted(result['seed']['labels']))
        self.assertNotIn('left', result['arms'])
        right = result['arms']['right']
        self.assertEqual(10, right['complete_to_bp'])
        for path in right['paths']:
            self.assertEqual(10, path['length_bp'])
            self.assertEqual('max_extension_bp', path['path_reason'])
        # the bases of each path, with the bound applied
        for path in right['paths']:
            flank = ''.join(right['segments'][s]['sequence'] for s in path['segments'])
            self.assertEqual(10, len(flank))

    def test_resolve(self):
        out = self._traverse({'sequence': self.X + self.P, 'discover': {},
                              'select': {'min_block_bp': len(self.X + self.P), 'max_seeds': 1}},
                             '--resolve')
        self.assertEqual(K, out['k'])
        names = {label['label']: label for label in out['labels']}
        self.assertEqual({'A', 'B'}, set(names))
        n = len(self.X + self.P) - K + 1
        self.assertEqual([[0, n]], names['A']['runs'])
        seeds = out['selection']['seeds']
        self.assertEqual(1, len(seeds))
        self.assertEqual(['A'], seeds[0]['labels'])
        self.assertEqual(self.X + self.P, seeds[0]['sequence'])

        # the selected seed is valid /traverse input, with its id
        traversed = self._traverse({'seeds': [{'seed_id': seeds[0]['seed_id'],
                                               'sequence': seeds[0]['sequence'],
                                               'labels': seeds[0]['labels']}]})
        self.assertFalse(traversed['results'][0]['seed']['seed_id_mismatch'])

    def test_invalid_requests(self):
        for request in [
            {},
            {'seeds': []},
            {'seeds': [{'sequence': self.X, 'labels': []}]},
            {'seeds': [{'sequence': self.X}], 'strategy': {'support': 'trace'}},
            {'seeds': [{'sequence': self.X, 'label': ['A']}]},
            {'seeds': [{'sequence': self.X, 'labels': ['C']}]},
            {'seeds': [{'sequence': self.X}], 'strategy': {'branching': {'max_label_branches': -1}}},
            {'seeds': [{'sequence': self.X}] * 11},
        ]:
            res = self._run(request)
            self.assertNotEqual(res.returncode, 0, request)
        self.assertNotEqual(0, self._run({'seeds': [{'sequence': self.X}]},
                                         '--traverse-max-seed-bp', '50').returncode)
        resolve = {'sequence': self.X, 'labels': ['A']}
        self.assertEqual(0, self._run(resolve, '--resolve').returncode)
        self.assertNotEqual(0, self._run(resolve, '--resolve',
                                         '--resolve-max-query-bp', '50').returncode)
        # a limit must be a whole number: 1e6 is refused, not read as 1
        self.assertNotEqual(0, self._run(resolve, '--resolve',
                                         '--resolve-max-query-bp', '1e6').returncode)


@parameterized_class(('mode',), input_values=[(mode,) for mode in GRAPH_MODES])
class TestTraverseBranchLimit(TraverseCLI, TraverseFixture):
    """One label C on both records, so that it reaches the fork after X."""

    @classmethod
    def _records(cls):
        return f'>C\n{cls.X + cls.P}\n>C\n{cls.X + cls.Q}\n'

    def test_branch_limit(self):
        # C is on both records, so it reaches the fork after X: without a branch allowance
        # it ends there, with one it follows both records
        out = self._traverse({'seeds': [{'sequence': self.X, 'labels': ['C']}],
                              'strategy': {'branching': {'max_label_branches': 0}}})
        right = out['results'][0]['arms']['right']
        self.assertEqual(1, len(right['paths']))
        self.assertEqual(['branch'], [e['reason'] for e in right['paths'][0]['end_labels']])
        for branches in [2, 'unlimited']:
            out = self._traverse({'seeds': [{'sequence': self.X, 'labels': ['C']}],
                                  'strategy': {'branching': {'max_label_branches': branches}}})
            right = out['results'][0]['arms']['right']
            flanks = {''.join(right['segments'][s]['sequence'] for s in p['segments'])
                      for p in right['paths']}
            self.assertEqual({self.P, self.Q}, flanks, branches)


class TestTraverseServer(TraverseFixture):
    mode = 'basic'

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.host = '127.0.0.1'
        os.environ['NO_PROXY'] = cls.host
        cls.servers = {}
        port = 4567
        for gated in [True, False]:
            flag = '--experimental-traverse' if gated else ''
            cls.servers[gated] = start_server(cls.graph, cls.annotation, cls.host, port, flag)
            port = cls.servers[gated][1] + 1
        # a server answers 503 while it loads: wait until each one answers
        deadline = time.time() + 60
        for gated, (_, port) in cls.servers.items():
            route = 'traverse/capabilities' if gated else 'stats'
            while True:
                try:
                    if requests.get(f'http://{cls.host}:{port}/{route}', timeout=1).status_code == 200:
                        break
                except requests.RequestException:
                    pass
                if time.time() > deadline:
                    raise RuntimeError(f'The server on port {port} did not get ready')
                time.sleep(0.05)

    @classmethod
    def tearDownClass(cls):
        for process, _ in cls.servers.values():
            process.kill()

    def _url(self, route, gated=True):
        return f'http://{self.host}:{self.servers[gated][1]}/{route}'

    def test_capabilities(self):
        res = requests.get(self._url('traverse/capabilities'))
        self.assertEqual(200, res.status_code, res.text)
        caps = res.json()
        self.assertEqual(1, caps['schema_version'])
        self.assertEqual(K, caps['k'])
        self.assertEqual('basic', caps['regime'])
        self.assertEqual(10, caps['limits']['max_seeds'])

    def test_traverse_and_resolve(self):
        res = requests.post(self._url('traverse'),
                            data=json.dumps({'seeds': [{'sequence': self.X, 'labels': ['A']}]}))
        self.assertEqual(200, res.status_code, res.text)
        paths = res.json()['results'][0]['arms']['right']['paths']
        self.assertEqual(1, len(paths))
        self.assertEqual(len(self.P), paths[0]['length_bp'])

        res = requests.post(self._url('resolve'),
                            data=json.dumps({'sequence': self.X, 'labels': ['A', 'B']}))
        self.assertEqual(200, res.status_code, res.text)
        self.assertEqual(2, len(res.json()['labels']))

    def test_invalid_request(self):
        # an unknown field is refused even in an otherwise valid request
        request = {'seeds': [{'sequence': self.X}]}
        self.assertEqual(200, requests.post(self._url('traverse'), data=json.dumps(request)).status_code)
        request['graph'] = 'x'
        self.assertEqual(400, requests.post(self._url('traverse'), data=json.dumps(request)).status_code)

    def test_python_client(self):
        client = GraphClient(self.host, self.servers[True][1])
        profile = client.resolve(self.X + self.P, discover={},
                                 select={'min_block_bp': len(self.X + self.P), 'max_seeds': 1})
        seed = profile['selection']['seeds'][0]
        self.assertEqual(['A'], seed['labels'])
        result = client.traverse([{'seed_id': seed['seed_id'], 'sequence': seed['sequence'],
                                   'labels': seed['labels']}],
                                 strategy={'direction': 'left'})['results'][0]
        self.assertFalse(result['seed']['seed_id_mismatch'])
        self.assertEqual(['left'], list(result['arms']))
        with self.assertRaises(RuntimeError):
            client.traverse([{'sequence': self.X}], strategy={'unknown': 1})

    def test_routes_are_gated(self):
        for route in ['traverse/capabilities']:
            self.assertEqual(404, requests.get(self._url(route, gated=False)).status_code)
        for route in ['traverse', 'resolve']:
            res = requests.post(self._url(route, gated=False), data=json.dumps({}))
            self.assertEqual(404, res.status_code)


if __name__ == '__main__':
    unittest.main()
