import json
import os
import random
import shlex
import socket
import subprocess
import time
import unittest

import requests
from parameterized import parameterized_class

from base import PROTEIN_MODE, TestingBase, METAGRAPH, start_server

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


class TraverseServer(TraverseFixture):
    """A server on the fixture's index; the tests post their requests to it."""
    server_flags = ''

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.host = '127.0.0.1'
        cls._prepare()
        os.environ['NO_PROXY'] = cls.host
        # the test classes run in parallel: each server gets a port that is free now
        with socket.socket() as s:
            s.bind((cls.host, 0))
            port = s.getsockname()[1]
        cls.server, cls.port = cls._start_server(port)
        # the server answers 503 while it loads: wait until it answers
        deadline = time.time() + 60
        while True:
            try:
                if requests.get(cls._url('stats'), timeout=1).status_code == 200:
                    break
            except requests.RequestException:
                pass
            if cls.server.poll() is not None or time.time() > deadline:
                cls.server.kill()
                raise RuntimeError(f'The server on port {cls.port} did not get ready')
            time.sleep(0.05)

    @classmethod
    def _prepare(cls):
        """Whatever the server needs, built before its port is chosen."""

    @classmethod
    def _start_server(cls, port):
        return start_server(cls.graph, cls.annotation, cls.host, port, cls.server_flags)

    @classmethod
    def tearDownClass(cls):
        cls.server.kill()
        super().tearDownClass()

    @classmethod
    def _url(cls, route):
        return f'http://{cls.host}:{cls.port}/{route}'

    def _post(self, route, request):
        return requests.post(self._url(route), data=json.dumps(request))

    def _traverse(self, request):
        res = self._post('traverse', request)
        self.assertEqual(200, res.status_code, res.text)
        return res.json()


@parameterized_class(('mode',), input_values=[(mode,) for mode in GRAPH_MODES])
class TestTraverse(TraverseServer):
    def test_fork(self):
        out = self._traverse({'seeds': [{'sequence': self.X, 'labels': ['A', 'B']}]})
        result = out['results'][0]
        self.assertEqual(['A', 'B'], result['seed']['labels'])
        self.assertFalse(result['seed']['labels_from_seed'])

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

    def test_search_labels_into_traverse(self):
        query = self.X + self.P
        res = self._post('search', {'FASTA': f'>q\n{query}', 'with_signature': True,
                                    'discovery_fraction': 0})
        self.assertEqual(200, res.status_code, res.text)
        results = res.json()[0]['results']
        self.assertEqual({'A', 'B'}, {r['sample'] for r in results})
        n = len(query) - K + 1

        # the query with the labels on all its k-mers is valid /traverse input
        whole = [r['sample'] for r in results if r['signature']['presence_mask'] == f'x{n}']
        self.assertEqual(['A'], whole)
        self.assertEqual([n], [r['kmer_count'] for r in results if r['sample'] == 'A'])
        traversed = self._traverse({'seeds': [{'sequence': query, 'labels': whole}]})
        self.assertEqual(['A'], traversed['results'][0]['seed']['labels'])

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
            {'seeds': [{'sequence': self.X}], 'graph': 'x'},
            {'seeds': [{'sequence': self.X}], 'graphs': ['x']},  # no list of graphs here
        ]:
            self.assertEqual(400, self._post('traverse', request).status_code, request)
        # there is no /resolve route
        res = self._post('resolve', {'sequence': self.X, 'labels': ['A']})
        self.assertEqual(404, res.status_code)


@parameterized_class(('mode',), input_values=[(mode,) for mode in GRAPH_MODES])
class TestTraverseBranchLimit(TraverseServer):
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


class TestTraverseGraphList(TraverseServer):
    """A server on a list of graphs: record A is in the graph `one`, record B in `two`."""
    mode = 'basic'

    @classmethod
    def _prepare(cls):
        rows = []
        for name, record, flank in [('one', 'A', cls.P), ('two', 'B', cls.Q)]:
            base = f'{cls.tempdir.name}/{name}'
            with open(base + '.fa', 'w') as f:
                f.write(f'>{record}\n{cls.X + flank}\n')
            cls._build_graph(base + '.fa', base + '.dbg', K, 'succinct', mode=cls.mode)
            cls._annotate_graph(base + '.fa', base + '.dbg', base, 'column')
            rows.append(f'{name},{base}.dbg,{base}.column.annodbg\n')
        with open(cls.tempdir.name + '/graphs.csv', 'w') as f:
            f.writelines(rows)

    @classmethod
    def _start_server(cls, port):
        # as start_server: a server that exits at once lost its port to another test
        for _ in range(100):
            command = [METAGRAPH, 'server_query', cls.tempdir.name + '/graphs.csv',
                       '--port', str(port), '--address', cls.host, '-p', '2']
            process = subprocess.Popen(command, stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL)
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                return process, port
            port += 1
        raise RuntimeError("Couldn't start server")

    def _flanks(self, walk):
        right = walk['arms']['right']
        return [''.join(right['segments'][s]['sequence'] for s in path['segments'])
                for path in right['paths']]

    def test_every_graph_is_walked(self):
        out = self._traverse({'seeds': [{'sequence': self.X, 'labels': ['A', 'B']},
                                        {'sequence': self.X}]})
        # each graph follows the labels it has
        for result in out['results']:
            walks = result['walks']
            self.assertEqual(['one', 'two'], [walk['graph'] for walk in walks])
            self.assertEqual([['A'], ['B']], [walk['seed']['labels'] for walk in walks])
            self.assertEqual([[self.P], [self.Q]], [self._flanks(walk) for walk in walks])

    def test_selected_graphs(self):
        seeds = [{'sequence': self.X}]
        walks = self._traverse({'seeds': seeds, 'graphs': ['two']})['results'][0]['walks']
        self.assertEqual(['two'], [walk['graph'] for walk in walks])
        for graphs in [['three'], 'two']:
            res = self._post('traverse', {'seeds': seeds, 'graphs': graphs})
            self.assertEqual(400, res.status_code, graphs)
        # a seed shorter than k in every graph
        self.assertEqual(400, self._post('traverse', {'seeds': [{'sequence': self.X[:K - 1]}]}).status_code)

    def test_graph_without_the_seed_or_its_labels(self):
        out = self._traverse({'seeds': [
            {'sequence': self.X, 'labels': ['A']},  # no graph but `one` has the label
            {'sequence': self.P[:40]},              # no graph but `one` has the seed
            {'sequence': self.Q[:20] + self.P[:20]},
        ]})
        self.assertEqual([['one'], ['one'], []],
                         [[walk['graph'] for walk in result['walks']] for result in out['results']])

    def test_label_in_no_selected_graph(self):
        res = self._post('traverse', {'seeds': [{'sequence': self.X, 'labels': ['A', 'C']}]})
        self.assertEqual(400, res.status_code)
        self.assertIn("'C'", res.text)


class TestTraverseLimits(TraverseServer):
    """The server's limits of a request."""
    mode = 'basic'
    server_flags = '--traverse-max-seeds 2'

    def test_limits(self):
        short = self.X[:50]
        self.assertEqual(200, self._post('traverse', {'seeds': [{'sequence': short}] * 2}).status_code)
        self.assertEqual(400, self._post('traverse', {'seeds': [{'sequence': short}] * 3}).status_code)
        # a seed over the fixed 10,000 bp limit
        self.assertEqual(400, self._post('traverse', {'seeds': [{'sequence': 'A' * 10001}]}).status_code)

    def test_malformed_limit(self):
        # a limit must be a whole number: 1e6 is refused, not read as 1
        command = (f'{METAGRAPH} server_query -i {self.graph} -a {self.annotation}'
                   f' --port {self.port + 1} --traverse-max-seeds 1e6')
        res = subprocess.run(shlex.split(command), stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, timeout=60)
        self.assertNotEqual(0, res.returncode)


if __name__ == '__main__':
    unittest.main()
