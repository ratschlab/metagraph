"""The /search answer budget (lane T5b-Y4, owner decision 2026-11-13; SPEC-multigraph-envelope.md
section 3): --search-max-answer-memory-mb (4,096 MiB by default, 0 = none) bounds what one POST
/search holds, each pair's rows (their per-position values) and the request's JSON answer and its
text, priced before they are built. Past it the whole request is refused 400
search_answer_too_large ({error, code, max_memory_mb}, error first), on a multi-graph server too,
the same request again and never cached; below it every answer is the one without a budget, byte
for byte. GET /capabilities states it as `search` (a single-graph server's top level, the
multigraph block on a multi-graph server).

Each test learns a request's charge from the server's log ("answer budget: values at most N B a
pair, answer M B"; the accounts only grow, so the bound a request needs is the larger of the two,
whatever the order its work charged in) and checks it at the bound of that many MiB rounded up
(answered, byte for byte as without a budget) and one MiB below (refused).
"""

import json
import math
import os
import random
import re
import shlex
import subprocess
import unittest
from tempfile import TemporaryDirectory

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH, graph_file_extension
from test_multigraph import Server

BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')
MIB = 1 << 20
ROW_DIFF_BRWT_COORD = '.row_diff_brwt_coord.annodbg'
MODES = {
    'plain': {},
    'abundance_sum': {'abundance_sum': True},
    'with_signature': {'with_signature': True},
    'query_counts': {'query_counts': True},
    'query_coords': {'query_coords': True},
}


def refusal_text(max_mb):
    return ('{"error":"search: the answer would hold more than --search-max-answer-memory-mb '
            f'{max_mb} MiB: nothing partial is sent; ask for fewer top_labels, a larger '
            'discovery_fraction, fewer or shorter FASTA records, or no per-position values '
            '(with_signature, query_counts, query_coords)","code":"search_answer_too_large",'
            f'"max_memory_mb":{max_mb}}}')


BLOCK = {'max_memory_mb': 4096,
         'rule': 'SPEC-multigraph-envelope.md section 3, the /search answer budget',
         'too_large': {'status': 400, 'code': 'search_answer_too_large',
                       'members': ['code', 'error', 'max_memory_mb']}}

CHARGE = re.compile(r'Request (\d+): answer budget: values at most (\d+) B a pair, answer (\d+) B, '
                    r'bound (\d+) MiB')


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestSearchAnswerBudget(TestingBase):
    k = 31
    record_length = 60
    records_per_pair = 1000
    started = 0

    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        d = cls.tempdir.name
        rng = random.Random(57)
        # two pairs of records of r237's form (a k=31 basic succinct graph, row_diff_brwt_coord
        # with its .seqs), a record's first base differing from the base after the previous
        # record's first k-mer, so that a query of first k-mers matches each record once
        cls.records = {}
        cls.pairs = {}
        last = None
        for name in ('pa', 'pb'):
            records = []
            for _ in range(cls.records_per_pair):
                record = ''.join(rng.choice('ACGT') for _ in range(cls.record_length))
                if last and record[0] == last[cls.k]:
                    record = ('A' if record[0] != 'A' else 'C') + record[1:]
                records.append(record)
                last = record
            cls.records[name] = records
            fasta = f'{d}/{name}.fa'
            with open(fasta, 'w') as f:
                for i, record in enumerate(records):
                    f.write(f'>{name}_rec{i}\n{record}\n')
            graph = f'{d}/{name}' + graph_file_extension['succinct']
            anno = f'{d}/{name}_anno'
            cls._build_graph(fasta, f'{d}/{name}', cls.k, 'succinct', 'basic')
            cls._annotate_graph(fasta, graph, anno, 'row_diff_brwt_coord', anno_type='filename')
            res = subprocess.run(shlex.split(METAGRAPH) + [
                'annotate', '--anno-filename', '--index-header-coords', '-i', graph, '-o', anno,
                fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert res.returncode == 0, res.stderr.decode()
            cls.pairs[name] = (graph, anno + ROW_DIFF_BRWT_COORD)
        cls.csv = f'{d}/graphs.csv'
        with open(cls.csv, 'w') as f:
            for name, (graph, anno) in cls.pairs.items():
                f.write(f'{name},{graph},{anno}\n')
        cls.servers = []

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers:
            server.stop()
        cls.tempdir.cleanup()

    def server(self, extra=(), multi=False, binary=None):
        # a log of its own (a stopped server leaves the list: its length repeats)
        TestSearchAnswerBudget.started += 1
        name = f'server_{TestSearchAnswerBudget.started}.log'
        args = [self.csv] if multi else ['-i', self.pairs['pa'][0], '-a', self.pairs['pa'][1]]
        server = Server(args + list(extra), os.path.join(self.tempdir.name, name), threads=2,
                        binary=binary)
        self.servers.append(server)
        return server

    def stop(self, server):
        server.stop()
        self.servers.remove(server)

    # -------------------------------------------------------------------- requests

    def first_kmers(self, name, n=None):
        return ''.join(r[:self.k] for r in self.records[name][:n])

    def body(self, mode, records, top_labels=10000, **more):
        fasta = ''.join(f'>q{i}\n{s}\n' for i, s in enumerate(records))
        return dict({'FASTA': fasta, 'top_labels': top_labels, 'discovery_fraction': 0.0},
                    **MODES[mode], **more)

    def charge(self, server, before):
        """the request's charge, the larger account, from the log line of the request after
        |before| lines"""
        lines = [m for m in CHARGE.finditer(server.text()) if m.start() >= before]
        self.assertEqual(1, len(lines), server.text()[-3000:])
        values, answer = int(lines[0].group(2)), int(lines[0].group(3))
        return values, answer

    def answered(self, server, body):
        before = len(server.text())
        r = server.post('search', body)
        self.assertEqual(200, r.status_code, r.text[:2000])
        return r.content, self.charge(server, before)

    def refused(self, server, body, max_mb):
        r = server.post('search', body)
        self.assertEqual(400, r.status_code, r.text[:2000])
        self.assertEqual(refusal_text(max_mb), r.text)
        self.assertEqual('application/json', r.headers['Content-Type'])
        self.assertNotIn('Content-Encoding', r.headers)
        self.assertEqual(['error', 'code', 'max_memory_mb'], list(json.loads(r.text)))
        return r

    # ------------------------------------------------------------------------ tests

    def test_capabilities(self):
        """GET /capabilities states the budget: at the top level of a single-graph server, in
        the multigraph block of a multi-graph server, null when 0"""
        single = self.server()
        self.assertEqual(BLOCK, single.get('capabilities').json()['search'])
        self.stop(single)
        off = self.server(['--search-max-answer-memory-mb', '0'])
        self.assertEqual(dict(BLOCK, max_memory_mb=None), off.get('capabilities').json()['search'])
        self.stop(off)
        multi = self.server(['--search-max-answer-memory-mb', '7'], multi=True)
        caps = multi.get('capabilities').json()
        self.assertNotIn('search', caps)
        self.assertEqual(dict(BLOCK, max_memory_mb=7), caps['multigraph']['search'])
        self.stop(multi)

    def test_flag_range(self):
        """the flag takes an integer in [0, 1048576]"""
        for value in ('1048577', '-1', '1.5', 'x'):
            res = subprocess.run(shlex.split(METAGRAPH) + [
                'server_query', '-i', self.pairs['pa'][0], '-a', self.pairs['pa'][1],
                '--search-max-answer-memory-mb', value, '--port', '1', '--address', '127.0.0.1'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
            self.assertNotEqual(0, res.returncode, value)
            self.assertIn('--search-max-answer-memory-mb', res.stderr.decode() + res.stdout.decode())

    def test_each_mode_just_under_and_over(self):
        """every mode: answered at the bound its charge needs, byte for byte as without a
        budget (and as the base binary, where one is given); refused one MiB below"""
        # 30 query records of every record's first k-mer: 30 x 1,000 rows
        body_records = [self.first_kmers('pa')] * 30
        unbounded = self.server(['--search-max-answer-memory-mb', '0'])
        default = self.server()
        needs = {}
        plain_answers = {}
        for mode in MODES:
            body = self.body(mode, body_records)
            r = unbounded.post('search', body)
            self.assertEqual(200, r.status_code, r.text[:2000])
            plain_answers[mode] = r.content
            got, (values, answer) = self.answered(default, body)
            self.assertEqual(r.content, got, mode)
            needs[mode] = max(values, answer)
            # the fixture passes one MiB: the bound below it refuses
            self.assertGreater(needs[mode], MIB, mode)
            print(f'[ {mode} ] values {values} B, answer {answer} B, text {len(got)} B',
                  flush=True)
        if BASE_BINARY:
            base = self.server(binary=BASE_BINARY)
            for mode in MODES:
                r = base.post('search', self.body(mode, body_records))
                self.assertEqual(plain_answers[mode], r.content, mode)
            self.stop(base)
        self.stop(unbounded)
        self.stop(default)
        for mode in MODES:
            bound = math.ceil(needs[mode] / MIB)
            under = self.server(['--search-max-answer-memory-mb', str(bound)])
            got, _ = self.answered(under, self.body(mode, body_records))
            self.assertEqual(plain_answers[mode], got, mode)
            self.stop(under)
            over = self.server(['--search-max-answer-memory-mb', str(bound - 1)])
            self.refused(over, self.body(mode, body_records), bound - 1)
            # the same request again, the same refusal; the server answers the next request
            self.refused(over, self.body(mode, body_records), bound - 1)
            small = self.body(mode, [self.records['pa'][0]])
            self.assertEqual(200, over.post('search', small).status_code)
            self.assertIn('search_answer_too_large', over.text())
            self.stop(over)

    def test_many_records(self):
        """a FASTA of thousands of records (top_labels 10,000), as the service forwards one:
        answered under its bound, refused over it"""
        body = self.body('plain', self.records['pa'] * 3, top_labels=10000)
        default = self.server()
        expected, (values, answer) = self.answered(default, body)
        self.stop(default)
        need = max(values, answer)
        self.assertGreater(need, MIB)
        bound = math.ceil(need / MIB)
        under = self.server(['--search-max-answer-memory-mb', str(bound)])
        self.assertEqual(expected, self.answered(under, body)[0])
        self.stop(under)
        over = self.server(['--search-max-answer-memory-mb', str(bound - 1)])
        self.refused(over, body, bound - 1)
        self.stop(over)

    def test_multi_graph(self):
        """a multi-graph server refuses the whole request when the pairs' answers together
        pass the bound, each pair alone within it; and when one pair's values pass it"""
        body_records = [self.first_kmers('pa'), self.first_kmers('pb')] * 10
        default = self.server(multi=True)
        single_needs = []
        for graphs in (['pa'], ['pb']):
            _, (values, answer) = self.answered(default, self.body('query_coords', body_records,
                                                                   graphs=graphs))
            single_needs.append(max(values, answer))
        both, (values_both, answer_both) = self.answered(
            default, self.body('query_coords', body_records, graphs=['pa', 'pb']))
        self.stop(default)
        bound = math.ceil(max(single_needs) / MIB)
        # the sum passes the bound each pair is within
        self.assertGreater(answer_both, bound * MIB)
        self.assertLessEqual(values_both, bound * MIB)
        server = self.server(['--search-max-answer-memory-mb', str(bound)], multi=True)
        for graphs in (['pa'], ['pb']):
            r = server.post('search', self.body('query_coords', body_records, graphs=graphs))
            self.assertEqual(200, r.status_code)
        self.refused(server, self.body('query_coords', body_records, graphs=['pa', 'pb']), bound)
        self.assertIn('the answer account passed', server.text())
        self.stop(server)
        # one pair's values past the bound: a long query of many positions on one pair
        long_query = [self.records['pa'][0] * 400]
        default = self.server(multi=True)
        _, (values, answer) = self.answered(
            default, self.body('query_coords', long_query, graphs=['pa', 'pb']))
        self.stop(default)
        if values > MIB and values > answer:
            bound = math.ceil(values / MIB) - 1
            server = self.server(['--search-max-answer-memory-mb', str(bound)], multi=True)
            self.refused(server, self.body('query_coords', long_query, graphs=['pa', 'pb']),
                         bound)
            self.assertRegex(server.text(), r'the values account of pair \(pa, ')
            self.stop(server)
        # the answer with both pairs at a bound within it is the unbounded answer
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        r = unbounded.post('search', self.body('query_coords', body_records, graphs=['pa', 'pb']))
        self.assertEqual(both, r.content)
        self.stop(unbounded)

    def test_in_ram_and_cache(self):
        """in_ram changes nothing (the same charge, the same refusal); a refused request is not
        cached and is refused again, an answered one is cached"""
        body = self.body('with_signature', [self.first_kmers('pa')] * 30, graphs=['pa'])
        # without the cache (on by default, 1 GB; its key does not tell in_ram apart: the second
        # request would be a hit)
        server = self.server(['--mmap', '--mem-cap-gb', '1', '--cache-size', '0'], multi=True)
        plain, (values, answer) = self.answered(server, body)
        in_ram, charge = self.answered(server, dict(body, in_ram=True))
        self.assertEqual(plain, in_ram)
        self.assertEqual((values, answer), charge)
        self.stop(server)
        # with the cache (--cache-size in GB) the answer account holds its copy of the text too:
        # the bound below the charge without it refuses
        bound = math.ceil(max(values, answer) / MIB) - 1
        self.assertGreaterEqual(bound, 1)
        server = self.server(['--mmap', '--mem-cap-gb', '1', '--cache-size', '1',
                              '--search-max-answer-memory-mb', str(bound)], multi=True)
        for request in (body, dict(body, in_ram=True), body):
            self.refused(server, request, bound)
        self.assertNotIn('cache HIT', server.text())
        self.assertNotIn('cache STORE', server.text())
        small = self.body('with_signature', [self.records['pa'][0]], graphs=['pa'])
        first = server.post('search', small)
        second = server.post('search', small)
        self.assertEqual((200, first.content), (second.status_code, second.content))
        self.assertIn('cache HIT', server.text())
        self.stop(server)


if __name__ == '__main__':
    unittest.main()
