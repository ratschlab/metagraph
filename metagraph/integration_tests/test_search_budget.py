"""The /search answer budget (lane T5b-Y4, owner decision 2026-11-13; per (record, pair) since lane
T7-P1; SPEC-multigraph-envelope.md section 3): --search-max-answer-memory-mb (4,096 MiB by default,
0 = none) bounds what one POST /search holds, each charge priced before what it prices is built:
the values of each unit (one FASTA record on one pair: its per-position values as they grow, its
rows as they are handed over), each unit's own JSON object (checked alone against the bound), and
the request's answer (every object, its text and the transport's and the cache's copies).

A unit past the bound refuses the whole request 400 search_answer_too_large ({error, account,
code, max_memory_mb}, error first; `account` "values" when some unit's values pass, else "answer")
-- or, with the request's `partial_refusal`, it is answered alone as a refused cut entry (no row;
returned 0, complete false, nth_score null, tied and tied_returned 0, nth_median null where the
cut has abundances, refused {account, code, max_memory_mb[, priced_mb]}) while every other unit
is answered byte for byte as without a budget. The request's answer account passing refuses the
whole request 400 "answer" either way. Below the bound every answer is the one without a budget,
byte for byte; a refusal is never cached. GET /capabilities states it as `search` (a single-graph
server's top level, the multigraph block on a multi-graph server) with the concurrency its model
total depends on (-p, --threads-each, --http-threads: every server here sets all three).

Each test learns a request's charge from the server's log ("answer budget: values at most N B a
unit, answer M B, bound B MiB, R of U units refused"; the accounts only grow, so the bound a
request needs is the larger of the two) or fixes a bound its fixture passes by far. Each refusal's
log line states the account whose charge passed the bound (a unit's with its record and pair), what
it held before that charge, and the charge; each unit refused alone has a line of its own, in the
answer's order (by record, then the pair's place among the record's cuts). `partial_refusal` is a
JSON boolean; any other value is refused 400.
"""

import ctypes
import json
import math
import os
import random
import re
import shlex
import struct
import subprocess
import sys
import threading
import time
import unittest
from tempfile import TemporaryDirectory

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
CODE = 'search_answer_too_large'


def refusal_text(max_mb, account):
    return ('{"error":"search: the answer would hold more than --search-max-answer-memory-mb '
            f'{max_mb} MiB: nothing partial is sent; ask for fewer top_labels, a larger '
            'discovery_fraction, fewer or shorter FASTA records, or no per-position values '
            f'(with_signature, query_counts, query_coords)","account":"{account}",'
            f'"code":"{CODE}","max_memory_mb":{max_mb}}}')


def block(max_mb=4096, pairs=2, records=2, handlers=2):
    return {'handler_threads': handlers,
            'max_memory_mb': max_mb,
            'pair_searches_at_once': pairs,
            'records_at_once': records,
            'rule': 'SPEC-multigraph-envelope.md section 3, the /search answer budget',
            'too_large': {'status': 400, 'code': CODE, 'accounts': ['answer', 'values'],
                          'members': ['account', 'code', 'error', 'max_memory_mb']},
            'unit_refused': {'field': 'partial_refusal', 'code': CODE,
                             'accounts': ['answer', 'values'],
                             'members': ['account', 'code', 'max_memory_mb', 'priced_mb']}}


CHARGE = re.compile(r'Request (\d+): answer budget: values at most (\d+) B a unit, answer (\d+) B, '
                    r'bound (\d+) MiB, (\d+) of (\d+) units refused\n')
REFUSAL = re.compile(r'Request (\d+): search_answer_too_large: (the (?:values|answer) account(?: of .*?)?) '
                     r'passed the bound of (\d+) MiB: held (\d+) B \((\d+\.\d) MiB\) before a charge of '
                     r'(\d+) B \((\d+\.\d) MiB\)\n')
UNIT = re.compile(r'Request (\d+): unit refused: record (\d+) on (.+?): the (values|answer) account '
                  r'passed the bound of (\d+) MiB: held (\d+) B \((\d+\.\d) MiB\) before a charge of '
                  r'(\d+) B \((\d+\.\d) MiB\)\n')


def compact(value):
    """|value|'s JSON without a space (its members in byte order)"""
    return json.dumps(value, separators=(',', ':'), sort_keys=True)


def styled(value):
    """the bytes the server writes for a /search answer |value| (jsoncpp's StreamWriterBuilder at
    its defaults: a tab a level, members in byte order, each member and each element of an array
    on a line of its own, a container that is a member's value on the line after its name; no
    float in these answers)"""
    out = []

    def write(v, level, own_line):
        if isinstance(v, (dict, list)) and v:
            out.append(('\n' + '\t' * level if own_line else '') + ('{' if isinstance(v, dict) else '['))
            items = sorted(v.items()) if isinstance(v, dict) else list(enumerate(v))
            for n, (key, child) in enumerate(items):
                out.append(',' if n else '')
                out.append('\n' + '\t' * (level + 1))
                if isinstance(v, dict):
                    out.append(json.dumps(key) + ' : ')
                write(child, level + 1, isinstance(v, dict))
            out.append('\n' + '\t' * level + ('}' if isinstance(v, dict) else ']'))
        elif isinstance(v, (dict, list)):
            out.append('{}' if isinstance(v, dict) else '[]')
        else:
            assert not isinstance(v, float), v
            out.append(json.dumps(v))

    write(value, 0, False)
    return ''.join(out)


def refused_entry(entry, account, max_mb, priced_mb=None):
    """the cut entry of a refused unit whose entry without a budget is |entry|: its tags and
    member set, no row"""
    out = dict(entry, returned=0, complete=False, nth_score=None, tied=0, tied_returned=0)
    if 'nth_median' in entry:
        out.update(nth_median=None, tied_same_median=0)
    out['refused'] = {'account': account, 'code': CODE, 'max_memory_mb': max_mb}
    if priced_mb is not None:
        out['refused']['priced_mb'] = priced_mb
    return out


def with_refused(answer, refused, max_mb):
    """|answer| (an answer without a budget, parsed) with the units of |refused| ({(record index,
    graph name or None on a single-graph server): (account, priced_mb or None)}) refused: their
    rows dropped, their cut entries refused, everything else as it is"""
    out = []
    for i, obj in enumerate(answer):
        obj = dict(obj, results=[row for row in obj['results']
                                 if (i, row.get('graph')) not in refused])
        cuts = []
        for cut in obj['cuts']:
            if (i, cut['graph']) in refused:
                account, priced_mb = refused[(i, cut['graph'])]
                cut = refused_entry(cut, account, max_mb, priced_mb)
            cuts.append(cut)
        obj['cuts'] = cuts
        out.append(obj)
    return out


def memory_readers():
    """(footprint, peak) of a process by pid, in bytes: on Linux the resident set and its high
    water mark (VmRSS, VmHWM), on macOS the physical footprint and its lifetime maximum
    (proc_pid_rusage, rusage_info_v4: 16 bytes of uuid, then ri_phys_footprint the 8th uint64
    and ri_lifetime_max_phys_footprint the 29th); None elsewhere"""
    if sys.platform.startswith('linux'):
        def status(pid, key):
            with open(f'/proc/{pid}/status') as f:
                for line in f:
                    if line.startswith(key + ':'):
                        return int(line.split()[1]) << 10
            raise AssertionError(key)
        return (lambda pid: status(pid, 'VmRSS')), (lambda pid: status(pid, 'VmHWM'))
    if sys.platform == 'darwin':
        libproc = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)

        def rusage(pid):
            buf = (ctypes.c_uint64 * 64)()
            assert libproc.proc_pid_rusage(pid, 4, ctypes.byref(buf)) == 0
            return buf
        return (lambda pid: rusage(pid)[2 + 7]), (lambda pid: rusage(pid)[2 + 28])
    return None


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
        # three pairs of records of r237's form (a k=31 basic succinct graph, row_diff_brwt_coord
        # with its .seqs), a record's first base differing from the base after the previous
        # record's first k-mer, so that a query of first k-mers matches each record once
        cls.records = {}
        cls.pairs = {}
        cls.fasta = {}
        last = None
        for name in ('pa', 'pb', 'pc'):
            records = []
            for _ in range(cls.records_per_pair):
                record = ''.join(rng.choice('ACGT') for _ in range(cls.record_length))
                if last and record[0] == last[cls.k]:
                    record = ('A' if record[0] != 'A' else 'C') + record[1:]
                records.append(record)
                last = record
            cls.records[name] = records
            cls.fasta[name] = f'{d}/{name}.fa'
            with open(cls.fasta[name], 'w') as f:
                for i, record in enumerate(records):
                    f.write(f'>{name}_rec{i}\n{record}\n')
            cls.pairs[name] = cls._coord_pair(d, name)
        # telomere-like (Y3's polya25k shape): 50 records each holding (CCCTAA)n, about 330
        # occurrences of each of its six k-mers, between flanks without a C, and one record of
        # unique sequence
        flank = lambda: ''.join(rng.choice('AGT') for _ in range(40))
        cls.unique = ''.join(rng.choice('ACGT') for _ in range(400))
        cls.fasta['pt'] = f'{d}/pt.fa'
        with open(cls.fasta['pt'], 'w') as f:
            for i in range(50):
                f.write(f'>t{i}\n{flank()}{"CCCTAA" * 333}{flank()}\n')
            f.write(f'>u0\n{cls.unique}\n')
        cls.pairs['pt'] = cls._coord_pair(d, 'pt')
        cls.repeat = 'CCCTAA' * 50
        # pa's graph with a column annotation of its headers, no coordinates and no .seqs: a
        # coordinate query is refused on it, and the other modes take the batched path
        cls._annotate_graph(cls.fasta['pa'], cls.pairs['pa'][0], f'{d}/pn_anno', 'column',
                            anno_type='header')
        cls.pairs['pn'] = (cls.pairs['pa'][0], f'{d}/pn_anno.column.annodbg')
        # pa also under the name pa2; the telomere pair under eleven names (as r237's eleven
        # pairs) beside pa under the name u
        cls.csv = f'{d}/graphs.csv'
        with open(cls.csv, 'w') as f:
            for name, pair in (('pa', 'pa'), ('pa2', 'pa'), ('pb', 'pb'), ('pc', 'pc'),
                               ('pn', 'pn')):
                f.write(f'{name},{cls.pairs[pair][0]},{cls.pairs[pair][1]}\n')
        cls.telomere_names = [f't{i:02d}' for i in range(11)]
        cls.telomere_csv = f'{d}/telomere.csv'
        with open(cls.telomere_csv, 'w') as f:
            for name in cls.telomere_names:
                f.write(f'{name},{cls.pairs["pt"][0]},{cls.pairs["pt"][1]}\n')
            f.write(f'u,{cls.pairs["pa"][0]},{cls.pairs["pa"][1]}\n')
        # polyA: a record holding the k-mer A^31 17,000 times and a record of unique sequence,
        # so that a query_coords row holds 17,000 coordinates at one position (test 24)
        cls.polya_unique = ''.join(rng.choice('ACGT') for _ in range(400))
        cls.fasta['pp'] = f'{d}/pp.fa'
        with open(cls.fasta['pp'], 'w') as f:
            f.write(f'>pp0\n{"A" * (17000 + cls.k - 1)}\n>pp1\n{cls.polya_unique}\n')
        cls.pairs['pp'] = cls._coord_pair(d, 'pp')
        # names whose unit lines sort otherwise than the answer's order (test 26): p with two
        # pairs listed against their paths' order (pb's, then pa's), and p+, whose lines sort
        # before p's ('+' before ',') while the name sorts after it
        cls.order_csv = f'{d}/order.csv'
        with open(cls.order_csv, 'w') as f:
            for name, pair in (('p', 'pb'), ('p', 'pa'), ('p+', 'pc')):
                f.write(f'{name},{cls.pairs[pair][0]},{cls.pairs[pair][1]}\n')
        cls.servers = []

    @classmethod
    def _coord_pair(cls, d, name):
        graph = f'{d}/{name}' + graph_file_extension['succinct']
        anno = f'{d}/{name}_anno'
        cls._build_graph(cls.fasta[name], f'{d}/{name}', cls.k, 'succinct', 'basic')
        cls._annotate_graph(cls.fasta[name], graph, anno, 'row_diff_brwt_coord',
                            anno_type='filename')
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'annotate', '--anno-filename', '--index-header-coords', '-i', graph, '-o', anno,
            cls.fasta[name]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert res.returncode == 0, res.stderr.decode()
        return graph, anno + ROW_DIFF_BRWT_COORD

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers:
            server.stop()
        cls.tempdir.cleanup()

    def server(self, extra=(), multi=False, binary=None, csv=None, pair='pa', p=2, e=2, h=2):
        """a server with -p |p|, --threads-each |e| and --http-threads |h| (None: not given)"""
        # a log of its own (a stopped server leaves the list: its length repeats)
        TestSearchAnswerBudget.started += 1
        name = f'server_{TestSearchAnswerBudget.started}.log'
        args = [csv or self.csv] if multi else ['-i', self.pairs[pair][0], '-a', self.pairs[pair][1]]
        args += ['--threads-each', str(e)] + (['--http-threads', str(h)] if h else [])
        server = Server(args + list(extra), os.path.join(self.tempdir.name, name), threads=p,
                        binary=binary)
        self.servers.append(server)
        return server

    def bounded(self, max_mb, **kwargs):
        return self.server(['--search-max-answer-memory-mb', str(max_mb)]
                           + list(kwargs.pop('extra', ())), **kwargs)

    def stop(self, server):
        server.stop()
        self.servers.remove(server)

    # -------------------------------------------------------------------- requests

    def first_kmers(self, name, n=None):
        return ''.join(r[:self.k] for r in self.records[name][:n])

    def whole(self, name, times=1):
        return ''.join(self.records[name]) * times

    def body(self, mode, records, top_labels=10000, **more):
        """|records|: sequences (named q<i>) or (name, sequence) pairs"""
        named = [r if isinstance(r, tuple) else (f'q{i}', r) for i, r in enumerate(records)]
        fasta = ''.join(f'>{name}\n{s}\n' for name, s in named)
        return dict({'FASTA': fasta, 'top_labels': top_labels, 'discovery_fraction': 0.0},
                    **MODES[mode], **more)

    def reference(self, server, body):
        """the answer of a server without a budget: its text, written as styled() writes it
        (the expected answers below are built from it and written by styled()), and parsed"""
        r = server.post('search', body)
        self.assertEqual(200, r.status_code, r.text[:2000])
        parsed = json.loads(r.text)
        self.assertEqual(r.text, styled(parsed))
        return r.text, parsed

    def charge(self, server, before):
        """the request's charge (values at most a unit, answer, refused units, units) from the
        log line of the request after |before| bytes"""
        lines = [m for m in CHARGE.finditer(server.text()) if m.start() >= before]
        self.assertEqual(1, len(lines), server.text()[-3000:])
        return tuple(int(lines[0].group(i)) for i in (2, 3, 5, 6))

    def units(self, server, before):
        """the units refused alone after |before| bytes of the log: (record, pair, account,
        held, charge), each checked against its bound"""
        out = []
        for m in UNIT.finditer(server.text()):
            if m.start() < before:
                continue
            bound, held, charge = (int(m.group(i)) for i in (5, 6, 8))
            if m.group(4) == 'values':
                self.assertLessEqual(held, bound * MIB)
            else:
                # the unit's object, checked alone: nothing held before it
                self.assertEqual(0, held)
            self.assertGreater(held + charge, bound * MIB)
            self.assertEqual((f'{held / MIB:.1f}', f'{charge / MIB:.1f}'),
                             (m.group(7), m.group(9)))
            out.append((int(m.group(2)), m.group(3), m.group(4), held, charge))
        return out

    def answered(self, server, body):
        before = len(server.text())
        r = server.post('search', body)
        self.assertEqual(200, r.status_code, r.text[:2000])
        values, answer, refused, _ = self.charge(server, before)
        self.assertEqual(0, refused)
        return r.content, (values, answer)

    def partial(self, server, body):
        """a 200 with units refused alone: (text, units refused, (values, answer, refused,
        units))"""
        before = len(server.text())
        r = server.post('search', dict(body, partial_refusal=True))
        self.assertEqual(200, r.status_code, r.text[:2000])
        charge = self.charge(server, before)
        units = self.units(server, before)
        self.assertEqual(charge[2], len(units))
        self.assertNotIn('search_answer_too_large: ', server.text()[before:])
        return r.text, units, charge

    def refused(self, server, body, max_mb, account):
        before = len(server.text())
        r = server.post('search', body)
        self.assertEqual(400, r.status_code, r.text[:2000])
        self.assertEqual(refusal_text(max_mb, account), r.text)
        self.assertEqual('application/json', r.headers['Content-Type'])
        self.assertNotIn('Content-Encoding', r.headers)
        self.assertEqual(['error', 'account', 'code', 'max_memory_mb'], list(json.loads(r.text)))
        # the log: the account's charge that passed the bound, on top of what it held
        lines = [m for m in REFUSAL.finditer(server.text()) if m.start() >= before]
        self.assertEqual(1, len(lines), server.text()[-3000:])
        self.assertTrue(lines[0].group(2).startswith(f'the {account} account'), lines[0].group(2))
        bound, held, charge = (int(lines[0].group(i)) for i in (3, 4, 6))
        self.assertEqual(max_mb, bound)
        self.assertLessEqual(held, max_mb * MIB)
        self.assertGreater(held + charge, max_mb * MIB)
        self.assertEqual((f'{held / MIB:.1f}', f'{charge / MIB:.1f}'),
                         (lines[0].group(5), lines[0].group(7)))
        self.assertNotIn('answer budget: values at most', server.text()[before:])
        return lines[0].group(2)

    def pair_name(self, name):
        """a pair as the budget's log names it, by the name it was selected under"""
        pair = 'pt' if name in self.telomere_names else {'pa2': 'pa', 'u': 'pa'}.get(name, name)
        graph, anno = self.pairs[pair]
        return f'pair ({name}, {graph}, {anno})'

    # ------------------------------------------------------------------------ tests

    def test_capabilities(self):
        """GET /capabilities states the budget (test 22): at the top level of a single-graph
        server, in the multigraph block of a multi-graph server, null when 0; the pair searches
        at once are -p on a multi-graph server and the handler threads on a single-graph one,
        the records at once --threads-each, the handler threads --http-threads (else -p)"""
        single = self.server()
        self.assertEqual(block(), single.get('capabilities').json()['search'])
        self.stop(single)
        off = self.server(['--search-max-answer-memory-mb', '0'], p=3, e=4, h=5)
        self.assertEqual(block(None, 5, 4, 5), off.get('capabilities').json()['search'])
        self.stop(off)
        default_h = self.server(p=3, e=1, h=None)
        self.assertEqual(block(4096, 3, 1, 3), default_h.get('capabilities').json()['search'])
        self.stop(default_h)
        multi = self.server(['--search-max-answer-memory-mb', '7'], multi=True, p=3, e=4, h=5)
        caps = multi.get('capabilities').json()
        self.assertNotIn('search', caps)
        self.assertEqual(block(7, 3, 4, 5), caps['multigraph']['search'])
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
        budget (and as the base binary, where one is given), with the field or without; refused
        one MiB below, on the account the larger charge is on"""
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
            got, _ = self.answered(default, dict(body, partial_refusal=True))
            self.assertEqual(r.content, got, mode)
            needs[mode] = (values, answer)
            # the fixture passes one MiB: the bound below it refuses
            self.assertGreater(max(values, answer), MIB, mode)
            print(f'[ {mode} ] values {values} B a unit, answer {answer} B, text {len(got)} B',
                  flush=True)
        if BASE_BINARY:
            base = self.server(binary=BASE_BINARY)
            for mode in MODES:
                r = base.post('search', self.body(mode, body_records))
                self.assertEqual(plain_answers[mode], r.content, mode)
            self.stop(base)
        self.stop(unbounded)
        self.stop(default)
        for mode, (values, answer) in needs.items():
            bound = math.ceil(max(values, answer) / MIB)
            under = self.bounded(bound)
            got, _ = self.answered(under, self.body(mode, body_records))
            self.assertEqual(plain_answers[mode], got, mode)
            self.stop(under)
            over = self.bounded(bound - 1)
            account = 'values' if values > (bound - 1) * MIB else 'answer'
            for field in (False, True):
                body = self.body(mode, body_records, **({'partial_refusal': True} if field else {}))
                if field and account == 'values':
                    continue
                self.refused(over, body, bound - 1, account)
                # the same request again, the same refusal; the server answers the next request
                self.refused(over, body, bound - 1, account)
            small = self.body(mode, [self.records['pa'][0]])
            self.assertEqual(200, over.post('search', small).status_code)
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
        under = self.bounded(bound)
        self.assertEqual(expected, self.answered(under, body)[0])
        self.stop(under)
        over = self.bounded(bound - 1)
        self.refused(over, body, bound - 1, 'values' if values > (bound - 1) * MIB else 'answer')
        self.stop(over)

    def test_a_unit_over_and_one_under_in_each_mode(self):
        """test 8: every mode, two records on a single-graph server, the first's values past the
        bound of 1 MiB (1,000 records of 60 bp, 20 times over), the second's within it. With
        `partial_refusal` a 200: the first record's object exactly (no row, its entry refused
        "values"), then the second's byte for byte as the answer without a budget of a body
        holding it alone. Without the field: 400 "values", the unit named in the log"""
        over, under = ('over', self.whole('pa', 20)), ('under', self.records['pa'][0])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'])
        bounded = self.bounded(1)
        for mode in MODES:
            with self.subTest(mode=mode):
                alone, alone_parsed = self.reference(unbounded, self.body(mode, [under]))
                _, both = self.reference(unbounded, self.body(mode, [over, under]))
                text, units, (_, _, refused, built) = self.partial(bounded,
                                                                    self.body(mode, [over, under]))
                self.assertEqual([(0, 'the graph', 'values')], [u[:3] for u in units])
                self.assertEqual((1, 2), (refused, built))
                expected = {'cuts': [refused_entry(both[0]['cuts'][0], 'values', 1)],
                            'results': [], 'seq_description': 'over'}
                self.assertEqual(styled([expected] + alone_parsed), text)
                # the second record's bytes, as the server wrote them alone
                self.assertTrue(text.endswith(',' + alone[1:]))
                where = self.refused(bounded, self.body(mode, [over, under]), 1, 'values')
                self.assertEqual('the values account of record 0 on the graph', where)

    def test_align_with_the_first_pair_refused(self):
        """test 9: with `align`, the record's unit on the first pair in pair order refused: the
        record-level members (seq_description and the first pair's alignment) are those of the
        answer without a budget, as are the second pair's rows and entry"""
        body = self.body('query_coords', [('rep', self.repeat)], top_labels=1000, align=True,
                         graphs=['t00', 'u'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True,
                                csv=self.telomere_csv)
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        self.assertIn('cigar', reference[0])
        bounded = self.bounded(16, multi=True, csv=self.telomere_csv)
        text, units, _ = self.partial(bounded, body)
        self.assertEqual([(0, self.pair_name('t00'), 'values')], [u[:3] for u in units])
        self.assertEqual(styled(with_refused(reference, {(0, 't00'): ('values', None)}, 16)),
                         text)
        self.stop(bounded)

    def test_a_diagonal_of_refused_units(self):
        """test 10: three records on three pairs, record i's values past the bound on pair i
        alone (its coordinates there): each answered unit's rows and entry byte for byte as
        without a budget, in pair order; without the field 400 "values" """
        records = [(f'r{i}', self.whole(name, 2)) for i, name in enumerate(('pa', 'pb', 'pc'))]
        body = self.body('query_coords', records, graphs=['pc', 'pa', 'pb'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        bounded = self.bounded(1, multi=True)
        text, units, (_, _, refused, built) = self.partial(bounded, body)
        diagonal = {(0, 'pa'): ('values', None), (1, 'pb'): ('values', None),
                    (2, 'pc'): ('values', None)}
        self.assertEqual([(i, self.pair_name(name), 'values') for (i, name) in sorted(diagonal)],
                         [u[:3] for u in units])
        self.assertEqual((3, 9), (refused, built))
        self.assertEqual(styled(with_refused(reference, diagonal, 1)), text)
        self.assertRegex(self.refused(bounded, body, 1, 'values'),
                         r'^the values account of record [012] on pair \(p[abc], ')
        self.stop(bounded)

    def test_the_answer_account_refuses_the_request(self):
        """test 11: units that each fit with an answer that does not: 400 "answer" with the
        field and without; a unit refused and the remaining answer past the bound: 400 "answer"
        with the field ("values" without it); refused stubs alone past it (many short repeats,
        each refused): 400 "answer" """
        many = [self.first_kmers('pa')] * 30
        default = self.server()
        _, (values, answer) = self.answered(default, self.body('plain', many))
        self.stop(default)
        self.assertLess(values, MIB)
        self.assertGreater(answer, MIB)
        bounded = self.bounded(1)
        for field in ({}, {'partial_refusal': True}):
            self.assertEqual('the answer account',
                             self.refused(bounded, self.body('plain', many, **field), 1, 'answer'))
        over = self.whole('pa', 20)
        body = self.body('plain', [over] + many)
        self.assertEqual('the answer account',
                         self.refused(bounded, dict(body, partial_refusal=True), 1, 'answer'))
        self.refused(bounded, body, 1, 'values')
        self.stop(bounded)
        stubs = self.body('query_coords', [self.repeat[:60]] * 3000)
        telomere = self.bounded(1, pair='pt')
        self.assertEqual('the answer account',
                         self.refused(telomere, dict(stubs, partial_refusal=True), 1, 'answer'))
        self.refused(telomere, stubs, 1, 'values')
        self.stop(telomere)

    def test_values_beside_an_answer_side_refusal(self):
        """test 12: without the field, a unit's values past the bound beside the answer account
        past it (the records before it): always 400 "values", over 20 runs, on a single-graph
        and on a multi-graph server"""
        many = [self.first_kmers('pa')] * 30
        single = self.bounded(1)
        body = self.body('plain', many + [self.whole('pa', 20)])
        for _ in range(20):
            self.assertEqual('the values account of record 30 on the graph',
                             self.refused(single, body, 1, 'values'))
        self.stop(single)
        multi = self.bounded(1, multi=True)
        body = self.body('query_coords', many + [self.whole('pb', 2)], graphs=['pa', 'pb'])
        for _ in range(20):
            self.assertEqual(f'the values account of record 30 on {self.pair_name("pb")}',
                             self.refused(multi, body, 1, 'values'))
        self.stop(multi)

    def test_the_object_check(self):
        """test 13: a verbose with_signature record (1,000 rows, each with a mask of 31,000
        characters) whose object alone passes the bound, its values within it: refused alone
        "answer" with priced_mb, the unit's own object in MiB rounded up (its log line's
        charge), the same at another bound and beside another record; the other record
        answered; without the field 400 "answer" """
        record = ('masks', self.first_kmers('pa'))
        other = ('other', self.records['pa'][1])
        unbounded = self.server(['--search-max-answer-memory-mb', '0', '--verbose-output'])
        alone, alone_parsed = self.reference(unbounded, self.body('with_signature', [other]))
        masks_text, both = self.reference(unbounded, self.body('with_signature', [record, other]))
        self.stop(unbounded)
        priced = set()
        for bound in (4, 8):
            server = self.bounded(bound, extra=['--verbose-output'])
            text, units, (values, _, _, _) = self.partial(
                    server, self.body('with_signature', [record, other]))
            self.assertLessEqual(values, bound * MIB)
            [(index, where, account, held, charge)] = units
            self.assertEqual((0, 'the graph', 'answer', 0), (index, where, account, held))
            priced_mb = math.ceil(charge / MIB)
            priced.add(priced_mb)
            self.assertGreater(priced_mb, bound)
            # its tree holds at least the object's text
            self.assertGreaterEqual(charge, len(compact(both[0])))
            expected = {'cuts': [refused_entry(both[0]['cuts'][0], 'answer', bound, priced_mb)],
                        'results': [], 'seq_description': 'masks'}
            self.assertEqual(styled([expected] + alone_parsed), text)
            self.assertTrue(text.endswith(',' + alone[1:]))
            text, units, _ = self.partial(server, self.body('with_signature', [record]))
            self.assertEqual([(0, 'the graph', 'answer', 0, charge)], units)
            self.assertEqual(f'the answer account of record 0 on the graph',
                             self.refused(server, self.body('with_signature', [record, other]),
                                          bound, 'answer'))
            self.stop(server)
        self.assertEqual(1, len(priced), priced)
        print(f'[ object check ] priced_mb {priced}, text {len(masks_text)} B', flush=True)

    def test_every_unit_refused(self):
        """test 14: every unit refused, with the field: 200, every entry refused, no row (two
        records on two pairs; two records on a single-graph server)"""
        record = self.whole('pa', 2) + self.whole('pb', 2)
        body = self.body('query_coords', [('r0', record), ('r1', record)], graphs=['pa', 'pb'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        every = {(i, name): ('values', None) for i in (0, 1) for name in ('pa', 'pb')}
        bounded = self.bounded(1, multi=True)
        text, units, (_, _, refused, built) = self.partial(bounded, body)
        self.assertEqual((4, 4), (refused, built))
        expected = with_refused(reference, every, 1)
        self.assertEqual([[], []], [obj['results'] for obj in expected])
        self.assertEqual(styled(expected), text)
        self.stop(bounded)
        body = self.body('plain', [self.whole('pa', 20)] * 2)
        unbounded = self.server(['--search-max-answer-memory-mb', '0'])
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        bounded = self.bounded(1)
        text, units, _ = self.partial(bounded, body)
        self.assertEqual([(0, 'the graph', 'values'), (1, 'the graph', 'values')],
                         [u[:3] for u in units])
        self.assertEqual(styled(with_refused(reference, {(0, None): ('values', None),
                                                          (1, None): ('values', None)}, 1)),
                         text)
        self.stop(bounded)

    def test_precedence_beside_another_pairs_error(self):
        """test 15 (section 3.3): a unit refused on one pair, another pair's error (a coordinate
        query on a pair without coordinates): with the field that pair's error, as without a
        budget; without it 400 "values" """
        body = self.body('query_coords', [self.whole('pa', 2)], graphs=['pa', 'pn'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        reference = unbounded.post('search', body)
        self.stop(unbounded)
        self.assertEqual(400, reference.status_code)
        self.assertIn('Annotation does not support k-mer coordinate queries', reference.text)
        bounded = self.bounded(1, multi=True)
        before = len(bounded.text())
        r = bounded.post('search', dict(body, partial_refusal=True))
        self.assertEqual((reference.status_code, reference.content), (r.status_code, r.content))
        self.assertEqual([(0, self.pair_name('pa'), 'values')],
                         [u[:3] for u in self.units(bounded, before)])
        self.refused(bounded, body, 1, 'values')
        self.stop(bounded)

    def test_top_labels_0_and_a_pair_under_two_names(self):
        """test 16: top_labels 0 beside a refused entry (the record's sequence alone passes);
        one pair under two names (`graphs` naming both): an entry under each name, refused or
        answered alike"""
        body = self.body('plain', [self.whole('pa', 20), self.records['pa'][0]], top_labels=0)
        unbounded = self.server(['--search-max-answer-memory-mb', '0'])
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        self.assertEqual([0, 0], [obj['cuts'][0]['returned'] for obj in reference])
        bounded = self.bounded(1)
        text, units, _ = self.partial(bounded, body)
        self.assertEqual([(0, 'the graph', 'values')], [u[:3] for u in units])
        self.assertEqual(styled(with_refused(reference, {(0, None): ('values', None)}, 1)), text)
        self.stop(bounded)
        body = self.body('query_coords', [self.whole('pa', 2), self.records['pa'][0]],
                         graphs=['pa2', 'pa'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        self.assertEqual(['pa', 'pa2'], [cut['graph'] for cut in reference[1]['cuts']])
        bounded = self.bounded(1, multi=True)
        text, units, _ = self.partial(bounded, body)
        self.assertEqual([(0, self.pair_name('pa'), 'values'), (0, self.pair_name('pa2'), 'values')],
                         [u[:3] for u in units])
        self.assertEqual(styled(with_refused(reference, {(0, 'pa'): ('values', None),
                                                          (0, 'pa2'): ('values', None)}, 1)),
                         text)
        self.stop(bounded)

    def measure(self, server, body, readers):
        """(status, answer bytes, wall seconds, peak memory above the server's ready footprint
        or None, the charge line or None)"""
        warm = server.post('search', self.body('plain', [self.records['pa'][0][:40]]))
        self.assertIn(warm.status_code, (200, 400))
        time.sleep(0.5)
        ready = readers[0](server.process.pid) if readers else None
        before = len(server.text())
        start = time.time()
        r = server.post('search', body)
        wall = time.time() - start
        above = readers[1](server.process.pid) - ready if readers else None
        lines = [m for m in CHARGE.finditer(server.text()) if m.start() >= before]
        charge = tuple(int(lines[0].group(i)) for i in (2, 3, 5, 6)) if lines else None
        return r.status_code, len(r.content), wall, above, charge

    def report(self, name, numbers):
        status, size, wall, above, charge = numbers
        print(f'[ {name} ] status {status}, answer {size} B, wall {wall:.2f} s, peak above ready '
              f'{"n/a" if above is None else f"{above / MIB:.1f} MiB"}, charge {charge}',
              flush=True)

    def test_telomere_like(self):
        """test 17: a (CCCTAA)n record and a unique record, query_coords, on a pair whose 50
        labels hold about 330 occurrences of each of the repeat's k-mers: with the field the
        repeat refused alone, the unique record answered byte for byte as alone; without it 400
        "values". Measured (printed): the answer, wall time and peak memory of the request of
        both records on eleven names of the pair (r237's eleven pairs), at the bound 0, the
        default and 16 MiB, and the base binary's where one is given"""
        body = self.body('query_coords', [('telomere', self.repeat), ('unique', self.unique[50:250])],
                         top_labels=1000)
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], pair='pt')
        alone, alone_parsed = self.reference(unbounded, self.body('query_coords',
                                                       [('unique', self.unique[50:250])],
                                                       top_labels=1000))
        _, both = self.reference(unbounded, body)
        self.stop(unbounded)
        self.assertEqual(50, both[0]['cuts'][0]['returned'])
        self.assertEqual(1, both[1]['cuts'][0]['returned'])
        bounded = self.bounded(16, pair='pt')
        text, units, _ = self.partial(bounded, body)
        self.assertEqual([(0, 'the graph', 'values')], [u[:3] for u in units])
        expected = {'cuts': [refused_entry(both[0]['cuts'][0], 'values', 16)], 'results': [],
                    'seq_description': 'telomere'}
        self.assertEqual(styled([expected] + alone_parsed), text)
        self.assertTrue(text.endswith(',' + alone[1:]))
        self.assertEqual('the values account of record 0 on the graph',
                         self.refused(bounded, body, 16, 'values'))
        self.stop(bounded)
        readers = memory_readers()
        eleven = dict(body, graphs=self.telomere_names, partial_refusal=True)
        runs = [('P1 bound 0', ['--search-max-answer-memory-mb', '0'], None),
                ('P1 default', [], None),
                ('P1 bound 16', ['--search-max-answer-memory-mb', '16'], None)]
        if BASE_BINARY:
            runs += [('base bound 0', ['--search-max-answer-memory-mb', '0'], BASE_BINARY),
                     ('base default', [], BASE_BINARY)]
        sizes = {}
        for name, extra, binary in runs:
            server = self.server(extra, multi=True, csv=self.telomere_csv, binary=binary)
            numbers = self.measure(server, eleven, readers)
            self.report(f'telomere x 11 pairs, {name}', numbers)
            sizes[name] = numbers[:2]
            self.stop(server)
        self.assertEqual(sizes['P1 bound 0'], sizes['P1 default'])
        self.assertEqual(200, sizes['P1 bound 16'][0])

    def test_determinism(self):
        """test 18: the same request 20 times on servers started with -p 1, 2 or 8 and
        --threads-each 1, 2 or 4: the same bytes (a partial 200 with the field, a 400 without
        it); two different requests at once, repeated: each answer its answer alone"""
        records = [(f'r{i}', self.whole(name, 2)) for i, name in enumerate(('pa', 'pb', 'pc'))]
        # and three records answered on their pair (100 rows each)
        records += [(f's{i}', self.first_kmers(name, 100))
                    for i, name in enumerate(('pa', 'pb', 'pc'))]
        body = self.body('query_coords', records, graphs=['pa', 'pb', 'pc'])
        texts = set()
        # without the result cache: every request is searched
        for p in (1, 2, 8):
            for e in (1, 2, 4):
                server = self.bounded(1, multi=True, p=p, e=e, extra=['--cache-size', '0'])
                for _ in range(20):
                    r = server.post('search', dict(body, partial_refusal=True))
                    texts.add((r.status_code, r.text))
                    r = server.post('search', body)
                    texts.add((r.status_code, r.text))
                self.stop(server)
        self.assertEqual({200, 400}, {status for status, _ in texts})
        self.assertEqual(2, len(texts), [text[:300] for _, text in texts])
        server = self.bounded(1, multi=True, extra=['--cache-size', '0'])
        a = dict(body, partial_refusal=True)
        b = self.body('with_signature', [self.whole('pb', 20), self.records['pb'][0]],
                      graphs=['pb'], partial_refusal=True)
        alone = {}
        for name, request in (('a', a), ('b', b)):
            r = server.post('search', request)
            self.assertEqual(200, r.status_code)
            alone[name] = r.content
        for _ in range(10):
            got = {}
            threads = [threading.Thread(target=lambda n=n, q=q: got.update(
                           {n: server.post('search', q).content})) for n, q in (('a', a), ('b', b))]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            self.assertEqual(alone, got)
        self.stop(server)

    def test_cache(self):
        """test 19 (--cache-size 1): a partly refused answer is stored and the same request is
        a hit with the same bytes; the same body without the field misses and its 400 is never
        stored; a refusal is never cached"""
        body = self.body('plain', [self.whole('pa', 20), self.records['pa'][0]])
        server = self.bounded(1, extra=['--cache-size', '1'])
        before = len(server.text())
        first = server.post('search', dict(body, partial_refusal=True))
        self.assertEqual(200, first.status_code)
        self.assertEqual(1, server.text()[before:].count('cache STORE'))
        second = server.post('search', dict(body, partial_refusal=True))
        self.assertEqual(first.content, second.content)
        self.assertEqual(1, server.text()[before:].count('cache HIT'))
        before = len(server.text())
        for _ in range(2):
            self.refused(server, body, 1, 'values')
        self.assertNotIn('cache HIT', server.text()[before:])
        self.assertNotIn('cache STORE', server.text()[before:])
        self.stop(server)

    def test_in_ram_bound_0_and_the_base_binary(self):
        """test 20: in_ram gives the bytes and the charge of mmap; the bound 0 refuses nothing
        and gives the base binary's bytes (METAGRAPH_BASE_BINARY), with the field or without;
        the default bound the same"""
        records = [(f'r{i}', self.whole(name, 2)) for i, name in enumerate(('pa', 'pb', 'pc'))]
        body = self.body('query_coords', records, graphs=['pa', 'pb', 'pc'], partial_refusal=True)
        server = self.bounded(1, multi=True, extra=['--mmap', '--mem-cap-gb', '1', '--cache-size', '0'])
        mapped, units, charge = self.partial(server, body)
        in_ram, units_in_ram, charge_in_ram = self.partial(server, dict(body, in_ram=True))
        self.assertEqual((mapped, units, charge), (in_ram, units_in_ram, charge_in_ram))
        self.stop(server)
        answers = set()
        for extra, binary in ((['--search-max-answer-memory-mb', '0'], None), ([], None),
                              (['--search-max-answer-memory-mb', '0'], BASE_BINARY or None),
                              ([], BASE_BINARY or None)):
            server = self.server(extra, multi=True, binary=binary)
            for request in (body, {k: v for k, v in body.items() if k != 'partial_refusal'}):
                r = server.post('search', request)
                self.assertEqual(200, r.status_code)
                answers.add(r.content)
            self.assertNotIn('unit refused', server.text())
            self.stop(server)
        self.assertEqual(1, len(answers))

    def test_peak_memory(self):
        """test 21: a pair's records whose values each fit (20 repeats of 300 bp on the
        telomere-like pair, about 53 MiB of coordinates each) hold about --threads-each records'
        values at once, not their sum: the peak above the ready footprint is at most 3 x the
        largest unit's values + 64 MiB (printed, and the base binary's, which held the pair's
        records together until the pair's answer was built)"""
        readers = memory_readers()
        if not readers:
            self.skipTest('no peak memory reading on this platform')
        body = self.body('query_coords', [self.repeat] * 20, top_labels=1000)
        server = self.server(pair='pt')
        status, size, wall, above, charge = numbers = self.measure(server, body, readers)
        self.report('20 repeats, one pair, P1 default', numbers)
        self.stop(server)
        self.assertEqual(200, status)
        values, answer, refused, built = charge
        self.assertEqual((0, 20), (refused, built))
        self.assertGreater(values, 32 * MIB)
        self.assertLessEqual(above, 3 * values + 64 * MIB)
        if BASE_BINARY:
            server = self.server(pair='pt', binary=BASE_BINARY)
            self.report('20 repeats, one pair, base default', self.measure(server, body, readers))
            self.stop(server)

    def test_l1_relaxation(self):
        """L1 (DECISIONS, phase 1b): a request the base refused only because a pair's records
        together passed the bound (here 20 repeats, each within 64 MiB, together about 1 GiB)
        is answered, without the field, byte for byte as without a budget"""
        body = self.body('query_coords', [self.repeat] * 20, top_labels=1000)
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], pair='pt')
        expected = unbounded.post('search', body).content
        self.stop(unbounded)
        server = self.bounded(64, pair='pt')
        got, (values, _) = self.answered(server, body)
        self.assertEqual(expected, got)
        self.assertLessEqual(values, 64 * MIB)
        self.stop(server)
        if BASE_BINARY:
            base = self.bounded(64, pair='pt', binary=BASE_BINARY)
            r = base.post('search', body)
            print(f'[ L1 ] base at 64 MiB: {r.status_code} {r.text[:200]}', flush=True)
            self.assertEqual(400, r.status_code)
            self.assertIn('the values account of the graph passed the bound', base.text())
            self.stop(base)

    def test_batched_path(self):
        """test 23: on a pair without .seqs (the batched path: one query graph for the batch,
        a record per loop iteration), plain: a unit refused alone inside a batch, the other
        record answered byte for byte as alone"""
        over, under = ('over', self.whole('pa', 20)), ('under', self.records['pa'][0])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], pair='pn')
        alone, alone_parsed = self.reference(unbounded, self.body('plain', [under]))
        _, both = self.reference(unbounded, self.body('plain', [over, under]))
        self.stop(unbounded)
        bounded = self.bounded(1, pair='pn')
        text, units, _ = self.partial(bounded, self.body('plain', [over, under]))
        self.assertEqual([(0, 'the graph', 'values')], [u[:3] for u in units])
        expected = {'cuts': [refused_entry(both[0]['cuts'][0], 'values', 1)], 'results': [],
                    'seq_description': 'over'}
        self.assertEqual(styled([expected] + alone_parsed), text)
        self.assertTrue(text.endswith(',' + alone[1:]))
        self.refused(bounded, self.body('plain', [over, under]), 1, 'values')
        self.stop(bounded)

    def test_multi_graph(self):
        """a multi-graph server refuses the whole request when the pairs' answers together
        pass the bound, each pair alone within it ("answer"); and when one unit's values pass
        it ("values", the record and its pair in the log)"""
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
        server = self.bounded(bound, multi=True)
        for graphs in (['pa'], ['pb']):
            r = server.post('search', self.body('query_coords', body_records, graphs=graphs))
            self.assertEqual(200, r.status_code)
        for field in ({}, {'partial_refusal': True}):
            self.assertEqual('the answer account', self.refused(
                server, self.body('query_coords', body_records, graphs=['pa', 'pb'], **field),
                bound, 'answer'))
        self.stop(server)
        # one unit's values past the bound: the coordinates of every record of one pair, twice
        # (the values come first, whatever the answer)
        long_query = [self.whole('pa', 2)]
        default = self.server(multi=True)
        _, (values, answer) = self.answered(
            default, self.body('query_coords', long_query, graphs=['pa', 'pb']))
        self.stop(default)
        self.assertGreater(values, MIB)
        bound = math.ceil(values / MIB) - 1
        server = self.bounded(bound, multi=True)
        self.assertRegex(self.refused(server, self.body('query_coords', long_query,
                                                        graphs=['pa', 'pb']), bound, 'values'),
                         r'^the values account of record 0 on pair \(pa, ')
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
        account = 'values' if values > bound * MIB else 'answer'
        server = self.bounded(bound, multi=True,
                              extra=['--mmap', '--mem-cap-gb', '1', '--cache-size', '1'])
        for request in (body, dict(body, in_ram=True), body):
            self.refused(server, request, bound, account)
        self.assertNotIn('cache HIT', server.text())
        self.assertNotIn('cache STORE', server.text())
        small = self.body('with_signature', [self.records['pa'][0]], graphs=['pa'])
        first = server.post('search', small)
        second = server.post('search', small)
        self.assertEqual((200, first.content), (second.status_code, second.content))
        self.assertIn('cache HIT', server.text())
        self.stop(server)

    def test_the_measuring_scratch_alone_past_the_bound(self):
        """test 24 (review m1): query_coords, a record whose row holds 17,000 coordinates at one
        position (polyA): its values within the bound of 1 MiB, the scratch of measuring its
        object (the open coordinate ranges, 2.5 buffers of 32,768 ranges of 24 B: 1,966,080 B)
        past it. With the field the unit is refused alone "answer" before it is measured,
        priced_mb that scratch in MiB rounded up (2), the other record answered byte for byte as
        alone; without it 400 "answer", the log naming the unit"""
        polya, unique = ('polya', 'A' * self.k), ('unique', self.polya_unique[50:250])
        body = self.body('query_coords', [polya, unique], top_labels=1000)
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], pair='pp')
        alone, alone_parsed = self.reference(unbounded, self.body('query_coords', [unique],
                                                                   top_labels=1000))
        _, both = self.reference(unbounded, body)
        self.stop(unbounded)
        self.assertEqual(1, both[0]['cuts'][0]['returned'])
        default = self.server(pair='pp')
        _, (values, _) = self.answered(default, body)
        self.stop(default)
        self.assertLess(values, MIB)
        scratch = 1966080
        bounded = self.bounded(1, pair='pp')
        text, units, (_, answer, refused, built) = self.partial(bounded, body)
        self.assertEqual([(0, 'the graph', 'answer', 0, scratch)], units)
        self.assertEqual((1, 2), (refused, built))
        self.assertLess(answer, MIB)
        expected = {'cuts': [refused_entry(both[0]['cuts'][0], 'answer', 1, 2)], 'results': [],
                    'seq_description': 'polya'}
        self.assertEqual(styled([expected] + alone_parsed), text)
        self.assertTrue(text.endswith(',' + alone[1:]))
        before = len(bounded.text())
        self.assertEqual('the answer account of record 0 on the graph',
                         self.refused(bounded, body, 1, 'answer'))
        [line] = [m for m in REFUSAL.finditer(bounded.text()) if m.start() >= before]
        self.assertEqual((0, scratch), (int(line.group(4)), int(line.group(6))))
        self.stop(bounded)

    def test_the_answer_account_beside_another_pairs_error(self):
        """test 25 (review m2, section 3.3): with the field, the request's answer account past
        the bound while the pairs work (the objects of 100 records of 100 k-mers on pa, each
        within it, no unit refused) beside another pair's error (a coordinate query on pn): 400
        "answer", not the pair's error, which a server without a budget answers. (A unit whose own
        object passes is refused alone and charges only its stub; an answer that passes only with
        its text, priced after every pair has answered, never meets the error: such requests get
        the pair's error.)"""
        many = [self.first_kmers('pa', 100)] * 100
        body = self.body('query_coords', many, graphs=['pa', 'pn'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        reference = unbounded.post('search', body)
        self.stop(unbounded)
        self.assertEqual(400, reference.status_code)
        self.assertIn('Annotation does not support k-mer coordinate queries', reference.text)
        default = self.server(multi=True)
        _, (values, answer) = self.answered(default, dict(body, graphs=['pa']))
        self.stop(default)
        self.assertLess(values, MIB)
        self.assertGreater(answer, MIB)
        bounded = self.bounded(1, multi=True)
        for field in ({'partial_refusal': True}, {}):
            before = len(bounded.text())
            self.assertEqual('the answer account',
                             self.refused(bounded, dict(body, **field), 1, 'answer'))
            self.assertEqual([], self.units(bounded, before))
        self.stop(bounded)

    def test_the_unit_lines_in_the_answers_order(self):
        """test 26 (review m3): the unit lines in the order of the cuts (by record, then the
        pair's place in the answer), not in their text's byte order: the name p with two pairs
        listed against their paths' order (pb's, then pa's) beside p+ (pc's), whose lines sort
        first; two records, each refused alone on every pair"""
        record = self.whole('pa', 2) + self.whole('pb', 2) + self.whole('pc', 2)
        body = self.body('query_coords', [('r0', record), ('r1', record)], graphs=['p+', 'p'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True,
                                csv=self.order_csv)
        _, reference = self.reference(unbounded, body)
        self.stop(unbounded)
        order = (('p', 'pb'), ('p', 'pa'), ('p+', 'pc'))
        self.assertEqual([(name, self.pairs[pair][0]) for name, pair in order],
                         [(cut['graph'], cut['graph_path']) for cut in reference[0]['cuts']])
        wheres = [f'pair ({name}, {self.pairs[pair][0]}, {self.pairs[pair][1]})'
                  for name, pair in order]
        self.assertEqual(wheres[::-1], sorted(wheres))
        bounded = self.bounded(1, multi=True, csv=self.order_csv)
        text, units, (_, _, refused, built) = self.partial(bounded, body)
        self.assertEqual((6, 6), (refused, built))
        self.assertEqual([(i, where, 'values') for i in (0, 1) for where in wheres],
                         [u[:3] for u in units])
        every = {(i, name): ('values', None) for i in (0, 1) for name in ('p', 'p+')}
        self.assertEqual(styled(with_refused(reference, every, 1)), text)
        self.stop(bounded)

    def test_partial_refusal_is_a_boolean(self):
        """test 27 (review m4): partial_refusal true or false, or absent; any other value (a
        number, a string, null, an array) refused 400 {"error": ...} (the generic refusal of a
        bad field) before any work, on a server with a bound and on one without"""
        body = self.body('plain', [self.records['pa'][0]])
        for server in (self.bounded(1), self.server(['--search-max-answer-memory-mb', '0'])):
            answers = set()
            for value in (True, False, None):
                r = server.post('search', body if value is None
                                else dict(body, partial_refusal=value))
                self.assertEqual(200, r.status_code, r.text[:2000])
                answers.add(r.content)
            self.assertEqual(1, len(answers))
            for value in (1, 0, 'true', None, []):
                before = len(server.text())
                r = server.post('search', dict(body, partial_refusal=value))
                self.assertEqual((400, {'error': 'partial_refusal should be a boolean (true or '
                                                 'false)'}), (r.status_code, r.json()))
                self.assertNotIn('answer budget', server.text()[before:])
            self.stop(server)

    def test_a_name_listed_twice(self):
        """test 28 (review, untested claim 7): a name listed twice in `graphs` is one name (the
        list is deduplicated before any pair is queried): one entry a record, one unit line"""
        body = self.body('query_coords', [self.whole('pa', 2), self.records['pa'][0]],
                         graphs=['pa', 'pa'])
        unbounded = self.server(['--search-max-answer-memory-mb', '0'], multi=True)
        once = unbounded.post('search', dict(body, graphs=['pa']))
        twice = unbounded.post('search', body)
        self.stop(unbounded)
        self.assertEqual((200, once.content), (twice.status_code, twice.content))
        self.assertEqual([['pa'], ['pa']],
                         [[cut['graph'] for cut in obj['cuts']] for obj in json.loads(twice.text)])
        bounded = self.bounded(1, multi=True)
        _, units, (_, _, refused, built) = self.partial(bounded, body)
        self.assertEqual([(0, self.pair_name('pa'), 'values')], [u[:3] for u in units])
        self.assertEqual((1, 2), (refused, built))
        self.stop(bounded)


def write_legacy_rows(path, num_labels, rows):
    """a row annotation (RowCompressed) in its legacy form: |num_labels| labels L0000000, L0000001,
    ... (seven digits, so that their names sort as their ids), the graph's rows |rows| (each a
    list of label ids), each row's ids + 1 and a 0 after them, 20 bits each"""
    with open(path, 'wb') as out:
        out.write(struct.pack('>Q', 0))       # obsolete string map keys
        out.write(struct.pack('<QB', 0, 64))  # obsolete integer map values
        out.write(struct.pack('>Q', num_labels))
        out.write(b''.join(bytes([8]) + f'L{i:07}'.encode() for i in range(num_labels)))
        out.write(struct.pack('>QQ', len(rows), num_labels))
        values = [v for row in rows for v in [j + 1 for j in row] + [0]]
        width = 20
        assert max(values) < 1 << width
        packed = sum(value << (width * i) for i, value in enumerate(values))
        out.write(struct.pack('<QB', width * len(values), width))
        out.write(packed.to_bytes(((width * len(values) + 63) // 64) * 8, 'little'))


@unittest.skipIf(PROTEIN_MODE, "the fixture is DNA")
class TestSearchColumnStoreGrowth(TestingBase):
    """A search of a row annotation of one million labels keeps its rows' values in a hash table
    keyed by their label ids, which grows past its reserve where the ids collide (GPT final
    review of trains 5a-7, M1): the 3,500 labels 256 apart on the k-mer A^31 double its 4,096
    buckets twice, 2,752,512 B at the second doubling, the old buckets beside the new. Each of
    its allocations is charged to the unit's values account before it is made, and a doubling's
    old buckets, once freed, are used again by what the unit charges after them, its values and
    its rows' handover (the Opus review of cf77380a0, MA-1: they stayed charged, and under 3 MiB
    2,800 and 3,000 rows, which the build before M1 answered, were refused). Under 3 MiB, 2,800
    and 3,000 rows are answered in full with and without `partial_refusal`, as without a
    budget; 3,500 rows pass the bound in the request's answer account, 400 "answer" either way.
    Under 1 MiB the table's first doubling passes it in the unit's values account, before its
    array is allocated: the account holds the reserve's 4,096 buckets, 458,752 B, and the doubled
    8,192 buckets, 917,504 B, are the charge (before M1, the grown table's 1,376,256 B beyond the
    reserve was charged once the table had grown to 1,835,008 B); with `partial_refusal` the unit
    is refused alone, without it 400 "values"."""

    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        d = cls.tempdir.name
        with open(f'{d}/graph.fa', 'w') as f:
            f.write('>a\n' + 'A' * 31 + '\n>c\n' + 'C' * 31 + '\n')
        # one thread: the first k-mer, A^31, is the first row
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'build', '-p', '1', '-k', '31', '--graph', 'hash', '-o', f'{d}/graph',
            f'{d}/graph.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert res.returncode == 0, res.stderr.decode()
        cls.graph = f'{d}/graph' + graph_file_extension['hash']
        cls.anno = f'{d}/anno.row.annodbg'
        write_legacy_rows(cls.anno, 1_000_000, [[256 * r for r in range(3500)], []])
        cls.servers = []

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers:
            server.stop()
        cls.tempdir.cleanup()

    def server(self, max_mb):
        log = os.path.join(self.tempdir.name, f'server_{len(self.servers)}.log')
        server = Server(['-i', self.graph, '-a', self.anno, '--threads-each', '1',
                         '--http-threads', '1', '--cache-size', '0',
                         '--search-max-answer-memory-mb', str(max_mb)], log, threads=1)
        self.servers.append(server)
        return server

    @staticmethod
    def body(top_labels, **more):
        return dict({'FASTA': '>q\n' + 'A' * 31 + '\n', 'with_signature': True,
                     'top_labels': top_labels, 'discovery_fraction': 0.0}, **more)

    def logged(self, server, before, pattern):
        return [m for m in pattern.finditer(server.text()) if m.start() >= before]

    def test_the_table_grown_and_its_rows_within_the_bound(self):
        unbounded = self.server(0)
        bounded = self.server(3)
        for top in (2800, 3000):
            reference = unbounded.post('search', self.body(top))
            self.assertEqual(200, reference.status_code, reference.text[:2000])
            self.assertEqual(top, len(json.loads(reference.text)[0]['results']))
            for field in ({}, {'partial_refusal': True}):
                with self.subTest(top=top, **field):
                    before = len(bounded.text())
                    r = bounded.post('search', self.body(top, **field))
                    self.assertEqual((200, reference.content), (r.status_code, r.content))
                    lines = self.logged(bounded, before, CHARGE)
                    self.assertEqual(1, len(lines), bounded.text()[-3000:])
                    values, refused = int(lines[0].group(2)), int(lines[0].group(5))
                    self.assertEqual(0, refused)
                    # the second doubling's peak, held at once: within the bound
                    self.assertLessEqual(2_752_512, values)
                    self.assertLessEqual(values, 3 * MIB)

    def test_3500_rows_pass_in_the_answer_account(self):
        bounded = self.server(3)
        for field in ({}, {'partial_refusal': True}):
            with self.subTest(**field):
                before = len(bounded.text())
                r = bounded.post('search', self.body(3500, **field))
                self.assertEqual((400, refusal_text(3, 'answer')), (r.status_code, r.text))
                lines = self.logged(bounded, before, REFUSAL)
                self.assertEqual(1, len(lines), bounded.text()[-3000:])
                self.assertEqual('the answer account', lines[0].group(2))
                self.assertEqual([], self.logged(bounded, before, UNIT))

    def test_the_growth_past_the_bound_refuses_the_unit(self):
        unbounded = self.server(0)
        bounded = self.server(1)
        reference = json.loads(unbounded.post('search', self.body(3500)).text)
        before = len(bounded.text())
        r = bounded.post('search', self.body(3500, partial_refusal=True))
        self.assertEqual(200, r.status_code, r.text[:2000])
        expected = dict(reference[0], results=[],
                        cuts=[refused_entry(reference[0]['cuts'][0], 'values', 1)])
        self.assertEqual(styled([expected]), r.text)
        units = self.logged(bounded, before, UNIT)
        self.assertEqual(1, len(units), bounded.text()[-3000:])
        self.assertEqual(('0', 'the graph', 'values', '1'), units[0].group(2, 3, 4, 5))
        # refused at the first doubling, before its array is allocated: the reserve held, the
        # doubled array the charge
        self.assertEqual((458_752, 917_504), (int(units[0].group(6)), int(units[0].group(8))))
        before = len(bounded.text())
        r = bounded.post('search', self.body(3500))
        self.assertEqual((400, refusal_text(1, 'values')), (r.status_code, r.text))
        lines = self.logged(bounded, before, REFUSAL)
        self.assertEqual(1, len(lines), bounded.text()[-3000:])
        self.assertEqual('the values account of record 0 on the graph', lines[0].group(2))
        self.assertEqual((458_752, 917_504), (int(lines[0].group(4)), int(lines[0].group(6))))


if __name__ == '__main__':
    unittest.main()
