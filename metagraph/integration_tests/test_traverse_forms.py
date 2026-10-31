"""The budget-charged annotation reads on every deployed annotation form (E14;
SPEC-labeled-traversal-core.md §6.8 and §10.3, DESIGN-traverse-graphlet.md §14.3): /traverse and
/resolve of a single-graph server on each form, without a budget, under a small memory budget and
under a work budget alone.

Two fixtures, one per kind of form:
- the row-diff forms: row_diff_brwt and row_diff_brwt_coord, charged before E14, as the controls,
  and row_diff_disk (metagut, metasub41), row_diff_sparse and row_diff_flat (Tara scaffolds) and
  row_diff_int_brwt (kingsford's counts). A 20,000 bp record whose row-diff paths run to its last
  k-mer (the anchors set at most 1,000,000 rows apart), so that a row read far from the end holds
  thousands of dependency rows: a walk to the left from near the end reaches a row whose read does
  not fit 1 MiB and stops there, in phase annotation_decode.
- the plain forms: column and brwt (a relaxed BRWT, the lab canary's form), read row by row
  without dependency rows, and flat, whose reads stay observed only. A k-mer carried by 100,000
  records, each a label of its own, whose row alone needs more than 1 MiB to decode.
"""

import json
import os
import random
import shlex
import shutil
import socket
import subprocess
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH


K = 31
# the anchors of the row-diff fixtures: no row-diff path is cut short by its length
MAX_PATH_LENGTH = 1000000
# what the decoder charges states (memory_bound_soft and a memory stop's message)
CHARGED = 'charged before they are held'
OBSERVED = 'are held before they can be charged'
DECODE_STOP = 'annotation reads are charged inside the decoder'


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def run(args, what):
    res = subprocess.run(shlex.split(METAGRAPH) + args, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    assert res.returncode == 0, (what, args, res.stderr.decode()[-3000:])


class Server:
    """A single-graph server_query on one (graph, annotation) pair"""

    def __init__(self, graph, anno, log):
        self.port = free_port()
        self.url = f'http://127.0.0.1:{self.port}'
        self.log = open(log, 'w')
        self.proc = subprocess.Popen(
            shlex.split(METAGRAPH) + ['server_query', '-i', graph, '-a', anno,
                                      '--address', '127.0.0.1', '--port', str(self.port),
                                      '-p', '2'],
            stdout=self.log, stderr=subprocess.STDOUT)
        for _ in range(1200):
            if self.proc.poll() is not None:
                break
            try:
                if requests.get(self.url + '/traverse/capabilities', timeout=2).ok:
                    return
            except requests.exceptions.RequestException:
                pass
            time.sleep(0.1)
        self.stop()
        raise AssertionError(f'the server on {anno} did not start (log: {log})')

    def post(self, route, body, status=200):
        ret = requests.post(f'{self.url}/{route}', data=json.dumps(body), timeout=300)
        assert ret.status_code == status, (route, ret.status_code, ret.text[:3000])
        return ret.json()

    def capabilities(self):
        return requests.get(self.url + '/traverse/capabilities', timeout=30).json()

    def stop(self):
        self.proc.kill()
        self.proc.wait()
        self.log.close()


def walks(result, counters=True):
    """What a walk found, whatever form its annotation was read from: the arms (segments, paths,
    their labels and ends, and unless |counters| is false their counters), the label dictionary
    and the outcome"""
    arms = json.loads(json.dumps(result.get('arms')))
    if not counters:
        for arm in (arms or {}).values():
            arm.pop('counters', None)
    return json.dumps({'outcome': result['outcome'], 'arms': arms,
                       'label_dict': result.get('label_dict')}, sort_keys=True)


def soft(result):
    (lim,) = [l for l in result.get('limitations', []) if l['kind'] == 'memory_bound_soft']
    return lim['effect']


class FormsBase(TestingBase):
    """Servers on one graph's annotation forms (FORMS), each built in a directory of its own
    from a copy of the graph (the row-diff transforms write sidecars beside it)"""

    FORMS = ()

    @classmethod
    def prepare(cls, root, fasta):
        """What the forms share, built once before them (none by default)"""

    @classmethod
    def build_form(cls, form, d, fasta):
        raise NotImplementedError

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.servers = {}
        root = cls.tempdir.name
        fasta = cls.write_fasta(root)
        run(['build', '-p', '2', '--mode', 'basic', '--graph', 'succinct', '-k', str(K),
             '-o', f'{root}/graph', fasta], 'build')
        try:
            cls.prepare(root, fasta)
            for form in cls.FORMS:
                d = os.path.join(root, form)
                os.makedirs(d)
                shutil.copyfile(f'{root}/graph.dbg', f'{d}/graph.dbg')
                anno = cls.build_form(form, d, fasta)
                assert os.path.exists(anno), anno
                cls.servers[form] = Server(f'{d}/graph.dbg', anno, f'{d}/server.log')
        except BaseException:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        for server in cls.servers.values():
            server.stop()
        cls.servers = {}
        super().tearDownClass()

    def traverse(self, form, seed, labels, direction, bounds, mode='constrain'):
        strategy = {'direction': direction, 'bounds': bounds,
                    'output': {'detail': 'summary', 'timing': False}}
        if mode == 'annotate':
            strategy['labels'] = {'mode': 'annotate'}
        seed = {'sequence': seed, **({'labels': labels} if labels else {})}
        out = self.servers[form].post('traverse', {'seeds': [seed], 'strategy': strategy})
        (result,) = out['results']
        return out, result


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
class TestTraverseRowDiffForms(FormsBase):
    """The row-diff forms (module docstring): the same walks and resolves as the controls, the
    counts read as presence, and a read that does not fit 1 MiB refused on every form"""

    FORMS = ('row_diff_brwt', 'row_diff_brwt_coord', 'row_diff_disk', 'row_diff_sparse',
             'row_diff_flat', 'row_diff_int_brwt')
    CONTROL = 'row_diff_brwt'
    L = 20000

    @classmethod
    def write_fasta(cls, root):
        rng = random.Random(14014)
        cls.long = ''.join(rng.choice('ACGT') for _ in range(cls.L))
        records = f'>long\n{cls.long}\n>tail\n{cls.long[-3000:]}\n>head\n{cls.long[:3000]}\n'
        with open(f'{root}/records.fa', 'w') as f:
            f.write(records)
        # the count annotation's: 'tail' twice, so that its k-mers count 2 where 'long' counts
        # 1 (the same k-mers: the graph is the same; a coordinate annotation refuses a header
        # repeated so)
        with open(f'{root}/counted.fa', 'w') as f:
            f.write(records + f'>tail\n{cls.long[-3000:]}\n')
        return f'{root}/records.fa'

    @classmethod
    def build_form(cls, form, d, fasta):
        coords = form.endswith('_coord')
        counts = '_int_' in form
        extra = (['--coordinates'] if coords else []) + (['--count-kmers'] if counts else [])
        if counts:
            fasta = os.path.join(os.path.dirname(fasta), 'counted.fa')
        run(['annotate', '-p', '2', '-i', f'{d}/graph.dbg', '--anno-header', '--anno-type',
             'column', '-o', f'{d}/anno', fasta] + extra, 'annotate')
        os.makedirs(f'{d}/rd')
        for stage in ('0', '1', '2'):
            run(['transform_anno', '-p', '2', '--anno-type', 'row_diff', '--row-diff-stage', stage,
                 '--max-path-length', str(MAX_PATH_LENGTH), '-i', f'{d}/graph.dbg',
                 '-o', f'{d}/rd/out', f'{d}/anno.column.annodbg'] + extra, f'row_diff {stage}')
        stored = f'{d}/rd/anno.{"column" if coords or counts else "row_diff"}.annodbg'
        run(['transform_anno', '-p', '2', '--anno-type', form, '--greedy', '-i',
             f'{d}/graph.dbg', '-o', f'{d}/anno', stored], form)
        anno = f'{d}/anno.{form}.annodbg'
        if form.endswith('brwt') or form.endswith('brwt_coord'):
            run(['relax_brwt', '-p', '2', '-o', f'{d}/anno', anno], 'relax_brwt')
        return anno

    def left_walk(self, form, bounds, mode='constrain'):
        """A walk to the left from 140 bp before the record's end: each level reads rows
        farther from the anchor at its end"""
        return self.traverse(form, self.long[-200:-140],
                             None if mode == 'annotate' else ['long'], 'left',
                             dict({'max_extension_bp': self.L}, **bounds), mode)

    def test_the_counts_are_stated_as_presence_only(self):
        """capabilities.counts (SPEC §10.3): "presence_only" on the count annotation, in the
        probe and in every /traverse and /resolve response; absent on the others"""
        for form, server in self.servers.items():
            expected = 'presence_only' if form == 'row_diff_int_brwt' else None
            with self.subTest(form=form):
                out, _ = self.left_walk(form, {'max_extension_bp': 10})
                resolved = server.post('resolve', {'sequence': self.long[:100],
                                                   'labels': ['long']})
                for caps in (server.capabilities(), out['capabilities'],
                             resolved['capabilities']):
                    self.assertEqual(expected is not None, 'counts' in caps)
                    self.assertEqual(expected, caps.get('counts'))

    def test_every_form_walks_as_the_control(self):
        """Without a budget (the default reads), under a memory budget large enough (the
        budget-aware reads, never refused) and under a work budget alone: the walks of every
        form are the control's, in both label modes. The count annotation's labels are its
        presence: 'tail' counts 2 where 'long' counts 1 on the same k-mers, so its row-diff
        diffs change counts without removing columns"""
        for mode in ('constrain', 'annotate'):
            for bounds in ({}, {'max_memory_mb': 64}, {'max_work_units': 10 ** 12}):
                _, control = self.left_walk(self.CONTROL, bounds, mode)
                self.assertEqual('complete', control['outcome']['walks'], (mode, bounds))
                self.assertEqual('complete', control['arms']['left']['status'])
                if mode == 'annotate':
                    self.assertEqual({'long', 'tail', 'head'},
                                     {l['name'] for l in control['label_dict']})
                for form in self.FORMS:
                    with self.subTest(mode=mode, bounds=bounds, form=form):
                        _, result = self.left_walk(form, bounds, mode)
                        # a budget-aware read of coordinates counts them as work too
                        self.assertEqual(walks(control, form != 'row_diff_brwt_coord'),
                                         walks(result, form != 'row_diff_brwt_coord'))
                        if 'max_memory_mb' in bounds:
                            self.assertIn(CHARGED, soft(result))

    def test_a_read_that_does_not_fit_is_refused_on_every_form(self):
        """Under max_memory_mb 1 the walk to the left reaches a row whose read alone, with its
        row-diff dependency rows, does not fit what the walk has left: refused whole, the walk
        stops in phase annotation_decode, partial, every walk up to complete_to_bp present (a
        prefix of the unbudgeted walk), stated as charged inside the decoder — on every form,
        the controls' and the new ones'"""
        for mode in ('constrain', 'annotate'):
            _, full = self.left_walk(self.CONTROL, {}, mode)
            for form in self.FORMS:
                with self.subTest(mode=mode, form=form):
                    _, result = self.left_walk(form, {'max_memory_mb': 1}, mode)
                    self.assertEqual('partial', result['outcome']['walks'])
                    stop = result['resource_stop']
                    self.assertEqual(('memory', 'annotation_decode'),
                                     (stop['resource'], stop['phase']))
                    self.assertIn('row-diff dependency rows', stop['message'])
                    self.assertIn('a read that does not fit is refused whole', stop['message'])
                    self.assertIn(CHARGED, soft(result))
                    arm = result['arms']['left']
                    self.assertEqual('truncated', arm['status'])
                    # a stop thousands of rows from the end, short of the whole record
                    self.assertGreater(arm['complete_to_bp'], 100)
                    self.assertLess(arm['complete_to_bp'],
                                    full['arms']['left']['complete_to_bp'])
                    (walk,) = [l for l in arm['limitations'] if l['kind'] == 'walk_domain']
                    self.assertGreater(walk['observed'], walk['limit'])

    def test_a_work_budget_charges_the_dependency_rows_on_every_form(self):
        """Under a work budget alone the row-diff forms read budget-aware too: a row's work
        counts its dependency rows, and the stop says so"""
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.left_walk(form, {'max_work_units': 20000})
                stop = result['resource_stop']
                self.assertEqual('work', stop['resource'])
                self.assertIn('(8 per key and per dependency row, 1 per entry and coordinate',
                              stop['message'])

    def test_resolve_on_every_form_is_the_controls(self):
        """/resolve reads the default decode on every form (it has no memory budget: a
        bounds.max_memory_mb is refused as an unknown field): the label profiles and the
        discovery of every form are the control's"""
        requests_ = [{'sequence': self.long[16000:17500], 'labels': ['long', 'tail', 'head']},
                     {'sequence': self.long[16000:17500],
                      'discover': {'max_labels': 10, 'kind': 'column'}},
                     {'sequence': self.long[1000:4000],
                      'discover': {'max_labels': 10, 'kind': 'column'}}]
        for request in requests_:
            control = self.servers[self.CONTROL].post('resolve', request)
            self.assertTrue(control['labels'])
            for form, server in self.servers.items():
                with self.subTest(form=form, request=list(request)):
                    out = server.post('resolve', request)
                    for key in ('labels', 'candidates', 'num_kmers'):
                        self.assertEqual(control.get(key), out.get(key), key)
        for form, server in self.servers.items():
            with self.subTest(form=form):
                ret = server.post('resolve', {'sequence': self.long[:100], 'labels': ['long'],
                                              'bounds': {'max_memory_mb': 1}}, status=400)
                self.assertIn('max_memory_mb', ret['error'])

    def test_pattern_reads_the_new_row_diff_forms_budgeted(self):
        """/pattern: output.labels "all" on every row-diff form is read budget-aware
        (annotation "budgeted"), the count annotation's labels being its presence"""
        request = {'patterns': [{'dna': self.long[17000:17020]}], 'output': {'labels': 'all'}}
        for form, server in self.servers.items():
            with self.subTest(form=form):
                out = server.post('pattern', request)
                (entry,) = out['patterns']
                self.assertEqual('budgeted', entry['annotation'])
                labels = {x['column'] for r in entry['results'] for x in r['labels']}
                self.assertEqual({'long', 'tail'}, labels)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
class TestTraversePlainForms(FormsBase):
    """The plain forms (module docstring): column and brwt read row by row under a memory
    budget, a row that does not fit refused; flat's reads observed only; under a work budget
    alone all three read the default way"""

    FORMS = ('column', 'brwt', 'flat')
    CHARGED_FORMS = ('column', 'brwt')
    W = 100000

    @classmethod
    def write_fasta(cls, root):
        rng = random.Random(14015)

        def seq(n):
            return ''.join(rng.choice('ACGT') for _ in range(n))
        cls.P, cls.wide, cls.Q = seq(120), seq(40), seq(60)
        path = f'{root}/records.fa'
        with open(path, 'w') as f:
            f.write(f'>path\n{cls.P}{cls.wide}{cls.Q}\n')
            for i in range(cls.W):
                f.write(f'>w{i}\n{cls.wide}\n')
        return path

    @classmethod
    def prepare(cls, root, fasta):
        # one column annotation, the others converted from it (no form writes graph sidecars)
        run(['annotate', '-p', '2', '-i', f'{root}/graph.dbg', '--anno-header', '--anno-type',
             'column', '-o', f'{root}/anno', fasta], 'annotate')
        cls.column = f'{root}/anno.column.annodbg'

    @classmethod
    def build_form(cls, form, d, fasta):
        if form == 'column':
            return cls.column
        run(['transform_anno', '-p', '2', '--anno-type', form, '-o', f'{d}/anno', cls.column],
            form)
        anno = f'{d}/anno.{form}.annodbg'
        if form == 'brwt':
            # relaxed, as a lab index's .brwt.annodbg
            run(['relax_brwt', '-p', '2', '-o', f'{d}/anno', anno], 'relax_brwt')
        return anno

    def right_walk(self, form, bounds):
        """From the start of 'path' to the right, into the k-mers the 100,000 records share"""
        return self.traverse(form, self.P[:60], ['path'], 'right',
                             dict({'max_extension_bp': 300}, **bounds))

    def test_no_counts(self):
        for form, server in self.servers.items():
            with self.subTest(form=form):
                self.assertNotIn('counts', server.capabilities())

    def test_every_form_walks_alike_where_it_fits(self):
        """Without a budget, under a memory budget the wide row fits and under a work budget
        alone the walks of the three forms are the same"""
        for bounds in ({}, {'max_memory_mb': 64}, {'max_work_units': 10 ** 12}):
            _, control = self.right_walk('flat', bounds)
            self.assertEqual('complete', control['outcome']['walks'], bounds)
            for form in self.FORMS:
                with self.subTest(bounds=bounds, form=form):
                    _, result = self.right_walk(form, bounds)
                    self.assertEqual(walks(control), walks(result))
                    if 'max_memory_mb' in bounds:
                        self.assertIn(CHARGED if form in self.CHARGED_FORMS else OBSERVED,
                                      soft(result))

    def test_a_wide_row_is_refused_on_the_plain_forms(self):
        """Under max_memory_mb 1 the walk reaches the k-mers of 100,000 labels: on column and
        brwt the read of such a row is charged as it is decoded and refused (phase
        annotation_decode, at the same depth on both: the same rows and demands), on flat it is
        read and observed"""
        stopped = {}
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.right_walk(form, {'max_memory_mb': 1})
                if form in self.CHARGED_FORMS:
                    self.assertEqual('partial', result['outcome']['walks'])
                    stop = result['resource_stop']
                    self.assertEqual(('memory', 'annotation_decode'),
                                     (stop['resource'], stop['phase']))
                    self.assertIn(CHARGED, soft(result))
                    # the walk reached the k-mers before the shared ones
                    stopped[form] = result['arms']['right']['complete_to_bp']
                    self.assertEqual(len(self.P) - 60 + K - 1, stopped[form])
                else:
                    self.assertNotEqual('annotation_decode',
                                        (result.get('resource_stop') or {}).get('phase'))
                    self.assertIn(OBSERVED, soft(result))
        self.assertEqual(stopped['column'], stopped['brwt'])

    def test_a_wide_seed_row_fails_the_seed_on_the_plain_forms(self):
        """A seed on the shared k-mers validated against one of their labels: on column and
        brwt the seed phase's read of a row is admitted by its demand and refused, the seed
        failed with the demand stated; on flat the row is read (observed) and the seed walks"""
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.traverse(form, self.wide[:35], ['w0'], 'right',
                                          {'max_extension_bp': 100, 'max_memory_mb': 1})
                if form in self.CHARGED_FORMS:
                    self.assertEqual('failed', result['outcome']['walks'])
                    self.assertIn('did not admit reading the annotation row of seed k-mer 0',
                                  result['error'])
                    self.assertEqual(('memory', 'annotation_decode'),
                                     (result['resource_stop']['resource'],
                                      result['resource_stop']['phase']))
                else:
                    self.assertNotIn('error', result)

    def test_a_work_budget_alone_keeps_the_default_reads(self):
        """Under a work budget alone the plain forms are read the default way (a plain row has
        no dependency rows to charge as work, DESIGN-traverse-graphlet.md §14.3): the same work
        stop on the three forms, its weights without dependency rows"""
        stops = {}
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.right_walk(form, {'max_work_units': 300})
                stop = result['resource_stop']
                self.assertEqual('work', stop['resource'])
                self.assertIn('(8 per key, 1 per entry and coordinate; near the budget',
                              stop['message'])
                stops[form] = (walks(result), stop['message'])
        self.assertEqual(stops['flat'], stops['column'])
        self.assertEqual(stops['flat'], stops['brwt'])

    def test_resolve_on_every_form_is_the_same(self):
        request = {'sequence': self.P + self.wide + self.Q, 'labels': ['path', 'w0', 'w99999']}
        control = self.servers['flat'].post('resolve', request)
        self.assertEqual(3, len(control['labels']))
        for form, server in self.servers.items():
            with self.subTest(form=form):
                out = server.post('resolve', request)
                for key in ('labels', 'candidates', 'num_kmers'):
                    self.assertEqual(control.get(key), out.get(key), key)

    def test_pattern_keeps_the_plain_forms_unbudgeted(self):
        """/pattern keeps its contract: a plain column or brwt annotation is "unbudgeted" there
        (SPEC-pattern-search.md §14.4), output.labels "all" refused without
        allow_unbudgeted_annotation"""
        request = {'patterns': [{'dna': self.P[10:30]}], 'output': {'labels': 'all'}}
        for form, server in self.servers.items():
            with self.subTest(form=form):
                ret = server.post('pattern', request, status=400)
                self.assertEqual('annotation_unbudgeted', ret['code'])


if __name__ == '__main__':
    unittest.main()
