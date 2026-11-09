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
- the plain forms: column and brwt (a relaxed BRWT, the lab canary's form), read under a memory
  budget without dependency rows -- a few labels cell by cell, as without a budget, more labels
  and annotate mode row by row -- and flat, whose reads stay observed only. A k-mer carried by
  100,000 records, each a label of its own, whose row alone needs more than 1 MiB to decode.
- the greedy fixture: four records whose columns a --greedy BRWT groups out of column order, so
  that its rows are read in the tree's order; annotate mode and its label cap on column, brwt,
  the greedy brwt and flat, without a budget and under budgets that do not bind.

With $METAGRAPH_BASE_BINARY (a build of release 10, before E14), each fixture's forms are also
served by that binary and the answers compared byte for byte apart from what E14 states by design.
"""

import json
import os
import random
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH


# a build before E14 to compare with (release 10), as test_traverse.py's regression panel takes it
BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')


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
    """A single-graph server_query on one (graph, annotation) pair, of this build or |binary|"""

    def __init__(self, graph, anno, log, binary=None):
        self.port = free_port()
        self.url = f'http://127.0.0.1:{self.port}'
        self.log = open(log, 'w')
        self.proc = subprocess.Popen(
            (shlex.split(METAGRAPH) if binary is None else [binary])
            + ['server_query', '-i', graph, '-a', anno,
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

    def post_text(self, route, body):
        ret = requests.post(f'{self.url}/{route}', data=json.dumps(body), timeout=300)
        return ret.status_code, ret.text

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
        cls.files = {}
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
                cls.files[form] = (f'{d}/graph.dbg', anno)
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

    @staticmethod
    def traverse_body(seed, labels, direction, bounds, mode='constrain', labels_strategy=None,
                      detail='summary'):
        strategy = {'direction': direction, 'bounds': bounds,
                    'output': {'detail': detail, 'timing': False}}
        if mode == 'annotate':
            strategy['labels'] = {'mode': 'annotate'}
        if labels_strategy:
            strategy['labels'] = dict(strategy.get('labels', {}), **labels_strategy)
        seed = {'sequence': seed, **({'labels': labels} if labels else {})}
        return {'seeds': [seed], 'strategy': strategy}

    def traverse(self, form, seed, labels, direction, bounds, mode='constrain',
                 labels_strategy=None):
        body = self.traverse_body(seed, labels, direction, bounds, mode, labels_strategy)
        out = self.servers[form].post('traverse', body)
        (result,) = out['results']
        return out, result


class RowDiffFixture:
    """The row-diff forms' fixture (module docstring): a 20,000 bp record and two parts of it"""

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
        return self.traverse(form, *self.left_walk_args(bounds, mode))

    def left_walk_args(self, bounds, mode='constrain'):
        return (self.long[-200:-140], None if mode == 'annotate' else ['long'], 'left',
                dict({'max_extension_bp': self.L}, **bounds), mode)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
class TestTraverseRowDiffForms(RowDiffFixture, FormsBase):
    """The row-diff forms (module docstring): the same walks and resolves as the controls, the
    counts read as presence, and a read that does not fit 1 MiB refused on every form"""

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
        """/resolve's explicit labels read the default decode on every form, and a discovery
        the budget-aware decode under the server's bound (SPEC-labeled-traversal-core.md §4.6;
        a request has no memory budget: a bounds.max_memory_mb is refused as an unknown
        field): the label profiles and the discovery of every form are the control's"""
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

    def test_pattern_memory_limit_on_every_form(self):
        """/pattern's memory account on the budgeted forms (SPEC-pattern-search.md §14.4): a
        pattern near the record's start, whose rows read with their row-diff dependency rows
        need more than max_memory_mb 1 holds, has them refused on every form -- in mode partial
        listed in rows_refused (phase label_discovery, reason max_memory, needed_bytes above
        available_bytes), its contexts labels_status "refused", the labels count at_least; in
        all_or_count the labels withheld (annotation_budget) --; with room (64 MiB) the same
        pattern's labels are read on every form, nothing refused"""
        pattern = {'dna': self.long[200:220]}
        for form, server in self.servers.items():
            with self.subTest(form=form):
                out = server.post('pattern', {'patterns': [pattern], 'mode': 'partial',
                                              'max_memory_mb': 64, 'output': {'labels': 'all'}})
                (entry,) = out['patterns']
                self.assertEqual('budgeted', entry['annotation'])
                self.assertEqual([], entry.get('rows_refused', []))
                self.assertEqual({'long', 'head'},
                                 {x['column'] for r in entry['results'] for x in r['labels']})
                out = server.post('pattern', {'patterns': [pattern], 'mode': 'partial',
                                              'max_memory_mb': 1, 'output': {'labels': 'all'}})
                (entry,) = out['patterns']
                self.assertTrue(entry['rows_refused'], entry)
                for row in entry['rows_refused']:
                    self.assertEqual(('label_discovery', 'max_memory'),
                                     (row['phase'], row['reason']))
                    self.assertGreater(row['needed_bytes'], row['available_bytes'])
                self.assertTrue(any(r['labels_status'] == 'refused' for r in entry['results']))
                self.assertEqual('at_least', entry['counts']['labels']['relation'])
                out = server.post('pattern', {'patterns': [pattern], 'mode': 'all_or_count',
                                              'max_memory_mb': 1, 'output': {'labels': 'all'}})
                (entry,) = out['patterns']
                self.assertTrue(entry['rows_refused'], entry)
                self.assertEqual('annotation_budget', entry['withheld']['reason'])


class PlainFixture:
    """The plain forms' fixture (module docstring): a record 'path' whose middle 40 bases 100,000
    records carry, each a label of its own"""

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

    def right_walk(self, form, bounds, wide=False):
        """From the start of 'path' to the right, into the k-mers the 100,000 records share.
        |wide|: 16 of those records' labels permitted beside 'path' (labels.extra, reachable at
        a constant cost), 17 labels: more than the 16 a plain annotation reads cell by cell, so
        that its rows are read whole -- the walk itself is the one of 'path' alone"""
        return self.traverse(form, *self.right_walk_args(bounds, wide))

    def right_walk_args(self, bounds, wide=False):
        extra = {'extra': [f'w{i}' for i in range(16)],
                 'change_cost': {'model': 'constant', 'value': 1},
                 'loss_budget': 1} if wide else None
        return (self.P[:60], ['path'], 'right', dict({'max_extension_bp': 300}, **bounds),
                'constrain', extra)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
class TestTraversePlainForms(PlainFixture, FormsBase):
    """The plain forms (module docstring): column and brwt read under a memory budget charged
    -- a few labels cell by cell, as without a budget, more labels by rows, a row that does not
    fit refused --; flat's reads observed only; under a work budget alone all three read the
    default way"""

    def test_no_counts(self):
        for form, server in self.servers.items():
            with self.subTest(form=form):
                self.assertNotIn('counts', server.capabilities())

    def test_every_form_walks_alike_where_it_fits(self):
        """Without a budget, under a memory budget the wide row fits and under a work budget
        alone the walks of the three forms are the same, with one label and with 17"""
        for wide in (False, True):
            for bounds in ({}, {'max_memory_mb': 64}, {'max_work_units': 10 ** 12}):
                _, control = self.right_walk('flat', bounds, wide)
                self.assertEqual('complete', control['outcome']['walks'], bounds)
                for form in self.FORMS:
                    with self.subTest(wide=wide, bounds=bounds, form=form):
                        _, result = self.right_walk(form, bounds, wide)
                        self.assertEqual(walks(control), walks(result))
                        if 'max_memory_mb' in bounds:
                            self.assertIn(CHARGED if form in self.CHARGED_FORMS else OBSERVED,
                                          soft(result))

    def test_a_few_labels_are_read_cell_by_cell_under_a_budget(self):
        """At most 16 labels on a plain annotation with direct access are read cell by cell
        under a memory budget as without one (access path direct; on column and brwt charged
        cell by cell, DESIGN-traverse-graphlet.md §14.3): the k-mers of 100,000 labels are
        never decoded whole, so the walk that max_memory_mb 1 stops when its rows are read
        whole (test_a_wide_row_is_refused_on_the_plain_forms) completes, with the walks, the
        access path and the direct reads of the unbudgeted walk"""
        for form in self.FORMS:
            with self.subTest(form=form):
                _, plain = self.right_walk(form, {})
                self.assertEqual('direct', plain['annotation']['access_path'])
                self.assertGreater(plain['annotation']['direct_reads'], 0)
                for bounds in ({'max_memory_mb': 1}, {'max_memory_mb': 64},
                               {'max_memory_mb': 1, 'max_work_units': 10 ** 12}):
                    _, result = self.right_walk(form, bounds)
                    self.assertEqual('complete', result['outcome']['walks'], bounds)
                    self.assertNotIn('resource_stop', result)
                    self.assertEqual(walks(plain, counters=False),
                                     walks(result, counters=False), bounds)
                    self.assertEqual(plain['annotation'], result['annotation'], bounds)
                    self.assertIn(CHARGED if form in self.CHARGED_FORMS else OBSERVED,
                                  soft(result))

    def test_a_wide_row_is_refused_on_the_plain_forms(self):
        """Under max_memory_mb 1 the walk with 17 permitted labels (read by rows) reaches the
        k-mers of 100,000 labels: on column and brwt the read of such a row is charged as it is
        decoded and refused (phase annotation_decode, at the same depth on both: the same rows
        and demands; the row read alone, without dependency rows), on flat it is read and
        observed"""
        stopped = {}
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.right_walk(form, {'max_memory_mb': 1}, wide=True)
                self.assertEqual('rows', result['annotation']['access_path'])
                if form in self.CHARGED_FORMS:
                    self.assertEqual('partial', result['outcome']['walks'])
                    stop = result['resource_stop']
                    self.assertEqual(('memory', 'annotation_decode'),
                                     (stop['resource'], stop['phase']))
                    self.assertIn('a row of the next level read alone needs more than',
                                  stop['message'])
                    self.assertNotIn('dependency row', stop['message'])
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
        """A seed on the shared k-mers validated against 17 of their labels (read by rows): on
        column and brwt the seed phase's read of a row is admitted by its demand and refused,
        the seed failed with the demand stated; on flat the row is read (observed) and the seed
        walks. Validated against one of their labels, the seed's cells are read one by one and
        it walks on every form"""
        for form in self.FORMS:
            with self.subTest(form=form):
                _, result = self.traverse(form, self.wide[:35], [f'w{i}' for i in range(17)],
                                          'right', {'max_extension_bp': 100, 'max_memory_mb': 1})
                if form in self.CHARGED_FORMS:
                    self.assertEqual('failed', result['outcome']['walks'])
                    self.assertIn('did not admit reading the annotation row of seed k-mer 0',
                                  result['error'])
                    self.assertNotIn('dependency row', result['error'])
                    self.assertEqual(('memory', 'annotation_decode'),
                                     (result['resource_stop']['resource'],
                                      result['resource_stop']['phase']))
                else:
                    self.assertNotIn('error', result)
                _, one = self.traverse(form, self.wide[:35], ['w0'], 'right',
                                       {'max_extension_bp': 100, 'max_memory_mb': 1})
                self.assertNotIn('error', one)
                self.assertEqual('direct', one['annotation']['access_path'])

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


class GreedyFixture:
    """Four records A, B, C, D (columns 0 to 3 in that order) sharing a segment S, A and C one
    stretch X after it, B and D another, Y, each with a tail of its own: a --greedy BRWT groups
    the similar columns, {A, C} and {B, D}, so it returns a row of S as A, C, B, D (its tree's
    order), not in column order. Forms: column, brwt (relaxed, column groups in order), the
    relaxed greedy brwt and flat"""

    FORMS = ('column', 'brwt', 'brwt_greedy', 'flat')
    CHARGED_FORMS = ('column', 'brwt', 'brwt_greedy')

    @classmethod
    def write_fasta(cls, root):
        rng = random.Random(14016)

        def seq(n):
            return ''.join(rng.choice('ACGT') for _ in range(n))
        cls.S, cls.X, cls.Y = seq(120), seq(300), seq(300)
        tails = [seq(60) for _ in range(4)]
        path = f'{root}/records.fa'
        with open(path, 'w') as f:
            for name, middle, tail in zip('ABCD', (cls.X, cls.Y, cls.X, cls.Y), tails):
                f.write(f'>{name}\n{cls.S}{middle}{tail}\n')
        return path

    @classmethod
    def prepare(cls, root, fasta):
        run(['annotate', '-p', '2', '-i', f'{root}/graph.dbg', '--anno-header', '--anno-type',
             'column', '-o', f'{root}/anno', fasta], 'annotate')
        cls.column = f'{root}/anno.column.annodbg'

    @classmethod
    def build_form(cls, form, d, fasta):
        if form == 'column':
            return cls.column
        anno_type = 'brwt' if form.startswith('brwt') else form
        run(['transform_anno', '-p', '2', '--anno-type', anno_type, '-o', f'{d}/anno']
            + (['--greedy'] if form == 'brwt_greedy' else []) + [cls.column], form)
        anno = f'{d}/anno.{anno_type}.annodbg'
        if anno_type == 'brwt':
            run(['relax_brwt', '-p', '2', '-o', f'{d}/anno', anno], 'relax_brwt')
        return anno

    def annotate_walk_args(self, bounds, cap=None):
        """An annotate walk from inside S to the right, into X and Y"""
        return (self.S[40:100], None, 'right', dict({'max_extension_bp': 200}, **bounds),
                'annotate', {'max_labels_per_node': cap} if cap else None)


def without_work_units(result):
    """|result|'s walks with their counters but work_units, which a budget adds"""
    arms = json.loads(json.dumps(result.get('arms')))
    for arm in (arms or {}).values():
        arm.get('counters', {}).pop('work_units', None)
    return json.dumps({'outcome': result['outcome'], 'arms': arms,
                       'label_dict': result.get('label_dict')}, sort_keys=True)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
class TestTraverseGreedyBrwt(GreedyFixture, FormsBase):
    """Annotate mode and its label cap on the plain forms: a row is read under a budget as
    without one -- on the greedy brwt in the tree's order --, so the labels a node keeps and the
    dictionary's numbering are the unbudgeted walk's"""

    BUDGETS = ({'max_memory_mb': 4096}, {'max_memory_mb': 4096, 'max_work_units': 10 ** 12},
               {'max_work_units': 10 ** 12})

    def test_the_greedy_brwt_reads_rows_out_of_column_order(self):
        """The fixture is what it is for: the unbudgeted walk on the greedy brwt names the
        labels of S's rows in another order than on the column annotation (the same labels)"""
        names = {}
        for form in ('column', 'brwt_greedy'):
            _, result = self.traverse(form, *self.annotate_walk_args({}))
            names[form] = [label['name'] for label in result['label_dict']]
        self.assertEqual(['A', 'B', 'C', 'D'], names['column'])
        self.assertEqual(sorted(names['column']), sorted(names['brwt_greedy']))
        self.assertNotEqual(names['column'], names['brwt_greedy'])

    def test_annotate_walks_are_the_unbudgeted_ones(self):
        """Under budgets that do not bind the annotate walk of every form is its unbudgeted
        walk byte for byte (the dictionary, its numbering, the labels each node keeps, the
        counters but the budget's work_units), with every label of a node kept and with
        max_labels_per_node 2, below the four of S's rows (the first two of the row read)"""
        for form in self.FORMS:
            for cap in (None, 2):
                _, plain = self.traverse(form, *self.annotate_walk_args({}, cap))
                self.assertEqual('complete', plain['outcome']['walks'], (form, cap))
                for bounds in self.BUDGETS:
                    with self.subTest(form=form, cap=cap, bounds=bounds):
                        _, result = self.traverse(form, *self.annotate_walk_args(bounds, cap))
                        self.assertEqual(without_work_units(plain), without_work_units(result))
                        self.assertEqual(plain['annotation'], result['annotation'])
                        if 'max_memory_mb' in bounds:
                            self.assertIn(CHARGED if form in self.CHARGED_FORMS else OBSERVED,
                                          soft(result))

    def test_constrain_walks_are_the_unbudgeted_ones(self):
        """A constrain walk on two labels (read cell by cell) and on all four: the unbudgeted
        walk under every budget that does not bind, on every form"""
        for form in self.FORMS:
            for labels in (['A', 'B'], ['A', 'B', 'C', 'D']):
                args = (self.S[40:100], labels, 'right', {'max_extension_bp': 200})
                _, plain = self.traverse(form, *args)
                for bounds in self.BUDGETS:
                    with self.subTest(form=form, labels=labels, bounds=bounds):
                        _, result = self.traverse(form, args[0], labels, 'right',
                                                  dict(args[3], **bounds))
                        self.assertEqual(without_work_units(plain), without_work_units(result))
                        self.assertEqual(plain['annotation'], result['annotation'])


def diff_paths(a, b, path=''):
    """The paths (.key, [index]) at which the JSON values |a| and |b| differ"""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f'{path}.{k}')
            else:
                out += diff_paths(a[k], b[k], f'{path}.{k}')
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff_paths(x, y, f'{path}[{i}]')
        return out
    return [] if a == b else [path]


class AgainstTheBase:
    """Each form of the fixture also served by $METAGRAPH_BASE_BINARY (a build before E14), on
    the same files: the answers compared byte for byte, apart from what E14 states by design
    (by_design)"""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.base_servers = {}
        try:
            for form, (graph, anno) in cls.files.items():
                cls.base_servers[form] = Server(
                        graph, anno, os.path.join(cls.tempdir.name, form, 'server_base.log'),
                        binary=BASE_BINARY)
        except BaseException:
            cls.tearDownClass()
            raise

    @classmethod
    def tearDownClass(cls):
        for server in getattr(cls, 'base_servers', {}).values():
            server.stop()
        cls.base_servers = {}
        super().tearDownClass()

    def differ(self, form, route, body):
        """Where this build's answer differs from the base's: their paths (none: byte for byte
        the same text)"""
        base = self.base_servers[form].post_text(route, body)
        new = self.servers[form].post_text(route, body)
        self.assertEqual(base[0], new[0], (form, route, base[1][:500], new[1][:500]))
        if base[1] == new[1]:
            return []
        return diff_paths(json.loads(base[1]), json.loads(new[1]))

    def by_design(self, form, path, answer, bounds):
        """What E14 changes by design in an answer of |form| under |bounds|: the counts
        statement of a count annotation (SPEC §10.3), and on a form whose reads E14 made
        budget-aware under these budgets, memory_bound_soft (its text says the reads are
        charged; what it observed held beyond the account) and, on a row-diff form, the work
        units with the rows' dependency rows (SPEC §6.8). And on every form read budget-aware
        under a memory budget, the controls included (their reads were budget-aware before
        E14), memory_bound_soft's effect text: the final_bytes fix (E5b follow-up, train 1;
        its name-copy bound reworded in train 2) restated what the account leaves out, so the
        effect is this build's statement and differs from a base before the fix"""
        if path == '.capabilities.counts':
            return form.endswith('int_brwt')
        memory = 'max_memory_mb' in bounds
        work = 'max_work_units' in bounds
        newly = form not in ('row_diff_brwt', 'row_diff_brwt_coord')
        row_diff = form.startswith('row_diff')
        m = re.match(r'\.results\[0\]\.limitations\[(\d+)\]\.', path)
        soft = m and answer['results'][0]['limitations'][int(m.group(1))]['kind'] \
            == 'memory_bound_soft'
        if soft and memory and path.endswith('.effect'):
            effect = answer['results'][0]['limitations'][int(m.group(1))]['effect']
            if CHARGED in effect and "at most twice the seed's labels x the longest" in effect:
                return True
        charged = newly and (memory or (work and row_diff))
        if not charged:
            return False
        if soft:
            return memory
        if re.match(r'\.results\[0\]\.arms\.\w+\.counters\.work_units$', path):
            return row_diff
        return False

    def assert_as_the_base(self, form, body, bounds):
        paths = self.differ(form, 'traverse', body)
        if not paths:
            return []
        answer = self.servers[form].post('traverse', body)
        unexpected = [p for p in paths if not self.by_design(form, p, answer, bounds)]
        self.assertEqual([], unexpected, (form, bounds))
        return paths

    @staticmethod
    def report(seen):
        """The differences by design seen, on stderr (form, bounds, the path without indexes)"""
        for item in sorted(seen):
            print('as the base, by design:', *item, file=sys.stderr)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
@unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                     "$METAGRAPH_BASE_BINARY (a build before E14) is not set")
class TestTraverseRowDiffFormsAsTheBase(RowDiffFixture, AgainstTheBase, FormsBase):
    """The row-diff forms against a build before E14: /traverse without a budget and under
    budgets that do not bind, in both label modes, and /resolve, byte for byte apart from the
    statements E14 adds by design"""

    def test_traverse_as_the_base(self):
        seen = set()
        for mode in ('constrain', 'annotate'):
            for bounds in ({}, {'max_memory_mb': 4096}, {'max_work_units': 10 ** 12},
                           {'max_memory_mb': 4096, 'max_work_units': 10 ** 12}):
                for form in self.FORMS:
                    with self.subTest(mode=mode, bounds=bounds, form=form):
                        body = self.traverse_body(*self.left_walk_args(bounds, mode))
                        paths = self.assert_as_the_base(form, body, bounds)
                        if not bounds and form != 'row_diff_int_brwt':
                            self.assertEqual([], paths)
                        seen |= {(form, json.dumps(bounds, sort_keys=True),
                                  re.sub(r'\[\d+\]', '[]', p)) for p in paths}
        self.report(seen)
        # the controls answer as before under every budget, apart from memory_bound_soft's
        # effect under a memory budget (the final_bytes fix's statement, by_design)
        self.assertEqual(set(), {p for p in seen
                                 if p[0] in ('row_diff_brwt', 'row_diff_brwt_coord')
                                 and not (p[2] == '.results[].limitations[].effect'
                                          and 'max_memory_mb' in p[1])})

    def test_resolve_as_the_base(self):
        for form in self.FORMS:
            for request in ({'sequence': self.long[16000:17500], 'labels': ['long', 'tail']},
                            {'sequence': self.long[1000:4000],
                             'discover': {'max_labels': 10, 'kind': 'column'}}):
                with self.subTest(form=form, request=list(request)):
                    # (a /resolve answer always states its timing: the run's own)
                    paths = [p for p in self.differ(form, 'resolve', request)
                             if not p.startswith('.timing.')]
                    self.assertEqual(['.capabilities.counts'] if form.endswith('int_brwt')
                                     else [], paths)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
@unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                     "$METAGRAPH_BASE_BINARY (a build before E14) is not set")
class TestTraversePlainFormsAsTheBase(PlainFixture, AgainstTheBase, FormsBase):
    """The plain forms against a build before E14: one label (read cell by cell) and 17 (read
    by rows), byte for byte apart from memory_bound_soft under a memory budget on column and
    brwt"""

    def test_traverse_as_the_base(self):
        seen = set()
        for wide in (False, True):
            for bounds in ({}, {'max_memory_mb': 64}, {'max_work_units': 10 ** 12},
                           {'max_memory_mb': 64, 'max_work_units': 10 ** 12}):
                for form in self.FORMS:
                    with self.subTest(wide=wide, bounds=bounds, form=form):
                        body = self.traverse_body(*self.right_walk_args(bounds, wide))
                        paths = self.assert_as_the_base(form, body, bounds)
                        if 'max_memory_mb' not in bounds or form == 'flat':
                            self.assertEqual([], paths)
                        seen |= {(form, json.dumps(bounds, sort_keys=True),
                                  re.sub(r'\[\d+\]', '[]', p)) for p in paths}
        self.report(seen)


@unittest.skipIf(PROTEIN_MODE, "the fixtures are DNA")
@unittest.skipUnless(BASE_BINARY and os.path.isfile(BASE_BINARY),
                     "$METAGRAPH_BASE_BINARY (a build before E14) is not set")
class TestTraverseGreedyBrwtAsTheBase(GreedyFixture, AgainstTheBase, FormsBase):
    """Annotate mode and its label cap on the plain forms (the greedy brwt among them) against
    a build before E14: byte for byte apart from memory_bound_soft under a memory budget on
    column and the brwts"""

    def test_traverse_as_the_base(self):
        seen = set()
        for cap in (None, 2):
            for bounds in ({}, {'max_memory_mb': 4096}, {'max_work_units': 10 ** 12},
                           {'max_memory_mb': 4096, 'max_work_units': 10 ** 12}):
                for form in self.FORMS:
                    with self.subTest(cap=cap, bounds=bounds, form=form):
                        body = self.traverse_body(*self.annotate_walk_args(bounds, cap))
                        paths = self.assert_as_the_base(form, body, bounds)
                        if 'max_memory_mb' not in bounds or form == 'flat':
                            self.assertEqual([], paths)
                        seen |= {(form, json.dumps(bounds, sort_keys=True),
                                  re.sub(r'\[\d+\]', '[]', p)) for p in paths}
        self.report(seen)


if __name__ == '__main__':
    unittest.main()
