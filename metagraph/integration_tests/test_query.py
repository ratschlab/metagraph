import unittest
from parameterized import parameterized, parameterized_class
import subprocess
import itertools
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import numpy as np
from helpers import get_test_class_name
from base import TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension
import hashlib


"""Test graph construction"""

DNA_MODE = os.readlink(METAGRAPH).endswith("_DNA")
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")

anno_file_extension = {'column': '.column.annodbg',
                       'column_coord': '.column_coord.annodbg',
                       'brwt_coord': '.brwt_coord.annodbg',
                       'row_diff_coord': '.row_diff_coord.annodbg',
                       'row_diff_brwt_coord': '.row_diff_brwt_coord.annodbg',
                       'row': '.row.annodbg',
                       'row_diff': '.row_diff.annodbg',
                       'row_sparse': '.row_sparse.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg',
                       'row_diff_sparse': '.row_diff_sparse.annodbg',
                       'row_diff_sparse_noswap': '.row_diff_sparse.annodbg',
                       'rb_brwt': '.rb_brwt.annodbg',
                       'brwt': '.brwt.annodbg',
                       'int_brwt': '.int_brwt.annodbg',
                       'row_diff_int_brwt': '.row_diff_int_brwt.annodbg',
                       'rbfish': '.rbfish.annodbg',
                       'flat': '.flat.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]
ANNO_TYPES = [anno_type for anno_type, _ in anno_file_extension.items()]

NUM_THREADS = 4

def product(graph_types, anno_types):
    result  = []
    for graph in graph_types:
        for anno in anno_types:
            if graph == 'succinct' or not anno.startswith('row_diff'):
                result.append((graph, anno))
    return result


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ) + product(['succinct_bloom', 'succinct_mask'], ['flat']),
    class_name_func=get_test_class_name
)
class TestQuery(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic',
                         '--mask-dummy' if cls.mask_dummy else '')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == out[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 46960' == out[1])
        assert('mode: basic' == out[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt
        )

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == out[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 46960' == out[1])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(f'representation: {cls.anno_repr}' == out[3])

    def test_query(self):
        query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_both(self):
        """query graph (fwd and reverse)"""
        query_command = '{exe} query --fwd-and-reverse -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 261390)

        query_command = '{exe} query --fwd-and-reverse --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 260215)

    def test_query_parallel(self):
        """query graph (multi-threaded)"""
        query_command = '{exe} query -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_both_parallel(self):
        """query graph (fwd and reverse, multi-threaded)"""
        query_command = '{exe} query --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 261390)

        query_command = '{exe} query --fwd-and-reverse --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 260215)

    def test_query_with_align(self):
        query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12248)
        else:
            self.assertEqual(len(res.stdout), 12244)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12354)
        else:
            self.assertEqual(len(res.stdout), 12350)

        # align to graph (multi-threaded)
        query_command = '{exe} query --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12248)
        else:
            self.assertEqual(len(res.stdout), 12244)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12354)
        else:
            self.assertEqual(len(res.stdout), 12350)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_with_align_both(self):
        """align to graph (fwd and reverse multi-threaded)"""
        query_command = '{exe} query --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 24565)

        query_command = '{exe} query --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 24777)

    def test_batch_query(self):
        query_command = '{exe} query --fast -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_both(self):
        """query graph (fwd and reverse)"""
        query_command = '{exe} query --fast --fwd-and-reverse -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 261390)

        query_command = '{exe} query --fast --fwd-and-reverse --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 260215)

    def test_batch_query_parallel(self):
        """query graph (multi-threaded)"""
        query_command = '{exe} query --fast -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_both_parallel(self):
        """query graph (fwd and reverse, multi-threaded)"""
        query_command = '{exe} query --fast --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 261390)

        query_command = '{exe} query --fast --fwd-and-reverse --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 260215)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12248)
        else:
            self.assertEqual(len(res.stdout), 12244)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12354)
        else:
            self.assertEqual(len(res.stdout), 12350)

        # align to graph (multi-threaded)
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12248)
        else:
            self.assertEqual(len(res.stdout), 12244)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        if DNA_MODE:
            self.assertEqual(len(res.stdout), 12354)
        else:
            self.assertEqual(len(res.stdout), 12350)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_with_align_both(self):
        """align to graph (fwd and reverse multi-threaded)"""
        query_command = '{exe} query --fast --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 24565)

        query_command = '{exe} query --fast --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 24777)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --fast --batch-size 100 -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --fast --batch-size 100 --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    def test_query_coordinates(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --query-coords \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 0.05 {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 139268)

        query_command = f'{METAGRAPH} query --query-coords \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 0.95 {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 31522)

    def test_query_coordinates_expanded(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --query-coords --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 0.05 {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 1619883)

        query_command = f'{METAGRAPH} query --query-coords --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 0.95 {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 492788)


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(['succinct'], ANNO_TYPES),
    class_name_func=get_test_class_name
)
class TestQueryTinyLinear(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.fasta_graph = cls.tempdir.name + '/file_graph.fa'
        with open(cls.fasta_graph, 'w') as f:
            f.write(f'>L\nAAAACCCCGGGGTTTT\n')

        cls.fasta_anno = cls.tempdir.name + '/file_anno.fa'
        with open(cls.fasta_anno, 'w') as f:
            f.write(f'>L1\nAAAACCCCGGGGTTTT\n')
            f.write(f'>L2\nAAAACCCCGGGGTTTT\n')
            f.write(f'>L3\nCCCCGGGG\n')

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(cls.fasta_graph, cls.tempdir.name + '/graph', 5, cls.graph_repr, 'basic', '--mask-dummy')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 5' == out[0])
        assert('nodes (k): 12' == out[1])
        assert('mode: basic' == out[2])

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(cls.fasta_anno,
                cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
                cls.tempdir.name + '/annotation', cls.anno_repr,
                separate, no_fork_opt, no_anchor_opt)

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  3' == out[0])
        assert('objects: 12' == out[1])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(f'representation: {cls.anno_repr}' == out[3])

    def test_query_coordinates(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --query-coords  --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 0.05 {self.fasta_graph}'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(set(res.stdout.decode().strip().split('\t')),
                         {'0', 'L',
                            '<L1>:0:1:2:3:4:5:6:7:8:9:10:11',
                            '<L2>:0:1:2:3:4:5:6:7:8:9:10:11',
                            '<L3>:::::0:1:2:3::::'})


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ) + product(['succinct_bloom', 'succinct_mask'], ['flat']),
    class_name_func=get_test_class_name
)
class TestQuery1Column(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic',
                         '--mask-dummy' if cls.mask_dummy else '')

        res = cls._get_stats(cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr])
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == out[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 46960' == out[1])
        assert('mode: basic' == out[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt,
            anno_type='label 1'
        )

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  1' == out[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 46960' == out[1])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert('representation: ' + cls.anno_repr == out[3])

    def test_query(self):
        query_command = f'{METAGRAPH} query \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 1.0 \
                            {TEST_DATA_DIR}/transcripts_1000.fa'
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(hashlib.sha224(res.stdout).hexdigest(), '254d173abb255a81a4ab8a685201a73de8dbad4546c378e0a645d454')

        query_command = f'{METAGRAPH} query --count-labels \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction 1.0 \
                            {TEST_DATA_DIR}/transcripts_1000.fa'
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(hashlib.sha224(res.stdout).hexdigest(), '1bd6c24373812064c3e17e73533de7b1e30baa3cca3a64b460e83cb4')


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        ['int_brwt', 'row_diff_int_brwt']
    ),
    class_name_func=get_test_class_name
)
class TestQueryCounts(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.kmer_counts_1 = {
            'AAA': 1,
            'AAC': 2,
            'ACC': 3,
            'CCC': 4,
            'CCG': 5,
            'CGG': 6,
            'GGG': 7,
            'GGT': 8,
            'GTT': 9,
            'TTT': 10,
            'TTA': 11,
            'TAA': 12,
        }
        cls.kmer_counts_2 = {
            'AAA': 11,
            'AAC': 12,
            'ACC': 13,
            'CCC': 14,
            'CCG': 15,
            'CGG': 16,
            'GGG': 17,
            'GGT': 18,
            'GTT': 19,
            'TTT': 20,
        }
        fasta_file = cls.tempdir.name + '/file.fa'
        with open(fasta_file, 'w') as f:
            for kmer, count in cls.kmer_counts_1.items():
                f.write(f'>L1\n{kmer}\n' * count)

            for kmer, count in cls.kmer_counts_2.items():
                f.write(f'>L2\n{kmer}\n' * count)

        cls.k = 3

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(fasta_file, cls.tempdir.name + '/graph',
                         cls.k, cls.graph_repr, 'basic', '--mask-dummy' if cls.mask_dummy else '')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 3' == out[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 12' == out[1])
        assert('mode: basic' == out[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(
            fasta_file,
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt
        )

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  2' == out[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 12' == out[1])
        assert('representation: ' + cls.anno_repr == out[3])

        cls.queries = [
            'AAA',
            'AAAA',
            'AAAAAAAAAAAAA',
            'CCC',
            'CCCC',
            'CCCCCCCCCCCCC',
            'TTT',
            'AAACCCGGGTTT',
            'AAACCCGGGTTTTTT',
            'AAACCCGGGTTTAAA',
            'TTTAAACCCGGG',
            'ACACACACACACATTTAAACCCGGG',
        ]

    def test_abundance_sum_query(self):
        query_file = self.tempdir.name + '/query.fa'
        for discovery_rate in np.linspace(0, 1, 5):
            expected_output = ''
            with open(query_file, 'w') as f:
                for i, s in enumerate(self.queries):
                    f.write(f'>s{i}\n{s}\n')
                    expected_output += f'{i}\ts{i}'
                    def get_count(d, kmer):
                        try:
                            return d[kmer]
                        except:
                            return 0

                    num_kmers = len(s) - self.k + 1

                    num_matches_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    count_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) for i in range(len(s) - self.k + 1)])

                    num_matches_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    count_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) for i in range(len(s) - self.k + 1)])

                    for (c, i, n) in [(count_1, 1, num_matches_1), (count_2, 0, num_matches_2)]:
                        if n >= discovery_rate * num_kmers:
                            expected_output += f'\t<L{2-i}>:{c}'

                    expected_output += '\n'

            query_command = f'{METAGRAPH} query --fast --count-kmers \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction {discovery_rate} {query_file}'

            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(res.stdout.decode(), expected_output)

    def test_count_query(self):
        query_file = self.tempdir.name + '/query.fa'
        for discovery_rate in np.linspace(0, 1, 5):
            expected_output = ''
            with open(query_file, 'w') as f:
                for i, s in enumerate(self.queries):
                    f.write(f'>s{i}\n{s}\n')
                    expected_output += f'{i}\ts{i}'
                    def get_count(d, kmer):
                        try:
                            return d[kmer]
                        except:
                            return 0

                    num_kmers = len(s) - self.k + 1

                    num_matches_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    counts_1 = [str(get_count(self.kmer_counts_1, s[i:i + self.k])) for i in range(len(s) - self.k + 1)]

                    num_matches_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    counts_2 = [str(get_count(self.kmer_counts_2, s[i:i + self.k])) for i in range(len(s) - self.k + 1)]

                    for (cs, i, n) in [(counts_1, 1, num_matches_1), (counts_2, 0, num_matches_2)]:
                        if n >= discovery_rate * num_kmers:
                            cs = ':'.join(cs)
                            expected_output += f'\t<L{2-i}>:{cs}'

                    expected_output += '\n'

            query_command = f'{METAGRAPH} query --fast --query-counts --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --discovery-fraction {discovery_rate} {query_file}'

            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self.assertEqual(res.stdout.decode(), expected_output)

        query_command = f'{METAGRAPH} query --fast --query-counts \
                        -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                        -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                        --discovery-fraction {discovery_rate} {query_file}'

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(res.stdout.decode(),
            "0\ts0\t<L1>:0=1\t<L2>:0=11\n"
            "1\ts1\t<L1>:0-1=1\t<L2>:0-1=11\n"
            "2\ts2\t<L1>:0-10=1\t<L2>:0-10=11\n"
            "3\ts3\t<L1>:0=4\t<L2>:0=14\n"
            "4\ts4\t<L1>:0-1=4\t<L2>:0-1=14\n"
            "5\ts5\t<L1>:0-10=4\t<L2>:0-10=14\n"
            "6\ts6\t<L1>:0=10\t<L2>:0=20\n"
            "7\ts7\t<L1>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9=10\t<L2>:0=11:1=12:2=13:3=14:4=15:5=16:6=17:7=18:8=19:9=20\n"
            "8\ts8\t<L1>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9-12=10\t<L2>:0=11:1=12:2=13:3=14:4=15:5=16:6=17:7=18:8=19:9-12=20\n"
            "9\ts9\t<L1>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9=10:10=11:11=12:12=1\n"
            "10\ts10\t<L1>:0=10:1=11:2=12:3=1:4=2:5=3:6=4:7=5:8=6:9=7\n"
            "11\ts11\n")

    def test_count_quantiles(self):
        query_file = self.tempdir.name + '/query.fa'
        quantiles = np.linspace(0, 1, 100)
        expected_output = ''
        with open(query_file, 'w') as f:
            for i, s in enumerate(self.queries):
                f.write(f'>s{i}\n{s}\n')
                expected_output += f'{i}\ts{i}\t<L1>'
                def get_count(d, kmer):
                    try:
                        return d[kmer]
                    except:
                        return 0
                for p in quantiles:
                    counts = [get_count(self.kmer_counts_1, s[i:i + self.k]) for i in range(len(s) - self.k + 1)]
                    expected_output += f':{np.quantile(counts, p, interpolation="lower")}'
                expected_output += f'\t<L2>'
                for p in quantiles:
                    counts = [get_count(self.kmer_counts_2, s[i:i + self.k]) for i in range(len(s) - self.k + 1)]
                    expected_output += f':{np.quantile(counts, p, interpolation="lower")}'
                expected_output += '\n'

        query_command = f'{METAGRAPH} query --fast --count-quantiles <SET_BELOW> \
                        -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                        -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                        --discovery-fraction 0.0 {query_file}'

        query_command = query_command.split()
        query_command[4] = ' '.join([str(p) for p in quantiles])
        res = subprocess.run(query_command, stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(res.stdout.decode(), expected_output)

        query_command = f'{METAGRAPH} query --fast --count-quantiles <SET_BELOW> \
                        -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                        -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                        --discovery-fraction 1.0 {query_file}'

        query_command = query_command.split()
        query_command[4] = ' '.join([str(p) for p in quantiles])
        res = subprocess.run(query_command, stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout.decode()), 5230)


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=(product(list(set(GRAPH_TYPES) - {'hashstr'}), ANNO_TYPES) +
                  product(['succinct_bloom', 'succinct_mask'], ['flat'])),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryCanonical(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'canonical',
                         '--mask-dummy' if cls.mask_dummy else '')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == out[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 91584' == out[1])
        assert('mode: canonical' == out[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr
        )

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == out[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 91584' == out[1])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert('representation: ' + cls.anno_repr == out[3])

    def test_query(self):
        query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_query_with_align(self):
        query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12969)

    def test_batch_query(self):
        query_command = '{exe} query --fast -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12969)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --fast --batch-size 100 -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --fast --batch-size 100 --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=(product(list(set(GRAPH_TYPES) - {'hashstr'}), ANNO_TYPES) +
                  product(['succinct_bloom', 'succinct_mask'], ['flat'])),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryPrimary(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'primary',
                         '--mask-dummy' if cls.mask_dummy else '')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == out[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 45792' == out[1])
        assert('mode: primary' == out[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr
        )

        # check annotation
        res = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == out[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 45792' == out[1])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert('representation: ' + cls.anno_repr == out[3])

    def test_query(self):
        query_command = '{exe} query -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_query_with_align(self):
        query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12969)

    def test_batch_query(self):
        query_command = '{exe} query --fast -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12969)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --fast --batch-size 100 -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --fast --batch-size 100 --count-labels -i {graph} -a {annotation} --discovery-fraction 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)


if __name__ == '__main__':
    unittest.main()
