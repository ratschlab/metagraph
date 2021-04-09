import unittest
from parameterized import parameterized, parameterized_class
import subprocess
import itertools
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
from helpers import get_test_class_name, build_annotation, graph_file_extension, anno_file_extension, GRAPH_TYPES, ANNO_TYPES, product


"""Test graph construction"""

METAGRAPH = './metagraph'
DNA_MODE = os.readlink(METAGRAPH).endswith("_DNA")
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

NUM_THREADS = 4


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES + ['succinct_bloom', 'succinct_mask'] if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ),
    class_name_func=get_test_class_name
)
class TestQuery(unittest.TestCase):
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

        construct_command = '{exe} build {mask_dummy} -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            mask_dummy='--mask-dummy' if cls.mask_dummy else '',
            num_threads=NUM_THREADS,
            repr=cls.graph_repr,
            outfile=cls.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        assert(res.returncode == 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == params_str[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 46960' == params_str[1])
        assert('mode: basic' == params_str[2])

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

        build_annotation(
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.anno_repr,
            cls.tempdir.name + '/annotation',
            separate,
            no_fork_opt,
            no_anchor_opt
        )

        # check annotation
        anno_stats_command = '{exe} stats -a {annotation}'.format(
            exe=METAGRAPH,
            annotation=cls.tempdir.name + '/annotation' + anno_file_extension[cls.anno_repr],
        )
        res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == params_str[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 46960' == params_str[1])
        assert('representation: ' + cls.anno_repr == params_str[3])

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
            self.assertEqual(len(res.stdout), 12241)
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
            self.assertEqual(len(res.stdout), 12347)
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
            self.assertEqual(len(res.stdout), 12241)
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
            self.assertEqual(len(res.stdout), 12347)
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
        self.assertEqual(len(res.stdout), 20522)

        query_command = '{exe} query --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 20636)

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
            self.assertEqual(len(res.stdout), 12241)
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
            self.assertEqual(len(res.stdout), 12347)
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
            self.assertEqual(len(res.stdout), 12241)
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
            self.assertEqual(len(res.stdout), 12347)
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
        self.assertEqual(len(res.stdout), 20522)

        query_command = '{exe} query --fast --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 20636)

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


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        list(set(GRAPH_TYPES) - {'hashstr'}) + ['succinct_bloom', 'succinct_mask'],
        ANNO_TYPES
    ),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryCanonical(unittest.TestCase):
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

        construct_command = '{exe} build {mask_dummy} --mode canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            mask_dummy='--mask-dummy' if cls.mask_dummy else '',
            num_threads=NUM_THREADS,
            repr=cls.graph_repr,
            outfile=cls.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        assert(res.returncode == 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == params_str[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 91584' == params_str[1])
        assert('mode: canonical' == params_str[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        build_annotation(
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.anno_repr,
            cls.tempdir.name + '/annotation'
        )

        # check annotation
        anno_stats_command = '{exe} stats -a {annotation}'.format(
            exe=METAGRAPH,
            annotation=cls.tempdir.name + '/annotation' + anno_file_extension[cls.anno_repr],
        )
        res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == params_str[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 91584' == params_str[1])
        assert('representation: ' + cls.anno_repr == params_str[3])

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
    input_values=product(
        list(set(GRAPH_TYPES) - {'hashstr'}) + ['succinct_bloom', 'succinct_mask'],
        ANNO_TYPES
    ),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryPrimary(unittest.TestCase):
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

        construct_command = '{exe} build {mask_dummy} --mode canonical -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            mask_dummy='--mask-dummy' if cls.mask_dummy else '',
            num_threads=NUM_THREADS,
            repr=cls.graph_repr,
            outfile=cls.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        assert(res.returncode == 0)

        transform_command = '{exe} transform --to-fasta --primary-kmers -p {num_threads} \
                -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=cls.tempdir.name + '/graph',
            input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
        )

        res = subprocess.run([transform_command], shell=True)
        assert(res.returncode == 0)

        construct_command = '{exe} build {mask_dummy} --mode primary -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            mask_dummy='--mask-dummy' if cls.mask_dummy else '',
            num_threads=NUM_THREADS,
            repr=cls.graph_repr,
            outfile=cls.tempdir.name + '/graph',
            input=cls.tempdir.name + '/graph.fasta.gz'
        )

        res = subprocess.run([construct_command], shell=True)
        assert(res.returncode == 0)

        stats_command = '{exe} stats {graph}'.format(
            exe=METAGRAPH,
            graph=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
        )
        res = subprocess.run(stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == params_str[0])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('nodes (k): 45792' == params_str[1])
        assert('mode: primary' == params_str[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        build_annotation(
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.anno_repr,
            cls.tempdir.name + '/annotation'
        )

        # check annotation
        anno_stats_command = '{exe} stats -a {annotation}'.format(
            exe=METAGRAPH,
            annotation=cls.tempdir.name + '/annotation' + anno_file_extension[cls.anno_repr],
        )
        res = subprocess.run(anno_stats_command.split(), stdout=PIPE)
        assert(res.returncode == 0)
        params_str = res.stdout.decode().split('\n')[2:]
        assert('labels:  100' == params_str[0])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert('objects: 45792' == params_str[1])
        assert('representation: ' + cls.anno_repr == params_str[3])

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
