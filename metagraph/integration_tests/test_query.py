import unittest
from parameterized import parameterized, parameterized_class
import subprocess
import itertools
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
from helpers import get_test_class_name


"""Test graph construction"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]

NUM_THREADS = 4


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=list(itertools.product(
        GRAPH_TYPES + ['succinct_bloom', 'succinct_mask'],
        ['row', 'column']
    )),
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
        assert('canonical mode: no' == params_str[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        # build annotation
        annotate_command = '{exe} annotate --anno-header -i {graph} \
                --anno-type {anno_repr} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            anno_repr=cls.anno_repr,
            outfile=cls.tempdir.name + '/annotation',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        assert(res.returncode == 0)

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

        # query graph (fwd and reverse)
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

        # query graph (multi-threaded)
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

        # query graph (fwd and reverse, multi-threaded)
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
        query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12241)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12347)

        # align to graph (multi-threaded)
        query_command = '{exe} query --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12241)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12347)

        # align to graph (fwd and reverse multi-threaded)
        query_command = '{exe} query --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 20522)

        query_command = '{exe} query --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
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

        # query graph (fwd and reverse)
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

        # query graph (multi-threaded)
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

        # query graph (fwd and reverse, multi-threaded)
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
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12241)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12347)

        # align to graph (multi-threaded)
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12241)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} -p {num_threads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12347)

        # align to graph (fwd and reverse multi-threaded)
        query_command = '{exe} query --fast --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 20522)

        query_command = '{exe} query --fast --fwd-and-reverse --align --count-labels -i {graph} -a {annotation} -p {num_theads} --discovery-fraction 0.0 {input}'.format(
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
    input_values=list(itertools.product(
        list(set(GRAPH_TYPES) - {'hashstr'}) + ['succinct_bloom', 'succinct_mask'],
        ['row', 'column']
    )),
    class_name_func=get_test_class_name
)
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

        construct_command = '{exe} build {mask_dummy} --canonical -p {num_threads} \
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
        assert('canonical mode: yes' == params_str[2])

        if cls.with_bloom:
            convert_command = '{exe} transform -o {outfile} --initialize-bloom {bloom_param} {input}'.format(
                exe=METAGRAPH,
                outfile=cls.tempdir.name + '/graph',
                bloom_param='--bloom-fpp 0.1',
                input=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            )
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        # build annotation
        annotate_command = '{exe} annotate --anno-header -i {graph} \
                --anno-type {anno_repr} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            anno_repr=cls.anno_repr,
            outfile=cls.tempdir.name + '/annotation',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        assert(res.returncode == 0)

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
        query_command = '{exe} query --align -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
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
        query_command = '{exe} query --align --fast -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        )
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12839)

        query_command = '{exe} query --align --fast --count-labels -i {graph} -a {annotation} --discovery-fraction 0.0 {input}'.format(
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
