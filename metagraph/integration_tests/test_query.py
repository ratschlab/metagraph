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
DNA_MODE = os.readlink(METAGRAPH).endswith("_DNA")
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg',
                       'row_diff': '.row_diff.annodbg',
                       'row_sparse': '.row_sparse.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg',
                       'row_diff_sparse': '.row_diff_sparse.annodbg',
                       'rb_brwt': '.rb_brwt.annodbg',
                       'brwt': '.brwt.annodbg',
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

def build_annotation(graph_filename, input_fasta, anno_repr, output_filename,
                     separate=False, no_fork_opt=False, no_anchor_opt=False):
    target_anno = anno_repr
    if anno_repr in {'rb_brwt', 'brwt', 'row_sparse'} or anno_repr.startswith('row_diff'):
        target_anno = anno_repr
        anno_repr = 'column'
    elif anno_repr in {'flat', 'rbfish'}:
        target_anno = anno_repr
        anno_repr = 'row'

    annotate_command = '{exe} annotate -p {num_threads} --anno-header -i {graph} \
            --anno-type {anno_repr} -o {outfile} {input}'.format(
        exe=METAGRAPH,
        num_threads=NUM_THREADS,
        graph=graph_filename,
        anno_repr=anno_repr,
        outfile=output_filename,
        input=input_fasta
    )
    res = subprocess.run([annotate_command], shell=True)
    assert(res.returncode == 0)

    if target_anno == anno_repr:
        return

    final_anno = target_anno
    if final_anno.startswith('row_diff'):
        target_anno = 'row_diff'

    annotate_command = '{exe} transform_anno -p {num_threads} \
            --anno-type {target_anno} -o {outfile} {input}'.format(
        exe=METAGRAPH,
        num_threads=NUM_THREADS,
        graph=graph_filename,
        target_anno=target_anno,
        outfile=output_filename,
        input=output_filename + anno_file_extension[anno_repr]
    )
    if target_anno == 'row_diff':
        annotate_command += ' -i ' + graph_filename

    if not no_fork_opt:
        print('-- Building RowDiff without fork optimization...')
        res = subprocess.run([annotate_command], shell=True)
        assert(res.returncode == 0)

    if target_anno == 'row_diff':
        without_input_anno = annotate_command.split(' ')
        without_input_anno.pop(-3)
        without_input_anno = ' '.join(without_input_anno)
        if not no_anchor_opt:
            if separate:
                print('-- Building RowDiff succ/pred...')
                res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 1'], shell=True)
                assert(res.returncode == 0)
            res = subprocess.run([annotate_command + ' --row-diff-stage 1'], shell=True)
            assert(res.returncode == 0)
        if separate:
            print('-- Assigning anchors...')
            res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 2'], shell=True)
            assert(res.returncode == 0)
        res = subprocess.run([annotate_command + ' --row-diff-stage 2'], shell=True)
        assert(res.returncode == 0)

        os.remove(output_filename + anno_file_extension[anno_repr])

        if final_anno != target_anno:
            annotate_command = f'{METAGRAPH} transform_anno --anno-type {final_anno} --greedy -o {output_filename} ' \
                               f'-i {graph_filename} -p {NUM_THREADS} {output_filename}.row_diff.annodbg'
            res = subprocess.run([annotate_command], shell=True)
            assert (res.returncode == 0)
            os.remove(output_filename + anno_file_extension['row_diff'])


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
