import unittest
import subprocess
from parameterized import parameterized, parameterized_class
from tempfile import TemporaryDirectory
import os
import gzip
import itertools
from helpers import get_test_class_name
from base import TestingBase, graph_file_extension
from test_query import anno_file_extension, GRAPH_TYPES, ANNO_TYPES, product


"""Test graph assemble"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

gfa_tests = {
    'compacted': {
        'fasta_path': TEST_DATA_DIR + '/transcripts_100.fa',
        'flag': '--compacted',
        'gfa_lines': 2887,
        'field_records': {
            'H':1,
            'S':1252,
            'L':1634,
        },
        'expected_P_lines':100,
    },
    'not_compacted': {
        'fasta_path': TEST_DATA_DIR + '/transcripts_100.fa',
        'flag': '',
        'gfa_lines': 183551,
        'field_records': {
            'H':1,
            'S':91584,
            'L':91966,
        },
        'expected_P_lines':100,
    }
}

GFAs = [name for name, _ in gfa_tests.items()]

NUM_THREADS = 4


class TestAnnotate(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    def generate_fasta_from_gfa(self, input: str, output: str):
        fasta_file = open(output, 'w')
        with open(input, 'r') as file:
            data = file.read()
        lines = data.rstrip("\n").split("\n")

        fasta_line_limit = 80
        num_segments = 0
        for line in lines:
            if line[0] != 'S':
                continue
            num_segments += 1
            segment = line.split("\t")[2]
            fasta_file.write(">" + str(num_segments) + "\n")
            act_segment_index = 0
            while act_segment_index < len(segment):
                if len(segment) - act_segment_index > fasta_line_limit:
                    fasta_file.write(segment[act_segment_index:act_segment_index + fasta_line_limit] + '\n')
                    act_segment_index += fasta_line_limit
                else:
                    fasta_file.write(segment[act_segment_index:] + '\n')
                    break
        fasta_file.close()

    @parameterized.expand(GFAs)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_assemble_gfa(self, gfa_test):
        k = 20
        construct_command = '{exe} build -p {num_threads} --mask-dummy \
                --mode canonical -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph',
            input=gfa_tests[gfa_test]['fasta_path']
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa {gfa_flag}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph.dbg',
            output_gfa=self.tempdir.name + '/assembled',
            gfa_flag=gfa_tests[gfa_test]['flag']
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        align_command = '{exe} align -i {graph_input} {input} \
                    -o {gfa_file} {gfa_compacted_flag}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph.dbg',
            input=gfa_tests[gfa_test]['fasta_path'],
            gfa_file=self.tempdir.name + '/assembled.gfa',
            gfa_compacted_flag=gfa_tests[gfa_test]['flag']
        )
        res = subprocess.run([align_command], shell=True)
        self.assertEqual(res.returncode, 0)

        with open(self.tempdir.name + '/assembled.gfa', 'r') as file:
            data = file.read()
        gfa_lines = data.rstrip("\n").split("\n")
        field_records = {}
        sequences = {}
        for line in gfa_lines:
            if line[0] in field_records:
                field_records[line[0]] += 1
            else:
                field_records[line[0]] = 1
            if line[0] == 'S':
                sequences[line.split('\t')[1]] = line.split('\t')[2]

        self.assertEqual(len(gfa_lines), gfa_tests[gfa_test]['gfa_lines'])
        self.assertEqual(field_records, gfa_tests[gfa_test]['field_records'])

        # Ensure valid links.
        for line in gfa_lines:
            if line[0] != 'L':
                continue
            self.assertEqual(line.split('\t')[2], "+")
            self.assertEqual(line.split('\t')[4], "+")
            self.assertEqual(line.split('\t')[5], str(k-1) + "M")
            self.assertEqual(
                sequences[line.split('\t')[1]][-(k-1):],
                sequences[line.split('\t')[3]][:(k-1)]
            )

        with open(self.tempdir.name + '/assembled.path.gfa', 'r') as file:
            data = file.read()
        gfa_lines = data.rstrip("\n").split("\n")
        self.assertEqual(len(gfa_lines), gfa_tests[gfa_test]['expected_P_lines'])

        # Ensure valid paths.
        for line in gfa_lines:
            self.assertEqual(line[0], 'P')
            path_nodes = line.split('\t')[2].split(',')
            for node_idx in range(len(path_nodes) - 1):
                cur_node = path_nodes[node_idx].rstrip('+')
                nxt_node = path_nodes[node_idx + 1].rstrip('+')
                self.assertEqual(sequences[cur_node][-(k-1):], sequences[nxt_node][:(k-1)])

    @parameterized.expand(GFAs)
    def test_round_robin_graph_size_via_gfa(self, gfa_test):
        k = 20
        construct_command = '{exe} build -p {num_threads} \
                -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph',
            input=gfa_tests[gfa_test]['fasta_path']
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa {gfa_flag}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph.dbg',
            output_gfa=self.tempdir.name + '/assembled',
            gfa_flag=gfa_tests[gfa_test]['flag']
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        self.generate_fasta_from_gfa(
            input=self.tempdir.name + '/assembled.gfa',
            output=self.tempdir.name + '/assembled.fa'
        )

        construct_command = '{exe} build -p {num_threads} \
                -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph_round_robin',
            input=self.tempdir.name + '/assembled.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        compare_command = '{exe} compare {fst_dbg} {scd_dbg}'.format(
            exe=METAGRAPH,
            fst_dbg=self.tempdir.name + '/graph.dbg',
            scd_dbg=self.tempdir.name + '/graph_round_robin.dbg'
        )
        res = subprocess.run(compare_command.split(), stdout=subprocess.PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual("Graphs are identical" in res.stdout.decode().split('\n')[2], True)

@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES + ['succinct_bloom', 'succinct_mask'] if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ),
    class_name_func=get_test_class_name
)
class TestDiffAssembly(TestingBase):
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

        k = 31
        construct_command = '{exe} build {mask_dummy} -p {num_threads} \
                --graph {repr} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            mask_dummy='--mask-dummy' if cls.mask_dummy else '',
            repr=cls.graph_repr,
            k=k,
            outfile=cls.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/metasub_fake_data.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        assert(res.returncode == 0)

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
            TEST_DATA_DIR + '/metasub_fake_data.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt
        )

    def test_diff_assembly(self):
        assemble_command = '{exe} assemble -p {num_threads} \
                -a {annotation} -o {outfile} \
                --label-mask-file {mask} {graph}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=self.tempdir.name + '/diff_contigs',
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            mask=TEST_DATA_DIR + '/example.diff.json'
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        results = dict()
        with gzip.open(self.tempdir.name + '/diff_contigs.fasta.gz', 'rt') as f:
            for head, seq in itertools.zip_longest(*[f]*2):
                head = head.rstrip()
                seq = seq.rstrip()
                if head not in results:
                    results[head] = [seq]
                else:
                    results[head].append(seq)

        self.assertEqual(len(results), 2)
        self.assertTrue('>metasub_other' in results)
        self.assertTrue('>metasub_by_kmer' in results)
        self.assertEqual(len(results['>metasub_other']), 1)
        self.assertEqual(len(results['>metasub_by_kmer']), 1)
        self.assertEqual(results['>metasub_other'][0], 'CTTGGATCACACTCTTCTCAGAGCCCAGGCCAGGGGCCCCCAAGAAAGGCTCTGGTGGAGAACCTGTGCATGAAGGCTGTCAACCAGTCCATAGGCAGGGCCATCAGGCACCAAAGGGATTCTGCCAGCATAGTGCTCCTGGACCAGTGATACACCCGGCACCCTGTCCTGGACATGCTGTTGGCCTGGATCTGAGCCCTCGTGGAGGTCAAAGCCACCTTTGGTTCTGCCATTGCTGCTGTGTGGAAGTTCACTCAAGTAGGCCTCTTCCTG')
        self.assertEqual(results['>metasub_by_kmer'][0], 'CTTGGATCACACTCTTCTCAGAGCCCAGGCCAGGGGCCCCCAAGAAAGGCTCTGGTGGAGAACCTGTGCATGAAGGCTGTCAACCAGTCCATAGGCAGGGCCATCAGGCACCAAAGGGATTCTGCCAGCATAGTGCTCCTGGACCAGTGATACACCCGGCACCCTGTCCTGGACATGCTGTTGGCCTGGATCTGAGCCCTCGTGGAGGTCAAAGCCACCTTTGGTTCTGCCATTGCTGCTGTGTGGAAGTTCACTCAAGTAGGCCTCTTCCTGACAGGCAGCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCTGCCCGCTGGAGACGGTGTTTGTCATGGGCCTGGTCTGCAGG')
