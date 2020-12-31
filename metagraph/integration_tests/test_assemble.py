import unittest
import subprocess
from parameterized import parameterized
from tempfile import TemporaryDirectory
import os


"""Test graph assemble"""

METAGRAPH = './metagraph'
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
        }
    },
    'not_compacted': {
        'fasta_path': TEST_DATA_DIR + '/transcripts_100.fa',
        'flag': '',
        'gfa_lines': 183551,
        'field_records': {
            'H':1,
            'S':91584,
            'L':91966,
        }
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
    def test_assemble_gfa(self, gfa_test):
        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=self.tempdir.name + '/graph',
            input=gfa_tests[gfa_test]['fasta_path']
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {graph} \
                    -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph.dbg',
            outfile=self.tempdir.name + '/annotation',
            input=gfa_tests[gfa_test]['fasta_path']
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input} {gfa_flag}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph.dbg',
            output_gfa=self.tempdir.name + '/assembled',
            anno_input=self.tempdir.name + '/annotation.column.annodbg',
            gfa_flag=gfa_tests[gfa_test]['flag']
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        with open(self.tempdir.name + '/assembled.gfa', 'r') as file:
            data = file.read()
        gfa_lines = data.rstrip("\n").split("\n")
        field_records = {}
        for line in gfa_lines:
            if line[0] in field_records:
                field_records[line[0]] += 1
            else:
                field_records[line[0]] = 1

        self.assertEqual(len(gfa_lines), gfa_tests[gfa_test]['gfa_lines'])
        self.assertEqual(field_records, gfa_tests[gfa_test]['field_records'])

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

        annotate_command = '{exe} annotate --anno-header -i {graph} \
                     -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph.dbg',
            outfile=self.tempdir.name + '/annotation',
            input=gfa_tests[gfa_test]['fasta_path']
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input} {gfa_flag}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph.dbg',
            output_gfa=self.tempdir.name + '/assembled',
            anno_input=self.tempdir.name + '/annotation.column.annodbg',
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
