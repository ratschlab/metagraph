import unittest
import subprocess
from tempfile import TemporaryDirectory
import os


"""Test graph assemble"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'gfa': '.gfa',
                        'annotation': '.column.annodbg',
                        'fasta': ".fa"}

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

    def test_assemble_compacted_gfa(self):
        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {graph} \
                    -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension['succinct'],
            outfile=self.tempdir.name + '/annotation',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input} --compacted'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            output_gfa=self.tempdir.name + '/assembled',
            anno_input=self.tempdir.name + '/annotation' + graph_file_extension["annotation"]
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        with open(self.tempdir.name + '/assembled' + graph_file_extension["gfa"], 'r') as file:
            data = file.read()
        gfa_lines = data.rstrip("\n").split("\n")
        fields_records = {}
        for line in gfa_lines:
            if line[0] in fields_records:
                fields_records[line[0]] += 1
            else:
                fields_records[line[0]] = 1

        self.assertEqual(len(gfa_lines), 2887)
        self.assertEqual(fields_records['H'], 1)
        self.assertEqual(fields_records['S'], 1252)
        self.assertEqual(fields_records['L'], 1634)

    def test_assemble_notcompacted_gfa(self):
        construct_command = '{exe} build --mask-dummy -p {num_threads} \
                --canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {graph} \
                    -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension['succinct'],
            outfile=self.tempdir.name + '/annotation',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            output_gfa=self.tempdir.name + '/assembled',
            anno_input=self.tempdir.name + '/annotation' + graph_file_extension["annotation"]
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        with open(self.tempdir.name + '/assembled' + graph_file_extension["gfa"], 'r') as file:
            data = file.read()
        gfa_lines = data.rstrip("\n").split("\n")
        fields_records = {}
        for line in gfa_lines:
            if line[0] in fields_records:
                fields_records[line[0]] += 1
            else:
                fields_records[line[0]] = 1

        self.assertEqual(len(gfa_lines), 183551)
        self.assertEqual(fields_records['H'], 1)
        self.assertEqual(fields_records['S'], 91584)
        self.assertEqual(fields_records['L'], 91966)

    def test_round_robin_graph_size_via_gfa(self):
        k = 20
        construct_command = '{exe} build -p {num_threads} \
                -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {graph} \
                     -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension['succinct'],
            outfile=self.tempdir.name + '/annotation',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input}'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            output_gfa=self.tempdir.name + '/assembled_not_compacted',
            anno_input=self.tempdir.name + '/annotation' + graph_file_extension["annotation"]
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        assemble_command = '{exe} assemble -v {graph_input} \
                    -o {output_gfa} --unitigs --to-gfa --annotator {anno_input} --compacted'.format(
            exe=METAGRAPH,
            graph_input=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            output_gfa=self.tempdir.name + '/assembled_compacted',
            anno_input=self.tempdir.name + '/annotation' + graph_file_extension["annotation"]
        )
        res = subprocess.run([assemble_command], shell=True)
        self.assertEqual(res.returncode, 0)

        self.generate_fasta_from_gfa(
            input=self.tempdir.name + '/assembled_not_compacted' + graph_file_extension["gfa"],
            output=self.tempdir.name + '/fasta_from_not_compacted' + graph_file_extension["fasta"]
        )
        self.generate_fasta_from_gfa(
            input=self.tempdir.name + '/assembled_compacted' + graph_file_extension["gfa"],
            output=self.tempdir.name + '/fasta_from_compacted' + graph_file_extension["fasta"]
        )

        construct_command = '{exe} build -p {num_threads} \
                -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph_round_robin_not_compacted',
            input=self.tempdir.name + '/fasta_from_not_compacted' + graph_file_extension["fasta"]
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        construct_command = '{exe} build -p {num_threads} \
                -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph_round_robin_compacted',
            input=self.tempdir.name + '/fasta_from_compacted' + graph_file_extension["fasta"]
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        compare_command = '{exe} compare {fst_dbg} {scd_dbg}'.format(
            exe=METAGRAPH,
            fst_dbg=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            scd_dbg=self.tempdir.name + '/graph_round_robin_not_compacted' + graph_file_extension["succinct"]
        )
        res = subprocess.run(compare_command.split(), stdout=subprocess.PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual("Graphs are identical" in res.stdout.decode().split('\n')[2], True)

        compare_command = '{exe} compare {fst_dbg} {scd_dbg}'.format(
            exe=METAGRAPH,
            fst_dbg=self.tempdir.name + '/graph' + graph_file_extension["succinct"],
            scd_dbg=self.tempdir.name + '/graph_round_robin_compacted' + graph_file_extension["succinct"]
        )
        res = subprocess.run(compare_command.split(), stdout=subprocess.PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual("Graphs are identical" in res.stdout.decode().split('\n')[2], True)
