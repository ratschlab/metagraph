import unittest
import subprocess
from tempfile import TemporaryDirectory
import os

import sys
sys.path.append('../external-libraries/gfapy')
import gfapy


"""Test graph assemble"""

METAGRAPH = './metagraph'
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

graph_file_extension = {'succinct': '.dbg',
                        'gfa': '.gfa',
                        'annotation': '.column.annodbg'}

NUM_THREADS = 4


class TestAnnotate(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

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
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension['succinct'],
            anno_repr='column',
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

        # vlevel=3 does a validation when reading at the highest level.
        gfa_obj = gfapy.Gfa(data.rstrip("\n"), vlevel=3, version='gfa1')
        self.assertEqual(len(gfa_obj.lines), 2887)
        self.assertEqual(len(gfa_obj._records["S"]), 1252)
        self.assertEqual(len(gfa_obj._records["L"]), 1634)
        self.assertEqual(len(gfa_obj._records["P"]), 0)

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
                    --anno-type {anno_repr} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension['succinct'],
            anno_repr='column',
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

        # vlevel=3 does a validation when reading at the highest level.
        gfa_obj = gfapy.Gfa(data.rstrip("\n"), vlevel=3, version='gfa1')
        self.assertEqual(len(gfa_obj.lines), 183551)
        self.assertEqual(len(gfa_obj._records["S"]), 91584)
        self.assertEqual(len(gfa_obj._records["L"]), 91966)
        self.assertEqual(len(gfa_obj._records["P"]), 0)