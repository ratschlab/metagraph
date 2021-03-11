import unittest
import subprocess
from subprocess import PIPE

from tempfile import TemporaryDirectory
import os


"""Test taxonomy classification workflow"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'
TAXO_DATA_DIR = TEST_DATA_DIR + "/taxo_data"

NUM_THREADS = 4

class TestTaxonomy(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy(self):
        k = 20
        construct_command = '{exe} build -p {num_threads} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            k=k,
            outfile=self.tempdir.name + '/graph',
            input=TAXO_DATA_DIR + '/taxo_input.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {dbg} -o {anno} -p {num_threads} {input_fasta}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            anno=self.tempdir.name + '/annotation',
            num_threads=NUM_THREADS,
            input_fasta=TAXO_DATA_DIR + '/taxo_input.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        transform_anno_tax_command = '{exe} transform_anno_tax --taxonomic-tree {tax_tree} \
                                     --lookup-table {lookup_table} -o {output} {anno}'.format(
            exe=METAGRAPH,
            tax_tree=TAXO_DATA_DIR + '/dumb_nodes.dmp',
            lookup_table=TAXO_DATA_DIR + '/dumb.accession2taxid',
            output=self.tempdir.name + '/taxoDB',
            anno=self.tempdir.name + '/annotation.column.annodbg'
        )
        res = subprocess.run([transform_anno_tax_command], shell=True)
        self.assertEqual(res.returncode, 0)

        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {taxoDB} \
                            --lca-coverage-threshold {lca_coverage}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAXO_DATA_DIR + '/taxo_query.fa',
            taxoDB=self.tempdir.name + '/taxoDB.taxo',
            lca_coverage=0.90
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        res_lines = res.stdout.decode().rstrip().split('\n')

        num_correct_predictions_tips = 0
        num_correct_predictions_internals = 0
        num_total_predictions_tips = 0
        num_total_predictions_internals = 0
        for line in res_lines:
            if line == "":
                continue
            query_expected = line.split(" ")[1].split("|")[3]
            query_prediction = line.split(" ")[7].split("'")[1]

            if line.split(" ")[1].split("|")[5] == "0":
                num_total_predictions_tips += 1
                if query_expected == query_prediction:
                    num_correct_predictions_tips += 1
                    print (num_correct_predictions_tips, ". line valid tip = ", line)
                else:
                    print (num_correct_predictions_tips, ". line wrong tip = ", line)
            else:
                num_total_predictions_internals += 1
                if query_expected == query_prediction:
                    num_correct_predictions_internals += 1
                    print (num_correct_predictions_tips, ". line valid inter = ", line)
                else:
                    print (num_correct_predictions_tips, ". line wrong inter = ", line)

        print ("num_correct_predictions_tips = ", num_correct_predictions_tips)
        print ("num_correct_predictions_internals = ", num_correct_predictions_internals)
        print ("num_total_predictions_tips = ", num_total_predictions_tips)
        print ("num_total_predictions_internals = ", num_total_predictions_internals)

        self.assertEqual(num_correct_predictions_tips, 547)
        self.assertEqual(num_correct_predictions_internals, 26)
        self.assertEqual(num_total_predictions_tips, 650)
        self.assertEqual(num_total_predictions_internals, 48)
