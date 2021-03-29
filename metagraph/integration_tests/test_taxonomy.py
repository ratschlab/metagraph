import unittest
import subprocess
from subprocess import PIPE
from parameterized import parameterized

from tempfile import TemporaryDirectory
import os


"""Test taxonomy classification workflow"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'
TAX_DATA_DIR = TEST_DATA_DIR + "/taxonomic_data"

tax_tests = {
    'one_thread': {
        'threads': 1,
    },
    '5_threads': {
        'threads': 5,
    },
    '16_threads': {
        'threads': 16,
    }
}

test_params = [name for name, _ in tax_tests.items()]

class TestTaxonomy(unittest.TestCase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()
        self.tax_parent = {}
        self.tax_root = -1
        tax_lines = open(TAX_DATA_DIR + '/dumb_nodes.dmp').readlines()
        for line in tax_lines:
            act_node = line.split('\t')[0].strip()
            act_parent = line.split('\t')[2].strip()
            self.tax_parent[act_node] = act_parent
            if act_node == act_parent:
                self.tax_root = act_node

    def is_descendant(self, node: str, query: str) -> bool:
        while query != self.tax_root:
            query = self.tax_parent[query]
            if query == node:
                return True
        return False

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy(self, tax_test):
        k = 20
        lca_coverage = 0.9
        construct_command = '{exe} build -p {num_threads} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=tax_tests[tax_test]['threads'],
            k=k,
            outfile=self.tempdir.name + '/graph',
            input=TAX_DATA_DIR + '/tax_input.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {dbg} -o {anno} -p {num_threads} {input_fasta}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            anno=self.tempdir.name + '/annotation',
            num_threads=tax_tests[tax_test]['threads'],
            input_fasta=TAX_DATA_DIR + '/tax_input.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        transform_anno_tax_command = '{exe} transform_anno_tax --taxonomic-tree {tax_tree} \
                                     --label-taxid-map {lookup_table} -o {output} -p {num_threads} {anno}'.format(
            exe=METAGRAPH,
            tax_tree=TAX_DATA_DIR + '/dumb_nodes.dmp',
            lookup_table=TAX_DATA_DIR + '/dumb.accession2taxid',
            output=self.tempdir.name + '/taxDB',
            num_threads=tax_tests[tax_test]['threads'],
            anno=self.tempdir.name + '/annotation.column.annodbg'
        )
        res = subprocess.run([transform_anno_tax_command], shell=True)
        self.assertEqual(res.returncode, 0)

        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {taxDB} \
                            --lca-coverage-threshold {lca_coverage} -p {num_threads}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            taxDB=self.tempdir.name + '/taxDB.taxdb',
            lca_coverage=lca_coverage,
            num_threads=tax_tests[tax_test]['threads'],
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        res_lines = res.stdout.decode().rstrip().split('\n')

        num_correct_predictions_tips = 0
        num_correct_predictions_internals = 0
        num_total_predictions_tips = 0
        num_total_predictions_internals = 0

        num_descendant_predictions_internals = 0
        num_ancestor_predictions_tips = 0
        num_ancestor_predictions_internals = 0

        num_wrong_predictions_tips = 0
        num_wrong_predictions_internals = 0
        for line in res_lines:
            if line == "":
                continue
            query_expected = line.split(" ")[1].split("|")[3].strip()
            query_prediction = line.split(" ")[7].split("'")[1].strip()

            if line.split(" ")[1].split("|")[5] == "0":
                # This taxid is a tip, has no children in the taxonomic tree.
                num_total_predictions_tips += 1
                if query_expected == query_prediction:
                    num_correct_predictions_tips += 1
                else:
                    if self.is_descendant(node=query_prediction, query=query_expected):
                        num_ancestor_predictions_tips += 1
                    else:
                        num_wrong_predictions_tips += 1
            else:
                # This taxid is an internal node.
                num_total_predictions_internals += 1
                if query_expected == query_prediction:
                    num_correct_predictions_internals += 1
                else:
                    if self.is_descendant(node=query_prediction, query=query_expected):
                        num_ancestor_predictions_internals += 1
                    elif self.is_descendant(node=query_expected, query=query_prediction):
                        num_descendant_predictions_internals += 1
                    else:
                        num_wrong_predictions_internals += 1
        self.assertEqual(num_total_predictions_tips, 650)
        self.assertEqual(num_total_predictions_internals, 48)

        self.assertEqual(num_correct_predictions_tips, 547)
        self.assertEqual(num_correct_predictions_internals, 24)

        self.assertEqual(num_ancestor_predictions_internals, 13)
        self.assertEqual(num_descendant_predictions_internals, 8)
        self.assertEqual(num_ancestor_predictions_tips, 98)

        self.assertEqual(num_wrong_predictions_internals, 3)
        self.assertEqual(num_wrong_predictions_tips, 5)
