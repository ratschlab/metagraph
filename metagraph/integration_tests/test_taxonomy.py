import unittest
import subprocess
from subprocess import PIPE
from parameterized import parameterized

from tempfile import TemporaryDirectory
import os


"""Test taxonomy classification framework"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein") # TODO - decide if we need to consider this "_Protein" case
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
        self.lca_coverage = 0.9
        self.k = 20
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

    def build_graph_and_annotation(self, num_threads: int):
        construct_command = '{exe} build -p {num_threads} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=num_threads,
            k=self.k,
            outfile=self.tempdir.name + '/graph',
            input=TAX_DATA_DIR + '/tax_input.fa'
        )
        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        annotate_command = '{exe} annotate --anno-header -i {dbg} -o {anno} -p {num_threads} {input_fasta}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            anno=self.tempdir.name + '/annotation',
            num_threads=num_threads,
            input_fasta=TAX_DATA_DIR + '/tax_input.fa'
        )
        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

    def get_prediction_statistics_from(self, res_lines: [str]) -> {}:
        result = {}
        result["num_correct_tips"] = 0
        result["num_correct_internals"] = 0
        result["num_total_tips"] = 0
        result["num_total_internals"] = 0

        result["num_descendant_prediction_internals"] = 0
        result["num_ancestor_prediction_tips"] = 0
        result["num_ancestor_prediction_internals"] = 0

        result["num_wrong_prediction_tips"] = 0
        result["num_wrong_prediction_internals"] = 0

        result["num_too_few_discovered_kmers"] = 0

        for line in res_lines:
            if line == "":
                continue
            query_expected = line.split(" ")[1].split("|")[1].strip()
            query_prediction = line.split(" ")[7].split("'")[1].strip()

            # TaxId 0 means that there were not enough discovered kmers in order to realize the tax classification.
            if query_prediction == "0":
                result["num_too_few_discovered_kmers"] += 1
                continue

            # All the tax nodes with ids {10001, 10002 .. 10008} represents internal nodes. TaxIds >= 10009 are reserved for the leaves.
            if int(line.split(" ")[1].split("|")[1]) >= 10009:
                # The current taxid is a tip, thus, it has no children in the taxonomic tree.
                result["num_total_tips"] += 1
                if query_expected == query_prediction:
                    result["num_correct_tips"] += 1
                else:
                    if self.is_descendant(node=query_prediction, query=query_expected):
                        result["num_ancestor_prediction_tips"] += 1
                    else:
                        result["num_wrong_prediction_tips"] += 1
            else:
                # The current taxid is an internal node.
                result["num_total_internals"] += 1
                if query_expected == query_prediction:
                    result["num_correct_internals"] += 1
                else:
                    if self.is_descendant(node=query_prediction, query=query_expected):
                        result["num_ancestor_prediction_internals"] += 1
                    elif self.is_descendant(node=query_expected, query=query_prediction):
                        result["num_descendant_prediction_internals"] += 1
                    else:
                        result["num_wrong_prediction_internals"] += 1
        return result

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets") # TODO - decide if this skipIf can be deleted.
    def test_taxonomy_taxdb(self, tax_test):
        self.build_graph_and_annotation(tax_tests[tax_test]['threads'])
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

        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-db {taxDB} \
                            --lca-coverage-fraction {lca_coverage} -p {num_threads}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            taxDB=self.tempdir.name + '/taxDB.taxdb',
            lca_coverage=self.lca_coverage,
            num_threads=tax_tests[tax_test]['threads'],
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res_lines = res.stdout.decode().rstrip().split('\n')
        statistics = self.get_prediction_statistics_from(res_lines)

        self.assertEqual(statistics["num_total_tips"], 118)
        self.assertEqual(statistics["num_total_internals"], 80)

        self.assertEqual(statistics["num_correct_tips"], 109)
        self.assertEqual(statistics["num_correct_internals"], 38)

        self.assertEqual(statistics["num_ancestor_prediction_internals"], 5)
        self.assertEqual(statistics["num_descendant_prediction_internals"], 34)
        self.assertEqual(statistics["num_ancestor_prediction_tips"], 9)

        self.assertEqual(statistics["num_wrong_prediction_internals"], 3)
        self.assertEqual(statistics["num_wrong_prediction_tips"], 0)

        self.assertEqual(statistics["num_too_few_discovered_kmers"], 2)

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy_getrows(self, tax_test):
        self.build_graph_and_annotation(tax_tests[tax_test]['threads'])
        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {tax_tree} \
                            --lca-coverage-fraction {lca_coverage} -p {num_threads} -a {anno}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            tax_tree=TAX_DATA_DIR + '/dumb_nodes.dmp',
            lca_coverage=self.lca_coverage,
            num_threads=tax_tests[tax_test]['threads'],
            anno=self.tempdir.name + '/annotation.column.annodbg',
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res_lines = res.stdout.decode().rstrip().split('\n')
        statistics = self.get_prediction_statistics_from(res_lines)

        self.assertEqual(statistics["num_total_tips"], 118)
        self.assertEqual(statistics["num_total_internals"], 80)

        self.assertEqual(statistics["num_correct_tips"], 109)
        self.assertEqual(statistics["num_correct_internals"], 38)

        self.assertEqual(statistics["num_ancestor_prediction_internals"], 5)
        self.assertEqual(statistics["num_descendant_prediction_internals"], 34)
        self.assertEqual(statistics["num_ancestor_prediction_tips"], 9)

        self.assertEqual(statistics["num_wrong_prediction_internals"], 3)
        self.assertEqual(statistics["num_wrong_prediction_tips"], 0)

        self.assertEqual(statistics["num_too_few_discovered_kmers"], 2)

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy_toplabels(self, tax_test):
        self.build_graph_and_annotation(tax_tests[tax_test]['threads'])
        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {tax_tree} \
                            --lca-coverage-fraction {lca_coverage} -p {num_threads} -a {anno} \
                            --top-label-fraction {top_label_fraction}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            tax_tree=TAX_DATA_DIR + '/dumb_nodes.dmp',
            lca_coverage=self.lca_coverage,
            num_threads=tax_tests[tax_test]['threads'],
            anno=self.tempdir.name + '/annotation.column.annodbg',
            top_label_fraction=0.7,
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res_lines = res.stdout.decode().rstrip().split('\n')
        statistics = self.get_prediction_statistics_from(res_lines)

        self.assertEqual(statistics["num_total_tips"], 118)
        self.assertEqual(statistics["num_total_internals"], 68)

        self.assertEqual(statistics["num_correct_tips"], 74)
        self.assertEqual(statistics["num_correct_internals"], 24)

        self.assertEqual(statistics["num_ancestor_prediction_internals"], 27)
        self.assertEqual(statistics["num_descendant_prediction_internals"], 15)
        self.assertEqual(statistics["num_ancestor_prediction_tips"], 44)

        self.assertEqual(statistics["num_wrong_prediction_internals"], 2)
        self.assertEqual(statistics["num_wrong_prediction_tips"], 0)

        self.assertEqual(statistics["num_too_few_discovered_kmers"], 14)
