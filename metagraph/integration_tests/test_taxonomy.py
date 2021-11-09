import unittest
import subprocess
from subprocess import PIPE
from parameterized import parameterized

from tempfile import TemporaryDirectory
import os


"""Test taxonomy classification framework"""

METAGRAPH = './metagraph'
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein") # TODO - decide if we need to consider this "_Protein" case
TAX_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data/taxonomic_data'

tax_tests = {
    'one_thread': {
        'threads': 1,
    },
    'nine_threads': {
        'threads': 9,
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

    def is_descendant(self, target: str, curr: str) -> bool:
        if curr == target:
            return True
        while curr != self.tax_root:
            curr = self.tax_parent[curr]
            if curr == target:
                return True
        return False

    def build_graph_and_anno_matrix(self, num_threads: int):
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

    def get_prediction_statistics_from_stdout(self, stdout_lines: [str]) -> {}:
        result = {"num_tip_hit": 0,
                  "num_internal_hit": 0,
                  "total_num_tips": 0,
                  "total_num_internals": 0,
                  "num_descendant_hit": 0,
                  "num_ancestor_hit_for_tips": 0,
                  "num_ancestor_hit_for_internals": 0,
                  "num_tip_misses": 0,
                  "num_internal_misses": 0,
                  "num_failed_classification": 0}

        for line in stdout_lines:
            if line == "":
                continue
            query_expected = line.split(" ")[1].split("|")[1].strip()
            query_prediction = line.split(" ")[7].split("'")[1].strip()

            # TaxId 0 is a wildcard for not enough discovered kmers to produce a confident classification.
            if query_prediction == "0":
                result["num_failed_classification"] += 1
                continue

            # All the tax nodes with ids {10001, 10002 .. 10008} represents internal nodes, while
            # taxIds >= 10009 are reserved for the leaves.
            if int(line.split(" ")[1].split("|")[1]) >= 10009:
                # The current taxid is a tip, thus, it has no children in the taxonomic tree.
                result["total_num_tips"] += 1
                if query_expected == query_prediction:
                    result["num_tip_hit"] += 1
                else:
                    if self.is_descendant(target=query_prediction, curr=query_expected):
                        result["num_ancestor_hit_for_tips"] += 1
                    else:
                        result["num_tip_misses"] += 1
            else:
                # The current taxid is an internal node.
                result["total_num_internals"] += 1
                if query_expected == query_prediction:
                    result["num_internal_hit"] += 1
                else:
                    if self.is_descendant(target=query_prediction, curr=query_expected):
                        result["num_ancestor_hit_for_internals"] += 1
                    elif self.is_descendant(target=query_expected, curr=query_prediction):
                        result["num_descendant_hit"] += 1
                    else:
                        result["num_internal_misses"] += 1
        return result

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy_getrows(self, tax_test):
        self.build_graph_and_anno_matrix(tax_tests[tax_test]['threads'])
        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {tax_tree} \
                            --min-lca-coverage {lca_coverage} --label-taxid-map {label_taxid_map} ' \
                            '-p {num_threads} -a {anno}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            tax_tree=TAX_DATA_DIR + '/dumb_nodes.dmp',
            lca_coverage=self.lca_coverage,
            label_taxid_map=TAX_DATA_DIR + '/dumb.accession2taxid',
            num_threads=tax_tests[tax_test]['threads'],
            anno=self.tempdir.name + '/annotation.column.annodbg',
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res_lines = res.stdout.decode().rstrip().split('\n')
        statistics = self.get_prediction_statistics_from_stdout(res_lines)

        self.assertEqual(statistics["total_num_tips"], 118)
        self.assertEqual(statistics["total_num_internals"], 80)

        self.assertEqual(statistics["num_tip_hit"], 109)
        self.assertEqual(statistics["num_internal_hit"], 38)

        self.assertEqual(statistics["num_ancestor_hit_for_internals"], 5)
        self.assertEqual(statistics["num_descendant_hit"], 34)
        self.assertEqual(statistics["num_ancestor_hit_for_tips"], 9)

        self.assertEqual(statistics["num_internal_misses"], 3)
        self.assertEqual(statistics["num_tip_misses"], 0)

        self.assertEqual(statistics["num_failed_classification"], 2)

    @parameterized.expand(test_params)
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_taxonomy_toplabels(self, tax_test):
        self.build_graph_and_anno_matrix(tax_tests[tax_test]['threads'])
        tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {tax_tree} \
                            --min-lca-coverage {lca_coverage} -p {num_threads} -a {anno} \
                            --label-taxid-map {label_taxid_map} \
                            --top-label-fraction {top_label_fraction}'.format(
            exe=METAGRAPH,
            dbg=self.tempdir.name + '/graph.dbg',
            fasta_queries=TAX_DATA_DIR + '/tax_query.fa',
            tax_tree=TAX_DATA_DIR + '/dumb_nodes.dmp',
            lca_coverage=self.lca_coverage,
            label_taxid_map=TAX_DATA_DIR + '/dumb.accession2taxid',
            num_threads=tax_tests[tax_test]['threads'],
            anno=self.tempdir.name + '/annotation.column.annodbg',
            top_label_fraction=0.7,
        )
        res = subprocess.run([tax_class_command], shell=True, stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res_lines = res.stdout.decode().rstrip().split('\n')
        statistics = self.get_prediction_statistics_from_stdout(res_lines)

        self.assertEqual(statistics["total_num_tips"], 118)
        self.assertEqual(statistics["total_num_internals"], 68)

        self.assertEqual(statistics["num_tip_hit"], 74)
        self.assertEqual(statistics["num_internal_hit"], 24)

        self.assertEqual(statistics["num_ancestor_hit_for_internals"], 27)
        self.assertEqual(statistics["num_descendant_hit"], 15)
        self.assertEqual(statistics["num_ancestor_hit_for_tips"], 44)

        self.assertEqual(statistics["num_internal_misses"], 2)
        self.assertEqual(statistics["num_tip_misses"], 0)

        self.assertEqual(statistics["num_failed_classification"], 14)
