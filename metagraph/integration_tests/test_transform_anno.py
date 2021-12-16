import unittest
from parameterized import parameterized, parameterized_class
import itertools
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import gzip
from base import TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension, anno_file_extension
from helpers import get_test_class_name

"""Test operations on annotation columns"""

DNA_MODE = os.readlink(METAGRAPH).endswith("_DNA")
PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")
TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../tests/data'

NUM_THREADS = 4


class TestColumnOperations(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        cls.graph_repr = 'succinct'
        cls.anno_repr = 'column'

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic', '--mask-dummy')

        res = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(res.returncode == 0)
        out = res.stdout.decode().split('\n')[2:]
        assert('k: 20' == out[0])
        assert('nodes (k): 46960' == out[1])
        assert('mode: basic' == out[2])

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            extra_params='--count-kmers'
        )

    def setUp(self):
        # use relative paths to get rid of dir names
        self.old_cwd = os.getcwd()
        os.chdir(self.tempdir.name)
        self.annotation = f'annotation{anno_file_extension[self.anno_repr]}';

        # check annotation
        res = self._get_stats(f'-a {self.annotation}')
        self.assertEqual(res.returncode, 0)
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('labels:  100', out[0])
        self.assertEqual('objects: 46960', out[1])
        self.assertEqual('density: 0.0185072', out[2])
        self.assertEqual(f'representation: {self.anno_repr}', out[3])

    def tearDown(self):
        os.chdir(self.old_cwd)

    def test_overlap(self):
        command = f'{METAGRAPH} transform_anno {self.annotation} \
            --compute-overlap {self.annotation} -o out -p {NUM_THREADS}'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(156421, len(res.stdout.decode()))

    def _check_aggregation_min(self, min_count, expected_density):
        command = f'{METAGRAPH} transform_anno {self.annotation} -p {NUM_THREADS} \
            --aggregate-columns --min-count {min_count} -o aggregated'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(res.returncode, 0)
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('labels:  1', out[0])
        self.assertEqual('objects: 46960', out[1])
        self.assertEqual(f'density: {expected_density}', out[2])
        self.assertEqual(f'representation: {self.anno_repr}', out[3])

    def test_aggregate_columns(self):
        self._check_aggregation_min(0, 1)
        self._check_aggregation_min(1, 1)
        self._check_aggregation_min(5, 0.0715077)
        self._check_aggregation_min(10, 0.00344974)
        self._check_aggregation_min(20, 0)

    def _check_aggregation_min_max_value(self, min_count, max_value, expected_density):
        command = f'{METAGRAPH} transform_anno {self.annotation} -p {NUM_THREADS} \
            --aggregate-columns --min-count {min_count} --max-value {max_value} -o aggregated'

        res = subprocess.run(command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)

        res = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(res.returncode, 0)
        out = res.stdout.decode().split('\n')[2:]
        self.assertEqual('labels:  1', out[0])
        self.assertEqual('objects: 46960', out[1])
        self.assertEqual(f'density: {expected_density}', out[2])
        self.assertEqual(f'representation: {self.anno_repr}', out[3])

    def test_aggregate_columns_filtered(self):
        self._check_aggregation_min_max_value(0, 0, 0)
        self._check_aggregation_min_max_value(1, 0, 0)
        self._check_aggregation_min_max_value(2, 0, 0)
        self._check_aggregation_min_max_value(3, 0, 0)
        self._check_aggregation_min_max_value(5, 0, 0)

        self._check_aggregation_min_max_value(0, 1, 0.99704)
        self._check_aggregation_min_max_value(1, 1, 0.99704)
        self._check_aggregation_min_max_value(2, 1, 0.392994)
        self._check_aggregation_min_max_value(3, 1, 0.183305)
        self._check_aggregation_min_max_value(5, 1, 0.0715077)

        self._check_aggregation_min_max_value(0, 2, 0.998807)
        self._check_aggregation_min_max_value(1, 2, 0.998807)
        self._check_aggregation_min_max_value(2, 2, 0.394825)
        self._check_aggregation_min_max_value(3, 2, 0.183986)
        self._check_aggregation_min_max_value(5, 2, 0.0715077)

        self._check_aggregation_min_max_value(0, 5, 0.998999)
        self._check_aggregation_min_max_value(1, 5, 0.998999)
        self._check_aggregation_min_max_value(2, 5, 0.395315)
        self._check_aggregation_min_max_value(3, 5, 0.184817)
        self._check_aggregation_min_max_value(5, 5, 0.0715077)

        self._check_aggregation_min_max_value(0, 1000, 1)
        self._check_aggregation_min_max_value(1, 1000, 1)
        self._check_aggregation_min_max_value(2, 1000, 0.395336)
        self._check_aggregation_min_max_value(3, 1000, 0.184817)
        self._check_aggregation_min_max_value(5, 1000, 0.0715077)


if __name__ == '__main__':
    unittest.main()
