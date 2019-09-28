import unittest
import subprocess
from tempfile import TemporaryDirectory
import glob
import os


"""Test graph construction"""


TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../data'


class TestBuild(unittest.TestCase):
    def test_simple_all_graphs(self):
        """
        Simple build test
        """

        tempdir = TemporaryDirectory()

        for representation in ['succinct', 'bitmap', 'hash', 'hashstr']:

            command = './metagraph build --graph {repr} -k 20 -o {outfile} {input}'.format(
                repr=representation,
                outfile=tempdir.name + '/graph',
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )

            res = subprocess.run([command], shell=True)
            self.assertEqual(res.returncode, 0)

    def test_build_tiny_k(self):
        """
        Simple build test
        """

        tempdir = TemporaryDirectory()

        args = ['./metagraph', 'build',
                    '-k', '2',
                    '-o', tempdir.name + '/graph',
                    TEST_DATA_DIR + '/transcripts_1000.fa']
        command = ' '.join(args)

        res = subprocess.run([command], shell=True)
        self.assertEqual(res.returncode, 0)


if __name__ == '__main__':
    unittest.main()
