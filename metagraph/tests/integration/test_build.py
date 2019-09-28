import unittest
import subprocess
from tempfile import NamedTemporaryFile
import glob
import os


"""Test graph construction"""


TEST_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../data'


def remove_files(prefix):
    for out_file in glob.glob(prefix + '*'):
        os.remove(out_file)


class TestBuild(unittest.TestCase):
    def setUp(self):
        self.__temp_filename = NamedTemporaryFile(delete=False)
        self.addCleanup(remove_files, self.__temp_filename.name)

    def test_simple_all_graphs(self):
        """
        Simple build test
        """

        for representation in ['succinct', 'bitmap', 'hash', 'hashstr']:

            command = './metagraph build --graph {repr} -k 20 -o {outfile} {input}'.format(
                repr=representation,
                outfile=self.__temp_filename.name,
                input=TEST_DATA_DIR + '/transcripts_1000.fa'
            )

            res = subprocess.run([command], shell=True)
            self.assertEqual(res.returncode, 0)

    def test_build_tiny_k(self):
        """
        Simple build test
        """

        args = ['./metagraph', 'build',
                    '-k', '2',
                    '-o', self.__temp_filename.name,
                    TEST_DATA_DIR + '/transcripts_1000.fa']
        command = ' '.join(args)
        print(command)
        res = subprocess.run([command], shell=True)
        self.assertEqual(res.returncode, 0)


if __name__ == '__main__':
    unittest.main()
