import os
import subprocess
import unittest
from subprocess import PIPE
from tempfile import TemporaryDirectory

script_path = os.path.dirname(os.path.realpath(__file__))

METAGRAPH = './metagraph'

TEST_DATA_DIR = os.path.join(script_path, '..', 'tests', 'data')

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        # 'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]

NUM_THREADS = 4


class TestingBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

    @staticmethod
    def _get_stats(graph_filename):
        stats_command = METAGRAPH + ' stats ' + graph_filename
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        return res

    @staticmethod
    def _build_graph(input, output, k, repr, canonical=False, primary=False, extra_params=''):
        construct_command = '{exe} build -p {num_threads} {canonical} {extra_params} \
                --graph {repr} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            extra_params=extra_params,
            k=k,
            repr=repr,
            canonical='--canonical' if canonical else '',
            outfile=output,
            input=input
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE,
                             stderr=PIPE)
        assert res.returncode == 0

        if primary:
            transform_command = '{exe} transform -p {num_threads} --to-fasta --primary-kmers \
                    -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                k=k,
                repr=repr,
                outfile='{}.fasta.gz'.format(output),
                input=output
            )

            res = subprocess.run([transform_command], shell=True, stdout=PIPE,
                                 stderr=PIPE)
            assert res.returncode == 0

            construct_command = '{exe} build -p {num_threads} {extra_params} \
                    --graph {repr} -k {k} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                extra_params=extra_params,
                k=k,
                repr=repr,
                outfile=output,
                input='{}.fasta.gz'.format(output)
            )

            res = subprocess.run([construct_command], shell=True, stdout=PIPE,
                                 stderr=PIPE)
            assert res.returncode == 0

    @staticmethod
    def _clean(graph, output, extra_params=''):
        clean_command = '{exe} clean -p {num_threads} \
                --to-fasta -o {outfile} {extra_params} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=output,
            extra_params=extra_params,
            input=graph
        )
        res = subprocess.run([clean_command], shell=True)
        assert res.returncode == 0

    @staticmethod
    def _annotate_graph(input, graph_path, output, anno_repr, primary=False):
        annotate_command = '{exe} annotate {fwd_and_rev} --anno-header -i {graph} \
                --anno-type {anno_repr} -o {outfile} -p {num_threads} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            fwd_and_rev='--canonical' if primary else '',
            graph=graph_path,
            anno_repr=anno_repr,
            outfile=output,
            input=input
        )
        res = subprocess.run([annotate_command], shell=True)
        assert res.returncode == 0
