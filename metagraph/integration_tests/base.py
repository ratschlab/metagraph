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
                        'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg',
                       'row_diff': '.row_diff.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]

NUM_THREADS = 4


class TestingBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

    @staticmethod
    def _get_stats(graph_path):
        stats_command = METAGRAPH + ' stats ' + graph_path
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        return res

    @staticmethod
    def _build_graph(input, output, k, repr, mode='basic', extra_params=''):
        if not output.endswith(graph_file_extension[repr]):
            output += graph_file_extension[repr]

        assert mode in ['basic', 'primary', 'canonical']

        construct_command = '{exe} build -p {num_threads} --mode {mode} {extra_params} \
                --graph {repr} -k {k} -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            extra_params=extra_params,
            k=k,
            repr=repr,
            mode='basic' if mode == 'basic' else 'canonical',
            outfile=output,
            input=input
        )

        res = subprocess.run([construct_command], shell=True, stdout=PIPE,
                             stderr=PIPE)
        assert res.returncode == 0

        if mode == 'primary':
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

            construct_command = '{exe} build --mode primary -p {num_threads} {extra_params} \
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
    def _annotate_graph(input, graph_path, output, anno_repr,
                        separate=False, no_fork_opt=False, no_anchor_opt=False):
        target_anno = anno_repr
        if (anno_repr in {'row_sparse', 'column_coord'} or
                anno_repr.endswith('brwt') or
                anno_repr.startswith('row_diff')):
            target_anno = anno_repr
            anno_repr = 'column'
        elif anno_repr in {'flat', 'rbfish'}:
            target_anno = anno_repr
            anno_repr = 'row'

        command = f'{METAGRAPH} annotate -p {NUM_THREADS} --anno-header \
                    -i {graph_path} --anno-type {anno_repr} \
                    -o {output} {input}'

        if target_anno.endswith('_coord'):
            command += ' --coordinates'

        with_counts = target_anno.endswith('int_brwt')
        if with_counts:
            command += ' --count-kmers'

        res = subprocess.run([command], shell=True)
        assert(res.returncode == 0)

        if target_anno == anno_repr:
            return

        final_anno = target_anno
        if final_anno.startswith('row_diff'):
            target_anno = 'row_diff'

        command = f'{METAGRAPH} transform_anno -p {NUM_THREADS} \
                    --anno-type {target_anno} -o {output} \
                    {output + anno_file_extension[anno_repr]}'

        other_args = ' --count-kmers' if with_counts else ''
        other_args += ' --coordinates' if final_anno.endswith('_coord') else ''

        if target_anno == 'row_diff':
            command += ' -i ' + graph_path

        if not no_fork_opt:
            if target_anno.startswith('row_diff'):
                print('-- Building RowDiff without fork optimization...')
            res = subprocess.run([command + other_args], shell=True)
            assert(res.returncode == 0)

        if target_anno == 'row_diff':
            without_input_anno = command.split(' ')
            without_input_anno.pop(-3)
            without_input_anno = ' '.join(without_input_anno)
            if not no_anchor_opt:
                if separate:
                    print('-- Building RowDiff succ/pred...')
                    res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 1'], shell=True)
                    assert(res.returncode == 0)
                res = subprocess.run([command + ' --row-diff-stage 1' + other_args], shell=True)
                assert(res.returncode == 0)
            if separate:
                print('-- Assigning anchors...')
                res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 2'], shell=True)
                assert(res.returncode == 0)
            res = subprocess.run([command + ' --row-diff-stage 2' + other_args], shell=True)
            assert(res.returncode == 0)

            if final_anno != target_anno:
                rd_type = 'column' if with_counts or final_anno.endswith('_coord') else 'row_diff'
                command = f'{METAGRAPH} transform_anno --anno-type {final_anno} --greedy -o {output} ' \
                                   f'-i {graph_path} -p {NUM_THREADS} {output}.{rd_type}.annodbg'
                res = subprocess.run([command], shell=True)
                assert (res.returncode == 0)
                os.remove(output + anno_file_extension[rd_type])
            else:
                os.remove(output + anno_file_extension[anno_repr])
