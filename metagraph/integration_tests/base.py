import os
import subprocess
import unittest
from subprocess import PIPE
from tempfile import TemporaryDirectory
import psutil
import shutil

script_path = os.path.dirname(os.path.realpath(__file__))

METAGRAPH_EXE = f'{os.getcwd()}/metagraph'
_target = os.readlink(METAGRAPH_EXE)
DNA_MODE = "DNA" in _target
PROTEIN_MODE = "Protein" in _target
METAGRAPH = METAGRAPH_EXE

def update_prefix(PREFIX):
    global METAGRAPH
    METAGRAPH = PREFIX + METAGRAPH_EXE

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

MEMORY_MAPPING = True
MMAP_FLAG = ' --mmap' if MEMORY_MAPPING else ''


class TestingBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

    @staticmethod
    def _get_stats(graph_path):
        stats_command = METAGRAPH + ' stats ' + graph_path + ' --mmap'
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        if res.returncode != 0:
            raise AssertionError(f"Command '{stats_command}' failed with return code {res.returncode} and error: {res.stderr.decode()}")
        stats_command = METAGRAPH + ' stats ' + graph_path + MMAP_FLAG
        res = subprocess.run(stats_command.split(), stdout=PIPE, stderr=PIPE)
        parsed = dict()
        parsed['returncode'] = res.returncode
        res = res.stdout.decode().split('\n')[2:]
        for line in res:
            if ': ' in line:
                x, y = map(str.strip, line.split(':', 1))
                assert(x not in parsed or parsed[x] == y)
                parsed[x] = y
        return parsed

    @staticmethod
    def _build_graph(input, output, k, repr, mode='basic', extra_params=''):
        if not isinstance(input, str):
            input = ' '.join(input)

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
        ) + MMAP_FLAG

        res = subprocess.run([construct_command], shell=True, stdout=PIPE, stderr=PIPE)
        if res.returncode != 0:
            raise AssertionError(f"Build command failed with return code {res.returncode}\nCommand: {construct_command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

        if mode == 'primary':
            transform_command = '{exe} transform -p {num_threads} --to-fasta --primary-kmers \
                    -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                k=k,
                repr=repr,
                outfile='{}.fasta.gz'.format(output),
                input=output
            ) + MMAP_FLAG

            res = subprocess.run([transform_command], shell=True, stdout=PIPE, stderr=PIPE)
            if res.returncode != 0:
                raise AssertionError(f"Transform command failed with return code {res.returncode}\nCommand: {transform_command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

            construct_command = '{exe} build --mode primary -p {num_threads} {extra_params} \
                    --graph {repr} -k {k} -o {outfile} {input}'.format(
                exe=METAGRAPH,
                num_threads=NUM_THREADS,
                extra_params=extra_params,
                k=k,
                repr=repr,
                outfile=output,
                input='{}.fasta.gz'.format(output)
            ) + MMAP_FLAG

            res = subprocess.run([construct_command], shell=True, stdout=PIPE, stderr=PIPE)
            if res.returncode != 0:
                raise AssertionError(f"Build primary command failed with return code {res.returncode}\nCommand: {construct_command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

    @staticmethod
    def _clean(graph, output, extra_params=''):
        clean_command = '{exe} clean -p {num_threads} \
                --to-fasta -o {outfile} {extra_params} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            outfile=output,
            extra_params=extra_params,
            input=graph
        ) + MMAP_FLAG
        res = subprocess.run([clean_command], shell=True, stdout=PIPE, stderr=PIPE)
        if res.returncode != 0:
            raise AssertionError(f"Clean command failed with return code {res.returncode}\nCommand: {clean_command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

    @staticmethod
    def _annotate_graph(input, graph_path, output, anno_repr,
                        separate=False, no_fork_opt=False, no_anchor_opt=False,
                        anno_type='header', extra_params='', num_threads=NUM_THREADS):
        if not isinstance(input, str):
            input = ' '.join(input)

        target_anno = anno_repr

        noswap = anno_repr.endswith('_noswap')
        if noswap:
            anno_repr = anno_repr[:-len('_noswap')]

        if (anno_repr in {'row_sparse', 'column_coord'} or
                anno_repr.endswith('_coord') or
                anno_repr.endswith('brwt') or
                anno_repr.startswith('row_diff')):
            target_anno = anno_repr
            anno_repr = 'column'
        elif anno_repr in {'flat', 'rbfish'}:
            target_anno = anno_repr
            anno_repr = 'row'

        command = f'{METAGRAPH} annotate -v -p {num_threads} --anno-{anno_type}\
                    -i {graph_path} --anno-type {anno_repr} {extra_params} \
                    -o {output} {input}'

        if target_anno.endswith('_coord'):
            command += ' --coordinates'

        with_counts = 'int_' in target_anno
        if with_counts:
            command += ' --count-kmers'

        # Monitor RAM and disk usage before command

        mem_before = psutil.virtual_memory()
        disk_before = shutil.disk_usage('/')
        sys_before_gb = (mem_before.total - mem_before.available) / 1024 / 1024 / 1024
        sys_total_gb = mem_before.total / 1024 / 1024 / 1024
        sys_before_pct = mem_before.percent
        disk_before_gb = disk_before.used / 1024 / 1024 / 1024
        disk_total_gb = disk_before.total / 1024 / 1024 / 1024
        disk_before_pct = (disk_before.used / disk_before.total) * 100

        print(f"\n\033[33m[RESOURCE BEFORE]\033[0m Annotation: RAM {sys_before_gb:.1f}/{sys_total_gb:.1f}GB ({sys_before_pct:.1f}%) | Disk {disk_before_gb:.1f}/{disk_total_gb:.1f}GB ({disk_before_pct:.1f}%)", flush=True)

        res = subprocess.run([command], shell=True, stdout=PIPE, stderr=PIPE)

        # Monitor RAM and disk usage after command
        mem_after = psutil.virtual_memory()
        disk_after = shutil.disk_usage('/')
        sys_after_gb = (mem_after.total - mem_after.available) / 1024 / 1024 / 1024
        sys_after_pct = mem_after.percent
        disk_after_gb = disk_after.used / 1024 / 1024 / 1024
        disk_after_pct = (disk_after.used / disk_after.total) * 100

        print(f"\033[33m[RESOURCE AFTER ]\033[0m Annotation: RAM {sys_after_gb:.1f}/{sys_total_gb:.1f}GB ({sys_after_pct:.1f}%) | Disk {disk_after_gb:.1f}/{disk_total_gb:.1f}GB ({disk_after_pct:.1f}%)", flush=True)
        if res.returncode != 0:
            raise AssertionError(f"Annotate command failed with return code {res.returncode}\nCommand: {command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

        if target_anno == anno_repr:
            return

        final_anno = target_anno
        if final_anno.startswith('row_diff'):
            target_anno = 'row_diff'

        command = f'{METAGRAPH} transform_anno {MMAP_FLAG} -p {num_threads} \
                    --anno-type {target_anno} -o {output} \
                    {output + anno_file_extension[anno_repr]}'

        other_args = ' --count-kmers' if with_counts else ''
        other_args += ' --coordinates' if final_anno.endswith('_coord') else ''
        other_args += ' --disk-swap \"\"' if noswap else ''

        if target_anno == 'row_diff':
            command += ' -i ' + graph_path

        if not no_fork_opt:
            if target_anno.startswith('row_diff'):
                print('-- Building RowDiff without fork optimization...')
            res = subprocess.run([command + other_args], shell=True, stdout=PIPE, stderr=PIPE)
            if res.returncode != 0:
                raise AssertionError(f"Transform_anno command failed with return code {res.returncode}\nCommand: {command + other_args}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")

        if target_anno == 'row_diff':
            without_input_anno = command.split(' ')
            without_input_anno.pop(-3)
            without_input_anno = ' '.join(without_input_anno)
            if not no_anchor_opt:
                if separate:
                    print('-- Building RowDiff succ/pred...')
                    res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 1'], shell=True, stdout=PIPE, stderr=PIPE)
                    if res.returncode != 0:
                        raise AssertionError(f"RowDiff stage 1 (separate) command failed with return code {res.returncode}\nCommand: echo \"\" | {without_input_anno} --row-diff-stage 1\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
                res = subprocess.run([command + ' --row-diff-stage 1' + other_args], shell=True, stdout=PIPE, stderr=PIPE)
                if res.returncode != 0:
                    raise AssertionError(f"RowDiff stage 1 command failed with return code {res.returncode}\nCommand: {command} --row-diff-stage 1{other_args}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
                subprocess.run([f'rm -f {output}.row_count'], shell=True)
            if separate:
                print('-- Assigning anchors...')
                res = subprocess.run(['echo \"\" | ' + without_input_anno + ' --row-diff-stage 2'], shell=True, stdout=PIPE, stderr=PIPE)
                if res.returncode != 0:
                    raise AssertionError(f"RowDiff stage 2 (separate) command failed with return code {res.returncode}\nCommand: echo \"\" | {without_input_anno} --row-diff-stage 2\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
            res = subprocess.run([command + ' --row-diff-stage 2' + other_args], shell=True, stdout=PIPE, stderr=PIPE)
            if res.returncode != 0:
                raise AssertionError(f"RowDiff stage 2 command failed with return code {res.returncode}\nCommand: {command} --row-diff-stage 2{other_args}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
            subprocess.run([f'rm -f {output}.row_reduction'], shell=True)

            if final_anno != target_anno:
                rd_type = 'column' if with_counts or final_anno.endswith('_coord') else 'row_diff'
                command = f'{METAGRAPH} transform_anno --anno-type {final_anno} --greedy -o {output} ' \
                                   f'-i {graph_path} -p {num_threads} {output}.{rd_type}.annodbg' + MMAP_FLAG
                res = subprocess.run([command], shell=True, stdout=PIPE, stderr=PIPE)
                if res.returncode != 0:
                    raise AssertionError(f"Transform_anno (final) command failed with return code {res.returncode}\nCommand: {command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
                subprocess.run([f'rm {output}{anno_file_extension[rd_type]}*'], shell=True)
            else:
                subprocess.run([f'rm {output}{anno_file_extension[anno_repr]}*'], shell=True)

        if final_anno.endswith('brwt') or final_anno.endswith('brwt_coord'):
            command = f'{METAGRAPH} relax_brwt -o {output} -p {num_threads} {output}.{final_anno}.annodbg' + MMAP_FLAG
            res = subprocess.run([command], shell=True, stdout=PIPE, stderr=PIPE)
            if res.returncode != 0:
                raise AssertionError(f"Relax_brwt command failed with return code {res.returncode}\nCommand: {command}\nStdout: {res.stdout.decode()}\nStderr: {res.stderr.decode()}")
