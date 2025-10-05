import unittest
from parameterized import parameterized, parameterized_class
import subprocess
import itertools
from subprocess import PIPE
from tempfile import TemporaryDirectory
import glob
import os
import numpy as np
from helpers import get_test_class_name
from base import PROTEIN_MODE, DNA_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension
import hashlib
import platform
import re
import psutil


"""Test graph construction"""

anno_file_extension = {'column': '.column.annodbg',
                       'column_coord': '.column_coord.annodbg',
                       'brwt_coord': '.brwt_coord.annodbg',
                       'row_diff_coord': '.row_diff_coord.annodbg',
                       'row_diff_brwt_coord': '.row_diff_brwt_coord.annodbg',
                       'row_diff_disk_coord': '.row_diff_disk_coord.annodbg',
                       'row': '.row.annodbg',
                       'row_diff': '.row_diff.annodbg',
                       'row_sparse': '.row_sparse.annodbg',
                       'row_diff_brwt': '.row_diff_brwt.annodbg',
                       'row_diff_disk': '.row_diff_disk.annodbg',
                       'row_diff_flat': '.row_diff_flat.annodbg',
                       'row_diff_sparse': '.row_diff_sparse.annodbg',
                       'row_diff_sparse_noswap': '.row_diff_sparse.annodbg',
                       'rb_brwt': '.rb_brwt.annodbg',
                       'brwt': '.brwt.annodbg',
                       'int_brwt': '.int_brwt.annodbg',
                       'row_diff_int_brwt': '.row_diff_int_brwt.annodbg',
                       'row_diff_int_disk': '.row_diff_int_disk.annodbg',
                       'rbfish': '.rbfish.annodbg',
                       'flat': '.flat.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]
ANNO_TYPES = [anno_type for anno_type, _ in anno_file_extension.items()]

NUM_THREADS = 4

MEMORY_MAPPING = True
MMAP_FLAG = ' --mmap' if MEMORY_MAPPING else ''


def product(graph_types, anno_types):
    result  = []
    for graph in graph_types:
        for anno in anno_types:
            if graph == 'succinct' or not anno.startswith('row_diff'):
                result.append((graph, anno))
    return result


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ) + product(['succinct_bloom', 'succinct_mask'], ['flat']),
    class_name_func=get_test_class_name
)
class TestQuery(TestingBase):
    def _run_with_memory_monitor(self, query_command, expected_stdout_len):
        """Run command with /usr/bin/time and report memory usage"""
        # Get system memory before
        mem_before = psutil.virtual_memory()
        sys_before_gb = mem_before.used / 1024 / 1024 / 1024
        sys_total_gb = mem_before.total / 1024 / 1024 / 1024
        sys_before_pct = mem_before.percent
        print(f"\n\033[33m[MEM BEFORE]\033[0m {self.graph_repr}+{self.anno_repr}: System before: {sys_before_gb:.1f}/{sys_total_gb:.1f}GB ({sys_before_pct:.1f}%)", flush=True)
        
        # Use /usr/bin/time to measure memory
        if platform.system() == 'Darwin':  # macOS
            time_cmd = '/usr/bin/time -l'
        else:  # Linux
            time_cmd = '/usr/bin/time -v'
        
        full_command = f"{time_cmd} {query_command}"
        res = subprocess.run(full_command, shell=True, stdout=PIPE, stderr=PIPE)
        
        # Get system memory after
        mem_after = psutil.virtual_memory()
        sys_after_gb = mem_after.used / 1024 / 1024 / 1024
        sys_after_pct = mem_after.percent
        
        # Extract memory info from stderr
        stderr_text = res.stderr.decode()
        if 'maximum resident set size' in stderr_text.lower():
            # Linux format: "	Maximum resident set size (kbytes): 3076"
            # macOS format: "  1234567  maximum resident set size"
            match = re.search(r'maximum resident set size \(kbytes\):\s*(\d+)', stderr_text, re.IGNORECASE)
            if match:
                mem_kb = int(match.group(1))
                mem_mb = mem_kb / 1024
            else:
                match = re.search(r'\s+(\d+)\s+maximum resident set size', stderr_text)
                if match:
                    mem_bytes = int(match.group(1))
                    mem_mb = mem_bytes / 1024 / 1024
                else:
                    mem_mb = 0
            
            if mem_mb > 0:
                print(f"\033[33m[ MEM AFTER]\033[0m {self.graph_repr}+{self.anno_repr}: Peak={mem_mb:.1f}MB | System after: {sys_after_gb:.1f}/{sys_total_gb:.1f}GB ({sys_after_pct:.1f}%)", flush=True)
        
        if res.returncode != 0:
            print(f"\n[ERROR] Command failed with return code {res.returncode}")
            print(f"STDERR: {stderr_text[:500]}")
        
        self.assertEqual(res.returncode, 0)
        if expected_stdout_len is not None:
            self.assertEqual(len(res.stdout), expected_stdout_len)
        
        return res

    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic',
                         '--mask-dummy' if cls.mask_dummy else '')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert('20' == stats_graph['k'])
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert('46960' == stats_graph['nodes (k)'])
        assert('basic' == stats_graph['mode'])

        if cls.with_bloom:
            convert_command = f'{METAGRAPH} transform -o {cls.tempdir.name}/graph \
                                --initialize-bloom --bloom-fpp 0.1 \
                                {cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}' + MMAP_FLAG
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt
        )

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert('100' == stats_annotation['labels'])
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert(stats_graph['max index (k)'] == stats_annotation['objects'])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(cls.anno_repr == stats_annotation['representation'])

    def test_query(self):
        query_command = '{exe} query --batch-size 0 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 137140)

        query_command = '{exe} query --batch-size 0 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_both(self):
        """query graph (fwd and reverse)"""
        query_command = '{exe} query --batch-size 0 --fwd-and-reverse -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 261390)

        query_command = '{exe} query --batch-size 0 --fwd-and-reverse --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 260215)

    def test_query_parallel(self):
        """query graph (multi-threaded)"""
        query_command = '{exe} query --batch-size 0 -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 137140)

        query_command = '{exe} query --batch-size 0 --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_both_parallel(self):
        """query graph (fwd and reverse, multi-threaded)"""
        query_command = '{exe} query --batch-size 0 --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 261390)

        query_command = '{exe} query --batch-size 0 --fwd-and-reverse --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 260215)

    def test_query_with_align(self):
        query_command = '{exe} query --batch-size 0 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        
        expected_len = 12249 if DNA_MODE else 12244
        self._run_with_memory_monitor(query_command, expected_len)

        query_command = '{exe} query --batch-size 0 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        
        expected_len = 12355 if DNA_MODE else 12350
        self._run_with_memory_monitor(query_command, expected_len)

        # align to graph (multi-threaded)
        query_command = '{exe} query --batch-size 0 --align -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        expected_len = 12249 if DNA_MODE else 12244
        self._run_with_memory_monitor(query_command, expected_len)

        query_command = '{exe} query --batch-size 0 --align --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        expected_len = 12355 if DNA_MODE else 12350
        self._run_with_memory_monitor(query_command, expected_len)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_query_with_align_both(self):
        """align to graph (fwd and reverse multi-threaded)"""
        query_command = '{exe} query --batch-size 0 --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 24567)

        query_command = '{exe} query --batch-size 0 --fwd-and-reverse --align --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 24779)

    def test_batch_query(self):
        query_command = '{exe} query --batch-size 100000000 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 137140)

        query_command = '{exe} query --batch-size 100000000 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_both(self):
        """query graph (fwd and reverse)"""
        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 261390)

        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 260215)

    def test_batch_query_parallel(self):
        """query graph (multi-threaded)"""
        query_command = '{exe} query --batch-size 100000000 -i {graph} -a {annotation} -p {num_threads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_threads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 137140)

        query_command = '{exe} query --batch-size 100000000 --query-mode matches -i {graph} -a {annotation} -p {num_threads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_threads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 136959)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_both_parallel(self):
        """query graph (fwd and reverse, multi-threaded)"""
        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG

        self._run_with_memory_monitor(query_command, 261390)

        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 260215)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --batch-size 100000000 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG

        expected_len = 12249 if DNA_MODE else 12244
        self._run_with_memory_monitor(query_command, expected_len)

        query_command = '{exe} query --batch-size 100000000 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG

        expected_len = 12355 if DNA_MODE else 12350
        self._run_with_memory_monitor(query_command, expected_len)

        # align to graph (multi-threaded)
        query_command = '{exe} query --batch-size 100000000 --align -i {graph} -a {annotation} -p {num_threads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        ) + MMAP_FLAG

        expected_len = 12249 if DNA_MODE else 12244
        self._run_with_memory_monitor(query_command, expected_len)

        query_command = '{exe} query --batch-size 100000000 --align --query-mode matches -i {graph} -a {annotation} -p {num_threads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_threads=NUM_THREADS
        ) + MMAP_FLAG

        expected_len = 12355 if DNA_MODE else 12350
        self._run_with_memory_monitor(query_command, expected_len)

    @unittest.skipIf(PROTEIN_MODE, "Reverse sequences for Protein alphabets are not defined")
    def test_batch_query_with_align_both(self):
        """align to graph (fwd and reverse multi-threaded)"""
        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse --align -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 24567)

        query_command = '{exe} query --batch-size 100000000 --fwd-and-reverse --align --query-mode matches -i {graph} -a {annotation} -p {num_theads} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa',
            num_theads=NUM_THREADS
        ) + MMAP_FLAG
        
        self._run_with_memory_monitor(query_command, 24779)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --batch-size 100000000 --batch-size 100 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137140)

        query_command = '{exe} query --batch-size 100000000 --batch-size 100 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 136959)

    def test_query_coordinates(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode coords \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 0.05 {TEST_DATA_DIR}/transcripts_100.fa' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 139268)

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode coords \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 0.95 {TEST_DATA_DIR}/transcripts_100.fa' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 31522)

    def test_query_coordinates_expanded(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode coords --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 0.05 {TEST_DATA_DIR}/transcripts_100.fa' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 1619883)

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode coords --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 0.95 {TEST_DATA_DIR}/transcripts_100.fa' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 492788)


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(['succinct'], ANNO_TYPES),
    class_name_func=get_test_class_name
)
class TestQueryTinyLinear(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.fasta_graph = cls.tempdir.name + '/file_graph.fa'
        with open(cls.fasta_graph, 'w') as f:
            f.write(f'>L\nAAAACCCCGGGGTTTT\n')

        cls.fasta_anno = cls.tempdir.name + '/file_anno.fa'
        with open(cls.fasta_anno, 'w') as f:
            f.write(f'>L1\nAAAACCCCGGGGTTTT\n')
            f.write(f'>L2\nAAAACCCCGGGGTTTT\n')
            f.write(f'>L3\nCCCCGGGG\n')

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(cls.fasta_graph, cls.tempdir.name + '/graph', 5, cls.graph_repr, 'basic', '--mask-dummy')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '5')
        assert(stats_graph['nodes (k)'] == '12')
        assert(stats_graph['mode'] == 'basic')

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(cls.fasta_anno,
                cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
                cls.tempdir.name + '/annotation', cls.anno_repr,
                separate, no_fork_opt, no_anchor_opt)

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert(stats_annotation['labels'] == '3')
        assert(stats_annotation['objects'] == stats_graph['max index (k)'])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(cls.anno_repr == stats_annotation['representation'])

    def test_query_coordinates(self):
        if not self.anno_repr.endswith('_coord'):
            self.skipTest('annotation does not support coordinates')

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode coords  --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 0.05 {self.fasta_graph}' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(set(res.stdout.decode().strip().split('\t')),
                         {'0', 'L',
                            '<L1>:0:1:2:3:4:5:6:7:8:9:10:11',
                            '<L2>:0:1:2:3:4:5:6:7:8:9:10:11',
                            '<L3>:::::0:1:2:3::::'})


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        ANNO_TYPES + ['row_diff_brwt_separate',
                      'row_diff_brwt_no_fork_opt',
                      'row_diff_brwt_no_anchor_opt']
    ) + product(['succinct_bloom', 'succinct_mask'], ['flat']),
    class_name_func=get_test_class_name
)
class TestQuery1Column(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic',
                         '--mask-dummy' if cls.mask_dummy else '')

        stats_graph = cls._get_stats(cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr])
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '20')
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert(stats_graph['nodes (k)'] == '46960')
        assert(stats_graph['mode'] == 'basic')

        if cls.with_bloom:
            convert_command = f'{METAGRAPH} transform -o {cls.tempdir.name}/graph \
                                --initialize-bloom --bloom-fpp 0.1 \
                                {cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}' + MMAP_FLAG
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        def check_suffix(anno_repr, suffix):
            match = anno_repr.endswith(suffix)
            if match:
                anno_repr = anno_repr[:-len(suffix)]
            return anno_repr, match

        cls.anno_repr, separate = check_suffix(cls.anno_repr, '_separate')
        cls.anno_repr, no_fork_opt = check_suffix(cls.anno_repr, '_no_fork_opt')
        cls.anno_repr, no_anchor_opt = check_suffix(cls.anno_repr, '_no_anchor_opt')

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            separate,
            no_fork_opt,
            no_anchor_opt,
            anno_type='label 1'
        )

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert(stats_annotation['labels'] == '1')
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert(stats_annotation['objects'] == stats_graph['max index (k)'])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(cls.anno_repr == stats_annotation['representation'])

    def test_query(self):
        query_command = f'{METAGRAPH} query --batch-size 0 \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 1.0 \
                            {TEST_DATA_DIR}/transcripts_1000.fa' + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(hashlib.sha224(res.stdout).hexdigest(), '254d173abb255a81a4ab8a685201a73de8dbad4546c378e0a645d454')

        query_command = f'{METAGRAPH} query --batch-size 0 --query-mode matches \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label 1.0 \
                            {TEST_DATA_DIR}/transcripts_1000.fa' + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(hashlib.sha224(res.stdout).hexdigest(), '1bd6c24373812064c3e17e73533de7b1e30baa3cca3a64b460e83cb4')


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=product(
        [repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)],
        [anno_type for anno_type in ANNO_TYPES if '_int_' in anno_type or anno_type.endswith('_coord')]
    ),
    class_name_func=get_test_class_name
)
class TestQueryCounts(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.kmer_counts_1 = {
            'AAA': 1,
            'AAC': 2,
            'ACC': 3,
            'CCC': 4,
            'CCG': 5,
            'CGG': 6,
            'GGG': 7,
            'GGT': 8,
            'GTT': 9,
            'TTT': 10,
            'TTA': 11,
            'TAA': 12,
        }
        cls.kmer_counts_2 = {
            'AAA': 11,
            'AAC': 12,
            'ACC': 13,
            'CCC': 14,
            'CCG': 15,
            'CGG': 16,
            'GGG': 17,
            'GGT': 18,
            'GTT': 19,
            'TTT': 20,
        }

        cls.fasta_file_1 = cls.tempdir.name + '/file_1.fa'
        with open(cls.fasta_file_1, 'w') as f:
            for kmer, count in cls.kmer_counts_1.items():
                f.write(f'>L1\n{kmer}\n' * count)

        cls.fasta_file_2 = cls.tempdir.name + '/file_2.fa'
        with open(cls.fasta_file_2, 'w') as f:
            for kmer, count in cls.kmer_counts_2.items():
                f.write(f'>L2\n{kmer}\n' * count)

        cls.k = 3

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph((cls.fasta_file_1, cls.fasta_file_2), cls.tempdir.name + '/graph',
                         cls.k, cls.graph_repr, 'basic', '--mask-dummy' if cls.mask_dummy else '')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '3')
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert(stats_graph['nodes (k)'] == '12')
        assert(stats_graph['mode'] == 'basic')

        if cls.with_bloom:
            convert_command = f'{METAGRAPH} transform -o {cls.tempdir.name}/graph \
                                --initialize-bloom --bloom-fpp 0.1 \
                                {cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}' + MMAP_FLAG
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        cls._annotate_graph(
            (cls.fasta_file_1, cls.fasta_file_2),
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr,
            anno_type='filename'
        )

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert(stats_annotation['labels'] == '2')
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert(stats_annotation['objects'] == stats_graph['max index (k)'])
        assert(stats_annotation['representation'] == cls.anno_repr)

        cls.queries = [
            'AAA',
            'AAAA',
            'AAAAAAAAAAAAA',
            'CCC',
            'CCCC',
            'CCCCCCCCCCCCC',
            'TTT',
            'AAACCCGGGTTT',
            'AAACCCGGGTTTTTT',
            'AAACCCGGGTTTAAA',
            'TTTAAACCCGGG',
            'ACACACACACACATTTAAACCCGGG',
        ]

    def _compare_unsorted_results(self, output, expected):
        self.assertEqual(len(output), len(expected))
        output_lines = output.split('\n')
        expected_lines = expected.split('\n')
        self.assertEqual(len(output_lines), len(expected_lines))
        for output_line, expected_line in zip(output_lines, expected_lines):
            self.assertCountEqual(output_line.split('\t'), expected_line.split('\t'))

    def test_abundance_sum_query(self):
        query_file = self.tempdir.name + '/query.fa'
        for discovery_rate in np.linspace(0, 1, 5):
            expected_output = ''
            with open(query_file, 'w') as f:
                for i, s in enumerate(self.queries):
                    f.write(f'>s{i}\n{s}\n')
                    expected_output += f'{i}\ts{i}'
                    def get_count(d, kmer):
                        try:
                            return d[kmer]
                        except:
                            return 0

                    num_kmers = len(s) - self.k + 1

                    num_matches_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    count_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) for i in range(len(s) - self.k + 1)])

                    num_matches_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    count_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) for i in range(len(s) - self.k + 1)])

                    for (c, label, n) in [(count_1, self.fasta_file_1, num_matches_1),
                                            (count_2, self.fasta_file_2, num_matches_2)]:
                        if n >= discovery_rate * num_kmers:
                            expected_output += f'\t<{label}>:{c}'

                    expected_output += '\n'

            query_command = f'{METAGRAPH} query --batch-size 100000000 --query-mode counts-sum \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label {discovery_rate} {query_file}' + MMAP_FLAG

            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self._compare_unsorted_results(res.stdout.decode(), expected_output)

    def test_count_query(self):
        query_file = self.tempdir.name + '/query.fa'
        for discovery_rate in np.linspace(0, 1, 5):
            expected_output = ''
            with open(query_file, 'w') as f:
                for i, s in enumerate(self.queries):
                    f.write(f'>s{i}\n{s}\n')
                    expected_output += f'{i}\ts{i}'
                    def get_count(d, kmer):
                        try:
                            return d[kmer]
                        except:
                            return 0

                    num_kmers = len(s) - self.k + 1

                    num_matches_1 = sum([get_count(self.kmer_counts_1, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    counts_1 = [str(get_count(self.kmer_counts_1, s[i:i + self.k])) for i in range(len(s) - self.k + 1)]

                    num_matches_2 = sum([get_count(self.kmer_counts_2, s[i:i + self.k]) > 0 for i in range(num_kmers)])
                    counts_2 = [str(get_count(self.kmer_counts_2, s[i:i + self.k])) for i in range(len(s) - self.k + 1)]

                    for (cs, label, n) in [(counts_1, self.fasta_file_1, num_matches_1),
                                            (counts_2, self.fasta_file_2, num_matches_2)]:
                        if n >= discovery_rate * num_kmers:
                            cs = ':'.join(cs)
                            expected_output += f'\t<{label}>:{cs}'

                    expected_output += '\n'

            query_command = f'{METAGRAPH} query --batch-size 100000000 --query-mode counts --verbose-output \
                            -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                            -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                            --min-kmers-fraction-label {discovery_rate} {query_file}' + MMAP_FLAG

            res = subprocess.run(query_command.split(), stdout=PIPE)
            self.assertEqual(res.returncode, 0)
            self._compare_unsorted_results(res.stdout.decode(), expected_output)

        query_command = f'{METAGRAPH} query --batch-size 100000000 --query-mode counts \
                        -i {self.tempdir.name}/graph{graph_file_extension[self.graph_repr]} \
                        -a {self.tempdir.name}/annotation{anno_file_extension[self.anno_repr]} \
                        --min-kmers-fraction-label {discovery_rate} {query_file}' + MMAP_FLAG

        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self._compare_unsorted_results(res.stdout.decode(),
            f"0\ts0\t<{self.fasta_file_1}>:0=1\t<{self.fasta_file_2}>:0=11\n"
            f"1\ts1\t<{self.fasta_file_1}>:0-1=1\t<{self.fasta_file_2}>:0-1=11\n"
            f"2\ts2\t<{self.fasta_file_1}>:0-10=1\t<{self.fasta_file_2}>:0-10=11\n"
            f"3\ts3\t<{self.fasta_file_1}>:0=4\t<{self.fasta_file_2}>:0=14\n"
            f"4\ts4\t<{self.fasta_file_1}>:0-1=4\t<{self.fasta_file_2}>:0-1=14\n"
            f"5\ts5\t<{self.fasta_file_1}>:0-10=4\t<{self.fasta_file_2}>:0-10=14\n"
            f"6\ts6\t<{self.fasta_file_1}>:0=10\t<{self.fasta_file_2}>:0=20\n"
            f"7\ts7\t<{self.fasta_file_1}>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9=10\t<{self.fasta_file_2}>:0=11:1=12:2=13:3=14:4=15:5=16:6=17:7=18:8=19:9=20\n"
            f"8\ts8\t<{self.fasta_file_1}>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9-12=10\t<{self.fasta_file_2}>:0=11:1=12:2=13:3=14:4=15:5=16:6=17:7=18:8=19:9-12=20\n"
            f"9\ts9\t<{self.fasta_file_1}>:0=1:1=2:2=3:3=4:4=5:5=6:6=7:7=8:8=9:9=10:10=11:11=12:12=1\n"
            f"10\ts10\t<{self.fasta_file_1}>:0=10:1=11:2=12:3=1:4=2:5=3:6=4:7=5:8=6:9=7\n"
            "11\ts11\n")


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=(product(list(set(GRAPH_TYPES) - {'hashstr'}), ANNO_TYPES) +
                  product(['succinct_bloom', 'succinct_mask'], ['flat'])),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryCanonical(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'canonical',
                         '--mask-dummy' if cls.mask_dummy else '')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '20')
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert(stats_graph['nodes (k)'] == '91584')
        assert(stats_graph['mode'] == 'canonical')

        if cls.with_bloom:
            convert_command = f'{METAGRAPH} transform -o {cls.tempdir.name}/graph \
                                --initialize-bloom --bloom-fpp 0.1 \
                                {cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}' + MMAP_FLAG
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr
        )

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert(stats_annotation['labels'] == '100')
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert(stats_annotation['objects'] == stats_graph['max index (k)'])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(cls.anno_repr == stats_annotation['representation'])

    def test_query(self):
        query_command = '{exe} query --batch-size 0 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 0 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_query_with_align(self):
        query_command = '{exe} query --batch-size 0 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12840)

        query_command = '{exe} query --batch-size 0 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12970)

    def test_batch_query(self):
        query_command = '{exe} query --batch-size 100000000 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 100000000 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --batch-size 100000000 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12840)

        query_command = '{exe} query --batch-size 100000000 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12970)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --batch-size 100000000 --batch-size 100 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 100000000 --batch-size 100 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)


@parameterized_class(('graph_repr', 'anno_repr'),
    input_values=(product(list(set(GRAPH_TYPES) - {'hashstr'}), ANNO_TYPES) +
                  product(['succinct_bloom', 'succinct_mask'], ['flat'])),
    class_name_func=get_test_class_name
)
@unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
class TestQueryPrimary(TestingBase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()

        cls.with_bloom = False
        if cls.graph_repr == 'succinct_bloom':
            cls.graph_repr = 'succinct'
            cls.with_bloom = True

        cls.mask_dummy = False
        if cls.graph_repr == 'succinct_mask':
            cls.graph_repr = 'succinct'
            cls.mask_dummy = True

        cls._build_graph(TEST_DATA_DIR + '/transcripts_100.fa',
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'primary',
                         '--mask-dummy' if cls.mask_dummy else '')

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '20')
        if cls.graph_repr != 'succinct' or cls.mask_dummy:
            assert(stats_graph['nodes (k)'] == '45792')
        assert(stats_graph['mode'] == 'primary')

        if cls.with_bloom:
            convert_command = f'{METAGRAPH} transform -o {cls.tempdir.name}/graph \
                                --initialize-bloom --bloom-fpp 0.1 \
                                {cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}' + MMAP_FLAG
            res = subprocess.run([convert_command], shell=True)
            assert(res.returncode == 0)

        cls._annotate_graph(
            TEST_DATA_DIR + '/transcripts_100.fa',
            cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
            cls.tempdir.name + '/annotation',
            cls.anno_repr
        )

        # check annotation
        stats_annotation = cls._get_stats(f'-a {cls.tempdir.name}/annotation{anno_file_extension[cls.anno_repr]}')
        assert(stats_annotation['returncode'] == 0)
        assert(stats_annotation['labels'] == '100')
        if cls.graph_repr != 'hashfast' and (cls.graph_repr != 'succinct' or cls.mask_dummy):
            assert(stats_annotation['objects'] == stats_graph['max index (k)'])

        if cls.anno_repr.endswith('_noswap'):
            cls.anno_repr = cls.anno_repr[:-len('_noswap')]

        assert(cls.anno_repr == stats_annotation['representation'])

    def test_query(self):
        query_command = '{exe} query --batch-size 0 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 0 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_query_with_align(self):
        query_command = '{exe} query --batch-size 0 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12840)

        query_command = '{exe} query --batch-size 0 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12970)

    def test_batch_query(self):
        query_command = '{exe} query --batch-size 100000000 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 100000000 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)

    def test_batch_query_with_align(self):
        query_command = '{exe} query --batch-size 100000000 --align -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12840)

        query_command = '{exe} query --batch-size 100000000 --align --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 0.0 --align-min-exact-match 0.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_100_tail10_snp.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 12970)

    def test_batch_query_with_tiny_batch(self):
        query_command = '{exe} query --batch-size 100000000 --batch-size 100 -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137269)

        query_command = '{exe} query --batch-size 100000000 --batch-size 100 --query-mode matches -i {graph} -a {annotation} --min-kmers-fraction-label 1.0 {input}'.format(
            exe=METAGRAPH,
            graph=self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr],
            annotation=self.tempdir.name + '/annotation' + anno_file_extension[self.anno_repr],
            input=TEST_DATA_DIR + '/transcripts_1000.fa'
        ) + MMAP_FLAG
        res = subprocess.run(query_command.split(), stdout=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertEqual(len(res.stdout), 137093)


if __name__ == '__main__':
    unittest.main()
