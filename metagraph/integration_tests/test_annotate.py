import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import filecmp
import glob
import os
from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, NUM_THREADS


"""Test graph annotation"""

graph_file_extension = {'succinct': '.dbg',
                        'bitmap': '.bitmapdbg',
                        'hash': '.orhashdbg',
                        # 'hashfast': '.hashfastdbg',
                        'hashstr': '.hashstrdbg'}

anno_file_extension = {'column': '.column.annodbg',
                       'row': '.row.annodbg'}

GRAPH_TYPES = [graph_type for graph_type, _ in graph_file_extension.items()]


class TestAnnotate(TestingBase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    @parameterized.expand([repr for repr in GRAPH_TYPES if not (repr == 'bitmap' and PROTEIN_MODE)])
    def test_simple_all_graphs(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --in-ram -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('46960', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} \
                                -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_100.fa'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('100', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.0185072 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_canonical(self, graph_repr):

        construct_command = '{exe} build --mask-dummy --in-ram -p {num_threads} \
                --graph {repr} --mode canonical -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('91584', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_100.fa'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('100', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.00948888 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_all_graphs_from_kmc(self, graph_repr):
        """
        Annotate non-canonical graph constructed from non-canonical KMC database
        """

        construct_command = f'{METAGRAPH} build --mask-dummy --in-ram -p {NUM_THREADS} \
                            --graph {graph_repr} -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('469983', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                1 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_simple_all_graphs_from_kmc_both(self, graph_repr):
        """
        Annotate non-canonical graph constructed from canonical KMC database
        """

        construct_command = f'{METAGRAPH} build --mask-dummy --in-ram -p {NUM_THREADS} \
                            --graph {graph_repr} -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_single \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.585342 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

            # both strands
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_both \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                1 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    # TODO: add 'hashstr' once the canonical mode is implemented for it
    @parameterized.expand(['succinct', 'bitmap', 'hash'])  # , 'hashstr']:
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_simple_all_graphs_from_kmc_both_canonical(self, graph_repr):
        """
        Annotate canonical graph with k-mers from KMC
        """

        construct_command = f'{METAGRAPH} build --mask-dummy --in-ram -p {NUM_THREADS} \
                            --graph {graph_repr} --mode canonical -k 11 \
                            -o {self.tempdir.name}/graph \
                            {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('11', stats_graph['k'])
        self.assertEqual('802920', stats_graph['nodes (k)'])
        self.assertEqual('canonical', stats_graph['mode'])

        for anno_repr in ['row', 'column']:
            # build annotation
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_single \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_single' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.5 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

            # both strands
            annotate_command = f'{METAGRAPH} annotate --anno-label LabelName -p {NUM_THREADS} \
                                -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                                --anno-type {anno_repr} -o {self.tempdir.name}/annotation_both \
                                {TEST_DATA_DIR}/transcripts_1000_kmc_counters_both_strands.kmc_suf'

            res = subprocess.run([annotate_command], shell=True)
            self.assertEqual(res.returncode, 0)

            # check annotation
            stats_annotation = self._get_stats('-a ' + self.tempdir.name + '/annotation_both' + anno_file_extension[anno_repr])
            self.assertEqual(stats_annotation['returncode'], 0)
            self.assertEqual('1', stats_annotation['labels'])
            self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
            self.assertAlmostEqual(
                0.5 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
                float(stats_annotation['density']),
                places=6)
            self.assertEqual(anno_repr, stats_annotation['representation'])

    def test_annotate_with_disk_swap(self):
        graph_repr = 'succinct'
        anno_repr = 'column'

        construct_command = '{exe} build --mask-dummy --in-ram -p {num_threads} \
                --graph {repr} -k 20 -o {outfile} {input}'.format(
            exe=METAGRAPH,
            num_threads=NUM_THREADS,
            repr=graph_repr,
            outfile=self.tempdir.name + '/graph',
            input=TEST_DATA_DIR + '/transcripts_100.fa'
        )

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        stats_graph = self._get_stats(self.tempdir.name + '/graph' + graph_file_extension[graph_repr])
        self.assertEqual(stats_graph['returncode'], 0)
        self.assertEqual('20', stats_graph['k'])
        self.assertEqual('46960', stats_graph['nodes (k)'])
        self.assertEqual('basic', stats_graph['mode'])

        # build annotation
        annotate_command = f'{METAGRAPH} annotate --anno-header \
                            --disk-swap {self.tempdir.name} --mem-cap-gb 1e-6 \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            --anno-type {anno_repr} -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # check annotation
        stats_annotation = self._get_stats('-a ' + f'{self.tempdir.name}/annotation{anno_file_extension[anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual('100', stats_annotation['labels'])
        self.assertEqual(stats_graph['max index (k)'], stats_annotation['objects'])
        self.assertAlmostEqual(
            0.0185072 * (int(stats_graph['nodes (k)']) / int(stats_graph['max index (k)'])),
            float(stats_annotation['density']),
            places=6)
        self.assertEqual(anno_repr, stats_annotation['representation'])

    @parameterized.expand(GRAPH_TYPES)
    def test_annotate_coordinates(self, graph_repr):

        construct_command = f'{METAGRAPH} build --mask-dummy --in-ram -p {NUM_THREADS} \
                              --graph {graph_repr} -k 11 \
                              -o {self.tempdir.name}/graph \
                              {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

    @parameterized.expand(GRAPH_TYPES)
    def test_annotate_coordinates_with_disk_swap(self, graph_repr):

        construct_command = f'{METAGRAPH} build --mask-dummy --in-ram -p {NUM_THREADS} \
                              --graph {graph_repr} -k 11 \
                              -o {self.tempdir.name}/graph \
                              {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([construct_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation with disk swap
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            --disk-swap {self.tempdir.name} --mem-cap-gb 1e-6 \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        # build annotation in RAM
        annotate_command = f'{METAGRAPH} annotate --anno-header -p {NUM_THREADS} --coordinates \
                            -i {self.tempdir.name}/graph{graph_file_extension[graph_repr]} \
                            -o {self.tempdir.name}/annotation_ram \
                            {TEST_DATA_DIR}/transcripts_100.fa'

        res = subprocess.run([annotate_command], shell=True)
        self.assertEqual(res.returncode, 0)

        self.assertTrue(filecmp.cmp(f'{self.tempdir.name}/annotation.column.annodbg',
                                    f'{self.tempdir.name}/annotation_ram.column.annodbg'))
        self.assertTrue(filecmp.cmp(f'{self.tempdir.name}/annotation.column.annodbg.coords',
                                    f'{self.tempdir.name}/annotation_ram.column.annodbg.coords'))


    def test_annotate_coordinates_hints_index_header_coords(self):
        """`metagraph annotate --coordinates --anno-filename` (without
        `--index-header-coords`) must point users at the second pass that
        builds the CoordToHeader (.seqs) mapping. This is the brand-new
        hint at the end of annotate_graph().

        The hint is emitted via `logger->info()`; the split_sink in
        common/logger.cpp routes info-level messages to stdout (warn+
        goes to stderr), so we look in stdout.
        """
        construct_command = (f'{METAGRAPH} build --mask-dummy --in-ram '
                             f'-p {NUM_THREADS} --graph succinct -k 11 '
                             f'-o {self.tempdir.name}/graph '
                             f'{TEST_DATA_DIR}/transcripts_100.fa')
        res = subprocess.run([construct_command], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0)

        annotate_cmd = (f'{METAGRAPH} annotate --anno-filename --coordinates '
                        f'-p {NUM_THREADS} '
                        f'-i {self.tempdir.name}/graph.dbg '
                        f'-o {self.tempdir.name}/annotation '
                        f'{TEST_DATA_DIR}/transcripts_100.fa')
        res = subprocess.run([annotate_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0)
        stdout = res.stdout.decode()
        self.assertIn('To enable per-sequence coordinate reporting', stdout)
        self.assertIn('--index-header-coords', stdout)
        self.assertIn(f'-i {self.tempdir.name}/graph.dbg', stdout)
        self.assertIn(f'-o {self.tempdir.name}/annotation', stdout)

        # When --index-header-coords is passed, the hint must not appear
        # (the user is already doing what the hint suggests).
        index_cmd = (f'{METAGRAPH} annotate --anno-filename '
                     f'--index-header-coords -p {NUM_THREADS} '
                     f'-i {self.tempdir.name}/graph.dbg '
                     f'-o {self.tempdir.name}/annotation '
                     f'{TEST_DATA_DIR}/transcripts_100.fa')
        res = subprocess.run([index_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0)
        self.assertNotIn('To enable per-sequence coordinate reporting',
                         res.stdout.decode())


class TestAnnotateHeaderChecks(TestingBase):
    """`metagraph annotate` refuses a header-derived name that is not valid UTF-8 (a
    --anno-header label, a --index-header-coords record header) and, with
    --index-header-coords, a record header repeating an earlier record's of the same file,
    naming the file, the 1-based record(s) and the bytes in hex (never raw), with nothing
    written; --allow-unchecked-headers skips both checks with a warning (the owner's decision
    of 2026-10-11). A record's reverse complement under --fwd-and-reverse is the same record,
    not a repeat."""
    K = 11
    SEQ_A = 'ACGTTGCAGGCTAGCTTACGGATCCAGTACGATGCCTTAGC'
    SEQ_B = 'TTGACCGTAGGCATCGATTCAAGGCTTACGACTGGATACCA'
    SEQ_C = 'GGATCCTTAAGCGCGTATACCAGTTGCAATGGCTCGATCAG'

    def setUp(self):
        self.tempdir = TemporaryDirectory()
        self.d = self.tempdir.name

    def _fasta(self, records, name='records.fa'):
        # headers are bytes: some are not UTF-8 on purpose
        path = os.path.join(self.d, name)
        with open(path, 'wb') as f:
            for header, seq in records:
                f.write(b'>' + header + b'\n' + seq.encode() + b'\n')
        return path

    def _graph(self, fasta):
        self._build_graph(fasta, os.path.join(self.d, 'graph'), self.K, 'succinct')
        return os.path.join(self.d, 'graph.dbg')

    def _annotate(self, graph, fasta, flags, out='anno'):
        cmd = (f'{METAGRAPH} annotate -p {NUM_THREADS} -i {graph} {flags} '
               f'-o {os.path.join(self.d, out)} {fasta}')
        return subprocess.run([cmd], shell=True, stdout=PIPE, stderr=PIPE)

    def _assert_refused(self, res, fasta, *fragments):
        self.assertNotEqual(res.returncode, 0)
        err = res.stderr.decode('utf-8', errors='replace')
        self.assertIn(fasta, err)
        for fragment in fragments:
            self.assertIn(fragment, err)
        # the offending bytes are named in hex only
        self.assertNotIn(b'\xff', res.stderr)
        self.assertEqual([], [f for f in os.listdir(self.d)
                              if f.startswith('anno')], 'nothing is written')

    def test_anno_header_refuses_a_label_that_is_not_utf8(self):
        fasta = self._fasta([(b'ok', self.SEQ_A), (b'bad\xff', self.SEQ_B),
                             (b'later', self.SEQ_C)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-header')
        self._assert_refused(res, fasta, 'record 2', 'not valid UTF-8', '626164ff',
                             '--allow-unchecked-headers')
        # the label is the name token: a comment that is not UTF-8 is dropped with the
        # comment, and refused once --header-comment-delim makes it part of the label
        fasta = self._fasta([(b'ok', self.SEQ_A), (b'name meta\xff', self.SEQ_B)], 'comment.fa')
        res = self._annotate(graph, fasta, '--anno-header', out='anno_name')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertTrue(os.path.exists(os.path.join(self.d, 'anno_name.column.annodbg')))
        res = self._annotate(graph, fasta, '--anno-header --header-comment-delim _', out='x')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('record 2', res.stderr.decode('utf-8', 'replace'))
        self.assertIn('6e616d655f6d657461ff', res.stderr.decode('utf-8', 'replace'))
        self.assertFalse(os.path.exists(os.path.join(self.d, 'x.column.annodbg')))

    def test_index_header_coords_refuses_a_record_header_that_is_not_utf8(self):
        fasta = self._fasta([(b'ok', self.SEQ_A), (b'bad\xff', self.SEQ_B)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords')
        self._assert_refused(res, fasta, 'record 2', 'not valid UTF-8', '626164ff')

    def test_index_header_coords_refuses_a_record_without_a_header_name(self):
        """A header of `>` alone (or a comment alone, without --header-comment-delim) yields
        no label: the record is refused by name, since a record of the mapping needs one."""
        fasta = self._fasta([(b'ok', self.SEQ_A), (b'', self.SEQ_B)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords')
        self._assert_refused(res, fasta, 'record 2', 'is empty')
        fasta = self._fasta([(b'ok', self.SEQ_A), (b' only-a-comment', self.SEQ_B)], 'comment.fa')
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords')
        self._assert_refused(res, fasta, 'record 2', 'is empty')
        # joined with the comment, the header is a name
        res = self._annotate(graph, fasta,
                             '--anno-filename --index-header-coords --header-comment-delim _',
                             out='joined')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertIn('sequences: 2 (ok\tonly-a-comment)\n', self._seqs_stats('joined'))

    def _seqs_stats(self, out):
        res = subprocess.run([f'{METAGRAPH} stats -v {os.path.join(self.d, out)}.seqs'],
                             shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        return res.stdout.decode()

    def _column_names(self, out):
        res = subprocess.run([f'{METAGRAPH} stats --print-col-names '
                              f'-a {os.path.join(self.d, out)}.column.annodbg'],
                             shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        lines = res.stdout.decode().splitlines()
        # the logger's info lines go to stdout too, before the count
        start = next(i for i, line in enumerate(lines) if line.startswith('Number of columns: '))
        names = lines[start + 1:]
        self.assertEqual('Number of columns: %d' % len(names), lines[start])
        return names

    def test_a_record_without_a_comment_is_named_without_one(self):
        """kseq keeps the comment buffer of the last record that had one: under
        --header-comment-delim a record without a comment is named by its name token alone,
        for a --anno-header label and a --index-header-coords record header, so `>x meta`
        followed by `>x` is two names, not a repeat."""
        fasta = self._fasta([(b'r1 meta', self.SEQ_A), (b'r2', self.SEQ_B),
                             (b'r3 other', self.SEQ_C)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-header --header-comment-delim _', out='labels')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertEqual(['r1_meta', 'r2', 'r3_other'], self._column_names('labels'))
        res = self._annotate(graph, fasta,
                             '--anno-filename --index-header-coords --header-comment-delim _',
                             out='records')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertIn('sequences: 3 (r1_meta\tr2\tr3_other)\n', self._seqs_stats('records'))
        fasta = self._fasta([(b'x meta', self.SEQ_A), (b'x', self.SEQ_B)], 'joined.fa')
        res = self._annotate(graph, fasta,
                             '--anno-filename --index-header-coords --header-comment-delim _',
                             out='joined')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertIn('sequences: 2 (x_meta\tx)\n', self._seqs_stats('joined'))
        # the same name token twice is a repeat, whatever the comments
        fasta = self._fasta([(b'x meta', self.SEQ_A), (b'x', self.SEQ_B)], 'repeat.fa')
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords', out='repeat')
        self._assert_refused(res, fasta, 'record 2', 'repeats the header of record 1')

    def test_count_kmers_from_headers_refuses_a_record_without_a_count(self):
        """--count-kmers without a .kmer_counts.gz file reads ka:f:N from each record's
        comment: a record without a comment after one with a count is inconsistent, not the
        previous record's count (kseq's comment buffer outlives its record)."""
        fasta = self._fasta([(b'r1 ka:f:3', self.SEQ_A), (b'r2', self.SEQ_B)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-header --count-kmers', out='counts')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('Inconsistent k-mer count headers',
                      res.stderr.decode('utf-8', 'replace'))
        fasta = self._fasta([(b'r1 ka:f:3', self.SEQ_A), (b'r2 ka:f:2', self.SEQ_B)],
                            'counted.fa')
        res = self._annotate(graph, fasta, '--anno-header --count-kmers', out='counted')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertEqual(['r1', 'r2'], self._column_names('counted'))

    def test_separately_writes_nothing_when_a_later_file_is_refused(self):
        """--separately annotates each file on its own and writes its annotation as soon as it
        is done: the headers of every file are checked before the first annotation starts, so
        a refused file leaves no annotation of another file behind (and no output
        directory)."""
        ok = self._fasta([(b'ok', self.SEQ_A), (b'fine', self.SEQ_C)], 'ok.fa')
        bad = self._fasta([(b'bad\xff', self.SEQ_B)], 'bad.fa')
        graph = self._graph(f'{ok} {bad}')
        cmd = (f'{METAGRAPH} annotate -p 1 -i {graph} --separately --anno-header '
               f'-o {os.path.join(self.d, "anno")} {ok} {bad}')
        res = subprocess.run([cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self._assert_refused(res, bad, 'record 1', 'not valid UTF-8', '626164ff')
        # the files in order pass, each annotated on its own
        good = self._fasta([(b'good', self.SEQ_B)], 'good.fa')
        cmd = (f'{METAGRAPH} annotate -p 1 -i {graph} --separately --anno-header '
               f'-o {os.path.join(self.d, "anno")} {ok} {good}')
        res = subprocess.run([cmd], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertEqual(['good.fa.column.annodbg', 'ok.fa.column.annodbg'],
                         sorted(os.listdir(os.path.join(self.d, 'anno'))))

    def test_index_header_coords_refuses_a_repeated_record_header(self):
        fasta = self._fasta([(b'dup', self.SEQ_A), (b'other', self.SEQ_B), (b'dup', self.SEQ_C)])
        graph = self._graph(fasta)
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords')
        self._assert_refused(res, fasta, 'record 3', 'repeats the header of record 1',
                             '647570', '--allow-unchecked-headers')
        # consecutive records too: a duplicate is a duplicate
        fasta = self._fasta([(b'dup', self.SEQ_A), (b'dup', self.SEQ_B)], 'consecutive.fa')
        res = self._annotate(graph, fasta, '--anno-filename --index-header-coords')
        self._assert_refused(res, fasta, 'record 2', 'repeats the header of record 1')
        # a repeat only counts within a file: two files may share a header
        other = self._fasta([(b'dup', self.SEQ_C)], 'other.fa')
        first = self._fasta([(b'dup', self.SEQ_A)], 'first.fa')
        res = self._annotate(graph, f'{first} {other}', '--anno-filename --index-header-coords')
        self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
        self.assertTrue(os.path.exists(os.path.join(self.d, 'anno.seqs')))

    def test_fwd_and_reverse_counts_a_records_reverse_complement_to_the_record(self):
        fasta = self._fasta([(b'a', self.SEQ_A), (b'b', self.SEQ_B)])
        graph = self._graph(fasta)
        kmers = 2 * (len(self.SEQ_A) - self.K + 1)
        for flags, out, total in (('', 'fwd', kmers), ('--fwd-and-reverse', 'both', 2 * kmers)):
            res = self._annotate(graph, fasta, f'--anno-filename --index-header-coords {flags}',
                                 out=out)
            self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
            stats = subprocess.run([f'{METAGRAPH} stats {os.path.join(self.d, out)}.seqs'],
                                   shell=True, stdout=PIPE, stderr=PIPE).stdout.decode()
            self.assertIn('total sequences: 2\n', stats, flags)
            self.assertIn(f'total k-mers: {total}\n', stats, flags)

    def test_allow_unchecked_headers_builds_with_a_warning(self):
        fasta = self._fasta([(b'dup', self.SEQ_A), (b'bad\xff', self.SEQ_B), (b'dup', self.SEQ_C)])
        graph = self._graph(fasta)
        for flags, written in (('--anno-header', 'labels.column.annodbg'),
                               ('--anno-filename --index-header-coords', 'records.seqs')):
            out = written.split('.')[0]
            res = self._annotate(graph, fasta, f'{flags} --allow-unchecked-headers', out=out)
            self.assertEqual(res.returncode, 0, res.stderr.decode('utf-8', 'replace'))
            self.assertTrue(os.path.exists(os.path.join(self.d, written)), written)
            err = res.stderr.decode('utf-8', 'replace')
            self.assertIn('[warning]', err)
            self.assertIn('Header names are not checked (--allow-unchecked-headers)', err)
        # unchecked, the mapping stores every record (a repeat is found at its first only)
        stats = subprocess.run([f'{METAGRAPH} stats {os.path.join(self.d, "records")}.seqs'],
                               shell=True, stdout=PIPE, stderr=PIPE).stdout.decode()
        self.assertIn('total sequences: 3\n', stats)


if __name__ == '__main__':
    unittest.main()
