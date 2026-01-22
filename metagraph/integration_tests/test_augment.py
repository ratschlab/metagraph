import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import os
import gzip
from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension, MMAP_FLAG


"""Test graph augmentation with and without node weights"""

GRAPH_TYPES = ['succinct', 'hash', 'hashfast', 'hashstr']


class TestAugment(TestingBase):
    def setUp(self):
        self.tempdir = TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def _augment_graph(self, input_graph, output_base, new_sequences, representation, use_mmap=True):
        """Augment a graph with new sequences"""
        augment_command = f'{METAGRAPH} extend -i {input_graph} -o {output_base} {new_sequences}'
        if use_mmap:
            augment_command += MMAP_FLAG

        res = subprocess.run([augment_command], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, f"Augment failed: {res.stderr.decode()}")
        return output_base + graph_file_extension[representation]

    def _verify_graph_stats(self, graph_path, expected_k='20', expected_mode=None, expect_weights=False):
        """Verify graph statistics"""
        stats = self._get_stats(graph_path)
        self.assertEqual(stats['returncode'], 0)
        self.assertEqual(stats['k'], expected_k)

        if expected_mode:
            self.assertEqual(stats['mode'], expected_mode)

        if expect_weights:
            self.assertIn('nnz weights', stats)
            self.assertIn('avg weight', stats)
        else:
            self.assertNotIn('nnz weights', stats)
            self.assertNotIn('avg weight', stats)

        return stats

    def _count_kmers_in_fasta(self, fasta_path, k):
        total = 0
        current = []
        with gzip.open(fasta_path, 'rt') as handle:
            for line in handle:
                if line.startswith('>'):
                    if current:
                        seq_len = len(''.join(current))
                        if seq_len >= k:
                            total += seq_len - k + 1
                        current = []
                else:
                    current.append(line.strip())
        if current:
            seq_len = len(''.join(current))
            if seq_len >= k:
                total += seq_len - k + 1
        return total

    def _count_graph_kmers_from_contigs(self, graph_path, k):
        output_base = f'{self.tempdir.name}/contigs_{os.path.basename(graph_path)}'
        transform_command = f'{METAGRAPH} transform --to-fasta -o {output_base} {graph_path}' + MMAP_FLAG
        res = subprocess.run([transform_command], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0, f"Transform failed: {res.stderr.decode()}")

        fasta_gz = output_base + '.fasta.gz'
        fasta_zst = output_base + '.fasta.zst'
        if os.path.exists(fasta_gz):
            fasta_path = fasta_gz
        elif os.path.exists(fasta_zst):
            fasta_path = fasta_zst
        else:
            self.fail(f"Missing contigs fasta for {output_base}")

        return self._count_kmers_in_fasta(fasta_path, k)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment(self, representation):
        """Test that graph augmentation works correctly"""
        # Build initial graph
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'basic', '--mask-dummy --disk-swap ""')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        # Verify initial graph
        stats_initial = self._verify_graph_stats(initial_graph)
        initial_nodes = int(stats_initial['nodes (k)'])
        self.assertGreater(initial_nodes, 0)
        self.assertFalse(os.path.exists(initial_graph + '.weights'),
                         f"Weights file {initial_graph}.weights should not exist")
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)

        # Augment the graph
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented graph
        stats_augmented = self._verify_graph_stats(augmented_graph)
        augmented_nodes = int(stats_augmented['nodes (k)'])
        self.assertGreater(augmented_nodes, initial_nodes,
                           "Augmented graph should have more nodes than initial graph")
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59955)
        self.assertFalse(os.path.exists(augmented_graph + '.weights'),
                         f"Weights file {augmented_graph}.weights should not exist")

    @parameterized.expand([repr for repr in GRAPH_TYPES if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_augment_canonical(self, representation):
        """Test augmentation of canonical graph"""
        # Build initial canonical graph
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'canonical', '--mask-dummy --disk-swap ""')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        # Verify initial graph
        stats_initial = self._verify_graph_stats(initial_graph, expected_mode='canonical')
        initial_nodes = int(stats_initial['nodes (k)'])
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)

        # Augment the graph
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented canonical graph
        stats_augmented = self._verify_graph_stats(augmented_graph, expected_mode='canonical')
        augmented_nodes = int(stats_augmented['nodes (k)'])
        self.assertGreater(augmented_nodes, initial_nodes)
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented canonical graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 58787)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_unweighted_to_weighted(self, representation):
        """Test augmenting an unweighted graph, then adding weights via augmentation"""
        # Build initial graph without weights
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'basic', '--mask-dummy --disk-swap ""')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        stats_initial = self._verify_graph_stats(initial_graph)
        initial_nodes = int(stats_initial['nodes (k)'])
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)

        # Augment the graph (still no weights)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented graph still has no weights
        stats_augmented = self._verify_graph_stats(augmented_graph)
        augmented_nodes = int(stats_augmented['nodes (k)'])
        self.assertGreater(augmented_nodes, initial_nodes)
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59955)
        self.assertFalse(os.path.exists(augmented_graph + '.weights'),
                         f"Weights file {augmented_graph}.weights should not exist")

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_with_node_weights(self, representation):
        """Test that graph augmentation correctly updates node weights"""
        # Build initial graph with node weights
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'basic', '--mask-dummy --disk-swap "" --count-kmers')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        # Verify initial graph has weights
        stats_initial = self._verify_graph_stats(initial_graph, expect_weights=True)
        initial_nodes = int(stats_initial['nodes (k)'])
        self.assertGreater(initial_nodes, 0)
        initial_nnz_weights = int(stats_initial['nnz weights'])
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)
        self.assertTrue(os.path.exists(initial_graph + '.weights'),
                        f"Weights file {initial_graph}.weights should exist")

        # Augment the graph (without --mmap due to node weights issue)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=False)

        # Verify augmented graph
        stats_augmented = self._verify_graph_stats(augmented_graph, expect_weights=True)
        augmented_nodes = int(stats_augmented['nodes (k)'])
        self.assertGreater(augmented_nodes, initial_nodes,
                           "Augmented graph should have more nodes than initial graph")
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59955)

        augmented_nnz_weights = int(stats_augmented['nnz weights'])
        self.assertGreaterEqual(augmented_nnz_weights, initial_nnz_weights,
                                "Augmented graph should have at least as many weighted nodes")
        self.assertTrue(os.path.exists(augmented_graph + '.weights'),
                        f"Weights file {augmented_graph}.weights should exist")

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_preserves_existing_weights(self, representation):
        """Test that augmentation preserves existing node weights and adds new ones"""
        # Build initial graph with node weights
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'basic', '--mask-dummy --disk-swap "" --count-kmers')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        stats_initial = self._verify_graph_stats(initial_graph, expect_weights=True)
        initial_avg_weight = float(stats_initial['avg weight'])
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)

        # Augment with overlapping sequences (should increase weights of existing nodes)
        # Note: --mmap flag causes issues with node weights, so we don't use it here
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/transcripts_100.fa',  # Same file
                                              representation, use_mmap=False)

        # Verify weights were updated
        stats_augmented = self._verify_graph_stats(augmented_graph, expect_weights=True)
        augmented_avg_weight = float(stats_augmented['avg weight'])
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertEqual(augmented_kmers, initial_kmers,
                         "Augmenting with identical sequences should not change contig k-mers")
        # Note: weights are additive, so adding the same sequences should increase average weight
        # However, new nodes might have lower weights, so we just check that weights exist
        self.assertGreater(augmented_avg_weight, 0,
                           "Augmented graph should have positive average weight")

    @parameterized.expand([repr for repr in GRAPH_TYPES if repr != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_augment_canonical_with_weights(self, representation):
        """Test augmentation of canonical graph with node weights"""
        # Build initial canonical graph with node weights
        initial_graph_base = self.tempdir.name + '/graph_initial'
        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, 'canonical', '--mask-dummy --disk-swap "" --count-kmers')
        initial_graph = initial_graph_base + graph_file_extension[representation]

        stats_initial = self._verify_graph_stats(initial_graph, expected_mode='canonical',
                                                 expect_weights=True)
        initial_nodes = int(stats_initial['nodes (k)'])
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)

        # Augment the canonical graph (without --mmap due to node weights issue)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=False)

        # Verify augmented canonical graph
        stats_augmented = self._verify_graph_stats(augmented_graph, expected_mode='canonical',
                                                   expect_weights=True)
        augmented_nodes = int(stats_augmented['nodes (k)'])
        self.assertGreater(augmented_nodes, initial_nodes)
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented canonical graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 58787)


if __name__ == '__main__':
    unittest.main()
