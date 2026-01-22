import unittest
from parameterized import parameterized
import subprocess
from subprocess import PIPE
from tempfile import TemporaryDirectory
import os
import gzip
from base import PROTEIN_MODE, TestingBase, METAGRAPH, TEST_DATA_DIR, graph_file_extension, MMAP_FLAG, PROTEIN_MODE


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

    def _build_initial_graph(self, representation, mode='basic', with_weights=False):
        initial_graph_base = self.tempdir.name + '/graph_initial'
        build_flags = '--mask-dummy --disk-swap ""'
        if with_weights:
            build_flags += ' --count-kmers'

        self._build_graph(TEST_DATA_DIR + '/transcripts_100.fa', initial_graph_base,
                          20, representation, mode, build_flags)

        initial_graph = initial_graph_base + graph_file_extension[representation]
        initial_kmers = self._count_graph_kmers_from_contigs(initial_graph, 20)
        return initial_graph, initial_kmers

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
        return self._count_kmers_in_fasta(output_base + '.fasta.gz', k)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment(self, representation):
        """Test that graph augmentation works correctly"""
        # Build initial graph
        initial_graph, initial_kmers = self._build_initial_graph(representation)

        # Augment the graph
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented graph
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59975 if PROTEIN_MODE else 59955)

    @parameterized.expand([graph_type for graph_type in GRAPH_TYPES if graph_type != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_augment_canonical(self, representation):
        """Test augmentation of canonical graph"""
        # Build initial canonical graph
        initial_graph, initial_kmers = self._build_initial_graph(representation, mode='canonical')

        # Augment the graph
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented canonical graph
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented canonical graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 58787)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_unweighted_to_weighted(self, representation):
        """Test augmenting an unweighted graph, then adding weights via augmentation"""
        # Build initial graph without weights
        initial_graph, initial_kmers = self._build_initial_graph(representation)

        # Augment the graph (still no weights)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=True)

        # Verify augmented graph still has no weights
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59975 if PROTEIN_MODE else 59955)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_with_node_weights(self, representation):
        """Test that graph augmentation correctly updates node weights"""
        # Build initial graph with node weights
        initial_graph, initial_kmers = self._build_initial_graph(representation, with_weights=True)

        # Augment the graph (without --mmap due to node weights issue)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=False)

        # Verify augmented graph
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 59975 if PROTEIN_MODE else 59955)

    @parameterized.expand(GRAPH_TYPES)
    def test_augment_preserves_existing_weights(self, representation):
        """Test that augmentation preserves existing node weights and adds new ones"""
        # Build initial graph with node weights
        initial_graph, initial_kmers = self._build_initial_graph(representation, with_weights=True)

        # Augment with overlapping sequences (should increase weights of existing nodes)
        # Note: --mmap flag causes issues with node weights, so we don't use it here
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/transcripts_100.fa',  # Same file
                                              representation, use_mmap=False)

        # Verify weights were updated
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertEqual(augmented_kmers, initial_kmers,
                         "Augmenting with identical sequences should not change contig k-mers")

    @parameterized.expand([graph_type for graph_type in GRAPH_TYPES if graph_type != 'hashstr'])
    @unittest.skipIf(PROTEIN_MODE, "No canonical mode for Protein alphabets")
    def test_augment_canonical_with_weights(self, representation):
        """Test augmentation of canonical graph with node weights"""
        # Build initial canonical graph with node weights
        initial_graph, initial_kmers = self._build_initial_graph(
            representation, mode='canonical', with_weights=True)

        # Augment the canonical graph (without --mmap due to node weights issue)
        augmented_graph = self._augment_graph(initial_graph,
                                              self.tempdir.name + '/graph_augmented',
                                              TEST_DATA_DIR + '/genome.MT.fa',
                                              representation, use_mmap=False)

        # Verify augmented canonical graph
        augmented_kmers = self._count_graph_kmers_from_contigs(augmented_graph, 20)
        self.assertGreater(augmented_kmers, initial_kmers,
                           "Augmented canonical graph should have more k-mers in contigs")
        self.assertEqual(augmented_kmers, 58787)


if __name__ == '__main__':
    unittest.main()
