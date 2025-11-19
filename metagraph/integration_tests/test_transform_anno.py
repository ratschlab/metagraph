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

        stats_graph = cls._get_stats(f'{cls.tempdir.name}/graph{graph_file_extension[cls.graph_repr]}')
        assert(stats_graph['returncode'] == 0)
        assert(stats_graph['k'] == '20')
        assert(stats_graph['nodes (k)'] == '46960')
        assert(stats_graph['mode'] == 'basic')

        cls.num_nodes = stats_graph['nodes (k)']
        cls.max_index = stats_graph['max index (k)']

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
        stats_annotation = self._get_stats(f'-a {self.annotation}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '100')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            0.0185072 * int(self.num_nodes) / int(self.max_index),
            places=6)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

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

        stats_annotation = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '1')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            float(expected_density) * int(self.num_nodes) / int(self.max_index),
            places=5)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

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

        stats_annotation = self._get_stats(f'-a aggregated{anno_file_extension[self.anno_repr]}')
        self.assertEqual(stats_annotation['returncode'], 0)
        self.assertEqual(stats_annotation['labels'], '1')
        self.assertEqual(stats_annotation['objects'], self.max_index)
        self.assertAlmostEqual(
            float(stats_annotation['density']),
            float(expected_density) * int(self.num_nodes) / int(self.max_index),
            places=5)
        self.assertEqual(stats_annotation['representation'], self.anno_repr)

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


class TestColumnOperationsWithCounts(TestingBase):
    """Test count aggregation by querying k-mers from the original FASTA files."""
    
    @classmethod
    def setUpClass(cls):
        cls.tempdir = TemporaryDirectory()
        cls.graph_repr = 'succinct'
        cls.anno_repr = 'column'
        
        # Define test datasets
        cls.datasets = ['logan_30', 'logan_30_alt']
        
        # Build graph from both files
        fasta_paths = [f'{TEST_DATA_DIR}/{dataset}.fa' for dataset in cls.datasets]
        cls._build_graph(fasta_paths,
                         cls.tempdir.name + '/graph',
                         20, cls.graph_repr, 'basic', '--mask-dummy')
        
        # Store test k-mers with their count arrays from each annotation
        # Format: kmer -> (logan_30_counts, logan_30_alt_counts)
        cls.test_kmers = {
            'AAAAAAAAAAAAAAAAAAAA': (
                [3,20,15,12,42,3,21,39,12,39,8,21,48,8,6,12,78,9,9,18,8,3,12,20,24,12,28,10,6,14],
                [9,44,65,18,183,30,72,177,18,162,12,60,200,8,9,44,612,36,33,39,24,24,24,44,72,39,136,18,14,32]
            ),
            'AAAAAAAAAAAAAAAAAAGC': (
                [2,2,3],
                [2,3,11]
            ),
            'AAAAAAAAAAAAAAAAAGCA': (
                [2],
                [2]
            ),
            'AAAAAAAAAAAAAAAAAAAC': (
                [1,5,3,4,14,1,7,13,14,5,3,7],
                [3,11,13,6,61,10,24,59,68,9,7,16]
            ),
            'AAAAAAAAAAAACCCAAAAT': (
                [14,5,3,7],
                [61,9,7,16]
            ),
            'AAAAAAAAAAAAGAGAGGGG': (
                [13],
                [54]
            ),
        }
        
        # Annotate each file with counts (use max count-width 16 for setup)
        for dataset in cls.datasets:
            fasta_path = f'{TEST_DATA_DIR}/{dataset}.fa'
            anno_base = f'{cls.tempdir.name}/{dataset}'
            
            cls._annotate_graph(
                fasta_path,
                cls.tempdir.name + '/graph' + graph_file_extension[cls.graph_repr],
                anno_base,
                cls.anno_repr,
                anno_type='header',
                extra_params='--count-kmers --count-width 16'
            )
        
        # Create query file with test k-mers once
        cls.test_query = f'{cls.tempdir.name}/test_kmers.fa'
        with open(cls.test_query, 'w') as f:
            for i, kmer in enumerate(cls.test_kmers.keys()):
                f.write(f'>kmer{i}\n{kmer}\n')
    
    @staticmethod
    def _convert_annotation(input_anno, output_base, output_repr):
        command = f'{METAGRAPH} transform_anno -o {output_base} --anno-type {output_repr} {input_anno}'
        res = subprocess.run(command.split(), stdout=PIPE, stderr=PIPE)
        assert res.returncode == 0, f"Conversion failed: {res.stderr.decode()}"
    
    def _query_counts_simple(self, anno_path, fasta_path):
        """Query k-mers and return {kmer: count} dict. For aggregated annotations only."""
        # Extract k-mers from query file
        kmers = []
        with open(fasta_path) as f:
            for line in f:
                if not line.startswith('>'):
                    kmers.append(line.strip().upper())
        
        # Query
        graph_path = self.tempdir.name + '/graph' + graph_file_extension[self.graph_repr]
        command = f'{METAGRAPH} query --query-mode counts --min-kmers-fraction-label 0 -i {graph_path} -a {anno_path} {fasta_path}'
        res = subprocess.run(command.split(), stdout=PIPE, stderr=PIPE)
        assert res.returncode == 0, f"Query failed: {res.stderr.decode()}"
        
        # Parse output: for aggregated anno, format is simple: "id\tname\t<mask>:0=count"
        counts = {}
        for line in res.stdout.decode().strip().split('\n'):
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                query_id = int(parts[0])
                results = parts[2]  # <mask>:0=count
                
                # Extract count from "<mask>:0=count"
                if ':' in results and '=' in results:
                    count_str = results.split('=')[-1]
                    count = int(count_str)
                    
                    if query_id < len(kmers):
                        kmer = kmers[query_id]
                        counts[kmer] = count
        
        return counts
    
    def _compute_expected_aggregated_count(self, logan_30_counts, logan_30_alt_counts,
                                           min_value=1, max_value=float('inf'),
                                           min_count=1, max_count=float('inf')):
        """
        Compute expected aggregated count based on transform_anno parameters.
        
        With --count-kmers:
          - First, filter counts by min_value and max_value
          - Then, sum all filtered counts
          - Finally, apply min_count and max_count to the aggregated sum
        
        Formula: min_count <= sum(filtered_counts) <= max_count
        """
        # Aggregate counts: sum all counts from both annotations
        all_counts = logan_30_counts + logan_30_alt_counts
        
        # Filter by min_value and max_value
        filtered_counts = [c for c in all_counts if min_value <= c <= max_value]
        
        # Sum filtered counts
        aggregated_sum = sum(filtered_counts)
        
        # Check if aggregated sum is in [min_count, max_count]
        if min_count <= aggregated_sum <= max_count:
            return aggregated_sum
        else:
            return 0
    
    def _aggregate_and_query(self, count_width=16, min_value=None, max_value=None,
                            min_count=None, max_count=None, output_name='aggregated'):
        """Helper to aggregate annotations and query test k-mers."""
        # Build command with parameters
        command = [METAGRAPH, 'transform_anno', '-p', str(NUM_THREADS),
                   '--count-kmers', '--count-width', str(count_width),
                   f'{self.tempdir.name}/{self.datasets[0]}{anno_file_extension[self.anno_repr]}',
                   f'{self.tempdir.name}/{self.datasets[1]}{anno_file_extension[self.anno_repr]}',
                   '--aggregate-columns']
        
        # Add filtering parameters if specified
        if min_value is not None:
            command.extend(['--min-value', str(min_value)])
        if max_value is not None:
            command.extend(['--max-value', str(max_value)])
        if min_count is not None:
            command.extend(['--min-count', str(min_count)])
        if max_count is not None:
            command.extend(['--max-count', str(max_count)])
        
        command.extend(['-o', f'{self.tempdir.name}/{output_name}'])
        
        res = subprocess.run(command, stdout=PIPE, stderr=PIPE)
        self.assertEqual(res.returncode, 0)
        
        # Convert aggregated to int_brwt for querying
        self._convert_annotation(
            self.tempdir.name + f'/{output_name}' + anno_file_extension[self.anno_repr],
            self.tempdir.name + f'/{output_name}_queryable',
            'int_brwt'
        )
        
        # Query aggregated annotation
        return self._query_counts_simple(
            self.tempdir.name + f'/{output_name}_queryable.int_brwt.annodbg',
            self.test_query
        )
    
    def test_aggregate_counts_basic(self):
        """Test basic count aggregation with default parameters (no filtering)."""
        agg_counts = self._aggregate_and_query(count_width=16, output_name='agg_basic')
        
        # Verify each test k-mer
        for kmer, (logan_30_counts, logan_30_alt_counts) in self.test_kmers.items():
            expected = self._compute_expected_aggregated_count(logan_30_counts, logan_30_alt_counts)
            actual = agg_counts.get(kmer, 0)
            
            self.assertEqual(actual, expected,
                f"K-mer {kmer}: expected {expected} (sum of {sum(logan_30_counts)}+{sum(logan_30_alt_counts)}), got {actual}")
    
    def test_aggregate_counts_with_value_filter(self):
        """Test count aggregation with min/max-value filtering."""
        # Filter: only count values in range [10, 100]
        agg_counts = self._aggregate_and_query(count_width=16, min_value=10, max_value=100,
                                               output_name='agg_value_filter')
        
        # Verify each test k-mer
        for kmer, (logan_30_counts, logan_30_alt_counts) in self.test_kmers.items():
            expected = self._compute_expected_aggregated_count(
                logan_30_counts, logan_30_alt_counts,
                min_value=10, max_value=100
            )
            actual = agg_counts.get(kmer, 0)
            
            self.assertEqual(actual, expected,
                f"K-mer {kmer}: expected {expected} with value filter [10,100], got {actual}")
    
    def test_aggregate_counts_with_count_filter(self):
        """Test count aggregation with min/max-count filtering."""
        # Filter: only include k-mers with aggregated count in [5, 100]
        agg_counts = self._aggregate_and_query(count_width=16, min_count=5, max_count=100,
                                               output_name='agg_count_filter')
        
        # Verify each test k-mer
        for kmer, (logan_30_counts, logan_30_alt_counts) in self.test_kmers.items():
            expected = self._compute_expected_aggregated_count(
                logan_30_counts, logan_30_alt_counts,
                min_count=5, max_count=100
            )
            actual = agg_counts.get(kmer, 0)
            
            self.assertEqual(actual, expected,
                f"K-mer {kmer}: expected {expected} with min_count=5, max_count=100, got {actual}")
    
    def test_aggregate_counts_with_width_saturation(self):
        """Test count aggregation with reduced count-width causing saturation."""
        # Use count-width 8 (max 255) - some counts will saturate
        agg_counts = self._aggregate_and_query(count_width=8, output_name='agg_width8')
        
        # Verify each test k-mer with saturation at 255
        for kmer, (logan_30_counts, logan_30_alt_counts) in self.test_kmers.items():
            # Compute expected with saturation: cap each count at 255, then sum
            saturated_30 = [min(c, 255) for c in logan_30_counts]
            saturated_alt = [min(c, 255) for c in logan_30_alt_counts]
            expected_raw = sum(saturated_30) + sum(saturated_alt)
            expected = min(expected_raw, 255)  # Final result also saturates at 255
            
            actual = agg_counts.get(kmer, 0)
            
            self.assertEqual(actual, expected,
                f"K-mer {kmer}: expected {expected} with count-width 8 (saturation at 255), got {actual}")


if __name__ == '__main__':
    unittest.main()
