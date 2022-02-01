import json
import os
import shlex
import time
import unittest
import subprocess
import socket
import requests

import pandas as pd

from metagraph.client import GraphClientJson, MultiGraphClient
from concurrent.futures import Future
from parameterized import parameterized, parameterized_class

from base import TestingBase, METAGRAPH, TEST_DATA_DIR

PROTEIN_MODE = os.readlink(METAGRAPH).endswith("_Protein")


class TestAPIBase(TestingBase):
    @classmethod
    def setUpClass(cls, fasta_path, mode='basic', anno_repr='column'):
        super().setUpClass()

        graph_path = cls.tempdir.name + '/graph.dbg'
        annotation_path_base = cls.tempdir.name + '/annotation'
        annotation_path = annotation_path_base + f'.{anno_repr}.annodbg'

        cls._build_graph(fasta_path, graph_path, 6, 'succinct', mode=mode)
        cls._annotate_graph(fasta_path, graph_path, annotation_path_base, anno_repr)

        cls.host = '127.0.0.1'
        os.environ['NO_PROXY'] = cls.host
        cls.port = 3456
        num_retries = 100
        while num_retries > 0:
            cls.server_process = cls._start_server(cls, graph_path, annotation_path)
            try:
                cls.server_process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                break
            cls.port += 1
            num_retries -= 1
        if num_retries == 0:
            raise "Couldn't start server"

        wait_time_sec = 1
        print("Waiting {} sec for the server (PID {}) to start up".format(
            wait_time_sec, cls.server_process.pid), flush=True)
        time.sleep(wait_time_sec)

    @classmethod
    def tearDownClass(cls):
        cls.server_process.kill()

    def _start_server(self, graph, annotation):
        construct_command = f'{METAGRAPH} server_query -i {graph} -a {annotation} \
                                            --port {self.port} --address {self.host} -p {2}'

        return subprocess.Popen(shlex.split(construct_command))


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIRaw(TestAPIBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode)

    def setUp(self) -> None:
        self.raw_post_request = lambda cmd, payload: requests.post(
                                    url=f'http://{self.host}:{self.port}/{cmd}',
                                    data=payload)

    def test_api_raw_incomplete_json(self):
        ret = self.raw_post_request('search',
            '{"FASTA": ">query\\nAATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA", "discovery_fraction": 0.1,')

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Bad json received:", ret.json()['error'])

    def test_api_raw_invalid_params(self):
        payload = json.dumps({
                    "num_labels": 'not_a_number',
                    "FASTA": "\n".join([">query",
                                        'AATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA',
                                        ]),
                    "discovery_fraction": 1 / 100,
                    })

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Value is not convertible to Int.", ret.json()['error'])

    def test_api_raw_missing_params(self):
        payload = json.dumps({
            "num_labels": 100,
            "discovery_fraction": 1 / 100
        })

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("No input sequences received from client", ret.json()['error'])

    def test_api_raw_invalid_discovery_fraction(self):
        payload = json.dumps({
            "FASTA": "\n".join([">query",
                                'TCGA',
                                ]),
            "discovery_fraction": 1.1,
            "num_labels": 1,
        })
        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)

    def test_api_raw_missing_num_labels(self):
        payload = json.dumps({
            "FASTA": "\n".join([">query",
                                'TCGA',
                                ]),
            "discovery_fraction": 0.1
        })
        ret = self.raw_post_request('search', payload)

        # server has uses some default value
        self.assertEqual(ret.status_code, 200)

    def test_api_raw_invalid_url(self):
        ret = self.raw_post_request('not_valid', {})
        self.assertEqual(ret.status_code, 404)

    def test_api_raw_no_sequence(self):
        payload = json.dumps({
                    "FASTA": "\n".join([">query",
                                        'SEQUENCE_NOT_IN_GRAPH',
                                        ]),
                    "discovery_fraction": 1 / 100,
                    "num_labels": 1,
                    })
        ret = self.raw_post_request('search', payload)
        json_ret = ret.json()

        self.assertEqual(ret.status_code, 200)
        self.assertEqual(len(json_ret), 1)
        self.assertEqual(json_ret[0]['results'], [])

    @parameterized.expand([(1,1), (3,1)])
    @unittest.skipIf(PROTEIN_MODE, "Protein graphs have different alignments")
    def test_api_raw_align_sequence(self, repetitions, dummy_arg):
        fasta_str = '\n'.join([ f">query{i}\nTCGATCGA" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str, "min_exact_match": 0})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'seq_description': 'query0',
                    'alignments': [{'score': 21, 'sequence': 'TCGATCAA', 'cigar': '6=1X1='}]}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )

    @parameterized.expand([(1,1), (3,1)])
    def test_api_raw_align_bad_sequence(self, repetitions, dummy_arg):
        fasta_str = '\n'.join([ f">query{i}\nNNNNNNNN" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str, "min_exact_match": 0})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'seq_description': 'query0', 'alignments': []}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )

    def test_api_raw_align_no_sequence(self):
        fasta_str = ">query\n"
        payload = json.dumps({"FASTA": fasta_str, "min_exact_match": 0})

        ret = self.raw_post_request('align', payload).json()

        self.assertEqual(ret[0]['alignments'], [])

    def test_api_raw_align_empty_fasta_desc(self):
        fasta_str = ">\nTCGATCGA"
        payload = json.dumps({"FASTA": fasta_str, "min_exact_match": 0})
        ret = self.raw_post_request('align', payload).json()

        self.assertEqual(ret[0]['seq_description'], '')

    def test_api_raw_search_empty_fasta_desc(self):
        fasta_str = ">\nCCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG"
        payload = json.dumps({"FASTA": fasta_str, 'num_labels': 5, 'min_exact_match': 0.1})
        ret = self.raw_post_request('search', payload).json()

        self.assertEqual(ret[0]['seq_description'], '')

    def test_api_raw_search_no_coordinate_support(self):
        fasta_str = ">query\nCCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG"
        payload = json.dumps({"FASTA": fasta_str, 'num_labels': 5, 'min_exact_match': 0.1,
                              'query_coords': True})

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Annotation does not support k-mer coordinate queries", ret.json()['error'])

    def test_api_raw_search_no_count_support(self):
        fasta_str = ">query\nCCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG"
        payload = json.dumps({"FASTA": fasta_str, 'num_labels': 5, 'min_exact_match': 0.1,
                              'abundance_sum': True})

        ret = self.raw_post_request('search', payload)

        self.assertEqual(ret.status_code, 400)
        self.assertIn("Annotation does not support k-mer count queries", ret.json()['error'])


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIClient(TestAPIBase):
    graph_name = 'test_graph'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode)

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

        # 'canonical' and 'primary' graphs represent more k-mers than 'basic', so
        # they get more matches
        cls.sample_query_expected_rows = 98 if cls.mode == 'basic' else 99
        cls.expected_matches = 840 if cls.mode == 'basic' else 1381

    def test_api_multiple_query_df(self):
        repetitions = 5
        ret = self.graph_client.search([self.sample_query] * repetitions, parallel=False,
                                       discovery_threshold=0.01)
        df = ret[self.graph_name]
        self.assertEqual((self.sample_query_expected_rows * repetitions, 3), df.shape)
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches * repetitions)

    def test_api_simple_query_df(self):
        ret = self.graph_client.search(self.sample_query, parallel=False,
                                       discovery_threshold=0.01)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)

    def test_api_simple_query_with_signature_df(self):
        ret = self.graph_client.search(self.sample_query, parallel=False,
                                       discovery_threshold=0.01, with_signature=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 4))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)

    def test_api_simple_query_align_df(self):
        ret = self.graph_client.search(self.sample_query, parallel=False,
                                       discovery_threshold=0.01, align=True, min_exact_match=0.01)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)

    def test_api_simple_query_align_no_result(self):
        # If aligned sequence does not have result when searched, make sure client doesn't fail
        sample_align_query = self.sample_query[:2] + 'GG' + self.sample_query[4:]
        ret = self.graph_client.search(sample_align_query, parallel=False,
                                       discovery_threshold=1.0, align=True)
        df = ret[self.graph_name]
        self.assertTrue(df.empty)

    def test_api_client_column_labels(self):
        ret = self.graph_client.column_labels()

        self.assertIn(self.graph_name, ret.keys())

        label_list = ret[self.graph_name]

        self.assertEqual(len(label_list), 100)
        self.assertTrue(all(l.startswith('ENST') for l in label_list))

    def test_api_align_df(self):
        repetitions = 4
        alignment_cnt = 3
        seq = ["TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"]
        ret = self.graph_client.align(seq * repetitions, parallel=False,
                                      max_alternative_alignments=alignment_cnt)

        align_res = ret[self.graph_name]
        self.assertIn('cigar', align_res.columns)
        # number of alignments returned per sequence is not necessarily equals max_alternative_alignments
        # but here it turns out to be the case
        self.assertEqual(len(align_res), repetitions *  alignment_cnt)

    def test_api_align_df_too_divergent(self):
        repetitions = 4
        alignment_cnt = 3
        seq = ["TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"]
        ret = self.graph_client.align(seq * repetitions, parallel=False,
                                      max_alternative_alignments=alignment_cnt, min_exact_match=1.0)

        align_res = ret[self.graph_name]
        self.assertIn('cigar', align_res.columns)
        self.assertEqual(len(align_res), 0)

    @unittest.expectedFailure
    def test_api_search_no_coordinate_support(self):
        ret = self.graph_client.search(self.sample_query, parallel=False,
                                       discovery_threshold=0.01, query_coords=True)

    @unittest.expectedFailure
    def test_api_search_no_count_support(self):
        ret = self.graph_client.search(self.sample_query, parallel=False,
                                       discovery_threshold=0.01, abundance_sum=True)


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIJson(TestAPIBase):
    graph_name = 'test_graph'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode)

        cls.graph_client = GraphClientJson(cls.host, cls.port)

        # 'canonical' and 'primary' graphs represent more k-mers than 'basic', so
        # they get more matches
        cls.sample_query_expected_rows = 98 if cls.mode == 'basic' else 99

    def setUp(self):
        if not self.graph_client.ready():
            self.fail("Server takes too long to initialize")

    def test_api_align_json(self):
        ret = self.graph_client.align("TCGATCGA")
        self.assertEqual(len(ret), 1)

    # do various queries
    def test_api_simple_query(self):
        res_list = self.graph_client.search(self.sample_query, discovery_threshold=0.01)
        self.assertEqual(len(res_list), 1)

        res_obj = res_list[0]['results']
        self.assertEqual(len(res_obj), self.sample_query_expected_rows)

        first_res = sorted(res_obj, key=lambda k: k['kmer_count'], reverse=True)[0]

        self.assertEqual(first_res['kmer_count'], 39)

        self.assertTrue(first_res['sample'].startswith('ENST00000456328.2|ENSG00000223972.5|'))
        self.assertTrue('properties' not in first_res.keys()) # doesn't have properties, so don't sent them

    def test_api_multiple_queries(self):
        repetitions = 4

        res_list = self.graph_client.search([self.sample_query] * repetitions)
        self.assertEqual(len(res_list), repetitions)

        # testing if the returned query indices range from 0 to n - 1
        self.assertEqual(sorted(range(0, repetitions)),
                         sorted([int(a['seq_description']) for a in res_list]))

    def test_api_stats(self):
        res = self.graph_client.stats()

        self.assertIn("graph", res.keys())
        graph_props = res['graph']
        self.assertEqual(graph_props["k"], 6)


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIClientWithProperties(TestAPIBase):
    """
    Testing whether properties encoded in sample name are properly processed
    """
    graph_name = 'test_graph_metasub'

    sample_query = 'GCCAGCATAGTGCTCCTGGACCAGTGATACACCCGGCACCCTGTCCTGGA'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/metasub_fake_data.fa', mode=cls.mode)

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

    def test_api_search_property_df(self):
        df = self.graph_client.search(self.sample_query, parallel=False)[self.graph_name]

        self.assertIn('kmer_count', df.columns)
        self.assertEqual(df['kmer_count'].dtype, int)
        self.assertIn('latitude', df.columns)
        self.assertEqual(df['latitude'].dtype, float)
        self.assertEqual(df.shape, (3, 10))

    def test_api_search_property_df_empty(self):
        df = self.graph_client.search("THISSEQUENCEDOESNOTEXIST", parallel=False)[self.graph_name]
        self.assertTrue(df.empty)


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIClientWithCoordinates(TestAPIBase):
    """
    Testing whether API works well given coordinate aware annotations
    """
    graph_name = 'test_client_coord'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode,
                           anno_repr='row_diff_brwt_coord')

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

        # 'canonical' and 'primary' graphs represent more k-mers than 'basic', so
        # they get more matches
        cls.sample_query_expected_rows = 98 if cls.mode == 'basic' else 99
        cls.expected_matches = 840 if cls.mode == 'basic' else 1381
        cls.expected_abundance_sum = 1176 if cls.mode == 'basic' else 2323

    def test_api_simple_query_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)

    def test_api_simple_query_abundance_sum_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, abundance_sum=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_abundance_sum)

    def test_api_simple_query_counts_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, query_counts=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 4))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)
        self.assertEqual(df['kmer_abundances'].size, self.sample_query_expected_rows)

    def test_api_simple_query_coords_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, query_coords=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 4))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)
        self.assertEqual(df['kmer_coords'].size, self.sample_query_expected_rows)


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIClientWithCounts(TestAPIBase):
    """
    Testing whether API works well given k-mer count aware annotations
    """
    graph_name = 'test_client_count'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode, anno_repr='row_diff_int_brwt')

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

        # 'canonical' and 'primary' graphs represent more k-mers than 'basic', so
        # they get more matches
        cls.sample_query_expected_rows = 98 if cls.mode == 'basic' else 99
        cls.expected_matches = 840 if cls.mode == 'basic' else 1381
        cls.expected_abundance_sum = 1176 if cls.mode == 'basic' else 2323

    def test_api_simple_query_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)

    def test_api_simple_query_abundance_sum_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, abundance_sum=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 3))
        self.assertEqual(df['kmer_count'].sum(), self.expected_abundance_sum)

    def test_api_simple_query_counts_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, query_counts=True)
        df = ret[self.graph_name]

        self.assertEqual(df.shape, (self.sample_query_expected_rows, 4))
        self.assertEqual(df['kmer_count'].sum(), self.expected_matches)
        self.assertEqual(df['kmer_abundances'].size, self.sample_query_expected_rows)

    @unittest.expectedFailure
    def test_api_search_no_coordinate_support(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01,
                                       parallel=False, query_coords=True)


# No canonical mode for Protein alphabets
@parameterized_class(('mode',),
    input_values=[('basic',)] + ([] if PROTEIN_MODE else [('canonical',), ('primary',)]))
class TestAPIClientParallel(TestAPIBase):
    """
    Testing whether or not parallel requests work
    """
    graph_names = ['test_graph', 'test_graph_2']

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'

    sample_align = ['TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA']

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa', mode=cls.mode)

        cls.graph_client = MultiGraphClient()

        # Here just query the same endpoint twice to test client
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_names[0])
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_names[1])

        # 'canonical' and 'primary' graphs represent more k-mers than 'basic', so
        # they get more matches
        cls.sample_query_expected_rows = 98 if cls.mode == 'basic' else 99

    def test_api_parallel_query_df(self):
        futures = self.graph_client.search(self.sample_query, parallel=True,
                                           discovery_threshold=0.01)

        self.assertEqual(len(futures), 2)

        for graph in self.graph_names:
            self.assertIn(graph, futures)
            self.assertIsInstance(futures[graph], Future)

        dfs = MultiGraphClient.wait_for_result(futures)

        for graph in self.graph_names:
            self.assertIn(graph, dfs)
            self.assertIsInstance(dfs[graph], pd.DataFrame)
            self.assertEqual((self.sample_query_expected_rows, 3), dfs[graph].shape)

    def test_api_parallel_align_df(self):
        repetitions = 4
        alignment_cnt = 3

        futures = self.graph_client.align(self.sample_align * repetitions, parallel=True,
                                          max_alternative_alignments=alignment_cnt)

        self.assertEqual(len(futures), 2)

        for graph in self.graph_names:
            self.assertIn(graph, futures)
            self.assertIsInstance(futures[graph], Future)

        dfs = MultiGraphClient.wait_for_result(futures)

        for graph in self.graph_names:
            self.assertIn(graph, dfs)
            self.assertIsInstance(dfs[graph], pd.DataFrame)
            self.assertIn('cigar', dfs[graph].columns)
            self.assertEqual(len(dfs[graph]), repetitions *  alignment_cnt)

    def test_api_parallel_query_error(self):
        futures = self.graph_client.search(self.sample_query, parallel=True,
                                           discovery_threshold=1.2)

        self.assertEqual(len(futures), 2)

        for graph in self.graph_names:
            self.assertIn(graph, futures)
            self.assertIsInstance(futures[graph], Future)

        dfs = MultiGraphClient.wait_for_result(futures)

        for graph in self.graph_names:
            self.assertIsInstance(dfs[graph], ValueError)
