import json
import os
import shlex
import time
from subprocess import Popen

import socket
import requests
from metagraph.client import GraphClientJson, MultiGraphClient
from parameterized import parameterized, parameterized_class

from base import TestingBase, METAGRAPH, TEST_DATA_DIR


class TestAPIBase(TestingBase):
    @classmethod
    def setUpClass(cls, fasta_path, canonical=False, primary=False):
        super().setUpClass()

        graph_path = cls.tempdir.name + '/graph.dbg'
        annotation_path = cls.tempdir.name + '/annotation.column.annodbg'

        cls._build_graph(cls, fasta_path, graph_path, 6, 'succinct',
                         canonical=canonical, primary=primary)
        cls._annotate_graph(cls, fasta_path, graph_path, annotation_path, 'column',
                            primary=primary)

        cls.host = socket.gethostbyname(socket.gethostname())
        cls.port = 3456
        os.environ['NO_PROXY'] = cls.host
        cls.server_process = cls._start_server(cls, graph_path, annotation_path)

        wait_time_sec = 1
        print("Waiting {} sec for the server (PID {}) to start up".format(wait_time_sec, cls.server_process.pid), flush=True)
        time.sleep(wait_time_sec)

    @classmethod
    def tearDownClass(cls):
        cls.server_process.kill()

    def _start_server(self, graph, annotation):
        construct_command = '{exe} server_query -i {graph} -a {annot} --port {port} --address {host} --align-local -p {threads}'.format(
            exe=METAGRAPH,
            graph=graph,
            annot=annotation,
            host=self.host,
            port=self.port,
            threads=2
        )

        return Popen(shlex.split(construct_command))


@parameterized_class(('mode',), input_values=[('canonical',), ('primary',)])
class TestAPIRaw(TestAPIBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa',
                           canonical=True, primary=(cls.mode == 'primary'))

    def setUp(self) -> None:
        self.raw_post_request = lambda cmd, payload: requests.post(url=f'http://{self.host}:{self.port}/{cmd}', data=payload)

    def test_api_raw_incomplete_json(self):
        ret = self.raw_post_request('search', '{"FASTA": ">query\\nAATAAAGGTGTGAGATAACCCCAGCGGTGCCAGGATCCGTGCA", "discovery_fraction": 0.1,')

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

        self.assertEqual(ret.status_code, 200)
        self.assertEqual(ret.json(), [])

    @parameterized.expand([(1,1), (3,1)])
    def test_api_raw_align_sequence(self, repetitions, dummy_arg):
        fasta_str = '\n'.join([ f">query{i}\nTCGATCGA" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str, "discovery_fraction": 0})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'seq_description': 'query0',
                    'alignments': [{'score': 12, 'sequence': 'TCGATC', 'cigar': '6=2S'}]}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )

    @parameterized.expand([(1,1), (3,1)])
    def test_api_raw_align_bad_sequence(self, repetitions, dummy_arg):
        fasta_str = '\n'.join([ f">query{i}\nNNNNNNNN" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str, "discovery_fraction": 0})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'seq_description': 'query0', 'alignments': []}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )

    def test_api_raw_align_empty_fasta_desc(self):
        fasta_str = ">\nTCGATCGA"
        payload = json.dumps({"FASTA": fasta_str, "discovery_fraction": 0})
        ret = self.raw_post_request('align', payload).json()

        self.assertEqual(ret[0]['seq_description'], '')

    def test_api_raw_search_empty_fasta_desc(self):
        fasta_str = ">\nCCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG"
        payload = json.dumps({"FASTA": fasta_str, 'num_labels': 5, 'discovery_fraction': 0.1})
        ret = self.raw_post_request('search', payload).json()

        self.assertEqual(ret[0]['seq_description'], '')

#@parameterized_class(('mode',), input_values=[('canonical',), ('primary',)])
@parameterized_class(('mode',), input_values=[('canonical',)])
class TestAPIRawChain(TestAPIBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa',
                           canonical=True, primary=(cls.mode == 'primary'))

    def setUp(self) -> None:
        self.raw_post_request = lambda cmd, payload: requests.post(url=f'http://{self.host}:{self.port}/{cmd}', data=payload)

    def _start_server(self, graph, annotation):
        construct_command = '{exe} server_query --align-min-seed-length 2 -i {graph} -a {annot} --port {port} --address {host} -p {threads}'.format(
            exe=METAGRAPH,
            graph=graph,
            annot=annotation,
            host=self.host,
            port=self.port,
            threads=2
        )

        return Popen(shlex.split(construct_command))

    @parameterized.expand([(1, 1), (3, 1)])
    def test_api_raw_align_sequence(self, repetitions, dummy_arg):
        fasta_str = '\n'.join([ f">query{i}\nTCGATCGA" for i in range(repetitions)])

        payload = json.dumps({"FASTA": fasta_str, "discovery_fraction": 0})

        ret = self.raw_post_request('align', payload)

        self.assertEqual(ret.status_code, 200)

        self.assertEqual(len(ret.json()), repetitions)
        expected = {'seq_description': 'query0',
                    'alignments': [{'score': 12, 'sequence': 'TCGATC', 'cigar': '6=2S'},
                                   {'score': 4, 'sequence': 'GA', 'cigar': '6S2='}]}
        self.assertDictEqual(ret.json()[0], expected)

        self.assertListEqual(
            [ret.json()[i]['seq_description'] for i in range(repetitions)],
            [f"query{i}" for i in range(repetitions)]
        )


@parameterized_class(('mode',), input_values=[('canonical',), ('primary',)])
class TestAPIClient(TestAPIBase):
    graph_name = 'test_graph'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'
    sample_query_expected_rows = 99

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa',
                           canonical=True, primary=(cls.mode == 'primary'))

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

    def test_api_multiple_query_df(self):
        repetitions = 5
        ret = self.graph_client.search([self.sample_query] * repetitions, discovery_threshold=0.01)
        df = ret[self.graph_name]
        self.assertEqual((self.sample_query_expected_rows * repetitions, 3), df.shape)

    def test_api_simple_query_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01)
        df = ret[self.graph_name]

        self.assertEqual((self.sample_query_expected_rows, 3), df.shape)

    def test_api_simple_query_align_df(self):
        ret = self.graph_client.search(self.sample_query, discovery_threshold=0.01, align=True)
        df = ret[self.graph_name]

        self.assertIn('cigar', df.columns)
        self.assertEqual((self.sample_query_expected_rows, 6), df.shape)

    def test_api_client_column_labels(self):
        ret = self.graph_client.column_labels()

        self.assertIn(self.graph_name, ret.keys())

        label_list = ret[self.graph_name]

        self.assertGreater(len(label_list), 0)
        self.assertTrue(all(l.startswith('ENST') for l in label_list))

    def test_api_align_df(self):
        repetitions = 4
        alignment_cnt = 3
        seq = ["TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"]
        ret = self.graph_client.align(seq * repetitions, max_alternative_alignments=alignment_cnt)

        align_res = ret[self.graph_name]
        self.assertIn('cigar', align_res.columns)
        # number of alignments returned per sequence is not necessarily equals max_alternative_alignments
        # but here it turns out to be the case
        self.assertEqual(len(align_res), repetitions *  alignment_cnt)

    def test_api_align_df_too_divergent(self):
        repetitions = 4
        alignment_cnt = 3
        seq = ["TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"]
        ret = self.graph_client.align(seq * repetitions, max_alternative_alignments=alignment_cnt, discovery_threshold=1.0)

        align_res = ret[self.graph_name]
        self.assertIn('cigar', align_res.columns)
        self.assertEqual(len(align_res), 0)


@parameterized_class(('mode',), input_values=[('canonical',), ('primary',)])
class TestAPIJson(TestAPIBase):
    graph_name = 'test_graph'

    sample_query = 'CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG'
    sample_query_expected_cols = 99

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/transcripts_100.fa',
                           canonical=True, primary=(cls.mode == 'primary'))

        cls.graph_client = GraphClientJson(cls.host, cls.port)

    def test_api_align_json(self):
        ret, _ = self.graph_client.align("TCGATCGA")
        self.assertEqual(len(ret), 1)

    # do various queries
    def test_api_simple_query(self):
        res_list, _ = self.graph_client.search(self.sample_query, discovery_threshold=0.01)

        self.assertEqual(len(res_list), 1)

        res_obj = res_list[0]['results']
        self.assertEqual(len(res_obj), self.sample_query_expected_cols)

        first_res = res_obj[0]

        self.assertEqual(first_res['kmer_count'], 39)

        self.assertTrue(first_res['sample'].startswith('ENST00000456328.2|ENSG00000223972.5|'))
        self.assertTrue('properties' not in first_res.keys()) # doesn't have properties, so don't sent them

    def test_api_multiple_queries(self):
        repetitions = 4

        res_list, _ = self.graph_client.search([self.sample_query] * repetitions)
        self.assertEqual(len(res_list), repetitions)

        # testing if the returned query indices range from 0 to n - 1
        self.assertEqual(sorted(range(0, repetitions)), sorted([int(a['seq_description']) for a in res_list]))

    def test_api_stats(self):
        res = self.graph_client.stats()[0]

        self.assertIn("graph", res.keys())
        graph_props = res['graph']
        self.assertEqual(graph_props["k"], 6)


@parameterized_class(('mode',), input_values=[('canonical',), ('primary',)])
class TestAPIClientWithProperties(TestAPIBase):
    """
    Testing whether properties encoded in sample name are properly processed
    """
    graph_name = 'test_graph_metasub'

    sample_query = 'GCCAGCATAGTGCTCCTGGACCAGTGATACACCCGGCACCCTGTCCTGGA'

    @classmethod
    def setUpClass(cls):
        super().setUpClass(TEST_DATA_DIR + '/metasub_fake_data.fa',
                           canonical=True, primary=(cls.mode == 'primary'))

        cls.graph_client = MultiGraphClient()
        cls.graph_client.add_graph(cls.host, cls.port, cls.graph_name)

    def test_api_search_property_df(self):
        df = self.graph_client.search(self.sample_query)[self.graph_name]

        self.assertIn('kmer_count', df.columns)
        self.assertEqual(df['kmer_count'].dtype, int)
        self.assertIn('latitude', df.columns)
        self.assertEqual(df['latitude'].dtype, float)
        self.assertEqual(df.shape, (3, 10))

    def test_api_search_property_df_empty(self):
        df = self.graph_client.search("THISSEQUENCEDOESNOTEXIST")[self.graph_name]
        self.assertTrue(df.empty)
