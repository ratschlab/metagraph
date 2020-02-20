import os
import shlex
import time
from subprocess import Popen

from metagraph.client import Client

from base import TestingBase, METAGRAPH, TEST_DATA_DIR


class TestApi(TestingBase):
    graph_name = 'test_graph'

    def setUp(self):
        super().setUp()

        fasta = TEST_DATA_DIR + '/transcripts_100.fa'

        graph_path = self.tempdir.name + '/graph.dbg'
        annotation_path = self.tempdir.name + '/annotation.column.annodbg'

        port = 3456
        self._build_graph(fasta, graph_path, 6, 'succinct')

        self._annotate_graph(fasta, graph_path, annotation_path, 'column')
        print(os.listdir(self.tempdir.name))

        self.server_process = self._start_server(graph_path, annotation_path, port)

        print("waiting for the server to start up on " + str(self.server_process.pid))
        time.sleep(1)

        self.graph_client = Client()

        self.graph_client.add_graph('127.0.0.1', port, self.graph_name)

    def tearDown(self) -> None:
        self.server_process.kill()


    def _start_server(self, graph, annotation, port):
        construct_command = '{exe} server_query -i {graph} -a {annot} --port {port} --address 127.0.0.1 '.format(
            exe=METAGRAPH,
            graph=graph,
            annot=annotation,
            port=port
        )

        return Popen(shlex.split(construct_command))

    # do various queries
    def test_api_client_query(self):
        ret = self.graph_client.search_json("CCTCTGTGGAATCCAATCTGTCTTCCATCCTGCGTGGCCGAGGG")

        self.assertIn(self.graph_name, ret.keys())

        json_obj = ret[self.graph_name][0][0]
        self.assertEqual(len(json_obj['results']), 98)

        first_ret = json_obj['results'][0]

        self.assertEqual(first_ret['sampleCount'], 39)

        self.assertTrue(first_ret['sampleName'].startswith('ENST00000456328.2|ENSG00000223972.5|'))
        self.assertTrue('properties' not in first_ret.keys()) # doesn't have properties, so don't sent them

    def test_api_client_column_labels(self):
        ret = self.graph_client.column_labels()

        self.assertIn(self.graph_name, ret.keys())

        label_list = ret[self.graph_name]

        self.assertGreater(len(label_list), 0)
        self.assertTrue(all(l.startswith('ENST') for l in label_list))
