# -*- coding: utf-8 -*-

import time
import sys
import os
import traceback
import json
import pandas as pd
from socket import error as socket_error
from .client_annotator import RemoteEngine
from .helpers import get_js_sample_list


"""Metagraph client."""


class Client:
    def __init__(self):
        self.servers = []

    def connect(self, host, port, server_name=""):
        print("Trying to connect to Metagraph server {}:{}".format(host, port))
        annotated_dbg = RemoteEngine(host, port)
        self.servers.append((host, port, server_name))

    def search(self, sequence='AACGCTAATGTAGATTGAT', discovery_threshold=1.0):
        if len(self.servers) == 0:
            print("Error: connect to a database first")

        results = []

        for (host, port, server_name) in self.servers:
            try:
                annotated_dbg = RemoteEngine(host, port)

                task_context = json.dumps({
                    "FASTA": "\n".join([ ">query",
                                            sequence,
                    ]),
                    "count_labels": True,
                    "discovery_fraction": discovery_threshold / 100,
                    "num_labels": 10000
                })

                annotations = annotated_dbg.compute(task_context)

                js_sample_list = get_js_sample_list(annotations)

                results.append(pd.DataFrame(json.loads(js_sample_list)))

            except socket_error as e:
                traceback.print_exc(limit=20, file=sys.stderr)
                output = "Error: Server does not respond"

        return results
