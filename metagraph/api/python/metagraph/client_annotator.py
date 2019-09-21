#!/usr/bin/env python
#
# coding: utf-8
#

import sys
import socket
import json
from .helpers import get_js_sample_list
from .connection import Connection


__author__ = 'Mikhail Karasikov'


class RemoteEngine():
    def __init__(self, host, port):
        self.engine_address = (host, port)

    def compute(self, message):
        """ send message -> receive result """
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        sock.connect(self.engine_address)

        connection = Connection(sock)

        connection.send_string(message)
        received = connection.receive_string()
        return received


if __name__ == '__main__':
    """ 1) start metagraph in the server mode
        2) run this script for testing
    """

    if len(sys.argv) != 3:
        print("Usage: {} <host> <port>".format(sys.argv[0]))
        exit(1)

    _, host, port = sys.argv

    sequence = 'CAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGG';

    annotated_dbg = RemoteEngine(host, int(port))

    message = json.dumps({
        "FASTA": "\n".join([ ">test",
                                sequence[:len(sequence) // 2 + 1],
                                sequence[len(sequence) // 2 + 1:] ]),
        "count_labels": True,
        "discovery_fraction": 0.75,
        "num_labels": 15
    })

    result = annotated_dbg.compute(message)
    js_sample_list = get_js_sample_list(result)
    print(js_sample_list)
