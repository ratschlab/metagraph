#!/usr/bin/env python2.7
#
# coding: utf-8
#

import os
import sys
import re
import subprocess
import numpy as np
from scipy.io import loadmat
from tempfile import NamedTemporaryFile
import socket
import sys
import struct
from helpers import get_js_sample_list
from connection import Connection
import json


__author__ = 'Mikhail Karasikov'


class ClientAnnotator():
    def __init__(self, host, port):
        self.host = host
        self.port = port

    def annotate(self, json_message):
        # send message
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        try:
            sock.connect((self.host, self.port))

            connection = Connection(sock)
            connection.send_string(json_message)
            received = connection.receive_string()
            return received, ''

        except subprocess.CalledProcessError as e:
            return None, e.output

        except Exception as e:
            return None, str(e)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: {} <host> <port>".format(sys.argv[0]))
        exit(1)

    _, host, port = sys.argv

    sequence = 'CAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGG';

    annotated_dbg = ClientAnnotator(host, int(port))

    message = json.dumps({
        "FASTA": "\n".join([">test", sequence[:len(sequence) / 2 + 1], sequence[len(sequence) / 2 + 1:]]),
        "perc_similarity": 0.75,
        "num_labels": 10
    })

    result, output = annotated_dbg.annotate(message)
    if result is not None:
        js_sample_list = get_js_sample_list(result)
        print(js_sample_list.encode('utf-8'))
    else:
        print('Error:')
        print(output)
