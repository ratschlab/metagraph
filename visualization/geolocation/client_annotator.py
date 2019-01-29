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


__author__ = 'Mikhail Karasikov'


BUFF_SIZE = 4096
MAX_STRING_SIZE = 10**8


class Connection:
    """ Class for communication through socket
    It guarantees correct transmission for TCP socket
    """

    def __init__(self, socket):
        self.socket = socket

    def send_string(self, string):
        self.socket.sendall(string.encode('utf-8'))

    def receive_string(self):
        received = []
        string_length = 0

        # Receive data by small chunks and reconstruct the string
        while True:
            received.append(self.socket.recv(BUFF_SIZE))
            string_length += len(received[-1])

            if string_length > MAX_STRING_SIZE:
                raise IOError("Received string length longer than maximum allowed" +
                                " (" + str(string_length) + " > " + str(MAX_STRING_SIZE) + ")")

            if len(received[-1]) < BUFF_SIZE:
                break

        return b"".join(received).decode("utf-8")


class ClientAnnotator():
    def __init__(self, host, port):
        self.host = host
        self.port = port

    def annotate(self, sequence):
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        try:
            sock.connect((self.host, self.port))

            connection = Connection(sock)
            connection.send_string(sequence)
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

    result, output = annotated_dbg.annotate(sequence)
    if result is not None:
        js_sample_list = get_js_sample_list(result)
        print(js_sample_list.encode('utf-8'))
    else:
        print('Error:')
        print(output)
