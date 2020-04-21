#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""TCP Connection Interface"""


import socket
import struct


MAX_STRING_SIZE = 10 ** 8
MAX_CHUNK_SIZE = 4096


class Connection:
    """ Socket communicator
    Correct and complete transmission for TCP sockets is guaranteed
    """

    def __init__(self, socket):
        self.socket = socket

    def send_string(self, string):
        string_length = len(string)
        self.socket.sendall(struct.pack('!L', string_length >> 32)
                                + struct.pack('!L', string_length & 0xFFFFFFFF)
                                + string.encode('utf-8'))

    def receive_string(self):
        data = b""
        while len(data) < 8:
            data += self.socket.recv(8 - len(data))

        string_length = ((struct.unpack('!L', data[:4])[0] << 32)
                            + struct.unpack('!L', data[4:])[0])

        if string_length > MAX_STRING_SIZE:
            raise IOError("Received string length longer than maximum allowed" +
                            " (" + str(string_length) + " > " + str(MAX_STRING_SIZE) + ")")
        # Receive data by small chunks and reconstruct the string
        received = []
        while string_length > 0:
            received.append(self.socket.recv(min(MAX_CHUNK_SIZE, string_length)))
            string_length -= len(received[-1])

        return b"".join(received).decode("utf-8")
