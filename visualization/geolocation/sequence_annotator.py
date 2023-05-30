#!/usr/bin/env python2.7
#
# coding: utf-8
#

import sys
import subprocess
from tempfile import NamedTemporaryFile
from helpers import get_js_sample_list


__author__ = 'Mikhail Karasikov'


class SequenceAnnotator():
    def __init__(self, executable, graph, annotation):
        self.executable = executable
        self.graph = graph
        self.annotation = annotation

    def annotate(self, sequence):
        try:
            f_in = NamedTemporaryFile('w', bufsize=100, suffix='.fa')
            f_in.write('>seq\n' + sequence)
            f_in.flush()

            result = subprocess.check_output([
                    self.executable, 'query',
                    '-i', self.graph,
                    '-a', self.annotation,
                    *('--query-mode matches'.split()),
                    f_in.name
                ],
                stderr=subprocess.STDOUT
            ).decode('utf-8')

            return result, ''

        except subprocess.CalledProcessError as e:
            return None, e.output

        except Exception as e:
            return None, str(e)


if __name__ == '__main__':
    executable = '../../metagraph/build/metagraph'
    if len(sys.argv) != 3:
        print("Usage: {} <graph> <annotation>".format(sys.argv[0]))
        exit(1)

    _, graph, annotation = sys.argv

    sequence = 'CAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGG';

    annotated_dbg = SequenceAnnotator(executable, graph, annotation)

    result, output = annotated_dbg.annotate(sequence)
    if result is not None:
        js_sample_list = get_js_sample_list(result)
        print(js_sample_list.encode('utf-8'))
    else:
        print('Error:')
        print(output)
