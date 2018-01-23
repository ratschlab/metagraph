# coding: utf-8

import os
import re
import subprocess
import numpy as np
from scipy.io import loadmat
from tempfile import NamedTemporaryFile


__author__ = 'Mikhail Karasikov'


class DBGVisualizer():
    def __init__(self, executable):
        self.executable = executable

    def visualize(self, k, sequences):
        output = []
        try:
            f_in = NamedTemporaryFile('w', bufsize=100, suffix='.fa')
            f_in.write('>seq\n' + '\n>seq\n'.join(sequences))
            f_in.flush()

            graph_filename, out_log = self.__build_graph(k, f_in.name)
            output.append(out_log)

            if graph_filename is None:
                return None, output

            state_repr, out_log = self.__get_state_representation(graph_filename)
            output.append(out_log)

            os.remove(graph_filename + '.W.dbg')
            os.remove(graph_filename + '.l.dbg')
            os.remove(graph_filename + '.F.dbg')

            return state_repr, '\n'.join([x for x in output if len(x)])

        except Exception as e:
            output.append(str(e))
            return None, '\n'.join([x for x in output if len(x)])

    def __build_graph(self, k, fasta_filename):
        try:
            graph_file = NamedTemporaryFile(suffix='.graph')
            graph_file.close()

            output = subprocess.check_output(
                [ self.executable, 'build', '-k', str(k), '-o', graph_file.name, fasta_filename ],
                stderr=subprocess.STDOUT
            ).decode()

            return graph_file.name, output

        except subprocess.CalledProcessError as e:
            return None, e.output

    def __get_state_representation(self, graph_filename):
        try:
            state_repr = subprocess.check_output(
                [ self.executable, 'stats', '--print', graph_filename ],
                stderr=subprocess.STDOUT
            ).decode()
            return state_repr, ''

        except subprocess.CalledProcessError as e:
            return None, e.output


def get_js_edge_list(state_repr):
    js_output = 'var links = [\n'

    try:
        start = -1
        for i, line in enumerate(state_repr.split('\n')):
            if line.startswith('Index'):
                start = i
                break

        for line in state_repr.split('\n')[start + 1:]:
            tokens = line.split('\t')
            if len(tokens) != 4:
                break
            source = tokens[2]
            target = source[1:] + tokens[3][0]

            js_output += '  {{source: "{}", target: "{}"}},\n'.format(source, target)
    except:
        pass

    js_output += '];\n'
    return js_output


def to_html(state_repr):
    html_table = ('<table id="graph_representation"'
                        ' class="table table-striped table-bordered centered-table"'
                        ' style="width:auto; font: 12px monospace;" align="center">\n')
    try:
        start = -1
        for i, line in enumerate(state_repr.split('\n')):
            if line.startswith('Index'):
                start = i
                break

        html_table += ("<tr>" +
                        "<th>#</th>" +
                        "<th>{}</th>" * 2 +
                        "<th style='text-align:left;'>{}</th>" +
                       "</tr>\n").format(*line.split('\t')[1:])

        for line in state_repr.split('\n')[start + 1:]:
            tokens = line.split('\t')
            if len(tokens) != 4:
                break

            html_table += ("<tr>" +
                                "<th>{}</th>" +
                                "<td>{}</td>" * 2 +
                                "<td style='text-align:left;'>{}</td>" +
                           "</tr>\n").format(*tokens)
    except:
        pass

    html_table += '</table>'
    return html_table


if __name__ == '__main__':
    executable = '../../metagraph/build/metagengraph'

    k = 3
    sequences = [
        'ACGACGGAGCAGCAGCGAGCAGCGAGCAGCG',
        'TTTTTTTTTTTTTTTTTTTTTTT',
        'ACAGCAGCGAGCG4324324234',
        '',
        '\n\n'
    ]

    graph = DBGVisualizer(executable)

    state_repr, output = graph.visualize(k, sequences)
    print('Output:')
    print(output)

    if state_repr is not None:
        state_repr_html_table = to_html(state_repr)
        print('Representation:')
        print(state_repr_html_table)

        js_edge_list = get_js_edge_list(state_repr)
        print('Edge list:')
        print(js_edge_list)
