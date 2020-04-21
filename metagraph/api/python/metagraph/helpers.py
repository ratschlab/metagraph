#!/usr/bin/env python2.7
#
# coding: utf-8
#


__author__ = 'Mikhail Karasikov'


def get_js_sample_list(annotations):
    js_output = []

    for x in annotations.split('\t'):
        try:
            x = x.split(':')
            if len(x) < 2:
                continue

            count = int(x[-1])

            x = ':'.join(x[:-1])

            if x.startswith('<') and x.endswith('>'):
                x = x[1:-1]
            else:
                continue
        except:
            continue

        try:
            labels = x.split(';')
            sample_name = labels[0]
            other_labels = labels[1:]

            tokens = [
                u'"Sample": "{}"'.format(sample_name),
                u'"Count": {}'.format(int(count)),
            ]

            for pair in other_labels:
                attribute, value = pair.split('=')
                try:
                    tokens.append(u'"{}": {}'.format(attribute, float(value)))
                except:
                    tokens.append(u'"{}": "{}"'.format(attribute, value))

            js_output.append('  { ' + ', '.join(tokens) + ' }')

        except:
            js_output.append(u'  {{ "Sample": "{}", "Count": {} }}'.format(x, count))

    return '[\n' + ',\n'.join(js_output) + '\n]'
