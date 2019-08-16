#!/usr/bin/env python2.7
#
# coding: utf-8
#


__author__ = 'Mikhail Karasikov'


def get_js_sample_list(annotations):
    js_output = 'var samples = [\n'

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
                u'"Count": "{}"'.format(count),
            ]

            for pair in other_labels:
                attribute, value = pair.split('=')
                try:
                    tokens.append(u'"{}": {}'.format(attribute, float(value)))
                except:
                    tokens.append(u'"{}": "{}"'.format(attribute, value))

            js_output += '  { ' + ', '.join(tokens) + ' },\n'

        except:
            js_output += u'  {{ "Sample": "{}", "Count": {} }},\n'.format(x, count)

    js_output += '];\n'
    return js_output
