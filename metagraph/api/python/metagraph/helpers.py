import uuid
import warnings

import pandas as pd


class SearchRefusedWarning(RuntimeWarning):
    """A /search answer in which the server refused some query sequences on some graphs (the
    request's `partial_refusal`, SPEC-multigraph-envelope.md section 3): their rows are absent
    from the data frame, which lists them in its attrs["refused"]."""


# master's n1 pair id (SPEC-multigraph-envelope.md section 1.2; src/common/configured_pair_id.hpp):
# "n1-" and the UUIDv5 in this namespace of the graph name and the two paths as the graph list
# spells them, each "<UTF-8 byte length>:<text>". A copy of the traverse library's
# configured_pair_id (train 7), until that library is merged here
_PAIR_ID_NAMESPACE = uuid.UUID('d3565c73-e7e9-5bca-b6a9-64cb5c4c9f93')


def _configured_pair_id(graph, graph_path, annotation_path):
    name = ''.join(f'{len(f.encode())}:{f}' for f in (graph, graph_path, annotation_path))
    return 'n1-' + str(uuid.uuid5(_PAIR_ID_NAMESPACE, name))


def _pair_id_of(cut):
    """A cut's pair_id: the one it states (a server of train 8), else the n1 id of the
    graph_path and annotation_path an older server states, None where it names no pair (a
    single-graph server)"""
    if cut.get('pair_id') is not None:
        return cut['pair_id']
    graph, graph_path, annotation_path = (cut.get(k) for k in
                                          ('graph', 'graph_path', 'annotation_path'))
    if all(isinstance(v, str) for v in (graph, graph_path, annotation_path)):
        return _configured_pair_id(graph, graph_path, annotation_path)
    return None


def refused_units(json_res):
    """The units (a query sequence on a graph) a /search answer refused alone: one dict per
    refused `cuts` entry, {"record": the query sequence's index in the answer,
    "seq_description", "graph", "pair_id" (the entry's pair: its name and its opaque id,
    null on a single-graph server; an older server's cut, which states the pair's paths in
    its place, gets the id of its paths), and the entry's `refused` members: "account",
    "code", "max_memory_mb", and with the answer account "priced_mb"}, in the answer's order.
    A refused unit never means that nothing matched: its rows were not sent."""
    refused = []
    for i, result in enumerate(json_res or []):
        for cut in result.get('cuts') or []:
            if isinstance(cut, dict) and isinstance(cut.get('refused'), dict):
                refused.append({'record': i,
                                'seq_description': result.get('seq_description'),
                                'graph': cut.get('graph'),
                                'pair_id': _pair_id_of(cut),
                                **cut['refused']})
    return refused


def df_from_search_result(json_res):
    def _build_dict(row, result):
        d = dict(row)
        # `properties` (optional): dictionary with metadata about the sample
        if 'properties' in d.keys():
            props = d.pop('properties')
        else:
            props = {}

        props['seq_description'] = result['seq_description']

        return {**d, **props}

    lst = [_build_dict(row, result) for result in json_res for row in
           result['results']]

    # columns of the table may vary on the graph, and inferred automatically
    df = pd.DataFrame(lst)
    # the units the server refused alone (the request's partial_refusal): the frame holds rows
    # only, so they are listed beside it, and warned about
    refused = refused_units(json_res)
    if refused:
        df.attrs['refused'] = refused
        named = '; '.join(
            f"record {u['record']} ({u['seq_description']!r}) on "
            + (f"graph {u['graph']!r}" if u['graph'] is not None else 'the graph')
            + (f" (pair {u['pair_id']})" if u['pair_id'] is not None else '')
            + f" ({u['account']} account)" for u in refused)
        warnings.warn(f"the server refused {len(refused)} query sequence(s) on a graph alone, "
                      f"past its answer budget ({refused[0]['code']}; their rows are absent, "
                      f"see the data frame's attrs['refused']): {named}",
                      SearchRefusedWarning, stacklevel=2)
    return df


def df_from_align_result(json_res):
    # flatten out json result
    lst = [(alignment['cigar'],
            alignment['score'],
            alignment['max_score'],
            alignment['sequence'],
            alignment['orientation'],
            result['seq_description'])
           for result in json_res for alignment in result['alignments']]

    df = pd.DataFrame(lst,
                      columns=['cigar', 'score', 'max_score', 'sequence', 'orientation', 'seq_description'])

    return df
