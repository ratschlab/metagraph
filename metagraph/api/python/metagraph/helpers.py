import warnings

import pandas as pd


class SearchRefusedWarning(RuntimeWarning):
    """A /search answer in which the server refused some query sequences on some graphs (the
    request's `partial_refusal`, SPEC-multigraph-envelope.md section 3): their rows are absent
    from the data frame, which lists them in its attrs["refused"]."""


def refused_units(json_res):
    """The units (a query sequence on a graph) a /search answer refused alone: one dict per
    refused `cuts` entry, {"record": the query sequence's index in the answer,
    "seq_description", "graph", "graph_path", "annotation_path", and the entry's `refused`
    members: "account", "code", "max_memory_mb", and with the answer account "priced_mb"}, in
    the answer's order. A refused unit never means that nothing matched: its rows were not
    sent."""
    refused = []
    for i, result in enumerate(json_res or []):
        for cut in result.get('cuts') or []:
            if isinstance(cut, dict) and isinstance(cut.get('refused'), dict):
                refused.append({'record': i,
                                'seq_description': result.get('seq_description'),
                                'graph': cut.get('graph'),
                                'graph_path': cut.get('graph_path'),
                                'annotation_path': cut.get('annotation_path'),
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
