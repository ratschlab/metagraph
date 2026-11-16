import os
from pathlib import Path
import json
from metagraph import helpers
import pytest


def _load_json_data(filename):
    path = Path(os.path.abspath(__file__)).parent / 'data' / filename

    with open(path, 'r') as f:
        return json.load(f)


@pytest.mark.parametrize("file_name,align,expected_shape", [
    ('search_response.json', False, (4, 15)),
    ('search_with_align_response.json', True, (354, 15))
])
def test_df_from_search_result(file_name, align, expected_shape):
    json_obj = _load_json_data(file_name)
    df = helpers.df_from_search_result(json_obj)

    assert df.shape == expected_shape

    expected_cols = ['kmer_count', 'sample', 'city', 'city_latitude', 'city_longitude',
                          'city_total_population', 'continent', 'latitude', 'longitude',
                          'metasub_name', 'num_reads', 'sample_type', 'station',
                          'surface_material', 'seq_description']

    assert list(df.columns) == expected_cols


@pytest.mark.parametrize("file_name,expected_shape", [
    ('align_response.json', (7, 4)),
    ('empty_align_response.json', (0, 4))
])
def test_df_from_align_result(file_name, expected_shape):
    json_obj = _load_json_data(file_name)
    df = helpers.df_from_align_result(json_obj)

    assert df.shape == expected_shape
    assert list(df.columns) == ['cigar', 'score', 'sequence', 'seq_description']


# Lane T7-P1 (test 24): a /search answer with units refused alone (the request's
# partial_refusal, SPEC-multigraph-envelope.md section 3): the frame holds the rows sent, a
# warning names each refused unit and the frame's attrs["refused"] lists them; an answer
# without one warns nothing
def _cut(graph, returned, refused=None):
    cut = {'annotation_path': f'/idx/{graph}.annodbg', 'complete': refused is None,
           'graph': graph, 'graph_path': f'/idx/{graph}.dbg', 'index_fp': None,
           'nth_score': None, 'returned': returned, 'tied': 0, 'tied_returned': 0}
    if refused is not None:
        cut['refused'] = refused
    return cut


def test_df_from_search_result_lists_the_refused_units():
    import warnings as w
    refused = {'account': 'values', 'code': 'search_answer_too_large', 'max_memory_mb': 4096}
    answer = [
        {'seq_description': 'fwd',
         'results': [{'sample': 'a1', 'kmer_count': 3, 'graph': 'a',
                      'graph_path': '/idx/a.dbg', 'annotation_path': '/idx/a.annodbg'}],
         'cuts': [_cut('a', 1), _cut('b', 0, refused)]},
        {'seq_description': 'rc',
         'results': [{'sample': 'b1', 'kmer_count': 2, 'graph': 'b',
                      'graph_path': '/idx/b.dbg', 'annotation_path': '/idx/b.annodbg'}],
         'cuts': [_cut('a', 0, dict(refused, account='answer', priced_mb=5000)), _cut('b', 1)]},
    ]
    with pytest.warns(helpers.SearchRefusedWarning) as caught:
        df = helpers.df_from_search_result(answer)
    assert df.shape[0] == 2
    assert list(df['sample']) == ['a1', 'b1']
    assert df.attrs['refused'] == [
        {'record': 0, 'seq_description': 'fwd', 'graph': 'b', 'graph_path': '/idx/b.dbg',
         'annotation_path': '/idx/b.annodbg', 'account': 'values',
         'code': 'search_answer_too_large', 'max_memory_mb': 4096},
        {'record': 1, 'seq_description': 'rc', 'graph': 'a', 'graph_path': '/idx/a.dbg',
         'annotation_path': '/idx/a.annodbg', 'account': 'answer',
         'code': 'search_answer_too_large', 'max_memory_mb': 4096, 'priced_mb': 5000}]
    message = str(caught[0].message)
    assert "record 0 ('fwd') on graph 'b' (values account)" in message
    assert "record 1 ('rc') on graph 'a' (answer account)" in message
    # a single-graph server's entry has null tags
    single = [{'seq_description': 'q', 'results': [],
               'cuts': [dict(_cut('x', 0, refused), graph=None, graph_path=None,
                             annotation_path=None)]}]
    with pytest.warns(helpers.SearchRefusedWarning, match=r"record 0 \('q'\) on the graph"):
        df = helpers.df_from_search_result(single)
    assert df.shape[0] == 0 and df.attrs['refused'][0]['graph'] is None
    # nothing refused: no warning, no attrs
    answer[0]['cuts'][1] = _cut('b', 0)
    answer[1]['cuts'][0] = _cut('a', 0)
    with w.catch_warnings():
        w.simplefilter('error')
        df = helpers.df_from_search_result(answer)
    assert 'refused' not in df.attrs
    assert helpers.refused_units(answer) == []
    # an answer of an older server (no cuts) too
    assert helpers.refused_units([{'seq_description': 'q', 'results': []}]) == []


def test_the_client_sends_partial_refusal_only_when_asked():
    from unittest import mock
    from metagraph.client import GraphClientJson, GraphClient
    sent = []

    class Response:
        ok = True
        status_code = 200

        def json(self):
            return [{'seq_description': '0', 'results': [], 'cuts': [
                _cut('g', 0, {'account': 'values', 'code': 'search_answer_too_large',
                              'max_memory_mb': 1})]}]

    def post(url, json, timeout):
        sent.append(json)
        return Response()

    with mock.patch('requests.post', post):
        GraphClientJson('localhost', 1).search('ACGT')
        GraphClientJson('localhost', 1).search('ACGT', partial_refusal=True)
        with pytest.warns(helpers.SearchRefusedWarning):
            df = GraphClient('localhost', 1).search('ACGT', partial_refusal=True)
    assert 'partial_refusal' not in sent[0]
    assert sent[1]['partial_refusal'] is True and sent[2]['partial_refusal'] is True
    assert df.attrs['refused'][0]['account'] == 'values'
