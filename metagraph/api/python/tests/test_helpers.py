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
