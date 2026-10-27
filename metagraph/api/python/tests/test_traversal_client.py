from unittest.mock import Mock, patch

import pytest

from metagraph.client import GraphApiError, GraphClient, GraphClientJson


@pytest.mark.parametrize('client_type', [GraphClientJson, GraphClient])
@pytest.mark.parametrize('labels', ['sample1', ['sample1']])
def test_resolve_preserves_label_names(client_type, labels):
    response = Mock(ok=True, status_code=200)
    response.json.return_value = {}
    with patch('metagraph.client.requests.post', return_value=response) as post:
        client_type('localhost', 5555).resolve('ACGTA', labels=labels)
    assert post.call_args[1]['json']['labels'] == ['sample1']


@pytest.mark.parametrize('client_type', [GraphClientJson, GraphClient])
def test_non_json_error_preserves_status_and_body(client_type):
    response = Mock(ok=False, status_code=503, text='temporarily unavailable')
    response.json.side_effect = ValueError('not JSON')
    with patch('metagraph.client.requests.post', return_value=response):
        with pytest.raises(GraphApiError) as error:
            client_type('localhost', 5555).traverse({})
    assert error.value.status_code == 503
    assert error.value.response == 'temporarily unavailable'


@pytest.mark.parametrize('client_type', [GraphClientJson, GraphClient])
def test_attempt_id_preserves_server_accepted_characters(client_type):
    response = Mock(ok=True, status_code=200)
    response.json.return_value = {}
    with patch('metagraph.client.requests.get', return_value=response) as get:
        client_type('localhost', 5555).traversal_attempt('run:1._-')
    assert get.call_args[1]['url'].endswith('/traverse/attempt/run:1._-')
