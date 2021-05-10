import pytest

from metagraph_workflows import rule_utils

@pytest.mark.parametrize("case, expected",
    [
        ('/my/path/sample.fasta', 'sample'),
        ('/my/path/sample.fasta.gz', 'sample'),
        ('/my/path/sample.txt', 'sample'),
        ('/my/path/sample', 'sample'),
        ('/my/path/sample/', 'sample'),
    ]
)
def test_get_sample_name(case, expected):
    assert rule_utils.get_sample_name(case) == expected