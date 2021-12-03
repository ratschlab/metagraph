
from metagraph_workflows import cli


def test_parse_additional_snakemake_args():
    assert cli._parse_additional_snakemake_args('') == {}
    assert cli._parse_additional_snakemake_args('dryrun=1') == {
        'dryrun': True}

    assert cli._parse_additional_snakemake_args(
        'some_param="hello world" another=123') == {
               'some_param': 'hello world', 'another': 123}
