import subprocess
import os
from unittest.mock import Mock

import pytest

from metagraph_workflows import cli, utils


def test_expand_executable_home_directory(monkeypatch):
    run = Mock(return_value=subprocess.CompletedProcess(
        [], 0, '', '-k --kmer-length [INT] length [31]\n'))
    monkeypatch.setattr(utils.subprocess, 'run', run)
    config = {'metagraph_cmd': '~/bin/metagraph'}
    utils.resolve_build_k(config)
    assert run.call_args.args[0][0] == os.path.expanduser('~/bin/metagraph')
    assert config['k'] == 31


@pytest.mark.parametrize('k', [3, 31, 21, 15, 12])
@pytest.mark.parametrize('stream', ['stdout', 'stderr'])
def test_resolve_selected_binary_default(monkeypatch, k, stream):
    result = subprocess.CompletedProcess([], 0, stdout='', stderr='')
    setattr(result, stream, f'\t-k --kmer-length [INT] \tlength of the k-mer to use [{k}]\n')
    run = Mock(return_value=result)
    monkeypatch.setattr(utils.subprocess, 'run', run)
    config = {'k': None, 'metagraph_cmd': 'env OMP_NUM_THREADS=2 "/path with spaces/metagraph"'}
    utils.resolve_build_k(config)
    assert config['k'] == k
    assert run.call_args.args[0] == [
        'env', 'OMP_NUM_THREADS=2', '/path with spaces/metagraph', 'build', '--help']
    # Reusing the resolved config must not query the binary again.
    utils.resolve_build_k(config)
    assert run.call_count == 1


@pytest.mark.parametrize('failure', [
    FileNotFoundError(), subprocess.CalledProcessError(1, 'metagraph'),
    subprocess.TimeoutExpired('metagraph', 10),
])
def test_unavailable_default_requires_explicit_k(monkeypatch, failure):
    monkeypatch.setattr(utils.subprocess, 'run', Mock(side_effect=failure))
    config = {'metagraph_cmd': 'metagraph'}
    with pytest.raises(ValueError, match='set -k explicitly'):
        utils.resolve_build_k(config)
    assert 'k' not in config


@pytest.mark.parametrize('help_text', [
    '', 'other option [31]', '-k --kmer-length [INT] length [0]',
    '-k --kmer-length [INT] length [31]\n-k --kmer-length [INT] length [12]',
])
def test_unrecognized_help_has_no_silent_fallback(monkeypatch, help_text):
    monkeypatch.setattr(utils.subprocess, 'run', Mock(return_value=
                       subprocess.CompletedProcess([], 0, help_text, '')))
    with pytest.raises(ValueError, match='set -k explicitly'):
        utils.resolve_build_k({'metagraph_cmd': 'metagraph'})


@pytest.mark.parametrize('k', [0, -1, True, 1.5, '31'])
def test_invalid_explicit_k(k):
    with pytest.raises(ValueError, match='positive integer'):
        utils.resolve_build_k({'k': k, 'metagraph_cmd': 'unused'})


@pytest.mark.parametrize('mode', ['default', 'explicit', 'existing_graph'])
def test_workflow_resolves_only_when_build_needs_default(monkeypatch, tmp_path, mode):
    samples = tmp_path / 'samples'
    samples.mkdir()
    graph = tmp_path / 'existing.dbg'
    graph.touch()
    run = Mock(return_value=subprocess.CompletedProcess(
        [], 0, '', '-k --kmer-length [INT] length [12]\n'))
    monkeypatch.setattr(utils.subprocess, 'run', run)
    invoke = Mock()
    monkeypatch.setattr(cli, '_invoke_snakemake', invoke)
    cli.run_workflow(tmp_path / 'out', samples, threads=1, dryrun=True,
                     metagraph_cmd='chosen-metagraph',
                     k=7 if mode == 'explicit' else None,
                     graph=graph if mode == 'existing_graph' else None)
    config = invoke.call_args.args[0]
    assert config['k'] == {'default': 12, 'explicit': 7, 'existing_graph': None}[mode]
    assert run.call_count == (1 if mode == 'default' else 0)


def test_cli_omitted_k_reaches_workflow_as_none(monkeypatch, tmp_path):
    workflow = Mock()
    monkeypatch.setattr(cli, 'run_workflow', workflow)
    cli.main(['build', str(tmp_path), '-o', str(tmp_path / 'out')])
    assert workflow.call_args.kwargs['k'] is None
