import pytest
import sys
import subprocess

import metagraph_workflows.utils

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
    assert metagraph_workflows.utils.get_sample_name(case) == expected


def test_get_time_wrapper_command_uses_python_wrapper():
    cmd = metagraph_workflows.utils.get_time_wrapper_command({})
    assert "metagraph_workflows.time_wrapper" in cmd


def test_time_wrapper_success_prints_timing():
    proc = subprocess.run(
        [sys.executable, "-m", "metagraph_workflows.time_wrapper", "--", sys.executable, "-c", "print('ok')"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert proc.returncode == 0
    assert "ok" in proc.stdout
    assert "[timing] wall_sec=" in proc.stderr


def test_time_wrapper_missing_command_returns_127():
    proc = subprocess.run(
        [sys.executable, "-m", "metagraph_workflows.time_wrapper", "--", "definitely_missing_binary_12345"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert proc.returncode == 127
    assert "[timing] failed_to_exec=" in proc.stderr