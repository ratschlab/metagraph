import os
import subprocess
from itertools import product
from pathlib import Path

import pytest

import metagraph_workflows
from metagraph_workflows import cli, utils
from metagraph_workflows.workflow_configs import AnnotationLabelsSource, \
    AnnotationFormats
COUNT_FORMATS = {
    AnnotationFormats.INT_BRWT,
    AnnotationFormats.ROW_DIFF_INT_BRWT,
    AnnotationFormats.ROW_DIFF_INT_DISK,
}


WORKFLOW_ROOT = Path(metagraph_workflows.__file__).parent / 'snakemake'


def run_wrapper(args_list):
    code_base = Path(os.path.realpath(__file__)).parent.parent

    process_args = ['python', '-m', 'metagraph_workflows.cli'] + args_list

    proc = subprocess.run([str(a) for a in process_args],
                          cwd=code_base, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    return proc


@pytest.fixture
def output_dir(tmpdir):
    return tmpdir / 'output'


@pytest.fixture
def sample_list_path(tmpdir):
    list_path = tmpdir / 'transcript_paths.txt'
    utils.create_transcript_path_list(WORKFLOW_ROOT / 'test_data', list_path)
    return list_path


@pytest.mark.parametrize('primary,annotation_format,annotation_label_src', list(product([False], [AnnotationFormats.ROW_DIFF_BRWT], [AnnotationLabelsSource.SEQUENCE_HEADERS])) +
    list(product([False, True], AnnotationFormats, [AnnotationLabelsSource.SEQUENCE_FILE_NAMES])))
def test_build_workflow(primary, annotation_format, annotation_label_src, sample_list_path, output_dir):

    base_args = ['build',
                 '--seqs-file-list-path', sample_list_path,
                 '-k', 5,
                 '--annotation-format', annotation_format.value,
                 '--annotation-labels-source', annotation_label_src.value]
    if annotation_format in COUNT_FORMATS:
        base_args += ['--with-counts']

    base_args += ['--build-primary-graph'] if primary else []

    ret = run_wrapper(base_args + [output_dir])

    if ret.returncode != 0:
        print("Workflow test was not successful:")
        print(ret.stdout.decode())

    assert ret.returncode == 0, ret.stderr

    assert len(output_dir.listdir()) > 1


def test_workflow_invocation_via_python(sample_list_path, output_dir):
    assert cli.run_build_workflow(output_dir, seqs_file_list_path=sample_list_path) is None


def test_workflow_invocation_additional_args(sample_list_path, output_dir):
    base_args = ['build',
                 '--seqs-file-list-path', sample_list_path,
                 '-k', 5,
                 '--additional-snakemake-args="summary=True"']

    proc = run_wrapper(base_args + [output_dir])

    assert proc.returncode == 0
    assert output_dir.exists()
    assert (output_dir / "config.yaml").exists()
    assert len([f for f in output_dir.listdir() if f.check(file=1)]) == 1


def test_with_counts_defaults_to_row_diff_int_brwt(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-counts',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "with_counts: true" in cfg
    assert "annotation_formats:" in cfg
    assert "- row_diff_int_brwt" in cfg


def test_with_counts_rejects_incompatible_annotation_format(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-counts',
        '--annotation-format', AnnotationFormats.BRWT.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "Count-aware mode is enabled" in proc.stdout.decode()


def test_with_counts_respects_explicit_annotation_format(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-counts',
        '--annotation-format', AnnotationFormats.INT_BRWT.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "- int_brwt" in cfg
    assert "- row_diff_int_brwt" not in cfg


def test_count_capable_format_auto_enables_with_counts(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_INT_DISK.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "with_counts: true" in cfg
    assert "- row_diff_int_disk" in cfg


def test_build_help_mentions_defaults():
    proc = run_wrapper(['build', '-h'])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "Default is relax.row_diff_brwt" in out
    assert "row_diff_int_brwt" in out


@pytest.mark.parametrize("count_width", [2, 12, 32])
def test_count_width_is_propagated_to_count_build_steps(sample_list_path, output_dir, count_width):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--count-width', str(count_width),
        '--dryrun',
        '--additional-snakemake-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0

    out = proc.stdout.decode()
    assert "--count-kmers" in out
    assert f"--count-width {count_width}" in out

    cfg = (output_dir / "config.yaml").read()
    assert f"count_width: {count_width}" in cfg
    assert "with_counts: true" in cfg


@pytest.mark.parametrize("invalid_count_width", [1, 64])
def test_count_width_out_of_range_fails(sample_list_path, output_dir, invalid_count_width):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--count-width', str(invalid_count_width),
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--count-width must be in range [2, 32]" in proc.stdout.decode()
