import os
import subprocess
import shutil
import re
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
COORD_FORMATS = {
    AnnotationFormats.BRWT_COORD,
    AnnotationFormats.ROW_DIFF_COORD,
    AnnotationFormats.ROW_DIFF_BRWT_COORD,
    AnnotationFormats.ROW_DIFF_DISK_COORD,
}


WORKFLOW_ROOT = Path(metagraph_workflows.__file__).parent / 'snakemake'


def _resolve_metagraph_cmd() -> str | None:
    """Best-effort resolution of the local `metagraph` binary for tests.

    Tests should normally run with a `metagraph` executable available on `PATH`.
    When that's not the case (e.g. fresh dev environments), we fall back to
    `../build/metagraph` if it exists.
    """
    if shutil.which("metagraph") is not None:
        return None

    code_base = Path(os.path.realpath(__file__)).parent.parent  # metagraph/workflows
    repo_root = code_base.parent  # metagraph/
    candidate = repo_root / "build" / "metagraph"
    if candidate.exists() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def run_wrapper(args_list):
    code_base = Path(os.path.realpath(__file__)).parent.parent
    normalized_args = list(args_list)
    if normalized_args and normalized_args[0] == "build":
        has_output_flag = any(arg in ("-o", "--output_dir") for arg in normalized_args)
        if not has_output_flag and len(normalized_args) > 1:
            last = normalized_args[-1]
            last_str = str(last)
            if not last_str.startswith("-"):
                normalized_args = normalized_args[:-1] + ["--output_dir", last]

    process_args = ['python', '-m', 'metagraph_workflows.cli'] + normalized_args

    # If tests are running without `metagraph` on PATH, inject `--metagraph-cmd`
    # pointing to the locally built binary (when available).
    if "--metagraph-cmd" not in normalized_args:
        metagraph_cmd = _resolve_metagraph_cmd()
        if metagraph_cmd is None:
            pytest.skip("metagraph executable not found in PATH and local build/metagraph missing")
        if not normalized_args:
            pytest.skip("empty args_list passed to run_wrapper")
        process_args = ['python', '-m', 'metagraph_workflows.cli'] + normalized_args + [
            "--metagraph-cmd", metagraph_cmd
        ]

    proc = subprocess.run(
        [str(a) for a in process_args],
        cwd=code_base,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

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
    list(product([False, True], AnnotationFormats, [AnnotationLabelsSource.FILE_NAMES])))
def test_build_workflow(primary, annotation_format, annotation_label_src, sample_list_path, output_dir):

    base_args = ['build',
                 '--seqs-file-list-path', sample_list_path,
                 '-k', 5,
                 '--annotation-format', annotation_format.value,
                 '--anno-source', annotation_label_src.value]
    if annotation_format in COUNT_FORMATS:
        base_args += ['--with-counts']
    if annotation_format in COORD_FORMATS:
        base_args += ['--with-coords']

    base_args += ['--primary'] if primary else []

    ret = run_wrapper(base_args + [output_dir])

    if ret.returncode != 0:
        print("Workflow test was not successful:")
        print(ret.stdout.decode())

    assert ret.returncode == 0, ret.stderr

    assert len(output_dir.listdir()) > 1


def test_workflow_invocation_via_python(sample_list_path, output_dir):
    metagraph_cmd = _resolve_metagraph_cmd()
    if metagraph_cmd is None and shutil.which("metagraph") is None:
        pytest.skip("metagraph executable not found in PATH and local build/metagraph missing")

    assert cli.run_build_workflow(
        output_dir,
        seqs_file_list_path=sample_list_path,
        metagraph_cmd=metagraph_cmd,
    ) is None


def test_workflow_invocation_additional_args(sample_list_path, output_dir):
    base_args = ['build',
                 '--seqs-file-list-path', sample_list_path,
                 '-k', 5,
                 '--extra-args="summary=True"']

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


def test_with_coordinates_defaults_to_row_diff_brwt_coord(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-coords',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "with_coordinates: true" in cfg
    assert "annotation_formats:" in cfg
    assert "- row_diff_brwt_coord" in cfg


def test_with_coordinates_rejects_incompatible_annotation_format(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-coords',
        '--annotation-format', AnnotationFormats.BRWT.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "Coordinate-aware mode is enabled" in proc.stdout.decode()


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


def test_coord_capable_format_auto_enables_with_coordinates(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_BRWT_COORD.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "with_coordinates: true" in cfg
    assert "- row_diff_brwt_coord" in cfg


def test_with_counts_and_with_coordinates_are_mutually_exclusive(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--with-counts',
        '--with-coords',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "mutually exclusive" in proc.stdout.decode()


def test_mixed_count_and_coord_formats_are_mutually_exclusive(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--annotation-format', AnnotationFormats.ROW_DIFF_BRWT_COORD.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "mutually exclusive" in proc.stdout.decode()


def test_build_help_mentions_defaults():
    proc = run_wrapper(['build', '-h'])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    # Argparse may insert newlines + indentation into long help strings;
    # normalize all whitespace so substring checks are stable.
    out_norm = re.sub(r'\s+', ' ', out).strip()
    assert "[relax.row_diff_brwt/row_diff_int_brwt/row_diff_brwt_coord]" in out_norm
    assert "row_diff_int_brwt" in out
    assert "row_diff_brwt_coord" in out


def test_dryrun_prints_summary(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "=== Workflow summary ===" in out
    assert "Mode: dry-run" in out


def test_missing_metagraph_executable_fails_fast(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--metagraph-cmd', 'definitely_missing_metagraph_binary_12345',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "was not found in PATH" in proc.stdout.decode()


def test_invalid_annotation_format_shows_suggestion(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', 'row_diff_int_brwt1',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    out = proc.stdout.decode()
    assert "Unsupported annotation format 'row_diff_int_brwt1'" in out
    assert "Did you mean 'row_diff_int_brwt'" in out


def test_invalid_coord_annotation_format_shows_suggestion(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', 'row_diff_brwt_coord1',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    out = proc.stdout.decode()
    assert "Unsupported annotation format 'row_diff_brwt_coord1'" in out
    assert "Did you mean 'row_diff_brwt_coord'" in out


@pytest.mark.parametrize("count_width", [2, 12, 32])
def test_count_width_is_propagated_to_count_build_steps(sample_list_path, output_dir, count_width):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotation-format', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--count-width', str(count_width),
        '--dryrun',
        '--extra-args=printshellcmds=True',
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


def test_annotate_threads_each_default_is_eight(sample_list_path, output_dir):
    # 16 threads / threads_each=8 -> parallel_cols=2, effective_each=8.
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--threads', '16',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "annotate_threads_each: 8" in cfg
    out = proc.stdout.decode()
    assert "--parallel 2" in out
    assert "--threads-each 8" in out


def test_annotate_threads_each_overrides_default(sample_list_path, output_dir):
    # 16 threads / threads_each=4 -> parallel_cols=4, effective_each=4.
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--threads', '16',
        '--annotate-threads-each', '4',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "annotate_threads_each: 4" in cfg
    out = proc.stdout.decode()
    assert "--parallel 4" in out
    assert "--threads-each 4" in out


def test_annotate_threads_each_redistributes_leftover(sample_list_path, output_dir):
    # 12 threads / threads_each=8 -> parallel_cols=ceil(12/8)=2,
    # effective_each=ceil(12/2)=6, total used = 12 (no waste).
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--threads', '12',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "--parallel 2" in out
    assert "--threads-each 6" in out


def test_annotate_threads_each_ceiling_overcommits_at_boundary(sample_list_path, output_dir):
    # 13 threads / threads_each=8 -> parallel_cols=ceil(13/8)=2,
    # effective_each=ceil(13/2)=7, total=14 (1-thread overcommit).
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--threads', '13',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "--parallel 2" in out
    assert "--threads-each 7" in out


def test_annotate_threads_each_must_be_positive(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        '--seqs-file-list-path', sample_list_path,
        '--annotate-threads-each', '0',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--annotate-threads-each must be >= 1" in proc.stdout.decode()
