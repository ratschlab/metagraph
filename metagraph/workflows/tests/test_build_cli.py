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
        if "-o" not in normalized_args and len(normalized_args) > 1:
            last = normalized_args[-1]
            last_str = str(last)
            if not last_str.startswith("-"):
                normalized_args = normalized_args[:-1] + ["-o", last]

    process_args = ['python', '-m', 'metagraph_workflows.cli'] + normalized_args

    # If tests are running without `metagraph` on PATH, inject `--metagraph-cmd`
    # pointing to the locally built binary (when available).
    if "--metagraph-cmd" not in normalized_args and shutil.which("metagraph") is None:
        metagraph_cmd = _resolve_metagraph_cmd()
        if metagraph_cmd is None:
            pytest.skip("metagraph executable not found in PATH and local build/metagraph missing")
        process_args = process_args + ["--metagraph-cmd", metagraph_cmd]

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


@pytest.mark.parametrize('primary,annotation_format,annotation_label_src', list(product([False], [AnnotationFormats.ROW_DIFF_BRWT], [AnnotationLabelsSource.HEADER])) +
    list(product([False, True], AnnotationFormats, [AnnotationLabelsSource.FILENAME])))
def test_build_workflow(primary, annotation_format, annotation_label_src, sample_list_path, output_dir):

    base_args = ['build',
                 sample_list_path,
                 '-k', 5,
                 '--anno-type', annotation_format.value,
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

    assert cli.run_workflow(
        output_dir,
        samples=sample_list_path,
        metagraph_cmd=metagraph_cmd,
    ) is None


def test_workflow_invocation_additional_args(sample_list_path, output_dir):
    base_args = ['build',
                 sample_list_path,
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
        sample_list_path,
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
        sample_list_path,
        '--with-counts',
        '--anno-type', AnnotationFormats.BRWT.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "Count-aware mode is enabled" in proc.stdout.decode()


def test_with_coordinates_defaults_to_row_diff_brwt_coord(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
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
        sample_list_path,
        '--with-coords',
        '--anno-type', AnnotationFormats.BRWT.value,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "Coordinate-aware mode is enabled" in proc.stdout.decode()


def test_with_counts_respects_explicit_annotation_format(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--with-counts',
        '--anno-type', AnnotationFormats.INT_BRWT.value,
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
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_INT_DISK.value,
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
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_BRWT_COORD.value,
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
        sample_list_path,
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
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--anno-type', AnnotationFormats.ROW_DIFF_BRWT_COORD.value,
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
    # Default formats are rendered bracketed in the list, e.g. `[relax.row_diff_brwt]`.
    assert "[relax.row_diff_brwt]" in out_norm
    assert "[row_diff_int_brwt]" in out_norm
    assert "[row_diff_brwt_coord]" in out_norm


def test_dryrun_prints_summary(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
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
        sample_list_path,
        '--metagraph-cmd', 'definitely_missing_metagraph_binary_12345',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "was not found in PATH" in proc.stdout.decode()


def test_invalid_annotation_format_shows_suggestion(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--anno-type', 'row_diff_int_brwt1',
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
        sample_list_path,
        '--anno-type', 'row_diff_brwt_coord1',
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
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
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
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_INT_BRWT.value,
        '--count-width', str(invalid_count_width),
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--count-width must be in range [2, 32]" in proc.stdout.decode()


def test_annotate_threads_each_default_is_eight(sample_list_path, output_dir):
    # 16 threads / threads_each=8 -> parallel_cols=2, effective_each=8.
    # threads_each defaults to 8 for binary mode (mode-derived).
    proc = run_wrapper([
        'build',
        sample_list_path,
        '-p', '16',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "-p 2" in out
    assert "--threads-each 8" in out


def test_annotate_threads_each_overrides_default(sample_list_path, output_dir):
    # 16 threads / threads_each=4 -> parallel_cols=4, effective_each=4.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '-p', '16',
        '--anno-threads-each', '4',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "annotate_threads_each: 4" in cfg
    out = proc.stdout.decode()
    assert "-p 4" in out
    assert "--threads-each 4" in out


def test_annotate_threads_each_redistributes_leftover(sample_list_path, output_dir):
    # 12 threads / threads_each=8 -> parallel_cols=ceil(12/8)=2,
    # effective_each=ceil(12/2)=6, total used = 12 (no waste).
    proc = run_wrapper([
        'build',
        sample_list_path,
        '-p', '12',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "-p 2" in out
    assert "--threads-each 6" in out


def test_annotate_threads_each_ceiling_overcommits_at_boundary(sample_list_path, output_dir):
    # 13 threads / threads_each=8 -> parallel_cols=ceil(13/8)=2,
    # effective_each=ceil(13/2)=7, total=14 (1-thread overcommit).
    proc = run_wrapper([
        'build',
        sample_list_path,
        '-p', '13',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "-p 2" in out
    assert "--threads-each 7" in out


def test_disk_swap_dir_propagates_to_metagraph_stages(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--disk-swap-dir', '/var/tmp/test-swap',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    assert "tmpdir: /var/tmp/test-swap" in cfg
    out = proc.stdout.decode()
    # The dir is quoted (so empty/whitespace paths round-trip safely too).
    assert '--disk-swap "/var/tmp/test-swap"' in out


def test_metagraph_always_runs_with_v(sample_list_path, output_dir):
    # `-v` is always passed to metagraph so log files capture every
    # trace; the workflow's --verbose only controls whether that output
    # also streams to the terminal.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "metagraph build  -v" in out


@pytest.mark.parametrize("flag", ["-v", "--verbose"])
def test_verbose_streams_metagraph_output_to_terminal(sample_list_path, output_dir, flag):
    # With --verbose, tee's stdout goes through; otherwise it's silenced
    # (> /dev/null) so metagraph output lives in {log} only.
    proc = run_wrapper([
        'build',
        sample_list_path,
        flag,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    assert "tee" in out
    assert "/dev/null" not in out


def test_default_silences_metagraph_terminal_output(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    # Without --verbose, every tee redirects stdout to /dev/null so
    # the user's terminal only sees Snakemake job-status lines.
    for line in out.splitlines():
        if 'tee ' not in line:
            continue
        assert '> /dev/null' in line, f"non-verbose tee missing silencer: {line}"


def test_disk_swap_dir_empty_string_disables_swap(sample_list_path, output_dir):
    # `--disk-swap-dir ""` is the explicit-off sentinel: every metagraph
    # invocation gets `--disk-swap ""` so nothing spills to disk.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--disk-swap-dir', '',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    for line in out.splitlines():
        if '--disk-swap' not in line:
            continue
        assert '--disk-swap ""' in line, f"Unexpected --disk-swap target: {line}"
    cfg = (output_dir / "config.yaml").read()
    assert "tmpdir:" not in cfg


def test_disk_swap_dir_defaults_to_output_temp(sample_list_path, output_dir):
    # When --disk-swap-dir is omitted, the workflow defaults to
    # <output_dir>/temp so metagraph transform_anno doesn't fall back to
    # its OUT_BASEDIR default and silently spill next to artifacts.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    expected = f'--disk-swap "{output_dir}/temp"'
    for line in out.splitlines():
        if '--disk-swap' not in line:
            continue
        assert expected in line, f"Unexpected --disk-swap target: {line}"
    cfg = (output_dir / "config.yaml").read()
    assert f"tmpdir: {output_dir}/temp" in cfg


def test_small_graph_step_runs(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "rule build_small_graph" in out
    # The transform step writes <base>_small.dbg via `-o <base>_small`.
    assert "metagraph transform" in out
    assert "--state small" in out


@pytest.mark.parametrize("fmt", [
    "brwt_coord", "row_diff_coord", "row_diff_brwt_coord", "row_diff_disk_coord",
])
def test_index_header_coords_fires_in_coords_filenames_mode(sample_list_path, output_dir, fmt):
    # `.seqs` only makes sense for coord-aware annotations indexed with
    # --anno-filename: the sidecar maps file-level coord ranges back to
    # original sequence headers.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--anno-type', fmt,
        '--anno-source', 'filename',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    assert "rule index_header_coords" in out
    # The loader (load_annotated_graph.cpp) strips the annotation's full
    # kExtension (e.g. `.row_diff_brwt_coord.annodbg`) and appends `.seqs`,
    # so it always looks for `<graph>.seqs`. A per-format file would never
    # be picked up.
    assert "/graph.seqs" in out
    assert f"/graph.{fmt}.seqs" not in out
    # Column order comes from the final annotation (not the input file list),
    # so BRWT-reordered columns line up.
    assert "stats --print-col-names" in out
    assert "--index-header-coords" in out


def test_with_coords_alone_fires_seqs_sidecar(sample_list_path, output_dir):
    # The default --anno-source is `filename`, so `--with-coords` alone is
    # enough to trigger the .seqs sidecar rule.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--with-coords',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    assert "rule index_header_coords" in out
    assert "--index-header-coords" in out
    # Loader looks for <graph>.seqs, not per-format.
    assert "/graph.seqs" in out


@pytest.mark.parametrize("flags,reason", [
    (["--with-coords", "--anno-source", "header"], "header mode doesn't need .seqs"),
    (["--anno-source", "filename"], "no --with-coords -> binary mode"),
    (["--with-counts", "--anno-source", "filename"], "counts + filename doesn't need .seqs"),
])
def test_index_header_coords_does_not_fire_outside_coords_filenames(
        sample_list_path, output_dir, flags, reason):
    proc = run_wrapper([
        'build',
        sample_list_path,
        *flags,
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    # Match the rule line specifically -- pytest's tmpdir path may contain
    # the substring "index_header_coords" when the test name does.
    assert "rule index_header_coords" not in proc.stdout.decode(), reason
    assert "--index-header-coords" not in proc.stdout.decode(), reason


def test_coords_mode_auto_picks_threads_each_16(sample_list_path, output_dir):
    # Snakefile derives annotate_threads_each from mode when unset; for
    # coords it becomes 16, so parallel_cols=ceil(16/16)=1.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--anno-type', AnnotationFormats.ROW_DIFF_BRWT_COORD.value,
        '-p', '16',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "-p 1" in out
    assert "--threads-each 16" in out


@pytest.mark.parametrize("fmt", ["row_diff_flat", "row_diff_sparse", "row_diff_disk"])
def test_row_diff_binary_formats_reach_shell(sample_list_path, output_dir, fmt):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--anno-type', fmt,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    # Resolved by the shared annotate_row_diff_binary rule, which passes
    # --anno-type <fmt> and -i graph.dbg (row-diff transform input).
    assert f"--anno-type {fmt}" in out
    assert "annotate_row_diff_binary" in out
    if fmt == "row_diff_disk":
        # disk variant gets --mem-cap-gb; the other two omit it.
        assert "--mem-cap-gb" in out


def test_brwt_parallel_nodes_default_is_10(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    out = proc.stdout.decode()
    assert "--parallel-nodes 10" in out


def test_brwt_subsample_default_and_override(sample_list_path, output_dir, tmpdir):
    # Default: 1000000 from default.yml; reaches the row_diff_brwt rule.
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0
    assert "--subsample 1000000" in proc.stdout.decode()

    # CLI override.
    other_out = tmpdir / 'other_out'
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--brwt-subsample', '200000',
        '--dryrun',
        '--extra-args=printshellcmds=True',
        other_out,
    ])
    assert proc.returncode == 0
    assert "--subsample 200000" in proc.stdout.decode()


@pytest.mark.parametrize("expr,expanded", [("1e6", 1000000), ("2.5e4", 25000)])
def test_brwt_subsample_accepts_scientific_notation(sample_list_path, output_dir, expr, expanded):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--brwt-subsample', expr,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    assert f"--subsample {expanded}" in proc.stdout.decode()


@pytest.mark.parametrize("bad", [0, 1, 999])
def test_brwt_subsample_below_1000_rejected(sample_list_path, output_dir, bad):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--brwt-subsample', str(bad),
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--brwt-subsample must be >= 1000" in proc.stdout.decode()


def test_mem_gb_sets_memory_mb(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--mem-gb', '12',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode == 0
    cfg = (output_dir / "config.yaml").read()
    # 12 GiB -> 12 * 1024 = 12288 MB.
    assert "memory_mb: 12288" in cfg


def test_mem_gb_must_be_positive(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--mem-gb', '0',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--mem-gb must be > 0" in proc.stdout.decode()


@pytest.fixture
def stub_graph_path(tmpdir):
    p = tmpdir / "graph_in.dbg"
    p.write("")  # zero-byte placeholder is enough for snakemake's existence check
    return p


def test_build_with_graph_skips_build_rules(sample_list_path, stub_graph_path, output_dir):
    proc = run_wrapper([
        'build',
        '--graph', stub_graph_path,
        sample_list_path,
        '--dryrun',
        '--extra-args=printshellcmds=True',
        '-o', output_dir,
    ])
    assert proc.returncode == 0, proc.stdout.decode()
    out = proc.stdout.decode()
    # The build pipeline must not appear when --graph is provided.
    for build_rule in (
        "build_joint_graph",
        "build_joint_primary",
        "build_canonical_graph_single_sample",
        "primarize_joint_graph",
        "primarize_canonical_graph_single_sample",
        "extract_kmer_counts",
    ):
        assert build_rule not in out, f"unexpected build rule in DAG: {build_rule}"
    # Annotate + transforms still run.
    for kept_rule in (
        "rule annotate:",
        "rule transform_rd_stage0",
        "rule transform_rd_stage1",
        "rule transform_rd_stage2",
        "rule annotate_row_diff_brwt",
    ):
        assert kept_rule in out, f"missing rule in DAG: {kept_rule}"
    # Small-state graph is intentionally skipped when --graph is provided.
    assert "rule build_small_graph" not in out

    cfg = (output_dir / "config.yaml").read()
    assert "external_graph: true" in cfg
    # The user's graph must be symlinked into the output dir.
    target = output_dir / "graph.dbg"
    assert target.exists()


def test_build_with_graph_requires_existing_graph(sample_list_path, output_dir, tmpdir):
    missing = tmpdir / "does_not_exist.dbg"
    proc = run_wrapper([
        'build',
        '--graph', missing,
        sample_list_path,
        '--dryrun',
        '-o', output_dir,
    ])
    assert proc.returncode != 0
    assert "Graph file not found" in proc.stdout.decode()


def test_annotate_threads_each_must_be_positive(sample_list_path, output_dir):
    proc = run_wrapper([
        'build',
        sample_list_path,
        '--anno-threads-each', '0',
        '--dryrun',
        output_dir,
    ])
    assert proc.returncode != 0
    assert "--anno-threads-each must be >= 1" in proc.stdout.decode()
