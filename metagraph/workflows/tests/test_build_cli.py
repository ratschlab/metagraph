import os
import subprocess
from itertools import product
from pathlib import Path

import pytest

from metagraph_workflows import utils
from metagraph_workflows.common import AnnotationFormats, AnnotationLabelsSource

WORKFLOW_ROOT = Path(__file__).parent.parent / 'workflows'


def run_wrapper(args_list):
    code_base = Path(os.path.realpath(__file__)).parent.parent

    process_args = ['python', '-m', 'metagraph_workflows.cli'] + args_list

    proc = subprocess.run([str(a) for a in process_args],
                          cwd=code_base, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    return proc


@pytest.mark.parametrize('primary,annotation_format,annotation_label_src', list(product([False], [AnnotationFormats.ROW_DIFF_BRWT], [AnnotationLabelsSource.SEQUENCE_HEADERS])) +
    list(product([False, True], AnnotationFormats, [AnnotationLabelsSource.SEQUENCE_FILE_NAMES])))
def test_build_workflow(primary, annotation_format, annotation_label_src, tmpdir):
    sample_list_path = tmpdir / 'transcript_paths.txt'
    utils.create_transcript_path_list(WORKFLOW_ROOT / 'example_data', sample_list_path)

    output_dir = tmpdir

    base_args = ['build',
                 '--seqs-file-list-path', sample_list_path,
                 '-k', 5,
                 '--annotation-format', annotation_format.value,
                 '--annotation-labels-source', annotation_label_src.value]

    base_args += ['--build-primary-graph'] if primary else []

    ret = run_wrapper(base_args + [output_dir])

    if ret.returncode != 0:
        print("Workflow test was not successful:")
        print(ret.stdout.decode())

    assert ret.returncode == 0, ret.stderr

    assert len(output_dir.listdir()) > 1


