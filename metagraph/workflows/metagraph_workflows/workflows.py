import importlib
import logging
import sys
from pathlib import Path
from typing import Iterable, Optional, Union, Dict, Any
import shlex

import snakemake

from .common import AnnotationLabelsSource, AnnotationFormats
from .constants import SEQS_FILE_LIST_PATH, SEQS_DIR_PATH

WORKFLOW_ROOT = Path(__file__).parent / 'snakemake'

LOGGING_FORMAT='%(asctime)s - %(levelname)s: %(message)s'

logging.basicConfig(format=LOGGING_FORMAT, level=logging.WARNING)


default_path = Path(WORKFLOW_ROOT / 'default.yml')

# TODO: use custom config object? fluent config?
def run_build_workflow(
        output_dir: Path,
        seqs_file_list_path: Optional[Path] = None,
        seqs_dir_path: Optional[Path] = None,
        k: Optional[int] = None,
        base_name: Optional[str] = None,
        build_primary_graph: bool = False,
        annotation_formats: Iterable[AnnotationFormats] = (),
        annotation_labels_source: Optional[AnnotationLabelsSource] = None,
        metagraph_cmd: Optional[str] = None,
        threads: Optional[int] = None,
        force: bool = False,
        verbose: bool = False,
        dryrun: bool = False,
        additional_snakemake_args: Optional[Dict[str, Any]] = None
) -> bool:
    # TODO: support str argumt?

    snakefile_path = Path(WORKFLOW_ROOT / 'Snakefile')

    config = snakemake.load_configfile(default_path)

    if not seqs_file_list_path and not seqs_dir_path:
        raise ValueError("seqs_file_list_path and seqs_dir_path cannot both be None")

    if seqs_file_list_path:
        config[SEQS_FILE_LIST_PATH] = str(seqs_file_list_path)
    if seqs_dir_path:
        config[SEQS_DIR_PATH] = str(seqs_dir_path)

    config['output_directory'] = str(output_dir)

    config['k'] = k if k else config['k']

    if annotation_labels_source:
        config['annotation_labels_source'] = annotation_labels_source.value

    config['base_name'] = base_name if base_name else config['base_name']
    config['build_primary_graph'] = build_primary_graph

    config['annotation_formats'] = [af.value for af in
                                    annotation_formats] if annotation_formats else config['annotation_formats']

    config['metagraph_cmd'] = metagraph_cmd if metagraph_cmd else config['metagraph_cmd']
    config['max_threads'] = threads if threads else snakemake.available_cpu_count()

    if verbose:
        importlib.reload(logging)
        logging.basicConfig(format=LOGGING_FORMAT, level=logging.INFO)
        logging.info("Dumping config:")
        for k, v in sorted(config.items(), key=lambda t: t[0]):
            logging.info(f"\t{k}: {v}")

    additional_args = additional_snakemake_args if additional_snakemake_args else {}

    was_successful = snakemake.snakemake(str(snakefile_path), config=config,
                                         scheduler='greedy',
                                         forceall=force,
                                         dryrun=dryrun,
                                         **additional_args
                                         )

    return was_successful


def setup_build_parser(parser):
    parser.add_argument('output_dir', type=Path)

    input_seq_group = parser.add_argument_group('input sequence paths', '')

    input_seq_group_xor = input_seq_group.add_mutually_exclusive_group(required=True)
    input_seq_group_xor.add_argument('--seqs-file-list-path',
                                     help='Path to text file containing paths of sequences files')
    input_seq_group_xor.add_argument('--seqs-dir-path',
                                     help="Path to directory containing sequence files")

    graph = parser.add_argument_group('graph', 'arguments for graph building')
    graph.add_argument('-k', type=int, default=None)
    graph.add_argument('--base-name', default=None)
    graph.add_argument('--build-primary-graph', default=False,
                       action='store_true')

    annotation = parser.add_argument_group('annotation',
                                           'arguments for annotations')
    annotation.add_argument('--annotation-format', action='append',
                            default=[],
                            help=f"Annotation format (can be used multiple times). "
                                 f"Possible values: {', '.join([v.value for v in AnnotationFormats])}")
    annotation.add_argument('--annotation-labels-source',
                            type=AnnotationLabelsSource,
                            default=AnnotationLabelsSource.SEQUENCE_HEADERS,
                            help=f"What should be used as column labels. Possible values: "
                                 f"{', '.join([v.value for v in AnnotationLabelsSource])}")

    workflow = parser.add_argument_group('workflow',
                                         'arguments for the workflow')
    workflow.add_argument('--threads', type=int, default=None)
    workflow.add_argument('--force', default=False, action='store_true')
    workflow.add_argument('--verbose', default=False, action='store_true')
    workflow.add_argument('--dryrun', default=False, action='store_true')
    workflow.add_argument('--metagraph-cmd', type=str, default=None)
    workflow.add_argument('--additional-snakemake-args', type=str, default='',
                          help='Additional arguments to pass to snakemake, e.g. --additional-snakemake-args="arg1=val1 arg2=val2"')

    parser.set_defaults(func=init_build)


def _convert_type(v: str) -> Any:
    if v.lower() == 'true' or v == '1':
        return True
    elif v.lower() == 'false' or v == '0':
        return False

    try:
        return float(v)
    except:
        pass

    try:
        return int(v)
    except:
        pass

    return v


def _parse_additional_snakemake_args(arg: str) -> Dict[str, Any]:
    ret = {}
    for a in shlex.split(arg):
        if '=' not in a:
            raise ValueError("ex")

        k, v = a.split('=')
        ret[k] = _convert_type(v)

    return ret


def init_build(args):
    was_successful = run_build_workflow(
        args.output_dir,
        seqs_file_list_path=args.seqs_file_list_path,
        seqs_dir_path=args.seqs_dir_path,
        k=args.k,
        base_name=args.base_name,
        build_primary_graph=args.build_primary_graph,
        annotation_formats=[AnnotationFormats(af) for af in args.annotation_format],
        annotation_labels_source=args.annotation_labels_source,
        metagraph_cmd=args.metagraph_cmd,
        threads=args.threads,
        force=args.force,
        verbose=args.verbose,
        dryrun=args.dryrun,
        additional_snakemake_args=_parse_additional_snakemake_args(args.additional_snakemake_args)
    )

    if not was_successful:
        logging.error("Workflow did not run successfully")
        sys.exit(1)
