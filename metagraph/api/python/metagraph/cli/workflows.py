import importlib
import logging
import sys
from pathlib import Path
from typing import Iterable, Optional

import snakemake

from metagraph.cli.common import AnnotationLabelsSource, AnnotationFormats

WORKFLOW_ROOT = Path(__file__).parent.parent.parent / 'workflows'

LOGGING_FORMAT='%(asctime)s - %(levelname)s: %(message)s'

logging.basicConfig(format=LOGGING_FORMAT, level=logging.WARNING)

# TODO: use custom config object? fluent config?
def run_build_workflow(
        sample_list_path: Path,
        output_dir: Path,
        k: Optional[int] = None,
        base_name: Optional[str] = None,
        build_primary_graph: bool = False,
        annotation_formats: Iterable[AnnotationFormats] = (),
        annotation_labels_source: Optional[AnnotationLabelsSource] = None,
        exec_cmd: str = None,
        threads: Optional[int] = None,
        force: bool = False,
        verbose: bool = False,
        dryrun: bool = False,
) -> bool:
    # TODO: support str argumt?

    snakefile_path = Path(WORKFLOW_ROOT / 'Snakefile')
    default_path = Path(WORKFLOW_ROOT / 'default.yml')

    config = snakemake.load_configfile(default_path)

    config['input_files_list_path'] = sample_list_path
    config['output_directory'] = output_dir

    config['k'] = k if k else config['k']

    config['annotation_labels_source'] = annotation_labels_source.value

    config['base_name'] = base_name if base_name else config['base_name']
    config['build_primary_graph'] = build_primary_graph

    config['annotation_formats'] = [af.value for af in
                                    annotation_formats] if annotation_formats else config['annotation_formats']

    config['exec_cmd'] = exec_cmd if exec_cmd else config['exec_cmd']
    config['max_threads'] = threads if threads else snakemake.available_cpu_count()

    config['verbose'] = verbose

    if config['verbose']:
        config['write_logs'] = True

        importlib.reload(logging)
        logging.basicConfig(format=LOGGING_FORMAT, level=logging.INFO)
        logging.info("Dumping config:")
        for k, v in sorted(config.items(), key=lambda t: t[0]):
            logging.info(f"\t{k}: {v}")

    was_successful = snakemake.snakemake(str(snakefile_path), config=config,
                                         scheduler='greedy',
                                         forceall=force,
                                         dryrun=dryrun
                                         )

    return was_successful


def setup_parser(parser):
    parser.add_argument('sample_list_path', type=Path)
    parser.add_argument('output_dir', type=Path)

    graph = parser.add_argument_group('graph', 'arguments for graph building')
    graph.add_argument('-k', type=int, default=None)
    graph.add_argument('--base-name', default=None)
    graph.add_argument('--build-primary-graph', default=False,
                       action='store_true')

    annotation = parser.add_argument_group('annotation',
                                           'arguments for annotations')
    annotation.add_argument('--annotation-format', action='append',
                            default=[])
    annotation.add_argument('--annotation-labels-source',
                            type=AnnotationLabelsSource,
                            default=None,
                            help=f"What should be used as column labels. Possible values: "
                                 f"{', '.join([v.value for v in AnnotationLabelsSource])}")

    workflow = parser.add_argument_group('workflow',
                                         'arguments for the workflow')
    workflow.add_argument('--threads', type=int, default=None)
    workflow.add_argument('--force', default=False, action='store_true')
    workflow.add_argument('--verbose', default=False, action='store_true')
    workflow.add_argument('--dryrun', default=False, action='store_true')
    workflow.add_argument('--exec-cmd', type=str, default=None)

    parser.set_defaults(func=init_build)


def init_build(args):
    was_successful = run_build_workflow(
        args.sample_list_path,
        args.output_dir,
        k=args.k,
        base_name=args.base_name,
        build_primary_graph=args.build_primary_graph,
        annotation_formats=[AnnotationFormats(af) for af in args.annotation_format],
        annotation_labels_source=AnnotationLabelsSource.SEQUENCE_HEADERS,
        exec_cmd=args.exec_cmd,
        threads=args.threads,
        force=args.force,
        verbose=args.verbose,
        dryrun=args.dryrun
    )

    if not was_successful:
        logging.error("Workflow did not run successfully")
        sys.exit(1)

