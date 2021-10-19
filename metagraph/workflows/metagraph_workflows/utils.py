import itertools
import logging
import re
import subprocess
from pathlib import Path
from typing import Union

from metagraph_workflows import constants, utils
from metagraph_workflows.constants import GNU_TIME_CMD, TMP_DIR, \
    RULE_CONFIGS_KEY, SEQS_FILE_LIST_PATH, SEQS_DIR_PATH

logger = logging.getLogger("metagraph_workflow")


def get_seqs_file_list_path(wdir, config):
    if SEQS_FILE_LIST_PATH in config:
        return config[SEQS_FILE_LIST_PATH]

    seqs_file_list_path = wdir/'sequence_file_list_path.txt'
    seqs_dir_path = config.get(SEQS_DIR_PATH, None)

    if not seqs_dir_path:
        raise ValueError(f"Neither {SEQS_FILE_LIST_PATH} nor {SEQS_DIR_PATH} parameter are set. Need either to proceed")

    utils.create_transcript_path_list(seqs_dir_path, seqs_file_list_path)
    return seqs_file_list_path


def take_value_or_default(key, default, config):
    return config[key] if (key in config.keys() and config[key]) else default


def create_transcript_path_list(path: Union[Path, str], transcript_path: Union[Path, str], suffix=''):
    paths = [str(p.absolute()) for p in Path(path).glob(f'*{suffix}')]

    with open(transcript_path, 'w') as f:
        f.write('\n'.join(paths))


def get_sample_name(l):
    file_name = Path(l.strip()).name

    m = re.compile(r'^([^.]*)\.(fasta|[a-zA-Z]{2,4})(\.gz)?$').match(file_name)
    if m:
        return m.groups()[0]

    return file_name


def derive_sample_dictionary(transcript_path_list_path: Union[Path, str]):
    with open(transcript_path_list_path) as f:
        ret = {get_sample_name(l): l.strip() for l in f}
    return ret


def get_build_single_sample_input(config, orig_samples_path, seq_ids_dict):
    def _sample_input(wildcards):
        sample_id = wildcards[0] # TODO:

        if config[constants.SAMPLE_IDS_PATH]:
            return orig_samples_path / f"{{sample_id}}{config[constants.SAMPLE_STAGING_FILE_ENDING]}"
        else:
            return seq_ids_dict[sample_id]

    return _sample_input


def get_build_joint_input(config, contigs_dir, seq_ids_dict, seqs_file_list_path):
    sample_ids = set()
    if constants.SAMPLE_IDS_PATH in config and config[constants.SAMPLE_IDS_PATH]:
        with open(config[constants.SAMPLE_IDS_PATH]) as f:
            sample_ids = {f"{l.strip()}" for l in f}

    def _get_build_graph_input(wildcards):
        if config[constants.PRIMARIZE_SAMPLES_SEPARATELY]:
            all_samples = sample_ids if sample_ids else seq_ids_dict.keys()
            return [contigs_dir/f"{sample_id}_primary.fasta.gz" for sample_id in all_samples]
        else:
            return seqs_file_list_path

    return _get_build_graph_input


def generate_col_paths(annotation_cols_path, seqs_file_list_path, config):
    sample_names = set()

    if constants.SAMPLE_IDS_PATH in config and config[constants.SAMPLE_IDS_PATH]:
        with open(config[constants.SAMPLE_IDS_PATH]) as f:
            sample_names = { f"{l.strip()}_primary.fasta.gz" for l in f}

    else:
        with open(seqs_file_list_path) as f:
            column_names = [f"{f.strip().rstrip('/').split('/')[-1]}" for f in
                            f.readlines()]

            duplicate_col_names = [grp_key for (grp_key, names_lst) in
                                   itertools.groupby(sorted(column_names)) if
                                   len(list(names_lst)) > 1]

            assert not duplicate_col_names, f"Found duplicate filenames: {', '.join(duplicate_col_names)}"

            if config[constants.PRIMARIZE_SAMPLES_SEPARATELY]:
                sample_names = {f"{get_sample_name(c)}_primary.fasta.gz" for c in column_names}
            else:
                sample_names = set(column_names)

    return [annotation_cols_path / f"{c}.column.annodbg" for c in
            sample_names]


def get_wdir(config):
    return Path(config['output_directory'])


def get_gnu_time_command(config):
    EMTPY_CMD = ''
    cmd = config.get(GNU_TIME_CMD, EMTPY_CMD)

    if cmd:
        test_cmd=[cmd, '--version']
        proc = subprocess.run(test_cmd, capture_output=True)
        if proc.returncode == 0:
            return f"{cmd} --verbose"
        else:
            logger.error(f"Command {' '.join(test_cmd)} for GNU time could not be executed successfully: {proc.stderr}")
    else:
        logger.warning("No GNU Time command provided.")

    return EMTPY_CMD


def get_log_path(rule_name, config, wildcards=None):
    log_dir = get_wdir(config)/'logs'

    if wildcards:
        wildcard_str = '_'.join([f"{{{w}}}" for w in wildcards])
        return f"{log_dir}/{rule_name}/{rule_name}_{wildcard_str}.log"
    else:
        return f"{log_dir}/{rule_name}.log"


def temp_dir_config(config):
    return f"--disk-swap {config[TMP_DIR]}" if TMP_DIR in config else '',


def get_rule_specific_config(rule, key, config):
    if RULE_CONFIGS_KEY in config and rule in config[
        RULE_CONFIGS_KEY] and key in config[RULE_CONFIGS_KEY][rule]:
        return config[RULE_CONFIGS_KEY][rule][key]
    return None