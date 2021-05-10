import itertools
import re
from pathlib import Path
from typing import Union

from metagraph_workflows import constants


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
