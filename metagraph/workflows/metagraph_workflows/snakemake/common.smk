import itertools

from metagraph_workflows import utils
from metagraph_workflows.constants import SEQS_FILE_LIST_PATH, SEQS_DIR_PATH


def take_value_or_default(key, default):
    return config[key] if (key in config.keys() and config[key]) else default

def get_log_opt(rule):
    return f" | tee {log_dir}/{graph}_{rule}.log" if write_logs else ''

def get_seqs_file_list_path():
    if SEQS_FILE_LIST_PATH in config:
        return config[SEQS_FILE_LIST_PATH]

    seqs_file_list_path = wdir/'sequence_file_list_path.txt'
    seqs_dir_path = config.get(SEQS_DIR_PATH, None)

    if not seqs_dir_path:
        raise ValueError(f"Neither {SEQS_FILE_LIST_PATH} nor {SEQS_DIR_PATH} parameter are set. Need either to proceed")

    utils.create_transcript_path_list(seqs_dir_path, seqs_file_list_path)
    return seqs_file_list_path


disk_cap_factor=0.9

def get_disk(key):
    return config[key] if key in config else config['default_disk_mb']

def get_disk_cap(wilcards, resources):
    return max(int(resources.disk_mb * disk_cap_factor/1024), 1)


def generate_col_paths(seq_file):
    with open(seq_file) as f:
        column_names = [f"{f.strip().split('/')[-1]}" for f in f.readlines()]

        duplicate_col_names = [grp_key for (grp_key, names_lst) in itertools.groupby(sorted(column_names)) if len(list(names_lst)) > 1]
        assert not duplicate_col_names, f"Found duplicate filenames: { ', '.join(duplicate_col_names)}"

        return [ annotation_cols_path/f"{c}.column.annodbg" for c in column_names]
