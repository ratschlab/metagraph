from metagraph_workflows.constants import RULE_CONFIGS_KEY, TMP_DIR
from pathlib import Path


def get_wdir(config):
    return Path(config['output_directory'])


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