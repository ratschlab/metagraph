from metagraph_workflows.constants import RULE_CONFIGS_KEY, TMP_DIR


def temp_dir_config(config):
    return f"--disk-swap {config[TMP_DIR]}" if TMP_DIR in config else '',


def get_rule_specific_config(rule, key, config):
    if RULE_CONFIGS_KEY in config and rule in config[
        RULE_CONFIGS_KEY] and key in config[RULE_CONFIGS_KEY][rule]:
        return config[RULE_CONFIGS_KEY][rule][key]
    return None