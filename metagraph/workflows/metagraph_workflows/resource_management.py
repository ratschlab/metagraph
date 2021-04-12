import logging
import os

from metagraph_workflows import constants

RULE_CONFIGS_KEY = 'rules'

MEM_MB_KEY = 'mem_mb'
DISK_MB_KEY = 'disk_mb'

MEM_CAP_MB_KEY = 'mem_cap_mb'

BASE_MEM = 2 * 1024
FALLBACK_MAX_MEM = 4 * 1024


def _get_rule_specific_config(rule, key, config):
    if RULE_CONFIGS_KEY in config and rule in config[
        RULE_CONFIGS_KEY] and key in config[RULE_CONFIGS_KEY][rule]:
        return config[RULE_CONFIGS_KEY][rule][key]
    return None


def _get_max_memory(config):
    return config.get(constants.MAX_MEMORY_MB, FALLBACK_MAX_MEM)


def columns_size_mb(columns_file):
    with open(columns_file) as f:
        col_file_size_bytes = sum(os.stat(l.strip()).st_size for l in f)
        return col_file_size_bytes / 1024 ** 2


class ResourceConfig():
    def __init__(self, rule_name):
        self.rule_name = rule_name

    def get_mem(self, config):
        def _get_mem():
            mem_mb = _get_rule_specific_config(self.rule_name, MEM_MB_KEY,
                                               config)
            if not mem_mb:
                mem_mb = _get_max_memory(config)
            return mem_mb

        return _get_mem

    def get_base_mem_estimate(self, wildcards, input, threads):
        return BASE_MEM

    def get_disk(self, config):
        def _get_disk():
            disk_mb = _get_rule_specific_config(self.rule_name, DISK_MB_KEY,
                                               config)
            if not disk_mb:
                disk_mb = _get_max_memory(config)
            return disk_mb

        return _get_disk

    def get_base_disk_estimate(self, wildcards, input, threads):
        return lkj


class SupportsMemoryCap(ResourceConfig):
    def get_mem(self, config):
        def _get_mem(wildcards, input, threads):
            estimated_base_memory = self.get_base_mem_estimate(wildcards, input,
                                                               threads)

            mem_mb = _get_rule_specific_config(self.rule_name, MEM_MB_KEY,
                                               config)

            if not mem_mb:
                mem_mb = _get_max_memory(config)

            mem_cap_mb = _get_rule_specific_config(self.rule_name,
                                                   MEM_CAP_MB_KEY, config)

            if mem_cap_mb:
                estimated_total_mem = mem_cap_mb + estimated_base_memory

                if estimated_total_mem > mem_mb:
                    logging.warning(
                        f"The estimated memory of {estimated_total_mem} MB "
                        f"is larger than the max memory {mem_mb}. Mem-cap: {mem_cap_mb}, "
                        f"estimated base memory {estimated_base_memory}")

                return estimated_total_mem

            return mem_mb

        return _get_mem

    def get_mem_cap(self, config):
        def _get_mem_cap(wildcards, input, threads, resources):
            mem_cap_mb = _get_rule_specific_config(self.rule_name,
                                                   MEM_CAP_MB_KEY, config)

            if not mem_cap_mb:
                avail_mem_mb = resources.get('mem_mb', _get_max_memory(config))

                mem_cap_mb = max((avail_mem_mb - self.get_base_mem_estimate(
                    wildcards, input, threads)), 1024)

            return float(mem_cap_mb / 1024)

        return _get_mem_cap


class BuildResources(SupportsMemoryCap):
    def __init__(self):
        self.rule_name = 'build'


class TransformRdStage0Resources(SupportsMemoryCap):
    def __init__(self):
        self.rule_name = 'transform_rd_stage0'

    def get_base_mem_estimate(self, wildcards, input, threads):
        return columns_size_mb(input.columns_file) + BASE_MEM


class TransformRdStage1Resources(SupportsMemoryCap):
    def __init__(self):
        self.rule_name = 'transform_rd_stage1'


class TransformRdStage2Resources(SupportsMemoryCap):
    def __init__(self):
        self.rule_name = 'transform_rd_stage2'
