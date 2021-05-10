import json
import math
import os
from pathlib import Path



from metagraph_workflows import constants
from metagraph_workflows.cfg_utils import get_rule_specific_config
from metagraph_workflows.constants import MEM_MB_KEY, DISK_MB_KEY, \
    MEM_BUFFER_MB_KEY, THREADS_KEY
from metagraph_workflows.utils import logger

BASE_MEM = 1 * 1024
FALLBACK_MAX_MEM = 4 * 1024
FALLBACK_MAX_DISK = 10 * 1024


# obviously wrong value to work around --dryrun issues when the resources requirements
# and other values depend on the input of rule (can be removed once is fixed
# https://github.com/snakemake/snakemake/issues/927)
TBD_VALUE = -1

def _get_max_memory(config):
    return config.get(constants.MAX_MEMORY_MB, FALLBACK_MAX_MEM)


def _get_max_disk(config):
    return config.get(constants.MAX_DISK_MB, FALLBACK_MAX_DISK)


def columns_size_mb(columns_file):
    with open(columns_file) as f:
        col_file_size_bytes = sum(os.stat(l.strip()).st_size for l in f)
        return col_file_size_bytes / 1024 ** 2


class ResourceConfig:
    def __init__(self, rule_name, config):
        self.rule_name = rule_name
        self.config = config

    def get_threads(self, max_threads) -> int:
        threads = get_rule_specific_config(self.rule_name, THREADS_KEY, self.config)

        if not threads:
            threads = max_threads
        return threads

    def get_mem(self):
        def _get_mem(wildcards, input, threads) -> int:
            mem_mb = get_rule_specific_config(self.rule_name, MEM_MB_KEY,
                                              self.config)
            if not mem_mb:
                mem_mb = self._get_mem_estimate(wildcards, input, threads)

                max_mem = _get_max_memory(self.config)
                if mem_mb > max_mem:
                    logger.warning(
                        f"The estimated memory of {mem_mb} MB "
                        f"is larger than the max memory {max_mem}.")

            return mem_mb

        return _get_mem

    def _get_mem_estimate(self, wildcards, input, threads):
        return _get_max_memory(self.config)

    def get_disk(self):
        def _get_disk(wildcards):
            disk_mb = get_rule_specific_config(self.rule_name, DISK_MB_KEY,
                                               self.config)
            if not disk_mb:
                disk_mb = _get_max_disk(self.config)
            return disk_mb

        return _get_disk


class SupportsMemBufferSize(ResourceConfig):
    MEM_OVERHEAD = BASE_MEM

    CAP_MEM_FRACTION = 0.85

    def get_mem_buffer_gib(self):
        """
        value for the `--mem-cap-gb` parameter (in GiB)
        """
        def _get_mem_buffer(wildcards, input, threads, resources):
            mem_cap_mb = get_rule_specific_config(self.rule_name,
                                                  MEM_BUFFER_MB_KEY, self.config)

            if not mem_cap_mb:
                mem_cap_mb = min(self._mem_buf_estimate(wildcards, resources, input, threads), self.config[constants.MAX_BUFFER_SIZE_MB])

                if mem_cap_mb == TBD_VALUE:
                    return TBD_VALUE

            return int(math.ceil(mem_cap_mb / 1024.0))

        return _get_mem_buffer

    def _mem_buf_estimate(self, wildcards, resources, input, threads):
        """
        Default estimation for mem cap: get a percentage of the available memory
        """
        avail_mem_mb = get_rule_specific_config(self.rule_name, MEM_MB_KEY,
                                          self.config)

        if not avail_mem_mb:
            avail_mem_mb = resources.get('mem_mb', _get_max_memory(self.config))

        return max(int(self.CAP_MEM_FRACTION * avail_mem_mb), 1024) # TODO: parametrize constant?


class SupportsMemBufferSizeWithEstimation(SupportsMemBufferSize):
    """
    Base class for cases where we have a heuristic to estimate the required mem cap.
    """
    def _mem_buf_estimate(self, wildcards, resources, input, threads) -> int:
        raise NotImplementedError("Mixing in SupportsMemoryCapWithEstimation requires reimplementing mem_cap_estimate")

    def _get_mem_estimate(self, wildcards, input, threads):
        mem_cap = self.get_mem_buffer_gib()(wildcards, input, threads, None)

        if mem_cap == TBD_VALUE:
            return TBD_VALUE

        mem_cap_mib = mem_cap*1024

        # adjusting memory s.t. mem_cap is CAP_MEM_FRACTION of the overall memory
        # (to be consistent with the default heuristic in SupportsMemoryCap
        return mem_cap_mib + max(int(mem_cap_mib/self.CAP_MEM_FRACTION * (1-self.CAP_MEM_FRACTION)), BASE_MEM)


class SupportsDiskCap(ResourceConfig):
    def get_disk_cap(self):
        return self.get_disk() # TODO: come up with a heuristic


class BuildGraphResources(SupportsMemBufferSize, SupportsDiskCap):
    pass


class BuildGraphResourcesWithKmerEstimates(SupportsMemBufferSizeWithEstimation, SupportsDiskCap):

    KMC_STATS_KEY = "Stats"
    KMC_UNIQUE_KMER_CNT = "#Unique_counted_k-mers"

    def _mem_buf_estimate(self, wildcards, resources, input, threads) -> int:
        kmc_json_path = Path(input['kmer'])

        if not kmc_json_path.exists():
            return TBD_VALUE

        with open(kmc_json_path, 'r') as f:
            kmc_data = json.load(f)

        unique_kmers = kmc_data[self.KMC_STATS_KEY][self.KMC_UNIQUE_KMER_CNT]

        bytes_per_kmer = 2.6
        kmer_count = 2.6 * unique_kmers  # 2x canonical+non-canonical +  ~30% for dummy kmers (typically it's 10%)
        required_ram = int(math.ceil(kmer_count * bytes_per_kmer / 1024**2))
        required_ram_mb = max(required_ram, 1024)

        return required_ram_mb


class PrimarizeCanonicalGraphSingleSampleResources(ResourceConfig):
    def __init__(self, config):
        super().__init__('primarize_canonical_graph_single_sample', config)

    def _get_mem_estimate(self, wildcards, input, threads):
        input_path = Path(input[0])

        if input_path.exists():
            file_size_mib = max(int(math.ceil(input_path.stat().st_size / 1024.0**2)), 1)
            logger.debug(f"File size of {input_path.name} is {file_size_mib}")

            # factor 2 is based on experiments on the mouse data set.
            # In most cases factor 1.3 to 1.5 would be enough, however, there are outliers
            return 2*file_size_mib

        return TBD_VALUE


class TransformRdStage0Resources(SupportsMemBufferSizeWithEstimation):
    def __init__(self, config):
        super().__init__('transform_rd_stage0', config)

    def _mem_buf_estimate(self, wildcards, resources, input, threads):
        if Path(input.columns_file).exists():
            return int(columns_size_mb(input.columns_file) + BASE_MEM)
        return TBD_VALUE


class TransformRdStage1Resources(SupportsMemBufferSize):
    def __init__(self, config):
        super().__init__('transform_rd_stage1', config)


class TransformRdStage2Resources(SupportsMemBufferSize):
    def __init__(self, config):
        super().__init__('transform_rd_stage2', config)
