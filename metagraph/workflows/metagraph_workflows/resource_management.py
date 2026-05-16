import json
import math
import os
from pathlib import Path
from snakemake.common.tbdstring import TBDString


from metagraph_workflows import workflow_configs
from metagraph_workflows.workflow_configs import MEM_MB_KEY, DISK_MB_KEY, \
    MEM_BUFFER_MB_KEY, THREADS_KEY
from metagraph_workflows.utils import logger, get_rule_specific_config

BASE_MEM = 1 * 1024
FALLBACK_MAX_MEM = 4 * 1024
FALLBACK_MAX_DISK = 10 * 1024



def _get_max_memory(config):
    return config.get(workflow_configs.MAX_MEMORY_MB, FALLBACK_MAX_MEM)


def _get_max_disk(config):
    return config.get(workflow_configs.MAX_DISK_MB, FALLBACK_MAX_DISK)


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

    # If True, clip the auto-derived buffer to `max_buffer_size_mb`.
    # build/annotate need this (their --mem-cap-gb is a preallocated
    # buffer that competes with other in-memory state). Row-diff
    # transforms scale with annotation size and routinely need
    # >> max_buffer_size_mb, so those subclasses opt out.
    APPLY_MAX_BUFFER_CAP = True

    def get_mem_buffer_gib(self):
        """
        value for the `--mem-cap-gb` parameter (in GiB)
        """
        def _get_mem_buffer(wildcards, input, threads, resources):
            mem_cap_mb = get_rule_specific_config(self.rule_name,
                                                  MEM_BUFFER_MB_KEY, self.config)

            if not mem_cap_mb:
                mem_cap_mb = self._mem_buf_estimate(wildcards, resources, input, threads)
                if mem_cap_mb == TBDString():
                    return TBDString()
                if self.APPLY_MAX_BUFFER_CAP:
                    mem_cap_mb = min(mem_cap_mb, self.config[workflow_configs.MAX_BUFFER_SIZE_MB])

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

        if avail_mem_mb == TBDString():
            return TBDString()
        return max(int(self.CAP_MEM_FRACTION * avail_mem_mb), 1024) # TODO: parametrize constant?


class SupportsMemBufferSizeWithEstimation(SupportsMemBufferSize):
    """
    Base class for cases where we have a heuristic to estimate the required mem cap.
    """
    def _mem_buf_estimate(self, wildcards, resources, input, threads) -> int:
        raise NotImplementedError("Mixing in SupportsMemoryCapWithEstimation requires reimplementing mem_cap_estimate")

    def _get_mem_estimate(self, wildcards, input, threads):
        mem_cap = self.get_mem_buffer_gib()(wildcards, input, threads, None)

        if mem_cap == TBDString():
            return TBDString()

        mem_cap_mib = mem_cap*1024

        # adjusting memory s.t. mem_cap is CAP_MEM_FRACTION of the overall memory
        # (to be consistent with the default heuristic in SupportsMemoryCap
        return mem_cap_mib + max(int(mem_cap_mib/self.CAP_MEM_FRACTION * (1-self.CAP_MEM_FRACTION)), BASE_MEM)


class BuildGraphResources(SupportsMemBufferSize):
    pass


class BuildGraphResourcesWithKmerEstimates(SupportsMemBufferSizeWithEstimation):

    KMC_STATS_KEY = "Stats"
    KMC_UNIQUE_KMER_CNT = "#Unique_counted_k-mers"

    def _mem_buf_estimate(self, wildcards, resources, input, threads) -> int:
        kmc_json_path = Path(input['kmer'])

        if not kmc_json_path.exists():
            return TBDString()

        with open(kmc_json_path, 'r') as f:
            kmc_data = json.load(f)

        unique_kmers = kmc_data[self.KMC_STATS_KEY][self.KMC_UNIQUE_KMER_CNT]

        # Two independent factors:
        #   - expansion factor 2.6: 2x for canonical+reverse, ~30% extra
        #     for dummy k-mers (typically 10-30%); ~= 2 * 1.3.
        #   - 2.6 bytes per stored k-mer during succinct-graph
        #     construction (peak working-set per k-mer, empirical).
        # Total peak RSS ~= 6.76 bytes per unique k-mer reported by KMC.
        # If this underestimates on a particular dataset, override via
        # `rules.<rule>.mem_buffer_mb` in the workflow config.
        EXPANSION_FACTOR = 2.6
        BYTES_PER_KMER = 2.6
        required_ram_bytes = unique_kmers * EXPANSION_FACTOR * BYTES_PER_KMER
        required_ram_mb = max(int(math.ceil(required_ram_bytes / 1024**2)), 1024)
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

        return TBDString()


class TransformRdStage0Resources(SupportsMemBufferSizeWithEstimation):
    APPLY_MAX_BUFFER_CAP = False

    def __init__(self, config):
        super().__init__('transform_rd_stage0', config)

    def _mem_buf_estimate(self, wildcards, resources, input, threads):
        if Path(input.columns_file).exists():
            return int(columns_size_mb(input.columns_file) + BASE_MEM)
        return TBDString()


class TransformRdStage1Resources(SupportsMemBufferSize):
    APPLY_MAX_BUFFER_CAP = False

    def __init__(self, config):
        super().__init__('transform_rd_stage1', config)


class TransformRdStage2Resources(SupportsMemBufferSize):
    APPLY_MAX_BUFFER_CAP = False

    def __init__(self, config):
        super().__init__('transform_rd_stage2', config)


class AnnotateResources(SupportsMemBufferSize):
    """
    `metagraph annotate --mem-cap-gb` is a *per-column* buffer, applied
    independently to each column constructed in parallel under --separately.
    Total memory used is approximately `mem_cap_gb * parallel_cols`, where
    `parallel_cols = max_threads // threads_each`. Override the buffer
    heuristic accordingly.
    """

    def __init__(self, config):
        super().__init__('annotate', config)
        threads_each = config.get(workflow_configs.ANNOTATE_THREADS_EACH, 1) or 1
        self.threads_each = max(int(threads_each), 1)

    def get_parallel_cols(self, threads):
        """Number of columns built in parallel.

        Ceiling division: when `threads` is not a multiple of
        `threads_each` we'd rather fit one extra column than leave
        Snakemake-reserved cores idle.
        """
        return max(math.ceil(int(threads) / self.threads_each), 1)

    def get_effective_threads_each(self, threads):
        """Threads passed to `--threads-each`.

        Ceiling division so no reserved core sits idle when `threads` is
        not divisible by `parallel_cols`. May overcommit by up to
        `parallel_cols - 1` threads at boundary inputs (e.g. 13/8 -> 2*7
        = 14), which the OS scheduler absorbs.
        """
        return max(math.ceil(int(threads) / self.get_parallel_cols(threads)), 1)

    def get_mem_buffer_gib(self):
        """
        Returns the per-column buffer in GiB for `--mem-cap-gb`. If
        `mem_buffer_mb` is set in the rule config, treat it as the
        per-column value (matches the CLI semantic). Otherwise split the
        global mem budget across columns built in parallel.
        """
        def _get_mem_buffer(wildcards, input, threads, resources):
            mem_cap_mb = get_rule_specific_config(self.rule_name,
                                                  MEM_BUFFER_MB_KEY, self.config)

            if not mem_cap_mb:
                avail_mem_mb = get_rule_specific_config(self.rule_name, MEM_MB_KEY,
                                                        self.config)
                if not avail_mem_mb:
                    avail_mem_mb = resources.get('mem_mb', _get_max_memory(self.config))

                if avail_mem_mb == TBDString():
                    return TBDString()

                total_buf_mb = max(int(self.CAP_MEM_FRACTION * avail_mem_mb), 1024)
                total_buf_mb = min(total_buf_mb, self.config[workflow_configs.MAX_BUFFER_SIZE_MB])

                parallel_cols = self.get_parallel_cols(threads)
                mem_cap_mb = max(total_buf_mb // parallel_cols, 1024)

            return int(math.ceil(mem_cap_mb / 1024.0))

        return _get_mem_buffer
