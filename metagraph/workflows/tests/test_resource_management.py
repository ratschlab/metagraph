import pytest
import math

from metagraph_workflows import resource_management as rm
from metagraph_workflows import workflow_configs

@pytest.fixture()
def config():
    return {
        workflow_configs.MAX_MEMORY_MB: 16000,
        workflow_configs.MAX_BUFFER_SIZE_MB: 50000
    }


def test_TransformRdStage1Resources(config):
    rule_name = 'transform_rd_stage1'
    inst = rm.TransformRdStage1Resources(config)

    # by default get max available memory
    assert inst.get_mem()(None, None, None) == 16000

    base_mem = 1024

    # now explicitly setting available memory for the rule
    mem = 8000
    config['rules'] = {rule_name: {'mem_mb': mem}}
    assert inst.get_mem()(None, None, None) == mem

    resources = {'mem_mb': mem}
    expected_mb = int(rm.SupportsMemBufferSize.CAP_MEM_FRACTION * mem)
    assert inst.get_mem_buffer_gib()(None, None, None, resources) == int(math.ceil(expected_mb / 1024.0))

    # now additionally setting mem cap explicitly
    mem_buffer = 2048
    config['rules'][rule_name]['mem_buffer_mb'] = mem_buffer
    assert inst.get_mem()(None, None, None) == mem
    assert inst.get_mem_buffer_gib()(None, None, None, resources) == int(math.ceil(mem_buffer / 1024))


def test_AnnotateResources_per_column_buffer(config):
    # annotate --mem-cap-gb is per-column, so the global budget must be
    # divided across columns built in parallel (= threads // threads_each).
    mem = 64000  # MB available to the rule
    config[workflow_configs.ANNOTATE_THREADS_EACH] = 8
    inst = rm.AnnotateResources(config)

    threads = 64
    parallel_cols = threads // 8  # 8 parallel columns

    assert inst.get_parallel_cols(threads) == parallel_cols

    resources = {'mem_mb': mem}
    expected_total = int(rm.SupportsMemBufferSize.CAP_MEM_FRACTION * mem)
    expected_total = min(expected_total, config[workflow_configs.MAX_BUFFER_SIZE_MB])
    expected_per_col_mb = max(expected_total // parallel_cols, 1024)
    expected_gib = int(math.ceil(expected_per_col_mb / 1024.0))
    assert inst.get_mem_buffer_gib()(None, None, threads, resources) == expected_gib

    # threads_each not set in config -> falls back to 1 (workflow default
    # comes from default.yml, see ANNOTATE_THREADS_EACH there).
    config2 = {
        workflow_configs.MAX_MEMORY_MB: 16000,
        workflow_configs.MAX_BUFFER_SIZE_MB: 50000,
    }
    inst2 = rm.AnnotateResources(config2)
    assert inst2.threads_each == 1
    assert inst2.get_parallel_cols(8) == 8

    # explicit mem_buffer_mb is treated as the per-column value
    config3 = {
        workflow_configs.MAX_MEMORY_MB: 16000,
        workflow_configs.MAX_BUFFER_SIZE_MB: 50000,
        workflow_configs.ANNOTATE_THREADS_EACH: 4,
        'rules': {'annotate': {'mem_buffer_mb': 3000}},
    }
    inst3 = rm.AnnotateResources(config3)
    assert inst3.get_mem_buffer_gib()(None, None, 16, {'mem_mb': 16000}) == int(math.ceil(3000 / 1024))
