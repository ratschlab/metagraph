import pytest
import math

from metagraph_workflows import resource_management as rm
from metagraph_workflows import constants

@pytest.fixture()
def config():
    return {
        constants.MAX_MEMORY_MB: 16000
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
    assert inst.get_mem_buffer_gib()(None, None, None, resources) == int(math.ceil(0.8 * mem / 1024))

    # now additionally setting mem cap explicitly
    mem_buffer = 2048
    config['rules'][rule_name]['mem_buffer_mb'] = mem_buffer
    assert inst.get_mem()(None, None, None) == mem
    assert inst.get_mem_buffer_gib()(None, None, None, resources) == int(math.ceil(mem_buffer / 1024))
