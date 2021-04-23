import pytest

from metagraph_workflows import resource_management as rm


@pytest.fixture()
def config():
    return {
        'max_memory_mb': 4100
    }


def test_TransformRdStage1Resources(config):
    rule_name = 'transform_rd_stage1'
    inst = rm.TransformRdStage1Resources()

    assert inst.get_mem(config)(None, None, None) == 4100

    base_mem = 1024
    mem = 2345
    config['rules'] = {rule_name: {'mem_mb': mem}}
    assert inst.get_mem(config)(None, None, None) == mem
    res = {'mem_mb': mem}
    assert inst.get_mem_cap(config)(None, None, None, res) == pytest.approx((mem-base_mem)/1024)

    config['rules'][rule_name]['mem_cap_mb'] = 1000
    assert inst.get_mem(config)(None, None, None) == 1000 + base_mem
    assert inst.get_mem_cap(config)(None, None, None, res) == pytest.approx(1000/1024)
