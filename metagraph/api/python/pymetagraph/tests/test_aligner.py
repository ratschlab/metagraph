import os

import pytest
from readfish.plugins.utils import Result

from pymetagraph import Aligner

FIXTURE_DIR = os.environ.get("PYMETAGRAPH_TEST_FIXTURE_DIR", "/tmp/pymetagraph_fixture")
GRAPH = os.path.join(FIXTURE_DIR, "idx.dbg")
ANNOTATOR = os.path.join(FIXTURE_DIR, "idx.column.annodbg")

MATCHING_SEQ = (
    "GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTC"
    "TCTTAGCCCAGACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGA"
)
EXPECTED_LABEL = (
    "ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|"
    "OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|"
)


@pytest.fixture
def aligner():
    a = Aligner(debug_log=None, input=GRAPH, annotator=ANNOTATOR)
    yield a
    a.disconnect()


@pytest.fixture
def align_aligner():
    a = Aligner(debug_log=None, input=GRAPH, annotator=ANNOTATOR, method="align")
    yield a
    a.disconnect()


def test_validate_missing_input():
    with pytest.raises(AttributeError):
        Aligner(debug_log=None)


def test_initialised(aligner):
    assert aligner.initialised is True


def test_describe(aligner):
    description = aligner.describe(regions=[], barcodes={})
    assert GRAPH in description


def test_map_reads(aligner):
    results = [
        Result(channel=1, read_id="matching", seq=MATCHING_SEQ),
        Result(channel=2, read_id="no_seq", seq=""),
        Result(channel=3, read_id="garbage", seq="N" * 100),
    ]

    by_id = {r.read_id: r for r in aligner.map_reads(iter(results))}

    assert len(by_id) == 3
    assert by_id["matching"].alignment_data
    assert by_id["matching"].alignment_data[0].ctg == EXPECTED_LABEL
    assert by_id["no_seq"].alignment_data == []
    assert by_id["garbage"].alignment_data == []


def test_map_reads_align_method(align_aligner):
    results = [
        Result(channel=1, read_id="matching", seq=MATCHING_SEQ),
        Result(channel=2, read_id="no_seq", seq=""),
        Result(channel=3, read_id="garbage", seq="N" * 100),
    ]

    by_id = {r.read_id: r for r in align_aligner.map_reads(iter(results))}

    assert len(by_id) == 3
    matching = by_id["matching"].alignment_data
    assert matching
    assert matching[0].ctg == EXPECTED_LABEL
    # Real coordinates, unlike method="query"'s dummy 0/len(seq)/+1.
    assert 0 <= matching[0].r_st < matching[0].r_en <= len(MATCHING_SEQ)
    assert matching[0].strand in (1, -1)
    assert by_id["no_seq"].alignment_data == []
    assert by_id["garbage"].alignment_data == []
