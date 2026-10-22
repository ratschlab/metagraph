"""Exercise the native alignment API against the CLI and reference sequences."""
import json
import os
import shutil
import subprocess
import struct
import sys

import pytest
import pymetagraph

METAGRAPH = shutil.which(os.environ.get("METAGRAPH") or "metagraph")
pytestmark = pytest.mark.skipif(METAGRAPH is None, reason="metagraph binary not found")
REFERENCE = "AACGTTGCAATGCCGATTCAGGCTAAC"
OTHER = "TTGATAGCCTTGACCGATAGCTTCCA"


def run(*args):
    return subprocess.run([METAGRAPH, *map(str, args)], check=True,
                          capture_output=True, text=True).stdout


def test_sshash_version_error_keeps_python_alive(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(f">reference\n{REFERENCE}\n")
    run("build", "--graph", "sshash", "--disk-swap", tmp_path, "-k", 5,
        "-o", tmp_path / "graph", fasta)
    graph = tmp_path / "graph.sshashdbg"
    data = bytearray(graph.read_bytes())
    # The native-format wrapper stores two size_t values and the mode enum
    # before the three SSHash version bytes.
    data[2 * struct.calcsize('N') + struct.calcsize('i')] = 5
    old_graph = tmp_path / "v5.sshashdbg"
    old_graph.write_bytes(data)
    code = """
import pymetagraph, sys
try:
    pymetagraph.AnnotatedGraph(sys.argv[1])
except RuntimeError as error:
    assert 'rebuild' in str(error)
else:
    raise AssertionError('An incompatible SSHash graph loaded')
assert pymetagraph.AnnotatedGraph(sys.argv[2]).k == 5
print('interpreter alive')
"""
    result = subprocess.run([sys.executable, "-c", code, str(old_graph), str(graph)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'interpreter alive' in result.stdout


@pytest.fixture(scope="module", params=["succinct", "sshash"])
def index(request, tmp_path_factory):
    folder = tmp_path_factory.mktemp("alignment_" + request.param)
    fasta = folder / "reference.fa"
    fasta.write_text(f">alpha\n{REFERENCE}\n>beta\n{OTHER}\n")
    run("build", "--graph", "succinct", "-k", 5, "-o", folder / "graph", fasta)
    if request.param == "sshash":
        # SSHash consumes a spectrum-preserving string set, not raw reads.
        run("transform", "--to-fasta", "-o", folder / "unitigs", folder / "graph.dbg")
        run("build", "--graph", "sshash", "--disk-swap", folder, "-k", 5,
            "-o", folder / "graph", folder / "unitigs.fasta.gz")
    extension = ".dbg" if request.param == "succinct" else ".sshashdbg"
    graph = folder / ("graph" + extension)
    run("annotate", "--anno-header", "--coordinates", "-i", graph,
        "-o", folder / "annotation", fasta)
    run("transform_anno", "--anno-type", "column_coord", "--coordinates",
        "-o", folder / "annotation", folder / "annotation.column.annodbg")
    return graph, folder / "annotation.column_coord.annodbg"


@pytest.fixture(scope="module")
def graph(index):
    return pymetagraph.AnnotatedGraph(*index, threads=2)


def test_alignment_matches_cli(graph, index, tmp_path):
    query = REFERENCE[4:22]
    fasta = tmp_path / "query.fa"
    fasta.write_text(f">query\n{query}\n")
    expected = [json.loads(line) for line in run("align", "--json", "--align-only-forwards",
                              "--align-min-seed-length", 5, "--align-min-exact-match", 0,
                              "-i", index[0], "-a", index[1], fasta).splitlines()]
    actual = graph.align(query, seed_length=5, min_exact_match=0, only_forwards=True)
    assert len(actual) == len(expected)
    assert [(a["score"], a["annotation"]) for a in actual] == [
        (e["score"], e["annotation"]) for e in expected]
    assert actual[0]["query_start"] == 0
    assert actual[0]["query_end"] == len(query)
    assert actual[0]["annotation"]["labels"][0]["nt_coords"] == "5-22"


def test_batches_and_empty_inputs(graph):
    queries = [REFERENCE, OTHER, "", "NNNNNNNN"]
    result = graph.align(queries, seed_length=5, only_forwards=True)
    assert result == [graph.align(q, seed_length=5, only_forwards=True) for q in queries]
    assert result[2:] == [[], []]
    assert graph.align([]) == []
    assert all(len(hits) <= 1 for hits in result)


def test_unannotated_reverse_alignment(index):
    graph = pymetagraph.AnnotatedGraph(index[0], threads=2)
    query = REFERENCE.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    result = graph.align(query, seed_length=5)
    assert len(result) == 1
    assert result[0]["strand"] == "-"
    assert result[0]["sequence"] == query
    assert result[0]["query_start"] == 0
    assert result[0]["query_end"] == len(query)
    assert result[0]["annotation"]["ref_sequence"] == REFERENCE
    assert graph.labels == []
    with pytest.raises(ValueError, match="annotation"):
        graph.query(query)


def test_chaining_options(graph):
    for connect in (False, True):
        for extend in (False, True):
            hits = graph.align(REFERENCE, seed_length=5, only_forwards=True,
                               connect_anchors=connect, extend_chains=extend)
            assert hits
            assert hits[0]["query_end"] > hits[0]["query_start"]


def test_search_limits(graph):
    query = REFERENCE[:8] + "C" + REFERENCE[9:]
    assert graph.align(REFERENCE, seed_length=5, max_num_seeds_per_locus=0) == []
    for options in ({"max_num_nodes_per_seq_char": 0}, {"max_ram_per_alignment": 0}):
        hits = graph.align(query, seed_length=5, min_exact_match=0,
                           only_forwards=True, **options)
        assert hits
        assert hits[0]["query_end"] - hits[0]["query_start"] < len(query)
    assert graph.align(query, seed_length=5, min_exact_match=0, only_forwards=True,
                       max_num_nodes_per_seq_char=1000, max_ram_per_alignment=1)[0]["annotation"]["cigar"] == "8=1X18="


def test_inexact_chain_errors_leave_workers_usable(graph):
    query = REFERENCE[:8] + "C" + REFERENCE[9:]
    # A raw chain has an unresolved graph gap and cannot use the alignment
    # JSON schema. The failed batch must drain and raise, not abort Python.
    with pytest.raises(ValueError, match="connect_anchors=True"):
        graph.align([query, REFERENCE], seed_length=5, min_exact_match=0,
                    only_forwards=True, connect_anchors=False, extend_chains=False)
    for extend in (False, True):
        hits = graph.align(query, seed_length=5, min_exact_match=0,
                           only_forwards=True, extend_chains=extend)
        assert hits
        assert hits[0]["annotation"]["ref_sequence"] == REFERENCE
        assert hits[0]["annotation"]["cigar"] == "8=1X18="
    assert graph.align(REFERENCE, seed_length=5, only_forwards=True)


@pytest.mark.parametrize("options", [
    {"xdrop": -1}, {"rel_score_cutoff": 1.1}, {"min_exact_match": float("nan")},
    {"max_alternative_alignments": 0}, {"seed_length": -1},
    {"rel_score_cutoff": float("inf")},
    {"max_num_nodes_per_seq_char": -1}, {"max_ram_per_alignment": float("nan")},
])
def test_invalid_alignment_options(graph, options):
    with pytest.raises((ValueError, TypeError)):
        graph.align(REFERENCE, **options)


def test_non_ascii_batch_leaves_workers_usable(graph):
    with pytest.raises(ValueError, match="ASCII"):
        graph.align([REFERENCE, "ACGéT", OTHER], seed_length=5)
    assert graph.align(REFERENCE, seed_length=5)
