"""Tests for pymetagraph.AnnotatedGraph.

The fixture is built from metagraph's own test data with the `metagraph`
binary (taken from $METAGRAPH, else from PATH); the tests are skipped if
neither is available. Query results are compared to `metagraph query`.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import pymetagraph

DATA = Path(__file__).resolve().parents[4] / "tests" / "data" / "metasub_fake_data.fa"
K = 21
METAGRAPH = shutil.which(os.environ.get("METAGRAPH") or "metagraph")

pytestmark = pytest.mark.skipif(METAGRAPH is None, reason="metagraph binary not found")


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        else:
            records[name] += line
    return records


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    out = tmp_path_factory.mktemp("index") / "graph"
    subprocess.run([METAGRAPH, "build", "-k", str(K), "-o", str(out), str(DATA)],
                   check=True, capture_output=True)
    subprocess.run([METAGRAPH, "annotate", "-i", f"{out}.dbg", "--anno-header",
                    "-o", str(out), str(DATA)], check=True, capture_output=True)
    return Path(f"{out}.dbg"), Path(f"{out}.column.annodbg")


@pytest.fixture(scope="module")
def graph(index):
    return pymetagraph.AnnotatedGraph(*index, threads=2)


def cli_query(index, fasta, discovery_fraction):
    """{sequence index: [(label, kmer_count)]} as reported by `metagraph query`."""
    result = subprocess.run(
        [METAGRAPH, "query", "-i", str(index[0]), "-a", str(index[1]),
         "--min-kmers-fraction-label", str(discovery_fraction), str(fasta)],
        check=True, capture_output=True, text=True)
    hits = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        # each hit is printed as <label>:<kmer_count>
        hits[int(fields[0])] = [(m.group(1), int(m.group(2)))
                                for m in (re.fullmatch(r"<(.*)>:(\d+)", f) for f in fields[2:] if f)]
    return hits


def test_properties(graph):
    assert graph.k == K
    assert sorted(graph.labels) == sorted(read_fasta(DATA))


def test_matches_cli(graph, index, tmp_path):
    seqs = list(read_fasta(DATA).values())
    queries = [seqs[0], seqs[1][:100], seqs[2][50:150], "ACGT" * 30]
    fasta = tmp_path / "queries.fa"
    fasta.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(queries)))

    for discovery_fraction in (0.0, 0.5, 1.0):
        expected = cli_query(index, fasta, discovery_fraction)
        result = graph.query(queries, num_top_labels=len(graph.labels),
                             discovery_fraction=discovery_fraction)
        for i, hits in enumerate(result):
            assert sorted(hits) == sorted(expected.get(i, []))


def test_counts_and_order(graph):
    seq = next(iter(read_fasta(DATA).values()))
    hits = graph.query(seq, num_top_labels=len(graph.labels))
    counts = [count for _, count in hits]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] == len(seq) - K + 1


def test_single_sequence_and_edge_cases(graph):
    seq = next(iter(read_fasta(DATA).values()))
    assert graph.query(seq) == graph.query([seq])[0]
    assert graph.query(["", "A" * (K - 1)]) == [[], []]
    assert graph.query([]) == []


def test_invalid_arguments(graph, index):
    with pytest.raises(FileNotFoundError):
        pymetagraph.AnnotatedGraph("missing.dbg", index[1])
    with pytest.raises(ValueError):
        pymetagraph.AnnotatedGraph(*index, threads=0)
    with pytest.raises(ValueError):
        graph.query(["ACGT"], discovery_fraction=1.5)


def test_concurrent_graph_construction(index, graph):
    from concurrent.futures import ThreadPoolExecutor
    query = next(iter(read_fasta(DATA).values()))[:100]
    expected = graph.query(query)

    def load_and_query(threads):
        loaded = pymetagraph.AnnotatedGraph(*index, threads=threads)
        return loaded.query(query)

    with ThreadPoolExecutor(max_workers=3) as executor:
        assert list(executor.map(load_and_query, (1, 2, 3, 3, 2, 1))) == [expected] * 6
