"""Tests for pymetagraph.AnnotatedGraph.

The fixture is built from metagraph's own test data with the `metagraph`
binary (taken from $METAGRAPH, else from PATH); the tests using it are skipped if
neither is available, unless PYMETAGRAPH_REQUIRE_METAGRAPH=1. Query results
are compared to `metagraph query`, alignments and mappings to `metagraph align`.
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
METAGRAPH = os.environ.get("METAGRAPH") or shutil.which("metagraph")

# CI sets PYMETAGRAPH_REQUIRE_METAGRAPH=1 so that a missing binary fails the
# tests instead of skipping them
if METAGRAPH is None and os.environ.get("PYMETAGRAPH_REQUIRE_METAGRAPH") == "1":
    raise RuntimeError("metagraph binary not found: set METAGRAPH or add it to PATH")


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
    # only the tests that need an index are skipped without the binary
    if METAGRAPH is None:
        pytest.skip("metagraph binary not found")
    out = tmp_path_factory.mktemp("index") / "graph"
    subprocess.run([METAGRAPH, "build", "-k", str(K), "-o", str(out), str(DATA)],
                   check=True, capture_output=True)
    subprocess.run([METAGRAPH, "annotate", "-i", f"{out}.dbg", "--anno-header",
                    "-o", str(out), str(DATA)], check=True, capture_output=True)
    return Path(f"{out}.dbg"), Path(f"{out}.column.annodbg")


@pytest.fixture(scope="module")
def graph(index):
    return pymetagraph.AnnotatedGraph(*index, threads=2)


def cli_query(index, fasta, discovery_fraction=None):
    """{sequence index: [(label, kmer_count)]} as reported by `metagraph query`,
    with its default options unless discovery_fraction is given."""
    options = [] if discovery_fraction is None else ["--min-kmers-fraction-label", str(discovery_fraction)]
    result = subprocess.run(
        [METAGRAPH, "query", "-i", str(index[0]), "-a", str(index[1]), *options, str(fasta)],
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


def write_queries(tmp_path):
    seqs = list(read_fasta(DATA).values())
    queries = [seqs[0], seqs[1][:100], seqs[2][50:150], seqs[2], "ACGT" * 30]
    fasta = tmp_path / "queries.fa"
    fasta.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(queries)))
    return queries, fasta


def test_defaults_match_cli(graph, index, tmp_path):
    queries, fasta = write_queries(tmp_path)
    expected = cli_query(index, fasta)
    result = graph.query(queries)
    assert any(len(hits) > 1 for hits in result)
    for i, hits in enumerate(result):
        assert sorted(hits) == sorted(expected.get(i, []))


def test_matches_cli(graph, index, tmp_path):
    queries, fasta = write_queries(tmp_path)

    for discovery_fraction in (0.0, 0.5, 1.0):
        expected = cli_query(index, fasta, discovery_fraction)
        result = graph.query(queries, discovery_fraction=discovery_fraction)
        for i, hits in enumerate(result):
            assert sorted(hits) == sorted(expected.get(i, []))
            # the CLI prints untruncated hits in label order, query sorts them
            counts = [count for _, count in hits]
            assert counts == sorted(counts, reverse=True)


def test_counts_and_order(graph):
    seq = next(iter(read_fasta(DATA).values()))
    hits = graph.query(seq)
    counts = [count for _, count in hits]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] == len(seq) - K + 1
    assert graph.query(seq, num_top_labels=1) == hits[:1]


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
    for fraction in (1.5, -0.1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            graph.query(["ACGT"], discovery_fraction=fraction)
        with pytest.raises(ValueError):
            graph.query("ACGT", presence_fraction=fraction)
        # the extension checks the same without the Python wrapper
        with pytest.raises(ValueError):
            graph._graph.query(["ACGT"], None, fraction, 0.0)


def mutate(seq, every=37):
    """seq with a substitution every |every| characters"""
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return "".join(swap.get(c, c) if i % every == every // 2 else c for i, c in enumerate(seq))


def write_align_queries(tmp_path):
    seqs = list(read_fasta(DATA).values())
    queries = [seqs[0], mutate(seqs[1]), seqs[2][30:230], "ACGT" * 30,
               mutate(seqs[0][:200], every=11)]
    # also the reverse complement of a reference sequence
    rc = {"A": "T", "C": "G", "G": "C", "T": "A"}
    queries.append("".join(rc.get(c, c) for c in reversed(seqs[1])))
    fasta = tmp_path / "align_queries.fa"
    fasta.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(queries)))
    return queries, fasta


def cli_align(index, fasta, labeled, options=()):
    """{sequence index: (strand, score, cigar, labels)} of the best alignment
    reported by `metagraph align`, None for unaligned sequences."""
    anno = ["-a", str(index[1])] if labeled else []
    result = subprocess.run(
        [METAGRAPH, "align", "-i", str(index[0]), *anno, *options, str(fasta)],
        check=True, capture_output=True, text=True)
    alignments = {}
    for line in result.stdout.splitlines():
        # header, query, then strand, path sequence, score, matches, cigar,
        # offset and, if labeled, labels
        fields = line.split("\t")
        if fields[2] == "*":
            alignments[int(fields[0])] = None
            continue
        # labels joined by ";" (which test labels contain as well)
        labels = fields[8] if labeled else ""
        alignments[int(fields[0])] = (fields[2], int(fields[4]), fields[6], labels)
    return alignments


@pytest.mark.parametrize("labeled", [True, False])
@pytest.mark.parametrize("min_exact_match", [None, 0.0, 0.3])
def test_align_matches_cli(graph, index, tmp_path, labeled, min_exact_match):
    queries, fasta = write_align_queries(tmp_path)
    options = [] if min_exact_match is None else ["--align-min-exact-match", str(min_exact_match)]
    expected = cli_align(index, fasta, labeled, options)
    result = graph.align(queries, labeled=labeled, min_exact_match=min_exact_match)
    assert len(result) == len(queries)
    assert sum(bool(alignments) for alignments in result) >= 4
    for i, alignments in enumerate(result):
        if expected[i] is None:
            assert alignments == []
            continue
        best = alignments[0]
        assert ("-" if best.orientation else "+", best.score, best.cigar,
                ";".join(best.labels)) == expected[i]
        assert bool(best.labels) == labeled
        assert 0 <= best.query_begin < best.query_end <= len(queries[i])
        assert best.score <= best.max_score


def test_align_masked_graph_matches_cli(tmp_path):
    """`metagraph align` drops the dummy k-mer mask before aligning, so a graph
    with a mask must align as without one."""
    if METAGRAPH is None:
        pytest.skip("metagraph binary not found")
    out = tmp_path / "masked"
    subprocess.run([METAGRAPH, "build", "-k", str(K), "--mask-dummy", "--in-ram",
                    "-o", str(out), str(DATA)], check=True, capture_output=True)
    assert Path(f"{out}.edgemask").exists()
    subprocess.run([METAGRAPH, "annotate", "-i", f"{out}.dbg", "--anno-header",
                    "-o", str(out), str(DATA)], check=True, capture_output=True)
    index = (Path(f"{out}.dbg"), Path(f"{out}.column.annodbg"))
    masked = pymetagraph.AnnotatedGraph(*index, threads=2)
    queries, fasta = write_align_queries(tmp_path)
    expected = cli_align(index, fasta, True, ["--align-min-exact-match", "0.0"])
    for i, alignments in enumerate(masked.align(queries, min_exact_match=0.0)):
        if expected[i] is None:
            assert alignments == []
            continue
        best = alignments[0]
        assert ("-" if best.orientation else "+", best.score, best.cigar,
                ";".join(best.labels)) == expected[i]
    # the label query is unaffected by dropping the mask
    queries, fasta = write_queries(tmp_path)
    expected = cli_query(index, fasta, 0.0)
    for i, hits in enumerate(masked.query(queries, discovery_fraction=0.0)):
        assert sorted(hits) == sorted(expected.get(i, []))


def test_align_alternatives_and_options(graph):
    seq = next(iter(read_fasta(DATA).values()))
    assert len(graph.align(seq, max_alternative_alignments=3)) >= 1
    exact = graph.align(seq)[0]
    assert exact.cigar == f"{len(seq)}="
    assert (exact.query_begin, exact.query_end) == (0, len(seq))
    assert exact.score == exact.max_score
    assert [a.cigar for a in graph.align(seq)] == [a.cigar for a in graph.align([seq])[0]]
    assert graph.align(["", "ACGT"]) == [[], []]
    assert graph.align([]) == []
    for kwargs in ({"min_exact_match": 1.5}, {"min_exact_match": float("nan")},
                   {"max_alternative_alignments": 0}, {"max_num_nodes_per_seq_char": 0}):
        with pytest.raises(ValueError):
            graph.align(seq, **kwargs)


def cli_map(index, fasta, options=()):
    """{sequence index: (strand, labels, matches, begin, end)} from
    `metagraph align --align-map-only`, None for unmapped sequences."""
    result = subprocess.run(
        [METAGRAPH, "align", "-i", str(index[0]), "-a", str(index[1]), "--align-map-only",
         *options, str(fasta)], check=True, capture_output=True, text=True)
    mappings = {}
    for line in result.stdout.splitlines():
        header, strand, labels, matches, begin, end = line.split("\t")
        mappings[int(header)] = None if strand == "*" else \
            (strand, labels, int(matches), int(begin), int(end))
    return mappings


@pytest.mark.parametrize("min_exact_match", [None, 0.0, 0.3])
def test_map_matches_cli(graph, index, tmp_path, min_exact_match):
    queries, fasta = write_align_queries(tmp_path)
    options = [] if min_exact_match is None else ["--align-min-exact-match", str(min_exact_match)]
    expected = cli_map(index, fasta, options)
    result = graph.map(queries, min_exact_match=min_exact_match)
    assert len(result) == len(queries)
    assert sum(m is not None for m in result) >= 4
    for i, mapping in enumerate(result):
        if expected[i] is None:
            assert mapping is None
            continue
        assert ("-" if mapping.orientation else "+", ";".join(mapping.labels),
                mapping.num_matches, mapping.query_begin, mapping.query_end) == expected[i]


def test_map_edge_cases(graph):
    seq = next(iter(read_fasta(DATA).values()))
    mapping = graph.map(seq)
    assert mapping.num_matches == len(seq)
    assert graph.map(["", "ACGT"]) == [None, None]
    assert graph.map([]) == []
    with pytest.raises(ValueError):
        graph.map(seq, min_exact_match=-1)


def test_worker_error_raises():
    with pytest.raises(RuntimeError, match="error in worker"):
        pymetagraph._pymetagraph_core._raise_in_worker(threads=2)
