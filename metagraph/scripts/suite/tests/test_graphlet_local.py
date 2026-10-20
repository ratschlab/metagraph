"""Offline tests of area G's local side: graphlet_doc (the suite's document digests) against the library's own
parser.body_digest, and run_graphlet_local (the W1/W2/W4/drift runner) on the library's fixtures: an op plan calling
every tool of GraphletTools with every paging mode and an error case per tool, paging assembly against the uncut
answer, hash-seed determinism across ways, a response's conversion against the same document from a file, the manifest
mode, the library export and the drift report. The tests import the library (metagraph/api/python of this checkout);
the modules under test never do (the runner's way processes do, by their own sys.path)."""
import copy
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))
DATA = os.path.join(LIB, "tests", "data", "traverse")
sys.path.insert(0, SUITE)

import graphlet_doc  # noqa: E402
import run_graphlet_local as rgl  # noqa: E402
import suite_common as sc  # noqa: E402

# the library (dataclass slots) needs Python 3.10+; the modules under test run on 3.8+
needs_library = unittest.skipIf(sys.version_info < (3, 10), "the traversal library needs Python 3.10 or later")


def library():
    """metagraph.traverse.parser of this checkout (tests only)."""
    if LIB not in sys.path:
        sys.path.insert(0, LIB)
    from metagraph.traverse import parser
    return parser


def read_json(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def responses():
    """(name, /traverse response of `detail: graphlet`) of every committed fixture that has one."""
    out = []
    for path in sorted(glob.glob(os.path.join(DATA, "documents", "*.graphlet.json"))
                       + glob.glob(os.path.join(DATA, "documents", "budgets", "*.graphlet.json"))
                       + glob.glob(os.path.join(DATA, "real", "mini_refseq", "*.graphlet.json.gz"))):
        resp = read_json(path)
        if any(isinstance(r, dict) and isinstance(r.get("graphlet"), str) for r in resp.get("results", [])):
            out.append((os.path.basename(path), resp))
    return out


def bodies():
    """(name, server body text) of every committed .mgt fixture and every response's graphlet."""
    out = []
    for path in sorted(glob.glob(os.path.join(DATA, "documents", "*.mgt"))
                       + glob.glob(os.path.join(DATA, "documents", "*", "*.mgt"))):
        with open(path, "r", encoding="utf-8", newline="") as fh:
            out.append((os.path.relpath(path, DATA), fh.read()))
    for name, resp in responses():
        for i, r in enumerate(resp["results"]):
            if isinstance(r.get("graphlet"), str):
                out.append(("%s[%d]" % (name, i), r["graphlet"]))
    return out


def with_seed_index(text, value):
    first = text.index("\n")
    head = text[:first].split(" ")
    head[graphlet_doc.H_SEED_INDEX] = str(value)
    return " ".join(head) + text[first:]


@needs_library
class GraphletDocTests(unittest.TestCase):
    def test_fixtures_present(self):
        self.assertGreaterEqual(len(responses()), 20)
        self.assertGreaterEqual(len(bodies()), 40)

    def test_unmasked_digest_is_the_librarys(self):
        parser = library()
        for name, body in bodies():
            with self.subTest(name):
                self.assertEqual(graphlet_doc.body_digest(body, mask_seed_index=False), parser.body_digest(body))

    def test_standalone_and_body_share_the_digest(self):
        parser = library()
        for name, resp in responses():
            for r in resp["results"]:
                if not isinstance(r.get("graphlet"), str):
                    continue
                with self.subTest(name):
                    text = parser.standalone_text(r["graphlet"], r, resp)
                    self.assertTrue(text.split("\n", 2)[1].startswith("J "))
                    self.assertEqual(graphlet_doc.body_digest(text, mask_seed_index=False),
                                     parser.body_digest(text))
                    self.assertEqual(graphlet_doc.body_digest(text), graphlet_doc.body_digest(r["graphlet"]))
                    d = graphlet_doc.digests(text)
                    self.assertEqual(d["seed_index"], 0)
                    self.assertEqual(d["bytes"], len(text.encode("utf-8")))
                    self.assertIsNotNone(d["envelope_digest"])

    def test_seed_index_masked(self):
        parser = library()
        name, body = bodies()[0]
        moved = with_seed_index(body, 3)
        self.assertNotEqual(parser.body_digest(moved), parser.body_digest(body))
        self.assertNotEqual(graphlet_doc.body_digest(moved, mask_seed_index=False),
                            graphlet_doc.body_digest(body, mask_seed_index=False))
        self.assertEqual(graphlet_doc.body_digest(moved), graphlet_doc.body_digest(body))
        self.assertEqual(graphlet_doc.frame(moved)["seed_index"], 3)

    def test_a_record_change_changes_the_digest(self):
        name, body = bodies()[0]
        lines = body.split("\n")
        lines[1] = lines[1] + "A"   # the S record's sequence, one base more
        self.assertNotEqual(graphlet_doc.body_digest("\n".join(lines)), graphlet_doc.body_digest(body))

    def test_envelope_digest(self):
        parser = library()
        name, resp = [x for x in responses() if x[0] == "fork.graphlet.json"][0]
        r = resp["results"][0]
        body = r["graphlet"]
        self.assertIsNone(graphlet_doc.envelope(body))
        self.assertIsNone(graphlet_doc.envelope_digest(body))
        base = graphlet_doc.envelope_digest(parser.standalone_text(body, r, resp))

        def digest_of(resp2, r2=None):
            r2 = r if r2 is None else r2
            return graphlet_doc.envelope_digest(parser.standalone_text(body, r2, resp2))
        # the per-request values: the response's timing and usage, the seed summary's timing
        self.assertEqual(digest_of(dict(resp, timing={"elapsed_ms": 12.5}, usage={"attempt_id": "x", "per_seed": [
            {"index": 0, "work": 3}, {"index": 1, "work": 4}]})), base)
        self.assertEqual(digest_of(resp, dict(r, timing={"serialize_ms": 1.0})), base)
        # E1: the annotation's lookahead counters differ between two engine runs of one request
        annotation = dict(r.get("annotation") or {}, direct_reads=7, keys_mapped=12345)
        self.assertEqual(digest_of(resp, dict(r, annotation=annotation)), base)
        self.assertNotEqual(digest_of(resp, dict(r, annotation=dict(annotation, rows_requested=1))), base)
        # VL16: the caller's seed name and its mismatch flag are request metadata; the validated id is not
        self.assertEqual(digest_of(resp, dict(r, seed=dict(r["seed"], seed_id="service-name"))), base)
        self.assertEqual(digest_of(resp, dict(r, seed=dict(r["seed"], seed_id="service-name",
                                                           seed_id_mismatch=True))), base)
        self.assertIn("seed_id_mismatch", r["seed"])
        self.assertNotEqual(digest_of(resp, dict(r, seed=dict(r["seed"], validated_seed_id="0" * 16))), base)
        strategy_ann = json.loads(json.dumps(resp["strategy"]))
        strategy_ann["annotation"] = {"batch_kmers": 7}
        self.assertNotEqual(digest_of(dict(resp, strategy=strategy_ann)), base)
        # everything else is compared: the strategy (its output.timing too), the capabilities, the summary
        strategy = json.loads(json.dumps(resp["strategy"]))
        strategy.setdefault("output", {})["timing"] = not strategy["output"].get("timing", False)
        self.assertNotEqual(digest_of(dict(resp, strategy=strategy)), base)
        caps = dict(resp["capabilities"], feature_level=99)
        self.assertNotEqual(digest_of(dict(resp, capabilities=caps)), base)
        self.assertNotEqual(digest_of(resp, dict(r, label_mode="annotate")), base)
        env = graphlet_doc.envelope(parser.standalone_text(body, r, resp))
        self.assertEqual(env["results"][0]["seed"], r["seed"])
        self.assertNotIn("graphlet", env["results"][0])

    def test_broken_frames_are_refused(self):
        name, body = bodies()[0]
        lines = body.split("\n")[:-1]
        cases = {
            "no final LF": body[:-1],
            "not MGT v1": body.replace("H mgt 1 ", "H mgt 2 ", 1),
            "H field count": body.replace("\n", " x\n", 1),
            "seed_index": with_seed_index(body, "01"),
            "Z count": "\n".join(lines[:-1] + ["Z %d" % (len(lines) + 1)]) + "\n",
            "no Z": "\n".join(lines[:-1]) + "\n",
            "J elsewhere": "\n".join(lines[:2] + ['J {}'] + lines[2:-1] + ["Z %d" % (len(lines) + 1)]) + "\n",
            "J not an object": "\n".join(lines[:1] + ["J []"] + lines[1:-1] + ["Z %d" % (len(lines) + 1)]) + "\n",
            "not UTF-8": body.encode("utf-8") + b"\xff",
            "not text": 5,
        }
        for case, text in cases.items():
            with self.subTest(case):
                with self.assertRaises(ValueError):
                    graphlet_doc.digests(text)
        # bytes are read as UTF-8
        self.assertEqual(graphlet_doc.body_digest(body.encode("utf-8")), graphlet_doc.body_digest(body))



# ============================================================================================== run_graphlet_local
REAL = os.path.join(DATA, "real", "mini_refseq")
DOCS = {
    "G-switch": os.path.join(REAL, "mini_ndm1_rc__switch1.graphlet.json.gz"),     # two arms, merges, cut lists
    "G-trace": os.path.join(REAL, "mini_ndm1__trace.graphlet.json.gz"),           # constrain; compare a
    "G-annot": os.path.join(REAL, "mini_ndm1__annotate_exh.graphlet.json.gz"),    # annotate, same seed; compare b
    "G-body": os.path.join(DATA, "documents", "fork.mgt"),                        # a bare body: no envelope
    "G-coords": os.path.join(REAL, "mini_rep_00__trace_coords.graphlet.json.gz"), # record coordinates
    "G-merge": os.path.join(DATA, "documents", "merge.graphlet.json"),            # one arm, two merges
}


def sel(name, **args):
    return {"select": name, "args": args}


def ref(op, path):
    return {"ref": op, "path": path}


def op(oid, tool, args, expect="ok", paging="none"):
    """An op of the plan; `expect` (the outcome: ok, an error code, or selector) is the test's, the runner ignores
    keys it does not read."""
    return {"id": oid, "tool": tool, "args": args, "paging": paging, "expect": expect}


def walks(arm, **kw):
    return dict({"arm": arm, "rank": "id", "n": 100, "min_bp": 0, "label": None, "spell": "none", "tail_bp": 60,
                 "route_consistent": True, "max_bytes": 262144}, **kw)


def labels(**kw):
    return dict({"arm": None, "rank": "direct_bp", "n": 100, "min_direct_bp": 0, "at_bp": None, "name": None,
                 "max_bytes": 262144, "coordinates": False}, **kw)


def seq(**kw):
    return dict({"arm": "left", "walk": None, "segment": None, "start": 0, "end": None, "orientation": "natural",
                 "with_seed": False}, **kw)


def export(**kw):
    return dict({"format": "fasta", "arm": None, "walks": None, "path": None, "what": None, "other": None}, **kw)


def fetch(**kw):
    return dict({"index": None, "seed": {"sequence": "ACGTACGTAC"}, "strategy": {}, "keep": True,
                 "max_graphlet_mb": 8, "replay": None, "allow_unverified_index": False, "max_bytes": 2048,
                 "budget": None, "coordinates": "auto"}, **kw)


def cont(**kw):
    return dict({"arm": "left", "walk": None, "overrides": None, "execute": False, "allow_unverified_index": False,
                 "max_bytes": 16384, "reset_branches": False, "budget": None}, **kw)


CHAIN = ["walk-first", "walk-3", "walk-9", "walk-last"]
MERGE = sel("walk_through_merge", walks=CHAIN)
CONT = sel("walk_with_continuation", walks=CHAIN, min_labels=1)
TOP = sel("top_label", labels="labels-tab", arm="left")
ENTERED = sel("label_entered_at_merge", walk="walk-merge")


def every_tool_plan():
    """Every tool of GraphletTools at least once, every paging mode, an error case per tool, every named selector on
    the library's real answers. Six documents: three enter as W1's (responses), G-body as a bare server body."""
    switch = [
        op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": 262144}),
        op("sum-cut", "graphlet_summary", {"arm": "left", "detail": "summary", "max_bytes": 900}),
        op("lims", "graphlet_summary", {"arm": None, "detail": "limitations", "max_bytes": 600}, paging="all"),
        op("lims-p1", "graphlet_summary", {"arm": None, "detail": "limitations", "max_bytes": 600},
           paging={"pages": 1}),
        op("sum-bad", "graphlet_summary", {"arm": None, "detail": "full", "max_bytes": 2048}, "bad_argument"),
        op("sum-arm", "graphlet_summary", {"arm": "up", "detail": "summary", "max_bytes": 2048}, "bad_arm"),
        op("labels-tab", "graphlet_labels", labels()),
        op("walks-id-l", "graphlet_walks", walks("left")),
        op("walks-len-l", "graphlet_walks", walks("left", rank="length", n=1)),
        op("walks-sup-paged", "graphlet_walks", walks("left", rank="support", spell="tail", max_bytes=1500),
           paging="all"),
        op("walks-p1", "graphlet_walks", walks("left", rank="support", spell="tail", max_bytes=1500),
           paging={"pages": 1}),
        op("walks-resume", "graphlet_walks", walks("left", rank="support", spell="tail", max_bytes=1500,
                                                   cursor={"cursor_of": "walks-p1"})),
        op("walks-reuse", "graphlet_walks", walks("left", rank="length", spell="tail", max_bytes=1500,
                                                  cursor={"cursor_of": "walks-p1"}), "bad_cursor"),
        op("walks-done", "graphlet_walks", walks("left", rank="length", cursor={"cursor_of": "walks-sup-paged"}),
           "selector"),
        op("walks-r", "graphlet_walks", walks("right", rank="loss", spell="full")),
        op("walks-label", "graphlet_walks", walks("left", label=TOP, route_consistent=False)),
        op("walks-bad-rank", "graphlet_walks", walks("left", rank="best"), "bad_argument"),
        op("walks-coords", "graphlet_walks", walks("left", coordinates=True), "no_coordinates"),
        op("walks-budget", "graphlet_walks", walks("left", budget={"work_units": 1000}), "bad_argument"),
        op("walk-first", "graphlet_walk", {"arm": "left", "walk": ref("walks-id-l", "/rows/0/walk"),
                                           "max_bytes": 262144}),
        op("walk-last", "graphlet_walk", {"arm": "left", "walk": ref("walks-id-l", "/rows/11/walk"),
                                          "max_bytes": 262144}),
        op("walk-longest", "graphlet_walk", {"arm": "left", "walk": ref("walks-len-l", "/rows/0/walk"),
                                             "max_bytes": 262144}),
        op("walk-3", "graphlet_walk", {"arm": "left", "walk": 3, "max_bytes": 262144}),
        op("walk-9", "graphlet_walk", {"arm": "left", "walk": 9, "max_bytes": 262144}),
        op("walk-merge", "graphlet_walk", {"arm": "left", "walk": MERGE, "max_bytes": 262144}),
        op("walk-cont", "graphlet_walk", {"arm": "left", "walk": CONT, "max_bytes": 262144}),
        op("walk-paged", "graphlet_walk", {"arm": "left", "walk": 9, "max_bytes": 1500}, paging="all"),
        op("walk-out", "graphlet_walk", {"arm": "left", "walk": ref("walks-id-l", "/total"), "max_bytes": 2048},
           "bad_argument"),
        op("walk-sel", "graphlet_walk", {"arm": "left", "walk": ref("walks-id-l", "/rows/40/walk"),
                                         "max_bytes": 2048}, "selector"),
        op("walk-sel-name", "graphlet_walk", {"arm": "left", "walk": {"select": ["walk_through_merge"]},
                                              "max_bytes": 2048}, "selector"),
        op("walks-label-merge", "graphlet_walks", walks("left", label=ENTERED)),
        op("support-merge", "graphlet_support", {"arm": "left", "walk": MERGE, "step": None, "max_bytes": 262144}),
        op("support-step", "graphlet_support", {"arm": "left", "walk": MERGE, "step": 50, "max_bytes": 1500},
           paging="all"),
        op("support-bad", "graphlet_support", {"arm": "middle", "walk": 0, "step": None, "max_bytes": 2048},
           "bad_arm"),
        op("labels-reach", "graphlet_labels", labels(arm="left", rank="reach_bp", max_bytes=1200), paging="all"),
        op("labels-at", "graphlet_labels", labels(arm="left", min_direct_bp=10, at_bp=100)),
        op("labels-name", "graphlet_labels", labels(name=TOP)),
        op("labels-name-merge", "graphlet_labels", labels(arm="left", name=ENTERED)),
        op("labels-coords", "graphlet_labels", labels(coordinates=True), "bad_argument"),
        op("labels-unknown", "graphlet_labels", labels(name={"name": "no such label"}), "unknown_label"),
        op("splits-l", "graphlet_splits", {"arm": "left", "n": 100, "min_labels_before": 0, "max_bytes": 262144}),
        op("splits-paged", "graphlet_splits", {"arm": "left", "n": 100, "min_labels_before": 0, "max_bytes": 1000},
           paging="all"),
        op("splits-bad", "graphlet_splits", {"arm": "left", "n": 0, "min_labels_before": 0, "max_bytes": 2048},
           "bad_argument"),
        op("claims-l", "graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": True, "labels": None,
                                           "max_bytes": 262144}),
        op("claims-paged", "graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": False,
                                               "labels": None, "max_bytes": 1500}, paging="all"),
        op("claims-labels", "graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": True,
                                                "labels": [TOP], "max_bytes": 262144}),
        op("claims-unknown", "graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": True,
                                                 "labels": [{"name": "no such label"}], "max_bytes": 2048},
           "unknown_label"),
        op("seq-walk", "graphlet_sequence", seq(walk=MERGE), paging="all"),
        op("seq-seed", "graphlet_sequence", seq(walk=MERGE, with_seed=True, orientation="walk"), paging="all"),
        op("seq-seg", "graphlet_sequence", seq(segment=3, orientation="walk")),
        op("seq-past-end", "graphlet_sequence", seq(
            walk=MERGE, start=sel("slice_past_end", length_of="seq-walk", part="start", keep=20),
            end=sel("slice_past_end", length_of="seq-walk", part="end", past=50))),
        op("seq-p2", "graphlet_sequence", seq(walk=MERGE, with_seed=True), paging={"pages": 2}),
        op("seq-both", "graphlet_sequence", seq(walk=0, segment=0), "bad_argument"),
        op("seq-seg-oor", "graphlet_sequence", seq(segment=999), "bad_argument"),
        op("seq-seed-seg", "graphlet_sequence", seq(segment=0, with_seed=True), "bad_argument"),
        op("seq-beyond", "graphlet_sequence", seq(walk=0, start=100000), "bad_argument"),
        op("exp-fasta", "graphlet_export", export()),
        op("exp-fasta-walks", "graphlet_export", export(arm="left", walks=[0, MERGE], path="G-switch.left.fa")),
        op("exp-gfa", "graphlet_export", export(format="gfa")),
        op("exp-json", "graphlet_export", export(format="json")),
        op("exp-json-walk", "graphlet_export", export(format="json", what="walk", arm="left", walks=[MERGE],
                                                      path="walk.json")),
        op("exp-json-compare", "graphlet_export", export(format="json", what="compare",
                                                         other={"graphlet": "G-switch"}, path="cmp.json")),
        op("exp-mgt", "graphlet_export", export(format="mgt")),
        op("exp-bad-what", "graphlet_export", export(what="walk"), "bad_argument"),
        op("exp-path", "graphlet_export", export(path="../escape.fa"), "path_not_allowed"),
        op("sub", "graphlet_subtrie", {"labels": [TOP], "arm": None, "mode": "any"}),
        op("sub-all", "graphlet_subtrie", {"labels": [TOP, sel("top_label", labels="labels-tab", arm="left",
                                                               index=1)], "arm": "left", "mode": "all"}),
        op("sub-walks", "graphlet_walks", walks("left", handle=ref("sub", "/handle"))),
        op("sub-summary", "graphlet_summary", {"handle": ref("sub", "/handle"), "arm": None, "detail": "summary",
                                               "max_bytes": 262144}),
        op("sub-cmp", "graphlet_compare", {"a": ref("sub", "/handle"), "b": {"graphlet": "G-switch"}, "arm": None,
                                           "mode": "claims", "max_bytes": 262144}, "view_unsupported"),
        op("sub-sub", "graphlet_subtrie", {"handle": ref("sub", "/handle"), "labels": [TOP], "arm": None,
                                           "mode": "any"}, "view_unsupported"),
        op("sub-bad", "graphlet_subtrie", {"labels": [], "arm": None, "mode": "any"}, "bad_argument"),
        op("cmp-self", "graphlet_compare", {"b": {"graphlet": "G-switch"}, "arm": None, "mode": "walks",
                                            "max_bytes": 262144}),
        op("save", "graphlet_save", {"path": "G-switch.saved.mgt"}),
        op("save-bad", "graphlet_save", {"path": "../saved.mgt"}, "path_not_allowed"),
        op("load", "graphlet_load", {"path": "G-switch.saved.mgt"}),
        op("load-summary", "graphlet_summary", {"handle": ref("load", "/handle"), "arm": "left",
                                                "detail": "summary", "max_bytes": 262144}),
        op("load-bad", "graphlet_load", {"path": "no-such-file.mgt"}, "io_error"),
        op("list", "graphlet_list", {"max_bytes": 262144}),
        op("list-paged", "graphlet_list", {"max_bytes": 900}, paging="all"),
        op("list-bad", "graphlet_list", {"cursor": "junk", "max_bytes": 2048}, "bad_cursor"),
        op("free", "graphlet_free", {"handle": ref("load", "/handle")}),
        op("free-again", "graphlet_free", {"handle": ref("load", "/handle")}, "unknown_handle"),
        op("caps", "traverse_capabilities", {"index": None, "max_bytes": 40000}, "unknown_index"),
        op("resolve", "traverse_resolve", {"index": None, "sequence": "ACGTACGT", "labels": None, "discover": None,
                                           "select": None, "limit": 10, "max_bytes": 2048}, "unknown_index"),
        op("fetch", "traverse_fetch", fetch(), "unknown_index"),
        op("fetch-bad", "traverse_fetch", fetch(seed=None), "bad_request"),
        op("cont", "traverse_continue", cont(walk=CONT)),
        op("cont-clash", "traverse_continue", cont(walk=CONT, overrides={"arm": "right"}), "bad_argument"),
        op("nope", "graphlet_nope", {}, "unknown_tool"),
    ]
    trace = [
        op("cmp-claims", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "claims",
                                              "max_bytes": 262144}),
        op("cmp-walks", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": "left", "mode": "walks",
                                             "max_bytes": 262144}),
        op("cmp-labels", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "labels",
                                              "max_bytes": 262144}),
        op("cmp-prefix", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "prefix_subset",
                                              "max_bytes": 262144}),
        op("cmp-paged", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "labels",
                                             "max_bytes": 1500}, paging="all"),
        op("cmp-bad", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "all",
                                           "max_bytes": 2048}, "bad_argument"),
        op("cmp-nodoc", "graphlet_compare", {"b": {"graphlet": "G-missing"}, "arm": None, "mode": "claims",
                                             "max_bytes": 2048}, "selector"),
    ]
    body = [
        op("body-sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": 2048}, "missing_envelope"),
        op("body-walks", "graphlet_walks", walks("left")),
    ]
    coords = [
        op("coords-tab", "graphlet_labels", labels(arm="left", n=5)),
        op("coords-walks", "graphlet_walks", walks("left", coordinates=True)),
        op("coords-name", "graphlet_labels", labels(arm="left", name=sel("top_label", labels="coords-tab"),
                                                    coordinates=True)),
    ]
    merge = [
        op("m-sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": 262144}),
        op("m-walk", "graphlet_walk", {"arm": None, "walk": 0, "max_bytes": 262144}),
        op("m-first", "graphlet_walks", walks(sel("first_arm", of="m-sum"))),
        op("m-second", "graphlet_labels", labels(arm=sel("second_arm_or_none", of="m-sum"))),
        op("m-entered", "graphlet_walks", walks(None, label=sel("label_entered_at_merge", walk="m-walk"))),
    ]
    return {"scenarios": [{"id": sid, "graphlet": {"ops": ops}} for sid, ops in (
        ("G-switch", switch), ("G-trace", trace), ("G-annot", []), ("G-body", body), ("G-coords", coords),
        ("G-merge", merge))]}


def write_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    return path


def run_offline(out, plan, docs=None, ways="W4", limits=None, hashseed=None, extra=()):
    """run_graphlet_local offline over |plan| and |docs| ({id: path}) into |out| -> {way: [records]}."""
    os.makedirs(out, exist_ok=True)
    argv = ["--plan", write_json(os.path.join(out, "plan.json"), plan), "--library-root", LIB, "--out", out,
            "--ways", ways, "--run-id", "test-run", "--quiet"]
    for sid, path in sorted((docs or DOCS).items()):
        argv += ["--doc", "%s=%s" % (sid, path)]
    if limits is not None:
        argv += ["--limits", write_json(os.path.join(out, "limits.json"), limits)]
    if hashseed:
        argv += ["--hashseed", hashseed]
    rc = rgl.main(argv + list(extra))
    assert rc == 0, rc
    return dict((w, records(out, w)) for w in ways.split(","))


def records(out, way, name="graphlet_ops.jsonl"):
    rows, bad = sc.read_jsonl(os.path.join(out, "arms", rgl.WAY_DIRS[way], name))
    assert bad == 0
    return rows


def outcome(rec):
    if rec["selector_error"]:
        return "selector"
    return rec["error"]["code"] if rec["error"] else "ok"


VOLATILE = ("elapsed_ms", "started_utc")


def comparable(rec, drop_way=False):
    """A record without its timing, hash seed (and way): what two processes must agree on. graphlet_list's rows carry
    the store's wall-clock created / accessed times, timing too."""
    r = copy.deepcopy(rec)
    for k in VOLATILE:
        r.pop(k, None)
    r["library"].pop("hashseed", None)
    if drop_way:
        r.pop("way", None)
    if r["tool"] == "graphlet_list" and r["answer"]:
        for row in r["answer"]["rows"]:
            row.pop("created", None)
            row.pop("accessed", None)
    return r


def library_tool_names():
    if LIB not in sys.path:
        sys.path.insert(0, LIB)
    from metagraph.traverse.mcp_tools import tool_names
    return set(tool_names())


@needs_library
class ToolTableTests(unittest.TestCase):
    def test_every_tool_is_placed(self):
        """Every library tool either reads a graphlet (HANDLE_ARG), makes a handle (MINTING), or neither (list, load,
        the backend's capabilities and resolve)."""
        self.assertTrue(set(rgl.HANDLE_ARG) | set(rgl.MINTING) | {"graphlet_list", "graphlet_load"}
                        >= library_tool_names() - {"traverse_capabilities", "traverse_resolve"})


@needs_library
class EveryToolTests(unittest.TestCase):
    """The every-tool plan on two ways with different hash seeds (W1, W2) in one invocation, sequence pages of 1 KiB."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        cls.plan = every_tool_plan()
        cls.out = os.path.join(cls.tmp, "run")
        # max_graphlet_mb: accepted, not applied (VL10)
        cls.recs = run_offline(cls.out, cls.plan, ways="W1,W2", limits={"sequence_max_bytes": 1024,
                                                                        "max_graphlet_mb": 8},
                               hashseed="W1=11,W2=12345")
        cls.by_op = dict((r["op"], r) for r in cls.recs["W1"])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_outcomes_as_planned(self):
        planned = dict((o["id"], o["expect"]) for s in self.plan["scenarios"] for o in s["graphlet"]["ops"])
        self.assertEqual(len(self.recs["W1"]), len(planned))
        for r in self.recs["W1"]:
            with self.subTest(r["op"]):
                self.assertEqual(outcome(r), planned[r["op"]], r["selector_error"] or r["error"])

    def test_every_tool_called_and_refused(self):
        names = library_tool_names()
        self.assertEqual(len(names), 19)
        called = set(r["tool"] for r in self.recs["W1"])
        self.assertLessEqual(names, called)
        answered = set(r["tool"] for r in self.recs["W1"] if r["answer"] is not None)
        refused = set(r["tool"] for r in self.recs["W1"] if r["error"] is not None)
        self.assertEqual(names - answered, {"traverse_capabilities", "traverse_resolve", "traverse_fetch"})
        self.assertLessEqual(names, refused)

    def test_paging_modes(self):
        by = self.by_op
        for name in ("lims", "walks-sup-paged", "walk-paged", "support-step", "labels-reach", "splits-paged",
                     "claims-paged", "list-paged", "seq-walk", "seq-seed"):
            with self.subTest(name):
                self.assertGreater(by[name]["pages"], 1)
                self.assertTrue(by[name]["paging_complete"])
        for name, n in (("lims-p1", 1), ("seq-p2", 2)):
            with self.subTest(name):
                self.assertEqual(by[name]["pages"], n)
                self.assertFalse(by[name]["paging_complete"])
        self.assertEqual(len(by["lims-p1"]["answer"]["rows"]) < by["lims-p1"]["answer"]["total"], True)
        self.assertTrue(by["seq-p2"]["answer"]["truncated"])
        self.assertEqual(by["lims"]["answer"]["rows"][:len(by["lims-p1"]["answer"]["rows"])],
                         by["lims-p1"]["answer"]["rows"])
        self.assertEqual(by["sum"]["pages"], 1)
        # a cursor an earlier op stopped at: selectors read it, args_resolved shows <cursor>
        n = len(by["walks-p1"]["answer"]["rows"])
        self.assertFalse(by["walks-p1"]["paging_complete"])
        self.assertEqual(by["walks-resume"]["args_resolved"]["cursor"], "<cursor>")
        self.assertEqual(by["walks-p1"]["answer"]["rows"] + by["walks-resume"]["answer"]["rows"],
                         by["walks-sup-paged"]["answer"]["rows"][:n + len(by["walks-resume"]["answer"]["rows"])])
        self.assertEqual(by["cmp-paged"]["pages"] > 1, True)
        for r in self.recs["W1"]:
            if r["answer"] is not None:
                self.assertNotIn("next_cursor", r["answer"])
                self.assertNotIn("_next_cursor", r["answer"])
                self.assertNotIn("next_from", r["answer"])
        # every page is kept raw, numbered from 0
        raw = os.path.join(self.out, "arms", "dl-lib", "raw")
        n = by["claims-paged"]["pages"]
        for i in range(n):
            page = sc.read_gz_json(os.path.join(raw, "G-switch.claims-paged.%d.json.gz" % i))
            self.assertEqual((page["page"], page["op"], page["reply"]["ok"]), (i, "claims-paged", True))
            self.assertEqual("cursor" in page["args"], i > 0)
        self.assertFalse(os.path.exists(os.path.join(raw, "G-switch.claims-paged.%d.json.gz" % n)))

    def test_selectors_on_real_answers(self):
        by = self.by_op
        merge_walk = by["walk-merge"]["args_resolved"]["walk"]
        self.assertTrue(any("merge_of" in row for row in by["walk-merge"]["answer"]["rows"]))
        self.assertEqual(merge_walk, 3)   # walk 3 is the lowest named walk through a merge (segment 11)
        self.assertEqual(by["walk-cont"]["args_resolved"]["walk"], 9)
        self.assertGreaterEqual(by["walk-cont"]["answer"]["continuation"]["labels"], 1)
        self.assertEqual(by["walk-first"]["args_resolved"]["walk"], 0)
        self.assertEqual(by["walk-out"]["args_resolved"]["walk"], 12)
        self.assertIn("index 40", by["walk-sel"]["selector_error"])
        self.assertIn("followed to its end", by["walks-done"]["selector_error"])
        self.assertTrue(by["walks-done"]["selector_error"].startswith("args.cursor: "))
        self.assertTrue(by["walk-sel"]["selector_error"].startswith("args.walk: "))
        self.assertIn("named by a string", by["walk-sel-name"]["selector_error"])   # VL8, through the runner
        top = by["labels-tab"]["answer"]["rows"]
        left = [r for r in top if r["arm"] == "left"]
        self.assertEqual(by["walks-label"]["args_resolved"]["label"], {"ref": left[0]["ref"]})
        self.assertEqual(by["sub-all"]["args_resolved"]["labels"], [{"ref": left[0]["ref"]}, {"ref": left[1]["ref"]}])
        entered = by["labels-name-merge"]["args_resolved"]["name"]
        self.assertEqual(entered, by["walks-label-merge"]["args_resolved"]["label"])
        self.assertEqual(by["labels-name-merge"]["answer"]["ref"], entered["ref"])
        length = by["seq-walk"]["answer"]["length"]
        past = by["seq-past-end"]
        self.assertEqual((past["args_resolved"]["start"], past["args_resolved"]["end"]), (length - 20, length + 50))
        self.assertEqual((past["answer"]["from"], past["answer"]["to"], len(past["answer"]["sequence"])),
                         (length - 20, length, 20))
        self.assertEqual(by["seq-walk"]["answer"]["sequence"][-20:], past["answer"]["sequence"])
        self.assertEqual(len(by["seq-walk"]["answer"]["sequence"]), length)
        merge = dict((r["op"], r) for r in self.recs["W1"] if r["scenario"] == "G-merge")
        self.assertEqual(merge["m-first"]["args_resolved"]["arm"], "right")
        self.assertIsNone(merge["m-second"]["args_resolved"]["arm"])
        self.assertEqual(merge["m-entered"]["args_resolved"]["label"], {"ref": "c:1"})

    def test_handles_as_symbols(self):
        by = self.by_op
        self.assertEqual(by["sum"]["args_resolved"]["handle"], "<handle>")
        self.assertEqual(by["sum"]["answer"]["handle"], "<handle>")
        self.assertEqual(by["sub"]["answer"]["handle"], "<handle:sub>")
        self.assertEqual(by["sub"]["answer"]["of"], "<handle>")
        self.assertEqual(by["sub-walks"]["args_resolved"]["handle"], "<handle:sub>")
        self.assertEqual(by["sub-walks"]["answer"]["evidence"]["view"]["of"], "<handle>")
        self.assertEqual(by["sub-cmp"]["error"], {"code": "view_unsupported", "fields": {"of": "<handle>"}})
        self.assertEqual(by["load-summary"]["args_resolved"]["handle"], "<handle:load>")
        self.assertEqual(by["free"]["answer"], {"freed": "<handle:load>"})
        self.assertEqual(by["cmp-self"]["args_resolved"]["b"], "<handle>")
        self.assertEqual(by["cmp-self"]["args_resolved"]["a"], "<handle>")
        self.assertEqual(by["exp-fasta"]["answer"]["path"], "<export_dir>/<handle>.fa")
        self.assertEqual(by["exp-json-compare"]["args_resolved"]["other"], "<handle>")
        cmp = dict((r["op"], r) for r in self.recs["W1"] if r["scenario"] == "G-trace")
        self.assertEqual(cmp["cmp-claims"]["args_resolved"]["b"], "<handle:G-annot>")
        self.assertEqual(cmp["cmp-claims"]["answer"]["b"], "<handle:G-annot>")
        self.assertIn("args.b", cmp["cmp-nodoc"]["selector_error"])
        # the store when G-switch lists it (its scenario runs first): its graphlet and what its ops made, oldest first
        self.assertEqual([row["handle"] for row in by["list"]["answer"]["rows"]],
                         ["<handle>", "<handle:sub>", "<handle:sub-all>", "<handle:load>"])
        with open(os.path.join(self.out, "arms", "dl-lib", "graphlet_ops.jsonl"), encoding="utf-8") as fh:
            self.assertIsNone(re.search(r"g_[0-9a-f]{12}", fh.read()))

    def test_files_kept(self):
        """FORMAT v2 B3: {name, bytes, sha256, raw}; the text kept at raw with handles as symbols."""
        names = {"exp-fasta": "<handle>.fa", "exp-fasta-walks": "G-switch.left.fa", "exp-gfa": "<handle>.gfa",
                 "exp-json": "<handle>.json", "exp-json-walk": "walk.json", "exp-json-compare": "cmp.json",
                 "exp-mgt": "<handle>.mgt", "save": "G-switch.saved.mgt"}
        for op_id, name in names.items():
            with self.subTest(op_id):
                f = self.by_op[op_id]["file"]
                self.assertEqual(set(f), {"name", "bytes", "sha256", "raw"})
                self.assertEqual(f["name"], name)
                self.assertTrue(f["raw"].startswith("raw/G-switch.%s.file." % op_id))
                with gzip.open(os.path.join(self.out, "arms", "dl-lib", f["raw"]), "rb") as fh:
                    data = fh.read()
                self.assertEqual(hashlib.sha256(data).hexdigest(), f["sha256"])
                self.assertEqual(f["bytes"], self.by_op[op_id]["answer"]["bytes"])
                self.assertIsNone(re.search(rb"g_[0-9a-f]{12}", data))
                if op_id != "exp-json-walk":   # the one file here that names a handle: kept as its symbol
                    self.assertEqual(len(data), f["bytes"])
        with gzip.open(os.path.join(self.out, "arms", "dl-lib", self.by_op["exp-json-walk"]["file"]["raw"])) as fh:
            self.assertEqual(json.loads(fh.read().decode("utf-8"))["handle"], "<handle>")
        with gzip.open(os.path.join(self.out, "arms", "dl-lib", self.by_op["exp-mgt"]["file"]["raw"])) as fh:
            mgt = fh.read().decode("utf-8")
        self.assertEqual(graphlet_doc.digests(mgt)["body_digest"], self.by_op["sum"]["doc"]["body_digest"])

    def test_record_shape(self):
        keys = {"run_id", "scenario", "way", "op", "tool", "args_resolved", "selector_error", "pages",
                "paging_complete", "answer", "error", "service_fields", "doc", "library", "elapsed_ms", "started_utc"}
        for r in self.recs["W1"]:
            self.assertEqual(set(r) - {"file"}, keys)
            self.assertEqual(r["way"], "W1")
            self.assertIsNone(r["service_fields"])
            self.assertEqual(set(r["library"]), {"commit", "python", "hashseed", "work_model", "limits"})
            self.assertEqual(r["library"]["hashseed"], 11)
            self.assertEqual(r["library"]["limits"]["sequence_max_bytes"], 1024)
            self.assertEqual(set(r["library"]["limits"]), {"max_bytes", "sequence_max_bytes", "local_limits"})
            self.assertTrue(r["library"]["commit"])
            self.assertEqual((r["answer"] is None), bool(r["error"] or r["selector_error"]))
            if r["error"]:
                self.assertEqual(set(r["error"]), {"code", "fields"})
                self.assertNotIn("message", r["error"]["fields"])
        docs = dict((d["scenario"], d) for d in records(self.out, "W1", "graphlet_docs.jsonl"))
        self.assertEqual(set(docs), set(DOCS))
        self.assertIsNone(docs["G-body"]["doc"]["envelope_digest"])
        for sid, d in docs.items():
            self.assertIsNone(d["error"], sid)
            self.assertEqual(d["load"]["handle"], "<handle>")
            with gzip.open(os.path.join(self.out, "arms", "dl-lib", "docs", sid + ".mgt.gz"), "rb") as fh:
                self.assertEqual(graphlet_doc.digests(fh.read())["body_digest"], d["doc"]["body_digest"])

    def test_hashseed_determinism(self):
        """W1 and W2 run the same documents in two processes with different PYTHONHASHSEED: their records are equal
        but for the way, the hash seed and the timing fields."""
        a, b = self.recs["W1"], self.recs["W2"]
        self.assertEqual(len(a), len(b))
        self.assertEqual([comparable(r, True) for r in a], [comparable(r, True) for r in b])
        self.assertEqual(set(r["library"]["hashseed"] for r in b), {12345})
        runs = [sc.read_json(os.path.join(self.out, "arms", d, "run.json")) for d in ("dl-lib", "dr-lib")]
        self.assertEqual([r["interpreter"]["hashseed"] for r in runs], ["11", "12345"])
        self.assertNotEqual(runs[0]["interpreter"]["hash_probe"], runs[1]["interpreter"]["hash_probe"])
        for run in runs:   # VL7: the probe equals a separate launch's with that seed, the environment was read
            self.assertEqual(run["interpreter"]["hash_probe"], run["interpreter"]["hash_probe_launch"])
            self.assertEqual(run["interpreter"]["ignore_environment"], 0)
        self.assertEqual(runs[0]["child"]["env_keys"].count("PYTHONHASHSEED"), 1)
        self.assertFalse([k for k in runs[0]["child"]["env_keys"] if k.startswith("PYTHON") and k != "PYTHONHASHSEED"])
        self.assertEqual(runs[0]["state"], "done")
        self.assertEqual(runs[0]["limits_not_applied"]["max_graphlet_mb"]["value"], 8)
        # no recorded capabilities, the table states no graphlet block: --limits / the defaults, recorded
        self.assertEqual(runs[0]["limits_source"]["missing"],
                         ["service_settings.graphlet.library (absent: --limits or the defaults used)"])
        for d in ("dl-lib", "dr-lib"):   # VL15: the library's scratch goes with the run (--keep-work keeps it)
            self.assertFalse(os.path.exists(os.path.join(self.out, "arms", d, "work")))
        for way in ("W1", "W2"):
            files = sorted(os.listdir(os.path.join(self.out, "arms", rgl.WAY_DIRS[way], "raw")))
            self.assertEqual(len(files), len(os.listdir(os.path.join(self.out, "arms", "dl-lib", "raw"))))

    def test_no_append(self):
        with self.assertRaises(SystemExit):
            run_offline(self.out, self.plan, ways="W1")


@needs_library
class DocumentErrorTests(unittest.TestCase):
    def test_missing_failed_and_broken_documents(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            failed = write_json(os.path.join(tmp, "failed.json"), {"results": [{"seed": {}, "error": "no labels"}]})
            with open(DOCS["G-body"], encoding="utf-8") as fh:
                cut = fh.read()[:-40]
            broken = os.path.join(tmp, "broken.mgt")
            with open(broken, "w", encoding="utf-8") as fh:
                fh.write(cut)
            # VL2: a direct arm's raw records of a transport failure, an HTTP error with an error body, a 200 answer
            # whose seed failed, and a 200 answer that is no traverse response
            seed_failed = {"results": [{"seed": {}, "error": "no labels"}]}
            raws = {"G-transport": {"http": 0, "body": None, "error": {"code": "timeout", "text": "timed out"}},
                    "G-http": {"http": 503, "body": {"error": "server_busy"},
                               "error": {"code": None, "text": "server_busy"}},
                    "G-seed": {"http": 200, "body": seed_failed, "error": None},
                    "G-empty": {"http": 200, "body": {"results": [{"seed": {}}]}, "error": None},
                    "G-noresults": {"http": 200, "body": {"status": "ok"}, "error": None}}
            for sid, answer in raws.items():
                sc.write_gz_json(os.path.join(tmp, sid + ".json.gz"), {"scenario": sid, "answer": answer})
            docs = {"G-nofile": os.path.join(tmp, "absent.mgt.gz"), "G-failed": failed, "G-broken": broken,
                    "G-annot": DOCS["G-annot"]}
            docs.update((sid, os.path.join(tmp, sid + ".json.gz")) for sid in raws)
            ops = [op("s", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": 2048}),
                   op("c", "graphlet_compare", {"a": {"graphlet": "G-annot"}, "b": "<handle>", "arm": None,
                                                "mode": "labels", "max_bytes": 2048}),
                   # ops whose arguments select from an earlier op's answer (async's review X4)
                   op("w", "graphlet_walks", walks("left", n=1)),
                   op("w0", "graphlet_walk", {"arm": "left", "walk": {"ref": "w", "path": "/rows/0/walk"},
                                              "max_bytes": 2048}),
                   op("cur", "graphlet_walks", walks("left", n=1, cursor={"cursor_of": "w"}))]
            plan = {"scenarios": [{"id": sid, "graphlet": {"ops": ops}} for sid in (
                "G-nofile", "G-failed", "G-broken", "G-nodoc") + tuple(sorted(raws)) + ("G-annot",)]}
            recs = run_offline(os.path.join(tmp, "run"), plan, docs=docs, ways="W4")["W4"]
            got = dict(((r["scenario"], r["op"]), r) for r in recs)
            self.assertEqual(got[("G-nofile", "s")]["error"]["code"], "no_document")
            self.assertEqual(got[("G-nodoc", "s")]["error"], {"code": "no_document", "fields": {"source": None}})
            self.assertEqual(got[("G-transport", "s")]["error"], {"code": "no_document", "fields": {
                "source": os.path.join(tmp, "G-transport.json.gz"), "stage": "fetch", "code": "timeout"}})
            self.assertEqual(got[("G-http", "s")]["error"], {"code": "document_refused", "fields": {
                "stage": "fetch", "http": 503, "code": None}})
            self.assertEqual(got[("G-seed", "s")]["error"], {"code": "no_graphlet", "fields": {}})
            self.assertEqual(got[("G-empty", "s")]["error"], {"code": "document_refused", "fields": {
                "stage": "convert", "reason": "no_graphlet_no_error"}})
            self.assertEqual(got[("G-noresults", "s")]["error"], {"code": "document_refused", "fields": {
                "stage": "convert", "reason": "no_results"}})
            self.assertEqual(got[("G-failed", "s")]["error"], {"code": "no_graphlet", "fields": {}})
            self.assertEqual(got[("G-broken", "s")]["error"]["code"], "document_refused")
            self.assertEqual(got[("G-broken", "s")]["error"]["fields"], {"stage": "load", "code": "format_error"})
            # X4: without a document a dependent op carries the document's error (no selector error), so a missing
            # MCP document reads as the document's failure on W4, never as a W1-W4 disagreement
            for sid in ("G-nofile", "G-failed", "G-broken", "G-nodoc") + tuple(sorted(raws)):
                for oid in ("w", "w0", "cur"):
                    self.assertIsNone(got[(sid, oid)]["selector_error"], (sid, oid))
                    self.assertEqual(got[(sid, oid)]["error"], got[(sid, "s")]["error"], (sid, oid))
            # with a document the selectors resolve (a real selector error stays one)
            self.assertIsNone(got[("G-annot", "w0")]["selector_error"])
            self.assertIsNone(got[("G-annot", "w0")]["error"])
            # a call on another graphlet still runs; <handle> of a scenario without one is left to the library
            self.assertEqual(got[("G-failed", "c")]["args_resolved"]["a"], "<handle:G-annot>")
            self.assertEqual(got[("G-failed", "c")]["error"]["code"], "unknown_handle")
            docs_lines = dict((d["scenario"], d) for d in records(os.path.join(tmp, "run"), "W4",
                                                                   "graphlet_docs.jsonl"))
            self.assertIn("no labels", docs_lines["G-failed"]["error_text"])
            self.assertIsNotNone(docs_lines["G-broken"]["doc_error"])
            self.assertIsNone(docs_lines["G-broken"]["doc"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


FAKE_LIB = os.path.join(HERE, "data", "fake_graphlet_library")


class SupervisionTests(unittest.TestCase):
    """VERIFY-L VL3: the way's process supervised, with a fake library whose graphlet_walks(mode=...) hangs, answers
    what is no JSON, ends the process, repeats a cursor, prints, or reads stdin: each op gets its record, the process
    is restarted after a hang or a death with its documents reloaded, and later ops and ways still run."""

    def test_misbehaving_library(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            def walks_op(oid, mode, expect="ok", **kw):
                o = op(oid, "graphlet_walks", {"arm": "left", "mode": mode, "max_bytes": 2048}, expect, **kw)
                return o
            ops = [walks_op("ok1", "ok"), dict(walks_op("sleep", "sleep", "timeout"), timeout_s=1.5),
                   walks_op("ok2", "ok"), walks_op("set", "set", "exception"), walks_op("ok3", "ok"),
                   walks_op("die", "die", "child_died"), walks_op("ok4", "ok"),
                   walks_op("loop", "loop", "paging_runaway", paging="all"), walks_op("print", "print"),
                   walks_op("stdin", "stdin"), walks_op("nan", "nan", "non_finite"),
                   walks_op("inf", "inf", "non_finite"),
                   op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": 2048})]
            plan = {"scenarios": [{"id": "G-a", "graphlet": {"ops": ops}}, {"id": "G-b", "graphlet": {"ops": [
                walks_op("b-ok", "ok")]}}]}
            docs = {"G-a": DOCS["G-body"], "G-b": DOCS["G-body"]}
            out = os.path.join(tmp, "run")
            recs = run_offline(out, plan, docs=docs, ways="W1,W2", extra=["--library-root", FAKE_LIB,
                                                                          "--library-commit", "fake"])
            planned = dict((o["id"], o["expect"]) for s in plan["scenarios"] for o in s["graphlet"]["ops"])
            for way in ("W1", "W2"):
                with self.subTest(way):
                    got = dict((r["op"], r) for r in recs[way])
                    self.assertEqual(set(got), set(planned))
                    self.assertEqual(dict((k, outcome(r)) for k, r in got.items()), planned)
                    self.assertEqual(got["sleep"]["error"], {"code": "timeout", "fields": {"timeout_s": 1.5,
                                                                                           "page": 0}})
                    self.assertEqual(got["die"]["error"], {"code": "child_died", "fields": {"exit": 3, "page": 0}})
                    self.assertEqual(got["set"]["error"], {"code": "exception", "fields": {"type": "TypeError"}})
                    self.assertEqual(got["loop"]["error"], {"code": "paging_runaway", "fields": {
                        "pages": 2, "reason": "repeats"}})
                    self.assertEqual(got["loop"]["pages"], 2)
                    self.assertEqual(got["nan"]["error"], {"code": "non_finite", "fields": {"pointer": "/rows/1/cost",
                                                                                            "page": 0}})
                    self.assertEqual(got["inf"]["error"]["fields"]["pointer"], "/limits/loss")
                    # VL9: the records are strict JSON (no NaN / Infinity)
                    ops_path = os.path.join(out, "arms", rgl.WAY_DIRS[way], "graphlet_ops.jsonl")
                    with open(ops_path, encoding="utf-8") as fh:
                        for line in fh:
                            json.loads(line, parse_constant=lambda c: self.fail("non-standard JSON %s" % c))
                    # after each restart the document is reloaded: the scenario's handle works again
                    for name in ("ok2", "ok4", "sum"):
                        self.assertEqual(got[name]["answer"]["handle"], "<handle>")
                    self.assertEqual(got["sum"]["answer"]["summary"], {"loaded": 1})   # the new process: G-a reloaded
                    run = sc.read_json(os.path.join(out, "arms", rgl.WAY_DIRS[way], "run.json"))
                    self.assertEqual([r["reason"] for r in run["child"]["restarts"]], ["timeout", "child_died"])
                    self.assertEqual(run["state"], "done")
                    with open(os.path.join(out, "arms", rgl.WAY_DIRS[way], "child.log"), encoding="utf-8") as fh:
                        log = fh.read()
                    self.assertIn("noise on stdout", log)
                    self.assertIn("noise on fd 1", log)
                    raw = sc.read_gz_json(os.path.join(out, "arms", rgl.WAY_DIRS[way], "raw", "G-a.sleep.0.json.gz"))
                    self.assertEqual(raw["reply"]["supervised"], "timeout")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_a_load_receipt_holding_nan(self):
        """VR2: a graphlet_load receipt holding a NaN refuses that document (stage load, its pointer); the other
        scenario runs."""
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            ops = [op("w", "graphlet_walks", {"arm": "left", "mode": "ok", "max_bytes": 2048})]
            plan = {"scenarios": [{"id": sid, "graphlet": {"ops": ops}} for sid in ("G-nanload", "G-a")]}
            out = os.path.join(tmp, "run")
            recs = run_offline(out, plan, docs={"G-nanload": DOCS["G-body"], "G-a": DOCS["G-body"]}, ways="W4",
                               extra=["--library-root", FAKE_LIB, "--library-commit", "fake"])["W4"]
            got = dict((r["scenario"], r) for r in recs)
            self.assertEqual(got["G-nanload"]["error"], {"code": "document_refused", "fields": {
                "stage": "load", "non_finite": "/summary/seed/length_bp"}})
            self.assertEqual(outcome(got["G-a"]), "ok")
            docs = dict((d["scenario"], d) for d in records(out, "W4", "graphlet_docs.jsonl"))
            self.assertIsNone(docs["G-nanload"]["load"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_a_way_that_restarts_too_often(self):
        """VR5: past MAX_RESTARTS (here 2) the way ends: the op in flight is recorded (child_died), run.json says
        failed and why, the next way runs, the exit code is 1."""
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        saved = rgl.MAX_RESTARTS
        try:
            rgl.MAX_RESTARTS = 2
            ops = [op("die%d" % i, "graphlet_walks", {"arm": "left", "mode": "die", "max_bytes": 2048})
                   for i in range(4)]
            plan = {"scenarios": [{"id": "G-a", "graphlet": {"ops": ops}}]}
            out = os.path.join(tmp, "run")
            argv = ["--plan", write_json(os.path.join(tmp, "plan.json"), plan), "--doc", "G-a=" + DOCS["G-body"],
                    "--library-root", FAKE_LIB, "--library-commit", "fake", "--out", out, "--ways", "W1,W2",
                    "--quiet"]
            self.assertEqual(rgl.main(argv), 1)
            for way in ("W1", "W2"):
                recs = records(out, way)
                self.assertEqual([(r["op"], outcome(r)) for r in recs],
                                 [("die0", "child_died"), ("die1", "child_died"), ("die2", "child_died")])
                run = sc.read_json(os.path.join(out, "arms", rgl.WAY_DIRS[way], "run.json"))
                self.assertEqual(run["state"], "failed")
                self.assertIn("restarted 2 times", run["failure"])
        finally:
            rgl.MAX_RESTARTS = saved
            shutil.rmtree(tmp, ignore_errors=True)

    def test_a_seed_that_does_not_take_effect(self):
        """VL7: a child launched with -E echoes the seed it was given, but ignores it: the way does not start."""
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        saved = rgl.CHILD_FLAGS
        try:
            rgl.CHILD_FLAGS = saved + ("-E",)
            plan = {"scenarios": [{"id": "G-a", "graphlet": {"ops": [
                op("w", "graphlet_walks", {"arm": "left", "mode": "ok", "max_bytes": 2048})]}}]}
            argv = ["--plan", write_json(os.path.join(tmp, "plan.json"), plan), "--doc", "G-a=" + DOCS["G-body"],
                    "--library-root", FAKE_LIB, "--library-commit", "fake", "--out", os.path.join(tmp, "run"),
                    "--ways", "W1", "--quiet"]
            self.assertEqual(rgl.main(argv), 1)
            self.assertFalse(os.path.exists(os.path.join(tmp, "run", "arms", "dl-lib", "graphlet_ops.jsonl")))
            rgl.CHILD_FLAGS = saved
            again = argv[:-5] + ["--out", os.path.join(tmp, "run2"), "--ways", "W1", "--quiet"]
            self.assertEqual(rgl.main(again), 0)
        finally:
            rgl.CHILD_FLAGS = saved
            shutil.rmtree(tmp, ignore_errors=True)

    @needs_library
    def test_a_way_that_fails_leaves_the_others(self):
        """A way whose library does not start (the fake refuses max_bytes 77) is reported and the run goes on with
        the next way (here drift, on the real library); the invocation then exits 1."""
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            out = os.path.join(tmp, "run")
            plan = {"scenarios": [{"id": "G-body", "graphlet": {"ops": [
                op("w", "graphlet_walks", walks("left", max_bytes=77))]}}]}
            argv = ["--plan", write_json(os.path.join(tmp, "plan.json"), plan), "--doc", "G-body=" + DOCS["G-body"],
                    "--library-root", FAKE_LIB, "--drift-library-root", LIB, "--out", out, "--ways", "W4,drift",
                    "--limits", write_json(os.path.join(tmp, "limits.json"), {"max_bytes": 77}), "--quiet"]
            self.assertEqual(rgl.main(argv), 1)
            self.assertFalse(os.path.exists(os.path.join(out, "arms", "mcp-lib", "graphlet_ops.jsonl")))
            self.assertEqual([r["op"] for r in records(out, "drift")], ["w"])
            # VR5: the failed way's run.json says so, and why; the way that ran is done with no failure
            failed = sc.read_json(os.path.join(out, "arms", "mcp-lib", "run.json"))
            self.assertEqual(failed["state"], "failed")
            self.assertIn("refuses max_bytes 77", failed["failure"])
            self.assertFalse(os.path.exists(os.path.join(out, "arms", "mcp-lib", "work")))
            ran = sc.read_json(os.path.join(out, "arms", "drift-lib", "run.json"))
            self.assertEqual((ran["state"], ran["failure"]), ("done", None))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@needs_library
class PagingTests(unittest.TestCase):
    """Paging assembly over small max_bytes equals the uncut one-page answer of a large max_bytes, for every list
    tool; the sequence's next_from loop equals one slice under a large ceiling."""

    LISTS = {
        "lims": ("G-switch", "graphlet_summary", {"arm": None, "detail": "limitations"}),
        "walks": ("G-switch", "graphlet_walks", dict(walks("left", rank="support", spell="tail"))),
        "walk": ("G-switch", "graphlet_walk", {"arm": "left", "walk": 9}),
        "support": ("G-switch", "graphlet_support", {"arm": "left", "walk": 9, "step": None}),
        "labels": ("G-switch", "graphlet_labels", labels()),
        "labels-name": ("G-switch", "graphlet_labels", labels(name={"ref": "h:5:6"})),
        "splits": ("G-switch", "graphlet_splits", {"arm": "left", "n": 100, "min_labels_before": 0}),
        "claims": ("G-switch", "graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": False,
                                                   "labels": None}),
        "compare": ("G-trace", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None, "mode": "claims"}),
        "compare-labels": ("G-trace", "graphlet_compare", {"b": {"graphlet": "G-annot"}, "arm": None,
                                                           "mode": "labels"}),
        "list": ("G-trace", "graphlet_list", {}),   # after G-trace's compares loaded G-annot: three entries
    }
    SIZES = (500, 700, 800, 900, 1000, 1500, 2000, 3000, 5000)

    def test_assembled_equals_uncut(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            ops = {"G-switch": [], "G-trace": []}
            for name, (sid, tool, args) in sorted(self.LISTS.items()):
                args = dict(args)
                args.pop("max_bytes", None)
                ops[sid].append(op(name + "@big", tool, dict(args, max_bytes=262144)))
                for size in self.SIZES:
                    ops[sid].append(op("%s@%d" % (name, size), tool, dict(args, max_bytes=size), paging="all"))
            plan = {"scenarios": [{"id": sid, "graphlet": {"ops": o}} for sid, o in sorted(ops.items())]}
            recs = run_offline(os.path.join(tmp, "run"), plan, ways="W4")["W4"]
            by = dict((r["op"], r) for r in recs)
            cut_rows = 0
            for name in sorted(self.LISTS):
                with self.subTest(name):
                    big = comparable(by[name + "@big"])["answer"]
                    self.assertEqual(by[name + "@big"]["pages"], 1)
                    self.assertGreater(big["total"], 1)
                    used = []
                    for size in self.SIZES:
                        r = by["%s@%d" % (name, size)]
                        if not r["error"] and r["answer"].get("row_truncated"):
                            cut_rows += self.check_cut_rows(comparable(r)["answer"], big)
                        if r["error"] or r["pages"] < 2 or r["answer"].get("row_truncated"):
                            continue
                        used.append(size)
                        self.assertTrue(r["paging_complete"])
                        small = comparable(r)["answer"]
                        # compare's base gives way per page (_page_cut): a field cut on EVERY page stays cut and
                        # named; the assembled answer is then the uncut one without exactly those optional fields
                        cut = small.pop("fields_cut", [])
                        self.assertLessEqual(set(cut), {"evidence", "export"})
                        self.assertEqual(small, dict((k, v) for k, v in big.items() if k not in cut),
                                         "max_bytes %d" % size)
                    self.assertTrue(used, "no size gave whole rows over several pages")
            cuts = [by["compare@%d" % size]["answer"].get("fields_cut") for size in self.SIZES
                    if by["compare@%d" % size]["answer"]]
            self.assertIn(["evidence"], cuts)
            self.assertGreater(cut_rows, 5)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def check_cut_rows(self, small, big):
        """FORMAT v2 A3: an answer with cut rows has every row of the uncut answer, in order; a row marked _cut_fields
        is the uncut row with each named field cut (a string's head, a list's head: mcp_tools._shrink), the rest
        equal; an unmarked row is the uncut row. -> the number of marked rows."""
        self.assertEqual(len(small["rows"]), len(big["rows"]))
        self.assertTrue(small["cut_fields"])
        marked = 0
        for row, whole in zip(small["rows"], big["rows"]):
            row = dict(row)
            fields = row.pop("_cut_fields", None)
            if fields is None:
                self.assertEqual(row, whole)
                continue
            marked += 1
            self.assertLessEqual(set(fields), set(small["cut_fields"]))
            for k, v in row.items():
                if k in fields:
                    self.assertEqual(v, whole[k][:len(v)], k)
                    self.assertLess(len(v), len(whole[k]), k)
                else:
                    self.assertEqual(v, whole[k], k)
        self.assertGreater(marked, 0)
        return marked

    def test_sequence_pages_equal_one_slice(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            ops = [op("s-%d" % i, "graphlet_sequence", seq(**kw), paging="all") for i, kw in enumerate((
                {"walk": 9}, {"walk": 9, "with_seed": True}, {"walk": 3, "orientation": "walk", "with_seed": True},
                {"walk": 9, "start": 100, "end": 1500}, {"segment": 3}))]
            plan = {"scenarios": [{"id": "G-switch", "graphlet": {"ops": ops}}]}
            small = run_offline(os.path.join(tmp, "small"), plan, limits={"sequence_max_bytes": 700})["W4"]
            large = run_offline(os.path.join(tmp, "large"), plan, limits={"sequence_max_bytes": 1 << 20})["W4"]
            self.assertTrue(all(r["pages"] == 1 for r in large))
            self.assertGreater(sum(r["pages"] > 1 for r in small), 3)
            for a, b in zip(small, large):
                with self.subTest(a["op"]):
                    a, b = comparable(a)["answer"], comparable(b)["answer"]
                    self.assertEqual((a.pop("ceiling_bytes"), b.pop("ceiling_bytes")), (700, 1 << 20))
                    self.assertEqual(a, b)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@needs_library
class BudgetTests(unittest.TestCase):
    def test_budgeted_pages_resume_to_the_whole_list(self):
        """Under local limits a small work budget stops a list after a few rows (complete false, a resume cursor);
        paging "all" resumes it page by page, each under its own budget, and the assembly takes total, complete, the
        filter counts, their hint and the local block from the last page: the WHOLE answer (but for the local block,
        whose usage is per call) is the unbudgeted call's. claims-rc is VERIFY-L VL1's reproduction: its first budget
        segment filters nothing, so `filtered` and `hint` appear on later pages only."""
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            lists = {"claims": ("graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": False,
                                                    "labels": None}, 800),
                     "claims-rc": ("graphlet_claims", {"arm": "left", "min_bp": 0, "route_consistent": True,
                                                       "labels": None}, 800),
                     "walks": ("graphlet_walks", walks("left", spell="full"), 1000),
                     "splits": ("graphlet_splits", {"arm": "left", "n": 100, "min_labels_before": 0}, 600),
                     "walk": ("graphlet_walk", {"arm": "left", "walk": 9}, 700)}
            ops = []
            for name, (tool, args, wu) in sorted(lists.items()):
                args = dict(args, max_bytes=262144)
                ops += [op(name, tool, args), op(name + "-b", tool, dict(args, budget={"work_units": wu}),
                                                 paging="all"),
                        op(name + "-b1", tool, dict(args, budget={"work_units": wu}), paging={"pages": 1})]
            plan = {"scenarios": [{"id": "G-switch", "graphlet": {"ops": ops}}]}
            by = dict((r["op"], r) for r in run_offline(os.path.join(tmp, "run"), plan, docs={
                "G-switch": DOCS["G-switch"]}, limits={"local_limits": {}})["W4"])
            for name in lists:
                with self.subTest(name):
                    whole, paged, first = by[name]["answer"], by[name + "-b"]["answer"], by[name + "-b1"]["answer"]
                    self.assertGreater(by[name + "-b"]["pages"], 1)
                    self.assertFalse(first["complete"])
                    self.assertFalse(by[name + "-b1"]["paging_complete"])
                    self.assertTrue(paged["complete"])
                    self.assertTrue(paged["local"]["complete"])
                    self.assertEqual(dict((k, v) for k, v in paged.items() if k != "local"),
                                     dict((k, v) for k, v in whole.items() if k != "local"))
                    self.assertEqual(first["rows"], whole["rows"][:len(first["rows"])])
                    self.assertEqual(by[name]["library"]["limits"]["local_limits"], {})
            # the reproduction: page 0 states neither filtered_so_far nor hint, the assembled answer both
            raw0 = sc.read_gz_json(os.path.join(tmp, "run", "arms", "mcp-lib", "raw", "G-switch.claims-rc-b.0.json.gz"))
            self.assertFalse({"filtered_so_far", "hint"} & set(raw0["reply"]["result"]["answer"]))
            self.assertEqual(set(by["claims-rc-b"]["answer"]) & {"filtered", "hint"}, {"filtered", "hint"})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@needs_library
class ReparsedSequenceTests(unittest.TestCase):
    def test_a_reparse_between_sequence_pages(self):
        """VERIFY-L VR1's reproduction, in process: local limits, a store that drops the body from RAM between two
        pages (ttl_ram_s 0.5, a clock one second on per page), a 1,200-byte sequence ceiling: page 1's `local` gains
        `parsed` (the body re-parsed). The pages assemble into the one-slice answer, `local` from the last page."""
        if LIB not in sys.path:
            sys.path.insert(0, LIB)
        import graphlet_paging as gp
        from metagraph.traverse import GraphletStore, parser
        from metagraph.traverse.mcp_tools import GraphletTools, ToolLimits
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            now = [1000.0]
            store = GraphletStore(os.path.join(tmp, "spool"), ttl_ram_s=0.5, clock=lambda: now[0])
            tools = GraphletTools(store, secret=b"k" * 32, export_dir=os.path.join(tmp, "x"),
                                  sequence_max_bytes=1200, local_limits=ToolLimits())
            resp = read_json(DOCS["G-switch"])
            r = resp["results"][0]
            path = os.path.join(tmp, "spool", "doc.mgt")
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write(parser.standalone_text(r["graphlet"], r, resp))
            handle = tools.graphlet_load(path=path)["handle"]
            args = dict(handle=handle, arm="left", walk=9)
            pages = [tools.graphlet_sequence(**args)]
            while gp.more("graphlet_sequence", "all", pages):
                now[0] += 1.0
                arg, value = gp.next_call("graphlet_sequence", pages[-1])
                pages.append(tools.graphlet_sequence(**dict(args, **{arg: value})))
            self.assertGreater(len(pages), 10)
            parsed = ["parsed" in (p.get("local") or {}) for p in pages]
            self.assertTrue(any(parsed[1:]), parsed)
            self.assertNotEqual(pages[0]["local"], pages[1]["local"])
            answer = gp.assemble("graphlet_sequence", pages, True)
            self.assertEqual(answer["local"], pages[-1]["local"])
            whole = GraphletTools(GraphletStore(os.path.join(tmp, "spool2")), secret=b"k" * 32,
                                  export_dir=os.path.join(tmp, "x2"), sequence_max_bytes=1 << 20)
            path2 = os.path.join(tmp, "spool2", "doc.mgt")   # graphlet_load reads its own spool or export dir
            shutil.copy(path, path2)
            one = whole.graphlet_sequence(**dict(args, handle=whole.graphlet_load(path=path2)["handle"]))
            self.assertEqual(answer["sequence"], one["sequence"])
            self.assertEqual((answer["from"], answer["to"], answer["length"]), (one["from"], one["to"], one["length"]))
            self.assertNotIn("truncated", answer)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@needs_library
class ConversionTests(unittest.TestCase):
    def test_response_equals_file(self):
        """A W1-style document (the response, converted in the way's process by standalone_text) and the same
        document saved by the library (parser.save of from_response, gzipped: W4's form) give the same records."""
        parser = library()
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            files = {}
            for sid in ("G-switch", "G-trace", "G-annot", "G-coords", "G-merge"):
                resp = read_json(DOCS[sid])
                g = parser.from_response(resp["results"][0], resp)
                plain = os.path.join(tmp, sid + ".mgt")
                parser.save(g, plain)
                with open(plain, "rb") as fh, gzip.open(plain + ".gz", "wb") as out:
                    out.write(fh.read())
                files[sid] = plain + ".gz"
            files["G-body"] = DOCS["G-body"]
            plan = every_tool_plan()
            a = run_offline(os.path.join(tmp, "resp"), plan, ways="W1")["W1"]
            b = run_offline(os.path.join(tmp, "file"), plan, docs=files, ways="W1")["W1"]
            self.assertEqual([comparable(r) for r in a], [comparable(r) for r in b])
            for sid in files:
                with gzip.open(os.path.join(tmp, "resp", "arms", "dl-lib", "docs", sid + ".mgt.gz")) as x, \
                        gzip.open(os.path.join(tmp, "file", "arms", "dl-lib", "docs", sid + ".mgt.gz")) as y:
                    self.assertEqual(x.read(), y.read(), sid)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def fake_run_dir(tmp, ids):
    """A run directory as the arms leave it: each scenario's /traverse response as the direct arms' raw records, and the
    MCP's export as its .mgt.gz (the library's standalone text of the same response)."""
    parser = library()
    run = os.path.join(tmp, "run")
    for sid, path in ids.items():
        resp = read_json(path)
        for arm in ("direct-local", "direct-remote"):
            rel = os.path.join(run, "arms", arm, sc.raw_name(sid, arm, 0))
            os.makedirs(os.path.dirname(rel), exist_ok=True)
            sc.write_gz_json(rel, {"scenario": sid, "arm": arm, "repeat": 0,
                                   "answer": {"http": 200, "body": resp, "error": None}})
        mgt = os.path.join(run, "arms", "mcp", "docs", sid + ".mgt.gz")   # FORMAT v2 C1
        os.makedirs(os.path.dirname(mgt), exist_ok=True)
        with gzip.open(mgt, "wb") as fh:
            r = resp["results"][0]
            fh.write(parser.standalone_text(r["graphlet"], r, resp).encode("utf-8"))
    return run


@needs_library
class ManifestModeTests(unittest.TestCase):
    def test_three_ways_agree(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            run = fake_run_dir(tmp, {"G-switch": DOCS["G-switch"], "G-trace": DOCS["G-trace"],
                                     "G-annot": DOCS["G-annot"]})
            plan = every_tool_plan()
            scenarios = [s for s in plan["scenarios"] if s["id"] in ("G-switch", "G-trace", "G-annot")]
            for s in scenarios:
                s["graphlet"]["ways"] = {"default": ["W1", "W2", "W3", "W4", "W5"]}
                if s["id"] == "G-trace":   # an op of the service only, and one of W1 alone
                    s["graphlet"]["ops"][0]["ways"] = ["W3"]
                    s["graphlet"]["ops"][1]["ways"] = ["W1"]
            manifest = write_json(os.path.join(tmp, "manifest.json"), {"manifest_version": "8-test",
                                                                       "scenarios": [{"id": "T-other"}] + scenarios})
            rc = rgl.main(["--run-dir", run, "--manifest", manifest, "--library-root", LIB, "--quiet"])
            self.assertEqual(rc, 0)
            recs = dict((w, records(run, w)) for w in ("W1", "W2", "W4"))
            trace = [r["op"] for r in recs["W4"] if r["scenario"] == "G-trace"]
            self.assertNotIn("cmp-claims", trace)
            self.assertNotIn("cmp-walks", trace)
            self.assertIn("cmp-walks", [r["op"] for r in recs["W1"] if r["scenario"] == "G-trace"])
            w1 = [comparable(r, True) for r in recs["W1"] if r["op"] != "cmp-walks"]
            self.assertEqual(w1, [comparable(r, True) for r in recs["W2"] if r["op"] != "cmp-walks"])
            self.assertEqual([r for r in w1 if r["scenario"] == "G-switch"],
                             [comparable(r, True) for r in recs["W4"] if r["scenario"] == "G-switch"])
            for way, source in (("W1", "arms/direct-local/raw/G-switch.direct-local.0.json.gz"),
                                ("W4", "arms/mcp/docs/G-switch.mgt.gz")):
                docs = dict((d["scenario"], d) for d in records(run, way, "graphlet_docs.jsonl"))
                self.assertEqual(docs["G-switch"]["source"], source)
                self.assertEqual(set(docs), {"G-switch", "G-trace", "G-annot"})   # G-annot: a compare's b
            self.assertEqual(sc.read_json(os.path.join(run, "arms", "mcp-lib", "run.json"))["run_id"], "run")
            # one scenario selected: its compare partner is loaded all the same
            rc = rgl.main(["--run-dir", run, "--manifest", manifest, "--library-root", LIB, "--quiet", "--out",
                           os.path.join(tmp, "one"), "--ways", "W4", "--scenarios", "G-trace"])
            self.assertEqual(rc, 0)
            one = records(os.path.join(tmp, "one"), "W4")
            self.assertEqual(set(r["scenario"] for r in one), {"G-trace"})
            self.assertEqual([comparable(r) for r in one],
                             [comparable(r) for r in recs["W4"] if r["scenario"] == "G-trace"])
            with self.assertRaises(SystemExit):
                rgl.main(["--run-dir", run, "--manifest", manifest, "--library-root", LIB, "--quiet", "--out",
                          os.path.join(tmp, "x"), "--scenarios", "G-nope"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_malformed_plans_refused(self):
        for bad in ({"ops": [{"id": "a", "tool": "graphlet_summary", "paging": "some"}]},
                    {"ops": [{"id": "a", "tool": "graphlet_summary"}, {"id": "a", "tool": "graphlet_walks"}]},
                    {"ops": [{"id": "a", "tool": 3}]}, {"ops": [{"id": "a", "tool": "x", "args": []}]},
                    {"ops": [{"id": "a", "tool": "x", "paging": {"pages": 0}}]}):
            with self.subTest(bad):
                problems = rgl.plan_problems(rgl.plan_scenarios(bad, ["G"]))
                self.assertTrue(problems)
        self.assertEqual(rgl.plan_problems(every_tool_plan()["scenarios"]), [])


@needs_library
class ExportAndDriftTests(unittest.TestCase):
    def test_export_library_and_drift_report(self):
        tmp = tempfile.mkdtemp(prefix="graphlet-local-")
        try:
            repo = subprocess.check_output(["git", "-C", SUITE, "rev-parse", "--show-toplevel"]).decode().strip()
            head = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"]).decode().strip()
            lib = os.path.join(tmp, "lib")
            self.assertEqual(rgl.export_library("HEAD", lib, repo, quiet=True), head)
            self.assertTrue(os.path.isfile(os.path.join(lib, "metagraph", "api", "python", "metagraph", "traverse",
                                                        "mcp_tools.py")))
            self.assertEqual(rgl.library_commit(lib), head)
            # VL15: a root given as the export's metagraph/api/python finds the export's LIBRARY_COMMIT three levels
            # up; a root that is no such directory looks there for none
            self.assertEqual(rgl.library_commit(os.path.join(lib, "metagraph", "api", "python")), head)
            self.assertEqual(rgl.library_commit(FAKE_LIB, "fake"), "fake")
            # the same commit again: left as it is
            self.assertEqual(rgl.export_library(head, lib, repo, quiet=True), head)
            with self.assertRaises(SystemExit):
                rgl.export_library("HEAD~1", lib, repo, quiet=True)
            plan = {"scenarios": [s for s in every_tool_plan()["scenarios"] if s["id"] in ("G-trace", "G-annot")]}
            out = os.path.join(tmp, "run")
            docs = {"G-trace": DOCS["G-trace"], "G-annot": DOCS["G-annot"]}
            run_offline(out, plan, docs=docs, ways="W4,drift", extra=["--drift-library-root", LIB])
            w4 = records(out, "W4")
            self.assertEqual(set(r["library"]["commit"] for r in w4), {rgl.library_commit(LIB)})
            report = rgl.drift_report(out, quiet=True)
            self.assertEqual((report["compared"], report["differ"]), (len(w4), []))
            # the exported copy runs too, and a drift record that differs is listed
            run_offline(os.path.join(tmp, "copy"), plan, docs=docs, ways="W4", extra=["--library-root", lib,
                                                                                    "--keep-work"])
            self.assertTrue(os.path.isdir(os.path.join(tmp, "copy", "arms", "mcp-lib", "work", "spool")))
            self.assertEqual(set(r["library"]["commit"] for r in records(os.path.join(tmp, "copy"), "W4")), {head})
            path = os.path.join(out, "arms", "drift-lib", "graphlet_ops.jsonl")
            rows = records(out, "drift")
            rows[0]["answer"]["equal"] = "changed"
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("".join(json.dumps(r) + "\n" for r in rows))
            report = rgl.drift_report(out, quiet=True)
            self.assertEqual(report["differ"], [{"scenario": "G-trace", "op": rows[0]["op"], "tool": "graphlet_compare",
                                                 "fields": ["answer"]}])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)



class LibraryLimitsFromTheCapabilities(unittest.TestCase):
    """The ways' GraphletTools from the capabilities' service_settings.graphlet.library block where it states a key
    (DESIGN-v8 §3), else --limits, else the defaults; a part absent or served null is not stated (the fallback, the
    missing limit recorded); which was used is recorded per key (run.json limits_source)."""

    BLOCK = {"max_bytes": 4096, "sequence_max_bytes": 32768, "local_limits": None, "max_graphlet_mb": None}

    def refs(self, library):
        doc = None if library is ... else {"service_settings": {"graphlet": {"library": library}}}
        return sc.References(None, {"service": doc}, {"service": "arms/mcp/run.json service_capabilities"})

    def limits_file(self, obj):
        d = tempfile.mkdtemp(prefix="rgl-limits-")
        self.addCleanup(shutil.rmtree, d, True)
        path = os.path.join(d, "limits.json")
        with open(path, "w") as fh:
            json.dump(obj, fh)
        return path

    def test_from_the_block(self):
        limits, ignored, rec = rgl.library_limits(refs=self.refs(dict(self.BLOCK)))
        self.assertEqual(limits, {"max_bytes": 4096, "sequence_max_bytes": 32768, "local_limits": None})
        self.assertEqual(set(rec["source"].values()), {"capabilities (live)"})
        self.assertEqual(rec["missing"], [])
        self.assertEqual(rec["document"], "arms/mcp/run.json service_capabilities")
        self.assertEqual(ignored, {})
        _l, ignored, _r = rgl.library_limits(refs=self.refs(dict(self.BLOCK, max_graphlet_mb=8)))
        self.assertEqual(ignored["max_graphlet_mb"]["value"], 8)
        self.assertEqual(ignored["max_graphlet_mb"]["from"], "capabilities (live)")

    def test_null_or_absent_parts_fall_back(self):
        path = self.limits_file({"max_bytes": 1024})
        limits, _i, rec = rgl.library_limits(path, refs=self.refs(dict(self.BLOCK, max_bytes=None)))
        self.assertEqual((limits["max_bytes"], limits["sequence_max_bytes"]), (1024, 32768))
        self.assertEqual(rec["source"]["max_bytes"], "--limits")
        self.assertEqual(rec["missing"], ["service_settings.graphlet.library.max_bytes (null; --limits used)"])
        block = dict(self.BLOCK)
        del block["sequence_max_bytes"]
        limits, _i, rec = rgl.library_limits(refs=self.refs(dict(block, max_bytes=10)))
        self.assertEqual(limits, rgl.DEFAULT_LIMITS)
        self.assertEqual(rec["source"], {"max_bytes": "default", "sequence_max_bytes": "default",
                                         "local_limits": "capabilities (live)"})
        self.assertEqual(rec["missing"], ["service_settings.graphlet.library.max_bytes (not a valid value 10; default "
                                          "used)", "service_settings.graphlet.library.sequence_max_bytes (absent)"])
        # the whole block null (a part the service failed to compute) or absent: --limits / the defaults
        for library in (None, ...):
            limits, _i, rec = rgl.library_limits(path, refs=self.refs(library))
            self.assertEqual(limits["max_bytes"], 1024)
            self.assertEqual(rec["missing"], ["service_settings.graphlet.library (absent: --limits or the defaults "
                                              "used)"])

    def test_a_run_directory_records_its_document(self):
        d = tempfile.mkdtemp(prefix="rgl-run-")
        self.addCleanup(shutil.rmtree, d, True)
        os.makedirs(os.path.join(d, "arms", "mcp"))
        with open(os.path.join(d, "arms", "mcp", "run.json"), "w") as fh:
            json.dump({"service_capabilities": {"document": {"service_settings": {"graphlet": {
                "library": dict(self.BLOCK)}}}}}, fh)
        limits, _i, rec = rgl.library_limits(None, d)
        self.assertEqual(limits["max_bytes"], 4096)
        self.assertIn("service_capabilities", rec["document"])
        limits, _i, rec = rgl.library_limits(None, os.path.join(d, "nothing"))
        self.assertEqual(limits, rgl.DEFAULT_LIMITS)

if __name__ == "__main__":
    unittest.main()
