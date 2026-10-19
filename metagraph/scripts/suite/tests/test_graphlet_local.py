"""Offline tests of area G's local side: graphlet_doc (the suite's document digests) against the library's own
parser.body_digest on the library's fixtures. The tests import the library (metagraph/api/python of this checkout);
the modules under test never do."""
import glob
import gzip
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))
DATA = os.path.join(LIB, "tests", "data", "traverse")
sys.path.insert(0, SUITE)

import graphlet_doc  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
