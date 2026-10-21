"""Offline tests of graphlet_oracle.py (W5, manifest v8 phase A) on recorded ENGINE answers only: the `detail: full`
and `detail: graphlet` responses of one request (api/python/tests/data/traverse: real/* snapshots of the live
mini_refseq, SRA and UHGG servers, coords/* of mini_refseq, documents/* written by the CLI; the two hand-made
documents, clipped_merge and limits, are left out). No library import: the library is tested against W5 in
test_graphlet_oracle_vs_library.py, never the other way round. Hand goldens are read off the fixtures' own JSON
(the walk, its runs, its continuation) and written here as literals."""
import ast
import copy
import gzip
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import checks  # noqa: E402
import graphlet_oracle as O  # noqa: E402

DATA = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python", "tests", "data", "traverse"))
HAND_MADE = {"clipped_merge", "limits"}
SUBDIRS = ("real/mini_refseq", "real/sra", "real/uhgg", "documents", "documents/budgets", "coords")
NOT_TWINS = {"real/sra/sra_16s_OK510341__time_tight", "coords/d3_memory_stop"}
# real pairs of the area G rehearsal 2 (refseq97, 2026-10-19, run 20261019T220627Z): the engine's `detail: graphlet`
# response body (the direct-local raw answer) and its full twin as v8 sends it, without bounds.max_memory_mb (VD4)
AREA_G = os.path.join(HERE, "data", "area_g")
REHEARSAL2 = ("G-cut-c1", "G-merge-c1")


def _load(path):
    with (gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")) as f:
        return json.load(f)


def engine_pairs():
    """(cell, seed index, graphlet response, full response) of every engine-made pair, every seed with a body."""
    out = []
    for sub in SUBDIRS:
        d = os.path.join(DATA, sub)
        for n in sorted(os.listdir(d)):
            for suf in (".graphlet.json.gz", ".graphlet.json"):
                cell = n[:-len(suf)]
                if not n.endswith(suf) or cell in HAND_MADE:
                    continue
                fp = os.path.join(d, cell + suf.replace("graphlet", "full"))
                if not os.path.exists(fp):
                    continue
                g, f = _load(os.path.join(d, n)), _load(fp)
                for i, rg in enumerate(g["results"]):
                    if "graphlet" in rg:
                        out.append((sub + "/" + cell, i, g, f))
    return out


def doc(name):
    """A documents/ pair: (graphlet response, full response, MGT body)."""
    g = _load(os.path.join(DATA, "documents", name + ".graphlet.json"))
    f = _load(os.path.join(DATA, "documents", name + ".full.json"))
    return g, f, g["results"][0]["graphlet"]


def rehearsal2(sid):
    """A rehearsal-2 pair (AREA_G): (graphlet response, full twin, MGT body)."""
    g = _load(os.path.join(AREA_G, "rehearsal2-20261019-%s.graphlet.json.gz" % sid))
    f = _load(os.path.join(AREA_G, "rehearsal2-20261019-%s.full.json.gz" % sid))
    return g, f, g["results"][0]["graphlet"]


def rec(tool, args, name, **ctx):
    g, f, mgt = doc(name)
    return O.answer(tool, args, f, mgt, ctx)


class EnginePairs(unittest.TestCase):
    """W5 on every engine pair: its derivations agree with what the engine itself states elsewhere (the graphlet's
    JSON summary, the A records, the G records), and the twin pre-check sorts the pairs."""

    @classmethod
    def setUpClass(cls):
        cls.pairs = engine_pairs()

    def test_corpus(self):
        self.assertGreaterEqual(len(self.pairs), 42)   # 40 twins, NOT_TWINS
        modes = {f["results"][i]["label_mode"] for _, i, _, f in self.pairs}
        self.assertEqual(modes, {"constrain", "annotate"})

    def test_twin_pre_check(self):
        """DESIGN-v8 §2: graphlet_arms of the body against the full twin; the time-stopped SRA cell and the coords
        cell stopped by a memory budget in full detail (SPEC §7.1: arrival order there) are no twins."""
        seen = set()
        for cell, i, g, f in self.pairs:
            why = O.twin_difference(f, g["results"][i]["graphlet"], i)
            if cell in NOT_TWINS:
                self.assertIsNotNone(why, cell)
                self.assertRegex(why, r"^(left|right) (segment|path) \d+: full ", cell)
                rec_ = O.answer("graphlet_summary", {}, f, g["results"][i]["graphlet"], {"seed_index": i})
                self.assertEqual(rec_["oracle"]["declared"], "twin_mismatch")
                self.assertEqual(O.differences(rec_, {"anything": 1}), [])
                seen.add(cell)
            else:
                self.assertIsNone(why, cell)
        self.assertEqual(seen, NOT_TWINS)

    def test_summary_is_the_engine_summary(self):
        """graphlet_summary from the twin equals the graphlet's own JSON summary (SPEC §7.5.5 counts, status,
        complete_to_bp, completeness_scope, outcome: computed by the engine for the other execution) and the A
        records' validation counts (SPEC §7.5.2)."""
        n = 0
        for cell, i, g, f in self.pairs:
            if cell in NOT_TWINS:
                continue
            rg = g["results"][i]
            r = O.answer("graphlet_summary", {}, f, rg["graphlet"], {"seed_index": i})
            s = r["answer"]["summary"]
            self.assertEqual(s["outcome"], rg["outcome"], cell)
            self.assertEqual(s["seed"]["length_bp"], rg["seed"]["length_bp"], cell)
            self.assertEqual(s["seed"]["mode"], rg["label_mode"], cell)
            arecs = {x.split(" ")[1]: x.split(" ") for x in rg["graphlet"].split("\n") if x.startswith("A ")}
            for side, a in rg["arms"].items():
                got, c = s["arms"][side], a["counts"]
                self.assertEqual((got["status"], got["complete_to_bp"], got["scope"]),
                                 (a["status"], a["complete_to_bp"], a["completeness_scope"]), cell)
                self.assertEqual((got["walks"], got["segments"], got["splits"], got["merges"], got["bases"],
                                  got["max_bp"]),
                                 (c["leaves"], c["segments"], c["splits"], c["merges"], c["bases"], c["max_bp"]),
                                 "%s %s" % (cell, side))
                A = arecs[side[0]]
                self.assertEqual([got["segments"], got["walks"], got["splits"], got["merges"], got["bases"]],
                                 [int(A[-6]), int(A[-4]), int(A[-3]), int(A[-2]), int(A[-1])], cell)
                n += 1
            self.assertEqual(r["answer"]["summary"]["index"]["verifiable"],
                             rg["graphlet"].split("\n", 1)[0].split(" ")[14] != "*", cell)
        self.assertGreater(n, 60)

    def test_segment_bases_in_both_orientations(self):
        """graphlet_sequence(segment=s): natural orientation from the twin's JSON sequences, walk orientation equal
        to the G record's bases, stored in walking order on both arms (SPEC §7.5.4): two engine sources."""
        n = 0
        for cell, i, g, f in self.pairs:
            if cell in NOT_TWINS:
                continue
            mgt = g["results"][i]["graphlet"]
            g_bases, side = {}, None
            for line in mgt.split("\n"):
                p = line.split(" ")
                if p[0] == "A":
                    side, k = ("left" if p[1] == "l" else "right"), 0
                elif p[0] == "G":
                    g_bases[(side, k)] = "" if p[10] == "." else p[10]
                    k += 1
            for (side, sid), bases in g_bases.items():
                if bases == "*":
                    continue
                r = O.answer("graphlet_sequence", {"arm": side, "segment": sid, "orientation": "walk"}, f, mgt,
                             {"seed_index": i})
                self.assertEqual(r["answer"]["sequence"], bases, "%s %s %d" % (cell, side, sid))
                nat = O.answer("graphlet_sequence", {"arm": side, "segment": sid}, f, mgt, {"seed_index": i})
                self.assertEqual(nat["answer"]["sequence"], bases[::-1] if side == "left" else bases)
                n += 1
        self.assertGreater(n, 1000)

    def test_molecules_and_continuations(self):
        """A walk with the seed (natural(left) + seed | seed + natural(right), SPEC §7.5.4) holds the engine's
        continuation where the C record derives it (SPEC §7.5.3: right the last n bases of seed + natural(right),
        left the first n of natural(left) + seed)."""
        n = 0
        for cell, i, g, f in self.pairs:
            if cell in NOT_TWINS:
                continue
            mgt = g["results"][i]["graphlet"]
            for side, a in f["results"][i]["arms"].items():
                for p in a["paths"]:
                    c = p.get("continuation")
                    if not c or not c.get("sequence") or "sequence" not in a["segments"][0]:
                        continue
                    m = O.answer("graphlet_sequence", {"arm": side, "walk": p["id"], "with_seed": True}, f, mgt,
                                 {"seed_index": i})["answer"]["sequence"]
                    want = m[:len(c["sequence"])] if side == "left" else m[len(m) - len(c["sequence"]):]
                    self.assertEqual(want, c["sequence"], "%s %s %d" % (cell, side, p["id"]))
                    n += 1
        self.assertGreater(n, 300)

    def test_every_record_names_its_clauses(self):
        """Every value W5 answers is claimed by a clause (the record builder refuses otherwise), every checked
        pointer resolves, and the record compares equal with its own answer (differences() is empty)."""
        for cell, i, g, f in self.pairs[::3]:
            if cell in NOT_TWINS:
                continue
            mgt = g["results"][i]["graphlet"]
            sides = [s for s in ("left", "right") if s in f["results"][i]["arms"]]
            for tool, args in [("graphlet_summary", {}), ("graphlet_summary", {"detail": "limitations"}),
                               ("graphlet_labels", {}), ("graphlet_export", {"format": "json"}),
                               ("graphlet_export", {"format": "mgt"}), ("graphlet_save", {})] + \
                    [(t, dict(a, arm=s)) for s in sides for t, a in
                     [("graphlet_walks", {"rank": "support"}), ("graphlet_walks", {"rank": "length", "spell": "full"}),
                      ("graphlet_walk", {"walk": 0}), ("graphlet_splits", {"min_labels_before": 0})]]:
                r = O.answer(tool, args, f, mgt, {"seed_index": i})
                self.assertIn("checked", r["oracle"])
                for ptr, clause in r["oracle"]["checked"].items():
                    self.assertTrue(clause.startswith(("SPEC", "DESIGN", "the ", "one ", "<side>", "more than",
                                                       "path id")), (tool, ptr, clause))
                    O.resolve_pointer(r["answer"], ptr)
                self.assertEqual(O.differences(r, r["answer"]), [], (cell, tool, args))


class HandGoldens(unittest.TestCase):
    """Answers read off small CLI documents by hand: `caps` (a switch hY -> hZ at 103 under a table cost, a
    max_steps leaf, no sequences), `merge` (two non-first-parent merges at 36 and 62), `fork` (both arms)."""

    def test_caps_walks(self):
        sup = rec("graphlet_walks", {"arm": "right", "rank": "support", "n": 10}, "caps")["answer"]
        # walk 0's only leaf label hZ was switched in at 103 (run 3): not full; walks 1 and 2 end hX, full
        self.assertEqual([(r["walk"], r["length_bp"], r["end"], r["n_alive"], r["n_full"]) for r in sup["rows"]],
                         [(1, 109, "edge_reuse", 1, 1), (2, 95, "loss_budget", 1, 1), (0, 126, "max_steps", 1, 0)])
        self.assertEqual(sup["rows"][0]["labels_full"], [{"name": "hX", "ref": "h:0:0"}])
        self.assertEqual(sup["total"], 3)
        ln = rec("graphlet_walks", {"arm": "right", "rank": "length"}, "caps")["answer"]
        self.assertEqual([r["walk"] for r in ln["rows"]], [0, 1, 2])
        mb = rec("graphlet_walks", {"arm": "right", "rank": "id", "min_bp": 100}, "caps")["answer"]
        self.assertEqual([r["walk"] for r in mb["rows"]], [0, 1])
        # hZ switched in on walk 0: route-consistent (no merge), kept; F2
        hz = rec("graphlet_walks", {"arm": "right", "rank": "id", "label": "hZ"}, "caps")["answer"]
        self.assertEqual(([r["walk"] for r in hz["rows"]], hz.get("filtered")), ([0], None))
        loss = rec("graphlet_walks", {"arm": "right", "rank": "loss"}, "caps")
        self.assertEqual((loss["answer"]["total"], loss["answer"]["rows"]), (3, []))
        self.assertIn("/rows/*", loss["oracle"]["consistency_only"])
        # no sequences in this retrieval: a spelled walk is consistency-only, a slice an error
        sp = rec("graphlet_walks", {"arm": "right", "spell": "full"}, "caps")
        self.assertIn("/rows/*/sequence", sp["oracle"]["consistency_only"])
        self.assertEqual(rec("graphlet_sequence", {"arm": "right", "walk": 0}, "caps")["oracle"]["declared"],
                         "expected_error")

    def test_caps_walk_labels_splits(self):
        w = rec("graphlet_walk", {"arm": "right", "walk": 0}, "caps")
        a = w["answer"]
        self.assertEqual([(r["segment"], r["from_bp"], r["to_bp"], r["labels_in"], r["labels_out"]) for r in a["rows"]],
                         [(0, 0, 20, 2, 2), (1, 20, 126, 1, 1)])
        self.assertEqual(a["rows"][1]["in"], [{"name": "hY", "ref": "h:0:1"}])
        self.assertEqual(a["rows"][1]["out"], [{"name": "hZ", "ref": "h:0:2"}])
        self.assertEqual(a["rows"][1]["ends"], [{"name": "hZ", "ref": "h:0:2", "at_bp": 126, "why": "max_steps"}])
        self.assertEqual(a["continuation"], {"length_bp": 23, "loss_used": 0.5, "labels": 1})
        self.assertIn("/rows/0/ends", w["oracle"]["absent"])
        self.assertIn("/rows/1/merge_of", w["oracle"]["absent"])
        self.assertNotIn("/continuation", w["oracle"]["absent"])
        self.assertIn("/continuation", rec("graphlet_walk", {"arm": "right", "walk": 1}, "caps")["oracle"]["absent"])
        lab = rec("graphlet_labels", {"rank": "direct_bp"}, "caps")["answer"]
        self.assertEqual([(r["name"], r["arm"], r["direct_bp"], r["reach_bp"], r["reentries"], r["runs"])
                          for r in lab["rows"]],
                         [("hX", "right", 109, 109, 0, 2), ("hY", "right", 103, 126, 0, 1), ("hZ", "right", 0, 0, 1, 1)])
        reach = rec("graphlet_labels", {"rank": "reach_bp"}, "caps")["answer"]
        self.assertEqual([r["name"] for r in reach["rows"]], ["hY", "hX", "hZ"])
        # at outward base 103 hY's run [0, 103) has ended and hZ's [103, 126) begun
        at = rec("graphlet_labels", {"at_bp": 103}, "caps")["answer"]
        self.assertEqual([r["name"] for r in at["rows"]], ["hX", "hZ"])
        self.assertEqual(rec("graphlet_labels", {"min_direct_bp": 104}, "caps")["answer"]["total"], 1)
        sp = rec("graphlet_splits", {"arm": "right", "min_labels_before": 0}, "caps")["answer"]
        self.assertEqual(sp["rows"], [
            {"at_bp": 20, "segment": 0, "kind": "divergence", "labels_before": 2,
             "branches": [{"char": "A", "labels": 1}, {"char": "C", "labels": 1}]},
            {"at_bp": 70, "segment": 2, "kind": "ambiguous", "labels_before": 1,
             "branches": [{"char": "C", "labels": 1}, {"char": "G", "labels": 1}]}])
        self.assertEqual(rec("graphlet_splits", {"arm": "right"}, "caps")["answer"]["total"], 1)

    def test_caps_continuation(self):
        """SPEC §7.1: the seed is the continuation's sequence and labels; the loss budget (1.0) lowered by the
        largest terminal loss of the continued labels (hZ: 0.5); max_label_branches unlimited stays unlimited. The
        loss_budget / branch_budget blocks are the tool's presentation (consistency-only); overrides make the
        strategy consistency-only (their merge is library-defined)."""
        r = rec("traverse_continue", {"arm": "right", "walk": 0, "execute": False}, "caps")
        a = r["answer"]
        self.assertEqual(a["request"]["seeds"], [{"sequence": "CCGAGTGATATATCACCTTGCTA", "labels": ["hZ"]}])
        self.assertEqual(a["request"]["strategy"], {"labels": {"loss_budget": 0.5},
                                                    "branching": {"max_label_branches": "unlimited"}})
        self.assertEqual(a["parent"], {"arm": "right", "walk": 0, "overlap_bp": 23})
        self.assertNotIn("loss_budget", a)
        self.assertIn("/loss_budget", r["oracle"]["consistency_only"])
        ov = rec("traverse_continue", {"arm": "right", "walk": 0, "execute": False,
                                       "overrides": {"labels": {"loss_budget": 3.0}}}, "caps")
        self.assertNotIn("strategy", ov["answer"]["request"])
        self.assertIn("/request/strategy/*", ov["oracle"]["consistency_only"])
        self.assertEqual(rec("traverse_continue", {"arm": "right", "walk": 0, "execute": False, "overrides": {"arm": 1}},
                             "caps")["oracle"]["expected_error"], "bad_argument")
        # fork: max_label_branches 0, a terminal branch count 0: 0; reset_branches keeps the original
        f0 = rec("traverse_continue", {"arm": "left", "walk": 0, "execute": False}, "fork")["answer"]["request"]
        self.assertEqual(f0["strategy"]["branching"], {"max_label_branches": 0})
        no = rec("traverse_continue", {"arm": "right", "walk": 1, "execute": False}, "caps")["oracle"]
        self.assertEqual((no["declared"], no["expected_error"]), ("expected_error", None))
        self.assertEqual(rec("traverse_continue", {"arm": "right", "walk": 0}, "caps")["oracle"]["declared"],
                         "not_applicable")

    def test_merge_displayed_evidence(self):
        """DESIGN §5.1: a.fa entered the displayed walk through a non-first parent at 62, b.fa at 36 (their end
        labels' route_bp): neither carries the whole walk; c.fa and both.fa do. These values are re-derived from
        DESIGN's formula (oracle.rederived, DESIGN-v8 §2)."""
        r = rec("graphlet_walks", {"arm": "right"}, "merge")
        self.assertEqual(r["oracle"]["rederived"], ["/rows/0/labels_full", "/rows/0/n_full", "/rows/0/walk"])
        by_length = rec("graphlet_walks", {"arm": "right", "rank": "length"}, "caps")["oracle"]["rederived"]
        self.assertNotIn("/rows/0/walk", by_length)      # the length order uses no DESIGN formula
        self.assertIn("/rows/0/n_full", by_length)
        w = r["answer"]["rows"][0]
        self.assertEqual((w["n_alive"], w["n_full"], w["labels_full"]),
                         (4, 2, [{"name": "c.fa", "ref": "c:2"}, {"name": "both.fa", "ref": "c:3"}]))
        a = rec("graphlet_walks", {"arm": "right", "label": {"name": "a.fa"}}, "merge")["answer"]
        self.assertEqual((a["total"], a["filtered"]), (0, {"merge_entered": 1}))
        a = rec("graphlet_walks", {"arm": "right", "label": {"name": "a.fa"}, "route_consistent": False},
                "merge")["answer"]
        self.assertEqual((a["total"], "filtered" in a), (1, False))
        rows = rec("graphlet_walk", {"arm": "right", "walk": 0}, "merge")["answer"]["rows"]
        self.assertEqual([(r["segment"], r.get("merge_of")) for r in rows],
                         [(0, None), (1, None), (3, [1, 2]), (5, None), (6, [5, 4])])
        s = rec("graphlet_summary", {}, "merge")["answer"]["summary"]["arms"]["right"]
        self.assertEqual((s["walks"], s["segments"], s["splits"], s["merges"], s["bases"], s["max_bp"]),
                         (1, 7, 2, 2, 132, 100))

    def test_fork_sequences_and_exports(self):
        g, f, mgt = doc("fork")
        seed = "CGTACCGTCGTAGCCATGCTGCTTCATTGCAGGTTCTATT"
        q = lambda **a: O.answer("graphlet_sequence", a, f, mgt)["answer"]   # noqa: E731
        self.assertEqual(q(arm="left", walk=0)["sequence"], "ATCCCTCACG")
        self.assertEqual(q(arm="left", walk=0, orientation="walk")["sequence"], "GCACTCCCTA")
        self.assertEqual(q(arm="left", walk=0, with_seed=True)["sequence"], "ATCCCTCACG" + seed)
        self.assertEqual(q(arm="right", walk=1, with_seed=True)["sequence"], seed + "TTACACCACG")
        self.assertEqual(q(arm="right", walk=1, start=2, end=5), {
            "arm": "right", "orientation": "natural", "with_seed": False, "from": 2, "length": 10, "sequence": "ACA",
            "to": 5, "evidence": q(arm="right", walk=1)["evidence"]})
        # SPEC §7.5.4 defines walking order for flanks only: the left arm's walk orientation with the seed is open
        r = O.answer("graphlet_sequence", {"arm": "left", "walk": 0, "orientation": "walk", "with_seed": True}, f, mgt)
        self.assertIn("/sequence", r["oracle"]["consistency_only"])
        self.assertEqual(O.answer("graphlet_sequence", {"arm": "left", "walk": 0, "start": 11}, f, mgt)["oracle"]
                         ["declared"], "expected_error")
        fa = O.answer("graphlet_export", {"format": "fasta"}, f, mgt)["answer"]
        self.assertEqual([x["id"] for x in fa["fasta"]], ["left_0", "left_1", "right_0", "right_1"])
        self.assertEqual(fa["fasta"][0]["sequence"], "ATCCCTCACG" + seed)
        self.assertEqual(fa["records"], 4)
        one = O.answer("graphlet_export", {"format": "fasta", "arm": "right", "walks": [1]}, f, mgt)["answer"]
        self.assertEqual(one["fasta"], [{"id": "right_1", "sequence": seed + "TTACACCACG"}])
        gfa = O.answer("graphlet_export", {"format": "gfa"}, f, mgt)["answer"]["gfa"]
        self.assertEqual(gfa["n_paths"], 4)
        self.assertIn(seed + "ATTGCTAAGA", gfa["paths"])

    def test_subtrie_view_receipt(self):
        r = rec("graphlet_subtrie", {"labels": ["acc2", {"ref": "h:0:2"}], "mode": "all", "arm": "left"}, "fork")
        a = r["answer"]
        self.assertEqual(a["view"], {"selectors": [{"ref": "h:0:1"}, {"ref": "h:0:2"}], "mode": "all", "arm": "left"})
        self.assertEqual(a["summary"]["labels"], [{"name": "acc2", "ref": "h:0:1"}, {"name": "acc3", "ref": "h:0:2"}])
        self.assertEqual(a["summary"]["arms"], {"left": {"complete_to_bp": 10, "scope": "united_history",
                                                         "qualified": "for the selected labels"}})
        self.assertIn("/summary/arms/*/walks", r["oracle"]["consistency_only"])
        self.assertEqual(rec("graphlet_subtrie", {"labels": ["nope"]}, "fork")["oracle"]["declared"], "expected_error")
        self.assertEqual(rec("graphlet_walks", {"arm": "left"}, "fork", view=True)["oracle"]["declared"],
                         "view_consistency_only")


class ClaimsGoldens(unittest.TestCase):
    """graphlet_claims (phase B): one claim per run end (SPEC §6.9), the engine's run fields, DESIGN §5.1's displayed
    evidence (re-derived), read off the CLI documents and DESIGN §5.1's own clipped-merge example."""

    def test_merge_evidence_and_filter(self):
        """merge: a.fa enters the displayed walk through the second parent at 62, b.fa at 36 (their leaf stamps);
        both.fa's runs closed by the merges at 36 and 62 (code m) carry their route; route_consistent leaves a.fa
        and b.fa out as merge_entered."""
        r = rec("graphlet_claims", {"arm": "right", "route_consistent": False}, "merge")
        rows = r["answer"]["rows"]
        self.assertEqual([(x["label"]["name"], x["from_bp"], x["to_bp"], x["evidence_from"]) for x in rows],
                         [("a.fa", 0, 100, 62), ("b.fa", 0, 100, 36), ("c.fa", 0, 100, 0), ("both.fa", 0, 62, 0),
                          ("both.fa", 0, 36, 0), ("both.fa", 0, 100, 0)])
        self.assertEqual([x.get("reason") for x in rows], ["max_extension_bp"] * 3 + [None, None, "max_extension_bp"])
        self.assertEqual([x.get("path_id") for x in rows], [0, 0, 0, None, None, 0])
        self.assertIn("/rows/0/evidence_from", r["oracle"]["rederived"])
        self.assertIn("/rows/3/reason", r["oracle"]["consistency_only"])     # a merge's close states no end reason
        on = rec("graphlet_claims", {"arm": "right"}, "merge")["answer"]
        self.assertEqual((on["total"], on["filtered"]), (4, {"merge_entered": 2}))
        one = rec("graphlet_claims", {"arm": "right", "labels": "both.fa", "route_consistent": False}, "merge")
        self.assertEqual(one["answer"]["total"], 3)
        mb = rec("graphlet_claims", {"arm": "right", "min_bp": 50, "route_consistent": False}, "merge")["answer"]
        # a.fa shows 100 - 62 = 38 < 50 and both.fa's run to 36 shows 36: out; b.fa shows 64, c.fa 100, both.fa 62, 100
        self.assertEqual([(x["label"]["name"], x["to_bp"]) for x in mb["rows"]],
                         [("b.fa", 100), ("c.fa", 100), ("both.fa", 62), ("both.fa", 100)])
        self.assertIn("/filtered", rec("graphlet_claims", {"arm": "right", "min_bp": 1}, "merge")["oracle"]
                      ["consistency_only"])                                          # the filters' order (PB3)

    def test_route_filter_on_stamps(self):
        """VO-B2: route_consistent filters on the claim's SPEC §7.1 merge stamp (T at a leaf, R inside), as phase A's
        walks filter does (F2), not on where the displayed evidence begins. merge, with c.fa's leaf stamps raised to
        10 (T and R): filtered, although its evidence begins at 0; an inner run whose R stamp and displayed stamp
        disagree (both.fa's run closed at 36: R raised to 5, displayed 0) leaves the selection open."""
        g, f, mgt = doc("merge")
        f = copy.deepcopy(f)
        arm = f["results"][0]["arms"]["right"]
        arm["runs"][2]["route_bp"] = 10
        for e in arm["paths"][0]["end_labels"]:
            if e["run"] == 2:
                e["route_bp"] = 10
        r = O.answer("graphlet_claims", {"arm": "right"}, f, mgt)["answer"]
        self.assertEqual((r["filtered"], [x["label"]["name"] for x in r["rows"]]),
                         ({"merge_entered": 3}, ["both.fa", "both.fa", "both.fa"]))
        arm["runs"][4]["route_bp"] = 5
        open_ = O.answer("graphlet_claims", {"arm": "right"}, f, mgt)
        # R2B-3: the rows before the undecided claim stay checked; it, every later row, the total and the counts
        # are consistency-only, and the list has no length check
        self.assertEqual([(x["label"]["name"], x["to_bp"]) for x in open_["answer"]["rows"]], [("both.fa", 62)])
        self.assertIn("/rows/0/to_bp", open_["oracle"]["checked"])
        self.assertTrue({"/total", "/filtered"} <= set(open_["oracle"]["consistency_only"]))
        self.assertNotIn("total", open_["answer"])
        self.assertNotIn("/rows/1", open_["oracle"].get("absent", {}))
        extra = copy.deepcopy(open_["answer"])
        extra["rows"].append({"label": {"name": "both.fa"}})
        self.assertEqual(O.differences(open_, extra), [])
        self.assertEqual(len(O.answer("graphlet_claims", {"arm": "right", "route_consistent": False}, f, mgt)
                             ["answer"]["rows"]), 6)          # without the filter nothing is open

    def test_switch_and_boundary(self):
        """caps: hZ entered by a switch from hY at 103 (cost 0.5) holds [103, 126); linear: acc2 present at the
        right boundary only, a zero-length run: a boundary claim without displayed bases (DESIGN §5.1, PB1)."""
        rows = rec("graphlet_claims", {"arm": "right"}, "caps")["answer"]["rows"]
        hz = rows[3]
        self.assertEqual((hz["label"]["name"], hz["entered_by"], hz["from_label"], hz["cost"], hz["evidence_from"],
                          hz["loss"], hz["qualifier"]),
                         ("hZ", "switch", {"name": "hY", "ref": "h:0:1"}, 0.5, 103, 0.5, None))
        self.assertEqual(rows[1]["reason"], "label_lost")
        lin = rec("graphlet_claims", {"arm": "right", "route_consistent": False}, "linear")
        b = lin["answer"]["rows"][1]
        self.assertEqual((b["label"]["name"], b["from_bp"], b["to_bp"], b["evidence_from"], b["kind"], b["reason"]),
                         ("acc2", 0, 0, None, "boundary", "branch"))
        self.assertIn("/rows/0/kind", lin["oracle"]["consistency_only"])            # no clause for the other kinds
        quorum = rec("graphlet_claims", {"arm": "right", "route_consistent": False}, "quorum")["answer"]["rows"]
        self.assertIn("minority", [x["qualifier"] for x in quorum])                 # R's Rm (SPEC §7.5.1)

    def test_design_clipped_merge_example(self):
        """DESIGN §5.1's clipped-merge case (documents/clipped_merge: HAND-MADE, its claims vector hand-derived from the
        walker's semantics, not the library): lineB enters the displayed path through the merge's second parent at
        14, so its displayed evidence starts there; lineA's at 0."""
        g, f, mgt = doc("clipped_merge")
        with open(os.path.join(DATA, "documents", "clipped_merge.claims.json"), encoding="utf-8") as fh:
            vec = json.load(fh)["claims"]["all"]
        rows = O.answer("graphlet_claims", {"arm": "right", "route_consistent": False}, f, mgt)["answer"]["rows"]
        self.assertEqual([(x["label"]["ref"], x["from_bp"], x["to_bp"], x["evidence_from"], x.get("route_bp"),
                           x.get("reason"), x["loss"], x["branches"]) for x in rows],
                         [(v["ref"], v["from_bp"], v["to_bp"], v["evidence_from"], v["route_bp"], v["reason"], v["loss"],
                           v["branches"]) for v in vec])

    def test_annotate_claims_are_consistency_only(self):
        r = rec("graphlet_claims", {"arm": "right"}, "annotate")
        self.assertEqual((r["oracle"]["checked"].keys() & {"/total"}, "/rows" in r["oracle"]["consistency_only"]),
                         (set(), True))


class LabelWalkGoldens(unittest.TestCase):
    """graphlet_labels(name=) (phase B): the named label's runs with their routes (DESIGN §5 label_walks, routes)."""

    def test_merge_routes_and_closes(self):
        """merge: both.fa's lineage splits at 20 (A and T) and its copies are closed by the merges at 36 and 62 (code
        m: merged into segments 3 and 6); the run that reaches the radius takes the displayed parents."""
        r = rec("graphlet_labels", {"name": "both.fa"}, "merge")
        rows = r["answer"]["rows"]
        self.assertEqual([(x["arm"], x["from_bp"], x["to_bp"], x["merged_into"], x["route"]) for x in rows],
                         [("right", 0, 62, 6, [0, 1, 3, 4]), ("right", 0, 36, 3, [0, 2]),
                          ("right", 0, 100, None, [0, 1, 3, 5, 6])])
        self.assertNotIn("walks_below", rows[0])            # DESIGN names leaves_below without defining it (VO-B3)
        self.assertIn("/rows/*/walks_below", r["oracle"]["consistency_only"])
        self.assertEqual(rows[2]["end"], "max_extension_bp")
        self.assertIn("/rows/0/end", r["oracle"]["consistency_only"])      # a merge's close: no label_end event
        self.assertEqual(set(r["oracle"]["rederived"]) >= {"/rows/0/route", "/rows/0/evidence_from"}, True)
        a = rec("graphlet_labels", {"name": {"ref": "c:0"}}, "merge")["answer"]["rows"]
        # a.fa enters segment 6 through its second parent (partitions [[1, 2, 3], [0]]): its own route passes 4
        self.assertEqual((a[0]["evidence_from"], a[0]["route"]), (62, [0, 1, 3, 4, 6]))

    def test_design_clipped_merge_routes(self):
        """DESIGN §5.1's clipped-merge case (HAND-MADE fixture and vector): lineB's own route passes the merge's second
        parent (segments 0, 2, 3), lineA's the first (0, 1, 3)."""
        g, f, mgt = doc("clipped_merge")
        with open(os.path.join(DATA, "documents", "clipped_merge.claims.json"), encoding="utf-8") as fh:
            routes = json.load(fh)["routes"]
        for ref, want in routes.items():
            rows = O.answer("graphlet_labels", {"name": {"ref": ref}}, f, mgt)["answer"]["rows"]
            self.assertEqual([x["route"] for x in rows], want, ref)

    def test_switch_lineages(self):
        """caps: hZ, switched in from hY at 103 on segment 1, its route through hY's segments; switch_chain: a silent
        switch-source end (Lw) has no SPEC word for its end (consistency-only); both arms of fork, left first."""
        hz = rec("graphlet_labels", {"name": "hZ"}, "caps")["answer"]["rows"]
        self.assertEqual([(x["from_bp"], x["to_bp"], x["evidence_from"], x["route"], x["end"]) for x in hz],
                         [(103, 126, 103, [0, 1], "max_steps")])
        sc = rec("graphlet_labels", {"name": {"id": 1}}, "switch_chain")
        self.assertIn("/rows/0/end", sc["oracle"]["consistency_only"])
        both = rec("graphlet_labels", {"name": "acc1"}, "fork")["answer"]
        self.assertEqual([x["arm"] for x in both["rows"]], ["left", "right"])
        self.assertEqual((both["name"], both["ref"], both["total"]), ("acc1", "h:0:0", 2))
        self.assertEqual(rec("graphlet_labels", {"name": "nope"}, "fork")["oracle"]["expected_error"], "unknown_label")
        # documents/same_name records two headers ACC1 from different columns: a bare name is ambiguous (DESIGN §5
        # Graphlet.label: select by ref), the refs are not
        self.assertEqual(rec("graphlet_labels", {"name": "ACC1"}, "same_name")["oracle"]["expected_error"],
                         "ambiguous_label")
        self.assertNotIn("declared", rec("graphlet_labels", {"name": {"ref": "h:0:0"}}, "same_name")["oracle"])
        self.assertEqual(rec("graphlet_labels", {"name": "acc1", "coordinates": True}, "fork")["oracle"]
                         ["expected_error"], "no_coordinates")


class SupportGoldens(unittest.TestCase):
    """graphlet_support (phase B): the displayed support profile along a walk (DESIGN §5, §6; re-derived from DESIGN
    §5.1 with the walk as the displayed path)."""

    def test_merge_profile(self):
        """merge, walk 0 (segments 0, 1, 3, 5, 6), as W5's claims of the same document say (DESIGN §5.1): a.fa's kept
        entry came through the second parent at 62, b.fa's at 36, so the displayed bases before those merges are not
        theirs; c.fa and both.fa (its runs closed at 36 and 62 kept where their chains are the walk) span it."""
        r = rec("graphlet_support", {"arm": "right", "walk": 0}, "merge")
        rows = r["answer"]["rows"]
        name = lambda xs: [x["name"] for x in xs]   # noqa: E731
        self.assertEqual([(x["from_bp"], x["to_bp"], x["n_labels"]) for x in rows], [(0, 36, 2), (36, 62, 3), (62, 100, 4)])
        self.assertEqual([(name(x.get("joined", [])), name(x.get("left", []))) for x in rows[1:]],
                         [(["b.fa"], []), (["a.fa"], [])])
        claims = {x["label"]["name"]: x["evidence_from"] for x in rec("graphlet_claims", {
            "arm": "right", "route_consistent": False}, "merge")["answer"]["rows"] if x["to_bp"] == 100}
        self.assertEqual((claims["a.fa"], claims["b.fa"]), (rows[2]["from_bp"], rows[1]["from_bp"]))
        self.assertIn("/rows/0/n_labels", r["oracle"]["rederived"])
        stepped = rec("graphlet_support", {"arm": "right", "walk": 0, "step": 10}, "merge")["oracle"]
        self.assertEqual(("/rows" in stepped["consistency_only"], stepped["checked"].keys() & {"/total"}),
                         (True, set()))                     # step=: the tool's merging rule (consistency-only)

    def test_switch_profile(self):
        """caps, walk 0 (segments 0, 1): hX takes C at the split at 20; hY switches to hZ at 103 on segment 1."""
        rows = rec("graphlet_support", {"arm": "right", "walk": 0}, "caps")["answer"]["rows"]
        self.assertEqual([(x["from_bp"], x["to_bp"], x["n_labels"], [y["name"] for y in x.get("joined", [])],
                           [y["name"] for y in x.get("left", [])]) for x in rows],
                         [(0, 20, 2, [], []), (20, 103, 1, [], ["hX"]), (103, 126, 1, ["hZ"], ["hY"])])

    def test_design_clipped_merge_profile(self):
        """DESIGN §5.1's clipped-merge case (HAND-MADE fixture and vector): the displayed profile."""
        g, f, mgt = doc("clipped_merge")
        with open(os.path.join(DATA, "documents", "clipped_merge.claims.json"), encoding="utf-8") as fh:
            want = json.load(fh)["support_profile_displayed"]
        rows = O.answer("graphlet_support", {"arm": "right", "walk": 0}, f, mgt)["answer"]["rows"]
        self.assertEqual([[x["from_bp"], x["to_bp"], x["n_labels"]] for x in rows], [[a, b, len(c)] for a, b, c in want])
        # which labels, not only how many (VO-B5): the per-base sets W5 derives, as refs, against the vector's
        d = O._Doc(f, mgt)
        sets = O._support_sets(d, "right", d.arms["right"]["paths"][0])
        got = []
        for b, cur in enumerate(sets):
            refs = sorted(d.labels[x]["ref"] for x in cur)
            if got and got[-1][2] == refs:
                got[-1][1] = b + 1
            else:
                got.append([b, b + 1, refs])
        self.assertEqual(got, [[a, b, sorted(c)] for a, b, c in want])
        self.assertEqual([x["ref"] for x in rows[1]["joined"]], sorted(set(want[1][2]) - set(want[0][2])))

    def test_errors_and_annotate(self):
        self.assertEqual(rec("graphlet_support", {"arm": "right", "walk": 3}, "caps")["oracle"]["expected_error"],
                         "bad_argument")
        # annotate: the labels the engine recorded at the walk's nodes (SPEC §7.4 label_sets), checked, not
        # re-derived; documents/annotate cuts lists at its cap, so who joined or left is not read across them
        g, f, mgt = doc("annotate")
        r = O.answer("graphlet_support", {"arm": "right", "walk": 0}, f, mgt)
        sets = [ls for sid in f["results"][0]["arms"]["right"]["paths"][0]["segments"]
                for ls in f["results"][0]["arms"]["right"]["segments"][sid].get("label_sets") or []]
        self.assertEqual(r["answer"]["rows"][0]["n_labels_total"], sets[0]["labels_total"])
        self.assertIn("/rows/0/n_labels", r["oracle"]["checked"])
        self.assertNotIn("rederived", r["oracle"])
        self.assertTrue(any(c.endswith("/joined") for c in r["oracle"]["consistency_only"] if c != "/rows/0/joined"))


class RecordSemantics(unittest.TestCase):
    """The record contract (FORMAT-v8 §4 and its extensions) and that W5 is not vacuous: a wrong answer differs."""

    def test_declared_records(self):
        g, f, mgt = doc("fork")
        for tool, args, why in [
                                ("graphlet_compare", {"mode": "walks"}, "properties_only"),
                                ("graphlet_list", {}, "not_applicable"), ("traverse_fetch", {}, "not_applicable"),
                                ("graphlet_walk", {"arm": "left", "walk": 2}, "expected_error"),
                                ("graphlet_walks", {"rank": "nope", "arm": "left"}, "expected_error"),
                                ("graphlet_walks", {}, "expected_error")]:
            r = O.answer(tool, args, f, mgt)
            self.assertEqual((r["answer"], r["oracle"]["declared"], r["oracle"]["checked"]), (None, why, {}), tool)
        err = O.answer("graphlet_walk", {"arm": "left", "walk": 2}, f, mgt)
        self.assertEqual(O.differences(err, {"error": "bad_argument", "message": "..."}), [])
        self.assertEqual([d["pointer"] for d in O.differences(err, {"walk": 2, "rows": []})], ["/error"])

    def test_wrong_answers_differ(self):
        """Mutations of W5's own answer, one checked value each: every one is named with its pointer and clause."""
        g, f, mgt = doc("merge")
        for tool, args, path, value in [
                ("graphlet_walks", {"arm": "right"}, ["rows", 0, "n_full"], 3),
                ("graphlet_walks", {"arm": "right"}, ["rows", 0, "labels_full", 0, "ref"], "c:0"),
                ("graphlet_walk", {"arm": "right", "walk": 0}, ["rows", 2, "merge_of"], [2, 1]),
                ("graphlet_walk", {"arm": "right", "walk": 0}, ["continuation", "labels"], 3),
                ("graphlet_labels", {}, ["rows", 0, "reach_bp"], 99),
                ("graphlet_splits", {"arm": "right", "min_labels_before": 0}, ["rows", 1, "kind"], "divergence"),
                ("graphlet_summary", {}, ["summary", "arms", "right", "merges"], 1),
                ("graphlet_summary", {}, ["evidence", "scope", "right"], "per_path"),
                ("graphlet_sequence", {"arm": "right", "walk": 0}, ["sequence"], "TTT"),
                ("graphlet_export", {"format": "json"}, ["json", "label_summary", 0, "right", "direct_bp"], 1)]:
            r = O.answer(tool, args, f, mgt)
            other = copy.deepcopy(r["answer"])
            node = other
            for k in path[:-1]:
                node = node[k]
            node[path[-1]] = value
            diffs = O.differences(r, other)
            self.assertEqual(len(diffs), 1, (tool, path, diffs))
            self.assertTrue(diffs[0]["clause"].startswith(("SPEC", "DESIGN")), diffs[0])
        # a missing key, a reordered list, a key that must be absent
        r = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, *doc("caps")[1:])
        other = copy.deepcopy(r["answer"])
        other["rows"].reverse()
        self.assertGreaterEqual(len(O.differences(r, other)), 4)
        del other["total"]
        self.assertIn("/total", [d["pointer"] for d in O.differences(r, other) if d.get("missing")])
        w = O.answer("graphlet_walk", {"arm": "right", "walk": 1}, *doc("caps")[1:])
        other = dict(copy.deepcopy(w["answer"]), continuation={"length_bp": 1, "loss_used": 0, "labels": 1})
        self.assertEqual([d["pointer"] for d in O.differences(w, other)], ["/continuation"])

    def test_compare_entry(self):
        """FORMAT v2 F2: compare(record_w5, answer) on a per-call record or a bare answer; F3: a fitting error is
        not_comparable; declarations are not_checked; an error where W5 answers disagrees; expected codes."""
        g, f, mgt = doc("caps")
        r = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, f, mgt)
        self.assertEqual(O.compare(r, r["answer"])["outcome"], "agree")
        line = {"way": "W1", "op": "w", "args_resolved": {}, "answer": copy.deepcopy(r["answer"]), "error": None}
        self.assertEqual(O.compare(r, line)["outcome"], "agree")
        line["answer"]["rows"][0]["length_bp"] = 1
        v = O.compare(r, line)
        self.assertEqual((v["outcome"], v["differences"][0]["pointer"]), ("disagree", "/rows/0/length_bp"))
        fit = {"way": "W3", "answer": None, "error": {"code": "result_too_large", "fields": {"max_bytes": 512}}}
        self.assertEqual(O.compare(r, fit)["outcome"], "not_comparable")
        self.assertEqual(O.compare(r, {"error": "op_busy"})["outcome"], "not_comparable")
        self.assertEqual(O.compare(r, {"error": "bad_arm"})["outcome"], "disagree")
        err = O.answer("graphlet_walks", {"arm": "right", "n": 0}, f, mgt)
        self.assertEqual(err["oracle"]["expected_error"], "bad_argument")
        self.assertEqual(O.compare(err, {"way": "W1", "answer": None, "error": {"code": "bad_argument", "fields": {}}})
                         ["outcome"], "agree")
        self.assertEqual(O.compare(err, {"error": "bad_arm"})["outcome"], "disagree")
        self.assertEqual(O.compare(err, r["answer"])["outcome"], "disagree")
        ph = O._declared("not_implemented_phase_b")
        self.assertEqual(O.compare(ph, {"rows": []})["outcome"], "not_checked")
        # no bases (sequences: false): the tool refuses a FASTA (library-defined, any code); a GFA may carry `*`
        # sequences, so W5 makes no claim there; a record checking nothing makes no claim
        fa = O.answer("graphlet_export", {"format": "fasta"}, f, mgt)
        self.assertEqual((fa["oracle"]["declared"], fa["oracle"]["expected_error"]), ("expected_error", None))
        self.assertIn("library-defined", fa["oracle"]["note"])
        gfa = O.answer("graphlet_export", {"format": "gfa"}, f, mgt)
        self.assertEqual((gfa["oracle"]["checked"], O.compare(gfa, {"error": "bad_argument"})["outcome"]),
                         ({}, "not_checked"))
        nothing = O.answer("graphlet_export", {"format": "gfa", "arm": "right"}, f, mgt)
        self.assertEqual(O.compare(nothing, {"error": "bad_argument"})["outcome"], "not_checked")

    def test_declarations_in_compare(self):
        """R2-5: a declaration by which W5 makes no claim is not_checked; one by which it could not judge is
        inconclusive (it blocks); so is anything that is no W5 record."""
        for why in O.NO_CLAIM:
            self.assertEqual(O.compare(O._declared(why), {"rows": []})["outcome"], "not_checked", why)
        for why in O.CANNOT_JUDGE:
            self.assertEqual(O.compare(O._declared(why), {"rows": []})["outcome"], "inconclusive", why)
        self.assertEqual(O.compare(O._declared("something_new"), {})["outcome"], "inconclusive")
        self.assertEqual(O.compare({"answer": None}, {})["outcome"], "inconclusive")
        # a twin mismatch from answer() itself
        g = _load(os.path.join(DATA, "real", "sra", "sra_16s_OK510341__time_tight.graphlet.json.gz"))
        f = _load(os.path.join(DATA, "real", "sra", "sra_16s_OK510341__time_tight.full.json.gz"))
        r = O.answer("graphlet_summary", {}, f, g["results"][0]["graphlet"])
        self.assertEqual(O.compare(r, {"summary": {}})["outcome"], "inconclusive")

    def test_failed_seed(self):
        """VR4: a seed the full twin states failed (SPEC §7.1: {seed, outcome walks failed, error}, no arms) has no
        graphlet: every op that reads one expects no_graphlet (lane L's code for the same answer), without any
        document; graphlet_load any code; backend tools stay not applicable."""
        f = _load(os.path.join(DATA, "real", "sra", "sra_batch3__batch_refuse.full.json.gz"))
        self.assertTrue(O.seed_failed(f, 1) and not O.seed_failed(f, 0))
        for tool, args in [("graphlet_walks", {"arm": "left"}), ("graphlet_summary", {}),
                           ("graphlet_export", {"format": "mgt"}), ("graphlet_claims", {"arm": "left"}),
                           ("graphlet_compare", {"mode": "walks"})]:
            r = O.answer(tool, args, f, None, {"seed_index": 1})
            self.assertEqual((r["oracle"]["declared"], r["oracle"]["expected_error"]), ("expected_error", "no_graphlet"))
            self.assertIn("SPEC §7.1", r["oracle"]["note"])
            w1 = {"way": "W1", "args_resolved": dict(args, handle="<handle>"), "answer": None,
                  "error": {"code": "no_graphlet", "fields": {}}}
            self.assertEqual(O.compare(r, w1)["outcome"], "agree", tool)
            self.assertEqual(O.compare(r, {"error": "bad_argument"})["outcome"], "disagree", tool)
        self.assertIsNone(O.answer("graphlet_load", {"path": "x.mgt"}, f, None, {"seed_index": 1})["oracle"]
                          ["expected_error"])
        self.assertEqual(O.answer("traverse_fetch", {}, f, None, {"seed_index": 1})["oracle"]["declared"],
                         "not_applicable")

    def test_refused_arguments(self):
        """F-D2: a path outside the export directory (path_not_allowed), a budget without local limits (bad_argument),
        a cursor of another listing (bad_cursor): expected errors, their codes library-defined; what W5 cannot
        judge (a subdirectory path, a budget under local limits, a resumed page) makes no claim."""
        g, f, mgt = doc("fork")
        code = lambda tool, args, **ctx: O.answer(tool, args, f, mgt, ctx)["oracle"]   # noqa: E731
        for tool, path in [("graphlet_export", "/etc/x.fa"), ("graphlet_export", "../x.fa"), ("graphlet_save", "/tmp/x.mgt"),
                           ("graphlet_load", "../../x.mgt"), ("graphlet_export", "<export_dir>/../x.fa")]:
            o = code(tool, {"format": "fasta", "path": path})
            self.assertEqual((o["declared"], o["expected_error"]), ("expected_error", "path_not_allowed"), path)
            self.assertIn("library-defined", o["note"])
        for tool, path in [("graphlet_export", "x.fa"), ("graphlet_export", "<export_dir>/x.fa"),
                           ("graphlet_load", "<spool_dir>/x.mgt"), ("graphlet_export", "sub/../a.fa"),
                           ("graphlet_export", "<export_dir>/sub/../a.fa")]:
            self.assertNotIn("declared", code(tool, {"format": "fasta", "path": path}), path)
        # R2B-2: resolved against the export directory the context names: inside is answered, outside refused
        self.assertNotIn("declared", code("graphlet_export", {"format": "fasta", "path": "/w/e/b.fa"}, export_dir="/w/e"))
        self.assertEqual(code("graphlet_export", {"format": "fasta", "path": "/w/e/../b.fa"}, export_dir="/w/e")
                         ["expected_error"], "path_not_allowed")
        self.assertEqual(code("graphlet_export", {"format": "fasta", "path": "/w/f/b.fa"}, export_dir="/w/e")
                         ["expected_error"], "path_not_allowed")
        for path in ("sub/x.fa", "~/c.fa", "", "<export_dir>/sub/x.fa", "<other>/x.fa"):
            self.assertEqual(code("graphlet_export", {"format": "fasta", "path": path})["declared"], "path_unknown", path)
        o = code("graphlet_walks", {"arm": "left", "budget": {"work_units": 10}})
        self.assertEqual((o["declared"], o["expected_error"]), ("expected_error", "bad_argument"))
        self.assertEqual(code("graphlet_walks", {"arm": "left", "budget": {"work_units": 10}},
                              local_limits={"view": {}})["declared"], "local_limits")
        self.assertNotIn("declared", code("graphlet_walks", {"arm": "left", "budget": None}))
        self.assertEqual(code("graphlet_walks", {"arm": "left", "cursor": "<cursor>"})["declared"], "cursor_page")
        o = code("graphlet_walks", {"arm": "left", "cursor": "<cursor>"}, cursor_bound_elsewhere=True)
        self.assertEqual((o["declared"], o["expected_error"]), ("expected_error", "bad_cursor"))
        # R2B-1: the cursor's source op decides on normalised arguments (arm short forms, labels resolved, DESIGN
        # §6's defaults): only another tool or another value is bad_cursor; the rest makes no claim
        src = {"tool": "graphlet_walks", "args": {"handle": "<handle>", "arm": "right", "rank": "id", "n": 1}}
        for extra, want in [({"arm": "r"}, "cursor_page"), ({"arm": "right", "spell": "none"}, "cursor_page"),
                            ({"arm": "right", "n": 3}, "cursor_page"), ({"arm": "right", "rank": "length"}, "bad_cursor"),
                            ({"arm": "right", "min_bp": 1}, "bad_cursor"), ({"arm": "left"}, "bad_cursor")]:
            args = dict({"handle": "<handle>", "rank": "id", "n": 1, "cursor": "<cursor>"}, **extra)
            o = code("graphlet_walks", args, cursor_source=src)
            self.assertEqual(o.get("expected_error") or o["declared"], want, extra)
        self.assertEqual(code("graphlet_labels", {"cursor": "<cursor>"}, cursor_source=src)["expected_error"], "bad_cursor")
        by_name = {"tool": "graphlet_walks", "args": {"arm": "right", "label": "acc1"}}
        self.assertEqual(code("graphlet_walks", {"arm": "right", "label": {"ref": "h:0:0"}, "cursor": "<cursor>"},
                              cursor_source=by_name)["declared"], "cursor_page")       # one label, named two ways
        no_arm = {"tool": "graphlet_walks", "args": {"rank": "id"}}
        self.assertEqual(code("graphlet_walks", {"arm": "right", "rank": "id", "cursor": "<cursor>"},
                              cursor_source=no_arm)["declared"], "cursor_page")      # arm has no documented default
        for why in ("cursor_page", "path_unknown", "local_limits"):
            self.assertEqual(O.compare(O._declared(why), {"rows": []})["outcome"], "not_checked")

    def test_patterns_and_position_order(self):
        """`*` patterns in checked, unordered and absent expand over the answers (FORMAT v2 F1); with
        ENDS_UNORDERED off (G3) a row's ends compare sorted by position only (SPEC §7.1)."""
        rec_ = {"answer": {"rows": [{"x": 1, "e": [{"at_bp": 2, "l": "a"}, {"at_bp": 2, "l": "b"}]},
                                    {"x": 2, "e": []}]},
                "oracle": {"checked": {"/rows/*/x": "c", "/rows/*/e": "c"}, "position_ordered": ["/rows/*/e"],
                           "absent": {"/rows/*/bad": "c"}}}
        same_set = {"rows": [{"x": 1, "e": [{"at_bp": 2, "l": "b"}, {"at_bp": 2, "l": "a"}]}, {"x": 2, "e": []}]}
        self.assertEqual(O.differences(rec_, same_set), [])
        self.assertEqual([d["pointer"] for d in O.differences(rec_, dict(same_set, rows=[
            dict(same_set["rows"][0], bad=1), same_set["rows"][1]]))], ["/rows/0/bad"])
        self.assertEqual(sorted(d["pointer"] for d in O.differences(rec_, {"rows": [{"x": 1}, {"x": 3, "e": []}]})),
                         ["/rows/0/e", "/rows/1/x"])
        g, f, mgt = doc("merge")
        old = O.ENDS_UNORDERED
        try:
            O.ENDS_UNORDERED = False
            r = O.answer("graphlet_walk", {"arm": "right", "walk": 0}, f, mgt)
        finally:
            O.ENDS_UNORDERED = old
        self.assertIn("/rows/4/ends", r["oracle"]["position_ordered"])
        other = copy.deepcopy(r["answer"])
        other["rows"][4]["ends"].reverse()      # all four at 100: ties unordered
        self.assertEqual(O.differences(r, other), [])
        other["rows"][4]["ends"][0] = dict(other["rows"][4]["ends"][0], at_bp=99)
        self.assertEqual([d["pointer"] for d in O.differences(r, other)], ["/rows/4/ends"])

    def test_list_lengths(self):
        """VO7: a whole list's length is checked (an absent pointer at its end): a repeated row, an extra branch or
        export record is a difference; rows of a listing W1 did not read to its end are not."""
        g, f, mgt = doc("merge")
        for tool, args, path in [("graphlet_walks", {"arm": "right"}, ["rows"]),
                                 ("graphlet_walk", {"arm": "right", "walk": 0}, ["rows"]),
                                 ("graphlet_labels", {}, ["rows"]),
                                 ("graphlet_splits", {"arm": "right", "min_labels_before": 0}, ["rows"]),
                                 ("graphlet_splits", {"arm": "right", "min_labels_before": 0}, ["rows", 0, "branches"]),
                                 ("graphlet_summary", {"detail": "limitations"}, ["rows"]),
                                 ("graphlet_export", {"format": "fasta"}, ["fasta"])]:
            r = O.answer(tool, args, f, mgt)
            other = copy.deepcopy(r["answer"])
            node = other
            for k in path:
                node = node[k]
            node.append(copy.deepcopy(node[-1]))
            self.assertEqual(len(O.differences(r, other)), 1, (tool, path))
        part = O.answer("graphlet_walks", {"arm": "right"}, f, mgt, {"paging": "none", "rows": 1})
        other = copy.deepcopy(part["answer"])
        other["rows"].append({"walk": 9})
        self.assertEqual(O.differences(part, other), [])

    def test_self_compare_rules(self):
        """VO6: without index_fp a graphlet with bases and depth compares unverifiable only; unknown is accepted where
        a rule of unknown also applies (no bases in a spelling mode: caps; complete_to_bp 0: no_bins)."""
        base = {"equal": None, "counts": {"only_in_a": 0, "only_in_b": 0, "differ": 0}, "rows": []}
        _, f, mgt = doc("fork")
        self.assertEqual(O.compare_self_properties(dict(base, mode="walks", comparable="unverifiable"), f, mgt), [])
        self.assertEqual(len(O.compare_self_properties(dict(base, mode="walks", comparable="unknown"), f, mgt)), 1)
        self.assertEqual(len(O.compare_self_properties(dict(base, mode="walks", comparable="maybe"), f, mgt)), 2)
        _, f, mgt = doc("caps")
        for c in ("unknown", "unverifiable"):
            self.assertEqual(O.compare_self_properties(dict(base, mode="walks", comparable=c), f, mgt), [])
        self.assertEqual(len(O.compare_self_properties(dict(base, mode="labels", comparable="unknown"), f, mgt)), 1)
        _, f, mgt = doc("no_bins")
        self.assertEqual(O.compare_self_properties(dict(base, mode="labels", comparable="unknown"), f, mgt), [])

    def test_unordered_lists(self):
        """The unordered reading (ENDS_UNORDERED on, the pre-G3 library's run order) compares ends as a set."""
        g, f, mgt = doc("merge")
        old = O.ENDS_UNORDERED
        try:
            O.ENDS_UNORDERED = True
            r = O.answer("graphlet_walk", {"arm": "right", "walk": 0}, f, mgt)
        finally:
            O.ENDS_UNORDERED = old
        self.assertIn("/rows/4/ends", r["oracle"]["unordered"])
        other = copy.deepcopy(r["answer"])
        other["rows"][4]["ends"].reverse()     # run order instead of (at_bp, label): the same set (F5)
        self.assertEqual(O.differences(r, other), [])
        other["rows"][4]["ends"].pop()
        self.assertEqual([d["pointer"] for d in O.differences(r, other)], ["/rows/4/ends"])

    def test_cut_rows(self):
        """A row cut to fit its page names the cut fields (DESIGN §6; FORMAT v2 A3: `_cut_fields` on each row of such
        a page): each must be a cut of W5's value in the suite's words (REVIEW D-R8: a prefix; a set field a
        subset); a top-level marker alone does not excuse a row."""
        g, f, mgt = doc("merge")
        r = O.answer("graphlet_walk", {"arm": "right", "walk": 0}, f, mgt)
        other = copy.deepcopy(r["answer"])
        row = other["rows"][4]
        row["_cut_fields"] = ["ends", "in"]
        row["in"] = row["in"][:2]
        row["ends"] = row["ends"][2:3]
        self.assertEqual(O.differences(r, other), [])
        row["in"] = row["in"][1:]
        self.assertEqual([d["pointer"] for d in O.differences(r, other)], ["/rows/4/in"])
        top = dict(copy.deepcopy(r["answer"]), row_truncated=True, cut_fields=["out"])
        top["rows"][0]["out"] = top["rows"][0]["out"][:1]
        self.assertEqual([d["pointer"] for d in O.differences(r, top)], ["/rows/0/out"])
        top["rows"][0]["_cut_fields"] = ["out"]
        self.assertEqual(O.differences(r, top), [])

    def test_public_cut_rule(self):
        """is_cut_of / cut_problem, the one cut rule lane C also uses (VERIFY-C VC1)."""
        full = {"segment": 3, "in": [1, 2, 3], "ends": [{"at_bp": 5, "l": 1}, {"at_bp": 2, "l": 2}], "seq": "ACGT"}
        cut = {"segment": 3, "in": [1, 2], "ends": [{"at_bp": 2, "l": 2}], "seq": "AC", "_cut_fields": ["in", "seq"]}
        self.assertTrue(O.is_cut_of(cut, full, ["in", "seq", "ends"], unordered=["ends"]))
        self.assertIsNone(O.cut_problem(cut, full, ["in", "seq", "ends"], position_ordered=["ends"]))
        self.assertEqual(O.cut_problem(cut, full, ["in", "seq"]), "ends: differs and is not named cut")
        self.assertEqual(O.cut_problem(dict(cut, segment=4), full, ["in", "seq", "ends"], unordered=["ends"]),
                         "segment: differs and is not named cut")
        self.assertIn("in: not a cut", O.cut_problem(dict(cut, **{"in": [2]}), full, ["in", "seq", "ends"],
                                                     unordered=["ends"]))
        self.assertEqual(O.cut_problem({"segment": 3}, full, ["in"]), "in: missing in the cut")
        self.assertTrue(O.is_cut_of("ACG", "ACGTT") and O.is_cut_of([1], [1, 2]) and not O.is_cut_of([2], [1, 2]))
        self.assertTrue(O.is_cut_of([2], [1, 2], unordered=[""]))
        # an empty cut is a cut: the library's row cut takes a string's tail with v[:max(0, len(v) - over - 1)] and
        # halves a list with v[:len(v) // 2], so a field can end as "" or [] (mcp_tools.py _shrink, :2713-2733,
        # called for a row larger than its page at :1068-1079); lane C follows this rule
        self.assertIsNone(O.cut_problem({"in": [], "seq": "", "segment": 3}, {"in": [7], "seq": "ACGT", "segment": 3},
                                        ["in", "seq"]))
        self.assertTrue(O.is_cut_of([], [{"at_bp": 1}], unordered=[""]) and O.is_cut_of("", "A"))
        self.assertTrue(O.is_cut_of([], [{"at_bp": 1}], position_ordered=[""]))
        self.assertEqual(O.cut_problem({"in": [], "segment": 3}, {"in": [7], "segment": 3}, []),
                         "in: differs and is not named cut")      # but only for a field named cut

    def test_paging_context(self):
        g, f, mgt = doc("caps")
        whole = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, f, mgt)
        two = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, f, mgt, {"paging": "none", "rows": 2})
        self.assertEqual(two["answer"]["rows"], whole["answer"]["rows"][:2])
        self.assertEqual(two["answer"]["total"], 3)
        blind = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, f, mgt, {"paging": "none"})
        self.assertEqual(blind["answer"]["rows"], [])
        self.assertIn("/rows/*", blind["oracle"]["consistency_only"])
        # a slice longer than a page: whole when every page is read, else consistency-only
        _, fk, mk = doc("fork")
        cut = O.answer("graphlet_sequence", {"arm": "left", "walk": 0, "with_seed": True}, fk, mk,
                       {"sequence_max_bytes": 64, "paging": "none"})
        self.assertNotIn("sequence", cut["answer"])
        self.assertIn("/next_from", cut["oracle"]["consistency_only"])
        allp = O.answer("graphlet_sequence", {"arm": "left", "walk": 0, "with_seed": True}, fk, mk,
                        {"sequence_max_bytes": 64})
        self.assertEqual(len(allp["answer"]["sequence"]), 50)
        self.assertNotIn("to", allp["answer"])

    def test_time_measured_values(self):
        """F4: a clock's values (a time limitation's observed, a time resource stop's amounts) are not compared
        with the twin (SPEC §6.8): coords/d3_kmer stopped at 1e-6 ms."""
        f = _load(os.path.join(DATA, "coords", "d3_kmer.full.json.gz"))
        g = _load(os.path.join(DATA, "coords", "d3_kmer.graphlet.json.gz"))
        mgt = g["results"][0]["graphlet"]
        r = O.answer("graphlet_summary", {"detail": "limitations"}, f, mgt)
        timed = [p for p in r["oracle"]["consistency_only"] if p.endswith("/observed")]
        self.assertEqual(len(timed), 3)
        j = O.answer("graphlet_export", {"format": "json"}, f, mgt)
        self.assertIn("/json/arms/left/limitations/0/observed", j["oracle"]["consistency_only"])
        self.assertIn("/json/arms/left/limitations/0/kind", j["oracle"]["checked"])
        other = copy.deepcopy(j["answer"])
        for side in ("left", "right"):
            other["json"]["arms"][side]["limitations"][0]["observed"] = 0.7
        self.assertEqual(O.differences(j, other), [])
        # a time resource stop (none among the twins): its scope, resource, phase, requested and actions are
        # compared field by field, its used, remaining, effective and message are not
        g, f, mgt = doc("fork")
        f["results"][0]["resource_stop"] = {"scope": "locus", "resource": "time", "phase": "traversal",
                                            "requested": 5.0, "effective": 5.0, "used": 5.3, "remaining": 0,
                                            "actions": ["raise_time_budget"], "message": "m"}
        r = O.answer("graphlet_summary", {"detail": "limitations"}, f, mgt)
        other = copy.deepcopy(r["answer"])
        other["resource_stop"] = dict(f["results"][0]["resource_stop"], used=9.9, remaining=1, message="x")
        self.assertEqual(O.differences(r, other), [])
        other["resource_stop"]["phase"] = "seed"
        self.assertEqual([d["pointer"] for d in O.differences(r, other)], ["/resource_stop/phase"])


class Normalisation(unittest.TestCase):
    def test_full_result_normalisation(self):
        """DESIGN §9 T37: events sorted within equal at_bp, needed_budgets sorted, timing removed; §2.5 growth
        label_ends null == {}; SPEC §8.3 annotation lookahead counters dropped."""
        r = {"timing": {"x": 1}, "annotation": {"access_path": "rows", "keys_mapped": 5, "direct_reads": 2,
                                                "rows_requested": 9},
             "arms": {"right": {"segments": [{"events": [{"at_bp": 5, "type": "switch"}, {"at_bp": 3, "type": "b"},
                                                         {"at_bp": 5, "type": "label_end"}]}],
                                "needed_budgets": [2.5, 0.5], "growth": [{"label_ends": None}, {"label_ends": {"R": 1}}],
                                "timing": {"y": 2}}}}
        n = O.normalise_full_result(r)
        self.assertEqual(n["annotation"], {"access_path": "rows", "rows_requested": 9})
        self.assertNotIn("timing", n)
        self.assertNotIn("timing", n["arms"]["right"])
        self.assertEqual([e["at_bp"] for e in n["arms"]["right"]["segments"][0]["events"]], [3, 5, 5])
        self.assertEqual(n["arms"]["right"]["needed_budgets"], [0.5, 2.5])
        self.assertEqual(n["arms"]["right"]["growth"][0]["label_ends"], {})
        self.assertEqual(r["arms"]["right"]["needed_budgets"], [2.5, 0.5])   # the input is left as it is
        # the engine's per-request counters (work.memory_bytes, the deploy of 2026-10-21): dropped
        r["work"] = {"memory_bytes": 2862, "steps": 75}
        n = O.normalise_full_result(r)
        self.assertNotIn("work", n)
        self.assertIn("work", r)

    def test_body_digest(self):
        """FORMAT-v8 §3: H seed_index masked, J and Z excluded."""
        _, _, body = doc("fork")
        lines = body.rstrip("\n").split("\n")
        h = lines[0].split(" ")
        h[11] = "3"
        moved = "\n".join([" ".join(h), 'J {"envelope": 1}'] + lines[1:-1] + ["Z %d" % (len(lines) + 1)]) + "\n"
        self.assertEqual(O.body_digest(body), O.body_digest(moved))
        self.assertNotEqual(O.body_digest(body), O.body_digest(body.replace("\nG ", "\nG  ", 1)))
        self.assertEqual(O.arms_digest(body), O.arms_digest(moved))

    def test_export_parsers(self):
        self.assertEqual(O.parse_fasta(">left_0 length_bp=4 seed=yes\nAC\nGT\n>right_1\n\nTT\n"),
                         [{"id": "left_0", "header": "left_0 length_bp=4 seed=yes", "sequence": "ACGT"},
                          {"id": "right_1", "header": "right_1", "sequence": "TT"}])
        # GFA 1.0 with k-1 overlaps on the L lines (k = 4), a reverse step and a P line with its own overlaps
        gfa = "\n".join(["H\tVN:Z:1.0", "S\ts\tACGTA", "S\ta\tGTAC", "S\tb\tTACG", "L\ts\t+\ta\t+\t3M",
                         "L\ta\t+\tb\t+\t3M", "P\tp0\ts+,a+,b+\t*", "S\tr\tCGTAC", "P\tp1\ts+,r-\t3M"])
        self.assertEqual(O.gfa_paths(gfa), ["ACGTACG", "ACGTACG"])
        self.assertEqual(O.export_view("gfa", {"records": 9}, gfa), {"records": 9, "gfa": {"n_paths": 2,
                                                                                         "paths": ["ACGTACG"] * 2}})


class CompareProperties(unittest.TestCase):
    """The compare properties (DESIGN-v8 §2) on engine documents, with compare answers written by hand: the body-only
    CLI pair of one locus at radii 60 and 100 (documents/compare, a body without twin) and the fork document."""

    @classmethod
    def setUpClass(cls):
        d = os.path.join(DATA, "documents", "compare")
        with open(os.path.join(d, "cli_bubbles_r60.mgt"), encoding="utf-8") as fh:
            cls.r60 = fh.read()
        with open(os.path.join(d, "cli_bubbles_r100.mgt"), encoding="utf-8") as fh:
            cls.r100 = fh.read()

    def test_self(self):
        g, f, mgt = doc("fork")      # no index_fp in H: a join is unverifiable (SPEC §7.5.7)
        ok = {"mode": "walks", "comparable": "unverifiable", "equal": None, "counts": None, "rows": []}
        self.assertEqual(O.compare_self_properties(ok, f, mgt), [])
        bad = dict(ok, comparable=True, equal=True)
        self.assertEqual(len(O.compare_self_properties(bad, f, mgt)), 1)
        self.assertEqual(len(O.compare_self_properties(dict(ok, rows=[{"side": "a"}]), f, mgt)), 1)
        g, f, mgt = doc("merge")
        h = mgt.split("\n", 1)
        stamped = h[0].rsplit(" ", 2)[0] + " " + "ab" * 32 + " " + h[0].rsplit(" ", 1)[1] + "\n" + h[1]
        self.assertEqual(O.compare_self_properties(dict(ok, comparable=True, equal=True,
                                                        counts={"only_in_a": 0, "only_in_b": 0, "differ": 0}),
                                                   f, stamped), [])
        self.assertEqual(len(O.compare_self_properties(dict(ok, comparable=False), f, stamped)), 1)

    def test_swap(self):
        ab = {"mode": "walks", "comparable": True, "equal": False, "depth_used": 60,
              "counts": {"only_in_a": 1, "only_in_b": 0, "differ": 1},
              "rows": [{"side": "a", "arm": "right", "walk": "AC", "labels": []},
                       {"side": "differ", "arm": "right", "walk": "GG", "a": ["x"], "b": ["y"]}]}
        ba = {"mode": "walks", "comparable": True, "equal": False, "depth_used": 60,
              "counts": {"only_in_a": 0, "only_in_b": 1, "differ": 1},
              "rows": [{"side": "differ", "arm": "right", "walk": "GG", "a": ["y"], "b": ["x"]},
                       {"side": "b", "arm": "right", "walk": "AC", "labels": []}]}
        self.assertEqual(O.compare_swap_properties(ab, ba), [])
        self.assertEqual(len(O.compare_swap_properties(ab, ab)), 2)

    def test_walks_and_prefix_subset_by_radius(self):
        """One locus at radii 60 and 100: within depth 60 one trie (SPEC §6.10), so walks mode has no difference
        and prefix_subset nothing to report; walks cut at a depth beyond 60 differ."""
        quiet = {"comparable": True, "equal": True, "depth_used": 60,
                 "counts": {"only_in_a": 0, "only_in_b": 0, "differ": 0}, "rows": []}
        self.assertEqual(O.compare_walks_properties(dict(quiet, mode="walks"), None, self.r60, None, self.r100), [])
        self.assertEqual(O.compare_prefix_subset_properties(dict(quiet, mode="prefix_subset"), None, self.r60, None,
                                                            self.r100), [])
        noisy = dict(quiet, mode="prefix_subset", rows=[{"side": "b", "arm": "right", "walk": "ACGT"}])
        self.assertEqual(len(O.compare_prefix_subset_properties(noisy, None, self.r60, None, self.r100)), 1)
        deep = dict(quiet, mode="prefix_subset", depth_used=80)
        self.assertTrue(O.compare_prefix_subset_properties(deep, None, self.r60, None, self.r100))
        # walks cut at 100 bp: the 60 bp walks are only in a, the 100 bp walks only in b (a set difference of the
        # spelled walks; a depth beyond complete_to_bp only exercises the set logic here)
        cut60, cut100 = (O._cut_walks(O._Doc(None, t), "right", 100) for t in (self.r60, self.r100))
        self.assertEqual((len(cut60), len(cut100), cut60 & cut100), (4, len(cut100), set()))
        rows = [{"side": "a", "arm": "right", "walk": w, "labels": []} for w in sorted(cut60)] + \
            [{"side": "b", "arm": "right", "walk": w, "labels": []} for w in sorted(cut100)]
        ans = dict(quiet, mode="walks", depth_used=100, equal=False, rows=rows)
        self.assertEqual(O.compare_walks_properties(ans, None, self.r60, None, self.r100), [])
        self.assertTrue(O.compare_walks_properties(dict(ans, rows=rows[1:]), None, self.r60, None, self.r100))


    def test_walks_property_on_one_arm(self):
        """VO2: the compare's arm comes from its arguments (the answer names none); rows of the left arm only are
        right for arm=left and short for both arms."""
        _, fa, ta = doc("fork")
        _, fb, tb = doc("linear")
        D = 10
        da, db = O._Doc(fa, ta), O._Doc(fb, tb)
        rows = []
        for side in ("left", "right"):
            wa, wb = O._cut_walks(da, side, D), O._cut_walks(db, side, D)
            rows.append([{"side": "a", "arm": side, "walk": w, "labels": []} for w in sorted(wa - wb)] +
                        [{"side": "b", "arm": side, "walk": w, "labels": []} for w in sorted(wb - wa)])
        self.assertTrue(rows[0] and rows[1])
        ans = {"mode": "walks", "comparable": "qualified", "equal": None, "depth_used": D, "rows": rows[0]}
        self.assertEqual(O.compare_walks_properties(ans, fa, ta, fb, tb, arm="left"), [])
        self.assertTrue(O.compare_walks_properties(ans, fa, ta, fb, tb))
        self.assertTrue(O.compare_walks_properties(dict(ans, rows=rows[0] + rows[1]), fa, ta, fb, tb, arm="left"))


class RunStep(unittest.TestCase):
    """FORMAT v2 F4: graphlet_oracle.py --run-dir RUN writes arms/oracle/graphlet_ops.jsonl from W1's records, the
    full twins (C2) and the engine's body (the direct-local raw answer; W1's document C1 otherwise)."""

    def make_run(self, root):
        import gzip as gz
        for sid, name in (("G-caps", "caps"), ("G-merge", "merge")):
            g, f, mgt = doc(name)
            os.makedirs(os.path.join(root, "arms", "direct-local", "raw"), exist_ok=True)
            os.makedirs(os.path.join(root, "arms", "direct-local", "docs"), exist_ok=True)
            os.makedirs(os.path.join(root, "arms", "dl-lib", "docs"), exist_ok=True)
            if sid == "G-caps":   # the engine's body from the raw answer
                with gz.open(os.path.join(root, "arms", "direct-local", "raw", "%s.direct-local.0.json.gz" % sid),
                             "wt", encoding="utf-8") as fh:
                    json.dump({"answer": {"body": g}}, fh)
            else:                 # W1's document only
                with gz.open(os.path.join(root, "arms", "dl-lib", "docs", "%s.mgt.gz" % sid), "wt",
                             encoding="utf-8") as fh:
                    fh.write(mgt)
            with gz.open(os.path.join(root, "arms", "direct-local", "docs", "%s.full.json.gz" % sid), "wt",
                         encoding="utf-8") as fh:
                json.dump(f, fh)
        # a failed seed (VR4): its full twin and its graphlet answer, neither with a graphlet, no W1 document
        g = _load(os.path.join(DATA, "real", "sra", "sra_batch3__batch_refuse.graphlet.json.gz"))
        f = _load(os.path.join(DATA, "real", "sra", "sra_batch3__batch_refuse.full.json.gz"))
        with gz.open(os.path.join(root, "arms", "direct-local", "raw", "G-fail.direct-local.0.json.gz"), "wt",
                     encoding="utf-8") as fh:
            json.dump({"answer": {"body": dict(g, results=[g["results"][1]])}}, fh)
        with gz.open(os.path.join(root, "arms", "direct-local", "docs", "G-fail.full.json.gz"), "wt",
                     encoding="utf-8") as fh:
            json.dump(dict(f, results=[f["results"][1]]), fh)
        caps_walks = O.answer("graphlet_walks", {"arm": "right", "rank": "length"}, *doc("caps")[1:])["answer"]
        w1 = [
            {"scenario": "G-caps", "op": "walks", "tool": "graphlet_walks", "paging_complete": True,
             "args_resolved": {"handle": "<handle>", "arm": "right", "rank": "length", "n": 10, "max_bytes": 2048},
             "answer": caps_walks, "library": {"limits": {"max_bytes": 2048, "sequence_max_bytes": 16384}}},
            {"scenario": "G-caps", "op": "walks-1", "tool": "graphlet_walks", "paging_complete": False,
             "args_resolved": {"handle": "<handle>", "arm": "right", "rank": "length", "n": 1},
             "answer": dict(caps_walks, rows=caps_walks["rows"][:1])},
            {"scenario": "G-caps", "op": "bad", "tool": "graphlet_walk", "selector_error": "no /rows/9/walk",
             "args_resolved": None},
            {"scenario": "G-caps", "op": "sub", "tool": "graphlet_subtrie",
             "args_resolved": {"handle": "<handle>", "labels": ["hX"]}},
            {"scenario": "G-caps", "op": "on-view", "tool": "graphlet_walks",
             "args_resolved": {"handle": "<handle:sub>", "arm": "right"}},
            {"scenario": "G-caps", "op": "load", "tool": "graphlet_load", "args_resolved": {"path": "x.mgt"}},
            {"scenario": "G-caps", "op": "on-load", "tool": "graphlet_summary",
             "args_resolved": {"handle": "<handle:load>", "max_bytes": 262144}},
            {"scenario": "G-merge", "op": "walk", "tool": "graphlet_walk", "paging_complete": True,
             "args_resolved": {"handle": "<handle>", "arm": "right", "walk": 0}},
            {"scenario": "G-merge", "op": "cmp", "tool": "graphlet_compare", "paging_complete": True,
             "args_resolved": {"a": "<handle>", "b": "<handle:G-caps>", "mode": "walks"}},
            {"scenario": "G-none", "op": "x", "tool": "graphlet_summary", "args_resolved": {"handle": "<handle>"}},
            # R2-5: symbols naming nothing, W1 without arguments
            {"scenario": "G-caps", "op": "q", "tool": "graphlet_walks", "args_resolved": {"handle": "<handle:?>",
                                                                                          "arm": "right"}},
            {"scenario": "G-caps", "op": "typo", "tool": "graphlet_walks",
             "args_resolved": {"handle": "<handle:G-cpas/sub>", "arm": "right"}},
            {"scenario": "G-caps", "op": "raw", "tool": "graphlet_walks",
             "args_resolved": {"handle": "g_0123456789ab", "arm": "right"}},
            {"scenario": "G-caps", "op": "noargs", "tool": "graphlet_walks", "args_resolved": None,
             "error": {"code": "no_document", "fields": {}}},
            {"scenario": "G-fail", "op": "fail-walks", "tool": "graphlet_walks", "paging_complete": None,
             "args_resolved": {"handle": "<handle>", "arm": "left"}, "answer": None,
             "error": {"code": "no_graphlet", "fields": {}}},
            # F-D2: a cursor from the one listing W1 did not read to its end (walks-1 above), with other arguments
            # (bad_cursor) and with its own (a resumed page: no claim)
            {"scenario": "G-caps", "op": "cur-other", "tool": "graphlet_walks",
             "args_resolved": {"handle": "<handle>", "arm": "right", "rank": "id", "n": 1, "cursor": "<cursor>"}},
            {"scenario": "G-caps", "op": "cur-same", "tool": "graphlet_walks",
             "args_resolved": {"handle": "<handle>", "arm": "right", "rank": "length", "n": 5, "cursor": "<cursor>"}},
            {"scenario": "G-caps", "op": "budget", "tool": "graphlet_walks", "library": {"limits": {"local_limits": None}},
             "args_resolved": {"handle": "<handle>", "arm": "right", "budget": {"work_units": 5}}}]
        with open(os.path.join(root, "arms", "dl-lib", "graphlet_ops.jsonl"), "w", encoding="utf-8") as fh:
            for r in w1:
                fh.write(json.dumps(dict(r, run_id="R", way="W1")) + "\n")
            # two torn lines (a kill mid-write): one still names its scenario and op, one names nothing
            fh.write('{"run_id": "R", "scenario": "G-merge", "way": "W1", "op": "torn", "tool": "graphlet_walks", '
                     '"answer": {"rows": [\n')
            fh.write('{"answer": {"rows": [{"walk": 0, "len\n')
        return w1

    def test_run_dir(self):
        import tempfile
        with tempfile.TemporaryDirectory() as root:
            w1 = self.make_run(root)
            path = O.write_run(root, quiet=True)
            with self.assertRaises(SystemExit):     # never appended
                O.write_run(root, quiet=True)
            with open(path, encoding="utf-8") as fh:
                lines = {r["op"]: r for r in map(json.loads, fh)}
            self.assertEqual(len(lines), len(w1) + 2)      # the two torn lines
            walks = lines["walks"]
            self.assertEqual((walks["way"], walks["args_resolved"], walks["library"]),
                             ("W5", w1[0]["args_resolved"], None))
            self.assertIn("the engine's body", walks["sources"]["mgt"])
            self.assertEqual(O.compare(walks, dict(w1[0], error=None))["outcome"], "agree")
            self.assertEqual(len(lines["walks-1"]["answer"]["rows"]), 1)       # W1 did not read to the end
            self.assertEqual(O.compare(lines["walks-1"], dict(w1[1], error=None))["outcome"], "agree")
            self.assertEqual(lines["bad"]["oracle"]["declared"], "selector_error")
            self.assertEqual(lines["on-view"]["oracle"]["declared"], "view_consistency_only")
            self.assertEqual(lines["on-load"]["answer"]["summary"]["arms"]["right"]["walks"], 3)  # the same graphlet
            self.assertIn("W1's document", lines["walk"]["sources"]["mgt"])
            self.assertEqual(lines["walk"]["answer"]["walk"], 0)
            self.assertEqual(lines["cmp"]["oracle"]["declared"], "properties_only")
            self.assertEqual(lines["x"]["oracle"]["declared"], "no_input")
            # R2-5: what W5 cannot judge is declared so, and compares inconclusive (it blocks)
            for op in ("q", "typo", "raw"):
                self.assertEqual(lines[op]["oracle"]["declared"], "unresolved_handle", op)
            self.assertEqual(lines["noargs"]["oracle"]["declared"], "no_input")
            o = lines["cur-other"]["oracle"]
            self.assertEqual((o["declared"], o["expected_error"]), ("expected_error", "bad_cursor"))
            self.assertEqual(lines["cur-same"]["oracle"]["declared"], "cursor_page")
            self.assertEqual(lines["budget"]["oracle"]["expected_error"], "bad_argument")
            fail = lines["fail-walks"]["oracle"]
            self.assertEqual((fail["declared"], fail["expected_error"]), ("expected_error", "no_graphlet"))
            self.assertEqual(O.compare(lines["fail-walks"], dict(w1[-1], error={"code": "no_graphlet", "fields": {}}))
                             ["outcome"], "agree")
            self.assertEqual((lines["torn"]["scenario"], lines["torn"]["tool"], lines["torn"]["oracle"]["declared"]),
                             ("G-merge", "graphlet_walks", "unreadable_w1"))
            run_level = lines[None]
            self.assertEqual((run_level["run_level"], run_level["scenario"], run_level["error"]["fields"]["lines"]),
                             (True, None, 1))
            for op in ("q", "typo", "raw", "noargs", "torn", "x", "bad", None):
                self.assertEqual(O.compare(lines[op], {"rows": []})["outcome"], "inconclusive", op)
            for op in ("on-view", "cmp"):
                self.assertEqual(O.compare(lines[op], {"rows": []})["outcome"], "not_checked", op)
            with open(os.path.join(root, "arms", "oracle", "run.json"), encoding="utf-8") as fh:
                self.assertEqual(json.load(fh)["unreadable_w1_lines"], 2)
            self.assertEqual(O.main(["--run-dir", root, "--force", "--quiet"]), 0)


class MemoryBudgetTwin(unittest.TestCase):
    """VERIFY-D VD4 (decision F-D4): v8 sends the full twin without bounds.max_memory_mb (SPEC §7.1 R21 (4): `full`
    detail under a memory budget keeps a merge's arrival order), so the graphlet answer states a memory_bound_soft
    (SPEC §7.0: stated exactly when the knob is set) the twin does not. W5 then reads the seed-level entries from the
    graphlet's own K records (SPEC §7.5.2 K; VERIFY-O-VD4 suggestion 1), checks the rest as before and names the
    difference in oracle.inputs. On the two real rehearsal-2 pairs."""

    @staticmethod
    def stated(g, f):
        """graphlet_summary(detail='limitations') as the graphlet states it: its seed-level entries (the response's
        summary), then each arm's (left first; the twin walked the same arms: twin pre-check), outcome and stop."""
        rg = g["results"][0]
        rows = [dict(e) for e in rg["limitations"]]
        for s in ("left", "right"):
            rows += [dict(e, arm=s) for e in (f["results"][0]["arms"].get(s) or {}).get("limitations") or []]
        return {"outcome": rg["outcome"], "resource_stop": rg.get("resource_stop"), "total": len(rows), "rows": rows}

    def test_the_pairs(self):
        for sid in REHEARSAL2:
            g, f, mgt = rehearsal2(sid)
            self.assertEqual(g["strategy"]["bounds"]["max_memory_mb"], 512)
            self.assertNotIn("max_memory_mb", f["strategy"]["bounds"])
            self.assertEqual([(e["kind"], e["knob"]) for e in g["results"][0]["limitations"]],
                             [("memory_bound_soft", "bounds.max_memory_mb")])
            self.assertEqual(f["results"][0]["limitations"], [])
            self.assertIsNone(O.twin_difference(f, mgt), sid)

    def test_k_records(self):
        """k_entries reads every K record as the JSON entry the graphlet's own summary states (SPEC §7.5.2 K), on
        every engine pair and both rehearsal-2 pairs: seed level, and each arm where the summary lists it."""
        pairs = engine_pairs() + [("r2/" + s, 0) + rehearsal2(s)[:2] for s in REHEARSAL2]
        n_arm = 0
        for cell, i, g, _ in pairs:
            rg = g["results"][i]
            ks = O.k_entries(rg["graphlet"])
            self.assertTrue(O.same([e for a, e in ks if a == "*"], rg.get("limitations") or []), cell)
            for a, side in (("l", "left"), ("r", "right")):
                summ = (rg.get("arms") or {}).get(side) or {}
                if "limitations" in summ:
                    n_arm += 1
                    self.assertTrue(O.same([e for x, e in ks if x == a], summ["limitations"]), (cell, side))
        self.assertGreater(n_arm, 60)
        # the typed values and the free text (SPEC §7.5.1)
        self.assertEqual(O.k_entries("K * x k.n u f:0.25 * a=s:p%25%2Cq%C3%A9,b=i:7 e%25f%0Ag"),
                         [("*", {"kind": "x", "knob": "k.n", "limit": "unlimited", "observed": 0.25, "a": "p%,qé",
                                 "b": 7, "effect": "e%f\ng"})])
        self.assertEqual(O.k_entries("K r y k.m i:3 * 12 . e"),
                         [("r", {"kind": "y", "knob": "k.m", "limit": 3, "complete_to_bp": 12, "effect": "e"})])
        for bad in ("K * x k i:01 * * . e", "K * x k q:1 * * . e", "K * x k i:1 * * a e", "K * x k i:1"):
            with self.assertRaises(ValueError):
                O.k_entries(bad)

    def test_limitations(self):
        for sid in REHEARSAL2:
            g, f, mgt = rehearsal2(sid)
            want = self.stated(g, f)
            self.assertGreater(len(want["rows"]), 2)
            r = O.answer("graphlet_summary", {"detail": "limitations"}, f, mgt, {"strategy": g["strategy"]})
            self.assertEqual(O.compare(r, want)["outcome"], "agree", sid)
            o = r["oracle"]
            # the seed-level block is the graphlet's K records, memory_bound_soft included, checked
            self.assertEqual(r["answer"]["rows"][0], g["results"][0]["limitations"][0])
            self.assertIn("SPEC §7.5.2 K (arm *)", o["checked"]["/rows/0/kind"])
            self.assertIn("/rows/0/observed", o["checked"])
            self.assertIn("/rows/1/kind", o["checked"])
            self.assertEqual(o["inputs"][:2], ["full", "mgt"])
            self.assertTrue(o["inputs"][2].startswith("twin without the memory budget: the graphlet request has "
                                                      "bounds.max_memory_mb 512"))
            self.assertIn("R21 (4)", o["inputs"][2])
            self.assertIn("§6.8", o["inputs"][2])
            other = copy.deepcopy(want)
            other["rows"][0]["observed"] = 9
            self.assertEqual([d["pointer"] for d in O.compare(r, other)["differences"]], ["/rows/0/observed"])
            # the rest is checked as before: each arm row in place, and the total, which counts the entry (§7.0)
            other = copy.deepcopy(want)
            other["rows"][1]["observed"] += 1
            self.assertEqual([d["pointer"] for d in O.compare(r, other)["differences"]], ["/rows/1/observed"])
            other = copy.deepcopy(want)
            del other["rows"][0]
            other["total"] -= 1
            self.assertEqual(O.compare(r, other)["outcome"], "disagree", sid)
            other = dict(want, outcome=dict(want["outcome"], walks="complete"))
            self.assertEqual(O.compare(r, other)["outcome"], "disagree", sid)
            # without the request's strategy W5 cannot tell (VD4 as measured: the rows shift by the entry)
            r0 = O.answer("graphlet_summary", {"detail": "limitations"}, f, mgt)
            self.assertEqual(r0["oracle"]["inputs"], ["full", "mgt"])
            self.assertEqual(O.compare(r0, want)["outcome"], "disagree", sid)
            # one arm's listing beside a seed-level entry is the tool's choice, as for a twin stating one
            ra = O.answer("graphlet_summary", {"detail": "limitations", "arm": "left"}, f, mgt,
                          {"strategy": g["strategy"]})
            self.assertIn("/rows/*", ra["oracle"]["consistency_only"])
            self.assertNotIn("total", ra["answer"])
            # a twin given as one result, without its envelope: the knob is read off its stated limitations (§7.0)
            r1 = O.answer("graphlet_summary", {"detail": "limitations"}, f["results"][0], mgt,
                          {"strategy": g["strategy"]})
            self.assertIn("(its result states no memory_bound_soft)", r1["oracle"]["inputs"][2])
            self.assertEqual(O.compare(r1, want)["outcome"], "agree", sid)
            self.assertIn("its strategy has no bounds.max_memory_mb and its result states no memory_bound_soft",
                          o["inputs"][2])
            # VERIFY-O-VD4 suggestion 2: a twin whose result states memory_bound_soft ran under the budget, whatever
            # its strategy echo says: no gap, the twin's entries as they are
            f2 = copy.deepcopy(f)
            f2["results"][0]["limitations"] = copy.deepcopy(g["results"][0]["limitations"])
            r2 = O.answer("graphlet_summary", {"detail": "limitations"}, f2, mgt, {"strategy": g["strategy"]})
            self.assertEqual(r2["oracle"]["inputs"], ["full", "mgt"])
            self.assertEqual(O.compare(r2, want)["outcome"], "agree", sid)
            # every record read from this twin names the difference
            w = O.answer("graphlet_walks", {"arm": "right"}, f, mgt, {"strategy": g["strategy"]})
            self.assertEqual(len(w["oracle"]["inputs"]), 3)

    def test_json_export(self):
        for sid in REHEARSAL2:
            g, f, mgt = rehearsal2(sid)
            j = O.answer("graphlet_export", {"format": "json"}, f, mgt, {"strategy": g["strategy"]})
            full = copy.deepcopy(f["results"][0])
            full["limitations"] = copy.deepcopy(g["results"][0]["limitations"])
            want = {"records": 1, "json": O.normalise_full_result(full)}
            self.assertEqual(O.compare(j, want)["outcome"], "agree", sid)
            self.assertIn("SPEC §7.5.2 K (arm *)", j["oracle"]["checked"]["/json/limitations"])
            self.assertIn("/json/arms/left/limitations", j["oracle"]["checked"])
            for change in (lambda x: x[0].update(observed=3), lambda x: x.append({"kind": "scope"}), list.clear):
                other = copy.deepcopy(want)
                change(other["json"]["limitations"])
                self.assertEqual([d["pointer"] for d in O.compare(j, other)["differences"]], ["/json/limitations"])
            other = copy.deepcopy(want)
            other["json"]["arms"]["left"]["limitations"] = []
            self.assertEqual(O.compare(j, other)["outcome"], "disagree", sid)
            j0 = O.answer("graphlet_export", {"format": "json"}, f, mgt)
            self.assertEqual(O.compare(j0, want)["outcome"], "disagree", sid)
            # the engine deploy's work block (work.memory_bytes) on the twin's result, none in the pinned library's
            # export: W5 answers from the twin, the export's view is read without it, they agree
            fw = copy.deepcopy(f)
            fw["results"][0]["work"] = {"memory_bytes": 2862, "steps": 75}
            jw = O.answer("graphlet_export", {"format": "json"}, fw, mgt, {"strategy": g["strategy"]})
            self.assertEqual(O.compare(jw, want)["outcome"], "agree", sid)
            self.assertNotIn("work", jw["answer"]["json"])

    def test_twin_with_the_budget(self):
        """documents/budgets/memory_soft: both executions under max_memory_mb 1 state memory_bound_soft; nothing
        changes: the entry is checked as the twin states it."""
        g = _load(os.path.join(DATA, "documents", "budgets", "memory_soft.graphlet.json"))
        f = _load(os.path.join(DATA, "documents", "budgets", "memory_soft.full.json"))
        mgt = g["results"][0]["graphlet"]
        for tool, args in (("graphlet_summary", {"detail": "limitations"}), ("graphlet_export", {"format": "json"})):
            r = O.answer(tool, args, f, mgt, {"strategy": g["strategy"]})
            self.assertEqual(r["oracle"]["inputs"], ["full", "mgt"])
            self.assertEqual(r, O.answer(tool, args, f, mgt))
        self.assertEqual(r["answer"]["json"]["limitations"][0]["kind"], "memory_bound_soft")
        self.assertIn("/json/limitations", r["oracle"]["checked"])

    def test_run_step(self):
        """The run step takes the request's strategy from the engine's response (the raw answer), else from W1's
        document's J record (SPEC §7.5.2: the response envelope)."""
        import tempfile
        g, f, mgt = rehearsal2("G-cut-c1")
        want = self.stated(g, f)
        rg = {k: v for k, v in g["results"][0].items() if k != "graphlet"}
        h, rest = mgt.split("\n", 1)
        with_j = h + "\nJ " + json.dumps(dict(g, results=[rg]), separators=(",", ":")) + "\n" + rest
        with tempfile.TemporaryDirectory() as root:
            for sub in ("raw", "docs"):
                os.makedirs(os.path.join(root, "arms", "direct-local", sub))
            os.makedirs(os.path.join(root, "arms", "dl-lib", "docs"))
            for sid in ("G-raw", "G-doc"):
                with gzip.open(os.path.join(root, "arms", "direct-local", "docs", "%s.full.json.gz" % sid), "wt",
                               encoding="utf-8") as fh:
                    json.dump(f, fh)
            with gzip.open(os.path.join(root, "arms", "direct-local", "raw", "G-raw.direct-local.0.json.gz"), "wt",
                           encoding="utf-8") as fh:
                json.dump({"answer": {"body": g}}, fh)
            with gzip.open(os.path.join(root, "arms", "dl-lib", "docs", "G-doc.mgt.gz"), "wt", encoding="utf-8") as fh:
                fh.write(with_j)
            w1 = [{"scenario": sid, "op": "lims", "tool": "graphlet_summary", "paging_complete": True,
                   "args_resolved": {"handle": "<handle>", "detail": "limitations", "max_bytes": 2048},
                   "answer": want, "error": None, "run_id": "R", "way": "W1"} for sid in ("G-raw", "G-doc")]
            with open(os.path.join(root, "arms", "dl-lib", "graphlet_ops.jsonl"), "w", encoding="utf-8") as fh:
                fh.write("".join(json.dumps(r) + "\n" for r in w1))
            with open(O.write_run(root, quiet=True), encoding="utf-8") as fh:
                lines = {r["scenario"]: r for r in map(json.loads, fh)}
        for sid, r in zip(("G-raw", "G-doc"), w1):
            self.assertEqual(O.compare(lines[sid], r)["outcome"], "agree", sid)
            self.assertTrue(lines[sid]["oracle"]["inputs"][2].startswith("twin without the memory budget"), sid)
        self.assertNotIn("strategy", lines["G-raw"]["sources"])
        self.assertIn(" J ", lines["G-doc"]["sources"]["strategy"])


class Independence(unittest.TestCase):
    """DESIGN-v8 §2 / REVIEW D-R5, mechanically: the oracle imports the standard library and the suite's checks and
    suite_common only (checks imports suite_common; both stdlib only), and names no library file or function."""

    STDLIB = {"__future__", "copy", "hashlib", "json", "typing", "gzip", "math", "re", "datetime", "os", "shlex",
              "sys", "zlib", "importlib", "importlib.util", "shutil", "time", "platform", "subprocess", "argparse",
              "posixpath"}

    def imports(self, path):
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        out = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                out |= {a.name for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                out.add(node.module)
        return out

    def test_imports(self):
        self.assertEqual(self.imports(os.path.join(SUITE, "graphlet_oracle.py")) - self.STDLIB, {"checks", "suite_common"})
        self.assertEqual(self.imports(os.path.join(SUITE, "checks.py")) - self.STDLIB, {"suite_common"})
        self.assertEqual(self.imports(os.path.join(SUITE, "suite_common.py")) - self.STDLIB, set())

    def test_no_library_citation(self):
        with open(os.path.join(SUITE, "graphlet_oracle.py"), encoding="utf-8") as fh:
            text = fh.read()
        for word in ("derive.py", "ops.py", "export.py", "mcp_tools", "model.py", "parser.py", "normalize_result",
                     "metagraph.traverse", "import metagraph", "from metagraph"):
            self.assertNotIn(word, text)


if __name__ == "__main__":
    unittest.main()
