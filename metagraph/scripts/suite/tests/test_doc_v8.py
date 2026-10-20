"""The v8 results document as a delta on v7 (2026-10-20): manifest v8 is v7's 725 scenarios unchanged, ten warm-ups of
T-walk-0024, area G (the graphlet tools compared five ways) and the stop scenarios; the unit mapping is v7's with
chapter R's warm-up unit extended and chapters G and S added, the prose's v8 statements are `manifest_version >= 8`
branches (a v5 / v6 / v7 build stays as it was). Offline; stdlib only."""

import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
DOC = os.path.join(SUITE, "doc")
for p in (SUITE, HERE, DOC):
    if p not in sys.path:
        sys.path.insert(0, p)

import suite_common as sc  # noqa: E402

import units as units_mod  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")


def prose(name):
    return sc.read_json(os.path.join(PROSE_DIR, "%s.json" % name))


class _Stub(object):
    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def v8_manifest_stub():
    """A manifest shaped like v8 for the mapping: v7's scenarios as they are plus, for every scenario the v8 mapping
    adds, a stub of its area (the real v8 manifest is frozen by the owner, so the test does not depend on it)."""
    v7 = sc.read_json(os.path.join(SUITE, "manifest", "suite_v7.json"))
    raw8 = sc.read_json(os.path.join(DOC, "units", "suite_v8.json"))
    ids7 = [s["id"] for s in v7["scenarios"]]
    scen = dict((s["id"], s) for s in v7["scenarios"])
    order = [k for k in raw8 if not k.startswith("_")]
    out = []
    for sid in order:
        out.append(scen.get(sid) or {"id": sid, "area": raw8[sid]["chapter"], "title": sid})
    return {"manifest_version": 8, "scenarios": out, "relations": v7.get("relations")}, ids7, raw8


class V8MappingTest(unittest.TestCase):
    def test_mapping_is_v7_plus_the_new_scenarios(self):
        m8, ids7, raw8 = v8_manifest_stub()
        raw7 = sc.read_json(os.path.join(DOC, "units", "suite_v7.json"))
        self.assertEqual(raw8["_manifest"]["version"], 8)
        self.assertEqual(units_mod.default_mapping_path({"manifest_version": 8, "content_sha256": "not frozen"}),
                         os.path.join(DOC, "units", "suite_v8.json"))
        # v7's entries, units and families as they are; R-warm gains the ten warm-ups of T-walk-0024
        for sid in ids7:
            self.assertEqual(raw8[sid], raw7[sid], sid)
        for uid, u in raw7["_units"].items():
            if uid == "R-warm":
                continue
            self.assertEqual(raw8["_units"][uid], u, uid)
        added = ["R-warm-%04d" % i for i in range(9, 19)]
        self.assertEqual(raw8["_units"]["R-warm"]["members"], raw7["_units"]["R-warm"]["members"] + added)
        self.assertEqual(raw8["_units"]["R-warm"]["n"], 18)
        for sid in added:
            self.assertRegex(raw8[sid]["variant"], r"\Aof T-walk-0024 \(\d+ bp\)\Z")
        for fid, f in raw7["_families"].items():
            self.assertEqual(raw8["_families"][fid], f, fid)
        # chapters G and S after T, the others in v7's order
        self.assertEqual(list(raw8["_chapters"]), ["R", "T", "G", "S"] + list(raw7["_chapters"])[2:])
        self.assertEqual(sorted(f for f, v in raw8["_families"].items() if v["chapter"] in ("G", "S")),
                         ["G-corpus", "G-pair", "G-warm", "S-base", "S-stop"])
        units = units_mod.resolve(m8, units_mod.load_mapping(os.path.join(DOC, "units", "suite_v8.json")))
        self.assertEqual(units_mod.validate(units, m8), [])
        self.assertEqual(len(units.chapter_scenarios("G")), 19)
        self.assertEqual(len(units.chapter_scenarios("S")), 4)


if __name__ == "__main__":
    unittest.main()
