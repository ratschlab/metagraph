#!/usr/bin/env python3
"""The unit mapping of a manifest v8 (doc/units/suite_v8.json) from the committed v7 mapping (doc/units/suite_v7.json):
manifest v8 is v7's 725 scenarios byte for byte, in v7's order, plus ten warm-ups of T-walk-0024 (R-warm-0009 ..
-0018, chapter R's warm-up unit), area G (the graphlet corpus, its compare partners and the long walk's warm-ups:
chapter G) and the stop scenarios with their base requests (chapter S). Every v7 entry, unit and family is kept as it
is, R-warm's unit gains its ten members; the mapping is keyed to the given v8 manifest's content sha (`_manifest`), so
the real build (build_real_v8.sh) makes it from the manifest it is given, whatever the freeze changed.

The premises, checked (the script refuses otherwise): the v7 mapping is keyed to the v7 manifest; v8 holds v7's
scenarios unchanged and in v7's order and v7's relations; every other v8 scenario is an R-warm warm-up, an area G or an
area S scenario. Pure stdlib; writes OUT only.

    python3 -I units_v8_build.py MANIFEST_V8 OUT.json [MANIFEST_V7 MAPPING_V7]

MANIFEST_V7 defaults to manifest/suite_v7.json, MAPPING_V7 to doc/units/suite_v7.json (this tree's).
"""
import os
import json
import re
import sys
from collections import OrderedDict

G_FAMILIES = OrderedDict([
    ("G-corpus", ("the corpus graphlets", "list", "the corpus graphlets, fetched and processed five ways",
                  "A graphlet of the corpus (DESIGN-v8 §4): fetched on each access path (the direct requests with "
                  "`detail: graphlet` and the `full` twin, `traverse_fetch` via MCP), then processed by every "
                  "graphlet tool its op plan names on the ways that apply; the comparison sections of this chapter "
                  "give the outcome per tool.")),
    ("G-pair", ("compare partners", "list", "the compare partners",
                "The partner graphlet of a `graphlet_compare` pair: fetched and processed as the corpus graphlets, "
                "its own op plan small; the pair's compare calls are counted under its corpus graphlet.")),
    ("G-warm", ("warm-ups of the long walk", "replicates", "warm-ups of the long walk",
                "A warm-up: a window of the long walk's region resolved on direct-local right before it, so that "
                "the walk meets rows the host has just read.")),
])
S_FAMILIES = OrderedDict([
    ("S-stop", ("stops of running requests", "list", "traverse_stop on a running request",
                "`traverse_stop` sent via MCP to a running walk or discovery: the stop's answer, the job's end, the "
                "host attempt's state and the release of its token, read from the service's release record.")),
    ("S-base", ("the stops' base requests", "list", "the stops' base requests run without a stop",
                "The request a stop interrupts, run without a stop twice on direct-local: the warm time (repeat 1) "
                "the stop must beat.")),
])


def g_family(sid):
    if sid.startswith("G-pair-"):
        return "G-pair"
    if sid.startswith("G-warm-"):
        return "G-warm"
    return "G-corpus"


def variant_of(s):
    t = str(s.get("title") or s["id"])
    t = re.sub(r"\Aarea G (pair )?[a-z_0-9 ]+: ", "", t)
    t = re.sub(r"\Awarm-up of (G-[a-z]+-\d{4}): ", r"of \1: ", t)
    if s["id"].startswith("G-warm-"):
        m = re.match(r"of (G-[a-z]+-\d{4}): ([\d,]+ bp)", t)
        t = "of %s (%s)" % m.groups() if m else t
    if s["id"].startswith("S-stop-"):
        t = "traverse_stop on a running %s" % ("walk" if "POST /traverse" == s.get("route") else "discovery")
    if s["id"].startswith("S-base-"):
        m = re.match(r"(S-stop-\d{4})'s base request", t)
        t = "of %s (%s)" % (m.group(1) if m else s["id"], "a walk" if s.get("route") == "POST /traverse" else
                            "a discovery")
    return t


def main(v7_path, v8_path, mapping_path, out_path):
    v7 = json.load(open(v7_path), object_pairs_hook=OrderedDict)
    v8 = json.load(open(v8_path), object_pairs_hook=OrderedDict)
    mapping = json.load(open(mapping_path), object_pairs_hook=OrderedDict)
    assert v7["manifest_version"] == 7 and v8["manifest_version"] == 8
    assert mapping["_manifest"]["content_sha256"] == v7["content_sha256"], "the mapping is not v7's"
    s7 = OrderedDict((s["id"], s) for s in v7["scenarios"])
    s8 = OrderedDict((s["id"], s) for s in v8["scenarios"])
    assert [x for x in s8 if x in s7] == list(s7), "v8 does not hold v7's ids in v7's order"
    assert all(s7[x] == s8[x] for x in s7), "v8 changes a v7 scenario"
    assert v7.get("relations") == v8.get("relations"), "the relations differ"
    new = [x for x in s8 if x not in s7]
    for x in new:
        a = s8[x].get("area")
        assert (x.startswith("R-warm-") and a == "R") or a in ("G", "S"), "an unexpected new scenario: %s" % x
    out = OrderedDict()
    out["_about"] = ("Unit mapping of manifest v8 (content_sha256 %s, %d scenarios) for the suite results document: "
                     "every scenario -> chapter, family, unit, unit_type, variant and optional group; _units carries "
                     "each unit's title, type and property line; _families each family's chapter and title. Built by "
                     "units_v8_build.py from the v7 mapping: v7's entries, units and families as they are; R-warm's "
                     "unit gains the ten warm-ups of T-walk-0024; chapter G (the graphlet tools: the corpus, its "
                     "compare partners, the long walk's warm-ups) and chapter S (the stops and their base requests) "
                     "are new." % (v8["content_sha256"], len(s8)))
    out["_manifest"] = OrderedDict([("version", 8), ("content_sha256", v8["content_sha256"]), ("scenarios", len(s8))])
    chapters = OrderedDict()
    for k, v in mapping["_chapters"].items():
        chapters[k] = v
        if k == "T":
            chapters["G"] = "Graphlet tools"
            chapters["S"] = "Stops"
    out["_chapters"] = chapters
    fams = OrderedDict(mapping["_families"])
    for fid, (title, _, _, _) in list(G_FAMILIES.items()):
        fams[fid] = OrderedDict([("chapter", "G"), ("title", title), ("units", [fid])])
    for fid, (title, _, _, _) in list(S_FAMILIES.items()):
        fams[fid] = OrderedDict([("chapter", "S"), ("title", title), ("units", [fid])])
    out["_families"] = fams
    units = OrderedDict((k, OrderedDict(v)) for k, v in mapping["_units"].items())
    for fid, (_, utype, utitle, prop) in list(G_FAMILIES.items()) + list(S_FAMILIES.items()):
        units[fid] = OrderedDict([("family", fid), ("unit_type", utype), ("unit_title", utitle), ("property", prop),
                                  ("n", 0), ("members", [])])
    entries = OrderedDict()
    for x in s8:
        if x in s7:
            entries[x] = mapping[x]
            continue
        s = s8[x]
        if x.startswith("R-warm-"):
            base = mapping["R-warm-0001"]
            e = OrderedDict(base)
            m = re.search(r"(\d[\d,]* bp)", s.get("title") or "")
            e["variant"] = "of %s (%s)" % (s.get("warms") or "?", m.group(1) if m else "?")
            entries[x] = e
            units["R-warm"]["members"].append(x)
            units["R-warm"]["n"] = len(units["R-warm"]["members"])
            continue
        fid = g_family(x) if s["area"] == "G" else ("S-stop" if x.startswith("S-stop-") else "S-base")
        u = units[fid]
        e = OrderedDict([("chapter", s["area"]), ("family", fid), ("unit_id", fid), ("unit_type", u["unit_type"]),
                         ("variant", variant_of(s)), ("unit_title", u["unit_title"])])
        entries[x] = e
        u["members"].append(x)
        u["n"] = len(u["members"])
    out["_units"] = units
    for k, v in entries.items():
        out[k] = v
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print("wrote %s: %d scenarios, %d new" % (out_path, len(entries), len(new)))


HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)


if __name__ == "__main__":
    if len(sys.argv) not in (3, 5):
        sys.stderr.write(__doc__)
        sys.exit(2)
    v7_manifest = sys.argv[3] if len(sys.argv) == 5 else os.path.join(SUITE, "manifest", "suite_v7.json")
    v7_mapping = sys.argv[4] if len(sys.argv) == 5 else os.path.join(HERE, "units", "suite_v7.json")
    main(v7_manifest, sys.argv[1], v7_mapping, sys.argv[2])
