#!/usr/bin/env python3
"""Build the suite's frozen scenario manifests (DESIGN §8.1): the full catalogue and the smoke selection.

    build_manifest.py --out manifest/suite_v2.json --smoke-out manifest/suite_smoke_v2.json [--seed N] [--no-fetch]

Network: E-utilities for record sequences not in truth/cache yet (truth.py); nothing else. The same seed and the
same records give the same `content_sha256` (the sha256 of the canonical JSON of `scenarios`) and the same
`relations_sha256` (of `relations`, the comparisons across scenarios); `created_utc` is outside them.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import areas  # noqa: E402
import checks  # noqa: E402
import truth  # noqa: E402
from manifest_lib import DB, Ctx, canonical_json, sha256_json  # noqa: E402

DEFAULT_SEED = 20261013
MANIFEST_VERSION = 2
CONVENTIONS = {
    "truth.start": "0-based position of the window's first base on the record's forward strand (whole-record coordinates)",
    "truth.window": "the record's bases on its forward strand (for a '-' source the reverse complement of the pattern)",
    "truth.flank": "record bases around the window, starting at truth.flank_start (0-based); clipped at the record's ends "
                   "and at the fetched region's ends",
    "nt_coords": "host placements are 1-based inclusive (SPEC-pattern-search.md §14.3)",
    "resolve runs": "0-based half-open k-mer intervals (SPEC-labeled-traversal-core.md §4.2)",
    "body_raw / headers": "optional scenario fields (lane S2's runner): send body_raw verbatim instead of the JSON "
                          "of request; headers override the runner's defaults (null removes one); the method is "
                          "the route's",
    "mcp.by_design": "case ids of the async session's service-by-design table (service_by_design.json)",
    "mcp.requires_service": "the service increment the MCP mapping needs (I12: max_mismatches)",
    "timing_class": "fast (< ~5 s), slow, long (> 100 s: direct-local only)",
}


def build(seed: int, fetch: bool, only_areas=None) -> dict:
    recs = truth.load_with_sequences(fetch=fetch)
    missing = [r.key for r in recs if r.fetch_error]
    ctx = Ctx(recs, seed)
    caps = {}
    for name in ("capabilities", "pattern_capabilities", "traverse_capabilities"):
        with open(os.path.join(truth.TRUTH_DIR, "staging_capabilities", name + ".json")) as f:
            caps[name] = json.load(f)
    gens = [("R", areas.gen_R), ("T", areas.gen_T), ("P", areas.gen_P), ("L", areas.gen_L), ("X", areas.gen_X),
            ("M", areas.gen_M), ("Q", areas.gen_Q), ("N", areas.gen_N)]
    scenarios = []
    for area, gen in gens:
        if only_areas and area not in only_areas:
            continue
        scenarios += gen(ctx)
    if not only_areas or "C" in only_areas:
        scenarios += areas.gen_C(ctx, caps)
    if not only_areas:
        scenarios += areas.gen_REL(ctx, scenarios, caps)   # the relations' twins (needs every area)
    if not only_areas or "D" in only_areas:
        scenarios += areas.gen_D(ctx, scenarios)
    tag_scenarios(scenarios, recs)
    relations = ctx.relations
    problems = validate(scenarios) + validate_relations(relations, scenarios)
    if problems:
        raise SystemExit("manifest invalid:\n  " + "\n  ".join(problems[:30]))
    tc = caps["traverse_capabilities"]
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "content_sha256": sha256_json(scenarios),
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seed": seed,
        "index": {"database": DB, "release_at_build": caps["capabilities"].get("release"),
                  "index_fp": tc.get("index_fp"), "index_ns": tc.get("index_ns"), "k": tc.get("k"),
                  "capabilities_snapshot": "truth/staging_capabilities/ (2026-10-13)",
                  "kmers_for_the_chance_model": areas.REFSEQ33M_KMERS},
        "ground_truth": {
            "genomes_source": truth.SOURCES["genomes"],
            "proteins_source": truth.SOURCES["manifest_v1"],
            "index_records_source": truth.SOURCES["index_records"],
            "last_record_scan_source": truth.SOURCES["last_record_scan"],
            "records": [r.describe() for r in recs],
            "records_not_fetched": missing,
        },
        "conventions": CONVENTIONS,
        "skipped": ctx.skipped,
        "counts": counts(scenarios, relations),
        "relations_sha256": sha256_json(relations),
        "relations": relations,
        "scenarios": scenarios,
    }
    return manifest


def _records_of(t: dict) -> list:
    out = []
    if isinstance(t, dict):
        if "record" in t and "genome" in t:
            out.append(t)
        for key in ("a", "b", "shuffled_from"):
            out += _records_of(t.get(key))
    return out


def _exact_text(e: dict) -> bool:
    a = e.get("args") or {}
    return (e["check"] in ("status", "pattern_slot_error") and any(x in a for x in ("text", "text_prefix", "text_contains")))


def tag_scenarios(scenarios, recs) -> None:
    """Tags read from the scenario and the index facts, after every area is built:
    - `last-record-verified`: a truth record that is the last of its taxid column (the 2026-10-09 defect); the
      truth carries last_in_column and last_record_scan (the scan's status, or "single-record column": 0 of 500
      single-record NC_ columns were found missing), so that the report can list a failure there apart;
    - `code-derived`: an exact error text (status, slot) that no smoke run exercised: read from the source and
      the fixtures (smoke scenarios are exercised by the smoke validation)."""
    by_acc = {r.accession: r for r in recs}
    for sc in scenarios:
        truths = []
        for t in [sc.get("truth")] + list(sc.get("truths") or []):
            truths += _records_of(t)
        last = False
        for t in truths:
            r = by_acc.get(t["record"])
            if r is None or not r.last_in_column:
                continue
            t["last_in_column"] = True
            t["last_record_scan"] = r.last_record_scan or ("single-record column" if r.records_in_column == 1 else None)
            last = True
        tags = set(sc["tags"])
        if last and "index-defect" not in tags:
            tags.add("last-record-verified")
        if "smoke" not in tags and any(_exact_text(e) for e in sc["expect"]):
            tags.add("code-derived")
        sc["tags"] = sorted(tags)


def counts(scenarios, relations=()) -> dict:
    by_area = {}
    for s in scenarios:
        by_area[s["area"]] = by_area.get(s["area"], 0) + 1
    by_type = {}
    for r in relations:
        by_type[r["type"]] = by_type.get(r["type"], 0) + 1
    return {"scenarios": len(scenarios), "executions_per_arm": sum(s["repeat"] for s in scenarios),
            "by_area": dict(sorted(by_area.items())),
            "mcp_not_implemented": sum(1 for s in scenarios if "not_implemented" in s["mcp"]),
            "relations": len(relations), "relations_by_type": dict(sorted(by_type.items())),
            "by_tag": {t: sum(1 for s in scenarios if t in s["tags"])
                       for t in ("last-record-verified", "code-derived", "machine-dependent", "smoke")}}


def validate_relations(relations, scenarios) -> list:
    """Unique ids, a known type, both ends scenarios of the manifest and distinct, the args each type reads."""
    problems = []
    ids = {s["id"]: s for s in scenarios}
    seen = set()
    need = {"prefix_equals": ("labels", "a_from_kmer", "b_from_kmer"), "zero_equals_absent": ("pattern", "count"),
            "subset": ("a_pattern", "b_pattern", "key"), "same_answer": ()}
    for r in relations:
        if r["id"] in seen:
            problems.append("duplicate relation id %s" % r["id"])
        seen.add(r["id"])
        if r["type"] not in need:
            problems.append("%s: type %s" % (r["id"], r["type"]))
            continue
        for end in ("a", "b"):
            if r[end] not in ids:
                problems.append("%s: %s %s is not a scenario" % (r["id"], end, r[end]))
        if r["a"] == r["b"]:
            problems.append("%s: a and b are one scenario" % r["id"])
        for f in need[r["type"]]:
            if f not in r["args"]:
                problems.append("%s: no args.%s" % (r["id"], f))
        if r["type"] == "zero_equals_absent" and r["a"] in ids and r["b"] in ids:
            a, b = ids[r["a"]], ids[r["b"]]
            if a["route"] != "POST /pattern" or b["route"] != "POST /resolve" or a["request"].get("strands") != "forward":
                problems.append("%s: needs a forward /pattern a and a /resolve b" % r["id"])
        if r["type"] == "prefix_equals" and r["a"] in ids and r["b"] in ids:
            a, b = ids[r["a"]], ids[r["b"]]
            sa, sb = a["request"]["sequence"], b["request"]["sequence"]
            fa, fb = r["args"]["a_from_kmer"], r["args"]["b_from_kmer"]
            # both queries are pieces of one sequence S: S[fa:] starts with a's, S[fb:] with b's
            lo = max(fa, fb)
            if sa[lo - fa:][:100] != sb[lo - fb:][:100]:
                problems.append("%s: the two queries are not pieces of one sequence" % r["id"])
    return problems


def validate(scenarios) -> list:
    """Unique ids, every check known to checks.py, the contract's fields present."""
    problems = []
    seen = set()
    for s in scenarios:
        if s["id"] in seen:
            problems.append("duplicate id %s" % s["id"])
        seen.add(s["id"])
        for f in ("id", "area", "title", "route", "request", "mcp", "truth", "expect", "timing_class", "repeat", "tags"):
            if f not in s:
                problems.append("%s: no %s" % (s.get("id"), f))
        if s["timing_class"] not in ("fast", "slow", "long"):
            problems.append("%s: timing_class %s" % (s["id"], s["timing_class"]))
        if not ("tool" in s["mcp"] or "not_implemented" in s["mcp"]):
            problems.append("%s: mcp without tool or not_implemented" % s["id"])
        if s["mcp"].get("tool") and "idempotency_key" in s["mcp"].get("arguments", {}):
            problems.append("%s: an idempotency_key in mcp.arguments" % s["id"])
        if not any(e["check"] == "status" for e in s["expect"]):
            problems.append("%s: no status expectation" % s["id"])
        for e in s["expect"]:
            problems += ["%s: %s" % (s["id"], m) for m in checks.validate_expectation(e)]
        req = s["request"]
        if s["route"] == "POST /pattern" and isinstance(req, dict) and "patterns" in req:
            for f in ("mode", "max_contexts", "time_budget_ms"):
                if f not in req:
                    problems.append("%s: /pattern request without %s" % (s["id"], f))
    return problems


def smoke_of(manifest: dict) -> dict:
    sel = [s for s in manifest["scenarios"] if "smoke" in s["tags"] and "needs-next-deploy" not in s["tags"]]
    m = dict(manifest)
    m["scenarios"] = sel
    m["content_sha256"] = sha256_json(sel)
    ids = {s["id"] for s in sel}
    m["relations"] = [r for r in manifest.get("relations", []) if r["a"] in ids and r["b"] in ids]
    m["relations_sha256"] = sha256_json(m["relations"])
    m["counts"] = counts(sel, m["relations"])
    m["selection"] = "the full manifest's scenarios tagged smoke (full content_sha256 %s)" % manifest["content_sha256"]
    return m


def write(path: str, manifest: dict) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        # one scenario per line keeps the frozen file reviewable in a diff
        head = {k: v for k, v in manifest.items() if k not in ("scenarios", "relations")}
        f.write("{\n")
        for k, v in head.items():
            f.write("%s: %s,\n" % (json.dumps(k), canonical_json(v)))
        rels = manifest.get("relations", [])
        f.write('"relations": [\n')
        for i, r in enumerate(rels):
            f.write(canonical_json(r) + (",\n" if i + 1 < len(rels) else "\n"))
        f.write("],\n")
        f.write('"scenarios": [\n')
        for i, s in enumerate(manifest["scenarios"]):
            f.write(canonical_json(s) + (",\n" if i + 1 < len(manifest["scenarios"]) else "\n"))
        f.write("]\n}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "manifest", "suite_v2.json"))
    ap.add_argument("--smoke-out", default=os.path.join(HERE, "manifest", "suite_smoke_v2.json"))
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--no-fetch", action="store_true", help="use the record cache only")
    ap.add_argument("--areas", help="comma-separated areas (development only; the frozen manifest has all)")
    args = ap.parse_args(argv)
    only = set(args.areas.split(",")) if args.areas else None
    m = build(args.seed, fetch=not args.no_fetch, only_areas=only)
    write(args.out, m)
    s = smoke_of(m)
    write(args.smoke_out, s)
    for sk in m["skipped"]:
        print("skipped: %(area)s %(scenario)s: %(reason)s" % sk, file=sys.stderr)
    print(json.dumps({"full": {"path": args.out, "content_sha256": m["content_sha256"], **m["counts"]},
                      "smoke": {"path": args.smoke_out, "content_sha256": s["content_sha256"], **s["counts"]}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
