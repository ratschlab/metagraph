#!/usr/bin/env python3
"""Build the suite's frozen scenario manifests (DESIGN §8.1): the full catalogue and the smoke selection.

    build_manifest.py --out manifest/suite_v4.json --smoke-out manifest/suite_smoke_v4.json [--seed N] [--no-fetch]
                      [--capabilities DIR]

Network: E-utilities for record sequences not in truth/cache yet (truth.py); nothing else. The same seed and the
same records give the same `content_sha256` (the sha256 of the canonical JSON of `scenarios`) and the same
`relations_sha256` (of `relations`, the comparisons across scenarios); `created_utc` is outside them.

The smoke manifest holds the scenarios of SMOKE_IDS (tagged `smoke-v4` in the full manifest): few enough for a
paired three-arm smoke run within 60 host requests, the MCP arm's included. The `smoke` tag keeps marking the
scenarios a smoke validation exercised (the earlier selections).

--capabilities DIR: the capabilities snapshot (capabilities.json, pattern_capabilities.json,
traverse_capabilities.json; default truth/staging_capabilities/): the C scenarios, the relations between the three
documents and the scenarios a host capability decides (the protein substitution cap: manifest_lib.HOST_MISMATCH_MAX)
are built from it, so the frozen manifest is built with a snapshot refreshed from the host it will run against.
The frozen v3 manifests stay as they are (the builder of v3 is commit 6f945db4).
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
import manifest_lib  # noqa: E402
import truth  # noqa: E402
from manifest_lib import DB, Ctx, canonical_json, sha256_json  # noqa: E402

DEFAULT_SEED = 20261013
MANIFEST_VERSION = 4
DEFAULT_CAPABILITIES = os.path.join(truth.TRUTH_DIR, "staging_capabilities")
# The v3 smoke selection (lane S3): the checks lane S3 changed, live, within 60 host requests in a paired three-arm
# run (each scenario MCP cold, direct-local, direct-remote, MCP warm; plus 9 capabilities GETs): the 503 past a
# budget (R-budget-0001), exact junction runs, a walk to a record's end, a 1 ms walk budget (a stop at depth 1),
# the telomeric seed's re-entry, an all_or_count hard region never withheld, a peptide, a count in suffix scope
# (properties the MCP may not expose), a 31-mer of a real repeat, and the max_mismatches 0 relation pair.
SMOKE_TAG = "smoke-v%d" % MANIFEST_VERSION
SMOKE_IDS = ("R-budget-0001", "R-junction-0001", "T-end-0001", "T-budget-0001", "N-hard-0022", "N-hard-0005",
             "X-peptide10-0009", "P-exact-0001", "N-realrep-0007", "M-refusal-0012", "M-twin-0001")
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


def load_capabilities(directory: str) -> dict:
    """The three capabilities documents of a snapshot directory."""
    caps = {}
    for name in ("capabilities", "pattern_capabilities", "traverse_capabilities"):
        with open(os.path.join(directory, name + ".json")) as f:
            caps[name] = json.load(f)
    return caps


def build(seed: int, fetch: bool, only_areas=None, capabilities_dir: str = DEFAULT_CAPABILITIES) -> dict:
    recs = truth.load_with_sequences(fetch=fetch)
    missing = [r.key for r in recs if r.fetch_error]
    ctx = Ctx(recs, seed)
    caps = load_capabilities(capabilities_dir)
    # the host's substitution caps per kind decide the M-protein scenarios and their by-design naming (G04, G05)
    manifest_lib.set_host_capabilities(caps["pattern_capabilities"])
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
                  "capabilities_snapshot": _snapshot_name(capabilities_dir),
                  "mismatches_max": dict(manifest_lib.HOST_MISMATCH_MAX),
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
        "not_applicable": not_applicable(caps),
        "counts": counts(scenarios, relations),
        "relations_sha256": sha256_json(relations),
        "relations": relations,
        "scenarios": scenarios,
    }
    return manifest


def _snapshot_name(directory: str) -> str:
    """The snapshot directory as the manifest names it: relative to this directory when inside it."""
    rel = os.path.relpath(os.path.abspath(directory), HERE)
    return (rel if not rel.startswith("..") else os.path.abspath(directory)) + "/"


def not_applicable(caps: dict) -> list:
    """DESIGN §3's areas that need a host property the snapshot does not state, left out with the reason (the
    integrator's decision 7: P-unmasked and M-unmasked only on an unmasked host)."""
    pc = caps["pattern_capabilities"]
    out = []
    if pc.get("counting") == "exact" and pc.get("mask") not in (None, "none"):
        out.append({"areas": ["P-unmasked", "M-unmasked"], "tag": "not_applicable_masked",
                    "reason": "the staging host is masked (capabilities counting %r, mask %r): counts are exact, "
                              "no upper_bound or estimate to test; build these on an unmasked host" % (
                                  pc.get("counting"), pc.get("mask"))})
    return out


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
        if sc.get("id") in SMOKE_IDS:
            tags.add(SMOKE_TAG)
        if not tags & {"smoke", SMOKE_TAG} and any(_exact_text(e) for e in sc["expect"]):
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
                       for t in ("last-record-verified", "code-derived", "machine-dependent", "smoke", SMOKE_TAG)}}


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
    sel = [s for s in manifest["scenarios"] if SMOKE_TAG in s["tags"] and "needs-next-deploy" not in s["tags"]]
    missing = sorted(set(SMOKE_IDS) - {s["id"] for s in sel})
    if missing and not manifest.get("partial_build"):
        raise SystemExit("smoke scenarios not in the manifest: %s" % ", ".join(missing))
    m = dict(manifest)
    m["scenarios"] = sel
    m["content_sha256"] = sha256_json(sel)
    ids = {s["id"] for s in sel}
    m["relations"] = [r for r in manifest.get("relations", []) if r["a"] in ids and r["b"] in ids]
    m["relations_sha256"] = sha256_json(m["relations"])
    m["counts"] = counts(sel, m["relations"])
    m["selection"] = "the full manifest's scenarios tagged %s (full content_sha256 %s)" % (
        SMOKE_TAG, manifest["content_sha256"])
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


# the committed manifests of earlier versions: never overwritten
FROZEN = ("suite_v1.json", "suite_smoke_v1.json", "suite_v2.json", "suite_smoke_v2.json", "suite_v3.json",
          "suite_smoke_v3.json")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "manifest", "suite_v%d.json" % MANIFEST_VERSION))
    ap.add_argument("--smoke-out", default=os.path.join(HERE, "manifest", "suite_smoke_v%d.json" % MANIFEST_VERSION))
    ap.add_argument("--capabilities", default=DEFAULT_CAPABILITIES,
                    help="the capabilities snapshot directory (default truth/staging_capabilities/)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--no-fetch", action="store_true", help="use the record cache only")
    ap.add_argument("--areas", help="comma-separated areas (development only; the frozen manifest has all)")
    args = ap.parse_args(argv)
    only = set(args.areas.split(",")) if args.areas else None
    for path in (args.out, args.smoke_out):
        frozen = os.path.basename(path) in FROZEN and os.path.abspath(os.path.dirname(path)) == \
            os.path.join(HERE, "manifest")
        if frozen and os.path.exists(path):
            raise SystemExit("%s is a frozen manifest: build into another file" % path)
    m = build(args.seed, fetch=not args.no_fetch, only_areas=only, capabilities_dir=args.capabilities)
    if only:
        m["partial_build"] = True
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
