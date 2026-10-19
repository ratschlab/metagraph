#!/usr/bin/env python3
"""Build the suite's frozen scenario manifests (DESIGN §8.1): the full catalogue and the smoke selection.

    build_manifest.py --out manifest/suite_v7.json --smoke-out manifest/suite_smoke_v7.json [--seed N] [--no-fetch]
                      [--capabilities DIR] [--manifest-version 6|7]

Network: E-utilities for record sequences not in truth/cache yet (truth.py); nothing else. The same seed and the
same records give the same `content_sha256` (the sha256 of the canonical JSON of `scenarios`) and the same
`relations_sha256` (of `relations`, the comparisons across scenarios); `created_utc` is outside them.

The smoke manifest holds the scenarios of SMOKE_IDS (tagged `smoke-v5` in the full manifest: the tag names the
smoke SELECTION, v5's, which v6 keeps, so that the scenarios of v6 the deploy gate selects stay byte-identical to
v5's): few enough for a paired three-arm smoke run within 60 host requests, the MCP arm's included. The `smoke` tag
keeps marking the scenarios a smoke validation exercised (the earlier selections).

v6 (the owner, 2026-10-15 / -16) is v5 with every /pattern request that named the default budget of 60,000 ms
naming 110,000 ms (manifest_lib.PATTERN_TIME_BUDGET_MS), on both arms, and, for gating POST /pattern (round 2): the
/pattern checks the MCP's answer cannot carry marked not implemented yet via MCP (areas.mark_pattern_mcp_fields), the
six heavy /pattern scenarios tagged budget_stress (areas.tag_pattern_budget_stress) and a warm-up right before each
/pattern scenario whose cold MCP call timed out in the v5 baseline (areas.pattern_warm_ups; the deploy gate's
executions kept by keep_gate_warm_repeats); no relation changed, and the gate's 237 resolve / traverse /
capabilities scenarios are byte-identical with their warm repeats, so the v5 acceptance carries over. The MCP
exposure baseline it reads is v4's with its request digests rebased (--exposure-rebase, below).

v7 (the owner, 2026-10-18) is v6 with exactly these changes: P-edge-0004's entry_fields check, which reads the
pattern text the MCP's /pattern answer leaves out by design (the service owner's statement, 2026-10-17), is classed
service-side by design at property level (mcp.by_design_checks under the table's pattern.lowercase,
areas.MCP_PATTERN_ENTRY_BY_DESIGN) instead of not implemented yet via MCP; v6 is built again, byte for byte, with
--manifest-version 6 (build(manifest_version=6)); and every /pattern scenario's budget_respected states the host's
finalisation reserve, reserve_ms, the snapshot's pattern_capabilities finalize_reserve_ms
(areas.state_pattern_reserve), which checks._pattern_work_time reads instead of the reserve the answer states; and
(the owner's scope addition of 2026-10-18, with the service bundle that carries them) the checks that read the
entry's absence_scope or strands are no longer marked not implemented via MCP (areas.MCP_PATTERN_ENTRY_CARRIED_V7:
the MCP's answer states them per database, checks.MCP_PATTERN_DATABASE_FIELDS), so P-scope-0004 / -0006 run on the
MCP arm again.

v8 (the owner's design of 2026-10-19, DESIGN-v8 / FORMAT-v8; not frozen) is v7's scenarios, unchanged and in v7's
order, followed by area G (areas.gen_G, area_g.py): the graphlet processing tools run five ways on a corpus of
graphlets, built from the corpus file the rehearsal writes (--corpus; provisional without one: the first candidate of
every slot, preconditions pending) with the library pin --library-commit; the slots the MCP cannot express on staging
(annotate, output.sequences false, trace with coordinates) are direct-only scenarios (W1, W2, W5; the owner's decision
of 2026-10-19; --no-direct-only leaves them out). The G block is appended
after every v7 layout step, so no v7 scenario moves and the deploy gate's warm repeats stay v7's; the relations are
v7's. A v8 build is refused into manifest/ without --freeze (the owner freezes v8); `graphlets` (outside the content
digest) states the corpus file, the pin and what the freeze still needs (pending). --rehearse G (rehearse_g.py) runs
the candidates on staging and writes the corpus file and the rehearsal report.

--capabilities DIR: the capabilities snapshot (capabilities.json, pattern_capabilities.json,
traverse_capabilities.json; default truth/staging_capabilities/): the C scenarios, the relations between the three
documents and the scenarios a host capability decides (the protein substitution cap: manifest_lib.HOST_MISMATCH_MAX)
are built from it, so the frozen manifest is built with a snapshot refreshed from the host it will run against.
The frozen v1-v4 manifests stay as they are (v3 was built at commit 6f945db4; 1387455f rebuilds v4 byte for
byte). Any suite_v*.json / suite_smoke_v*.json that exists in manifest/ or that git tracks is refused as an output, by
whatever path, and --areas (a partial build) needs outputs of its own; the manifest's index records the sha256 of each
capabilities document beside the snapshot's directory.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
from typing import Optional

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import area_g  # noqa: E402
import areas  # noqa: E402
import checks  # noqa: E402
import manifest_lib  # noqa: E402
import truth  # noqa: E402
from manifest_lib import DB, Ctx, canonical_json, sha256_json  # noqa: E402

DEFAULT_SEED = 20261013
MANIFEST_VERSION = 7
# the versions this builder makes: v6 and v7 (frozen; rebuilt byte for byte) and v8 (area G; not frozen)
BUILDABLE_VERSIONS = (6, 7, 8)
GRAPHLET_VERSION = 8
DEFAULT_CORPUS = os.path.join(truth.TRUTH_DIR, "area_g", "corpus.json")
# The smoke selection's version (SMOKE_TAG): v5's selection, kept by v6 (the tag is part of every scenario's tags, and
# the gate's scenarios must stay byte-identical to v5's)
SMOKE_SELECTION_VERSION = 5
DEFAULT_CAPABILITIES = os.path.join(truth.TRUTH_DIR, "staging_capabilities")
# The v3 smoke selection (lane S3): the checks lane S3 changed, live, within 60 host requests in a paired three-arm
# run (each scenario MCP cold, direct-local, direct-remote, MCP warm; plus 9 capabilities GETs): the 503 past a
# budget (R-budget-0001), exact junction runs, a walk to a record's end, a 1 ms walk budget (a stop at depth 1),
# the telomeric seed's re-entry, an all_or_count hard region never withheld, a peptide, a count in suffix scope
# (properties the MCP may not expose), a 31-mer of a real repeat, and the max_mismatches 0 relation pair.
SMOKE_TAG = "smoke-v%d" % SMOKE_SELECTION_VERSION
# v5 adds a forced cap via MCP (T-caps-0002); the warm-ups of the smoke scenarios that have one (R-twin-0002 before
# R-budget-0001, R-warm-0006 before N-hard-0022) travel with them (tag_scenarios: `warms`), so that a smoke run meets
# them as the full run does: 14 scenarios, 28 direct requests of a paired smoke run with the closing GETs
SMOKE_IDS = ("R-budget-0001", "R-junction-0001", "T-end-0001", "T-budget-0001", "N-hard-0022", "N-hard-0005",
             "X-peptide10-0009", "P-exact-0001", "N-realrep-0007", "M-refusal-0012", "M-twin-0001", "T-caps-0002")
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
    "warms": "a warm-up's target (lane V5): the id of the scenario whose annotation rows it reads right before it; "
             "the builder keeps the pair adjacent (the warm-up first) and every selection that takes the target takes "
             "its warm-up (suite_common.select_scenarios)",
    "run_unique": "top-level request fields the direct runner makes unique per run and arm (lane S5b, A-10; lane "
                  "S5c, A2-06): it appends '-', the run id (the run directory's name, which the arms of a run share), "
                  "'-' and the arm to each, cut to the id pattern's 128 characters, and records the request it sent; "
                  "the manifest keeps the base value",
}


CAPABILITY_DOCUMENTS = ("capabilities", "pattern_capabilities", "traverse_capabilities")
# The checks a baseline run shows the MCP surface exposes, per scenario (mcp.exposes; compare.exposure_drops reads a
# listed check that reads not_exposed as a disagreement): written by --exposure-from from a recorded run
# v6: the v4 baseline with its request digests rebased to the budget-free rule of mcp_request_sha256 (--exposure-rebase:
# the exposes lists are v4-20261014T061053Z.json's, byte for byte), so that the /pattern requests of v6, which differ
# from v4's in their time_budget_ms alone, keep their lists
DEFAULT_EXPOSURE = os.path.join(truth.TRUTH_DIR, "mcp_exposure", "v4-20261014T061053Z-budget-free.json")
EXPOSURE_RULE = ("the checks decided (pass, fail or finding) on EVERY recorded MCP execution of the scenario whose "
                 "status check passed (an execution with another status decides nothing, so a scenario with one has "
                 "no list); written into a scenario's mcp.exposes only where its MCP request is the baseline's "
                 "(a time_budget_ms aside, traverse_resolve's and pattern_search's: a budget shapes no field of the "
                 "answer), restricted to the checks the scenario expects")


def load_capabilities(directory: str) -> dict:
    """The three capabilities documents of a snapshot directory."""
    caps = {}
    for name in CAPABILITY_DOCUMENTS:
        with open(os.path.join(directory, name + ".json")) as f:
            caps[name] = json.load(f)
    return caps


def snapshot_digests(directory: str) -> dict:
    """{document name: sha256 of its file} of a capabilities snapshot (A-09: the manifest records which snapshot
    built it, beside the directory's name)."""
    out = {}
    for name in CAPABILITY_DOCUMENTS:
        with open(os.path.join(directory, name + ".json"), "rb") as f:
            out[name + ".json"] = hashlib.sha256(f.read()).hexdigest()
    return out


def build(seed: int, fetch: bool, only_areas=None, capabilities_dir: str = DEFAULT_CAPABILITIES,
          exposure_file: Optional[str] = DEFAULT_EXPOSURE, manifest_version: int = MANIFEST_VERSION,
          corpus_file: Optional[str] = None, library_commit: Optional[str] = None,
          include_direct_only: bool = True) -> dict:
    if manifest_version not in BUILDABLE_VERSIONS:
        raise SystemExit("manifest version %r: this builder makes %s" % (
            manifest_version, " and ".join("v%d" % v for v in BUILDABLE_VERSIONS)))
    recs = truth.load_with_sequences(fetch=fetch)
    missing = [r.key for r in recs if r.fetch_error]
    ctx = Ctx(recs, seed)
    caps = load_capabilities(capabilities_dir)
    # the host's substitution caps per kind decide the M-protein scenarios and their by-design naming (G04, G05); the
    # routes' time-budget caps the effective budget budget_respected states above them (A2-03)
    manifest_lib.set_host_capabilities(caps["pattern_capabilities"], caps["traverse_capabilities"])
    cap = manifest_lib.HOST_BUDGET_CAP_MS.get("POST /resolve")
    if cap is not None and cap < manifest_lib.RESOLVE_TIME_BUDGET_MS:
        raise SystemExit("the snapshot's /resolve time-budget cap %d ms is below the builder's explicit-label budget "
                         "%d ms: the host would lower every resolve's budget" % (cap, manifest_lib.RESOLVE_TIME_BUDGET_MS))
    gens = [("R", areas.gen_R), ("T", areas.gen_T), ("P", areas.gen_P), ("L", areas.gen_L), ("X", areas.gen_X),
            ("M", areas.gen_M), ("Q", areas.gen_Q), ("N", areas.gen_N)]
    scenarios = []
    for area, gen in gens:
        if only_areas and area not in only_areas:
            continue
        scenarios += gen(ctx)
        if area == "T" and manifest_version >= GRAPHLET_VERSION:
            # v8 (the owner, 2026-10-19): T-walk-0024's ten warm-ups (areas.T_WALK_0024_REGION), right before it
            scenarios += areas.t_walk_0024_warm_ups(ctx, scenarios)
    if not only_areas or "C" in only_areas:
        scenarios += areas.gen_C(ctx, caps)
    if not only_areas:
        scenarios += areas.gen_REL(ctx, scenarios, caps)   # the relations' twins (needs every area)
    if not only_areas or "D" in only_areas:
        scenarios += areas.gen_D(ctx, scenarios)
    areas.state_chance_bounds(scenarios)
    # v6: what the MCP's /pattern answer cannot carry, at property level; v7: what it leaves out by design, classed so
    # v7: absence_scope and strands carried per database (areas.MCP_PATTERN_ENTRY_CARRIED_V7), no longer marked
    areas.mark_pattern_mcp_fields(scenarios, by_design=manifest_version >= 7,
                                  carried=areas.MCP_PATTERN_ENTRY_CARRIED_V7 if manifest_version >= 7 else ())
    scenario_level_marks(scenarios)
    areas.tag_pattern_budget_stress(scenarios)   # v6: the heavy /pattern scenarios outside the gate
    # v6's /pattern warm-ups must not move the deploy gate's executions: the warm-repeat status of every scenario on
    # its routes is the one of the layout without them (v5's), kept by keep_gate_warm_repeats
    reference = place_warm_ups(list(scenarios))
    keep_off_warm_repeat(reference)
    reference = gate_warm_status(reference)
    scenarios = place_warm_ups(scenarios + areas.pattern_warm_ups(scenarios))
    keep_off_warm_repeat(scenarios)
    keep_gate_warm_repeats(scenarios, reference)
    if manifest_version >= 7:
        # v7: the host's /pattern finalisation reserve in budget_respected's arguments (the warm-ups' included)
        areas.state_pattern_reserve(scenarios, caps["pattern_capabilities"])
    exposure = None
    if exposure_file:
        with open(exposure_file) as f:
            exposure = json.load(f)
        apply_exposure(scenarios, exposure)
    graphlets = None
    g_problems: list = []
    if manifest_version >= GRAPHLET_VERSION and (not only_areas or "G" in only_areas):
        # v8: area G after every v7 layout step (no v7 scenario moves; the gate's warm repeats are v7's)
        corpus = area_g.load_corpus(corpus_file)
        seeds_from = scenarios if any(s["area"] == "T" for s in scenarios) else areas.gen_T(Ctx(recs, seed))
        g = areas.gen_G(ctx, seeds_from, corpus, library_commit, include_direct_only)
        scenarios += g
        g_problems, pending = area_g.validate_g(g, _table_cases())
        graphlets = {"corpus": _corpus_identity(corpus_file, corpus), "library_commit": library_commit,
                     "include_direct_only": include_direct_only, "pending": pending, "summary": area_g.summary(g)}
    tag_scenarios(scenarios, recs)
    relations = ctx.relations
    problems = validate(scenarios) + validate_relations(relations, scenarios) + warm_repeat_problems(scenarios) + \
        gate_warm_problems(scenarios, reference) + g_problems
    if problems:
        raise SystemExit("manifest invalid:\n  " + "\n  ".join(problems[:30]))
    tc = caps["traverse_capabilities"]
    manifest = {
        "manifest_version": manifest_version,
        "content_sha256": sha256_json(scenarios),
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seed": seed,
        "index": {"database": DB, "release_at_build": caps["capabilities"].get("release"),
                  "index_fp": tc.get("index_fp"), "index_ns": tc.get("index_ns"), "k": tc.get("k"),
                  "capabilities_snapshot": _snapshot_name(capabilities_dir),
                  "capabilities_sha256": snapshot_digests(capabilities_dir),
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
        "mcp_exposure_baseline": None if exposure is None else {
            "file": os.path.relpath(os.path.abspath(exposure_file), HERE), "source": exposure.get("source"),
            "rule": exposure.get("rule"), "sha256": sha256_json(exposure),
            "scenarios_with_exposes": sum(1 for sc in scenarios if (sc.get("mcp") or {}).get("exposes"))},
        "skipped": ctx.skipped,
        "not_applicable": not_applicable(caps),
        "counts": counts(scenarios, relations),
        "relations_sha256": sha256_json(relations),
        "relations": relations,
        "scenarios": scenarios,
    }
    if graphlets is not None:
        manifest = dict((k, v) for k, v in manifest.items() if k not in ("relations", "scenarios"))
        manifest["graphlets"] = graphlets
        manifest["relations"] = relations
        manifest["scenarios"] = scenarios
    return manifest


def _corpus_identity(path: Optional[str], corpus: dict) -> dict:
    """The corpus file a v8 build read (relative to this directory when inside it, its sha256), or the provisional
    corpus."""
    if not path:
        return {"file": None, "provisional": True,
                "rule": "no rehearsal yet: the first candidate of every slot (direct-only ones by default), "
                        "preconditions pending"}
    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    rel = os.path.relpath(os.path.abspath(path), HERE)
    return {"file": rel if not rel.startswith("..") else os.path.abspath(path), "sha256": digest,
            "provisional": bool(corpus.get("provisional")), "rehearsal": corpus.get("rehearsal")}


# the tools whose time_budget_ms the request digest leaves out (EXPOSURE_RULE): traverse_resolve's since v5,
# pattern_search's since v6 (the budget 60,000 -> 110,000 ms would otherwise drop every /pattern scenario's list).
# EXPOSURE_RULE's "a budget shapes no field of the answer" is the short form, embedded in the frozen manifests: a
# pattern_search budget decides how far the answer gets (complete, or limited by time: SPEC-pattern-search §7.9),
# never which fields it carries, and a larger budget lets no fewer checks decide, so a list recorded under one budget
# is a lower bound under a larger one (v6's 110,000 ms against v4's 60,000 ms)
BUDGET_FREE_DIGEST_TOOLS = ("traverse_resolve", "pattern_search")


def mcp_request_sha256(mcp: dict) -> Optional[str]:
    """The digest of a scenario's MCP request (tool and arguments; the time_budget_ms of BUDGET_FREE_DIGEST_TOOLS left
    out), None for a scenario the MCP does not run."""
    if not isinstance(mcp, dict) or not mcp.get("tool"):
        return None
    args = dict(mcp.get("arguments") or {})
    if mcp["tool"] in BUDGET_FREE_DIGEST_TOOLS:
        args.pop("time_budget_ms", None)
    return sha256_json({"tool": mcp["tool"], "arguments": args})


def rebase_exposure(baseline: dict, source_manifest: dict, baseline_file: str) -> dict:
    """|baseline| with every scenario's mcp_request_sha256 recomputed by the current mcp_request_sha256 from the MCP
    request of |source_manifest|, the manifest the baseline was recorded against (its source.manifest_content_sha256
    must be that manifest's content_sha256, and every scenario of the baseline must be in it); the exposes lists
    and the source stay, the rule is the current one, and source.rebased_from names the file rebased (its sha256).
    A ValueError names what does not fit."""
    want = (baseline.get("source") or {}).get("manifest_content_sha256")
    if want != source_manifest.get("content_sha256"):
        raise ValueError("the baseline was recorded against the manifest content %s, not the given one (%s)" % (
            want, source_manifest.get("content_sha256")))
    scen = dict((s.get("id"), s) for s in source_manifest.get("scenarios") or [])
    missing = sorted(sid for sid in (baseline.get("scenarios") or {}) if sid not in scen)
    if missing:
        raise ValueError("baseline scenarios not in the manifest: %s" % ", ".join(missing[:10]))
    out = {}
    for sid, entry in sorted((baseline.get("scenarios") or {}).items()):
        out[sid] = {"exposes": list(entry.get("exposes") or []),
                    "mcp_request_sha256": mcp_request_sha256(scen[sid].get("mcp"))}
    with open(baseline_file, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    source = dict(baseline.get("source") or {})
    source["rebased_from"] = {"file": os.path.relpath(os.path.abspath(baseline_file), HERE), "sha256": digest,
                              "rule": baseline.get("rule"),
                              "digest_left_out": "time_budget_ms of %s" % ", ".join(BUDGET_FREE_DIGEST_TOOLS)}
    return {"source": source, "rule": EXPOSURE_RULE, "scenarios": out}


def exposure_baseline(rows: list, manifest: dict, source: dict) -> dict:
    """The exposure baseline of a recorded run (EXPOSURE_RULE): {source, rule, scenarios: {id: {exposes,
    mcp_request_sha256}}} from the run's execution lines |rows| (the MCP arm's) and its manifest."""
    scen = dict((sc.get("id"), sc) for sc in manifest.get("scenarios") or [])
    per: dict = {}
    for r in rows:
        if r.get("arm") != "mcp" or r.get("scenario") not in scen:
            continue
        checks_ = r.get("checks") or []
        status_ok = any(c.get("check") == "status" and c.get("outcome") == "pass" for c in checks_)
        decided = set(c.get("check") for c in checks_ if c.get("outcome") in ("pass", "fail", "finding")) \
            if status_ok else set()
        per.setdefault(r["scenario"], []).append(decided)
    out = {}
    for sid, sets in sorted(per.items()):
        common = set.intersection(*sets) if sets else set()
        out[sid] = {"exposes": sorted(common), "mcp_request_sha256": mcp_request_sha256(scen[sid].get("mcp"))}
    return {"source": source, "rule": EXPOSURE_RULE, "scenarios": out}


def apply_exposure(scenarios: list, baseline: Optional[dict]) -> int:
    """Fill mcp.exposes from |baseline| (EXPOSURE_RULE); returns how many scenarios got a list."""
    n = 0
    for sc in scenarios:
        b = ((baseline or {}).get("scenarios") or {}).get(sc["id"])
        sha = mcp_request_sha256(sc.get("mcp"))
        if not b or sha is None or sha != b.get("mcp_request_sha256"):
            continue
        expected = set(e["check"] for e in sc["expect"])
        exposes = [c for c in b.get("exposes") or [] if c in expected]
        if exposes:
            sc["mcp"]["exposes"] = exposes
            n += 1
    return n


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


def _marked_checks(mcp: dict) -> list:
    """The checks a scenario's property-level MCP marks name (mcp.not_implemented_checks, mcp.by_design_checks)."""
    return [e.get("check") for key in ("not_implemented_checks", "by_design_checks") for e in mcp.get(key) or []
            if isinstance(e, dict)]


def covers_every_property(sc: dict) -> bool:
    """Whether the property-level MCP marks of |sc| name every property check it expects (checks.is_property)."""
    marked = set(_marked_checks(sc.get("mcp") or {}))
    props = set(e["check"] for e in sc["expect"] if checks.is_property(e["check"]))
    return bool(marked) and bool(props) and props <= marked


def scenario_level_marks(scenarios: list) -> list:
    """The owner's rule (2026-10-15): a missing MCP feature is marked at the narrowest scope; when the property-level
    marks would cover every property check of a scenario, that scope is the scenario. Its MCP mapping becomes
    mcp.not_implemented (the manifest's whole-scenario form; it has no whole-scenario by-design form) with the marks'
    reasons, a service-side case stated as "service-side by design: ...", and keeps mcp.by_design. Otherwise the
    scenario would pass via MCP on its generic checks alone. Returns the ids changed."""
    out = []
    for sc in scenarios:
        mcp = sc.get("mcp") or {}
        if "tool" not in mcp or not covers_every_property(sc):
            continue
        why = ["service-side by design: the table's case %s; the MCP answer cannot carry %s" % (e["case"], e["check"])
               for e in mcp.get("by_design_checks") or []]
        why += ["not implemented yet via MCP: %s (%s)" % (e["reason"], e["check"])
                for e in mcp.get("not_implemented_checks") or []]
        sc["mcp"] = dict({"not_implemented": "; ".join(why)},
                         **({"by_design": mcp["by_design"]} if mcp.get("by_design") else {}))
        out.append(sc["id"])
    return out


def place_warm_ups(scenarios: list) -> list:
    """|scenarios| with each warm-up (tag warm-up, `warms`: the id of the scenario it warms) moved to just before the
    scenario it warms, several in the order built. A runner runs an arm's scenarios in manifest order (the paired
    loop runs each scenario's executions, MCP cold first, before the next scenario's), so the warmed scenario's first
    execution, the MCP arm's cold one, meets the rows the warm-up read (triage-v4 G22 / G24: cold reads under the
    load decided the budget). A warm-up whose target is not in the manifest stays where it is."""
    ids = set(s["id"] for s in scenarios)
    warm = {}
    rest = []
    for s in scenarios:
        target = s.get("warms") if "warm-up" in s.get("tags", []) else None
        if target in ids:
            warm.setdefault(target, []).append(s)
        else:
            rest.append(s)
    out = []
    for s in rest:
        out += warm.get(s["id"], [])
        out.append(s)
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
        if sc.get("id") in SMOKE_IDS or sc.get("warms") in SMOKE_IDS:   # a warm-up travels with its target
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
                       for t in ("last-record-verified", "code-derived", "machine-dependent", "index-empirical",
                                 "budget_stress", "warm-up", "smoke", SMOKE_TAG)}}


def area_list(text):
    """The areas of --areas ('P, T' as 'P,T', empties left out; SD-369: 'A, B' dropped B without a word), None for
    all."""
    out = set(a.strip() for a in (text or "").split(",") if a.strip())
    return out or None


# a step of a same_answer args.ignore path that names no object key: a list index (negative too) or '*'
_INDEX_STEP = re.compile(r"(-?[0-9]+|\*)\Z")
_NO_KEY_CHARS = re.compile(r"[\[\]*\s]")


def _ignore_problems(ignore) -> list:
    """Why a same_answer args.ignore is not a list of dotted object-key paths (checks._drop and compare.relation_drop:
    a path names object keys and looks through lists): a list index, '*', a bracket or white space in a step matches
    nothing in either reader (the owner's decision on SD-011: refused, not given a meaning), nor does an empty step."""
    if not isinstance(ignore, list):
        return ["args.ignore %s is not a list of paths" % json.dumps(ignore)]
    out = []
    for path in ignore:
        if not isinstance(path, str):
            out.append("args.ignore %s is not a dotted path" % json.dumps(path))
            continue
        bad = [step for step in path.split(".")
               if step == "" or _INDEX_STEP.match(step) or _NO_KEY_CHARS.search(step)]
        if bad:
            out.append("args.ignore %r: the step %r is no object key (an index, '[', ']', '*', white space or an "
                       "empty step matches nothing; a path looks through lists)" % (path, bad[0]))
    return out


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
        if r["type"] == "same_answer" and "ignore" in r["args"]:
            problems += ["%s: %s" % (r["id"], m) for m in _ignore_problems(r["args"]["ignore"])]
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


# The paired plan adds a warm repeat to every WARM_EVERY-th scenario of the manifest (suite_common.repeat_count, by
# the scenario's index; run_suite.sh's and gate_policy.py's default --warm-every, the async runner's too)
WARM_EVERY = 10


def _wants_200(s) -> bool:
    return 200 in manifest_lib.allowed_statuses(s.get("expect"))


def single_execution(s) -> str:
    """Why one execution of |s| must decide alone ('' when a warm repeat may follow): a run-unique scenario (run_unique:
    its id is the run's and the arm's) that expects the attempt's 200, whose repeat would resend the id and meet the
    retained attempt's 409 (T-attempts-0002; VERIFY-V5A; one that expects the 409, T-attempts-0003, meets it again); a
    stated-stop design (truth.stop_design: R-budget), whose every repeat must stop for the gate (VERIFY-V5A's
    re-verification: R-budget-0002 at index 50 had two)."""
    if s.get("run_unique") and _wants_200(s):
        return "run-unique"
    if isinstance(s.get("truth"), dict) and s["truth"].get("stop_design"):
        return "a stated-stop design"
    return ""


def warm_repeat_problems(scenarios, warm_every: int = WARM_EVERY) -> list:
    """A single-execution scenario (single_execution) must not get the paired plan's warm repeat (every
    |warm_every|-th scenario by manifest index: suite_common.repeat_count, which async's runner mirrors). Checked on
    every build, the frozen one included."""
    return ["%s: %s at manifest index %d, which the warm repeat (every %dth scenario) would repeat" % (
                s["id"], single_execution(s), i, warm_every)
            for i, s in enumerate(scenarios)
            if single_execution(s) and warm_every and i % warm_every == 0]


def keep_off_warm_repeat(scenarios: list, warm_every: int = WARM_EVERY) -> list:
    """Moves |scenarios| (in place, after place_warm_ups) so that no single-execution scenario sits on a warm-repeat
    index. A scenario may move unless it is a warm-up, a scenario a warm-up warms, run-unique or single-execution; the
    scenarios around one on that index that may not move form its unit (its warm-ups, a run-unique pair such as
    T-attempts-0002 / -0003), and the nearest scenario before the unit in its area that may move goes right after the
    unit, which shifts the unit one index down; the moved one gets the warm repeat. Draws, ids and every pair stay as
    built. Returns [(moved id, the id it now follows)]."""
    moved = []
    for _ in range(len(scenarios)):
        bad = [i for i, s in enumerate(scenarios) if single_execution(s) and warm_every and i % warm_every == 0]
        if not bad:
            return moved
        i = bad[0]
        warmed = set(s.get("warms") for s in scenarios if s.get("warms"))

        def free(s):
            return not (s.get("warms") or s.get("id") in warmed or s.get("run_unique") or single_execution(s)
                        or "warm-up" in (s.get("tags") or []))
        first, last = i, i
        while first > 0 and not free(scenarios[first - 1]):
            first -= 1
        while last + 1 < len(scenarios) and not free(scenarios[last + 1]):
            last += 1
        j = next((j for j in range(first - 1, -1, -1) if scenarios[j]["area"] == scenarios[i]["area"]
                  and free(scenarios[j])), None)
        if j is None:
            raise SystemExit("%s: on the warm-repeat index %d, and nothing before it in its area can move" % (
                scenarios[i]["id"], i))
        x = scenarios.pop(j)
        scenarios.insert(last, x)
        moved.append((x["id"], scenarios[last - 1]["id"]))
    raise SystemExit("keep_off_warm_repeat: no placement found")


# The deploy gate's routes (gate_policy.DEFAULT_ROUTES, pinned equal by test_builder_v6): the scenarios whose executions
# a manifest change must not move
GATE_ROUTES = ("POST /resolve", "POST /traverse", "GET /traverse/capabilities")


def _warm_repeated(i: int, s: dict, warm_every: int) -> bool:
    """Whether the paired plan gives the scenario at index |i| the warm repeat (suite_common.repeat_count: a
    scenario of one repeat on every |warm_every|-th index)."""
    return s.get("repeat", 1) == 1 and bool(warm_every) and i % warm_every == 0


def gate_warm_status(scenarios: list, warm_every: int = WARM_EVERY) -> dict:
    """{id: whether it gets the paired plan's warm repeat} of the scenarios on the deploy gate's routes (GATE_ROUTES)."""
    return dict((s["id"], _warm_repeated(i, s, warm_every)) for i, s in enumerate(scenarios)
                if s.get("route") in GATE_ROUTES)


def keep_gate_warm_repeats(scenarios: list, reference: dict, warm_every: int = WARM_EVERY) -> list:
    """Moves |scenarios| (in place, after keep_off_warm_repeat) so that every scenario on the deploy gate's routes has
    the warm-repeat status |reference| (gate_warm_status of the layout without a version's added warm-ups) states:
    scenarios added before it shift its index, and with it the gate's required executions. As keep_off_warm_repeat:
    the scenarios around one whose status differs that may not move form its unit (the gate's own scenarios, warm-ups
    and their targets, run-unique or single-execution ones), and the nearest scenario before the unit in its area that
    may move (a /pattern scenario, say) goes right after the unit, which shifts the unit one index down; the gate's
    scenarios, their content and their order stay. Returns [(moved id, the id it now follows)]."""
    moved = []
    warmed = set(s.get("warms") for s in scenarios if s.get("warms"))

    def free(s):
        return not (s.get("route") in GATE_ROUTES or s.get("warms") or s.get("id") in warmed or s.get("run_unique")
                    or single_execution(s) or "warm-up" in (s.get("tags") or []))
    for _ in range(10 * len(scenarios)):
        bad = [i for i, s in enumerate(scenarios) if s["id"] in reference
               and _warm_repeated(i, s, warm_every) != reference[s["id"]]]
        if not bad:
            return moved
        i = bad[0]
        first, last = i, i
        while first > 0 and not free(scenarios[first - 1]):
            first -= 1
        while last + 1 < len(scenarios) and not free(scenarios[last + 1]):
            last += 1
        # in its area, else (a scenario placed in another area's block, a warm-up, as R-warm-0006 among N's) in the
        # area of the unit's first scenario, else any
        j = None
        for area in (scenarios[i]["area"], scenarios[first]["area"], None):
            j = next((j for j in range(first - 1, -1, -1) if area in (None, scenarios[j]["area"])
                      and free(scenarios[j])), None)
            if j is not None:
                break
        if j is None:
            raise SystemExit("%s: its warm repeat would change (index %d), and nothing before it can move" % (
                scenarios[i]["id"], i))
        x = scenarios.pop(j)
        scenarios.insert(last, x)
        moved.append((x["id"], scenarios[last - 1]["id"]))
    raise SystemExit("keep_gate_warm_repeats: no placement found")


def gate_warm_problems(scenarios: list, reference: dict, warm_every: int = WARM_EVERY) -> list:
    """A scenario on the deploy gate's routes whose warm-repeat status is not the reference's (keep_gate_warm_repeats)."""
    return ["%s: on the deploy gate's routes, its warm repeat %s by the index %d" % (
                s["id"], "added" if _warm_repeated(i, s, warm_every) else "removed", i)
            for i, s in enumerate(scenarios)
            if s["id"] in reference and _warm_repeated(i, s, warm_every) != reference[s["id"]]]


def _table_cases() -> set:
    """The case ids of the service-by-design table beside this script."""
    with open(os.path.join(HERE, "service_by_design.json")) as f:
        return set(c.get("id") for c in json.load(f).get("cases") or [])


def validate(scenarios, cases=None) -> list:
    """Unique ids, every check known to checks.py, the contract's fields present; every mcp.by_design id a case of the
    table |cases| (the bundled service_by_design.json: a case the table lacks explains nothing, so the scenario's
    MCP answer would read unexplained); mcp.not_implemented_checks and mcp.by_design_checks name checks the scenario
    expects (with a reason, a case of the table)."""
    problems = []
    seen = set()
    cases = _table_cases() if cases is None else cases
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
        n_status = len([e for e in s["expect"] if e["check"] == "status"])
        if not n_status:
            problems.append("%s: no status expectation" % s["id"])
        elif n_status > 1:
            # SD-094: one, so that every reader of the allowed statuses (first, last or every entry) reads the same
            problems.append("%s: %d status expectations, not one" % (s["id"], n_status))
        expected = set(e["check"] for e in s["expect"])
        by_design = s["mcp"].get("by_design") or []
        if not isinstance(by_design, list):
            # SD-030: never iterated character by character; every manifest states a list
            problems.append("%s: mcp.by_design %s is not a list of case ids" % (s["id"], json.dumps(by_design)))
            by_design = []
        unknown = [c for c in by_design if c not in cases]
        if unknown:
            problems.append("%s: mcp.by_design names %s, not a case of the bundled by-design table" % (
                s["id"], ", ".join(unknown)))
        for e in s["mcp"].get("not_implemented_checks") or []:
            if not isinstance(e, dict) or e.get("check") not in expected or e.get("check") not in checks.CHECKS \
                    or not str(e.get("reason") or "").strip():
                problems.append("%s: mcp.not_implemented_checks %r: a check it expects, with a reason" % (s["id"], e))
        for e in s["mcp"].get("by_design_checks") or []:
            if not isinstance(e, dict) or e.get("check") not in expected or e.get("case") not in cases:
                problems.append("%s: mcp.by_design_checks %r: a check it expects and a case of the table" % (s["id"], e))
        if (s["mcp"].get("not_implemented_checks") or s["mcp"].get("by_design_checks")) and "tool" not in s["mcp"]:
            problems.append("%s: property-level MCP classes on a scenario the MCP does not run" % s["id"])
        elif covers_every_property(s):
            problems.append("%s: property-level MCP marks cover every property check: a scenario-level mark "
                            "(mcp.not_implemented, scenario_level_marks)" % s["id"])
        for e in s["expect"]:
            if not isinstance(e.get("args", {}), dict):
                problems.append("%s: the expectation %s's args %s is not an object" % (
                    s["id"], e.get("check"), json.dumps(e["args"])))
                continue
            problems += ["%s: %s" % (s["id"], m) for m in checks.validate_expectation(e)]
        disc = s["request"].get("discover") if s["route"] == "POST /resolve" and isinstance(s["request"], dict) else None
        if isinstance(disc, dict) and _status_ok(s):
            # one question on both arms: the MCP sends its limit and kind as they are (50 when none, refused above
            # 1,000), the host's default is 1,000 (the async session's constraint for v5)
            ml = disc.get("max_labels")
            by = (s.get("mcp") or {}).get("by_design") or []
            cap = manifest_lib.SERVICE_DISCOVER_MAX_LABELS
            if "max_labels" not in disc or "kind" not in disc:
                problems.append("%s: a discovery that does not name max_labels and kind" % s["id"])
            elif isinstance(ml, int) and ml > cap and "resolve.discover_limit_above_cap_host_answers" not in by:
                problems.append("%s: discover.max_labels %d above the MCP's %d outside its by-design case" % (
                    s["id"], ml, cap))
            elif "tool" in s["mcp"] and (s["mcp"]["arguments"].get("limit") != ml
                                         or s["mcp"]["arguments"].get("kind") != disc.get("kind")):
                problems.append("%s: the MCP arguments do not state the discovery's limit and kind" % s["id"])
        cb = (s.get("truth") or {}).get("chance_bound") if isinstance(s.get("truth"), dict) else None
        if cb and cb.get("p", 0) > areas.CHANCE_BOUND_MAX:
            problems.append("%s: a negative control whose chance bound %g exceeds %g" % (
                s["id"], cb["p"], areas.CHANCE_BOUND_MAX))
        req = s["request"]
        if s["route"] == "POST /pattern" and isinstance(req, dict) and "patterns" in req:
            for f in ("mode", "max_contexts", "time_budget_ms"):
                if f not in req:
                    problems.append("%s: /pattern request without %s" % (s["id"], f))
    return problems


def _status_ok(s: dict) -> bool:
    """Whether the scenario's status expectation admits a 200 (a request answered, not a refusal)."""
    return 200 in manifest_lib.allowed_statuses(s.get("expect"))


def smoke_of(manifest: dict) -> dict:
    sel = [s for s in manifest["scenarios"] if SMOKE_TAG in s["tags"] and "needs-next-deploy" not in s["tags"]]
    missing = sorted(set(SMOKE_IDS) - {s["id"] for s in sel})
    if missing and not manifest.get("partial_build"):
        raise SystemExit("smoke scenarios not in the manifest: %s" % ", ".join(missing))
    m = dict(manifest)
    m.pop("graphlets", None)   # v8: area G is not in the smoke selection
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


# a manifest file of manifest/ (any version): never written over once it exists or git tracks it
MANIFEST_NAME = re.compile(r"^suite(?:_smoke)?_v[0-9]+\.json$")
MANIFEST_DIR = os.path.join(HERE, "manifest")


def _tracked(path: str) -> bool:
    """Whether git tracks |path| in this tree."""
    try:
        p = subprocess.run(["git", "ls-files", "--error-unmatch", os.path.basename(path)], cwd=os.path.dirname(path),
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
        return p.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def frozen_target(path: str) -> str:
    """Why |path| must not be written ('' when it may): an existing or git-tracked suite_v*.json / suite_smoke_v*.json
    of manifest/, reached by any path (realpath, and samefile against each of manifest/'s manifests: a symbolic link
    or another spelling is the same file; A-08)."""
    real = os.path.realpath(path)
    candidates = [real] + [p for p in glob.glob(os.path.join(MANIFEST_DIR, "suite*_v*.json"))
                           if os.path.exists(path) and os.path.samefile(p, path)]
    for c in candidates:
        if os.path.dirname(c) != os.path.realpath(MANIFEST_DIR) or not MANIFEST_NAME.match(os.path.basename(c)):
            continue
        if os.path.exists(c):
            return "%s is a manifest of manifest/ (%s): build into another file" % (path, os.path.basename(c))
        if _tracked(c):
            return "%s is a manifest git tracks: build into another file" % path
    return ""


def manifest_path(path: str) -> bool:
    """Whether |path|, by any spelling (realpath: relative, `..`, a symbolic link to manifest/ or to a file there),
    names a manifest file of manifest/ (MANIFEST_NAME), existing or not (A2-05 of the review of lane S5b)."""
    real = os.path.realpath(path)
    return os.path.dirname(real) == os.path.realpath(MANIFEST_DIR) and bool(MANIFEST_NAME.match(os.path.basename(real)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=None, help="default manifest/suite_v<version>.json")
    ap.add_argument("--smoke-out", default=None, help="default manifest/suite_smoke_v<version>.json")
    ap.add_argument("--manifest-version", type=int, default=MANIFEST_VERSION, choices=BUILDABLE_VERSIONS,
                    help="the version to build (default %d; 6 rebuilds the frozen v6)" % MANIFEST_VERSION)
    ap.add_argument("--capabilities", default=DEFAULT_CAPABILITIES,
                    help="the capabilities snapshot directory (default truth/staging_capabilities/)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--no-fetch", action="store_true", help="use the record cache only")
    ap.add_argument("--areas", help="comma-separated areas (development only; the frozen manifest has all)")
    ap.add_argument("--exposure", default=DEFAULT_EXPOSURE,
                    help="the MCP exposure baseline (default %s; '' for none)" % os.path.relpath(DEFAULT_EXPOSURE, HERE))
    ap.add_argument("--exposure-from", metavar="RUN_DIR",
                    help="write the exposure baseline of a recorded run (its executions.reevaluated.jsonl, else "
                         "executions.jsonl, and manifest.json) to --exposure-out, then exit")
    ap.add_argument("--exposure-out", help="where --exposure-from / --exposure-rebase writes")
    ap.add_argument("--exposure-rebase", metavar="BASELINE",
                    help="write BASELINE with its request digests recomputed by the current rule from the manifest it "
                         "was recorded against (--exposure-source-manifest) to --exposure-out, then exit")
    ap.add_argument("--exposure-source-manifest", metavar="MANIFEST",
                    help="the manifest --exposure-rebase's baseline was recorded against (manifest/suite_v4.json)")
    ap.add_argument("--corpus", default=None,
                    help="v8: the area G corpus file the rehearsal wrote (default truth/area_g/corpus.json when it "
                         "exists, else the provisional corpus)")
    ap.add_argument("--library-commit", default=None,
                    help="v8: the library commit the manifest is frozen against (graphlet.library_commit)")
    ap.add_argument("--include-direct-only", action="store_true", default=True,
                    help="v8: the area G slots the MCP cannot express as direct-only scenarios (W1, W2, W5): the "
                         "default (the owner, 2026-10-19)")
    ap.add_argument("--no-direct-only", dest="include_direct_only", action="store_false",
                    help="v8: leave the direct-only area G slots out")
    ap.add_argument("--freeze", action="store_true",
                    help="v8: allow the default outputs manifest/suite_v8.json / suite_smoke_v8.json (the owner's "
                         "freeze)")
    ap.add_argument("--rehearse", metavar="AREA", default=None,
                    help="G: the area G rehearsal (rehearse_g.py; its own arguments follow --)")
    args, rest = ap.parse_known_args(argv)
    if args.rehearse:
        if args.rehearse != "G":
            raise SystemExit("--rehearse: only area G has a rehearsal")
        import rehearse_g
        return rehearse_g.main([x for x in rest if x != "--"])
    if rest:
        ap.error("unrecognized arguments: %s" % " ".join(rest))
    args.out = args.out or os.path.join(HERE, "manifest", "suite_v%d.json" % args.manifest_version)
    args.smoke_out = args.smoke_out or os.path.join(HERE, "manifest", "suite_smoke_v%d.json" % args.manifest_version)
    if args.exposure_rebase:
        if not args.exposure_out or not args.exposure_source_manifest:
            raise SystemExit("--exposure-rebase needs --exposure-source-manifest and --exposure-out")
        with open(args.exposure_rebase) as f:
            baseline = json.load(f)
        with open(args.exposure_source_manifest) as f:
            source_manifest = json.load(f)
        try:
            doc = rebase_exposure(baseline, source_manifest, args.exposure_rebase)
        except ValueError as e:
            raise SystemExit("--exposure-rebase: %s" % e)
        os.makedirs(os.path.dirname(os.path.abspath(args.exposure_out)), exist_ok=True)
        with open(args.exposure_out, "w") as f:
            json.dump(doc, f, indent=1, sort_keys=True)
            f.write("\n")
        print("%d scenarios, %d with exposed checks, digests rebased" % (
            len(doc["scenarios"]), sum(1 for v in doc["scenarios"].values() if v["exposes"])))
        return 0
    if args.exposure_from:
        if not args.exposure_out:
            raise SystemExit("--exposure-from needs --exposure-out")
        run = args.exposure_from
        lines = os.path.join(run, "executions.reevaluated.jsonl")
        if not os.path.exists(lines):
            lines = os.path.join(run, "executions.jsonl")
        with open(lines) as f:
            rows = [json.loads(x) for x in f if x.strip()]
        with open(os.path.join(run, "manifest.json")) as f:
            man = json.load(f)
        checks_sha = sorted(set((r.get("reevaluation") or {}).get("checks_sha256") for r in rows
                                if r.get("arm") == "mcp") - {None})
        run_ids = sorted(set(r.get("run_id") for r in rows if r.get("arm") == "mcp") - {None})
        source = {"run": run_ids or [os.path.basename(os.path.normpath(run))], "lines": os.path.basename(lines),
                  "manifest_content_sha256": man.get("content_sha256"), "checks_sha256": checks_sha}
        doc = exposure_baseline(rows, man, source)
        os.makedirs(os.path.dirname(os.path.abspath(args.exposure_out)), exist_ok=True)
        with open(args.exposure_out, "w") as f:
            json.dump(doc, f, indent=1, sort_keys=True)
            f.write("\n")
        print("%d scenarios, %d with exposed checks" % (len(doc["scenarios"]),
                                                       sum(1 for v in doc["scenarios"].values() if v["exposes"])))
        return 0
    only = area_list(args.areas)
    if args.manifest_version >= GRAPHLET_VERSION and not args.freeze and (
            manifest_path(args.out) or manifest_path(args.smoke_out)):
        raise SystemExit("v%d is not frozen: build into files outside manifest/ (--out, --smoke-out), or pass --freeze "
                         "(the owner's freeze)" % args.manifest_version)
    if args.corpus is None and args.manifest_version >= GRAPHLET_VERSION and os.path.exists(DEFAULT_CORPUS):
        args.corpus = DEFAULT_CORPUS
    if only and (manifest_path(args.out) or manifest_path(args.smoke_out)):
        raise SystemExit("--areas builds a partial manifest: give --out and --smoke-out outside manifest/ (not a "
                         "suite_v*.json / suite_smoke_v*.json there, by any spelling)")
    for path in (args.out, args.smoke_out):
        why = frozen_target(path)
        if why:
            raise SystemExit(why)
    m = build(args.seed, fetch=not args.no_fetch, only_areas=only, capabilities_dir=args.capabilities,
              exposure_file=args.exposure or None, manifest_version=args.manifest_version, corpus_file=args.corpus,
              library_commit=args.library_commit, include_direct_only=args.include_direct_only)
    if only:
        m["partial_build"] = True
    write(args.out, m)
    s = smoke_of(m)
    write(args.smoke_out, s)
    for sk in m["skipped"]:
        print("skipped: %(area)s %(scenario)s: %(reason)s" % sk, file=sys.stderr)
    for p in (m.get("graphlets") or {}).get("pending") or []:
        print("pending (area G, before the freeze): %s" % p, file=sys.stderr)
    print(json.dumps({"full": {"path": args.out, "content_sha256": m["content_sha256"], **m["counts"]},
                      "smoke": {"path": args.smoke_out, "content_sha256": s["content_sha256"], **s["counts"]}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
