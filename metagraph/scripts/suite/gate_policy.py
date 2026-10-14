#!/usr/bin/env python3
"""Writes a deploy-gate policy for gate.py from a manifest, and prints its counts.

    gate_policy.py MANIFEST [--routes R1,R2] [--arms mcp,direct-local] [--by-design PATH] [--advanced-path]
                   [--warm-every 10] [--require-agreement] [--declared FILE] [--name NAME] [--out POLICY]
                   [--ids-out FILE]

The selection is by ROUTE (default: POST /resolve, POST /traverse and GET /traverse/capabilities: every scenario the
resolve / traverse service answers, the N and D scenarios on those routes included, not only areas R and T). For the
MCP arm the scenarios it cannot express are explicit exclusions, each derived from the manifest and the by-design
table, never from an answer a run received:
  manifest_not_implemented  the manifest's mcp.not_implemented (no MCP call maps the request);
  table_not_expressible     a case the scenario names (mcp.by_design) whose table agreement is not_implemented_mcp
                            (e.g. resolve.not_expressible.bounds_max_query_bp: the MCP sends another request), unless
                            the case cannot change this request (not_expressible_moot: a query within the
                            service's own cap);
  advanced_path             traverse_fetch arguments on the advanced strategy path (an `overrides` key outside the
                            MCP's basic list, a preset other than walk / trie, a `strategy` object), which a deployment
                            without TRAVERSAL_ADVANCED_ENABLED refuses; --advanced-path for one that enables it;
  needs_direct_arm          a policy without a direct arm only: a named answered / clamped case that explains checks
                            (the service answers what the host refuses, e.g. resolve.duplicate_labels): compare
                            verifies such a pair only against a direct execution.
The deploy gate runs two arms, mcp and direct-local (paired): direct-local verifies the answered pairs and lets the
gate catch an MCP vs direct disagreement. Owner-declared entries come from --declared FILE, a versioned JSON file
beside the policy: {declared_version 1, manifest {content_sha256}, exclusions [{scenario, arms, repeats (optional:
only these repeats), reason}], not_required [{scenario, arms, check, reason}], relations_excluded [{relation, arms,
reason}]}; an exclusion has the rule `declared`, a not_required entry exempts a property check the scenario expects on
those arms (the MCP never exposes it: it stays required on direct-local); each entry needs a reason, names a selected
scenario (a manifest relation) and (not_required) a check the scenario expects, and the file names the manifest it
was written for. No property check is exempt otherwise.
A relation is required on an arm when both of its scenarios are selected and neither is excluded there; the others are
listed with the reason. The policy pins the manifest's content_sha256, relations_sha256 and index_fp and
the by-design table's sha256 (gate.py fails on another), and lists the executions each arm must record (the paired
plan: suite_common.mcp_repeats for the MCP arm, plan_executions for a direct arm, with --warm-every). --ids-out writes
the ids the policy's arms run (the selection less what is excluded on all of them), one per line in manifest order:
run_suite.sh --ids FILE, or async's runner standalone (--ids, comma-separated: gate.py policy-ids POLICY).
"""

import argparse
import hashlib
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as cmp  # noqa: E402
import gate  # noqa: E402
import suite_common as sc  # noqa: E402

DEFAULT_ROUTES = ("POST /resolve", "POST /traverse", "GET /traverse/capabilities")
# The traverse_fetch overrides the MCP passes on while the advanced strategy path is off (strategy key -> the sub-keys
# allowed, None = the whole key): the service's _TRAVERSE_BASIC_OVERRIDES (app/mcp_server.py, async's runner branch
# at e334715b); traverse_fetch's documented presets are walk and trie.
MCP_BASIC_OVERRIDES = {"bounds": None, "direction": None, "exhaustive": None,
                       "labels": ("max_labels_per_node", "seed_label_kind")}
MCP_BASIC_PRESETS = (None, "walk", "trie")
# The service always sends its own /resolve query cap as bounds.max_query_bp (app/traversal/resolve.py; app/settings.py
# TRAVERSAL_RESOLVE_MAX_QUERY_BP 2000, TRAVERSAL_RESOLVE_DISCOVER_MAX_QUERY_BP 500 for a discovery): the case stating
# bounds.max_query_bp not expressible changes nothing for a query no longer than that cap whose direct request states
# no bound of its own.
SERVICE_RESOLVE_MAX_QUERY_BP = {"labels": 2000, "discover": 500}


def not_expressible_moot(case, scenario):
    """Whether a not-expressible case cannot change a scenario's MCP request: today only the resolve query cap
    (SERVICE_RESOLVE_MAX_QUERY_BP) on a query within it; every other not-expressible case applies."""
    if (case.get("service_answer") or {}).get("missing_argument") != "bounds.max_query_bp":
        return False
    args = (scenario.get("mcp") or {}).get("arguments") or {}
    request = scenario.get("request") if isinstance(scenario.get("request"), dict) else {}
    seq = args.get("sequence")
    if not isinstance(seq, str) or "max_query_bp" in (request.get("bounds") or {}):
        return False
    mode = "discover" if args.get("discover") else "labels"
    return len(seq) <= SERVICE_RESOLVE_MAX_QUERY_BP[mode]


def advanced_keys(arguments):
    """The traverse_fetch arguments on the advanced strategy path, as dotted names ([] for none)."""
    out = []
    ov = (arguments or {}).get("overrides")
    if isinstance(ov, dict):
        for key, value in ov.items():
            allowed = MCP_BASIC_OVERRIDES.get(key, "absent")
            if allowed is None:
                continue
            if allowed == "absent":
                out.append("overrides.%s" % key)
                continue
            for sub in (value if isinstance(value, dict) else {}):
                if sub not in allowed:
                    out.append("overrides.%s.%s" % (key, sub))
    if (arguments or {}).get("preset") not in MCP_BASIC_PRESETS:
        out.append("preset %s" % arguments.get("preset"))
    if "strategy" in (arguments or {}):
        out.append("strategy")
    return out


def explaining_answered_cases(scenario, cases):
    """The answered / clamped cases the scenario names that explain checks (compare verifies their pair only
    against a direct execution)."""
    out = []
    for cid in cmp.case_ids(scenario):
        case = cases.get(cid) or {}
        if cmp.judged_by_checks(case) and cmp.explains(case):
            out.append(cid)
    return out


def mcp_exclusion(scenario, cases, advanced_path=False, direct_arm=True):
    """The MCP arm's exclusion of a scenario ({scenario, arms, rules, reason}) or None; |direct_arm|: the policy runs
    a direct arm (else needs_direct_arm applies)."""
    mcp = scenario.get("mcp") or {}
    rules, reasons = [], []
    if mcp.get("not_implemented"):
        rules.append("manifest_not_implemented")
        reasons.append("the manifest marks it not implemented via MCP (%s)" % mcp["not_implemented"])
    for cid in cmp.case_ids(scenario):
        case = cases.get(cid) or {}
        if case.get("agreement") == "not_implemented_mcp" and not not_expressible_moot(case, scenario):
            rules.append("table_not_expressible")
            reasons.append("it names %s, which the by-design table states not expressible via MCP" % cid)
    if mcp.get("tool") == "traverse_fetch" and not advanced_path:
        keys = advanced_keys(mcp.get("arguments"))
        if keys:
            rules.append("advanced_path")
            reasons.append("its traverse_fetch arguments need the advanced strategy path (%s), which the deployment "
                           "does not enable" % ", ".join(keys))
    answered = [] if direct_arm else explaining_answered_cases(scenario, cases)
    if answered:
        rules.append("needs_direct_arm")
        reasons.append("it names %s, answered by the service where the host refuses: compare verifies that pair only "
                       "against a direct execution, which this policy does not run" % ", ".join(answered))
    if not rules:
        return None
    return {"scenario": scenario.get("id"), "arms": ["mcp"], "rules": sorted(set(rules), key=rules.index),
            "reason": "; ".join(reasons)}


DECLARED_VERSIONS = (1,)


def check_declared(declared, manifest, selected):
    """The owner-declared entries (module docstring), validated: a ValueError names every problem."""
    problems = []
    if not isinstance(declared, dict):
        raise ValueError("the declared file is not a JSON object")
    if declared.get("declared_version") not in DECLARED_VERSIONS:
        problems.append("declared_version %r is not one of %s" % (declared.get("declared_version"), DECLARED_VERSIONS))
    want = (declared.get("manifest") or {}).get("content_sha256")
    if want != sc.manifest_identity(manifest)["content_sha256"]:
        problems.append("it was written for the manifest content_sha256 %s, not this one" % want)
    scen = dict((s.get("id"), s) for s in manifest.get("scenarios", []))
    rels = set(r.get("id") for r in manifest.get("relations", []) or [])
    for kind in ("exclusions", "not_required", "relations_excluded"):
        for i, e in enumerate(declared.get(kind) or []):
            where = "%s[%d]" % (kind, i)
            if not isinstance(e, dict):
                problems.append("%s is not an object" % where)
                continue
            if not str(e.get("reason") or "").strip():
                problems.append("%s states no reason" % where)
            if kind == "relations_excluded":
                if e.get("relation") not in rels:
                    problems.append("%s: %s is not a relation of the manifest" % (where, e.get("relation")))
            elif e.get("scenario") not in selected:
                problems.append("%s: %s is not a selected scenario" % (where, e.get("scenario")))
            reps = e.get("repeats")
            if kind == "exclusions" and reps is not None and (
                    not isinstance(reps, list) or not reps
                    or any(not isinstance(r, int) or isinstance(r, bool) or r < 0 for r in reps)):
                problems.append("%s: repeats %r (a list of repeat indices)" % (where, reps))
            arms = e.get("arms")
            if not isinstance(arms, list) or not arms or any(a not in sc.ARMS for a in arms):
                problems.append("%s: arms %r (one or more of %s)" % (where, arms, ", ".join(sc.ARMS)))
            if kind == "not_required":
                expected = [x.get("check") for x in (scen.get(e.get("scenario")) or {}).get("expect") or []
                            if isinstance(x, dict)]
                if e.get("check") not in expected:
                    problems.append("%s: %s expects no check %r" % (where, e.get("scenario"), e.get("check")))
    if problems:
        raise ValueError("; ".join(problems))


def build_policy(manifest, routes=DEFAULT_ROUTES, arms=("mcp",), table_path=cmp.DEFAULT_BY_DESIGN,
                 advanced_path=False, warm_every=10, require_agreement=False, not_required=None, name=None,
                 manifest_file=None, table=None, declared=None, declared_file=None):
    """The policy (module docstring); |table| (cases by id, info) in place of reading |table_path|; |declared| the
    owner's entries (check_declared), |not_required| more property exemptions given directly."""
    cases, info = table if table is not None else cmp.load_by_design(table_path)
    routes = list(routes)
    arms = [a for a in sc.ARMS if a in arms]
    direct_arm = any(a in sc.DIRECT_ARMS for a in arms)
    chosen = [s.get("id") for s in manifest.get("scenarios", []) if s.get("route") in routes]
    if declared is not None:
        check_declared(declared, manifest, set(chosen))
    declared_ex = {}
    for e in (declared or {}).get("exclusions") or []:
        declared_ex.setdefault(e["scenario"], []).append(e)
    selected, exclusions, required = [], [], []
    for index, scen in enumerate(manifest.get("scenarios", [])):
        if scen.get("route") not in routes:
            continue
        sid = scen.get("id")
        selected.append(sid)
        ex = mcp_exclusion(scen, cases, advanced_path, direct_arm) if "mcp" in arms else None
        if ex:
            exclusions.append(ex)
        for d in declared_ex.get(sid) or []:
            on = [a for a in arms if a in d["arms"] and not (a == "mcp" and ex)]
            if on:
                entry = {"scenario": sid, "arms": on, "rules": ["declared"], "reason": d["reason"]}
                if d.get("repeats") is not None:
                    entry["repeats"] = sorted(set(d["repeats"]))
                exclusions.append(entry)
        for arm in arms:
            mine = [e for e in exclusions if e["scenario"] == sid and arm in e["arms"]]
            if any("repeats" not in e for e in mine):
                continue
            off = set(r for e in mine for r in e.get("repeats") or [])
            reps = [r for _, _, r in gate.planned_executions(manifest, None, None, [sid], warm_every, [arm])
                    if r not in off]
            if reps:
                required.append({"scenario": sid, "arm": arm, "repeats": reps})
    excluded_on = {}
    for ex in exclusions:
        if "repeats" not in ex:   # a repeat-level exclusion leaves the scenario's first repeat to the relations
            excluded_on.setdefault(ex["scenario"], set()).update(ex["arms"])
    declared_rel = {}
    for e in (declared or {}).get("relations_excluded") or []:
        declared_rel.setdefault(e["relation"], []).append(e)
    chosen = set(selected)
    relations, rel_excluded = [], []
    for rel in manifest.get("relations", []) or []:
        a, b = rel.get("a"), rel.get("b")
        outside = [x for x in (a, b) if x not in chosen]
        if outside:
            rel_excluded.append({"relation": rel.get("id"), "arms": arms,
                                 "reason": "%s not in the selection" % ", ".join(outside)})
            continue
        on = [arm for arm in arms if not any(arm in excluded_on.get(x, ()) for x in (a, b))]
        off = [arm for arm in arms if arm not in on]
        if off:
            rel_excluded.append({"relation": rel.get("id"), "arms": off, "reason": "%s excluded on %s" % (
                ", ".join(x for x in (a, b) if set(off) & excluded_on.get(x, set())), ", ".join(off))})
        for d in declared_rel.get(rel.get("id")) or []:
            gone = [arm for arm in on if arm in d["arms"]]
            if gone:
                on = [arm for arm in on if arm not in gone]
                rel_excluded.append({"relation": rel.get("id"), "arms": gone, "rules": ["declared"],
                                     "reason": d["reason"]})
        if on:
            relations.append({"relation": rel.get("id"), "type": rel.get("type"), "a": a, "b": b, "arms": on})
    index = manifest.get("index") or {}
    ident = sc.manifest_identity(manifest)
    policy = {
        "policy": "suite deploy gate (gate.py)", "policy_version": 1,
        "name": name or "%s-%s" % (os.path.splitext(os.path.basename(manifest_file or "manifest"))[0],
                                   "-".join(arms)),
        "generated_utc": sc.iso_utc(time.time()),
        "generated_by": "gate_policy.py",
        "manifest": {"file": manifest_file, "manifest_version": manifest.get("manifest_version"),
                     "content_sha256": ident["content_sha256"], "relations_sha256": ident["relations_sha256"],
                     "index_fp": index.get("index_fp")},
        "by_design_table": dict((k, (info or {}).get(k)) for k in ("sha256", "version", "service_commit", "cases"))
        if info else None,
        "selection": {"rule": "route", "routes": routes, "scenarios": selected},
        "arms": arms, "warm_every": warm_every, "require_agreement": bool(require_agreement),
        "advanced_path": bool(advanced_path),
        "exclusions": exclusions,
        "required": required,
        "relations": relations, "relations_excluded": rel_excluded,
        "properties": {"required": "every property check the scenario expects (checks.is_property)",
                       "not_required": [dict(e, rules=["declared"]) for e in (declared or {}).get("not_required") or []]
                       + list(not_required or [])},
        "declared": None if declared is None else {"file": declared_file,
                                                   "declared_version": declared.get("declared_version"),
                                                   "sha256": hashlib.sha256(sc.canonical_json(declared)).hexdigest()},
    }
    policy["counts"] = counts(manifest, policy)
    return policy


def counts(manifest, policy):
    scen = dict((s.get("id"), s) for s in manifest.get("scenarios", []))
    by_route, by_area = {}, {}
    for sid in policy["selection"]["scenarios"]:
        s = scen.get(sid) or {}
        by_route[s.get("route")] = by_route.get(s.get("route"), 0) + 1
        by_area[s.get("area")] = by_area.get(s.get("area"), 0) + 1
    rules = {}
    for ex in policy["exclusions"]:
        for r in ex["rules"]:
            rules[r] = rules.get(r, 0) + 1
    per_arm = {}
    for item in policy["required"]:
        a = per_arm.setdefault(item["arm"], {"scenarios": 0, "executions": 0})
        a["scenarios"] += 1
        a["executions"] += len(item["repeats"])
    return {"selected": len(policy["selection"]["scenarios"]), "by_route": by_route, "by_area": by_area,
            "excluded": len(policy["exclusions"]), "exclusions_by_rule": rules, "required": per_arm,
            "relations_required": len(policy["relations"]), "relations_excluded": len(policy["relations_excluded"]),
            "properties_not_required": len(policy["properties"]["not_required"])}


def describe(policy):
    c = policy["counts"]
    out = ["policy %s (manifest %s, content %s, relations %s)" % (
        policy["name"], policy["manifest"].get("file"), policy["manifest"]["content_sha256"][:12],
        policy["manifest"]["relations_sha256"][:12])]
    out.append("  selected by route: %d (%s); by area %s" % (
        c["selected"], ", ".join("%s %d" % kv for kv in sorted(c["by_route"].items())),
        ", ".join("%s %d" % kv for kv in sorted(c["by_area"].items()))))
    out.append("  excluded: %d (%s)" % (
        c["excluded"], ", ".join("%s %d" % kv for kv in sorted(c["exclusions_by_rule"].items())) or "none"))
    for ex in policy["exclusions"]:
        if "declared" in ex["rules"]:
            out.append("    declared: %s%s excluded on %s: %s" % (
                ex["scenario"], "".join(" r%d" % r for r in ex.get("repeats") or []), ", ".join(ex["arms"]),
                ex["reason"]))
        else:
            out.append("    %s on %s: %s" % (ex["scenario"], ", ".join(ex["arms"]), ex["reason"]))
    for arm, v in sorted(c["required"].items()):
        out.append("  required on %s: %d scenario(s), %d execution(s)" % (arm, v["scenarios"], v["executions"]))
    out.append("  relations: %d required, %d excluded" % (c["relations_required"], c["relations_excluded"]))
    for r in policy["relations"]:
        out.append("    required %s (%s %s, %s) on %s" % (r["relation"], r["type"], r["a"], r["b"], ", ".join(r["arms"])))
    for r in policy["relations_excluded"]:
        out.append("    %s %s excluded on %s: %s" % ("declared:" if "declared" in (r.get("rules") or ()) else "",
                                                    r["relation"], ", ".join(r["arms"]), r["reason"]))
    out.append("  property checks not required: %d" % c["properties_not_required"])
    for nr in policy["properties"]["not_required"]:
        out.append("    %s%s %s not required on %s: %s" % ("declared: " if "declared" in (nr.get("rules") or ()) else "",
                                                       nr.get("scenario"), nr.get("check"),
                                                       ", ".join(nr.get("arms") or sc.ARMS), nr.get("reason")))
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("manifest")
    ap.add_argument("--routes", default=",".join(DEFAULT_ROUTES))
    ap.add_argument("--arms", default="mcp")
    ap.add_argument("--by-design", default=cmp.DEFAULT_BY_DESIGN)
    ap.add_argument("--advanced-path", action="store_true",
                    help="the deployment enables the advanced strategy path: no advanced_path exclusion")
    ap.add_argument("--warm-every", type=int, default=10)
    ap.add_argument("--require-agreement", action="store_true",
                    help="the MCP arm must agree with a direct arm on every required property")
    ap.add_argument("--declared", default=None, help="the owner's declared exclusions and not-required properties")
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", default=None, help="write the policy here (default: print the counts only)")
    ap.add_argument("--ids-out", default=None, help="write the ids the policy's arms run here, one per line")
    a = ap.parse_args(argv)
    manifest = sc.read_json(a.manifest)
    if not isinstance(manifest, dict):
        sys.stderr.write("gate_policy.py: %s is not a manifest\n" % a.manifest)
        return 2
    problems = sc.manifest_problems(manifest)
    if problems:
        sys.stderr.write("gate_policy.py: %s does not verify: %s\n" % (a.manifest, "; ".join(problems)))
        return 2
    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    bad = [x for x in arms if x not in sc.ARMS]
    if bad or not arms:
        sys.stderr.write("gate_policy.py: --arms takes %s\n" % ", ".join(sc.ARMS))
        return 2
    declared = None
    if a.declared:
        declared = sc.read_json(a.declared)
        if declared is None:
            sys.stderr.write("gate_policy.py: --declared %s is not a JSON file\n" % a.declared)
            return 2
    rel = os.path.relpath(os.path.abspath(a.manifest), HERE)
    try:
        policy = build_policy(manifest, [x.strip() for x in a.routes.split(",") if x.strip()], arms, a.by_design,
                              a.advanced_path, a.warm_every, a.require_agreement, None, a.name,
                              rel if not rel.startswith("..") else os.path.abspath(a.manifest),
                              declared=declared, declared_file=a.declared)
    except ValueError as e:
        sys.stderr.write("gate_policy.py: --declared %s: %s\n" % (a.declared, e))
        return 2
    sys.stdout.write(describe(policy))
    if a.out:
        sc.write_json_atomic(a.out, policy)
        with open(a.out, "rb") as fh:
            sys.stdout.write("wrote %s (sha256 %s)\n" % (a.out, hashlib.sha256(fh.read()).hexdigest()[:12]))
    if a.ids_out:
        ids = gate.policy_ids(policy)
        with open(a.ids_out, "w") as fh:
            fh.write("# %s: the scenarios the arms %s run (%d), by route less the exclusions\n" % (
                policy["name"], ", ".join(policy["arms"]), len(ids)))
            fh.write("".join("%s\n" % i for i in ids))
        sys.stdout.write("wrote %s (%d id(s))\n" % (a.ids_out, len(ids)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
