#!/usr/bin/env python3
"""Writes a deploy-gate policy for gate.py from a manifest, and prints its counts.

    gate_policy.py MANIFEST [--routes R1,R2] [--arms mcp] [--by-design PATH] [--advanced-path]
                   [--warm-every 10] [--require-agreement] [--not-required FILE] [--name NAME] [--out POLICY]

The selection is by ROUTE (default: POST /resolve, POST /traverse and GET /traverse/capabilities: every scenario the
resolve / traverse service answers, the N and D scenarios on those routes included, not only areas R and T). For the
MCP arm the scenarios it cannot express are explicit exclusions, each derived from the manifest and the by-design
table, never from an answer a run received:
  manifest_not_implemented  the manifest's mcp.not_implemented (no MCP call maps the request);
  table_not_expressible     a case the scenario names (mcp.by_design) whose table agreement is not_implemented_mcp
                            (e.g. resolve.not_expressible.bounds_max_query_bp: the MCP sends another request);
  advanced_path             traverse_fetch arguments on the advanced strategy path (an `overrides` key outside the
                            MCP's basic list, a preset other than walk / trie, a `strategy` object), which a deployment
                            without TRAVERSAL_ADVANCED_ENABLED refuses; --advanced-path for one that enables it.
A relation is required on an arm when both of its scenarios are selected and neither is excluded there; the others are
listed with the reason. No property check is exempt unless --not-required FILE (a JSON list of {scenario (an id or
"*"), arms, check, reason}) says so. The policy pins the manifest's content_sha256, relations_sha256 and index_fp and
the by-design table's sha256 (gate.py fails on another), and lists the executions each arm must record (the paired
plan: suite_common.mcp_repeats for the MCP arm, plan_executions for a direct arm, with --warm-every).
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


def mcp_exclusion(scenario, cases, advanced_path=False):
    """The MCP arm's exclusion of a scenario ({scenario, arms, rules, reason}) or None."""
    mcp = scenario.get("mcp") or {}
    rules, reasons = [], []
    if mcp.get("not_implemented"):
        rules.append("manifest_not_implemented")
        reasons.append("the manifest marks it not implemented via MCP (%s)" % mcp["not_implemented"])
    for cid in cmp.case_ids(scenario):
        case = cases.get(cid) or {}
        if case.get("agreement") == "not_implemented_mcp":
            rules.append("table_not_expressible")
            reasons.append("it names %s, which the by-design table states not expressible via MCP" % cid)
    if mcp.get("tool") == "traverse_fetch" and not advanced_path:
        keys = advanced_keys(mcp.get("arguments"))
        if keys:
            rules.append("advanced_path")
            reasons.append("its traverse_fetch arguments need the advanced strategy path (%s), which the deployment "
                           "does not enable" % ", ".join(keys))
    if not rules:
        return None
    return {"scenario": scenario.get("id"), "arms": ["mcp"], "rules": sorted(set(rules), key=rules.index),
            "reason": "; ".join(reasons)}


def build_policy(manifest, routes=DEFAULT_ROUTES, arms=("mcp",), table_path=cmp.DEFAULT_BY_DESIGN,
                 advanced_path=False, warm_every=10, require_agreement=False, not_required=None, name=None,
                 manifest_file=None, table=None):
    """The policy (module docstring); |table| (cases by id, info) in place of reading |table_path|."""
    cases, info = table if table is not None else cmp.load_by_design(table_path)
    routes = list(routes)
    arms = [a for a in sc.ARMS if a in arms]
    selected, exclusions, required = [], [], []
    for index, scen in enumerate(manifest.get("scenarios", [])):
        if scen.get("route") not in routes:
            continue
        sid = scen.get("id")
        selected.append(sid)
        ex = mcp_exclusion(scen, cases, advanced_path) if "mcp" in arms else None
        if ex:
            exclusions.append(ex)
        for arm in arms:
            if arm == "mcp" and ex:
                continue
            reps = gate.planned_executions(manifest, None, None, [sid], warm_every, [arm])
            if reps:
                required.append({"scenario": sid, "arm": arm, "repeats": [r for _, _, r in reps]})
    excluded_on = dict((ex["scenario"], set(ex["arms"])) for ex in exclusions)
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
        if on:
            relations.append({"relation": rel.get("id"), "type": rel.get("type"), "a": a, "b": b, "arms": on})
        if off:
            rel_excluded.append({"relation": rel.get("id"), "arms": off, "reason": "%s excluded on %s" % (
                ", ".join(x for x in (a, b) if set(off) & excluded_on.get(x, set())), ", ".join(off))})
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
                       "not_required": list(not_required or [])},
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
    out.append("  excluded on the MCP arm: %d (%s)" % (
        c["excluded"], ", ".join("%s %d" % kv for kv in sorted(c["exclusions_by_rule"].items())) or "none"))
    for ex in policy["exclusions"]:
        out.append("    %s: %s" % (ex["scenario"], ex["reason"]))
    for arm, v in sorted(c["required"].items()):
        out.append("  required on %s: %d scenario(s), %d execution(s)" % (arm, v["scenarios"], v["executions"]))
    out.append("  relations: %d required, %d excluded" % (c["relations_required"], c["relations_excluded"]))
    for r in policy["relations"]:
        out.append("    required %s (%s %s, %s) on %s" % (r["relation"], r["type"], r["a"], r["b"], ", ".join(r["arms"])))
    for r in policy["relations_excluded"]:
        out.append("    excluded %s on %s: %s" % (r["relation"], ", ".join(r["arms"]), r["reason"]))
    out.append("  property checks not required: %d" % c["properties_not_required"])
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
    ap.add_argument("--not-required", default=None, help="a JSON list of {scenario, arms, check, reason}")
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", default=None, help="write the policy here (default: print the counts only)")
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
    not_required = sc.read_json(a.not_required) if a.not_required else []
    if not isinstance(not_required, list):
        sys.stderr.write("gate_policy.py: --not-required is not a JSON list\n")
        return 2
    rel = os.path.relpath(os.path.abspath(a.manifest), HERE)
    policy = build_policy(manifest, [x.strip() for x in a.routes.split(",") if x.strip()], arms, a.by_design,
                          a.advanced_path, a.warm_every, a.require_agreement, not_required, a.name,
                          rel if not rel.startswith("..") else os.path.abspath(a.manifest))
    sys.stdout.write(describe(policy))
    if a.out:
        sc.write_json_atomic(a.out, policy)
        with open(a.out, "rb") as fh:
            sys.stdout.write("wrote %s (sha256 %s)\n" % (a.out, hashlib.sha256(fh.read()).hexdigest()[:12]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
