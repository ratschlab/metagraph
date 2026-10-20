#!/usr/bin/env python3
"""The area G rehearsal (DESIGN-v8 §3, §4, §8 phase P2; build_manifest.py --rehearse G): every candidate graphlet of
area_g.candidates run on staging, judged, and the corpus file the v8 build reads written from what it showed.

    rehearse_g.py --out DIR [--plan]                    the rehearsal manifest and the plan (no request): default
    rehearse_g.py --out DIR --execute [--mcp-cmd T]     MCP cold once per candidate, then direct-local x3 (+ twin)
    rehearse_g.py --out DIR --analyse --dl-run RUN --mcp-dir MCP [--library-root LIB] [--mcp-graphlet-run R]
    options: --library-commit SHA (the pin written into the corpus), --skip-conditional (no over-8-MiB candidate),
             --slots a,b (only these slots and their pairs)

Per candidate (DESIGN-v8 §3 and the lane D brief):
  MCP           one traverse_fetch (preset + overrides) through the async session's runner, cold, BEFORE the direct
                arm (the baseline's paired order): the job's graphlet body, traverse_status.overrides (the echo of the
                strategy the service sent) and the engine's elapsed time;
  direct-local  the `detail: graphlet` request three times (repeat 0, 1, 2) and the `detail: full` twin once
                (run_direct.py, area G: arms/direct-local/docs/<id>.full.json.gz);
  verdicts      echo (traverse_status.overrides == request.strategy, numbers by value), bodies (graphlet_doc
                body_digest equal over the three direct answers and the MCP's), envelope (the engine envelope without
                its per-request values, graphlet_doc.envelope_view, equal), twin (graphlet_precheck.twin_difference:
                the graphlet's arms are the full twin's), headroom (every engine time <= 25 % of the effective time
                budget, min(the strategy's time_budget_ms, the host's 60,000)), slot (the shape the slot needs: merges,
                cut lists, a resource stop by work, complete_to_bp 0, one linear walk, no graphlet, a walk with its seed
                over 16 kbp, a body over 8 MiB), plan (the op plan run by the local library on the direct document,
                W1, every op answering what area_g.expected_outcome says; selector flags the plan reads adjusted and the
                plan run again once); fit (every whole list W3 runs at max_bytes 262,144 fits it uncut in W1's
                answer, with room for W3's longer handle and its own fields: DESIGN-v8 §5);
  facts         from the document by the suite's own reader (checks.graphlet_arms): arms, walks, segments, merges,
                labels, walk lengths, bases, cut lists, labels per node, body bytes, continuations (the twin's paths);
                frozen as graphlet.preconditions; the service-check facts as graphlet.facts;
  timeouts      per op: min(lane, max(TIMEOUT_FLOOR_S, TIMEOUT_FACTOR x the median elapsed)) over the local runs and
                W3's (--mcp-graphlet-run, the async runner's area G records), written only when W3's are measured
                (else the lane: a local median says nothing of an MCP call); the medians always recorded; an op whose
                median exceeds a third of its lane is named (DESIGN-v8 §3);
  W3 headroom   (DESIGN-v8 §3, refined 2026-10-20) every op W3 runs within a third of its service lane: interactive
                ops a third of lanes.interactive.alarm_s (3 s -> 1 s), heavy ops (lanes.heavy_ops) a third of
                lanes.heavy.timeout_s (120 s -> 40 s); an interactive op whose record states it escalated (never
                inferred from its time: VERIFY-D VD10): allowed, recorded (escalations), judged as heavy. The values
                from the run's recorded capabilities, service_settings.graphlet.lanes (graphlet_lanes: a part absent or null is the
                fallback, recorded missing); a miss fails op_headroom.
A candidate that fails a verdict leaves the corpus (DESIGN-v8 §3: no op_timeout, no time stop is an expected outcome).
Per slot the first candidate (area_g order) that passes every verdict is chosen; the pairs follow their slot's chosen
candidate (its partner must pass too); the three chosen graphlets with the most segments are the richest. A
direct-only candidate (a strategy the MCP cannot express) is judged on the direct arm alone and enters the corpus only
with the owner's decision (--include-direct-only at the build).

Writes DIR/rehearsal_manifest.json, DIR/plan.json and DIR/PLAN.md (before), DIR/rehearsal.json, DIR/corpus.json and
DIR/REPORT.md (the analysis). Default DIR: ~/.cache/metagraph-traverse-plans/<UTC date>/v8/rehearsal/<UTC>/.
"""

from __future__ import annotations

import argparse
import copy
import datetime
import glob
import gzip
import json
import os
import re
import shlex
import statistics
import subprocess
import sys
from typing import Any, Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import area_g  # noqa: E402
import checks  # noqa: E402
import graphlet_doc  # noqa: E402
import graphlet_precheck as gp  # noqa: E402
import manifest_lib as ml  # noqa: E402
import suite_common as sc  # noqa: E402

V7 = os.path.join(HERE, "manifest", "suite_v7.json")
DL_REPEATS = 3
TIMEOUT_FLOOR_S = 5.0
TIMEOUT_FACTOR = 10.0
SEQUENCE_CEILING = 16384           # graphlet_sequence's fixed ceiling (the service's and the local ways')
IN_SERVICE_BYTES = 8 * 1024 * 1024
LIB_DEFAULT = os.path.normpath(os.path.join(HERE, "..", "..", "api", "python"))
# the engine envelope's per-request values (graphlet_doc.envelope_view drops them from a J line)
RESULT_BODY_KEYS = ("graphlet", "graphlet_bytes", "graphlet_lines")
# selector flags of the op plan the shadow run may turn off (op id stem -> facts key)
FLAG_OF_OP = (("walk-merge-", "merge_walk"), ("walk-cont-", "continuation"), ("labels-entered-", "entered"),
              ("walks-entered-", "entered"))


def utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def default_out() -> str:
    now = datetime.datetime.now(datetime.timezone.utc)
    return os.path.expanduser(os.path.join("~", ".cache", "metagraph-traverse-plans", now.strftime("%Y-%m-%d"), "v8",
                                           "rehearsal", now.strftime("%Y%m%dT%H%M%SZ")))


# ===================================================================================================== candidates
def build_context():
    """(ctx, v7 scenarios): the builder's context over the record cache (None without it: the long-walk candidates,
    drawn from the records, are left out) and the frozen v7 scenarios whose T seeds the candidates reuse."""
    import truth
    v7 = sc.read_json(V7)["scenarios"]
    try:
        recs = truth.load_with_sequences(fetch=False)
        if any(r.fetch_error for r in recs):
            raise ValueError("records missing from the cache")
        ctx = ml.Ctx(recs, 20261013)
    except Exception as e:   # no cache: a context without records
        ctx = _NoRecords(str(e))
    return ctx, v7


class _NoRecords(object):
    def __init__(self, why: str):
        self.why = why
        self.skipped: List[dict] = []

    def rec(self, key):
        raise KeyError("%s (no record cache: %s)" % (key, self.why))

    def rng(self, *parts):
        import random
        return random.Random(0)

    def skip(self, area, what, reason):
        self.skipped.append({"area": area, "scenario": what, "reason": reason})


def candidate_ids(cands: List[dict]) -> Dict[str, str]:
    """{cid: rehearsal scenario id}: G-<stem>-c<n> (the slot's stem, the candidate's number)."""
    stems = dict((s[0], s[1]) for s in area_g.SLOTS)
    out = {}
    for c in cands:
        n = c["cid"].split(".")[1]
        out[c["cid"]] = "%s-c%s" % (stems[c["slot"]], n.replace("+", "-"))
    return out


def partner_candidates(cands: List[dict]) -> List[Tuple[dict, dict, str]]:
    """[(partner candidate, its corpus candidate, pair)] for every candidate of a slot that has pairs (area_g.PAIRS):
    the partner the build would make if that candidate were chosen."""
    out = []
    ids = dict((s[0], "%s-0001" % s[1]) for s in area_g.SLOTS)
    for c in cands:
        specs = area_g.pair_plan({c["slot"]: c}, ids, {})
        for ps in specs:
            if ps.get("partner_candidate") is not None:
                out.append((ps["partner_candidate"], c, ps["pair"]))
    return out


def rehearsal_scenarios(ctx: Any, v7: List[dict], skip_conditional: bool = False,
                        slots: Optional[List[str]] = None) -> Tuple[List[dict], Dict[str, dict]]:
    """(the rehearsal manifest's scenarios, {scenario id: {cid, slot, pair?, of?, candidate}}): every candidate as an
    area G scenario (area_g.g_scenario with its expected facts: the plan the shadow run tests), repeat 3 on the
    direct arm, each pair's partner right after its corpus candidate."""
    cands = [c for c in area_g.candidates(ctx, v7)
             if (not slots or c["slot"] in slots) and not (skip_conditional and c["conditional"])]
    ids = candidate_ids(cands)
    partners = partner_candidates(cands)
    out, meta = [], {}
    corpus = area_g.default_corpus()
    for c in cands:
        sid = ids[c["cid"]]
        facts = area_g.facts_of(c, corpus)
        rich = c["slot"] in area_g.DEFAULT_RICH
        ops = area_g.op_plan(sid, c["slot"], facts, rich, [(sid, False)] if c["slot"] == "two_arms" else [],
                             direct_only=not area_g.expressible(c),
                             extras={"other": ids[cands[0]["cid"]] if c["expected"].get("failed") else sid})
        s = area_g.g_scenario(sid, "rehearsal: %s candidate %s" % (c["slot"], c["cid"]), c, corpus, None, ops,
                              ["rehearsal"])
        s["repeat"] = DL_REPEATS
        out.append(s)
        meta[sid] = {"cid": c["cid"], "slot": c["slot"], "candidate": c}
        for pc, of, pair in partners:
            if of["cid"] != c["cid"]:
                continue
            pid = "G-pair-%s-c%s" % (pair.replace("_", ""), c["cid"].split(".")[1])
            pfacts = area_g.facts_of(pc, corpus)
            pops = area_g.op_plan(pid, pc["slot"], pfacts, False, [(sid, pair in ("radius", "two_seeds"))],
                                  extras={}, partner_of=sid)
            ps = area_g.g_scenario(pid, "rehearsal: pair %s of %s" % (pair, c["cid"]), pc, corpus, None, pops,
                                   ["rehearsal", "pair"])
            ps["repeat"] = DL_REPEATS
            out.append(ps)
            meta[pid] = {"cid": pc["cid"], "slot": c["slot"], "pair": pair, "of": c["cid"], "candidate": pc}
    return out, meta


# The second rehearsal (the coordinator, 2026-10-19), short and direct-local mostly: the twins without
# max_memory_mb (F-D4) of the candidates whose twin failed (SPEC-core R21 (4)), direct-local only, one graphlet
# request and the twin each; and the long walk with its warm-ups (the owner's decision), on fresh windows (the first
# rehearsal walked the old ones), its MCP call cold after the warm-ups, three direct repeats and the twin.
REHEARSAL2_TWIN = ("two_arms.1", "merges.1", "cut_lists.1", "cut_lists.2", "cut_lists.3", "no_bases.1", "no_bases.2")
REHEARSAL2_TWIN_PARTNERS = ("two_arms.1+shift10", "merges.1+cap2")
REHEARSAL2_LONG = ("long_walk.4", "long_walk.5")
REHEARSAL2_DIRECT = ("rehearsal 2 (the coordinator, 2026-10-19): direct-local only, the twin's re-check without "
                     "bounds.max_memory_mb (F-D4)")


def rehearsal2_scenarios(ctx: Any, v7: List[dict]) -> Tuple[List[dict], Dict[str, dict]]:
    """(the second rehearsal's scenarios, meta): the twin re-checks (REHEARSAL2_TWIN and their partners: repeat 1,
    the MCP mapping not implemented), then each long walk of REHEARSAL2_LONG right after its two warm-ups
    (area_g.long_walk_warm_ups: direct-local only, the long class) with its MCP mapping and repeat 3. Ids and op plans
    as in the first rehearsal (a candidate keeps its G-<stem>-c<n>)."""
    cands = area_g.candidates(ctx, v7)
    by = dict((c["cid"], c) for c in cands)
    ids = candidate_ids(cands)
    corpus = area_g.default_corpus()
    out, meta = [], {}
    partners = dict((pc["cid"], (pc, of, pair)) for pc, of, pair in partner_candidates(cands))
    for cid in REHEARSAL2_TWIN:
        c = by[cid]
        sid = ids[cid]
        ops = area_g.op_plan(sid, c["slot"], area_g.facts_of(c, corpus), c["slot"] in area_g.DEFAULT_RICH,
                             [(sid, False)] if c["slot"] == "two_arms" else [], direct_only=True,
                             extras={"other": sid})
        sc_ = area_g.g_scenario(sid, "rehearsal 2: the twin of %s" % cid, c, corpus, None, ops,
                                ["rehearsal", "rehearsal-2"])
        _direct_local_only(sc_)
        out.append(sc_)
        meta[sid] = {"cid": cid, "slot": c["slot"], "candidate": c}
        for pcid in REHEARSAL2_TWIN_PARTNERS:
            pc, of, pair = partners[pcid]
            if of["cid"] != cid:
                continue
            pid = "G-pair-%s-c%s" % (pair.replace("_", ""), cid.split(".")[1])
            pops = area_g.op_plan(pid, pc["slot"], area_g.facts_of(pc, corpus), False,
                                  [(sid, pair in ("radius", "two_seeds"))], extras={}, partner_of=sid)
            ps = area_g.g_scenario(pid, "rehearsal 2: the twin of %s" % pcid, pc, corpus, None, pops,
                                   ["rehearsal", "rehearsal-2", "pair"])
            _direct_local_only(ps)
            out.append(ps)
            meta[pid] = {"cid": pcid, "slot": c["slot"], "pair": pair, "of": cid, "candidate": pc}
    for cid in REHEARSAL2_LONG:
        c = by.get(cid)
        if c is None:
            raise SystemExit("rehearsal 2: no candidate %s (the record cache: SUITE_TRUTH_CACHE)" % cid)
        sid = ids[cid]
        n = cid.split(".")[1]
        warm = area_g.long_walk_warm_ups(ctx, c, sid, ("G-warm-c%s-1" % n, "G-warm-c%s-2" % n))
        for w in warm:
            out.append(w)
            meta[w["id"]] = {"warm_up_of": sid}
        ops = area_g.op_plan(sid, c["slot"], area_g.facts_of(c, corpus), False, [], extras={"other": sid,
                                                                                            "row_cut": True,
                                                                                            "sequence_cut": True})
        sc_ = area_g.g_scenario(sid, "rehearsal 2: %s with its warm-ups" % cid, c, corpus, None, ops,
                                ["rehearsal", "rehearsal-2"], timing_class="slow")
        sc_["repeat"] = DL_REPEATS
        out.append(sc_)
        meta[sid] = {"cid": cid, "slot": c["slot"], "candidate": c}
    return out, meta


def _direct_local_only(scen: dict) -> None:
    """A twin re-check of the second rehearsal: one graphlet request (repeat 1) and its twin, no MCP job (the
    first rehearsal recorded the MCP's), so a direct-only scenario of the rehearsal's manifest."""
    scen["repeat"] = 1
    scen["mcp"] = {"not_implemented": REHEARSAL2_DIRECT}
    g = scen["graphlet"]
    g["ways"]["default"] = list(area_g.DIRECT_ONLY_WAYS)
    g["not_exposed"]["*"] = REHEARSAL2_DIRECT
    for o in g["ops"]:
        o["ways"] = [w for w in o.get("ways") or area_g.DIRECT_ONLY_WAYS if w in area_g.DIRECT_ONLY_WAYS]
        o["not_exposed"] = {"W3": REHEARSAL2_DIRECT, "W4": REHEARSAL2_DIRECT}


def plan2(scenarios: List[dict], meta: Dict[str, dict]) -> dict:
    """What the second rehearsal sends (run_suite.sh --paired --skip-remote --warm-every 0: per scenario in order the
    MCP arm cold, direct-local, the MCP arm's warm repeats), and its expected duration."""
    rows, dl, mcp_jobs, seconds = [], 0, 0, 0.0
    for s_ in scenarios:
        m = meta[s_["id"]]
        if m.get("warm_up_of"):
            n_dl, jobs, t = 1, 0, 20.0   # a cold resolve of ~8.3 kbp under a single-record column (estimate)
        elif "tool" in s_["mcp"]:
            # the paired plan's MCP executions (suite_common.mcp_repeats: an area G scenario runs once on MCP since
            # 2026-10-20; rehearsal 2 of 2026-10-19 still ran its long walks' two MCP warm repeats)
            n_dl, jobs = DL_REPEATS + 1, len(sc.mcp_repeats(0, s_, 0))
            t = jobs * (1.0 + MCP_OVERHEAD_S) + n_dl * (0.5 + DL_PAUSE_S)
        else:
            n_dl, jobs = s_.get("repeat", 1) + 1, 0
            t = n_dl * (0.5 + DL_PAUSE_S)
        dl += n_dl
        mcp_jobs += jobs
        seconds += t
        rows.append({"scenario": s_["id"], "candidate": m.get("cid") or "warm-up of %s" % m["warm_up_of"],
                     "mcp_jobs": jobs, "direct_local_requests": n_dl, "expected_s": round(t, 1)})
    dl += 4
    return {"candidates": sum(1 for s_ in scenarios if not meta[s_["id"]].get("warm_up_of")),
            "warm_ups": sum(1 for s_ in scenarios if meta[s_["id"]].get("warm_up_of")),
            "mcp_jobs": mcp_jobs, "direct_local_requests": dl, "direct_remote_requests": 0,
            "expected_minutes": round(seconds / 60.0, 1), "rows": rows,
            "load": "one request in flight at a time, in manifest order (run_suite.sh --paired --skip-remote): the "
                    "twin re-checks on direct-local only (one graphlet request and its twin each; the MCP runner "
                    "only records their not-implemented mark), then per long walk its two warm-up resolves on "
                    "direct-local, its MCP call cold, three direct-local repeats and the twin (the MCP arm's warm "
                    "repeats only where suite_common.mcp_repeats plans them: none for area G since 2026-10-20)"}


def execute2(out: str, manifest_path: str, mcp_cmd: str) -> int:
    """The second rehearsal in one paired run (staging: only with the coordinator's go)."""
    env = dict(os.environ, SUITE_RUNS=os.path.join(out, "runs"))
    return subprocess.call([os.path.join(HERE, "run_suite.sh"), "--manifest", manifest_path, "--paired",
                            "--skip-remote", "--mcp-cmd", mcp_cmd, "--warm-every", "0", "--no-gate", "--no-pdf"],
                           env=env, stdin=subprocess.DEVNULL)


# ===================================================================================================== rehearsal 3 (S)
# The stop rehearsal (VERIFY-D VD6; VERIFY-D-R2 R2-3, the coordinator 2026-10-20: after async's bundle with the
# release record is deployed), in production's warm state: first the scenarios that warm the host before the stops
# in a real v8 run (STOP_WARMERS: as the manifest has them, MCP and direct-local), then each stop scenario's base
# request once to completion on direct-local (S-base-c<n>: its warm engine time recorded), then each stop scenario
# STOP_TRIES times, interleaved, MCP only (direct_not_run). The entry condition: every try a taken stop (stop.took,
# no race) with every stop check passing. With it, the direct-local rehearsal of the constrain / annotate pair (VD8,
# R2-5: G-annot-c1 and its constrain partner G-pair-annot-c1, three graphlet requests and the twin each), whose
# compare rows and cmp-paged's pages the analysis reads from the built plans.
STOP_TRIES = 3
REHEARSAL3_PAIR = "annotate.1+constrain"
STOP_WARMERS = ("T-walk-0028", "R-warm-0005", "T-labels-0003")
REHEARSAL3_DIRECT = ("rehearsal S (VERIFY-D-R2 R2-3): the base request run once to completion on direct-local, to "
                     "warm the host as a real run does before the stop")


def rehearsal3_scenarios(ctx: Any, v7: List[dict]) -> Tuple[List[dict], Dict[str, dict]]:
    """(the stop rehearsal's scenarios, meta): STOP_WARMERS (v7's scenarios), the base runs S-base-c<n>, try t of
    S-stop-000n as S-stop-c<n>-t<t> (a copy of the built scenario, the tries of the two interleaved), then G-annot-c1
    and G-pair-annot-c1 (direct-local)."""
    import areas
    stops = areas.gen_stop(ctx, v7)
    if len(stops) != len(areas.STOP_IDS):
        raise SystemExit("rehearsal S: the stop scenarios do not build (%s)" % [x["id"] for x in stops])
    by = dict((x["id"], x) for x in v7)
    out, meta = [], {}
    for wid in STOP_WARMERS:
        if wid not in by:
            raise SystemExit("rehearsal S: no %s in v7" % wid)
        x = copy.deepcopy(by[wid])
        x.pop("warms", None)   # each runs on its own here, in this order
        out.append(x)
        meta[wid] = {"warmer": True}
    for n, (route, req) in enumerate(areas.stop_base_requests(v7), 1):
        sid = "S-base-c%d" % n
        x = ml.scenario(sid, "S", "rehearsal S: %s's base request once to completion (warms the host)"
                        % areas.STOP_IDS[n - 1], route, copy.deepcopy(req), {}, [{"check": "status",
                                                                                  "args": {"http": 200}}],
                        timing_class="slow", tags=["rehearsal", "rehearsal-s", "stop-base"],
                        mcp={"not_implemented": REHEARSAL3_DIRECT})
        out.append(x)
        meta[sid] = {"base_of": areas.STOP_IDS[n - 1]}
    for t in range(1, STOP_TRIES + 1):
        for n, base in enumerate(stops, 1):
            x = copy.deepcopy(base)
            x["id"] = "S-stop-c%d-t%d" % (n, t)
            x["tags"] = sorted(set(x["tags"] + ["rehearsal", "rehearsal-s"]))
            out.append(x)
            meta[x["id"]] = {"stop_of": base["id"], "try": t}
    cands = area_g.candidates(ctx, v7)
    partners = dict((pc["cid"], (pc, of, pair)) for pc, of, pair in partner_candidates(cands))
    if REHEARSAL3_PAIR not in partners:
        raise SystemExit("rehearsal S: no partner %s" % REHEARSAL3_PAIR)
    pc, of, pair = partners[REHEARSAL3_PAIR]
    corpus = area_g.default_corpus()
    of_sid = candidate_ids([of])[of["cid"]]
    pid = "G-pair-annot-c1"
    aops = area_g.op_plan(of_sid, of["slot"], area_g.facts_of(of, corpus), False, [], direct_only=True,
                          extras={"other": of_sid})
    annot = area_g.g_scenario(of_sid, "rehearsal S: %s, the annotate side of the pair (R2-5)" % of["cid"], of, corpus,
                              None, aops, ["rehearsal", "rehearsal-s"])
    pops = area_g.op_plan(pid, pc["slot"], area_g.facts_of(pc, corpus), False, [(of_sid, False)],
                          direct_only=True, extras={"no_compare": True}, partner_of=of_sid,
                          reason=area_g.direct_reason(pc))
    partner = area_g.g_scenario(pid, "rehearsal S: the constrain partner of %s (VD8)" % of["cid"], pc, corpus, None,
                                pops, ["rehearsal", "rehearsal-s", "pair"])
    for x, m in ((partner, {"cid": pc["cid"], "slot": pc["slot"], "pair": pair, "of": of["cid"], "candidate": pc}),
                 (annot, {"cid": of["cid"], "slot": of["slot"], "candidate": of})):
        x["repeat"] = DL_REPEATS
        out.append(x)
        meta[x["id"]] = m
    return out, meta


def plan3(scenarios: List[dict], meta: Dict[str, dict]) -> dict:
    """What the stop rehearsal sends, and its expected duration: a stop try is one MCP job (submit, running, delay,
    traverse_stop, terminal, then the admin reads until the release record, at most the lease + 35 s and 90 s after
    the end); the partner three direct-local requests and its twin."""
    rows, mcp_jobs, dl, seconds = [], 0, 0, 0.0
    for x in scenarios:
        m = meta[x["id"]]
        if m.get("stop_of"):
            jobs, n_dl, t = 1, 0, 2.0 + STOP_FLOW_DELAY_S_DEFAULT + 10.0 + 60.0
        elif m.get("warmer"):
            jobs = len(sc.mcp_repeats(0, x, 0))
            n_dl = 0 if x.get("direct_not_run") else int(x.get("repeat") or 1)
            t = jobs * (1.0 + MCP_OVERHEAD_S) + n_dl * (30.0 + DL_PAUSE_S)
        elif m.get("base_of"):
            jobs, n_dl, t = 0, 1, 20.0
        else:
            jobs, n_dl = 0, DL_REPEATS + 1
            t = n_dl * (0.5 + DL_PAUSE_S)
        mcp_jobs += jobs
        dl += n_dl
        seconds += t
        what = ("%s try %d" % (m["stop_of"], m["try"]) if m.get("stop_of") else
                "warms the host as in a v8 run" if m.get("warmer") else
                "%s's base request, once to completion" % m["base_of"] if m.get("base_of") else
                "the pair: %s" % m.get("cid"))
        rows.append({"scenario": x["id"], "what": what, "mcp_jobs": jobs, "direct_local_requests": n_dl,
                     "expected_s": round(t, 1)})
    dl += 4
    return {"stop_tries": sum(1 for x in scenarios if meta[x["id"]].get("stop_of")), "mcp_jobs": mcp_jobs,
            "direct_local_requests": dl, "direct_remote_requests": 0, "expected_minutes": round(seconds / 60.0, 1),
            "rows": rows,
            "load": "one request in flight at a time, in manifest order (run_suite.sh --paired --skip-remote): "
                    "the warmers as in a v8 run (MCP cold, direct-local), each stop's base request once on "
                    "direct-local, then per stop try one MCP job (traverse_fetch of the repeat-free 10 kbp sau "
                    "walk, or traverse_resolve of the 500 bp 16S discovery), stopped %g s after it runs (the direct "
                    "arms send nothing for them: direct_not_run); the pair's two graphlets three direct-local "
                    "requests and a twin each" % STOP_FLOW_DELAY_S_DEFAULT}


STOP_FLOW_DELAY_S_DEFAULT = 0.0


def plan3_markdown(p: dict, out: str, mpath: str, mcp_cmd: str) -> str:
    lines = ["# Stop rehearsal (rehearsal S) plan", "",
             "Stop tries %d (MCP jobs %d), direct-local requests %d, direct-remote 0; expected about %.0f minutes "
             "of staging time (most of it the runner waiting for each attempt's release record)." % (
                 p["stop_tries"], p["mcp_jobs"], p["direct_local_requests"], p["expected_minutes"]), "",
             "Load: %s." % p["load"], "",
             "Precondition: async's bundle with the release record (feat/stop-release-record) deployed on staging, "
             "and the MCP runner at a52bd3d6 or later (the stop flow, the release wait).", "",
             "Command (the coordinator's GO first): `METAGRAPH_ADMIN_KEY_FILE=<key file> rehearse_g.py --rehearsal 3 "
             "--out %s --execute` = `run_suite.sh --manifest %s --paired --skip-remote --mcp-cmd '<%s>' --warm-every 0 "
             "--no-gate --no-pdf` (SUITE_RUNS=%s/runs)" % (out, mpath, mcp_cmd.split(" ")[1] if " " in mcp_cmd
                                                         else mcp_cmd, out), "",
             "Then offline: `rehearse_g.py --rehearsal 3 --out %s --analyse --dl-run <the run> --prior <rehearsal 2's "
             "directory>`: per stop try taken / race and the four stop checks, the entry verdict (3 of 3 taken stops, "
             "every check passing), and the partner's corpus entry merged into a copy of truth/area_g/corpus.json "
             "(OUT/corpus.json) with the built plans re-checked." % out, "",
             "| scenario | what | MCP jobs | direct-local | expected s |", "|---|---|---|---|---|"]
    for r in p["rows"]:
        lines.append("| %s | %s | %d | %d | %.1f |" % (r["scenario"], r["what"], r["mcp_jobs"],
                                                        r["direct_local_requests"], r["expected_s"]))
    return "\n".join(lines) + "\n"


def stop_try_verdicts(manifest: dict, meta: Dict[str, dict], run_dir: str) -> dict:
    """Per stop try of a recorded stop rehearsal: {taken, race, checks {name: (outcome, detail)}, outcome}; per stop
    scenario the entry verdict: every try taken (stop.took, no race) and every check passing."""
    import checks
    tries, by_stop = {}, {}
    for x in manifest["scenarios"]:
        m = meta.get(x["id"]) or {}
        if not m.get("stop_of"):
            continue
        raw = _read(os.path.join(run_dir, "arms", "mcp", sc.raw_name(x["id"], "mcp", 0))) or {}
        ans = raw.get("answer") if isinstance(raw.get("answer"), dict) else None
        body = (ans or {}).get("body") if isinstance((ans or {}).get("body"), dict) else {}
        flow = raw.get("stop_flow") if isinstance(raw.get("stop_flow"), dict) else body.get("stop_flow")
        if ans is None:   # the execution ended without an answer: the raw record's stop_flow alone
            ans = {"http": None, "body": {"stop_flow": flow} if flow else None}
        res = checks.evaluate(x, ans, "mcp")
        stop = (flow or {}).get("stop") or {}
        t = {"taken": bool(stop.get("took")) and not (flow or {}).get("race"), "race": (flow or {}).get("race"),
             "flow_error": (flow or {}).get("error"), "recorded": flow is not None,
             "checks": dict((r["check"], (r["outcome"], r["detail"])) for r in res),
             "outcome": checks.scenario_outcome(res)}
        tries[x["id"]] = t
        by_stop.setdefault(m["stop_of"], []).append(x["id"])
    entry = {}
    for stop_id, ids in sorted(by_stop.items()):
        taken = sum(1 for i in ids if tries[i]["taken"])
        passed = sum(1 for i in ids if tries[i]["outcome"] == "pass")
        entry[stop_id] = {"tries": len(ids), "taken": taken, "passed": passed,
                          "verdict": "pass" if taken == passed == len(ids) == STOP_TRIES else
                          "fail: %d of %d taken, %d passed" % (taken, len(ids), passed)}
    return {"tries": tries, "entry": entry}


def analyse3(out: str, manifest: dict, meta: Dict[str, dict], dl_run: str, library_root: str,
             library_commit: Optional[str], prior: Optional[str]) -> dict:
    """The stop rehearsal's judgement: stop_try_verdicts, the partner analysed as any candidate (analyse), its corpus
    entry merged into a copy of the committed corpus (OUT/corpus.json; the committed file is not touched), the built
    plans of that corpus re-checked on both rehearsals' and this one's documents. Writes OUT/REPORT-S.md and
    OUT/rehearsal_s.json."""
    stops = stop_try_verdicts(manifest, meta, dl_run)
    gmeta = dict((k, v) for k, v in meta.items() if v.get("cid"))
    res = analyse(os.path.join(out, "partner"), {"scenarios": [x for x in manifest["scenarios"] if x["id"] in gmeta]},
                  gmeta, dl_run, os.path.join(dl_run, "arms", "mcp"), library_root, library_commit) \
        if gmeta else {"results": {}, "corpus": {"candidates": {}}}
    committed = sc.read_json(os.path.join(HERE, "truth", "area_g", "corpus.json"))
    merged = copy.deepcopy(committed)
    for cid, entry in (res["corpus"].get("candidates") or {}).items():
        old = (merged.get("candidates") or {}).get(cid) or {}
        entry = copy.deepcopy(entry)
        if "page_bytes" in (old.get("facts") or {}) and "page_bytes" not in (entry.get("facts") or {}):
            entry.setdefault("facts", {})["page_bytes"] = old["facts"]["page_bytes"]   # measured, kept
        merged.setdefault("candidates", {})[cid] = entry
    merged["rehearsal"] = dict(merged.get("rehearsal") or {}, stop_rehearsal=utc_now())
    with open(os.path.join(out, "corpus.json"), "w") as fh:
        json.dump(merged, fh, indent=1, sort_keys=True)
    built = None
    priors = [x for x in (prior or "").split(",") if x]
    if priors:
        # the first prior's rehearsal.json (rehearsal 2's holds rehearsal 1's results merged); the documents of
        # every prior's runs (rehearsal 1's graphlets live in its own run)
        prev = sc.read_json(os.path.join(priors[0], "rehearsal.json")) or {}
        pmeta = dict((prev.get("meta") or {}), **dict((k, dict((kk, vv) for kk, vv in v.items() if kk != "candidate"))
                                                     for k, v in gmeta.items()))
        presults = dict((prev.get("results") or {}), **res["results"])
        dl_runs = [dl_run] + [d for p_ in priors for d in sorted(glob.glob(os.path.join(p_, "runs", "*")))]
        built = built_plan_check(out, merged, pmeta, presults, dl_runs, library_root,
                                 merged.get("library_commit") or library_commit)
    bases = {}
    rows, _ = sc.read_jsonl(os.path.join(dl_run, "arms", "direct-local", "executions.jsonl"))
    for r in rows:
        m = meta.get(r.get("scenario")) or {}
        if m.get("base_of"):
            bases[m["base_of"]] = {"scenario": r.get("scenario"), "http": r.get("http"),
                                   "engine_ms": (r.get("server") or {}).get("elapsed_ms"), "wall_ms": r.get("wall_ms")}
    pair = {}
    if built and not built.get("error"):
        cmp_ = (built.get("compares") or {}).get("G-annot-0001") or {}
        rows_total = sum(sum((v.get("counts") or {}).values()) for v in cmp_.values()
                         if isinstance(v.get("counts"), dict))
        paged = (cmp_.get("cmp-paged") or {}).get("pages")
        pair = {"compares": cmp_, "rows": rows_total, "cmp_paged_pages": paged,
                "verdict": "pass" if rows_total > 0 and (paged or 0) > 1 else
                "fail: rows %d, cmp-paged %s page(s)" % (rows_total, paged)}
    record = {"stops": stops, "partner": res["results"], "built_plans": built, "bases": bases, "pair": pair}
    with open(os.path.join(out, "rehearsal_s.json"), "w") as fh:
        json.dump(record, fh, indent=1, sort_keys=True, default=str)
    lines = ["# Stop rehearsal (rehearsal S) report", "", "## Entry verdict (%d of %d taken stops, every check passing)"
             % (STOP_TRIES, STOP_TRIES), "", "| stop scenario | tries | taken | passed | verdict |", "|---|---|---|---|---|"]
    for sid_, e in stops["entry"].items():
        lines.append("| %s | %d | %d | %d | %s |" % (sid_, e["tries"], e["taken"], e["passed"], e["verdict"]))
    lines += ["", "## Tries", "", "| try | taken | race | outcome | checks not passing |", "|---|---|---|---|---|"]
    for tid, t in stops["tries"].items():
        bad = "; ".join("%s %s: %s" % (k, v[0], str(v[1])[:120]) for k, v in t["checks"].items() if v[0] != "pass")
        lines.append("| %s | %s | %s | %s | %s |" % (tid, t["taken"], t["race"] or "-", t["outcome"],
                                                     bad.replace("|", "/") or "-"))
    lines += ["", "## The base requests, once to completion on direct-local (the warm engine time the stop must beat)",
              ""]
    for stop_id, b_ in sorted(bases.items()):
        lines.append("- %s (%s): http %s, engine %s ms, wall %s ms" % (stop_id, b_["scenario"], b_["http"],
                                                                      b_["engine_ms"], b_["wall_ms"]))
    if pair:
        lines += ["", "## The constrain / annotate pair (R2-5): %s" % pair["verdict"], "",
                  "compare rows (all modes) %d; cmp-paged %s page(s)" % (pair["rows"], pair["cmp_paged_pages"])]
    for sid_, r in res["results"].items():
        lines += ["", "## The rehearsed graphlet %s: %s" % (sid_, "passed" if r.get("passed") else "not passed: %s" % "; ".join(
            "%s %s" % kv for kv in r.get("verdicts", {}).items() if not str(kv[1]).startswith("pass")))]
    if built:
        lines += ["", "## The built plans of OUT/corpus.json: %d scenario(s), %d op(s), problems in %d, not rehearsed "
                  "%s" % (built.get("checked", 0), built.get("ops", 0), len(built.get("problems") or {}),
                          ", ".join(built.get("not_rehearsed") or []) or "none")]
    with open(os.path.join(out, "REPORT-S.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return record


def rehearsal_manifest(scenarios: List[dict], skipped: List[dict]) -> dict:
    return {"manifest_version": 8, "rehearsal": "area G (rehearse_g.py)", "partial_build": True,
            "content_sha256": ml.sha256_json(scenarios), "created_utc": utc_now(), "index": {"database": ml.DB},
            "skipped": skipped, "relations_sha256": ml.sha256_json([]),
            "relations": [], "scenarios": scenarios}


# ===================================================================================================== the plan
# the engine time the plan expects per candidate on a cold host (v7's evidence where it has one, else the slot's)
EXPECTED_COLD_MS = {"two_arms": 2500, "merges": 3000, "one_arm": 2000, "cut_lists": 3000, "resource_stop": 800,
                    "depth_zero": 500, "linear": 1500, "failed_seed": 500, "long_walk": 20000, "trace": 4000,
                    "annotate": 1500, "no_bases": 3000, "over_8mib": 30000}
MCP_OVERHEAD_S = 6.0               # submit, queue, terminal lag, results, export, download (v7: ~2-5 s)
DL_PAUSE_S = 1.0                   # run_direct's --pause between requests


def plan(scenarios: List[dict], meta: Dict[str, dict]) -> dict:
    """What the rehearsal sends: per candidate the MCP job (expressible ones) and DL_REPEATS + 1 direct-local
    requests (+ the run's three capabilities GETs and the closing one), and the expected duration (cold engine times
    from v7's evidence, warm direct repeats at a tenth of it)."""
    rows = []
    mcp_jobs = dl = 0
    seconds = 0.0
    for s in scenarios:
        m = meta[s["id"]]
        cold = EXPECTED_COLD_MS.get(m["slot"], 2000) / 1000.0
        expressible = "tool" in s["mcp"]
        mcp = 1 if expressible else 0
        mcp_jobs += mcp
        n_dl = DL_REPEATS + 1
        dl += n_dl
        t = (cold + MCP_OVERHEAD_S if expressible else 0.0) + (cold if not expressible else cold / 10.0) + \
            (DL_REPEATS - 1) * cold / 10.0 + cold + n_dl * DL_PAUSE_S
        seconds += t
        rows.append({"scenario": s["id"], "candidate": m["cid"], "slot": m["slot"], "pair": m.get("pair"),
                     "mcp_jobs": mcp, "direct_local_requests": n_dl, "expected_s": round(t, 1),
                     "time_budget_ms": s["request"]["strategy"]["bounds"]["time_budget_ms"],
                     "conditional": m["candidate"]["conditional"], "direct_only": not expressible})
    dl += 4
    return {"candidates": len(scenarios), "mcp_jobs": mcp_jobs,
            "mcp_tool_calls_about": mcp_jobs * 6, "direct_local_requests": dl, "direct_remote_requests": 0,
            "expected_minutes": round(seconds / 60.0, 1), "rows": rows,
            "load": "one request in flight at a time: the MCP jobs first (cold, one per candidate, the async runner's "
                    "queue check), then the direct-local arm on mex (run_suite.sh --sequential --skip-remote "
                    "--no-mcp: three repeats and one twin per candidate); no direct-remote request"}


def plan2_markdown(p: dict, out: str, mpath: str, mcp_cmd: str) -> str:
    lines = ["# Area G rehearsal 2 plan", "",
             "Candidates %d and %d warm-ups: MCP jobs %d, direct-local requests %d, direct-remote 0; expected about "
             "%.0f minutes of staging time." % (p["candidates"], p["warm_ups"], p["mcp_jobs"],
                                                p["direct_local_requests"], p["expected_minutes"]), "",
             "Load: %s." % p["load"], "",
             "Command: `rehearse_g.py --rehearsal 2 --out %s --execute` = `run_suite.sh --manifest %s --paired "
             "--skip-remote --mcp-cmd '<%s>' --warm-every 0 --no-gate --no-pdf` (SUITE_RUNS=%s/runs)" % (
                 out, mpath, mcp_cmd.split(" ")[1] if " " in mcp_cmd else mcp_cmd, out), "",
             "Then offline: `rehearse_g.py --rehearsal 2 --out %s --analyse --dl-run <the run> --prior <rehearsal 1's "
             "directory>`." % out, "",
             "| scenario | candidate | MCP jobs | direct-local | expected s |", "|---|---|---|---|---|"]
    for r in p["rows"]:
        lines.append("| %s | %s | %d | %d | %.1f |" % (r["scenario"], r["candidate"], r["mcp_jobs"],
                                                        r["direct_local_requests"], r["expected_s"]))
    return "\n".join(lines) + "\n"


def plan_markdown(p: dict, out: str, mcp_cmd: str) -> str:
    lines = ["# Area G rehearsal plan", "",
             "Candidates %d (%d MCP jobs, about %d MCP tool calls), direct-local requests %d, direct-remote 0; "
             "expected about %.0f minutes of staging time." % (p["candidates"], p["mcp_jobs"],
                                                               p["mcp_tool_calls_about"], p["direct_local_requests"],
                                                               p["expected_minutes"]),
             "", "Load: %s." % p["load"], "",
             "Commands (`rehearse_g.py --out %s --execute` runs them in this order):" % out, "",
             "1. MCP, cold, one job per expressible candidate: `%s` (ids: the rehearsal manifest's, repeat 0)"
             % mcp_cmd,
             "2. direct-local x3 + twin: `run_suite.sh --manifest %s --sequential --skip-remote --no-mcp "
             "--warm-every 0 --no-gate --no-pdf`" % os.path.join(out, "rehearsal_manifest.json"),
             "3. offline: `rehearse_g.py --out %s --analyse --dl-run <the run dir of 2> --mcp-dir %s`" % (
                 out, os.path.join(out, "mcp")), "",
             "| scenario | candidate | slot | MCP | direct-local | budget ms | expected s | notes |",
             "|---|---|---|---|---|---|---|---|"]
    for r in p["rows"]:
        notes = ", ".join(x for x in (("pair %s" % r["pair"]) if r["pair"] else "",
                                      "conditional" if r["conditional"] else "",
                                      "direct-only" if r["direct_only"] else "") if x)
        lines.append("| %s | %s | %s | %d | %d | %s | %.1f | %s |" % (r["scenario"], r["candidate"], r["slot"],
                                                                     r["mcp_jobs"], r["direct_local_requests"],
                                                                     r["time_budget_ms"], r["expected_s"], notes))
    return "\n".join(lines) + "\n"


# ===================================================================================================== execution
def async_mcp_template() -> str:
    """run_suite.sh's ASYNC_MCP_CMD (the async session's MCP runner), read from the script so the two stay one."""
    with open(os.path.join(HERE, "run_suite.sh")) as fh:
        for line in fh:
            m = re.match(r'ASYNC_MCP_CMD="(.*)"\s*$', line.strip())
            if m:
                return m.group(1)
    raise SystemExit("run_suite.sh states no ASYNC_MCP_CMD")


def expand(template: str, **values: str) -> str:
    out = template
    for k, v in values.items():
        out = out.replace("{%s}" % k, shlex.quote(v))
    return out


def execute(out: str, manifest_path: str, scenarios: List[dict], mcp_cmd: str) -> int:
    """The MCP step (cold, once per expressible candidate), then the direct-local arm (run_suite.sh). Never call
    this before the coordinator's go: it loads staging."""
    ids = [s["id"] for s in scenarios if "tool" in s["mcp"]]
    mcp_dir = os.path.join(out, "mcp")
    os.makedirs(mcp_dir, exist_ok=True)
    cmd = expand(mcp_cmd, ids=",".join(ids), repeat="0", repeats="1", out=mcp_dir, manifest=manifest_path,
                 checks_dir=HERE)
    sys.stderr.write("rehearse_g: MCP step (%d jobs)\n" % len(ids))
    rc = subprocess.call(["bash", "-c", cmd], stdin=subprocess.DEVNULL)
    if rc != 0:
        sys.stderr.write("rehearse_g: the MCP step exited %d (see %s)\n" % (rc, mcp_dir))
    env = dict(os.environ, SUITE_RUNS=os.path.join(out, "runs"))
    sys.stderr.write("rehearse_g: direct-local step\n")
    rc2 = subprocess.call([os.path.join(HERE, "run_suite.sh"), "--manifest", manifest_path, "--sequential",
                           "--skip-remote", "--no-mcp", "--warm-every", "0", "--no-gate", "--no-pdf"], env=env,
                          stdin=subprocess.DEVNULL)
    return rc or rc2


# ===================================================================================================== facts
def document_facts(response: dict, full: Optional[dict] = None) -> dict:
    """The structural facts of a /traverse `detail: graphlet` response (one seed) by the suite's own reader
    (checks.graphlet_arms, graphlet_label_dict), with the full twin's continuations: {arms, walks, segments, merges,
    labels, label_names, walk_bp, bases, cut_lists, max_labels_per_node, body_bytes, complete_to_bp, continuation,
    failed, stop}."""
    r = ((response or {}).get("results") or [{}])[0] or {}
    text = r.get("graphlet")
    if not isinstance(text, str):
        return {"failed": True, "error": r.get("error")}
    arms = checks.graphlet_arms(text)
    names = [x.get("name") for x in checks.graphlet_label_dict(text)]
    segs = dict((a, arms[a]["segments"]) for a in arms)
    f: Dict[str, Any] = {
        "arms": [a for a in ("left", "right") if a in arms],
        "walks": dict((a, len(arms[a]["paths"])) for a in arms),
        "segments": dict((a, len(segs[a])) for a in arms),
        "merges": dict((a, sum(1 for s in segs[a] if len(s.get("parents") or []) > 1)) for a in arms),
        "labels": len(names), "label_names": names,
        "walk_bp": dict((a, max([p["length_bp"] for p in arms[a]["paths"]] or [0])) for a in arms),
        "bases": all("sequence" in s for a in segs for s in segs[a]) if any(segs.values()) else True,
        "body_bytes": len(text.encode("utf-8")),
        "seed_bp": len(((response.get("results") or [{}])[0].get("seed") or {}).get("sequence") or "") or None,
    }
    all_segs = [s for a in segs for s in segs[a]]
    if all_segs and all("labels_total" in s for s in all_segs):
        f["max_labels_per_node"] = max(s.get("labels_total") or 0 for s in all_segs)
    # cut label lists: the engine states them (a label_lists limitation, label_evidence lower_bound); the segments'
    # entry lists (labels_truncated) can all be whole while other recorded lists (end sets, events) are cut
    lims = [x for x in (r.get("limitations") or []) if isinstance(x, dict)]
    f["cut_lists"] = (any(x.get("kind") == "label_lists" for x in lims)
                      or (r.get("outcome") or {}).get("label_evidence") == "lower_bound"
                      or any(s.get("labels_truncated") for s in all_segs))
    rarms = r.get("arms") or {}
    f["complete_to_bp"] = dict((a, (rarms.get(a) or {}).get("complete_to_bp")) for a in f["arms"])
    rs = r.get("resource_stop")
    f["stop"] = (rs or {}).get("resource") if isinstance(rs, dict) else None
    if full:
        fr = ((full.get("results") or [{}])[0] or {}).get("arms") or {}
        f["continuation"] = dict((a, any(isinstance(p.get("continuation"), dict)
                                         for p in (fr.get(a) or {}).get("paths") or [])) for a in f["arms"])
    return f


def preconditions_of(f: dict) -> dict:
    """graphlet.preconditions from measured facts (graphlet_precheck.PRECONDITION_FACTS): the arms, at least the walks,
    segments and merges seen per arm, the top three labels present, the longest walk per arm, the bases."""
    if f.get("failed"):
        return {"no_graphlet": True}
    pre: Dict[str, Any] = {"arms": list(f["arms"]), "walks_min": dict(f["walks"]),
                           "segments_min": dict(f["segments"]), "bases": bool(f["bases"])}
    if any(f["merges"].values()):
        pre["merges_min"] = dict(f["merges"])
    if f.get("label_names"):
        pre["labels_present"] = [n for n in f["label_names"][:3] if isinstance(n, str) and n]
    if any(f["walk_bp"].values()):
        pre["walk_bp_min"] = dict((a, v) for a, v in f["walk_bp"].items())
    return pre


def service_facts(f: dict) -> dict:
    """graphlet.facts for the service checks (graphlet_service_checks.scenario_facts)."""
    out = {}
    if "cut_lists" in f:
        out["cut_lists"] = f["cut_lists"]
    if "max_labels_per_node" in f:
        out["max_labels_per_node"] = f["max_labels_per_node"]
    if f.get("body_bytes") is not None:
        out["body_over_8mib"] = f["body_bytes"] > IN_SERVICE_BYTES
    return out


SLOT_RULES = {
    "two_arms": "two arms, at least one branching (more than one walk)",
    "merges": "a merge on an arm",
    "one_arm": "one arm with more than one walk",
    "cut_lists": "cut label lists (a segment's recorded labels fewer than its labels_total)",
    "resource_stop": "a resource stop by work units",
    "depth_zero": "complete_to_bp 0 on an arm",
    "linear": "one arm, one walk, one segment",
    "failed_seed": "no graphlet (the seed failed) on every execution",
    "long_walk": "a walk whose length with the seed exceeds the 16,384-byte sequence ceiling",
    "trace": "a graphlet (support trace)",
    "annotate": "a graphlet (annotate)",
    "no_bases": "no bases",
    "over_8mib": "a body over 8 MiB",
}


def slot_problem(slot: str, f: dict, seed_bp: int) -> str:
    """Why the facts do not show the shape |slot| needs ('' when they do: SLOT_RULES)."""
    if slot == "failed_seed":
        return "" if f.get("failed") else "the seed has a graphlet"
    if f.get("failed"):
        return "the seed failed"
    ok = {
        "two_arms": len(f["arms"]) == 2 and max(f["walks"].values()) > 1,
        "merges": any(f["merges"].values()),
        "one_arm": len(f["arms"]) == 1 and max(f["walks"].values()) > 1,
        "cut_lists": bool(f.get("cut_lists")),
        "resource_stop": f.get("stop") == "work",
        "depth_zero": any(v == 0 for v in f["complete_to_bp"].values()),
        "linear": len(f["arms"]) == 1 and list(f["walks"].values()) == [1] and list(f["segments"].values()) == [1],
        "long_walk": max(f["walk_bp"].values() or [0]) + seed_bp > SEQUENCE_CEILING,
        "trace": True, "annotate": True,
        "no_bases": f.get("bases") is False,
        "over_8mib": (f.get("body_bytes") or 0) > IN_SERVICE_BYTES,
    }.get(slot, True)
    return "" if ok else "not the slot's shape: %s" % SLOT_RULES.get(slot, slot)


# ===================================================================================================== analysis
def _read(path: str) -> Any:
    if not os.path.exists(path):
        return None
    try:
        with (gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")) as f:
            return json.load(f)
    except (OSError, ValueError, EOFError):
        return None


def envelope_digest(response: dict) -> Optional[str]:
    """The engine envelope of a one-seed response without its per-request values (graphlet_doc.envelope_view over a
    J-like envelope: the response with results cut to the seed's result without its body)."""
    if not isinstance(response, dict) or not isinstance(response.get("results"), list) or not response["results"]:
        return None
    env = dict(response)
    r0 = dict((k, v) for k, v in (response["results"][0] or {}).items() if k not in RESULT_BODY_KEYS)
    env["results"] = [r0]
    return ml.sha256_json(graphlet_doc.envelope_view(env))


def body_digest(response: dict) -> Optional[str]:
    text = (((response or {}).get("results") or [{}])[0] or {}).get("graphlet")
    if not isinstance(text, str):
        return None
    try:
        return graphlet_doc.body_digest(text)
    except ValueError:
        return "unreadable"


def engine_ms(response: dict) -> Optional[float]:
    """The engine's elapsed time of a response (its timing.elapsed_ms, else its seed result's), None without one."""
    for t in (((response or {}).get("timing") or {}).get("elapsed_ms"),
              ((((response or {}).get("results") or [{}])[0] or {}).get("timing") or {}).get("elapsed_ms")):
        if isinstance(t, (int, float)) and not isinstance(t, bool):
            return float(t)
    return None


def mcp_execution(mcp_dir: str, sid: str) -> dict:
    """The MCP execution of a candidate (the async runner's raw record, repeat 0): {status, response, echo, error}."""
    raw = _read(os.path.join(mcp_dir, sc.raw_name(sid, "mcp", 0)))
    if not isinstance(raw, dict):
        return {"status": "missing"}
    echo = None
    for c in raw.get("calls") or []:
        if c.get("tool") == "traverse_status" and isinstance((c.get("payload") or {}).get("overrides"), dict):
            echo = c["payload"]["overrides"]
    ans = raw.get("answer") or {}
    body = ans.get("body") if isinstance(ans.get("body"), dict) else None
    err = ans.get("error") or raw.get("runner_error")
    return {"status": "answered" if body else "error", "response": body, "echo": echo, "error": err,
            "not_exposed": raw.get("not_exposed"), "outcome": raw.get("outcome")}


def dl_executions(dl_run: str, sid: str, repeats: int = DL_REPEATS) -> List[dict]:
    out = []
    for rep in range(repeats):
        raw = _read(os.path.join(dl_run, "arms", "direct-local", sc.raw_name(sid, "direct-local", rep)))
        ans = (raw or {}).get("answer") or {}
        out.append({"repeat": rep, "http": ans.get("http"), "response": ans.get("body")
                    if isinstance(ans.get("body"), dict) else None})
    return out


def judge(s: dict, meta: dict, dl: List[dict], mcp: dict, twin: Optional[dict]) -> dict:
    """The verdicts of one candidate (module docstring) and its measured facts."""
    v: Dict[str, Any] = {}
    strategy = s["request"]["strategy"]
    budget = min(float(strategy["bounds"]["time_budget_ms"]), float(area_g.HOST_MAX_TIME_MS))
    expressible = "tool" in s["mcp"]
    responses = [d["response"] for d in dl] + ([mcp.get("response")] if expressible else [])
    names = ["dl%d" % d["repeat"] for d in dl] + (["mcp"] if expressible else [])
    if any(d["http"] != 200 for d in dl):
        v["direct"] = "fail: http %s" % [d["http"] for d in dl]
    if expressible:
        if mcp.get("status") != "answered" and not (meta["slot"] == "failed_seed"):
            v["mcp"] = "fail: %s %s" % (mcp.get("status"), json.dumps(mcp.get("error"))[:200])
        echo = mcp.get("echo")
        if echo is None:
            v["echo"] = "fail: no traverse_status.overrides recorded"
        else:
            d = gp.first_difference(echo, strategy)
            v["echo"] = "pass" if d is None else "fail: %s echo %s, request %s" % (gp._dotted(d[0]), gp.compact(d[1]),
                                                                                  gp.compact(d[2]))
    else:
        v["mcp"] = "not expressible (direct-only)"
    bodies = dict((n, body_digest(r)) for n, r in zip(names, responses))
    envs = dict((n, envelope_digest(r)) for n, r in zip(names, responses))
    if meta["slot"] == "failed_seed":
        v["bodies"] = "pass" if all(b is None for b in bodies.values()) else "fail: a graphlet %s" % bodies
    else:
        vals = set(bodies.values())
        v["bodies"] = "pass" if len(vals) == 1 and None not in vals and "unreadable" not in vals else \
            "fail: %s" % dict((k, (x or "none")[:12]) for k, x in bodies.items())
    if meta["slot"] == "failed_seed":
        # the service keeps no body for a failed seed (its record's summary only): the error and outcome compared
        fails = dict((n, (((r or {}).get("results") or [{}])[0] or {})) for n, r in zip(names, responses))
        keys = set(json.dumps([x.get("error"), (x.get("outcome") or {}).get("walks")], sort_keys=True)
                   for x in fails.values())
        v["envelope"] = "pass (failed seed: error and outcome equal)" if len(keys) == 1 else \
            "fail: the failed seed's error or outcome differs: %s" % sorted(keys)
    else:
        ev = set(envs.values())
        v["envelope"] = "pass" if len(ev) == 1 and None not in ev else "fail: %s" % dict(
            (k, (x or "none")[:12]) for k, x in envs.items())
    times = dict((n, engine_ms(r)) for n, r in zip(names, responses))
    worst = max([t for t in times.values() if t is not None] or [None]) if any(
        t is not None for t in times.values()) else None
    limit = area_g.HEADROOM_FRACTION * budget
    v["headroom"] = "pass" if worst is not None and worst <= limit else \
        "fail: engine %s ms against %d ms (25 %% of %d ms)" % (worst, limit, budget)
    resp0 = dl[0]["response"] if dl else None
    facts = document_facts(resp0, twin) if resp0 else {"failed": True}
    flags = area_g.facts_of(meta["candidate"], area_g.default_corpus())
    for k in ("coordinates", "annotate"):
        if flags.get(k) and k not in facts:
            facts[k] = flags[k]
    if twin is None:
        v["twin"] = "fail: no twin kept"
    elif facts.get("failed"):
        v["twin"] = "pass" if not ((twin.get("results") or [{}])[0] or {}).get("arms") else "fail: the twin walked"
    else:
        diff = gp.twin_difference(resp0["results"][0]["graphlet"], twin, checks)
        v["twin"] = "pass" if diff is None else "fail: %s%s" % (gp.compact(diff), twin_cause(resp0, twin))
    seed_bp = len(s["request"]["seeds"][0]["sequence"])
    sp = slot_problem(meta["slot"], facts, seed_bp) if not meta.get("pair") else ""
    v["slot"] = "pass" if not sp else "fail: " + sp
    return {"verdicts": v, "facts": facts, "times_ms": times, "bodies": bodies, "echo": mcp.get("echo")}


def twin_cause(response: dict, twin: dict) -> str:
    """Why a twin differs, where the reading can say: a merge whose parents are one set in another order means
    SPEC-core §7.1 (R21 (4), feature level 6): `graphlet` detail shows a merge's majority parent first, `full` detail
    under a memory budget (the presets' max_memory_mb) keeps the arrival order, so the first-parent paths differ."""
    try:
        arms = checks.graphlet_arms(response["results"][0]["graphlet"])
        farms = ((twin.get("results") or [{}])[0] or {}).get("arms") or {}
        for a, arm in arms.items():
            fs = (farms.get(a) or {}).get("segments") or []
            for i, seg in enumerate(arm["segments"]):
                fp = fs[i].get("parents") if i < len(fs) else None
                if len(seg["parents"]) > 1 and fp is not None and seg["parents"] != fp and \
                        sorted(seg["parents"]) == sorted(fp):
                    return (" (by design: SPEC-core §7.1 R21 (4), a merge's majority parent first in graphlet detail, "
                            "the arrival order in full detail under the memory budget max_memory_mb %s)"
                            % ((twin.get("strategy") or {}).get("bounds") or {}).get("max_memory_mb"))
    except Exception:
        pass
    return ""


def passed(verdicts: dict) -> bool:
    return all(str(x).startswith(("pass", "not expressible")) for x in verdicts.values())


def shadow_run(manifest: dict, sids: List[str], dl_run: str, library_root: str, work: str, repeat: int = 0) -> Dict:
    """The op plans of |sids| run by the local library (run_graphlet_local, W1) on the direct-local documents of
    repeat |repeat|: {(scenario, op): record}."""
    os.makedirs(work, exist_ok=True)
    sub = {"manifest_version": 8, "scenarios": [s for s in manifest["scenarios"] if s["id"] in sids]}
    sub["content_sha256"] = ml.sha256_json(sub["scenarios"])
    mpath = os.path.join(work, "manifest.json")
    with open(mpath, "w") as fh:
        json.dump(sub, fh)
    out = os.path.join(work, "out")
    import shutil
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(os.path.join(out, "arms"), exist_ok=True)
    os.symlink(os.path.join(os.path.abspath(dl_run), "arms", "direct-local"), os.path.join(out, "arms",
                                                                                          "direct-local"))
    import run_graphlet_local as rgl
    rgl.main(["--run-dir", out, "--manifest", mpath, "--library-root", library_root, "--ways", "W1",
              "--repeat", str(repeat), "--quiet"])
    rows, _ = sc.read_jsonl(os.path.join(out, "arms", "dl-lib", "graphlet_ops.jsonl"))
    return dict(((r.get("scenario"), r.get("op")), r) for r in rows)


def plan_problems(s: dict, facts: dict, recs: Dict) -> List[dict]:
    """The ops of |s| whose W1 outcome is not area_g.expected_outcome's."""
    out = []
    for o in s["graphlet"]["ops"]:
        if "W1" not in (o.get("ways") or s["graphlet"]["ways"]["default"]):
            continue
        r = recs.get((s["id"], o["id"]))
        got = "missing" if r is None else "selector" if r.get("selector_error") else \
            ((r.get("error") or {}).get("code") or "ok") if r.get("error") else "ok"
        want = area_g.expected_outcome(o, facts)
        if got != want:
            out.append({"op": o["id"], "got": got, "want": want,
                        "detail": gp.compact((r or {}).get("selector_error") or ((r or {}).get("error") or {}).get(
                            "fields"))})
    return out


# what W3's answer adds to the library's beside the rows (its 36-character handle, walk_display, label_evidence,
# qualified_cause, hint): a whole list designed to fit 262,144 bytes keeps this much room in W1's answer
W3_ROOM_BYTES = 4096


def fit_problems(s: dict, recs: Dict) -> List[str]:
    """The ops of |s| W3 runs at the service's ceiling (max_bytes 262,144) whose W1 answer would not fit it uncut
    (DESIGN-v8 §5): a row cut (row_truncated), or an answer within W3_ROOM_BYTES of the ceiling."""
    out = []
    for o in s["graphlet"]["ops"]:
        ways = o.get("ways") or s["graphlet"]["ways"]["default"]
        if "W3" not in ways or o["args"].get("max_bytes") != area_g.MAXB:
            continue
        r = recs.get((s["id"], o["id"])) or {}
        ans = r.get("answer")
        if not isinstance(ans, dict):
            continue
        size = len(json.dumps(ans, separators=(",", ":"), ensure_ascii=False).encode("utf-8"))
        if ans.get("row_truncated") or size > area_g.MAXB - W3_ROOM_BYTES:
            out.append("%s: %d bytes%s" % (o["id"], size, ", a row cut" if ans.get("row_truncated") else ""))
    return out


def flags_from(problems: List[dict], facts: dict) -> bool:
    """Turn off the plan's selector flags an op's selector failure shows (FLAG_OF_OP); whether any changed."""
    changed = False
    for p in problems:
        if p["got"] != "selector":
            continue
        for stem, key in FLAG_OF_OP:
            if p["op"].startswith(stem):
                arm = p["op"][len(stem):].split("-")[0]
                if (facts.get(key) or {}).get(arm) is not False:
                    facts.setdefault(key, {})[arm] = False
                    changed = True
    return changed


# W3's headroom (DESIGN-v8 §3, refined by the coordinator 2026-10-20): an op W3 runs keeps a third of its lane's
# limit: an interactive op a third of the service's alarm (lanes.interactive.alarm_s: 3 s -> 1 s; past the alarm the
# service escalates the op ONCE to the heavy lane, same digest, same answer), a heavy op a third of the heavy lane's
# timeout (lanes.heavy.timeout_s: 120 s -> 40 s). An escalation is allowed and recorded: the escalated op is judged
# as a heavy one. Only a stated escalation counts (the record's op meta or the X-Graphlet-Op-Escalated header), never
# one inferred from the time (VERIFY-D VD10: monotonic, an unstated interactive op over the third fails however long).
# The values come from the capabilities' service_settings.graphlet.lanes (async's limits block: the run's recorded
# live document, else the by-design table's service_settings block); a part absent or served null
# (a part the service failed to compute) is "not stated": its fallback below is used and the missing limit recorded.
W3_HEADROOM_DIVISOR = 3.0
GRAPHLET_LANES = "service_settings.graphlet.lanes"
LANES_FALLBACK = {"interactive.alarm_s": 3.0, "interactive.timeout_s": 20.0, "heavy.timeout_s": 120.0,
                  "heavy_ops": ("graphlet_compare", "graphlet_export:fasta|gfa|json", "graphlet_labels:name")}


def graphlet_lanes(refs: Optional[Any] = None) -> dict:
    """The lanes W3's headroom reads: {alarm_s, interactive_timeout_s, heavy_timeout_s, heavy_ops (a frozenset: the
    served list's order follows the service's tuples), sources {part: live | table | fallback}, missing [what the
    capabilities did not state]}. |refs| a suite_common.References (default: the by-design table's blocks). Never
    raises: a document it cannot read is a missing limit like any other."""
    refs = refs if refs is not None else sc.default_references()
    got: Dict[str, Any] = {}
    sources: Dict[str, str] = {}
    missing: List[str] = []
    for part, fallback in LANES_FALLBACK.items():
        path = "%s.%s" % (GRAPHLET_LANES, part)
        try:
            value, source, _mismatch = refs.resolve("service", path)
        except Exception as e:   # a document of an unexpected shape: not stated
            value, source = None, "unknown (%s)" % type(e).__name__
        raw = value
        if part == "heavy_ops":
            ok = isinstance(raw, list) and bool(raw) and all(isinstance(x, str) and x for x in raw)
            value = frozenset(raw) if ok else frozenset(fallback)
        else:
            ok = isinstance(raw, (int, float)) and not isinstance(raw, bool) and raw > 0
            value = float(raw) if ok else float(fallback)
        if not ok:
            why = "absent" if source.startswith("unknown") else ("null" if raw is None else "not a valid value")
            if source.startswith("unknown"):
                why = _null_ancestor(getattr(refs, "live", {}).get("service"), path) or why
            missing.append("%s (%s; fallback %s)" % (path, why, sorted(fallback) if part == "heavy_ops" else
                                                     fallback))
        got[part] = value
        sources[part] = source if ok else "fallback"
    return {"alarm_s": got["interactive.alarm_s"], "interactive_timeout_s": got["interactive.timeout_s"],
            "heavy_timeout_s": got["heavy.timeout_s"], "heavy_ops": got["heavy_ops"], "sources": sources,
            "missing": missing}


def _null_ancestor(doc: Any, path: str) -> Optional[str]:
    """"null at P" when the document serves an ancestor P of |path| as null (a block the service failed to compute),
    else None."""
    keys = path.split(".")
    cur = doc
    for i, k in enumerate(keys[:-1]):
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
        if cur is None:
            return "null at %s" % ".".join(keys[:i + 1])
    return None


def w3_lane(o: dict, lanes: dict) -> str:
    """The lane an op starts on at the service: heavy when an entry of lanes.heavy_ops names it ("tool", or
    "tool:qualifier": a format alternative "fasta|gfa|json" of graphlet_export, else an argument the op sends
    non-null, e.g. graphlet_labels:name), else interactive."""
    args = o.get("args") or {}
    for entry in lanes["heavy_ops"]:
        tool, _, qual = entry.partition(":")
        if tool != o["tool"]:
            continue
        if not qual:
            return "heavy"
        if tool == "graphlet_export":
            if args.get("format") in qual.split("|"):
                return "heavy"
        elif args.get(qual) is not None:
            return "heavy"
    return "interactive"


def stated_escalation(rec: dict) -> Optional[bool]:
    """Whether the service states it escalated the op (its X-Graphlet-Op-Escalated header or op meta, where the
    runner recorded them; async's records keep neither at d876dc80), None when unstated."""
    for meta in (rec.get("op_meta"), rec.get("service_meta")):
        if isinstance(meta, dict) and isinstance(meta.get("escalated"), bool):
            return meta["escalated"]
    headers = rec.get("headers")
    if isinstance(headers, dict):
        for k, v in headers.items():
            if str(k).lower() == "x-graphlet-op-escalated":
                return str(v).lower() == "true"
    return None


def w3_headroom(ops: List[dict], w3: Dict[str, dict], lanes: dict) -> dict:
    """W3's headroom over one scenario's ops (|w3|: {op: W3 record}, one execution each): {slow: [...], escalations:
    [...], timed: n}. An interactive op whose record states it escalated: recorded, judged against the heavy lane's
    third; an op stating none is interactive, judged against a third of the alarm (no inference from its time: the
    rule stays monotonic, VERIFY-D VD10)."""
    slow, escalations, timed = [], [], 0
    interactive = lanes["alarm_s"] / W3_HEADROOM_DIVISOR
    heavy = lanes["heavy_timeout_s"] / W3_HEADROOM_DIVISOR
    for o in ops:
        r = w3.get(o["id"])
        ms = (r or {}).get("elapsed_ms")
        if not isinstance(ms, (int, float)) or isinstance(ms, bool):
            continue
        timed += 1
        t = ms / 1000.0
        lane = w3_lane(o, lanes)
        limit = heavy if lane == "heavy" else interactive
        if lane == "interactive" and stated_escalation(r):
            escalations.append("%s: %.2f s, escalated to heavy (stated)" % (o["id"], t))
            lane, limit = "heavy (escalated)", heavy
        if t > limit:
            slow.append("%s: W3 %.2f s, %s lane's third %.3g s" % (o["id"], t, lane, limit))
    return {"slow": slow, "escalations": escalations, "timed": timed}


def op_timeouts(records: List[Dict], ops: List[dict], lanes: Optional[dict] = None) -> Tuple[dict, List[str]]:
    """({op: timeout_s}, [ops over a third of their lane]) from the elapsed times the records state (every way's);
    the lanes' timeouts from |lanes| (graphlet_lanes), else area_g.LANE_S."""
    out, slow, medians = {}, [], {}
    lane_s = dict(area_g.LANE_S)
    if lanes:
        lane_s = {"interactive": lanes["interactive_timeout_s"], "heavy": lanes["heavy_timeout_s"]}
    for o in ops:
        ms = [r[o["id"]]["elapsed_ms"] for r in records if o["id"] in r and
              isinstance(r[o["id"]].get("elapsed_ms"), (int, float))]
        if not ms:
            continue
        med = statistics.median(ms) / 1000.0
        medians[o["id"]] = round(med * 1000.0, 3)
        lane = lane_s[area_g.lane_of(o)]
        out[o["id"]] = round(min(lane, max(TIMEOUT_FLOOR_S, TIMEOUT_FACTOR * med)), 1)
        if med > lane / 3.0:
            slow.append("%s: median %.1f s, lane %.0f s" % (o["id"], med, lane))
    op_timeouts.medians = medians
    return out, slow


def analyse(out: str, manifest: dict, meta: Dict[str, dict], dl_run: str, mcp_dir: str, library_root: str,
            library_commit: Optional[str] = None, mcp_graphlet_run: Optional[str] = None,
            prior: Optional[str] = None) -> dict:
    """The rehearsal's judgement (module docstring): rehearsal.json, corpus.json, REPORT.md under |out|."""
    results = {}
    warm_ups = {}
    for s in manifest["scenarios"]:
        if meta.get(s["id"], {}).get("warm_up_of"):
            raw = _read(os.path.join(dl_run, "arms", "direct-local", sc.raw_name(s["id"], "direct-local", 0)))
            ans = (raw or {}).get("answer") or {}
            warm_ups[s["id"]] = {"of": meta[s["id"]]["warm_up_of"], "http": ans.get("http"),
                                 "engine_ms": engine_ms(ans.get("body") if isinstance(ans.get("body"), dict) else {})}
    scen = dict((s["id"], s) for s in manifest["scenarios"] if isinstance(s.get("graphlet"), dict))
    first_of_seed: Dict[str, str] = {}
    for sid, s in scen.items():
        dl = dl_executions(dl_run, sid, int(s.get("repeat") or 1))
        mcp = mcp_execution(mcp_dir, sid) if "tool" in s["mcp"] else {}
        twin = _read(os.path.join(dl_run, "arms", "direct-local", "docs", sid + ".full.json.gz"))
        results[sid] = judge(s, meta[sid], dl, mcp, twin)
        shared_seed_headroom(sid, s, meta[sid], results, first_of_seed)
    # the shadow run: op plans rebuilt from the measured facts, run by the local library on the direct documents
    graphed = [sid for sid, r in results.items()
               if not r["facts"].get("failed") and r["verdicts"].get("direct") is None]
    shadow_manifest = {"scenarios": []}
    other = next((sid for sid in scen if not results[sid]["facts"].get("failed")), None)
    for sid in sorted(scen):
        shadow_manifest["scenarios"].append(replan(scen[sid], meta[sid], results[sid]["facts"], other))
    recs_by_rep = []
    for attempt in range(2):
        recs = shadow_run(shadow_manifest, list(scen), dl_run, library_root, os.path.join(out, "shadow", str(attempt)))
        changed = False
        for s in shadow_manifest["scenarios"]:
            probs = plan_problems(s, results[s["id"]]["facts"], recs)
            results[s["id"]]["plan_problems"] = probs
            changed |= flags_from(probs, results[s["id"]]["facts"])
        if not changed:
            break
        shadow_manifest["scenarios"] = [replan(scen[s["id"]], meta[s["id"]], results[s["id"]]["facts"], other)
                                        for s in shadow_manifest["scenarios"]]
    recs_by_rep.append(recs)
    for rep in range(1, DL_REPEATS):
        more = [sid for sid in graphed if int(scen[sid].get("repeat") or 1) > rep]
        if more:
            recs_by_rep.append(shadow_run(shadow_manifest, more, dl_run, library_root,
                                          os.path.join(out, "shadow", "repeat%d" % rep), repeat=rep))
    w3 = {}
    if mcp_graphlet_run:
        rows, _ = sc.read_jsonl(os.path.join(mcp_graphlet_run, "arms", "mcp", "graphlet_ops.jsonl"))
        w3 = dict(((r.get("scenario"), r.get("op")), r) for r in rows)
    lanes = lanes_of_runs([d for d in (mcp_graphlet_run, mcp_dir) if d])
    for s in shadow_manifest["scenarios"]:
        r = results[s["id"]]
        r["verdicts"]["plan"] = "pass" if not r["plan_problems"] else "fail: %d op(s): %s" % (
            len(r["plan_problems"]), ", ".join("%s %s (want %s)" % (p["op"], p["got"], p["want"])
                                               for p in r["plan_problems"][:6]))
        per = [dict((op, rec) for (sid, op), rec in recs.items() if sid == s["id"]) for recs in recs_by_rep]
        per.append(dict((op, rec) for (sid, op), rec in w3.items() if sid == s["id"]))
        timeouts, r["slow_ops"] = op_timeouts(per, s["graphlet"]["ops"], lanes)
        r["op_medians_ms"] = dict(op_timeouts.medians)
        hr = w3_headroom(s["graphlet"]["ops"], per[-1], lanes)
        r["w3_headroom"] = hr
        if hr["escalations"]:
            r["escalations"] = hr["escalations"]
        r["slow_ops"] = r["slow_ops"] + hr["slow"]
        # timeout_s goes to every runner, W3's too: without W3's own times (the async runner's area G records,
        # --mcp-graphlet-run) the local library's milliseconds would set a floor no MCP call was measured against,
        # so the ops keep their lanes (area_g.set_timeouts) and only the medians are recorded
        r["timeouts"] = timeouts if w3 else {}
        r["fit_problems"] = fit_problems(s, recs_by_rep[0])
        r["verdicts"]["fit"] = "pass" if not r["fit_problems"] else "fail: %s" % "; ".join(r["fit_problems"][:4])
        if r["slow_ops"]:
            r["verdicts"]["op_headroom"] = "fail: %s" % "; ".join(r["slow_ops"][:4])
        r["passed"] = passed(r["verdicts"])
        if "sum" in per[0] and isinstance(per[0]["sum"].get("answer"), dict):
            uncut = len(json.dumps(per[0]["sum"]["answer"], separators=(",", ":")))
            r["facts"]["summary_cut_bytes"] = max(1100, uncut * 2 // 3)
    meta = dict((sid, meta[sid]) for sid in scen)
    if prior:
        results, meta, scen = merge_prior(prior, results, meta, scen)
    corpus = choose_corpus(results, meta, scen, library_commit)
    dl_runs = [dl_run] + (sorted(glob.glob(os.path.join(prior, "runs", "*"))) if prior else [])
    built = built_plan_check(out, corpus, meta, results, dl_runs, library_root, library_commit)
    write_outputs(out, results, meta, scen, corpus, warm_ups, lanes, built)
    return {"results": results, "corpus": corpus, "warm_ups": warm_ups, "lanes": lanes, "built": built}


def _paged_once(o: dict, r: Optional[dict]) -> bool:
    """An op meant to page (its id names a paging pass, paging all) whose record took one page: it does not page."""
    return bool(r) and "-paged" in o["id"] and o.get("paging") == "all" and not r.get("error") and \
        (r.get("pages") or 0) <= 1


def built_plan_check(out: str, corpus: dict, meta: Dict[str, dict], results: Dict, dl_runs: List[str],
                     library_root: str, library_commit: Optional[str]) -> dict:
    """VERIFY-D VD9: the plans the build makes from the chosen corpus (area_g.gen_G: the chosen rich set, every
    corpus-wide extra, the pairs), shadow-run by the local library (W1) on the rehearsals' own documents: each built
    scenario reads the document of the rehearsal scenario of its candidate (the ids renamed in its plan). -> {checked,
    problems {built id: [...]}, unpaged {built id: [op]}, not_rehearsed [built id], skipped_ops {built id: [op]}}.
    A built scenario whose candidate no rehearsal ran (not_rehearsed) is left out, and so are the ops that name it."""
    ctx, v7 = build_context()
    try:
        built = [x for x in area_g.gen_G(ctx, v7, corpus, library_commit) if isinstance(x.get("graphlet"), dict)]
    except SystemExit as e:
        return {"error": str(e)}
    by_cid = {}
    for sid_, m in meta.items():
        if m.get("cid") and sid_ in results:
            by_cid.setdefault(m["cid"], sid_)
    rename = dict((b["id"], by_cid.get(b["graphlet"]["candidate"])) for b in built)
    absent = sorted(b for b, r in rename.items() if r is None)
    text = json.dumps(built)
    for b, r in rename.items():
        if r:
            text = text.replace('"%s"' % b, '"%s"' % r)
    scen = [x for x in json.loads(text) if x["id"] in set(r for r in rename.values() if r)]
    skipped = {}
    for x in scen:
        keep = []
        for o in x["graphlet"]["ops"]:
            named = [b for b in absent if '"%s"' % b in json.dumps(o["args"])]
            if named:
                skipped.setdefault(x["id"], []).append(o["id"])
            else:
                keep.append(o)
        x["graphlet"]["ops"] = keep
    # one direct-local arm holding every rehearsal document the built plans read (the rehearsals' runs)
    combo = os.path.join(out, "shadow", "built-docs")
    raw = os.path.join(combo, "arms", "direct-local", "raw")
    os.makedirs(raw, exist_ok=True)
    for x in scen:
        name = sc.raw_name(x["id"], "direct-local", 0).split("/", 1)[1]
        for d in dl_runs:
            src = os.path.join(d, "arms", "direct-local", "raw", name)
            if os.path.exists(src):
                dst = os.path.join(raw, name)
                if not os.path.exists(dst):
                    os.symlink(os.path.abspath(src), dst)
                break
    recs = shadow_run({"scenarios": scen}, [x["id"] for x in scen], combo, library_root,
                      os.path.join(out, "shadow", "built"))
    back = dict((r, b) for b, r in rename.items() if r)
    problems, unpaged, compares, unfit = {}, {}, {}, {}
    for x in scen:
        for o in x["graphlet"]["ops"]:
            r = recs.get((x["id"], o["id"]))
            if o["tool"] == "graphlet_compare" and r and isinstance(r.get("answer"), dict):
                ans = r["answer"]
                compares.setdefault(back[x["id"]], {})[o["id"]] = {
                    "comparable": ans.get("comparable"), "counts": ans.get("counts"), "total": ans.get("total"),
                    "pages": r.get("pages")}
        probs = plan_problems(x, (results.get(x["id"]) or {}).get("facts") or {}, recs)
        if probs:
            problems[back[x["id"]]] = probs
        once = [o["id"] for o in x["graphlet"]["ops"] if _paged_once(o, recs.get((x["id"], o["id"])))]
        if once:
            unpaged[back[x["id"]]] = once
        fit = fit_problems(x, recs)   # W3's whole-list calls at 262,144 must stay uncut
        if fit:
            unfit[back[x["id"]]] = fit
    return {"checked": len(scen), "ops": sum(len(x["graphlet"]["ops"]) for x in scen), "problems": problems,
            "unpaged": unpaged, "unfit": unfit, "not_rehearsed": absent, "compares": compares,
            "skipped_ops": dict((back[k], v) for k, v in skipped.items())}


def lanes_of_runs(run_dirs: List[str]) -> dict:
    """graphlet_lanes of the first run directory whose MCP arm recorded the service's capabilities document, else of
    the by-design table alone (every part it lacks: the fallback, recorded missing)."""
    for d in run_dirs:
        try:
            refs = sc.run_references(d)
        except Exception:
            continue
        if refs.live.get("service") is not None:
            lanes = graphlet_lanes(refs)
            lanes["document"] = refs.origins.get("service") or d
            return lanes
    lanes = graphlet_lanes()
    lanes["document"] = None
    return lanes


# candidates renamed since the first rehearsal (the owner's decision of 2026-10-19: trace with coordinates its own slot)
RENAMED = {"trace.3": ("coordinates.1", "coordinates")}


def merge_prior(prior: str, results: Dict, meta: Dict, scen: Dict) -> Tuple[Dict, Dict, Dict]:
    """The first rehearsal's results (|prior|: its directory) with this one's: a candidate judged again takes this
    rehearsal's verdicts; the order is the first rehearsal's, then the candidates new here."""
    ctx, v7 = build_context()
    cands = dict((c["cid"], c) for c in area_g.candidates(ctx, v7))
    for _seed_cid, pc_ in [(pc["cid"], pc) for pc, _of, _pair in partner_candidates(list(cands.values()))]:
        cands[pc_["cid"]] = pc_
    old = sc.read_json(os.path.join(prior, "rehearsal.json")) or {}
    old_scen = [x for x in (sc.read_json(os.path.join(prior, "rehearsal_manifest.json")) or {}).get("scenarios") or []
                if isinstance(x.get("graphlet"), dict)]
    m_results, m_meta, m_scen = {}, {}, {}
    for s_ in old_scen:
        sid = s_["id"]
        m = dict((old.get("meta") or {}).get(sid) or {})
        if m.get("cid") in RENAMED:
            m["cid"], m["slot"] = RENAMED[m["cid"]]
        if m.get("cid") not in cands or sid not in (old.get("results") or {}):
            continue
        m["candidate"] = cands[m["cid"]]
        m_results[sid], m_meta[sid], m_scen[sid] = old["results"][sid], m, s_
    for sid, s_ in scen.items():
        if sid in m_results and (s_.get("mcp") or {}).get("not_implemented") == REHEARSAL2_DIRECT:
            # a twin re-check (direct-local only): the first rehearsal's verdicts with this one's twin, plan, fit and
            # slot readings (the MCP's, the echo and the cold engine time stay the first rehearsal's)
            r = copy.deepcopy(m_results[sid])
            r2 = results[sid]
            for k in ("twin", "plan", "fit", "slot"):
                if k in r2["verdicts"]:
                    r["verdicts"][k] = r2["verdicts"][k]
            r["facts"] = r2["facts"]
            r["twin_rechecked"] = "rehearsal 2 (the twin without bounds.max_memory_mb, F-D4)"
            r["passed"] = passed(r["verdicts"])
            m_results[sid] = r
            continue
        m_results[sid], m_meta[sid], m_scen[sid] = results[sid], meta[sid], s_
    return m_results, m_meta, m_scen


def shared_seed_headroom(sid: str, s: dict, meta: dict, results: Dict, first_of_seed: Dict[str, str]) -> None:
    """A candidate (not a pair partner: a partner runs after its graphlet in the baseline too) whose seed an earlier
    candidate of the rehearsal already walked met warm rows: its headroom is judged against the cold engine time of
    the first candidate with that seed as well (the baseline runs only the chosen one, cold)."""
    key = s["request"]["seeds"][0]["sequence"]
    first = first_of_seed.setdefault(key, sid)
    if first == sid or meta.get("pair"):
        return
    cold = [t for t in results[first]["times_ms"].values() if t is not None]
    if not cold:
        return
    budget = min(float(s["request"]["strategy"]["bounds"]["time_budget_ms"]), float(area_g.HOST_MAX_TIME_MS))
    worst = max(cold)
    if worst > area_g.HEADROOM_FRACTION * budget:
        results[sid]["verdicts"]["headroom"] = (
            "fail: the seed's cold engine time %.0f ms (%s ran it first; this candidate met warm rows) against %d ms "
            "(25 %% of %d ms)" % (worst, first, area_g.HEADROOM_FRACTION * budget, budget))


def replan(s: dict, meta: dict, facts: dict, other: Optional[str] = None) -> dict:
    """|s| with its op plan rebuilt from the measured |facts| (as the build would make it); |other| the graphlet a
    failed seed's compare names."""
    s = copy.deepcopy(s)
    c = meta["candidate"]
    f = dict(area_g.facts_of(c, area_g.default_corpus()))
    f.update(dict((k, v) for k, v in facts.items() if v is not None))
    if meta.get("pair"):
        try:
            b, extras = _partner_b(s), {}
        except ValueError:   # a partner whose compares run on the other graphlet (the constrain / annotate pair)
            b, extras = s["id"], {"no_compare": True}
        ops = area_g.op_plan(s["id"], c["slot"], f, False, [(b, meta["pair"] in ("radius", "two_seeds"))],
                             partner_of=b, extras=extras, direct_only=not area_g.expressible(c),
                             reason=area_g.direct_reason(c))
    else:
        ops = area_g.op_plan(s["id"], c["slot"], f, c["slot"] in area_g.DEFAULT_RICH,
                             [(s["id"], False)] if c["slot"] == "two_arms" else [],
                             direct_only=not area_g.expressible(c), extras={"other": other or s["id"]})
    area_g.set_timeouts(ops)
    s["graphlet"]["ops"] = ops
    return s


def _partner_b(s: dict) -> str:
    for o in s["graphlet"]["ops"]:
        b = o["args"].get("b")
        if isinstance(b, dict) and "graphlet" in b:
            return b["graphlet"]
    raise ValueError("%s: a pair without a compare" % s["id"])


def choose_corpus(results: Dict, meta: Dict, scen: Dict, library_commit: Optional[str]) -> dict:
    """corpus.json: per slot the first passing candidate (area_g order), its facts, preconditions, echo, timeouts,
    reference digests; the pairs whose partner of the chosen candidate passed; the three richest chosen graphlets."""
    order = [m["cid"] for sid, m in sorted(meta.items(), key=lambda kv: list(scen).index(kv[0]))]
    by_cid = dict((m["cid"], sid) for sid, m in meta.items())
    slots, cands = {}, {}
    for slot, _stem, _title, _cond in area_g.SLOTS:
        tried = [cid for cid in order if meta[by_cid[cid]]["slot"] == slot and not meta[by_cid[cid]].get("pair")]
        pick, why = None, []
        for cid in tried:
            r = results[by_cid[cid]]
            if r["passed"]:
                pick = cid
                break
            why.append("%s: %s" % (cid, "; ".join("%s %s" % (k, v) for k, v in r["verdicts"].items()
                                                  if not str(v).startswith(("pass", "not expressible")))))
        slots[slot] = {"candidate": pick,
                       "why": ("passed every verdict" + ("; before it %s" % "; ".join(why) if why else ""))
                       if pick else ("; ".join(why) if why else "no candidate rehearsed"),
                       "tried": tried}
    for sid, r in results.items():
        m = meta[sid]
        entry = {"verdicts": r["verdicts"], "passed": r["passed"], "facts": r["facts"], "times_ms": r["times_ms"],
                 "timeouts": r.get("timeouts") or {}, "plan_problems": r.get("plan_problems") or []}
        if r["passed"]:
            entry["preconditions"] = preconditions_of(r["facts"])
            entry["service_facts"] = service_facts(r["facts"])
            if r.get("echo") is not None:
                entry["strategy"] = r["echo"]
            body = next((b for b in r["bodies"].values() if b), None)
            entry["reference"] = {"body_digest": body, "engine_ms": r["times_ms"]}
        cands[m["cid"]] = entry
    chosen = [(slot, slots[slot]["candidate"]) for slot in slots if slots[slot]["candidate"]]
    direct_only = set(m["cid"] for m in meta.values() if not area_g.expressible(m["candidate"]))
    # the richest: most segments, the graphlets the MCP can express first (a direct-only one enters the corpus only
    # with the owner's decision), one per seed (VERIFY-D-R2 R2-5: merges.1 and cut_lists.1 are one seed and graph)
    ranked = sorted(chosen, key=lambda x: (x[1] in direct_only,
                                           -sum((cands[x[1]]["facts"].get("segments") or {}).values())))
    rich, seeds = [], set()
    for slot, cid in ranked:
        seed = (meta[by_cid[cid]]["candidate"].get("seed") or {}).get("sequence")
        if seed in seeds:
            continue
        seeds.add(seed)
        rich.append((slot, cid))
    rich = rich[:area_g.RICHEST]
    return {"version": 1, "provisional": False,
            "rehearsal": {"utc": utc_now()},
            "library_commit": library_commit, "slots": slots, "candidates": cands,
            "rich": [slot for slot, _ in rich]}


def lanes_record(lanes: Optional[dict]) -> Optional[dict]:
    if not lanes:
        return None
    return dict(lanes, heavy_ops=sorted(lanes["heavy_ops"]),
                headroom_s={"interactive": lanes["alarm_s"] / W3_HEADROOM_DIVISOR,
                            "heavy": lanes["heavy_timeout_s"] / W3_HEADROOM_DIVISOR})


def write_outputs(out: str, results: Dict, meta: Dict, scen: Dict, corpus: dict,
                  warm_ups: Optional[Dict] = None, lanes: Optional[dict] = None, built: Optional[dict] = None) -> None:
    with open(os.path.join(out, "rehearsal.json"), "w") as fh:
        json.dump({"results": results, "meta": dict((k, dict((kk, vv) for kk, vv in v.items() if kk != "candidate"))
                                                    for k, v in meta.items()), "warm_ups": warm_ups or {},
                   "w3_lanes": lanes_record(lanes), "built_plans": built},
                  fh, indent=1, sort_keys=True, default=str)
    with open(os.path.join(out, "corpus.json"), "w") as fh:
        json.dump(corpus, fh, indent=1, sort_keys=True)
    lines = ["# Area G rehearsal report", "", "## Corpus", "", "| slot | chosen | why |", "|---|---|---|"]
    for slot, e in corpus["slots"].items():
        lines.append("| %s | %s | %s |" % (slot, e["candidate"] or "-", e["why"].replace("|", "/")[:400]))
    lines += ["", "Richest: %s" % ", ".join(corpus["rich"]), "", "## Candidates", "",
              "| scenario | candidate | passed | engine ms (dl0 dl1 dl2 mcp) | verdicts not passed |",
              "|---|---|---|---|---|"]
    for sid in scen:
        r = results[sid]
        t = " ".join("%s" % (int(x) if isinstance(x, float) else x) for x in r["times_ms"].values())
        bad = "; ".join("%s: %s" % (k, v) for k, v in r["verdicts"].items()
                        if not str(v).startswith(("pass", "not expressible")))
        lines.append("| %s | %s | %s | %s | %s |" % (sid, meta[sid]["cid"], "yes" if r["passed"] else "no", t,
                                                     bad.replace("|", "/")[:300]))
    if warm_ups:
        lines += ["", "## Warm-ups (direct-local)", "", "| warm-up | of | http | engine ms |", "|---|---|---|---|"]
        for wid, w in warm_ups.items():
            lines.append("| %s | %s | %s | %s |" % (wid, w["of"], w["http"], w["engine_ms"]))
    lines += ["", "## Per-op medians (the local library, W1, over the three direct documents)", "",
              "| scenario | ops timed | median of medians ms | slowest ops (median ms) |", "|---|---|---|---|"]
    for sid in scen:
        med = results[sid].get("op_medians_ms") or {}
        if not med:
            continue
        top = sorted(med.items(), key=lambda kv: -kv[1])[:3]
        lines.append("| %s | %d | %.2f | %s |" % (sid, len(med), statistics.median(med.values()),
                                                  ", ".join("%s %.1f" % kv for kv in top)))
    if built:
        lines += ["", "## The built plans (VERIFY-D VD9: the corpus's own plans, shadow-run on the rehearsals' documents)",
                  ""]
        if built.get("error"):
            lines.append("not built: %s" % built["error"])
        else:
            lines.append("%d scenario(s), %d op(s) run by the local library (W1); %d with an op off its expected "
                         "outcome; not rehearsed: %s." % (built["checked"], built["ops"], len(built["problems"]),
                                                         ", ".join(built["not_rehearsed"]) or "none"))
            for bid, probs in sorted(built["problems"].items()):
                lines.append("- %s: %s" % (bid, "; ".join("%s %s (want %s)" % (p_["op"], p_["got"], p_["want"])
                                                          for p_ in probs[:8])))
            for bid, ops in sorted(built["unpaged"].items()):
                lines.append("- %s: one page only: %s" % (bid, ", ".join(ops)))
    lr = lanes_record(lanes)
    if lr:
        lines += ["", "## W3 headroom (DESIGN-v8 §3, refined 2026-10-20)", "",
                  "Every op W3 runs keeps a third of its lane: interactive %.3g s (a third of the %.3g s alarm), "
                  "heavy %.3g s (a third of the %.3g s heavy timeout); heavy ops: %s. An interactive op past the alarm "
                  "escalates once to heavy: allowed, recorded, judged as heavy. Read from %s (%s)%s." % (
                      lr["headroom_s"]["interactive"], lr["alarm_s"], lr["headroom_s"]["heavy"],
                      lr["heavy_timeout_s"], ", ".join(lr["heavy_ops"]), lr.get("document") or "no recorded document",
                      ", ".join("%s %s" % kv for kv in sorted(lr["sources"].items())),
                      "; missing: " + "; ".join(lr["missing"]) if lr["missing"] else ""), "",
                  "| scenario | W3 ops timed | over their third | escalations |", "|---|---|---|---|"]
        for sid in scen:
            hr = results[sid].get("w3_headroom") or {}
            if hr.get("timed") or hr.get("escalations"):
                lines.append("| %s | %d | %s | %s |" % (sid, hr.get("timed") or 0, "; ".join(hr.get("slow") or [])
                                                        or "-", "; ".join(hr.get("escalations") or []) or "-"))
    with open(os.path.join(out, "REPORT.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ===================================================================================================== main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=None, help="the rehearsal's directory (default under the v8 plans folder)")
    ap.add_argument("--plan", action="store_true", help="write the manifest and the plan only (the default)")
    ap.add_argument("--execute", action="store_true", help="run the rehearsal on staging (after the coordinator's go)")
    ap.add_argument("--analyse", action="store_true", help="judge a recorded rehearsal (offline)")
    ap.add_argument("--dl-run", help="--analyse: the direct-local run directory (run_suite.sh's)")
    ap.add_argument("--mcp-dir", help="--analyse: the MCP step's directory (default OUT/mcp)")
    ap.add_argument("--mcp-graphlet-run", help="--analyse: a run directory with the async runner's area G records")
    ap.add_argument("--mcp-cmd", default=None, help="the MCP runner template (default run_suite.sh's ASYNC_MCP_CMD)")
    ap.add_argument("--library-root", default=LIB_DEFAULT, help="the library for the shadow run (default this tree's)")
    ap.add_argument("--library-commit", default=None, help="the pin written into the corpus file")
    ap.add_argument("--skip-conditional", action="store_true", help="leave out the conditional (over 8 MiB) slot")
    ap.add_argument("--slots", default=None, help="comma-separated slots (default all)")
    ap.add_argument("--rehearsal", type=int, choices=(1, 2, 3), default=1,
                    help="2: the second rehearsal (the twin re-checks, the long walk with its warm-ups); 3: the stop "
                         "rehearsal (S: each stop scenario 3 times, the constrain / annotate partner)")
    ap.add_argument("--prior", default=None, help="--analyse: the first rehearsal's directory, merged into the corpus "
                                              "(--rehearsal 3: comma-separated, the latest first; every one's runs "
                                              "lend their documents)")
    a = ap.parse_args(argv)
    out = os.path.abspath(a.out or default_out())
    os.makedirs(out, exist_ok=True)
    ctx, v7 = build_context()
    slots = [x.strip() for x in a.slots.split(",")] if a.slots else None
    if a.rehearsal == 3:
        scenarios, meta = rehearsal3_scenarios(ctx, v7)
    elif a.rehearsal == 2:
        scenarios, meta = rehearsal2_scenarios(ctx, v7)
    else:
        scenarios, meta = rehearsal_scenarios(ctx, v7, a.skip_conditional, slots)
    skipped = list(getattr(ctx, "skipped", []) or [])
    manifest = rehearsal_manifest(scenarios, skipped)
    mpath = os.path.join(out, "rehearsal_manifest.json")
    if not os.path.exists(mpath) or a.plan or not (a.execute or a.analyse):
        with open(mpath, "w") as fh:
            json.dump(manifest, fh, sort_keys=True)
    else:
        manifest = sc.read_json(mpath)
    mcp_cmd = a.mcp_cmd or async_mcp_template()
    p = plan3(manifest["scenarios"], meta) if a.rehearsal == 3 else plan2(manifest["scenarios"], meta) \
        if a.rehearsal == 2 else plan(manifest["scenarios"], meta)
    with open(os.path.join(out, "plan.json"), "w") as fh:
        json.dump(p, fh, indent=1, sort_keys=True)
    with open(os.path.join(out, "PLAN.md"), "w") as fh:
        fh.write(plan3_markdown(p, out, mpath, mcp_cmd) if a.rehearsal == 3 else plan2_markdown(p, out, mpath, mcp_cmd)
                 if a.rehearsal == 2 else plan_markdown(p, out, mcp_cmd))
    if a.execute:
        if a.rehearsal in (2, 3):
            return execute2(out, mpath, mcp_cmd)
        return execute(out, mpath, manifest["scenarios"], mcp_cmd)
    if a.analyse:
        if not a.dl_run:
            ap.error("--analyse needs --dl-run")
        if a.rehearsal == 3:
            rec = analyse3(out, manifest, meta, a.dl_run, a.library_root, a.library_commit, a.prior)
            print(json.dumps(rec["stops"]["entry"], indent=1))
            return 0
        default_mcp = os.path.join(a.dl_run, "arms", "mcp") if a.rehearsal == 2 else os.path.join(out, "mcp")
        res = analyse(out, manifest, meta, a.dl_run, a.mcp_dir or default_mcp, a.library_root,
                      a.library_commit, a.mcp_graphlet_run, a.prior)
        n = sum(1 for e in res["corpus"]["slots"].values() if e["candidate"])
        print("rehearsal: %d of %d slots chosen; %s" % (n, len(res["corpus"]["slots"]), os.path.join(out, "REPORT.md")))
        return 0
    if a.rehearsal == 3:
        print(json.dumps({"out": out, "stop_tries": p["stop_tries"], "mcp_jobs": p["mcp_jobs"],
                          "direct_local_requests": p["direct_local_requests"],
                          "expected_minutes": p["expected_minutes"]}, indent=1))
        return 0
    print(json.dumps({"out": out, "candidates": p["candidates"], "mcp_jobs": p["mcp_jobs"],
                      "direct_local_requests": p["direct_local_requests"], "expected_minutes": p["expected_minutes"],
                      "skipped": skipped}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
