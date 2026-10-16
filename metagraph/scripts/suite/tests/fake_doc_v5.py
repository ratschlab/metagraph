"""A FAKE run shaped like a v5 run, for the results document's tests and scaffold: every scenario of the frozen v5
manifest (manifest/suite_v5.json, or the manifest given) executed on the three arms with synthetic answers (the
status the scenario expects, every check passing), the MCP arm with its cold and warm repeats and the manifest's
marks (mcp.not_implemented: no call; a property-level mark: the marked check read not_exposed), one warm-up run again
for a resume, and samples; merged by merge_run.py and compared with the live by-design table. Nothing of it is
evidence about any host: it proves that the generator reads v5 ids, units and fields. Offline, deterministic."""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
for p in (SUITE, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import line, sample, write_arm  # noqa: E402

V5_MANIFEST = os.path.join(SUITE, "manifest", "suite_v5.json")
T0 = 1792000000.0
# the warm-up run again for a resume on direct-local, right before its target (suite_common.resume_rewarms)
RERUN_WARM_UP = "R-warm-0001"


def expected_http(scenario):
    for e in scenario.get("expect") or []:
        if isinstance(e, dict) and e.get("check") == "status":
            http = (e.get("args") or {}).get("http")
            if isinstance(http, list) and http:
                return http[0]
            if isinstance(http, int):
                return http
    return 200


def checks_for(scenario, repeat, arm):
    """Every check the scenario expects passing (determinism only on a later repeat); on the MCP arm a check a
    property-level mark names reads not_exposed (suite_common.mcp_check_marks)."""
    marks = sc.mcp_check_marks(scenario) if arm == "mcp" else {}
    out = []
    for e in scenario.get("expect") or []:
        name = e.get("check") if isinstance(e, dict) else None
        if not name:
            continue
        if name in sc.LATER_REPEATS_ONLY and repeat == 0:
            continue
        if name in marks:
            out.append({"check": name, "outcome": "not_exposed", "detail": "the MCP answer carries no such field"})
        else:
            out.append({"check": name, "outcome": "pass", "detail": "ok"})
    return out


def build_v5_run(run_dir, manifest_path=None):
    """Writes the arms of the fake run under |run_dir|, merges them and writes compare.json; returns the manifest."""
    manifest = sc.read_json(manifest_path or V5_MANIFEST)
    scenarios = [s for s in manifest.get("scenarios") or [] if isinstance(s, dict) and s.get("id")]
    for arm in sc.ARMS:
        os.makedirs(os.path.join(run_dir, "arms", arm), exist_ok=True)
        sc.write_json_atomic(os.path.join(run_dir, "arms", arm, "manifest.json"), manifest)
    local, remote, mcp = [], [], []
    t = T0
    ph = {"submit_ms": 40.0, "queue_ms": 150.0, "backend_ms": 200.0, "terminal_lag_ms": 60.0, "paging_ms": 20.0}
    for s in scenarios:
        sid = s["id"]
        http = expected_http(s)
        repeats = s.get("repeat") if isinstance(s.get("repeat"), int) and s.get("repeat") > 0 else 1
        ni = (s.get("mcp") or {}).get("not_implemented")
        # the paired loop's order: MCP cold, direct-local, direct-remote, MCP warm
        if not ni:
            mcp.append(line(sid, "mcp", 0, t, 600.0, http=http, checks=checks_for(s, 0, "mcp"), digest="d-%s" % sid,
                            phases=ph, server_ms=5.0))
        t += 1.0
        for r in range(repeats):
            local.append(line(sid, "direct-local", r, t, 12.0, http=http, checks=checks_for(s, r, "direct-local"),
                              digest="d-%s" % sid, memory=1000, server_ms=4.0))
            t += 0.5
        for r in range(repeats):
            remote.append(line(sid, "direct-remote", r, t, 60.0, http=http, checks=checks_for(s, r, "direct-remote"),
                               digest="d-%s" % sid, server_ms=4.0))
            t += 0.5
        if ni:
            mcp.append(line(sid, "mcp", 0, t, 1.0, outcome="not_implemented_mcp", http=None,
                            checks=[], digest="none"))
        else:
            for r in range(1, max(2, repeats)):
                mcp.append(line(sid, "mcp", r, t, 400.0, http=http, checks=checks_for(s, r, "mcp"),
                                digest="d-%s" % sid, phases=ph, server_ms=4.0))
                t += 0.5
        t += 1.0
    # the warm-up's re-run for a resume: right before its target's direct-local execution
    by_id = dict((s["id"], s) for s in scenarios)
    warm = by_id.get(RERUN_WARM_UP)
    target = warm.get("warms") if warm else None
    if target and any(l["scenario"] == target for l in local):
        t_target = min(sc.parse_utc(l["started_utc"]) for l in local if l["scenario"] == target)
        rerun = line(RERUN_WARM_UP, "direct-local", 0, t_target - 0.3, 10.0, http=expected_http(warm),
                     checks=checks_for(warm, 0, "direct-local"), digest="d-%s" % RERUN_WARM_UP, server_ms=4.0)
        rerun["rerun_for_resume"] = True
        rerun["rerun_target"] = target
        rerun["raw"] = sc.rerun_raw_name(RERUN_WARM_UP, "direct-local", 0, t_target - 0.3)
        local.append(rerun)
    samples = [sample(T0 + i * 60.0, 50.0 + i, 10 << 30, 3000.0, 40.0) for i in range(int((t - T0) // 60) + 1)]
    write_arm(run_dir, "direct-local", local, samples=samples,
              run_extra={"manifest": {"computed_sha256": manifest.get("content_sha256"),
                                      "manifest_version": manifest.get("manifest_version"),
                                      "scenarios_total": len(scenarios)}})
    write_arm(run_dir, "direct-remote", remote)
    write_arm(run_dir, "mcp", mcp, run_extra={"git_head": "fake-mcp-runner"})
    merge_run.merge(run_dir)
    table = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
    sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, table))
    return manifest


if __name__ == "__main__":
    build_v5_run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    print("built", sys.argv[1])
