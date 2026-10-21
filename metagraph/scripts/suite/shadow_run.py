#!/usr/bin/env python3
"""shadow_run.py: one full manifest v8 run end to end, the SHADOW RUN (the coordinator, round 6; nothing is frozen).

    shadow_run.py [--out DIR] [--dry-run] [--from STEP] [--run-dir RUN] [--mcp-cmd TEMPLATE]
                  [--library-commit SHA] [--expect-service-commit SHA] [--skip-remote]

Steps, in order (each a subprocess; --dry-run prints them and the expected load, sends nothing):
  build    build_manifest.py --manifest-version 8 --stop-scenarios into DIR (outside manifest/: not frozen)
  declared the owner's v7 declared entries (gate/declared_v7_rt.json) re-pinned to the shadow manifest into DIR (v7's
           725 scenarios are byte-identical inside v8, so the entries name the same scenarios; nothing in gate/ moves)
  policy   gate_policy.py on that manifest: the v8 policy, v7's rtp parameters (routes resolve / traverse /
           capabilities / pattern, arms direct-local + mcp, the re-pinned declared file), with its graphlets block
  library  the deployed library (--library-commit, default the corpus's pin) as a tree for the local ways
           (git archive of metagraph/api/python into DIR/library-<commit>)
  run      run_suite.sh --manifest <it> --paired, all three arms (direct-local, direct-remote, mcp: async's runner, the
           MCP arm and W3's graphlet ops), --no-gate (the gate comes last, after area G's ways)
  local    run_graphlet_local.py --run-dir RUN: W1, W2, W4 at the deployed library, drift at this tree's
  oracle   graphlet_oracle.py --run-dir RUN: W5
  compare  compare_graphlets.py RUN: the pre-checks, the five-way comparison, the service checks, the by-design
           pairs (RUN/graphlet_compare.json, RUN/arms/precheck/)
  report   report.py RUN --no-pdf --reuse: the report with chapter G
  gate     gate.py RUN DIR/policy_v8_shadow.json [--expect-service-commit]: the v1-v8 judgement with area G
--from STEP resumes at a step (RUN from --run-dir, else the newest run under DIR/runs): local and oracle first set
the records they wrote before aside (RUN/arms/<way dir>.prev-<UTC>, logged), the library is reused only when its
extraction is whole (its mark, written last, names the commit; logged), and every step's log line names the code's
commit (+dirty when the tree has edits). Staging is touched by the run step only; the coordinator's GO comes first."""

import argparse
import datetime
import json
import os
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
STEPS = ("build", "declared", "policy", "library", "run", "local", "oracle", "compare", "report", "gate")
POLICY_ROUTES = ("POST /resolve", "POST /traverse", "GET /traverse/capabilities", "POST /pattern")
POLICY_ARMS = ("direct-local", "mcp")
DECLARED = os.path.join(HERE, "gate", "declared_v7_rt.json")
TRUTH_CACHE = os.path.expanduser("~/.cache/metagraph-traverse-plans/2026-10-14/lane-s5/truth-cache")
ADMIN_KEY_FILE = "/Users/raetsch/tmp/mgss/secrets/staging_admin_key"
# a paired v7 run of 725 scenarios on staging: 3,015 executions in 3 h 51 min (20261019T005230Z)
V7_SECONDS_PER_EXECUTION = (3 * 3600 + 51 * 60) / 3015.0
W3_SECONDS_PER_CALL = 1.5        # an area G MCP graphlet call (submit, wait, page; the rehearsals' MCP arm)
LOCAL_SECONDS_PER_CALL = 0.05    # the local library's call, W1 / W2 / W4 / drift / W5


def default_out() -> str:
    now = datetime.datetime.now(datetime.timezone.utc)
    return os.path.expanduser(os.path.join("~", ".cache", "metagraph-traverse-plans", now.strftime("%Y-%m-%d"), "v8",
                                           "shadow", now.strftime("%Y%m%dT%H%M%SZ")))


def corpus_pin() -> str:
    import build_manifest as bm
    with open(bm.DEFAULT_CORPUS) as fh:
        return json.load(fh).get("library_commit") or ""


def paths(out: str, library_commit: str) -> dict:
    return {"manifest": os.path.join(out, "suite_v8_shadow.json"),
            "smoke": os.path.join(out, "suite_smoke_v8_shadow.json"),
            "policy": os.path.join(out, "policy_v8_shadow.json"), "runs": os.path.join(out, "runs"),
            "declared": os.path.join(out, "declared_v8_shadow.json"),
            "library": os.path.join(out, "library-%s" % library_commit[:12])}


def steps(a: argparse.Namespace, run_dir: str) -> list:
    """[(name, argv, env additions)] of the pipeline (module docstring)."""
    p = paths(a.out, a.library_commit)
    py = sys.executable
    env = {"SUITE_TRUTH_CACHE": TRUTH_CACHE}
    run_env = dict(env, SUITE_RUNS=p["runs"], METAGRAPH_ADMIN_KEY_FILE=ADMIN_KEY_FILE)
    run_argv = [os.path.join(HERE, "run_suite.sh"), "--manifest", p["manifest"], "--paired", "--no-gate"]
    if a.skip_remote:
        run_argv.append("--skip-remote")
    run_argv += ["--mcp-cmd", a.mcp_cmd or "default"]   # default: run_suite.sh's ASYNC_MCP_CMD (async's runner)
    gate_argv = [py, os.path.join(HERE, "gate.py"), run_dir, p["policy"]]
    if a.expect_service_commit:
        gate_argv += ["--expect-service-commit", a.expect_service_commit]
    return [
        ("build", [py, os.path.join(HERE, "build_manifest.py"), "--no-fetch", "--manifest-version", "8",
                   "--stop-scenarios", "--out", p["manifest"], "--smoke-out", p["smoke"]], env),
        ("declared", [py, os.path.abspath(__file__), "--derive-declared", p["manifest"], p["declared"]], env),
        ("policy", [py, os.path.join(HERE, "gate_policy.py"), p["manifest"], "--routes", ",".join(POLICY_ROUTES),
                    "--arms", ",".join(POLICY_ARMS), "--declared", p["declared"], "--name", "v8-shadow-rtp-mcp",
                    "--out", p["policy"]], env),
        ("library", ["sh", "-c", "rm -rf %s && mkdir -p %s && git -C %s archive %s metagraph/api/python | tar -x -C %s "
                     "&& echo %s > %s" % (
                         shlex.quote(p["library"]), shlex.quote(p["library"]), shlex.quote(REPO),
                         shlex.quote(a.library_commit), shlex.quote(p["library"]), shlex.quote(a.library_commit),
                         shlex.quote(os.path.join(p["library"], LIBRARY_MARK)))],
         {}),
        ("run", run_argv, run_env),
        ("local", [py, os.path.join(HERE, "run_graphlet_local.py"), "--run-dir", run_dir, "--manifest",
                   os.path.join(run_dir, "manifest.json"), "--ways", "W1,W2,W4,drift", "--library-root", p["library"],
                   "--library-commit", a.library_commit, "--drift-library-root", REPO], env),
        ("oracle", [py, os.path.join(HERE, "graphlet_oracle.py"), "--run-dir", run_dir], env),
        ("compare", [py, os.path.join(HERE, "compare_graphlets.py"), run_dir], env),
        ("report", [py, os.path.join(HERE, "report.py"), run_dir, "--no-pdf"], env),
        ("gate", gate_argv, env),
    ]


def derive_declared(manifest_path: str, out_path: str, source: str = DECLARED) -> dict:
    """The owner's v7 declared file re-pinned to |manifest_path| (a v8 manifest holding v7's scenarios unchanged):
    the same entries, the manifest identity the shadow's, the derivation stated."""
    import suite_common as sc
    with open(source) as fh:
        d = json.load(fh)
    with open(manifest_path) as fh:
        m = json.load(fh)
    d = dict(d, manifest={"file": manifest_path, "content_sha256": sc.manifest_identity(m)["content_sha256"]},
             about="the shadow run's re-pin of %s (shadow_run.py): its entries, for a v8 manifest whose v7 "
                   "scenarios are byte-identical to v7's; v7's file: %s" % (os.path.relpath(source, HERE),
                                                                            d.get("about", "")[:200]))
    with open(out_path, "w") as fh:
        json.dump(d, fh, indent=1, sort_keys=True)
        fh.write("\n")
    return d


def load(manifest: dict) -> dict:
    """The expected load of a paired run of |manifest| (all arms unless skip_remote is handled by the caller):
    executions per arm, area G's W3 calls and the local ways' calls, and the expected duration."""
    import area_g
    import suite_common as sc
    scen = manifest["scenarios"]
    sel = list(enumerate(scen))
    dl = len(sc.plan_executions(sel, "direct-local", 10)[0])
    dr = len(sc.plan_executions(sel, "direct-remote", 10)[0])
    mcp = sum(len(sc.mcp_repeats(i, s, 10)) for i, s in sel)
    g = area_g.summary([s for s in scen if isinstance(s.get("graphlet"), dict)])
    w3 = (g.get("calls_per_way") or {}).get("W3", 0)
    local = sum(v for k, v in (g.get("calls_per_way") or {}).items() if k != "W3")
    seconds = (dl + dr + mcp) * V7_SECONDS_PER_EXECUTION + w3 * W3_SECONDS_PER_CALL + \
        (local + (g.get("calls_per_way") or {}).get("W1", 0)) * LOCAL_SECONDS_PER_CALL
    return {"scenarios": len(scen), "direct_local": dl, "direct_remote": dr, "mcp": mcp, "area_g": g,
            "w3_calls": w3, "local_calls": local, "expected_hours": round(seconds / 3600.0, 1)}


def newest_run(runs: str) -> str:
    if not os.path.isdir(runs):
        return ""
    dirs = sorted(d for d in os.listdir(runs) if os.path.isdir(os.path.join(runs, d)))
    return os.path.join(runs, dirs[-1]) if dirs else ""


def code_identity() -> str:
    """The commit of the tree whose code the steps run (VERIFY-D-R4 R4-4), "+dirty" when it has edits."""
    try:
        head = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return head + ("+dirty" if dirty else "")
    except OSError:
        return "unknown"


def log_line(log: str, text: str) -> None:
    with open(log, "a") as fh:
        fh.write("== %s %s\n" % (datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), text))


def run_step(name: str, argv: list, env: dict, log: str, runner=subprocess.call) -> int:
    log_line(log, "%s (code %s): %s" % (name, code_identity(), " ".join(shlex.quote(x) for x in argv)))
    return runner(argv, env=dict(os.environ, **env), stdin=subprocess.DEVNULL)


# the records a post-run step writes, set aside when the step runs again (--from): way dirs under RUN/arms
STEP_OUTPUTS = {"local": ("dl-lib", "dr-lib", "mcp-lib", "drift-lib"), "oracle": ("oracle",)}
LIBRARY_MARK = ".shadow-library-commit"


def set_aside(run_dir: str, name: str, log: str) -> list:
    """Moves the records |name| wrote before (STEP_OUTPUTS) to <dir>.prev-<UTC>, logged; -> the moves."""
    moved = []
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    for d in STEP_OUTPUTS.get(name, ()):
        src = os.path.join(run_dir, "arms", d)
        if os.path.isdir(src):
            dst = "%s.prev-%s" % (src, stamp)
            os.rename(src, dst)
            moved.append((src, dst))
            log_line(log, "%s: the earlier records set aside: %s -> %s" % (name, src, dst))
    return moved


def library_complete(path: str, commit: str) -> bool:
    """Whether DIR/library-<commit> is a whole extraction of |commit| (its mark, written last, names it)."""
    mark = os.path.join(path, LIBRARY_MARK)
    if not os.path.isfile(mark) or not os.path.isdir(os.path.join(path, "metagraph", "api", "python")):
        return False
    with open(mark) as fh:
        return fh.read().strip() == commit


def main(argv=None, runner=subprocess.call) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=None, help="the shadow run's directory (default under the v8 plans folder)")
    ap.add_argument("--dry-run", action="store_true", help="print the steps and the expected load; send nothing")
    ap.add_argument("--from", dest="start", choices=STEPS, default="build", help="resume at this step")
    ap.add_argument("--to", dest="stop", choices=STEPS, default="gate", help="stop after this step")
    ap.add_argument("--run-dir", default=None, help="the run directory of the post-run steps (default: the newest)")
    ap.add_argument("--mcp-cmd", default=None, help="the MCP runner template (default run_suite.sh's ASYNC_MCP_CMD)")
    ap.add_argument("--library-commit", default=None, help="the deployed library (default the corpus's pin)")
    ap.add_argument("--expect-service-commit", default=None, help="the deployed service commit, for gate.py")
    ap.add_argument("--skip-remote", action="store_true", help="leave the direct-remote arm out (not the shadow)")
    ap.add_argument("--derive-declared", nargs=2, metavar=("MANIFEST", "OUT"), help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.derive_declared:
        derive_declared(*a.derive_declared)
        return 0
    a.out = os.path.abspath(a.out or default_out())
    a.library_commit = a.library_commit or corpus_pin()
    p = paths(a.out, a.library_commit)
    run_dir = a.run_dir or newest_run(p["runs"]) or os.path.join(p["runs"], "<the run run_suite.sh makes>")
    todo = steps(a, run_dir)
    names = [n for n, _, _ in todo]
    chosen = todo[names.index(a.start):names.index(a.stop) + 1]
    if a.dry_run:
        print("shadow run: %s (library %s)" % (a.out, a.library_commit))
        for name, argv_, env in chosen:
            print("  %-8s %s%s" % (name, " ".join("%s=%s" % kv for kv in sorted(env.items()) if kv[0] != "SUITE_RUNS"
                                                    and kv[0] != "SUITE_TRUTH_CACHE") + (" " if env else ""),
                                   " ".join(shlex.quote(x) for x in argv_)))
        if os.path.exists(p["manifest"]):
            with open(p["manifest"]) as fh:
                print(json.dumps({"expected_load": load(json.load(fh))}, indent=1))
        return 0
    os.makedirs(a.out, exist_ok=True)
    log = os.path.join(a.out, "shadow.log")
    for name in [n for n, _, _ in chosen]:
        if name in ("local", "oracle", "compare", "report", "gate"):
            # the run step made the run directory: the post-run steps read it
            run_dir = a.run_dir or newest_run(p["runs"])
            if not run_dir:
                print("shadow run: no run under %s for step %s" % (p["runs"], name), file=sys.stderr)
                return 2
        if name == "library" and library_complete(p["library"], a.library_commit):
            log_line(log, "library: reused %s (complete: its mark names %s)" % (p["library"], a.library_commit))
            continue
        if name in STEP_OUTPUTS:
            set_aside(run_dir, name, log)
        _n, argv_, env = next(x for x in steps(a, run_dir) if x[0] == name)
        rc = run_step(name, argv_, env, log, runner)
        with open(log, "a") as fh:
            fh.write("== %s exit %s\n" % (name, rc))
        if name == "gate":
            return rc
        if rc != 0:
            print("shadow run: step %s exited %s (see %s); resume with --from %s" % (name, rc, log, name),
                  file=sys.stderr)
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
