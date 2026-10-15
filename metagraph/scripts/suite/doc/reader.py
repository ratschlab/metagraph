"""Reading a recorded run for the results document (DESIGN-DOC §9.1): one `Run` object holding the manifest, the
executions (the re-evaluated set when the run has one, else the recorded set), compare.json (computed with
compare.py when absent, as report.py does), run.json, reevaluated.json, the samples, the gate verdict, the
by-design table, the owner-declared file and the checks module whose docstrings feed the check catalogue.

Everything is read from the run directory and the suite directory; nothing is written into the run.
"""

import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE_DIR = os.path.dirname(HERE)
if SUITE_DIR not in sys.path:
    sys.path.insert(0, SUITE_DIR)

import compare as cmp  # noqa: E402
import gate as gate_mod  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402

ARMS = ("direct-local", "direct-remote", "mcp")


class RunError(Exception):
    """A run directory the generator cannot read (usage / missing inputs: exit 2)."""


def _repeat_key(line):
    r = line.get("repeat")
    return r if isinstance(r, int) else -1


def _index_lines(path):
    """{(scenario, arm, repeat): 1-based physical line number} of a JSON-lines file; the last line of a key wins,
    as the merge's reading does. A missing or torn line is skipped, never raised."""
    out = {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for i, raw in enumerate(fh, 1):
                raw = raw.strip()
                if not raw:
                    continue
                try:
                    row = json.loads(raw)
                except ValueError:
                    continue
                if isinstance(row, dict):
                    out[(row.get("scenario"), row.get("arm"), row.get("repeat"))] = i
    except OSError:
        pass
    return out


def sha256_of(path):
    """The sha256 of a file, or None when it cannot be read."""
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return None


def git_commit_of(path, cwd):
    """The short commit that last touched `path` in the checkout at `cwd` (the version pin of a README
    reference), or None when git or the checkout is unavailable. Deterministic for one checkout; no network."""
    try:
        proc = subprocess.run(["git", "log", "-1", "--format=%h", "--", path], cwd=cwd, stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return None
    text = proc.stdout.decode("utf-8", "replace").strip() if proc.returncode == 0 else ""
    return text or None


def excluded_repeats(entry):
    """The MCP repeats compare left out of a scenario's agreement, each with its reason: [(repeat, reason)].
    `mcp_agreement.host_state_repeats` names the executions inconclusive by a host state (a resolve_timeout, a
    relayed deadline), which read neither for nor against the direct arm (compare.mcp_agreement); the owner's
    second review, point 2: an agreement is stated over the repeats actually compared."""
    ma = (entry or {}).get("mcp_agreement") if isinstance(entry, dict) else None
    out = []
    for r in (ma or {}).get("host_state_repeats") or []:
        out.append((r, "an inconclusive host-state outcome"))
    return out


def depth_over_compared_repeats(comp, rows, is_property):
    """report.agreement_depth judged over the MCP repeats compare compared: the executions of an excluded repeat
    (excluded_repeats) are left out before the depth is read, and every entry names `repeats_compared` and
    `repeats_excluded` ([{repeat, reason}]), so that 'fully checked' holds over the compared repeats and the
    detail block, the chapter line and the prose can name the excluded one."""
    scen = (comp or {}).get("scenarios") or {}
    excl = dict((sid, excluded_repeats(e)) for sid, e in scen.items())
    left_out = set((sid, r) for sid, items in excl.items() for r, _ in items)
    kept = [r for r in rows if not (r.get("arm") == "mcp" and (r.get("scenario"), r.get("repeat")) in left_out)]
    depth = report.agreement_depth(comp, kept, is_property)
    for sid, d in depth.items():
        reps = set(r.get("repeat") for r in kept if r.get("scenario") == sid and r.get("arm") == "mcp")
        d["repeats_compared"] = sorted(reps, key=lambda k: (not isinstance(k, int), str(k)))
        d["repeats_excluded"] = [{"repeat": r, "reason": why} for r, why in excl.get(sid) or []]
    return depth


class Run(object):
    def __init__(self, run_dir, recorded=False, by_design=None, checks_dir=None, declared=None, policy=None):
        self.dir = os.path.abspath(run_dir)
        if not os.path.isdir(self.dir):
            raise RunError("not a directory: %s" % run_dir)
        self.manifest = sc.read_json(os.path.join(self.dir, "manifest.json"))
        if not isinstance(self.manifest, dict) or not isinstance(self.manifest.get("scenarios"), list):
            raise RunError("%s/manifest.json: missing or without scenarios" % run_dir)
        self.scenarios = {}
        self.order = []
        for s in self.manifest["scenarios"]:
            if isinstance(s, dict) and s.get("id"):
                self.scenarios[s["id"]] = s
                self.order.append(s["id"])
        self.relations = [r for r in (self.manifest.get("relations") or []) if isinstance(r, dict)]
        reevaluated = (not recorded) and os.path.exists(os.path.join(self.dir, "executions.reevaluated.jsonl"))
        self.reevaluated = reevaluated
        self.run = (sc.read_json(os.path.join(self.dir, sc.reevaluated_name("run.json", reevaluated)))
                    or sc.read_json(os.path.join(self.dir, "run.json")) or {})
        rows, bad = sc.read_jsonl(os.path.join(self.dir, sc.reevaluated_name("executions.jsonl", reevaluated)))
        if not rows and not os.path.exists(os.path.join(self.dir, sc.reevaluated_name("executions.jsonl",
                                                                                     reevaluated))):
            raise RunError("%s: no executions file (merge the run first)" % run_dir)
        self.rows = rows
        self.unparsable = bad
        self.comp = report.load_compare(self.dir, by_design, reevaluated)
        self.reeval = sc.read_json(os.path.join(self.dir, "reevaluated.json")) if reevaluated else None
        self.samples, _ = sc.read_jsonl(os.path.join(self.dir, "samples.jsonl"))
        self.load = report.load_figures(self.samples)
        self.gate = (sc.read_json(os.path.join(self.dir, sc.reevaluated_name("gate.json", reevaluated)))
                     or sc.read_json(os.path.join(self.dir, "gate.json")))
        self.mcp_run = sc.read_json(os.path.join(self.dir, "arms", "mcp", "run.json")) or {}
        self.orchestrator = sc.read_json(os.path.join(self.dir, "orchestrator.json")) or {}
        judged = cmp.judging_checks(self.dir, reevaluated) or {}
        self.judged = judged
        self.checks_dir = checks_dir or judged.get("dir") or SUITE_DIR
        self.checks = cmp.load_checks_module(self.checks_dir)
        self.is_property = gate_mod._is_property_fn(self.checks_dir)
        self.depth = depth_over_compared_repeats(self.comp, self.rows, self.is_property)
        self.missing, self.plan_source = gate_mod.missing_planned(self.dir, self.manifest, self.rows)
        self.missing_set = set((k[0], k[1]) for k in self.missing)
        self.table, self.table_meta = self._load_table(by_design, reevaluated)
        self.cases = report.report_cases(self.comp, self.dir) or self.table
        self.declared = self._load_declared(declared)
        self.policy, self.policy_path = self._load_policy(policy)
        self.gate_offline = False
        self.gate_error = None
        if not isinstance(self.gate, dict) and isinstance(self.policy, dict):
            # no verdict recorded in the run: the gate is evaluated offline from the run directory (gate.judge
            # reads the merged lines, compare.json, the manifest and the plan; it writes nothing)
            try:
                self.gate = gate_mod.judge(self.dir, self.policy, reevaluated)
                self.gate["policy"]["path"] = self.policy_path
                self.gate_offline = True
            except Exception as e:  # a run the gate cannot judge: the document says why, never a crash
                self.gate = None
                self.gate_error = "%s: %s" % (type(e).__name__, e)
        self.lines = {}
        for r in self.rows:
            self.lines.setdefault((r.get("scenario"), r.get("arm")), []).append(r)
        for key in self.lines:
            self.lines[key].sort(key=_repeat_key)
        # the physical line of every execution in the executions file, for the detail blocks' artifact lines
        self.executions_file = sc.reevaluated_name("executions.jsonl", reevaluated)
        self.line_numbers = _index_lines(os.path.join(self.dir, self.executions_file))
        verdict, facts = report.verdict_text(self.run, self.comp, self.rows, reevaluated, (self.missing,
                                                                                            self.plan_source),
                                             self.manifest, self.cases)
        self.verdict = verdict
        self.verdict_facts = facts

    # -- the by-design table and the declared file ------------------------------------------------------------

    def _load_table(self, by_design, reevaluated):
        """(cases by id, meta) of the by-design table the document prints: the table the re-evaluation replayed
        (reevaluated.json by_design: path and sha256), found by its sha among the replayed path, the suite
        directory's copy and the run's snapshot, so that a live table edited after the run never changes the
        document (VERIFY-R2 H2); else the comparator's choice with the sha-mismatch warning kept in the meta
        (`warnings`, printed in I.4 and appendix F)."""
        try:
            path, warnings = cmp.table_for(self.dir, by_design, reevaluated)
        except Exception as e:  # a table the comparator cannot place: the document shows none
            path, warnings = None, ["no by-design table could be placed: %s" % e]
        used = ((self.reeval or {}).get("by_design") if isinstance((self.reeval or {}).get("by_design"), dict)
                else None) if reevaluated else None
        replayed_sha = (used or {}).get("sha256")
        if by_design is None and replayed_sha:
            candidates = [path, (used or {}).get("path"), os.path.join(SUITE_DIR, "service_by_design.json")]
            try:
                candidates.append(cmp.snapshot_table(self.dir))
            except Exception:
                pass
            for cand in candidates:
                if cand and os.path.exists(cand) and sha256_of(cand) == replayed_sha:
                    if cand != path:
                        warnings = []  # the replayed table itself: nothing to warn about
                        path = cand
                    break
        if not path or not os.path.exists(path):
            return {}, {"warnings": warnings, "replayed_sha256": replayed_sha,
                        "replayed_cases": (used or {}).get("cases")}
        data = sc.read_json(path) or {}
        cases = {}
        for c in data.get("cases") or []:
            if isinstance(c, dict) and c.get("id"):
                cases[c["id"]] = c
        sha = sha256_of(path)
        meta = {"path": path, "version": data.get("version"), "service_commit": data.get("service_commit"),
                "generated_utc": data.get("generated_utc"), "cases": len(cases), "sha256": sha,
                "replayed_sha256": replayed_sha, "replayed_cases": (used or {}).get("cases"),
                "replayed": bool(replayed_sha) and sha == replayed_sha, "warnings": list(warnings or [])}
        return cases, meta

    def _load_declared(self, declared):
        if declared:
            return sc.read_json(declared)
        sha = self.manifest.get("content_sha256")
        gate_dir = os.path.join(SUITE_DIR, "gate")
        if not sha or not os.path.isdir(gate_dir):
            return None
        for name in sorted(os.listdir(gate_dir)):
            if not (name.startswith("declared_") and name.endswith(".json")):
                continue
            data = sc.read_json(os.path.join(gate_dir, name))
            if isinstance(data, dict) and ((data.get("manifest") or {}).get("content_sha256") == sha):
                data = dict(data)
                data["_file"] = name
                return data
        return None

    def _load_policy(self, policy):
        """(policy, path): the given policy file, else the first policy of the suite's gate directory pinning the
        run's manifest (gate/policy_*.json with manifest.content_sha256 equal to the manifest's), else None."""
        if policy:
            try:
                return gate_mod.load_policy(policy), os.path.abspath(policy)
            except (OSError, ValueError):
                return sc.read_json(policy), os.path.abspath(policy)
        sha = self.manifest.get("content_sha256")
        gate_dir = os.path.join(SUITE_DIR, "gate")
        if not sha or not os.path.isdir(gate_dir):
            return None, None
        for name in sorted(os.listdir(gate_dir)):
            if not (name.startswith("policy_") and name.endswith(".json")):
                continue
            path = os.path.join(gate_dir, name)
            try:
                data = gate_mod.load_policy(path)
            except (OSError, ValueError):
                continue
            if (data.get("manifest") or {}).get("content_sha256") == sha:
                return data, path
        return None, None

    # -- accessors ----------------------------------------------------------------------------------------------

    def scenario(self, sid):
        return self.scenarios.get(sid) or {}

    def arm_lines(self, sid, arm):
        """The execution lines of a scenario on an arm, in repeat order."""
        return self.lines.get((sid, arm)) or []

    def comp_entry(self, sid):
        return (self.comp.get("scenarios") or {}).get(sid) or {}

    def comp_arm(self, sid, arm):
        return (self.comp_entry(sid).get("arms") or {}).get(arm) or {}

    def planned_missing(self, sid, arm):
        return (sid, arm) in self.missing_set

    def comp_relations(self):
        """compare's relation rows by id."""
        out = {}
        for r in self.comp.get("relations") or []:
            if isinstance(r, dict) and r.get("id"):
                out[r["id"]] = r
        return out

    def summary(self):
        return self.comp.get("summary") or {}

    def run_id(self):
        return self.run.get("run_id") or self.comp.get("run_id") or os.path.basename(self.dir)

    def line_no(self, line):
        """The physical line of an execution in the executions file (1-based), or None."""
        return self.line_numbers.get((line.get("scenario"), line.get("arm"), line.get("repeat")))

    def records(self):
        """The manifest's ground-truth records (key, label, accession, taxid, kingdom, ...), in manifest order."""
        gt = self.manifest.get("ground_truth") or {}
        return [r for r in (gt.get("records") or []) if isinstance(r, dict) and r.get("key")]

    def budget_ms(self, sid):
        """The time budget a scenario's request names (bounds.time_budget_ms, else a top-level time_budget_ms),
        or None."""
        req = self.scenario(sid).get("request") or {}
        if not isinstance(req, dict):
            return None
        bounds = req.get("bounds") if isinstance(req.get("bounds"), dict) else {}
        for v in (bounds.get("time_budget_ms"), req.get("time_budget_ms"), req.get("budget_ms")):
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                return v
        return None

    def check_doc(self, name):
        """The docstring of a check kind from the judging checks module, or None."""
        mod = self.checks
        table = getattr(mod, "CHECKS", None) if mod is not None else None
        fn = table.get(name) if isinstance(table, dict) else None
        doc = getattr(fn, "__doc__", None) if fn is not None else None
        return " ".join(doc.split()) if doc else None

    def session(self):
        return report.orchestrator_session(self.run)
