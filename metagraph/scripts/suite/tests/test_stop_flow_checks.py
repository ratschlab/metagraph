"""The stop scenarios' checks (checks.py, stop flows; manifest v8 --stop-scenarios) on synthetic stop_flow records of
async's runner (suite_stop_flow.py) with the admin job record's attempt rows and their token release record
(RELEASE-RECORD.md, 2026-10-20): a /traverse stop released on evidence (finished, tombstone, would_release with the
switch off) or a hold that fits the stop's evidence, a /resolve stop's 20 s grace or its lease cap, the lease bound
with async's margin (or token_lease_ms), a null release on a terminal job inconclusive, a race inconclusive, the direct
arms not exposed."""
import copy
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import areas  # noqa: E402
import checks  # noqa: E402

STOP_SENT, STOP_ANSWERED = "2026-10-20T10:00:02.000000+00:00", "2026-10-20T10:00:02.400Z"
LEASE = "2026-10-20T10:01:00+00:00"


def attempt(release=None, **kw):
    a = {"attempt_id": "a1", "task_id": "t1", "n": 1, "sent": True, "host_slug": "mex", "server_instance": "s1",
         "expect_server_instance": "s1", "dispatched_at": "2026-10-20T10:00:00+00:00", "lease_until": LEASE,
         "settled_at": "2026-10-20T10:00:03+00:00", "settled_by": "cancel", "error_class": None, "http_status": 200,
         "bound_ms": 60000, "not_after_ms": 61000, "has_semaphore_token": False, "release_stated": True,
         "release": release}
    a.update(kw)
    return a


def record(tool="traverse_fetch", attempts=None, stop_answer=None, status="stopped", finished_at=None, **kw):
    ans = stop_answer if stop_answer is not None else {"ok": True, "status": "stopping", "running": 1,
                                                       "cancel_supported": True, "bound_seconds": 58}
    rec = {"kind": "stop", "tool": tool, "job_id": "j1", "after": "running", "delay_s": 2.0, "race": None,
           "running_seen": {"utc": "2026-10-20T10:00:00.5Z", "reads": 2, "attempt_id": "a1"},
           "stop": {"sent_utc": STOP_SENT, "answered_utc": STOP_ANSWERED, "answer": ans, "error": None},
           "terminal": {"utc": "2026-10-20T10:00:05Z", "status": status, "terminal": True,
                        "finished_at": finished_at or "2026-10-20T10:00:04+00:00"},
           "attempts": attempts if attempts is not None else [attempt(released("finished", detail="finished"))],
           "tasks": [], "release_reads": [], "admin": "read", "error": None}
    rec.update(kw)
    return rec


_NULL = object()


def released(reason, site="stop", decision=_NULL, detail=None, at="2026-10-20T10:00:04+00:00"):
    """A release record; its decision null for completed and grace (async's rule), else released, unless given."""
    if decision is _NULL:
        decision = None if reason in ("completed", "grace") else "released"
    return {"at": at, "reason": reason, "site": site, "decision": decision, "detail": detail,
            "decided_at": at}


def scenario(route="traverse"):
    if route == "resolve":
        return {"id": "S-stop-0002", "route": "POST /resolve", "mcp": {"tool": "traverse_resolve"},
                "expect": copy.deepcopy(areas.STOP_CHECKS_RESOLVE)}
    return {"id": "S-stop-0001", "route": "POST /traverse", "mcp": {"tool": "traverse_fetch"},
            "expect": copy.deepcopy(areas.STOP_CHECKS_TRAVERSE)}


def run(rec, route="traverse", arm="mcp"):
    res = checks.evaluate(scenario(route), {"http": 200, "body": {"stop_flow": rec}}, arm)
    return dict((r["check"], (r["outcome"], r["detail"])) for r in res)


class TraverseStop(unittest.TestCase):
    def release_outcome(self, rel, **att):
        return run(record(attempts=[attempt(rel, **att)]))["token_release"]

    def test_released_on_evidence(self):
        out = run(record())
        self.assertEqual(dict((k, v[0]) for k, v in out.items()),
                         {"stop_answer": "pass", "job_terminal": "pass", "attempt_state": "pass",
                          "token_release": "pass"}, out)
        self.assertIn("finished (site stop)", out["token_release"][1])
        res = checks.evaluate(scenario(), {"http": 200, "body": {"stop_flow": record()}}, "mcp")
        self.assertEqual(checks.scenario_outcome(res), "pass")
        self.assertEqual(self.release_outcome(released("tombstone", site="worker"))[0], "pass")
        # R3-2: on a healthy host the worker's completed on the host's answer, decision null (rehearsal S2's record)
        completed = released("completed", site="worker", decision=None, detail="response")
        self.assertEqual(self.release_outcome(completed)[0], "pass")

    def test_holds_and_the_sweep_are_inconclusive(self):
        # R3-2: a regression to "always held to the lease" blocks (inconclusive, never a pass)
        for rel in (released("lease", decision="held", detail="stopping", at=LEASE),
                    released("lease", decision="would_release", detail="tombstone_would_release_ok", at=LEASE),
                    released("lease", decision="after_lease", detail="h4_finished", at="2026-10-20T10:01:01+00:00"),
                    released("swept", site="sweep", detail="clock_fallback"),
                    released("swept", site="sweep", decision=None, detail="finished")):
            out = self.release_outcome(rel)
            self.assertEqual(out[0], "inconclusive", (rel, out))
        self.assertEqual(self.release_outcome(released("grace", decision=None, detail="grace_20s"))[0], "fail")
        self.assertEqual(self.release_outcome(released("finished", decision="held"))[0], "fail")
        self.assertEqual(self.release_outcome(released("whatever"))[0], "fail")

    def test_asyncs_null_decision_rule(self):
        # RELEASE-RECORD.md "Settled": decision null exactly for completed and grace; another value there fails
        out = self.release_outcome(released("completed", site="worker", decision="released", detail="response"))
        self.assertEqual(out[0], "fail")
        self.assertIn("async's rule: null", out[1])

    def test_no_later_than_the_lease_with_asyncs_margin(self):
        late = released("finished", detail="finished", at="2026-10-20T10:01:20+00:00")   # lease + 20 s: within 30 s
        self.assertEqual(self.release_outcome(late)[0], "pass")
        later = released("finished", detail="finished", at="2026-10-20T10:01:31+00:00")
        out = self.release_outcome(later)
        self.assertEqual(out[0], "fail")
        self.assertIn("after the token's lease (lease_until + 30 s", out[1])
        # the row's token_lease_ms, where stated, is the bound
        out = self.release_outcome(late, token_lease_ms=70000)
        self.assertEqual(out[0], "fail", out)
        self.assertIn("dispatched_at + token_lease_ms 70000", out[1])
        self.assertEqual(self.release_outcome(released("finished", detail="finished", at=LEASE),
                                              token_lease_ms=70000)[0], "pass")
        self.assertEqual(self.release_outcome(dict(released("finished"), decision="error",
                                                   detail="release_failed"))[0], "fail")

    def test_a_null_release_is_inconclusive(self):
        out = self.release_outcome(None)
        self.assertEqual(out[0], "inconclusive")
        self.assertIn("release null on a terminal job", out[1])
        out = self.release_outcome(None, release_stated=False)
        self.assertEqual(out[0], "inconclusive")
        self.assertIn("older than the release record", out[1])
        # one attempt fails, another is null: the fail wins
        out = run(record(attempts=[attempt(None), attempt(released("grace", detail="grace_20s"), attempt_id="a2")]))
        self.assertEqual(out["token_release"][0], "fail")
        res = checks.evaluate(scenario(), {"http": 200, "body": {"stop_flow": record(attempts=[attempt(None)])}},
                              "mcp")
        self.assertEqual(checks.scenario_outcome(res), "inconclusive")
        out = run(record(admin="not read: no admin key"))
        self.assertEqual((out["token_release"][0], out["attempt_state"][0]), ("inconclusive", "inconclusive"))
        self.assertEqual(run(record(attempts=[]))["token_release"][0], "inconclusive")

    def test_a_race_is_inconclusive_on_every_check(self):
        for rec in (record(race="terminal_before_running"), record(race="stop_already_terminal"),
                    record(finished_at="2026-10-20T10:00:01+00:00")):
            out = run(rec)
            self.assertEqual(set(v[0] for v in out.values() if v), {"inconclusive"}, out)
            self.assertIn("race", " ".join(v[1] for v in out.values()))

    def test_the_stop_answer_the_end_and_the_attempt(self):
        out = run(record(stop_answer={"running": 1, "cancel_supported": False, "bound_seconds": 58}))
        self.assertEqual(out["stop_answer"][0], "fail")
        self.assertIn("cancel_supported False, expected True", out["stop_answer"][1])
        self.assertEqual(run(record(stop={"sent_utc": STOP_SENT, "answered_utc": STOP_ANSWERED, "answer": None,
                                          "error": "forbidden"}))["stop_answer"][0], "fail")
        out = run(record(finished_at="2026-10-20T10:01:30+00:00", stop_answer={
            "running": 1, "cancel_supported": True, "bound_seconds": 10}))
        self.assertEqual(out["job_terminal"][0], "fail")
        self.assertIn("its bound 10 s", out["job_terminal"][1])
        self.assertEqual(run(record(status="failed"))["job_terminal"][0], "fail")
        self.assertEqual(run(record(attempts=[attempt(released("finished"), settled_at=None)]))["attempt_state"][0],
                         "fail")
        out = run(record(attempts=[attempt(released("finished"), server_instance="s2")]))
        self.assertIn("ran on s2, reserved for s1", out["attempt_state"][1])

    def test_the_direct_arms_do_not_expose_them(self):
        res = checks.evaluate(scenario(), {"http": 200, "body": {"walks": []}}, "direct-local")
        got = dict((r["check"], r["outcome"]) for r in res)
        for name in checks.STOP_FLOW_CHECKS:
            self.assertEqual(got[name], "not_exposed", name)
        # on the MCP arm an answer without the record: not exposed (checks.Missing)
        res = checks.evaluate(scenario(), {"http": 200, "body": {}}, "mcp")
        self.assertEqual(set(r["outcome"] for r in res if r["check"] in checks.STOP_FLOW_CHECKS), {"not_exposed"})


class ResolveStop(unittest.TestCase):
    ANSWER = {"ok": True, "running": 1, "cancel_supported": True, "cancel": "disconnect", "bound_seconds": 58}

    def rec(self, rel, **att):
        return record(tool="traverse_resolve", stop_answer=dict(self.ANSWER),
                      attempts=[attempt(rel, settled_by="disconnect", **att)])

    def test_grace_20s_after_the_close(self):
        out = run(self.rec(released("grace", decision=None, detail="grace_20s", at="2026-10-20T10:00:22.6+00:00")),
                  "resolve")
        self.assertEqual(dict((k, v[0]) for k, v in out.items() if k in checks.STOP_FLOW_CHECKS),
                         {"stop_answer": "pass", "job_terminal": "pass", "attempt_state": "pass",
                          "token_release": "pass"}, out)
        out = run(self.rec(released("grace", decision=None, detail="grace_20s", at="2026-10-20T10:00:45+00:00")),
                  "resolve")
        self.assertEqual(out["token_release"][0], "fail")
        self.assertIn("expected the close + 20 s", out["token_release"][1])
        out = run(self.rec(released("grace", decision=None, detail="grace_20s", at="2026-10-20T10:00:05+00:00")),
                  "resolve")
        self.assertEqual(out["token_release"][0], "fail")
        out = run(self.rec(released("grace", decision=None, detail="grace_0", at="2026-10-20T10:00:03+00:00")),
                  "resolve")
        self.assertIn("expected 'grace_20s'", out["token_release"][1])
        out = run(self.rec(released("grace", decision="released", detail="grace_20s",
                                    at="2026-10-20T10:00:22.6+00:00")), "resolve")
        self.assertIn("async's rule: null", out["token_release"][1])   # a decision on grace fails

    def test_the_lease_first_and_a_completed_resolve(self):
        lease = "2026-10-20T10:00:15+00:00"   # before the close + 20 s: capped
        rel = released("lease", decision="held", detail="grace_capped", at=lease)
        self.assertEqual(run(self.rec(rel, lease_until=lease), "resolve")["token_release"][0], "pass")
        out = run(self.rec(dict(rel, at=LEASE), lease_until=LEASE), "resolve")   # the lease ends long after
        self.assertEqual(out["token_release"][0], "fail")
        self.assertIn("grace_capped although the lease", out["token_release"][1])
        # R3-2: a completed /resolve ended before the stop (a race): inconclusive; the sweep likewise
        self.assertEqual(run(self.rec(released("completed", site="worker", decision=None, detail="answer")),
                             "resolve")["token_release"][0], "inconclusive")
        self.assertEqual(run(self.rec(released("swept", site="sweep", detail="finished")),
                             "resolve")["token_release"][0], "inconclusive")
        self.assertEqual(run(self.rec(released("finished", detail="finished")), "resolve")["token_release"][0],
                         "fail")

    def test_the_resolve_stop_answer(self):
        out = run(record(tool="traverse_resolve", stop_answer=dict(self.ANSWER, cancel_supported=False)), "resolve")
        self.assertEqual(out["stop_answer"][0], "fail")


class BuiltScenarios(unittest.TestCase):
    """VERIFY-D VD5: the BUILT stop scenarios' expect lists (areas.gen_stop on v7's scenarios) evaluated on a stop_flow
    body in the runner's recorded shape (suite_stop_flow.py a52bd3d6: {"http": 200, "body": {"stop_flow": ...}}): no
    check judges the stop body as a traversal or a resolve answer; a correct stop passes, a race is inconclusive."""

    @classmethod
    def setUpClass(cls):
        import suite_common as sc
        v7 = sc.read_json(os.path.join(os.path.dirname(HERE), "manifest", "suite_v7.json"))["scenarios"]

        class Ctx(object):
            def skip(self, *a):
                raise AssertionError("skipped %s" % (a,))
        cls.built = dict((s["id"], s) for s in areas.gen_stop(Ctx(), v7))

    def evaluate(self, sid, rec):
        res = checks.evaluate(self.built[sid], {"http": 200, "body": {"stop_flow": rec}}, "mcp")
        return res, checks.scenario_outcome(res)

    def test_a_correct_stop_passes_on_every_built_check(self):
        res, outcome = self.evaluate("S-stop-0001", record())
        self.assertEqual(outcome, "pass", res)
        self.assertEqual([r["check"] for r in res], ["status", "stop_answer", "job_terminal", "attempt_state",
                                                     "token_release"])
        rel = released("grace", detail="grace_20s", at="2026-10-20T10:00:22.6+00:00")
        rec = record(tool="traverse_resolve", stop_answer=dict(ResolveStop.ANSWER),
                     attempts=[attempt(rel, settled_by="disconnect")])
        res, outcome = self.evaluate("S-stop-0002", rec)
        self.assertEqual(outcome, "pass", res)
        self.assertFalse([r for r in res if r["check"] in ("invariants", "answer_fields", "budget_respected")])

    def test_a_race_or_a_null_release_is_inconclusive_never_a_fail(self):
        for sid in ("S-stop-0001", "S-stop-0002"):
            self.assertEqual(self.evaluate(sid, record(race="terminal_before_running"))[1], "inconclusive", sid)
        self.assertEqual(self.evaluate("S-stop-0001", record(attempts=[attempt(None)]))[1], "inconclusive")


class TakenAndTimed(unittest.TestCase):
    """VERIFY-D VD7: already_stopping (an earlier call of the same stop took: async's STOP_TOOK) counts as taken; a
    stop never answered is inconclusive. Async's bundle review (2026-10-20): L1 the token's hold past lease_until is
    at most 30 s on a contract host, 25 s legacy, 0 for /resolve; M1 a /traverse release may be written well after the
    job read stopped (within the lease bound it is on time)."""

    def test_already_stopping_counts_as_taken(self):
        stop = {"sent_utc": STOP_SENT, "answered_utc": STOP_ANSWERED, "error": "already_stopping", "calls": 2,
                "took": True, "answer": {"error": {"code": "already_stopping", "message": "a stop was requested"}}}
        out = run(record(stop=stop))
        self.assertEqual(out["stop_answer"][0], "pass", out)
        self.assertIn("an earlier call of this stop took (2 call(s))", out["stop_answer"][1])
        # R2-4: only what the answer states is checked; the rest is named unseen
        self.assertIn("not seen: bound_seconds_gt, cancel_supported, running", out["stop_answer"][1])
        stated = dict(stop, answer={"error": {"code": "already_stopping"}, "running": 0})
        self.assertEqual(run(record(stop=stated))["stop_answer"][0], "fail")    # a stated running 0, expected 1
        # its bound unseen: the job's end is judged against the attempts' lease (LEASE 10:01:00)
        self.assertEqual(out["job_terminal"][0], "pass", out)
        late = run(record(stop=stop, finished_at="2026-10-20T10:01:30+00:00"))
        self.assertEqual(late["job_terminal"][0], "fail", late)
        # R2-4: a stop not taken is an error, whatever the reason
        for err in ("temporarily_unavailable", "timeout", "forbidden"):
            out = run(record(stop=dict(stop, error=err, took=False, answer=None)))
            self.assertEqual(out["stop_answer"][0], "fail", (err, out))
            self.assertIn("the stop was not taken: %s" % err, out["stop_answer"][1])
        out = run(record(stop=dict(stop, error=None, took=False, answer={"running": 1})))
        self.assertEqual(out["stop_answer"][0], "fail")

    def test_the_tokens_hold_past_the_lease(self):
        at = "2026-10-20T10:01:27+00:00"     # lease_until + 27 s
        rel = released("finished", detail="finished", at=at)
        self.assertEqual(run(record(attempts=[attempt(rel)]))["token_release"][0], "pass")   # contract: 30 s
        out = run(record(attempts=[attempt(rel, bound_ms=None, not_after_ms=None)]))        # legacy: 25 s
        self.assertEqual(out["token_release"][0], "fail", out)
        self.assertIn("lease_until + 25 s (legacy)", out["token_release"][1])
        # /resolve records before publishing: no hold past its lease
        lease = "2026-10-20T10:00:21+00:00"
        rel = released("grace", detail="grace_20s", at="2026-10-20T10:00:22.4+00:00")
        out = run(record(tool="traverse_resolve", stop_answer=dict(ResolveStop.ANSWER),
                         attempts=[attempt(rel, lease_until=lease)]), "resolve")
        self.assertEqual(out["token_release"][0], "fail", out)
        self.assertIn("lease_until + 0 s (resolve)", out["token_release"][1])

    def test_a_release_written_long_after_the_end_is_on_time_within_the_lease(self):
        # M1: the job read stopped at 10:00:04; the stop site wrote the release 40 s later, before lease + 30 s
        rel = released("finished", detail="finished", at="2026-10-20T10:00:44+00:00")
        self.assertEqual(run(record(attempts=[attempt(rel)]))["token_release"][0], "pass")


class StopScenarioShape(unittest.TestCase):
    def test_the_flow_states_only_asyncs_keys(self):
        flow_keys = {"kind", "after", "delay_s", "observe", "poll_s", "release_wait_s"}   # suite_stop_flow.FLOW_KEYS
        for checks_list in (areas.STOP_CHECKS_TRAVERSE, areas.STOP_CHECKS_RESOLVE):
            self.assertEqual([e["check"] for e in checks_list], list(checks.STOP_FLOW_CHECKS))
        import inspect
        src = inspect.getsource(areas.gen_stop)
        self.assertNotIn('"expect":', src)
        self.assertTrue(flow_keys >= {"kind", "after", "delay_s", "observe"})


if __name__ == "__main__":
    unittest.main()
