"""Synthetic area G runs (manifest v8, FORMAT-v8 §1-§5) for the tests of lane C: a run directory with a manifest of
graphlet scenarios, the five ways' per-call records, the documents the ways kept, the `detail: full` twin, the MCP
echo and the runner's identity records. The documents are real: the mini_refseq graphlets of the library's fixtures
(engine answers recorded once), spliced into the standalone form (H, J, body, Z) without the library.

The digests are graphlet_doc's (lane L's reader, cherry-picked); install_graphlet_doc_stub falls back to a stub of
FORMAT-v8 §3's definitions (body_digest: the sha256 of the body records, H's seed_index masked, J and Z left out;
envelope_digest: the J line's canonical JSON without its top-level timing and usage) where that module is absent."""

import copy
import gzip
import hashlib
import json
import os
import shutil
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
MINI = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python", "tests", "data", "traverse", "real",
                                     "mini_refseq"))
LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))

PIN = "a387dfdd0123456789abcdef0123456789abcdef"
SERVICE = "26a52cb63efd95752a109b5e09b1163bfffe60ca"
WAY_DIRS = {"W1": "dl-lib", "W2": "dr-lib", "W3": "mcp", "W4": "mcp-lib", "W5": "oracle", "drift": "drift-lib"}


def _stub_body_digest(text):
    import graphlet_precheck as gp
    return hashlib.sha256("\n".join(gp.body_records(text)).encode("utf-8")).hexdigest()


def _stub_envelope_digest(text):
    import graphlet_precheck as gp
    env = gp.envelope_of(text)
    if env is None:
        return None
    env = dict((k, v) for k, v in env.items() if k not in ("timing", "usage"))
    return hashlib.sha256(json.dumps(env, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def install_graphlet_doc_stub():
    """Put a graphlet_doc stub in sys.modules (and into graphlet_precheck) unless lane L's module loads; returns the
    module in use."""
    try:
        import graphlet_doc  # noqa: F401
        if not getattr(graphlet_doc, "_lane_c_stub", False) and hasattr(graphlet_doc, "body_digest"):
            mod = graphlet_doc
        else:
            raise ImportError
    except ImportError:
        mod = types.ModuleType("graphlet_doc")
        mod.body_digest = _stub_body_digest
        mod.envelope_digest = _stub_envelope_digest
        mod._lane_c_stub = True
        sys.modules["graphlet_doc"] = mod
    import graphlet_precheck as gp
    gp.graphlet_doc = mod
    return mod


class OracleStub(object):
    """A stand-in for lane O's graphlet_oracle.compare (FORMAT v2 F2) until its commit lands: a test double of its
    documented semantics on the cases the tests use (checked pointers equal as JSON values, a checked pointer the
    other answer lacks, expected_error requiring an error), recording every call it gets."""

    def __init__(self, result=None):
        self.calls = []
        self.result = result   # a fixed result to return instead (a list, or {outcome, differences})

    def compare(self, record, other):
        import graphlet_precheck as gp
        self.calls.append((record, other))
        if self.result is not None:
            if isinstance(self.result, Exception):
                raise self.result
            return self.result
        o = record.get("oracle") or {}
        if o.get("declared") == "expected_error":
            return [] if isinstance(other, dict) and "error" in other else [
                {"pointer": "/error", "clause": "an error expected", "oracle": "an error", "missing": True}]
        out = []
        for ptr, clause in sorted((o.get("checked") or {}).items()):
            _, mine = gp.pointer_get(record.get("answer"), ptr)
            found, theirs = gp.pointer_get(other, ptr)
            if not found:
                out.append({"pointer": ptr, "clause": clause, "oracle": mine, "missing": True})
            elif gp.first_difference(mine, theirs):
                out.append({"pointer": ptr, "clause": clause, "oracle": mine, "other": theirs})
        return out


def install_oracle_stub(result=None):
    """Point compare_graphlets at an OracleStub (returned) unless lane O's graphlet_oracle with compare() loads."""
    import compare_graphlets as cg
    real = sys.modules.get("graphlet_oracle")
    if real is not None and callable(getattr(real, "compare", None)) and result is None:
        cg.graphlet_oracle = real
        return real
    stub = OracleStub(result)
    cg.graphlet_oracle = stub
    return stub


def load_fixture(cell, kind):
    with gzip.open(os.path.join(MINI, "%s.%s.json.gz" % (cell, kind))) as fh:
        return json.load(fh)


def standalone(body, envelope, seed_index=None):
    """An MGT document in the standalone form: H (its seed_index set to |seed_index| when given), the J line of
    |envelope| (None: none), the body records, Z counting the lines."""
    lines = body.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    head, rest = lines[0], lines[1:]
    if rest and rest[-1].startswith("Z "):
        rest = rest[:-1]
    if seed_index is not None:
        f = head.split(" ")
        f[11] = str(seed_index)
        head = " ".join(f)
    out = [head]
    if envelope is not None:
        out.append("J " + json.dumps(envelope, sort_keys=True, separators=(",", ":"), ensure_ascii=False))
    out.extend(rest)
    out.append("Z %d" % (len(out) + 1))
    return "\n".join(out) + "\n"


class Fixture(object):
    """One graphlet of the mini_refseq fixtures: its body, the engine's envelope, the full twin."""

    def __init__(self, cell="mini_ndm1__trace"):
        self.cell = cell
        self.graphlet = load_fixture(cell, "graphlet")
        self.full = load_fixture(cell, "full")
        self.request = load_fixture(cell, "request")
        r = self.graphlet["results"][0]
        self.body = r["graphlet"]
        env = dict((k, v) for k, v in self.graphlet.items() if k != "results")
        env["results"] = [dict((k, v) for k, v in r.items() if k not in ("graphlet", "graphlet_bytes",
                                                                         "graphlet_lines"))]
        env["timing"] = {"elapsed_ms": 12.5}
        self.envelope = env

    def doc(self, seed_index=None, envelope=True, **env_changes):
        env = copy.deepcopy(self.envelope)
        for k, v in env_changes.items():
            env[k] = v
        return standalone(self.body, env if envelope else None, seed_index)


def strategy(detail="graphlet"):
    return {"bounds": {"max_extension_bp": 1000, "time_budget_ms": 15000, "max_work_units": 100000000},
            "labels": {"mode": "constrain", "max_labels_per_node": 64}, "direction": "both",
            "output": {"detail": detail, "timing": True}}


def scenario(sid="G-merge-0001", ops=None, ways=None, pre=None, pin=PIN, **graphlet_extra):
    g = {"twin": {"strategy": strategy("full")},
         "preconditions": pre if pre is not None else {"arms": ["left", "right"]},
         "ops": ops if ops is not None else [
             {"id": "walks-len", "tool": "graphlet_walks",
              "args": {"arm": "left", "rank": "length", "n": 1, "max_bytes": 2048}, "paging": "all"}],
         "ways": {"default": ways or ["W1", "W2", "W3", "W4", "W5"]}}
    if pin:
        g["library_commit"] = pin
    g.update(graphlet_extra)
    return {"id": sid, "area": "G", "route": "POST /traverse", "title": "graphlet %s" % sid,
            "request": {"seeds": [{"sequence": "ACGT"}], "strategy": strategy()},
            "mcp": {"tool": "traverse_fetch", "preset": "walk", "overrides": {"bounds": {"max_extension_bp": 1000}}},
            "expect": [{"check": "status", "http": 200}], "graphlet": g}


def manifest(scenarios, relations=None, **top):
    import suite_common as sc
    m = {"manifest_version": 8, "scenarios": scenarios, "relations": relations or [],
         "index": {"database": "refseq33m-experimental", "index_fp": "fp-test"}}
    m.update(top)
    ident = sc.manifest_identity(m)
    m["content_sha256"] = ident["content_sha256"]
    m["relations_sha256"] = ident["relations_sha256"]
    return m


class GRun(object):
    """A synthetic run directory for area G (call cleanup() when done)."""

    def __init__(self, scenarios=None, fixture=None, **manifest_top):
        self.dir = tempfile.mkdtemp(prefix="lane_c_run_")
        self.fixture = fixture or Fixture()
        self.scenarios = scenarios if scenarios is not None else [scenario()]
        self.manifest = manifest(self.scenarios, **manifest_top)
        self.write_json("manifest.json", self.manifest)

    @classmethod
    def attach(cls, run_dir, manifest, fixture=None):
        """A GRun writing into an existing run directory whose manifest.json is |manifest| (not rewritten)."""
        self = cls.__new__(cls)
        self.dir = run_dir
        self.fixture = fixture or Fixture()
        self.manifest = manifest
        self.scenarios = [s for s in manifest["scenarios"] if isinstance(s.get("graphlet"), dict)]
        return self

    def cleanup(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def path(self, *parts):
        p = os.path.join(self.dir, *parts)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        return p

    def write_json(self, rel, obj):
        with open(self.path(rel), "w") as fh:
            json.dump(obj, fh, sort_keys=True)

    def append(self, rel, obj):
        with open(self.path(rel), "a") as fh:
            fh.write(json.dumps(obj, sort_keys=True) + "\n")

    def write_doc(self, way, sid, text):
        with gzip.open(self.path("arms", WAY_DIRS[way], "docs", sid + ".mgt.gz"), "wb") as fh:
            fh.write(text.encode("utf-8"))

    def write_twin(self, sid, body=None):
        with gzip.open(self.path("arms", "direct-local", "docs", sid + ".full.json.gz"), "wt") as fh:
            json.dump(body if body is not None else self.fixture.full, fh)

    def write_echo(self, sid, overrides):
        self.append(os.path.join("arms", "mcp", "graphlet_jobs.jsonl"), {"scenario": sid, "overrides": overrides})

    def write_mcp_runs(self, library_commit=PIN, service_commit=SERVICE, n=1):
        for _ in range(n):
            self.append(os.path.join("arms", "mcp", "mcp_runs.jsonl"),
                        {"service_commit": service_commit, "identity_changed": False,
                         "library_commit": library_commit})

    def write_file(self, way, raw, text):
        """Keeps a file's text where a record's file.raw names it (FORMAT v2.1 B3: relative to the way's directory,
        "raw/<name>") and returns its `file` block."""
        path = self.path("arms", WAY_DIRS[way], raw)
        if raw.endswith(".gz"):
            with gzip.open(path, "wb") as fh:
                fh.write(text.encode("utf-8"))
        else:
            with open(path, "w") as fh:
                fh.write(text)
        return {"name": os.path.basename(raw), "bytes": len(text.encode("utf-8")),
                "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(), "raw": raw}

    def record(self, way, sid, op, tool, args, answer=None, error=None, service_fields=None, doc=None,
               library=None, selector_error=None, oracle=None, pages=1, file=None):
        """Appends one per-call record (FORMAT-v8 §3) to arms/<way>/graphlet_ops.jsonl."""
        if library is None:
            library = None if way == "W5" else ({"commit": PIN} if way == "W3" else
                                                {"commit": PIN, "python": "3.13.1", "hashseed": 17, "work_model": 2})
        rec = {"run_id": "test", "scenario": sid, "way": way, "op": op, "tool": tool,
               "args_resolved": dict(args, handle="<handle>"), "selector_error": selector_error, "pages": pages,
               "answer": answer if error is None else None, "error": error,
               "service_fields": service_fields if way == "W3" else None, "doc": doc, "library": library,
               "elapsed_ms": 1.0, "started_utc": "2026-10-19T00:00:00Z"}
        if oracle is not None:
            rec["oracle"] = oracle
        if file is not None:
            rec["file"] = file
        self.append(os.path.join("arms", WAY_DIRS[way], "graphlet_ops.jsonl"), rec)
        return rec

    def documents(self, sid, texts=None):
        """Writes W1/W2/W4's documents (|texts| {way: text}, default the fixture's with seed_index 0 / 3 / 0) and
        returns {way: {body_digest, envelope_digest}} as the stub computes them."""
        mod = install_graphlet_doc_stub()
        texts = texts or {"W1": self.fixture.doc(0), "W2": self.fixture.doc(3), "W4": self.fixture.doc(0)}
        out = {}
        for way, text in texts.items():
            self.write_doc(way, sid, text)
            out[way] = {"body_digest": mod.body_digest(text), "envelope_digest": mod.envelope_digest(text)}
        return out

    def complete_scenario(self, sid="G-merge-0001", ops=None, docs=None, echo=True, twin=True, mcp_runs=True):
        """Every input of a scenario that passes its pre-checks: documents, twin, echo, identity, and a record per
        op and way with the same answer (|ops| [(op id, tool, args, answer)], default the scenario's single walks
        op with a small answer)."""
        digests = self.documents(sid, docs)
        if twin:
            self.write_twin(sid)
        scen = next(s for s in self.scenarios if s["id"] == sid)
        if echo:
            self.write_echo(sid, scen["request"]["strategy"])
        if mcp_runs:
            self.write_mcp_runs()
        ops = ops if ops is not None else [(op["id"], op["tool"], op["args"], walks_answer())
                                           for op in scen["graphlet"]["ops"]]
        for op_id, tool, args, answer in ops:
            for way in ("W1", "W2", "W3", "W4"):
                d = digests.get(way if way != "W3" else "W4")
                ans = copy.deepcopy(answer)
                sf = None
                if way == "W3":
                    ans = dict(ans, handle="0d5b9b1a-1111-2222-3333-444455556666")
                    if tool in ("graphlet_walks", "graphlet_walk"):
                        sf = {"walk_display": "first_parent"}
                self.record(way, sid, op_id, tool, args, answer=ans, doc=d, service_fields=sf)
        return digests


def walks_answer(handle="g_0123456789ab", rows=None):
    return {"handle": handle, "arm": "left",
            "evidence": {"complete_to_bp": {"left": 100}, "exact": True, "support": "kmer"},
            "total": 2, "rows": rows if rows is not None else [
                {"walk": 1, "length_bp": 120, "end": "dead_end", "n_alive": 3, "n_full": 2,
                 "labels_full": [{"name": "NC_000962.3", "ref": "c:0"}], "complete": True, "beyond_certified_bp": 0},
                {"walk": 0, "length_bp": 80, "end": "dead_end", "n_alive": 1, "n_full": 1,
                 "labels_full": [], "complete": True, "beyond_certified_bp": 0}]}
