"""The fake GraphletTools: graphlet_walks(mode=...) misbehaves as its mode says (see the package docstring)."""
import os
import sys
import time

_NAMES = ["graphlet_load", "graphlet_summary", "graphlet_walks"]


def tool_names():
    return list(_NAMES)


class ToolLimits(object):
    def __init__(self, **kw):
        pass


class GraphletTools(object):
    def __init__(self, store, secret=None, export_dir=None, max_bytes=2048, sequence_max_bytes=16384,
                 local_limits=None):
        if max_bytes == 77:   # a library that does not start (the tests' way that fails)
            raise RuntimeError("this fake library refuses max_bytes 77")
        self.loaded = 0
        self.handles = set()

    def _unknown(self, handle):
        return {"error": "unknown_handle", "message": "no handle %s" % handle, "replayable": False}

    def graphlet_load(self, path, budget=None):
        self.loaded += 1
        handle = "g_%012x" % (0xfa4e000 + self.loaded)
        self.handles.add(handle)
        if "nanload" in os.path.basename(path):   # a receipt holding a NaN (VERIFY-L VR2)
            return {"handle": handle, "summary": {"seed": {"length_bp": float("nan")}}}
        return {"handle": handle, "summary": None}

    def graphlet_summary(self, handle, arm=None, detail="summary", cursor=None, max_bytes=None, budget=None):
        if handle not in self.handles:
            return self._unknown(handle)
        return {"handle": handle, "summary": {"loaded": self.loaded}}

    def graphlet_walks(self, handle, arm=None, mode="ok", cursor=None, max_bytes=None):
        if handle not in self.handles:
            return self._unknown(handle)
        if mode == "sleep":
            time.sleep(600)
        if mode == "set":
            return {"handle": handle, "rows": [{"labels": {1, 2}}], "total": 1}
        if mode == "die":
            os._exit(3)
        if mode == "nan":
            return {"handle": handle, "rows": [{"cost": 1.5}, {"cost": float("nan")}], "total": 2}
        if mode == "inf":
            return {"handle": handle, "rows": [], "total": 0, "limits": {"loss": float("inf")}}
        if mode == "loop":
            return {"handle": handle, "rows": [1], "total": 2, "next_cursor": "the-same-cursor"}
        if mode == "print":
            print("fake library noise on stdout")
            sys.stdout.flush()
            os.write(1, b"fake library noise on fd 1\n")
        if mode == "stdin":
            assert sys.stdin.read() == "" and os.read(0, 16) == b""
        return {"handle": handle, "rows": [1], "total": 1}
