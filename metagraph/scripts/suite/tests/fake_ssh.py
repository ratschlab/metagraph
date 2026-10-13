#!/usr/bin/env python3
"""A stand-in for ssh, for the orchestrator's offline tests of the ssh path (SUITE_MEX=user@host, with a
directory holding this script as `ssh` first on PATH). It runs the remote command with `bash -c` in
FAKE_SSH_HOME (mex's home directory), as ssh runs it in the login shell, and appends one JSON line per call
to FAKE_SSH_LOG: {pid, op: run | exit | check | down, control_path, command, killed}.

Multiplexing as far as the tests need it: a run call registers itself under its ControlPath (the literal
option value) while it runs; `-O exit` on a path ends every call registered under it with exit 255, as the
master's end cuts its sessions, while the remote command keeps running (it meets a closed pipe at its next
write). While the file FAKE_SSH_HOME/.ssh-down exists every call fails at once with exit 255 (no route to
mex); a call already running ends with 255 as soon as the file appears (a dropped connection). rsync's
remote side (`rsync --server ...`) is executed directly, with stdin and stdout passed through.
"""

import json
import os
import signal
import subprocess
import sys
import threading
import time

OPTS_WITH_ARG = set("bcDEeFIiJLlmOoPpQRSWw")


def parse(argv):
    """(options {letter: [values]}, -o values {name: value}, host, command words)."""
    opts, o_opts, i = {}, {}, 0
    while i < len(argv):
        a = argv[i]
        if not a.startswith("-") or a == "-":
            break
        letter = a[1:2]
        if letter in OPTS_WITH_ARG:
            if len(a) > 2:
                value = a[2:]
            else:
                i += 1
                value = argv[i] if i < len(argv) else ""
            opts.setdefault(letter, []).append(value)
            if letter == "o" and "=" in value:
                name, v = value.split("=", 1)
                o_opts[name.strip().lower()] = v.strip()
        else:
            opts.setdefault(letter, []).append(True)
        i += 1
    host = argv[i] if i < len(argv) else ""
    return opts, o_opts, host, argv[i + 1:]


def state_dir(home, control_path):
    safe = "".join(c if c.isalnum() else "_" for c in control_path or "none")
    d = os.path.join(home, ".fake-ssh", safe)
    os.makedirs(d, exist_ok=True)
    return d


def log(entry):
    path = os.environ.get("FAKE_SSH_LOG")
    if path:
        entry = dict(entry, pid=os.getpid(), t=time.time())
        with open(path, "a") as fh:
            fh.write(json.dumps(entry) + "\n")


def main(argv):
    home = os.environ["FAKE_SSH_HOME"]
    down = os.path.join(home, ".ssh-down")
    opts, o_opts, host, words = parse(argv)
    control_path = o_opts.get("controlpath")
    command = " ".join(words)
    if os.path.exists(down):
        log({"op": "down", "control_path": control_path, "command": command})
        sys.stderr.write("ssh: connect to host %s port 22: Network is unreachable\n" % host)
        return 255
    if "O" in opts:
        op = opts["O"][-1]
        killed = []
        if op == "exit":
            d = state_dir(home, control_path)
            for name in os.listdir(d):
                try:
                    os.kill(int(name), signal.SIGTERM)
                    killed.append(int(name))
                except (ValueError, ProcessLookupError):
                    pass
        log({"op": op, "control_path": control_path, "command": "", "killed": killed})
        return 0
    if words[:2] == ["rsync", "--server"]:
        log({"op": "run", "control_path": control_path, "command": command})
        os.chdir(home)
        os.execvp("bash", ["bash", "-c", command])
    reg = os.path.join(state_dir(home, control_path), str(os.getpid()))
    open(reg, "w").close()

    def cut(*_):
        try:
            os.unlink(reg)
        except OSError:
            pass
        os._exit(255)
    signal.signal(signal.SIGTERM, cut)
    log({"op": "run", "control_path": control_path, "command": command})
    try:
        child = subprocess.Popen(["bash", "-c", command], cwd=home, stdout=subprocess.PIPE,
                                 start_new_session=True)

        def watch_down():
            while child.poll() is None:
                if os.path.exists(down):
                    log({"op": "dropped", "control_path": control_path, "command": command})
                    cut()
                time.sleep(0.05)
        threading.Thread(target=watch_down, daemon=True).start()
        out = sys.stdout.buffer
        while True:
            chunk = child.stdout.read1(65536) if hasattr(child.stdout, "read1") else child.stdout.read(65536)
            if not chunk:
                break
            out.write(chunk)
            out.flush()
        return child.wait()
    finally:
        try:
            os.unlink(reg)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
