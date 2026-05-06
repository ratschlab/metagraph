"""Portable subprocess timing wrapper for workflow shell commands."""

import argparse
import resource
import subprocess
import sys
import time


def _rss_kb(ru_maxrss: int) -> int:
    # On Linux ru_maxrss is in KB, on macOS it's in bytes.
    if sys.platform == "darwin":
        return int(ru_maxrss / 1024)
    return int(ru_maxrss)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    ns = parser.parse_args(argv)

    cmd = ns.command
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if not cmd:
        print("time_wrapper: missing command", file=sys.stderr)
        return 2

    start = time.perf_counter()
    try:
        child = subprocess.Popen(cmd)
    except OSError as e:
        print(f"[timing] failed_to_exec={e}", file=sys.stderr)
        return 127
    exit_code = child.wait()
    elapsed = time.perf_counter() - start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)

    # Print concise, portable metrics footer.
    print(
        "[timing] wall_sec={:.3f} user_sec={:.3f} sys_sec={:.3f} max_rss_kb={} exit_code={}".format(
            elapsed,
            usage.ru_utime,
            usage.ru_stime,
            _rss_kb(usage.ru_maxrss),
            exit_code,
        ),
        file=sys.stderr,
    )
    return exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
