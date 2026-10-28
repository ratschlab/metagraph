#!/usr/bin/env python3
"""Run one engine on one read set, resuming past reads that hit the per-read timeout.

usage: run_arm.py BENCH ENGINE TIMEOUT_SEC READS.fa OUT_PREFIX -- [metagraph align options]

Writes OUT_PREFIX.tsv (one line per read, see align_engine_bench.cpp; a read that
crashed the process gets status CRASH) and
OUT_PREFIX.time.json: wall seconds and max RSS summed/maxed over all segments
(each segment reloads the graph and annotation), number of segments.
"""
import json, re, subprocess, sys, time
bench, engine, timeout, reads, out = sys.argv[1:6]
opts = sys.argv[sys.argv.index('--') + 1:]
nreads = sum(1 for l in open(reads) if l.startswith('>'))
skip, segments, wall, rss, crashes, timeouts = 0, 0, 0.0, 0, 0, 0
# stop an arm after this many timed-out reads; the rest are reported SKIPPED
MAX_TIMEOUTS = 25
with open(out + '.tsv', 'w') as f:
    pass
while skip < nreads:
    segments += 1
    cmd = ['/usr/bin/time', '-v', bench, engine, timeout, str(skip), 'align'] + opts + [reads]
    with open(out + '.tsv', 'a') as f, open(out + '.stderr', 'a') as e:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        e.write(p.stderr)
    m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', p.stderr)
    parts = [float(x) for x in m.group(1).split(':')]
    wall += sum(v * 60 ** i for i, v in enumerate(reversed(parts)))
    rss = max(rss, int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', p.stderr).group(1)))
    done = sum(1 for _ in open(out + '.tsv'))
    if p.returncode not in (0, 3):
        # a crash (e.g. SIGSEGV): output is flushed per read, so read number
        # |done| crashed; record it and resume after it
        names = [l[1:].split()[0] for l in open(reads) if l.startswith('>')]
        with open(out + '.tsv', 'a') as f:
            f.write(f'{names[done]}\tCRASH\t0\t*\t*\t*\t*\t*\t*\texit {p.returncode}\n')
        done += 1
        crashes += 1
    if p.returncode == 0 and done < nreads:
        sys.exit(f'{engine} {reads}: only {done}/{nreads} reads')
    if p.returncode == 3:
        timeouts += 1
    if timeouts >= MAX_TIMEOUTS and done < nreads:
        names = [l[1:].split()[0] for l in open(reads) if l.startswith('>')]
        with open(out + '.tsv', 'a') as f:
            for name in names[done:]:
                f.write(f'{name}\tSKIPPED\t0\t*\t*\t*\t*\t*\t*\t*\n')
        break
    skip = done
json.dump({'wall_s': wall, 'max_rss_kb': rss, 'segments': segments, 'reads': nreads, 'crashes': crashes},
          open(out + '.time.json', 'w'))
