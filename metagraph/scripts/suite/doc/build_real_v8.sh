#!/bin/sh
# build_real_v8.sh RUN_DIR MANIFEST POLICY [ALSO_POLICY ...]
#
# The v8 results document of a real run, in one step, with the delivery bundle the v7 edition shipped
# (doc.pdf, evidence.json, gate_problems.json, the identity input when given, SHA256SUMS) and the build's checks.
#
#   RUN_DIR      the run directory (its manifest.json, records, re-evaluation if any)
#   MANIFEST     the frozen manifest file the run ran (manifest/suite_v8.json after the owner's freeze): its content
#                sha256 must be the run's; the unit mapping is built from it (doc/units_v8_build.py), so a freeze
#                that changed the manifest (a by-design case leaving, scenarios added) needs no edit here
#   POLICY       the deploy gate's policy (gate/policy_v8_rtp_mcp.json); ALSO_POLICY: further policies judged
#                offline and printed beside it (gate/policy_v8_rt_mcp.json); every policy must pin MANIFEST
#
# Environment (all optional):
#   OUT        the output folder (default ~/.cache/metagraph-traverse-plans/<UTC date>/doc-v8-real/<run id>); it
#              must not exist yet (a new folder per build: nothing is overwritten)
#   IDENTITY   the deployed service's identity input (make_doc --identity; copied into the bundle)
#   PREVIOUS   the previous edition's run (default: the v7 baseline run 20261019T005230Z_refseq97-1009f77c)
#   RECORDED=1 the recorded set instead of the re-evaluated one (make_doc --recorded)
#   XELATEX    default /Library/TeX/texbin/xelatex
#   ALLOW_DIRTY=1  build although the suite's tree has uncommitted changes (evidence.json's generator_commit is
#              the last commit of doc/, so a dirty tree is refused by default)
#
# The area G by-design table is the run's own snapshot (scripts/service_by_design_g.json, else the tree's), read by
# the comparison as the gate reads it: never a constant here. The declared file is the deploy policy's own.
# Long: ~5 min (the comparison is recomputed in memory, xelatex runs two or three times); run it in the background.
# Exit 0 built and every check clean, 1 built with a check to look at, 2 refused (inputs), 3 the build failed.
set -u
SUITE=$(cd "$(dirname "$0")/.." && pwd)
if [ $# -lt 3 ]; then
  sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
abspath() { case "$1" in /*) echo "$1";; *) echo "$(pwd)/$1";; esac; }
RUN=$(abspath "$1"); MAN=$(abspath "$2"); shift 2
POL=$(abspath "$1"); shift
ALSO=""
for p in "$@"; do ALSO="$ALSO $(abspath "$p")"; done
XELATEX=${XELATEX:-/Library/TeX/texbin/xelatex}
PY="nice -n 19 python3 -I -B"

[ -f "$RUN/manifest.json" ] || { echo "refused: $RUN has no manifest.json"; exit 2; }
for f in "$MAN" "$POL" $ALSO ${IDENTITY:+"$IDENTITY"}; do
  [ -f "$f" ] || { echo "refused: no file $f"; exit 2; }
done

# the inputs agree: the run ran MANIFEST, every policy pins it
python3 -I - "$RUN/manifest.json" "$MAN" "$POL" $ALSO <<'EOF' || exit 2
import json, sys
run_m, man = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
sha = man.get("content_sha256")
bad = []
if man.get("manifest_version") != 8:
    bad.append("%s: manifest version %s, not 8" % (sys.argv[2], man.get("manifest_version")))
if run_m.get("content_sha256") != sha:
    bad.append("the run's manifest is %s, the given manifest %s" % (run_m.get("content_sha256"), sha))
for p in sys.argv[3:]:
    pin = (json.load(open(p)).get("manifest") or {}).get("content_sha256")
    if pin != sha:
        bad.append("%s pins manifest %s, not %s" % (p, pin, sha))
if bad:
    print("refused:\n  " + "\n  ".join(bad))
    sys.exit(1)
print("manifest v8 content %s, %d scenarios: the run's and every policy's" % (sha, len(man.get("scenarios") or [])))
EOF

RUN_ID=$(python3 -I -c 'import json,os,sys; r=json.load(open(os.path.join(sys.argv[1],"run.json"))); print(r.get("run_id") or os.path.basename(sys.argv[1].rstrip("/")))' "$RUN")
OUT=$(abspath "${OUT:-$HOME/.cache/metagraph-traverse-plans/$(date -u +%Y-%m-%d)/doc-v8-real/$RUN_ID}")
[ -e "$OUT" ] && { echo "refused: $OUT exists (a new folder per build)"; exit 2; }

# the tree: evidence.json names the last commit of doc/ as the generator; a dirty tree would not be that commit
cd "$SUITE" || exit 2
TOP=$(git rev-parse --show-toplevel)
DIRTY=$(git status --porcelain -- . | grep -v '^?? ' || true)
if [ -n "$DIRTY" ] && [ "${ALLOW_DIRTY:-}" != 1 ]; then
  echo "refused: the suite's tree has uncommitted changes (ALLOW_DIRTY=1 to build anyway):"; echo "$DIRTY" | head -10
  exit 2
fi
if [ -z "${PREVIOUS:-}" ]; then
  V7=runs/20261019T005230Z_refseq97-1009f77c
  PREVIOUS=$SUITE/$V7
  [ -d "$PREVIOUS" ] || PREVIOUS=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")/metagraph/scripts/suite/$V7
fi
[ -d "$PREVIOUS" ] || { echo "refused: no previous run $PREVIOUS (PREVIOUS=...)"; exit 2; }

mkdir -p "$OUT/inputs" || exit 2
echo "run      $RUN_ID ($RUN)"
echo "tree     $(git rev-parse --short HEAD) ($TOP)${DIRTY:+ DIRTY}; doc/ last commit $(git log -1 --format=%h -- doc)"
echo "out      $OUT"

# the unit mapping keyed to MANIFEST (the committed one when it is the same file)
$PY doc/units_v8_build.py "$MAN" "$OUT/inputs/units_v8.json" > "$OUT/inputs/units.log" 2>&1 || {
  cat "$OUT/inputs/units.log"; echo "refused: the unit mapping cannot be built from $MAN"; exit 2; }
if cmp -s "$OUT/inputs/units_v8.json" doc/units/suite_v8.json; then
  UNITS=$SUITE/doc/units/suite_v8.json
  echo "units    doc/units/suite_v8.json (keyed to the given manifest)"
else
  UNITS=$OUT/inputs/units_v8.json
  echo "units    built from the given manifest: $UNITS (doc/units/suite_v8.json is keyed to another content:"
  echo "         commit this file as doc/units/suite_v8.json; the document is the same either way)"
fi

# the declared file the deploy policy pins (relative to the suite's tree, as gate/ names it)
DECL=$(python3 -I -c 'import json,sys; print((json.load(open(sys.argv[1])).get("declared") or {}).get("file") or "")' "$POL")
case "$DECL" in ""|/*) ;; *) DECL=$SUITE/$DECL;; esac
DECL_ARG=""
if [ -n "$DECL" ]; then
  [ -f "$DECL" ] || { echo "refused: the policy's declared file $DECL is missing"; exit 2; }
  DECL_ARG="--declared $DECL"
fi

set -- --out "$OUT/build" --policy "$POL" --units "$UNITS" --previous "$PREVIOUS" --xelatex "$XELATEX"
for p in $ALSO; do set -- "$@" --also-policy "$p"; done
[ -n "$DECL_ARG" ] && set -- "$@" $DECL_ARG
[ -n "${IDENTITY:-}" ] && set -- "$@" --identity "$(abspath "$IDENTITY")"
[ "${RECORDED:-}" = 1 ] && set -- "$@" --recorded
echo "make_doc $RUN $*" | sed "s#$HOME#~#g"
echo "build    started $(date -u +%H:%M:%SZ)"
$PY doc/make_doc.py "$RUN" "$@" > "$OUT/build.log" 2>&1
RC=$?
echo "build    ended $(date -u +%H:%M:%SZ), make_doc exit $RC"
if [ $RC -ne 0 ] || [ ! -f "$OUT/build/doc.pdf" ]; then
  tail -20 "$OUT/build.log"; exit 3
fi

# the bundle (as the v7 delivery): doc.pdf, evidence.json, gate_problems.json, the identity input, SHA256SUMS
cp -p "$OUT/build/doc.pdf" "$OUT/build/evidence.json" "$OUT/build/gate_problems.json" "$OUT/" || exit 3
FILES="doc.pdf evidence.json gate_problems.json"
if [ -n "${IDENTITY:-}" ]; then cp -p "$IDENTITY" "$OUT/" && FILES="$FILES $(basename "$IDENTITY")"; fi
(cd "$OUT" && shasum -a 256 $FILES > SHA256SUMS)

# the checks
B=$OUT/build
pdftotext -layout "$B/doc.pdf" "$B/doc.txt"
PROBLEMS=0
note() { echo "  CHECK: $*"; PROBLEMS=$((PROBLEMS + 1)); }
echo "checks"
PAGES=$(pdfinfo "$B/doc.pdf" | awk '/^Pages/{print $2}')
echo "  pages: $PAGES"
N=$(grep -c 'Missing character' "$B/doc.log"); echo "  missing characters: $N"; [ "$N" = 0 ] || note "missing characters"
N=$(grep -c -i 'undefined' "$B/doc.log"); echo "  undefined references: $N"; [ "$N" = 0 ] || note "undefined references"
N=$(grep -c 'multiply defined' "$B/doc.log"); echo "  multiply defined labels: $N"; [ "$N" = 0 ] || note "multiply defined labels"
N=$(grep -c '^Overfull' "$B/doc.log"); echo "  overfull boxes: $N"; grep '^Overfull' "$B/doc.log" | sed 's/^/    /'
# the known 1.1 pt box (the cover's title line) is the only one a clean build has
BIG=$(grep '^Overfull' "$B/doc.log" | awk -F'[(p]' '{ if ($2 + 0 > 1.2) n++ } END { print n + 0 }')
[ "$BIG" = 0 ] || note "$BIG overfull box(es) above 1.2 pt"
N=$(grep -c -E 'UNRESOLVED|to be written|\{[a-z_]+(\[[^]]*\])*(:[a-z]+)?\}' "$B/doc.txt")
echo "  placeholders / unresolved in the PDF text: $N"; [ "$N" = 0 ] || note "placeholders in the PDF text"
grep -n -E 'UNRESOLVED|to be written|\{[a-z_]+(\[[^]]*\])*(:[a-z]+)?\}' "$B/doc.txt" | head -5 | sed 's/^/    /'
N=$(grep -c 'WARNING' "$OUT/build.log"); echo "  make_doc warnings: $N"; [ "$N" = 0 ] || { note "make_doc warnings"; grep 'WARNING' "$OUT/build.log" | cut -c1-300 | sed 's/^/    /'; }
N=$(grep -c -E '/Users/|/home/|~/|/private/|/var/folders|raetsch' "$B/doc.txt")
echo "  local paths in the PDF text: $N"; [ "$N" = 0 ] || note "local paths in the PDF text"
python3 -I - "$RUN" "$OUT" "$MAN" "$UNITS" "$POL" $ALSO <<'EOF' || note "evidence, pins or verdicts"
import hashlib, json, os, sys
run, out, man, units = sys.argv[1:5]
policies = sys.argv[5:]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
e = json.load(open(os.path.join(out, "evidence.json")))
bad = 0
for i in e["inputs"]:
    p = os.path.join(run, i["path"])
    if i.get("sha256"):
        ok = os.path.exists(p) and sha(p) == i["sha256"]
    else:   # arms/mcp/raw/: the number of records the document read (reader.mcp_tool_calls' *.json, *.json.gz)
        ok = os.path.isdir(p) and i.get("files") == sum(1 for f in os.listdir(p) if f.endswith((".json", ".json.gz")))
    bad += not ok
    if not ok:
        print("    evidence input MISMATCH: %s" % i["path"])
n = len(e["inputs"])
print("  evidence inputs: %d of %d match the run directory" % (n - bad, n))
rec = os.path.join(out, e["gate"]["record"])
rec_ok = sha(rec) == e["gate"]["record_sha256"]
print("  gate record sha256: %s" % ("ok" if rec_ok else "MISMATCH"))
g = e["gate"]
print("  deploy gate: %s (%s)%s" % (g.get("verdict"), os.path.basename((g.get("policy") or {}).get("path") or ""),
      "".join("; %s: %s" % (os.path.basename(a.get("path") or ""), a.get("verdict")) for a in g.get("also_policies") or [])))
t = e.get("by_design_table_g") or {}
pins = []
for p in policies:
    pin = ((json.load(open(p)).get("graphlets") or {}).get("by_design_g_table") or {}).get("sha256")
    pins.append(pin == t.get("sha256"))
    if pin != t.get("sha256"):
        print("    %s pins area G's table %s, the run's is %s" % (os.path.basename(p), (pin or "")[:12],
                                                                   (t.get("sha256") or "")[:12]))
print("  area G by-design table: %s, %s cases, sha256 %s (%s)" % (
    t.get("path"), t.get("cases"), (t.get("sha256") or "")[:12],
    "every policy pins it" if all(pins) else "NOT pinned by every policy"))
verdicts = [g.get("verdict") == "GATE PASS"] + [a.get("verdict") == "PASS" for a in g.get("also_policies") or []]
print("  chapter G inputs: %d digests (%s)" % (len(e.get("graphlet_inputs") or {}), ", ".join(sorted(e.get("graphlet_inputs") or {}))))
m = json.load(open(man))
u = json.load(open(units))["_manifest"]
uk = u.get("content_sha256") == m.get("content_sha256") and u.get("scenarios") == len(m.get("scenarios") or [])
print("  units _manifest: %s, %s scenarios (%s)" % (u.get("content_sha256", "")[:12], u.get("scenarios"),
                                                    "the manifest's" if uk else "NOT the manifest's"))
if not all(verdicts):
    print("    a gate verdict is not PASS: the document states it (I.1, I.4, appendix F); look before sharing")
sys.exit(0 if not bad and rec_ok and uk and all(pins) and all(verdicts) and e.get("graphlet_inputs") else 1)
EOF
echo "bundle   $OUT: $(cd "$OUT" && echo $FILES) SHA256SUMS"
if [ $PROBLEMS = 0 ]; then echo "CHECKS CLEAN"; exit 0; fi
echo "CHECKS: $PROBLEMS to look at"; exit 1
