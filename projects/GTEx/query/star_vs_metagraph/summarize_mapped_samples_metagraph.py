import sys
import gzip

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <result.txt.gz>\n' % sys.argv[0])
    sys.exit(1)
infname = sys.argv[1]

edits = dict()
unaligned = 0
READLEN = 76
for l,line in enumerate(gzip.open(infname, 'r')):
    if l > 0 and l % 100000 == 0:
        sys.stderr.write('.')
        if l % 1000000 == 0:
            sys.stderr.write('%i\n' % l)
        sys.stderr.flush()
    line = line.decode('utf-8').strip()
    if int(line) == 0:
        unaligned += 1
        continue
    e = READLEN - 40 - int(line)
    try:
        edits[e] += 1
    except KeyError:
        edits[e] = 1
edits[READLEN] = unaligned

for e in sorted(edits):
    print('%i\t%i' % (e, edits[e]))
    
