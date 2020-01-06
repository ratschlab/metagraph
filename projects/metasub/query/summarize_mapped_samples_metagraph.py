import sys
import os
import gzip

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <K> <result.txt.gz>\n' % sys.argv[0])
    sys.exit(1)
K = int(sys.argv[1])
infname = sys.argv[2]

edits = dict()
unaligned = 0
READLEN = 76
sample = '<' + os.path.basename(infname).split('.')[1] + '>'
for l,line in enumerate(gzip.open(infname, 'r')):
    if l > 0 and l % 1000 == 0:
        sys.stderr.write('.')
        if l % 10000 == 0:
            sys.stderr.write('%i\n' % l)
        sys.stderr.flush()
    line = line.decode('utf-8').strip()
    sl = line.strip().split('\t')
    ### skip logging etx
    if len(sl) < 2 or not sl[1].startswith('SRR'):
        continue
    ### unmapped
    if len(sl) == 2:
        unaligned += 1
        continue
    ### construct map dict (expensive)
    #mapped = dict([_.split(':') for _ in sl[2:]])
    #if sample in mapped:
    #    del mapped[sample]
    #matches = max([int(_) for _ in mapped.values()])
    if sl[2].split(':')[0] != sample:
        matches = int(sl[2].split(':')[1])
    elif len(sl) > 3:
        matches = int(sl[3].split(':')[1])
    else:
        unaligned += 1
        continue
    e = READLEN - K - matches + 1
    try:
        edits[e] += 1
    except KeyError:
        edits[e] = 1
edits[READLEN] = unaligned

for e in sorted(edits):
    print('%i\t%i' % (e, edits[e]))
    
