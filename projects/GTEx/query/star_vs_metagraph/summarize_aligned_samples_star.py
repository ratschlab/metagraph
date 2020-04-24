import sys
import gzip
import re

edits = dict()
unaligned = 0
READLEN = 76
for l,line in enumerate(sys.stdin):
#for l,line in enumerate(open(sys.argv[1], 'r')):
    if l > 0 and l % 100000 == 0:
        sys.stderr.write('.')
        if l % 1000000 == 0:
            sys.stderr.write('%i\n' % l)
        sys.stderr.flush()
    sl = line.strip().split('\t')
    if sl[2] == '*':
        unaligned += 1
        continue
    e = READLEN - sum([int(_) for _ in re.split('M', re.sub(r'[0-9]+[SHPIDN]', '', sl[5]))[:-1]])
    #if e > 0:
    #    c = sum([int(_) for _ in re.split('[^0-9]', re.sub('[0-9]+=', '', sl[6]))[:-1]])
    try:
        edits[e] += 1
    except KeyError:
        edits[e] = 1
edits[READLEN] = unaligned

for e in sorted(edits):
    print('%i\t%i' % (e, edits[e]))
    
