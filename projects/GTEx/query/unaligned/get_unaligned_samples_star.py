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
        if int(sl[1]) & 64 == 64:
            print(sl[0])
        else:
            print(sl[0][:-1] + '2')
        unaligned += 1
