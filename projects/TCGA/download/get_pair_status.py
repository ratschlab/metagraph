import sys
import os

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <header.txt>\n' % sys.argv[0])
    sys.exit(1)
fname = sys.argv[1]

for line in open(fname, 'r'):
    if not line.startswith('@PG'):
        continue
    sl = line.strip().split('\t')
    ssl = list(filter(lambda x: len(x) > 0, sl[4].split(' ')))
    idx = ssl.index('--readFilesIn')
    if ssl[idx + 2].startswith('/'):
        print('paired')
    else:
        print('single')


