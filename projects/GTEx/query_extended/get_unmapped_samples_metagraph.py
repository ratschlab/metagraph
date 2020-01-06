import sys
import os
import gzip

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <result.txt.gz>\n' % sys.argv[0])
    sys.exit(1)
infname = sys.argv[1]

edits = dict()
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
        print(sl[1].split(':')[0])
        continue
    
    ### sort counts in case of batched realign
    if os.path.basename(infname).startswith('realign'):
        sample_counts = sorted(sl[2:], key=lambda x: int(x.split(':')[1]), reverse=True)
    else:
        sample_counts = sl[2:]

    if sample_counts[0].split(':')[0] != sample or len(sample_counts) > 1:
        pass
    else:
        print(sl[1].split(':')[0])
        continue

