import sys
import os
import glob
import pdb

basedir = '/cluster/work/grlab/projects/metagenome/raw_data/tcga/data'

flist = glob.glob(os.path.join(basedir, '*', '*', '*.txt'))
print('found %i files in %s' % (len(flist), basedir))

pair_status = dict()
for f,fname in enumerate(flist):
    if f > 0 and f % 100 == 0:
         sys.stderr.write('.')
         if f % 1000 == 0:
            sys.stderr.write('%i\n' % f)
         sys.stderr.flush()

    uid = fname.split('/')[-2]
    if uid in pair_status:
        continue
    for line in open(fname, 'r'):
        if not line.startswith('@PG'):
            continue
        sl = line.strip().split('\t')
        ssl = list(filter(lambda x: len(x) > 0, sl[4].split(' ')))
        idx = ssl.index('--readFilesIn')
        if ssl[idx + 2].startswith('/'):
            pair_status[uid] = 'paired'
        else:
            pair_status[uid] = 'single'

out = open('pair_status.tsv', 'w')
for uid in pair_status:
    out.write('%s\t%s\n' % (uid, pair_status[uid]))
out.close()
