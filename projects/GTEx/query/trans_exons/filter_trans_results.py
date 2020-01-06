import sys
import os
import pickle
import re

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <trans_result.txt>\n' % sys.argv[0])
    sys.exit(1)
infname = sys.argv[1]

### load transcriptome alignment
tr_align = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all.mapped_transcriptome.sam'
tr_align_pickle = tr_align + '.pickle'
if not os.path.exists(tr_align_pickle):
    transcriptome_hits = set()
    for line in open(tr_align, 'r'):
        if line[0] == '@':
            continue
        sl = line.strip().split('\t')
        if len(sl) < 5:
            continue
        if sl[5] == '80M':
            transcriptome_hits.add(sl[0])
    pickle.dump(transcriptome_hits, open(tr_align_pickle, 'wb'))
else:
    transcriptome_hits = pickle.load(open(tr_align_pickle, 'rb'))

### load genome alignments
ge_align = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.trans_exons_all.mapped_genome.Aligned.out.sam'
ge_align_pickle = ge_align + '.pickle'
if not os.path.exists(ge_align_pickle):
    genome_hits = set()
    for line in open(ge_align, 'r'):
        if line[0] == '@':
            continue
        sl = line.strip().split('\t')
        if len(sl) < 5:
            continue
        matches = sum([int(_) for _ in re.split(r'M', re.sub(r'[0-9]*[DIHSN]', '', sl[5]))[:-1]])
        if matches >= 70:
            genome_hits.add(sl[0])
    pickle.dump(genome_hits, open(ge_align_pickle, 'wb'))
else:
    genome_hits = pickle.load(open(ge_align_pickle, 'rb'))


hits = dict()
for line in open(infname, 'r'):
    if line[0].startswith('['):
        continue
    sl = line.strip().split('\t')
    if sl[2] == '*':
        continue
    if sl[5] != '80':
        continue
    #print('> ' + sl[0])
    #print(sl[1])
    try:
        hits[sl[1]].append(sl[0])
    except KeyError:
        hits[sl[1]] = [sl[0]]

for hit in hits:
    broken = False
    for h in hits[hit]:
        if h in transcriptome_hits or h in genome_hits:
            broken = True
            break
    if broken:
        continue
    print('>' + ';'.join(hits[hit]))
    print(hit)
