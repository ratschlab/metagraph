import sys
import os
import scipy as sp
import h5py
import gzip

basedir = '/cluster/work/grlab/projects/metagenome/data/gtex/'
result_file = os.path.join(basedir, 'queries', 'gencode.v30.all_junctions.result.txt.gz')
hdf5_out_file = os.path.join(basedir, 'queries', 'gencode.v30.all_junctions.result.hdf5')

matching = sp.zeros((350000, 10000), dtype='int16')
expression = sp.zeros((350000, 10000), dtype='int')

strain2idx = dict()
trans2idx = dict()
transcripts = []
strains = []
for i,_line in enumerate(gzip.open(result_file, 'r')):
    line = _line.decode('utf-8')
    if i > 0 and i % 10 == 0: 
        sys.stderr.write('.')
        if i % 100 == 0:
            sys.stderr.write('%i (data for %i exons in %i strains)\n' % (i, len(trans2idx), len(strain2idx)))
        sys.stderr.flush()
    if line[0] in ['P', '#', '\n']:
        continue
    key = line.strip().split('\t')
    ### add index of current transcripts
    if not key[1] in trans2idx:
        trans2idx[key[1]] = len(trans2idx)
        transcripts.append(key[1])
    ### add index of strains
    for match in key[2:]:
        try:
            ssl = match.split(':')
            strain = ssl[0].split('.')[0]
            exp = float('.'.join(ssl[0].split('.')[5:7]))
            cnt = int(ssl[1])
            if not strain in strain2idx:
               strain2idx[strain] = len(strain2idx)
               strains.append(strain)
            matching[trans2idx[key[1]], strain2idx[strain]] += cnt
            expression[trans2idx[key[1]], strain2idx[strain]] += int(exp * cnt)
        except:
            pass

### remove unused rows / cols
matching = matching[:len(transcripts)][:, :len(strains)]
expression = expression[:len(transcripts)][:, :len(strains)]

OUT = h5py.File(hdf5_out_file, 'w')
OUT.create_dataset(name='matching', data=matching, compression='gzip')
OUT.create_dataset(name='expression', data=expression, compression='gzip')
OUT.create_dataset(name='exons', data=sp.array(transcripts).view(sp.chararray).encode('utf-8'), compression='gzip')
OUT.create_dataset(name='strains', data=sp.array(strains).view(sp.chararray).encode('utf-8'), compression='gzip')
OUT.close()

