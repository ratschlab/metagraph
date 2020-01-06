import sys
import scipy as sp
import h5py

result_file = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.abundances.tsv'
hdf5_out_file = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.abundances.hdf5'

matching = sp.zeros((210000, 10000), dtype='int16')
expression = sp.zeros((210000, 10000), dtype='int')

strain2idx = dict()
trans2idx = dict()
transcripts = []
strains = []
for i,line in enumerate(open(result_file, 'r')):
    if i > 0 and i % 10 == 0: 
        sys.stderr.write('.')
        if i % 100 == 0:
            sys.stderr.write('%i (data for %i transcripts in %i strains)\n' % (i, len(trans2idx), len(strain2idx)))
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
            strain = ssl[0].split('.')[0][1:]
            exp = int(ssl[0][-2])
            cnt = int(ssl[1])
            if not strain in strain2idx:
               strain2idx[strain] = len(strain2idx)
               strains.append(strain)
            matching[trans2idx[key[1]], strain2idx[strain]] += cnt
            expression[trans2idx[key[1]], strain2idx[strain]] += (exp + 1) * cnt
        except:
            pass

### remove unused rows / cols
matching = matching[:len(transcripts)][:, :len(strains)]
expression = expression[:len(transcripts)][:, :len(strains)]

OUT = h5py.File(hdf5_out_file, 'w')
OUT.create_dataset(name='matching', data=matching, compression='gzip')
OUT.create_dataset(name='expression', data=expression, compression='gzip')
OUT.create_dataset(name='transcripts', data=sp.array(transcripts).view(sp.chararray).encode('utf-8'), compression='gzip')
OUT.create_dataset(name='strains', data=sp.array(strains).view(sp.chararray).encode('utf-8'), compression='gzip')
OUT.close()

