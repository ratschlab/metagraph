import sys
import scipy as sp
import os
import cPickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pylab import cm 

### generate GTEx tissue dictionary from metadata
mfile = '/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt'
metadata = sp.loadtxt(mfile, dtype='str', delimiter='\t')
metadata = metadata[1:, :]
metadata = metadata[:, [16, 23]]
metadict = dict([(x[0], x[1]) for x in metadata])

tissueID = dict([(x, i) for i, x in enumerate(sp.unique(metadict.values()))])
_, normalizers = sp.unique(metadict.values(), return_counts=True)

counts_all = []
dfile = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.tsv'
dpickle = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.pickle'

gene_names = []

if not os.path.exists(dpickle):
    for l, line in enumerate(open(dfile, 'r')):
        if l > 0 and l % 1000 == 0:
            sys.stdout.write('.')
            if l % 10000 == 0:
                sys.stdout.write('%i\n' % l)
            sys.stdout.flush()
        if not '\t' in line:
            continue
        sl = line.strip('\n').split('\t')
        if len(sl[-1]) == 0:
            continue
        hits = [_.split('/')[-1].split('.')[0] for _ in sl[-1].split(':')]
        hits = [metadict[_] for _ in hits]
        uhits, count = sp.unique(hits, return_counts=True)
        uhits = [tissueID[_] for _ in uhits]
        counts = sp.zeros((len(tissueID),), dtype='int')
        counts[uhits] = count
        gene = sl[1].split('|')[1]
        gene_names.append(gene)
        counts_all.append(counts) 
    cPickle.dump((counts_all, normalizers, gene_names), open(dpickle, 'w'), -1)
else:
    counts_all, normalizers, gene_names = cPickle.load(open(dpickle, 'r'))
counts_all = sp.array(counts_all, dtype='float') / normalizers
gene_names = sp.array(gene_names)

### filter counts
kidx = sp.sum(counts_all > 0.2, axis=1) < 10
counts_all = counts_all[kidx, :]
gene_names = gene_names[kidx]

kidx = counts_all.max(axis=1) > 0.8
counts_all = counts_all[kidx, :]
gene_names = gene_names[kidx]

ucounts, ccounts = sp.unique(gene_names, return_counts=True)
ucounts = ucounts[ccounts > 1]

kidx = sp.in1d(gene_names, ucounts)
counts_all = counts_all[kidx, :]
gene_names = gene_names[kidx]

sidx = sp.argsort(gene_names)
gene_names = gene_names[sidx]
counts_all = counts_all[sidx, :]

### find interesting genes
gene_dict = dict()
for i, gene in enumerate(gene_names):
        try:
            gene_dict[gene].append(i)
        except KeyError:
            gene_dict[gene] = [i]
interesting = []
for gene in gene_dict:
    counts = counts_all[gene_dict[gene], :]
    if max(counts.max(axis=0) - counts.min(axis=0)) > 0.95:
        interesting.append(gene)

kidx = sp.in1d(gene_names, interesting)
counts_all = counts_all[kidx, :]
gene_names = gene_names[kidx]

### plotting
_, icounts = sp.unique(gene_names, return_inverse=True)
which_gene = (icounts % 2 == 1).astype('float')
counts_all = sp.c_[which_gene, counts_all]
_, iidx = sp.unique(gene_names, return_index=True)

fig = plt.figure(figsize=(10, 20), dpi=200)
ax = fig.add_subplot(111)
ax.matshow(counts_all, aspect='auto', cmap=cm.Blues)
ax.set_xticks(sp.arange(counts_all.shape[1]))
ax.set_xticklabels(sp.unique(metadict.values()), rotation=90)
ax.set_yticks(iidx)
ax.set_yticklabels(ucounts)
plt.savefig('transcripts_gtex.heatmap.pdf', format='pdf', bbox_inches='tight')

