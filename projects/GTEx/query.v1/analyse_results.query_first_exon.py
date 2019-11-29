import sys
import scipy as sp
import os
import pickle
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pylab import cm 

### constants
K = 41

### generate GTEx tissue dictionary from metadata
mfile = '/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt'
metadata = sp.loadtxt(mfile, dtype='str', delimiter='\t')
metadata = metadata[1:, :]
metadata = metadata[:, [16, 23]]
metadict = dict([(x[0], x[1]) for x in metadata])

tissueID = dict([(x, i) for i, x in enumerate(sp.unique(metadict.values()))])
_, normalizers = sp.unique(metadict.values(), return_counts=True)

basedir = '/cluster/work/grlab/projects/metagenome/data/gtex'
plotdir = os.path.join(basedir, 'plots', 'queries_first_exon')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

### load aggregated results
IN = h5py.File(os.path.join(basedir, 'queries', 'gencode.v30.first_exons.result.hdf5'), 'r')
strains = IN['strains'][:].view(sp.chararray).decode('utf-8')
exons = IN['exons'][:].view(sp.chararray).decode('utf-8')

tissues = sp.array([metadict[_] for _ in strains])
tissues_u = sp.unique(tissues)

### load query lengths
tlen_file = os.path.join(basedir, 'queries', 'gencode.v30.first_exons.len.tsv')
exons2len = dict()
exons2genes = dict()
genes2exons = dict()
for line in open(tlen_file, 'r'):
    sl = line.strip().split('\t')
    exons2len[sl[0]] = int(sl[1])
    gene = sl[2].split(';')[0].rsplit(':', 2)[0]
    exons2genes[sl[0]] = gene
    try:
        genes2exons[gene].append(sl[0])
    except KeyError:
        genes2exons[gene] = [sl[0]]
exon_lens = sp.array([exons2len[_] for _ in exons])
genes = sp.array([exons2genes[_] for _ in exons])
gene_names = sp.array([_.split(':')[1] for _ in genes])

### get mean tissue expression and matching per tissue
expression = IN['expression'][:]

### remove exons that are
### 1: unexpressed
### 2: are in a gene with only one transcript annotated
is_expressed = expression > (exon_lens - K + 1)[:, sp.newaxis]
kidx1 = sp.sum(is_expressed, axis=1) > 0
kidx2 = sp.array([len(genes2exons[_]) > 1 for _ in genes]) 
kidx = kidx1 & kidx2
expression = expression[kidx, :]
is_expressed = is_expressed[kidx]
exons = exons[kidx]
exon_lens = exon_lens[kidx]
genes = genes[kidx]
gene_names = gene_names[kidx]
genes_u = sp.unique(genes)

mean_expression = sp.zeros((expression.shape[0], tissues_u.shape[0]), dtype='float')
is_expressed_fraction = sp.zeros((expression.shape[0], tissues_u.shape[0]), dtype='float')
for i, tissue in enumerate(tissues_u):
    t_idx = sp.where(tissues == tissue)[0]
    e_idx = is_expressed[:, t_idx]
    mean_expression[:, i] = sp.mean(expression[:, t_idx] * e_idx, axis=1) /  (exon_lens - K + 1)
    is_expressed_fraction[:, i] = sp.sum(e_idx, axis=1) / t_idx.shape[0]
del expression

### find best matching rows (more then 90%)
matching = IN['matching'][:]
matching = matching[kidx, :]
mean_matching = sp.zeros((matching.shape[0], tissues_u.shape[0]), dtype='float')
for i, tissue in enumerate(tissues_u):
    t_idx = sp.where(tissues == tissue)[0]
    mean_matching[:, i] = sp.mean(matching[:, t_idx], axis=1) / (exon_lens - K + 1)
best_matching = matching.max(axis=1) / (exon_lens - K + 1)
min_matching = mean_matching.min(axis=1)
max_matching = mean_matching.max(axis=1)
    
### get genes of interest
rmsd_max = []
coord_max = []
for gene in genes_u:
    gidx = sp.where(genes == gene)[0]
    curr = sp.around(mean_expression[gidx, :] * mean_matching[gidx, :])
    rmsd = [0]
    coords = [[0, 0]]
    for i in range(curr.shape[0]):
        for j in range(i + 1, curr.shape[0]):
            if min(curr[i, :].max(), curr[j, :].max()) < 4:
                continue
            z = curr[i, :] - curr[j, :]
            rmsd.append(sp.sqrt(sp.dot(z, z.T) / z.shape[0]))
            coords.append([i, j])
    r = sp.argmax(sp.array(rmsd))
    rmsd_max.append(rmsd[r])
    coord_max.append(coords[r])

TOP = 20
pmatrix = []
glabels = []
elabels = []
s_idx = sp.argsort(rmsd_max)[::-1]
for i in s_idx[:TOP]:
    gidx = sp.where(genes == genes_u[i])[0]
    cmax = coord_max[i]
    pmatrix.append(sp.around(mean_expression[gidx, :] * mean_matching[gidx, :])[cmax])
    glabels.append(genes[gidx][cmax])
    elabels.append(exons[gidx][cmax])
pmatrix = sp.vstack(pmatrix)
elabels = sp.hstack(elabels)
glabels = sp.hstack(glabels)

fig = plt.figure(figsize=(10, 20), dpi=200)
ax = fig.add_subplot(111)
ax.matshow(pmatrix.astype('float'), aspect='auto', cmap=cm.Blues)
ax.set_xticks(sp.arange(pmatrix.shape[1]))
ax.set_xticklabels(tissues_u, rotation=90)
ax.set_yticks(sp.arange(pmatrix.shape[0]))
ylabels = sp.array(['-'.join([glabels[_].split(':')[-1], elabels[_]]) for _ in range(elabels.shape[0])])
ax.set_yticklabels(ylabels)
plt.savefig(os.path.join(plotdir, 'diff_first_exons.heatmap.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(plotdir, 'diff_first_exons.heatmap.png'), format='png', bbox_inches='tight')

     
import pdb
pdb.set_trace()

#
#if not os.path.exists(dpickle):
#    for l, line in enumerate(open(dfile, 'r')):
#        if l > 0 and l % 1000 == 0:
#            sys.stdout.write('.')
#            if l % 10000 == 0:
#                sys.stdout.write('%i\n' % l)
#            sys.stdout.flush()
#        if not '\t' in line:
#            continue
#        sl = line.strip('\n').split('\t')
#        if len(sl[-1]) == 0:
#            continue
#        hits = [_.split('/')[-1].split('.')[0] for _ in sl[-1].split(':')]
#        hits = [metadict[_] for _ in hits]
#        uhits, count = sp.unique(hits, return_counts=True)
#        uhits = [tissueID[_] for _ in uhits]
#        counts = sp.zeros((len(tissueID),), dtype='int')
#        counts[uhits] = count
#        gene = sl[1].split('|')[1]
#        gene_names.append(gene)
#        counts_all.append(counts) 
#    cPickle.dump((counts_all, normalizers, gene_names), open(dpickle, 'w'), -1)
#else:
#    counts_all, normalizers, gene_names = cPickle.load(open(dpickle, 'r'))
#counts_all = sp.array(counts_all, dtype='float') / normalizers
#gene_names = sp.array(gene_names)
#
#### filter counts
#kidx = sp.sum(counts_all > 0.2, axis=1) < 10
#counts_all = counts_all[kidx, :]
#gene_names = gene_names[kidx]
#
#kidx = counts_all.max(axis=1) > 0.8
#counts_all = counts_all[kidx, :]
#gene_names = gene_names[kidx]
#
#ucounts, ccounts = sp.unique(gene_names, return_counts=True)
#ucounts = ucounts[ccounts > 1]
#
#kidx = sp.in1d(gene_names, ucounts)
#counts_all = counts_all[kidx, :]
#gene_names = gene_names[kidx]
#
#sidx = sp.argsort(gene_names)
#gene_names = gene_names[sidx]
#counts_all = counts_all[sidx, :]
#
#### find interesting genes
#gene_dict = dict()
#for i, gene in enumerate(gene_names):
#        try:
#            gene_dict[gene].append(i)
#        except KeyError:
#            gene_dict[gene] = [i]
#interesting = []
#for gene in gene_dict:
#    counts = counts_all[gene_dict[gene], :]
#    if max(counts.max(axis=0) - counts.min(axis=0)) > 0.95:
#        interesting.append(gene)
#
#kidx = sp.in1d(gene_names, interesting)
#counts_all = counts_all[kidx, :]
#gene_names = gene_names[kidx]
#
#### plotting
#_, icounts = sp.unique(gene_names, return_inverse=True)
#which_gene = (icounts % 2 == 1).astype('float')
#counts_all = sp.c_[which_gene, counts_all]
#_, iidx = sp.unique(gene_names, return_index=True)
#
#fig = plt.figure(figsize=(10, 20), dpi=200)
#ax = fig.add_subplot(111)
#ax.matshow(counts_all, aspect='auto', cmap=cm.Blues)
#ax.set_xticks(sp.arange(counts_all.shape[1]))
#ax.set_xticklabels(sp.unique(metadict.values()), rotation=90)
#ax.set_yticks(iidx)
#ax.set_yticklabels(ucounts)
#plt.savefig('exons_gtex.heatmap.pdf', format='pdf', bbox_inches='tight')
#
