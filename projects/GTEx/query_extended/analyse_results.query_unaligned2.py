import sys
import os
import scipy as sp
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/tools/python')
import viz.axes as axs
import gzip

basedir = '/cluster/work/grlab/projects/metagenome/data/gtex/'
plotdir = os.path.join(basedir, 'plots', 'queries_extended_unmapped')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
 
if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <result.txt.gz>\n' % sys.argv[0])
    sys.exit(1)
result_file = sys.argv[1]
#result_file = os.path.join(basedir, 'queries', 'query_unlabeled.SRR2135314.result.txt')
sample_id = result_file.split('/')[-1].split('.')[1]
strategy = result_file.split('/')[-1].split('.')[0]

matching = sp.zeros((20000, 10000), dtype='bool')

reads_per_strain = dict()
strains_per_read = dict()
strain2idx = dict()
strains = []
unmapped = 0
total_reads = 0
for i,line in enumerate(gzip.open(result_file, 'rb')):
    line = line.decode('utf-8')
    if i > 0 and i % 1000 == 0: 
        sys.stderr.write('.')
        if i % 10000 == 0:
            sys.stderr.write('%i (data for %i reads in %i strains)\n' % (i, len(read2idx), len(strain2idx)))
        sys.stderr.flush()
    ### skip empty and text lines
    if line[0] in ['#', '[', '\n']:
        continue
    key = line.strip().split('\t')
    ### add index of current read
    total_reads += 1
    if len(key) < 3:
        unmapped += 1
        continue
    ### add index of strains
    matches = [_.split(':')[0][1:-1] for _ in key[2:]]
    quals = [_.split(':')[1] for _ in key[2:]]
    handled = False
    for m,match in enumerate(matches):
        ### ignore if we are only aligning against ourselves
        if match == sample_id:
            continue
        if not match.startswith('SRR'):
            match = 'ref'
        ### ignore too low quality matches 
        if int(quals[m]) < 20:
            continue
        handled = True
        if not match in strain2idx:
            strain2idx[match] = len(strain2idx)
            strains.append(match)
        matching[total_reads, strain2idx[match]] = 1
        try:
            reads_per_strain[match] += 1
        except KeyError:
            reads_per_strain[match] = 1
    if not handled:
        unmapped += 1
    cnt = sum([int(_.startswith('SRR')) for _ in matches])
    try:
        strains_per_read[cnt] += 1
    except KeyError:
        strains_per_read[cnt] = 1

kidx = sp.where(sp.sum(matching, axis=1) > 0)[0]
matching = matching[kidx, :]

### generate GTEx tissue dictionary from metadata
mfile = '/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt'
metadata = sp.loadtxt(mfile, dtype='str', delimiter='\t')
metadata = metadata[1:, :]
metadata = metadata[:, [16, 23, 32]]
dict_tissue = dict([(x[0], x[1]) for x in metadata])
dict_tissue['ref'] = 'reference'
dict_donor = dict([(x[0], x[2]) for x in metadata])
dict_donor['ref'] = 'reference'

donors = sp.array([dict_donor[_] for _ in strains])
donors_u = sp.unique(donors)
matching_per_donor = sp.zeros((matching.shape[0], donors_u.shape[0]), dtype='bool')
print('\nsummarizing per donor')
for i, d in enumerate(donors_u):
    if i > 0 and i % 10 == 0:
        sys.stderr.write('.')
        if i % 100 == 0:
            sys.stderr.write('%i\n' % i)
    d_idx = sp.where(donors == d)[0]
    matching_per_donor[:, i] = matching[:, d_idx].max(axis=1)

tissues = sp.array([dict_tissue[_] for _ in strains])
tissues_u = sp.unique(tissues)
matching_per_tissue = sp.zeros((matching.shape[0], tissues_u.shape[0]), dtype='bool')
print('\nsummarizing per tissue')
for i, d in enumerate(tissues_u):
    d_idx = sp.where(tissues == d)[0]
    matching_per_tissue[:, i] = matching[:, d_idx].max(axis=1)

gs = gridspec.GridSpec(3, 1)

### plot distribution of matching donors
fig = plt.figure(figsize=(8, 8), dpi=200)
fig.suptitle(sample_id, fontsize=14), 
ax = fig.add_subplot(gs[0, 0])
bmax = sp.sum(matching_per_donor > 0, axis=1).max()
brange = sp.arange(520, bmax + 2, 2)
brange[0] = 1
#brange[-1] = sp.sum(matching_per_donor > 0, axis=1).max()
match_hist, _ = sp.histogram(matching_per_donor.sum(axis=1), brange)
ax.bar(sp.arange(len(match_hist)) + 1.5, match_hist, width=1.0)
ax.bar(0.5, unmapped, width=1.0, color='grey')
#match_hist[0] = unmapped
match_hist = sp.r_[[unmapped], match_hist]
ax.fill_between(sp.arange(match_hist.shape[0]) + 0.5, 0, sp.cumsum(match_hist), color='grey', alpha=0.2)
xticks = sp.arange(0, match_hist.shape[0])
ax.set_xticks(xticks + 0.5)
xticklabels = sp.r_[['0'], brange.astype('str')]
xticklabels[1] = '<' + xticklabels[2]
ax.set_xticklabels(xticklabels, rotation=90)
ax.set_xlabel('donors')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

### plot distribution of matching tissues
ax = fig.add_subplot(gs[1, 0])
bmax = sp.sum(matching_per_tissue > 0, axis=1).max()
brange = sp.arange(bmax+1)
brange[0] = 1
#brange[-1] = sp.sum(matching_per_tissue > 0, axis=1).max()
match_hist, _ = sp.histogram(matching_per_tissue.sum(axis=1), brange)
ax.bar(sp.arange(len(match_hist)) + 1.5, match_hist, width=1.0, color='green')
ax.bar(0.5, unmapped, width=1.0, color='grey')
#match_hist[0] = unmapped
match_hist = sp.r_[[unmapped], match_hist]
ax.fill_between(sp.arange(match_hist.shape[0]) + 0.5, 0, sp.cumsum(match_hist), color='grey', alpha=0.2)
xticks = sp.arange(0, match_hist.shape[0])
ax.set_xticks(xticks + 0.5)
xticklabels = sp.r_[['0'], brange.astype('str')]
ax.set_xticklabels(xticklabels, rotation=90)
ax.set_xlabel('tissues')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

### plot distribution of matching samples
ax = fig.add_subplot(gs[2, 0])
bmax = sp.sum(matching > 0, axis=1).max()
brange = sp.arange(8500, bmax, 50)
brange[0] = 1
match_hist, _ = sp.histogram(matching.sum(axis=1), brange)
ax.bar(sp.arange(len(match_hist)) + 1.5, match_hist, width=1.0, color='gold')
ax.bar(0.5, unmapped, width=1.0, color='grey')
match_hist = sp.r_[[unmapped], match_hist]
ax.fill_between(sp.arange(match_hist.shape[0]) + 0.5, 0, sp.cumsum(match_hist), color='grey', alpha=0.2)
xticks = sp.arange(0, match_hist.shape[0])
ax.set_xticks(xticks + 0.5)
xticklabels = sp.r_[['0'], brange.astype('str')]
xticklabels[1] = '<' + xticklabels[2]
ax.set_xticklabels(xticklabels, rotation=90)
ax.set_xlabel('samples')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

plt.tight_layout()
plt.savefig(os.path.join(plotdir, strategy + '.' + sample_id + '.unaligned_sample_dist.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(plotdir, strategy + '.' + sample_id + '.unaligned_sample_dist.png'), format='png', bbox_inches='tight')
plt.close(fig)

