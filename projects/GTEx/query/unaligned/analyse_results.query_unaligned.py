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
plotdir = os.path.join(basedir, 'plots', 'queries_unaligned')
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
 
if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <result.txt.gz>\n' % sys.argv[0])
    sys.exit(1)
result_file = sys.argv[1]
#result_file = os.path.join(basedir, 'queries', 'query_unlabeled.SRR2135314.result.txt')
sample_id = result_file.split('/')[-1].split('.')[-3]

matching = sp.zeros((210000, 10000), dtype='bool')

reads_per_strain = dict()
strains_per_read = dict()
unmapped = 0
total_reads = 0
for i,line in enumerate(open(result_file, 'r')):
    if i > 0 and i % 1000 == 0: 
        sys.stderr.write('.')
        if i % 10000 == 0:
            sys.stderr.write('%i (data for %i reads in %i strains)\n' % (i, len(read2idx), len(strain2idx)))
        sys.stderr.flush()
    ### skip empty and text lines
    if line[0] in ['P', '#', '\n']:
        continue
    key = line.strip().split('\t')
    ### add index of current read
    total_reads += 1
    if len(key) < 3:
        unmapped += 1
        continue
    ### add index of strains
    try:
        matches = [_.split('/')[-1].split('.')[0] for _ in key[2].split(':')]
        for match in matches:
            try:
                reads_per_strain[match] += 1
            except KeyError:
                reads_per_strain[match] = 1
        cnt = len(matches)
        try:
            strains_per_read[cnt] += 1
        except KeyError:
            strains_per_read[cnt] = 1
    except:
        pass

### remove samples that have too complete columns
#kidx = sp.where(sp.sum(matching > 0, axis=0) < 60000)[0]
#matching = matching[:, kidx]
#strains = strains[kidx]
#nidx = sp.where(sp.sum(matching, axis=1) == 0)[0]
#kidx = sp.where(sp.sum(matching, axis=1) > 0)[0]
#unmapped = sp.r_[unmapped, reads[nidx]]
#matching = matching[kidx, :]
#reads = reads[kidx]

### generate GTEx tissue dictionary from metadata
mfile = '/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt'
metadata = sp.loadtxt(mfile, dtype='str', delimiter='\t')
metadata = metadata[1:, :]
metadata = metadata[:, [16, 23, 32]]
dict_tissue = dict([(x[0], x[1]) for x in metadata])
dict_donor = dict([(x[0], x[2]) for x in metadata])

donors = sp.array([dict_donor[_] for _ in reads_per_strain])
donors_u = sp.unique(donors)
reads_per_donor = dict()
print('\nsummarizing per donor')
for i, d in enumerate(donors_u):
    if i > 0 and i % 10 == 0:
        sys.stderr.write('.')
        if i % 100 == 0:
            sys.stderr.write('%i\n' % i)
    d_idx = sp.where(donors == d)[0]
    reads_per_donor[d] = 0
    for dd in d_idx:
        reads_per_donor[d] += reads_per_strain[donors[dd]]

gs = gridspec.GridSpec(1, 2)

### plot distribution of matching donors
fig = plt.figure(figsize=(10, 10), dpi=200)
fig.suptitle(sample_id, fontsize=14), 
ax = fig.add_subplot(gs[0, 0])
brange = sp.arange(20)
#brange[-1] = sp.sum(matching_per_donor > 0, axis=1).max()
match_hist, _ = sp.histogram(reads_per_donor.values(), brange)
ax.bar(range(len(match_hist)), match_hist, width=1.0)
ax.bar(0.5, unmapped, width=1.0, color='grey')
match_hist[0] = unmapped
ax.fill_between(sp.arange(match_hist.shape[0]) + 0.5, 0, sp.cumsum(match_hist), color='grey', alpha=0.2)
xticks = sp.arange(0, 20, 2)
ax.set_xticks(xticks + 0.5)
xticks = xticks.astype('str')
xticks[-1] = '>' + xticks[-1]
ax.set_xticklabels(xticks)
ax.set_xlabel('donors')
ax.set_title('per donor')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
#plt.savefig(os.path.join(plotdir, sample_id + '.unmapped_donor_dist.pdf'), format='pdf', bbox_inches='tight')
#plt.close(fig)

### plot distribution of matching samples
#fig = plt.figure(figsize=(10, 10), dpi=200)
ax = fig.add_subplot(gs[0, 1])
brange = sp.arange(20)
brange[-1] = sp.sum(matching > 0, axis=1).max()
match_hist, _ = sp.histogram(reads_per_donor.values(), brange)
ax.bar(range(len(match_hist)), match_hist, width=1.0, color='gold')
ax.bar(0.5, unmapped, width=1.0, color='grey')
match_hist[0] = unmapped
ax.fill_between(sp.arange(match_hist.shape[0]) + 0.5, 0, sp.cumsum(match_hist), color='grey', alpha=0.2)
xticks = sp.arange(0, 20, 2)
ax.set_xticks(xticks + 0.5)
xticks = xticks.astype('str')
xticks[-1] = '>' + xticks[-1]
ax.set_xticklabels(xticks)
ax.set_xlabel('samples')
ax.set_title('per sample')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
plt.savefig(os.path.join(plotdir, sample_id + '.unmapped_sample_dist.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(plotdir, sample_id + '.unmapped_sample_dist.png'), format='png', bbox_inches='tight')
plt.close(fig)

