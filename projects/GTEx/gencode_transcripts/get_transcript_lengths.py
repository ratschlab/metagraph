import sys
import pickle
import scipy as sp

transcript_fa = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/gencode.v29.transcripts.fa'
transcript_lens = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/gencode.v29.transcripts.length.tsv'

curr_id = ''
seq = 0
out = open(transcript_lens, 'w')
for line in open(transcript_fa, 'r'):
    if line[0] == '>':
        if seq > 0:
            out.write('%s\t%i\n' % (curr_id, seq))
        seq = 0
        curr_id = line.strip()[1:]
        continue
    seq += len(line.strip())
if seq > 0:
    out.write('%s\t%i\n' % (curr_id, seq))
out.close()
