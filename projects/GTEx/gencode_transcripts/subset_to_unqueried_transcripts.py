import sys
import pickle
import scipy as sp

transcript_fa = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/gencode.v29.transcripts.fa'
transcript_fa_to_complete = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/gencode.v29.transcripts_to_complete.fa'
result_completed = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.abundances.completed'
result_pickle = '/cluster/work/grlab/projects/metagenome/data/gtex/gencode.v29.result.abundances.pickle'

### collect completed queries
completed = set(sp.loadtxt(result_completed, dtype='str', delimiter='\t'))

def write_seq(curr_id, seq, out, w=80):
    out.write('>' + curr_id + '\n')
    for i in range(0, len(seq), w):
        out.write(seq[i:min(len(seq), i+w)] + '\n')

curr_id = ''
seq = []
out = open(transcript_fa_to_complete, 'w')
for line in open(transcript_fa, 'r'):
    if line[0] == '>':
        if len(seq) > 0 and not curr_id in completed:
            write_seq(curr_id, ''.join(seq), out)
        seq = []
        curr_id = line.strip()[1:]
        continue
    seq.append(line.strip())
if len(seq) > 0 and not curr_id in completed:
    write_seq(curr_id, ''.join(seq), out)
out.close()
