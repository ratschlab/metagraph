import sys
import scipy as sp

genome_file = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/GRCh38.p12.genome.fa'
annotation = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v29_transcriptome/gencode.v30.annotation.gtf' 

first_exons_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.first_exons.fa'
first_exons_len_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.first_exons.len.tsv'
all_exons_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.all_exons.fa'
all_exons_len_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/gencode.v30.all_exons.len.tsv'

def read_tags(s):
    tags = []
    for _ in s.strip(';').split(';'):
        tag = _.strip(' ').split(' ')
        tags.append((tag[0], tag[1].strip('"')))
    return dict(tags)


REV_DICT = {'A':'T',
            'T':'A',
            'C':'G',
            'G':'C',
            'N':'N'}
def rev_comp(s):
    return ''.join([REV_DICT[_] for _ in s][::-1])
    

def write_fasta(out, seqid, seq, w=80):
    out.write('>' + seqid + '\n')
    for i in range(0, len(seq), 80):
        out.write(seq[i:min(i+80, len(seq))] + '\n')

### extract transcripts
transcripts = dict()
for line in open(annotation, 'r'):
    if line[0] == '#':
        continue
    sl = line.strip().split('\t')
    if sl[2] != 'exon':
        continue
    tags = read_tags(sl[8])
    key = ':'.join([tags['gene_id'], tags['gene_name'], tags['transcript_id'], sl[6]])
    value = [sl[0], sl[3], sl[4]]
    try:
        transcripts[key].append(value)
    except KeyError:
        transcripts[key] = [value]

### parse genome
seq = []
curr_id = ''
genome = dict()
for line in open(genome_file, 'r'):
    if line[0] == '>':
        if len(seq) > 0:
            genome[curr_id] = ''.join(seq)
        curr_id = line.strip().split()[0][1:]
        sys.stderr.write('Parsing %s\n' % curr_id)
        seq = []
        continue
    seq.append(line.strip())
if len(seq) > 0:
    genome[curr_id] = ''.join(seq)
seq = []
        

### extract exons and first exons
first_exons = dict([])
exons = dict([])
donors = dict()
acceptors = dict()
for t in transcripts:
    ### sort exons by coordinate
    curr = sp.array(transcripts[t])
    s_idx = sp.argsort(curr[:, 1].astype('int'))
    strand = t.split(':')[-1]
    if strand == '-':
        s_idx = s_idx[::-1]
    ### get first exons
    key = ':'.join(transcripts[t][s_idx[0]]) + ':' + strand
    try:
        first_exons[key].append(t)
    except KeyError:
        first_exons[key] = [t]
    try:
        exons[key].append(t)
    except KeyError:
        exons[key] = [t]

    ### get all exons
    j = s_idx[0]
    for i in s_idx[1:]:
        chrm = transcripts[t][i][0]
        s0 = int(transcripts[t][j][1])
        e0 = int(transcripts[t][j][2])
        s1 = int(transcripts[t][i][1])
        e1 = int(transcripts[t][i][2])
        if strand == '-':
            donor = rev_comp(genome[chrm][s1-3:s1-1])
            acceptor = rev_comp(genome[chrm][e0:e0+2])
        else:
            donor = genome[chrm][e0:e0+2]
            acceptor = genome[chrm][s1-3:s1-1]
        try:
            donors[donor] += 1
        except KeyError:
            donors[donor] = 1
        try:
            acceptors[acceptor] += 1
        except KeyError:
            acceptors[acceptor] = 1
        try:
            exons[key].append(t)
        except KeyError:
            exons[key] = [t]
        j = i

### write out sequences for first exons
out = open(first_exons_out, 'w')
out_len = open(first_exons_len_out, 'w')
for exon in first_exons:
    strand = exon.split(':')[3]
    chrm = exon.split(':')[0]
    start = int(exon.split(':')[1]) - 1
    end = int(exon.split(':')[2])
    seq = genome[chrm][start:end]
    if strand == '-':
        seq = rev_comp(seq)
    write_fasta(out, exon + '\t' + ';'.join(first_exons[exon]), seq)
    out_len.write('\t'.join([exon, str(len(seq)), ';'.join(first_exons[exon])]) + '\n')
out.close()
out_len.close()

### write out sequences for all exons
out = open(all_exons_out, 'w')
out_len = open(all_exons_len_out, 'w')
for exon in exons:
    strand = exon.split(':')[3]
    chrm = exon.split(':')[0]
    start = int(exon.split(':')[1]) - 1
    end = int(exon.split(':')[2])
    seq = genome[chrm][start:end]
    if strand == '-':
        seq = rev_comp(seq)
    write_fasta(out, exon, seq)
    out_len.write('\t'.join([exon, str(len(seq)), ';'.join(first_exons[exon])]) + '\n')
out.close()
out_len.close()

import pdb
pdb.set_trace()
