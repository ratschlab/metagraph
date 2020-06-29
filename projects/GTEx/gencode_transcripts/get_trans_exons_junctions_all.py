import sys
import scipy as sp

genome_file = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v30/GRCh38.p12.genome.fa'
annotation = '/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v30/gencode.v30.annotation.gtf' 

trans_exons_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/trans_exons/gencode.v30.trans_exons_all.fa'
trans_exons_len_out = '/cluster/work/grlab/projects/metagenome/data/gtex/queries/trans_exons/gencode.v30.trans_exons_all.len.tsv'

K=41

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
        

### extract relevant junctions
junctions = dict([])
donors = dict()
acceptors = dict()
for t in transcripts:
    ### sort exons by coordinate
    curr = sp.array(transcripts[t])
    s_idx = sp.argsort(curr[:, 1].astype('int'))
    strand = t.split(':')[-1]
    if strand == '-':
        s_idx = s_idx[::-1]

    ### get all junctions
    j = s_idx[0]
    for k1 in range(len(s_idx) - 1):
        j = s_idx[k1]
        for k2 in range(k1 + 1, len(s_idx)):
            i = s_idx[k2]   
            chrm = transcripts[t][i][0]
            s0 = int(transcripts[t][j][1])
            e0 = int(transcripts[t][j][2])
            s1 = int(transcripts[t][i][1])
            e1 = int(transcripts[t][i][2])

            ### check length constraints for flanking exons
            ### see visual in https://genome.cshlp.org/content/24/1/25/F1.large.jpg
            if (e0 - s0 + 1) < K - 1 or (e1 - s1 + 1) < K - 1:
                continue
            if strand == '-':
                #key = ':'.join([chrm, strand, str(e1 - K + 1), str(e1), str(s0 - 1), str(s0 - 1 + K - 1)])
                key = ':'.join([chrm, strand, str(s1 - 1), str(s1 - 1 + K - 1), str(e0 - K + 1), str(e0)])
            else:
                #key = ':'.join([chrm, strand, str(e0 - K + 1), str(e0), str(s1 - 1), str(s1 - 1 + K - 1)])
                key = ':'.join([chrm, strand, str(s0 - 1), str(s0 - 1 + K - 1), str(e1 - K + 1), str(e1)])
            try:
                junctions[key].append(t)
            except KeyError:
                junctions[key] = [t]

### write out sequences for junctions
out = open(trans_exons_out, 'w')
out_len = open(trans_exons_len_out, 'w')
for junction in junctions:
    chrm, strand, s0, e0, s1, e1 = junction.split(':')
    #seq = genome[chrm][int(s0):int(e0)] + genome[chrm][int(s1):int(e1) + K - 1]
    seq = rev_comp(genome[chrm][int(s0):int(e0)]) + rev_comp(genome[chrm][int(s1):int(e1)])
    if strand == '-':
        seq = rev_comp(seq)
    write_fasta(out, junction + '\t' + ';'.join(junctions[junction]), seq)
    out_len.write('\t'.join([junction, str(int(s1) - int(e0)), ';'.join(junctions[junction])]) + '\n')
out.close()
out_len.close()

