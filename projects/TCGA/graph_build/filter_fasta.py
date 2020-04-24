import sys
import gzip

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <infile> <K>\n' % sys.argv[0])
    sys.exit(1)
fname = sys.argv[1]
K = int(sys.argv[2])

cnt = 0
for line in gzip.open(fname, 'rb'):
    line = line.decode('utf-8')
    if line.startswith('>'): 
        cnt += 1
    elif len(line.strip()) >= 2 * K:
        print('>seq_%i' % cnt + '\n' + line.strip())
        
