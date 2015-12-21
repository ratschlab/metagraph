import sys

upper = int(sys.argv[1])
kmer = int(sys.argv[2])
text = []
for i, line in enumerate(sys.stdin):
    if line[0] == '>':
        continue
    text.append(line.strip())
text = ''.join(text)

for i in range(upper):
    print text[i:i+kmer]
    
