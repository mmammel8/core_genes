import re
import sys

filename1 = "EC11.fasta"
fileOut = "EC11error.txt"

seq = dict()
targs = []
data_file = open(filename1, 'r')
data = ''.join(data_file.readlines())
data_file.close()
lines = data.split('\n')
ok = 0
for line in lines:
    if len(line) > 2:
        if line[0] == '>':
            targeta = line[1:].strip()
            targs.append(targeta)
        else:
            seq[targeta] = line.strip()

dim = len(targs)
dist = [[0 for _ in range(dim)] for _ in range(dim)]
for targeta in range(dim):
    fasta = seq[targs[targeta]]
    for targetb in range(dim):
        fastq = seq[targs[targetb]]
        ens = 0
        mm = 0
        for index in range(len(fasta)):
            if fasta[index] == 'N' or fasta[index] == '-':
                ens += 1
            elif fasta[index] != fastq[index] and fastq[index] != '-' and fastq[index] != 'N':
                #line = str(index) + ',' + snp_list[index] + ',' + fasta[index] + ',' + fastq[index]  + '\n'
                #out_file.write(line)
                mm += 1
        if targeta == targetb:
            dist[targeta][targetb] = ens
        else:
            dist[targeta][targetb] = mm

out_file = open(fileOut, 'w')
out = ",".join(targs) + '\n'    
out_file.write(out)
for line in dist:
    out = str(line) + '\n'
    out_file.write(out)
out_file.close()
