import os
import sys
import re
from subprocess import Popen, PIPE

#directory of alignment files to save
mypath = "W:/Mark_backup/ECSHSNPs2/aln18/"
#mypath = "W:/Mark_backup/ECSHSNPs2/alne1/"

#list of accession numbers to add to alignments
flist = "list_2.txt"

#list of genes to BLAST
#format is: locus	start	end	strand	length	locus
acc = "U00096_3"
fgenes = acc + ".txt"

#BLAST databases to use
mydb = 'WGS1'

#minimum match length
MINLEN = 0.9

#number of processors for BLAST
NUM_THREADS = 8

#must include fasta file of all accessions used in fgenes file
#saved in same directory, and named accession.fna

#globals
blin = "seq2.txt"
blout = "result2.txt"
orgs = []
minperc = []
tally1 = 0 #number of ok genes not saved because no mm
tally2 = 0 #number of ok genes with mismatch in alignment
base_val = {'0':0, 'A':1, 'C':2, 'G':4, 'T':8, 'N':15, 'X':0, 'R':5, 'Y':10, 'S':6, 'W':9, 'K':12, 'M':3, 'B':14, 'D':13, 'H':11, 'V':7, '-':15, '?':15}

class Alignment:
    """
    Class for holding alignment of a gene
    """
    def __init__(self, qlen, locus):
        """
        start with query gene
        """
        self.qlen = qlen
        self.locus = locus
        self.alignment = ""
        self.count = [0 for _ in range(len(orgs))]
        self.mm = 0
        self.sseq = ""
        self.qseq = ""
        self.org1 = None
        self.qst = 0
        self.qen = 0
        self.perc = 0
        self.usedset = set()

    def save_alignment(self):
        """
        append alignment of locus to text file in directory
        """
        fname = mypath + self.locus + ".txt"
        data_file = open(fname, 'a')
        data_file.write(self.alignment)
        data_file.close()   
  
    def process_match(self):
        """
        add to alignment if ok
        """ 
        idx = -1
        for ii in range(len(orgs)):
            if self.org1.startswith(orgs[ii]):
                idx = ii
        if idx > -1 and orgs[idx] not in self.usedset:
            if self.qen < self.qst:
                len1 = self.qst - self.qen + 1
            else:
                len1 = self.qen - self.qst + 1
            self.qseq = self.qseq.upper()
            self.sseq = self.sseq.upper()
            if self.perc >= minperc[idx] and len1 >= MINLEN * self.qlen:
                aseq = ""
                for pos in range(self.qst - 1):
                    aseq += '-'
                for pos in range(len(self.qseq)):
                    b1 = self.qseq[pos:pos + 1]
                    if b1 > '@':
                        b2 = self.sseq[pos:pos + 1]
                        aseq += b2
                        if base_val[b1] & base_val[b2] == 0:
                            self.mm += 1
                for pos in range(self.qen, self.qlen):
                    aseq += '-'
                self.count[idx] += 1
                if self.count[idx] == 1:
                    self.alignment += ">" + orgs[idx] + "\n" + aseq + "\n" 

    def capture(self):
        """
        Read BLAST result file and return aligning sequences
        """
        global tally1, tally2

        #first check for accessions already present
        fname = mypath + self.locus + ".txt"
        data_file = open(fname, 'r')
        data = ''.join(data_file.readlines())
        data_file.close()
        lines = data.split('\n')
        for line in lines:
            if len(line) > 1 and line[0] == '>':
                self.usedset.add(line[1:])

        returnval = 0
        bool = True
        data_file = open(blout, 'r')
        data = ''.join(data_file.readlines())
        data_file.close()
        lines = data.split('\n')
        # Sbjct  21809  GAGTTACGCAAACTCAACCAGTTTCGTACGTTTGCTCGAGGTTTTGATAATGTCCGCCAG  21750
        progq = re.compile("Query\s+(\d+)\s+(\S+)\s+(\d+)")
        progs = re.compile("Sbjct\s+\d+\s+(\S+)\s+\d+")
        # Identities = 677/687 (99%), Gaps = 0/687 (0%)
        progi = re.compile("\sIdentities\s=\s(\d+)/(\d+)")
        progh = re.compile(">(\S+)")

        for line in lines:
            if len(line) > 1:
    	        result = progq.match(line)
                if  result != None:
                    self.qseq += result.group(2)
                    if self.qst == 0: 
                        self.qst = int(result.group(1))
                    self.qen = int(result.group(3))
    	        result = progs.match(line)
                if  result != None:
                    self.sseq += result.group(1)
                result = progi.match(line)
                if result != None:
                    idN = float(result.group(1))
                    idD = float(result.group(2))
                    self.perc = 100.0 * idN / idD     
                result = progh.match(line)
                if result != None:
                    if self.org1 != None:
                        self.process_match()
                    self.org1 = result.group(1)
                    #clear
                    self.qst = 0
                    self.qen = 0
                    self.sseq = ""
                    self.qseq = "" 
                elif re.match(" Score =", line) != None:
                    if self.qst > 0 and self.org1 != None:
                        self.process_match()
                    #clear
                    self.qst = 0
                    self.qen = 0
                    self.sseq = ""
                    self.qseq = ""
                elif re.match("Lambda", line) != None and bool:
                    #end of section
                    bool = False
                    if self.org1 != None:
                        self.process_match()
                    #clear
                    self.qst = 0
                    self.qen = 0
                    self.sseq = ""
                    self.qseq = ""
        self.save_alignment()

def blast_it():
    """
    perform BLAST command
    """
    blast_args = ['blastn', '-task', 'blastn', '-db', mydb, '-query', blin, '-out', blout, '-num_descriptions', '1', '-num_alignments', '1600', '-num_threads', str(NUM_THREADS)]
    process = Popen(blast_args, stdout=PIPE)
    (stdout, stderr) = process.communicate()

def load(acc):
    """
    load genome sequence
    """
    seq = ""
    data_file = open(acc + ".fna", 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        if len(line) >= 1 and line[0] != '>':
            seq += line.strip()
    return seq

def rev_comp(seq):
    """
    return reverse complement
    """
    complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    return "".join([complement[base] for base in reversed(seq)])

def save_q(seq, strand, locus):
    """
    save seq of locus to file for BLAST
    rev compl if strand is minus
    """
    if strand == '-':
        seq = rev_comp(seq)
    data_file = open(blin, 'w')
    out = ">" + locus + "\n"
    data_file.write(out)
    out = seq + "\n"
    data_file.write(out)
    data_file.close()  

if __name__ == '__main__':
    #get list of strains to include
    data_file = open(flist, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    for line in lines:
        row = line.split('\t')
        if len(row) > 1:
            orgs.append(row[0])
            minperc.append(float(row[3]))

    #read gene file, BLAST each gene and get result
    prev_end = 0
    data_file = open(fgenes, 'r')
    data = ''.join(data_file.readlines())
    data_file.close()
    lines = data.split('\n')
    gseq = ''
    genome_seq = load(acc)
    gacc = acc
    lg = len(genome_seq)
    assert(lg > 1)
    for line in lines:
        row = line.split('\t')
        if len(row) > 4:
            locus, s, e, strand = row[:4]
            start = int(s)
            end = int(e)
            if start <= prev_end:
                start = prev_end + 1
            prev_end = end
            if start > 0 and start < end and end <= lg:
                fname = mypath + locus + ".txt"
                if os.path.isfile(fname): #only add to genes in directory
                    save_q(genome_seq[start - 1:end], strand, locus)
                    print locus
                    blast_it()
                    al1 = Alignment(end - start + 1, locus)
                    al1.capture()
