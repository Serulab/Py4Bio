from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline

INPUT_SEQUENCE = open('../../samples/hsc1.fasta')
OUTPUT_SEQUENCE = 'primer.txt'
sfile = open('../../samples/hsc1.fasta')
myseq = SeqIO.read(sfile, 'fasta')
title = myseq.id
seq = str(myseq.seq).upper()
win_size = 45
i = 0
number_l = []
while i <= (len(seq) - win_size):
    number_l.append(seq[i:i + win_size].count('AAT'))
    i += 1 # This is the same as i = i+1
pos = number_l.index(max(number_l))
pr_cl = Primer3Commandline(sequence=INPUT_SEQUENCE, auto=True)
pr_cl.outfile = OUTPUT_SEQUENCE
pr_cl.osize = 18
pr_cl.maxsize = 20
pr_cl.minsize = 15
pr_cl.explainflag = 1
pr_cl.target = (pos, win_size)
pr_cl.prange = (win_size, len(seq))
primer_cl()
