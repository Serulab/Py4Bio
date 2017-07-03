
from Bio import SeqIO

out_dir = '../samples/bioinfo/seqs/'

with open('../samples/vectorssmall.fasta') as fh:
    generator = SeqIO.parse(fh, 'fasta')
    for seq in generator:
        name = seq.name.split('|')[1]
        with open(out_dir + name + '.fasta', 'w') as fw:
            SeqIO.write(seq, fw, 'fasta')
