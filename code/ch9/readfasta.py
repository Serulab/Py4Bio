from Bio import SeqIO

FILE_IN = '../../samples/3seqs.fas'

with open(FILE_IN) as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        id_ = record.id
        seq = record.seq
        print('Name: {0}, size: {1}'.format(id_, len(seq)))
