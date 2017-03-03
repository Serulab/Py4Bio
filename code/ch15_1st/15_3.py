from Bio import SeqIO

INPUT_FILE = '../../samples/fasta22.fas'
OUTPUT_FILE = 'fasta22_out2.fas'

def retseq(seq_fh):
    """
    Parse a fasta file and returns non empty records
    :seq_fh: File handle of the input sequence
    :return: Non empty sequences
    """
    for record in SeqIO.parse(seq_fh, 'fasta'):
        if len(record.seq):
            yield record

with open(INPUT_FILE) as in_fh:
    with open(OUTPUT_FILE, 'w') as out_fh:
        SeqIO.write(retseq(in_fh), out_fh, 'fasta')
