from Bio import SeqIO

INPUT_FILE = '../../samples/fasta22.fas'
OUTPUT_FILE = 'fasta22_out.fas'

def retseq(seq_fh):
    """
    Parse a fasta file and store non empty records
    into the fullseqs list.
    :seq_fh: File handle of the input sequence
    :return: A list with non empty sequences
    """
    fullseqs = []
    for record in SeqIO.parse(seq_fh,'fasta'):
        if len(record.seq):
            fullseqs.append(record)
    return fullseqs

with open(INPUT_FILE) as in_fh:
    with open(OUTPUT_FILE, 'w') as out_fh:
        SeqIO.write(retseq(in_fh), out_fh, 'fasta')
