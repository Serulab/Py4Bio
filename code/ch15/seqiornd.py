import random

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

TOTAL_SEQUENCES = 500
MIN_SIZE = 400
MAX_SIZE = 1500

def new_rnd_seq(seq_len):
    """
    Generate a random DNA sequence with a sequence length
    of "sl" (int).
    return: A string with a DNA sequence.
    """
    s = ''
    while len(s) < seq_len:
        s += random.choice('ATCG')
    return s

with open('randomseqs.txt','w') as new_fh:
    for i in range(1, TOTAL_SEQUENCES + 1):
        # Select a random number between MIN_SIZE and MAX_SIZE
        rsl = random.randint(MIN_SIZE, MAX_SIZE)
        # Generate the random sequence
        rawseq = new_rnd_seq(rsl)
        # Generate a correlative name
        seqname = 'Sequence_number_{0}'.format(i)
        rec = SeqRecord(Seq(rawseq), id=seqname, description='')
        SeqIO.write([rec], new_fh, 'fasta')
