from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
with open('../../samples/NC2033.txt') as fh:
    with open('NC2033.fasta','w') as f_out:
        rawseq = fh.read().replace('\n', '')
        record = (SeqRecord(Seq(rawseq), 'NC2033.txt', '', ''),)
        SeqIO.write(record, f_out, 'fasta')
