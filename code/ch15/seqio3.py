from Bio import SeqIO
INPUT_FILE = '../../samples/fasta22.fas'
OUTPUT_FILE = 'fasta22_out3.fas'
with open(INPUT_FILE) as in_fh:
   with open(OUTPUT_FILE, 'w') as out_fh:
       for record in SeqIO.parse(in_fh, 'fasta'):
           if len(record.seq):
               SeqIO.write([record], out_fh, 'fasta')
