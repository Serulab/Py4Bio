from Bio import SeqIO

INPUT_FILE = 'fasta22_out.fas'
OUTPUT_FILE = 'fasta33.fas'

with open(INPUT_FILE) as in_fh:
    with open(OUTPUT_FILE, 'w') as out_fh:
        for record in SeqIO.parse(in_fh,'fasta'):
            # Modify description
            record.description += '[Rattus norvegicus]'
            SeqIO.write([record], out_fh, 'fasta')
