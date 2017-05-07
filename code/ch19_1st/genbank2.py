from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

GB_FILE = '../../samples/NC_006581.gb'
OUT_FILE = 'tg.fasta'
with open(GB_FILE) as gb_fh:
    record = SeqIO.read(gb_fh, 'genbank')
seqs_for_fasta = []
tg = (['cox2'],['atp6'],['atp9'],['cob'])
for feature in record.features:
    if feature.qualifiers.get('gene') in tg and feature.type=='gene':
        # Get the name of the gene
        genename = feature.qualifiers.get('gene')
        # Get the start position
        startpos = feature.location.start.position
        # Get the required slice
        newfrag = record.seq[startpos-1000: startpos]
        # Build a SeqRecord object
        newrec = SeqRecord(newfrag, genename[0] + ' 1000bp upstream',
                           '','')
        seqs_for_fasta.append(newrec)
with open(OUT_FILE,'w') as outf:
    # Write all the sequences as a FASTA file.
    SeqIO.write(seqs_for_fasta, outf, 'fasta')
