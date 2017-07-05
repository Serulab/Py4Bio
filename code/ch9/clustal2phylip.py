from Bio import AlignIO
fi = open('../../samples/example.aln')
with open('example.phy', 'w') as fo:
    align = AlignIO.read(fi, 'clustal')
    AlignIO.write([align], fo, 'phylip')
